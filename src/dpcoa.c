// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Double principal coordinate analysis (DPCoA) of Pavoine, Dufour &
 * Chessel (2004), J. Theor. Biol. 228, 523-537: written from the paper
 * and checked against the outputs of ade4's dpcoa(), whose code was not
 * read (see plans/dpcoa.md).
 *
 * Notation. Sample i has total count n_i, weight w_i = n_i / sum(n), and
 * relative abundances p_i (one per taxon, summing to 1). The taxa have
 * weights pbar = sum_i w_i p_i, and dissimilarities delta (Euclidean).
 * DPCoA places the taxa at points x_k of a Euclidean space, with
 * |x_k - x_l|^2 = delta_kl^2, centred so that sum_k pbar_k x_k = 0, and
 * each sample at y_i = sum_k p_ik x_k. Its axes are the principal axes of
 * the samples, weighted by w: the eigenvectors of Y' W Y, which has the
 * same non-zero eigenvalues as the n x n matrix
 *
 *   G = W^1/2 Y Y' W^1/2 = W^1/2 Q' K Q W^1/2,  Q = P - pbar 1',
 *
 * for any kernel K whose induced distances are delta^2 (K_kk + K_ll -
 * 2 K_kl = delta_kl^2): Q's columns sum to 0, so the parts of K that
 * aren't inner products cancel. Two kernels:
 *
 *   Tree:    delta^2 = the patristic distance, and K_kl = the length of
 *            the path the root shares with k and l. So K q, for any q,
 *            is a clade sum up the tree (for each node, q summed over
 *            the tips below it, times the node's branch length) and then
 *            a path sum down it: O(D), and no D x D matrix.
 *   taxa_dm: K = -1/2 J delta^2 J, a product with the packed distances
 *            (gower.c, as pcoa() uses), after checking that delta is
 *            Euclidean (K has no negative eigenvalue below rounding).
 *
 * The top k eigenpairs of G come from block Lanczos (lanczos.c), with
 * products G V computed from the sparse table: u = W^1/2 v; q = Q u =
 * P u - pbar sum(u) (taxon-major, one taxon per row); r = K q, less its
 * pbar-weighted mean; G v = W^1/2 P' r (sample-major, one sample per
 * row). Basis vectors are kept orthogonal to sqrt(w), which spans G's
 * null space from the centring.
 *
 * Results, with G's unit eigenvectors u_a and eigenvalues lambda_a:
 *   samples:  u_a sqrt(lambda_a) / sqrt(w)  (ade4's $li);
 *   taxa:     (K Q W^1/2 u_a, less its pbar-weighted mean) /
 *             sqrt(lambda_a)  (ade4's $dls: x_k projected on axis a);
 *   rao:      each sample's Rao quadratic entropy, 1/2 sum_kl p_ik p_il
 *             delta_kl^2; with a tree, sum over branches of length L,
 *             with a = p_i summed over the tips below, of L a (1 - a);
 *   trace:    total inertia, sum_i w_i |y_i|^2 (the among-sample part of
 *             Rao's diversity); with a tree, sum over branches of L
 *             times the w-weighted sum of (a - abar)^2; with taxa_dm,
 *             Rao's diversity of pbar less the w-weighted mean of rao.
 *
 * Determinism: every sum runs in a fixed order, and each value is
 * computed by one thread, so results are identical for any number of
 * threads.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON


// Frees new_ecomatrix()'s allocations, then raises an R error.
#define STOP(...) do { free_all(); error(__VA_ARGS__); } while (0)

// Vectors per pass over the table, and a macro that
// repeats F(v) for v = 0 to DP_W - 1.
#define DP_W 4
#define DP_REP(F) F(0) F(1) F(2) F(3)

// Rows (samples or taxa) per parallel task. Fixed, so
// that results never depend on the thread count.
#define DP_ROWS 64


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output rows that each task owns.
//------------------------------------------------------
static int           dp_n;      // Samples.
static int           dp_D;      // Taxa.
static const int    *dp_rpos;   // Sample-major: sample i's entries are
static const int    *dp_ridx;   //   dp_ridx/dp_rval[dp_rpos[i] ..
static const double *dp_rval;   //   dp_rpos[i + 1]), values p_ik.
static const int    *dp_cpos;   // Taxon-major: taxon k's entries are
static const int    *dp_cidx;   //   dp_cidx/dp_cval[dp_cpos[k] ..
static const double *dp_cval;   //   dp_cpos[k + 1]), in sample order.
static const double *dp_x;      // Input, row-major, DP_W columns.
static double       *dp_y;      // Output, row-major, DP_W columns.

// Tree: nodes 0 .. D - 1 are the tips (taxa), D the root,
// then the other internal nodes; edges in preorder.
static int           dp_ne;     // Edges (2D - 2).
static const int    *dp_epar;   // Parent and child of each edge.
static const int    *dp_echd;
static const double *dp_nlen;   // Each node's branch length (root: 0).
static const double *dp_abar;   // Each node's pbar summed over its tips.

// Per-sample results (rao, and its share of the trace).
static double       *dp_rao;
static double       *dp_tr;
static double       *dp_scratch;  // One (2D - 1)-vector per thread.

// taxa_dm: the packed lower triangle of delta.
static const double *dp_dist;


//------------------------------------------------------
// y = P x: one row per taxon, x n x DP_W (row-major).
//------------------------------------------------------
static void *mult_taxa(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n_tasks   = (dp_D + DP_ROWS - 1) / DP_ROWS;

  #define ZERO(v) s[v] = 0;
  #define DOT(v)  s[v] += a * xi[v];
  #define OUT(v)  yk[v] = s[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int k0 = task * DP_ROWS;
    int k1 = (k0 + DP_ROWS < dp_D) ? k0 + DP_ROWS : dp_D;
    for (int k = k0; k < k1; k++) {
      double s[DP_W];
      DP_REP(ZERO)
      for (int e = dp_cpos[k]; e < dp_cpos[k + 1]; e++) {
        double        a  = dp_cval[e];
        const double *xi = dp_x + (size_t) dp_cidx[e] * DP_W;
        DP_REP(DOT)
      }
      double *yk = dp_y + (size_t) k * DP_W;
      DP_REP(OUT)
    }
  }

  #undef ZERO
  #undef DOT
  #undef OUT
  return NULL;
}


//------------------------------------------------------
// y = P' x: one row per sample, x D x DP_W (row-major).
//------------------------------------------------------
static void *mult_samples(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n_tasks   = (dp_n + DP_ROWS - 1) / DP_ROWS;

  #define ZERO(v) s[v] = 0;
  #define DOT(v)  s[v] += a * xk[v];
  #define OUT(v)  yi[v] = s[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i0 = task * DP_ROWS;
    int i1 = (i0 + DP_ROWS < dp_n) ? i0 + DP_ROWS : dp_n;
    for (int i = i0; i < i1; i++) {
      double s[DP_W];
      DP_REP(ZERO)
      for (int e = dp_rpos[i]; e < dp_rpos[i + 1]; e++) {
        double        a  = dp_rval[e];
        const double *xk = dp_x + (size_t) dp_ridx[e] * DP_W;
        DP_REP(DOT)
      }
      double *yi = dp_y + (size_t) i * DP_W;
      DP_REP(OUT)
    }
  }

  #undef ZERO
  #undef DOT
  #undef OUT
  return NULL;
}


//------------------------------------------------------
// Rao's quadratic entropy of each sample, and its share
// of the trace, with a tree: a = p_i summed over the
// tips below each node, then sums over the branches.
//------------------------------------------------------
static void *rao_tree(void *arg) {

  int     thread_i  = ((worker_t *) arg)->i;
  int     n_threads = ((worker_t *) arg)->n;
  int     D         = dp_D;
  int     n_nodes   = 2 * D - 1;
  int     n_tasks   = (dp_n + DP_ROWS - 1) / DP_ROWS;
  double *a         = dp_scratch + (size_t) thread_i * n_nodes;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i0 = task * DP_ROWS;
    int i1 = (i0 + DP_ROWS < dp_n) ? i0 + DP_ROWS : dp_n;
    for (int i = i0; i < i1; i++) {

      memset(a, 0, (size_t) n_nodes * sizeof(double));
      for (int e = dp_rpos[i]; e < dp_rpos[i + 1]; e++)
        a[dp_ridx[e]] += dp_rval[e];
      for (int e = dp_ne - 1; e >= 0; e--)
        a[dp_epar[e]] += a[dp_echd[e]];

      double rao = 0, tr = 0;
      for (int e = 0; e < dp_ne; e++) {
        int    c = dp_echd[e];
        double d = a[c] - dp_abar[c];
        rao += dp_nlen[c] * a[c] * (1 - a[c]);
        tr  += dp_nlen[c] * d * d;
      }
      dp_rao[i] = rao;
      dp_tr[i]  = tr;
    }
  }

  return NULL;
}


//------------------------------------------------------
// Rao's quadratic entropy of each sample with taxa_dm:
// sum over pairs k < l of p_ik p_il delta_kl^2, over the
// sample's non-zero taxa.
//------------------------------------------------------
static void *rao_dist(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D         = dp_D;
  int n_tasks   = (dp_n + DP_ROWS - 1) / DP_ROWS;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i0 = task * DP_ROWS;
    int i1 = (i0 + DP_ROWS < dp_n) ? i0 + DP_ROWS : dp_n;
    for (int i = i0; i < i1; i++) {
      double rao = 0;
      for (int e = dp_rpos[i]; e < dp_rpos[i + 1]; e++) {
        int    k  = dp_ridx[e];
        double pk = dp_rval[e], s = 0;
        for (int f = dp_rpos[i]; f < e; f++) {
          int l = dp_ridx[f];
          if (l == k) continue;
          int    lo = (k < l) ? k : l, hi = (k < l) ? l : k;
          double d  = dp_dist[dist_col_start(lo, D) + (size_t) (hi - lo - 1)];
          s += dp_rval[f] * d * d;
        }
        rao += pk * s;
      }
      dp_rao[i] = rao;
    }
  }

  return NULL;
}


//------------------------------------------------------
// For taxa_dm: column k's part of Rao's diversity of
// pbar, pbar_k * sum_{l > k} pbar_l delta_kl^2.
//------------------------------------------------------
static const double *dp_pbar;
static double       *dp_colpart;

static void *rao_pooled(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int D         = dp_D;
  for (int k = thread_i; k < D - 1; k += n_threads) {
    const double *col = dp_dist + dist_col_start(k, D);
    double s = 0;
    for (int l = k + 1; l < D; l++) {
      double d = col[l - k - 1];
      s += dp_pbar[l] * d * d;
    }
    dp_colpart[k] = dp_pbar[k] * s;
  }
  return NULL;
}


//------------------------------------------------------
// The products, for lanczos_op_t and the taxa scores.
//------------------------------------------------------
typedef struct {
  int           n, D;
  const double *sqw;       // sqrt(w), per sample.
  const double *pbar;      // Taxon weights.
  int           tree;      // 1: tree kernel; 0: gower.
  gower_t      *g;         // For taxa_dm.
  double       *xn, *yn;   // n x DP_W scratch (row-major).
  double       *xd, *yd;   // D x DP_W scratch.
  double       *node;      // (2D - 1) x DP_W scratch (tree).
  double       *qc, *rc;   // D x DP_W scratch, column-major (gower).
  int           n_threads;
  int           nt_n, nt_D; // Threads worth using per pass.
} dp_ctx_t;


// r = K q for the DP_W columns of ctx->yd (row-major),
// less each column's pbar-weighted mean, into ctx->xd.
static void kernel_apply(dp_ctx_t *c) {

  int D = c->D;
  double *q = c->yd, *r = c->xd;

  if (c->tree) {

    int     nn  = 2 * D - 1;
    double *acc = c->node;

    // Clade sums of q, children before parents.
    memset(acc, 0, (size_t) nn * DP_W * sizeof(double));
    memcpy(acc, q, (size_t) D * DP_W * sizeof(double));
    for (int e = dp_ne - 1; e >= 0; e--) {
      double *ap = acc + (size_t) dp_epar[e] * DP_W;
      double *ac = acc + (size_t) dp_echd[e] * DP_W;
      for (int v = 0; v < DP_W; v++) ap[v] += ac[v];
    }

    // Times branch lengths, then path sums from the root
    // (in place: each node's clade sum is read before its
    // path sum is written, parents first).
    double *root = acc + (size_t) D * DP_W;
    for (int v = 0; v < DP_W; v++) root[v] = 0;
    for (int e = 0; e < dp_ne; e++) {
      int     ch = dp_echd[e];
      double  L  = dp_nlen[ch];
      double *ap = acc + (size_t) dp_epar[e] * DP_W;
      double *ac = acc + (size_t) ch * DP_W;
      for (int v = 0; v < DP_W; v++) ac[v] = ap[v] + L * ac[v];
    }
    memcpy(r, acc, (size_t) D * DP_W * sizeof(double));
  }

  else {

    // Column-major for gower_mult(), centred (q sums to 0
    // already, up to rounding).
    for (int v = 0; v < DP_W; v++)
      for (int k = 0; k < D; k++)
        c->qc[(size_t) v * D + k] = q[(size_t) k * DP_W + v];
    for (int v = 0; v < DP_W; v++) lz_center(c->qc + (size_t) v * D, D);
    gower_mult(c->qc, c->rc, DP_W, c->g);
    for (int v = 0; v < DP_W; v++)
      for (int k = 0; k < D; k++)
        r[(size_t) k * DP_W + v] = c->rc[(size_t) v * D + k];
  }

  for (int v = 0; v < DP_W; v++) {
    double s = 0;
    for (int k = 0; k < D; k++) s += c->pbar[k] * r[(size_t) k * DP_W + v];
    for (int k = 0; k < D; k++) r[(size_t) k * DP_W + v] -= s;
  }
}


// Taxa side for up to DP_W columns of V (n x ., column-
// major, from column v0): ctx->xd = K Q W^1/2 V, less
// each column's pbar-weighted mean (row-major).
static void taxa_side(dp_ctx_t *c, const double *V, int v0, int mc) {

  int n = c->n, D = c->D;

  double S[DP_W];
  for (int v = 0; v < DP_W; v++) {
    double s = 0;
    for (int i = 0; i < n; i++) {
      double a = (v < mc) ? c->sqw[i] * V[(size_t) (v0 + v) * n + i] : 0;
      c->xn[(size_t) i * DP_W + v] = a;
      s += a;
    }
    S[v] = s;
  }

  dp_x = c->xn;
  dp_y = c->yd;
  run_parallel(mult_taxa, c->nt_D, D);

  for (int k = 0; k < D; k++)
    for (int v = 0; v < DP_W; v++)
      c->yd[(size_t) k * DP_W + v] -= c->pbar[k] * S[v];

  kernel_apply(c);
}


// W = G V: a lanczos_op_t.
static void dp_op(const double *V, double *W, int m, void *ctx) {

  dp_ctx_t *c = (dp_ctx_t *) ctx;
  int       n = c->n;

  for (int v0 = 0; v0 < m; v0 += DP_W) {
    int mc = (m - v0 < DP_W) ? m - v0 : DP_W;
    taxa_side(c, V, v0, mc);
    dp_x = c->xd;
    dp_y = c->yn;
    run_parallel(mult_samples, c->nt_n, n);
    for (int v = 0; v < mc; v++)
      for (int i = 0; i < n; i++)
        W[(size_t) (v0 + v) * n + i] = c->sqw[i] * c->yn[(size_t) i * DP_W + v];
  }
}


// -K V, for K's smallest eigenvalue (taxa_dm's check).
static void neg_gower(const double *V, double *W, int m, void *ctx) {
  gower_mult(V, W, m, ctx);
  size_t len = (size_t) m * ((gower_t *) ctx)->n;
  for (size_t i = 0; i < len; i++) W[i] = -W[i];
}



//======================================================
// R interface.
//
// sexp_counts, sexp_margin: the count table (taxa in the
//   order of the tree's tip ids / taxa_dm's rows).
// sexp_epar, sexp_echd, sexp_elen: the tree's edges in
//   preorder, from C_philr_tree() (0-based ids: tips are
//   the taxa, the root is D), or NULL.
// sexp_dist: taxa_dm's packed lower triangle, or NULL.
// sexp_k: axes wanted (0 to min(n, D) - 1).
//
// Returns list(points, taxa, eig, rao, trace, steps,
// converged, empty, neg). If some samples are empty,
// only `empty` (their 1-based indices) is set. If taxa_dm
// is not Euclidean, only `neg` (K's smallest eigenvalue,
// relative to trace(K)) is set.
//======================================================
SEXP C_dpcoa(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_epar, SEXP sexp_echd,
    SEXP sexp_elen, SEXP sexp_dist, SEXP sexp_k, SEXP sexp_n_threads ) {

  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_n_threads);
  int tree      = !isNull(sexp_epar);

  const char *names[] = {
    "points", "taxa", "eig", "rao", "trace", "steps", "converged",
    "empty", "neg", "" };
  SEXP result = PROTECT(mkNamed(VECSXP, names));

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples, D = em->n_otus, nnz = em->nnz;
  const int    *pos = em->pos_vec;
  const int    *idx = em->otu_vec;
  const double *val = em->val_vec;

  if (k < 0 || k > n - 1 || k > D - 1 || n_threads < 1)
    STOP("Invalid arguments to C_dpcoa()."); // # nocov
  if (tree && (XLENGTH(sexp_epar) != 2 * (R_xlen_t) D - 2 ||
               XLENGTH(sexp_echd) != XLENGTH(sexp_epar) ||
               XLENGTH(sexp_elen) != XLENGTH(sexp_epar)))
    STOP("Invalid arguments to C_dpcoa()."); // # nocov
  if (!tree && XLENGTH(sexp_dist) != (R_xlen_t) D * (D - 1) / 2)
    STOP("Invalid arguments to C_dpcoa()."); // # nocov

  for (int e = 0; e < nnz; e++) {
    double x = val[e];
    if (!R_FINITE(x)) STOP("`counts` must not contain NA, NaN, or infinite values.");
    if (x < 0)        STOP("`counts` must not contain negative values.");
  }


  //--------------------------------------------------
  // Weights and relative abundances. Empty samples go
  // back to R for the message.
  //--------------------------------------------------
  double *tot = (double *) R_alloc(n, sizeof(double));
  double  T   = 0;
  int     n_empty = 0;
  for (int i = 0; i < n; i++) {
    double s = 0;
    for (int e = pos[i]; e < pos[i + 1]; e++) s += val[e];
    tot[i] = s;
    T     += s;
    if (!(s > 0)) n_empty++;
  }

  if (n_empty) {
    SEXP sexp_empty = PROTECT(allocVector(INTSXP, n_empty));
    for (int i = 0, j = 0; i < n; i++)
      if (!(tot[i] > 0)) INTEGER(sexp_empty)[j++] = i + 1;
    SET_VECTOR_ELT(result, 7, sexp_empty);
    free_all();
    UNPROTECT(2);
    return result;
  }

  // Non-zero entries only, as p_ik.
  int    *rpos = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  int    *ridx = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *rval = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  double *sqw  = (double *) R_alloc(n, sizeof(double));
  double *pbar = (double *) R_alloc(D, sizeof(double));
  int    *cpos = (int *)    R_alloc((size_t) D + 1, sizeof(int));
  memset(cpos, 0, ((size_t) D + 1) * sizeof(int));
  memset(pbar, 0, (size_t) D * sizeof(double));

  int m_nz = 0;
  rpos[0] = 0;
  for (int i = 0; i < n; i++) {
    for (int e = pos[i]; e < pos[i + 1]; e++) {
      if (val[e] == 0) continue;
      ridx[m_nz] = idx[e];
      rval[m_nz] = val[e] / tot[i];
      pbar[idx[e]] += val[e];
      cpos[idx[e] + 1]++;
      m_nz++;
    }
    rpos[i + 1] = m_nz;
    sqw[i] = sqrt(tot[i] / T);
  }

  // Taxon-major copy, in sample order; pbar summed over
  // samples in order.
  for (int j = 0; j < D; j++) cpos[j + 1] += cpos[j];
  int    *cidx = (int *)    R_alloc((size_t) m_nz + 1, sizeof(int));
  double *cval = (double *) R_alloc((size_t) m_nz + 1, sizeof(double));
  {
    int *fill = (int *) R_alloc((size_t) D + 1, sizeof(int));
    memcpy(fill, cpos, (size_t) D * sizeof(int));
    for (int i = 0; i < n; i++)
      for (int e = rpos[i]; e < rpos[i + 1]; e++) {
        int f = fill[ridx[e]]++;
        cidx[f] = i;
        cval[f] = rval[e];
      }
  }
  free_all(); // The ecomatrix isn't needed any more.

  // Taxon weights: column totals (summed in sample order
  // above) over the grand total. The basis that Lanczos
  // keeps vectors orthogonal to: sqrt(w), normalized.
  for (int j = 0; j < D; j++) pbar[j] /= T;
  double *f1 = (double *) R_alloc(n, sizeof(double));
  {
    double s = 0;
    for (int i = 0; i < n; i++) s += sqw[i] * sqw[i];
    s = sqrt(s);
    for (int i = 0; i < n; i++) f1[i] = sqw[i] / s;
  }

  dp_n = n; dp_D = D;
  dp_rpos = rpos; dp_ridx = ridx; dp_rval = rval;
  dp_cpos = cpos; dp_cidx = cidx; dp_cval = cval;

  int nt_n = threads_for(n_threads, (double) m_nz * DP_W + n, 2e5);
  int nt_D = threads_for(n_threads, (double) m_nz * DP_W + D, 2e5);


  //--------------------------------------------------
  // The kernel: tree, or taxa_dm (checked Euclidean).
  //--------------------------------------------------
  double *rao   = (double *) R_alloc(n, sizeof(double));
  double  trace = 0, pooled = 0;
  gower_t g;
  int     steps = 0, converged = 1;

  if (tree) {

    int nn = 2 * D - 1, ne = 2 * D - 2;
    const int    *epar = INTEGER(sexp_epar);
    const int    *echd = INTEGER(sexp_echd);
    const double *elen = REAL(sexp_elen);

    double *nlen = (double *) R_alloc(nn, sizeof(double));
    double *abar = (double *) R_alloc(nn, sizeof(double));
    memset(nlen, 0, (size_t) nn * sizeof(double));
    memset(abar, 0, (size_t) nn * sizeof(double));
    for (int e = 0; e < ne; e++) {
      if (epar[e] < D || epar[e] >= nn || echd[e] < 0 || echd[e] >= nn)
        error("Invalid arguments to C_dpcoa()."); // # nocov
      nlen[echd[e]] = elen[e];
    }
    memcpy(abar, pbar, (size_t) D * sizeof(double));
    for (int e = ne - 1; e >= 0; e--) abar[epar[e]] += abar[echd[e]];

    dp_ne = ne; dp_epar = epar; dp_echd = echd; dp_nlen = nlen; dp_abar = abar;

    double *tr = (double *) R_alloc(n, sizeof(double));
    int     nt = threads_for(n_threads, (double) n * nn * 3, 2e5);
    dp_rao     = rao;
    dp_tr      = tr;
    dp_scratch = (double *) R_alloc((size_t) nn * nt, sizeof(double));
    run_parallel(rao_tree, nt, n);

    for (int i = 0; i < n; i++) trace += (tot[i] / T) * tr[i];

    for (int e = 0; e < ne; e++) {
      double a = abar[echd[e]];
      pooled += nlen[echd[e]] * a * (1 - a);
    }
  }

  else {

    const double *dist = REAL(sexp_dist);
    double trace_k = gower_init(&g, dist, D, 16, n_threads);

    // K's smallest eigenvalue, as -K's largest; stop as
    // soon as one is proven below -sqrt(eps) trace(K).
    if (trace_k > 0 && D > 2) {
      double  val = 0, res = 0;
      double *vec = (double *) R_alloc(D, sizeof(double));
      double  lim = sqrt(DBL_EPSILON) * trace_k;
      lanczos_t lz = {
        .n = D, .k = 1, .dims = D - 1, .deflate = 1, .scale = trace_k,
        .op = neg_gower, .ctx = &g, .n_threads = n_threads,
        .vecs = vec, .vals = &val, .resid = &res,
        .early_stop = 1, .stop_above = lim, .stop_tol = INFINITY };
      lanczos(&lz);
      if (val > lim) {
        SET_VECTOR_ELT(result, 8, ScalarReal(-val / trace_k));
        UNPROTECT(1);
        return result;
      }
    }

    dp_dist = dist;
    dp_rao  = rao;
    run_parallel(rao_dist, threads_for(n_threads, (double) m_nz * 20, 2e5), n);

    dp_pbar    = pbar;
    dp_colpart = (double *) R_alloc(D, sizeof(double));
    dp_colpart[D - 1] = 0;
    run_parallel(rao_pooled, threads_for(n_threads, 0.5 * D * (D - 1.0), 1e6), D);

    double within = 0;
    for (int j = 0; j < D; j++) pooled += dp_colpart[j];
    for (int i = 0; i < n; i++) within += (tot[i] / T) * rao[i];
    trace = pooled - within;
  }

  // No variation among the samples (identical
  // compositions): an among-sample diversity at
  // rounding-error level, relative to the total.
  if (!(trace > 1e-12 * pooled)) trace = 0;


  //--------------------------------------------------
  // The top k eigenpairs of G.
  //--------------------------------------------------
  dp_ctx_t ctx = {
    .n = n, .D = D, .sqw = sqw, .pbar = pbar, .tree = tree, .g = &g,
    .n_threads = n_threads, .nt_n = nt_n, .nt_D = nt_D };
  ctx.xn = (double *) R_alloc((size_t) n * DP_W, sizeof(double));
  ctx.yn = (double *) R_alloc((size_t) n * DP_W, sizeof(double));
  ctx.xd = (double *) R_alloc((size_t) D * DP_W, sizeof(double));
  ctx.yd = (double *) R_alloc((size_t) D * DP_W, sizeof(double));
  if (tree) {
    ctx.node = (double *) R_alloc((size_t) (2 * D - 1) * DP_W, sizeof(double));
  }
  else {
    ctx.qc = (double *) R_alloc((size_t) D * DP_W, sizeof(double));
    ctx.rc = (double *) R_alloc((size_t) D * DP_W, sizeof(double));
  }

  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_taxa   = PROTECT(allocMatrix(REALSXP, D, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_rao    = PROTECT(allocVector(REALSXP, n));
  double *points = REAL(sexp_points);
  double *taxa   = REAL(sexp_taxa);
  double *eig    = REAL(sexp_eig);
  memcpy(REAL(sexp_rao), rao, (size_t) n * sizeof(double));

  if (k > 0 && trace > 0) {
    double *res = (double *) R_alloc(k, sizeof(double));
    lanczos_t lz = {
      .n = n, .k = k, .dims = n - 1, .deflate = 0, .basis = f1, .n_basis = 1,
      .scale = trace, .op = dp_op, .ctx = &ctx, .n_threads = n_threads,
      .vecs = points, .vals = eig, .resid = res };
    lanczos(&lz);
    steps     = lz.steps;
    converged = lz.converged;
  }
  else {
    memset(points, 0, (size_t) n * k * sizeof(double));
    memset(eig,    0, (size_t) k * sizeof(double));
  }


  //--------------------------------------------------
  // Taxa scores; scale and fix each axis's sign.
  //--------------------------------------------------
  for (int v0 = 0; v0 < k; v0 += DP_W) {
    int mc = (k - v0 < DP_W) ? k - v0 : DP_W;
    taxa_side(&ctx, points, v0, mc);
    for (int v = 0; v < mc; v++)
      for (int j = 0; j < D; j++)
        taxa[(size_t) (v0 + v) * D + j] = ctx.xd[(size_t) j * DP_W + v];
  }

  for (int a = 0; a < k; a++) {

    double *x = points + (size_t) n * a;
    double *t = taxa   + (size_t) D * a;
    double  s = (eig[a] > 0) ? sqrt(eig[a]) : 0;

    for (int i = 0; i < n; i++) x[i] = (s > 0) ? x[i] * s / sqw[i] : 0;

    int    imax = lz_sign_index(x, n);
    double tm   = (s > 0) ? 1 / s : 0;
    if (x[imax] < 0) {
      tm = -tm;
      for (int i = 0; i < n; i++) x[i] = -x[i];
    }
    for (int j = 0; j < D; j++) t[j] *= tm;
  }

  SET_VECTOR_ELT(result, 0, sexp_points);
  SET_VECTOR_ELT(result, 1, sexp_taxa);
  SET_VECTOR_ELT(result, 2, sexp_eig);
  SET_VECTOR_ELT(result, 3, sexp_rao);
  SET_VECTOR_ELT(result, 4, ScalarReal(trace));
  SET_VECTOR_ELT(result, 5, ScalarInteger(steps));
  SET_VECTOR_ELT(result, 6, ScalarLogical(converged));

  UNPROTECT(5);
  return result;
}
