// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Principal Component Analysis (PCA) of a count table.
 *
 * With Z the n x p table after normalization (samples in
 * rows), J = I - 11'/n the centering matrix, and D the
 * diagonal matrix of feature standard deviations (or I when
 * not scaling), PCA is the SVD of
 *
 *     A = J * Z * D^-1  =  U * S * V'
 *
 * Scores are U * S, loadings V, and each axis's variance is
 * s^2 / (n - 1), as in stats::prcomp().
 *
 * Only the top k singular triplets are wanted, so A is never
 * formed. Z is a sparse matrix plus a rank-1 term:
 *
 *     Z = S + z0 * 1'
 *
 * where z0[i] is the value that every zero of sample i takes
 * after normalization (log(pseudocount) - mean log for CLR;
 * 0 otherwise) and S holds the non-zeros minus z0. So
 *
 *     Z  * x = S  * x + z0 * sum(x)
 *     Z' * y = S' * y + 1  * (z0' * y)
 *
 * cost O(nnz + n + p) per vector, and J is just subtracting
 * the mean.
 *
 * Precision: a feature with no zeros is stored minus its
 * mean (shift_vec), as T = S - 1 * shift', which keeps it
 * sparse. Z * x then lacks 1 * (shift' * x), and Z' * y
 * lacks shift * (1' * y), but every product here is either
 * centered afterwards (Z * x) or made from a centered y
 * (Z' * y), so both are exactly zero. Without the shift, a
 * feature whose mean is large next to its spread (e.g.
 * values near 1e9) loses most of its digits to
 * cancellation. Block Lanczos (lanczos.c) then finds the top k
 * eigenpairs of the smaller of
 *
 *     A * A'  (n x n; when n - 1 <= p)  ->  U, s^2
 *     A' * A  (p x p; otherwise)        ->  V, s^2
 *
 * and the other factor follows from one more product:
 * V = A' * U / s, or U * S = A * V.
 *
 * Determinism: S * X is computed in fixed blocks of PCA_ROWS
 * samples from a by-sample (CSR) copy of S, and S' * Y in
 * fixed blocks of PCA_ROWS features from a by-feature (CSC)
 * copy, so every output element is summed by one thread in a
 * fixed order. Everything else runs on the main thread.
 *
 * Sign convention: as in pcoa(), each axis is flipped so that
 * its score with the largest absolute value is positive (the
 * first such sample, if tied).
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// Rows (samples or features) per parallel task. Fixed, so
// that summation order never depends on the thread count.
#define PCA_ROWS 64

// Vectors per pass over S, and a macro that repeats F(v)
// for v = 0 to PCA_W - 1.
#define PCA_W 4
#define PCA_REP(F) F(0) F(1) F(2) F(3)


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output rows that each task owns.
//------------------------------------------------------
static int           n_sams;   // Samples (rows of Z).
static int           n_feat;   // Features (columns of Z).

static const int    *r_pos;    // CSR: sample i's non-zeros are
static const int    *r_idx;    //   r_idx/r_val[r_pos[i]..r_pos[i+1]).
static const double *r_val;    //   Values of Z; subtract z0[i] for S.
static const double *z0_vec;   // Zeros' values per sample, or NULL.

static int          *c_pos;    // CSC: feature j's non-zeros are
static int          *c_idx;    //   c_idx/c_val[c_pos[j]..c_pos[j+1]).
static double       *c_val;    //   Values of S, minus shift_vec[j].
static double       *shift_vec; // Per-feature shift (0 if any zeros).

static const double *x_mtx;    // Input, row-major, PCA_W columns.
static double       *y_mtx;    // Output, row-major, PCA_W columns.
static double        x_add[PCA_W]; // Rank-1 term per column.

static double       *mu_vec;   // Per-feature mean of Z.
static double       *ss_vec;   // Per-feature centered sum of squares.
static double       *amax_vec; // Per-feature largest |Z|.
static double        z0_sum;   // sum(z0)
static double        z0_ctr;   // sum((z0 - mean(z0))^2)
static int           z0_const; // All z0 equal (or no z0)?
static double        z0_amax;  // max(abs(z0))


//------------------------------------------------------
// Y = Z * X, for X p x PCA_W and Y n x PCA_W (row-major).
// x_add[v] must hold sum(X[, v]).
//------------------------------------------------------
static void *mult_z(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = n_sams;
  int n_tasks   = (n + PCA_ROWS - 1) / PCA_ROWS;

  #define ZERO(v) s[v] = 0;
  #define DOT(v)  s[v] += a * xj[v];
  #define OUT(v)  yi[v] = s[v] + z0 * x_add[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * PCA_ROWS;
    int i1 = (i0 + PCA_ROWS < n) ? i0 + PCA_ROWS : n;

    for (int i = i0; i < i1; i++) {

      double z0 = z0_vec ? z0_vec[i] : 0;
      double s[PCA_W];
      PCA_REP(ZERO)

      for (int e = r_pos[i]; e < r_pos[i + 1]; e++) {
        int           j  = r_idx[e];
        double        a  = r_val[e] - z0 - shift_vec[j];
        const double *xj = x_mtx + (size_t) j * PCA_W;
        PCA_REP(DOT)
      }

      double *yi = y_mtx + (size_t) i * PCA_W;
      PCA_REP(OUT)
    }
  }

  #undef ZERO
  #undef DOT
  #undef OUT

  return NULL;
}


//------------------------------------------------------
// Y = Z' * X, for X n x PCA_W and Y p x PCA_W (row-major).
// x_add[v] must hold z0' * X[, v].
//------------------------------------------------------
static void *mult_zt(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int p         = n_feat;
  int n_tasks   = (p + PCA_ROWS - 1) / PCA_ROWS;

  #define ZERO(v) s[v] = 0;
  #define DOT(v)  s[v] += a * xi[v];
  #define OUT(v)  yj[v] = s[v] + x_add[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int j0 = task * PCA_ROWS;
    int j1 = (j0 + PCA_ROWS < p) ? j0 + PCA_ROWS : p;

    for (int j = j0; j < j1; j++) {

      double s[PCA_W];
      PCA_REP(ZERO)

      for (int e = c_pos[j]; e < c_pos[j + 1]; e++) {
        double        a  = c_val[e];
        const double *xi = x_mtx + (size_t) c_idx[e] * PCA_W;
        PCA_REP(DOT)
      }

      double *yj = y_mtx + (size_t) j * PCA_W;
      PCA_REP(OUT)
    }
  }

  #undef ZERO
  #undef DOT
  #undef OUT

  return NULL;
}


//------------------------------------------------------
// Each feature's mean, centered sum of squares, and
// largest absolute value, from the CSC copy of S.
//------------------------------------------------------
static void *col_stats(void *arg) {

  int    thread_i  = ((worker_t *) arg)->i;
  int    n_threads = ((worker_t *) arg)->n;
  int    n         = n_sams;
  double zbar      = z0_sum / n;

  for (int j = thread_i; j < n_feat; j += n_threads) {

    int e0 = c_pos[j], e1 = c_pos[j + 1], nz = e1 - e0;

    double cs = 0;
    for (int e = e0; e < e1; e++) cs += c_val[e];
    double m = (cs + z0_sum) / n;

    // Non-zeros: Z = S + z0.
    double ss = 0, amax = 0, ss0 = 0;
    for (int e = e0; e < e1; e++) {
      double z0 = z0_vec ? z0_vec[c_idx[e]] : 0;
      double z  = c_val[e] + z0;
      ss += (z - m) * (z - m);
      if (fabs(z) > amax) amax = fabs(z);
      if (!z0_const) ss0 += (z0 - m) * (z0 - m);
    }

    // Zeros: each sample's z0. Over all samples, sum((z0 -
    // m)^2) = z0_ctr + n * (zbar - m)^2; take away the
    // samples where feature j is non-zero.
    if (nz < n) {
      if (z0_const) {
        double z0 = z0_vec ? z0_vec[0] : 0;
        ss += (n - nz) * (z0 - m) * (z0 - m);
        if (fabs(z0) > amax) amax = fabs(z0);
      }
      else {
        double zero_ss = z0_ctr + n * (zbar - m) * (zbar - m) - ss0;
        if (zero_ss > 0) ss += zero_ss;
        if (z0_amax > amax) amax = z0_amax;
      }
    }

    mu_vec[j]   = m;
    ss_vec[j]   = ss;
    amax_vec[j] = amax;
  }

  return NULL;
}


//------------------------------------------------------
// Row-major helpers for m <= PCA_W columns of length len.
//------------------------------------------------------

// Columns of `V` (column-major, length len) into the
// row-major scratch `x`, each element times w[row] (if
// w), zero-padded to PCA_W columns.
static void pack(const double *V, int len, int m, const double *w, double *x) {
  for (int i = 0; i < len; i++)
    for (int v = 0; v < PCA_W; v++) {
      double a = (v < m) ? V[(size_t) v * len + i] : 0;
      x[(size_t) i * PCA_W + v] = w ? a * w[i] : a;
    }
}

// The reverse: m columns of row-major `y` into `W`, each
// element times w[row] (if w).
static void unpack(const double *y, int len, int m, const double *w, double *W) {
  for (int v = 0; v < m; v++)
    for (int i = 0; i < len; i++) {
      double a = y[(size_t) i * PCA_W + v];
      W[(size_t) v * len + i] = w ? a * w[i] : a;
    }
}

// Subtract each column's mean.
static void center_rows(double *y, int len) {
  for (int v = 0; v < PCA_W; v++) {
    double s = 0;
    for (int i = 0; i < len; i++) s += y[(size_t) i * PCA_W + v];
    s /= len;
    for (int i = 0; i < len; i++) y[(size_t) i * PCA_W + v] -= s;
  }
}


//------------------------------------------------------
// Y = Z * X (n x PCA_W) and Y = Z' * X (p x PCA_W).
//------------------------------------------------------
typedef struct {
  double       *xn, *yn;   // n x PCA_W scratch.
  double       *xp, *yp;   // p x PCA_W scratch.
  const double *w_feat;    // 1 / sd per feature, or NULL.
  const double *w2_feat;   // 1 / var per feature, or NULL.
  int           n_threads;
} pca_ctx_t;

static void run_z(const double *x, double *y, int n_threads) {
  for (int v = 0; v < PCA_W; v++) {
    double s = 0;
    for (int j = 0; j < n_feat; j++) s += x[(size_t) j * PCA_W + v];
    x_add[v] = s;
  }
  x_mtx = x;
  y_mtx = y;
  run_parallel(mult_z, n_threads, n_sams);
}

static void run_zt(const double *x, double *y, int n_threads) {
  for (int v = 0; v < PCA_W; v++) {
    double s = 0;
    if (z0_vec)
      for (int i = 0; i < n_sams; i++) s += z0_vec[i] * x[(size_t) i * PCA_W + v];
    x_add[v] = s;
  }
  x_mtx = x;
  y_mtx = y;
  run_parallel(mult_zt, n_threads, n_feat);
}


//------------------------------------------------------
// lanczos_op_t's. With A = J * Z * D^-1:
//   sample side:  W = A * A' * V = J * Z * D^-2 * Z' * J * V
//   feature side: W = A' * A * V = D^-1 * Z' * J * Z * D^-1 * V
//------------------------------------------------------
static void op_samples(const double *V, double *W, int m, void *ctx) {

  pca_ctx_t *c = (pca_ctx_t *) ctx;
  int        n = n_sams, p = n_feat;

  for (int v0 = 0; v0 < m; v0 += PCA_W) {
    int mc = (m - v0 < PCA_W) ? m - v0 : PCA_W;
    pack(V + (size_t) v0 * n, n, mc, NULL, c->xn);
    center_rows(c->xn, n);
    run_zt(c->xn, c->yp, c->n_threads);
    if (c->w2_feat)
      for (int j = 0; j < p; j++)
        for (int v = 0; v < PCA_W; v++)
          c->yp[(size_t) j * PCA_W + v] *= c->w2_feat[j];
    run_z(c->yp, c->yn, c->n_threads);
    center_rows(c->yn, n);
    unpack(c->yn, n, mc, NULL, W + (size_t) v0 * n);
  }
}

static void op_features(const double *V, double *W, int m, void *ctx) {

  pca_ctx_t *c = (pca_ctx_t *) ctx;
  int        n = n_sams, p = n_feat;

  for (int v0 = 0; v0 < m; v0 += PCA_W) {
    int mc = (m - v0 < PCA_W) ? m - v0 : PCA_W;
    pack(V + (size_t) v0 * p, p, mc, c->w_feat, c->xp);
    run_z(c->xp, c->yn, c->n_threads);
    center_rows(c->yn, n);
    run_zt(c->yn, c->yp, c->n_threads);
    unpack(c->yp, p, mc, c->w_feat, W + (size_t) v0 * p);
  }
}


//======================================================
// R interface. Returns list(scores, loadings, eig, total,
// steps, converged, resid): scores n x k (U * S), loadings
// p x k (V), eig the k values s^2, and total the sum of
// all s^2 (squared Frobenius norm of A). Axes with s^2 <= 0
// have scores and loadings of 0.
//======================================================
SEXP C_pca(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_k, SEXP sexp_scale,
    SEXP sexp_n_threads ) {

  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    k           = asInteger(sexp_k);
  int    scale       = asLogical(sexp_scale);
  int    n_threads   = asInteger(sexp_n_threads);

  init_n_ptrs(10);

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  // Before normalizing: 'binary' would make NA a presence.
  for (int e = 0; e < em->nnz; e++)
    if (!R_FINITE(em->val_vec[e])) {
      free_all();
      error("`counts` must not contain NA, NaN, or infinite values.");
    }

  if (norm) normalize(em, norm, n_threads, pseudocount);

  int n   = em->n_samples;
  int p   = em->n_otus;
  int nnz = em->nnz;

  if (n < 2 || p < 1 || k < 1 || k > n - 1 || k > p || n_threads < 1) {
    free_all();                                       // # nocov
    error("Invalid arguments to C_pca().");           // # nocov
  }

  n_sams = n;
  n_feat = p;
  r_pos  = em->pos_vec;
  r_idx  = em->otu_vec;
  r_val  = em->val_vec;
  z0_vec = em->clr_vec; // Non-NULL only for CLR.


  //--------------------------------------------------
  // CSC copy of S, by counting sort (keeps each
  // feature's samples in increasing order).
  //--------------------------------------------------
  c_pos = (int *)    R_alloc((size_t) p + 1, sizeof(int));
  c_idx = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  c_val = (double *) R_alloc((size_t) nnz + 1, sizeof(double));

  memset(c_pos, 0, ((size_t) p + 1) * sizeof(int));
  for (int e = 0; e < nnz; e++) c_pos[r_idx[e] + 1]++;
  for (int j = 0; j < p; j++)   c_pos[j + 1] += c_pos[j];

  int *fill = (int *) R_alloc((size_t) p + 1, sizeof(int));
  memcpy(fill, c_pos, ((size_t) p + 1) * sizeof(int));

  for (int i = 0; i < n; i++) {
    double z0 = z0_vec ? z0_vec[i] : 0;
    for (int e = r_pos[i]; e < r_pos[i + 1]; e++) {
      int dst    = fill[r_idx[e]]++;
      c_idx[dst] = i;
      c_val[dst] = r_val[e] - z0;
    }
  }


  //--------------------------------------------------
  // Feature means and variances.
  //--------------------------------------------------
  z0_sum = 0; z0_ctr = 0; z0_const = 1; z0_amax = 0;
  if (z0_vec) {
    for (int i = 0; i < n; i++) z0_sum += z0_vec[i];
    double zbar = z0_sum / n;
    for (int i = 0; i < n; i++) {
      z0_ctr += (z0_vec[i] - zbar) * (z0_vec[i] - zbar);
      if (z0_vec[i] != z0_vec[0]) z0_const = 0;
      if (fabs(z0_vec[i]) > z0_amax) z0_amax = fabs(z0_vec[i]);
    }
  }

  mu_vec   = (double *) R_alloc(p, sizeof(double));
  ss_vec   = (double *) R_alloc(p, sizeof(double));
  amax_vec = (double *) R_alloc(p, sizeof(double));
  run_parallel(col_stats, n_threads, p);

  // Features with no zeros: store S minus its column mean.
  shift_vec = (double *) R_alloc(p, sizeof(double));
  for (int j = 0; j < p; j++) {
    double s = 0;
    if (c_pos[j + 1] - c_pos[j] == n) {
      for (int e = c_pos[j]; e < c_pos[j + 1]; e++) s += c_val[e];
      s /= n;
      for (int e = c_pos[j]; e < c_pos[j + 1]; e++) c_val[e] -= s;
    }
    shift_vec[j] = s;
  }

  // Feature weights: 1 / sd (and 1 / var) when scaling.
  // Constant features get 0 and drop out, where prcomp()
  // would stop with an error.
  double *w_feat  = NULL;
  double *w2_feat = NULL;
  if (scale) {
    w_feat  = (double *) R_alloc(p, sizeof(double));
    w2_feat = (double *) R_alloc(p, sizeof(double));
  }

  double total = 0;
  for (int j = 0; j < p; j++) {
    double ss = ss_vec[j];
    if (scale) {
      double tiny = n * (16 * DBL_EPSILON * amax_vec[j]) * (16 * DBL_EPSILON * amax_vec[j]);
      if (ss > tiny) {
        w_feat[j]  = 1 / sqrt(ss / (n - 1));
        w2_feat[j] = 1 / (ss / (n - 1));
        total     += n - 1;
      }
      else {
        w_feat[j]  = 0;
        w2_feat[j] = 0;
      }
    }
    else {
      total += ss;
    }
  }

  if (!R_FINITE(total)) {
    free_all();
    error("`counts` must not contain NA, NaN, or infinite values.");
  }


  //--------------------------------------------------
  // Top k eigenpairs of A * A' or A' * A.
  //--------------------------------------------------
  int side_n = (n - 1 <= p); // Sample side?
  int len    = side_n ? n : p;

  pca_ctx_t ctx = {
    .xn        = (double *) R_alloc((size_t) n * PCA_W, sizeof(double)),
    .yn        = (double *) R_alloc((size_t) n * PCA_W, sizeof(double)),
    .xp        = (double *) R_alloc((size_t) p * PCA_W, sizeof(double)),
    .yp        = (double *) R_alloc((size_t) p * PCA_W, sizeof(double)),
    .w_feat    = w_feat,
    .w2_feat   = w2_feat,
    .n_threads = threads_for(n_threads, ((double) nnz + n + p) * PCA_W, 1e6) };

  SEXP sexp_scores   = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_loadings = PROTECT(allocMatrix(REALSXP, p, k));
  SEXP sexp_eig      = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_res      = PROTECT(allocVector(REALSXP, k));
  double *scores     = REAL(sexp_scores);
  double *loadings   = REAL(sexp_loadings);
  double *eig        = REAL(sexp_eig);

  lanczos_t lz = {
    .n = len, .k = k, .dims = side_n ? n - 1 : p, .deflate = side_n,
    .scale = total, .op = side_n ? op_samples : op_features, .ctx = &ctx,
    .vecs = side_n ? scores : loadings, .vals = eig, .resid = REAL(sexp_res) };

  // No variation at all (e.g., every feature constant): no
  // axes, rather than eigenvectors of rounding noise.
  if (total > 0) {
    lanczos(&lz);
  }
  else {
    memset(lz.vecs, 0, (size_t) len * k * sizeof(double));
    memset(eig,     0, (size_t) k * sizeof(double));
    memset(REAL(sexp_res), 0, (size_t) k * sizeof(double));
    lz.steps     = 0;
    lz.converged = 1;
  }


  //--------------------------------------------------
  // The other factor: loadings = A' * U / s, or
  // scores = A * V.
  //--------------------------------------------------
  for (int v0 = 0; v0 < k; v0 += PCA_W) {

    int mc = (k - v0 < PCA_W) ? k - v0 : PCA_W;

    if (side_n) {
      pack(scores + (size_t) v0 * n, n, mc, NULL, ctx.xn);
      center_rows(ctx.xn, n);
      run_zt(ctx.xn, ctx.yp, n_threads);
      unpack(ctx.yp, p, mc, w_feat, loadings + (size_t) v0 * p);
    }
    else {
      pack(loadings + (size_t) v0 * p, p, mc, w_feat, ctx.xp);
      run_z(ctx.xp, ctx.yn, n_threads);
      center_rows(ctx.yn, n);
      unpack(ctx.yn, n, mc, NULL, scores + (size_t) v0 * n);
    }
  }


  //--------------------------------------------------
  // Scale, and fix each axis's sign.
  //--------------------------------------------------
  for (int j = 0; j < k; j++) {

    double *x = scores   + (size_t) n * j;
    double *l = loadings + (size_t) p * j;
    double  s = (eig[j] > 0) ? sqrt(eig[j]) : 0;

    // Sample side: x = U, l = A' * U = V * s.
    // Feature side: x = A * V = U * s, l = V.
    double xmult = side_n ? s : 1;
    double lmult = side_n ? ((s > 0) ? 1 / s : 0) : ((s > 0) ? 1 : 0);
    if (s == 0) xmult = 0;

    int imax = lz_sign_index(x, n);

    if (x[imax] < 0) { xmult = -xmult; lmult = -lmult; }

    for (int i = 0; i < n; i++) x[i] *= xmult;
    for (int i = 0; i < p; i++) l[i] *= lmult;

    // Features left out by `scale` have no loading. (On the
    // feature side, the starting vectors leave a trace.)
    if (w_feat)
      for (int i = 0; i < p; i++)
        if (w_feat[i] == 0) l[i] = 0;
  }

  free_all();

  const char *names[] = {
    "scores", "loadings", "eig", "total", "steps", "converged", "resid", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_scores);
  SET_VECTOR_ELT(sexp_result, 1, sexp_loadings);
  SET_VECTOR_ELT(sexp_result, 2, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 3, ScalarReal(total));
  SET_VECTOR_ELT(sexp_result, 4, ScalarInteger(lz.steps));
  SET_VECTOR_ELT(sexp_result, 5, ScalarLogical(lz.converged));
  SET_VECTOR_ELT(sexp_result, 6, sexp_res);

  UNPROTECT(5);
  return sexp_result;
}
