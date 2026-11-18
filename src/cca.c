// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * The chi-square table of correspondence analysis, for
 * cca() (and the planned ca()), as in ter Braak (1986) and
 * Legendre & Legendre (2012, sections 9.2 and 11.2).
 *
 * With X the n x p table of non-negative abundances, N its
 * grand total, P = X / N, and row and column weights
 * r = P 1 and c = P' 1, CA analyzes the matrix of
 * standardised residuals
 *
 *     Qbar = Dr^-1/2 (P - r c') Dc^-1/2
 *          = A - sqrt(r) sqrt(c)',   A = Dr^-1/2 P Dc^-1/2
 *
 * A is as sparse as X, so products with Qbar are a pass
 * over X's non-zeros plus a rank-one correction; Qbar is
 * never formed. Every column of Qbar is orthogonal to
 * sqrt(r), which takes the place of the ones vector of the
 * unweighted methods. Its squared Frobenius norm, the total
 * inertia (Pearson's chi-square / N), is sum(A^2) - 1.
 *
 * Features with no observations have c = 0: their columns
 * of A and Qbar are 0, so they drop out (loadings 0).
 * Samples with no observations are an error (r = 0).
 *
 * cca_perm_inertia() gives the constrained inertia of
 * Qbar on a design whose rows are permuted, for vegan's
 * permutation-based adjusted R-squared (RsquareAdj() on a
 * cca result, after Peres-Neto et al. 2006);
 * cca_perm_inertia_fixed() the same for partial models,
 * permuting Qbar's rows against a fixed design, as
 * RsquareAdj() does for them.
 *
 * Determinism: each output element is summed by one thread
 * in a fixed order (features for A' V, samples for A V,
 * permutations for the inertia), so results are identical
 * for any number of threads.
 */

#include "ecodive.h"

#define CCA_ROWS 64   // Rows per parallel task.
#define CCA_QR_TOL 1e-7 // As DBRDA_QR_TOL in dbrda.c.


static int           cc_n, cc_p, cc_nnz;
static const int    *r_pos, *r_idx;  // CSR (samples)
static double       *r_val;          // ... values of A
static int          *c_pos, *c_idx;  // CSC (features)
static double       *c_val;          // ... values of A
static double       *sr, *sc;        // sqrt(r), sqrt(c)
static const double *raw_val;        // P (CSR), for permutations

static const double *op_x;  // n x m input
static double       *op_t;  // p x m
static double       *op_y;  // n x m output
static const double *op_sx; // m: sqrt(r)' X  (for T)
static const double *op_st; // m: sqrt(c)' T  (for Y)
static int           op_m;


// T = A' X - sqrt(c) (sqrt(r)' X): one feature per row.
static void *mult_t(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = cc_n, p = cc_p, m = op_m;
  int n_tasks = (p + CCA_ROWS - 1) / CCA_ROWS;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int j1 = (task + 1) * CCA_ROWS;
    if (j1 > p) j1 = p;
    for (int j = task * CCA_ROWS; j < j1; j++)
      for (int v = 0; v < m; v++) {
        const double *x = op_x + (size_t) v * n;
        double        s = 0;
        for (int e = c_pos[j]; e < c_pos[j + 1]; e++) s += c_val[e] * x[c_idx[e]];
        op_t[(size_t) v * p + j] = s - sc[j] * op_sx[v];
      }
  }
  return NULL;
}


// Y = A T - sqrt(r) (sqrt(c)' T): one sample per row.
static void *mult_y(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = cc_n, p = cc_p, m = op_m;
  int n_tasks = (n + CCA_ROWS - 1) / CCA_ROWS;

  for (int task = thread_i; task < n_tasks; task += n_threads) {
    int i1 = (task + 1) * CCA_ROWS;
    if (i1 > n) i1 = n;
    for (int i = task * CCA_ROWS; i < i1; i++)
      for (int v = 0; v < m; v++) {
        const double *t = op_t + (size_t) v * p;
        double        s = 0;
        for (int e = r_pos[i]; e < r_pos[i + 1]; e++) s += r_val[e] * t[r_idx[e]];
        op_y[(size_t) v * n + i] = s - sr[i] * op_st[v];
      }
  }
  return NULL;
}


static int     cc_threads;
static double *buf_sx, *buf_st, *buf_t; // Scratch, for buf_m vectors.
static int     buf_m;

static void grow_buffers(int m) {
  if (m <= buf_m) return;
  buf_sx = (double *) R_alloc(m, sizeof(double));
  buf_st = (double *) R_alloc(m, sizeof(double));
  buf_t  = (double *) R_alloc((size_t) cc_p * m, sizeof(double));
  buf_m  = m;
}

// L = Qbar' X (p x m), for X n x m.
static void cca_t(const double *X, double *L, int m) {

  int     n  = cc_n;
  grow_buffers(m);
  double *sx = buf_sx;
  for (int v = 0; v < m; v++) {
    double s = 0;
    for (int i = 0; i < n; i++) s += sr[i] * X[(size_t) v * n + i];
    sx[v] = s;
  }
  op_x = X; op_t = L; op_sx = sx; op_m = m;
  run_parallel(mult_t, threads_for(cc_threads, (double) (cc_nnz + cc_p) * m, 1e6), cc_p);
}


//======================================================
// Reads the table, checks it, and sets up the products.
// Returns the total inertia, or -1 if some samples have
// no observations: then empty[i] is 1 for each of them
// (empty must hold n ints; *n_out is set either way).
//======================================================
double cca_table_init(
    SEXP sexp_counts, SEXP sexp_margin, int n_threads, int *n_out,
    int *p_out, int **empty ) {

  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples, p = em->n_otus, nnz = em->nnz;

  for (int e = 0; e < nnz; e++) {
    double x = em->val_vec[e];
    if (!R_FINITE(x)) {
      free_all();
      error("`counts` must not contain NA, NaN, or infinite values.");
    }
    if (x < 0) {
      free_all();
      error("`counts` must not contain negative values.");
    }
  }

  *n_out = n;
  *p_out = p;

  // Row sums, in order; the grand total from them.
  double *rs = (double *) R_alloc(n + 1, sizeof(double));
  double *cs = (double *) R_alloc(p + 1, sizeof(double));
  memset(cs, 0, (size_t) (p + 1) * sizeof(double));

  int    n_empty = 0;
  double total_n = 0;
  *empty = (int *) R_alloc(n + 1, sizeof(int));
  for (int i = 0; i < n; i++) {
    double s = 0;
    for (int e = em->pos_vec[i]; e < em->pos_vec[i + 1]; e++) {
      s += em->val_vec[e];
      cs[em->otu_vec[e]] += em->val_vec[e];
    }
    rs[i] = s;
    total_n += s;
    (*empty)[i] = !(s > 0);
    n_empty += (*empty)[i];
  }
  if (n_empty) return -1;

  cc_n = n; cc_p = p; cc_nnz = nnz; cc_threads = n_threads;
  buf_m = 0;

  sr = (double *) R_alloc(n + 1, sizeof(double));
  sc = (double *) R_alloc(p + 1, sizeof(double));
  for (int i = 0; i < n; i++) sr[i] = sqrt(rs[i] / total_n);
  for (int j = 0; j < p; j++) sc[j] = sqrt(cs[j] / total_n);

  // CSR copy with A's values, and P for permutations.
  int    *rp  = (int *)    R_alloc((size_t) n + 1,   sizeof(int));
  int    *ri  = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  double *rv  = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  double *pv  = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  double  ss  = 0;
  int     dst = 0;
  rp[0] = 0;
  for (int i = 0; i < n; i++) {
    for (int e = em->pos_vec[i]; e < em->pos_vec[i + 1]; e++) {
      double x = em->val_vec[e];
      if (x == 0) continue;
      int j   = em->otu_vec[e];
      ri[dst] = j;
      pv[dst] = x / total_n;
      rv[dst] = pv[dst] / (sr[i] * sc[j]);
      ss     += rv[dst] * rv[dst];
      dst++;
    }
    rp[i + 1] = dst;
  }
  nnz = cc_nnz = dst;
  r_pos = rp; r_idx = ri; r_val = rv; raw_val = pv;

  // CSC copy, by counting sort.
  c_pos = (int *)    R_alloc((size_t) p + 1,   sizeof(int));
  c_idx = (int *)    R_alloc((size_t) nnz + 1, sizeof(int));
  c_val = (double *) R_alloc((size_t) nnz + 1, sizeof(double));
  memset(c_pos, 0, ((size_t) p + 1) * sizeof(int));
  for (int e = 0; e < nnz; e++) c_pos[ri[e] + 1]++;
  for (int j = 0; j < p; j++)   c_pos[j + 1] += c_pos[j];
  int *fill = (int *) R_alloc((size_t) p + 1, sizeof(int));
  memcpy(fill, c_pos, ((size_t) p + 1) * sizeof(int));
  for (int i = 0; i < n; i++)
    for (int e = rp[i]; e < rp[i + 1]; e++) {
      int d = fill[ri[e]]++;
      c_idx[d] = i;
      c_val[d] = rv[e];
    }

  free_all(); // The ecomatrix: everything kept is R_alloc()'d.

  double total = ss - 1;
  if (total < 0) total = 0;
  return total;
}


const double *cca_sqrt_r(void) { return sr; }
const double *cca_sqrt_c(void) { return sc; }


// W = Qbar Qbar' V (n x m each): a lanczos_op_t.
void cca_gram_mult(const double *V, double *W, int m, void *ctx) {

  (void) ctx;
  int p = cc_p;
  grow_buffers(m);
  double *T = buf_t, *st = buf_st;

  cca_t(V, T, m);

  for (int v = 0; v < m; v++) {
    double s = 0;
    for (int j = 0; j < p; j++) s += sc[j] * T[(size_t) v * p + j];
    st[v] = s;
  }

  op_t = T; op_y = W; op_st = st; op_m = m;
  run_parallel(mult_y, threads_for(cc_threads, (double) (cc_nnz + cc_n) * m, 1e6), cc_n);
}


// L = Dc^-1/2 Qbar' U (p x m): vegan's species scores in
// scaling 2 for unit eigenvectors U. Empty features get 0.
void cca_cross(const double *U, int m, double *L) {
  int p = cc_p;
  cca_t(U, L, m);
  for (int v = 0; v < m; v++)
    for (int j = 0; j < p; j++) {
      double *l = L + (size_t) v * p + j;
      *l = (sc[j] > 0) ? *l / sc[j] : 0;
    }
}



//======================================================
// Permuted constrained inertia. For permutation k, sample
// p[a] takes design row a (vegan's community rows Y[p, ]
// against fixed constraints, with the weights going with
// the community). The design's columns are centred with
// the weights r, scaled by sqrt(r), orthonormalized
// (dropping aliased ones as cord_design() does), and the
// inertia is ||Qbar' Q||^2 = sum_j (sum_i P_ij q_i /
// sqrt(r_i))^2 / c_j, since sqrt(r)' q = 0.
//======================================================
static const double  *pi_x;     // n x px design (constraints)
static int            pi_px;
static perm_source_t *pi_ps;
static int            pi_nperm;
static double        *pi_out;   // n_perm
static double        *pi_buf;   // per thread: n * px + p + n ints' worth

static void *perm_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = cc_n, p = cc_p, px = pi_px;

  size_t  per = (size_t) n * px + p + n;
  double *Q   = pi_buf + per * thread_i;
  double *t   = Q + (size_t) n * px;
  int    *pm  = (int *) (t + p);  // n ints fit in n doubles

  for (int k = thread_i; k < pi_nperm; k += n_threads) {

    perm_source_get(pi_ps, k + 1, pm);

    int    rq  = 0;
    double sum = 0;

    for (int c = 0; c < px; c++) {

      double       *q  = Q + (size_t) rq * n;
      const double *xc = pi_x + (size_t) c * n;

      for (int a = 0; a < n; a++) q[pm[a]] = xc[a];

      double mu = 0;
      for (int i = 0; i < n; i++) mu += sr[i] * sr[i] * q[i];
      for (int i = 0; i < n; i++) q[i] = sr[i] * (q[i] - mu);

      double norm0 = 0;
      for (int i = 0; i < n; i++) norm0 += q[i] * q[i];
      norm0 = sqrt(norm0);

      for (int pass = 0; pass < 2; pass++)
        for (int b = 0; b < rq; b++) {
          const double *qb = Q + (size_t) b * n;
          double        d  = 0;
          for (int i = 0; i < n; i++) d += qb[i] * q[i];
          for (int i = 0; i < n; i++) q[i] -= d * qb[i];
        }

      double norm = 0;
      for (int i = 0; i < n; i++) norm += q[i] * q[i];
      norm = sqrt(norm);

      if (!(norm0 > 0) || norm < CCA_QR_TOL * norm0) continue;
      for (int i = 0; i < n; i++) q[i] /= norm;

      // t = P' (q / sqrt(r)), then the inertia.
      memset(t, 0, (size_t) p * sizeof(double));
      for (int i = 0; i < n; i++) {
        double y = q[i] / sr[i];
        for (int e = r_pos[i]; e < r_pos[i + 1]; e++) t[r_idx[e]] += raw_val[e] * y;
      }
      for (int j = 0; j < p; j++)
        if (sc[j] > 0) sum += t[j] * t[j] / (sc[j] * sc[j]);

      rq++;
    }

    pi_out[k] = sum;
  }

  return NULL;
}


// The constrained inertia for each of n_perm permutations
// (out, length n_perm) of the n x px design x.
void cca_perm_inertia(
    const double *x, int px, perm_source_t *ps, int n_perm, double *out ) {

  int n = cc_n, p = cc_p;
  int n_threads = threads_for(
    cc_threads, ((double) cc_nnz + (double) n * px) * px * n_perm, 1e6 );
  if (n_threads > n_perm) n_threads = n_perm;
  if (n_threads < 1) n_threads = 1;

  size_t per = (size_t) n * px + p + n;
  pi_buf   = (double *) R_alloc(per * n_threads + 1, sizeof(double));
  pi_x     = x;
  pi_px    = px;
  pi_ps    = ps;
  pi_nperm = n_perm;
  pi_out   = out;

  run_parallel(perm_worker, n_threads, n_perm);
}



//======================================================
// Partial models, as vegan's RsquareAdj() for a cca with
// conditions: here the rows of Qbar are permuted against
// a fixed orthonormal design Q (n x rq; the first rz
// columns span the conditions), built with the unpermuted
// weights. For permutation k, with q~[p[a]] = q[a], the
// fitted inertia of each column is ||Qbar' q~||^2, where
//   Qbar' q~ = Dc^-1/2 P' Dr^-1/2 q~ - sqrt(c) (sqrt(r)' q~)
// (q~ is no longer orthogonal to sqrt(r)). out_all gets
// the sum over all rq columns, out_z over the first rz.
//======================================================
static const double *pf_q;
static int           pf_rq, pf_rz;
static double       *pf_z;

static void *perm_fixed_worker(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = cc_n, p = cc_p;

  size_t  per = (size_t) n + p + n;
  double *q   = pi_buf + per * thread_i;
  double *t   = q + n;
  int    *pm  = (int *) (t + p);

  for (int k = thread_i; k < pi_nperm; k += n_threads) {

    perm_source_get(pi_ps, k + 1, pm);

    double sum = 0, sum_z = 0;

    for (int c = 0; c < pf_rq; c++) {

      const double *qc = pf_q + (size_t) c * n;
      for (int a = 0; a < n; a++) q[pm[a]] = qc[a];

      double s = 0;
      for (int i = 0; i < n; i++) s += sr[i] * q[i];

      memset(t, 0, (size_t) p * sizeof(double));
      for (int i = 0; i < n; i++) {
        double y = q[i] / sr[i];
        for (int e = r_pos[i]; e < r_pos[i + 1]; e++) t[r_idx[e]] += raw_val[e] * y;
      }

      double ss = 0;
      for (int j = 0; j < p; j++) {
        if (!(sc[j] > 0)) continue;
        double v = t[j] / sc[j] - sc[j] * s;
        ss += v * v;
      }

      sum += ss;
      if (c < pf_rz) sum_z += ss;
    }

    pi_out[k] = sum;
    pf_z[k]   = sum_z;
  }

  return NULL;
}


void cca_perm_inertia_fixed(
    const double *Q, int rq, int rz, perm_source_t *ps, int n_perm,
    double *out_all, double *out_z ) {

  int n = cc_n, p = cc_p;
  int n_threads = threads_for(
    cc_threads, ((double) cc_nnz + n + p) * rq * n_perm, 1e6 );
  if (n_threads > n_perm) n_threads = n_perm;
  if (n_threads < 1) n_threads = 1;

  size_t per = (size_t) n + p + n;
  pi_buf   = (double *) R_alloc(per * n_threads + 1, sizeof(double));
  pf_q     = Q;
  pf_rq    = rq;
  pf_rz    = rz;
  pi_ps    = ps;
  pi_nperm = n_perm;
  pi_out   = out_all;
  pf_z     = out_z;

  run_parallel(perm_fixed_worker, n_threads, n_perm);
}
