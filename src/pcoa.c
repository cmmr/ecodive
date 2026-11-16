// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Principal Coordinates Analysis (PCoA), a.k.a. classical
 * multidimensional scaling.
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, the
 * coordinates come from the top eigenpairs of
 *
 *     B = -1/2 * J * A * J
 *
 * Only the k largest (algebraic) eigenpairs are wanted, with
 * k small. So instead of forming B (n x n) and decomposing it
 * in full like cmdscale(), they are found by block Lanczos
 * with thick restarts (lanczos.c), which only needs products
 * B * X. Those are computed straight from the packed lower
 * triangle of D that R's `dist` objects hold.
 *
 * Every basis vector is kept orthogonal to the ones vector, so
 * J * X = X, and B * X = -1/2 * J * (A * X): one product with
 * A, then subtracting each column's mean. (This is the same as
 * using A's row means and grand mean, since the means of A * X
 * are those row means times X.)
 *
 * Determinism: results are bit-identical for any number of
 * threads. A * X is computed in fixed blocks of PCOA_ROWS
 * rows; each block's rows are summed by one thread, in a fixed
 * order, with no shared accumulators. Everything else runs on
 * the main thread. The starting vectors come from a fixed PCG
 * stream.
 *
 * Sign convention: each axis is flipped so that its coordinate
 * with the largest absolute value is positive (the first such
 * sample, if tied).
 */

#include "ecodive.h"

// Rows of A * X per parallel task. Fixed, so that summation
// order never depends on the number of threads.
#define PCOA_ROWS 64


// Vectors per pass over D in mult_a(), and a macro that
// repeats F(v) for v = 0 to PCOA_W - 1.
#define PCOA_W 4
#define PCOA_REP(F) F(0) F(1) F(2) F(3)


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output rows that each task owns.
//------------------------------------------------------
static int           n_sams;   // Samples (rows of B).
static const double *d_vec;    // Packed lower triangle of D.
static const double *x_mtx;    // n_sams x PCOA_W, row-major.
static double       *y_mtx;    // n_sams x PCOA_W, row-major.
static double       *ss_vec;   // Per-column sums of D^2.


//------------------------------------------------------
// Offset of element (j + 1, j) in a dist vector. Column
// j's elements, rows j + 1 to n - 1, follow it in order.
//------------------------------------------------------
static inline size_t col_start(int j, int n) {
  return (size_t) j * (size_t) n - (size_t) j * (size_t) (j + 1) / 2;
}


//------------------------------------------------------
// Sum of squares of each column of the lower triangle.
// Their total, divided by n, is trace(B).
//------------------------------------------------------
static void *sum_squares(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = n_sams;

  for (int j = thread_i; j < n - 1; j += n_threads) {
    const double *d  = d_vec + col_start(j, n);
    int           nd = n - j - 1;
    double        ss = 0;
    for (int i = 0; i < nd; i++)
      ss += d[i] * d[i];
    ss_vec[j] = ss;
  }

  return NULL;
}


//------------------------------------------------------
// The runs that mult_a() reads from successive columns
// are far apart, which hardware prefetchers don't
// foresee. Ask for column j's run, rows i0 to i0 + nr - 1,
// ahead of time.
//------------------------------------------------------
#define PCOA_AHEAD 8

static inline void prefetch_run(int j, int i0, int nr, int n) {
#if defined(__GNUC__) || defined(__clang__)
  if (j < i0) {
    const char *p = (const char *) (d_vec + col_start(j, n) + (i0 - j - 1));
    for (size_t b = 0; b < (size_t) nr * sizeof(double); b += 64)
      __builtin_prefetch(p + b, 0, 0);
  }
#else
  (void) j; (void) i0; (void) nr; (void) n;
#endif
}


//------------------------------------------------------
// Y = A * X, where A = D^2 and X has PCOA_W columns.
// Each task computes every element of PCOA_ROWS rows of
// Y, reading the parts of D's packed lower triangle that
// those rows need as contiguous runs.
//------------------------------------------------------
static void *mult_a(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = n_sams;
  int n_tasks   = (n + PCOA_ROWS - 1) / PCOA_ROWS;

  const double *x = x_mtx;
  double       *y = y_mtx;

  // Unrolled loops over the PCOA_W vectors, which lets the
  // compiler keep xj[] and s[] in registers.
  #define X_J(v)  xj[v] = xc[v];
  #define AXPY(v) yr[v] += a * xj[v];
  #define DOT(v)  s[v]  += a * xi[v];
  #define BOTH(v) s[v]  += a * xi[v]; yi[v] += a * xj[v];
  #define ZERO(v) s[v]  = 0;
  #define ADD(v)  yc[v] += s[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * PCOA_ROWS;
    int i1 = (i0 + PCOA_ROWS < n) ? i0 + PCOA_ROWS : n;
    int nr = i1 - i0;

    double *yb = y + (size_t) i0 * PCOA_W;
    memset(yb, 0, (size_t) nr * PCOA_W * sizeof(double));

    double xj[PCOA_W], s[PCOA_W];

    // Columns left of the block: rows i0..i1-1 of column j
    // are contiguous in the dist vector.
    for (int j = 0; j < i0; j++) {
      const double *d  = d_vec + col_start(j, n) + (i0 - j - 1);
      const double *xc = x + (size_t) j * PCOA_W;
      prefetch_run(j + PCOA_AHEAD, i0, nr, n);
      PCOA_REP(X_J)
      for (int r = 0; r < nr; r++) {
        double  a  = d[r] * d[r];
        double *yr = yb + (size_t) r * PCOA_W;
        PCOA_REP(AXPY)
      }
    }

    // The block's own columns: below-diagonal elements
    // within the block feed both of their rows; the rest
    // of each column feeds row j as a dot product.
    for (int j = i0; j < i1; j++) {

      const double *d  = d_vec + col_start(j, n) - (j + 1);
      const double *xc = x + (size_t) j * PCOA_W;
      PCOA_REP(X_J)
      PCOA_REP(ZERO)

      for (int i = j + 1; i < i1; i++) {
        double        a  = d[i] * d[i];
        const double *xi = x + (size_t) i * PCOA_W;
        double       *yi = y + (size_t) i * PCOA_W;
        PCOA_REP(BOTH)
      }

      for (int i = i1; i < n; i++) {
        double        a  = d[i] * d[i];
        const double *xi = x + (size_t) i * PCOA_W;
        PCOA_REP(DOT)
      }

      double *yc = y + (size_t) j * PCOA_W;
      PCOA_REP(ADD)
    }
  }

  #undef X_J
  #undef AXPY
  #undef DOT
  #undef BOTH
  #undef ZERO
  #undef ADD

  return NULL;
}
//------------------------------------------------------
// W = B * V for `m` columns of length n (column-major),
// using xr and yr (n x PCOA_W each) as row-major
// scratch. A lanczos_op_t.
//------------------------------------------------------
typedef struct {
  double *xr;
  double *yr;
  int     n_threads;
} pcoa_ctx_t;

static void mult_b(const double *V, double *W, int m, void *ctx) {

  int     n         = n_sams;
  double *xr        = ((pcoa_ctx_t *) ctx)->xr;
  double *yr        = ((pcoa_ctx_t *) ctx)->yr;
  int     n_threads = ((pcoa_ctx_t *) ctx)->n_threads;

  // A column's products don't depend on the columns beside
  // it, so padding with zeros changes nothing.
  for (int v0 = 0; v0 < m; v0 += PCOA_W) {

    int mc = (m - v0 < PCOA_W) ? m - v0 : PCOA_W;

    for (int i = 0; i < n; i++)
      for (int v = 0; v < PCOA_W; v++)
        xr[(size_t) i * PCOA_W + v] = (v < mc) ? V[(size_t) (v0 + v) * n + i] : 0;

    x_mtx = xr;
    y_mtx = yr;
    run_parallel(mult_a, n_threads, n);

    for (int v = 0; v < mc; v++) {
      double *w = W + (size_t) (v0 + v) * n;
      for (int i = 0; i < n; i++)
        w[i] = yr[(size_t) i * PCOA_W + v];
      lz_center(w, n);
      for (int i = 0; i < n; i++)
        w[i] *= -0.5;
    }
  }
}


//======================================================
// R interface. `sexp_dist` is a double vector holding a
// dist object's packed lower triangle for `sexp_n`
// samples (n >= 2). Returns list(points, eig, trace,
// steps, converged, resid), with points an n x k matrix of
// eigenvectors scaled by sqrt(eig), or 0 where eig <= 0.
//======================================================
SEXP C_pcoa(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_n_threads) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_n_threads);
  int dims      = n - 1; // B's rank is at most n - 1.

  if (n < 2 || k < 1 || k > dims || n_threads < 1 ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_pcoa()."); // # nocov

  n_sams = n;
  d_vec  = REAL(sexp_dist);


  //--------------------------------------------------
  // trace(B) = (1/n) * sum_{i<j} d_ij^2, summed per
  // column in parallel, then over columns in order.
  //--------------------------------------------------
  ss_vec = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, n_threads, n);

  double total = 0;
  for (int j = 0; j < n - 1; j++) total += ss_vec[j];
  if (!R_FINITE(total))
    error("`dm` must not contain NA, NaN, or infinite values.");
  double trace = total / n;


  //--------------------------------------------------
  // The top k eigenpairs of B.
  //--------------------------------------------------
  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_res    = PROTECT(allocVector(REALSXP, k));
  double *points   = REAL(sexp_points);
  double *eig      = REAL(sexp_eig);

  pcoa_ctx_t ctx = {
    .xr        = (double *) R_alloc((size_t) n * PCOA_W, sizeof(double)),
    .yr        = (double *) R_alloc((size_t) n * PCOA_W, sizeof(double)),
    .n_threads = n_threads };

  lanczos_t lz = {
    .n = n, .k = k, .dims = dims, .deflate = 1, .scale = trace,
    .op = mult_b, .ctx = &ctx,
    .vecs = points, .vals = eig, .resid = REAL(sexp_res) };

  lanczos(&lz);


  //--------------------------------------------------
  // Coordinates, with each axis's sign fixed.
  //--------------------------------------------------
  for (int j = 0; j < k; j++) {

    double *x = points + (size_t) n * j;

    int imax = 0;
    for (int i = 1; i < n; i++)
      if (fabs(x[i]) > fabs(x[imax])) imax = i;

    double mult = (eig[j] > 0) ? sqrt(eig[j]) : 0;
    if (x[imax] < 0) mult = -mult;

    for (int i = 0; i < n; i++) x[i] *= mult;
  }

  const char *names[] = { "points", "eig", "trace", "steps", "converged", "resid", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(trace));
  SET_VECTOR_ELT(sexp_result, 3, ScalarInteger(lz.steps));
  SET_VECTOR_ELT(sexp_result, 4, ScalarLogical(lz.converged));
  SET_VECTOR_ELT(sexp_result, 5, sexp_res);

  UNPROTECT(4);
  return sexp_result;
}
