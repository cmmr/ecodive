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
 * with thick restarts, which only needs products B * X. Those
 * are computed straight from the packed lower triangle of D
 * that R's `dist` objects hold.
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

#include <float.h> // DBL_EPSILON

// Rows of A * X per parallel task. Fixed, so that summation
// order never depends on the number of threads.
#define PCOA_ROWS 64

// Default block size: at least this many vectors, and at
// least k.
#define PCOA_BLOCK 4

// Largest Krylov basis kept before a thick restart.
#define PCOA_BASIS 60

// Converged when each wanted Ritz pair's residual norm,
// || B u - theta u ||, is at most PCOA_TOL * scale, where
// scale is the larger of trace(B) and the largest |theta|.
// Rounding in B * X is on the order of DBL_EPSILON * scale
// (|A|'s largest eigenvalue is about 2 * trace(B)), so
// residuals can stall above PCOA_TOL. A residual at most
// PCOA_TOL_STALL * scale that hasn't halved in PCOA_STALL
// steps is accepted as converged.
#define PCOA_TOL       1e-14
#define PCOA_TOL_STALL 1e-10
#define PCOA_STALL     10

// Give up after this many products with blocks of vectors.
#define PCOA_MAX_STEPS 10000

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
// Small dense helpers. Vectors have length n; matrices
// are column-major.
//------------------------------------------------------
static double dot(const double *x, const double *y, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) s += x[i] * y[i];
  return s;
}

static void center(double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) s += x[i];
  s /= n;
  for (int i = 0; i < n; i++) x[i] -= s;
}


//------------------------------------------------------
// W = B * V for `m` columns of length n (column-major),
// using xr and yr (n x PCOA_W each) as row-major
// scratch.
//------------------------------------------------------
static void mult_b(
    const double *V, double *W, int n, int m,
    double *xr, double *yr, int n_threads ) {

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
      center(w, n);
      for (int i = 0; i < n; i++)
        w[i] *= -0.5;
    }
  }
}


//------------------------------------------------------
// Uniform random values in [-1, 1) from a PCG stream.
//------------------------------------------------------
static void fill_random(double *x, size_t n, pcg32_random_t *rng) {
  for (size_t i = 0; i < n; i++)
    x[i] = pcg32_random_r(rng) / 2147483648.0 - 1.0;
}


//------------------------------------------------------
// Orthonormalize the `m` columns of R against the `s`
// columns of V, the ones vector, and each other, with
// two passes of classical Gram-Schmidt. A column that is
// (numerically) in their span is replaced by a random
// one; if that is too, the space is exhausted and the
// column is dropped. Returns how many columns are kept,
// moved to the front of R.
//------------------------------------------------------
static int orthonormalize(
    const double *V, int s, double *R, int m, int n,
    pcg32_random_t *rng ) {

  int kept = 0;

  for (int c = 0; c < m; c++) {

    double *r = R + (size_t) kept * n;
    if (c != kept) memcpy(r, R + (size_t) c * n, (size_t) n * sizeof(double));

    for (int attempt = 0; attempt < 2; attempt++) {

      if (attempt) fill_random(r, n, rng);

      double norm0 = sqrt(dot(r, r, n));

      for (int pass = 0; pass < 2; pass++) {
        center(r, n);
        for (int j = 0; j < s; j++) {
          const double *q = V + (size_t) j * n;
          double        h = dot(q, r, n);
          for (int i = 0; i < n; i++) r[i] -= h * q[i];
        }
        for (int j = 0; j < kept; j++) {
          const double *q = R + (size_t) j * n;
          double        h = dot(q, r, n);
          for (int i = 0; i < n; i++) r[i] -= h * q[i];
        }
      }

      double norm = sqrt(dot(r, r, n));

      if (norm > 1e-8 * norm0) {
        for (int i = 0; i < n; i++) r[i] /= norm;
        kept++;
        break;
      }
    }
  }

  return kept;
}


//------------------------------------------------------
// Eigen-decomposition of a symmetric s x s matrix A
// (column-major; only its lower triangle is read): a
// Householder reduction to tridiagonal form, then the QL
// algorithm with implicit shifts. After the public-domain
// JAMA routines tred2() and tql2(), which follow EISPACK.
// Eigenvalues go to `vals` in decreasing order, and the
// matching unit eigenvectors to the columns of `vecs`.
// `tmp` needs room for 2 * s doubles.
//------------------------------------------------------
static void sym_eigen(const double *A, int s, double *vals, double *vecs, double *tmp) {

  #define V(i, j) vecs[(size_t) (j) * s + (i)]

  double *d = vals;
  double *e = tmp;

  for (int j = 0; j < s; j++)
    for (int i = j; i < s; i++)
      V(i, j) = V(j, i) = A[(size_t) j * s + i];


  //--------------------------------------------------
  // Householder tridiagonalization (tred2).
  //--------------------------------------------------
  for (int j = 0; j < s; j++) d[j] = V(s - 1, j);

  for (int i = s - 1; i > 0; i--) {

    double scale = 0, h = 0;
    for (int k = 0; k < i; k++) scale += fabs(d[k]);

    if (scale == 0) {
      e[i] = d[i - 1];
      for (int j = 0; j < i; j++) {
        d[j]    = V(i - 1, j);
        V(i, j) = 0;
        V(j, i) = 0;
      }
    }
    else {
      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h    += d[k] * d[k];
      }
      double f = d[i - 1];
      double g = sqrt(h);
      if (f > 0) g = -g;
      e[i]     = scale * g;
      h        = h - f * g;
      d[i - 1] = f - g;
      for (int j = 0; j < i; j++) e[j] = 0;

      for (int j = 0; j < i; j++) {
        f       = d[j];
        V(j, i) = f;
        g       = e[j] + V(j, j) * f;
        for (int k = j + 1; k <= i - 1; k++) {
          g    += V(k, j) * d[k];
          e[k] += V(k, j) * f;
        }
        e[j] = g;
      }

      f = 0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f    += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) e[j] -= hh * d[j];

      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i - 1; k++)
          V(k, j) -= (f * e[k] + g * d[k]);
        d[j]    = V(i - 1, j);
        V(i, j) = 0;
      }
    }
    d[i] = h;
  }

  // Accumulate the transformations.
  for (int i = 0; i < s - 1; i++) {
    V(s - 1, i) = V(i, i);
    V(i, i)     = 1;
    double h    = d[i + 1];
    if (h != 0) {
      for (int k = 0; k <= i; k++) d[k] = V(k, i + 1) / h;
      for (int j = 0; j <= i; j++) {
        double g = 0;
        for (int k = 0; k <= i; k++) g += V(k, i + 1) * V(k, j);
        for (int k = 0; k <= i; k++) V(k, j) -= g * d[k];
      }
    }
    for (int k = 0; k <= i; k++) V(k, i + 1) = 0;
  }
  for (int j = 0; j < s; j++) {
    d[j]        = V(s - 1, j);
    V(s - 1, j) = 0;
  }
  V(s - 1, s - 1) = 1;
  e[0]            = 0;


  //--------------------------------------------------
  // QL with implicit shifts (tql2).
  //--------------------------------------------------
  for (int i = 1; i < s; i++) e[i - 1] = e[i];
  e[s - 1] = 0;

  double f = 0, tst1 = 0;

  for (int l = 0; l < s; l++) {

    // Find a small subdiagonal element.
    double t = fabs(d[l]) + fabs(e[l]);
    if (t > tst1) tst1 = t;
    int m = l;
    while (m < s - 1 && fabs(e[m]) > DBL_EPSILON * tst1) m++;

    // If m == l, d[l] is already an eigenvalue.
    if (m > l) {
      for (int iter = 0; iter < 100; iter++) {

        // Compute the implicit shift.
        double g = d[l];
        double p = (d[l + 1] - g) / (2 * e[l]);
        double r = hypot(p, 1);
        if (p < 0) r = -r;
        d[l]       = e[l] / (p + r);
        d[l + 1]   = e[l] * (p + r);
        double dl1 = d[l + 1];
        double h   = g - d[l];
        for (int i = l + 2; i < s; i++) d[i] -= h;
        f += h;

        // Implicit QL transformation.
        p = d[m];
        double c = 1, c2 = 1, c3 = 1, sn = 0, s2 = 0;
        double el1 = e[l + 1];
        for (int i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = sn;
          g  = c * e[i];
          h  = c * p;
          r  = hypot(p, e[i]);
          e[i + 1] = sn * r;
          sn = e[i] / r;
          c  = p / r;
          p  = c * d[i] - sn * g;
          d[i + 1] = h + sn * (c * g + sn * d[i]);
          for (int k = 0; k < s; k++) {
            h           = V(k, i + 1);
            V(k, i + 1) = sn * V(k, i) + c * h;
            V(k, i)     = c * V(k, i) - sn * h;
          }
        }
        p    = -sn * s2 * c3 * el1 * e[l] / dl1;
        e[l] = sn * p;
        d[l] = c * p;

        if (fabs(e[l]) <= DBL_EPSILON * tst1) break;
      }
    }
    d[l] = d[l] + f;
    e[l] = 0;
  }


  //--------------------------------------------------
  // Sort by decreasing eigenvalue (selection sort).
  //--------------------------------------------------
  for (int i = 0; i < s - 1; i++) {
    int k = i;
    for (int j = i + 1; j < s; j++)
      if (d[j] > d[k]) k = j;
    if (k != i) {
      double t = d[k]; d[k] = d[i]; d[i] = t;
      for (int r = 0; r < s; r++) {
        t = V(r, i); V(r, i) = V(r, k); V(r, k) = t;
      }
    }
  }

  #undef V
}


//------------------------------------------------------
// X = V * Y[, 0:m), for V n x s and Y s x s.
//------------------------------------------------------
static void combine(const double *V, int n, int s, const double *Y, int m, double *X) {
  for (int c = 0; c < m; c++) {
    double       *x = X + (size_t) c * n;
    const double *y = Y + (size_t) c * s;
    memset(x, 0, (size_t) n * sizeof(double));
    for (int j = 0; j < s; j++) {
      const double *v = V + (size_t) j * n;
      double        a = y[j];
      for (int i = 0; i < n; i++) x[i] += a * v[i];
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
  // Sizes: block p >= k; basis of at most M columns;
  // q Ritz vectors kept at each restart.
  //--------------------------------------------------
  int p = (k > PCOA_BLOCK) ? k : PCOA_BLOCK;
  int M = PCOA_BASIS;
  if (M < 2 * k + 2 * p) M = 2 * k + 2 * p;
  if (p > dims) p = dims;
  if (M > dims) M = dims;
  int q = (M - p + k) / 2;
  if (q < k) q = k;

  double *V   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *T   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *R   = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *xr  = (double *) R_alloc((size_t) n * PCOA_W, sizeof(double));
  double *yr  = (double *) R_alloc((size_t) n * PCOA_W, sizeof(double));
  double *H   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Hc  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Y   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *tmp = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *eig = (double *) R_alloc(M, sizeof(double));
  double *res = (double *) R_alloc(k, sizeof(double));
  double *u   = (double *) R_alloc(n, sizeof(double));
  double *bu  = (double *) R_alloc(n, sizeof(double));

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20260116u, 54u);


  //--------------------------------------------------
  // Starting block.
  //--------------------------------------------------
  fill_random(R, (size_t) n * p, &rng);
  int s_v = orthonormalize(V, 0, R, p, n, &rng); // columns in V
  int s_w = 0;                                   // ... with W = B * V
  memcpy(V, R, (size_t) n * s_v * sizeof(double));

  int    steps = 0, converged = 0, stalled = 0;
  double best  = INFINITY;

  while (1) {

    //------------------------------------------------
    // Multiply the newest block by B, and extend
    // H = V' * W to the new columns.
    //------------------------------------------------
    int m = s_v - s_w;
    mult_b(V + (size_t) n * s_w, W + (size_t) n * s_w, n, m, xr, yr, n_threads);
    steps++;

    for (int b = s_w; b < s_v; b++)
      for (int a = 0; a < s_v; a++)
        H[(size_t) b * M + a] = dot(V + (size_t) n * a, W + (size_t) n * b, n);

    for (int b = s_w; b < s_v; b++) {
      for (int a = 0; a < s_w; a++)
        H[(size_t) a * M + b] = H[(size_t) b * M + a];
      for (int a = s_w; a < b; a++) {
        double avg = (H[(size_t) b * M + a] + H[(size_t) a * M + b]) / 2;
        H[(size_t) b * M + a] = avg;
        H[(size_t) a * M + b] = avg;
      }
    }
    s_w = s_v;


    //------------------------------------------------
    // Rayleigh-Ritz.
    //------------------------------------------------
    for (int j = 0; j < s_w; j++)
      memcpy(Hc + (size_t) j * s_w, H + (size_t) j * M, (size_t) s_w * sizeof(double));
    sym_eigen(Hc, s_w, eig, Y, tmp);

    double scale = trace;
    for (int j = 0; j < s_w; j++)
      if (fabs(eig[j]) > scale) scale = fabs(eig[j]);

    double worst = 0;
    for (int j = 0; j < k; j++) {
      combine(V, n, s_w, Y + (size_t) j * s_w, 1, u);
      combine(W, n, s_w, Y + (size_t) j * s_w, 1, bu);
      for (int i = 0; i < n; i++) bu[i] -= eig[j] * u[i];
      res[j] = sqrt(dot(bu, bu, n)) / scale;
      if (res[j] > worst) worst = res[j];
    }

    if (worst < best / 2) { best = worst; stalled = 0; }
    else                  { stalled++;                }

    converged = worst <= PCOA_TOL ||
      (worst <= PCOA_TOL_STALL && stalled >= PCOA_STALL);

    if (converged || s_w == dims || steps >= PCOA_MAX_STEPS) break;

    R_CheckUserInterrupt();


    //------------------------------------------------
    // Next block: the newest products, orthonormalized
    // against the basis.
    //------------------------------------------------
    memcpy(R, W + (size_t) n * (s_w - m), (size_t) n * m * sizeof(double));
    int r = orthonormalize(V, s_v, R, m, n, &rng);
    if (r == 0) break; // # nocov (only if the space is exhausted)


    //------------------------------------------------
    // Thick restart: keep the top q Ritz vectors.
    //------------------------------------------------
    if (s_v + r > M) {

      int keep = (q + r <= M) ? q : M - r;

      combine(V, n, s_w, Y, keep, T);
      memcpy(V, T, (size_t) n * keep * sizeof(double));
      combine(W, n, s_w, Y, keep, T);
      memcpy(W, T, (size_t) n * keep * sizeof(double));

      for (int b = 0; b < keep; b++)
        for (int a = 0; a <= b; a++) {
          double h = (a == b)
            ? dot(V + (size_t) n * a, W + (size_t) n * b, n)
            : (dot(V + (size_t) n * a, W + (size_t) n * b, n) +
               dot(V + (size_t) n * b, W + (size_t) n * a, n)) / 2;
          H[(size_t) b * M + a] = h;
          H[(size_t) a * M + b] = h;
        }

      s_v = s_w = keep;
    }

    memcpy(V + (size_t) n * s_v, R, (size_t) n * r * sizeof(double));
    s_v += r;
  }


  //--------------------------------------------------
  // Coordinates, with each axis's sign fixed.
  //--------------------------------------------------
  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  double *points   = REAL(sexp_points);

  combine(V, n, s_w, Y, k, points);

  for (int j = 0; j < k; j++) {

    double *x = points + (size_t) n * j;

    int imax = 0;
    for (int i = 1; i < n; i++)
      if (fabs(x[i]) > fabs(x[imax])) imax = i;

    double mult = (eig[j] > 0) ? sqrt(eig[j]) : 0;
    if (x[imax] < 0) mult = -mult;

    for (int i = 0; i < n; i++) x[i] *= mult;

    REAL(sexp_eig)[j] = eig[j];
  }

  SEXP sexp_res = PROTECT(allocVector(REALSXP, k));
  memcpy(REAL(sexp_res), res, (size_t) k * sizeof(double));

  const char *names[] = { "points", "eig", "trace", "steps", "converged", "resid", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(trace));
  SET_VECTOR_ELT(sexp_result, 3, ScalarInteger(steps));
  SET_VECTOR_ELT(sexp_result, 4, ScalarLogical(converged));
  SET_VECTOR_ELT(sexp_result, 5, sexp_res);

  UNPROTECT(4);
  return sexp_result;
}
