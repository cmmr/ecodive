// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Distance-based redundancy analysis (db-RDA), as in
 * vegan::dbrda() (McArdle & Anderson 2001).
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, Gower's
 * centred matrix is
 *
 *     G = -1/2 * J * A * J
 *
 * The design's columns (centred) are orthonormalized by
 * Gram-Schmidt into Q = [Qz, Qx]: Qz spans the conditioning
 * variables, and Qx the constraints after removing Qz. Then
 *
 *   total inertia        trace(G)
 *   conditional inertia  trace(Qz' G Qz)
 *   constrained axes     eigenpairs (lambda, y) of Qx' G Qx,
 *                        with u = Qx y
 *   unconstrained axes   eigenpairs of (I-H) G (I-H), with
 *                        H = 11'/n + Q Q'
 *
 * vegan builds G and H G H in full (n x n) and decomposes
 * them. Here the constrained part only needs G * Q (a few
 * products with G, computed straight from the packed lower
 * triangle of D) and a p x p eigenproblem. Unconstrained axes
 * are only computed when there are fewer positive
 * constrained axes than requested, by block Lanczos with
 * thick restarts on (I-H) G (I-H), keeping every basis
 * vector orthogonal to Q and the ones vector.
 *
 * Site scores for constrained axes are vegan's "wa" scores,
 * (I - Qz Qz') G u / lambda; for unconstrained axes they are
 * the unit eigenvectors (vegan's CA$u). The R side applies
 * vegan's scaling constant.
 *
 * The products with G, the small dense helpers and the
 * Lanczos iteration are adapted from src/pcoa.c, which
 * isn't to be edited while sibling methods are written.
 * TODO: consolidate with pcoa.c (and pca's Lanczos) into a
 * shared file.
 *
 * Determinism: results are bit-identical for any number of
 * threads. G * X is computed in fixed blocks of DBRDA_ROWS
 * rows, each summed by one thread in a fixed order. All else
 * runs on the main thread. Lanczos starting vectors come from
 * a fixed PCG stream.
 *
 * Sign convention: each axis is flipped so that its score
 * with the largest absolute value is positive (the first such
 * sample, if tied).
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// Rows of A * X per parallel task. Fixed, so that summation
// order never depends on the number of threads.
#define DBRDA_ROWS 64

// Lanczos: default block size, largest basis before a thick
// restart, convergence tolerances and step limit. As in
// pcoa.c.
#define DBRDA_BLOCK     4
#define DBRDA_BASIS     60
#define DBRDA_TOL       1e-14
#define DBRDA_TOL_STALL 1e-10
#define DBRDA_STALL     10
#define DBRDA_MAX_STEPS 10000

// A design column is aliased (redundant) when Gram-Schmidt
// leaves less than this fraction of its centred norm: the
// same rule and tolerance as R's qr(), which vegan uses.
#define DBRDA_QR_TOL 1e-7

// Vectors per pass over D in mult_a(), and a macro that
// repeats F(v) for v = 0 to DBRDA_W - 1.
#define DBRDA_W 4
#define DBRDA_REP(F) F(0) F(1) F(2) F(3)


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output rows that each task owns.
//------------------------------------------------------
static int           n_sams;   // Samples (rows of G).
static const double *d_vec;    // Packed lower triangle of D.
static const double *x_mtx;    // n_sams x DBRDA_W, row-major.
static double       *y_mtx;    // n_sams x DBRDA_W, row-major.
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
// Their total, divided by n, is trace(G).
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
// Prefetch column j's run, rows i0 to i0 + nr - 1, ahead
// of time (see pcoa.c).
//------------------------------------------------------
#define DBRDA_AHEAD 8

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
// Y = A * X, where A = D^2 and X has DBRDA_W columns.
// Each task computes every element of DBRDA_ROWS rows of
// Y, reading the parts of D's packed lower triangle that
// those rows need as contiguous runs.
//------------------------------------------------------
static void *mult_a(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = n_sams;
  int n_tasks   = (n + DBRDA_ROWS - 1) / DBRDA_ROWS;

  const double *x = x_mtx;
  double       *y = y_mtx;

  #define X_J(v)  xj[v] = xc[v];
  #define AXPY(v) yr[v] += a * xj[v];
  #define DOT(v)  s[v]  += a * xi[v];
  #define BOTH(v) s[v]  += a * xi[v]; yi[v] += a * xj[v];
  #define ZERO(v) s[v]  = 0;
  #define ADD(v)  yc[v] += s[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * DBRDA_ROWS;
    int i1 = (i0 + DBRDA_ROWS < n) ? i0 + DBRDA_ROWS : n;
    int nr = i1 - i0;

    double *yb = y + (size_t) i0 * DBRDA_W;
    memset(yb, 0, (size_t) nr * DBRDA_W * sizeof(double));

    double xj[DBRDA_W], s[DBRDA_W];

    for (int j = 0; j < i0; j++) {
      const double *d  = d_vec + col_start(j, n) + (i0 - j - 1);
      const double *xc = x + (size_t) j * DBRDA_W;
      prefetch_run(j + DBRDA_AHEAD, i0, nr, n);
      DBRDA_REP(X_J)
      for (int r = 0; r < nr; r++) {
        double  a  = d[r] * d[r];
        double *yr = yb + (size_t) r * DBRDA_W;
        DBRDA_REP(AXPY)
      }
    }

    for (int j = i0; j < i1; j++) {

      const double *d  = d_vec + col_start(j, n) - (j + 1);
      const double *xc = x + (size_t) j * DBRDA_W;
      DBRDA_REP(X_J)
      DBRDA_REP(ZERO)

      for (int i = j + 1; i < i1; i++) {
        double        a  = d[i] * d[i];
        const double *xi = x + (size_t) i * DBRDA_W;
        double       *yi = y + (size_t) i * DBRDA_W;
        DBRDA_REP(BOTH)
      }

      for (int i = i1; i < n; i++) {
        double        a  = d[i] * d[i];
        const double *xi = x + (size_t) i * DBRDA_W;
        DBRDA_REP(DOT)
      }

      double *yc = y + (size_t) j * DBRDA_W;
      DBRDA_REP(ADD)
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

// x -= F * (F' x), for the f orthonormal columns of F.
static void project_out(const double *F, int f, double *x, int n) {
  for (int j = 0; j < f; j++) {
    const double *q = F + (size_t) j * n;
    double        h = dot(q, x, n);
    for (int i = 0; i < n; i++) x[i] -= h * q[i];
  }
}


//------------------------------------------------------
// W = G * V for `m` columns of length n (column-major),
// using xr and yr (n x DBRDA_W each) as row-major
// scratch. If f > 0, each column of W then has its part
// in the span of F's (orthonormal) columns removed:
// W = (I - FF') G V.
//------------------------------------------------------
static void mult_g(
    const double *V, double *W, int n, int m,
    const double *F, int f,
    double *xr, double *yr, int n_threads ) {

  for (int v0 = 0; v0 < m; v0 += DBRDA_W) {

    int mc = (m - v0 < DBRDA_W) ? m - v0 : DBRDA_W;

    for (int i = 0; i < n; i++)
      for (int v = 0; v < DBRDA_W; v++)
        xr[(size_t) i * DBRDA_W + v] = (v < mc) ? V[(size_t) (v0 + v) * n + i] : 0;

    x_mtx = xr;
    y_mtx = yr;
    run_parallel(mult_a, n_threads, n);

    for (int v = 0; v < mc; v++) {
      double *w = W + (size_t) (v0 + v) * n;
      for (int i = 0; i < n; i++)
        w[i] = yr[(size_t) i * DBRDA_W + v];
      center(w, n);
      for (int i = 0; i < n; i++)
        w[i] *= -0.5;
      project_out(F, f, w, n);
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
// Orthonormalize the `m` columns of R against the ones
// vector, the `f` columns of F, the `s` columns of V, and
// each other, with two passes of classical Gram-Schmidt.
// A column that is (numerically) in their span is replaced
// by a random one; if that is too, the column is dropped.
// Returns how many columns are kept, moved to the front.
//------------------------------------------------------
static int orthonormalize(
    const double *F, int f, const double *V, int s,
    double *R, int m, int n, pcg32_random_t *rng ) {

  int kept = 0;

  for (int c = 0; c < m; c++) {

    double *r = R + (size_t) kept * n;
    if (c != kept) memcpy(r, R + (size_t) c * n, (size_t) n * sizeof(double));

    for (int attempt = 0; attempt < 2; attempt++) {

      if (attempt) fill_random(r, n, rng);

      double norm0 = sqrt(dot(r, r, n));

      for (int pass = 0; pass < 2; pass++) {
        center(r, n);
        project_out(F, f, r, n);
        project_out(V, s, r, n);
        project_out(R, kept, r, n);
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
// (column-major; only its lower triangle is read), by
// Householder tridiagonalization and QL with implicit
// shifts (JAMA's tred2() and tql2(), after EISPACK), as in
// pcoa.c. Eigenvalues go to `vals` in decreasing order,
// unit eigenvectors to the columns of `vecs`. `tmp` needs
// room for 2 * s doubles.
//------------------------------------------------------
static void sym_eigen(const double *A, int s, double *vals, double *vecs, double *tmp) {

  #define V(i, j) vecs[(size_t) (j) * s + (i)]

  double *d = vals;
  double *e = tmp;

  for (int j = 0; j < s; j++)
    for (int i = j; i < s; i++)
      V(i, j) = V(j, i) = A[(size_t) j * s + i];

  // Householder tridiagonalization (tred2).
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

  // QL with implicit shifts (tql2).
  for (int i = 1; i < s; i++) e[i - 1] = e[i];
  e[s - 1] = 0;

  double f = 0, tst1 = 0;

  for (int l = 0; l < s; l++) {

    double t = fabs(d[l]) + fabs(e[l]);
    if (t > tst1) tst1 = t;
    int m = l;
    while (m < s - 1 && fabs(e[m]) > DBL_EPSILON * tst1) m++;

    if (m > l) {
      for (int iter = 0; iter < 100; iter++) {

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

  // Sort by decreasing eigenvalue (selection sort).
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
// X = V * Y[, 0:m), for V n x s and Y s x s (or any Y
// with leading dimension s).
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


//------------------------------------------------------
// The k largest eigenpairs of (I - FF') G (I - FF'),
// restricted to vectors orthogonal to the ones vector,
// by block Lanczos with thick restarts (as in pcoa.c).
// `dims` is the dimension of that space. Writes unit
// eigenvectors to `vecs` (n x k) and eigenvalues to
// `vals`. `scale` (trace(G)) sets the convergence
// tolerance. Returns 1 if converged.
//------------------------------------------------------
static int lanczos(
    int n, int k, int dims, const double *F, int f, double scale0,
    int n_threads, double *vecs, double *vals, int *n_steps ) {

  int p = (k > DBRDA_BLOCK) ? k : DBRDA_BLOCK;
  int M = DBRDA_BASIS;
  if (M < 2 * k + 2 * p) M = 2 * k + 2 * p;
  if (p > dims) p = dims;
  if (M > dims) M = dims;
  int q = (M - p + k) / 2;
  if (q < k) q = k;

  double *V   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *T   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *R   = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *xr  = (double *) R_alloc((size_t) n * DBRDA_W, sizeof(double));
  double *yr  = (double *) R_alloc((size_t) n * DBRDA_W, sizeof(double));
  double *H   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Hc  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Y   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *tmp = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *eig = (double *) R_alloc(M, sizeof(double));
  double *u   = (double *) R_alloc(n, sizeof(double));
  double *bu  = (double *) R_alloc(n, sizeof(double));

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20260116u, 54u);

  fill_random(R, (size_t) n * p, &rng);
  int s_v = orthonormalize(F, f, V, 0, R, p, n, &rng);
  int s_w = 0;
  memcpy(V, R, (size_t) n * s_v * sizeof(double));

  int    steps = 0, converged = 0, stalled = 0;
  double best  = INFINITY;

  while (1) {

    // Multiply the newest block, and extend H = V' * W.
    int m = s_v - s_w;
    mult_g(V + (size_t) n * s_w, W + (size_t) n * s_w, n, m, F, f, xr, yr, n_threads);
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

    // Rayleigh-Ritz.
    for (int j = 0; j < s_w; j++)
      memcpy(Hc + (size_t) j * s_w, H + (size_t) j * M, (size_t) s_w * sizeof(double));
    sym_eigen(Hc, s_w, eig, Y, tmp);

    double scale = scale0;
    for (int j = 0; j < s_w; j++)
      if (fabs(eig[j]) > scale) scale = fabs(eig[j]);

    double worst = 0;
    for (int j = 0; j < k; j++) {
      combine(V, n, s_w, Y + (size_t) j * s_w, 1, u);
      combine(W, n, s_w, Y + (size_t) j * s_w, 1, bu);
      for (int i = 0; i < n; i++) bu[i] -= eig[j] * u[i];
      double res = sqrt(dot(bu, bu, n)) / scale;
      if (res > worst) worst = res;
    }

    if (worst < best / 2) { best = worst; stalled = 0; }
    else                  { stalled++;                }

    converged = worst <= DBRDA_TOL ||
      (worst <= DBRDA_TOL_STALL && stalled >= DBRDA_STALL);

    if (converged || s_w == dims || steps >= DBRDA_MAX_STEPS) break;

    R_CheckUserInterrupt();

    // Next block: the newest products, orthonormalized.
    memcpy(R, W + (size_t) n * (s_w - m), (size_t) n * m * sizeof(double));
    int r = orthonormalize(F, f, V, s_v, R, m, n, &rng);
    if (r == 0) break; // # nocov (only if the space is exhausted)

    // Thick restart: keep the top q Ritz vectors.
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

  // A converged run (or exhausted space) has s_w >= k.
  if (s_w < k) k = s_w; // # nocov

  combine(V, n, s_w, Y, k, vecs);
  memcpy(vals, eig, (size_t) k * sizeof(double));

  *n_steps = steps;
  return converged || s_w == dims;
}


//------------------------------------------------------
// Flip x (and y, if given) so that x's element with the
// largest absolute value is positive.
//------------------------------------------------------
static void fix_sign(double *x, int n, double *y, int ny) {
  int imax = 0;
  for (int i = 1; i < n; i++)
    if (fabs(x[i]) > fabs(x[imax])) imax = i;
  if (x[imax] < 0) {
    for (int i = 0; i < n;  i++) x[i] = -x[i];
    for (int i = 0; i < ny; i++) y[i] = -y[i];
  }
}


//======================================================
// R interface.
//
//   sexp_dist   Packed lower triangle of D (double),
//               for n >= 2 samples.
//   sexp_x      n x p design matrix (double): the first
//               `nz` columns are conditioning variables,
//               the rest constraints. Centred here.
//   sexp_k      Axes wanted (>= 1).
//   sexp_adj    Multiplier for G (vegan's 1/(n-1) for
//               "mean squared distances", else 1).
//
// Returns list(
//   scores        n x m: wa scores for constrained axes,
//                 unit eigenvectors for unconstrained ones,
//   eig           m eigenvalues, times adj,
//   n_con         how many of the m axes are constrained,
//   con_eig       all non-zero constrained eigenvalues,
//   total, conditional, constrained   inertia,
//   kept          per design column, 1 if not aliased,
//   biplot        (kept constraint columns) x n_con,
//   u             n x (positive constrained axes) unit
//                 eigenvectors (lc scores, unscaled),
//   steps, converged )
//======================================================
SEXP C_dbrda(
    SEXP sexp_dist, SEXP sexp_n, SEXP sexp_x, SEXP sexp_nz,
    SEXP sexp_k, SEXP sexp_adj, SEXP sexp_n_threads ) {

  int    n         = asInteger(sexp_n);
  int    nz        = asInteger(sexp_nz);
  int    k         = asInteger(sexp_k);
  double adj       = asReal(sexp_adj);
  int    n_threads = asInteger(sexp_n_threads);
  int    p         = (n > 0) ? (int) (XLENGTH(sexp_x) / n) : 0;

  if (n < 2 || k < 1 || n_threads < 1 || nz < 0 || nz > p ||
      !isReal(sexp_x) || XLENGTH(sexp_x) != (R_xlen_t) n * p ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2 || !(adj > 0))
    error("Invalid arguments to C_dbrda()."); // # nocov

  const double ZERO = sqrt(DBL_EPSILON); // vegan's threshold

  n_sams = n;
  d_vec  = REAL(sexp_dist);


  //--------------------------------------------------
  // trace(G) = (1/n) * sum_{i<j} d_ij^2.
  //--------------------------------------------------
  ss_vec = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, n_threads, n);

  double ss = 0;
  for (int j = 0; j < n - 1; j++) ss += ss_vec[j];
  if (!R_FINITE(ss))
    error("`dm` must not contain NA, NaN, or infinite values.");
  double trace_g = ss / n;      // unscaled
  double total   = trace_g * adj;


  //--------------------------------------------------
  // Centre the design's columns, and orthonormalize them
  // in order by Gram-Schmidt (two passes), dropping the
  // aliased ones. Xc keeps the centred columns.
  //--------------------------------------------------
  double *Xc   = (double *) R_alloc((size_t) n * (p + 1), sizeof(double));
  double *Q    = (double *) R_alloc((size_t) n * (p + 1), sizeof(double));
  int    *kept = (int *)    R_alloc(p + 1, sizeof(int));
  int     rz = 0, rx = 0;

  memcpy(Xc, REAL(sexp_x), (size_t) n * p * sizeof(double));

  for (int c = 0; c < p; c++) {

    double *xc = Xc + (size_t) c * n;
    double *q  = Q  + (size_t) (rz + rx) * n;

    for (int i = 0; i < n; i++)
      if (!R_FINITE(xc[i]))
        error("The design must not contain NA, NaN, or infinite values.");

    center(xc, n);
    memcpy(q, xc, (size_t) n * sizeof(double));

    double norm0 = sqrt(dot(q, q, n));
    project_out(Q, rz + rx, q, n);
    project_out(Q, rz + rx, q, n);
    double norm = sqrt(dot(q, q, n));

    kept[c] = norm0 > 0 && norm >= DBRDA_QR_TOL * norm0;
    if (!kept[c]) continue;

    for (int i = 0; i < n; i++) q[i] /= norm;
    if (c < nz) rz++; else rx++;
  }

  const double *Qz = Q;
  const double *Qx = Q + (size_t) n * rz;
  int           rq = rz + rx;


  //--------------------------------------------------
  // GQ = G * Q, then the conditional inertia and the
  // constrained matrix C = Qx' G Qx.
  //--------------------------------------------------
  double *xr = (double *) R_alloc((size_t) n * DBRDA_W, sizeof(double));
  double *yr = (double *) R_alloc((size_t) n * DBRDA_W, sizeof(double));
  double *GQ = (double *) R_alloc((size_t) n * (rq + 1), sizeof(double));

  mult_g(Q, GQ, n, rq, NULL, 0, xr, yr, n_threads);
  const double *GQx = GQ + (size_t) n * rz;

  double conditional = 0;
  for (int j = 0; j < rz; j++)
    conditional += dot(Qz + (size_t) j * n, GQ + (size_t) j * n, n);
  conditional *= adj;
  if (conditional < ZERO) conditional = 0;

  double *C    = (double *) R_alloc((size_t) rx * rx + 1, sizeof(double));
  double *Y    = (double *) R_alloc((size_t) rx * rx + 1, sizeof(double));
  double *lam  = (double *) R_alloc(rx + 1, sizeof(double));
  double *tmp  = (double *) R_alloc(2 * (size_t) rx + 2, sizeof(double));

  for (int b = 0; b < rx; b++)
    for (int a = b; a < rx; a++) {
      double h = (a == b)
        ? dot(Qx + (size_t) n * a, GQx + (size_t) n * b, n)
        : (dot(Qx + (size_t) n * a, GQx + (size_t) n * b, n) +
           dot(Qx + (size_t) n * b, GQx + (size_t) n * a, n)) / 2;
      C[(size_t) b * rx + a] = h;
      C[(size_t) a * rx + b] = h;
    }

  if (rx > 0) sym_eigen(C, rx, lam, Y, tmp);


  //--------------------------------------------------
  // vegan's rule for zero eigenvalues: |lambda| below
  // max(ZERO, ZERO * largest). The rest, positive and
  // negative, make up the constrained inertia.
  //--------------------------------------------------
  int    n_nonzero = 0, n_pos = 0;
  double constrained = 0;
  double con_zero    = (rx > 0 && ZERO * lam[0] * adj > ZERO) ? ZERO * lam[0] * adj : ZERO;

  int *con_idx = (int *) R_alloc(rx + 1, sizeof(int)); // nonzero eigenvalues' indices
  for (int j = 0; j < rx; j++) {
    double l = lam[j] * adj;
    if (fabs(l) < con_zero) continue;
    con_idx[n_nonzero++] = j;
    constrained += l;
    if (l > 0) n_pos++;
  }

  int kc = (k < n_pos) ? k : n_pos;


  //--------------------------------------------------
  // Unconstrained axes, if more are wanted: the top
  // eigenpairs of (I-H) G (I-H), in the space orthogonal
  // to Q and the ones vector.
  //--------------------------------------------------
  int     kr = 0, steps = 0, converged = 1;
  int     dims_r = n - 1 - rq;
  double *r_vecs = NULL, *r_vals = NULL;

  if (kc < k && dims_r > 0) {

    kr     = (k - kc < dims_r) ? k - kc : dims_r;
    r_vecs = (double *) R_alloc((size_t) n * kr, sizeof(double));
    r_vals = (double *) R_alloc(kr, sizeof(double));

    converged = lanczos(n, kr, dims_r, Q, rq, trace_g, n_threads, r_vecs, r_vals, &steps);

    // Keep the positive, non-zero ones (vegan's ordResid).
    double res_zero = (ZERO * r_vals[0] * adj > ZERO) ? ZERO * r_vals[0] * adj : ZERO;
    int    keep     = 0;
    while (keep < kr && r_vals[keep] * adj >= res_zero) keep++;
    kr = keep;
  }


  //--------------------------------------------------
  // Output.
  //--------------------------------------------------
  int m = kc + kr;

  int n_bp = 0;
  for (int c = nz; c < p; c++) n_bp += kept[c];

  SEXP sexp_scores  = PROTECT(allocMatrix(REALSXP, n, m));
  SEXP sexp_eig     = PROTECT(allocVector(REALSXP, m));
  SEXP sexp_con_eig = PROTECT(allocVector(REALSXP, n_nonzero));
  SEXP sexp_kept    = PROTECT(allocVector(LGLSXP,  p));
  SEXP sexp_biplot  = PROTECT(allocMatrix(REALSXP, n_bp, kc));
  SEXP sexp_u       = PROTECT(allocMatrix(REALSXP, n, n_pos));

  double *scores = REAL(sexp_scores);
  double *eig    = REAL(sexp_eig);
  double *biplot = REAL(sexp_biplot);
  double *u      = REAL(sexp_u);

  for (int j = 0; j < n_nonzero; j++)
    REAL(sexp_con_eig)[j] = lam[con_idx[j]] * adj;

  for (int c = 0; c < p; c++)
    LOGICAL(sexp_kept)[c] = kept[c];

  // Constrained axes: u = Qx y; wa = (I - Qz Qz') G u / lambda.
  // u is returned for every positive axis (vegan's centroids
  // use them all), wa only for the first kc.
  for (int j = 0; j < n_pos; j++) {

    int           e  = con_idx[j];
    const double *y  = Y + (size_t) e * rx;
    double       *uj = u + (size_t) j * n;

    combine(Qx, n, rx, y, 1, uj);
    if (j >= kc) continue;

    double *wa = scores + (size_t) j * n;
    combine(GQx, n, rx, y, 1, wa);
    project_out(Qz, rz, wa, n);
    for (int i = 0; i < n; i++) wa[i] /= lam[e];

    fix_sign(wa, n, uj, n);
    eig[j] = lam[e] * adj;

    // Biplot scores: correlations of the centred constraint
    // columns with u.
    int row = 0;
    for (int c = nz; c < p; c++) {
      if (!kept[c]) continue;
      const double *xc = Xc + (size_t) c * n;
      biplot[(size_t) j * n_bp + row++] = dot(xc, uj, n) / sqrt(dot(xc, xc, n));
    }
  }

  // Unconstrained axes.
  for (int j = 0; j < kr; j++) {
    double *x = scores + (size_t) (kc + j) * n;
    memcpy(x, r_vecs + (size_t) j * n, (size_t) n * sizeof(double));
    fix_sign(x, n, NULL, 0);
    eig[kc + j] = r_vals[j] * adj;
  }

  const char *names[] = {
    "scores", "eig", "n_con", "con_eig", "total", "conditional",
    "constrained", "kept", "biplot", "u", "steps", "converged", "" };

  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result,  0, sexp_scores);
  SET_VECTOR_ELT(sexp_result,  1, sexp_eig);
  SET_VECTOR_ELT(sexp_result,  2, ScalarInteger(kc));
  SET_VECTOR_ELT(sexp_result,  3, sexp_con_eig);
  SET_VECTOR_ELT(sexp_result,  4, ScalarReal(total));
  SET_VECTOR_ELT(sexp_result,  5, ScalarReal(conditional));
  SET_VECTOR_ELT(sexp_result,  6, ScalarReal(constrained));
  SET_VECTOR_ELT(sexp_result,  7, sexp_kept);
  SET_VECTOR_ELT(sexp_result,  8, sexp_biplot);
  SET_VECTOR_ELT(sexp_result,  9, sexp_u);
  SET_VECTOR_ELT(sexp_result, 10, ScalarInteger(steps));
  SET_VECTOR_ELT(sexp_result, 11, ScalarLogical(converged));

  UNPROTECT(7);
  return sexp_result;
}
