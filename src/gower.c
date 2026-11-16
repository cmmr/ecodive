// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Products with Gower's centered matrix, for pcoa() and
 * dbrda(). With A = D^2 (elementwise) and J = I - 11'/n,
 *
 *     G = -1/2 * J * A * J
 *
 * gower_mult() computes W = (I - F F') * G * V, where F is an
 * optional orthonormal basis (dbrda's constraints) and V's
 * columns are orthogonal to the ones vector (so J * V = V):
 * one product with A, straight from the packed lower
 * triangle of D that R's `dist` objects hold, then
 * subtracting each column's mean, halving, and projecting
 * out F.
 *
 * Determinism: results are bit-identical for any number of
 * threads. A * X is computed in fixed blocks of GOWER_ROWS
 * rows; each block's rows are summed by one thread, in a
 * fixed order, with no shared accumulators. Each element of
 * A * X is summed in the same order whatever the number of
 * vectors per pass. Everything else runs on the main thread.
 */

#include "ecodive.h"

// Rows of A * X per parallel task. Fixed, so that summation
// order never depends on the number of threads.
#define GOWER_ROWS 64

// Vectors per unrolled group, and a macro that repeats F(v)
// for v = 0 to GOWER_W - 1.
#define GOWER_W 4
#define GOWER_REP(F) F(0) F(1) F(2) F(3)

// How many columns ahead to prefetch.
#define GOWER_AHEAD 8


//------------------------------------------------------
// Shared with the worker threads. Read-only, except for
// the output rows that each task owns.
//------------------------------------------------------
static int           g_n;      // Samples.
static const double *g_dist;   // Packed lower triangle of D.
static const double *g_x;      // n x (GOWER_W * g_ng), row-major.
static double       *g_y;      // n x (GOWER_W * g_ng), row-major.
static int           g_ng;     // Groups of GOWER_W vectors.
static double       *g_ss;     // Per-column sums of D^2.


//------------------------------------------------------
// Sum of squares of each column of the lower triangle.
//------------------------------------------------------
static void *sum_squares(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = g_n;

  for (int j = thread_i; j < n - 1; j += n_threads) {
    const double *d  = g_dist + dist_col_start(j, n);
    int           nd = n - j - 1;
    double        ss = 0;
    for (int i = 0; i < nd; i++)
      ss += d[i] * d[i];
    g_ss[j] = ss;
  }

  return NULL;
}


//------------------------------------------------------
// The runs that mult_a() reads from successive columns
// are far apart, which hardware prefetchers don't
// foresee. Ask for column j's run, rows i0 to i0 + nr - 1,
// ahead of time.
//------------------------------------------------------
static inline void prefetch_run(int j, int i0, int nr, int n) {
#if defined(__GNUC__) || defined(__clang__)
  if (j < i0) {
    const char *p = (const char *) (g_dist + dist_col_start(j, n) + (i0 - j - 1));
    for (size_t b = 0; b < (size_t) nr * sizeof(double); b += 64)
      __builtin_prefetch(p + b, 0, 0);
  }
#else
  (void) j; (void) i0; (void) nr; (void) n;
#endif
}


//------------------------------------------------------
// Y = A * X, where A = D^2. Each task computes every
// element of GOWER_ROWS rows of Y, reading the parts of
// D's packed lower triangle that those rows need as
// contiguous runs, once for all g_ng groups of vectors.
//------------------------------------------------------
static void *mult_a(void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = g_n;
  int ng        = g_ng;
  int width     = ng * GOWER_W;
  int n_tasks   = (n + GOWER_ROWS - 1) / GOWER_ROWS;

  const double *x = g_x;
  double       *y = g_y;

  // Unrolled loops over a group's GOWER_W vectors, which
  // lets the compiler keep xj[] and s[] in registers.
  #define X_J(v)  xj[v] = xc[v];
  #define AXPY(v) yr[v] += a * xj[v];
  #define DOT(v)  s[v]  += a * xi[v];
  #define BOTH(v) s[v]  += a * xi[v]; yi[v] += a * xj[v];
  #define ZERO(v) s[v]  = 0;
  #define ADD(v)  yc[v] += s[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * GOWER_ROWS;
    int i1 = (i0 + GOWER_ROWS < n) ? i0 + GOWER_ROWS : n;
    int nr = i1 - i0;

    double *yb = y + (size_t) i0 * width;
    memset(yb, 0, (size_t) nr * width * sizeof(double));

    double xj[GOWER_W], s[GOWER_W];

    // Columns left of the block: rows i0..i1-1 of column j
    // are contiguous in the dist vector.
    for (int j = 0; j < i0; j++) {
      const double *d = g_dist + dist_col_start(j, n) + (i0 - j - 1);
      prefetch_run(j + GOWER_AHEAD, i0, nr, n);
      for (int g = 0; g < ng; g++) {
        const double *xc = x + (size_t) j * width + g * GOWER_W;
        GOWER_REP(X_J)
        for (int r = 0; r < nr; r++) {
          double  a  = d[r] * d[r];
          double *yr = yb + (size_t) r * width + g * GOWER_W;
          GOWER_REP(AXPY)
        }
      }
    }

    // The block's own columns: below-diagonal elements
    // within the block feed both of their rows; the rest
    // of each column feeds row j as a dot product. Element
    // (i, j), i > j, is d[i - j - 1].
    for (int j = i0; j < i1; j++) {

      const double *d = g_dist + dist_col_start(j, n);

      for (int g = 0; g < ng; g++) {

        const double *xc = x + (size_t) j * width + g * GOWER_W;
        GOWER_REP(X_J)
        GOWER_REP(ZERO)

        for (int i = j + 1; i < i1; i++) {
          double        a  = d[i - j - 1] * d[i - j - 1];
          const double *xi = x + (size_t) i * width + g * GOWER_W;
          double       *yi = y + (size_t) i * width + g * GOWER_W;
          GOWER_REP(BOTH)
        }

        for (int i = i1; i < n; i++) {
          double        a  = d[i - j - 1] * d[i - j - 1];
          const double *xi = x + (size_t) i * width + g * GOWER_W;
          GOWER_REP(DOT)
        }

        double *yc = y + (size_t) j * width + g * GOWER_W;
        GOWER_REP(ADD)
      }
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
// Sets up `g` for products with the n x n Gower matrix of
// `dist`, up to `max_m` vectors per pass over it (scratch
// from R_alloc), and returns trace(G) = sum_{i<j} d^2 / n.
//------------------------------------------------------
double gower_init(gower_t *g, const double *dist, int n, int max_m, int n_threads) {

  if (max_m < 1) max_m = 1;
  int ng = (max_m + GOWER_W - 1) / GOWER_W;

  g->n         = n;
  g->dist      = dist;
  g->F         = NULL;
  g->f         = 0;
  g->max_m     = ng * GOWER_W;
  g->xr        = (double *) R_alloc((size_t) n * g->max_m, sizeof(double));
  g->yr        = (double *) R_alloc((size_t) n * g->max_m, sizeof(double));
  g->n_threads = threads_for(n_threads, 0.5 * n * (n - 1.0) * GOWER_W, 1e6);

  // Summed per column in parallel, then over columns in order.
  g_n    = n;
  g_dist = dist;
  g_ss   = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, threads_for(n_threads, 0.5 * n * (n - 1.0), 1e6), n);

  double total = 0;
  for (int j = 0; j < n - 1; j++) total += g_ss[j];
  if (!R_FINITE(total))
    error("`dm` must not contain NA, NaN, or infinite values.");

  return total / n;
}


//------------------------------------------------------
// W = (I - F F') * G * V for `m` columns of length n
// (column-major), with V's columns orthogonal to the ones
// vector. A lanczos_op_t: `ctx` is a gower_t.
//------------------------------------------------------
void gower_mult(const double *V, double *W, int m, void *ctx) {

  gower_t *g = (gower_t *) ctx;
  int      n = g->n;

  for (int v0 = 0; v0 < m; v0 += g->max_m) {

    int mc    = (m - v0 < g->max_m) ? m - v0 : g->max_m;
    int ng    = (mc + GOWER_W - 1) / GOWER_W;
    int width = ng * GOWER_W;

    // A column's products don't depend on the columns
    // beside it, so padding with zeros changes nothing.
    for (int i = 0; i < n; i++)
      for (int v = 0; v < width; v++)
        g->xr[(size_t) i * width + v] = (v < mc) ? V[(size_t) (v0 + v) * n + i] : 0;

    g_n    = n;
    g_dist = g->dist;
    g_x    = g->xr;
    g_y    = g->yr;
    g_ng   = ng;
    run_parallel(mult_a, g->n_threads, n);

    for (int v = 0; v < mc; v++) {
      double *w = W + (size_t) (v0 + v) * n;
      for (int i = 0; i < n; i++)
        w[i] = g->yr[(size_t) i * width + v];
      lz_center(w, n);
      for (int i = 0; i < n; i++)
        w[i] *= -0.5;
      lz_project_out(g->F, g->f, w, n);
    }
  }
}
