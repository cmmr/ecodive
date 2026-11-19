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
 * For pcoa()'s negative eigenvalue corrections, distances
 * can be read as d + shift (gower_set_shift()), and
 * gower_mult2() computes products with G and with
 * G2 = -1/2 * J * D * J together, from one pass over D.
 * Without a shift, gower_mult() runs exactly as before.
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
static double       *g_mx;     // Per-column largest distance, or NULL.
static double        g_shift;  // Distances are read as d + g_shift.
static double       *g_y2;     // As g_y, for products with D (mult_a_both).


//------------------------------------------------------
// Sum of squares of each column of the lower triangle,
// with distances read as d + g_shift, and (if g_mx isn't
// NULL) each column's largest distance as given.
//------------------------------------------------------
static void *sum_squares(void *arg) {

  int    thread_i  = ((worker_t *) arg)->i;
  int    n_threads = ((worker_t *) arg)->n;
  int    n         = g_n;
  double sh        = g_shift;

  for (int j = thread_i; j < n - 1; j += n_threads) {
    const double *d  = g_dist + dist_col_start(j, n);
    int           nd = n - j - 1;
    double        ss = 0;
    if (sh == 0) {
      for (int i = 0; i < nd; i++)
        ss += d[i] * d[i];
    }
    else {
      for (int i = 0; i < nd; i++) {
        double t = d[i] + sh;
        ss += t * t;
      }
    }
    g_ss[j] = ss;
    if (g_mx) {
      double mx = d[0];
      for (int i = 1; i < nd; i++)
        if (d[i] > mx) mx = d[i];
      g_mx[j] = mx;
    }
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
// As mult_a(), with distances read as t = d + g_shift:
// Y = A * X with A's elements t^2, and, if `both`, also
// Y2 = D * X with D's elements t (to g_y2), from the same
// pass over the distances. For pcoa()'s Cailliez
// correction. (mult_a() itself is left as it is, so that
// products without a shift stay exactly as they were.)
//------------------------------------------------------
static inline void mult_a_gen(void *arg, int both) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n         = g_n;
  int ng        = g_ng;
  int width     = ng * GOWER_W;
  int n_tasks   = (n + GOWER_ROWS - 1) / GOWER_ROWS;
  double sh     = g_shift;

  const double *x  = g_x;
  double       *y  = g_y;
  double       *y2 = g_y2;

  #define X_J(v)  xj[v] = xc[v];
  #define AXPY(v) yr[v] += a * xj[v]; if (both) yr2[v] += t * xj[v];
  #define DOT(v)  s[v]  += a * xi[v]; if (both) s2[v]  += t * xi[v];
  #define BOTH(v) s[v]  += a * xi[v]; yi[v] += a * xj[v]; \
                  if (both) { s2[v] += t * xi[v]; yi2[v] += t * xj[v]; }
  #define ZERO(v) s[v]  = 0; s2[v] = 0;
  #define ADD(v)  yc[v] += s[v]; if (both) yc2[v] += s2[v];

  for (int task = thread_i; task < n_tasks; task += n_threads) {

    int i0 = task * GOWER_ROWS;
    int i1 = (i0 + GOWER_ROWS < n) ? i0 + GOWER_ROWS : n;
    int nr = i1 - i0;

    double *yb  = y + (size_t) i0 * width;
    double *yb2 = both ? y2 + (size_t) i0 * width : NULL;
    memset(yb, 0, (size_t) nr * width * sizeof(double));
    if (both) memset(yb2, 0, (size_t) nr * width * sizeof(double));

    double xj[GOWER_W], s[GOWER_W], s2[GOWER_W];

    for (int j = 0; j < i0; j++) {
      const double *d = g_dist + dist_col_start(j, n) + (i0 - j - 1);
      prefetch_run(j + GOWER_AHEAD, i0, nr, n);
      for (int g = 0; g < ng; g++) {
        const double *xc = x + (size_t) j * width + g * GOWER_W;
        GOWER_REP(X_J)
        for (int r = 0; r < nr; r++) {
          double  t   = d[r] + sh;
          double  a   = t * t;
          double *yr  = yb + (size_t) r * width + g * GOWER_W;
          double *yr2 = both ? yb2 + (size_t) r * width + g * GOWER_W : NULL;
          GOWER_REP(AXPY)
        }
      }
    }

    for (int j = i0; j < i1; j++) {

      const double *d = g_dist + dist_col_start(j, n);

      for (int g = 0; g < ng; g++) {

        const double *xc = x + (size_t) j * width + g * GOWER_W;
        GOWER_REP(X_J)
        GOWER_REP(ZERO)

        for (int i = j + 1; i < i1; i++) {
          double        t   = d[i - j - 1] + sh;
          double        a   = t * t;
          const double *xi  = x + (size_t) i * width + g * GOWER_W;
          double       *yi  = y + (size_t) i * width + g * GOWER_W;
          double       *yi2 = both ? y2 + (size_t) i * width + g * GOWER_W : NULL;
          GOWER_REP(BOTH)
        }

        for (int i = i1; i < n; i++) {
          double        t  = d[i - j - 1] + sh;
          double        a  = t * t;
          const double *xi = x + (size_t) i * width + g * GOWER_W;
          GOWER_REP(DOT)
        }

        double *yc  = y + (size_t) j * width + g * GOWER_W;
        double *yc2 = both ? y2 + (size_t) j * width + g * GOWER_W : NULL;
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
}

static void *mult_a_shift(void *arg) { mult_a_gen(arg, 0); return NULL; }
static void *mult_a_both(void *arg)  { mult_a_gen(arg, 1); return NULL; }


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
  g->n_req     = n_threads;
  g->shift     = 0;
  g->max_d     = 0;
  g->yr2       = NULL;

  // Summed per column in parallel, then over columns in order.
  g_n     = n;
  g_dist  = dist;
  g_shift = 0;
  g_ss    = (double *) R_alloc(n, sizeof(double));
  g_mx    = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, threads_for(n_threads, 0.5 * n * (n - 1.0), 1e6), n);

  double total = 0;
  for (int j = 0; j < n - 1; j++) total += g_ss[j];
  if (!R_FINITE(total))
    error("`dm` must not contain NA, NaN, or infinite values.");

  for (int j = 0; j < n - 1; j++)
    if (g_mx[j] > g->max_d) g->max_d = g_mx[j];
  g_mx = NULL;

  return total / n;
}


//------------------------------------------------------
// From now on, read each distance d (between different
// samples) as d + shift, and return the trace of the
// Gower matrix of those distances.
//------------------------------------------------------
double gower_set_shift(gower_t *g, double shift) {

  int n = g->n;
  g->shift = shift;

  g_n     = n;
  g_dist  = g->dist;
  g_shift = shift;
  g_mx    = NULL;
  g_ss    = (double *) R_alloc(n, sizeof(double));
  run_parallel(sum_squares, threads_for(g->n_threads, 0.5 * n * (n - 1.0), 1e6), n);

  double total = 0;
  for (int j = 0; j < n - 1; j++) total += g_ss[j];

  return total / n;
}


//------------------------------------------------------
// Copies columns v0 to v0 + mc - 1 of V (n x ., column-
// major) to g->xr, row-major and padded with zeros to a
// multiple of GOWER_W columns, and sets the shared state
// for a pass. A column's products don't depend on the
// columns beside it, so padding changes nothing. Returns
// the padded width.
//------------------------------------------------------
static int gower_load(gower_t *g, const double *V, int v0, int mc) {

  int n     = g->n;
  int ng    = (mc + GOWER_W - 1) / GOWER_W;
  int width = ng * GOWER_W;

  for (int i = 0; i < n; i++)
    for (int v = 0; v < width; v++)
      g->xr[(size_t) i * width + v] = (v < mc) ? V[(size_t) (v0 + v) * n + i] : 0;

  g_n     = n;
  g_dist  = g->dist;
  g_x     = g->xr;
  g_y     = g->yr;
  g_y2    = g->yr2;
  g_ng    = ng;
  g_shift = g->shift;

  return width;
}


//------------------------------------------------------
// w = (I - F F') * (-1/2) * J * y for column v of the
// row-major n x width result y.
//------------------------------------------------------
static void gower_finish(const gower_t *g, const double *y, int width, int v, double *w) {
  int n = g->n;
  for (int i = 0; i < n; i++)
    w[i] = y[(size_t) i * width + v];
  lz_center(w, n);
  for (int i = 0; i < n; i++)
    w[i] *= -0.5;
  lz_project_out(g->F, g->f, w, n);
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
    int width = gower_load(g, V, v0, mc);

    // Threads for this pass's work: `width` vectors, not
    // just GOWER_W. (Results are the same for any number.)
    int nt = threads_for(g->n_req, 0.5 * n * (n - 1.0) * width, 1e6);
    if (nt < g->n_threads) nt = g->n_threads;
    run_parallel((g->shift == 0) ? mult_a : mult_a_shift, nt, n);

    for (int v = 0; v < mc; v++)
      gower_finish(g, g->yr, width, v, W + (size_t) (v0 + v) * n);
  }
}


//------------------------------------------------------
// As gower_mult(), but also W2 = (I - F F') * G2 * V,
// where G2 = -1/2 * J * D * J is the Gower matrix of the
// distances themselves rather than their squares (with
// the same shift), from the same passes over `dist`.
//------------------------------------------------------
void gower_mult2(const double *V, double *W, double *W2, int m, gower_t *g) {

  int n = g->n;

  if (!g->yr2)
    g->yr2 = (double *) R_alloc((size_t) n * g->max_m, sizeof(double));

  for (int v0 = 0; v0 < m; v0 += g->max_m) {

    int mc    = (m - v0 < g->max_m) ? m - v0 : g->max_m;
    int width = gower_load(g, V, v0, mc);

    run_parallel(mult_a_both, g->n_threads, n);

    for (int v = 0; v < mc; v++) {
      gower_finish(g, g->yr,  width, v, W  + (size_t) (v0 + v) * n);
      gower_finish(g, g->yr2, width, v, W2 + (size_t) (v0 + v) * n);
    }
  }
}
