// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Block Lanczos with thick restarts, for the eigenpairs
 * with the k largest algebraic eigenvalues of a symmetric,
 * possibly indefinite operator B of size n x n (pcoa()'s B
 * has negative eigenvalues for non-Euclidean distances),
 * given only a function that computes products B * V for
 * a few vectors at a time. Shared by pcoa() and pca().
 *
 * With `deflate` set, every basis vector is kept orthogonal
 * to the ones vector. pcoa() and the sample-side of pca()
 * use it: their B is doubly centered, so the ones vector is
 * in its null space.
 *
 * With `basis` set, every basis vector is also kept
 * orthogonal to its `n_basis` orthonormal columns. dbrda()
 * uses it for the unconstrained axes, which lie in the
 * space orthogonal to the constraints.
 *
 * Optional fields of lanczos_t (see ecodive.h) set a larger
 * basis, a Rayleigh-Ritz step only every few products, a
 * looser tolerance, an early stop once the largest
 * eigenvalue is known to exceed a value, another starting
 * stream, and output of the final Ritz vectors. pcoa()'s
 * negative eigenvalue corrections use them. Left at 0,
 * lanczos() runs exactly as without them.
 *
 * lz_sym_eigen() and lz_gen_eigenvalues() follow the JAMA
 * routines tred2()/tql2() and orthes()/hqr2(), which are in
 * the public domain (JAMA: A Java Matrix Package, by The
 * MathWorks and the National Institute of Standards and
 * Technology), and themselves follow EISPACK.
 *
 * Determinism: the starting vectors come from a fixed PCG
 * stream, and every value is computed in a fixed order.
 * With lanczos_t.n_threads > 1, the steps on blocks of
 * vectors (products of the basis with small matrices, and
 * dot products) are split among threads by rows or by dot
 * product, but each element is still computed by one
 * thread, in the same order as with one: the results are
 * bit-identical for any number of threads. The product
 * function must be deterministic too (see pcoa.c and
 * pca.c).
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON

// Default block size: at least this many vectors, and at
// least k.
#define LZ_BLOCK 4

// Largest Krylov basis kept before a thick restart.
#define LZ_BASIS 60

// Converged when each wanted Ritz pair's residual norm,
// || B u - theta u ||, is at most LZ_TOL * scale, where
// scale is the larger of the caller's `scale` (trace(B))
// and the largest |theta|. Rounding in B * X is on the
// order of DBL_EPSILON * scale, so residuals can stall
// above LZ_TOL. A residual at most LZ_TOL_STALL * scale
// that hasn't halved in LZ_STALL steps is accepted as
// converged.
#define LZ_TOL       1e-14
#define LZ_TOL_STALL 1e-10
#define LZ_STALL     10

// Give up after this many products with blocks of vectors.
#define LZ_MAX_STEPS 10000


//------------------------------------------------------
// Small dense helpers. Vectors have length n; matrices
// are column-major.
//------------------------------------------------------
double lz_dot(const double *x, const double *y, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) s += x[i] * y[i];
  return s;
}

void lz_center(double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) s += x[i];
  s /= n;
  for (int i = 0; i < n; i++) x[i] -= s;
}

#define dot       lz_dot
#define center    lz_center
#define sym_eigen lz_sym_eigen
#define combine   lz_combine
#define fill_random     lz_fill_random
#define orthonormalize  lz_orthonormalize



//------------------------------------------------------
// Uniform random values in [-1, 1) from a PCG stream.
//------------------------------------------------------
void lz_fill_random(double *x, size_t n, pcg32_random_t *rng) {
  for (size_t i = 0; i < n; i++)
    x[i] = pcg32_random_r(rng) / 2147483648.0 - 1.0;
}


//------------------------------------------------------
// Orthonormalize the `m` columns of R against the `s`
// columns of V, the ones vector (if `deflate`), the `f`
// columns of F, and each other, with two passes of
// classical Gram-Schmidt. A
// column that is (numerically) in their span is replaced
// by a random one; if that is too, the space is exhausted
// and the column is dropped. Returns how many columns are
// kept, moved to the front of R.
//------------------------------------------------------
int lz_orthonormalize(
    const double *V, int s, double *R, int m, int n,
    int deflate, const double *F, int f, pcg32_random_t *rng ) {

  int kept = 0;

  for (int c = 0; c < m; c++) {

    double *r = R + (size_t) kept * n;
    if (c != kept) memcpy(r, R + (size_t) c * n, (size_t) n * sizeof(double));

    for (int attempt = 0; attempt < 2; attempt++) {

      if (attempt) fill_random(r, n, rng);

      double norm0 = sqrt(dot(r, r, n));

      for (int pass = 0; pass < 2; pass++) {
        if (deflate) center(r, n);
        for (int j = 0; j < f; j++) {
          const double *q = F + (size_t) j * n;
          double        h = dot(q, r, n);
          for (int i = 0; i < n; i++) r[i] -= h * q[i];
        }
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
// As lz_orthonormalize(), for lanczos()'s blocks of more
// than LZ_BLOCK vectors: block classical Gram-Schmidt,
// applied twice (CGS2, as stable as the two passes of the
// one-column-at-a-time version). All m columns are
// projected off the ones vector, F and V at once (two
// passes), with dot products and products with V split
// among threads (see lz_nt); then off each other, one
// column at a time. Results are the same for any number
// of threads. It runs for blocks wider than LZ_BLOCK
// (k > LZ_BLOCK): which version runs depends only on the
// block, never on the number of threads, and k <= LZ_BLOCK
// runs exactly as before.
//------------------------------------------------------
static void combine_ld(
    const double *V, int n, int s, const double *Y, int ldy, int m,
    double *X, int sub );
static void dots(const double **a, const double **b, double *out, int nd, int n);
static void cross_tn(const double *Q, int nq, const double *R, int m, int n, double *H);

// H, a and b need room for max(s, f, 1) * m values,
// norm0 for m.
static int orthonormalize_block(
    const double *V, int s, double *R, int m, int n,
    int deflate, const double *F, int f, pcg32_random_t *rng,
    double *H, double *norm0, const double **a, const double **b ) {

  // Each column's norm before projection.
  for (int c = 0; c < m; c++) a[c] = b[c] = R + (size_t) c * n;
  dots(a, b, norm0, m, n);
  for (int c = 0; c < m; c++) norm0[c] = sqrt(norm0[c]);

  for (int pass = 0; pass < 2; pass++) {
    if (deflate)
      for (int c = 0; c < m; c++) center(R + (size_t) c * n, n);
    for (int z = 0; z < 2; z++) {
      const double *Q  = z ? V : F;
      int           nq = z ? s : f;
      if (nq == 0 || !Q) continue;
      // H = Q' R (nq x m), then R -= Q H.
      cross_tn(Q, nq, R, m, n, H);
      combine_ld(Q, n, nq, H, nq, m, R, 1);
    }
  }

  // Then off each other.
  int kept = 0;
  for (int c = 0; c < m; c++) {

    double *r = R + (size_t) kept * n;
    if (c != kept) memcpy(r, R + (size_t) c * n, (size_t) n * sizeof(double));

    for (int pass = 0; pass < 2 && kept > 0; pass++) {
      cross_tn(R, kept, r, 1, n, H);
      combine_ld(R, n, kept, H, kept, 1, r, 1);
    }

    double norm = sqrt(dot(r, r, n));

    if (norm > 1e-8 * norm0[c]) {
      for (int i = 0; i < n; i++) r[i] /= norm;
      kept++;
      continue;
    }

    // In the span: a random column instead, projected as
    // lz_orthonormalize() does (rare; on this thread).
    fill_random(r, n, rng);
    double n0 = sqrt(dot(r, r, n));
    for (int pass = 0; pass < 2; pass++) {
      if (deflate) center(r, n);
      for (int j = 0; j < f; j++) {
        const double *q = F + (size_t) j * n;
        double        h = dot(q, r, n);
        for (int i = 0; i < n; i++) r[i] -= h * q[i];
      }
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
    norm = sqrt(dot(r, r, n));
    if (norm > 1e-8 * n0) {
      for (int i = 0; i < n; i++) r[i] /= norm;
      kept++;
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
// `tmp` needs room for 2 * s doubles. Returns 0, or 1 if
// the QL iteration failed to converge for an eigenvalue
// (after 100 iterations; LAPACK reports this as info > 0).
//------------------------------------------------------
int sym_eigen(const double *A, int s, double *vals, double *vecs, double *tmp) {

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
  int    failed = 0;

  for (int l = 0; l < s; l++) {

    // Find a small subdiagonal element.
    double t = fabs(d[l]) + fabs(e[l]);
    if (t > tst1) tst1 = t;
    int m = l;
    while (m < s - 1 && fabs(e[m]) > DBL_EPSILON * tst1) m++;

    // If m == l, d[l] is already an eigenvalue.
    if (m > l) {
      int done = 0;
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

        if (fabs(e[l]) <= DBL_EPSILON * tst1) { done = 1; break; }
      }
      if (!done) failed = 1; // # nocov
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
  return failed;
}


//------------------------------------------------------
// Threads for the steps on blocks of vectors. lanczos()
// sets lz_nt from lanczos_t.n_threads for its run, and
// resets it to 1; the exported helpers called from
// elsewhere run on the calling thread.
//
// Work is split into pieces that don't depend on the
// number of threads: blocks of LZ_ROWS rows, or whole dot
// products. Each element is computed by one thread, in the
// same order as on one thread.
//------------------------------------------------------
#define LZ_ROWS 256

static int lz_nt = 1;

// How many threads for `work` multiply-adds: at least
// ~1e6 each (as threads_for()), else one.
static int lz_threads(double work) {
  return (lz_nt > 1) ? threads_for(lz_nt, work, 1e6) : 1;
}

// X[, c] = V * Y[, c] (sub = 0) or X[, c] -= V * Y[, c]
// (sub = 1), for c < m, for rows r0 .. r1 - 1. Y has
// leading dimension ldy.
static void combine_rows(
    const double *V, int n, int s, const double *Y, int ldy, int m,
    double *X, int sub, int r0, int r1 ) {

  for (int c = 0; c < m; c++) {
    double       *x = X + (size_t) c * n;
    const double *y = Y + (size_t) c * ldy;
    if (!sub) memset(x + r0, 0, (size_t) (r1 - r0) * sizeof(double));
    for (int j = 0; j < s; j++) {
      const double *v = V + (size_t) j * n;
      double        a = y[j];
      if (sub) for (int i = r0; i < r1; i++) x[i] -= a * v[i];
      else     for (int i = r0; i < r1; i++) x[i] += a * v[i];
    }
  }
}

static struct {
  const double *V, *Y;
  double       *X;
  int           n, s, ldy, m, sub;
} cb;

static void *combine_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n_blocks  = (cb.n + LZ_ROWS - 1) / LZ_ROWS;
  for (int b = thread_i; b < n_blocks; b += n_threads) {
    int r0 = b * LZ_ROWS;
    int r1 = (r0 + LZ_ROWS < cb.n) ? r0 + LZ_ROWS : cb.n;
    combine_rows(cb.V, cb.n, cb.s, cb.Y, cb.ldy, cb.m, cb.X, cb.sub, r0, r1);
  }
  return NULL;
}

static void combine_ld(
    const double *V, int n, int s, const double *Y, int ldy, int m,
    double *X, int sub ) {

  if (s == 0 || m == 0) return;
  int nt = lz_threads((double) n * s * m);
  if (nt <= 1) { combine_rows(V, n, s, Y, ldy, m, X, sub, 0, n); return; }

  cb.V = V; cb.n = n; cb.s = s; cb.Y = Y; cb.ldy = ldy; cb.m = m;
  cb.X = X; cb.sub = sub;
  run_parallel(combine_worker, nt, 100);
}


//------------------------------------------------------
// X = V * Y[, 0:m), for V n x s and Y s x s.
//------------------------------------------------------
void combine(const double *V, int n, int s, const double *Y, int m, double *X) {
  combine_ld(V, n, s, Y, s, m, X, 0);
}


//------------------------------------------------------
// out[t] = dot(a[t], b[t], n) for t < nd, each by one
// thread, in lz_dot()'s order.
//------------------------------------------------------
static struct {
  const double **a, **b;
  double        *out;
  int            n, nd;
} dt;

static void *dots_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  for (int t = thread_i; t < dt.nd; t += n_threads)
    dt.out[t] = dot(dt.a[t], dt.b[t], dt.n);
  return NULL;
}

static void dots(const double **a, const double **b, double *out, int nd, int n) {
  int nt = lz_threads((double) n * nd);
  if (nt > nd) nt = nd;
  if (nt <= 1) { for (int t = 0; t < nd; t++) out[t] = dot(a[t], b[t], n); return; }
  dt.a = a; dt.b = b; dt.out = out; dt.n = n; dt.nd = nd;
  run_parallel(dots_worker, nt, 100);
}

// Room for the dot products of one step (set by lanczos()).
static const double **dt_a, **dt_b;
static double        *dt_out;
static int            dt_max = 0;


//------------------------------------------------------
// H = Q' R (nq x m, column-major), for Q n x nq and R
// n x m: each block of LZ_ROWS rows' part by one thread
// (to cx_part), then added up in block order. So the
// result depends only on n, not on the threads, and each
// block of Q is read once for all m columns of R.
//------------------------------------------------------
static double *cx_part = NULL;   // n_blocks x (nq x m), set by lanczos()
static size_t  cx_max  = 0;      //   room per block

static struct {
  const double *Q, *R;
  int           n, nq, m;
} cx;

static void cross_rows(int r0, int r1, double *P) {
  int n = cx.n, nq = cx.nq, m = cx.m, len = r1 - r0;
  for (int c = 0; c < m; c++) {
    const double *rc = cx.R + (size_t) c * n + r0;
    int j = 0;
    for (; j + 3 < nq; j += 4) {
      const double *q0 = cx.Q + (size_t) j * n + r0;
      const double *q1 = q0 + n, *q2 = q1 + n, *q3 = q2 + n;
      double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
      for (int i = 0; i < len; i++) {
        double x = rc[i];
        s0 += q0[i] * x; s1 += q1[i] * x; s2 += q2[i] * x; s3 += q3[i] * x;
      }
      P[(size_t) c * nq + j]     = s0;
      P[(size_t) c * nq + j + 1] = s1;
      P[(size_t) c * nq + j + 2] = s2;
      P[(size_t) c * nq + j + 3] = s3;
    }
    for (; j < nq; j++) {
      const double *q = cx.Q + (size_t) j * n + r0;
      double s0 = 0;
      for (int i = 0; i < len; i++) s0 += q[i] * rc[i];
      P[(size_t) c * nq + j] = s0;
    }
  }
}

static void *cross_worker(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n_blocks  = (cx.n + LZ_ROWS - 1) / LZ_ROWS;
  size_t sz     = (size_t) cx.nq * cx.m;
  for (int b = thread_i; b < n_blocks; b += n_threads) {
    int r0 = b * LZ_ROWS;
    int r1 = (r0 + LZ_ROWS < cx.n) ? r0 + LZ_ROWS : cx.n;
    cross_rows(r0, r1, cx_part + (size_t) b * sz);
  }
  return NULL;
}

static void cross_tn(const double *Q, int nq, const double *R, int m, int n, double *H) {

  size_t sz       = (size_t) nq * m;
  int    n_blocks = (n + LZ_ROWS - 1) / LZ_ROWS;
  if (sz == 0) return;
  if (sz > cx_max) error("lanczos: cross_tn() scratch too small."); // # nocov

  cx.Q = Q; cx.R = R; cx.n = n; cx.nq = nq; cx.m = m;

  int nt = lz_threads((double) n * sz);
  if (nt > n_blocks) nt = n_blocks;
  if (nt > 1) run_parallel(cross_worker, nt, 100);
  else        { worker_t w = { .i = 0, .n = 1 }; cross_worker(&w); }

  memcpy(H, cx_part, sz * sizeof(double));
  for (int b = 1; b < n_blocks; b++) {
    const double *P = cx_part + (size_t) b * sz;
    for (size_t t = 0; t < sz; t++) H[t] += P[t];
  }
}


//------------------------------------------------------
// Extend H = V' * W (M x M, column-major) to columns s_w
// to s_v - 1 of V and W, symmetrizing.
//------------------------------------------------------
void lz_extend_h(
    const double *V, const double *W, double *H,
    int M, int s_w, int s_v, int n ) {

  int nd = (s_v - s_w) * s_v;
  if (lz_nt > 1 && nd <= dt_max) {
    int t = 0;
    for (int b = s_w; b < s_v; b++)
      for (int a = 0; a < s_v; a++, t++) {
        dt_a[t] = V + (size_t) n * a;
        dt_b[t] = W + (size_t) n * b;
      }
    dots(dt_a, dt_b, dt_out, nd, n);
    t = 0;
    for (int b = s_w; b < s_v; b++)
      for (int a = 0; a < s_v; a++, t++)
        H[(size_t) b * M + a] = dt_out[t];
  }
  else {
    for (int b = s_w; b < s_v; b++)
      for (int a = 0; a < s_v; a++)
        H[(size_t) b * M + a] = dot(V + (size_t) n * a, W + (size_t) n * b, n);
  }

  for (int b = s_w; b < s_v; b++) {
    for (int a = 0; a < s_w; a++)
      H[(size_t) a * M + b] = H[(size_t) b * M + a];
    for (int a = s_w; a < b; a++) {
      double avg = (H[(size_t) b * M + a] + H[(size_t) a * M + b]) / 2;
      H[(size_t) b * M + a] = avg;
      H[(size_t) a * M + b] = avg;
    }
  }
}


//------------------------------------------------------
// H = V' * W for the first `keep` columns, symmetrized,
// after a thick restart.
//------------------------------------------------------
void lz_restart_h(const double *V, const double *W, double *H, int M, int keep, int n) {

  // V'W in full, then symmetrized as below.
  if (lz_nt > 1 && keep * keep <= dt_max) {
    int t = 0;
    for (int b = 0; b < keep; b++)
      for (int a = 0; a < keep; a++, t++) {
        dt_a[t] = V + (size_t) n * a;
        dt_b[t] = W + (size_t) n * b;
      }
    dots(dt_a, dt_b, dt_out, keep * keep, n);
    for (int b = 0; b < keep; b++)
      for (int a = 0; a <= b; a++) {
        double h = (a == b) ? dt_out[(size_t) b * keep + a]
          : (dt_out[(size_t) b * keep + a] + dt_out[(size_t) a * keep + b]) / 2;
        H[(size_t) b * M + a] = h;
        H[(size_t) a * M + b] = h;
      }
    return;
  }

  for (int b = 0; b < keep; b++)
    for (int a = 0; a <= b; a++) {
      double h = (a == b)
        ? dot(V + (size_t) n * a, W + (size_t) n * b, n)
        : (dot(V + (size_t) n * a, W + (size_t) n * b, n) +
           dot(V + (size_t) n * b, W + (size_t) n * a, n)) / 2;
      H[(size_t) b * M + a] = h;
      H[(size_t) a * M + b] = h;
    }
}


//------------------------------------------------------
// Rayleigh-Ritz for the basis V (n x s) with W = B * V and
// H = V' * W (leading dimension M): Ritz values to eig
// (decreasing), their coordinates in V to the columns of
// Y, and the first k pairs' residual norms, relative to
// the larger of `trace` and the largest |Ritz value|, to
// res. Returns the largest of those; sets *ql_failed.
//------------------------------------------------------
static double rayleigh_ritz(
    const double *V, const double *W, const double *H, int M, int s, int n,
    int k, double trace, double *Hc, double *eig, double *Y, double *tmp,
    double *u, double *bu, double *res, int *ql_failed ) {

  for (int j = 0; j < s; j++)
    memcpy(Hc + (size_t) j * s, H + (size_t) j * M, (size_t) s * sizeof(double));
  *ql_failed = sym_eigen(Hc, s, eig, Y, tmp);

  double scale = trace;
  for (int j = 0; j < s; j++)
    if (fabs(eig[j]) > scale) scale = fabs(eig[j]);

  double worst = 0;

  if (lz_nt > 1 && k <= dt_max) {
    // All k at once: U = V Y, BU = W Y (u and bu have room
    // for n x k), then each residual's norm.
    combine(V, n, s, Y, k, u);
    combine(W, n, s, Y, k, bu);
    for (int j = 0; j < k; j++) {
      double       *r = bu + (size_t) j * n;
      const double *x = u  + (size_t) j * n;
      for (int i = 0; i < n; i++) r[i] -= eig[j] * x[i];
      dt_a[j] = dt_b[j] = r;
    }
    dots(dt_a, dt_b, dt_out, k, n);
    for (int j = 0; j < k; j++) {
      res[j] = sqrt(dt_out[j]) / scale;
      if (res[j] > worst) worst = res[j];
    }
    return worst;
  }

  for (int j = 0; j < k; j++) {
    combine(V, n, s, Y + (size_t) j * s, 1, u);
    combine(W, n, s, Y + (size_t) j * s, 1, bu);
    for (int i = 0; i < n; i++) bu[i] -= eig[j] * u[i];
    res[j] = sqrt(dot(bu, bu, n)) / scale;
    if (res[j] > worst) worst = res[j];
  }

  return worst;
}


//======================================================
// The k largest eigenpairs of the n x n operator that
// lz->op applies. See lanczos_t in ecodive.h. B's range
// must have at most lz->dims dimensions (n - 1 or fewer
// with lz->deflate, and n_basis fewer again with a
// basis), and k <= dims.
//======================================================
void lanczos(lanczos_t *lz) {

  int    n       = lz->n;
  int    k       = lz->k;
  int    dims    = lz->dims;
  int    deflate = lz->deflate;
  double trace   = lz->scale;

  const double *F = lz->basis;
  int           f = F ? lz->n_basis : 0;

  // Optional settings (see ecodive.h).
  int    every     = (lz->rr_every > 1) ? lz->rr_every : 1;
  double tol       = (lz->tol > 0) ? lz->tol : LZ_TOL;
  double tol_stall = (lz->tol > 0) ? lz->tol : LZ_TOL_STALL;
  int    limit     = (every == 1) ? LZ_STALL : (LZ_STALL + every - 1) / every;
  int    max_steps = (lz->max_steps > 0) ? lz->max_steps : LZ_MAX_STEPS;

  //--------------------------------------------------
  // Sizes: block p >= k; basis of at most M columns;
  // q Ritz vectors kept at each restart.
  //--------------------------------------------------
  int p = (k > LZ_BLOCK) ? k : LZ_BLOCK;
  int M = (lz->max_basis > 0) ? lz->max_basis : LZ_BASIS;
  if (M < 2 * k + 2 * p) M = 2 * k + 2 * p;
  if (p > dims) p = dims;
  if (M > dims) M = dims;
  int q = (M - p + k) / 2;
  if (q < k) q = k;

  double *V   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *T   = (double *) R_alloc((size_t) n * q, sizeof(double)); // restarts keep <= q
  double *R   = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *H   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Hc  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Y   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *tmp = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *eig = (double *) R_alloc(M, sizeof(double));
  double *res = lz->resid;
  // Threads for the steps on blocks (see lz_nt). With more
  // than one, u and bu hold all k Ritz vectors at once.
  lz_nt  = (lz->n_threads > 1) ? lz->n_threads : 1;
  dt_max = 0;
  int ku = (lz_nt > 1) ? k : 1;
  double *u   = (double *) R_alloc((size_t) n * ku, sizeof(double));
  double *bu  = (double *) R_alloc((size_t) n * ku, sizeof(double));
  if (lz_nt > 1) {
    dt_max = M * M;
    dt_a   = (const double **) R_alloc(dt_max, sizeof(double *));
    dt_b   = (const double **) R_alloc(dt_max, sizeof(double *));
    dt_out = (double *)        R_alloc(dt_max, sizeof(double));
  }

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20260116u, lz->stream ? lz->stream : 54u);


  //--------------------------------------------------
  // Starting block.
  //--------------------------------------------------
  // Scratch for orthonormalize_block(), for blocks of up to
  // p columns against up to max(M, f) others.
  size_t  nb_max = (size_t) ((M > f) ? M : f) * p;
  double *Hb     = NULL, *nb0 = NULL;
  const double **ab = NULL, **bb = NULL;
  if (p > LZ_BLOCK) {
    cx_max  = nb_max;
    cx_part = (double *) R_alloc(nb_max * (size_t) ((n + LZ_ROWS - 1) / LZ_ROWS), sizeof(double));
    Hb  = (double *) R_alloc(nb_max, sizeof(double));
    nb0 = (double *) R_alloc(p,      sizeof(double));
    ab  = (const double **) R_alloc(nb_max, sizeof(double *));
    bb  = (const double **) R_alloc(nb_max, sizeof(double *));
  }
  #define ORTHO(Vb, sb, Rb, mb) (((mb) > LZ_BLOCK) \
    ? orthonormalize_block(Vb, sb, Rb, mb, n, deflate, F, f, &rng, Hb, nb0, ab, bb) \
    : orthonormalize(Vb, sb, Rb, mb, n, deflate, F, f, &rng))

  fill_random(R, (size_t) n * p, &rng);
  int s_v = ORTHO(V, 0, R, p); // columns in V
  int s_w = 0;                                            // ... with W = B * V
  memcpy(V, R, (size_t) n * s_v * sizeof(double));

  int    steps = 0, converged = 0, stalled = 0, ql_failed = 0, fresh = 0;
  double best  = INFINITY;

  while (1) {

    //------------------------------------------------
    // Multiply the newest block by B, and extend
    // H = V' * W to the new columns.
    //------------------------------------------------
    int m = s_v - s_w;
    lz->op(V + (size_t) n * s_w, W + (size_t) n * s_w, m, lz->ctx);
    steps++;

    lz_extend_h(V, W, H, M, s_w, s_v, n);
    s_w = s_v;


    //------------------------------------------------
    // Rayleigh-Ritz: every `every` steps, and whenever
    // this may be the last step or a restart is due
    // (which needs current Ritz vectors).
    //------------------------------------------------
    fresh = (every == 1 || steps % every == 0 || s_w == dims ||
             steps >= max_steps || s_v + m > M);

    if (fresh) {

      double worst = rayleigh_ritz(
        V, W, H, M, s_w, n, k, trace, Hc, eig, Y, tmp, u, bu, res, &ql_failed );

      if (worst < best / 2) { best = worst; stalled = 0; }
      else                  { stalled++;                }

      // A basis of the whole space makes Rayleigh-Ritz exact.
      converged = !ql_failed && (s_w == dims || worst <= tol ||
        (worst <= tol_stall && stalled >= limit));

      if (lz->early_stop && !ql_failed && eig[0] > lz->stop_above && res[0] <= lz->stop_tol)
        converged = 1;

      if (converged || s_w == dims || steps >= max_steps) break;
    }

    // pca() holds malloc'd buffers (memory.c): free them first.
    if (user_interrupted()) { lz_nt = 1; dt_max = 0; cx_max = 0; cx_part = NULL; free_all(); error("Interrupted."); } // # nocov


    //------------------------------------------------
    // Next block: the newest products, orthonormalized
    // against the basis.
    //------------------------------------------------
    memcpy(R, W + (size_t) n * (s_w - m), (size_t) n * m * sizeof(double));
    int r = ORTHO(V, s_v, R, m);
    if (r == 0) { // # nocov start: only if the space is exhausted,
      if (!fresh) // which makes Rayleigh-Ritz exact.
        rayleigh_ritz(V, W, H, M, s_w, n, k, trace, Hc, eig, Y, tmp, u, bu, res, &ql_failed);
      converged = !ql_failed;
      break;
    } // # nocov end


    //------------------------------------------------
    // Thick restart: keep the top q Ritz vectors. (Y is
    // current: r <= m, so s_v + m > M above.)
    //------------------------------------------------
    if (s_v + r > M) {

      int keep = (q + r <= M) ? q : M - r;

      combine(V, n, s_w, Y, keep, T);
      memcpy(V, T, (size_t) n * keep * sizeof(double));
      combine(W, n, s_w, Y, keep, T);
      memcpy(W, T, (size_t) n * keep * sizeof(double));

      lz_restart_h(V, W, H, M, keep, n);

      s_v = s_w = keep;
    }

    memcpy(V + (size_t) n * s_v, R, (size_t) n * r * sizeof(double));
    s_v += r;
  }


  //--------------------------------------------------
  // The wanted Ritz pairs, and any other Ritz vectors
  // asked for.
  //--------------------------------------------------
  combine(V, n, s_w, Y, k, lz->vecs);
  memcpy(lz->vals, eig, (size_t) k * sizeof(double));

  lz->n_ritz_out = 0;
  if (lz->ritz && lz->n_ritz > 0) {
    int nr = (lz->n_ritz < s_w) ? lz->n_ritz : s_w;
    combine(V, n, s_w, Y, nr, lz->ritz);
    lz->n_ritz_out = nr;
  }

  lz->steps     = steps;
  lz->converged = converged;

  lz_nt   = 1;
  dt_max  = 0;
  cx_max  = 0;
  cx_part = NULL;
}


//------------------------------------------------------
// Where to read an axis's sign from: the first element
// (lowest index) whose magnitude is within a relative
// 1e-8 of the largest. An axis is flipped so that this
// element is positive. Allowing for rounding keeps a tie
// (e.g., two samples symmetric about the center) from
// being broken by the last bits of a computation, which
// differ between compilers and platforms.
//------------------------------------------------------
int lz_sign_index(const double *x, int n) {
  double amax = 0;
  for (int i = 0; i < n; i++)
    if (fabs(x[i]) > amax) amax = fabs(x[i]);
  double floor_ = amax * (1 - 1e-8);
  for (int i = 0; i < n; i++)
    if (fabs(x[i]) >= floor_) return i;
  return 0; // # nocov
}


//------------------------------------------------------
// x -= F * (F' x), for the f orthonormal columns of F
// (n x f, column-major; nothing if f is 0).
//------------------------------------------------------
void lz_project_out(const double *F, int f, double *x, int n) {
  for (int j = 0; j < f; j++) {
    const double *q = F + (size_t) j * n;
    double        h = lz_dot(q, x, n);
    for (int i = 0; i < n; i++) x[i] -= h * q[i];
  }
}


//------------------------------------------------------
// Eigenvalues of a real general s x s matrix A (column-
// major; overwritten): reduction to upper Hessenberg form
// by Householder reflections, then the Francis double-
// shift QR algorithm. After the public-domain JAMA
// routines orthes() and hqr2(), which follow EISPACK,
// without the eigenvectors. Real parts go to wr and
// imaginary parts to wi. A real eigenvalue found alone,
// or in a 2 x 2 block with real roots, has wi exactly 0,
// but a repeated real eigenvalue can come out as a
// complex pair with a tiny wi. `ort` needs room for s
// doubles. There is no balancing step, so A's entries
// should be of similar magnitude (see pcoa.c). Returns 0,
// or 1 if the QR iteration failed to converge.
//------------------------------------------------------
int lz_gen_eigenvalues(double *A, int s, double *wr, double *wi, double *ort) {

  #define H(i, j) A[(size_t) (j) * s + (i)]

  //--------------------------------------------------
  // Householder reduction to Hessenberg form (orthes).
  //--------------------------------------------------
  int low = 0, high = s - 1;

  for (int m = low + 1; m <= high - 1; m++) {

    double scale = 0;
    for (int i = m; i <= high; i++) scale += fabs(H(i, m - 1));
    if (scale == 0) continue;

    double h = 0;
    for (int i = high; i >= m; i--) {
      ort[i] = H(i, m - 1) / scale;
      h     += ort[i] * ort[i];
    }
    double g = sqrt(h);
    if (ort[m] > 0) g = -g;
    h      = h - ort[m] * g;
    ort[m] = ort[m] - g;

    for (int j = m; j < s; j++) {
      double f = 0;
      for (int i = high; i >= m; i--) f += ort[i] * H(i, j);
      f = f / h;
      for (int i = m; i <= high; i++) H(i, j) -= f * ort[i];
    }

    for (int i = 0; i <= high; i++) {
      double f = 0;
      for (int j = high; j >= m; j--) f += ort[j] * H(i, j);
      f = f / h;
      for (int j = m; j <= high; j++) H(i, j) -= f * ort[j];
    }

    ort[m]      = scale * ort[m];
    H(m, m - 1) = scale * g;
  }

  for (int j = 0; j < s; j++)
    for (int i = j + 2; i < s; i++)
      H(i, j) = 0;


  //--------------------------------------------------
  // Francis double-shift QR (hqr2, eigenvalues only).
  //--------------------------------------------------
  int    nn = s, n = nn - 1, iter = 0, total = 0;
  double eps = DBL_EPSILON, exshift = 0;
  double p = 0, q = 0, r = 0, ss = 0, w = 0, x = 0, y = 0, z = 0;

  double norm = 0;
  for (int i = 0; i < nn; i++)
    for (int j = (i > 0) ? i - 1 : 0; j < nn; j++)
      norm += fabs(H(i, j));

  while (n >= low) {

    // Look for a single small subdiagonal element.
    int l = n;
    while (l > low) {
      ss = fabs(H(l - 1, l - 1)) + fabs(H(l, l));
      if (ss == 0) ss = norm;
      if (fabs(H(l, l - 1)) < eps * ss) break;
      l--;
    }

    if (l == n) { // One root found.

      wr[n] = H(n, n) + exshift;
      wi[n] = 0;
      n--;
      iter = 0;
    }

    else if (l == n - 1) { // Two roots found.

      w = H(n, n - 1) * H(n - 1, n);
      p = (H(n - 1, n - 1) - H(n, n)) / 2.0;
      q = p * p + w;
      z = sqrt(fabs(q));
      x = H(n, n) + exshift;

      if (q >= 0) { // Real pair.
        z         = (p >= 0) ? p + z : p - z;
        wr[n - 1] = x + z;
        wr[n]     = (z != 0) ? x - w / z : x + z;
        wi[n - 1] = 0;
        wi[n]     = 0;
      }
      else {        // Complex pair.
        wr[n - 1] = x + p;
        wr[n]     = x + p;
        wi[n - 1] = z;
        wi[n]     = -z;
      }

      n   -= 2;
      iter = 0;
    }

    else { // No convergence yet.

      if (++total > 50 * nn) return 1; // # nocov

      // Form the shift.
      x = H(n, n);
      y = 0;
      w = 0;
      if (l < n) {
        y = H(n - 1, n - 1);
        w = H(n, n - 1) * H(n - 1, n);
      }

      // Wilkinson's original ad hoc shift.
      if (iter == 10) { // # nocov start
        exshift += x;
        for (int i = low; i <= n; i++) H(i, i) -= x;
        ss = fabs(H(n, n - 1)) + fabs(H(n - 1, n - 2));
        x  = y = 0.75 * ss;
        w  = -0.4375 * ss * ss;
      }

      // MATLAB's ad hoc shift.
      if (iter == 30) {
        ss = (y - x) / 2.0;
        ss = ss * ss + w;
        if (ss > 0) {
          ss = sqrt(ss);
          if (y < x) ss = -ss;
          ss = x - w / ((y - x) / 2.0 + ss);
          for (int i = low; i <= n; i++) H(i, i) -= ss;
          exshift += ss;
          x = y = w = 0.964;
        }
      } // # nocov end

      iter++;

      // Look for two consecutive small subdiagonal elements.
      int m = n - 2;
      while (m >= l) {
        z  = H(m, m);
        r  = x - z;
        ss = y - z;
        p  = (r * ss - w) / H(m + 1, m) + H(m, m + 1);
        q  = H(m + 1, m + 1) - z - r - ss;
        r  = H(m + 2, m + 1);
        ss = fabs(p) + fabs(q) + fabs(r);
        p /= ss;
        q /= ss;
        r /= ss;
        if (m == l) break;
        if (fabs(H(m, m - 1)) * (fabs(q) + fabs(r)) <
            eps * (fabs(p) * (fabs(H(m - 1, m - 1)) + fabs(z) + fabs(H(m + 1, m + 1)))))
          break;
        m--;
      }

      for (int i = m + 2; i <= n; i++) {
        H(i, i - 2) = 0;
        if (i > m + 2) H(i, i - 3) = 0;
      }

      // Double QR step on rows l to n and columns m to n.
      for (int k = m; k <= n - 1; k++) {

        int notlast = (k != n - 1);

        if (k != m) {
          p = H(k, k - 1);
          q = H(k + 1, k - 1);
          r = notlast ? H(k + 2, k - 1) : 0;
          x = fabs(p) + fabs(q) + fabs(r);
          if (x == 0) continue; // # nocov
          p /= x;
          q /= x;
          r /= x;
        }

        ss = sqrt(p * p + q * q + r * r);
        if (p < 0) ss = -ss;
        if (ss == 0) continue; // # nocov

        if (k != m)      H(k, k - 1) = -ss * x;
        else if (l != m) H(k, k - 1) = -H(k, k - 1);

        p = p + ss;
        x = p / ss;
        y = q / ss;
        z = r / ss;
        q = q / p;
        r = r / p;

        // Row modification.
        for (int j = k; j < nn; j++) {
          p = H(k, j) + q * H(k + 1, j);
          if (notlast) {
            p           += r * H(k + 2, j);
            H(k + 2, j) -= p * z;
          }
          H(k, j)     -= p * x;
          H(k + 1, j) -= p * y;
        }

        // Column modification.
        int imax = (n < k + 3) ? n : k + 3;
        for (int i = 0; i <= imax; i++) {
          p = x * H(i, k) + y * H(i, k + 1);
          if (notlast) {
            p           += z * H(i, k + 2);
            H(i, k + 2) -= p * r;
          }
          H(i, k)     -= p;
          H(i, k + 1) -= p * q;
        }
      }
    }
  }

  #undef H

  return 0;
}
