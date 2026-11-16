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
 * Determinism: everything here runs on the main thread, in
 * a fixed order, and the starting vectors come from a fixed
 * PCG stream. The product function must be deterministic
 * too (see pcoa.c and pca.c).
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

#define dot    lz_dot
#define center lz_center



//------------------------------------------------------
// Uniform random values in [-1, 1) from a PCG stream.
//------------------------------------------------------
static void fill_random(double *x, size_t n, pcg32_random_t *rng) {
  for (size_t i = 0; i < n; i++)
    x[i] = pcg32_random_r(rng) / 2147483648.0 - 1.0;
}


//------------------------------------------------------
// Orthonormalize the `m` columns of R against the `s`
// columns of V, the ones vector (if `deflate`), and each
// other, with two passes of classical Gram-Schmidt. A
// column that is (numerically) in their span is replaced
// by a random one; if that is too, the space is exhausted
// and the column is dropped. Returns how many columns are
// kept, moved to the front of R.
//------------------------------------------------------
static int orthonormalize(
    const double *V, int s, double *R, int m, int n,
    int deflate, pcg32_random_t *rng ) {

  int kept = 0;

  for (int c = 0; c < m; c++) {

    double *r = R + (size_t) kept * n;
    if (c != kept) memcpy(r, R + (size_t) c * n, (size_t) n * sizeof(double));

    for (int attempt = 0; attempt < 2; attempt++) {

      if (attempt) fill_random(r, n, rng);

      double norm0 = sqrt(dot(r, r, n));

      for (int pass = 0; pass < 2; pass++) {
        if (deflate) center(r, n);
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
// The k largest eigenpairs of the n x n operator that
// lz->op applies. See lanczos_t in ecodive.h. B's range
// must have at most lz->dims dimensions (n - 1 or fewer
// with lz->deflate), and k <= dims.
//======================================================
void lanczos(lanczos_t *lz) {

  int    n       = lz->n;
  int    k       = lz->k;
  int    dims    = lz->dims;
  int    deflate = lz->deflate;
  double trace   = lz->scale;

  //--------------------------------------------------
  // Sizes: block p >= k; basis of at most M columns;
  // q Ritz vectors kept at each restart.
  //--------------------------------------------------
  int p = (k > LZ_BLOCK) ? k : LZ_BLOCK;
  int M = LZ_BASIS;
  if (M < 2 * k + 2 * p) M = 2 * k + 2 * p;
  if (p > dims) p = dims;
  if (M > dims) M = dims;
  int q = (M - p + k) / 2;
  if (q < k) q = k;

  double *V   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *W   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *T   = (double *) R_alloc((size_t) n * M, sizeof(double));
  double *R   = (double *) R_alloc((size_t) n * p, sizeof(double));
  double *H   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Hc  = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *Y   = (double *) R_alloc((size_t) M * M, sizeof(double));
  double *tmp = (double *) R_alloc((size_t) 2 * M, sizeof(double));
  double *eig = (double *) R_alloc(M, sizeof(double));
  double *res = lz->resid;
  double *u   = (double *) R_alloc(n, sizeof(double));
  double *bu  = (double *) R_alloc(n, sizeof(double));

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 20260116u, 54u);


  //--------------------------------------------------
  // Starting block.
  //--------------------------------------------------
  fill_random(R, (size_t) n * p, &rng);
  int s_v = orthonormalize(V, 0, R, p, n, deflate, &rng); // columns in V
  int s_w = 0;                                            // ... with W = B * V
  memcpy(V, R, (size_t) n * s_v * sizeof(double));

  int    steps = 0, converged = 0, stalled = 0;
  double best  = INFINITY;

  while (1) {

    //------------------------------------------------
    // Multiply the newest block by B, and extend
    // H = V' * W to the new columns.
    //------------------------------------------------
    int m = s_v - s_w;
    lz->op(V + (size_t) n * s_w, W + (size_t) n * s_w, m, lz->ctx);
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

    converged = worst <= LZ_TOL ||
      (worst <= LZ_TOL_STALL && stalled >= LZ_STALL);

    if (converged || s_w == dims || steps >= LZ_MAX_STEPS) break;

    R_CheckUserInterrupt();


    //------------------------------------------------
    // Next block: the newest products, orthonormalized
    // against the basis.
    //------------------------------------------------
    memcpy(R, W + (size_t) n * (s_w - m), (size_t) n * m * sizeof(double));
    int r = orthonormalize(V, s_v, R, m, n, deflate, &rng);
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
  // The wanted Ritz pairs.
  //--------------------------------------------------
  combine(V, n, s_w, Y, k, lz->vecs);
  memcpy(lz->vals, eig, (size_t) k * sizeof(double));

  lz->steps     = steps;
  lz->converged = converged;
}
