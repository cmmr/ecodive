// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// Ported from gemelli 0.0.13 (gemelli/rpca.py, preprocessing.py,
// matrix_completion.py, and optspace.py), which carries this notice:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2019, Cameron Martino
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived from
//     this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
 * Robust Aitchison PCA (RPCA): Martino et al. (2019), mSystems 4(1),
 * e00016-19, as implemented by gemelli's rpca().
 *
 *   1. Filter the table (gemelli's rpca_table_processing()).
 *   2. Robust clr: log(x) - mean(log(x)) over each sample's non-zero
 *      entries (normalize.c's rclr, on the filtered table). Zeros are
 *      missing values.
 *   3. OptSpace matrix completion (Keshavan, Montanari & Oh 2010) of the
 *      n x m rclr matrix M, with rank r:
 *        - rescale M so that ||M||_F^2 = nnz * r;
 *        - U0, V0: M's top r singular subspaces, scaled by sqrt(n) and
 *          sqrt(m);
 *        - S = argmin_S ||P(U S V' - M)||_F (P keeps observed entries);
 *        - gradient descent on the Grassmann manifold with gemelli's line
 *          search, re-solving S after each step.
 *   4. PCA of the double-centered completed matrix J_n U S V' J_m.
 *
 * gemelli does all of this with dense n x m arrays. Here, nothing of size
 * n x m (or n x n) is ever formed: every step works from the observed
 * entries (CSR and CSC copies) and the n x r and m x r factors.
 *
 *   - The table comes from new_ecomatrix() (ecomatrix.c), so every input
 *     class that pca() takes works, with repeated sparse entries summed
 *     and explicit zeros dropped.
 *   - The initial SVD is lanczos.c's block Lanczos on M'M or MM' (the
 *     smaller side), whose products are two sparse passes. OptSpace's
 *     iterations are equivariant under rotations U -> UQ, V -> VR (S
 *     absorbs them, and every other term depends only on U S V' and on
 *     the rows' norms), so only the subspaces of U0 and V0 matter: any
 *     orthonormal bases of them give the same completed matrix.
 *   - S solves the r^2 x r^2 normal equations A s = c, with
 *       A = sum_k (U_k U_k') (x) (sum_{l in row k} V_l V_l'),
 *     summed over rows, which is gemelli's matrix in singular_values().
 *     Both factors of each term are symmetric, so only their unique
 *     products are accumulated. The minimum-norm solution, with
 *     eigenvalues below 1e-12 times the largest dropped (lz_sym_eigen()),
 *     is lstsq(A, c, rcond = 1e-12)'s.
 *   - The gradients and the distortion are sums over observed entries.
 *   - The line search's distortion along U + tdU, V + tdV is a quartic in
 *     t whose five coefficients take one pass over the observed entries;
 *     each trial step then costs O((n + m) r) instead of a pass.
 *   - J_n U S V' J_m = Qu (Ru S Rv') Qv', with Qu Ru and Qv Rv thin QR
 *     factorizations of the column-centered U and V, so the final SVD is
 *     that of an r x r matrix. (Neither a QR nor an SVD is among the shared
 *     helpers, so thin_qr() and small_svd() are here. An eigensolution of
 *     K'K would lose the small singular values' relative accuracy.)
 *
 * Threads: run_parallel(). Each pass over the observed entries is split
 * into fixed blocks of rows or of columns, which depend only on the data,
 * not on the number of threads. Each block is handled by one thread,
 * writes only its own outputs and its own partial sums, and the main
 * thread adds the partial sums in block order. So results are
 * bit-identical for any number of threads. Everything else (the Lanczos
 * basis, the r x r and r^2 x r^2 problems) runs on the main thread.
 */

#include "ecodive.h"

#include <float.h>  // DBL_EPSILON, DBL_MIN

// normalize.c's code for the robust clr.
#define RPCA_NORM_RCLR  5

// Largest rank. Re-fitting S takes time in proportion to r^4 per sample.
#define RPCA_MAX_K      10

// Blocks of observed entries: at most RPCA_MAX_BLOCKS, each with at least
// RPCA_MIN_BLOCK entries (except the last).
#define RPCA_MAX_BLOCKS 256
#define RPCA_MIN_BLOCK  4096

// The least work (in multiply-adds) worth giving each thread.
#define RPCA_MIN_WORK   65536

// gemelli's OptSpace() defaults, which its rpca() doesn't expose.
static double step_size; // OptSpace's step_size: 10000, from R.
#define RPCA_RES_LIMIT  20
#define RPCA_RCOND      1e-12


//------------------------------------------------------
// Shared with the worker threads.
//------------------------------------------------------
static int     n_rows, n_cols, rnk;     // n samples, m features, rank r.
static int    *r_pos, *r_idx;           // CSR: row k is r_pos[k] .. r_pos[k+1]-1
static double *r_val;
static int    *c_pos, *c_idx;           // CSC: column l likewise.
static double *c_val;
static int     n_rb, n_cb;              // Number of row / column blocks.
static int    *rb, *cb;                 // Block boundaries (rows / columns).

static int           p_w;               // Columns of p_x and p_y.
static const double *p_x;               // Row-major inputs.
static double       *p_y;               // Row-major outputs.
static const double *p_U, *p_V, *p_US, *p_VS, *p_dUS, *p_dV;
static double       *p_part;            // Per-block partial sums.

static const double *p_U_dm;            // Distances between rows of this.
static double       *p_dist;
static int           p_j0, p_j1;        // Columns of the triangle to fill.


//======================================================
// Worker functions. Each processes blocks thread_i,
// thread_i + n_threads, ...
//======================================================

#define FOREACH_BLOCK(nb)                                       \
  int thread_i  = ((worker_t *) arg)->i;                        \
  int n_threads = ((worker_t *) arg)->n;                        \
  for (int blk = thread_i; blk < (nb); blk += n_threads)


//------------------------------------------------------
// p_y (n x w) = M * p_x (m x w).
//------------------------------------------------------
static void *w_mult_rows(void *arg) {
  int w = p_w;
  FOREACH_BLOCK(n_rb) {
    for (int k = rb[blk]; k < rb[blk + 1]; k++) {
      double *y = p_y + (size_t) k * w;
      for (int v = 0; v < w; v++) y[v] = 0;
      for (int e = r_pos[k]; e < r_pos[k + 1]; e++) {
        const double *x = p_x + (size_t) r_idx[e] * w;
        double        a = r_val[e];
        for (int v = 0; v < w; v++) y[v] += a * x[v];
      }
    }
  }
  return NULL;
}


//------------------------------------------------------
// p_y (m x w) = M' * p_x (n x w).
//------------------------------------------------------
static void *w_mult_cols(void *arg) {
  int w = p_w;
  FOREACH_BLOCK(n_cb) {
    for (int l = cb[blk]; l < cb[blk + 1]; l++) {
      double *y = p_y + (size_t) l * w;
      for (int v = 0; v < w; v++) y[v] = 0;
      for (int e = c_pos[l]; e < c_pos[l + 1]; e++) {
        const double *x = p_x + (size_t) c_idx[e] * w;
        double        a = c_val[e];
        for (int v = 0; v < w; v++) y[v] += a * x[v];
      }
    }
  }
  return NULL;
}


//------------------------------------------------------
// With R = P(M - U S V'): p_y (n x r) = R * VS, and
// p_part[blk] = the block's sum of R^2.
//------------------------------------------------------
static void *w_grad_rows(void *arg) {
  int r = rnk;
  FOREACH_BLOCK(n_rb) {
    double ss = 0;
    for (int k = rb[blk]; k < rb[blk + 1]; k++) {
      const double *us = p_US + (size_t) k * r;
      double       *y  = p_y  + (size_t) k * r;
      for (int a = 0; a < r; a++) y[a] = 0;
      for (int e = r_pos[k]; e < r_pos[k + 1]; e++) {
        const double *v  = p_V  + (size_t) r_idx[e] * r;
        const double *vs = p_VS + (size_t) r_idx[e] * r;
        double x = 0;
        for (int a = 0; a < r; a++) x += us[a] * v[a];
        double res = r_val[e] - x;
        ss += res * res;
        for (int a = 0; a < r; a++) y[a] += res * vs[a];
      }
    }
    p_part[blk] = ss;
  }
  return NULL;
}


//------------------------------------------------------
// p_y (m x r) = R' * US.
//------------------------------------------------------
static void *w_grad_cols(void *arg) {
  int r = rnk;
  FOREACH_BLOCK(n_cb) {
    for (int l = cb[blk]; l < cb[blk + 1]; l++) {
      const double *v = p_V + (size_t) l * r;
      double       *y = p_y + (size_t) l * r;
      for (int a = 0; a < r; a++) y[a] = 0;
      for (int e = c_pos[l]; e < c_pos[l + 1]; e++) {
        const double *us = p_US + (size_t) c_idx[e] * r;
        double x = 0;
        for (int a = 0; a < r; a++) x += us[a] * v[a];
        double res = c_val[e] - x;
        for (int a = 0; a < r; a++) y[a] += res * us[a];
      }
    }
  }
  return NULL;
}


//------------------------------------------------------
// Coefficients of the distortion along the search line.
// At (U + t dU) S (V + t dV)', each observed residual is
// e0 + b t + c t^2. p_part[4 blk + 0..3] get the block's
// sums of e0 b, b^2 + 2 e0 c, b c, and c^2.
//------------------------------------------------------
static void *w_line_rows(void *arg) {
  int r = rnk;
  FOREACH_BLOCK(n_rb) {
    double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    for (int k = rb[blk]; k < rb[blk + 1]; k++) {
      const double *us  = p_US  + (size_t) k * r;
      const double *dus = p_dUS + (size_t) k * r;
      for (int e = r_pos[k]; e < r_pos[k + 1]; e++) {
        const double *v  = p_V  + (size_t) r_idx[e] * r;
        const double *dv = p_dV + (size_t) r_idx[e] * r;
        double x = 0, b = 0, c = 0;
        for (int a = 0; a < r; a++) {
          x += us[a]  * v[a];
          b += dus[a] * v[a] + us[a] * dv[a];
          c += dus[a] * dv[a];
        }
        double e0 = x - r_val[e];
        s1 += e0 * b;
        s2 += b * b + 2 * e0 * c;
        s3 += b * c;
        s4 += c * c;
      }
    }
    double *part = p_part + (size_t) blk * 4;
    part[0] = s1; part[1] = s2; part[2] = s3; part[3] = s4;
  }
  return NULL;
}


//------------------------------------------------------
// Normal equations for S, per block of rows:
//   A[(a,b),(i,j)] = sum U_ka U_ki V_lj V_lb
//   C[(a,b)]       = sum U_ka M_kl V_lb
// over observed (k, l), with (a,b) at a * r + b.
//
// A[(a,b),(i,j)] = uu(a,i) * G(b,j), summed over rows k,
// with uu = U_k U_k' and G = sum_l V_l V_l', both
// symmetric. So only T[(a,i),(b,j)] for a <= i, b <= j
// is accumulated (np x np, np = r (r + 1) / 2); A is
// expanded from it after the blocks are added up. Each
// block's T, then C (r^2), go to p_part.
//------------------------------------------------------
#define RPCA_NP (RPCA_MAX_K * (RPCA_MAX_K + 1) / 2)

static void *w_sv_rows(void *arg) {

  int    r  = rnk, np = r * (r + 1) / 2;
  size_t sz = (size_t) np * np + (size_t) r * r;
  double G[RPCA_NP], c[RPCA_MAX_K], uu[RPCA_NP];

  FOREACH_BLOCK(n_rb) {

    double *T = p_part + (size_t) blk * sz;
    double *C = T + (size_t) np * np;
    memset(T, 0, sz * sizeof(double));

    for (int k = rb[blk]; k < rb[blk + 1]; k++) {

      const double *u = p_U + (size_t) k * r;

      for (int i = 0; i < np; i++) G[i] = 0;
      for (int i = 0; i < r;  i++) c[i] = 0;

      // G(b,j), b <= j, at j * (j + 1) / 2 + b.
      for (int e = r_pos[k]; e < r_pos[k + 1]; e++) {
        const double *v = p_V + (size_t) r_idx[e] * r;
        double        x = r_val[e];
        for (int j = 0, q = 0; j < r; j++) {
          c[j] += x * v[j];
          for (int b = 0; b <= j; b++, q++)
            G[q] += v[j] * v[b];
        }
      }

      for (int i = 0, q = 0; i < r; i++)
        for (int a = 0; a <= i; a++, q++)
          uu[q] = u[a] * u[i];

      for (int a = 0; a < r; a++)
        for (int b = 0; b < r; b++)
          C[a * r + b] += u[a] * c[b];

      for (int p = 0; p < np; p++) {
        double *Trow = T + (size_t) p * np;
        double  x    = uu[p];
        for (int q = 0; q < np; q++)
          Trow[q] += x * G[q];
      }
    }
  }
  return NULL;
}


//------------------------------------------------------
// Euclidean distances between the rows of p_U_dm
// (n x r), into columns p_j0 .. p_j1-1 of a dist
// object's packed lower triangle. Each column is
// written by one thread.
//------------------------------------------------------
static void *w_dist(void *arg) {
  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;
  int n = n_rows, r = rnk;
  for (int j = p_j0 + thread_i; j < p_j1; j += n_threads) {
    double       *d  = p_dist + dist_col_start(j, n);
    const double *uj = p_U_dm + (size_t) j * r;
    for (int i = j + 1; i < n; i++) {
      const double *ui = p_U_dm + (size_t) i * r;
      double ss = 0;
      for (int a = 0; a < r; a++) {
        double t = ui[a] - uj[a];
        ss += t * t;
      }
      d[i - j - 1] = sqrt(ss);
    }
  }
  return NULL;
}


//======================================================
// Running passes over the observed entries.
//======================================================

static int nt_rows, nt_cols; // Threads for row / column passes.
static int n_tasks;          // Observed entries, for run_parallel().

static void run_rows(pthread_func_t f) { run_parallel(f, nt_rows, n_tasks); }
static void run_cols(pthread_func_t f) { run_parallel(f, nt_cols, n_tasks); }


//------------------------------------------------------
// Split `count` rows (or columns), with entry offsets
// `pos`, into blocks of at least `target` entries (the
// last can be smaller). Returns the number of blocks;
// their boundaries go to `bounds`.
//------------------------------------------------------
static int make_blocks(const int *pos, int count, int target, int *bounds) {
  int nb = 0;
  bounds[0] = 0;
  for (int i = 0; i < count; ) {
    int start = pos[i];
    while (i < count && pos[i] - start < target) i++;
    bounds[++nb] = i;
  }
  return nb;
}


//======================================================
// Small dense helpers. All run on the main thread.
//======================================================

//------------------------------------------------------
// Thin QR of A (n x r, column-major, overwritten by Q)
// by classical Gram-Schmidt with two passes. R (r x r,
// column-major) is upper triangular. A column with
// nothing left after projection becomes zero, with
// R[j,j] = 0; A = QR still holds.
//------------------------------------------------------
static void thin_qr(double *A, int n, int r, double *R) {

  for (int i = 0; i < r * r; i++) R[i] = 0;

  for (int j = 0; j < r; j++) {
    double *x = A + (size_t) j * n;
    for (int pass = 0; pass < 2; pass++)
      for (int i = 0; i < j; i++) {
        const double *q = A + (size_t) i * n;
        double        h = lz_dot(q, x, n);
        for (int t = 0; t < n; t++) x[t] -= h * q[t];
        R[j * r + i] += h;
      }
    double norm = sqrt(lz_dot(x, x, n));
    R[j * r + j] = norm;
    for (int t = 0; t < n; t++) x[t] = (norm > 0) ? x[t] / norm : 0;
  }
}


//------------------------------------------------------
// The rotation that zeros the off-diagonal element of
// [[alpha, gamma], [gamma, beta]]: tangent t, cosine c,
// sine s, as in Golub & Van Loan's sym.schur2.
//------------------------------------------------------
static void jacobi_rotation(double alpha, double beta, double gamma,
                            double *c, double *s, double *t) {
  double tau = (beta - alpha) / (2 * gamma);
  double tt  = (fabs(tau) > 1e150)
    ? 0.5 / tau
    : ((tau >= 0) ? 1 : -1) / (fabs(tau) + sqrt(1 + tau * tau));
  *t = tt;
  *c = 1 / sqrt(1 + tt * tt);
  *s = tt * *c;
}


//------------------------------------------------------
// SVD of an r x r matrix K (column-major) by one-sided
// Jacobi (Hestenes): K Z = G, with Z orthogonal and G's
// columns orthogonal. G overwrites K. Columns are then
// sorted by decreasing norm (the singular values).
//------------------------------------------------------
static void small_svd(double *K, int r, double *Z, double *sv) {

  #define G_(i, j) K[(size_t) (j) * r + (i)]
  #define Z_(i, j) Z[(size_t) (j) * r + (i)]

  for (int j = 0; j < r; j++)
    for (int i = 0; i < r; i++)
      Z_(i, j) = (i == j);

  for (int sweep = 0; sweep < 100; sweep++) {

    int rotated = 0;

    for (int p = 0; p < r - 1; p++)
      for (int q = p + 1; q < r; q++) {

        double alpha = 0, beta = 0, gamma = 0;
        for (int i = 0; i < r; i++) {
          alpha += G_(i, p) * G_(i, p);
          beta  += G_(i, q) * G_(i, q);
          gamma += G_(i, p) * G_(i, q);
        }

        if (fabs(gamma) < DBL_MIN ||
            fabs(gamma) <= DBL_EPSILON * sqrt(alpha) * sqrt(beta))
          continue;

        double c, sn, t;
        jacobi_rotation(alpha, beta, gamma, &c, &sn, &t);

        for (int i = 0; i < r; i++) {
          double gp = G_(i, p), gq = G_(i, q);
          G_(i, p) = c * gp - sn * gq;
          G_(i, q) = sn * gp + c * gq;
          double zp = Z_(i, p), zq = Z_(i, q);
          Z_(i, p) = c * zp - sn * zq;
          Z_(i, q) = sn * zp + c * zq;
        }

        rotated = 1;
      }

    if (!rotated) break;
  }

  for (int j = 0; j < r; j++) {
    double ss = 0;
    for (int i = 0; i < r; i++) ss += G_(i, j) * G_(i, j);
    sv[j] = sqrt(ss);
  }

  for (int i = 0; i < r - 1; i++) {
    int k = i;
    for (int j = i + 1; j < r; j++)
      if (sv[j] > sv[k]) k = j;
    if (k != i) {
      double t = sv[k]; sv[k] = sv[i]; sv[i] = t;
      for (int row = 0; row < r; row++) {
        t = G_(row, i); G_(row, i) = G_(row, k); G_(row, k) = t;
        t = Z_(row, i); Z_(row, i) = Z_(row, k); Z_(row, k) = t;
      }
    }
  }

  #undef G_
  #undef Z_
}


//======================================================
// Initial SVD: lanczos.c on B = M'M (when m <= n) or
// B = MM' (otherwise), whose dimension d is the smaller
// of m and n. Each product B X is two passes over the
// observed entries, RPCA_MAX_K vectors at a time.
//======================================================

typedef struct {
  int     right;          // 1: B = M'M (d = m); 0: B = MM' (d = n).
  double *xr, *tr, *yr;   // Row-major scratch, RPCA_MAX_K wide.
} op_ctx_t;

static void op_b(const double *X, double *Y, int w, void *ctx) {

  op_ctx_t *c = (op_ctx_t *) ctx;
  int       d = c->right ? n_cols : n_rows;

  for (int v0 = 0; v0 < w; v0 += RPCA_MAX_K) {

    int           mc = (w - v0 < RPCA_MAX_K) ? w - v0 : RPCA_MAX_K;
    const double *Xc = X + (size_t) v0 * d;
    double       *Yc = Y + (size_t) v0 * d;

    for (int i = 0; i < d; i++)
      for (int v = 0; v < mc; v++)
        c->xr[(size_t) i * mc + v] = Xc[(size_t) v * d + i];

    p_w = mc;

    if (c->right) {
      p_x = c->xr; p_y = c->tr; run_rows(w_mult_rows);
      p_x = c->tr; p_y = c->yr; run_cols(w_mult_cols);
    }
    else {
      p_x = c->xr; p_y = c->tr; run_cols(w_mult_cols);
      p_x = c->tr; p_y = c->yr; run_rows(w_mult_rows);
    }

    for (int i = 0; i < d; i++)
      for (int v = 0; v < mc; v++)
        Yc[(size_t) v * d + i] = c->yr[(size_t) i * mc + v];
  }
}


//======================================================
// OptSpace pieces.
//======================================================

// Work space for solve_s().
static double *sv_part, *sv_T, *sv_A, *sv_vals, *sv_vecs, *sv_tmp;

//------------------------------------------------------
// S = argmin ||P(U S V' - M)||_F, as gemelli's
// singular_values(): the minimum-norm solution of the
// normal equations, dropping eigenvalues of A below
// RPCA_RCOND times the largest.
//------------------------------------------------------
static void solve_s(const double *U, const double *V, double *S) {

  int    r = rnk, r2 = r * r, np = r * (r + 1) / 2;
  size_t nt = (size_t) np * np, sz = nt + r2;
  double C[RPCA_MAX_K * RPCA_MAX_K], x[RPCA_MAX_K * RPCA_MAX_K];

  p_U = U; p_V = V; p_part = sv_part;
  run_rows(w_sv_rows);

  double *T = sv_T, *A = sv_A, *vals = sv_vals, *vecs = sv_vecs;

  for (size_t i = 0; i < nt; i++) T[i] = 0;
  for (int i = 0; i < r2; i++) C[i] = 0;
  for (int blk = 0; blk < n_rb; blk++) {
    const double *pp = sv_part + (size_t) blk * sz;
    for (size_t i = 0; i < nt; i++) T[i] += pp[i];
    for (int i = 0; i < r2; i++) C[i] += pp[nt + i];
  }

  // A[(a,b),(i,j)] = T[(min(a,i), max(a,i)), (min(b,j), max(b,j))],
  // with pair (x, y), x <= y, at y * (y + 1) / 2 + x. Symmetric, so
  // row-major and column-major are the same.
  #define PAIR(x, y) (((x) <= (y)) ? (y) * ((y) + 1) / 2 + (x) : (x) * ((x) + 1) / 2 + (y))
  for (int a = 0; a < r; a++)
    for (int b = 0; b < r; b++)
      for (int i = 0; i < r; i++)
        for (int j = 0; j < r; j++)
          A[(size_t) (a * r + b) * r2 + (i * r + j)] =
            T[(size_t) PAIR(a, i) * np + PAIR(b, j)];
  #undef PAIR

  // A QL failure (after 100 iterations for one eigenvalue) leaves
  // usable, if less accurate, eigenpairs.
  (void) lz_sym_eigen(A, r2, vals, vecs, sv_tmp);

  double lim = 0;
  for (int i = 0; i < r2; i++)
    if (fabs(vals[i]) > lim) lim = fabs(vals[i]);
  lim *= RPCA_RCOND;

  for (int i = 0; i < r2; i++) x[i] = 0;
  for (int e = 0; e < r2; e++) {
    if (!(fabs(vals[e]) > lim)) continue;
    const double *q = vecs + (size_t) e * r2;
    double        h = lz_dot(q, C, r2) / vals[e];
    for (int i = 0; i < r2; i++) x[i] += h * q[i];
  }

  for (int i = 0; i < r2; i++) S[i] = x[i]; // S[i * r + j] = S_ij
}

// US = U S (rows x r), for U row-major.
static void mult_s(const double *U, const double *S, int rows, int transpose, double *US) {
  int r = rnk;
  for (int k = 0; k < rows; k++) {
    const double *u = U  + (size_t) k * r;
    double       *y = US + (size_t) k * r;
    for (int j = 0; j < r; j++) {
      double s = 0;
      for (int i = 0; i < r; i++)
        s += u[i] * (transpose ? S[j * r + i] : S[i * r + j]);
      y[j] = s;
    }
  }
}

//------------------------------------------------------
// gemelli's grassmann_manifold_one(): the penalty on
// rows of U (rows x r) whose squared norm is at least
// 2 * step_size * r, at U + t dU.
//------------------------------------------------------
static double penalty(const double *U, const double *dU, double t, int rows) {
  int    r     = rnk;
  double denom = 2 * step_size * r, sum = 0;
  for (int k = 0; k < rows; k++) {
    double ss = 0;
    for (int a = 0; a < r; a++) {
      double x = U[(size_t) k * r + a] + t * dU[(size_t) k * r + a];
      ss += x * x;
    }
    double step = ss / denom;
    if (step < 1) continue;
    double g = exp((step - 1) * (step - 1)) - 1;
    if (g != INFINITY) sum += g;
  }
  return sum;
}

//------------------------------------------------------
// One factor's descent direction, as gemelli's
// gradient_decent():
//   dU = -R VS + U Qu + rho * grassmann_manifold_two(U),
//   Qu = U' R VS / rows.
// RVS (rows x r) holds R VS on entry.
//------------------------------------------------------
static void direction(const double *U, const double *RVS, int rows,
                      double rho, double *dU) {
  int    r = rnk;
  double Q[RPCA_MAX_K * RPCA_MAX_K];

  for (int a = 0; a < r; a++)
    for (int b = 0; b < r; b++) {
      double s = 0;
      for (int k = 0; k < rows; k++)
        s += U[(size_t) k * r + a] * RVS[(size_t) k * r + b];
      Q[a * r + b] = s / rows;
    }

  double denom = 2 * step_size * r;

  for (int k = 0; k < rows; k++) {
    const double *u  = U   + (size_t) k * r;
    const double *g  = RVS + (size_t) k * r;
    double       *du = dU  + (size_t) k * r;

    double step = lz_dot(u, u, r) / denom;
    double reg  = 2 * exp((step - 1) * (step - 1)) * (step - 1);
    if (reg < 0) reg = 0;
    reg /= step_size * r;

    for (int b = 0; b < r; b++) {
      double s = 0;
      for (int a = 0; a < r; a++) s += u[a] * Q[a * r + b];
      du[b] = -g[b] + s + rho * u[b] * reg;
    }
  }
}


//======================================================
// Filters, as gemelli's rpca_table_processing():
//   features: total > mfc, and present (> 0) in more
//     than mff percent of samples;
//   samples: total over the remaining features > msc;
//   then features left without any non-zero count.
// With `drop_single`, also samples with a single
// non-zero feature, whose rclr carries no information
// (it is 0); they go to `single` (0-based), and their
// number to *n_single.
//
// The table is new_ecomatrix()'s: N0 samples, P0
// features, sample i's entries at tpos[i] ..
// tpos[i+1]-1. On return, s_new[i] (f_new[j]) is the
// sample's (feature's) index among those kept, or -1;
// *n and *m are the numbers kept. Scratch is R_alloc'd.
//======================================================
void rpca_filter(
    int N0, int P0, const int *tpos, const int *tidx, const double *tval,
    double msc, double mfc, double mff, int drop_single,
    int *s_new, int *f_new, int *single, int *n_single, int *n, int *m ) {

  double *fsum = (double *) R_alloc(P0 > 0 ? P0 : 1, sizeof(double));
  int    *fcnt = (int *)    R_alloc(P0 > 0 ? P0 : 1, sizeof(int));

  for (int j = 0; j < P0; j++) { fsum[j] = 0; fcnt[j] = 0; }
  for (int e = 0; e < tpos[N0]; e++) {
    fsum[tidx[e]] += tval[e];
    fcnt[tidx[e]]++;
  }

  for (int j = 0; j < P0; j++)
    f_new[j] = (fsum[j] > mfc && (double) fcnt[j] / N0 > mff / 100) ? 0 : -1;

  int ns = 0;
  for (int i = 0; i < N0; i++) {
    double ss = 0;
    int    nz = 0;
    for (int e = tpos[i]; e < tpos[i + 1]; e++)
      if (f_new[tidx[e]] == 0) {
        ss += tval[e];
        nz++;
      }
    s_new[i] = (ss > msc && (nz >= 2 || !drop_single)) ? 0 : -1;
    if (drop_single && ss > msc && nz == 1) single[ns++] = i;
  }

  // Features present in a kept sample.
  for (int j = 0; j < P0; j++) fcnt[j] = 0;
  for (int i = 0; i < N0; i++)
    if (s_new[i] == 0)
      for (int e = tpos[i]; e < tpos[i + 1]; e++)
        fcnt[tidx[e]] = 1;

  int nn = 0, mm = 0;
  for (int j = 0; j < P0; j++) f_new[j] = (f_new[j] == 0 && fcnt[j]) ? mm++ : -1;
  for (int i = 0; i < N0; i++) s_new[i] = (s_new[i] == 0) ? nn++ : -1;

  *n = nn; *m = mm; *n_single = ns;
}


//======================================================
// A copy of the named list `x` with `k` more elements,
// `vals`, named `names`.
//======================================================
SEXP rpca_list_append(SEXP x, const char **names, SEXP *vals, int k) {

  int  n0  = length(x);
  SEXP res = PROTECT(allocVector(VECSXP, n0 + k));
  SEXP nms = PROTECT(allocVector(STRSXP, n0 + k));
  SEXP old = getAttrib(x, R_NamesSymbol);

  for (int i = 0; i < n0; i++) {
    SET_VECTOR_ELT(res, i, VECTOR_ELT(x, i));
    SET_STRING_ELT(nms, i, STRING_ELT(old, i));
  }
  for (int i = 0; i < k; i++) {
    SET_VECTOR_ELT(res, n0 + i, vals[i]);
    SET_STRING_ELT(nms, n0 + i, mkChar(names[i]));
  }
  setAttrib(res, R_NamesSymbol, nms);

  UNPROTECT(2);
  return res;
}


//======================================================
// OptSpace and PCA (steps 3 and 4 above) of an n x m
// table of robust clr values, given as CSR: row k's
// observed entries are rpos[k] .. rpos[k+1]-1, at
// columns ridx[] (increasing) with values rval[]. All
// three are R_alloc'd by the caller; rval is rescaled
// in place. Rows and columns may have no observed
// entries (rpca_phylo()'s all-missing columns and
// rows): m counts every column, and they all take part
// in OptSpace's scaling, its updates of U and V, and
// the final centering, as in gemelli. r is the rank,
// at most min(n, m) - 1. `no_var` is the error message
// for a table whose observed values are all 0.
//
// `step` is OptSpace's step_size (gemelli: 10000;
// tests use a smaller value to exercise its penalty).
//
// Returns list(scores, loadings, eig, k, cost, n_iters,
// dist, u0, v0, u, s, v, ss, lanczos_steps,
// lanczos_converged, threads): scores n x r (u * s),
// loadings m x r (v), eig the r values s^2 (largest
// first), the rank r, the fit error after each
// iteration, and the robust clr values' sum of squares.
// dist, u0, v0, u, s and v are NULL unless asked for.
//======================================================
SEXP rpca_core(
    int n, int m, int *rpos, int *ridx, double *rval,
    int r, int max_iter, double tol, int want_dm, int debug,
    double step, int n_threads, const char *no_var ) {

  if (r < 1 || r > RPCA_MAX_K || r > n - 1 || r > m - 1 || max_iter < 1 ||
      n_threads < 1 || !(step > 0))
    error("Invalid arguments to rpca_core()."); // # nocov

  step_size    = step;
  int    nnz   = rpos[n];
  size_t nnz_z = (size_t) nnz;

  n_rows = n; n_cols = m; rnk = r;
  r_pos = rpos; r_idx = ridx; r_val = rval;

  double ss = 0;
  for (int e = 0; e < nnz; e++) ss += rval[e] * rval[e];

  if (!(ss > 0)) error("%s", no_var);


  //--------------------------------------------------
  // Blocks, and threads: enough work for each, and no
  // more than there are blocks.
  //--------------------------------------------------
  size_t target_z = (nnz_z + RPCA_MAX_BLOCKS - 2) / (RPCA_MAX_BLOCKS - 1);
  int    target   = (target_z < RPCA_MIN_BLOCK) ? RPCA_MIN_BLOCK : (int) target_z;

  rb   = (int *) R_alloc((size_t) n + 1, sizeof(int));
  n_rb = make_blocks(rpos, n, target, rb);

  int *cpos = (int *) R_alloc((size_t) m + 1, sizeof(int));
  for (int l = 0; l <= m; l++) cpos[l] = 0;
  for (int e = 0; e < nnz; e++) cpos[ridx[e] + 1]++;
  for (int l = 0; l < m; l++) cpos[l + 1] += cpos[l];
  cb   = (int *) R_alloc((size_t) m + 1, sizeof(int));
  n_cb = make_blocks(cpos, m, target, cb);

  int nt  = threads_for(n_threads, (double) nnz * (r + 1), RPCA_MIN_WORK);
  nt_rows = (n_rb < nt) ? n_rb : nt;
  nt_cols = (n_cb < nt) ? n_cb : nt;
  n_tasks = nnz;
  int nt_max = (nt_rows > nt_cols) ? nt_rows : nt_cols;

  int nprot = 0;

  double rescal = sqrt((double) nnz * r / ss);
  double eps    = nnz / sqrt((double) m * n);
  double rho    = eps * n;

  for (int e = 0; e < nnz; e++) rval[e] *= rescal;


  //--------------------------------------------------
  // CSC copy, rows in increasing order in each column.
  //--------------------------------------------------
  int    *cidx = (int *)    R_alloc(nnz_z, sizeof(int));
  double *cval = (double *) R_alloc(nnz_z, sizeof(double));
  int    *fill = (int *)    R_alloc(m, sizeof(int));

  for (int l = 0; l < m; l++) fill[l] = cpos[l];
  for (int k = 0; k < n; k++)
    for (int e = rpos[k]; e < rpos[k + 1]; e++) {
      int at = fill[ridx[e]]++;
      cidx[at] = k;
      cval[at] = rval[e];
    }

  c_pos = cpos; c_idx = cidx; c_val = cval;


  //--------------------------------------------------
  // Work space for solve_s(), the gradient, and the
  // line search.
  //--------------------------------------------------
  int    np  = r * (r + 1) / 2, r2 = r * r;
  size_t part_len = (size_t) n_rb * ((size_t) np * np + (size_t) r2);
  if (part_len < (size_t) n_rb * 4) part_len = (size_t) n_rb * 4;
  sv_part = (double *) R_alloc(part_len, sizeof(double));
  sv_T    = (double *) R_alloc((size_t) np * np, sizeof(double));
  sv_A    = (double *) R_alloc((size_t) r2 * r2, sizeof(double));
  sv_vecs = (double *) R_alloc((size_t) r2 * r2, sizeof(double));
  sv_vals = (double *) R_alloc(r2, sizeof(double));
  sv_tmp  = (double *) R_alloc((size_t) 2 * r2, sizeof(double));


  //--------------------------------------------------
  // Initial factors: bases of the top singular
  // subspaces, scaled by sqrt(n) and sqrt(m).
  //--------------------------------------------------
  size_t nr_ = (size_t) n * r, mr_ = (size_t) m * r;

  double *U    = (double *) R_alloc(nr_, sizeof(double)); // row-major
  double *V    = (double *) R_alloc(mr_, sizeof(double));
  double *US   = (double *) R_alloc(nr_, sizeof(double));
  double *VS   = (double *) R_alloc(mr_, sizeof(double));
  double *RVS  = (double *) R_alloc(nr_, sizeof(double));
  double *RUS  = (double *) R_alloc(mr_, sizeof(double));
  double *dU   = (double *) R_alloc(nr_, sizeof(double));
  double *dV   = (double *) R_alloc(mr_, sizeof(double));
  double *dUS  = (double *) R_alloc(nr_, sizeof(double));
  double  S[RPCA_MAX_K * RPCA_MAX_K];

  int lz_conv, lz_steps;
  {
    int right = (m <= n);
    int d     = right ? m : n;
    int other = right ? n : m;

    op_ctx_t octx = {
      .right = right,
      .xr    = (double *) R_alloc((size_t) d     * RPCA_MAX_K, sizeof(double)),
      .tr    = (double *) R_alloc((size_t) other * RPCA_MAX_K, sizeof(double)),
      .yr    = (double *) R_alloc((size_t) d     * RPCA_MAX_K, sizeof(double)) };

    double *top = (double *) R_alloc((size_t) d * r, sizeof(double));
    double *oth = (double *) R_alloc((size_t) other * r, sizeof(double));

    // trace(B) = ||M||_F^2 = nnz * r after the rescaling. A basis of
    // at least 8 blocks (of max(r, 4) vectors): for larger r, fewer
    // restarts save more products with B (each two passes over the
    // observed entries) than the larger basis costs.
    int basis = 8 * ((r > 4) ? r : 4);
    lanczos_t lz = {
      .n = d, .k = r, .dims = d, .deflate = 0, .scale = (double) nnz * r,
      .max_basis = (basis > 60) ? basis : 0,
      .op = op_b, .ctx = &octx, .n_threads = n_threads, .vecs = top,
      .vals  = (double *) R_alloc(r, sizeof(double)),
      .resid = (double *) R_alloc(r, sizeof(double)) };
    lanczos(&lz);
    lz_conv  = lz.converged;
    lz_steps = lz.steps;

    // The other side: M Y or M' Y, orthonormalized.
    for (int i = 0; i < d; i++)
      for (int a = 0; a < r; a++)
        octx.xr[(size_t) i * r + a] = top[(size_t) a * d + i];
    p_w = r; p_x = octx.xr; p_y = octx.tr;
    if (right) run_rows(w_mult_rows); else run_cols(w_mult_cols);
    for (int i = 0; i < other; i++)
      for (int a = 0; a < r; a++)
        oth[(size_t) a * other + i] = octx.tr[(size_t) i * r + a];

    // Keeps all r columns: other > r, so the random
    // replacement for a dependent column fits.
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, 20261117u, 52u);
    int kept = lz_orthonormalize(NULL, 0, oth, r, other, 0, NULL, 0, &rng);
    if (kept < r)
      error("Robust PCA could not initialize OptSpace."); // # nocov

    const double *cu = right ? oth : top;
    const double *cv = right ? top : oth;
    double        fu = sqrt((double) n), fv = sqrt((double) m);
    for (int k = 0; k < n; k++)
      for (int a = 0; a < r; a++)
        U[(size_t) k * r + a] = cu[(size_t) a * n + k] * fu;
    for (int l = 0; l < m; l++)
      for (int a = 0; a < r; a++)
        V[(size_t) l * r + a] = cv[(size_t) a * m + l] * fv;
  }

  SEXP sexp_u0 = R_NilValue, sexp_v0 = R_NilValue;
  if (debug) {
    sexp_u0 = PROTECT(allocMatrix(REALSXP, n, r)); nprot++;
    sexp_v0 = PROTECT(allocMatrix(REALSXP, m, r)); nprot++;
    for (int k = 0; k < n; k++) for (int a = 0; a < r; a++)
      REAL(sexp_u0)[(size_t) a * n + k] = U[(size_t) k * r + a];
    for (int l = 0; l < m; l++) for (int a = 0; a < r; a++)
      REAL(sexp_v0)[(size_t) a * m + l] = V[(size_t) l * r + a];
  }


  //--------------------------------------------------
  // OptSpace's gradient descent. Only R_alloc'd memory
  // is in use from here on, so interrupts are safe.
  //--------------------------------------------------
  double *cost = (double *) R_alloc((size_t) max_iter + 1, sizeof(double));
  int     iter = 0;

  #define GRADIENT() do {                                          \
    mult_s(U, S, n, 0, US);                                        \
    mult_s(V, S, m, 1, VS);                                        \
    p_US = US; p_V = V; p_VS = VS; p_y = RVS; p_part = sv_part;    \
    run_rows(w_grad_rows);                                         \
    double sq = 0;                                                 \
    for (int blk = 0; blk < n_rb; blk++) sq += sv_part[blk];       \
    cost[iter] = sqrt(sq) / sqrt((double) nnz);                    \
    if (!R_FINITE(cost[iter]))                                     \
      error("Robust PCA diverged (non-finite values in OptSpace)."); \
  } while (0)

  solve_s(U, V, S);
  GRADIENT();

  for (iter = 1; iter <= max_iter; ) {

    R_CheckUserInterrupt();

    p_US = US; p_y = RUS;
    run_cols(w_grad_cols);

    direction(U, RVS, n, rho, dU);
    direction(V, RUS, m, rho, dV);

    double norm_update = 0;
    for (size_t i = 0; i < nr_; i++) norm_update += dU[i] * dU[i];
    {
      double t = 0;
      for (size_t i = 0; i < mr_; i++) t += dV[i] * dV[i];
      norm_update += t;
    }

    // Line search, as gemelli's line_search().
    mult_s(dU, S, n, 0, dUS);
    p_US = US; p_dUS = dUS; p_V = V; p_dV = dV; p_part = sv_part;
    run_rows(w_line_rows);

    double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    for (int blk = 0; blk < n_rb; blk++) {
      c1 += sv_part[4 * blk + 0];
      c2 += sv_part[4 * blk + 1];
      c3 += sv_part[4 * blk + 2];
      c4 += sv_part[4 * blk + 3];
    }

    double g0   = penalty(U, dU, 0, n) + penalty(V, dV, 0, m);
    double line = -1e-1;
    for (int i = 0; i < RPCA_RES_LIMIT; i++) {
      double t  = line;
      double dc = t * c1 + 0.5 * t * t * c2 + t * t * t * c3 + 0.5 * t * t * t * t * c4 +
        rho * (penalty(U, dU, t, n) + penalty(V, dV, t, m) - g0);
      if (dc <= 0.5 * line * norm_update) break;
      line /= 2;
    }

    for (size_t i = 0; i < nr_; i++) U[i] += line * dU[i];
    for (size_t i = 0; i < mr_; i++) V[i] += line * dV[i];

    solve_s(U, V, S);
    GRADIENT();

    if (cost[iter] < tol) { iter++; break; }
    iter++;
  }
  int n_iters = iter - 1; // Iterations run; cost[0 .. n_iters].

  #undef GRADIENT


  //--------------------------------------------------
  // Undo the rescaling, and sort by S's diagonal, as
  // gemelli's svd_sort().
  //--------------------------------------------------
  for (int i = 0; i < r2; i++) S[i] /= rescal;

  int ord[RPCA_MAX_K];
  for (int a = 0; a < r; a++) ord[a] = a;
  for (int a = 0; a < r - 1; a++) {
    int best = a;
    for (int b = a + 1; b < r; b++)
      if (S[ord[b] * r + ord[b]] > S[ord[best] * r + ord[best]]) best = b;
    int t = ord[a]; ord[a] = ord[best]; ord[best] = t;
  }
  {
    double S2[RPCA_MAX_K * RPCA_MAX_K], row[RPCA_MAX_K];
    for (int a = 0; a < r; a++)
      for (int b = 0; b < r; b++)
        S2[a * r + b] = S[ord[a] * r + ord[b]];
    memcpy(S, S2, sizeof(double) * r2);
    for (int k = 0; k < n; k++) {
      for (int a = 0; a < r; a++) row[a] = U[(size_t) k * r + ord[a]];
      memcpy(U + (size_t) k * r, row, sizeof(double) * r);
    }
    for (int l = 0; l < m; l++) {
      for (int a = 0; a < r; a++) row[a] = V[(size_t) l * r + ord[a]];
      memcpy(V + (size_t) l * r, row, sizeof(double) * r);
    }
  }


  //--------------------------------------------------
  // PCA of J_n U S V' J_m = Qu (Ru S Rv') Qv'.
  //--------------------------------------------------
  double *Qu = (double *) R_alloc(nr_, sizeof(double)); // column-major
  double *Qv = (double *) R_alloc(mr_, sizeof(double));
  double  Ru[RPCA_MAX_K * RPCA_MAX_K], Rv[RPCA_MAX_K * RPCA_MAX_K];
  double  K[RPCA_MAX_K * RPCA_MAX_K],  T1[RPCA_MAX_K * RPCA_MAX_K];
  double  Z[RPCA_MAX_K * RPCA_MAX_K],  sv[RPCA_MAX_K];

  for (int a = 0; a < r; a++) {
    double mu = 0;
    for (int k = 0; k < n; k++) mu += U[(size_t) k * r + a];
    mu /= n;
    for (int k = 0; k < n; k++) Qu[(size_t) a * n + k] = U[(size_t) k * r + a] - mu;
    double mv = 0;
    for (int l = 0; l < m; l++) mv += V[(size_t) l * r + a];
    mv /= m;
    for (int l = 0; l < m; l++) Qv[(size_t) a * m + l] = V[(size_t) l * r + a] - mv;
  }
  thin_qr(Qu, n, r, Ru);
  thin_qr(Qv, m, r, Rv);

  // K = Ru S Rv' (all column-major; S is row-major).
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r; j++) {
      double s = 0;
      for (int a = 0; a < r; a++) s += Ru[a * r + i] * S[a * r + j];
      T1[j * r + i] = s; // (Ru S)[i, j]
    }
  for (int i = 0; i < r; i++)
    for (int j = 0; j < r; j++) {
      double s = 0;
      for (int a = 0; a < r; a++) s += T1[a * r + i] * Rv[a * r + j];
      K[j * r + i] = s;
    }

  small_svd(K, r, Z, sv); // K Z = G (in K)

  SEXP sexp_scores   = PROTECT(allocMatrix(REALSXP, n, r)); nprot++;
  SEXP sexp_loadings = PROTECT(allocMatrix(REALSXP, m, r)); nprot++;
  SEXP sexp_eig      = PROTECT(allocVector(REALSXP, r));    nprot++;
  double *scores   = REAL(sexp_scores);
  double *loadings = REAL(sexp_loadings);

  // scores = Qu G = u * s; loadings = Qv Z = v.
  for (int j = 0; j < r; j++) {
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int a = 0; a < r; a++) s += Qu[(size_t) a * n + i] * K[j * r + a];
      scores[(size_t) j * n + i] = s;
    }
    for (int i = 0; i < m; i++) {
      double s = 0;
      for (int a = 0; a < r; a++) s += Qv[(size_t) a * m + i] * Z[j * r + a];
      loadings[(size_t) j * m + i] = s;
    }
  }

  // Each axis's largest score (the first, if tied) is positive.
  for (int j = 0; j < r; j++) {

    double *x = scores   + (size_t) n * j;
    double *y = loadings + (size_t) m * j;

    if (x[lz_sign_index(x, n)] < 0) {
      for (int i = 0; i < n; i++) x[i] = -x[i];
      for (int i = 0; i < m; i++) y[i] = -y[i];
    }

    REAL(sexp_eig)[j] = sv[j] * sv[j];
  }


  //--------------------------------------------------
  // gemelli's distances: between rows of OptSpace's U,
  // in chunks of about 2^24, checking for interrupts
  // between them.
  //--------------------------------------------------
  SEXP sexp_dist = R_NilValue;
  if (want_dm) {
    sexp_dist = PROTECT(allocVector(REALSXP, (R_xlen_t) n * (n - 1) / 2)); nprot++;
    p_dist = REAL(sexp_dist);
    p_U_dm = U;
    p_j0 = 0;
    while (p_j0 < n - 1) {
      double left = 0;
      p_j1 = p_j0;
      while (p_j1 < n - 1 && left < 16777216.0) left += n - 1 - p_j1++;
      int ntd = threads_for(n_threads, left * r, RPCA_MIN_WORK);
      run_parallel(w_dist, ntd, p_j1 - p_j0);
      p_j0 = p_j1;
      R_CheckUserInterrupt();
    }
  }


  //--------------------------------------------------
  // Results.
  //--------------------------------------------------
  SEXP sexp_cost = PROTECT(allocVector(REALSXP, n_iters + 1)); nprot++;
  memcpy(REAL(sexp_cost), cost, sizeof(double) * (n_iters + 1));

  SEXP sexp_u = R_NilValue, sexp_s = R_NilValue, sexp_v = R_NilValue;
  if (debug) {
    sexp_u = PROTECT(allocMatrix(REALSXP, n, r)); nprot++;
    sexp_s = PROTECT(allocMatrix(REALSXP, r, r)); nprot++;
    sexp_v = PROTECT(allocMatrix(REALSXP, m, r)); nprot++;
    for (int k = 0; k < n; k++) for (int a = 0; a < r; a++)
      REAL(sexp_u)[(size_t) a * n + k] = U[(size_t) k * r + a];
    for (int l = 0; l < m; l++) for (int a = 0; a < r; a++)
      REAL(sexp_v)[(size_t) a * m + l] = V[(size_t) l * r + a];
    for (int a = 0; a < r; a++) for (int b = 0; b < r; b++)
      REAL(sexp_s)[b * r + a] = S[a * r + b];
  }

  const char *names[] = {
    "scores", "loadings", "eig", "k", "cost", "n_iters", "dist",
    "u0", "v0", "u", "s", "v", "ss", "lanczos_steps", "lanczos_converged",
    "threads", "" };
  SEXP res = PROTECT(mkNamed(VECSXP, names)); nprot++;
  SET_VECTOR_ELT(res,  0, sexp_scores);
  SET_VECTOR_ELT(res,  1, sexp_loadings);
  SET_VECTOR_ELT(res,  2, sexp_eig);
  SET_VECTOR_ELT(res,  3, ScalarInteger(r));
  SET_VECTOR_ELT(res,  4, sexp_cost);
  SET_VECTOR_ELT(res,  5, ScalarInteger(n_iters));
  SET_VECTOR_ELT(res,  6, sexp_dist);
  SET_VECTOR_ELT(res,  7, sexp_u0);
  SET_VECTOR_ELT(res,  8, sexp_v0);
  SET_VECTOR_ELT(res,  9, sexp_u);
  SET_VECTOR_ELT(res, 10, sexp_s);
  SET_VECTOR_ELT(res, 11, sexp_v);
  SET_VECTOR_ELT(res, 12, ScalarReal(ss));
  SET_VECTOR_ELT(res, 13, ScalarInteger(lz_steps));
  SET_VECTOR_ELT(res, 14, ScalarLogical(lz_conv));
  SET_VECTOR_ELT(res, 15, ScalarInteger(nt_max));

  UNPROTECT(nprot);
  return res;
}



// Frees new_ecomatrix()'s allocations, then raises an R error.
#define STOP(...) do { free_all(); error(__VA_ARGS__); } while (0)


//======================================================
// R interface for rpca().
//
// `sexp_counts`/`sexp_margin`: anything new_ecomatrix()
// takes. `sexp_k`: rank. `sexp_iters`, `sexp_tol`:
// OptSpace's maximum iterations and tol. `sexp_msc`,
// `sexp_mfc`, `sexp_mff`: gemelli's table filters.
// `sexp_dm`: also return distances. `sexp_debug`: also
// return the intermediate results that the tests check
// against gemelli's. `sexp_step`: OptSpace's step_size.
//
// Returns rpca_core()'s list, followed by samples,
// features, single, rclr: the kept samples and features
// (1-based), the samples dropped for having a single
// non-zero feature, and (with `sexp_debug`) the robust
// clr values in the order of the non-zero counts by
// sample.
//======================================================
SEXP C_rpca(
    SEXP sexp_counts, SEXP sexp_margin,  SEXP sexp_k,
    SEXP sexp_iters,  SEXP sexp_tol,     SEXP sexp_msc,
    SEXP sexp_mfc,    SEXP sexp_mff,     SEXP sexp_dm,
    SEXP sexp_debug,  SEXP sexp_step,    SEXP sexp_n_threads ) {

  int    k_req     = asInteger(sexp_k);
  int    max_iter  = asInteger(sexp_iters);
  double tol       = asReal(sexp_tol);
  double msc       = asReal(sexp_msc);
  double mfc       = asReal(sexp_mfc);
  double mff       = asReal(sexp_mff);
  int    want_dm   = asLogical(sexp_dm);
  int    debug     = asLogical(sexp_debug);
  int    n_threads = asInteger(sexp_n_threads);
  double step      = asReal(sexp_step);

  if (k_req < 1 || k_req > RPCA_MAX_K || max_iter < 1 || n_threads < 1 || !(step > 0))
    error("Invalid arguments to C_rpca()."); // # nocov


  //--------------------------------------------------
  // Read the table: each sample's non-zero entries, in
  // order of feature, repeats summed.
  //--------------------------------------------------
  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int           N0   = em->n_samples, P0 = em->n_otus, e1 = em->nnz;
  const int    *tpos = em->pos_vec, *tidx = em->otu_vec;
  const double *tval = em->val_vec;

  for (int e = 0; e < e1; e++) {
    if (!R_FINITE(tval[e])) STOP("`counts` must not contain NA, NaN, or infinite values.");
    if (tval[e] < 0)        STOP("`counts` must not contain negative values.");
  }


  //--------------------------------------------------
  // Filters.
  //--------------------------------------------------
  int *f_new  = (int *) R_alloc(P0 > 0 ? P0 : 1, sizeof(int));
  int *s_new  = (int *) R_alloc(N0 > 0 ? N0 : 1, sizeof(int));
  int *single = (int *) R_alloc(N0 > 0 ? N0 : 1, sizeof(int));
  int  n, m, n_single;

  rpca_filter(N0, P0, tpos, tidx, tval, msc, mfc, mff, 1,
              s_new, f_new, single, &n_single, &n, &m);

  if (n < 2 || m < 2) {
    if (n_single > 0)
      STOP(
        "After filtering, %d sample%s and %d feature%s remain (%d sample%s "
        "with a single non-zero feature %s removed); at least 2 of each are needed.",
        n, (n == 1) ? "" : "s", m, (m == 1) ? "" : "s",
        n_single, (n_single == 1) ? "" : "s", (n_single == 1) ? "was" : "were" );
    STOP(
      "After filtering, %d sample%s and %d feature%s remain; at least 2 of each are needed.",
      n, (n == 1) ? "" : "s", m, (m == 1) ? "" : "s" );
  }

  int r = (k_req < n - 1) ? k_req : n - 1;
  if (r > m - 1) r = m - 1;


  //--------------------------------------------------
  // CSR copy of the kept table (rows = samples).
  //--------------------------------------------------
  int *s_idx = (int *) R_alloc(n, sizeof(int));
  int *f_idx = (int *) R_alloc(m, sizeof(int));
  for (int i = 0; i < N0; i++) if (s_new[i] >= 0) s_idx[s_new[i]] = i;
  for (int j = 0; j < P0; j++) if (f_new[j] >= 0) f_idx[f_new[j]] = j;

  int nnz = 0;
  for (int k = 0; k < n; k++)
    for (int e = tpos[s_idx[k]]; e < tpos[s_idx[k] + 1]; e++)
      if (f_new[tidx[e]] >= 0) nnz++;
  size_t nnz_z = (size_t) nnz;

  int    *rpos = (int *)    R_alloc((size_t) n + 1, sizeof(int));
  int    *ridx = (int *)    R_alloc(nnz_z, sizeof(int));
  double *rval = (double *) R_alloc(nnz_z, sizeof(double));

  nnz = 0;
  for (int k = 0; k < n; k++) {
    int i = s_idx[k];
    rpos[k] = nnz;
    for (int e = tpos[i]; e < tpos[i + 1]; e++)
      if (f_new[tidx[e]] >= 0) {
        ridx[nnz] = f_new[tidx[e]];
        rval[nnz] = tval[e];
        nnz++;
      }
  }
  rpos[n] = nnz;


  //--------------------------------------------------
  // Robust clr (normalize.c), on a copy that
  // normalize() makes; then new_ecomatrix()'s memory
  // can go. Each sample is done by one thread, so the
  // number of threads doesn't change the result.
  //--------------------------------------------------
  {
    ecomatrix_t fem = {
      .n_samples = n, .n_otus = m, .nnz = nnz,
      .pos_vec = rpos, .otu_vec = ridx, .val_vec = rval };
    int nt = threads_for(n_threads, (double) nnz * (r + 1), RPCA_MIN_WORK);
    normalize(&fem, RPCA_NORM_RCLR, nt, 0);
    memcpy(rval, fem.val_vec, nnz_z * sizeof(double));
  }
  free_all();

  int nprot = 0;

  SEXP sexp_rclr = R_NilValue;
  if (debug) {
    sexp_rclr = PROTECT(allocVector(REALSXP, nnz)); nprot++;
    memcpy(REAL(sexp_rclr), rval, nnz_z * sizeof(double));
  }

  SEXP core = PROTECT(rpca_core(
    n, m, rpos, ridx, rval, r, max_iter, tol, want_dm, debug, step, n_threads,
    "There is no variation to ordinate: every sample's non-zero counts are equal." ));
  nprot++;

  SEXP sexp_sam    = PROTECT(allocVector(INTSXP, n));        nprot++;
  SEXP sexp_fea    = PROTECT(allocVector(INTSXP, m));        nprot++;
  SEXP sexp_single = PROTECT(allocVector(INTSXP, n_single)); nprot++;
  for (int k = 0; k < n; k++)        INTEGER(sexp_sam)[k]    = s_idx[k] + 1;
  for (int l = 0; l < m; l++)        INTEGER(sexp_fea)[l]    = f_idx[l] + 1;
  for (int i = 0; i < n_single; i++) INTEGER(sexp_single)[i] = single[i] + 1;

  const char *names[] = { "samples", "features", "single", "rclr" };
  SEXP        vals[]  = { sexp_sam, sexp_fea, sexp_single, sexp_rclr };
  SEXP        res     = rpca_list_append(core, names, vals, 4);

  UNPROTECT(nprot);
  return res;
}
