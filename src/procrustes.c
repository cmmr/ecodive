// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The statistic follows scipy.spatial.procrustes() and
// scipy.linalg.orthogonal_procrustes() from SciPy 1.18.1
// (scipy/spatial/_procrustes.py, scipy/linalg/_procrustes.py):
// both configurations centred and scaled to unit sum of squares,
// and the fit given by the sum of the singular values of their
// cross-product. SciPy's licence:
//
//   Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//
//   3. Neither the name of the copyright holder nor the names of its
//      contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ecodive.h"

#include <float.h>  // DBL_EPSILON, DBL_MIN


//======================================================
// PROTEST (procrustes_test()): the symmetric Procrustes
// correlation r of two configurations, X (n x k, the
// target) and Y (n x k, the one permuted), for the
// observed samples and for each permutation of Y's rows.
//
// X and Y arrive centred, scaled to unit sum of squares,
// and zero-padded to the same k columns. Then
//   r = sum of the singular values of X' Y,
// which is SciPy's `scale` from orthogonal_procrustes(),
// and sqrt(1 - m^2) for the symmetric m^2 (its
// procrustes() disparity). Permutation p takes row i of
// Y to be Y[p[i], ], as vegan's Y[take, ]: the centring
// and scaling are unchanged by it.
//
// Each permutation forms the k x k cross-product, O(n k^2),
// and its singular values by one-sided Jacobi (Hestenes),
// which is accurate to rounding for small k. Each is done
// by one thread, start to finish, so results are
// identical for any thread count.
//======================================================

static int           pr_n;
static int           pr_k;
static int           pr_total;   // Permutations + 1.
static const double *pr_x;       // n x k, row-major.
static const double *pr_y;       // n x k, row-major.
static double       *pr_out;     // pr_total values of r.
static perm_source_t pr_ps;

static int          *pr_iwork;   // n ints per thread.
static double       *pr_dwork;   // 2 k^2 doubles per thread.
static int           pr_k0;      // This round's permutations:
static int           pr_k1;      //   k0 .. k1 - 1.


//------------------------------------------------------
// The sum of the singular values of the k x k matrix G
// (column-major; overwritten), by one-sided Jacobi:
// rotate pairs of columns until all are orthogonal; the
// singular values are then the columns' norms.
//------------------------------------------------------
static double nuclear_norm(double *G, int k) {

  #define G_(i, j) G[(size_t) (j) * k + (i)]

  for (int sweep = 0; sweep < 100; sweep++) {

    int rotated = 0;

    for (int p = 0; p < k - 1; p++)
      for (int q = p + 1; q < k; q++) {

        double alpha = 0, beta = 0, gamma = 0;
        for (int i = 0; i < k; i++) {
          alpha += G_(i, p) * G_(i, p);
          beta  += G_(i, q) * G_(i, q);
          gamma += G_(i, p) * G_(i, q);
        }

        if (fabs(gamma) < DBL_MIN ||
            fabs(gamma) <= DBL_EPSILON * sqrt(alpha) * sqrt(beta))
          continue;

        // The rotation that zeros gamma (Golub & Van Loan's
        // sym.schur2).
        double tau = (beta - alpha) / (2 * gamma);
        double t   = (fabs(tau) > 1e150)
          ? 0.5 / tau
          : ((tau >= 0) ? 1 : -1) / (fabs(tau) + sqrt(1 + tau * tau));
        double c   = 1 / sqrt(1 + t * t);
        double s   = t * c;

        for (int i = 0; i < k; i++) {
          double gp = G_(i, p), gq = G_(i, q);
          G_(i, p) = c * gp - s * gq;
          G_(i, q) = s * gp + c * gq;
        }

        rotated = 1;
      }

    if (!rotated) break;
  }

  double sum = 0;
  for (int j = 0; j < k; j++) {
    double ss = 0;
    for (int i = 0; i < k; i++) ss += G_(i, j) * G_(i, j);
    sum += sqrt(ss);
  }

  #undef G_

  return sum;
}


//------------------------------------------------------
// r for permutation p: the cross-product X' Y[p, ], then
// its nuclear norm. G is k^2 doubles of workspace.
//------------------------------------------------------
static double protest_perm(const int *p, double *G) {

  int n = pr_n, k = pr_k;

  for (int j = 0; j < k * k; j++) G[j] = 0;

  // G[a, b] = sum_i X[i, a] Y[p[i], b], summed over i in
  // order.
  for (int i = 0; i < n; i++) {
    const double *xi = pr_x + (size_t) i    * k;
    const double *yi = pr_y + (size_t) p[i] * k;
    for (int b = 0; b < k; b++) {
      double  v  = yi[b];
      double *gb = G + (size_t) b * k;
      for (int a = 0; a < k; a++) gb[a] += xi[a] * v;
    }
  }

  return nuclear_norm(G, k);
}


static void *protest_worker(void *arg) {

  worker_t *w = (worker_t *) arg;
  int      *p = pr_iwork + (size_t) w->i * pr_n;
  double   *G = pr_dwork + (size_t) w->i * pr_k * pr_k;

  for (int t = pr_k0 + w->i; t < pr_k1; t += w->n) {
    perm_source_get(&pr_ps, t, p);
    pr_out[t] = protest_perm(p, G);
  }

  return NULL;
}


//======================================================
// R interface.
//
// x, y:    k x n matrices (so each sample's coordinates
//          are contiguous): the target and the
//          configuration permuted, centred, scaled to unit
//          sum of squares, and zero-padded to k rows.
// perms:   n x n_perm integer matrix of 0-based
//          permutations, or NULL for random ones.
// strata:  Each sample's stratum, 0 .. n_strata - 1, or NULL.
// opts:    c(n_perm, seed, n_threads).
//
// Returns the n_perm + 1 values of r, the first for the
// observed samples.
//======================================================
SEXP C_protest(
    SEXP sexp_x, SEXP sexp_y, SEXP sexp_perms, SEXP sexp_strata,
    SEXP sexp_opts ) {

  int k         = nrows(sexp_x);
  int n         = ncols(sexp_x);
  int n_perm    = INTEGER(sexp_opts)[0];
  int seed      = INTEGER(sexp_opts)[1];
  int n_threads = INTEGER(sexp_opts)[2];

  if (n < 1 || k < 1 || n_perm < 0 || n_threads < 1 ||
      nrows(sexp_y) != k || ncols(sexp_y) != n)
    error("Invalid arguments to C_protest()."); // # nocov

  pr_n     = n;
  pr_k     = k;
  pr_total = n_perm + 1;
  pr_x     = REAL(sexp_x);
  pr_y     = REAL(sexp_y);

  perm_source_init(&pr_ps, n, sexp_perms, sexp_strata, seed);

  if (n_threads > pr_total) n_threads = pr_total;

  pr_iwork = (int *)    R_alloc((size_t) n * n_threads, sizeof(int));
  pr_dwork = (double *) R_alloc((size_t) k * k * n_threads, sizeof(double));

  SEXP sexp_r = PROTECT(allocVector(REALSXP, pr_total));
  pr_out = REAL(sexp_r);


  //--------------------------------------------------
  // Rounds of permutations, checking for interrupts
  // between them: at least ~1e8 operations per thread,
  // or every remaining one.
  //--------------------------------------------------
  double per_perm = (double) n * (k * k + 2) + 10.0 * k * k * k;

  int per_round = (int) fmin(1e8 * n_threads / per_perm, (double) pr_total);
  if (per_round < n_threads) per_round = n_threads;

  for (pr_k0 = 0; pr_k0 < pr_total; pr_k0 = pr_k1) {

    pr_k1 = (pr_total - pr_k0 > per_round) ? pr_k0 + per_round : pr_total;

    int nt = threads_for(n_threads, per_perm * (pr_k1 - pr_k0), 1e6);
    if (nt > pr_k1 - pr_k0) nt = pr_k1 - pr_k0;
    run_parallel(protest_worker, nt, (nt > 1) ? 100 : 1);

    if (pr_k1 < pr_total) R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return sexp_r;
}
