// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The orthogonal projections to latent structures (O-PLS) fit in this file
// (the weight vector w = X'y / |X'y|, and each orthogonal component's
// scores t, loadings p, orthogonal weights w_o, scores t_o and loadings p_o,
// with the deflation X <- X - t_o p_o') is ported from pyopls 20.3.post1:
// `OPLS.fit()` and `_center_scale_xy()` in pyopls/opls.py. That file is
// covered by this license:
//
//   MIT License
//
//   Copyright (c) 2019 Wright State University
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


/*
 * Orthogonal partial least squares discriminant analysis (oplsda()):
 * O-PLS (Trygg & Wold 2002) of one centred (and scaled) response y, the
 * indicator of the second of two groups, on a count table.
 *
 * With A = J * Z * D^-1 the centred (and scaled) table of pca.c (n x p,
 * never formed) and X_a = A - T_o P_o' what is left after the orthogonal
 * components so far (plsda.c's pls_x_t):
 *
 *     w   = A' y / |A' y|          the predictive weights (fixed)
 *
 * and, for each orthogonal component,
 *
 *     t   = X_a w,  p = X_a' t / (t't)
 *     w_o = p - (w'p) w,  w_o /= |w_o|      orthogonal weights
 *     t_o = X_a w_o                         orthogonal scores
 *     p_o = X_a' t_o / (t_o't_o)            orthogonal loadings
 *     X_{a+1} = X_a - t_o p_o'
 *
 * (pyopls's OPLS.fit()). The predictive component is then the one-component
 * PLS of y on the filtered table X_k (as pyopls's README fits it with
 * scikit-learn's PLSRegression(1)):
 *
 *     t_p = X_k w,  p_p = X_k' t_p / (t_p't_p),  q_p = y' t_p / (t_p't_p).
 *
 * X_k' y = A' y (each t_o is orthogonal to y), so its weights are w again.
 * All the scores are orthogonal, so the variances |t|^2 |p|^2 / (n - 1)
 * that the components explain add up.
 *
 * If |A'y| is at most sqrt(DBL_EPSILON) |A| |y| (|A| the Frobenius norm),
 * the table has no covariance with y beyond rounding error, and there is
 * no predictive component: C_oplsda() returns NULL.
 *
 * Fitting of orthogonal components stops early, with fewer of them, once
 * |w_o| (before normalizing) falls to sqrt(DBL_EPSILON) times |p|: what is
 * left of the table then has no variation orthogonal to y, beyond rounding
 * error.
 *
 * Signs: each component is flipped so that its score with the largest
 * absolute value is positive (lz_sign_index(), as for every ecodive
 * ordination), with its weights and loadings (and q_p).
 *
 * Determinism: the products are pca.c's, summed in fixed blocks; all else
 * runs on the main thread in a fixed order.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON



//======================================================
// R interface: O-PLS of the n-vector y (centred and, if
// `scale`, scaled, by the caller) on the table, as
// pca_table_init() reads it, with up to sexp_k
// orthogonal components.
//
// Returns list(scores (n x (1 + k)), weights, loadings
// (p x (1 + k)), y_loading, eig (1 + k), total), where
// k is the number of orthogonal components fitted and
// column 1 is the predictive component; or NULL if
// A'y is 0, beyond rounding error (no predictive
// component).
//======================================================
SEXP C_oplsda(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_scale, SEXP sexp_y,
    SEXP sexp_k, SEXP sexp_n_threads ) {

  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    scale       = asLogical(sexp_scale) == TRUE;
  int    max_k       = asInteger(sexp_k);
  int    n_threads   = asInteger(sexp_n_threads);

  if (max_k < 0 || n_threads < 1 || !isReal(sexp_y))
    error("Invalid arguments to C_oplsda()."); // # nocov

  init_n_ptrs(10);

  pca_table_t tbl;
  double total = pca_table_init(
    sexp_counts, sexp_margin, norm, pseudocount, scale, n_threads, &tbl );

  int n = tbl.n, p = tbl.p;

  if (length(sexp_y) != n) {
    free_all();                                 // # nocov
    error("Invalid arguments to C_oplsda().");  // # nocov
  }

  int m = max_k + 1; // Room for the predictive component too.

  pls_x_t x = {
    .n = n, .p = p, .k = 0, .max_k = max_k, .n_threads = n_threads,
    .T = (double *) R_alloc((size_t) n * (max_k > 0 ? max_k : 1), sizeof(double)),
    .P = (double *) R_alloc((size_t) p * (max_k > 0 ? max_k : 1), sizeof(double)) };

  const double *y = REAL(sexp_y);

  double *w  = (double *) R_alloc(p, sizeof(double));
  double *M  = (double *) R_alloc(p, sizeof(double));
  double *Wo = (double *) R_alloc((size_t) p * m, sizeof(double));
  double *t  = (double *) R_alloc(n, sizeof(double));
  double *pl = (double *) R_alloc(p, sizeof(double));
  double *tp = (double *) R_alloc(n, sizeof(double));
  double *pp = (double *) R_alloc(p, sizeof(double));

  // The predictive weights: w = A'y / |A'y|. None if |A'y| is 0 beyond
  // rounding error: at most sqrt(DBL_EPSILON) |A| |y|.
  double s  = pls_weights(&x, y, 1, M, w);
  double ny = sqrt(lz_dot(y, y, n));
  if (!(s > sqrt(DBL_EPSILON) * sqrt(total) * ny)) {
    free_all();
    return R_NilValue;
  }

  int k = 0;

  for (int a = 0; a < max_k; a++) {

    // t = X_a w, p = X_a' t / t't.
    pls_x_mult(&x, w, t);
    double tt = lz_dot(t, t, n);
    if (!(tt > 0)) break; // # nocov
    pls_xt_mult(&x, t, 1, pl);
    for (int j = 0; j < p; j++) pl[j] /= tt;

    // w_o = p - (w'p / w'w) w, normalized.
    double *wo  = Wo + (size_t) a * p;
    double  wp  = lz_dot(w, pl, p) / lz_dot(w, w, p);
    for (int j = 0; j < p; j++) wo[j] = pl[j] - wp * w[j];
    double nwo = sqrt(lz_dot(wo, wo, p));
    double npl = sqrt(lz_dot(pl, pl, p));
    if (!(nwo > sqrt(DBL_EPSILON) * npl)) break;
    for (int j = 0; j < p; j++) wo[j] /= nwo;

    // t_o = X_a w_o, p_o = X_a' t_o / t_o't_o.
    pls_x_mult(&x, wo, t);
    tt = lz_dot(t, t, n);
    if (!(tt > 0)) break; // # nocov
    pls_xt_mult(&x, t, 1, pl);
    for (int j = 0; j < p; j++) pl[j] /= tt;

    pls_x_deflate(&x, t, pl);
    k = a + 1;

    if (user_interrupted()) { free_all(); error("Interrupted."); } // # nocov
  }

  // The predictive component, from the filtered table.
  pls_x_mult(&x, w, tp);
  double tt = lz_dot(tp, tp, n);
  pls_xt_mult(&x, tp, 1, pp);
  for (int j = 0; j < p; j++) pp[j] /= tt;
  double q = lz_dot(y, tp, n) / tt;


  //--------------------------------------------------
  // Results, each component flipped so that its largest
  // score is positive.
  //--------------------------------------------------
  int kk = k + 1;

  SEXP sexp_scores = PROTECT(allocMatrix(REALSXP, n, kk));
  SEXP sexp_w      = PROTECT(allocMatrix(REALSXP, p, kk));
  SEXP sexp_p      = PROTECT(allocMatrix(REALSXP, p, kk));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, kk));

  for (int a = 0; a < kk; a++) {

    const double *ta = (a == 0) ? tp : x.T + (size_t) (a - 1) * n;
    const double *wa = (a == 0) ? w  : Wo  + (size_t) (a - 1) * p;
    const double *pa = (a == 0) ? pp : x.P + (size_t) (a - 1) * p;

    double sgn = (ta[lz_sign_index(ta, n)] < 0) ? -1 : 1;
    if (a == 0) q *= sgn;

    for (int i = 0; i < n; i++) REAL(sexp_scores)[(size_t) a * n + i] = sgn * ta[i];
    for (int j = 0; j < p; j++) {
      REAL(sexp_w)[(size_t) a * p + j] = sgn * wa[j];
      REAL(sexp_p)[(size_t) a * p + j] = sgn * pa[j];
    }

    REAL(sexp_eig)[a] = lz_dot(ta, ta, n) * lz_dot(pa, pa, p) / (n - 1);
  }

  free_all();

  const char *names[] = {
    "scores", "weights", "loadings", "y_loading", "eig", "total", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_scores);
  SET_VECTOR_ELT(sexp_result, 1, sexp_w);
  SET_VECTOR_ELT(sexp_result, 2, sexp_p);
  SET_VECTOR_ELT(sexp_result, 3, ScalarReal(q));
  SET_VECTOR_ELT(sexp_result, 4, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 5, ScalarReal(total));

  UNPROTECT(5);
  return sexp_result;
}
