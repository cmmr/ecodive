// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The partial least squares fit in this file (centring and scaling, the
// weights, scores and loadings of each component, the regression-mode
// deflation of X and Y, and the sign rule applied while fitting) is ported
// from scikit-learn 1.9.1: `_PLS.fit()`, `_center_scale_xy()`,
// `_get_first_singular_vectors_svd()` and `_svd_flip_1d()` in
// sklearn/cross_decomposition/_pls.py, as used by `PLSRegression`. That
// file is covered by this license:
//
//   BSD 3-Clause License
//
//   Copyright (c) 2007-2026 The scikit-learn developers.
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
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


/*
 * Partial least squares discriminant analysis (plsda()): PLS regression
 * (scikit-learn's PLSRegression; mode A, regression deflation) of a dummy
 * matrix of the groups on a count table.
 *
 * With A = J * Z * D^-1 the centred (and scaled) table of pca.c (n x p,
 * never formed) and Y the centred (and scaled) n x g dummy matrix,
 * component a (a = 1, 2, ...) works on the deflated matrices
 *
 *     X_a = A - T P'       (T, P: the scores and loadings so far)
 *     Y_a = Y - T Q'       (Q: the Y loadings so far)
 *
 * and finds
 *
 *     w = the top left singular vector of M = X_a' Y_a  (p x g)
 *     t = X_a w                  scores
 *     p = X_a' t / (t't)         X loadings
 *     q = Y_a' t / (t't)         Y loadings (= scikit-learn's y_weights)
 *
 * then deflates: X_{a+1} = X_a - t p', Y_{a+1} = Y_a - t q'.
 *
 * scikit-learn's PLSRegression finds w by NIPALS, a power iteration that
 * converges to that singular vector (it stops once the squared change in
 * w is below `tol`, 1e-6 by default). Here it is found directly, as
 * scikit-learn's `algorithm = "svd"` does: the top eigenvector c of the
 * small g x g matrix M'M, and w = M c / |M c|. So the result is NIPALS'
 * fixed point, without its stopping error. With one column of Y (as for
 * OPLS-DA), w = M / |M|, which is also NIPALS' one step.
 *
 * X_a is never formed: X_a v = A v - T (P' v), X_a' u = A' u - P (T' u),
 * through pca.c's products, which read only the non-zero values of a
 * sparse table. Y_a is n x g and is deflated in place.
 *
 * A component explains the variance |t|^2 |p|^2 / (n - 1) of X (the
 * deflations are orthogonal projections, so these add up to what the
 * components explain together), and |t|^2 |q|^2 of Y's sum of squares.
 *
 * Fitting stops early, with fewer components, once the top singular value
 * of M falls to sqrt(DBL_EPSILON) times the first component's: what is
 * left of X then has no covariance with Y, beyond rounding error (e.g.,
 * once the components span every direction of X).
 *
 * Signs: each component is flipped so that its score with the largest
 * absolute value is positive (lz_sign_index(), as for every ecodive
 * ordination); w, t, p and q flip together, and the deflations don't
 * depend on the sign. (scikit-learn makes w's largest element positive.)
 *
 * Determinism: the products are pca.c's, summed in fixed blocks; all else
 * runs on the main thread in a fixed order.
 */

#include "ecodive.h"

#include <float.h> // DBL_EPSILON


//======================================================
// The deflated table X_a = A - T P', for plsda() (and,
// later, oplsda()). pca_table_init() must have set up A.
//======================================================

// out = X_a v (n), for v of length p.
void pls_x_mult(const pls_x_t *x, const double *v, double *out) {

  int n = x->n, p = x->p;

  pca_mult(v, 1, out, x->n_threads);

  for (int b = 0; b < x->k; b++) {
    double pv = lz_dot(x->P + (size_t) b * p, v, p);
    const double *t = x->T + (size_t) b * n;
    for (int i = 0; i < n; i++) out[i] -= t[i] * pv;
  }
}


// L = X_a' U (p x m), for U n x m.
void pls_xt_mult(const pls_x_t *x, const double *U, int m, double *L) {

  int n = x->n, p = x->p;

  pca_cross(U, m, L, x->n_threads);

  for (int v = 0; v < m; v++) {
    double *l = L + (size_t) v * p;
    for (int b = 0; b < x->k; b++) {
      double tu = lz_dot(x->T + (size_t) b * n, U + (size_t) v * n, n);
      const double *pb = x->P + (size_t) b * p;
      for (int j = 0; j < p; j++) l[j] -= pb[j] * tu;
    }
  }
}


// Deflate by score t (n) and loading p (p): append them
// as the next column of T and P. Room for x->max_k.
void pls_x_deflate(pls_x_t *x, const double *t, const double *p) {
  if (x->k >= x->max_k) { free_all(); error("Too many PLS components."); } // # nocov
  memcpy(x->T + (size_t) x->k * x->n, t, (size_t) x->n * sizeof(double));
  memcpy(x->P + (size_t) x->k * x->p, p, (size_t) x->p * sizeof(double));
  x->k++;
}


// The X weights of the next component: w (p), the unit
// top left singular vector of M = X_a' Y (p x g), for Y
// n x g. Returns M's top singular value (0 if M is 0;
// w is then 0). `M` needs p * g doubles of scratch.
double pls_weights(const pls_x_t *x, const double *Y, int g, double *M, double *w) {

  int p = x->p;

  pls_xt_mult(x, Y, g, M);

  double *G    = (double *) R_alloc((size_t) g * g, sizeof(double));
  double *vals = (double *) R_alloc(g, sizeof(double));
  double *vecs = (double *) R_alloc((size_t) g * g, sizeof(double));
  double *tmp  = (double *) R_alloc((size_t) 2 * g, sizeof(double));

  for (int a = 0; a < g; a++)
    for (int b = a; b < g; b++)
      G[(size_t) a * g + b] = G[(size_t) b * g + a] =
        lz_dot(M + (size_t) a * p, M + (size_t) b * p, p);

  if (g == 1) { vals[0] = G[0]; vecs[0] = 1; }
  else if (lz_sym_eigen(G, g, vals, vecs, tmp)) {
    free_all();                                                    // # nocov
    error("PLS weights: eigen-decomposition failed to converge."); // # nocov
  }

  // w = M c / |M c|.
  for (int j = 0; j < p; j++) {
    double s = 0;
    for (int a = 0; a < g; a++) s += M[(size_t) a * p + j] * vecs[a];
    w[j] = s;
  }
  double nw = sqrt(lz_dot(w, w, p));
  if (!(nw > 0)) { memset(w, 0, (size_t) p * sizeof(double)); return 0; }
  for (int j = 0; j < p; j++) w[j] /= nw;

  // |M c| is the top singular value; more accurate than
  // sqrt(vals[0]).
  return nw;
}



//======================================================
// R interface: PLS regression of the n x g matrix Y
// (centred and, if `scale`, scaled, by the caller) on
// the table, as pca_table_init() reads it.
//
// Returns list(scores (n x k), weights (p x k),
// loadings (p x k), y_loadings (g x k), eig (k), total,
// sv (k)), where k is the number of components fitted
// (at most sexp_k): eig is each component's variance of
// X, |t|^2 |p|^2 / (n - 1); total the squared Frobenius
// norm of A (so trace = total / (n - 1)); sv the top
// singular value of X_a' Y_a.
//======================================================
SEXP C_plsda(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_norm,
    SEXP sexp_pseudocount, SEXP sexp_scale, SEXP sexp_y,
    SEXP sexp_k, SEXP sexp_n_threads ) {

  int    norm        = asInteger(sexp_norm);
  double pseudocount = asReal(sexp_pseudocount);
  int    scale       = asLogical(sexp_scale) == TRUE;
  int    max_k       = asInteger(sexp_k);
  int    n_threads   = asInteger(sexp_n_threads);
  int    g           = ncols(sexp_y);

  if (max_k < 1 || n_threads < 1 || g < 1 || !isReal(sexp_y))
    error("Invalid arguments to C_plsda()."); // # nocov

  init_n_ptrs(10);

  pca_table_t tbl;
  double total = pca_table_init(
    sexp_counts, sexp_margin, norm, pseudocount, scale, n_threads, &tbl );

  int n = tbl.n, p = tbl.p;

  if (nrows(sexp_y) != n) {
    free_all();                                 // # nocov
    error("Invalid arguments to C_plsda().");   // # nocov
  }

  pls_x_t x = {
    .n = n, .p = p, .k = 0, .max_k = max_k, .n_threads = n_threads,
    .T = (double *) R_alloc((size_t) n * max_k, sizeof(double)),
    .P = (double *) R_alloc((size_t) p * max_k, sizeof(double)) };

  double *Y  = (double *) R_alloc((size_t) n * g, sizeof(double));
  double *M  = (double *) R_alloc((size_t) p * g, sizeof(double));
  double *W  = (double *) R_alloc((size_t) p * max_k, sizeof(double));
  double *Q  = (double *) R_alloc((size_t) g * max_k, sizeof(double));
  double *ev = (double *) R_alloc(max_k, sizeof(double));
  double *sv = (double *) R_alloc(max_k, sizeof(double));
  double *t  = (double *) R_alloc(n, sizeof(double));
  double *pl = (double *) R_alloc(p, sizeof(double));
  memcpy(Y, REAL(sexp_y), (size_t) n * g * sizeof(double));

  double y_eps = DBL_EPSILON;
  double sv1   = 0;
  int    k     = 0;

  for (int a = 0; a < max_k; a++) {

    // Columns of Y that are all but 0 become 0 (scikit-learn).
    for (int v = 0; v < g; v++) {
      double *y = Y + (size_t) v * n;
      int     z = 1;
      for (int i = 0; i < n && z; i++) if (!(fabs(y[i]) < 10 * y_eps)) z = 0;
      if (z) memset(y, 0, (size_t) n * sizeof(double));
    }

    double *w = W + (size_t) a * p;
    double  s = pls_weights(&x, Y, g, M, w);

    if (a == 0) sv1 = s;
    if (!(s > 0) || s <= sqrt(DBL_EPSILON) * sv1) break;

    // Scores, X loadings, Y loadings.
    pls_x_mult(&x, w, t);
    double tt = lz_dot(t, t, n);
    if (!(tt > 0)) break; // # nocov

    pls_xt_mult(&x, t, 1, pl);
    for (int j = 0; j < p; j++) pl[j] /= tt;

    double *q = Q + (size_t) a * g;
    for (int v = 0; v < g; v++) q[v] = lz_dot(Y + (size_t) v * n, t, n) / tt;

    // Deflate.
    for (int v = 0; v < g; v++) {
      double *y = Y + (size_t) v * n;
      for (int i = 0; i < n; i++) y[i] -= t[i] * q[v];
    }
    pls_x_deflate(&x, t, pl);

    ev[a] = tt * lz_dot(pl, pl, p) / (n - 1);
    sv[a] = s;
    k     = a + 1;

    if (user_interrupted()) { free_all(); error("Interrupted."); } // # nocov
  }


  //--------------------------------------------------
  // Results, each component flipped so that its largest
  // score is positive.
  //--------------------------------------------------
  SEXP sexp_scores = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_w      = PROTECT(allocMatrix(REALSXP, p, k));
  SEXP sexp_p      = PROTECT(allocMatrix(REALSXP, p, k));
  SEXP sexp_q      = PROTECT(allocMatrix(REALSXP, g, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_sv     = PROTECT(allocVector(REALSXP, k));

  for (int a = 0; a < k; a++) {

    const double *ta = x.T + (size_t) a * n;
    double sgn = (ta[lz_sign_index(ta, n)] < 0) ? -1 : 1;

    for (int i = 0; i < n; i++) REAL(sexp_scores)[(size_t) a * n + i] = sgn * ta[i];
    for (int j = 0; j < p; j++) {
      REAL(sexp_w)[(size_t) a * p + j] = sgn * W[(size_t) a * p + j];
      REAL(sexp_p)[(size_t) a * p + j] = sgn * x.P[(size_t) a * p + j];
    }
    for (int v = 0; v < g; v++) REAL(sexp_q)[(size_t) a * g + v] = sgn * Q[(size_t) a * g + v];

    REAL(sexp_eig)[a] = ev[a];
    REAL(sexp_sv)[a]  = sv[a];
  }

  free_all();

  const char *names[] = {
    "scores", "weights", "loadings", "y_loadings", "eig", "total", "sv", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_scores);
  SET_VECTOR_ELT(sexp_result, 1, sexp_w);
  SET_VECTOR_ELT(sexp_result, 2, sexp_p);
  SET_VECTOR_ELT(sexp_result, 3, sexp_q);
  SET_VECTOR_ELT(sexp_result, 4, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 5, ScalarReal(total));
  SET_VECTOR_ELT(sexp_result, 6, sexp_sv);

  UNPROTECT(7);
  return sexp_result;
}
