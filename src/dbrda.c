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
 * The Lanczos iteration and the small dense helpers are
 * shared with pcoa() and pca() (lanczos.c), and the products
 * with G with pcoa() (gower.c).
 *
 * Determinism: results are bit-identical for any number of
 * threads. G * X is computed (gower.c) in fixed blocks of
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

// A design column is aliased (redundant) when Gram-Schmidt
// leaves less than this fraction of its centred norm: the
// same rule and tolerance as R's qr(), which vegan uses.
#define DBRDA_QR_TOL 1e-7


//------------------------------------------------------
// Small dense helpers, shared with lanczos.c. Vectors
// have length n; matrices are column-major.
//------------------------------------------------------
#define dot       lz_dot
#define center    lz_center
#define sym_eigen lz_sym_eigen
#define combine   lz_combine

#define project_out lz_project_out



//------------------------------------------------------
// Flip x (and y, if given) so that x's element with the
// largest absolute value is positive.
//------------------------------------------------------
static void fix_sign(double *x, int n, double *y, int ny) {
  int imax = lz_sign_index(x, n);
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

  //--------------------------------------------------
  // trace(G) = (1/n) * sum_{i<j} d_ij^2, and products
  // with G (gower.c, shared with pcoa()).
  //--------------------------------------------------
  gower_t g;
  double trace_g = gower_init(&g, REAL(sexp_dist), n, 16, n_threads); // unscaled
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
  double *GQ = (double *) R_alloc((size_t) n * (rq + 1), sizeof(double));

  gower_mult(Q, GQ, rq, &g);
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

  if (rx > 0 && sym_eigen(C, rx, lam, Y, tmp))
    error("The constrained eigenproblem did not converge."); // # nocov


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

    g.F = Q;
    g.f = rq;

    lanczos_t lz = {
      .n = n, .k = kr, .dims = dims_r, .deflate = 1,
      .basis = Q, .n_basis = rq, .scale = trace_g,
      .op = gower_mult, .ctx = &g,
      .vecs = r_vecs, .vals = r_vals,
      .resid = (double *) R_alloc(kr, sizeof(double)) };

    lanczos(&lz);
    steps     = lz.steps;
    converged = lz.converged;

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
