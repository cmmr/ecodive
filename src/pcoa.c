// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit

/*
 * Principal Coordinates Analysis (PCoA), a.k.a. classical
 * multidimensional scaling.
 *
 * With A = D^2 (elementwise) and J = I - 11'/n, the
 * coordinates come from the top eigenpairs of
 *
 *     B = -1/2 * J * A * J
 *
 * Only the k largest (algebraic) eigenpairs are wanted, with
 * k small. So instead of forming B (n x n) and decomposing it
 * in full like cmdscale(), they are found by block Lanczos
 * with thick restarts (lanczos.c), which only needs products
 * B * X. Those are computed straight from the packed lower
 * triangle of D that R's `dist` objects hold (gower.c, shared
 * with dbrda()).
 *
 * Every basis vector is kept orthogonal to the ones vector, so
 * J * X = X, and B * X = -1/2 * J * (A * X): one product with
 * A, then subtracting each column's mean.
 *
 * Determinism: results are bit-identical for any number of
 * threads (see gower.c). The starting vectors come from a
 * fixed PCG stream.
 *
 * Sign convention: each axis is flipped so that its coordinate
 * with the largest absolute value is positive (the first such
 * sample, allowing for rounding; see lz_sign_index()).
 */

#include "ecodive.h"


//======================================================
// R interface. `sexp_dist` is a double vector holding a
// dist object's packed lower triangle for `sexp_n`
// samples (n >= 2). Returns list(points, eig, trace,
// steps, converged, resid), with points an n x k matrix of
// eigenvectors scaled by sqrt(eig), or 0 where eig <= 0.
//======================================================
SEXP C_pcoa(SEXP sexp_dist, SEXP sexp_n, SEXP sexp_k, SEXP sexp_n_threads) {

  int n         = asInteger(sexp_n);
  int k         = asInteger(sexp_k);
  int n_threads = asInteger(sexp_n_threads);
  int dims      = n - 1; // B's rank is at most n - 1.

  if (n < 2 || k < 1 || k > dims || n_threads < 1 ||
      XLENGTH(sexp_dist) != (R_xlen_t) n * (n - 1) / 2)
    error("Invalid arguments to C_pcoa()."); // # nocov


  //--------------------------------------------------
  // trace(B) = (1/n) * sum_{i<j} d_ij^2, and products
  // with B for up to Lanczos's block size of vectors at
  // a time.
  //--------------------------------------------------
  gower_t g;
  double trace = gower_init(&g, REAL(sexp_dist), n, (k > 16) ? k : 16, n_threads);


  //--------------------------------------------------
  // The top k eigenpairs of B.
  //--------------------------------------------------
  SEXP sexp_points = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP sexp_eig    = PROTECT(allocVector(REALSXP, k));
  SEXP sexp_res    = PROTECT(allocVector(REALSXP, k));
  double *points   = REAL(sexp_points);
  double *eig      = REAL(sexp_eig);

  lanczos_t lz = {
    .n = n, .k = k, .dims = dims, .deflate = 1, .scale = trace,
    .op = gower_mult, .ctx = &g,
    .vecs = points, .vals = eig, .resid = REAL(sexp_res) };

  lanczos(&lz);


  //--------------------------------------------------
  // Coordinates, with each axis's sign fixed.
  //--------------------------------------------------
  for (int j = 0; j < k; j++) {

    double *x = points + (size_t) n * j;

    int imax = lz_sign_index(x, n);

    double mult = (eig[j] > 0) ? sqrt(eig[j]) : 0;
    if (x[imax] < 0) mult = -mult;

    for (int i = 0; i < n; i++) x[i] *= mult;
  }

  const char *names[] = { "points", "eig", "trace", "steps", "converged", "resid", "" };
  SEXP sexp_result = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(sexp_result, 0, sexp_points);
  SET_VECTOR_ELT(sexp_result, 1, sexp_eig);
  SET_VECTOR_ELT(sexp_result, 2, ScalarReal(trace));
  SET_VECTOR_ELT(sexp_result, 3, ScalarInteger(lz.steps));
  SET_VECTOR_ELT(sexp_result, 4, ScalarLogical(lz.converged));
  SET_VECTOR_ELT(sexp_result, 5, sexp_res);

  UNPROTECT(4);
  return sexp_result;
}
