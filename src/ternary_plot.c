// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The sums behind ternary_plot(): each sample's three part sums (samples
// mode), or each feature's mean relative abundance in each of three sample
// groups (features mode). Read sparsely, through new_ecomatrix().
// Single-threaded, with every sum in the samples' and the stored values'
// order, so the results do not depend on anything but the input. Written
// from the formulas in plans/ternary_plot.md (clean room).

#include "ecodive.h"


//======================================================
// R interface: each sample's sum over the features of
// each of the three parts.
//
// sexp_part: for each feature, its 0-based part (0, 1 or
//            2), or -1 for none.
//
// Returns an n x 3 matrix.
//======================================================
SEXP C_ternary_sums (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_part) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int n = em->n_samples;
  const int *part = INTEGER(sexp_part);

  SEXP sexp_result = PROTECT(allocMatrix(REALSXP, n, 3));
  double *out = REAL(sexp_result);

  for (int s = 0; s < n; s++) {
    double sum[3] = {0, 0, 0};
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
      int k = part[em->otu_vec[e]];
      if (k >= 0) sum[k] += em->val_vec[e];
    }
    for (int k = 0; k < 3; k++) out[s + (size_t) n * k] = sum[k];
  }

  free_all();
  UNPROTECT(1);
  return sexp_result;
}



//======================================================
// R interface: each feature's mean relative abundance in
// each of three groups of samples, and over all of them.
// Each sample is closed over every feature first; samples
// whose total is not above 0 are left out.
//
// sexp_group: for each sample, its 0-based group (0, 1 or
//             2), or -1 for none.
//
// Returns list(means (p x 3), overall (p), totals (n),
// sizes (3): the samples used in each group).
//======================================================
SEXP C_ternary_means (SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_group) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);
  int n = em->n_samples, p = em->n_otus;
  const int *group = INTEGER(sexp_group);

  SEXP sexp_means   = PROTECT(allocMatrix(REALSXP, p, 3));
  SEXP sexp_overall = PROTECT(allocVector(REALSXP, p));
  SEXP sexp_totals  = PROTECT(allocVector(REALSXP, n));
  SEXP sexp_sizes   = PROTECT(allocVector(INTSXP, 3));
  double *means   = REAL(sexp_means);
  double *overall = REAL(sexp_overall);
  double *totals  = REAL(sexp_totals);
  int    *sizes   = INTEGER(sexp_sizes);

  for (size_t i = 0; i < (size_t) p * 3; i++) means[i] = 0;
  for (int j = 0; j < p; j++) overall[j] = 0;
  for (int k = 0; k < 3; k++) sizes[k] = 0;

  for (int s = 0; s < n; s++) {
    double total = 0;
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) total += em->val_vec[e];
    totals[s] = total;
    int k = group[s];
    if (k < 0 || !(total > 0)) continue;
    sizes[k]++;
    double *col = means + (size_t) p * k;
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++) {
      double v = em->val_vec[e] / total;
      col[em->otu_vec[e]] += v;
      overall[em->otu_vec[e]] += v;
    }
  }

  int used = sizes[0] + sizes[1] + sizes[2];
  for (int k = 0; k < 3; k++) {
    double *col = means + (size_t) p * k;
    for (int j = 0; j < p; j++) col[j] = sizes[k] ? col[j] / sizes[k] : R_NaN;
  }
  for (int j = 0; j < p; j++) overall[j] = used ? overall[j] / used : R_NaN;

  free_all();

  SEXP sexp_result = PROTECT(allocVector(VECSXP, 4));
  SEXP sexp_names  = PROTECT(allocVector(STRSXP, 4));
  SET_VECTOR_ELT(sexp_result, 0, sexp_means);   SET_STRING_ELT(sexp_names, 0, mkChar("means"));
  SET_VECTOR_ELT(sexp_result, 1, sexp_overall); SET_STRING_ELT(sexp_names, 1, mkChar("overall"));
  SET_VECTOR_ELT(sexp_result, 2, sexp_totals);  SET_STRING_ELT(sexp_names, 2, mkChar("totals"));
  SET_VECTOR_ELT(sexp_result, 3, sexp_sizes);   SET_STRING_ELT(sexp_names, 3, mkChar("sizes"));
  setAttrib(sexp_result, R_NamesSymbol, sexp_names);

  UNPROTECT(6);
  return sexp_result;
}
