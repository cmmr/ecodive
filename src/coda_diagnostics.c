// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// zero_pattern(): which features are zero in each sample, grouped into
// patterns. A sample's pattern is the set of its zero features, or
// equally the list of its non-zero ones, which the sparse table already
// holds (each sample's entries in feature order, zeros dropped). Samples
// are sorted by that list (its length, then its features, then the sample
// itself, a total order), so equal lists are adjacent; each run of equal
// lists is one pattern. Nothing is densified but the patterns themselves.

#include "ecodive.h"
#include <stdlib.h>   // qsort

static const int *zp_pos;   // per sample: start of its non-zero features
static const int *zp_otu;   // the features, increasing within a sample


// Compares two samples' lists of non-zero features.
static int zp_list_cmp(int a, int b) {
  int la = zp_pos[a + 1] - zp_pos[a];
  int lb = zp_pos[b + 1] - zp_pos[b];
  if (la != lb) return la < lb ? -1 : 1;
  const int *x = zp_otu + zp_pos[a];
  const int *y = zp_otu + zp_pos[b];
  for (int k = 0; k < la; k++)
    if (x[k] != y[k]) return x[k] < y[k] ? -1 : 1;
  return 0;
}


static int zp_cmp(const void *p, const void *q) {
  int a = *(const int *) p, b = *(const int *) q;
  int c = zp_list_cmp(a, b);
  if (c) return c;
  return a < b ? -1 : a > b;
}


//======================================================
// R interface for zero_pattern().
//
// Returns a list:
//   group     per sample (n): its pattern, 1-based, the
//             patterns numbered by their first sample.
//   n_nonzero per feature (D): samples where it is not 0.
//   patterns  logical, patterns x features: TRUE = zero.
//======================================================
SEXP C_zero_pattern(SEXP sexp_counts, SEXP sexp_margin) {

  init_n_ptrs(10);
  ecomatrix_t *em = new_ecomatrix(sexp_counts, sexp_margin);

  int n = em->n_samples;
  int D = em->n_otus;

  zp_pos = em->pos_vec;
  zp_otu = em->otu_vec;

  int *idx  = (int *) R_alloc((size_t) n + 1, sizeof(int));
  int *run  = (int *) R_alloc((size_t) n + 1, sizeof(int));   // run of each sample
  int *rmin = (int *) R_alloc((size_t) n + 1, sizeof(int));   // a run's first sample
  for (int s = 0; s < n; s++) idx[s] = s;
  if (n > 1) qsort(idx, (size_t) n, sizeof(int), zp_cmp);

  // Runs of equal lists, in sorted order; each run's first (smallest)
  // sample is its first element, as ties are broken by sample.
  int n_runs = 0;
  for (int k = 0; k < n; k++) {
    if (k == 0 || zp_list_cmp(idx[k - 1], idx[k]) != 0) rmin[n_runs++] = idx[k];
    run[idx[k]] = n_runs - 1;
  }

  // Number the patterns by their first sample.
  int *id = (int *) R_alloc((size_t) n_runs + 1, sizeof(int));
  for (int r = 0; r < n_runs; r++) id[r] = -1;
  int next = 0;
  int *first = (int *) R_alloc((size_t) n_runs + 1, sizeof(int));

  SEXP sexp_group = PROTECT(allocVector(INTSXP, n));
  int *group = INTEGER(sexp_group);
  for (int s = 0; s < n; s++) {
    int r = run[s];
    if (id[r] < 0) { id[r] = next; first[next] = s; next++; }
    group[s] = id[r] + 1;
  }

  SEXP sexp_nnz = PROTECT(allocVector(INTSXP, D));
  int *nnz = INTEGER(sexp_nnz);
  for (int f = 0; f < D; f++) nnz[f] = 0;
  for (int e = 0; e < em->pos_vec[n]; e++) nnz[em->otu_vec[e]]++;

  SEXP sexp_pat = PROTECT(allocMatrix(LGLSXP, n_runs, D));
  int *pat = LOGICAL(sexp_pat);
  for (size_t k = 0; k < (size_t) n_runs * D; k++) pat[k] = 1;
  for (int p = 0; p < n_runs; p++) {
    int s = first[p];
    for (int e = em->pos_vec[s]; e < em->pos_vec[s + 1]; e++)
      pat[p + (size_t) n_runs * em->otu_vec[e]] = 0;
  }

  SEXP result = PROTECT(allocVector(VECSXP, 3));
  SEXP names  = PROTECT(allocVector(STRSXP, 3));
  SET_VECTOR_ELT(result, 0, sexp_group);
  SET_VECTOR_ELT(result, 1, sexp_nnz);
  SET_VECTOR_ELT(result, 2, sexp_pat);
  SET_STRING_ELT(names, 0, mkChar("group"));
  SET_STRING_ELT(names, 1, mkChar("n_nonzero"));
  SET_STRING_ELT(names, 2, mkChar("patterns"));
  setAttrib(result, R_NamesSymbol, names);

  free_all();
  UNPROTECT(5);
  return result;
}
