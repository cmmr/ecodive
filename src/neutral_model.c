// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit
//
// The per-feature sums behind neutral_model() (R/neutral_model.r): for each
// group of samples and each feature, its total count and the number of the
// group's samples it is present in. Everything else (the fit of the
// immigration probability m, the intervals and the alternative models) is
// done in R over these K-length vectors.
//
// One pass over the table's stored values, in feature order (da_table_init()
// in da_tests.c, with no normalization). Each feature's sums run over its
// samples in order, so the results are the same for any number of threads.

#include "ecodive.h"

#define NORM_NONE 0


static const da_table_t *nm_tab;
static const int        *nm_group;  // per sample: 0 .. G - 1, or -1 (left out)
static int               nm_G;
static double           *nm_total;  // G x D
static int              *nm_occ;    // G x D


static void *nm_worker (void *arg) {

  int thread_i  = ((worker_t *) arg)->i;
  int n_threads = ((worker_t *) arg)->n;

  const da_table_t *tab = nm_tab;
  int D = tab->D, G = nm_G;

  for (int f = thread_i; f < D; f += n_threads) {

    size_t  o     = (size_t) f * G;
    double *total = nm_total + o;
    int    *occ   = nm_occ   + o;

    for (int g = 0; g < G; g++) { total[g] = 0; occ[g] = 0; }

    for (int e = tab->fpos[f]; e < tab->fpos[f + 1]; e++) {
      int    g = nm_group[tab->fsam[e]];
      double v = tab->fval[e];
      if (g < 0) continue;
      total[g] += v;
      if (v > 0) occ[g]++;
    }
  }

  return NULL;
}



//======================================================
// R interface. For a count table of n samples and D
// features (`margin` as in ecomatrix.c), `group` (one
// integer per sample: 0 .. G - 1, or -1 to leave the
// sample out) and G:
//
// Returns list(total, occ): G x D matrices of each
// feature's total count, and the number of samples it is
// present in (count > 0), in each group.
//======================================================
SEXP C_neutral_stats(
    SEXP sexp_counts, SEXP sexp_margin, SEXP sexp_group, SEXP sexp_n_groups,
    SEXP sexp_n_threads ) {

  int G     = asInteger(sexp_n_groups);
  int n_thr = asInteger(sexp_n_threads);
  if (n_thr < 1) n_thr = 1;
  if (G < 1) error("No groups.");  // # nocov

  da_table_t tab;
  da_table_init(sexp_counts, sexp_margin, NORM_NONE, 0, n_thr, &tab);
  int n = tab.n, D = tab.D;

  if (LENGTH(sexp_group) != n)
    error("`group` must have one value per sample.");  // # nocov

  SEXP sexp_total = PROTECT(allocMatrix(REALSXP, G, D));
  SEXP sexp_occ   = PROTECT(allocMatrix(INTSXP,  G, D));

  // Threads only for enough work: about 10 ns per stored value.
  int nnz = tab.fpos[D];
  int nt  = threads_for(n_thr, ((double) nnz + D) * 10, 1e6);

  nm_tab   = &tab;
  nm_group = INTEGER(sexp_group);
  nm_G     = G;
  nm_total = REAL(sexp_total);
  nm_occ   = INTEGER(sexp_occ);

  run_parallel(nm_worker, nt, D < 100 ? 100 : D);

  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, sexp_total);
  SET_VECTOR_ELT(result, 1, sexp_occ);

  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("total"));
  SET_STRING_ELT(names, 1, mkChar("occ"));
  setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(4);
  return result;
}
