// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include "ecodive.h"


extern SEXP C_alpha_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_blas_settable(void);
extern SEXP C_blas_threads(SEXP);
extern SEXP C_dbrda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nmds_prep(SEXP, SEXP);
extern SEXP C_nmds_run(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nmds_scale(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pcoa(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_permanova(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pthreads(void);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_umap_graph(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_rand(SEXP, SEXP, SEXP);
extern SEXP C_umap_sgd(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_spectral(SEXP, SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern void init_blas_threads(void);


static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 6},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  8},
  {"C_blas_settable", (DL_FUNC) &C_blas_settable, 0},
  {"C_blas_threads",  (DL_FUNC) &C_blas_threads,  1},
  {"C_dbrda",     (DL_FUNC) &C_dbrda,     7},
  {"C_nmds_prep", (DL_FUNC) &C_nmds_prep, 2},
  {"C_nmds_run",  (DL_FUNC) &C_nmds_run,  6},
  {"C_nmds_scale", (DL_FUNC) &C_nmds_scale, 4},
  {"C_pca",       (DL_FUNC) &C_pca,       7},
  {"C_pcoa",      (DL_FUNC) &C_pcoa,      4},
  {"C_permanova", (DL_FUNC) &C_permanova, 9},
  {"C_pthreads",  (DL_FUNC) &C_pthreads,  0},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    5},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_umap_graph",    (DL_FUNC) &C_umap_graph,    5},
  {"C_umap_rand",     (DL_FUNC) &C_umap_rand,     3},
  {"C_umap_sgd",      (DL_FUNC) &C_umap_sgd,      5},
  {"C_umap_spectral", (DL_FUNC) &C_umap_spectral, 3},
  {"C_unifrac",   (DL_FUNC) &C_unifrac,   7},
  {NULL, NULL, 0}
};


// # nocov start
void R_init_ecodive(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  init_blas_threads();
}
// # nocov end
