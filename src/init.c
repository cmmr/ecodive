// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include "ecodive.h"


extern SEXP C_alpha_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_anosim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_blas_settable(void);
extern SEXP C_blas_threads(SEXP);
extern SEXP C_cca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dbrda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_envfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_isomap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mantel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mrpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nmds(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_normalize(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_plsda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pcoa(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pcoa_correct(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pcoa_biplot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_permanova(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_permdisp(SEXP, SEXP, SEXP);
extern SEXP C_protest(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_permdisp_perm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pair_seed(SEXP, SEXP);
extern SEXP C_pthreads(void);
extern SEXP C_stop_thread_pool(void);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rpca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_simper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_tsne_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_tsne_run(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_ab(SEXP, SEXP);
extern SEXP C_umap_graph(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_rand(SEXP, SEXP, SEXP);
extern SEXP C_umap_sgd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_spectral(SEXP, SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern void init_blas_threads(void);


static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 6},
  {"C_anosim",    (DL_FUNC) &C_anosim,    6},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  8},
  {"C_blas_settable", (DL_FUNC) &C_blas_settable, 0},
  {"C_blas_threads",  (DL_FUNC) &C_blas_threads,  1},
  {"C_cca",       (DL_FUNC) &C_cca,       9},
  {"C_dbrda",     (DL_FUNC) &C_dbrda,     8},
  {"C_envfit",    (DL_FUNC) &C_envfit,   12},
  {"C_isomap",    (DL_FUNC) &C_isomap,    7},
  {"C_mantel",    (DL_FUNC) &C_mantel,    7},
  {"C_mrpp",      (DL_FUNC) &C_mrpp,      7},
  {"C_nmds",      (DL_FUNC) &C_nmds,      6},
  {"C_normalize", (DL_FUNC) &C_normalize, 6},
  {"C_pca",       (DL_FUNC) &C_pca,       7},
  {"C_plsda",     (DL_FUNC) &C_plsda,     8},
  {"C_rda",       (DL_FUNC) &C_rda,       9},
  {"C_pcoa",      (DL_FUNC) &C_pcoa,      4},
  {"C_pcoa_correct", (DL_FUNC) &C_pcoa_correct, 5},
  {"C_pcoa_biplot",  (DL_FUNC) &C_pcoa_biplot,  7},
  {"C_permanova", (DL_FUNC) &C_permanova, 9},
  {"C_permdisp",  (DL_FUNC) &C_permdisp,  3},
  {"C_permdisp_perm", (DL_FUNC) &C_permdisp_perm, 6},
  {"C_protest",   (DL_FUNC) &C_protest,   5},
  {"C_pair_seed",     (DL_FUNC) &C_pair_seed,     2},
  {"C_pthreads",  (DL_FUNC) &C_pthreads,  0},
  {"C_stop_thread_pool", (DL_FUNC) &C_stop_thread_pool, 0},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    5},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_rpca",      (DL_FUNC) &C_rpca,     12},
  {"C_simper",    (DL_FUNC) &C_simper,    7},
  {"C_tsne_p",    (DL_FUNC) &C_tsne_p,    7},
  {"C_tsne_run",  (DL_FUNC) &C_tsne_run,  4},
  {"C_umap_ab",       (DL_FUNC) &C_umap_ab,       2},
  {"C_umap_graph",    (DL_FUNC) &C_umap_graph,    5},
  {"C_umap_rand",     (DL_FUNC) &C_umap_rand,     3},
  {"C_umap_sgd",      (DL_FUNC) &C_umap_sgd,      6},
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
