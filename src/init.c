// Copyright (c) 2026 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include "ecodive.h"


extern SEXP C_alpha_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancombc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancombc_global(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancombc_dunnett(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancombc_trend(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ancombc_nonzero(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_anosim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_balance(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_bioenv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_blas_settable(void);
extern SEXP C_blas_threads(SEXP);
extern SEXP C_capda_trace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cpu_info(SEXP);
extern SEXP C_cord_cross(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dbrda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dpcoa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_decontam(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_da_lm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_da_tests(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_feature_lm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_feature_lmm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dirmult_ttest(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_envfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_isomap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_helmert_basis(SEXP);
extern SEXP C_lda(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_linda_prep(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_linda_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_linda_mode(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_lda_loo(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_log_ratio(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_log_ratio_coords(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cross_assoc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mantel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mrpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mult_default_delta(SEXP);
extern SEXP C_nmds(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_normalize(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_calinski(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pam(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_silhouette(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_philr_coords(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_philr_stats(SEXP, SEXP, SEXP);
extern SEXP C_philr_tree(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_proportionality(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_oplsda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
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
extern SEXP C_qmp_scale(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rarefaction_curve(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_coverage_depth(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rpca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rpca_phylo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_simper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sparcc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_network_paths(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_network_greedy(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_network_barrat(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_network_eigen(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_perm_indices(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_medoid_accum(SEXP, SEXP, SEXP);
extern SEXP C_row_medians(SEXP, SEXP);
extern SEXP C_replace_zeros_mult(SEXP, SEXP, SEXP);
extern SEXP C_rng_dist(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_st_collapse(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_st_cond_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sourcetracker(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_tsne_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_tsne_run(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_ab(SEXP, SEXP);
extern SEXP C_umap_graph(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_rand(SEXP, SEXP, SEXP);
extern SEXP C_umap_sgd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_umap_spectral(SEXP, SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_zero_counts(SEXP, SEXP);

extern void init_blas_threads(void);
extern void init_cpu_info(void);


static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 6},
  {"C_ancom",     (DL_FUNC) &C_ancom,     10},
  {"C_ancombc",   (DL_FUNC) &C_ancombc,   12},
  {"C_ancombc_global",  (DL_FUNC) &C_ancombc_global,   4},
  {"C_ancombc_dunnett", (DL_FUNC) &C_ancombc_dunnett,  6},
  {"C_ancombc_trend",   (DL_FUNC) &C_ancombc_trend,   10},
  {"C_ancombc_nonzero", (DL_FUNC) &C_ancombc_nonzero,  4},
  {"C_anosim",    (DL_FUNC) &C_anosim,    6},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  8},
  {"C_balance",   (DL_FUNC) &C_balance,   5},
  {"C_bioenv",    (DL_FUNC) &C_bioenv,    5},
  {"C_blas_settable", (DL_FUNC) &C_blas_settable, 0},
  {"C_blas_threads",  (DL_FUNC) &C_blas_threads,  1},
  {"C_capda_trace", (DL_FUNC) &C_capda_trace, 6},
  {"C_cca",       (DL_FUNC) &C_cca,      10},
  {"C_cpu_info",  (DL_FUNC) &C_cpu_info,  1},
  {"C_cord_cross", (DL_FUNC) &C_cord_cross, 8},
  {"C_dbrda",     (DL_FUNC) &C_dbrda,     8},
  {"C_dpcoa",     (DL_FUNC) &C_dpcoa,     8},
  {"C_decontam",  (DL_FUNC) &C_decontam,  8},
  {"C_da_lm",         (DL_FUNC) &C_da_lm,         8},
  {"C_da_tests",      (DL_FUNC) &C_da_tests,      7},
  {"C_feature_lm",    (DL_FUNC) &C_feature_lm,    4},
  {"C_feature_lmm",   (DL_FUNC) &C_feature_lmm,   4},
  {"C_dirmult_ttest", (DL_FUNC) &C_dirmult_ttest, 11},
  {"C_envfit",    (DL_FUNC) &C_envfit,   12},
  {"C_isomap",    (DL_FUNC) &C_isomap,    7},
  {"C_lda",       (DL_FUNC) &C_lda,       5},
  {"C_linda_prep", (DL_FUNC) &C_linda_prep, 5},
  {"C_linda_fit",  (DL_FUNC) &C_linda_fit, 12},
  {"C_linda_mode", (DL_FUNC) &C_linda_mode, 5},
  {"C_helmert_basis", (DL_FUNC) &C_helmert_basis, 1},
  {"C_lda_loo",   (DL_FUNC) &C_lda_loo,   5},
  {"C_log_ratio", (DL_FUNC) &C_log_ratio, 5},
  {"C_log_ratio_coords", (DL_FUNC) &C_log_ratio_coords, 6},
  {"C_cross_assoc", (DL_FUNC) &C_cross_assoc, 6},
  {"C_mantel",    (DL_FUNC) &C_mantel,    7},
  {"C_mrpp",      (DL_FUNC) &C_mrpp,      7},
  {"C_mult_default_delta", (DL_FUNC) &C_mult_default_delta, 1},
  {"C_nmds",      (DL_FUNC) &C_nmds,      6},
  {"C_normalize", (DL_FUNC) &C_normalize, 6},
  {"C_calinski",  (DL_FUNC) &C_calinski,  6},
  {"C_pam",       (DL_FUNC) &C_pam,       4},
  {"C_silhouette", (DL_FUNC) &C_silhouette, 5},
  {"C_pca",       (DL_FUNC) &C_pca,       7},
  {"C_philr_coords", (DL_FUNC) &C_philr_coords, 9},
  {"C_philr_stats",  (DL_FUNC) &C_philr_stats,  3},
  {"C_philr_tree",   (DL_FUNC) &C_philr_tree,   4},
  {"C_proportionality", (DL_FUNC) &C_proportionality, 8},
  {"C_oplsda",    (DL_FUNC) &C_oplsda,    8},
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
  {"C_qmp_scale", (DL_FUNC) &C_qmp_scale, 5},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    5},
  {"C_rarefaction_curve", (DL_FUNC) &C_rarefaction_curve, 5},
  {"C_coverage_depth",    (DL_FUNC) &C_coverage_depth,    4},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_medoid_accum", (DL_FUNC) &C_medoid_accum, 3},
  {"C_row_medians",  (DL_FUNC) &C_row_medians,  2},
  {"C_replace_zeros_mult", (DL_FUNC) &C_replace_zeros_mult, 3},
  {"C_rng_dist",  (DL_FUNC) &C_rng_dist,  5},
  {"C_st_collapse",   (DL_FUNC) &C_st_collapse,   4},
  {"C_st_cond_prob",  (DL_FUNC) &C_st_cond_prob,  9},
  {"C_sourcetracker", (DL_FUNC) &C_sourcetracker, 16},
  {"C_rpca",      (DL_FUNC) &C_rpca,     12},
  {"C_rpca_phylo", (DL_FUNC) &C_rpca_phylo, 18},
  {"C_simper",    (DL_FUNC) &C_simper,    7},
  {"C_sparcc",    (DL_FUNC) &C_sparcc,   11},
  {"C_network_paths",  (DL_FUNC) &C_network_paths,  5},
  {"C_network_greedy", (DL_FUNC) &C_network_greedy, 4},
  {"C_network_barrat", (DL_FUNC) &C_network_barrat, 4},
  {"C_network_eigen",  (DL_FUNC) &C_network_eigen,  5},
  {"C_perm_indices",   (DL_FUNC) &C_perm_indices,   4},
  {"C_tsne_p",    (DL_FUNC) &C_tsne_p,    7},
  {"C_tsne_run",  (DL_FUNC) &C_tsne_run,  4},
  {"C_umap_ab",       (DL_FUNC) &C_umap_ab,       2},
  {"C_umap_graph",    (DL_FUNC) &C_umap_graph,    5},
  {"C_umap_rand",     (DL_FUNC) &C_umap_rand,     3},
  {"C_umap_sgd",      (DL_FUNC) &C_umap_sgd,      6},
  {"C_umap_spectral", (DL_FUNC) &C_umap_spectral, 3},
  {"C_unifrac",   (DL_FUNC) &C_unifrac,   7},
  {"C_zero_counts", (DL_FUNC) &C_zero_counts, 2},
  {NULL, NULL, 0}
};


// # nocov start
void R_init_ecodive(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  init_blas_threads();
  init_cpu_info();
}
// # nocov end
