# To-do: new methods, by priority

Work top to bottom. Each item links to its plan; read [README.md](README.md)
first. Tick an item, and add the commit, when its method is merged with tests
and docs.

**Order (renumbered 2026-11-20, at the user's request):**

1. **Completed items come first,** in the order they were built (items
   1–36).
2. **Open items follow, ranked by ease of implementation and demand in the
   field.** Each open item has a demand score and an effort score (1–5; see
   "Ranking inputs" below). Items are ranked by demand minus effort; ties go
   to the higher demand, then the lower effort.
   - **Demand** comes from OpenAlex citations 2022–26 of the method's paper
     and Europe PMC microbiome full-text mentions 2022–26, retrieved
     2026-11-19 (earlier plans) and 2026-11-20 (the new ones).
   - **Effort** counts clean-room work as harder than a port of permissive
     code (see the 2026-11-18 rules below), and counts licence care
     (non-commercial references) as extra effort.
3. **Dependencies win over score:** an item never comes before an item it
   builds on. For example, `alpha_metrics` comes before
   `rarefaction_curve()`, which uses its `q`, and `moran_i()` before
   `dbmem()`.
4. **Superseded methods go last** (`dca()`, `largevis()`).

The old numbers still appear in commit messages ("plans: tick item 40"),
so "Renumbering map" at the bottom gives each old number's new one. Every
plan's **Priority** line and every "item N" cross-reference in `plans/` was
updated to the new numbers. `new_methods.md` is the exception: it is
historical and uses the survey's own numbering.

| # | Do | Plan | Tier | Why here |
|--:|:---|:-----|:-----|:---------|
| 1 | [x] Split `mantel_partial()` out of `mantel(given =)` (3074a20) | [mantel_partial.md](mantel_partial.md) | housekeeping | It is a distinct named test (vegan's `mantel.partial()`), and the API is still unreleased, so splitting it now is cheap. |
| 2 | [x] Pairwise tests: `pairwise` and `p_adjust` on `permanova()`, `anosim()`, `mrpp()` and `permdisp()` (575da24) | [pairwise_tests.md](pairwise_tests.md) | port (q2-diversity, BSD-3) | In nearly every microbiome paper; simple, built on existing tests. |
| 3 | [x] `rda()`, including partial RDA (46db60b) | [rda.md](rda.md) | port (skbio, BSD-3) + ecodive's dbrda core | High demand (tb-RDA, phyloseq "RDA"). Several later items build on its core. |
| 4 | [x] `cca()` (e33c7e1) | [cca_ca.md](cca_ca.md) | port (skbio) | High demand in ecology; phyloseq "CCA". |
| 5 | [x] `ca()` (80c8a86) | [cca_ca.md](cca_ca.md) | port (skbio) | Cheap once `cca()`'s weighted core exists. |
| 6 | [x] `procrustes()` and `procrustes_test()` (a1c2a66) | [procrustes.md](procrustes.md) | port (SciPy, BSD-3) | In QIIME 2; simple. |
| 7 | [x] `plsda()` (b3c81a4) | [plsda_oplsda.md](plsda_oplsda.md) | port (scikit-learn, BSD-3) | Popular in microbiome work through mixOmics. |
| 8 | [x] `simper()` (dbf0c2b) | [simper.md](simper.md) | port (ecopy, MIT) + Clarke's permutation p-values | High demand in community ecology; simple. |
| 9 | [x] `bioenv()` (96381a6) | [bioenv.md](bioenv.md) | port (skbio) | In QIIME 2; simple. |
| 10 | [x] `lda()` (7323057) | [lda_capda.md](lda_capda.md) | port (scikit-learn) | Simple; `capda()` needs it. |
| 11 | [x] `oplsda()` (4c7e646) | [plsda_oplsda.md](plsda_oplsda.md) | port (pyopls, MIT) | Metabolomics demand; small on top of `plsda()`. |
| 12 | [x] `varpart()` (eb06b58) | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | High demand in ecology; formulas only. |
| 13 | [x] `forward_select()` (d0b6339) | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | Usually used alongside varpart. |
| 14 | [x] `alr()`, `ilr()` (51fd7d0) | [log_ratio_transforms.md](log_ratio_transforms.md) | port (scikit-bio, BSD-3) | Simplest compositional gap; `balance()` and `philr()` build on the ILR code. |
| 15 | [x] `replace_zeros()` wrapper (+ `list_zero_replacements()`, `match_zero_replacement()`) and `replace_zeros_mult()` (7d84f63) | [zero_replacement.md](zero_replacement.md) | port (scikit-bio, BSD-3) | Principled zeros for every log-ratio method; the wrapper is ready for the Bayesian methods (item 55). |
| 16 | [x] `log_ratio()`, `balance()` (c56ea2a) | [log_ratio.md](log_ratio.md) | port (Qurro / scikit-bio, BSD-3) | Trivial; log ratios are what current compositional guidance recommends analysing. |
| 17 | [x] `vlr()`, `rho_p()`, `phi_s()` (60d304e) | [proportionality.md](proportionality.md) | port (scikit-bio `vlr`) + formula (propr output-only) | High demand (propr); all-pairs C over features is ecodive's strength. |
| 18 | [x] `qmp_scale()`, `qmp_rarefy()`; per-sample `rarefy(depth =)` (068453a) | [qmp.md](qmp.md) | port (raeslab QMP.R, BSD-3) + formula | Simple, over `rarefy()`'s C; the absolute-abundance side of the debate. |
| 19 | [x] `ancom()` (63822a8) | [ancom.md](ancom.md) | port (scikit-bio, BSD-3) | In QIIME 2 and widely used; D² independent tests suit threaded C. |
| 20 | [x] `dirmult_ttest()` (2291a31), then its scale model (c88b95b) | [dirmult_ttest.md](dirmult_ttest.md) | port (scikit-bio, BSD-3; scale model from ALDEx3, MIT) | ALDEx2-style effect sizes, high demand; needs a PCG Gamma sampler. |
| 21 | [x] `rpca_phylo()` (4b259ba) | [rpca_phylo.md](rpca_phylo.md) | port (gemelli, BSD-3) on `rpca()` | QIIME 2's phylogenetic RPCA; reuses `rpca()` and the UniFrac tree walk. Phylogeny-aware. |
| 22 | [x] `ancombc()`, `ancombc2()` (49a6e59) | [ancombc.md](ancombc.md) | port (scikit-bio, BSD-3, derived from Artistic-2.0 ANCOMBC: keep both notices) | Among the most used DA methods now; large (~1,000 lines of core). |
| 23 | [x] `ancombc_global()`, `ancombc_pairwise()`, `ancombc_dunnett()`, `ancombc_trend()`, `structural_zeros()` (7d9b859) | [ancombc.md](ancombc.md) | port (scikit-bio + ANCOMBC notices), as above | Post-hoc tests on a fitted `ancombc()` / `ancombc2()`; `ancombc_trend()` needs a small QP solver. |
| 24 | [x] `capda()` (06f2fbf) | [lda_capda.md](lda_capda.md) | derived (`pcoa()` + `lda()`) | Anderson's CAP, popular in marine ecology. |
| 25 | [x] `philr()`, `philr_distance()` (5f8d7bb) | [philr.md](philr.md) | **clean room** (philr GPL-3; formulas from the open-access paper) | First among the clean-room methods: phylogeny-aware, high demand, and an O(D) formula per sample. |
| 26 | [x] `dpcoa()` (6ed9799) | [dpcoa.md](dpcoa.md) | **clean room** | Second among the clean-room methods, after `philr()`: phyloseq "DPCoA", phylogeny-aware. |
| 27 | [x] Documentation: Schloss (2026) and depth in the RPCA family; `scale_sd` sensitivity; ML caveats; predicted loads in `qmp_scale()`; `log_base` on DA results (60b71df) | [compositional_docs.md](compositional_docs.md) | docs | Almost no code (one `log_base` attribute), and it corrects how existing functions are presented: Schloss's 2026 rebuttal of robust Aitchison PCA, scale-assumption FDR results. |
| 28 | [x] Effect (MW), MA and forest plots for DA results (`plot(type = )`), with new dispersion and mean-abundance columns in `dirmult_ttest()` (48cf7ad) | [da_plots.md](da_plots.md) | port (ALDEx3's effect statistics, MIT; ALDEx2 output-only) | Low effort on results users already have; Gloor et al. (2016) show volcano plots, ecodive's only DA display, mislead. |
| 29 | [x] `da_test()` wrapper (+ `list_da_tests()`, `match_da_test()`); `da_wilcoxon()`, `da_ttest()`, `da_lm()` (4d0226c) | [da_test.md](da_test.md) | formula | Very low effort; the plain tests are what Wirbel (2024) and Pelto (2025) found best at controlling false discoveries; four DA methods now justify the wrapper. |
| 30 | [x] `beta_div_rarefied()`, `alpha_div_rarefied()` (963647f) | [rarefied_diversity.md](rarefied_diversity.md) | derived (on `rarefy()`) | Low effort; rarefaction is the most mentioned normalization (~3,600 papers a year) and Schloss (2024, 2026) argues for it; lets users compare with RPCA. |
| 31 | [x] `linda()` (f09370e) | [linda.md](linda.md) | clean room (LinDA GPL) | Short estimator sharing `ancombc()`'s code; 494 citations and fast-rising; most robust under strong compositional effects. |
| 32 | [x] `sparcc()` (21c3da7) | [sparcc.md](sparcc.md) | clean room (CC-BY paper) | The most-used compositional network method (1,848 citations 2022–26); builds on `vlr()` and the existing Gamma sampler. |
| 33 | [x] `assoc_network()`, `assoc_network_compare()`, `assoc_network_diff()` (a4d81fd; with `feature_assoc()`) | [assoc_network.md](assoc_network.md) | clean room (NetCoMi GPL) | Textbook formulas over `rho_p()` / `sparcc()`; NetCoMi uptake grew 3.1× from 2021–22 to 2024–25. |
| 34 | [x] `da_sensitivity()` (8d25396) | [da_sensitivity.md](da_sensitivity.md) | derived (on `da_test()`) | Reruns and plots DA calls over pseudocounts, zero treatments and `scale_sd`; the benchmarks' central concern. Needs item 29. |
| 35 | [x] `enterotype()` (20f28e1) | [community_typing.md](community_typing.md) | clean room (PAM; `cluster` is GPL) | Low effort on `jsd()`; enterotyping is among the most cited community analyses (Arumugam 2011: 3,182 citations 2022–26). |
| 36 | [x] `sourcetracker()` (12034fd) | [source_tracking.md](source_tracking.md) | port (SourceTracker2, BSD-3) | Source tracking is in steady high demand (Knights 2011: ~200 citations a year); simple; SourceTracker2's per-read Python loop is slow and its `multiprocessing` isn't reproducible, so threaded C with per-sink PCG streams is a clear win. |
| 37 | [x] `pielou()`, `goods_coverage()`, `hill()`, `renyi()`, `tsallis()` and other scikit-bio alpha metrics; `alpha_div(q =)` (9f1add2; several metrics at once dropped at the user's request) | [alpha_metrics.md](alpha_metrics.md) | port (scikit-bio, BSD-3) | Pielou in 2,001 and Good's coverage in 1,236 microbiome papers; completes QIIME 2's core alpha set; each metric is a one-line formula in the existing alpha C switch. First because items 38, 48 and 64 use it. |
| 38 | [x] `rarefaction_curve()` with `plot()` (0d65d46; also `rarefy(coverage =)`) | [rarefaction_curves.md](rarefaction_curves.md) | derived (on `rarefy()` and `alpha_div_rarefied()`'s loops) + formula (Hurlbert, Heck) | Rarefaction curves appear in 8,700 microbiome papers; closed forms for richness and coverage, resampling for other metrics; low-to-medium effort. |
| 39 | [x] `cross_assoc()`: associations between the features of two tables, with FDR and a starred heatmap (70f49b4) | [cross_assoc.md](cross_assoc.md) | formula + port (SciPy p-values, BSD-3) | Taxa-vs-environment or -metabolite Spearman heatmaps appear in 7,120 papers; a threaded p × q kernel on existing rank and Kendall code; items 56 and 84 build on it. |
| 40 | [x] `decontam()`, `decontam_not()` (5f01091) | [decontam.md](decontam.md) | port (decontam, Artistic-2.0: keep its notices) | Rising demand (3,160 citations, 478 → 751 a year; 1,617 mentions) and the Fierer et al. (2025) low-biomass guidelines; a small port of per-feature two-sample scores. |
| 41 | [x] `lefse()` (joins `da_test()`) (18c8d9d) | [lefse.md](lefse.md) | port (LEfSe, MIT; lefser Artistic-2.0) | By far the most used DA method (19,755 mentions, 8.9× MaAsLin2); its rank tests and LDA exist in ecodive; medium effort. Its weak FDR (Nearing 2022) is documented, pointing to `da_sensitivity()`. |
| 42 | [x] `dirmult_lm()` (6fd2b9e) | [dirmult_lm.md](dirmult_lm.md) | port (ALDEx3, MIT) | Extends `dirmult_ttest()` beyond two groups with covariates and load-informed scale priors; low-to-moderate effort; ALDEx2 is widely used (1,178 mentions), ALDEx3 is new. |
| 43 | [ ] Normalization `norm` choices: CSS, TMM, RLE, GMPR | [normalization_factors.md](normalization_factors.md) | port (CSS, Artistic-2.0) + clean room | Low effort each; CSS and TMM are mentioned ~120–145 times a year; they can later centre scale models. RLE is shared with item 61 (whichever lands first writes it; PyDESeq2, MIT, can settle RLE's open points). |
| 44 | [x] Biplot scalings (`plot(ord, scaling = )`) with links (811adb2) | [biplot_scaling.md](biplot_scaling.md) | formula | Low effort; the compositional biplot Gloor (2017) recommends; Aitchison & Greenacre 2002 has 830 citations. |
| 45 | [x] `valencia()` (582e456) | [community_typing.md](community_typing.md) | port (VALENCIA, MIT) | Low effort; the standard vaginal community-state-type classifier, rising (499 citations, 121 in 2026). |
| 46 | [x] Rank plots and `plot()` of `log_ratio()` / `balance()` results (Qurro-style), `log_ratio_features()` (32f2e17) | [log_ratio_plots.md](log_ratio_plots.md) | port (Qurro, BSD) | Low effort on ported code; turns DA effects and loadings into interpretable log ratios (reference frames: 668 citations). Uses item 44's `scaling` and `links`. |
| 47 | [x] `composition_barplot()` (68b12bd) | [composition_barplot.md](composition_barplot.md) | port (microshades idea, MIT) | Low effort; stacked bars are in nearly every paper (though not a compositional analysis); a better-ordered version of a ubiquitous figure. |
| 48 | [ ] `alpha_test()` wrapper (+ `ALPHA_TESTS`, `list_alpha_tests()`, `match_alpha_test()`); `alpha_wilcoxon()`, `alpha_kruskal()`, `alpha_ttest()`, `alpha_lm()` (with `(1 \| subject)`), `alpha_cor()` | [alpha_test.md](alpha_test.md) | formula + derived (on the `da_*` kernels and `linda()`'s mixed model) | The most common statistics in microbiome papers (alpha diversity with Wilcoxon 15,920, Kruskal-Wallis 12,841, mixed models 2,852 mentions); medium-to-high effort for the mixed-model contrasts; lets rbiom drop emmeans and lme4 for alpha diversity. |
| 49 | [x] `maaslin3()` (7bd94a3) (scope approved by the user, 2026-11-20) | [maaslin3.md](maaslin3.md) | port (maaslin3, MIT) | The highest-uptake DA family (MaAsLin2 2,850 citations; MaAsLin3 rising steeply); moderate effort; widens ecodive into general per-feature regression, as a one-stop shop needs. |
| 50 | [x] `indicator_species()` (9737ada) (IndVal.g, r.g; multilevel pattern analysis) | [indicator_species.md](indicator_species.md) | **clean room** (indicspecies GPL) | Solid demand (Dufrêne & Legendre 2,642 and De Cáceres & Legendre 2,494 citations; 1,487 mentions); reuses ecodive's permutation machinery and `strata`; the formulas and permutation statistic are pinned down from outputs. |
| 51 | [x] `neutral_model()` (0490b63) (Sloan's neutral community model) | [neutral_model.md](neutral_model.md) | port (neufit, BSD-2) + formula | Widely used (Sloan 2006: 1,470 citations; 531 mentions) with no standard documented package; a 1-D least-squares fit in base R plus one small C pass; low-to-medium effort. |
| 52 | [x] `beta_partition()`, `beta_partition_multi()` and turnover / nestedness beta metrics (8d5eb46) | [beta_partition.md](beta_partition.md) | formula + port (abdiv, MIT); clean room (betapart GPL) for multiple-site | Heavily cited (Baselga 2010: 2,447; betapart 1,623); short published formulas on the existing beta kernels; low-to-medium effort. |
| 53 | [x] `batch_mmuphin()`, `batch_combat()` (442ed7c) (phase 1 of batch correction) | [batch_correction.md](batch_correction.md) | port (MMUPHin, MIT; sva, Artistic-2.0) | Multi-study pooling needs batch correction (490 mentions; MMUPHin 192 citations, rising); a deterministic MIT port with ComBat's empirical-Bayes core; low-to-medium effort. Its DA meta-analysis waits for item 49. |
| 54 | [ ] `ternary_plot()` | [ternary_plot.md](ternary_plot.md) | clean room (ggtern, Ternary GPL) | The most used classic CoDa display (331 microbiome papers 2022–26); low-to-medium effort. |
| 55 | [ ] `replace_zeros_czm()`, `replace_zeros_gbm()`, `replace_zeros_sq()`, `replace_zeros_bl()` | [zero_replacement.md](zero_replacement.md) | **clean room** (zCompositions GPL; the defining paper is closed access, so formulas are confirmed from outputs) | Simple once the formulas are pinned down; `cmultRepl()`'s GBM is the field's default. |
| 56 | [x] `first_differences()`, `first_distances()`, `volatility()`, then `nmit()` (phases 1 and 2) (b701ec5) | [longitudinal.md](longitudinal.md) | port (q2-longitudinal, BSD-3) | q2-longitudinal (244 citations) has no R equivalent; phase 1 is pair bookkeeping on existing inputs, low effort; `nmit()` reuses item 39's correlation kernels and `permanova()`. |
| 57 | [ ] `principal_response_curves()` | [principal_response_curves.md](principal_response_curves.md) | derived (partial `rda()`) | Experiments over time; low effort. |
| 58 | [x] `mantel_correlogram()` (97ddfb1) | [mantel_correlogram.md](mantel_correlogram.md) | derived (`mantel()` on distance classes) | Distance decay, biogeography; low effort. |
| 59 | [ ] `moran_i()` | [moran_i.md](moran_i.md) | port (PySAL esda, BSD-3) | Simple; `dbmem()` (item 78) needs it. |
| 60 | [ ] `variation_array()`, `zero_pattern()` | [coda_diagnostic_plots.md](coda_diagnostic_plots.md) | clean room | Low effort diagnostics (Aitchison's variation array; zero patterns as in zCompositions); modest uptake. |
| 61 | [ ] `deseq2()` (joins `da_test()`; phased: Wald, then LRT, then shrinkage) | [deseq2.md](deseq2.md) | port (PyDESeq2, MIT) + clean room (LRT, R's `lowess()`; DESeq2 LGPL output-only) | Very high demand (10,234 microbiome mentions) but high effort: a large per-feature NB GLM pipeline. Benchmarks show inflated FDR on microbiome data, so it comes after the cheaper items. |
| 62 | [ ] `random_forest()` with CV, OOB, importance, ROC (added by the user, 2026-11-20) | [random_forest.md](random_forest.md) | port (scikit-learn trees, q2-sample-classifier; BSD-3) | The machine learning microbiome papers use (10,679 mentions, beside DESeq2's 10,234); an exact, threaded, seed-deterministic port of scikit-learn's trees; high effort. Whichever of items 62 and 86 lands first writes the shared MT19937 test code. |
| 63 | [ ] `feature_volatility()` (phase 3 of longitudinal) | [longitudinal.md](longitudinal.md) | port (q2-longitudinal, BSD-3) | Small once item 62 exists: q2's feature volatility is a random forest on time plus volatility charts of its top features. |
| 64 | [ ] `hill_extrapolation()`, `hill_standardize()` (iNEXT-style) | [hill_extrapolation.md](hill_extrapolation.md) | **clean room** (iNEXT GPL; formulas and bootstrap from Chao et al. 2014 and its Appendix G) | iNEXT has 3,689 citations (Chao 2014: 3,087) but is slow at microbiome depths; threaded deterministic bootstraps; coverage standardization feeds item 48; medium-to-high effort. |
| 65 | [ ] `spiec_easi_mb()`, `spiec_easi_glasso()` | [spiec_easi.md](spiec_easi.md) | clean room (SpiecEasi GPL) | High demand (1,335 citations 2022–26) but medium-to-high effort: graphical lasso and StARS selection. |
| 66 | [ ] `dmm()` | [community_typing.md](community_typing.md) | clean room (DirichletMultinomial LGPL; mothur GPL) | Medium effort (EM with model selection); steady use (506 citations 2022–26). |
| 67 | [ ] `splsda()` | [splsda.md](splsda.md) | clean room | mixOmics users; sparse PLS-DA is its own algorithm; medium effort. |
| 68 | [ ] `feast()` | [source_tracking.md](source_tracking.md) | **clean room** (FEAST is CC BY-NC-SA: non-commercial; even its outputs as fixtures need the user's OK) | FEAST is the other standard source-tracking method, rising (~140 citations a year); an EM fully stated in the paper, but its initialisation must be settled first. |
| 69 | [ ] `mantel_regression()` (MRM) | [mantel_regression.md](mantel_regression.md) | clean room | Simple algebra; biogeography. |
| 70 | [ ] `diffusion_map()` | [diffusion_map.md](diffusion_map.md) | port (pydiffmap MIT / scanpy BSD-3) | Growing use; uses neighbour and eigensolver code ecodive already has. |
| 71 | [ ] `pacmap()` | [pacmap.md](pacmap.md) | port (PaCMAP, Apache-2.0) | Growing use; keep the Apache notices. |
| 72 | [ ] `mpd()`, `mntd()`, `mpd_ses()`, `mntd_ses()` (NRI / NTI) | [phylo_structure.md](phylo_structure.md) | clean room (picante GPL) + port (phylocom, BSD-3) | The standard phylogenetic-structure summary (picante 3,272 citations); 999-draw null models are slow in R; tree-based threaded C avoids the cophenetic matrix; medium-to-high effort. Item 73 shares its tree code. |
| 73 | [ ] `beta_mntd()`, `beta_mpd()`, `beta_nti()`, `raup_crick_bray()`, `assembly_processes()` (Stegen framework; iCAMP later) | [assembly_processes.md](assembly_processes.md) | **clean room** (Stegen's scripts unlicensed; iCAMP, NST, picante GPL) | Surging (Stegen 2013: 2,024 citations, 317 → 479 a year; iCAMP 136 → 420); 999 randomizations per sample pair is where threaded deterministic C wins; high effort. Needs item 72. |
| 74 | [ ] `coda_lasso()`, then `coda_lasso_pairwise()` | [coda_lasso.md](coda_lasso.md) | port (zeroSum MIT; coda4microbiome MIT; glmnet output-only) | Compositional ML: log-contrast lasso / elastic net for prediction and taxon selection (coda4microbiome, selbal); moderate effort; CV folds thread naturally. Item 81's ConQuR reuses its lasso. |
| 75 | [ ] `secom_pearson()`, `secom_distance()` | [secom.md](secom.md) | port (ANCOMBC, Artistic-2.0) | Bias-corrected correlations (its own closed-form sample-bias estimate, not ANCOM-BC's EM); rising (195 citations); medium effort; three GPL helpers it calls must be written clean-room. |
| 76 | [ ] `mirkat()` | [mirkat.md](mirkat.md) | clean room | In demand in microbiome methods work; needs Davies or Liu p-values written from the papers. |
| 77 | [ ] `coinertia()` | [coinertia.md](coinertia.md) | clean room (RV from hoggorm, BSD-2) | Multi-omics; medium effort. |
| 78 | [ ] `dbmem()` | [dbmem.md](dbmem.md) | clean room | Spatial ecology. Needs item 59. |
| 79 | [ ] `kulczynski()`, `mountford()`, `raup_crick()`, `binomial_deviance()`, `chao_jaccard()`, `cao()`, `mahalanobis_distance()`, `alt_gower()`, `chisq_distance()` | [vegdist_metrics.md](vegdist_metrics.md) | formula; clean room (vegan GPL) | Completes `vegdist()` parity so vegan users find every metric, threaded; each is a small C kernel; low demand ("Kulczynski": 32 mentions). |
| 80 | [ ] `fastancom()` | [fastancom.md](fastancom.md) | clean room (GPL) | Low effort; the only microbiome-specific method Wirbel (2024) found to control FDR, but little use (22 citations). |
| 81 | [ ] `batch_combat_seq()`, `batch_conqur()`, then the `batch_correct()` wrapper (phase 2 of batch correction) | [batch_correction.md](batch_correction.md) | port (sva, Artistic-2.0) + clean room (edgeR parts GPL; ConQuR GPL-3) | ComBat-seq (1,786 citations, all fields) and ConQuR (167, rising) need edgeR's dispersions rebuilt and a quantile-regression solver plus lasso (item 74); high effort. |
| 82 | [ ] `rpca_joint()` | [rpca_joint.md](rpca_joint.md) | port (gemelli, BSD-3) | Multi-omics; a large port (~1,000 C lines) on `rpca()`'s OptSpace code. |
| 83 | [ ] `tree_plot()` | [tree_plot.md](tree_plot.md) | first principles (ape GPL) | Medium-to-high effort (a base-R tree drawer); shows `rpca_phylo()` clade loadings and `philr()` balances on the tree. Item 85 reuses its tree walk. |
| 84 | [ ] `halla()` (hierarchical all-against-all association testing) | [halla.md](halla.md) | port (HAllA, MIT; SciPy linkage, BSD-3) | Small but rising use (88 citations, 5 → 28 a year); turns item 39's pairwise table into blocks of associated features; medium effort. Needs item 39. |
| 85 | [ ] `coda_dendrogram()` | [coda_dendrogram.md](coda_dendrogram.md) | clean room (compositions GPL) | Medium effort; well cited in CoDa but rarely used in microbiome papers. |
| 86 | [ ] `ctf()` | [ctf.md](ctf.md) | port (gemelli, BSD-3) | Longitudinal designs; the largest gemelli port, with reference quirks to settle with the user first. |
| 87 | [ ] `phate()` | [phate.md](phate.md) | clean room; the reference's license is non-commercial | Single-cell crossover; legal care needed even for fixtures. |
| 88 | [ ] `rademu()` | [rademu.md](rademu.md) | port (radEmu, MIT) | High effort and little uptake so far (≤15 citations); methodologically strong, so kept for later. |
| 89 | [ ] `dca()` | [dca.md](dca.md) | clean room; **superseded** (criticised; NMDS preferred) | Kept for users who expect it, written to match vegan. |
| 90 | [ ] `largevis()` | [largevis.md](largevis.md) | port (uwot's BSD-2 headers); **superseded** by UMAP | Cheap, since `umap()` already uses these headers. |


## Ranking inputs

Scores for the open items (1 = lowest, 5 = highest), as set on 2026-11-20.
Rank is by demand − effort, then demand, then effort; dependencies (last
column) override the score. Re-rank with this table when demand or a plan
changes.

| # | Item | Demand | Effort | D − E | After |
|--:|:-----|-------:|-------:|------:|:------|
| 37 | alpha metrics | 4 | 2 | 2 | (moved up: 38, 48 and 64 use it) |
| 38 | `rarefaction_curve()` | 5 | 2 | 3 | 37 |
| 39 | `cross_assoc()` | 5 | 2 | 3 | |
| 40 | `decontam()` | 4 | 1 | 3 | |
| 41 | `lefse()` | 5 | 3 | 2 | |
| 42 | `dirmult_lm()` | 4 | 2 | 2 | |
| 43 | normalization factors | 3 | 1 | 2 | |
| 44 | biplot scalings | 3 | 1 | 2 | |
| 45 | `valencia()` | 3 | 1 | 2 | |
| 46 | log-ratio plots | 3 | 1 | 2 | 44 |
| 47 | `composition_barplot()` | 3 | 1 | 2 | |
| 48 | `alpha_test()` | 5 | 4 | 1 | 37 |
| 49 | `maaslin3()` | 4 | 3 | 1 | |
| 50 | `indicator_species()` | 4 | 3 | 1 | |
| 51 | `neutral_model()` | 3 | 2 | 1 | |
| 52 | `beta_partition()` | 3 | 2 | 1 | |
| 53 | batch correction, phase 1 | 3 | 2 | 1 | |
| 54 | `ternary_plot()` | 3 | 2 | 1 | |
| 55 | Bayesian zero replacements | 3 | 2 | 1 | |
| 56 | longitudinal, phases 1–2 | 2 | 1 | 1 | 39 |
| 57 | `principal_response_curves()` | 2 | 1 | 1 | |
| 58 | `mantel_correlogram()` | 2 | 1 | 1 | |
| 59 | `moran_i()` | 2 | 1 | 1 | |
| 60 | `variation_array()`, `zero_pattern()` | 2 | 1 | 1 | |
| 61 | `deseq2()` | 5 | 5 | 0 | 43 (shares RLE) |
| 62 | `random_forest()` | 5 | 5 | 0 | |
| 63 | `feature_volatility()` | 2 | 1 | 1 | 62 |
| 64 | `hill_extrapolation()` | 4 | 4 | 0 | 38 (coverage helpers) |
| 65 | `spiec_easi_*()` | 4 | 4 | 0 | |
| 66 | `dmm()` | 3 | 3 | 0 | |
| 67 | `splsda()` | 3 | 3 | 0 | |
| 68 | `feast()` | 3 | 3 | 0 | (licence: ask the user) |
| 69 | `mantel_regression()` | 2 | 2 | 0 | |
| 70 | `diffusion_map()` | 2 | 2 | 0 | |
| 71 | `pacmap()` | 2 | 2 | 0 | |
| 72 | `mpd()` ... `mntd_ses()` | 3 | 4 | −1 | |
| 73 | assembly processes | 4 | 5 | −1 | 72 |
| 74 | `coda_lasso()` | 2 | 3 | −1 | |
| 75 | `secom_*()` | 2 | 3 | −1 | |
| 76 | `mirkat()` | 2 | 3 | −1 | |
| 77 | `coinertia()` | 2 | 3 | −1 | |
| 78 | `dbmem()` | 2 | 3 | −1 | 59 |
| 79 | vegdist metrics | 1 | 2 | −1 | |
| 80 | `fastancom()` | 1 | 2 | −1 | |
| 81 | batch correction, phase 2 | 3 | 5 | −2 | 53, 74 |
| 82 | `rpca_joint()` | 2 | 4 | −2 | |
| 83 | `tree_plot()` | 2 | 4 | −2 | |
| 84 | `halla()` | 1 | 3 | −2 | 39 |
| 85 | `coda_dendrogram()` | 1 | 3 | −2 | 60, 83 |
| 86 | `ctf()` | 2 | 5 | −3 | |
| 87 | `phate()` | 1 | 4 | −3 | (licence) |
| 88 | `rademu()` | 1 | 5 | −4 | |
| 89 | `dca()` | — | — | — | superseded: last |
| 90 | `largevis()` | — | — | — | superseded: last |

Demand bands, roughly: 5 = over ~7,000 microbiome mentions (or the
field's default figure or test); 4 = over ~2,000 recent citations or
~1,000 mentions; 3 = ~300–1,000 mentions or ~1,000–2,000 citations;
2 = a few hundred citations; 1 = under ~100 citations or mentions.
Effort bands follow each plan's own estimate: 1 low, 2 low-to-medium,
3 medium, 4 medium-to-high, 5 high.


## History of the order

**2026-11-18 rules (the user's).** These still decide the effort scores:

1. Favour methods that are **simple**, that have an **obvious permissively
   licensed reference** to port, and that **microbiome researchers ask for
   now**.
2. **Clean-room** methods are harder. "Clean room" means there is no
   permissive code to port, so a substantial algorithm must be written
   from the paper and checked against GPL outputs.
3. **Superseded** methods go last.

Methods marked *derived* have no permissive reference code, but are thin
layers over ecodive's own functions plus a published formula (varpart is
adjusted R² from partial RDAs, for example). They aren't clean-room in a
meaningful sense.

**Compositional methods (added 2026-11-19).** The plans in
[compositional.md](compositional.md) were inserted by the same rules.
`rho_p()`, `phi_s()` and `qmp_scale()` are formulas on top of ported code,
so they count as derived, like varpart. `philr()` and the Bayesian zero
replacements have no permissive code, so they count as clean room.

**Source tracking and compositional regression (added 2026-11-19).**
The user asked whether machine learning belongs in ecodive. The answer
was no for neural networks and general-purpose learners: microbiome tables
are too small for deep learning, and ranger and XGBoost are already fast,
multithreaded and reproducible. The answer was yes for costly methods
specific to microbiome data, which don't yet have a fast, reproducible
implementation: `sourcetracker()` (now 36), `coda_lasso()` (now 74) and
`feast()` (now 68). SparCC was unplanned then; it was added later (item
32, done).

**Second compositional batch (2026-11-20).** The plans of the 2026-11-19
literature review were added and the open items re-ranked by **ease of
implementation and demand in the field** (citations and full-text
mentions in 2022–2026), with clean-room work counting as harder than a
port.

**One-stop-shop batch (2026-11-20).** The user wants ecodive to be the
only package a researcher needs for a modern, comprehensive microbiome
analysis. A gap review against vegan, scikit-bio / QIIME 2, picante,
iNEXT and recent papers found the gaps at the two ends of the workflow:
descriptive and QC steps (alpha metrics and tests, rarefaction curves,
coverage-based extrapolation, decontamination) and ecological
interpretation (indicator species, community assembly, beta partitioning,
the neutral model), plus LEfSe, DESeq2, two-table associations, batch
correction, q2-longitudinal and random forests. The user approved all of
them (items 37–41, 48, 50–53, 56, 61–64, 72, 73, 79, 81, 84) and
`maaslin3()`'s widened scope (item 49). This reverses the 2026-11-19 "no"
for one general-purpose learner: **random forests are in** (item 62),
because they are the machine learning microbiome papers actually use;
neural networks and boosting stay out. Then the whole list was renumbered:
completed items first, open items by the scores above.

Dropped: `ordisurf()` (it would need a GAM written in C; the user dropped
it on 2026-11-18).

**Still out of scope** (reasons in the gap review): sequence processing
(DADA2), taxonomy assignment, functional prediction (PICRUSt2 and
Tax4Fun need reference databases), shotgun profiling, compositional
mediation, neural networks and boosting.


## Renumbering map (2026-11-20)

Old number → new number. Items 1–24 kept their numbers. New plans have
no old number.

| Old | New | | Old | New | | Old | New | | Old | New |
|--:|--:|:-|--:|--:|:-|--:|--:|:-|--:|--:|
| 25 | 57 | | 37 | 30 | | 49 | 54 | | 61 | 47 |
| 26 | 58 | | 38 | 31 | | 50 | 42 | | 62 | 85 |
| 27 | 59 | | 39 | 32 | | 51 | 43 | | 63 | 83 |
| 28 | 70 | | 40 | 33 | | 52 | 45 | | 64 | 67 |
| 29 | 71 | | 41 | 34 | | 53 | 60 | | 65 | 77 |
| 30 | 82 | | 42 | 49 | | 54 | 80 | | 66 | 78 |
| 31 | 86 | | 43 | 44 | | 55 | 66 | | 67 | 88 |
| 32 | 25 | | 44 | 46 | | 56 | 65 | | 68 | 87 |
| 33 | 26 | | 45 | 35 | | 57 | 75 | | 69 | 89 |
| 34 | 27 | | 46 | 55 | | 58 | 68 | | 70 | 90 |
| 35 | 28 | | 47 | 36 | | 59 | 69 | | | |
| 36 | 29 | | 48 | 74 | | 60 | 76 | | | |


## Existing functions checked for splitting

2026-11-19: checked again for the compositional plans (`normalize(norm)`,
`pseudocount`, `rarefy()`). Nothing needs splitting; see
[compositional.md](compositional.md).

### 2026-11-18

The rule is one function per named method variant (see README).

**Split:**
- `mantel(given =)` becomes `mantel_partial()` (item 1).

**Kept as arguments,** because every reference treats them as options of
one method, not as separate methods:

| Argument | Why it stays an argument |
|:---------|:-------------------------|
| `pcoa(correction)` | ape, vegan and `cmdscale()` all treat Lingoes and Cailliez as options of PCoA. |
| `permdisp(type)` | vegan and Anderson's PERMDISP2 both treat centroid and median as options. |
| `permdisp(bias_adjust)` | An option in vegan too. |
| `mantel(method)` | The choice of correlation coefficient. |
| `mantel(alternative)` | A test option. |
| `nmds(ties)` | An option of SMACOF in scikit-learn and vegan. |
| `tsne(gradient)` | A computation route (exact or Barnes-Hut), not a different method. |
| `tsne(init)`, `umap(init)` | Starting layouts. |
| `pca(scale)` | `prcomp()`'s `scale.`. |
| `pca(norm)` | The transformation applied first, as in every function with `norm`. |
| `permanova(by)` | Which ANOVA table to compute. |
| `mrpp(weights)` | Group weights. |
| `rpca(distances)` | An extra output. |
| `envfit(data / counts)` | What is fitted, not how. |

**Already split:** `mds()` / `nmds()`, `dbrda()` / `capscale()`, and the
UniFrac family.
