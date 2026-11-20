# To-do: new methods, by priority

Work top to bottom. Each item links to its plan; read [README.md](README.md)
first. Tick an item, and add the commit, when its method is merged with tests
and docs.

The order follows three rules the user set (2026-11-18):

1. Favour methods that are **simple**, that have an **obvious permissively
   licensed reference** to port, and that **microbiome researchers ask for
   now**.
2. **Clean-room** methods go to the bottom, with DPCoA first among them.
   "Clean room" means there is no permissive code to port, so a substantial
   algorithm must be written from the paper and checked against GPL outputs.
3. **Superseded** methods go last.

Methods marked *derived* have no permissive reference code, but are thin
layers over ecodive's own functions plus a published formula (varpart is
adjusted R² from partial RDAs, for example). They aren't clean-room in a
meaningful sense, so they sit just after the methods they build on.

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
| 15 | [x] `replace_zeros()` wrapper (+ `list_zero_replacements()`, `match_zero_replacement()`) and `replace_zeros_mult()` (7d84f63) | [zero_replacement.md](zero_replacement.md) | port (scikit-bio, BSD-3) | Principled zeros for every log-ratio method; the wrapper is ready for the Bayesian methods (later item). |
| 16 | [x] `log_ratio()`, `balance()` (c56ea2a) | [log_ratio.md](log_ratio.md) | port (Qurro / scikit-bio, BSD-3) | Trivial; log ratios are what current compositional guidance recommends analysing. |
| 17 | [x] `vlr()`, `rho_p()`, `phi_s()` (60d304e) | [proportionality.md](proportionality.md) | port (scikit-bio `vlr`) + formula (propr output-only) | High demand (propr); all-pairs C over features is ecodive's strength. |
| 18 | [x] `qmp_scale()`, `qmp_rarefy()`; per-sample `rarefy(depth =)` (068453a) | [qmp.md](qmp.md) | port (raeslab QMP.R, BSD-3) + formula | Simple, over `rarefy()`'s C; the absolute-abundance side of the debate. |
| 19 | [x] `ancom()` (63822a8) | [ancom.md](ancom.md) | port (scikit-bio, BSD-3) | In QIIME 2 and widely used; D² independent tests suit threaded C. |
| 20 | [x] `dirmult_ttest()` (2291a31), then its scale model (c88b95b) | [dirmult_ttest.md](dirmult_ttest.md) | port (scikit-bio, BSD-3; scale model from ALDEx3, MIT) | ALDEx2-style effect sizes, high demand; needs a PCG Gamma sampler. |
| 21 | [x] `rpca_phylo()` (4b259ba) | [rpca_phylo.md](rpca_phylo.md) | port (gemelli, BSD-3) on `rpca()` | QIIME 2's phylogenetic RPCA; reuses `rpca()` and the UniFrac tree walk. Phylogeny-aware. |
| 22 | [x] `ancombc()`, `ancombc2()` (49a6e59) | [ancombc.md](ancombc.md) | port (scikit-bio, BSD-3, derived from Artistic-2.0 ANCOMBC: keep both notices) | Among the most used DA methods now; large (~1,000 lines of core). |
| 23 | [x] `ancombc_global()`, `ancombc_pairwise()`, `ancombc_dunnett()`, `ancombc_trend()`, `structural_zeros()` (7d9b859) | [ancombc.md](ancombc.md) | port (scikit-bio + ANCOMBC notices), as above | Post-hoc tests on a fitted `ancombc()` / `ancombc2()`; `ancombc_trend()` needs a small QP solver. |
| 24 | [x] `capda()` (06f2fbf) | [lda_capda.md](lda_capda.md) | derived (`pcoa()` + `lda()`) | Anderson's CAP, popular in marine ecology. |
| 25 | [ ] `principal_response_curves()` | [principal_response_curves.md](principal_response_curves.md) | derived (partial `rda()`) | Experiments over time. |
| 26 | [ ] `mantel_correlogram()` | [mantel_correlogram.md](mantel_correlogram.md) | derived (`mantel()` on distance classes) | Distance decay, biogeography. |
| 27 | [ ] `moran_i()` | [moran_i.md](moran_i.md) | port (PySAL esda, BSD-3) | Simple; `dbmem()` needs it. |
| 28 | [ ] `diffusion_map()` | [diffusion_map.md](diffusion_map.md) | port (pydiffmap MIT / scanpy BSD-3) | Growing use; uses neighbour and eigensolver code ecodive already has. |
| 29 | [ ] `pacmap()` | [pacmap.md](pacmap.md) | port (PaCMAP, Apache-2.0) | Growing use; keep the Apache notices. |
| 30 | [ ] `rpca_joint()` | [rpca_joint.md](rpca_joint.md) | port (gemelli, BSD-3) | Multi-omics; a large port (~1,000 C lines) on `rpca()`'s OptSpace code. |
| 31 | [ ] `ctf()` | [ctf.md](ctf.md) | port (gemelli, BSD-3) | Longitudinal designs; the largest gemelli port, with reference quirks to settle with the user first. |
| 32 | [x] `philr()`, `philr_distance()` (5f8d7bb) | [philr.md](philr.md) | **clean room** (philr GPL-3; formulas from the open-access paper) | First among the clean-room methods: phylogeny-aware, high demand, and an O(D) formula per sample. |
| 33 | [x] `dpcoa()` (6ed9799) | [dpcoa.md](dpcoa.md) | **clean room** | Second among the clean-room methods, after `philr()`: phyloseq "DPCoA", phylogeny-aware. |
| 34 | [x] Documentation: Schloss (2026) and depth in the RPCA family; `scale_sd` sensitivity; ML caveats; predicted loads in `qmp_scale()`; `log_base` on DA results (60b71df) | [compositional_docs.md](compositional_docs.md) | docs | Almost no code (one `log_base` attribute), and it corrects how existing functions are presented: Schloss's 2026 rebuttal of robust Aitchison PCA, scale-assumption FDR results. |
| 35 | [x] Effect (MW), MA and forest plots for DA results (`plot(type = )`), with new dispersion and mean-abundance columns in `dirmult_ttest()` (48cf7ad) | [da_plots.md](da_plots.md) | port (ALDEx3's effect statistics, MIT; ALDEx2 output-only) | Low effort on results users already have; Gloor et al. (2016) show volcano plots, ecodive's only DA display, mislead. |
| 36 | [x] `da_test()` wrapper (+ `list_da_tests()`, `match_da_test()`); `da_wilcoxon()`, `da_ttest()`, `da_lm()` (4d0226c) | [da_test.md](da_test.md) | formula | Very low effort; the plain tests are what Wirbel (2024) and Pelto (2025) found best at controlling false discoveries; four DA methods now justify the wrapper. |
| 37 | [x] `beta_div_rarefied()`, `alpha_div_rarefied()` (963647f) | [rarefied_diversity.md](rarefied_diversity.md) | derived (on `rarefy()`) | Low effort; rarefaction is the most mentioned normalization (~3,600 papers a year) and Schloss (2024, 2026) argues for it; lets users compare with RPCA. |
| 38 | [x] `linda()` (f09370e) | [linda.md](linda.md) | clean room (LinDA GPL) | Short estimator sharing `ancombc()`'s code; 494 citations and fast-rising; most robust under strong compositional effects. |
| 39 | [x] `sparcc()` (21c3da7) | [sparcc.md](sparcc.md) | clean room (CC-BY paper) | The most-used compositional network method (1,848 citations 2022–26); builds on `vlr()` and the existing Gamma sampler. |
| 40 | [x] `assoc_network()`, `assoc_network_compare()`, `assoc_network_diff()` (a4d81fd; with `feature_assoc()`) | [assoc_network.md](assoc_network.md) | clean room (NetCoMi GPL) | Textbook formulas over `rho_p()` / `sparcc()`; NetCoMi uptake grew 3.1× from 2021–22 to 2024–25. |
| 41 | [x] `da_sensitivity()` (8d25396) | [da_sensitivity.md](da_sensitivity.md) | derived (on `da_test()`) | Reruns and plots DA calls over pseudocounts, zero treatments and `scale_sd`; the benchmarks' central concern. Needs item 36. |
| 42 | [ ] `maaslin3()` — **scope decision: ask the user first** | [maaslin3.md](maaslin3.md) | port (maaslin3, MIT) | The highest-uptake DA family (MaAsLin2 2,850 citations; MaAsLin3 rising steeply); moderate effort; widens ecodive into general per-feature regression. |
| 43 | [ ] Biplot scalings (`plot(ord, scaling = )`) with links | [biplot_scaling.md](biplot_scaling.md) | formula | Low effort; the compositional biplot Gloor (2017) recommends; Aitchison & Greenacre 2002 has 830 citations. |
| 44 | [ ] Rank plots and `plot()` of `log_ratio()` / `balance()` results (Qurro-style) | [log_ratio_plots.md](log_ratio_plots.md) | port (Qurro, BSD) | Low effort on ported code; turns DA effects and loadings into interpretable log ratios (reference frames: 668 citations). |
| 45 | [x] `enterotype()` (20f28e1) | [community_typing.md](community_typing.md) | clean room (PAM; `cluster` is GPL) | Low effort on `jsd()`; enterotyping is among the most cited community analyses (Arumugam 2011: 3,182 citations 2022–26). |
| 46 | [ ] `replace_zeros_czm()`, `replace_zeros_gbm()`, `replace_zeros_sq()`, `replace_zeros_bl()` | [zero_replacement.md](zero_replacement.md) | **clean room** (zCompositions GPL; the defining paper is closed access, so formulas are confirmed from outputs) | Simple once the formulas are pinned down; `cmultRepl()`'s GBM is the field's default. |
| 47 | [x] `sourcetracker()` (12034fd) | [source_tracking.md](source_tracking.md) | port (SourceTracker2, BSD-3) | Source tracking is in steady high demand (Knights 2011: ~200 citations a year); simple; SourceTracker2's per-read Python loop is slow and its `multiprocessing` isn't reproducible, so threaded C with per-sink PCG streams is a clear win. |
| 48 | [ ] `coda_lasso()`, then `coda_lasso_pairwise()` | [coda_lasso.md](coda_lasso.md) | port (zeroSum MIT; coda4microbiome MIT; glmnet output-only) | Compositional ML: log-contrast lasso / elastic net for prediction and taxon selection, in active use (coda4microbiome, selbal); builds on the log-ratio conventions; CV folds thread naturally. |
| 49 | [ ] `ternary_plot()` | [ternary_plot.md](ternary_plot.md) | clean room (ggtern, Ternary GPL) | The most used classic CoDa display (331 microbiome papers 2022–26); low-to-medium effort. |
| 50 | [ ] `dirmult_lm()` | [dirmult_lm.md](dirmult_lm.md) | port (ALDEx3, MIT) | Extends `dirmult_ttest()` beyond two groups with covariates and load-informed scale priors; low-to-moderate effort; ALDEx2 is widely used, ALDEx3 is new. |
| 51 | [ ] Normalization `norm` choices: CSS, TMM, RLE, GMPR | [normalization_factors.md](normalization_factors.md) | port (CSS, Artistic-2.0) + clean room | Low effort each; CSS and TMM are mentioned ~120–145 times a year; they can later centre scale models. |
| 52 | [ ] `valencia()` | [community_typing.md](community_typing.md) | port (VALENCIA, MIT) | Low effort; the standard vaginal community-state-type classifier, rising (499 citations, 121 in 2026). |
| 53 | [ ] `variation_array()`, `zero_pattern()` | [coda_diagnostic_plots.md](coda_diagnostic_plots.md) | clean room | Low effort diagnostics (Aitchison's variation array; zero patterns as in zCompositions); modest uptake. |
| 54 | [ ] `fastancom()` | [fastancom.md](fastancom.md) | clean room (GPL) | Low effort; the only microbiome-specific method Wirbel (2024) found to control FDR, but little use (22 citations). |
| 55 | [ ] `dmm()` | [community_typing.md](community_typing.md) | clean room (DirichletMultinomial LGPL; mothur GPL) | Medium effort (EM with model selection); steady use (506 citations 2022–26). |
| 56 | [ ] `spiec_easi_mb()`, `spiec_easi_glasso()` | [spiec_easi.md](spiec_easi.md) | clean room (SpiecEasi GPL) | High demand (1,335 citations 2022–26) but medium-to-high effort: graphical lasso and StARS selection. |
| 57 | [ ] `secom_pearson()`, `secom_distance()` | [secom.md](secom.md) | port (ANCOMBC, Artistic-2.0) | Bias-corrected correlations (its own closed-form sample-bias estimate, not ANCOM-BC's EM); rising (195 citations); medium effort; three GPL helpers it calls must be written clean-room. |
| 58 | [ ] `feast()` | [source_tracking.md](source_tracking.md) | **clean room** (FEAST is CC BY-NC-SA: non-commercial; even its outputs as fixtures need the user's OK) | FEAST is the other standard source-tracking method, rising (~140 citations a year); an EM fully stated in the paper, but its initialisation must be settled first. |
| 59 | [ ] `mantel_regression()` (MRM) | [mantel_regression.md](mantel_regression.md) | clean room | Simple algebra; biogeography. |
| 60 | [ ] `mirkat()` | [mirkat.md](mirkat.md) | clean room | In demand in microbiome methods work; needs Davies or Liu p-values written from the papers. |
| 61 | [ ] `composition_barplot()` | [composition_barplot.md](composition_barplot.md) | port (microshades idea, MIT) | Low effort, but stacked bars are not a compositional analysis; a better-ordered version of a ubiquitous figure. |
| 62 | [ ] `coda_dendrogram()` | [coda_dendrogram.md](coda_dendrogram.md) | clean room (compositions GPL) | Medium effort; well cited in CoDa but rarely used in microbiome papers. |
| 63 | [ ] `tree_plot()` | [tree_plot.md](tree_plot.md) | first principles (ape GPL) | Medium-to-high effort (a base-R tree drawer); shows `rpca_phylo()` clade loadings and `philr()` balances on the tree. |
| 64 | [ ] `splsda()` | [splsda.md](splsda.md) | clean room | mixOmics users; sparse PLS-DA is its own algorithm. |
| 65 | [ ] `coinertia()` | [coinertia.md](coinertia.md) | clean room (RV from hoggorm, BSD-2) | Multi-omics. |
| 66 | [ ] `dbmem()` | [dbmem.md](dbmem.md) | clean room | Spatial ecology. |
| 67 | [ ] `rademu()` | [rademu.md](rademu.md) | port (radEmu, MIT) | High effort and little uptake so far (≤15 citations); methodologically strong, so kept for later. |
| 68 | [ ] `phate()` | [phate.md](phate.md) | clean room; the reference's license is non-commercial | Single-cell crossover; legal care needed even for fixtures. |
| 69 | [ ] `dca()` | [dca.md](dca.md) | clean room; **superseded** (criticised; NMDS preferred) | Kept for users who expect it, written to match vegan. |
| 70 | [ ] `largevis()` | [largevis.md](largevis.md) | port (uwot's BSD-2 headers); **superseded** by UMAP | Cheap, since `umap()` already uses these headers. |

**Second compositional batch (2026-11-20).** Items 34 onward were
re-ranked when the plans of the 2026-11-19 literature review were added
(see [compositional.md](compositional.md)). The user set the rule: rank
by **ease of implementation and demand in the field** (citations and
full-text mentions in 2022–2026). Clean-room work counts as harder than a
port, and superseded methods have little demand, so they stay last. Items
before 34 keep their numbers.

**Compositional methods (added 2026-11-19).** The plans in
[compositional.md](compositional.md) were inserted by the same rules. The
permissive ports come first because they are simple and in high demand in
microbiome work now. `rho_p()`, `phi_s()` and `qmp_scale()` are formulas
on top of ported code, so they count as derived, like varpart. `philr()`
and the Bayesian zero replacements have no permissive code, so they go
with the clean-room methods. `philr()` leads that group as the simpler
of its two phylogeny-aware methods.

**Source tracking and compositional regression (added 2026-11-19).**
The user asked whether machine learning belongs in ecodive. The answer
was no for neural networks and general-purpose learners: microbiome tables
are too small for deep learning, and ranger and XGBoost are already fast,
multithreaded and reproducible. The answer was yes for costly methods
specific to microbiome data, which don't yet have a fast, reproducible
implementation. Two of those weren't planned yet. The user had them
placed after item 34, so that no item that is done or in progress is
renumbered (`dpcoa()`, item 33, was being implemented when they were
added):
- `sourcetracker()` (35) and `coda_lasso()` (36) are permissive ports in
  high demand;
- `feast()` (37) is clean room.

Items from `mantel_regression()` on moved down by three. SparCC stays
unplanned, for the reasons in [compositional.md](compositional.md).

Dropped: `ordisurf()` (it would need a GAM written in C; the user dropped
it on 2026-11-18).


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
