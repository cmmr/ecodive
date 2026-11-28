# Parameter names: glossary and proposed renames (2026-11-18)

**Status:** all renames in §2 and §3 were approved by the user and applied on 2026-11-18.

I compared every exported signature (dumped from the installed package)
and every signature proposed in the method plans, looking for two things:
the same name meaning different things, and different names meaning the
same thing.

- §1 is the glossary that existing and new functions should follow.
- §2 lists proposed changes to **existing** functions. They need the
  user's OK. All of them are in functions new in the unreleased 2.4.0.
- §3 lists changes to the **planned** signatures.
- §4 lists inconsistencies left as they are on purpose.


## 1. Glossary: one name, one meaning

| Name | Meaning | Used by |
|:-----|:--------|:--------|
| `counts` | a count table (first argument) | metrics, `pca`, `rpca`, `normalize`, `rarefy`, `count_ord`, and all planned count-table methods |
| `dm` | beta diversity distances (first argument) | `pcoa`, ..., all tests, `dist_ord`, `dist_test` |
| `ord` | an ordination result (first argument) | `envfit`, `pcoa_biplot`; planned `procrustes` |
| `x`, `y` | two inputs of equal standing, or a first argument that can be either `dm` or `counts` | planned `coinertia(x, y)`, `varpart(x, ...)`, `forward_select(x, ...)`; `enterotype(x, ...)` (with `counts` as a later, optional argument when `x` is a dist); `mantel(dm, y)` |
| `n_axes` | **the number of axes returned** (proposed rename of `k`; see §2.0) | every ordination |
| `axes` | **which** axes of an input ordination to use (indices) | `envfit`; planned `procrustes`; `log_ratio_features` (one axis, 2026-11-24, item 46) |
| `groups`, `data`, `formula` | the grouping variable or model, and the sample metadata it comes from. DA model methods (2026-11-20, `da_test()`): `formula` is the model and `groups` its term of interest; `groups` alone is `formula = ~ groups` | tests, `dbrda`, `capscale`, DA methods, `neutral_model` (one fit per group, 2026-11-24), and planned constrained and discriminant methods |
| `given` | distances (or variables) whose effect is partialled out | `mantel` (moving to `mantel_partial`); planned `bioenv`; `cross_assoc` (covariates: a formula over `data`, column names, or a table) |
| `margin`, `norm`, `pseudocount` | the count table's orientation and transformation (`norm = 'log'`, log-TSS, added 2026-11-20; `norm = 'chisq'`, the chi-square transformation, added 2026-11-27, item 79) | everything that takes `counts` |
| `scale` | standardise features to unit variance | `pca`; planned `rda`, `plsda`, `oplsda`, `splsda` |
| `scaling` | how a biplot scales samples against features: `'form'` (1), `'covariance'` (2), `'contribution'`; vegan's 1, 2, 3 for `rda`, `cca`, `ca` (item 44, 2026-11-24); vegan's 1, 2, 3 (`'sites'`, `'species'`, `'symmetric'`) for the coefficients and weights of `principal_response_curves` (item 57) | `plot.ecodive_ord`, `principal_response_curves` |
| `links` | pairs of features joined in a covariance biplot (Aitchison & Greenacre's term) | `plot.ecodive_ord` |
| `method` | a correlation coefficient (`'pearson'`, `'spearman'`, `'kendall'`), or the method name in `dist_ord` and `count_ord` | `mantel`; planned `mantel_partial`, `mantel_correlogram`, `bioenv`; `mantel_regression` takes only `'pearson'` (least squares) and `'spearman'` (least squares on ranks: ecodist's `mrank = TRUE`) (2026-11-26, item 69) |
| `metric` | a beta diversity metric name | `beta_div`, `alpha_div`; planned `bioenv` |
| `alternative` | `'greater'`, `'less'` or `'two.sided'`; `mantel_correlogram()` adds `'sign'` (its default, vegan's: each class tested in the direction of its statistic's sign; 2026-11-25, item 58) | `mantel`, `mantel_partial`, `mantel_correlogram`, `moran_i` (for its analytic p-values too) |
| `permutations`, `strata`, `seed` | the permutation test (a count or a matrix), restricted permutations, and the RNG seed | all tests, `envfit` |
| `pairwise`, `p_adjust` | per-pair tests, and their p-value adjustment | planned on all group tests; `p_adjust` also on `mantel_correlogram` and optionally `simper` |
| `breaks`, `n_classes`, `cutoff`, `progressive` | **(new, 2026-11-25, item 58)** distance classes' boundaries, or their number (equal widths; Sturges' rule by default); whether to test only the classes vegan's `cutoff` keeps; and a progressive multiple-testing adjustment (vegan's `progressive`), alongside `p_adjust` (which also takes vegan's `'sidak'` there) | `mantel_correlogram` |
| `n_<unit>` | **an exact count** | `n_neighbors`, `n_starts`, `n_epochs`; planned `n_iters` (pacmap), `n_segments` (dca), `n_classes`, `n_landmarks`, `n_orthogonal` |
| `max_<unit>` | **a cap** (it may stop earlier) | `max_iters` |
| `tol` | a convergence tolerance | `mds`, `nmds`, `rpca`, `coda_lasso` (zeroSum's `precision`), `coda_lasso_pairwise` (glmnet's `thresh`); planned `splsda` |
| `init`, `learning_rate` | starting layout, step size | `tsne`, `umap`; planned `pacmap`, `largevis` |
| `space` | **(new)** sample positions: coordinates, geographic distances, or spatial weights. In `moran_i()` (2026-11-25, item 59): a square matrix is spatial weights, a `dist` is distances and other tables (or column names of `data`) are coordinates; distances and coordinates become inverse-distance weights. In `dbmem()` (2026-11-26, item 78): coordinates or a `dist` only (a square matrix of weights is an error) | `moran_i`, `dbmem` |
| `autocor` | **(new, 2026-11-26, item 78)** which Moran's eigenvector maps to keep, by the sign of their Moran's I: `'positive'`, `'negative'`, `'non-null'`, `'all'` (adespatial's `MEM.autocor`; replaces the plan's `moran_positive`) | `dbmem` |
| `style` | **(new, 2026-11-25, item 59)** how spatial weights are transformed: `'row'` (row-standardised), `'binary'`, `'variance'` (variance-stabilising), `'none'` (spdep's `style`, esda's `transformation`; not `standardize` or `transform`, which mean other things in `linda()` and `maaslin3()`) | `moran_i` |
| `numerator`, `denominator` | **(new, 2026-11-19)** features on the top and bottom of a log ratio (names, indices or a logical vector); `denominator = NULL` means the clr's geometric mean | `alr` (`denominator` only), `log_ratio`, `balance`, `rho_p`, `phi_s` (`denominator` only: one feature for an alr, or a set); the rank plots (`type = 'rank'`, item 46) of the DA plot methods and `plot.ecodive_ord()`, the features to highlight; `log_ratio_features()` returns them as a list |
| `basis` | **(new)** an ILR basis or a sequential binary partition | `ilr` |
| `load` | **(new, 2026-11-19)** each sample's measured total microbial load: one finite value > 0 per sample (matched by name, else by position), or a column name of `data`, by the same `per_sample()` rule as `groups` | `qmp_scale`, `qmp_rarefy`, `dirmult_lm` (as its scale model) |
| `contrast` | **(new)** the samples compared, `c(treatment, reference)`; `reference = NA` pools all other samples | `dirmult_ttest`, `da_wilcoxon`, `da_ttest` |
| `conf_level` | **(new, 2026-11-19)** the confidence level of an interval, strictly between 0 and 1 (as `t.test(conf.level)`, snake_case) | `dirmult_ttest`, `da_ttest`, `da_lm`; `plot.ecodive_ancombc` (its forest plot's Wald intervals); `neutral_model` (the interval for m and the Wilson band); `ternary_plot` (its ellipses, item 54) |
| `pseudocounts`, `zeros`, `scale_sds`, `crossed` | **(new, 2026-11-20)** a sensitivity analysis's grids, plural so that they don't capture the method's own `pseudocount` and `scale_sd` (they follow `...`, so they match only exactly; no DA method may take these names): pseudocounts to try, zero replacements (`replace_zeros()` ids) to try, scale uncertainties to try (`TRUE` for each grid's default values), and `crossed` (`TRUE`: every combination, with pseudocounts and zero replacements as alternatives) | `da_sensitivity` |
| `zero_handling`, `corr_cut`, `outlier_pct`, `standardize` | **(new, 2026-11-20)** how zeros of counts are replaced before logs (`'adaptive'`, `'pseudocount'`, `'imputation'`; not `zeros`, which `da_sensitivity()` takes), the p-value cutoff of the adaptive choice, the proportion of each feature's high values winsorized, and whether numeric covariates are centred and scaled first (LinDA's `zero.handling`, `corr.cut`, `outlier.pct`, and its scaling of `meta.dat`) | `linda` |
| `frac`, `adjust` | **(new, 2026-11-27, item 55)** zCompositions' names: the fraction of a detection limit that a zero becomes, strictly between 0 and 1 (`czm`'s `threshold / n`, and with `adjust` a feature's smallest observed proportion); and whether imputed values above that smallest proportion are lowered to `frac` of it. `delta` (mult) is the replaced value itself, as a proportion | `replace_zeros_czm`, `replace_zeros_gbm`, `replace_zeros_sq`, `replace_zeros_bl` |
| `var_equal` | **(new, 2026-11-20)** `TRUE` for Student's t-test with a pooled variance, `FALSE` (the default) for Welch's (as `t.test(var.equal)`, snake_case) | `da_ttest` |
| `part_weights`, `ilr_weights` | **(new)** PhILR's taxon and balance weights | planned `philr` |
| `subjects`, `states` | **(new)** per-sample subject ids and time points / conditions (numeric time points for the longitudinal functions), by the same `per_sample()` rule as `groups` | `first_differences`, `first_distances`, `volatility`, `nmit` (2026-11-24, item 56), `feature_volatility` (item 63); planned `ctf` |
| `treatment`, `time`, `control`, `plots` | **(new, 2026-11-25, item 57)** each sample's treatment and time point, by `per_sample()`'s rule (time points are categories, so not `states`, which the longitudinal functions take as numbers); the control treatment (a level of `treatment`); and the experimental units whose series of samples are permuted whole (permute's `Plots(strata = )`) | `principal_response_curves` |
| `baseline`, `replicates` | **(new, 2026-11-24)** the state that every other is compared with (`NULL`: consecutive states), and what to do with a subject's several samples at one state (`'error'`, `'drop'`, `'random'`; q2-longitudinal's `replicate_handling`) | `first_differences`, `first_distances` |
| `matrix_norm` | **(new, 2026-11-24)** the norm of a matrix difference (`'frobenius'`, `'nuclear'`); not `norm`, which normalizes counts | `nmit` |
| `n_resamples` | **(new, 2026-11-24)** an exact count of bootstrap resamples (`n_<unit>` rule) | `volatility`, `feature_volatility` |
| `test_size` | **(new, 2026-11-26, item 63)** the fraction of samples held out of a model to measure its accuracy on (scikit-learn's and q2's `test_size`; `0`: none) | `feature_volatility` |
| `rfe`, `rfe_step`, `rfe_folds` | **(new, 2026-11-26, item 63)** recursive feature elimination with cross-validation (scikit-learn's RFECV; q2's `optimize_feature_selection`): whether to run it, how many features to drop each round (a fraction below 1 or a count; q2's `step`), and its number of folds (q2's `cv`; not `n_folds`, which is the forest's own cross-validation) | `feature_volatility` |
| `backward` | **(new, 2026-11-26, item 69 follow-up)** whether to simplify a regression model by backward elimination of terms (Lichstein 2007), with `alpha` as the level a term's p-value must reach to stay (as `forward_select()`'s `alpha` for adding one) | `mantel_regression` |
| `fit` | **(new)** a fitted result given to a post-hoc test | planned `ancombc_*` |
| `family`, `l1_ratio`, `lambda`, `n_lambdas`, `lambda_min_ratio` | **(new, 2026-11-19; exist 2026-11-27, item 74)** the outcome distribution (`'gaussian'`, `'binomial'`; `NULL` chooses from `y`), elastic-net mixing (1 = lasso), and the penalty path | `coda_lasso`, `coda_lasso_pairwise` |
| `n_folds`, `folds`, `select` | **(new)** cross-validation: a fold count or a per-sample fold vector (overrides `seed`), and which λ to report (`'min'`, `'1se'`) | `random_forest` (`n_folds`, `folds`; exists, 2026-11-25); `coda_lasso`, `coda_lasso_pairwise` (all three; 2026-11-27) |
| `y` | **(item 62, 2026-11-25)** what a supervised model predicts: one value per sample or a column name of `data` (a factor, character or logical `y` → classification; numeric → regression) | `random_forest`; `coda_lasso`, `coda_lasso_pairwise` (two classes → binomial, numeric → gaussian) |
| `standardize` (penalized fits) | **(item 74, 2026-11-27)** penalize each coefficient in proportion to its column's sd, as if the columns were scaled (glmnet's and zeroSum's `standardize`); coefficients stay on the original scale. (In `linda()`, `standardize` scales the covariates; both are "scale before fitting".) | `coda_lasso`, `coda_lasso_pairwise` |
| `n_trees`, `max_features`, `min_samples_leaf`, `min_samples_split`, `max_depth`, `criterion`, `class_weight`, `max_samples` | **(item 62, 2026-11-25)** a random forest's settings, scikit-learn's names (`n_trees` for its `n_estimators`, by the `n_<unit>` rule): trees per forest; features per split (a count, a fraction below 1, `'sqrt'`, `'log2'`, `'all'`; a vector is a grid); fewest samples per leaf (a vector is a grid) and per split; depth limit; impurity (`'gini'`, `'entropy'`, `'squared_error'`); class weights (`'balanced'`, `'balanced_subsample'`, or named weights); bootstrap size (a count or a fraction below 1) | `random_forest`; `feature_volatility` (`n_trees`, and the others through `...`) |
| `importance`, `keep_forest` | **(item 62, 2026-11-25)** which feature importances to compute (`'impurity'`, `'permutation'`, `'both'`), and whether to keep the fitted model in the result for `predict()` | `random_forest`; `feature_volatility` (`'impurity'` or `'permutation'`: what ranks the features) |
| `sources`, `sink_depth` | **(new)** per-sample source-environment labels (`NA` marks a sink), and the sinks' rarefaction depth (`depth` is the sources', as in `rarefy`) | `sourcetracker` (exists); `sources` also planned `feast` |
| `n_restarts`, `n_burnin`, `thin`, `loo`, `features` | **(new)** the Gibbs schedule (independent chains, passes before the first draw, passes between draws), leave-one-out over environments, and per-feature output. `n_burnin` and `thin` replace SourceTracker's `burnin` and `delay` (the user's choice, 2026-11-20): `n_burnin` is an exact count (`n_<unit>` rule) and `thin` is the usual MCMC term | `sourcetracker` (exists) |
| `n_types` | **(new, 2026-11-20)** the candidate numbers of community types (a single value fixes it) | `enterotype`, `dmm` |
| `ch` | **(new, 2026-11-20)** the form of the Calinski-Harabasz index: `'medoids'` (clusterSim's, as in Arumugam et al.'s script) or `'centroids'` (the distances' pseudo-F) | `enterotype` |
| `centroids` | **(new)** reference type profiles to classify against: a types × features matrix, types named in the row names | `valencia` |
| `pool` | **(new, 2026-11-24, item 51)** a count table defining the source metacommunity (Burns et al.'s `pool`): each feature's mean relative abundance is its share of the pool's counts; matched to `counts`' features by name | `neutral_model` |
| `totals` | **(new, 2026-11-24)** each sample's total read count, by which its counts are divided instead of its sum (VALENCIA's `read_count`); one per sample, by name or position | `valencia` |
| `eta`, `nu` | **(new, 2026-11-27, item 66)** the shape and rate of a Gamma prior on every Dirichlet parameter (DirichletMultinomial's values, 0.1 and 0.1, by default) | `dmm` |
| `criterion` (community types) | **(2026-11-27, item 66)** which model-selection criterion chooses the number of types, the smallest winning: `'laplace'`, `'aic'` or `'bic'` (random_forest's `criterion` is its impurity, above; the meanings never meet in one function) | `dmm` |
| `q` | **(new, 2026-11-24)** the order of a Hill number or of Rényi's or Tsallis's entropy: a single number ≥ 0, `Inf` allowed; default 2, as scikit-bio's `order` (vegan's `scales`) | `hill`, `renyi`, `tsallis`, `alpha_div`, `alpha_div_rarefied`, `rarefaction_curve` |
| `bias_corrected` | **(new, 2026-11-24)** Chao1's bias-corrected form (`TRUE`, the default, as scikit-bio) or the classic one (`FALSE`), as scikit-bio's `bias_corrected` | `chao1`, `alpha_div`, `alpha_div_rarefied` |
| `conc` | **(new, 2026-11-24)** each sample's DNA concentration before sequencing (any quantity proportional to the input DNA; only its inverse relation with a contaminant's frequency matters): finite and > 0, matched by name, else by position, or a column name of `data` (`per_sample()`'s rule). Kept as decontam's name rather than reusing `load` (a measured total microbial load used as a scale), the user's choice | `decontam`, `decontam_frequency_plot` |
| `neg` | **(new, 2026-11-24)** which samples are negative controls: logical, `TRUE` a control, `NA` left out of the tests; by `per_sample()`'s rule | `decontam`, `decontam_not`, `decontam_frequency_plot` |
| `batch`, `batch_combine` | **(new, 2026-11-24)** per-sample batches (sequencing run, extraction plate) scored separately, by `per_sample()`'s rule; and how a feature's per-batch scores are combined (`'minimum'`, `'product'`, `'fisher'`: decontam's `batch.combine`, snake case). In batch correction (item 53), the batches whose effects are removed: at least two, a factor's levels in order, else the sorted values | `decontam`, `batch_mmuphin`, `batch_combat`; planned `batch_combat_seq`, `batch_conqur`, `batch_correct` |
| `zero_inflation` | **(new, 2026-11-24, item 53)** leave zeros out of every fit, so a feature is adjusted only in the batches it is present in, and zeros stay zero (MMUPHin's `control$zero_inflation`) | `batch_mmuphin` |
| `parametric`, `mean_only`, `ref_batch` | **(new, 2026-11-24, item 53)** ComBat's parametric (`TRUE`) or non-parametric empirical Bayes priors (sva's `par.prior`); adjust only the batches' means (`mean.only`); the batch left unchanged that the others are adjusted to (`ref.batch`; ConQuR's `batch_ref`; `reference` is avoided, §5) | `batch_combat`; `ref_batch` also planned `batch_conqur` |
| `threshold` | a cutoff. In the networks, an association cutoff; in `decontam()` and `decontam_not()` (2026-11-24), a **score** cutoff with decontam's strict `<` (a feature is called when its score is below it; the frequency score is a heuristic, not a p-value, so not `alpha`), two values for `method = 'either'` / `'both'`; in `dbmem()` (2026-11-26), the truncation **distance** (pairs within it are neighbours; `NULL` for the minimum spanning tree's longest edge); in `replace_zeros_czm()` (2026-11-27, item 55, zCompositions' name), a zero's upper limit **in reads** (above 0, at most 1; the limit is `threshold / n`) | `assoc_network`, `assoc_pairs`, `decontam`, `decontam_not`, `dbmem`, `replace_zeros_czm` |
| `subgroups`, `same_subgroups`, `min_subgroup` | **(new, 2026-11-24)** subclasses within the `groups` classes, by `per_sample()`'s rule (LEfSe's `-s`); whether to compare only same-named subgroups (`-e`); the size under which a subgroup is compared by its median alone (`--min_c`) | `lefse` |
| `total` | **(new, 2026-11-24)** each sample is rescaled to sum to this (`NULL`: as given); `norm` takes names, not a target sum | `lefse` |
| `wilcoxon`, `wilcoxon_alpha` | **(new, 2026-11-24)** whether to run LEfSe's Wilcoxon step, and its level (LEfSe's `--wilc`, `-w`) | `lefse` |
| `lda_threshold`, `n_boots`, `boot_fraction`, `one_against_one` | **(new, 2026-11-24)** LEfSe's LDA-score cutoff (strict `>`, `-l`), bootstrap count and subsample fraction (`-b`, `-f`), and the multi-class rule (`-y`) | `lefse` |
| `features` | in `decontam_frequency_plot()` (2026-11-24), **which** features to draw (names, indices or a logical vector, as `numerator`); not `sourcetracker(features = TRUE)`'s per-feature output switch | `decontam_frequency_plot` |
| `depths`, `n_depths`, `max_depth` | **(new, 2026-11-24, item 38)** the depths of rarefaction curves (whole numbers ≥ 1; plural, since `depth` is one depth per sample), the number of points in the default grid (`n_<unit>` rule; q2's `steps`), and the largest depth computed (q2's name) | `rarefaction_curve` |
| `analytic` | **(new, item 38)** which metrics to compute in closed form: `NULL` (those that have one), `TRUE` (all, an error if one has none) or `FALSE` (resample all) | `rarefaction_curve` |
| `coverage` | **(new, item 38)** a target sample coverage in (0, 1], for coverage-based rarefaction (Chao & Jost 2012); exclusive with `depth`. Also planned `hill_standardize(coverage =)` | `rarefy` |
| `n` | **(new, 2026-11-24, item 46)** in `log_ratio_features()`, how many features to take from **each** end of a ranking (a whole number ≥ 1; Qurro's autoselection count). Not `n_features`, which the DA plots use for how many to show | `log_ratio_features` |
| `index`, `n_sites` | **(new, 2026-11-24, item 52)** the family of a partition of beta diversity (`'sorensen'`, `'jaccard'`, `'bray'`, `'ruzicka'`: betapart's `index.family`; not `family`, which is an outcome distribution), and how many samples each random draw takes (betapart's `sites`; with `n_draws`, betapart's `samples`) | `beta_partition`, `beta_partition_multi` |
| `index` (association) | **(item 50, 2026-11-24)** the association index of `indicator_species()`: `'indval_g'`, `'indval'`, `'r'`, `'r_g'` (indicspecies' `func` in snake case; not `func`, which the registries use for a function, nor `method`). The same sense as `beta_partition(index)`: which index to compute | `indicator_species` |
| `min_order`, `max_order`, `combinations` | **(new, item 50)** the smallest and largest number of groups in a combination of groups (multipatt's `min.order`, `max.order`; `max_order = 1` is its `duleg = TRUE`), and the combinations considered, by name (`'a+b'`) or id (multipatt's `restcomb`) | `indicator_species` |
| `negative` | **(new, item 50)** whether a negative association can be a feature's best (multipatt's `allow.negative`; correlation indices only) | `indicator_species` |
| `type` | the kind of plot of a `plot()` method; `'rank'` (item 46) on DA results and ordinations. On `plot.ecodive_ord()`, any value but `'samples'` and `'rank'` goes on to `plot.default()` (`type = 'n'`) | `plot.*` |
| `feature_groups`, `n_feature_groups` | **(new, 2026-11-24, item 47)** each feature's higher-level group (e.g. its phylum; the feature-side counterpart of `groups`): one value per feature, matched by name when named by every feature, else by position, `NA` = `"Unassigned"`; and how many groups to draw (`n_<unit>` rule) | `composition_barplot` |
| `parts` | **(new, 2026-11-26, item 54)** the three parts of a ternary plot (the CoDa term): a selector of exactly 3 features (names, indices or a logical vector, `feature_set()`'s forms), or a list of 3 selectors whose features are summed (amalgamation), one of which may be `NA` for every other feature | `ternary_plot` |
| `center` | **(new, item 54)** centre compositions by perturbation with the inverse of their closed geometric mean (R's spelling, as `scale(center =)`) | `ternary_plot` |
| `n_patterns` | **(new, 2026-11-27, item 60)** how many zero patterns a plot draws, the most frequent first (`n_<unit>` rule; with `n_features`, the features drawn) | `plot.ecodive_zero_pattern` |
| `ellipses` | **(new, item 54)** which regions to draw around each group in ilr coordinates: `'none'`, `'confidence'` (Hotelling's region for the mean) or `'prediction'`, at `conf_level` | `ternary_plot` |
| `grid` | **(new, item 54)** whether to draw gridlines (`TRUE`/`FALSE`) | `ternary_plot` |
| `heatmap` | **(new, item 47)** also draw a heatmap (of centred clr values) below a plot | `composition_barplot` |
| `clust` | how to order samples from their distances: an `hclust()` linkage, `'pcoa'`, or `'none'`/`FALSE`; `composition_barplot()` (item 47) orders by Aitchison distances with it, with the same values | `beta_div` and the metric functions, `vlr`, `rho_p`, `phi_s`, `sparcc`, `composition_barplot` |
| `cpus` | threads (always last) | everything threaded |


## 2. Proposed changes to existing functions (need the user's OK)

0. **`k` → `n_axes` in every ordination** (`pcoa`, `mds`, `nmds`,
   `isomap`, `tsne`, `umap`, `dbrda`, `capscale`, `pca`, `rpca`), and in
   every planned one. The user suggested this; recommended.
   - `k` breaks the `n_<unit>` rule.
   - It is overloaded in this domain: `vegan::isomap(k)` is a neighbour
     count, and "k" means k-nearest neighbours in the t-SNE, UMAP,
     Isomap and PHATE literature, beside ecodive's own `n_neighbors`.
   - The references don't agree anyway: `cmdscale`/vegan use `k`,
     scikit-learn/uwot `n_components`, Rtsne `dims`, mixOmics `ncomp`,
     ade4 `nf`.
   - `n_axes` (how many) pairs with `envfit()`'s `axes` (which).
   - All the affected functions are unreleased, so a plain rename is
     enough (an old `k =` call gets R's "unused argument" error).
   - Scope: 10 R files, their Rd pages, about 200 `k =` uses in tests
     (some are fixture fields or local variables, so each needs
     checking), `vignettes/ordination.Rmd`, NEWS, the benchmark script,
     and the plans.
1. **`isomap(k = 2L)` → `k = 3L`.** See the `k` default rule in §3.1:
   isomap's axes have eigenvalues, like pcoa's. Its default of 2 comes
   from scikit-learn's `n_components`.
2. **Drop `input` from `list_dist_ords()` and `match_dist_ord()`, and the
   `input` column from `ORDINATIONS`.**
   - It was added before `count_ord()` existed. It's always `"dist"` now
     that count-table methods have their own registry (`COUNT_ORDS`), and
     `list_count_ords()` has no such argument.
   - Removing it makes the two `list_*` / `match_*` pairs identical in
     shape.
   - It touches tests in `test-dist_ord.r` (`input = 'dist'`, and the
     "Options are: NULL, dist" error).
3. **`rarefy(seed = 0)` → `seed = 0L`,** as every other function writes
   it. This is cosmetic; behaviour doesn't change.

Nothing else in the existing signatures conflicts with the glossary.


## 3. Changes to planned signatures

### 3.1 A rule for the default number of axes (`n_axes`, formerly `k`)

Today ecodive has `k = 3` for pcoa, pca, rpca, dbrda and capscale, and
`k = 2` for mds, nmds, isomap, tsne and umap. The proposed rule:
**3 when the axes have eigenvalues (`eig` is defined), 2 when they
don't** (layouts optimised in a fixed number of dimensions).

| Method | `n_axes` default |
|:-------|:------------|
| `rda`, `cca`, `ca`, `dpcoa`, `coinertia`, `diffusion_map` | 3 |
| `dca` | **3, not 4** (DCA traditionally reports 4 axes; `n_axes = 4` still works) |
| `lda`, `capda` | `min(3, number of groups − 1)` |
| `plsda`, `oplsda`, `splsda` | 2 (components have no eigenvalues; matches mixOmics' `ncomp = 2`) |
| `phate`, `pacmap`, `largevis` | 2 |
| `principal_response_curves` | 1 (the method's convention; documented as the exception) |
| `isomap` (existing) | 3, if §2.1 is accepted |

### 3.2 Renames in the plans

| Plan | Now | Proposed | Why |
|:-----|:----|:---------|:----|
| [lda_capda.md](lda_capda.md) | `capda(n_axes =)` | `n_pcoa_axes` | `n_axes` is the number of axes *returned*; this is the number of PCoA axes *used*. |
| [procrustes.md](procrustes.md) | `procrustes(x, y, k, ...)` | `procrustes(ord, target, axes = NULL, ...)` | It's asymmetric (`ord` is rotated onto `target`), so `x`/`y` hide which is which. Choosing input axes is `axes` (as in `envfit`), not a count. Same for `procrustes_test()`. |
| [procrustes.md](procrustes.md) | `scale = TRUE` | `dilate = TRUE` | `scale` means "standardise features" everywhere else; here it means allowing a uniform dilation of `ord` (Gower's term). |
| [largevis.md](largevis.md) | `gamma = 7` | `repulsion_strength = 7` | `umap()` already calls this parameter `repulsion_strength`, and so does uwot's `lvish()`. `gamma` also means something else in `phate()`. |
| [moran_i.md](moran_i.md) | `moran_i(x, w, ...)` | `moran_i(x, space, ...)` | `weights` already means group weights in `mrpp()`, and `w` is cryptic. One name for sample positions across the spatial functions. |
| [dbmem.md](dbmem.md) | `dbmem(x, ...)` | `dbmem(space, ...)` | As above. |
| [mantel_regression.md](mantel_regression.md) | `method = 'linear' / 'rank'` | `method = 'pearson' / 'spearman'` | The same vocabulary as the rest of the Mantel family. A rank regression is the Spearman analogue; document the mapping to ecodist's `mrank`. |
| [mantel_correlogram.md](mantel_correlogram.md) | (none) | add `alternative = 'greater'` | Every other Mantel function has it. |
| [bioenv.md](bioenv.md) | `vars = NULL` | `formula = NULL` (one-sided, naming the candidates; default all numeric columns) | Variables come from `data` through `formula` everywhere else. |
| [bioenv.md](bioenv.md) | `upto = NULL` | `max_vars = NULL` | A cap, so `max_*` (vegan's `upto`). |
| [bioenv.md](bioenv.md) | `partial = NULL` | `given = NULL` | The same concept as `mantel_partial()`'s `given`. |
| [varpart_forward_select.md](varpart_forward_select.md) | `forward_select(x, scope, ...)` | `forward_select(x, formula, ...)` | As for bioenv: candidates are named in a formula over `data` (vegan's name for it is `scope`). |
| [splsda.md](splsda.md) | `keep = 10L` | `n_features = 10L` | An exact count per component (mixOmics' `keepX`). |
| [dca.md](dca.md) | `rescale = 4L` | `n_rescales = 4L` | An exact count of rescaling cycles (vegan's `iresc`). |
| [mirkat.md](mirkat.md) | `outcome = 'continuous'` | drop it; infer it from `groups` (numeric → continuous; two categories → binary; more → error) | No other test asks for the variable's type; `dist_test()` already checks it from the values. |


## 4. Inconsistencies kept on purpose

- **The order of `margin` and other arguments in the alpha and beta metric
  functions** varies, e.g. `fisher(counts, digits, margin)` and
  `generalized_unifrac(counts, tree, alpha, margin)`. These functions are
  released on CRAN, and reordering would break positional calls; the
  names are already consistent.
- **`alpha`** means UniFrac's α in `generalized_unifrac()` and
  `beta_div()`, and would mean Coifman–Lafon's density normalisation in
  `diffusion_map()`. Both are the standard names in their literatures, the
  functions are unrelated, and `diffusion_map()` never passes through
  `beta_div()`. The docs say which α is meant. In the differential
  abundance tests (`ancom()`, `dirmult_ttest()`) `alpha` is the
  significance level, as in statistics generally; those functions never
  pass through `beta_div()` either.
- **`sourcetracker(alpha1, alpha2, beta)`** keeps SourceTracker's prior
  names (2026-11-19), which users know from its documentation. None is a
  significance level, and `beta` is not beta diversity; the docs say so.
  The elastic-net mixing in `coda_lasso()` is `l1_ratio`, not glmnet's
  `alpha`, because `alpha` is already a significance level in `ancom()`,
  `forward_select()` and `mantel_regression()` (backward elimination).
- **`envfit(na.rm =)`** keeps the dot, against the snake_case rule,
  because `na.rm` is base R's universal idiom.
- **Wrapper selector names differ by design:** `beta_div(metric)`,
  `dist_ord(method)`, `count_ord(method)`, `dist_test(test)`. `test` keeps
  tests free to take a `method` argument (the correlation), as
  [README.md](README.md) explains.
- **`umap(n_epochs)` against `tsne(max_iters)`:** uwot runs exactly
  `n_epochs` (a count), while t-SNE may stop early (a cap). They follow the
  `n_` / `max_` rule.
- **Method-specific `n_neighbors` defaults** (isomap 5, umap 15, pacmap
  10, phate 5, diffusion_map 15) follow each reference.


## 5. Compositional plans (2026-11-19)

The plans in [compositional.md](compositional.md) follow §1 and add the
rows marked "(new)". Two choices worth recording:

- **`reference` is not used as an argument name.** It would mean a
  feature in `alr()` (scikit-bio's `ref_idx`) and a group in
  `dirmult_ttest()` (scikit-bio's `reference`). They become `denominator`
  and `contrast`.
- **`n_draws`, `max_power_iters`, `n_test_samples`** follow the
  `n_<unit>` / `max_<unit>` rule, renaming scikit-bio's `draws` and
  gemelli's `max_iterations_rptm` / `n_test_samples`. `ctf()` takes
  `n_starts` for gemelli's `n_initializations` (as in `mds()` / `nmds()`).


## 6. Second compositional batch (2026-11-20)

The plans added on 2026-11-20 (listed in [compositional.md](compositional.md))
reuse §1's names and add these:

| Name | Meaning | Used by (planned) |
|:-----|:--------|:------------------|
| `scaling` | the biplot scaling (`'form'`, `'covariance'`, `'contribution'`, and vegan's 1, 2, 3 for `rda`, `cca`, `ca`); used by `plot.ecodive_ord()` since item 44 | `plot.ecodive_ord` |
| `threshold` | the association cutoff that makes an edge | `assoc_network`, `assoc_pairs` (exist) and the association functions' sparsification |
| `average`, `sd` | how repeated rarefactions are combined (`'mean'`, `'median'`, `'medoid'`), and whether to return their spread | `beta_div_rarefied`, `alpha_div_rarefied` (with `rarefy()`'s `times`) |
| `pseudocounts`, `zeros`, `scale_sds` | the settings grids of a sensitivity analysis. Plural, and after `...`, so they can't capture the method's own `pseudocount` / `scale_sd` (arguments after `...` match only exactly) | `da_sensitivity` |
| `load_sd` | the uncertainty of a measured load in a scale model: the standard deviation of each log2 load's error, one value or one per sample (or a column of `data`) | `dirmult_lm` (exists, item 42; with `qmp_scale()`'s `load`) |
| `scale_model`, `se_type` | **(added 2026-11-24, item 42)** the scale model without measured loads (`'clr'`, `'tss'`: ALDEx3's `clr.sm()`, `tss.sm()`; not `scale`, which standardises features), and a regression's standard errors (`'hc3'`, `'hc0'`, `'classical'`: ALDEx3's `test`) | `dirmult_lm` (exists) |
| `n_draws` | also used by `dirmult_lm` and `sparcc` | |
| `max_exclusions`, `exclusion_threshold` | **(added 2026-11-20, item 32)** SparCC's cap on exclusion iterations (`max_<unit>` rule; `0L` for none) and the \|r\| above which a pair is excluded from the variance estimate. Not `threshold`, which is `assoc_network()`'s edge cutoff | `sparcc` (exists) |
| `sign`, `clusters`, `hubs`, `hub_quantile`, `weighted_degree` | **(added 2026-11-20, item 33)** how an association becomes an edge length and weight (`'signed'`, `'unsigned'`, `'positive'`: NetCoMi's `dissFunc`), the clustering method (`'greedy'`, `'none'`), which centralities define a hub and the quantile it must be above (NetCoMi's `hubPar`, `hubQuant`), and whether degree sums the edge weights (`weightDeg`) | `assoc_network` (exists) |
| `quantile` | **(item 33)** the quantile a "most central" feature is above (NetCoMi's `jaccQuant`) | `assoc_network_compare` (exists) |
| `assoc_args`, `network_args` | **(item 33)** lists of the arguments that made two networks' associations and networks, passed again when the group-permutation test re-estimates them (`<step>_args`) | `assoc_network_compare` (exists) |
| `n_samples` | **(item 33)** the numbers of samples behind two correlation matrices, as the `n_samples` attribute of an association result | `assoc_network_diff` (exists) |
| `y_norm` | **(added 2026-11-24, item 39; the user chose it over `norm_y`)** the transformation of the second table `y`, with `norm`'s choices; `norm` stays the first table's | `cross_assoc` (exists); planned `halla` |
| `min_samples` | **(item 39)** the fewest samples (after leaving out `NA`s) a pair needs for a correlation, at least 3 | `cross_assoc` (exists) |
| `transform`, `augment`, `standardize` | **(added 2026-11-24, item 49)** maaslin3's names: the transform after `norm` (`'log'` = log2, `'none'`), the prevalence model's data augmentation, and the scaling of numeric covariates (as `linda()`'s `standardize`) | `maaslin3` (exists) |
| `median_comparison_abundance`, `median_comparison_prevalence`, `median_comparison_abundance_threshold`, `median_comparison_prevalence_threshold`, `subtract_median`, `warn_prevalence`, `small_random_effects` | **(item 49)** maaslin3's own names, kept snake case and all because its users type them: whether each model is tested against the median coefficient, the radius around it given p = 1, whether reported coefficients have it subtracted, the flag for prevalence associations induced by abundance ones, and the grouping as a fixed effect in the prevalence model | `maaslin3` (exists) |
| `n_sims` | **(item 49)** the number of Monte Carlo draws of the median comparison (maaslin3's fixed 10,000). Not `n_draws`, which counts posterior draws that each give a full result | `maaslin3` (exists) |

Open naming questions for the user are flagged "ask" in the plans:
`bias_counts` ([secom.md](secom.md)), and the `da_test()` groups/formula
convention ([da_test.md](da_test.md)).
