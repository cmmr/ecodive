# `random_forest()`: random forests for classification and regression

**Priority:** 62. **Tier:** port.
- **scikit-learn 1.9.1** (**BSD-3**: `COPYING` "Copyright (c) 2007-2026 The
  scikit-learn developers"; PyPI licence expression `BSD-3-Clause`; checked
  2026-11-20): the trees and forest. Files, all read at tag `1.9.1`:
  `sklearn/tree/_tree.pyx` (tree builder, `Tree`, importances),
  `_splitter.pyx` (best splitter), `_partitioner.pyx` (dense and sparse
  partitioners, the sort), `_criterion.pyx` (Gini, MSE), `_utils.pyx`,
  `sklearn/utils/_random.pxd` (`our_rand_r`), `sklearn/tree/_classes.py`,
  `sklearn/ensemble/_forest.py`, `_bootstrap.py`, `_base.py`
  (`_set_random_states`), `sklearn/model_selection/_split.py`
  (`KFold`, `StratifiedKFold`), `sklearn/metrics/_ranking.py`
  (`roc_curve`, `auc`).
- **q2-sample-classifier 2026.10.0** (**BSD-3**: `LICENSE` "Copyright (c)
  2017-2026, QIIME 2 development team", checked 2026-11-20): the workflow
  (nested cross-validation, mean importances over folds, the accuracy
  summary, ROC averaging).
- **numpy** (BSD-3) legacy `RandomState`: MT19937, `randint`'s masked
  rejection and `random_sample`, only for the parity mode (see "Same
  results"), as [ctf.md](ctf.md) ports MT19937 for gemelli.
- **Output-only:** randomForest 4.7-1.2 (**GPL (>= 2)**) and ranger 0.18.0
  (**GPL-3**) (CRAN / GitHub `DESCRIPTION`, checked 2026-11-20). Their docs
  (Rd, vignettes) and outputs may be used; their source may not.

Read [README.md](README.md) first, then [coda_lasso.md](coda_lasso.md),
whose CV-fold conventions (`n_folds`, `folds`) this plan shares.


## Context: why ecodive has this at all

On 2026-11-19 the user decided against general-purpose learners in ecodive
(ranger and XGBoost are fast, multithreaded and reproducible, and
microbiome tables are small; see TODO.md, "Source tracking and
compositional regression"). **On 2026-11-20 the user reversed this for
random forests**, to make ecodive a one-stop shop for microbiome analysis:
random forests are the machine learning microbiome papers actually use, and
users should not need a second toolkit for the most common supervised
analysis. The reversal covers random forests only: **no boosting, no neural
networks, and no hyperparameter search beyond a small grid scored by
out-of-bag error.**


## What it is

A random forest (Breiman 2001) averages many decision trees (CART;
Breiman et al. 1984), each grown on a bootstrap sample of the samples and
choosing each split among a random subset of the features. In microbiome
work it predicts a sample's class (disease vs healthy, body site) or a
number (age, a clinical score) from its feature table, ranks the features
by importance, and reports cross-validated accuracy and ROC/AUC.
QIIME 2's q2-sample-classifier wraps scikit-learn's forests for exactly
this; R users use randomForest or ranger.

Demand (retrieved 2026-11-20): Breiman (2001) has **98,118** OpenAlex
citations 2022–26 (all fields); "random forest" has **10,679** Europe PMC
microbiome full-text mentions 2022–26, about as many as DESeq2 (10,234) and
far more than MaAsLin2 (2,232) (LEfSe: 19,755).


## Names

- **`random_forest()`**: spelled out (rule 4: `rf` is cryptic and
  `forest` is vague). One function for both classification and
  regression, chosen from `y` (a factor, character or logical `y` →
  classification; numeric → regression), as randomForest, ranger and
  [coda_lasso.md](coda_lasso.md) do. scikit-learn's two classes are one
  method with two criteria, not two named variants.
- Not `*da` (rule 3): that suffix is for discriminant *ordinations*; a
  forest gives no axes.
- S3 methods: `predict.ecodive_random_forest()`,
  `plot.ecodive_random_forest()`, `print.ecodive_random_forest()`.
- **Same names in other packages:** none (`randomForest::randomForest()`,
  `ranger::ranger()`).


## What "same results as the reference" means here

A forest is random: scikit-learn with a fixed `random_state` gives one
forest, ranger with a fixed `seed` another, and ecodive's PCG streams a
third. Three options:

1. **Match scikit-learn bit for bit** by reproducing its random draws:
   numpy's MT19937 seeded by `random_state`, `randint()` per tree for the
   tree seeds, `RandomState(tree_seed).randint(0, n, n)` for the bootstrap,
   `RandomState(tree_seed).randint(0, 2^31 - 1)` to seed the splitter's
   xorshift `our_rand_r()`, and `rand_int()` for every feature draw.
2. **Match on distributions only:** accuracy, OOB error and importance
   ranks over many seeds agree with scikit-learn's over many
   `random_state`s.
3. Both.

**Recommendation: 3.** The user-facing `seed` uses PCG (README's rule:
one stream per tree, so the forest is `identical()` for any `cpus`), and
differs from scikit-learn only as a different seed would (allowed:
"differences that amount only to a different random seed are fine"). But
the **tests run the same code with numpy's draws** (an internal
`random_forest_details(..., numpy_random_state = 0L)`) and require
**identical trees** to scikit-learn's: the same split feature and
threshold at every node, the same leaf values, impurities and sample
counts, the same `feature_importances_`, `predict_proba()` and OOB
predictions (doubles at 1e-12, integers exactly). Distributional checks
alone cannot catch a wrong tie rule, an off-by-one in `min_samples_leaf`
or a criterion bug that shifts accuracy by 1%; bit-for-bit trees can. The
cost is small, because the draws are few and simple:
- MT19937 (as in [ctf.md](ctf.md); verified there against numpy over
  1,000 draws), numpy's legacy `randint` for int64 (for a range below 2^32:
  `next_uint32() & mask` with rejection while above the range; the
  "masked" method), and `random_sample` (53-bit
  `((a >> 5) * 2^26 + (b >> 6)) / 2^53`) for `choice(p = )`;
- sklearn's `our_rand_r()` (a 32-bit xorshift, `<< 13, >> 17, << 5`, a
  zero state replaced by 1, result mod 2^31) and `rand_int(low, high) =
  low + our_rand_r() % (high - low)`.

The parity mode exists only for tests and fixtures, like `ctf()`'s MT19937,
unless the user wants it exposed (**ask**: e.g. `seed = numpy(42)` or an
`rng = 'numpy'` argument, which would let q2-sample-classifier users
reproduce a QIIME 2 run exactly). Note README's PCG rule would need the
same exception as ctf.md's.

**Floating point:** scikit-learn converts `X` to **float32** before
fitting and predicting (`DTYPE = np.float32`; feature values, the
constant-feature test and thresholds `v[p-1]/2.0 + v[p]/2.0` all use the
float32 values). The port does the same (it halves `X`'s memory too).
Values that differ only beyond float32 precision (relative 6e-8) become
ties; for relative abundances and counts this never matters. Say so in
"Differences" if the user agrees (**ask**).


## Design

### Signature

```r
random_forest(counts, y, data = NULL, margin = 1L, norm = 'percent',
              pseudocount = NULL, n_trees = 500L, max_features = NULL,
              min_samples_leaf = 1L, min_samples_split = 2L,
              max_depth = NULL, criterion = NULL, class_weight = NULL,
              max_samples = NULL, importance = 'impurity',
              n_folds = 5L, folds = NULL, keep_forest = TRUE,
              seed = 0L, cpus = n_cpus())
```

| ecodive | scikit-learn | randomForest / ranger | Default here |
|:--------|:-------------|:----------------------|:-------------|
| `n_trees` | `n_estimators` (100) | `ntree` / `num.trees` (500) | **500 (ask)** |
| `max_features` | `max_features` (`'sqrt'` classifier; 1.0 regressor) | `mtry` (√p; p/3 regression) / `mtry` (√p both) | `NULL`: sklearn's (`max(1, floor(sqrt(p)))`; p for regression) **(ask)** |
| `min_samples_leaf` | `min_samples_leaf` (1) | `nodesize` (1; 5 regression) / `min.node.size` (1; 5) | 1 |
| `min_samples_split` | `min_samples_split` (2) | — | 2 |
| `max_depth` | `max_depth` (None) | `maxnodes` / `max.depth` | `NULL` (none) |
| `criterion` | `'gini'` / `'squared_error'` | Gini / variance | `NULL`: by type; `'entropy'` too (cheap) |
| `class_weight` | `class_weight` | `classwt` / `class.weights` | `NULL`; `'balanced'`, `'balanced_subsample'`, or a named numeric vector |
| `max_samples` | `max_samples` | `sampsize` / `sample.fraction` | `NULL` (n); a count or a fraction |
| `n_folds`, `folds` | q2's `cv = 5` | — | 5 |

- **`counts`, `margin`, `norm`, `pseudocount`:** any table
  `validate_counts()` reads; `norm = 'percent'` by default (q2 passes
  whatever table it is given, usually counts; relative abundances make
  samples comparable). Trees are invariant to monotone per-feature
  transforms, so `'clr'` vs `'log'` matters only through the per-sample
  centring. Any numeric table works with `norm = 'none'`. **Ask** about the
  default.
- **`y`:** a vector with one value per sample or a column name of `data`
  (`per_sample()`, as [coda_lasso.md](coda_lasso.md)). `NA` in `y`: those
  samples are dropped with a warning (q2's `missing_samples = 'ignore'`
  analogue) **or** an error (q2's default `'error'`): **ask**.
- **`max_features`:** an integer count, a fraction of p, `'sqrt'`,
  `'log2'`, or a **vector of candidates** (the grid, below).
  `min_samples_leaf` may also be a vector; the grid is their
  `expand.grid()`.
- **`importance`:** `'impurity'` (mean decrease in impurity, sklearn's
  `feature_importances_`; free), `'permutation'` (Breiman's out-of-bag
  permutation importance; costs about `n_trees × p_used × n_oob` tree
  walks), or `'both'`.
- **`keep_forest`:** store the full-data forest for `predict()`
  (memory: below).
- Bootstrap is always on (OOB error and permutation importance need it).
  sklearn's `bootstrap = FALSE`, `oob_score` (always computed here),
  `warm_start`, `max_leaf_nodes` (best-first growth), `min_impurity_decrease`,
  `min_weight_fraction_leaf`, `ccp_alpha`, `monotonic_cst` and missing-value
  support are **not ported** (`NA` in the table is an error). Add any
  later if asked.

### Algorithm: the trees (port of scikit-learn 1.9.1)

1. **Per-tree draws.** PCG mode: tree t of training set f uses
   `rng_seed(seed, f * 2^32 + t)` for its bootstrap and its splitter, so a
   tree's stream depends on neither `n_trees` nor the grid (common random
   numbers across grid points make their OOB comparison less noisy). numpy mode: forest seeds
   `MT(random_state).randint(2^31 - 1)` once per tree, in order
   (`_set_random_states()`); per tree, a fresh `MT(tree_seed)` for the
   bootstrap and another fresh `MT(tree_seed)` whose first
   `randint(0, 2^31 - 1)` seeds `our_rand_r()`.
2. **Bootstrap** (`_parallel_build_trees()`): `n_bootstrap` =
   `max_samples` (or n) draws with replacement; the tree is grown on all
   samples with **weights = draw counts** (`np.bincount`), samples with
   weight 0 left out of the tree (`Splitter.init()`). **scikit-learn 1.9.1
   changed class weights:** with `class_weight` (not
   `'balanced_subsample'`), the combined sample weights now drive a
   **weighted bootstrap** (`random_instance.choice(n, n_bootstrap, p = w /
   sum(w))`: cumulative sum, `random_sample()`, `searchsorted(side =
   'right')`), and the tree gets plain draw counts; with
   `'balanced_subsample'`, uniform draws and counts × balanced weights
   computed on the drawn sample. Earlier versions multiplied counts by the
   class weights. Pin the fixtures to 1.9.1 and say so.
   `'balanced'` weights are \eqn{n / (k \cdot n_c)}.
3. **Builder** (`DepthFirstTreeBuilder.build()`): a stack; the root
   first; each node is a leaf if depth ≥ `max_depth`, n < `min_samples_split`
   (raised to at least 2 × `min_samples_leaf`), n < 2 × `min_samples_leaf`,
   or impurity ≤ machine epsilon; otherwise split, and it is still a leaf
   if the split's position is the end or its improvement + ε < 0
   (`min_impurity_decrease = 0`). Children are pushed **right first, then
   left**, so the left subtree is grown first: the order of draws from the
   tree's generator depends on it. Nodes are numbered in creation order.
4. **Best split** (`node_split_best()`): sample up to `max_features`
   features **without replacement** with the Fisher–Yates scheme that also
   tracks constant features (known from ancestors, found at this node,
   drawn): `f_j = rand_int(n_drawn_constants, f_i - n_found_constants)`;
   known constants are swapped forward and don't count; a feature whose
   values span ≤ `FEATURE_THRESHOLD` (1e-7) at this node is newly constant
   and doesn't count either; the loop continues past `max_features` until
   at least one non-constant feature is evaluated. For each evaluated
   feature: sort the node's (float32 value, sample) pairs
   (`sort()`: introsort with median-of-three and heapsort fallback, **not
   stable**: port it exactly, since the order of tied samples changes
   floating-point sums in regression); scan positions, skipping within
   ties (`next_p()`), rejecting positions that leave fewer than
   `min_samples_leaf` samples on a side; the first position with a
   strictly larger proxy improvement wins (ties go to the earlier feature
   drawn and the smaller position); threshold = `v[p-1]/2.0 + v[p]/2.0`.
   Then partition the samples in place (the partition order feeds the
   children's sorts).
5. **Sparse input** (`SparsePartitioner`): CSC by feature; the node's
   non-zeros are split into negative / zero / positive parts and only the
   non-zeros are sorted (counts have no negatives). scikit-learn's tests
   require sparse and dense trees to be equal; port both partitioners and
   test both against the fixtures. Never densify a sparse `counts`.
6. **Criteria** (`_criterion.pyx`): Gini (classification; weighted class
   counts per node, left and right), entropy optionally, MSE
   (regression: weighted sums of y and y²). `update(new_pos)` adds samples
   to the left from `pos` or subtracts from the right from `end`,
   whichever is fewer (`reverse_reset()`): port exactly, again for the
   sums' order. `proxy_impurity_improvement()` chooses the split;
   `impurity_improvement()` (weighted, relative to the root's weighted n)
   is stored for importances.
7. **Leaf values:** classification, the weighted class fractions; regression,
   the weighted mean.

### Forest outputs

- **Predictions:** classification probabilities = the mean over trees of
  the leaf fractions (soft voting, sklearn's `predict_proba()`), summed in
  tree order; class = the first maximum in class order. randomForest votes
  by majority instead: say so. Regression: the mean of the leaf means.
- **OOB** (`_compute_oob_predictions()`): each sample averaged over the
  trees whose bootstrap missed it; OOB accuracy (classification) or R²
  (regression). A sample never out of bag: sklearn warns and divides by 1
  (prediction 0); ecodive warns and gives `NA` (a documented difference).
- **Impurity importance** (`Tree.compute_feature_importances()`,
  `BaseForest.feature_importances_`): per tree, \eqn{\sum_{nodes}
  (w_n i_n - w_l i_l - w_r i_r)} per feature over its splits, divided by the
  root's weighted n, normalised to sum 1; averaged over trees with more
  than one node; renormalised to sum 1. randomForest's MeanDecreaseGini and
  ranger's `'impurity'` are not normalised: say so. MDI favours features
  with many distinct values (Strobl et al. 2007), which is why
  `'permutation'` exists.
- **Permutation importance** (Breiman 2001; randomForest's `type = 1`,
  `scale = FALSE`; ranger's `importance = 'permutation'`; written from
  their papers and docs): per tree and per feature **used by that tree**,
  permute the feature's values among the tree's OOB samples (PCG stream per
  (tree, feature): `rng_seed(seed, 2^61 + t * p + j)`), and record the
  increase in OOB error (misclassification rate, or MSE); a feature the tree
  doesn't use contributes 0. Report the mean over all trees and its SD.
  randomForest's default scaled importance is mean / (SD / √n_trees): the
  docs show how to get it. Only the trees' walks for the OOB samples are
  repeated, and only where the permuted feature is on the path.

### Cross-validation and the grid (q2-sample-classifier's "ncv")

q2-sample-classifier's `classify_samples_ncv` / `regress_samples_ncv`
(`nested_cross_validation()`, `_fit_and_predict_cv()`): `StratifiedKFold`
(classification) or `KFold` (regression), `n_splits = cv = 5`, shuffled
with `random_state`; per fold, fit on the training part and predict the
held-out part; report the out-of-fold predictions and probabilities for
every sample, accuracy per fold (mean ± `np.std`, ddof 0), and the mean of
the folds' importances (`_mean_feature_importance()`).

- **Folds:** `folds` (one fold id per sample) overrides `seed`, as in
  [coda_lasso.md](coda_lasso.md); fixtures pass scikit-learn's folds.
  Otherwise port `StratifiedKFold._make_test_folds()`: encode classes in
  order of first appearance; `allocation[i] = bincount(sort(y)[i::K])`;
  for each class, the fold labels `repeat(0..K-1, allocation[, k])`
  shuffled (PCG Fisher–Yates on a reserved stream, `2^62`) and assigned to
  that class's samples in order. `KFold`: shuffle indices, then K
  contiguous chunks, the first `n %% K` one larger. A class with fewer
  samples than `n_folds` warns (sklearn's warning); every class smaller is
  an error. `n_folds = 0` or `1`: no CV (OOB only).
- **Also a full-data forest** (fold index K) for `predict()`, the OOB
  results and the reported importances; the CV importances are kept per
  fold for error bars.
- **The grid** (`max_features` and/or `min_samples_leaf` vectors): inside
  every training set (each CV fold and the full data), grow one forest per
  grid point and keep the one with the lowest OOB error (ties: the first
  grid point). This is nested and honest (the held-out fold never chooses
  a setting) and needs no inner CV, since OOB error is free. q2's
  `parameter_tuning` (a 20-draw `RandomizedSearchCV` over many settings,
  with inner CV) is not ported: too costly for its gain. The full-data
  choice is reported.

### Threads, memory, determinism

- **Tasks are trees:** one task per (training set, grid point, tree),
  `(n_folds + 1) × n_grid × n_trees` of them, dealt by `run_parallel()`.
  Each task has its own PCG stream (or numpy-mode seeds) and its own
  workspace; nothing it computes depends on another task, so results are
  `identical()` for any `cpus`.
- Shared, read-only: `X` as float32, dense feature-major (n × p, 4 bytes
  per value) or CSC; `y` (double); the training-set masks.
- Per thread: `samples` (n ints), `feature_values` (n floats), `features`
  and constant flags (p ints), the criterion's sums (3 × n_classes
  doubles), the stack. About 4(n + 2p) bytes plus the node arrays.
- Trees are grown into thread-local growable arrays, then copied into
  per-task slots. **Every reduction** (predictions, OOB sums, importances)
  runs after the trees, **in tree order**, per sample or per feature
  (threaded over samples / features, never over trees).
- Stored forest (`keep_forest`): per node a feature (int), threshold
  (double), left/right (int) and n_classes values (double): about
  (24 + 8k) bytes per node, ≤ 2n nodes per tree. 500 trees on 500 samples,
  2 classes: ≤ 20 MB. The CV forests are not kept.
- Cost: O(n_trees × n log n × max_features × depth) per training set;
  `threads_for()` decides how many threads small forests use.

### Result

A list of class `ecodive_random_forest`:

- `importance`: data frame `feature`, `impurity` (full-data forest),
  `impurity_cv_mean`, `impurity_cv_sd` (over folds, q2's reported value is
  the mean), and `permutation`, `permutation_sd` when computed; sorted by
  the requested importance, decreasing.
- `predictions`: data frame, one row per sample: `sample`, `observed`,
  `predicted` (out-of-fold), `fold`, `oob_predicted`, and for
  classification one probability column per class (`prob_<class>`,
  out-of-fold).
- `performance`: classification: `accuracy` (out-of-fold, overall),
  `accuracy_fold_mean`, `accuracy_fold_sd` (ddof 0, q2's), `balanced_accuracy`,
  `baseline_accuracy` (the largest class's share, q2's
  `_calculate_baseline_accuracy()`), `accuracy_ratio`, `auc` (binary) or
  `auc_macro`, `auc_micro`, `oob_accuracy`. Regression: `mse`, `r_squared`
  (out-of-fold), `mse_fold_mean`, `mse_fold_sd`, `slope`, `intercept`,
  `r` (q2's `scipy.stats.linregress` of predicted on observed), `oob_r_squared`.
- `confusion`: observed × predicted counts (classification).
- `roc`: data frame `class` (each class one-vs-rest, plus `'micro'` and
  `'macro'` for multiclass), `fpr`, `tpr`, `threshold`, from the
  out-of-fold probabilities: sklearn's `roc_curve(drop_intermediate =
  TRUE)` per class; micro = all classes binarised and pooled; macro =
  q2's `_roc_macro_average()` (every class's TPR interpolated on the union
  of FPRs and averaged); AUC by the trapezoid rule (binary AUC equals the
  Mann–Whitney statistic, ties counted ½).
- `tuning`: data frame of grid point × training set OOB errors, when there
  is a grid.
- `forest` (when `keep_forest`): the node arrays and offsets per tree,
  `classes`, `features`, and the transformation (`norm`, `pseudocount`,
  `margin`) so `predict()` can repeat it.
- Attributes: `method = 'random_forest'`, `type` (`'classification'` /
  `'regression'`), `n_trees`, `max_features` and `min_samples_leaf` (as
  used, or chosen by the grid), `criterion`, `class_weight`, `n_folds`
  (they qualify the results; rule 1).

`print()`: type, n samples × p features, the chosen settings, CV and OOB
performance, the top 10 features.

`predict.ecodive_random_forest(object, counts, type = c('class',
'prob'), margin = 1L, cpus = n_cpus())` (regression returns numbers):
features matched **by name**; features missing from `counts` are taken as
0 with a warning naming how many (count tables often lack rare taxa);
extra features are ignored. Threaded over blocks of samples; each
sample's sum over trees in tree order.


## Reference implementations and fixtures

- `data-raw/expected/random_forest.py` (scikit-learn 1.9.1, numpy ≥ 2, so
  `randint`'s default integer is int64 on every platform; record both
  versions): for each case, every tree's `tree_` arrays (`feature`,
  `threshold`, `children_left`, `children_right`, `value`, `impurity`,
  `n_node_samples`, `weighted_n_node_samples`), `estimators_samples_`,
  `feature_importances_`, `predict_proba()` / `predict()`,
  `oob_decision_function_` / `oob_prediction_`; `StratifiedKFold` and
  `KFold` folds; `roc_curve()` / `auc()`; 1,000 values each of MT19937
  `random_sample()`, `randint(0, n)` for several n, and `randint(2^31 - 1)`;
  and 1,000 values of `our_rand_r()` from a five-line Python copy of the
  BSD-3 function (it is pure; the tree-parity tests then confirm it in
  place). Text outputs under `data-raw/expected/random_forest/`,
  folded into `tests/tinytests/expected/random_forest.rds` by
  `data-raw/expected/random_forest.r`.
- Cases: binary classification (60 samples × 40 taxa, counts with many
  zeros, ties, a constant feature), three classes, regression (age-like
  y); dense and CSC; `class_weight` `'balanced'` and `'balanced_subsample'`;
  `max_samples = 0.7`; `min_samples_leaf = 3`; `max_depth = 3`;
  `max_features` 1, `'log2'`, 0.5; 20–50 trees each (enough to exercise
  the draws, small fixtures).
- **q2-sample-classifier 2026.10.0** `classify_samples_ncv` /
  `regress_samples_ncv` on one case with fixed `random_state` and the
  resulting folds: the out-of-fold predictions, probabilities and mean
  importances, which the numpy mode reproduces exactly.
- **randomForest / ranger** (output-only, in a scratch library):
  distributional comparisons only (OOB error over 20 seeds; permutation
  importance ranks by Spearman ρ between implementations' mean
  importances).


## Plotting

`plot.ecodive_random_forest(x, type = c('importance', 'roc', 'confusion',
'predictions', 'tuning'), n_features = 20, ...)`, base graphics, reading only
the result:

- `'importance'` (default): horizontal bars of the top `n_features`, the
  largest at the top, by the importance requested (`'impurity'`: the
  full-data forest's bars with the CV folds' SD as whiskers;
  `'permutation'`: mean with SD / √n_trees whiskers). Title "Random Forest
  Feature Importance"; subtitle which importance.
- `'roc'` (classification): one curve per class (one-vs-rest), plus macro
  (and micro) averages for multiclass, the diagonal, AUCs in the legend,
  as q2's `_roc_plot()`; binary: one curve. Title "Random Forest ROC
  (cross-validated)".
- `'confusion'`: a heatmap of row-normalised proportions (q2's
  `normalize = True`), counts written in the cells, observed classes as
  rows; subtitle: overall accuracy, baseline accuracy and their ratio.
- `'predictions'` (regression): out-of-fold predicted (y) against observed
  (x), the 1:1 line; subtitle R² and MSE (q2's regplot).
- `'tuning'`: OOB error per grid point, one line per training set.
- Colours: `plot_palette()` for classes, `PLOT_HIGHLIGHT` for single
  series. Returns `x` invisibly; `par()` restored; no RNG.


## Tests

- **Parity (numpy mode):** identical trees for every fixture case (node
  arrays: integers exactly, doubles at 1e-12), importances, probabilities,
  OOB predictions; sparse equals dense; the RNG streams equal the stored
  numpy / `our_rand_r()` values.
- **q2 parity:** the q2-sample-classifier case's out-of-fold predictions
  and mean importances, given its folds.
- **PCG mode:** `identical()` for `cpus` 1, 2, 3, 4, 8 (and with a grid);
  the same `seed` gives the same forest; tree t of a 30-tree forest equals
  tree t of a 50-tree forest with the same seed (streams per tree).
- **Distributional:** over 20 seeds, mean OOB accuracy and mean impurity
  importances within 3 standard errors of scikit-learn's over 20
  `random_state`s; a planted signal (two informative features among 50)
  ranks first by both importances.
- `predict()` on new data equals a base-R walk of the stored node arrays;
  feature matching by name (shuffled columns give the same predictions;
  missing features warn).
- ROC/AUC equal sklearn's `roc_curve()` / `auc()` on stored
  probabilities; macro and micro averages equal q2's.
- Edge cases: one class (error); a class smaller than `n_folds` (warning)
  or every class smaller (error); `y` with `NA` (per the decision); p = 1;
  n = 2 with `n_folds = 0`; a constant table (single-node trees: all
  importances 0, as sklearn); `max_features` > p (error); `folds` with an
  empty fold (error); a factor `y` with an unused level (dropped, with a
  warning).
- `plot()` of every type leaves `.Random.seed` and `par()` unchanged.
- Smoke lines: classification and regression, 10 trees, `n_folds = 2`.


## Docs

- Rd page in a new pkgdown section "Prediction" (with `coda_lasso()`).
- What a forest estimates; CV versus OOB; that importances are
  associations with the prediction, not effects; MDI's bias toward
  features with many distinct values; that features should not be
  filtered on the outcome before CV (leakage); that compositional effects
  still apply (link the ML caveats of item 27 in
  [compositional_docs.md](compositional_docs.md) and `coda_lasso()`).
- Recipes: q2-sample-classifier settings (`n_trees = 100`, sklearn's
  defaults, `n_folds = 5`); randomForest's scaled permutation importance;
  the microbiota maturity index ([longitudinal.md](longitudinal.md)).
- "Differences from scikit-learn, q2-sample-classifier, randomForest and
  ranger": RNG (PCG; numpy parity only in tests unless exposed); defaults;
  soft voting; normalised MDI; `NA` for never-OOB samples; no RFE
  (`optimize_feature_selection`) and no randomized search; float32.
- NEWS; parameter_names.md: `n_trees`, `max_features`,
  `min_samples_leaf`, `min_samples_split`, `max_depth`, `criterion`,
  `class_weight`, `max_samples`, `importance`, `keep_forest`.
- `inst/COPYRIGHTS` stanzas (scikit-learn BSD-3; q2-sample-classifier
  BSD-3; numpy BSD-3 for the legacy RNG; MT19937's own notice) and
  `Authors@R` credit for the scikit-learn developers if not already there
  (`src/nmds.c` already ports scikit-learn).


## Decide before writing code (ask the user)

1. Defaults: `n_trees = 500` (randomForest, ranger, and most microbiome
   papers) or 100 (scikit-learn, q2); regression `max_features`: all
   features (sklearn ≥ 1.1) or p/3 (randomForest).
2. Expose the numpy parity mode to users, or keep it internal.
3. float32 features, as scikit-learn.
4. `norm` default (`'percent'` or `'none'`).
5. `NA` in `y`: drop with a warning, or error.
6. Include `criterion = 'entropy'` (cheap) and `'balanced_subsample'`.


## Effort

High: about 2,000–2,500 lines of C (builder and `Tree` ~500, splitter and
two partitioners ~600, criteria ~300, forest, OOB, permutation importance,
folds and the grid ~500, RNGs ~150), 600 lines of R, 300 of plotting, and a
large fixture script. The port is mechanical but must be exact, so budget
time for the parity tests first.


## References

- Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32.
  doi:10.1023/A:1010933404324
- Breiman, L., Friedman, J. H., Olshen, R. A., & Stone, C. J. (1984).
  *Classification and Regression Trees*. Wadsworth.
- Pedregosa, F., et al. (2011). Scikit-learn: machine learning in Python.
  *Journal of Machine Learning Research*, 12, 2825–2830.
- Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley,
  G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning
  tools for microbiome classification and regression. *Journal of Open
  Source Software*, 3, 934. doi:10.21105/joss.00934
- Liaw, A., & Wiener, M. (2002). Classification and regression by
  randomForest. *R News*, 2(3), 18–22.
- Wright, M. N., & Ziegler, A. (2017). ranger: a fast implementation of
  random forests for high dimensional data in C++ and R. *Journal of
  Statistical Software*, 77(1), 1–17. doi:10.18637/jss.v077.i01
- Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias in
  random forest variable importance measures. *BMC Bioinformatics*, 8, 25.
  doi:10.1186/1471-2105-8-25
- Marsaglia, G. (2003). Xorshift RNGs. *Journal of Statistical Software*,
  8(14), 1–6. doi:10.18637/jss.v008.i14
- Matsumoto, M., & Nishimura, T. (1998). Mersenne twister. *ACM
  Transactions on Modeling and Computer Simulation*, 8, 3–30.
  doi:10.1145/272991.272995
- Fawcett, T. (2006). An introduction to ROC analysis. *Pattern
  Recognition Letters*, 27, 861–874. doi:10.1016/j.patrec.2005.10.010
