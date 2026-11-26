# Random forests

Grows a random forest (Breiman 2001) of classification or regression
trees on a feature table: to predict each sample's class (disease or
health, body site) or a number (age, a clinical score), to rank the
features by importance, and to estimate how well that works, by
cross-validation and out of bag. The trees are an exact port of
scikit-learn's, and the workflow is q2-sample-classifier's.

## Usage

``` r
random_forest(
  counts,
  y,
  data = NULL,
  margin = 1L,
  norm = "percent",
  pseudocount = NULL,
  n_trees = 500L,
  max_features = NULL,
  min_samples_leaf = 1L,
  min_samples_split = 2L,
  max_depth = NULL,
  criterion = NULL,
  class_weight = NULL,
  max_samples = NULL,
  importance = "impurity",
  n_folds = 5L,
  folds = NULL,
  keep_forest = TRUE,
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of feature values (samples \\\times\\ features), or
  any of the table types listed under "Input Types" (sparse tables stay
  sparse). Typically counts.

- y:

  What to predict: a vector with one value per sample, or the name of a
  column of `data`. A factor, character or logical `y` makes a
  classification forest (the classes are its values, or a factor's
  levels); a numeric `y`, a regression forest. Samples whose `y` is `NA`
  are left out, with a warning.

- data:

  A data frame of sample metadata, with one row per sample (matched to
  `counts` by row names when both have them), when `y` names one of its
  columns. Default: `NULL`

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  How to normalize each sample before growing the trees: `'percent'`
  (relative abundances, the default), `'none'` (the values as given, for
  any numeric table), or one of
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)'s
  other choices. Trees don't change under monotone transforms of a
  feature, so only transforms across features (percent, clr, chord, ...)
  matter. Default: `'percent'`

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- n_trees:

  The number of trees in each forest. Default: `500`

- max_features:

  How many features each split chooses among: a count, a fraction of the
  features (below 1), `'sqrt'`, `'log2'`, or `'all'`; `NULL` for
  scikit-learn's defaults (`'sqrt'` for classification, `'all'` for
  regression). A vector of these makes a grid (see Details). Default:
  `NULL`

- min_samples_leaf:

  The fewest (distinct) samples a leaf may hold. A vector of values
  makes a grid (see Details). Default: `1`

- min_samples_split:

  The fewest (distinct) samples a node needs to be split; at least twice
  `min_samples_leaf` is used. Default: `2`

- max_depth:

  The deepest a tree may grow, or `NULL` for no limit. Default: `NULL`

- criterion:

  The impurity to split on: `'gini'` (the default for classification) or
  `'entropy'`; `'squared_error'` for regression. Default: `NULL`

- class_weight:

  Weights for the classes (classification only): `NULL` (equal),
  `'balanced'` (\\n / (k n_c)\\ for a class of \\n_c\\ of the \\n\\
  training samples, \\k\\ classes), `'balanced_subsample'` (the same,
  computed on each tree's bootstrap sample), or a named numeric vector
  of weights (classes not named get 1). As in scikit-learn 1.9, weights
  other than `'balanced_subsample'` drive a weighted bootstrap (samples
  are drawn in proportion to their class weight). Default: `NULL`

- max_samples:

  The size of each tree's bootstrap sample: `NULL` (as many draws as
  training samples), a count, or a fraction (below 1) of the training
  samples (of their total class weight, with `class_weight`). Default:
  `NULL`

- importance:

  Which feature importances to compute: `'impurity'` (mean decrease in
  impurity, free), `'permutation'` (out-of-bag permutation importance),
  or `'both'`. Default: `'impurity'`

- n_folds:

  The number of cross-validation folds; `0` or `1` for none (out-of-bag
  estimates only). Default: `5`

- folds:

  A fold for each sample instead (any labels; each label is a fold),
  which overrides `n_folds` and `seed` for the folds. Matched to the
  samples like `y`. Default: `NULL`

- keep_forest:

  Keep the full-data forest in the result, for
  [predict()](https://cmmr.github.io/ecodive/reference/predict.ecodive_random_forest.md).
  Default: `TRUE`

- seed:

  An integer seed for the bootstraps, the trees' feature draws, the
  folds and the permutations. The same `seed` gives the same forest for
  any `cpus`. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A list of class `ecodive_random_forest`:

- `importance`:

  A data frame, one row per feature, sorted by the requested importance
  (largest first): `feature`, `impurity` (the full-data forest's),
  `impurity_cv_mean` and `impurity_cv_sd` (over the cross-validation
  folds' forests; q2 reports the mean), and, with
  `importance = 'permutation'` or `'both'`, `permutation` and
  `permutation_sd`.

- `predictions`:

  A data frame, one row per sample: `sample`, `observed`, `predicted`
  (cross-validated: by the forest of the fold that held the sample out),
  `fold`, `oob_predicted`, and for classification one cross-validated
  probability per class (`prob_<class>`).

- `performance`:

  A named numeric vector. Classification: `accuracy` (cross-validated),
  `accuracy_fold_mean`, `accuracy_fold_sd` (over folds, `sd` with
  divisor n, as q2), `balanced_accuracy`, `baseline_accuracy` (the
  largest class's share), `accuracy_ratio`, `auc` (two classes) or
  `auc_macro` and `auc_micro`, and `oob_accuracy`. Regression: `mse`,
  `r_squared`, `mse_fold_mean`, `mse_fold_sd`, `slope`, `intercept` and
  `r` (of the predicted on the observed values, as q2), and
  `oob_r_squared`. Without cross-validation, only the OOB values are not
  `NA`.

- `confusion`:

  Classification: the observed (rows) by predicted (columns) counts,
  cross-validated (out of bag without cross-validation).

- `roc`:

  Classification: a data frame `class`, `fpr`, `tpr`, `threshold`:
  scikit-learn's `roc_curve()` of each class against the rest (only the
  second class with two classes), and for more than two classes q2's
  `'micro'` and `'macro'` averages; from the cross-validated
  probabilities (out of bag without cross-validation).

- `tuning`:

  With a grid: a data frame of `training_set` (`'full'` or
  `'fold <k>'`), `max_features`, `min_samples_leaf`, `oob_error` and
  `chosen`.

- `forest`:

  With `keep_forest`: the full-data forest, for
  [predict()](https://cmmr.github.io/ecodive/reference/predict.ecodive_random_forest.md).

Attributes: `method` (`'random_forest'`), `type` (`'classification'` or
`'regression'`), `n_trees`, `max_features` and `min_samples_leaf` (of
the full-data forest: as chosen, with a grid), `criterion`,
`class_weight` and `n_folds`.

## Details

### The trees

Each tree is grown on a bootstrap sample of the training samples (drawn
with replacement; the number of draws of each sample is its weight),
splitting each node at the threshold, over a random subset of
`max_features` features, that most reduces the weighted impurity: Gini
or entropy for classes, the squared error for numbers. Trees grow until
their leaves are pure or reach `min_samples_leaf`, `min_samples_split`
or `max_depth`. A tree predicts the (weighted) class fractions or mean
of the leaf a sample falls in, and the forest averages its trees: class
probabilities (soft voting, as scikit-learn), whose largest gives the
predicted class (the first class on ties), or the mean of the trees'
predictions.

This is scikit-learn 1.9.1's `RandomForestClassifier` and
`RandomForestRegressor`, ported to C: given the same random draws, the
trees are identical (the tests check every node against scikit-learn's).
The draws themselves come from ecodive's PCG generator, one stream per
tree, so a forest is the same on any number of threads, and differs from
scikit-learn's as a different seed would.

### Out of bag and cross-validated

About a third of the samples are left out of each tree's bootstrap.
Predicting each sample from only the trees that never saw it gives
**out-of-bag** (OOB) predictions and accuracy (or \\R^2\\) for free.
**Cross-validation** grows a separate forest on each fold's training
samples and predicts its held-out samples (q2-sample-classifier's
`classify-samples-ncv` and `regress-samples-ncv`): with `n_folds` folds,
stratified by class for classification (as scikit-learn's
`StratifiedKFold`), or plain (`KFold`) for regression. Both are honest
estimates of how well the forest predicts new samples; cross-validation
also gives ROC curves, a confusion matrix and fold-to-fold variation. A
full-data forest is also grown: it gives the OOB results, the reported
importances and
[predict()](https://cmmr.github.io/ecodive/reference/predict.ecodive_random_forest.md).

### Choosing settings

A vector of `max_features` or `min_samples_leaf` values makes a grid
(all their combinations). Within each training set (each fold's, and the
full data), a forest is grown for every grid point and the one with the
lowest out-of-bag error is kept (the first on ties). The held-out
samples never choose a setting, so the cross-validated estimates stay
honest. The `tuning` element lists every grid point's OOB error.

### Importance

`'impurity'` importance (mean decrease in impurity, scikit-learn's
`feature_importances_`) is each feature's share of the trees' total
impurity decrease, averaged over trees and scaled to sum to 1. It is
biased toward features with many distinct values (Strobl et al. 2007).
`'permutation'` importance (Breiman 2001) permutes each feature among
each tree's out-of-bag samples and records how much the tree's
out-of-bag error (misclassification rate, or mean squared error) grows;
it is the mean over all trees (a tree that doesn't split on the feature
adds 0), with its standard deviation over trees. Importances measure
association with the prediction, not effects; correlated features share
it.

Features should not be filtered by their association with `y` before
cross-validation: that leaks the held-out samples into the model and
inflates the estimates. Relative abundances are compositional, so a
feature's importance can reflect changes in others; see "Machine
learning" below on transforms (Yerke et al. 2024 found relative
abundances did best with random forests).

### Recipes

- **q2-sample-classifier's settings:** `n_trees = 100`, `n_folds = 5`,
  `norm = 'none'` (q2 passes the table it's given).

- **randomForest's scaled permutation importance** (its default
  `importance(fit)`): `importance = 'permutation'`, then
  `permutation / (permutation_sd / sqrt(n_trees))`.

- **A microbiota maturity index** (Subramanian et al. 2014): a
  regression forest of age on healthy children's samples, then
  [predict()](https://cmmr.github.io/ecodive/reference/predict.ecodive_random_forest.md)
  for the others; their relative maturity is the predicted minus a
  spline of predicted age on age in the healthy children.

## Differences from scikit-learn, q2-sample-classifier, randomForest and ranger

- **Random draws:** PCG streams, not numpy's MT19937, so forests differ
  from scikit-learn's as with a different `random_state`. (A parity mode
  that reproduces scikit-learn's draws is used by the tests.)

- **Defaults:** 500 trees (randomForest's and ranger's; scikit-learn's
  and q2's 100), with scikit-learn's `max_features` (all features for
  regression, where randomForest uses a third) and leaf sizes (1 for
  regression, where randomForest and ranger use 5). Tables are relative
  abundances by default (`norm = 'percent'`).

- **Soft voting:** classes are predicted from the mean of the trees'
  class probabilities, as scikit-learn; randomForest counts votes.

- **Importances:** impurity importances sum to 1, as scikit-learn's;
  randomForest's `MeanDecreaseGini` and ranger's `'impurity'` are not
  scaled. Permutation importance is unscaled (randomForest's
  `scale = FALSE`, ranger's).

- **Never out of bag:** a sample that every bootstrap drew gets `NA` as
  its OOB prediction, with a warning, and is left out of the OOB scores;
  scikit-learn predicts 0 for it.

- **Folds:** as scikit-learn's `StratifiedKFold` (classification) and
  `KFold` (regression), shuffled by PCG instead of numpy. All features
  are kept in every fold (q2's DictVectorizer drops features absent from
  a fold's training samples).

- **Float32:** values are rounded to single precision before growing and
  predicting, as scikit-learn does; values that differ only beyond that
  (a relative 6e-8) are ties.

- **Not ported:** q2's recursive feature elimination
  (`optimize_feature_selection`) and randomized parameter search, and
  scikit-learn's `bootstrap = FALSE`, `max_leaf_nodes`,
  `min_impurity_decrease`, `min_weight_fraction_leaf`, `ccp_alpha`,
  `monotonic_cst` and missing values (`NA` in `counts` is an error).

## Machine learning

Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
and compositional geometry; don't expect them to improve prediction.
With random forests predicting 65 metadata variables in four datasets,
Yerke et al. (2024) found that alr, clr and ilr (PhILR) "performed
generally slightly worse or only as well as compositionally naïve
transformations", and relative abundances did best, by a small margin.
Compare with `norm = 'percent'` (see
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md))
before choosing a log-ratio predictor table.

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## Pseudocount

The `pseudocount` parameter is only relevant when `norm = 'clr'` or
`norm = 'log'`.

Zeros are undefined in the centered log-ratio (CLR) transformation and
in the log of proportions. If `norm` is `'clr'` or `'log'`,
`pseudocount` is `NULL` (the default), and zeros are detected, the
function uses half the minimum non-zero value (`min(x[x>0]) / 2`) and
issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`), or
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
whose result needs no pseudocount.

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

See [`aitchison`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5-32.
[doi:10.1023/A:1010933404324](https://doi.org/10.1023/A%3A1010933404324)

Breiman, L., Friedman, J. H., Olshen, R. A., & Stone, C. J. (1984).
*Classification and Regression Trees*. Wadsworth.

Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
*Journal of Machine Learning Research*, 12, 2825-2830.

Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley, G.
A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning
tools for microbiome classification and regression. *Journal of Open
Source Software*, 3, 934.
[doi:10.21105/joss.00934](https://doi.org/10.21105/joss.00934)

Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias
in random forest variable importance measures. *BMC Bioinformatics*, 8,
25. [doi:10.1186/1471-2105-8-25](https://doi.org/10.1186/1471-2105-8-25)

Subramanian, S., et al. (2014). Persistent gut microbiota immaturity in
malnourished Bangladeshi children. *Nature*, 510, 417-421.
[doi:10.1038/nature13421](https://doi.org/10.1038/nature13421)

Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
normalizations outperform compositional data transformations in machine
learning applications. *Microbiome*, 12, 45.
[doi:10.1186/s40168-023-01747-z](https://doi.org/10.1186/s40168-023-01747-z)

Fawcett, T. (2006). An introduction to ROC analysis. *Pattern
Recognition Letters*, 27, 861-874.
[doi:10.1016/j.patrec.2005.10.010](https://doi.org/10.1016/j.patrec.2005.10.010)

## See also

[`predict.ecodive_random_forest()`](https://cmmr.github.io/ecodive/reference/predict.ecodive_random_forest.md),
[`plot.ecodive_random_forest()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_random_forest.md),
[`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md),
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md)

## Examples

``` r
    # Which body site? 20 samples from 4 sites, 5-fold cross-validation.
    counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
    site   <- rep(rownames(ex_counts), each = 5)
    fit    <- random_forest(counts, site, n_trees = 50, cpus = 1)
    fit
#> # Random forest classification: 20 samples x 6 features
#> # 50 trees; max_features = 2, min_samples_leaf = 1, criterion = gini
#> # Classes: Gums, Nose, Saliva, Stool
#> # 5-fold CV accuracy: 1.00 (folds 1.00 +/- 0; baseline 0.250, ratio 4.00)
#> # CV AUC: macro 1.00, micro 1.00
#> # Out-of-bag accuracy: 1.00
#> # Top 6 features by impurity importance:
#>            feature impurity
#>      Streptococcus   0.3300
#>        Haemophilus   0.1920
#>  Propionibacterium   0.1510
#>     Staphylococcus   0.1480
#>        Bacteroides   0.1090
#>    Corynebacterium   0.0702

    # Predict new samples.
    predict(fit, ex_counts, cpus = 1)
#> Saliva   Gums   Nose  Stool 
#> Saliva   Gums   Nose  Stool 
#> Levels: Gums Nose Saliva Stool
```
