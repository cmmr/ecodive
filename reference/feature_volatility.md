# Feature Volatility

Which features change most predictably over time, and how:
q2-longitudinal's `feature-volatility` pipeline. A random-forest
regression of each sample's state (time point) on its features ranks the
features by importance; the most important ones are summarized (their
mean change between states, global mean, spread and so on), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
importances beside the volatility chart of a chosen feature.

## Usage

``` r
feature_volatility(
  counts,
  subjects,
  states,
  groups = NULL,
  data = NULL,
  n_features = 100L,
  test_size = 0.2,
  importance = "impurity",
  margin = 1L,
  norm = "percent",
  pseudocount = NULL,
  n_trees = 100L,
  conf_level = 0.95,
  n_resamples = 1000L,
  seed = 0L,
  cpus = n_cpus(),
  ...
)
```

## Arguments

- counts:

  A count table: samples \\\times\\ features, or any of the table types
  listed under **Input Types** (sparse tables stay sparse).

- subjects:

  Optional: the subject (individual) of each sample, for the chart's
  trajectories (q2's `individual_id_column`): one value per sample, or a
  column name of `data`. No `NA`s. `NULL` for none.

- states:

  The state (time point) of each sample, the forest's target: numeric,
  such as days or visit numbers (q2's `state_column`), one value per
  sample or a column name of `data`. No `NA`s; at least two distinct
  values.

- groups:

  Optional: the group of each sample, a mean line per group in the
  chart, given as `states` is. No `NA`s. Default: `NULL`, one group of
  every sample.

- data:

  Optional sample metadata: a data frame with one row per sample,
  holding the columns named in `subjects`, `states` and `groups`. Rows
  are matched to the samples in `counts` by row name, and rows for other
  samples are ignored; without row names, they must be in the samples'
  order. Default: `NULL`.

- n_features:

  How many of the most important features to keep (q2's
  `feature_count`): a whole number, or `Inf` for all of them. Default:
  `100`, as q2.

- test_size:

  The fraction of the samples held out of the forest, to measure its
  accuracy on (q2's `test_size`): \\\lceil\\ `test_size` \\\times n
  \rceil\\ samples, drawn at random. `0` trains on every sample
  (accuracy then only out of bag). Default: `0.2`, as q2.

- importance:

  What ranks the features: `'impurity'` (the forest's mean decrease in
  impurity, scikit-learn's `feature_importances_`, as q2; the default)
  or `'permutation'` (out-of-bag permutation importance, slower). See
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  How the forest's table is normalized: `'percent'` (relative
  abundances, the default), `'none'` (the values as given, as q2 does),
  or any other `norm` of
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md).
  The statistics and the chart always use relative abundances, as q2's
  do.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- n_trees:

  The number of trees. Default: `100`, as q2.

- conf_level:

  The confidence level of the chart's intervals: a number strictly
  between 0 and 1. Default: `0.95`, as q2's chart.

- n_resamples:

  The number of bootstrap resamples for each interval. Default: `1000`,
  as q2's chart.

- seed:

  An integer seed for the split, the forest and the bootstrap intervals.
  Default: `0`.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

- ...:

  More settings of the forest, passed to
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md):
  `max_features`, `min_samples_leaf`, `min_samples_split`, `max_depth`,
  `criterion`, `max_samples`, `n_folds` (default `0`: no
  cross-validation, as q2's fit; the held-out samples measure the
  accuracy) and `keep_forest`.

## Value

A list of class `ecodive_feature_volatility`:

- `features`:

  A data frame of the `n_features` most important features, most
  important first (q2's `feature_metadata.tsv`): `feature`,
  `importance`, `cumulative_avg_decrease` and `cumulative_avg_increase`
  (the sums of the negative and of the positive differences between the
  mean relative abundances at consecutive states), `net_avg_change`
  (their sum: the last state's mean minus the first's), and over all the
  samples `global_variance`, `global_mean`, `global_median`, `global_sd`
  and `global_cv` (SD / mean).

- `means`:

  A data frame, one row per feature, group and state: `feature`,
  `group`, `state`, `mean` (of the relative abundances), `lower` and
  `upper` (the bootstrap interval of the mean, as
  [`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)'s)
  and `n`.

- `abundances`:

  The relative abundances of those features: a samples \\\times\\
  features matrix.

- `samples`:

  A data frame, one row per sample, in the input order: `sample`,
  `subject` (with `subjects`), `group`, `state` and `set` (`"train"` or
  `"test"`).

- `predictions`:

  A data frame of the held-out samples (in the input order): `sample`,
  `observed` (the state) and `predicted` (by the forest). `NULL` with
  `test_size = 0`.

- `performance`:

  A named numeric vector: of the held-out samples' predictions, `mse`,
  `r_squared` (\\1 - SS\_{res}/SS\_{tot}\\), `slope`, `intercept` and
  `r` (of the predicted on the observed states, as q2's scatter plot);
  and `oob_r_squared`, the forest's out of bag. The held-out values are
  `NA` with `test_size = 0`.

- `forest`:

  The
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
  result, trained on the training samples: every feature's importances,
  out-of-bag predictions, and (with `keep_forest`, the default) the
  forest for [`predict()`](https://rdrr.io/r/stats/predict.html).

Attributes `conf_level` and `n_resamples` describe the intervals, and
`test_size` the held-out samples. Its title (see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md))
is "Feature Volatility".

## Details

**The pipeline**, as q2's:

1.  Hold out a random `test_size` of the samples (scikit-learn's
    `train_test_split()`, not stratified).

2.  Grow a random forest of `n_trees` regression trees of the state on
    the features of the other samples, and predict the held-out samples'
    states, to see how well the features track time.

3.  Rank the features by the forest's importance and keep the
    `n_features` most important.

4.  Summarize those features' relative abundances (of the whole table,
    every sample): q2's `_summarize_feature_stats()`. The mean at each
    state is over every sample at that state, whatever its group.

5.  For the chart, each feature's mean and bootstrap interval for each
    group at each state, as
    [`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)
    computes them.

**Differences from q2-longitudinal**:

- **No recursive feature elimination.** q2 first runs scikit-learn's
  `RFECV` (cross-validated recursive feature elimination) and keeps only
  the features it selects, often a handful, so its importances are those
  of a forest on those features alone; ecodive ranks every feature by a
  forest on all of them. Its top features are much the same, and
  `n_features` sets how many to keep.

- **Random draws** (the split, the forest and the intervals) are
  ecodive's own seeded PCG streams, not numpy's or JavaScript's, so they
  differ from q2's as one `random_state` does from another. Given q2's
  split and `random_state`, the forest and its importances are the same
  as q2's (see
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)).

- **The forest's table** is relative abundances by default; q2 grows it
  on the table it is given (usually counts): use `norm = 'none'` for
  that.

- **Features absent from every training sample** are kept, with
  importance 0; q2 drops them. Ties in importance are ranked in the
  table's feature order.

- **`global_cv`** is the standard deviation over the mean, as q2's
  "Global CV (%)" column, which despite its name is not multiplied by
  100.

- q2's chart lets the viewer pick the feature and the grouping
  interactively; here
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
  feature per call, and the grouping is set by `groups`.

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

Bokulich, N. A., Dillon, M. R., Zhang, Y., Rideout, J. R., Bolyen, E.,
Li, H., Albert, P. S., & Caporaso, J. G. (2018). q2-longitudinal:
longitudinal and paired-sample analyses of microbiome data. *mSystems*,
3, e00219-18.
[doi:10.1128/mSystems.00219-18](https://doi.org/10.1128/mSystems.00219-18)

Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley, G.
A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning
tools for microbiome classification and regression. *Journal of Open
Source Software*, 3, 934.
[doi:10.21105/joss.00934](https://doi.org/10.21105/joss.00934)

Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5-32.
[doi:10.1023/A:1010933404324](https://doi.org/10.1023/A%3A1010933404324)

## See also

[`plot.ecodive_feature_volatility()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_feature_volatility.md),
[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md),
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)

Other longitudinal:
[`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md),
[`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md),
[`nmit()`](https://cmmr.github.io/ecodive/reference/nmit.md),
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)

## Examples

``` r
    # 8 subjects sampled on days 0, 7, 14 and 21; feature 1 rises.
    meta <- data.frame(
      subject = rep(sprintf('s%d', 1:8), each = 4),
      day     = rep(c(0, 7, 14, 21), 8),
      diet    = rep(c('x', 'y'), each = 16) )
    counts <- outer(seq_len(32), 1:6, function (i, j) (i * j) %% 11 + 5)
    counts[, 1] <- counts[, 1] + 3 * meta$day
    colnames(counts) <- sprintf('F%d', 1:6)

    fv <- feature_volatility(counts, 'subject', 'day', 'diet', data = meta,
                             n_trees = 20, cpus = 1)
    fv
#> Feature Volatility
#> Random forest of 20 trees, trained on 25 of 32 samples; on the 7 held out: MSE 4.6, R^2 0.93, r 0.984; OOB R^2 0.844
#> The 6 most important of 6 features, with relative abundances over time:
#> # A data.frame: 6 × 10
#>   feature importance cumulative_avg_decrease cumulative_avg_increase
#>   <chr>        <dbl>                   <dbl>                   <dbl>
#> 1 F1         0.806                    0                     0.423
#> 2 F6         0.0613                  -0.0986                0
#> 3 F3         0.0566                  -0.107                 0
#> 4 F4         0.0527                  -0.0708                0
#> 5 F5         0.016                   -0.0715                0.000586
#> 6 F2         0.00694                 -0.0756                0
#> # … 6 more columns: net_avg_change <dbl>, global_variance <dbl>,
#> #   global_mean <dbl>, global_median <dbl>, global_sd <dbl>, global_cv <dbl>
#> Means by state with 95% bootstrap intervals in $means; trajectories in $abundances
    plot(fv)
```
