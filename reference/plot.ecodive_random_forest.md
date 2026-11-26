# Plot a random forest

Five views of a
[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
result, drawn from the result alone:

## Usage

``` r
# S3 method for class 'ecodive_random_forest'
plot(
  x,
  type = c("importance", "roc", "confusion", "predictions", "tuning"),
  n_features = 20,
  importance = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
  result.

- type:

  Which view: `'importance'`, `'roc'`, `'confusion'`, `'predictions'` or
  `'tuning'`. Default: `'importance'`

- n_features:

  For `'importance'`: how many features to show. Default: `20`

- importance:

  For `'importance'`: `'impurity'` or `'permutation'`; `NULL` for the
  one the result was sorted by. Default: `NULL`

- ...:

  Passed on to the base plotting function (e.g. `main`, `xlab`, `col`).

## Value

`x`, invisibly.

## Details

- `'importance'` (the default): horizontal bars of the `n_features` most
  important features, largest at the top. With impurity importance, the
  full-data forest's values, with whiskers of the cross-validation
  folds' standard deviation; with permutation importance, the mean over
  trees, with whiskers of its standard error (SD / sqrt(trees)).

- `'roc'` (classification): the ROC curve of each class against the rest
  (one curve for two classes), and for more classes the micro and macro
  averages, with AUCs in the legend, as q2-sample-classifier draws them.
  Cross-validated (out of bag without cross-validation).

- `'confusion'` (classification): a heatmap of the confusion matrix,
  each observed class's row as proportions (shading) and counts
  (written).

- `'predictions'` (regression): cross-validated (or out-of-bag)
  predictions against the observed values, with the 1:1 line.

- `'tuning'` (with a grid): each grid point's out-of-bag error, one line
  per training set.

## See also

[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)

## Examples

``` r
    counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
    site   <- rep(rownames(ex_counts), each = 5)
    fit    <- random_forest(counts, site, n_trees = 50, cpus = 1)
    plot(fit)

    plot(fit, 'roc')

    plot(fit, 'confusion')
```
