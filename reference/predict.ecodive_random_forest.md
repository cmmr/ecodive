# Predict with a random forest

Predicts the class (or class probabilities), or the value, of new
samples from the full-data forest of
[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md).

## Usage

``` r
# S3 method for class 'ecodive_random_forest'
predict(
  object,
  counts,
  type = c("class", "prob"),
  margin = 1L,
  cpus = n_cpus(),
  ...
)
```

## Arguments

- object:

  A
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
  result, grown with `keep_forest = TRUE`.

- counts:

  The new samples' feature table, of any type
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
  takes.

- type:

  For classification: `'class'` (the predicted classes, a factor) or
  `'prob'` (a samples x classes matrix of probabilities). Ignored for
  regression, which returns numbers. Default: `'class'`

- margin:

  As for
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md):
  `1` if samples are rows of `counts`, `2` if columns. Default: `1`

- cpus:

  The number of CPU threads. Default:
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md)

- ...:

  Not used.

## Value

A factor of classes, a matrix of class probabilities, or a numeric
vector, named by sample.

## Details

Features are matched by name. Features of the forest missing from
`counts` are taken as 0, with a warning saying how many (count tables
often lack rare features); features the forest doesn't know are dropped.
The table is then normalized as the forest's was (same `norm` and
`pseudocount`), over the forest's features only.

## See also

[`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)

## Examples

``` r
    counts <- ex_counts[rep(1:4, each = 5), ] + outer(1:20, 1:6) %% 7
    site   <- rep(rownames(ex_counts), each = 5)
    fit    <- random_forest(counts, site, n_trees = 50, n_folds = 0, cpus = 1)
    predict(fit, ex_counts, type = 'prob', cpus = 1)
#>        Gums Nose Saliva Stool
#> Saliva 0.38 0.08   0.50  0.04
#> Gums   0.98 0.02   0.00  0.00
#> Nose   0.00 1.00   0.00  0.00
#> Stool  0.06 0.02   0.02  0.90
```
