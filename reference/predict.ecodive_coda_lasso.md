# Predictions from a compositional lasso

The fitted model's linear predictor (the signature score) or
probabilities for new samples.

## Usage

``` r
# S3 method for class 'ecodive_coda_lasso'
predict(
  object,
  counts,
  data = NULL,
  given = NULL,
  margin = 1L,
  type = c("link", "response"),
  ...
)

# S3 method for class 'ecodive_coda_lasso_pairwise'
predict(
  object,
  counts,
  data = NULL,
  given = NULL,
  margin = 1L,
  type = c("link", "response"),
  ...
)
```

## Arguments

- object:

  The result of
  [`coda_lasso()`](https://cmmr.github.io/ecodive/reference/coda_lasso.md)
  or
  [`coda_lasso_pairwise()`](https://cmmr.github.io/ecodive/reference/coda_lasso.md).

- counts:

  A feature table of the new samples, with (at least) the fit's selected
  features, matched by name (by position when the table has no feature
  names and as many features as the fit).

- data:

  Metadata of the new samples, for a fit with a `given` formula or
  column names. Default: `NULL`

- given:

  For a fit whose `given` was a data frame: the new samples' covariates,
  as a data frame with the same columns. Default: `NULL`

- margin:

  `1` if the new samples are rows, `2` if columns. Default: `1`

- type:

  `'link'`, the linear predictor (the log odds for a binomial fit), or
  `'response'` (the probability of the second class for a binomial fit;
  the same as `'link'` for gaussian). Default: `'link'`

- ...:

  Not used.

## Value

A numeric vector, one value per new sample.

## Details

The fit's `pseudocount` is added before taking logs (with none, zeros in
the selected features are an error). Only the selected features enter,
so the others needn't be in `counts`; because the coefficients sum to
zero, each sample's total doesn't matter.

## Examples

``` r
    set.seed(1)
    counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
    y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
    fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
    predict(fit, counts[1:3, ])
#> [1] -0.2530881  0.3797884  0.1201738
```
