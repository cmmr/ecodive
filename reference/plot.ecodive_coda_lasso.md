# Plot a compositional lasso

`type = 'signature'` draws coda4microbiome's two standard figures side
by side: the selected features' coefficients as bars (numerator,
positive, in red; denominator, negative, in blue), and the signature
score (the linear predictor) by outcome: a box per class for binomial
`y`, against the observed value for gaussian `y`. `type = 'path'` draws
the coefficients along the \\\lambda\\ path, and `type = 'cv'` the
cross-validated error with its standard errors; both mark the chosen
\\\lambda\\.

## Usage

``` r
# S3 method for class 'ecodive_coda_lasso'
plot(x, type = c("signature", "path", "cv"), n_features = 30, ...)

# S3 method for class 'ecodive_coda_lasso_pairwise'
plot(x, type = c("signature", "path", "cv"), n_features = 30, ...)
```

## Arguments

- x:

  The result of
  [`coda_lasso()`](https://cmmr.github.io/ecodive/reference/coda_lasso.md)
  or
  [`coda_lasso_pairwise()`](https://cmmr.github.io/ecodive/reference/coda_lasso.md).

- type:

  `'signature'`, `'path'` or `'cv'`. Default: `'signature'`

- n_features:

  The most features to show in the signature's bars (the largest
  \\\|\beta\|\\). Default: `30`

- ...:

  Passed to the plotting functions (`main`, `xlab`, `col`, ...).

## Value

`x`, invisibly.

## Details

coda4microbiome rescales the coefficients to \\\sum\|\beta\| = 2\\
before drawing them; that changes neither their signs nor their order.

## Examples

``` r
    set.seed(1)
    counts <- matrix(rpois(240, 50), 30, 8, dimnames = list(NULL, LETTERS[1:8]))
    y      <- log(counts[, 'A'] / counts[, 'B']) + rnorm(30, sd = 0.1)
    fit    <- coda_lasso(counts, y, n_folds = 5, cpus = 1)
    plot(fit)

    plot(fit, type = 'path')

    plot(fit, type = 'cv')
```
