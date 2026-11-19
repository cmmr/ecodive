# Plot a Dirichlet-multinomial t-test

A volcano plot of
[`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)'s
results: each feature's log2 fold change against its adjusted p-value,
with the features that differ (those with `reject = TRUE`) highlighted,
and the most significant of them named.

## Usage

``` r
# S3 method for class 'ecodive_dirmult_ttest'
plot(x, labels = 10, ...)
```

## Arguments

- x:

  The result of
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md).

- labels:

  How many of the differentially abundant features to name, those with
  the smallest adjusted p-values. Default `10`; `0` for none.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `xlim` or `cex`.

## Value

`x`, invisibly.

## Details

The y axis is \\-\log\_{10}\\ of `p_adjusted` (an adjusted p-value of 0
is drawn at the smallest positive double). The dashed line marks
`alpha`; a feature also needs its confidence interval for the log2 fold
change to exclude 0 to count as differentially abundant, so some
features above the line may not be highlighted.

## Examples

``` r
    res <- dirmult_ttest(round(iris[, 1:4] * 10), 'Species', iris)
    plot(res)
```
