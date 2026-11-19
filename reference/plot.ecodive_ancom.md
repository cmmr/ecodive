# Plot an ANCOM result

A dot chart of
[`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md)'s \\W\\
statistic for the features with the largest, with the features declared
differentially abundant (`reject = TRUE`) highlighted, and a line at the
cutoff that separates them.

## Usage

``` r
# S3 method for class 'ecodive_ancom'
plot(x, n_features = 30, ...)
```

## Arguments

- x:

  The result of
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md).

- n_features:

  How many features to show, those with the largest \\W\\. Default `30`;
  `Inf` for all.

- ...:

  Passed on to
  [`graphics::dotchart()`](https://rdrr.io/r/graphics/dotchart.html),
  such as `main`, `xlab`, `cex` or `col` (its `color`, replacing the
  highlighting); and `sub`, to
  [`graphics::title()`](https://rdrr.io/r/graphics/title.html).

## Value

`x`, invisibly.

## Details

A feature's \\W\\ is the number of other features whose log ratio with
it differs significantly between the groups, out of \\D - 1\\ for \\D\\
features. ANCOM rejects the features with \\W\\ above `cutoff` \\\times
D\\ (the dashed line). There is no line when the largest \\W\\ was below
`theta`, so that nothing was rejected.

## Examples

``` r
    set.seed(1)
    counts <- matrix(rpois(200, 20), nrow = 10)
    counts[1:5, 1:3] <- counts[1:5, 1:3] * 8
    colnames(counts) <- paste0('F', 1:20)
    res <- ancom(counts, rep(c('a', 'b'), each = 5))
    plot(res)
```
