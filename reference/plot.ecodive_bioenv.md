# Plot a BIO-ENV result

The correlation of the best subset of variables of each size with the
distances, from
[`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md), with
the best subset overall circled and named in the subtitle.

## Usage

``` r
# S3 method for class 'ecodive_bioenv'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md).

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `ylab`, `ylim` or `col`.

## Value

`x`, invisibly.

## Examples

``` r
    res <- bioenv(bray(iris[, 1:4]), iris[, 1:4])
    plot(res)
```
