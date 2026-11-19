# Plot a forward selection

The cumulative adjusted \\R^2\\ as
[`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md)
added each variable, with a dashed line at that of the global model (all
candidates together), and the next best variable, which was not added,
as an open point.

## Usage

``` r
# S3 method for class 'ecodive_forward_select'
plot(x, ...)
```

## Arguments

- x:

  The result of
  [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md).

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `ylab`, `ylim` or `col`.

## Value

`x`, invisibly.

## Examples

``` r
    fs <- forward_select(
      iris[, 3:4], ~ Species + Sepal.Length + Sepal.Width, iris,
      norm = 'hellinger', permutations = 99 )
    plot(fs)
```
