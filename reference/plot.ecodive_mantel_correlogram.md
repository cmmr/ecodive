# Plot a Mantel correlogram

Draws the Mantel statistic of each tested distance class against the
class's midpoint, as
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of
[`vegan::mantel.correlog()`](https://vegandevs.github.io/vegan/reference/mantel.correlog.html)
does: filled squares for the classes whose adjusted p-value is at most
`alpha`, and open squares for the others.

## Usage

``` r
# S3 method for class 'ecodive_mantel_correlogram'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  A
  [`mantel_correlogram()`](https://cmmr.github.io/ecodive/reference/mantel_correlogram.md)
  result.

- alpha:

  The significance level of the filled squares, compared with
  `p_adjusted` (with `p_value` when there are no adjusted p-values).
  Default: `0.05`.

- ...:

  Passed on to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html):
  `main`, `xlab`, `ylab`, `col`, `xlim`, `ylim`, and so on.

## Value

`x`, invisibly.

## Examples

``` r
    set.seed(1)
    xy <- matrix(runif(60, 0, 10), 30, 2)
    counts <- matrix(rpois(30 * 10, 5), 30, 10)
    counts[, 1:3] <- counts[, 1:3] + round(2 * xy[, 1])

    res <- mantel_correlogram(bray(counts), xy)
    plot(res)

    plot(res, alpha = 0.01, main = 'Distance decay')

```
