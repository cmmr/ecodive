# Plot a variation partitioning

A Venn diagram of
[`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md)'s
result: a circle (or, for four sets, an ellipse) for each set of
variables, with the adjusted \\R^2\\ of each fraction in its region, and
that of the residuals below.

## Usage

``` r
# S3 method for class 'ecodive_varpart'
plot(x, digits = 2L, col = NULL, ...)
```

## Arguments

- x:

  The result of
  [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md).

- digits:

  The number of decimal places of the fractions. Default `2`.

- col:

  The sets' colours, one per set. Defaults to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"Dark 3"` palette.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub` or `cex`.

## Value

`x`, invisibly.

## Details

The regions are all the same size, whatever their fractions: their
values are written in them, not drawn to scale. Fractions shared by
several sets can be negative, when the sets explain more together than
the sum of what they explain apart. The fractions are the
`adj_r_squared` column of `x`, and its row names (`[a]`, `[b]`, ...)
follow the regions as in
[`vegan::showvarparts()`](https://vegandevs.github.io/vegan/reference/varpart.html):
`[a]` is the first set's own fraction, `[b]` the second's, and so on.

## Examples

``` r
    vp <- varpart(
      iris[, 3:4], species = ~ Species, sepal = ~ Sepal.Length + Sepal.Width,
      data = iris, norm = 'hellinger' )
    plot(vp)


    vp <- varpart(
      iris[, 3:4], species = ~ Species, length = ~ Sepal.Length,
      width = ~ Sepal.Width, data = iris, norm = 'hellinger' )
    plot(vp)
```
