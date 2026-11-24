# Plot Associations Between Two Tables

The heatmap of a
[`cross_assoc()`](https://cmmr.github.io/ecodive/reference/cross_assoc.md)
result: one row per feature of `x`, one column per variable of `y`,
coloured by the correlation, with stars on the cells whose adjusted
p-values are small.

## Usage

``` r
# S3 method for class 'ecodive_cross_assoc'
plot(
  x,
  n_x = 30,
  n_y = 30,
  clust = "average",
  stars = c(0.05, 0.01, 0.001),
  col = NULL,
  key = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`cross_assoc()`](https://cmmr.github.io/ecodive/reference/cross_assoc.md)
  result.

- n_x, n_y:

  The most features of `x` (rows) and variables of `y` (columns) to
  show: those with the smallest adjusted p-value (ties broken by the
  largest \\\|r\|\\). Default: `30` each.

- clust:

  `'average'` (the default) orders the rows and the columns by
  average-linkage clustering
  ([`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)) on \\1 -
  r\\ between their profiles of correlations in the shown block; any
  other method of
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) also works.
  `'none'` keeps the input order.

- stars:

  Thresholds of the adjusted p-values for `*`, `**`, `***` (any number
  of them, largest first). Default: `c(0.05, 0.01, 0.001)`. `NULL` marks
  the cells with `reject` by a dot instead.

- col:

  Colours for the scale from -1 to 1. Default: the diverging
  `hcl.colors(65, "Blue-Red 3")`.

- key:

  Whether to draw the colour key. Default: `TRUE`.

- ...:

  Passed on to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html), such as
  `main`, `xlab` and `ylab`.

## Value

`x`, invisibly.

## Examples

``` r
    set.seed(1)
    counts <- matrix(rpois(30 * 8, 20), 30, 8,
      dimnames = list(NULL, paste0('Taxon', 1:8)))
    env <- data.frame(pH = counts[, 1] / 10 + rnorm(30), temp = rnorm(30))
    plot(cross_assoc(counts, env))

```
