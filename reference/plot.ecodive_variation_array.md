# Plot a variation array

Draws a
[`variation_array()`](https://cmmr.github.io/ecodive/reference/variation_array.md)
as two heatmaps on one grid: the variances of the log ratios above the
diagonal, their means below, and beside each row the feature's share of
the total variance.

## Usage

``` r
# S3 method for class 'ecodive_variation_array'
plot(x, digits = 2L, key = TRUE, ...)
```

## Arguments

- x:

  A variation array, from
  [`variation_array()`](https://cmmr.github.io/ecodive/reference/variation_array.md).

- digits:

  The number of decimal places of the values written in the cells, when
  there are at most 15 features. Default `2`.

- key:

  Whether to draw the two colour keys on the right. Default `TRUE`.

- ...:

  Passed on to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html), such as
  `main`, `sub`, or `useRaster`.

## Value

`x`, invisibly.

## Details

The features run from the first at the top to the last at the bottom,
and from left to right in the same order. Their order is the one `x` was
made with (its `clust`), or else average-linkage clustering of the
distances \\\sqrt{\tau\_{ij}}\\, so that nearly proportional features
sit together. `x` itself isn't changed.

Cells above the diagonal are coloured by \\\tau\_{ij}\\, from light
(proportional) to dark red; cells below by the mean log ratio of the
row's feature over the column's, blue for negative and red for positive,
symmetric about 0. The bars on the right are each feature's clr variance
as a percentage of the total variance ("% clr var", as in
Pawlowsky-Glahn & Egozcue 2011).

## See also

[`variation_array()`](https://cmmr.github.io/ecodive/reference/variation_array.md),
[`plot.ecodive_dist()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_dist.md)
for [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) alone.

## Examples

``` r
    x <- replace_zeros(ex_counts, 'mult')
    plot(variation_array(x))

```
