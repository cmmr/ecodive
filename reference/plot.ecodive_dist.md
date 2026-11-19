# Plot a distance matrix

A heatmap of a distance matrix from
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md) or
any of its metric functions (such as
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) or
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md)),
or of the proportionality of features from
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) or
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md): one row
and one column per sample (or feature), each cell coloured by the
distance between its row's and its column's. Nothing is recomputed: the
samples are drawn in their order, or sorted by `groups`, and not
clustered.

## Usage

``` r
# S3 method for class 'ecodive_dist'
plot(x, groups = NULL, data = NULL, col = NULL, key = TRUE, ...)
```

## Arguments

- x:

  A distance matrix of class `ecodive_dist`, as ecodive's distance
  functions return. (`plot(as.dist(m))` of another matrix would draw
  base R's plot of the distances as numbers; give it the class with
  `class(d) <- c("ecodive_dist", "dist")`.)

- groups:

  Sorts the samples by group: the name of a column of `data`, or a
  vector with one value per sample (matched by name when both are named,
  as for
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)).
  Factor, character or logical groups are drawn as blocks, in the order
  of a factor's levels, separated by white lines and labelled with the
  group; numeric ones sort the samples from smallest to largest.

- data:

  A data frame of sample metadata, with one row per sample, when
  `groups` names one of its columns.

- col:

  The colour scale, from the smallest value to the largest. Defaults to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"YlOrRd"` palette, light for near and dark for far; or, when some
  values are negative (as from
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md)), its
  `"Blue-Red 3"` palette, symmetric about 0.

- key:

  Whether to draw the colour key on the right. Default `TRUE`.

- ...:

  Passed on to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html), such as
  `main`, `sub`, or `useRaster`.

## Value

`x`, invisibly.

## Details

The rows run from the first sample at the top to the last at the bottom,
and the columns from left to right in the same order, so the matrix is
symmetric about the diagonal from the top left; only the rows are
labelled. The diagonal, which a `dist` object doesn't store, is left
blank.

The heatmap has one cell per pair of samples, so it needs memory for an
\\n \times n\\ matrix: subset or summarize thousands of samples first.
Large heatmaps are drawn as raster images, where the device can.

The class `ecodive_dist` marks distances that ecodive computed. It comes
before `"dist"`, so every function for `dist` objects works on them as
before.

## Examples

``` r
    dm <- bray(iris[, 1:4])
    plot(dm)


    # Sorted by species: three blocks along the diagonal.
    plot(dm, 'Species', iris)


    # Proportionality, with negative values.
    plot(rho_p(ex_counts + 1))
```
