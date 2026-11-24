# Plot a partition of beta diversity

For
[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md):
each pair of samples as a point, its turnover (or balanced variation)
against its nestedness (or abundance gradient), with dashed lines where
the total dissimilarity is 0.25, 0.5, 0.75 and

1.  For
    [`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md):
    with draws, the density of each component's values over the draws,
    as betapart's vignette plots them; without, one stacked bar of the
    two components of all the samples.

## Usage

``` r
# S3 method for class 'ecodive_beta_partition'
plot(x, groups = NULL, col = NULL, legend = "topright", ...)

# S3 method for class 'ecodive_beta_partition_multi'
plot(x, col = NULL, legend = "topright", ...)
```

## Arguments

- x:

  The result of
  [`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md)
  or
  [`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md).

- groups:

  Optional, for
  [`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md):
  a vector or factor with each sample's group, in the input order of the
  samples. Pairs within a group are coloured by group; pairs between
  groups are grey.

- col:

  Colours: of the groups, or of the three components (for
  [`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md)).
  Defaults to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"Dark 3"` palette.

- legend:

  Where to put the legend: a position keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html)
  (`"topright"`, ...), or `NULL` for none.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `pch` or `cex`.

## Value

`x`, invisibly.

## Examples

``` r
    bp <- beta_partition(ex_counts)
    plot(bp, groups = c('a', 'a', 'b', 'b'))


    plot(beta_partition_multi(ex_counts, n_sites = 3, n_draws = 50))
```
