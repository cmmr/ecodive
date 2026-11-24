# Plot an Indicator Species Analysis

The features found by
[`indicator_species()`](https://cmmr.github.io/ecodive/reference/indicator_species.md)
(`reject`): their index as bars, or their groups' mean abundances as a
heatmap.

## Usage

``` r
# S3 method for class 'ecodive_indicator_species'
plot(x, type = c("bar", "heatmap"), n_features = 30, ...)
```

## Arguments

- x:

  The result of
  [`indicator_species()`](https://cmmr.github.io/ecodive/reference/indicator_species.md).

- type:

  `'bar'` (the default): one horizontal bar per feature, its length the
  index (`stat`), grouped and coloured by `combination`, the largest
  first within each. `'heatmap'`: the features (rows) by the groups
  (columns), each tile shaded by the group's mean abundance as a
  proportion of the feature's largest group mean, with the groups of the
  feature's `combination` outlined.

- n_features:

  How many features to show, those with the largest (absolute) index.
  Default: `30`.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab` or `sub`.

## Value

`x`, invisibly.

## Examples

``` r
    set.seed(1)
    grp <- rep(c('a', 'b', 'c'), c(6, 8, 10))
    x <- matrix(rpois(24 * 6, 2), 24, 6, dimnames = list(NULL, paste0('f', 1:6)))
    x[grp != 'a', 1:2] <- 0
    x[grp == 'a', 3] <- 0
    res <- indicator_species(x, grp, p_adjust = 'none', cpus = 1)
    plot(res)

    plot(res, 'heatmap')
```
