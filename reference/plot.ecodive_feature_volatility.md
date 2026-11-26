# Plot a Feature Volatility Result

The importances of the most important features as bars, beside the
volatility chart of one of them (as
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)'s
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)): its relative
abundance over time, each group's mean with a bootstrap interval, each
subject's trajectory, and the global mean and control limits.

## Usage

``` r
# S3 method for class 'ecodive_feature_volatility'
plot(
  x,
  feature = 1L,
  n_bars = 20L,
  spaghetti = TRUE,
  limits = TRUE,
  legend = "topright",
  col = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`feature_volatility()`](https://cmmr.github.io/ecodive/reference/feature_volatility.md)
  result.

- feature:

  The feature to chart: its name, or its rank (`1`, the default, is the
  most important).

- n_bars:

  How many of the most important features to draw as bars. Default:
  `20`.

- spaghetti:

  Draw each subject's trajectory. Default: `TRUE`.

- limits:

  Draw the global mean and control limits. Default: `TRUE`.

- legend:

  Where to put the legend of groups (a keyword of
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html), such
  as `"topright"`), or `NULL` for none. Only drawn for two or more
  groups.

- col:

  The groups' colours. Default:
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  "Dark 3" palette.

- ...:

  Passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)
  for the chart, such as `main`, `xlab`, `ylab` or `ylim`.

## Value

`x`, invisibly.

## Examples

``` r
    meta <- data.frame(
      subject = rep(sprintf('s%d', 1:8), each = 4),
      day     = rep(c(0, 7, 14, 21), 8) )
    counts <- outer(seq_len(32), 1:6, function (i, j) (i * j) %% 11 + 5)
    counts[, 1] <- counts[, 1] + 3 * meta$day
    colnames(counts) <- sprintf('F%d', 1:6)
    fv <- feature_volatility(counts, 'subject', 'day', data = meta,
                             n_trees = 20, cpus = 1)

    plot(fv)

    plot(fv, feature = fv$features$feature[[1]], spaghetti = FALSE, limits = FALSE)
```
