# Plot Longitudinal Results

Charts of
[`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md),
[`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md)
and
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)
results.

## Usage

``` r
# S3 method for class 'ecodive_first_differences'
plot(x, feature = NULL, ...)

# S3 method for class 'ecodive_first_distances'
plot(x, ...)

# S3 method for class 'ecodive_volatility'
plot(x, spaghetti = TRUE, limits = TRUE, legend = "topright", col = NULL, ...)
```

## Arguments

- x:

  A result of
  [`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md),
  [`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md)
  or
  [`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md).

- feature:

  For a
  [`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md)
  result of a count table: the feature to draw (a column name). Default:
  the first.

- ...:

  Passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `ylab` or `ylim`.

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

## Value

`x`, invisibly.

## Details

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
[`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md)
or
[`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md)
result draws the values against the (later) state: a box of them at each
state, and each subject's values joined by a faint line; a dashed line
at 0 for differences.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
[`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)
result draws q2-longitudinal's control chart: the state on the x axis;
each group's means, joined by a line, with their confidence intervals as
error bars, one colour per group; each subject's trajectory as a thin,
faint line in its group's colour (the "spaghetti"); the global mean as a
dashed line and the \\\pm\\ 2 and 3 SD control limits as dotted lines.

## Examples

``` r
    meta <- data.frame(
      subject = rep(c('a', 'b', 'c', 'd'), each = 3),
      day     = rep(c(0, 7, 14), 4),
      diet    = rep(c('x', 'y'), each = 6),
      shannon = c(2.1, 2.4, 2.2, 1.8, 2.0, 2.3, 2.5, 2.2, 2.6, 1.9, 2.1, 2.0) )

    plot(volatility('shannon', 'subject', 'day', 'diet', data = meta))

    plot(first_differences('shannon', 'subject', 'day', meta))
```
