# Plot a fastANCOM Result

Views of a
[`fastancom()`](https://cmmr.github.io/ecodive/reference/fastancom.md)
result for one comparison, with the features declared differentially
abundant (`reject = TRUE`) highlighted:

## Usage

``` r
# S3 method for class 'ecodive_fastancom'
plot(
  x,
  covariate = NULL,
  type = c("w", "volcano", "ma", "forest", "rank"),
  labels = 10,
  n_features = 30,
  conf_level = 0.95,
  numerator = NULL,
  denominator = NULL,
  ...
)
```

## Arguments

- x:

  The result of
  [`fastancom()`](https://cmmr.github.io/ecodive/reference/fastancom.md).

- covariate:

  Which comparison to plot, as named in the result's `covariate` column.
  Defaults to the global test of the term of interest when there is one,
  else to its first column.

- type:

  The plot: `"w"` (the default), `"volcano"`, `"ma"`, `"forest"`, or
  `"rank"`. A global test has \\W\\ but no log2 fold change, so only
  `"w"`.

- labels:

  How many of the differentially abundant features to name (in the
  volcano, effect and MA plots), those with the smallest adjusted
  p-values. In the rank plot, the highlighted features with the largest
  effects. Default `10`; `0` for none.

- n_features:

  For `"w"`: how many features to show, those with the largest \\W\\;
  for `"forest"`, those with the smallest adjusted p-values. Default
  `30`; `Inf` for all.

- conf_level:

  For `"forest"`: the confidence level of the intervals. Default `0.95`.

- numerator, denominator:

  For `"rank"`: the features to highlight as the numerator (red) and
  denominator (blue) of a log ratio, such as those that
  [`log_ratio_features()`](https://cmmr.github.io/ecodive/reference/log_ratio_features.md)
  picks: their names, their indices among the result's features, or a
  logical vector with one value per feature, as for
  [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md).
  They must not overlap.

- ...:

  Passed on to
  [`graphics::dotchart()`](https://rdrr.io/r/graphics/dotchart.html)
  (`"w"`) or
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub`, `xlab` or `cex`.

## Value

`x`, invisibly.

## Details

- `"w"` (the default): a dot chart of the features' \\W\\, those with
  the largest at the top, with a dashed line at the smallest \\W\\ that
  is rejected (the `cutoff` attribute).

- `"volcano"`: the log2 fold change against the adjusted p-value.

- `"ma"`: the log2 fold change against the mean log2 abundance.

- `"forest"`: the log2 fold changes and their confidence intervals of
  the `n_features` most significant features.

- `"rank"`: every feature's log2 fold change, sorted.

Log fold changes are in log2, against the reference set (see
[`fastancom()`](https://cmmr.github.io/ecodive/reference/fastancom.md)).
In the volcano plot, the y axis is \\-\log\_{10}\\ of `p_adjusted`, with
a dashed line at `alpha`; fastANCOM's decisions come from \\W\\, so the
highlighted features needn't be those above the line. The forest plot's
bars are \\\mathrm{log2\\fc} \pm q\\\mathrm{se}\\, with \\q\\ the \\(1 +
\\`conf_level`\\)/2\\ quantile of the t distribution with the result's
`dof`, not adjusted for testing many features.

## Rank plots

A rank plot (`type = "rank"`) draws every feature as a bar from 0 to its
effect, sorted from the smallest at the left to the largest at the right
(ties in the result's order), as Qurro's rank plot does (Fedarko et al.
2020,
[doi:10.1093/nargab/lqaa023](https://doi.org/10.1093/nargab/lqaa023) ).
With relative abundances, every feature's effect is known only up to a
shift shared by all of them, the unknown change in the total (Morton et
al. 2019,
[doi:10.1038/s41467-019-10656-5](https://doi.org/10.1038/s41467-019-10656-5)
): the ranking is what the data show, and the log ratio of features from
its two ends can be read without knowing the total. Pick them with
[`log_ratio_features()`](https://cmmr.github.io/ecodive/reference/log_ratio_features.md),
show where they sit with `numerator` and `denominator`, and plot their
log ratio sample by sample with
[`plot.ecodive_log_ratio()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_log_ratio.md).

Without `numerator` or `denominator`, the differentially abundant
features (`reject = TRUE`) are drawn dark and the others light; with
them, the numerator's bars are red and the denominator's blue, with a
legend. The `labels` highlighted features with the largest effects are
named. Features with no estimate (`NaN`) are left out and counted in the
subtitle; for a model result, the bars are `covariate`'s effects.

## Examples

``` r
    set.seed(1)
    x <- matrix(rpois(20 * 40, 40), 20, 40)
    x[11:20, 1:4] <- x[11:20, 1:4] * 5
    res <- fastancom(x, groups = rep(c('a', 'b'), each = 10))
    plot(res)

    plot(res, type = 'volcano')

    plot(res, type = 'forest', n_features = 10)
```
