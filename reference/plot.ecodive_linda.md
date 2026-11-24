# Plot a LinDA Result

Views of a
[`linda()`](https://cmmr.github.io/ecodive/reference/linda.md) result
for one covariate, each feature a point (or a row), with the features
declared differentially abundant (`reject = TRUE`) highlighted and the
most significant of them named:

## Usage

``` r
# S3 method for class 'ecodive_linda'
plot(
  x,
  covariate = NULL,
  type = c("volcano", "effect", "ma", "forest", "rank"),
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
  [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md).

- covariate:

  Which covariate's log2 fold changes to plot, as named in the result's
  `covariate` column (such as `"statusmoderate"` or `"age"`). Defaults
  to the first column of the term of interest (`groups`), or else to the
  first covariate.

- type:

  The plot: `"volcano"` (the default), `"effect"`, `"ma"`, `"forest"`,
  or `"rank"`.

- labels:

  How many of the differentially abundant features to name (in the
  volcano, effect and MA plots), those with the smallest adjusted
  p-values. In the rank plot, the highlighted features with the largest
  effects. Default `10`; `0` for none.

- n_features:

  For `"forest"`: how many features to show, those with the smallest
  adjusted p-values (ties in the table's order). Default `30`; `Inf` for
  all.

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
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub`, `xlab`, `xlim` or `cex`.

## Value

`x`, invisibly.

## Details

- `"volcano"` (the default): the log2 fold change against the adjusted
  p-value.

- `"effect"`: the log2 fold change against the residual standard
  deviation (for a two-level comparison).

- `"ma"`: the log2 fold change against the mean clr abundance.

- `"forest"`: the log2 fold changes and their confidence intervals of
  the `n_features` most significant features.

- `"rank"`: every feature's log2 fold change, sorted (see **Rank
  plots**).

Log fold changes are in log2, as LinDA estimates them (the result's
`log_base`), bias-corrected.

- **Volcano plot.** The y axis is \\-\log\_{10}\\ of `p_adjusted` (an
  adjusted p-value of 0 is drawn at the smallest positive double), with
  a dashed line at `alpha`.

- **Effect plot.** The x axis is the feature's residual standard
  deviation (the `resid_sd` attribute, log2); the dashed lines are \\y =
  \pm x\\. Only for a covariate whose design column is 0 or 1 (a
  factor's level against its reference).

- **MA plot.** The x axis is the feature's mean centred log2 ratio (the
  `mean` attribute).

- **Forest plot.** Bars over \\\mathrm{log2\\fc} \pm q\\\mathrm{se}\\,
  with \\q\\ the \\(1 + \\`conf_level`\\)/2\\ quantile of the t
  distribution with `dof` degrees of freedom, as the p-values; they
  aren't adjusted for testing many features.

Features whose values are `NaN` are left out. The subtitle names the
covariate and counts the features rejected.

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
    x <- matrix(rpois(20 * 30, 40), 20, 30)
    x[11:20, 1:3] <- x[11:20, 1:3] * 5
    res <- linda(x, groups = rep(c('a', 'b'), each = 10))
    plot(res)

    plot(res, type = 'effect')

    plot(res, type = 'forest', n_features = 10)

    plot(res, type = 'rank', numerator = 1:3)
```
