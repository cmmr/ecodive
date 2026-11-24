# Plot a Dirichlet-multinomial Linear Model

Views of a
[`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md)
result for one covariate, each feature a point (or a row), with the
features declared differentially abundant (`reject = TRUE`) highlighted
and the most significant of them named:

## Usage

``` r
# S3 method for class 'ecodive_dirmult_lm'
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
  [`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md).

- covariate:

  Which covariate's estimates to plot, as named in the result's
  `covariate` column (such as `"statusmoderate"` or `"age"`). Defaults
  to the first column of the term of interest (`groups`), or else to the
  first covariate that isn't the intercept.

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

- `"volcano"` (the default): the estimate against the adjusted p-value.

- `"effect"`: the estimate against the residual standard deviation (for
  a two-level comparison).

- `"ma"`: the estimate against the mean clr abundance.

- `"forest"`: the estimates and approximate confidence intervals of the
  `n_features` most significant features.

- `"rank"`: every feature's estimate, sorted (see **Rank plots**).

Estimates are in log2 units (the result's `log_base`). The subtitle
names the covariate, counts the features rejected, and gives the scale
model.

- **Volcano plot.** The y axis is \\-\log\_{10}\\ of `p_adjusted` (an
  adjusted p-value of 0 is drawn at the smallest positive double), with
  a dashed line at `alpha`.

- **Effect plot.** The x axis is the feature's residual standard
  deviation (the `resid_sd` attribute, log2), which for a model of one
  two-level factor is the pooled within-group standard deviation; the
  dashed lines are \\y = \pm x\\. Only for a covariate whose design
  column is 0 or 1.

- **MA plot.** The x axis is the feature's mean clr value (the `mean`
  attribute, log2).

- **Forest plot.** Bars over \\\mathrm{estimate} \pm q\\\mathrm{se}\\,
  with \\q\\ the \\(1 + \\`conf_level`\\)/2\\ quantile of the t
  distribution with `dof` degrees of freedom: approximate, as the
  estimate and standard error are averages over the draws, and the
  p-values come from the draws instead. They aren't adjusted for testing
  many features.

Features whose values are `NaN` are left out.

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
    x <- matrix(rpois(16 * 20, 50), 16, 20)
    x[9:16, 1:3] <- x[9:16, 1:3] * 4
    res <- dirmult_lm(x, groups = rep(c('a', 'b'), each = 8), n_draws = 200)
    plot(res)

    plot(res, type = 'effect')

    plot(res, type = 'forest', n_features = 10)
```
