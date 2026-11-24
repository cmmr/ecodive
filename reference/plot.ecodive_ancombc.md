# Plot an ANCOM-BC result

Plots the results of
[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
[`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
of their post-hoc tests, from what each holds:

## Usage

``` r
# S3 method for class 'ecodive_ancombc'
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

# S3 method for class 'ecodive_ancombc_pairwise'
plot(x, n_features = 30, digits = 2L, ...)

# S3 method for class 'ecodive_ancombc_dunnett'
plot(x, n_features = 30, digits = 2L, ...)

# S3 method for class 'ecodive_ancombc_global'
plot(x, n_features = 30, ...)

# S3 method for class 'ecodive_ancombc_trend'
plot(x, n_features = 30, ...)
```

## Arguments

- x:

  The result of
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  or
  [`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md).

- covariate:

  For [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md):
  which covariate's log fold changes to plot, as named in the result's
  `covariate` column (such as `"statusmoderate"` or `"age"`). Defaults
  to the first covariate that isn't the intercept.

- type:

  For [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md):
  the plot, `"volcano"` (the default), `"effect"`, `"ma"`, `"forest"`,
  or `"rank"`.

- labels:

  For [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md):
  how many of the differentially abundant features to name (in the
  volcano, effect and MA plots), those with the smallest adjusted
  p-values. In the rank plot, the highlighted features with the largest
  effects. Default `10`; `0` for none.

- n_features:

  How many features to show: for the forest plot, those with the
  smallest adjusted p-value (ties in the table's order); for the
  heatmap, those with the smallest adjusted p-value in any comparison;
  for the dot chart, those with the largest statistic. Default `30`;
  `Inf` for all.

- conf_level:

  For the forest plot: the confidence level of the intervals. Default
  `0.95`.

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
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)
  for the volcano, effect, MA and forest plots,
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html) for the
  heatmap (where `col` is the colour scale, from most negative to most
  positive), and
  [`graphics::dotchart()`](https://rdrr.io/r/graphics/dotchart.html) for
  the dot chart (`col` as its `color`), such as `main`, `sub`, `xlab` or
  `cex`.

- digits:

  For the heatmap: the number of decimal places of the log fold changes
  written in its cells. Default `2`.

## Value

`x`, invisibly.

## Details

- [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md):
  one covariate's bias-corrected log fold changes, as one of four plots
  (`type`): `"volcano"` (the default), against their adjusted p-values;
  `"effect"`, against each feature's residual standard deviation (Gloor
  et al. 2016; ALDEx2's `"MW"` plot), for a two-level comparison;
  `"ma"`, against each feature's mean clr abundance (ALDEx2's `"MA"`
  plot); or `"forest"`, with their confidence intervals, for the
  `n_features` most significant features; or `"rank"`, all of them
  sorted by their log fold change (see **Rank plots**).

- [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md):
  a heatmap of the log fold change of each feature (rows) in each
  comparison (columns), with the significant ones in bold.

- [`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
  [`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md):
  a dot chart of each feature's test statistic, largest first.

Features declared differentially abundant (`reject = TRUE`) are
highlighted in all of them. When the fit was made with a sensitivity
analysis (`sensitivity`, see
[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)), the
robust decisions (`reject_robust`) are drawn filled (in bold in the
heatmap) and named, and the others that are significant hollow (in
italics).

Log fold changes are natural logs, as ANCOM-BC estimates them; the axis
label and the heatmap's subtitle say so from the result's `log_base`
attribute. Coefficients that couldn't be estimated (`NaN`) are left out.

- **Volcano plot.** The y axis is \\-\log\_{10}\\ of `p_adjusted`, with
  a dashed line at `alpha`.

- **Effect plot.** The x axis is the feature's residual standard
  deviation (the `resid_sd` attribute, natural logs); the dashed lines
  are \\y = \pm x\\, so features beyond them differ by more than they
  vary within the groups. Only for a covariate whose design column is 0
  or 1 (a factor's level against its reference): for a model of that
  factor alone, `resid_sd` is the pooled within-group standard
  deviation; with other covariates, it is the spread left after them.

- **MA plot.** The x axis is the feature's mean clr abundance (the
  `mean` attribute; the rclr with zeros as missing values).

- **Forest plot.** Bars over Wald intervals, \\\mathrm{log\\fc} \pm
  q\\\mathrm{se}\\, with \\q\\ the \\(1 + \\`conf_level`\\)/2\\ quantile
  of the standard normal
  ([`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md))
  or of the t distribution with the feature's `dof`
  ([`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)),
  as the p-values; they aren't adjusted for testing many features, so a
  bar can miss 0 while the feature isn't rejected. The ANCOM-BC2
  vignette's bar charts of log fold changes and standard errors show the
  same.

In the heatmap, red cells are increases and blue cells decreases, on a
scale symmetric about 0, and the value of each cell is written in it, in
bold when that comparison is significant for that feature. For
[`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md),
each feature's label gives the pattern with the largest statistic.

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
    x <- rbind(
      c(1, 10, 16, 0, 0, 12, 0, 7), c(0, 10, 45, 2, 1, 25, 1, 9),
      c(4,  8, 22, 1, 5, 14, 1, 5), c(1, 11, 22, 1, 4, 14, 1, 3),
      c(4,  7, 36, 0, 0, 19, 0, 9), c(2, 30, 26, 0, 0, 12, 1, 0),
      c(1, 36, 10, 1, 5, 15, 2, 9), c(1, 25, 11, 2, 2,  8, 0, 5),
      c(1, 36, 23, 0, 2, 11, 1, 2), c(3, 60, 33, 3, 4, 14, 1, 6),
      c(3, 57, 30, 3, 1,  3, 0, 4), c(0, 65, 23, 0, 2,  2, 0, 8),
      c(2, 46, 19, 0, 1,  2, 0, 1), c(0, 41, 26, 0, 2,  0, 0, 3),
      c(3, 46, 20, 0, 0,  3, 0, 0) )
    dimnames(x) <- list(paste0('S', 1:15), paste0('F', 1:8))
    meta <- data.frame(
      status = rep(c('mild', 'moderate', 'severe'), each = 5),
      age    = c(25, 48, 23, 35, 52, 51, 18, 29, 40, 29, 44, 39, 26, 37, 46),
      row.names = rownames(x) )
    fit <- ancombc2(x, ~ status + age, meta, groups = 'status')

    # Moderate against mild, then age.
    plot(fit)

    plot(fit, 'age')


    # Each feature's difference against its spread, its abundance, and
    # with its interval.
    plot(fit, type = 'effect')

    plot(fit, type = 'ma')

    plot(fit, 'age', type = 'forest')


    # The features ranked, with the two at each end as a log ratio.
    sel <- log_ratio_features(fit, n = 2)
    plot(fit, type = 'rank', numerator = sel$numerator, denominator = sel$denominator)

    plot(do.call(log_ratio, c(list(x), sel)), 'status', meta)
#> Warning: 1 sample ('S14') has a numerator or denominator that sums to 0, so its log ratio is NA. A positive `pseudocount` gives it a value.


    # Every pair of groups.
    plot(ancombc_pairwise(fit))


    # Any difference between the groups.
    plot(ancombc_global(fit))
```
