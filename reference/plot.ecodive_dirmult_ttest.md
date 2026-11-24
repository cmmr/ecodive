# Plot a Dirichlet-multinomial t-test

Four views of
[`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)'s
results, each feature a point (or a row), with the features that differ
(`reject = TRUE`) highlighted and the most significant of them named:

## Usage

``` r
# S3 method for class 'ecodive_dirmult_ttest'
plot(
  x,
  type = c("volcano", "effect", "ma", "forest", "rank"),
  labels = 10,
  n_features = 30,
  numerator = NULL,
  denominator = NULL,
  ...
)
```

## Arguments

- x:

  The result of
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md).

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

- `"effect"`: the log2 fold change against the pooled within-group
  standard deviation (Gloor et al. 2016; ALDEx2's `"MW"` plot).

- `"ma"`: the log2 fold change against the mean clr abundance (ALDEx2's
  `"MA"` plot).

- `"forest"`: the log2 fold changes and their intervals of the
  `n_features` most significant features, largest at the top.

- `"rank"`: every feature's log2 fold change, sorted (see **Rank
  plots**).

**Volcano plot.** The y axis is \\-\log\_{10}\\ of `p_adjusted` (an
adjusted p-value of 0 is drawn at the smallest positive double). The
dashed line marks `alpha`; a feature also needs its interval for the
log2 fold change to exclude 0 to count as differentially abundant, so
some features above the line may not be highlighted. Gloor et al. (2016)
argue that volcano plots mislead: with many samples, tiny differences
reach small p-values, so read them with the effect plot.

**Effect plot.** Each feature's difference between the groups (y,
`log2_fc`) against its dispersion within them (x, `pooled_sd`), both in
log2 units. The dashed lines are \\y = \pm x\\: features beyond them
differ between the groups by more than they vary within them (an average
`effect` of about 1 or more). Features near the y axis vary little, so a
small difference can be significant; features far to the right need a
large difference. This is ALDEx2's `aldex.plot(type = "MW")`, with
ALDEx3's pooled standard deviation in place of ALDEx2's `diff.win` (see
[`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)'s
**Effect sizes**).

**MA plot.** The log2 fold change against the feature's mean clr value
(`mean`): how the differences depend on abundance. Rare features (on the
left) vary most, from sampling alone.

**Forest plot.** A point at each feature's `log2_fc` and a bar from
`ci_lower` to `ci_upper`, which are the envelope of the draws' intervals
(the lowest lower and highest upper bound), not an interval of one draw:
wider than a single confidence interval.

Features whose value is `NaN` are left out of each plot. The subtitle
names the contrast and counts the features rejected.

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

## References

Gloor, G. B., Macklaim, J. M., & Fernandes, A. D. (2016). Displaying
variation in large datasets: plotting a visual summary of effect sizes.
*Journal of Computational and Graphical Statistics*, 25, 971-979.
[doi:10.1080/10618600.2015.1131161](https://doi.org/10.1080/10618600.2015.1131161)

## Examples

``` r
    res <- dirmult_ttest(round(iris[, 1:4] * 10), 'Species', iris)
    plot(res)

    plot(res, 'effect')

    plot(res, 'ma')

    plot(res, 'forest')

    plot(res, 'rank', numerator = 'Petal.Width', denominator = 'Sepal.Width')
```
