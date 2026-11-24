# Plot a MaAsLin 3 Result

Views of a
[`maaslin3()`](https://cmmr.github.io/ecodive/reference/maaslin3.md)
result for one covariate and one of its two models, each feature a point
(or a row), with the features that the model declares associated
(`reject_abundance` or `reject_prevalence`) highlighted and the most
significant of them named:

## Usage

``` r
# S3 method for class 'ecodive_maaslin3'
plot(
  x,
  covariate = NULL,
  model = c("abundance", "prevalence"),
  type = c("volcano", "ma", "forest", "rank"),
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
  [`maaslin3()`](https://cmmr.github.io/ecodive/reference/maaslin3.md).

- covariate:

  Which covariate to plot, as named in the result's `covariate` column
  (such as `"statusmoderate"` or `"age"`). Defaults to the first column
  of the term of interest (`groups`), or else to the first covariate.

- model:

  `"abundance"` (the default) or `"prevalence"`.

- type:

  The plot: `"volcano"` (the default), `"ma"`, `"forest"`, or `"rank"`.

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

- `"volcano"` (the default): the coefficient against its adjusted
  p-value.

- `"ma"`: the coefficient against the feature's mean transformed
  abundance where present.

- `"forest"`: the coefficients and their confidence intervals of the
  `n_features` most significant features.

- `"rank"`: every feature's coefficient, sorted (see **Rank plots**).

Abundance coefficients are in the units of the result's `norm` and
`transform` (log2 by default, its `log_base`); prevalence coefficients
are log odds ratios.

- **Volcano plot.** The y axis is \\-\log\_{10}\\ of the model's
  adjusted p-value (`q_abundance` or `q_prevalence`; 0 is drawn at the
  smallest positive double), with a dashed line at `alpha`.

- **MA plot.** The x axis is the feature's mean transformed abundance
  where present (the `mean` attribute).

- **Forest plot.** Bars over \\\mathrm{coef} \pm q\\\mathrm{se}\\, with
  \\q\\ the \\(1 + \\`conf_level`\\)/2\\ quantile of the t distribution
  with `dof_abundance` degrees of freedom (abundance), or of the normal
  distribution (prevalence); they aren't adjusted for testing many
  features. With a median comparison, the coefficients were tested
  against the median (named in the subtitle), not against 0.

Features whose values are `NaN` are left out. The subtitle names the
covariate and the model, and counts the features rejected.

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
    x <- matrix(rpois(30 * 25, 10) * rbinom(30 * 25, 1, 0.7), 30, 25)
    x[16:30, 1:3] <- x[16:30, 1:3] * 6
    res <- maaslin3(x, groups = rep(c('a', 'b'), each = 15))
    plot(res)

    plot(res, model = 'prevalence')

    plot(res, type = 'forest', n_features = 10)

    plot(res, type = 'rank', numerator = 1:3)
```
