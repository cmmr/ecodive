# Plot a Per-feature Test Result

Views of the results of
[`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
or [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
each feature a point (or a row), with the features declared
differentially abundant (`reject = TRUE`) highlighted and the most
significant of them named:

## Usage

``` r
# S3 method for class 'ecodive_da_ttest'
plot(
  x,
  type = c("volcano", "effect", "ma", "forest"),
  labels = 10,
  n_features = 30,
  ...
)

# S3 method for class 'ecodive_da_wilcoxon'
plot(x, type = c("volcano", "ma"), labels = 10, ...)

# S3 method for class 'ecodive_da_lm'
plot(
  x,
  covariate = NULL,
  type = c("volcano", "effect", "ma", "forest"),
  labels = 10,
  n_features = 30,
  ...
)
```

## Arguments

- x:

  The result of
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  or
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md).

- type:

  The plot: `"volcano"` (the default), `"effect"`, `"ma"`, or `"forest"`
  (for
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  `"volcano"` or `"ma"`).

- labels:

  How many of the differentially abundant features to name (in the
  volcano, effect and MA plots), those with the smallest adjusted
  p-values. Default `10`; `0` for none.

- n_features:

  For `"forest"`: how many features to show, those with the smallest
  adjusted p-values (ties in the table's order). Default `30`; `Inf` for
  all.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub`, `xlab`, `xlim` or `cex`.

- covariate:

  For
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md):
  which covariate's coefficients to plot, as named in the result's
  `covariate` column (such as `"statusmoderate"` or `"age"`). Defaults
  to the first column of the term of interest (`groups`), or else to the
  first covariate that isn't the intercept.

## Value

`x`, invisibly.

## Details

- `"volcano"` (the default): the effect against the adjusted p-value.

- `"effect"`
  ([`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)):
  the effect against the within-group standard deviation (Gloor et al.
  2016; ALDEx2's `"MW"` plot).

- `"ma"`: the effect against the mean abundance (ALDEx2's `"MA"` plot).

- `"forest"`
  ([`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)):
  the effects and their confidence intervals of the `n_features` most
  significant features, largest at the top.

The effect is the difference in medians
([`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
in means
([`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
or the coefficient
([`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
in the units of the result's `norm`: a log fold change for `'clr'`,
`'rclr'` and `'log'` (natural logs, from the `log_base` attribute), a
difference in proportions for `'percent'`.

- **Volcano plot.** The y axis is \\-\log\_{10}\\ of `p_adjusted` (an
  adjusted p-value of 0 is drawn at the smallest positive double), with
  a dashed line at `alpha`.

- **Effect plot.** The x axis is the within-group dispersion:
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)'s
  `pooled_sd`, or
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)'s
  `resid_sd` attribute (the residual standard deviation, which for a
  model of one two-level factor is the pooled within-group standard
  deviation). The dashed lines are \\y = \pm x\\: features beyond them
  differ between the groups by more than they vary within them. For
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  only a covariate whose design column is 0 or 1 (a factor's level
  against its reference) has one; use `"ma"` or `"forest"` for others.

- **MA plot.** The x axis is the feature's mean transformed abundance:
  the `mean` column
  ([`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md))
  or attribute
  ([`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)).

- **Forest plot.** Points at the effects, bars over their confidence
  intervals (`ci_lower`, `ci_upper`, at the result's `conf_level`, not
  adjusted for testing many features).

Features whose values are `NaN` are left out. The subtitle names the
comparison and counts the features rejected.

## References

Gloor, G. B., Macklaim, J. M., & Fernandes, A. D. (2016). Displaying
variation in large datasets: plotting a visual summary of effect sizes.
*Journal of Computational and Graphical Statistics*, 25, 971-979.
[doi:10.1080/10618600.2015.1131161](https://doi.org/10.1080/10618600.2015.1131161)

## Examples

``` r
    x   <- ex_counts + 1
    grp <- c('oral', 'oral', 'other', 'other')
    res <- da_ttest(x, grp)
    plot(res)

    plot(res, 'effect')

    plot(res, 'forest', n_features = 10)

    plot(da_lm(x, groups = grp), type = 'ma')
```
