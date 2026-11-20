# Plot a Per-feature Test Result

A volcano plot of the results of
[`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
or [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md):
each feature's effect against its adjusted p-value, with the features
declared differentially abundant (`reject = TRUE`) highlighted and the
most significant of them named.

## Usage

``` r
# S3 method for class 'ecodive_da_ttest'
plot(x, labels = 10, ...)

# S3 method for class 'ecodive_da_wilcoxon'
plot(x, labels = 10, ...)

# S3 method for class 'ecodive_da_lm'
plot(x, covariate = NULL, labels = 10, ...)
```

## Arguments

- x:

  The result of
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  or
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md).

- labels:

  How many of the differentially abundant features to name, those with
  the smallest adjusted p-values. Default `10`; `0` for none.

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

The x axis is the effect: the difference in medians
([`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
in means
([`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
or the coefficient
([`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)),
in the units of the result's `norm` (natural logs for `'clr'`, `'rclr'`
and `'log'`). The y axis is \\-\log\_{10}\\ of `p_adjusted` (an adjusted
p-value of 0 is drawn at the smallest positive double), with a dashed
line at `alpha`. Features whose effect or p-value is `NaN` are left out.

## Examples

``` r
    x   <- ex_counts + 1
    grp <- c('oral', 'oral', 'other', 'other')
    plot(da_ttest(x, grp))

    plot(da_lm(x, groups = grp))
```
