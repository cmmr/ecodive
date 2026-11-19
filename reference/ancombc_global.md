# Post-hoc tests for ANCOM-BC and ANCOM-BC2

Multigroup tests on a fit from
[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) or
[`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md), for
a categorical variable with three or more levels (Lin & Peddada 2024),
as scikit-bio's `ANCOMBCResult` methods:

## Usage

``` r
ancombc_global(
  fit,
  groups = attr(fit, "groups"),
  alpha = attr(fit, "alpha"),
  p_adjust = attr(fit, "p_adjust"),
  cpus = n_cpus()
)

ancombc_pairwise(
  fit,
  groups = attr(fit, "groups"),
  alpha = attr(fit, "alpha"),
  p_adjust = attr(fit, "p_adjust"),
  cpus = n_cpus()
)

ancombc_dunnett(
  fit,
  groups = attr(fit, "groups"),
  alpha = attr(fit, "alpha"),
  p_adjust = attr(fit, "p_adjust"),
  n_bootstraps = 100L,
  seed = 0L,
  cpus = n_cpus()
)

ancombc_trend(
  fit,
  groups = attr(fit, "groups"),
  alpha = attr(fit, "alpha"),
  p_adjust = attr(fit, "p_adjust"),
  patterns = NULL,
  nodes = NULL,
  n_bootstraps = 100L,
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- fit:

  A result of
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) or
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md).

- groups:

  The categorical variable to compare: the name of a main-effect term of
  the fit's formula with at least three levels. Default: the fit's own
  `groups`. (Any such term can be given, since the fit keeps every
  coefficient's covariances.)

- alpha:

  The significance level of the adjusted p-values. Default: the fit's.

- p_adjust:

  The adjustment of the p-values (any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)): across
  the features for `ancombc_global()` and `ancombc_trend()`; within each
  feature, as the family-wise error control of mdFDR, for
  `ancombc_pairwise()` and `ancombc_dunnett()`. Default: the fit's.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

- n_bootstraps:

  The number of bootstrap draws of the null statistics, per feature.
  Default: `100`, as scikit-bio and ANCOMBC. Each p-value of the trend
  test is a proportion of these draws, so it needs many more (about
  10,000) to be stable.

- seed:

  An integer seed for the random draws. The same `seed` gives the same
  results, for any number of `cpus`.

- patterns:

  `ancombc_trend()`: the ordered patterns to test, as a named list of
  constraint matrices, each with one column per coefficient of `groups`
  (in the design's order) and one row per constraint `C %*% beta >= 0`.
  Default: `increasing` and `decreasing`, the coefficients (against the
  reference level, 0) rising, or falling, from level to level.

- nodes:

  `ancombc_trend()`: for each pattern, the coefficient (an index into
  the `groups` coefficients) whose contrast with the reference and with
  the last level gives the test statistic. Default: the last, for every
  pattern.

## Value

A data frame of class `c('ecodive_ancombc_<test>', 'data.frame')`.

`ancombc_global()` and `ancombc_trend()`: one row per feature, with
`feature`, (`ancombc_trend()`: `pattern`, the pattern with the largest
statistic,) `w` (the statistic), `p_value`, `p_adjusted` and `reject`.

`ancombc_pairwise()` and `ancombc_dunnett()`: one row per feature and
comparison, with `feature`, `comparison` (`'b - a'`: level b against
level a; with no intercept, a level alone is its own coefficient),
`log_fc`, `se`, `w`, `p_value`, `p_adjusted` and `reject`. As in
scikit-bio and ANCOMBC, `p_value` is set to 1 for the features that the
screening global test didn't find (and wherever it would be 0 or `NaN`),
and `p_adjusted` is the mdFDR-adjusted value.

Attributes: `method` (the fit's), `groups`, `alpha`, `p_adjust`, and for
the bootstrapped tests `n_bootstraps` and `seed`.

## Details

- `ancombc_global()`: does the feature differ between any of the levels?
  A Wald-type test of all the variable's coefficients at once.

- `ancombc_pairwise()`: every level against the reference, and every
  pair of the other levels, with the mixed directional false discovery
  rate (mdFDR) screened by the global test.

- `ancombc_dunnett()`: every level against the reference only, with
  mdFDR screened by a bootstrapped global test (Dunnett's type).

- `ancombc_trend()`: does the feature follow an ordered pattern across
  the levels (by default increasing or decreasing)? Constrained
  estimates, with a bootstrapped p-value.

With \\\hat\beta_f\\ the feature's \\k\\ coefficients of `groups` and
\\V_f\\ their covariance matrix (from the fit):

- **Global test:** \\W = \hat\beta' V^+ \hat\beta\\, referred to
  \\\chi^2_k\\
  ([`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md))
  or \\F\_{k, \nu}\\
  ([`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  with the feature's residual degrees of freedom), with the two-sided
  p-value \\2 \min(P, 1 - P)\\, as in ANCOMBC and scikit-bio.

- **Pairwise and Dunnett's tests:** each comparison's \\W = \hat\delta /
  \mathrm{se}\\, against the normal or t distribution as in the fit.
  mdFDR (Guo, Sarkar & Peddada 2010): a global screen at `alpha`, with
  Benjamini-Hochberg across features (the global test for pairwise; for
  Dunnett's, the largest \\\|W\|\\ against `n_bootstraps` draws of the
  largest of \\k\\ null t or normal values); then each screened
  feature's comparisons are adjusted with `p_adjust` for a family of \\m
  D / R\\ tests (\\m\\ comparisons, \\D\\ features, \\R\\ passing the
  screen).

- **Trend test:** for each pattern, the estimate closest to
  \\\hat\beta\\ in the metric \\V^{-1}\\ that satisfies its constraints
  (a small quadratic programme); its statistic is
  \\\max(\|\beta\_{node}\|, \|\beta\_{node} - \beta_k\|)\\, and \\W\\
  the largest over the patterns. The p-value is the proportion of
  `n_bootstraps` null draws (standard normal coefficients, projected the
  same way and scaled by the standard errors) whose statistic exceeds
  \\W\\; then `p_adjust` across features.

Features whose coefficients can't be estimated get `NaN` statistics and
a p-value of 1.

**Random numbers.** The bootstraps use ecodive's own generator (PCG),
with `seed`: one stream per feature, so results are identical for any
`cpus` and a feature's results don't depend on the others.
`ancombc_global()` and `ancombc_pairwise()` draw nothing.

## Differences from scikit-bio and ANCOMBC

The statistics and deterministic p-values are scikit-bio's (0.7.4), to
about `1e-10` given the same fit; the bootstrapped p-values agree with
scikit-bio's in distribution (they agree exactly given the same null
draws), not draw for draw, since the random numbers differ. The trend
test's constrained estimates are the exact solutions of its quadratic
programmes, which scikit-bio's SLSQP and ANCOMBC's `quadprog`
approximate. Comparisons are labelled with level names (scikit-bio:
`status[T.severe]_status[T.moderate]`), and `groups` can differ from the
fit's.

## References

Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
microbiomes with covariate adjustments and repeated measures. *Nature
Methods*, 21, 83–91.

Guo, W., Sarkar, S. K., & Peddada, S. D. (2010). Controlling false
discoveries in multidimensional directional decisions, with applications
to gene expression data on ordered categories. *Biometrics*, 66,
485–492.

## See also

[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
[`structural_zeros()`](https://cmmr.github.io/ecodive/reference/structural_zeros.md).

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

    ancombc_global(fit)
#>   feature           w      p_value   p_adjusted reject
#> 1      F1  1.78597500 4.568396e-01 1.000000e+00  FALSE
#> 2      F2 60.28742380 2.361777e-06 1.889422e-05   TRUE
#> 3      F3  0.03998274 7.811002e-02 4.686601e-01  FALSE
#> 4      F4  4.78034917 2.334486e-01 1.000000e+00  FALSE
#> 5      F5  1.38517980 6.225672e-01 1.000000e+00  FALSE
#> 6      F6 25.90625478 2.216400e-04 1.551480e-03   TRUE
#> 7      F7  0.86510213 9.745627e-01 1.000000e+00  FALSE
#> 8      F8  1.58903521 5.128828e-01 1.000000e+00  FALSE
    res <- ancombc_pairwise(fit)
    res[res$reject, ]
#>    feature        comparison    log_fc        se         w      p_value
#> 4       F2   moderate - mild  1.578050 0.2867837  5.502581 1.855143e-04
#> 5       F2     severe - mild  1.880123 0.2997796  6.271683 6.073015e-05
#> 17      F6     severe - mild -1.728554 0.3428363 -5.041923 5.051992e-04
#> 18      F6 severe - moderate -1.653628 0.4016451 -4.117137 2.086614e-03
#>      p_adjusted reject
#> 4  0.0020406570   TRUE
#> 5  0.0007287618   TRUE
#> 17 0.0060623903   TRUE
#> 18 0.0229527490   TRUE
    ancombc_dunnett(fit)
#>    feature      comparison       log_fc        se           w      p_value
#> 1       F1 moderate - mild -0.130278285 0.4277787 -0.30454595 1.000000e+00
#> 2       F1   severe - mild  0.449835345 0.4455405  1.00963962 1.000000e+00
#> 3       F2 moderate - mild  1.578050339 0.2867837  5.50258086 1.855143e-04
#> 4       F2   severe - mild  1.880122789 0.2997796  6.27168297 6.073015e-05
#> 5       F3 moderate - mild -0.054194121 0.2713215 -0.19974133 1.000000e+00
#> 6       F3   severe - mild -0.004017009 0.2968972 -0.01352997 1.000000e+00
#> 7       F4 moderate - mild  0.932813583 0.4266513  2.18636071 1.000000e+00
#> 8       F4   severe - mild  0.600334750 0.5156519  1.16422490 1.000000e+00
#> 9       F5 moderate - mild  0.069268610 0.3884956  0.17829960 1.000000e+00
#> 10      F5   severe - mild -0.448051616 0.4668154 -0.95980468 1.000000e+00
#> 11      F6 moderate - mild -0.074926233 0.2834519 -0.26433494 7.968898e-01
#> 12      F6   severe - mild -1.728554271 0.3428363 -5.04192317 5.051992e-04
#> 13      F7 moderate - mild  0.347893536 0.4396022  0.79138265 1.000000e+00
#> 14      F7   severe - mild  0.000000000 0.5163252  0.00000000 1.000000e+00
#> 15      F8 moderate - mild  0.028229795 0.4110188  0.06868249 1.000000e+00
#> 16      F8   severe - mild -0.530294725 0.4443275 -1.19347709 1.000000e+00
#>      p_adjusted reject
#> 1  1.0000000000  FALSE
#> 2  1.0000000000  FALSE
#> 3  0.0012985999   TRUE
#> 4  0.0004858412   TRUE
#> 5  1.0000000000  FALSE
#> 6  1.0000000000  FALSE
#> 7  1.0000000000  FALSE
#> 8  1.0000000000  FALSE
#> 9  1.0000000000  FALSE
#> 10 1.0000000000  FALSE
#> 11 1.0000000000  FALSE
#> 12 0.0040415935   TRUE
#> 13 1.0000000000  FALSE
#> 14 1.0000000000  FALSE
#> 15 1.0000000000  FALSE
#> 16 1.0000000000  FALSE
    ancombc_trend(fit, n_bootstraps = 1000)
#>   feature    pattern          w p_value p_adjusted reject
#> 1      F1 increasing 0.51303279   0.278      1.000  FALSE
#> 2      F2 increasing 1.88012279   0.000      0.000   TRUE
#> 3      F3 decreasing 0.03140829   0.988      1.000  FALSE
#> 4      F4 increasing 0.83253144   0.119      0.714  FALSE
#> 5      F5 decreasing 0.48538297   0.355      1.000  FALSE
#> 6      F6 decreasing 1.72855427   0.000      0.000   TRUE
#> 7      F7 increasing 0.23158725   0.789      1.000  FALSE
#> 8      F8 decreasing 0.53852924   0.228      1.000  FALSE
```
