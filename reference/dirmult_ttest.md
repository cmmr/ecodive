# Dirichlet-multinomial Welch t-tests (ALDEx2-style)

Finds features whose abundance differs between two sets of samples,
carrying the uncertainty of low counts and zeros through the test
(Fernandes et al. 2013, 2014). Each sample's proportions are drawn many
times from their Dirichlet posterior given its counts; every draw is
clr-transformed and each feature is compared with Welch's t-test; and
the results are averaged over the draws. This is the ALDEx2 approach, as
scikit-bio implements it, optionally with ALDEx2's and ALDEx3's scale
model, which allows for uncertainty in the samples' total loads (Nixon
et al. 2025).

## Usage

``` r
dirmult_ttest(
  counts,
  groups,
  data = NULL,
  contrast = NULL,
  margin = 1L,
  pseudocount = 0.5,
  scale_sd = 0,
  n_draws = 128L,
  alpha = 0.05,
  conf_level = 0.95,
  p_adjust = "holm",
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Raw
  counts are recommended: the Dirichlet posterior's spread depends on
  the counts' size, so proportions or rescaled values give wrong
  uncertainties.

- groups:

  Each sample's group: a vector or factor with one value per sample,
  matched to the samples in `counts` by name when both are named, else
  in the same order. Or, when `data` is given, the name of one of its
  columns (or a vector with one value per row of `data`). Every distinct
  value is a group, numbers included; there must be at least two. A
  factor's unused levels are dropped.

- data:

  Sample metadata: a data frame with one row per sample, holding the
  columns named in `groups` or `strata`. Rows are matched to the samples
  in `counts` by row name, and rows for other samples are ignored;
  without row names, they must be in the same order as the samples.
  Default: `NULL`.

- contrast:

  Which samples to compare, as `c(treatment, reference)`, two values of
  `groups`. A positive `log2_fc` means more abundant in the treatment
  group. Use `NA` as the reference to compare the treatment group with
  all other samples pooled; with two named groups, samples in other
  groups are left out. With `NULL` (the default), the treatment is the
  **first** group (in sorted order, or a factor's first level) and the
  reference is all other samples: with two groups, the first is compared
  with the second, the reverse of R's usual "second level versus first".

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  The Dirichlet prior: a positive number added to every value, zero or
  not, before drawing. Default: `0.5`, as scikit-bio and ALDEx2. This is
  not the zero patch that other functions' `pseudocount` is (no `NULL`
  rule): every value gets it.

- scale_sd:

  The scale model's uncertainty: the standard deviation, in log2 units,
  of a random difference in total load between the two sides, drawn
  afresh for each draw (see **Details**). `0` (the default) assumes the
  clr's reference (the geometric mean) is the same on both sides, as
  scikit-bio's `dirmult_ttest()` and ALDEx2's default do. Nixon et
  al. (2025) and ALDEx3 suggest `0.5`; use it with more draws
  (`n_draws = 1000` or more) and judge by `p_adjusted` (see
  **Details**).

- n_draws:

  The number of Monte Carlo draws from the Dirichlet posterior
  (scikit-bio's `draws`). Default: `128`.

- alpha:

  The significance level for `reject`, applied to the adjusted p-values.
  Strictly between 0 and 1. Default: `0.05`, as scikit-bio (which fixes
  it).

- conf_level:

  The confidence level of each draw's interval for the difference.
  Strictly between 0 and 1. Default: `0.95`, as scikit-bio (which fixes
  it).

- p_adjust:

  How to adjust the averaged p-values for testing every feature: one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html)
  (`'holm'`, the default, `'hochberg'`, `'hommel'`, `'bonferroni'`,
  `'BH'`, `'BY'`, `'fdr'`, or `'none'`).

- seed:

  An integer seed for the random draws. The same `seed` gives the same
  results, for any number of `cpus`.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A data frame of class `c('ecodive_dirmult_ttest', 'data.frame')`, one
row per feature, in the order of `counts`, with columns:

- `feature`: the feature's name (its index, as text, when `counts` has
  no feature names).

- `t`: Welch's t statistic, averaged over the draws (scikit-bio's
  `T-statistic`).

- `log2_fc`: the difference between the groups' mean clr values,
  treatment minus reference, in log2 units, averaged over the draws
  (`Log2(FC)`).

- `ci_lower`, `ci_upper`: the lowest lower bound and the highest upper
  bound of the draws' intervals for that difference, in log2 units
  (`CI(2.5)`, `CI(97.5)` for the default `conf_level`).

- `p_value`: Welch's two-sided p-value, averaged over the draws
  (`pvalue`).

- `p_adjusted`: `p_value` adjusted across the features with `p_adjust`
  (`qvalue`).

- `reject`: `TRUE` when `p_adjusted <= alpha` and the interval
  (`ci_lower`, `ci_upper`) excludes 0 (`Signif`).

Attributes: `contrast` (the treatment and reference groups; the
reference is `NA` when all other samples were pooled), `n` (the number
of samples on each side), `n_draws`, `pseudocount`, `scale_sd`, `alpha`,
`conf_level`, `p_adjust` and `seed`.

## Details

With \\x\_{sf}\\ the abundance of feature \\f\\ in sample \\s\\ and
\\c\\ the `pseudocount`, each of the `n_draws` draws:

1.  Draws \\G\_{sf} \sim \mathrm{Gamma}(x\_{sf} + c, 1)\\,
    independently. Each sample's \\G\_{s\cdot}\\, divided by its sum, is
    a draw from the Dirichlet posterior of its proportions; the clr
    makes the division unnecessary.

2.  Takes the clr: \\z\_{sf} = \log G\_{sf} - \frac{1}{D}\sum_g \log
    G\_{sg}\\.

3.  For each feature, compares the treatment's \\z\\ with the
    reference's by Welch's t-test: the difference of means \\\Delta\\,
    its standard error \\se = \sqrt{s_1^2/n_1 + s_2^2/n_2}\\, the
    Welch-Satterthwaite degrees of freedom \\\nu\\, \\t = \Delta / se\\,
    the two-sided p-value, and the interval \\\Delta \pm t\_{\nu}\\
    se\\, with \\t\_\nu\\ the \\(1 + \\`conf_level`\\)/2\\ quantile of
    Student's t distribution. With `scale_sd` \\\> 0\\, a random
    \\\Lambda \log 2\\ is first added to \\\Delta\\ (see **Scale
    uncertainty** below).

`t`, `log2_fc` and `p_value` are the averages over the draws (`log2_fc`
divided by \\\log 2\\). The interval is **not** a quantile of the draws:
it runs from the lowest lower bound to the highest upper bound of any
draw, as scikit-bio's, so a feature is only rejected when no draw's
interval includes 0. This conservative envelope guards against false
positives among rare features, whose draws vary most.

**Scale uncertainty.** A clr difference is a log fold change of the
features' absolute abundances only if the groups' geometric means are
the same, i.e. if the typical feature doesn't change; when the total
load differs between the groups, every difference is off by the same
unknown amount, and treating it as known gives false positives (Nixon et
al. 2025). The scale model makes it uncertain instead: in each draw, the
log2 difference in scale between the treatment and the reference is
\\\Lambda \sim N(0,\\ `scale_sd`\\^2)\\, shared by every sample on a
side and every feature, and \\\Lambda \log 2\\ is added to every
\\\Delta\\ (the variances, within sides, don't change). The draws'
intervals then spread, so a feature is rejected only when its change is
clear of plausible changes in load. This is ALDEx3's `clr.sm()`
(ALDEx2's `aldex.clr(gamma =)`) for the design "treatment versus
reference": ALDEx3 also draws a shift for the intercept, which cancels
from every difference. `scale_sd = 0` assumes no change in load (the clr
as it is); larger values express less confidence in that. `log2_fc` has
the same expected value either way: the shifts average out of `log2_fc`
and `t` as the draws increase, but not out of the p-values or the
intervals. Because the interval is the envelope of the draws' intervals,
with `scale_sd > 0` it keeps widening as `n_draws` grows (it reaches the
most extreme shifts drawn), so `reject` becomes stricter with more
draws; the averaged p-values converge instead, though slowly, being
driven by the rare large shifts (ALDEx3 uses 2,000 draws). So with a
scale model, use a larger `n_draws` and, for decisions that don't depend
on it, filter on `p_adjusted` alone (as ALDEx2 and ALDEx3 do). The shift
is drawn for the group that comes first (in sorted order, or a factor's
level order) against the other, so reversing a `contrast` of two groups
gives exactly the mirror image, with or without a scale model.
[`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md) tests log
ratios without the assumption, but gives no effect sizes.

A side with a single sample has no variance: its features' `t`,
`p_value`, intervals and `p_adjusted` are `NaN`, with a warning, and
none are rejected (as scikit-bio's, which gives them silently).

**Matching scikit-bio.** With `scale_sd = 0` and the same Gamma draws,
the results are those of scikit-bio's `dirmult_ttest()` (0.7.4) to about
`1e-13`: the tests feed both the same draws. The draws themselves come
from ecodive's own Gamma sampler (Marsaglia & Tsang 2000, on PCG), not
NumPy's, so for the same data the results differ from scikit-bio's by
Monte Carlo error only. The t distribution is R's
[`stats::pt()`](https://rdrr.io/r/stats/TDist.html) and
[`stats::qt()`](https://rdrr.io/r/stats/TDist.html) (SciPy's `stdtr()`
and `stdtrit()` in scikit-bio), and `p_adjust` is
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), with `NaN`
p-values left out of the family, as scikit-bio. Differences from
scikit-bio: the optional scale model (scikit-bio has none); `alpha` and
`conf_level` are arguments (scikit-bio fixes them at 0.05 and 0.95);
`pseudocount` must be positive; and `seed` defaults to `0`, not to a
random seed.

**ALDEx2** (Bioconductor) follows the same idea, but differs in its
details: it reports both Welch's and Wilcoxon's tests, averages
Benjamini-Hochberg adjusted p-values over the draws, and gives its own
effect size; its scale model (`gamma`) is off by default. ALDEx3 fits
linear models to each draw instead of t-tests. Their results agree with
this function's in direction and rough size, not exactly.

**Computation.** The draws run in parallel C, one draw at a time, with
memory for one draw (8 bytes per value of the samples in the contrast).
Every sample of every draw has its own random stream, so results are
identical for any number of `cpus`. 500 samples \\\times\\ 5,000
features, 128 draws, take a few seconds on 4 cores.

## Same names in other packages

scikit-bio (Python) has `dirmult_ttest()`, which this follows. ALDEx2's
`aldex()` with `test = "t"` is the original method (see above).

## Input Types

The `counts` parameter is designed to accept a simple numeric matrix,
but seamlessly supports objects from the following biological data
packages:

- `phyloseq`

- `rbiom`

- `SummarizedExperiment`

- `TreeSummarizedExperiment`

For large datasets, standard matrix operations may be slow. See
[`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
for details on using optimized formats (e.g. sparse matrices) and
parallel processing.

## References

Fernandes, A. D., Macklaim, J. M., Linn, T. G., Reid, G., & Gloor, G. B.
(2013). ANOVA-like differential expression (ALDEx) analysis for mixed
population RNA-Seq. *PLoS ONE*, 8, e67019.

Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A.,
Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of
high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA
gene sequencing and selective growth experiments by compositional data
analysis. *Microbiome*, 2, 15.

Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating
gamma variables. *ACM Transactions on Mathematical Software*, 26,
363-372.

Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating
scale uncertainty into differential expression analyses using ALDEx2.
*Genome Biology*, 26, 139.
[doi:10.1186/s13059-025-03609-3](https://doi.org/10.1186/s13059-025-03609-3)

## See also

[`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md).

## Examples

``` r
    # scikit-bio's example: b1 and b2 are less abundant under treatment.
    x <- rbind(
      c( 20, 110, 100, 101, 100, 103, 104),
      c( 33, 110, 120, 100, 101, 100, 102),
      c( 12, 110, 100, 110, 100,  50,  90),
      c(202, 201,   9,  10,  10,  11,  11),
      c(200, 202,  10,  10,  13,  10,  10),
      c(203, 201,  14,  10,  10,  13,  12) )
    colnames(x) <- paste0('b', 1:7)
    grp <- rep(c('treatment', 'placebo'), each = 3)

    res <- dirmult_ttest(x, grp, contrast = c('treatment', 'placebo'))
    res
#> # A data.frame: 7 × 8
#> # Dirichlet-multinomial t-test
#>   feature      t log2_fc ci_lower ci_upper p_value p_adjusted reject
#>   <chr>    <dbl>   <dbl>    <dbl>    <dbl>   <dbl>      <dbl>  <lgl>
#> 1 b1      -15.4    -5      -7.82     -2.3  0.00311    0.0186    TRUE
#> 2 b2      -17      -2.52   -3.72     -1.33 0.00102    0.00717   TRUE
#> 3 b3        6.9     1.61   -1.24      4.97 0.0236     0.0651   FALSE
#> 4 b4        6.75    1.7    -0.389     4.04 0.012      0.0602   FALSE
#> 5 b5        6.26    1.52   -1.3       4.55 0.0193     0.0651   FALSE
#> 6 b6        3.57    1.15   -1.25      3.29 0.0487     0.0651   FALSE
#> 7 b7        8.32    1.53   -1.11      3.81 0.0163     0.0651   FALSE
#> # Attributes: contrast, n, n_draws, pseudocount, scale_sd, alpha, conf_level,
#> #   p_adjust, seed
    plot(res)


    # With a scale model, the intervals widen and fewer features are
    # rejected.
    dirmult_ttest(x, grp, contrast = c('treatment', 'placebo'), scale_sd = 0.5)
#> # A data.frame: 7 × 8
#> # Dirichlet-multinomial t-test
#>   feature      t log2_fc ci_lower ci_upper p_value p_adjusted reject
#>   <chr>    <dbl>   <dbl>    <dbl>    <dbl>   <dbl>      <dbl>  <lgl>
#> 1 b1      -15.4    -4.97   -8.05    -1.83  0.00331     0.0199   TRUE
#> 2 b2      -16.8    -2.49   -4.52    -0.387 0.00154     0.0108   TRUE
#> 3 b3        7.08    1.64   -1.65     4.99  0.0402      0.139   FALSE
#> 4 b4        6.86    1.73   -1.18     4.31  0.0278      0.139   FALSE
#> 5 b5        6.44    1.55   -2.41     4.85  0.0397      0.139   FALSE
#> 6 b6        3.68    1.18   -1.81     3.6   0.102       0.139   FALSE
#> 7 b7        8.37    1.56   -0.774    4.57  0.0279      0.139   FALSE
#> # Attributes: contrast, n, n_draws, pseudocount, scale_sd, alpha, conf_level,
#> #   p_adjust, seed

    # The default contrast: the first group ('placebo', in sorted order)
    # against the rest.
    dirmult_ttest(x, grp, n_draws = 32)$log2_fc
#> [1]  4.986264  2.544825 -1.607468 -1.736224 -1.538447 -1.163018 -1.485933
```
