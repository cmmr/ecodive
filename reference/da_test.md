# Differential Abundance Testing

Runs any of ecodive's differential-abundance (DA) methods by name, as
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
does for ordinations, with one convention for the groups or model to
compare.

## Usage

``` r
da_test(counts, method, groups = NULL, data = NULL, formula = NULL, ...)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**.

- method:

  The name of a DA method: `"ancom"`, `"ancombc"`, `"ancombc2"`,
  `"da_lm"`, `"da_ttest"`, `"da_wilcoxon"`, `"dirmult_lm"`,
  `"dirmult_ttest"`, `"lefse"`, or `"linda"`. Flexible matching is
  supported (see below). Programmatic access via
  [`list_da_tests()`](https://cmmr.github.io/ecodive/reference/da_tests.md).

- groups:

  The grouping to compare: the name of a column of `data`, or a vector
  with one value per sample. For a model method, the term of interest:
  alone, it is the model (`formula = ~ groups`); with `formula`, it
  names one of its terms (see below).

- data:

  A data frame of sample metadata, with a row per sample, matched to the
  samples of `counts` by name when both have names, else by position.

- formula:

  For a model method (`design = "formula"` in
  [`list_da_tests()`](https://cmmr.github.io/ecodive/reference/da_tests.md)),
  the model: a one-sided formula over the columns of `data`, such as
  `~ status + age`. A group method takes a formula of one variable,
  `~ status`, as `groups = "status"`.

- ...:

  Further arguments for the method's function, such as `margin`,
  `alpha`, `p_adjust`, `cpus`, `contrast`
  ([`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)),
  `norm` and `pseudocount`. Each method keeps its own defaults; see
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md),
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
  [`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md), and
  [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md). `seed`
  is ignored by the methods that have no random parts (all but
  [`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md),
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  and [`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md)).

## Value

The method's result, unchanged: `da_test(counts, 'ancombc2', ...)` is
identical to `ancombc2(counts, ...)`, attributes included. Every result
is a data frame with a row per feature (and covariate, for the model
methods) and a `method` attribute, the method's id.

## Details

**List of DA Methods**

|  |  |  |
|----|----|----|
| Option / Function Name | Method Name | Design |
| `ancom` | Analysis of Composition of Microbiomes (ANCOM) | groups |
| `ancombc` | ANCOM with bias correction (ANCOM-BC) | formula |
| `ancombc2` | ANCOM-BC2 | formula |
| `da_lm` | Linear model per feature | formula |
| `da_ttest` | Welch t-test per feature | groups |
| `da_wilcoxon` | Wilcoxon rank-sum test per feature | groups |
| `dirmult_lm` | Dirichlet-multinomial linear model with a scale model (ALDEx3) | formula |
| `dirmult_ttest` | Dirichlet-multinomial t-test (ALDEx2-style) | groups |
| `lefse` | LDA effect size (LEfSe) | groups |
| `linda` | Linear models for differential abundance (LinDA) | formula |

**Groups and models.** *Group methods* compare the samples of one
categorical variable: `groups`, a column name of `data` or one value per
sample. *Model methods* fit a linear model, `formula`, of any number of
variables; `groups` names its term of interest (the one their plots and
post-hoc tests use by default), and `groups` alone is the same as
`formula = ~ groups`. So `da_test(counts, m, 'status', meta)` compares
the samples by `status` for every method `m`, and covariates go in
`formula`: `da_test(counts, 'da_lm', 'status', meta, ~ status + age)`. A
one-variable `formula` works for the group methods too; more than one
variable is an error for them.

The methods' assumptions differ, and so do their results: see each
method's page. In benchmarks on real data (Nearing et al. 2022; Wirbel
et al. 2024; Pelto et al. 2025), the plain per-feature tests
([`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md))
controlled false discoveries well, and the methods often disagree: a
feature found by several methods is a firmer finding than one found by
one.

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters
(digits count only to tell `'ANCOM-BC2'` from `'ANCOM-BC'`). E.g.
`method = 'Wilcoxon'`, `'Mann-Whitney'` or `'rank sum'` selects
`da_wilcoxon`; `'t-test'` or `'Welch'` selects `da_ttest`; `'lm'` or
`'linear model'` selects `da_lm`; and `'ALDEx2'` selects
`dirmult_ttest`, an ALDEx2-style method (not ALDEx2 itself, whose
results agree with it in direction and rough size only); `'ALDEx3'`
selects `dirmult_lm`, a port of ALDEx3's linear models.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

**Same names in other packages**

None known.

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

Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
differential abundance methods produce different results across 38
datasets. *Nature Communications*, 13, 342.
[doi:10.1038/s41467-022-28034-z](https://doi.org/10.1038/s41467-022-28034-z)

Wirbel, J., Essex, M., Forslund, S. K., & Zeller, G. (2024). A realistic
benchmark for differential abundance testing and confounder adjustment
in human microbiome studies. *Genome Biology*, 25, 247.
[doi:10.1186/s13059-024-03390-9](https://doi.org/10.1186/s13059-024-03390-9)

Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary
methods provide more replicable results in microbial differential
abundance analysis. *Briefings in Bioinformatics*, 26, bbaf130.
[doi:10.1093/bib/bbaf130](https://doi.org/10.1093/bib/bbaf130)

## See also

[`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
[`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
[`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md),
[`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
[`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md),
[`linda()`](https://cmmr.github.io/ecodive/reference/linda.md),
[`list_da_tests()`](https://cmmr.github.io/ecodive/reference/da_tests.md);
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
for the zeros;
[`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md)
to rerun a method over pseudocounts, zero replacements and scale
uncertainties.

## Examples

``` r
    grp <- c('oral', 'oral', 'other', 'other')
    da_test(ex_counts, 'da_wilcoxon', grp, pseudocount = 1)
#> # A data.frame: 6 × 8
#> # Wilcoxon Rank-Sum Test
#>   feature               w   auc median_diff   mean p_value p_adjusted reject
#>   <chr>             <dbl> <dbl>       <dbl>  <dbl>   <dbl>      <dbl>  <lgl>
#> 1 Streptococcus         4  1           4.43  1.48    0.333          1  FALSE
#> 2 Bacteroides           2  0.5        -1.95  0.141   1              1  FALSE
#> 3 Corynebacterium       0  0          -2.99 -0.689   0.333          1  FALSE
#> 4 Haemophilus           4  1           4.4   0.452   0.333          1  FALSE
#> 5 Propionibacterium     1  0.25       -1.61 -0.686   0.667          1  FALSE
#> 6 Staphylococcus        0  0          -2.28 -0.702   0.333          1  FALSE
#> # Attributes: method, contrast, n, norm, pseudocount, alpha, p_adjust,
#> #   log_base

    # The same comparison by every method, with a model method's covariate.
    meta <- data.frame(site = grp, depth = rowSums(ex_counts))
    da_test(ex_counts + 1, 'da_lm', 'site', meta, ~ site + depth)
#> # A data.frame: 18 × 10
#> # Linear Model
#>    feature           covariate   estimate       se       t  ci_lower ci_upper
#>    <chr>             <chr>          <dbl>    <dbl>   <dbl>     <dbl>    <dbl>
#>  1 Streptococcus     (Intercept)  2.66    0.582      4.58   -4.73     10.1
#>  2 Streptococcus     siteother   -4.76    0.427    -11.2   -10.2       0.659
#>  3 Streptococcus     depth        0.00168 0.000831   2.03   -0.00887   0.0122
#>  4 Bacteroides       (Intercept)  3.14    6.47       0.485 -79        85.3
#>  5 Bacteroides       siteother    3.22    4.74       0.679 -57        63.4
#>  6 Bacteroides       depth       -0.00645 0.00923   -0.699  -0.124     0.111
#>  7 Corynebacterium   (Intercept) -3.72    2.99      -1.24  -41.7      34.3
#>  8 Corynebacterium   siteother    2.5     2.19       1.14  -25.3      30.4
#>  9 Corynebacterium   depth        0.00249 0.00427    0.583  -0.0517    0.0567
#> 10 Haemophilus       (Intercept)  4.43    0.877      5.05   -6.72     15.6
#> 11 Haemophilus       siteother   -3.83    0.643     -5.96  -12         4.34
#> 12 Haemophilus       depth       -0.00289 0.00125   -2.31   -0.0188    0.013
#> 13 Propionibacterium (Intercept) -3.03    3         -1.01  -41.1      35.1
#> 14 Propionibacterium siteother    1.12    2.2        0.51  -26.8      29
#> 15 Propionibacterium depth        0.0025  0.00428    0.583  -0.0519    0.0569
#> 16 Staphylococcus    (Intercept) -3.48    1.94      -1.8   -28.1      21.1
#> 17 Staphylococcus    siteother    1.75    1.42       1.23  -16.3      19.8
#> 18 Staphylococcus    depth        0.00267 0.00276    0.965  -0.0324    0.0377
#> # … 3 more columns: p_value <dbl>, p_adjusted <dbl>, reject <lgl>
#> # Attributes: method, formula, groups, norm, pseudocount, conf_level, alpha,
#> #   p_adjust, log_base, design, dof, mean, resid_sd

    # Which features each method finds.
    x <- rbind(
      c(12, 11, 10, 10, 10, 10, 10), c( 9, 11, 12, 10, 10, 10, 10),
      c( 1, 11, 10, 11, 10,  5,  9), c(22, 21,  9, 10, 10, 10, 10),
      c(20, 22, 10, 10, 13, 10, 10), c(23, 21, 14, 10, 10, 10, 10) )
    colnames(x) <- paste0('b', 1:7)
    meta <- data.frame(trt = rep(c('treatment', 'placebo'), each = 3))
    ids  <- c('ancom', 'da_ttest', 'da_wilcoxon', 'da_lm')
    sapply(ids, function (m) {
      res <- da_test(x, m, 'trt', meta)
      if ('covariate' %in% names(res)) res <- res[res$covariate != '(Intercept)', ]
      res$reject })
#>      ancom da_ttest da_wilcoxon da_lm
#> [1,] FALSE    FALSE       FALSE FALSE
#> [2,]  TRUE    FALSE       FALSE FALSE
#> [3,] FALSE    FALSE       FALSE FALSE
#> [4,] FALSE    FALSE       FALSE FALSE
#> [5,] FALSE    FALSE       FALSE FALSE
#> [6,] FALSE    FALSE       FALSE FALSE
#> [7,] FALSE    FALSE       FALSE FALSE
```
