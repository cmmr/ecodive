# Distance Test Wrapper Function

Runs any of ecodive's statistical tests of a distance matrix by name, as
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
does for ordination methods: does one variable, such as treatment
groups, explain the differences between samples?

## Usage

``` r
dist_test(dm, test, groups, data = NULL, ...)
```

## Arguments

- dm:

  The distances to test: a `dist` object such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  or a symmetric matrix of distances.

- test:

  The name of a test: `"anosim"`, `"mantel"`, `"mrpp"`, `"permanova"`,
  or `"permdisp"`. Flexible matching is supported (see below).
  Programmatic access via
  [`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md).

- groups:

  The variable to test: a vector with one value per sample, matched to
  the samples in `dm` by name when both are named, else in the same
  order as `dm`. Or, when `data` is given, the name of one of its
  columns (or a vector with one value per row of `data`). Each method
  takes categorical groups, numeric values, or both: see
  [`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md).

- data:

  Optional sample metadata: a data frame with one row per sample,
  holding `groups` (and `strata`) as columns. Rows are matched to the
  samples in `dm` by row name. Default: `NULL`.

- ...:

  Further arguments for the method's function, such as `permutations`,
  `strata`, `seed`, or `cpus`, which every method takes, or a method's
  own, such as `weights` for
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) or
  `method` (the correlation) for
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md). Each
  method keeps its own defaults. For `mantel`, `groups` is passed as
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)'s
  `y`.

## Value

The method's result: a list of class `"ecodive_test"` whose first fields
are `statistic`, `p_value`, `z`, and `n_perm`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) turns it
into one row with the same columns for every method, so results can be
[`rbind()`](https://rdrr.io/r/base/cbind.html)-ed.
`dist_test(dm, 'mrpp', groups, data, ...)` is the same as
`mrpp(dm, groups, data, ...)`.

## Details

**List of Tests**

|  |  |  |
|----|----|----|
| Option / Function Name | Test Name | Groups |
| `anosim` | Analysis of Similarities (ANOSIM) | categorical |
| `mantel` | Mantel Test | numeric |
| `mrpp` | Multi-Response Permutation Procedure (MRPP) | categorical |
| `permanova` | Permutational Multivariate Analysis of Variance | categorical or numeric |
| `permdisp` | Homogeneity of Multivariate Dispersions (PERMDISP) | categorical |

Categorical groups are factors, character, or logical vectors; a method
that takes only categorical groups gives an error for numbers, and the
reverse. `mantel` correlates `dm` with the differences between samples'
values; to compare two distance matrices, call
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)
directly. To test several variables at once, call
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
with a `formula`, either directly or as
`dist_test(dm, 'permanova', data = meta, formula = ~ a + b)`.

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `test = 'PERMANOVA'` or `'perma'` selects `permanova`. Common other
names also work: `'adonis'` or `'adonis2'` for `permanova`;
`'betadisper'` or `'PERMDISP2'` for `permdisp`; and
`'Analysis of Similarities'` for `anosim`.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

**Same names in other packages**

[`vegan::anosim()`](https://vegandevs.github.io/vegan/reference/anosim.html),
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html),
and
[`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html)
share names with ecodive's functions, and whichever package was attached
last masks the other's. Write
[`ecodive::anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md)
and so on, or use `dist_test()`, which always runs ecodive's, to be
sure.

## See also

[`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md),
[`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md),
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
[`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)

## Examples

``` r
    # 30 samples in three treatment groups. The 'high' group has more
    # of the first five genera.
    set.seed(1)
    counts <- matrix(rpois(30 * 20, 5), 30, 20)
    rownames(counts) <- paste0('S', 1:30)
    meta <- data.frame(
      row.names = rownames(counts),
      treatment = rep(c('ctrl', 'low', 'high'), each = 10),
      ph        = round(rnorm(30, 7, 0.3), 1) )
    high <- meta$treatment == 'high'
    counts[high, 1:5] <- counts[high, 1:5] + 5

    dm <- bray(counts)

    dist_test(dm, 'permanova', 'treatment', meta)
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F = 3.442, p = 0.001, z = 8.275
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.2032. Model tested as a whole:
#> 
#>          Df SumOfSqs     R2     F Pr(>F)    
#> Model     2   0.1813 0.2032 3.442  0.001 ***
#> Residual 27   0.7111 0.7968                 
#> Total    29   0.8924 1.0000                 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

    # Every test of treatment, as one table.
    tests <- list_dist_tests(val = 'id', groups = 'categorical')
    do.call(rbind, lapply(tests, function (m)
      as.data.frame(dist_test(dm, m, 'treatment', meta, cpus = 1)) ))
#>        test stat_name  statistic p_value         z n_perm
#> 1    anosim         R  0.3216543   0.001  6.741598    999
#> 2      mrpp     delta  0.2243903   0.001 -7.321666    999
#> 3 permanova         F  3.4422193   0.001  8.274905    999
#> 4  permdisp         F 15.6818804   0.001 13.423920    999

    # A numeric variable.
    dist_test(dm, 'adonis2', 'ph', meta)
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F = 0.3849, p = 0.972, z = -1.488
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.01356. Model tested as a whole:
#> 
#>          Df SumOfSqs     R2     F Pr(>F)
#> Model     1   0.0121 0.0136 0.385  0.972
#> Residual 28   0.8803 0.9864             
#> Total    29   0.8924 1.0000             
    dist_test(dm, 'mantel', 'ph', meta, method = 'spearman')
#> Mantel test (Spearman)
#> 
#> r = -0.0818, p = 0.863, z = -1.089
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30
```
