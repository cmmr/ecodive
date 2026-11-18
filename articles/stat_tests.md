# Statistical Tests

## Introduction

An ordination plot can suggest that groups of samples differ. ecodive’s
permutation tests check whether the difference is real: does a variable,
such as treatment or body site, explain the distances between samples?

Each test takes a distance matrix (`dm`) from any beta diversity metric,
plus the variable to test: either a vector with one value per sample, or
the name of a column in a data frame of sample metadata (`data`).

| Function | Test | Variable | Question |
|:---|:---|:---|:---|
| [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md) | Permutational Multivariate Analysis of Variance | either | Do group centroids differ? How much variance (R²) does it explain? |
| [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md) | Homogeneity of Multivariate Dispersions | categorical | Do groups differ in their spread? |
| [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md) | Analysis of Similarities | categorical | Are distances within groups smaller than between them (ranks)? |
| [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) | Multi-Response Permutation Procedure | categorical | Are distances within groups smaller than expected? |
| [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) | Mantel Test | numeric | Do distances correlate with differences in a numeric variable, or with another distance matrix? |

They give the same statistics as `vegan`’s `adonis2()`, `betadisper()`,
[`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
[`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md), and
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md), with
permutations run in parallel C, and the p-values are identical for any
number of threads.

The examples use R’s built-in `iris` data, as in
[`vignette('ordination')`](https://cmmr.github.io/ecodive/articles/ordination.md):
150 flowers, four measurements each, and three species.

``` r

library(ecodive)

dm <- bray(iris[, 1:4])
```

## PERMANOVA

[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
is the most widely used test. Give the grouping variable as a column
name of `data`, or as a vector.

``` r

permanova(dm, 'Species', iris)
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F = 532.7, p = 0.001, z = 577.6
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.8788. Model tested as a whole:
#> 
#>           Df SumOfSqs     R2     F Pr(>F)    
#> Model      2   2.3173 0.8788 532.7  0.001 ***
#> Residual 147   0.3197 0.1212                 
#> Total    149   2.6370 1.0000                 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Species explain most of the variation between flowers (`R2`), and no
permutation of the labels did as well (`p`). The p-value cannot be
smaller than 1 / (`permutations` + 1), so use more permutations for a
more precise one.

A `formula` tests several variables at once. `by = 'terms'` tests each
in turn, as `adonis2()` does, and `strata` restricts permutations to
within blocks (such as subjects in a repeated-measures design).

``` r

permanova(dm, data = iris, formula = ~ Species + Sepal.Width, by = 'terms')
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F =   825, p = 0.001, z =  1105
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.9443. Terms added sequentially (first to last):
#> 
#>              Df SumOfSqs     R2      F Pr(>F)    
#> Species       2   2.3173 0.8788 1151.7  0.001 ***
#> Sepal.Width   1   0.1728 0.0655  171.8  0.001 ***
#> Residual    146   0.1469 0.0557                  
#> Total       149   2.6370 1.0000                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Checking the spread

PERMANOVA can also give a small p-value when groups have the same center
but different spreads.
[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
tests the spread on its own; a significant result means a PERMANOVA
difference may come partly or wholly from dispersion.

``` r

permdisp(dm, 'Species', iris)
#> Warning: Some squared distances to group centres are negative, and were set to
#> 0.
#> Homogeneity of multivariate dispersions (PERMDISP)
#> 
#> F = 0.5676, p = 0.573, z = -0.4401
#> Permutations: 999 (free)
#> 
#> Distances to group centroids:
#>     setosa versicolor  virginica 
#>    0.03781    0.04248    0.04103 
#> 
#>            Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
#> Groups      2 0.00057 0.0002852 0.568    999  0.573
#> Residuals 147 0.07386 0.0005024
```

Here the three species are equally spread, so the PERMANOVA result above
reflects a real difference in their centers. The warning comes from
Bray-Curtis distances not being Euclidean: a flower can then sit at a
negative squared distance from its species’ center, which
[`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
sets to zero (Anderson 2006).

## ANOSIM and MRPP

[`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md)
compares the ranks of distances within and between groups. `R` near 1
means every group is distinct, and near 0 means no difference.
[`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) compares
the average distance within groups with what random groups would give.

``` r

anosim(dm, 'Species', iris)
#> Analysis of similarities (ANOSIM)
#> 
#> R = 0.8576, p = 0.001, z =  79.4
#> Permutations: 999 (free)
#> 
#> Mean rank of the distances, between and within groups:
#>            pairs mean_rank
#> Between     7500   7163.90
#> setosa      1225   2239.66
#> versicolor  1225   2457.85
#> virginica   1225   2418.15

mrpp(dm, 'Species', iris)
#> Multi-response permutation procedure (MRPP)
#> 
#> delta = 0.05795, p = 0.001, z = -78.9
#> Permutations: 999 (free)
#> 
#> Mean within-group distances (weights: n):
#>        setosa versicolor virginica
#> delta 0.05502     0.0602   0.05861
#> n          50         50        50
#> 
#> Expected delta: 0.158
#> Chance-corrected within-group agreement A: 0.6332
```

## Mantel test

[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)
correlates `dm` with a numeric variable, or with a second distance
matrix: here, whether flowers with similar petals also have similar
sepals. `method` can be `'pearson'`, `'spearman'`, or `'kendall'`.

``` r

sepals <- bray(iris[, 1:2])
petals <- bray(iris[, 3:4])

mantel(sepals, petals, method = 'spearman')
#> Mantel test (Spearman)
#> 
#> r =  0.69, p = 0.001, z = 35.65
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 150
```

## Running tests by name

[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
runs any test by name, as
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
does for ordinations. It checks that the variable suits the test
(categorical or numeric), and names are matched forgivingly: `'adonis2'`
finds `permanova` and `'betadisper'` finds `permdisp`.

``` r

dist_test(dm, 'adonis2', 'Species', iris)
#> Permutational multivariate analysis of variance (PERMANOVA)
#> 
#> F = 532.7, p = 0.001, z = 577.6
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.8788. Model tested as a whole:
#> 
#>           Df SumOfSqs     R2     F Pr(>F)    
#> Model      2   2.3173 0.8788 532.7  0.001 ***
#> Residual 147   0.3197 0.1212                 
#> Total    149   2.6370 1.0000                 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Each result converts to a one-row data frame with the same columns, so
several tests combine into one table.

``` r

tests <- list_dist_tests(val = 'id', groups = 'categorical')

do.call(rbind, lapply(tests, function (test)
  as.data.frame(dist_test(dm, test, 'Species', iris)) ))
#> Warning: Some squared distances to group centres are negative, and were set to
#> 0.
#>        test stat_name    statistic p_value           z n_perm
#> 1    anosim         R   0.85763011   0.001  79.3963367    999
#> 2      mrpp     delta   0.05794548   0.001 -78.9020283    999
#> 3 permanova         F 532.74262193   0.001 577.5590509    999
#> 4  permdisp         F   0.56763808   0.573  -0.4401028    999
```

## Reproducibility

Permutations come from ecodive’s own random number generator, set by
`seed` (default `0`), so the same call always gives the same p-value, on
any computer and with any number of `cpus`. R’s
[`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect.

``` r

identical(
  permanova(dm, 'Species', iris, seed = 1, cpus = 1),
  permanova(dm, 'Species', iris, seed = 1, cpus = 2) )
#> [1] TRUE
```
