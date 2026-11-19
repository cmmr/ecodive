# Structural zeros

Finds the features that are absent, or nearly absent, from a group of
samples (Kaul et al. 2017), as scikit-bio's `struc_zero()` and ANCOMBC's
`struc_zero` option. A feature that is a structural zero in a group is
differentially (less) abundant there, whatever a test of its abundance
says, and should be left out of analyses of that group; one that is a
structural zero in every group should be removed.

## Usage

``` r
structural_zeros(counts, groups, data = NULL, margin = 1L, neg_lb = FALSE)
```

## Arguments

- counts:

  A numeric matrix of abundances (samples \\\times\\ features), or any
  of the objects listed under **Input Types**: the table **before** any
  pseudocount or zero replacement.

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

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- neg_lb:

  `FALSE` (the default): a feature is a structural zero in a group when
  it is zero in every sample of that group. `TRUE`: when the lower bound
  of the 95% interval of its proportion of non-zero samples, \\\hat p -
  1.96 \sqrt{\hat p (1 - \hat p) / n}\\, is at most 0; recommended for
  large groups.

## Value

A logical matrix, features \\\times\\ groups, `TRUE` where the feature
is a structural zero in the group. Groups are in sorted order
(`sort(method = 'radix')`, as scikit-bio's), or in a factor's level
order.

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

Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
microbiome data in the presence of excess zeros. *Frontiers in
Microbiology*, 8, 2114.

## See also

[`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
[`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md).

## Examples

``` r
    x <- rbind(
      c( 7, 1,  0, 11,  3,  1), c( 1, 1, 0, 13, 13,  0),
      c(11, 5,  0,  1,  4,  1), c( 2, 2, 0, 16,  4,  0),
      c( 0, 1,  0,  0,  6,  0), c(14, 8, 7,  9,  0,  5),
      c( 0, 7,  4,  1,  0, 26), c( 8, 1, 4, 28,  0, 10),
      c( 2, 2,  2,  4,  0,  5), c( 6, 4, 10,  1,  0,  9) )
    colnames(x) <- paste0('f', 0:5)
    grp <- rep(c('treatment', 'placebo'), each = 5)
    structural_zeros(x, grp)
#>    placebo treatment
#> f0   FALSE     FALSE
#> f1   FALSE     FALSE
#> f2   FALSE      TRUE
#> f3   FALSE     FALSE
#> f4    TRUE     FALSE
#> f5   FALSE     FALSE
    structural_zeros(x, grp, neg_lb = TRUE)
#>    placebo treatment
#> f0   FALSE     FALSE
#> f1   FALSE     FALSE
#> f2   FALSE      TRUE
#> f3   FALSE     FALSE
#> f4    TRUE     FALSE
#> f5   FALSE      TRUE
```
