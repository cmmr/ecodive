# Partial Mantel Test

Tests whether two distance matrices are correlated after removing the
effect of a third: for instance, whether samples with more different
communities are also more different in pH, beyond what geographic
distance between them explains. Gives the same statistic as
`vegan::mantel.partial()` (Pearson, Spearman, or Kendall partial
correlation), with permutation p-values, computed in parallel C straight
from the distances.
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) is the
simple Mantel test.

## Usage

``` r
mantel_partial(
  dm,
  y,
  given,
  data = NULL,
  method = "pearson",
  permutations = 999L,
  strata = NULL,
  alternative = "greater",
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- dm:

  A `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  or [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), or a
  symmetric numeric matrix of distances. Must not contain `NA`s or
  negative values. Its samples are the ones permuted.

- y:

  The distances to correlate with `dm`: a `dist` object or a symmetric
  numeric matrix, like `dm`, matched to `dm`'s samples by name when both
  have names, otherwise in the same order. Or numeric values, such as
  age or pH, whose absolute differences (as `stats::dist(y)`) are the
  distances: with `data`, the name of one of its numeric columns or a
  vector with one value per row of `data`; without `data`, a vector with
  one value per sample, matched to the samples by name when it has
  names, otherwise in the same order.

- given:

  The distances whose effect is removed from the correlation of `dm` and
  `y`, as `vegan::mantel.partial()`'s `zdis`: such as geographic
  distances between the samples. Given in any of the forms `y` takes,
  including as a column name of `data`. Required.

- data:

  Optional sample metadata: a data frame with one row per sample, whose
  columns `y`, `given`, and `strata` can name. Rows are matched to the
  samples in `dm` by row name, and rows for other samples are ignored;
  without row names, they must be in the same order as `dm`. Default:
  `NULL`.

- method:

  The correlation coefficient: `'pearson'`, `'spearman'`, or `'kendall'`
  (tau-b), as in [`stats::cor()`](https://rdrr.io/r/stats/cor.html).
  Default: `'pearson'`.

- permutations:

  The number of permutations. Default: `999`. Alternatively, a matrix
  with one permutation of `1:n` per row, as from
  `permute::shuffleSet()`, to use instead of random ones.

- strata:

  Permute samples only within these groups (blocks), for designs such as
  repeated measures. With `data`, the name of one of its columns, or a
  vector with one value per row of `data`. Without `data`, a vector with
  one value per sample, matched to the samples by name when it has
  names, otherwise in the same order. Default: `NULL` (free
  permutation).

- alternative:

  Which correlations count as more extreme than the observed one:
  `'greater'` (more positive, as in `vegan::mantel()`), `'less'` (more
  negative), or `'two.sided'` (larger in absolute value). Default:
  `'greater'`.

- seed:

  An integer seed for the random permutations. The same seed gives the
  same p-values. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A list of class `c('ecodive_mantel', 'ecodive_test')`, as from
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md), with:

- `statistic`: The partial Mantel statistic, `c(r = ...)`: the partial
  correlation of `dm` and `y` given `given`. `vegan::mantel.partial()`'s
  `statistic`.

- `p_value`: The permutation p-value.

- `z`: The standardized effect size, `(r - mean(r_perm)) / sd(r_perm)`.

- `n_perm`: The number of permutations.

- `method`: The correlation coefficient: `'pearson'`, `'spearman'`, or
  `'kendall'`.

- `alternative`: `'greater'`, `'less'`, or `'two.sided'`.

- `partial`: `TRUE`.

- `n`: The number of samples.

- `permutation`: How samples were permuted: `'free'`, `'within strata'`,
  or `'supplied matrix'`.

Attribute `r_perm` is a vector of the statistics of the permutations, as
`vegan::mantel.partial()`'s `perm`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
main results as one row, with `test = 'mantel_partial'`.

## Details

**Method**

The partial Mantel test (Smouse et al. 1986; Legendre & Legendre 2012)
correlates the \\n(n-1)/2\\ distances of `dm` (\\x\\) and `y` after
removing the effect of `given` (\\z\\), as the partial correlation
\$\$r\_{xy.z} = \frac{r\_{xy} - r\_{xz} r\_{yz}}{\sqrt{1 - r\_{xz}^2}
\sqrt{1 - r\_{yz}^2}}\$\$ where each \\r\\ is the Pearson, Spearman, or
Kendall correlation of two sets of distances. As in
`vegan::mantel.partial()`, the p-value comes from permutations of the
samples of `dm`, with `y` and `given` held fixed: \\r\_{xy}\\ and
\\r\_{xz}\\ are recomputed for each permutation, and \\r\_{yz}\\ is
computed once. The p-value is the proportion of permutations, counting
the observed one, with a partial correlation at least as extreme, as in
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md).

Its p-values can be too small when the distances are spatially
autocorrelated (Legendre & Fortin 2010; Guillot & Rousset 2013).

When `y` and `given`, or `dm` and `given`, are perfectly correlated,
there is no partial correlation: the statistic and p-value are `NA`,
with a warning.

**Speed, memory, and reproducibility**

As for [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md):
each permutation is one pass over the distances (two correlations,
\\r\_{xy}\\ and \\r\_{xz}\\, in the same pass), with no \\n \times n\\
matrix, and permutations run in parallel. They come from `seed`, so
p-values are identical from run to run and for any number of `cpus`, and
agree with `vegan::mantel.partial()`'s exactly when the same
`permutations` matrix is given to both.

**Differences from vegan::mantel.partial()**

- `y` and `given` may also be numeric vectors, one value per sample, or
  names of columns of `data`.

- Samples are matched by name between `dm`, `y`, `given`, `data`, and
  `strata`.

- `alternative` adds lower-tail and two-sided tests.

- There's no `na.rm`: distances must not be `NA`.

## References

Guillot, G., & Rousset, F. (2013). Dismantling the Mantel tests.
*Methods in Ecology and Evolution*, 4(4), 336-344.
[doi:10.1111/2041-210x.12018](https://doi.org/10.1111/2041-210x.12018)

Legendre, P., & Fortin, M.-J. (2010). Comparison of the Mantel test and
alternative approaches for detecting complex multivariate relationships
in the spatial analysis of genetic data. *Molecular Ecology Resources*,
10(5), 831-844.
[doi:10.1111/j.1755-0998.2010.02866.x](https://doi.org/10.1111/j.1755-0998.2010.02866.x)

Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
section 10.5. Elsevier.

Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression
and correlation extensions of the Mantel test of matrix correspondence.
*Systematic Zoology*, 35(4), 627-632.
[doi:10.2307/2413122](https://doi.org/10.2307/2413122)

## See also

[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md),
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)

## Examples

``` r
    # 30 samples along a pH gradient, which shifts the first five genera,
    # at random locations.
    set.seed(1)
    counts <- matrix(rpois(30 * 20, 5), 30, 20)
    rownames(counts) <- paste0('S', 1:30)
    meta <- data.frame(
      row.names = rownames(counts),
      ph  = round(runif(30, 6, 8), 1),
      lat = runif(30, 30, 35),
      lon = runif(30, -100, -95) )
    counts[, 1:5] <- counts[, 1:5] + round(4 * (meta$ph - 6))

    dm  <- bray(counts)
    geo <- dist(meta[, c('lat', 'lon')])

    # Are samples that differ more in pH also more different in
    # community composition, beyond the effect of geographic distance?
    mantel_partial(dm, 'ph', given = geo, data = meta)
#> Partial Mantel test (Pearson) of Bray-Curtis Dissimilarity
#> 
#> r = 0.5931, p = 0.001, z = 9.759
#> Permutations: 999 (free)
#> 
#> Correlation: Pearson's product-moment correlation
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30

    # With rank correlation.
    mantel_partial(dm, 'ph', given = geo, data = meta, method = 'spearman')
#> Partial Mantel test (Spearman) of Bray-Curtis Dissimilarity
#> 
#> r = 0.5528, p = 0.001, z = 9.969
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30

    # The same, through dist_test().
    dist_test(dm, 'mantel_partial', 'ph', meta, given = geo, method = 'spearman')
#> Partial Mantel test (Spearman) of Bray-Curtis Dissimilarity
#> 
#> r = 0.5528, p = 0.001, z = 9.969
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30
```
