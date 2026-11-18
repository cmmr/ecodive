# Mantel Test

Tests whether two distance matrices are correlated: for instance,
whether samples with more different communities are also further apart,
or more different in pH. Gives the same statistic as
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)
(Pearson, Spearman, or Kendall correlation), and as
[`vegan::mantel.partial()`](https://vegandevs.github.io/vegan/reference/mantel.html)
with `given`, with permutation p-values. Computes them in parallel C,
straight from the distances, without the \\n \times n\\ matrices that
`mantel()` builds for every permutation.

## Usage

``` r
mantel(
  dm,
  y,
  data = NULL,
  given = NULL,
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

- data:

  Optional sample metadata: a data frame with one row per sample, whose
  columns `y`, `given`, and `strata` can name. Rows are matched to the
  samples in `dm` by row name, and rows for other samples are ignored;
  without row names, they must be in the same order as `dm`. Default:
  `NULL`.

- given:

  For a partial Mantel test: distances whose effect is removed from the
  correlation of `dm` and `y`, as in
  [`vegan::mantel.partial()`](https://vegandevs.github.io/vegan/reference/mantel.html)'s
  `zdis`. Given like `y`, including as a column name of `data`. Default:
  `NULL` (a simple Mantel test).

- method:

  The correlation coefficient: `'pearson'`, `'spearman'`, or `'kendall'`
  (tau-b), as in [`stats::cor()`](https://rdrr.io/r/stats/cor.html).
  Default: `'pearson'`.

- permutations:

  The number of permutations. Default: `999`. Alternatively, a matrix
  with one permutation of `1:n` per row, as from
  [`permute::shuffleSet()`](https://gavinsimpson.github.io/permute/reference/shuffleSet.html),
  to use instead of random ones.

- strata:

  Permute samples only within these groups (blocks), for designs such as
  repeated measures. With `data`, the name of one of its columns, or a
  vector with one value per row of `data`. Without `data`, a vector with
  one value per sample, matched to the samples by name when it has
  names, otherwise in the same order. Default: `NULL` (free
  permutation).

- alternative:

  Which correlations count as more extreme than the observed one:
  `'greater'` (more positive, as in
  [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)),
  `'less'` (more negative), or `'two.sided'` (larger in absolute value).
  Default: `'greater'`.

- seed:

  An integer seed for the random permutations. The same seed gives the
  same p-values. Default: `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A list of class `c('ecodive_mantel', 'ecodive_test')`, with:

- `statistic`: The Mantel statistic, `c(r = ...)`: the correlation of
  `dm` and `y`, or their partial correlation with `given`.
  [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)'s
  `statistic`.

- `p_value`: The permutation p-value.

- `z`: The standardized effect size, `(r - mean(r_perm)) / sd(r_perm)`.

- `n_perm`: The number of permutations.

- `method`: The correlation coefficient: `'pearson'`, `'spearman'`, or
  `'kendall'`.

- `alternative`: `'greater'`, `'less'`, or `'two.sided'`.

- `partial`: Whether it is a partial Mantel test (with `given`).

- `n`: The number of samples.

- `permutation`: How samples were permuted: `'free'`, `'within strata'`,
  or `'supplied matrix'`.

Attribute `r_perm` is a vector of the statistics of the permutations, as
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)'s
`perm`. [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives the main results as one row.

## Details

**Method**

The Mantel test (Mantel 1967; Legendre & Legendre 2012) correlates the
\\n(n-1)/2\\ distances of `dm` with the matching distances of `y`. The
distances aren't independent of one another, so the p-value comes from
permutations of the samples: the rows and columns of `dm` are shuffled
together, `y` is held fixed, and the correlation is recomputed. The
p-value is the proportion of permutations, counting the observed one,
with a correlation at least as extreme: \\(1 + \\\\r\_{perm} \ge r\\) /
(1 + permutations)\\ for `alternative = 'greater'`. As in
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html),
permuted statistics within `sqrt(.Machine$double.eps)` of the observed
one count as ties.

The partial Mantel test (Smouse et al. 1986) correlates `dm` and `y`
after removing the effect of `given`, as the partial correlation
\$\$r\_{xy.z} = \frac{r\_{xy} - r\_{xz} r\_{yz}}{\sqrt{1 - r\_{xz}^2}
\sqrt{1 - r\_{yz}^2}}\$\$ with the samples of `dm` permuted, as in
[`vegan::mantel.partial()`](https://vegandevs.github.io/vegan/reference/mantel.html).
Its p-values can be too small when the distances are spatially
autocorrelated (Legendre & Fortin 2010; Guillot & Rousset 2013).

**Speed and memory**

`dm` and `y` are centred and scaled once (after ranking them, for
Spearman), so each permuted correlation is one pass over the distances,
reading `dm` through the permutation: \\O(n^2)\\ per permutation, with
no \\n \times n\\ matrix. Kendall's tau-b is counted with Knight's
(1966) merge sort, in \\O(n^2 \log n)\\ per permutation, where
[`stats::cor()`](https://rdrr.io/r/stats/cor.html) takes \\O(n^4)\\.
Permutations run in parallel.

**Reproducibility**

The permutations come from `seed`, not from R's random number generator,
so [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect, and
each permutation has a random number stream of its own. P-values are
identical from run to run and for any number of `cpus`. The random
permutations differ from
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)'s,
so p-values agree with it only up to sampling error, unless the same
`permutations` matrix is given to both.

**Differences from vegan::mantel()**

- One function does both the simple and the partial test.

- `y` and `given` may also be numeric vectors, one value per sample, or
  names of columns of `data`.

- Samples are matched by name between `dm`, `y`, `given`, `data`, and
  `strata`.

- `alternative` adds lower-tail and two-sided tests.

- There's no `na.rm`: distances must not be `NA`.

- When there are fewer possible permutations than requested,
  [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)
  (via `permute`) runs each of them once instead. `mantel()` always
  draws random permutations.

## References

Guillot, G., & Rousset, F. (2013). Dismantling the Mantel tests.
*Methods in Ecology and Evolution*, 4(4), 336-344.
[doi:10.1111/2041-210x.12018](https://doi.org/10.1111/2041-210x.12018)

Knight, W. R. (1966). A computer method for calculating Kendall's tau
with ungrouped data. *Journal of the American Statistical Association*,
61(314), 436-439.
[doi:10.1080/01621459.1966.10480879](https://doi.org/10.1080/01621459.1966.10480879)

Legendre, P., & Fortin, M.-J. (2010). Comparison of the Mantel test and
alternative approaches for detecting complex multivariate relationships
in the spatial analysis of genetic data. *Molecular Ecology Resources*,
10(5), 831-844.
[doi:10.1111/j.1755-0998.2010.02866.x](https://doi.org/10.1111/j.1755-0998.2010.02866.x)

Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
section 10.5. Elsevier.

Mantel, N. (1967). The detection of disease clustering and a generalized
regression approach. *Cancer Research*, 27(2), 209-220.

Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression
and correlation extensions of the Mantel test of matrix correspondence.
*Systematic Zoology*, 35(4), 627-632.
[doi:10.2307/2413122](https://doi.org/10.2307/2413122)

## See also

[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)

## Examples

``` r
    # 30 samples along a pH gradient, which shifts the first five genera.
    set.seed(1)
    counts <- matrix(rpois(30 * 20, 5), 30, 20)
    rownames(counts) <- paste0('S', 1:30)
    ph  <- setNames(round(runif(30, 6, 8), 1), rownames(counts))
    geo <- dist(cbind(lat = runif(30, 30, 35), lon = runif(30, -100, -95)))
    counts[, 1:5] <- counts[, 1:5] + round(4 * (ph - 6))

    dm <- bray(counts)

    # Are samples that differ more in pH also more different in
    # community composition?
    mantel(dm, ph)
#> Mantel test (Pearson)
#> 
#> r = 0.593, p = 0.001, z = 9.786
#> Permutations: 999 (free)
#> 
#> Correlation: Pearson's product-moment correlation
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30

    # With rank correlation.
    mantel(dm, ph, method = 'spearman')
#> Mantel test (Spearman)
#> 
#> r = 0.5528, p = 0.001, z =    10
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30

    # After removing the effect of geographic distance.
    mantel(dm, ph, given = geo)
#> Partial Mantel test (Pearson)
#> 
#> r = 0.5931, p = 0.001, z = 9.759
#> Permutations: 999 (free)
#> 
#> Correlation: Pearson's product-moment correlation
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30

    # The same as a column of sample metadata.
    meta <- data.frame(row.names = rownames(counts), ph = ph)
    mantel(dm, 'ph', data = meta, method = 'spearman')
#> Mantel test (Spearman)
#> 
#> r = 0.5528, p = 0.001, z =    10
#> Permutations: 999 (free)
#> 
#> Correlation: Spearman's rank correlation rho
#> Alternative: greater (r is more positive than by chance)
#> Samples: 30
```
