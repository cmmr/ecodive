# Balanced variation and abundance gradients of abundance dissimilarity

Baselga's partition of the Bray-Curtis and Ruzicka (Soergel)
dissimilarities between two samples into **balanced variation** in
abundance (some features gain what others lose) and **abundance
gradients** (every feature is less abundant in one sample).

## Usage

``` r
bray_balanced(
  counts,
  margin = 1L,
  norm = "none",
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

bray_gradient(
  counts,
  margin = 1L,
  norm = "none",
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

soergel_balanced(
  counts,
  margin = 1L,
  norm = "none",
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

soergel_gradient(
  counts,
  margin = 1L,
  norm = "none",
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Values
  must not be negative.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- norm:

  Normalize the incoming counts first: `'none'` (the default, as
  betapart and abdiv), `'percent'` (relative abundances), `'chord'`,
  `'hellinger'` or `'binary'`. `'clr'`, `'rclr'` and `'log'` give
  negative values and are refused.

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- clust:

  Whether and how to cluster the samples (rows and columns) before
  returning the distance matrix: `NULL` (the default) or `FALSE` to keep
  them in their input order, or one of
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)'s linkage
  methods, `"average"` (UPGMA), `"complete"`, `"single"`, `"ward.D"`,
  `"ward.D2"`, `"mcquitty"` (WPGMA), `"median"` or `"centroid"`, or
  `"pcoa"` to order them along the first axis of
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).
  Clustering **reorders the samples** of the result, with their labels
  (samples without names are labelled with their positions in the
  input): its `order` attribute gives each one's position in the input,
  so that `meta[attr(d, "order"), ]` lines metadata up with them.
  ecodive's ordinations and tests put the samples back in their input
  order before analysing them, so give them metadata in its usual order.
  Other packages see the reordered matrix: give them metadata reordered
  by `order`, or matched by name. The matrix's `display` list keeps
  `clust` (`FALSE` as `"none"`; see
  [ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)),
  which
  [`plot.ecodive_dist()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_dist.md)
  follows: with `NULL` it clusters the samples for the heatmap
  (`"average"`), and with `"none"` it keeps their order. `clust` needs
  every distance, so it cannot be used with `pairs`.
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md)'s
  similarities are clustered as \\1 - \rho_p\\.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A distance matrix: a `dist` object of class `c("ecodive_dist", "dist")`,
as [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
returns.

## Details

With abundances \\x_k\\ and \\y_k\\ of each feature in the two samples,
\\A = \sum_k \min(x_k, y_k)\\ (the abundance they share), \\B = \sum_k
\max(x_k - y_k, 0)\\ and \\C = \sum_k \max(y_k - x_k, 0)\\ (what only
the first and only the second has); Baselga (2013, 2017):

|  |  |  |
|----|----|----|
| Function | Component | Formula |
| [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) | Bray-Curtis dissimilarity (total) | \\\frac{B + C}{2A + B + C}\\ |
| `bray_balanced()` | balanced variation, \\d\_{BC-bal}\\ | \\\frac{\min(B, C)}{A + \min(B, C)}\\ |
| `bray_gradient()` | abundance gradient, \\d\_{BC-gra}\\ | \\d\_{BC} - d\_{BC-bal}\\ |
| [`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md) | Ruzicka dissimilarity (total) | \\\frac{B + C}{A + B + C}\\ |
| `soergel_balanced()` | balanced variation, \\d\_{RUZ-bal}\\ | \\\frac{2\min(B, C)}{A + 2\min(B, C)}\\ |
| `soergel_gradient()` | abundance gradient, \\d\_{RUZ-gra}\\ | \\d\_{RUZ} - d\_{RUZ-bal}\\ |

[`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md) is
the Ruzicka dissimilarity (vegan's quantitative Jaccard), so its
components are named after it; `"ruzicka_balanced"` and
`"ruzicka_gradient"` select them in
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).
The two components add up to the total,
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) or
[`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md). On
presence/absence data (`norm = 'binary'`) they are the turnover and
nestedness of
[`sorensen_turnover()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md).

When a sample is empty (all zeros), both components are undefined
(`NaN`), as in abdiv and betapart. Neither is a metric.

**Same names in other packages:** abdiv's functions are named
`bray_curtis_balanced()`, `bray_curtis_gradient()`, `ruzicka_balanced()`
and `ruzicka_gradient()`, taking two vectors; they give the same values.
betapart's `beta.pair.abund()` returns all three matrices of a family,
as
[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md)
does.

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

Baselga, A. (2013). Separating the two components of abundance-based
dissimilarity: balanced changes in abundance vs. abundance gradients.
*Methods in Ecology and Evolution*, 4, 552-557.
[doi:10.1111/2041-210X.12029](https://doi.org/10.1111/2041-210X.12029)

Baselga, A. (2017). Partitioning abundance-based multiple-site
dissimilarity into components: balanced variation in abundance and
abundance gradients. *Methods in Ecology and Evolution*, 8, 799-808.
[doi:10.1111/2041-210X.12693](https://doi.org/10.1111/2041-210X.12693)

## See also

[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md)
for the three matrices of a family at once,
[`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md)
for whole sets of samples,
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).

Other Partitioned metrics:
[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md),
[`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md),
[`sorensen_turnover()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md)

## Examples

``` r
    bray_balanced(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Balanced Variation
#>       Saliva  Gums Nose
#> Gums   0.27
#> Nose   0.922 0.968
#> Stool  0.988 0.989 0.99
#> # Attributes: method
    bray_gradient(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Abundance Gradient
#>        Saliva    Gums    Nose
#> Gums  0.321
#> Nose  0.0384  0.00208
#> Stool 0.00326 0.00205 0.00238
#> # Attributes: method

    # Of relative abundances.
    soergel_balanced(ex_counts, norm = 'percent')
#> # A dist: 4 samples (lower triangle)
#> # Ruzicka Balanced Variation
#>       Saliva  Gums  Nose
#> Gums   0.598
#> Nose   0.985 0.986
#> Stool  0.995 0.996 0.996
#> # Attributes: method
```
