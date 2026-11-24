# Turnover and nestedness of presence/absence dissimilarity

Baselga's partition of the Sorensen and Jaccard dissimilarities between
two samples into **turnover** (features replaced by others) and
**nestedness-resultant** dissimilarity (one sample holding a subset of
the other's features).

## Usage

``` r
sorensen_turnover(
  counts,
  margin = 1L,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

sorensen_nestedness(
  counts,
  margin = 1L,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

jaccard_turnover(
  counts,
  margin = 1L,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)

jaccard_nestedness(
  counts,
  margin = 1L,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Counts
  are converted to presence/absence. Values must not be negative.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

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

With \\a\\ features in both samples, and \\b\\ and \\c\\ only in the
first and only in the second (Baselga 2010, 2012):

|  |  |  |
|----|----|----|
| Function | Component | Formula |
| [`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md) | Sorensen dissimilarity (total) | \\\frac{b + c}{2a + b + c}\\ |
| `sorensen_turnover()` | Simpson dissimilarity, \\\beta\_{sim}\\ | \\\frac{\min(b, c)}{a + \min(b, c)}\\ |
| `sorensen_nestedness()` | nestedness-resultant, \\\beta\_{sne}\\ | \\\beta\_{sor} - \beta\_{sim}\\ |
| [`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md) | Jaccard dissimilarity (total) | \\\frac{b + c}{a + b + c}\\ |
| `jaccard_turnover()` | \\\beta\_{jtu}\\ | \\\frac{2\min(b, c)}{a + 2\min(b, c)}\\ |
| `jaccard_nestedness()` | \\\beta\_{jne}\\ | \\\beta\_{jac} - \beta\_{jtu}\\ |

Turnover is 0 when one sample's features are a subset of the other's,
and 1 when the samples share none; nestedness is 0 when they have the
same number of features. Turnover and nestedness add up to the total,
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md) or
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md).

When a sample has no features, turnover and nestedness are undefined
(`NaN`), as in abdiv and betapart, although the total is 1 (or `NaN`
when both samples are empty).

Neither component is a metric: a sample nested in two others has
turnover 0 to both, while they can have turnover 1 to each other.

The partition is debated: Podani & Schmera (2011) and Carvalho et al.
(2012) split the same totals differently, into replacement and richness
difference (adespatial's `beta.div.comp()`). See Baselga & Leprieur
(2015) for a comparison.

**Same names in other packages:** abdiv's functions are named
`sorenson_turnover()` and `sorenson_nestedness()` (with an "o"), and
`jaccard_turnover()` and `jaccard_nestedness()`, taking two vectors;
they give the same values. betapart's `beta.pair()` returns all three
matrices of a family, as
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

Baselga, A. (2010). Partitioning the turnover and nestedness components
of beta diversity. *Global Ecology and Biogeography*, 19, 134-143.
[doi:10.1111/j.1466-8238.2009.00490.x](https://doi.org/10.1111/j.1466-8238.2009.00490.x)

Baselga, A. (2012). The relationship between species replacement,
dissimilarity derived from nestedness, and nestedness. *Global Ecology
and Biogeography*, 21, 1223-1232.
[doi:10.1111/j.1466-8238.2011.00756.x](https://doi.org/10.1111/j.1466-8238.2011.00756.x)

Lennon, J. J., Koleff, P., Greenwood, J. J. D., & Gaston, K. J. (2001).
The geographical structure of British bird distributions: diversity,
spatial turnover and scale. *Journal of Animal Ecology*, 70, 966-979.
[doi:10.1046/j.0021-8790.2001.00563.x](https://doi.org/10.1046/j.0021-8790.2001.00563.x)

Podani, J., & Schmera, D. (2011). A new conceptual and methodological
framework for exploring and explaining pattern in presence-absence data.
*Oikos*, 120, 1625-1638.
[doi:10.1111/j.1600-0706.2011.19451.x](https://doi.org/10.1111/j.1600-0706.2011.19451.x)

Carvalho, J. C., Cardoso, P., & Gomes, P. (2012). Determining the
relative roles of species replacement and species richness differences
in generating beta-diversity patterns. *Global Ecology and
Biogeography*, 21, 760-771.
[doi:10.1111/j.1466-8238.2011.00694.x](https://doi.org/10.1111/j.1466-8238.2011.00694.x)

Baselga, A., & Leprieur, F. (2015). Comparing methods to separate
components of beta diversity. *Methods in Ecology and Evolution*, 6,
1069-1079.
[doi:10.1111/2041-210X.12388](https://doi.org/10.1111/2041-210X.12388)

## See also

[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md)
for the three matrices of a family at once,
[`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md)
for whole sets of samples,
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).

Other Partitioned metrics:
[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md),
[`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md),
[`bray_balanced()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md)

## Examples

``` r
    sorensen_turnover(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Simpson Dissimilarity (Sorensen Turnover)
#>       Saliva Gums Nose
#> Gums    0
#> Nose    0     0
#> Stool   0.25  0.2    0
#> # Attributes: method
    sorensen_nestedness(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Sorensen Nestedness
#>       Saliva      Gums   Nose
#> Gums  0.111
#> Nose  0.2       0.0909
#> Stool 0.0833 -5.55e-17 0.0909
#> # Attributes: method

    # Turnover + nestedness = total.
    all.equal(
      as.vector(jaccard_turnover(ex_counts) + jaccard_nestedness(ex_counts)),
      as.vector(jaccard(ex_counts)) )
#> [1] TRUE
```
