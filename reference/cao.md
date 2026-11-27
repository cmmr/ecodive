# Cao dissimilarity (CYd)

The Cao (CY) index of dissimilarity (Cao et al. 1997), for counts.

## Usage

``` r
cao(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

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

## Details

The Cao dissimilarity is defined as:
\$\$\frac{1}{S}\sum\_{i=1}^{S}\left\[\ln\frac{X_i + Y_i}{2} -
\frac{X_i\ln Y_i + Y_i\ln X_i}{X_i + Y_i}\right\]\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute counts of the \\i\\-th feature, each
  raised to at least 0.1 (so zeros count as 0.1).

- \\S\\ : The number of features found in either sample (the sum is over
  those).

It is 0 for identical samples, has no upper limit, and uses natural
logarithms. It is intended for counts: the 0.1 that stands in for zeros
assumes that the smallest non-zero value is 1. ecodive doesn't check for
whole numbers (vegan warns).

It equals vegan's `vegdist(x, "cao")`. abdiv's `cy_dissimilarity()` is a
different scaling (base-10 logarithms).

Base R Equivalent:

    x <- ex_counts[1,]
    y <- ex_counts[2,]
    k <- x > 0 | y > 0
    x <- pmax(x[k], 0.1); y <- pmax(y[k], 0.1); n <- x + y
    mean(log(n / 2) - (x * log(y) + y * log(x)) / n)

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

Cao, Y., Williams, W. P., & Bark, A. W. (1997). Similarity measure bias
in river benthic Aufwuchs community analysis. *Water Environment
Research*, 69, 95-106.
[doi:10.2175/106143097X125227](https://doi.org/10.2175/106143097X125227)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
`vignette('bdiv_guide')`

Other Abundance metrics:
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alt_gower()`](https://cmmr.github.io/ecodive/reference/alt_gower.md),
[`bhattacharyya()`](https://cmmr.github.io/ecodive/reference/bhattacharyya.md),
[`binomial_deviance()`](https://cmmr.github.io/ecodive/reference/binomial_deviance.md),
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md),
[`canberra()`](https://cmmr.github.io/ecodive/reference/canberra.md),
[`chao_jaccard()`](https://cmmr.github.io/ecodive/reference/chao_jaccard.md),
[`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md),
[`chisq_distance()`](https://cmmr.github.io/ecodive/reference/chisq_distance.md),
[`chord()`](https://cmmr.github.io/ecodive/reference/chord.md),
[`clark()`](https://cmmr.github.io/ecodive/reference/clark.md),
[`divergence()`](https://cmmr.github.io/ecodive/reference/divergence.md),
[`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md),
[`gower()`](https://cmmr.github.io/ecodive/reference/gower.md),
[`hellinger()`](https://cmmr.github.io/ecodive/reference/hellinger.md),
[`horn()`](https://cmmr.github.io/ecodive/reference/horn.md),
[`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md),
[`jsd()`](https://cmmr.github.io/ecodive/reference/jsd.md),
[`kulczynski()`](https://cmmr.github.io/ecodive/reference/kulczynski.md),
[`lorentzian()`](https://cmmr.github.io/ecodive/reference/lorentzian.md),
[`mahalanobis_distance()`](https://cmmr.github.io/ecodive/reference/mahalanobis_distance.md),
[`manhattan()`](https://cmmr.github.io/ecodive/reference/manhattan.md),
[`matusita()`](https://cmmr.github.io/ecodive/reference/matusita.md),
[`minkowski()`](https://cmmr.github.io/ecodive/reference/minkowski.md),
[`morisita()`](https://cmmr.github.io/ecodive/reference/morisita.md),
[`motyka()`](https://cmmr.github.io/ecodive/reference/motyka.md),
[`psym_chisq()`](https://cmmr.github.io/ecodive/reference/psym_chisq.md),
[`robust_aitchison()`](https://cmmr.github.io/ecodive/reference/robust_aitchison.md),
[`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md),
[`squared_chisq()`](https://cmmr.github.io/ecodive/reference/squared_chisq.md),
[`squared_chord()`](https://cmmr.github.io/ecodive/reference/squared_chord.md),
[`squared_euclidean()`](https://cmmr.github.io/ecodive/reference/squared_euclidean.md),
[`topsoe()`](https://cmmr.github.io/ecodive/reference/topsoe.md),
[`wave_hedges()`](https://cmmr.github.io/ecodive/reference/wave_hedges.md)

## Examples

``` r
    cao(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Cao Dissimilarity (CYd)
#>       Saliva Gums Nose
#> Gums   0.535
#> Nose   4.11  3.9
#> Stool  3.06  2.83 4.15
#> # Attributes: method
```
