# Chao-Jaccard dissimilarity

Chao et al.'s (2005) abundance-based Jaccard dissimilarity, which
estimates the shares of individuals that belong to shared species,
including shared species not seen.

## Usage

``` r
chao_jaccard(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
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

The Chao-Jaccard dissimilarity is defined as: \$\$1 - \frac{UV}{U + V -
UV}\$\$

with \$\$U = \min\left(1,\\ \frac{C_X}{N_X} + \frac{N_Y -
1}{N_Y}\cdot\frac{f\_{1}^{Y}}{2f\_{2}^{Y}}\cdot\frac{S_X}{N_X}\right)\$\$
and \\V\\ the same with \\X\\ and \\Y\\ swapped.

Where:

- \\N_X\\, \\N_Y\\ : The total counts of each sample.

- \\C_X\\ : The counts in \\X\\ of the features shared with \\Y\\.

- \\f\_{1}^{Y}\\, \\f\_{2}^{Y}\\ : The number of shared features with 1
  and with 2 counts in \\Y\\. When \\f\_{2}^{Y} = 0\\, 1 is used.

- \\S_X\\ : The counts in \\X\\ of the shared features with 1 count in
  \\Y\\.

It ranges from 0 to 1: 1 when the samples share no features (and `NaN`
when both are empty). It is not a metric.

It equals vegan's `vegdist(x, "chao")`. It is the Jaccard-type index;
vegan's `chaodist()` offers the Sorensen type and others. Named
`chao_jaccard()`, not `chao()`, to keep it apart from the alpha metric
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md);
`beta_div(x, 'chao')` finds it.

Base R Equivalent:

    x <- ex_counts[1,]
    y <- ex_counts[2,]
    u <- function (x, y) {
      sh <- x > 0 & y > 0
      f1 <- sum(sh & y == 1); f2 <- max(1, sum(sh & y == 2))
      min(1, sum(x[sh]) / sum(x) +
        (sum(y) - 1) / sum(y) * f1 / (2 * f2) * sum(x[sh & y == 1]) / sum(x))
    }
    U <- u(x, y); V <- u(y, x)
    1 - U * V / (U + V - U * V)

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

Chao, A., Chazdon, R. L., Colwell, R. K., & Shen, T.-J. (2005). A new
statistical approach for assessing similarity of species composition
with incidence and abundance data. *Ecology Letters*, 8, 148-159.
[doi:10.1111/j.1461-0248.2004.00707.x](https://doi.org/10.1111/j.1461-0248.2004.00707.x)

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
[`cao()`](https://cmmr.github.io/ecodive/reference/cao.md),
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
    chao_jaccard(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Chao-Jaccard Dissimilarity
#>         Saliva     Gums Nose
#> Gums  0.000566
#> Nose  0.602    0.0114
#> Stool 0.00325  0.000814    0
#> # Attributes: method
```
