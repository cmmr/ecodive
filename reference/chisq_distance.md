# Chi-square distance

The Euclidean distance after the chi-square transformation: the distance
that correspondence analysis
([`ca()`](https://cmmr.github.io/ecodive/reference/ca.md)) preserves.

## Usage

``` r
chisq_distance(
  counts,
  margin = 1L,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)
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

The chi-square distance is defined as:
\$\$\sqrt{T\sum\_{i=1}^{n}\frac{1}{c_i}\left(\frac{X_i}{X_T} -
\frac{Y_i}{Y_T}\right)^2}\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute counts of the \\i\\-th feature.

- \\X_T\\, \\Y_T\\ : Total counts in each sample.

- \\c_i\\ : The total count of the \\i\\-th feature over all the samples
  of `counts` (features with \\c_i = 0\\ are left out).

- \\T\\ : The total count of the table.

- \\n\\ : The number of features.

It is a metric. It depends on the whole table, through \\c_i\\ and
\\T\\. It is `euclidean(counts, norm = 'chisq')` (see
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)).
It is not
[`squared_chisq()`](https://cmmr.github.io/ecodive/reference/squared_chisq.md)
or
[`psym_chisq()`](https://cmmr.github.io/ecodive/reference/psym_chisq.md),
which compare each pair of samples on its own.

It equals vegan's `vegdist(x, "chisq")`, which is Euclidean distance
after `decostand(x, "chi.square")` (Legendre & Gallagher 2001).
Differences from vegan: an empty sample is all zeros after the
transformation, so its distances are finite (vegan gives `NA`).

Base R Equivalent:

    x <- ex_counts
    x <- x[, colSums(x) > 0]
    y <- sqrt(sum(x)) * sweep(x / rowSums(x), 2, sqrt(colSums(x)), '/')
    dist(y)

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

Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
transformations for ordination of species data. *Oecologia*, 129,
271-280.
[doi:10.1007/s004420100716](https://doi.org/10.1007/s004420100716)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md),
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
[`chao_jaccard()`](https://cmmr.github.io/ecodive/reference/chao_jaccard.md),
[`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md),
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
    chisq_distance(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Chi-Square Distance
#>       Saliva Gums Nose
#> Gums    1.56
#> Nose    2.48 2.25
#> Stool   2.83 2.64 2.69
#> # Attributes: method
```
