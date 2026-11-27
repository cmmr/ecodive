# Mahalanobis distance

The Euclidean distance after the features are centred, scaled to unit
variance, and decorrelated: distances that account for the features'
covariance.

## Usage

``` r
mahalanobis_distance(
  counts,
  margin = 1L,
  norm = "none",
  pseudocount = NULL,
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

- norm:

  Normalize the incoming counts. Options are:

  - `'none'`: No transformation.

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'hellinger'`: Square root of relative abundance (the Hellinger
    transformation).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  - `'rclr'`: Robust centered log ratio.

  - `'log'`: Log relative abundance (log-TSS), with `pseudocount` added
    to every value first.

  - `'chisq'`: The chi-square transformation, which uses the whole
    table's feature totals (see
    [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)).

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

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

The Mahalanobis distance is defined as: \$\$\sqrt{(X -
Y)^{T}\Sigma^{-1}(X - Y)}\$\$

Where:

- \\X\\, \\Y\\ : The two samples' vectors of counts (after `norm`).

- \\\Sigma\\ : The features' sample covariance matrix (with \\n - 1\\),
  over all the samples of `counts`.

It is a metric. It depends on the whole table: removing samples changes
the covariance and so every distance. ecodive computes it without
inverting \\\Sigma\\: from a QR decomposition of the centred table \\X_c
= QR\\, the whitened table \\\sqrt{n - 1}\\X_c R^{-1}\\, whose Euclidean
distances are the Mahalanobis distances.

It needs more samples than features, and no feature that is constant or
a linear combination of others (a diagonal element of \\R\\ below
\\10^{-9}\\ times the largest column norm counts as zero). Otherwise it
stops with an error: reduce the features first, e.g. to the leading axes
of [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md).
`norm = 'clr'`, `'percent'` and `'chisq'` make the features collinear
(each sample's values have a fixed sum), so they always stop it.

It equals vegan's `vegdist(x, "mahalanobis")`. Differences from vegan:
where \\\Sigma\\ is singular, vegan returns nearly constant distances
with a warning, and ecodive an error. Named `mahalanobis_distance()`,
not [`mahalanobis()`](https://rdrr.io/r/stats/mahalanobis.html), to
leave [`stats::mahalanobis()`](https://rdrr.io/r/stats/mahalanobis.html)
unmasked; `beta_div(x, 'mahalanobis')` finds it. Unlike the other
metrics, it accepts negative values.

Base R Equivalent:

    x <- as.matrix(iris[, 1:4])
    dist(scale(x, scale = FALSE) %*% solve(chol(cov(x))))

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

## Pseudocount

The `pseudocount` parameter is only relevant when `norm = 'clr'` or
`norm = 'log'`.

Zeros are undefined in the centered log-ratio (CLR) transformation and
in the log of proportions. If `norm` is `'clr'` or `'log'`,
`pseudocount` is `NULL` (the default), and zeros are detected, the
function uses half the minimum non-zero value (`min(x[x>0]) / 2`) and
issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`), or
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
whose result needs no pseudocount.

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

See [`aitchison`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Mahalanobis, P. C. (1936). On the generalized distance in statistics.
*Proceedings of the National Institute of Sciences of India*, 2, 49-55.

Mardia, K. V., Kent, J. T., & Bibby, J. M. (1979). *Multivariate
Analysis*. Academic Press.

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
    # It needs more samples than features: 150 flowers, 4 measurements.
    dm <- mahalanobis_distance(iris[, 1:4])
    dm
#> # A dist: 150 samples (lower triangle)
#> # Mahalanobis Distance
#>        1    2     3     4     5     6    7     8    9   10
#> 2  1.35
#> 3  0.969 0.97
#> 4  1.41  1.45 0.717
#> 5  0.59  1.81 1.13  1.33
#> 6  1.14  2.45 1.92  2.24  0.945
#> 7  1.36  1.94 0.972 0.951 1.04  1.74
#> 8  0.443 1.34 0.779 1.02  0.532 1.37  1.13
#> 9  1.84  1.42 0.93  0.62  1.87  2.75  1.31 1.52
#> 10 1.13  1.02 0.924 0.928 1.36  2.19  1.63 0.847 1.24
#> 11 0.697 1.89 1.66  2.03  0.932 0.903 1.93 1.02  2.52 1.63
#> # … 139 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method
```
