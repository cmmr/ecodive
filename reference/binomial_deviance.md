# Binomial deviance dissimilarity

The deviance of a binomial model of each feature's split between two
samples (Anderson & Millar 2004).

## Usage

``` r
binomial_deviance(
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

The binomial deviance dissimilarity is defined as:
\$\$\sum\_{i=1}^{n}\frac{1}{N_i}\left\[X_i\ln\frac{X_i}{N_i} +
Y_i\ln\frac{Y_i}{N_i} - N_i\ln\frac{1}{2}\right\]\$\$

Where:

- \\X_i\\, \\Y_i\\ : Absolute counts of the \\i\\-th feature.

- \\N_i = X_i + Y_i\\, over the features with \\N_i \> 0\\, and \\0 \ln
  0 = 0\\.

- \\n\\ : The number of features.

It is 0 for identical samples (`NaN` for two empty samples) and has no
upper limit: each feature found in only one sample adds \\\ln 2\\. With
`norm = 'binary'` it is \\\ln(2) (A + B - 2J)\\, as vegan's
`vegdist(x, "binomial", binary = TRUE)`. The values depend on the scale
of the counts: proportions give other values than counts.

It equals vegan's `vegdist(x, "binomial")` and abdiv's
`binomial_deviance()`. Named `binomial_deviance()`, not
[`binomial()`](https://rdrr.io/r/stats/family.html), to leave
[`stats::binomial()`](https://rdrr.io/r/stats/family.html) unmasked;
`beta_div(x, 'binomial')` finds it.

Base R Equivalent:

    x <- ex_counts[1,]
    y <- ex_counts[2,]
    k <- x + y > 0; x <- x[k]; y <- y[k]; n <- x + y
    xl <- ifelse(x > 0, x * log(x / n), 0)
    yl <- ifelse(y > 0, y * log(y / n), 0)
    sum((xl + yl - n * log(1/2)) / n)

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

Anderson, M. J., & Millar, R. B. (2004). Spatial variation and effects
of habitat on temperate reef fish assemblages in northeastern New
Zealand. *Journal of Experimental Marine Biology and Ecology*, 305,
191-221.
[doi:10.1016/j.jembe.2003.12.011](https://doi.org/10.1016/j.jembe.2003.12.011)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
`vignette('bdiv_guide')`

Other Abundance metrics:
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alt_gower()`](https://cmmr.github.io/ecodive/reference/alt_gower.md),
[`bhattacharyya()`](https://cmmr.github.io/ecodive/reference/bhattacharyya.md),
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
    binomial_deviance(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Binomial Deviance Dissimilarity
#>       Saliva Gums Nose
#> Gums    1.05
#> Nose    3.01 3.24
#> Stool   4.07 3.35 3.28
#> # Attributes: method
```
