# Mountford dissimilarity

Mountford's (1962) index of similarity, from Fisher's log-series, as a
dissimilarity.

## Usage

``` r
mountford(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
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

The Mountford dissimilarity is \\\ln 2 - M\\, where \\M\\ is the
positive root of \$\$e^{AM} + e^{BM} = 1 + e^{(A + B - J)M}\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

It ranges from 0 to \\\ln 2\\: \\\ln 2\\ when the samples share no
features (\\J = 0\\), and 0 when one sample's features are all in the
other (\\J = A\\ or \\J = B\\), so not only for identical samples. It is
`NaN` when either sample is empty. It is not a metric.

\\M\\ is solved by Newton's method, safeguarded by bisection, to full
double precision, in a form that can't overflow. It equals vegan's
`vegdist(x, "mountford")` to about \\10^{-10}\\ while \\A\\ and \\B\\
are at most about 50. Differences from vegan: for richer samples,
vegan's iterations can stop far from the root (by up to 0.65 at a few
hundred features), and ecodive's value is the root.

Base R Equivalent:

    x <- ex_counts[1,] > 0
    y <- ex_counts[2,] > 0
    A <- sum(x); B <- sum(y); J <- sum(x & y)
    g <- function (M) expm1(-(B - J) * M) + expm1(-(A - J) * M) - expm1(-(A + B - J) * M)
    log(2) - uniroot(g, c(1e-9, 10), tol = 1e-12)$root  # for 0 < J < min(A, B)

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

Mountford, M. D. (1962). An index of similarity and its application to
classification problems. In P. W. Murphy (Ed.), *Progress in Soil
Zoology* (pp. 43-50). Butterworths.

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
`vignette('bdiv_guide')`

Other Presence/Absence metrics:
[`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md),
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
[`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md),
[`raup_crick()`](https://cmmr.github.io/ecodive/reference/raup_crick.md),
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)

## Examples

``` r
    mountford(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Mountford Dissimilarity
#>       Saliva   Gums Nose
#> Gums   0
#> Nose   0     0
#> Stool  0.265 0.0172    0
#> # Attributes: method
```
