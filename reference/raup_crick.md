# Raup-Crick dissimilarity

The probability, under random placement of features, of the two samples
sharing at least as many features as they do (Raup & Crick 1979): the
analytic, presence/absence form of vegan's `vegdist()`.

## Usage

``` r
raup_crick(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
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

The Raup-Crick dissimilarity is: \$\$P(X \ge J), \quad X \sim
\mathrm{Hypergeometric}(S, A, B)\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

- \\S\\ : The number of features (columns) of `counts`, **including
  features absent from every sample**: adding an empty feature changes
  every value.

- \\X\\ : The number of features two random samples of \\A\\ and \\B\\
  of the \\S\\ features have in common, with all features equally
  likely.

It ranges from 0 to 1, and is not a metric. Samples that share no
features have a dissimilarity of exactly 1, while identical samples can
be slightly above 0. The smaller tail of the hypergeometric distribution
is summed directly, so small values keep their precision.

It equals vegan's `vegdist(x, "raup")`, which computes \\1 - P(X \< J)\\
and so returns 0 where the true value is below about \\10^{-16}\\.

This is not vegan's `raupcrick()` (Chase et al. 2011), which simulates
null communities with features drawn in proportion to their frequencies,
nor Stegen et al.'s abundance-based Raup-Crick (\\RC\_{bray}\\).

Base R Equivalent:

    x <- ex_counts[1,] > 0
    y <- ex_counts[2,] > 0
    S <- length(x); A <- sum(x); B <- sum(y); J <- sum(x & y)
    phyper(J - 1, A, S - A, B, lower.tail = FALSE)

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

Raup, D. M., & Crick, R. E. (1979). Measurement of faunal similarity in
paleontology. *Journal of Paleontology*, 53, 1213-1227.

Chase, J. M., Kraft, N. J. B., Smith, K. G., Vellend, M., & Inouye, B.
D. (2011). Using null models to disentangle variation in community
dissimilarity from variation in alpha-diversity. *Ecosphere*, 2, art24.
[doi:10.1890/ES10-00117.1](https://doi.org/10.1890/ES10-00117.1)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
`vignette('bdiv_guide')`

Other Presence/Absence metrics:
[`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md),
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
[`mountford()`](https://cmmr.github.io/ecodive/reference/mountford.md),
[`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md),
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)

## Examples

``` r
    raup_crick(ex_counts)
#> # A dist: 4 samples (lower triangle)
#> # Raup-Crick Dissimilarity
#>       Saliva Gums Nose
#> Gums   0.333
#> Nose   1        1
#> Stool  1        1    1
#> # Attributes: method
```
