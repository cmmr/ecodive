# Proportionality Correlation

The proportionality coefficient \\\rho_p\\ of Erb & Notredame (2016),
propr's default `"rho"`: a correlation-like measure of how close each
pair of features is to keeping a constant ratio. `1` means perfectly
proportional.

## Usage

``` r
rho_p(
  counts,
  margin = 1L,
  denominator = NULL,
  pseudocount = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts,
  proportions, and the output of
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  are all accepted: each sample's scale cancels out of the log ratios.
  Pairs are taken between the **features** (the columns when
  `margin = 1`), and the variance is over the samples.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- denominator:

  The reference that each feature's log abundance is taken relative to,
  for the variances \\\mathrm{var}(z_i)\\. `NULL` (the default) uses the
  centered log ratio (clr): each sample's geometric mean of all
  features, as propr's default `ivar = "clr"`. Otherwise features, as
  names, (1-based) indices, or a logical vector with one value per
  feature: their geometric mean. One feature gives the additive log
  ratio (alr, as `alr(denominator =)`), and that feature is then left
  out of the result (its alr is 0 in every sample); with several, every
  feature is kept.

- pseudocount:

  Added to every value before taking logs. `NULL` (the default) adds
  nothing when `counts` has no zeros, and otherwise uses half the
  smallest non-zero value, with a warning. `0` adds nothing, and is an
  error when `counts` has zeros. See **Zeros**.

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
  every distance, so it cannot be used with `pairs`. `rho_p()`'s
  similarities are clustered as \\1 - \rho_p\\.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A `dist` object over the features, as
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) returns, with
`method` `"rho_p"`.

\\\rho_p\\ is a **similarity** in a `dist` container:
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) puts 0 on the
diagonal, where \\\rho_p\\ is 1, and it must not be given as `dm` to an
ordination or test. Threshold it to find proportional pairs
(`which(rho > 0.7)`, or `as.matrix(rho) > 0.7` for an adjacency matrix),
or use `1 - rho_p(x)` (which is \\\mathrm{vlr}\_{ij} /
(\mathrm{var}(z_i) + \mathrm{var}(z_j))\\) for a dissimilarity.

## Details

With \\z\\ the log-ratio-transformed abundances (clr or alr),
\$\$\rho\_{p,ij} = 1 - \frac{\mathrm{var}(z_i -
z_j)}{\mathrm{var}(z_i) + \mathrm{var}(z_j)} =
\frac{2\\\mathrm{cov}(z_i, z_j)}{\mathrm{var}(z_i) +
\mathrm{var}(z_j)}\$\$ where \\\mathrm{var}(z_i - z_j)\\ is
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md). It lies in
\\\[-1, 1\]\\: 1 for proportional features, 0 for none, -1 for
reciprocal ones (\\z_i = -z_j\\). It is Lin's concordance correlation
coefficient of \\z_i\\ and \\z_j\\. Unlike
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md), it depends
on the `denominator`. A pair of features that are both constant relative
to the reference has no defined \\\rho_p\\ (`NaN`).

This is the coefficient that propr computes with
`propr(counts, metric = "rho")` (see **Zeros** for its zero handling).
propr also estimates false discovery rates by permutation
(`updateCutoffs()`), and offers `ivar = "iqlr"`, a clr over the features
of mid-range variance; neither is offered here.

## Zeros

Logs of zero are undefined, so zeros need replacing first. With
`pseudocount = NULL` (the default), no pseudocount is used when `counts`
has no zeros; otherwise half the smallest non-zero value is added to
every value, with a warning. The recommended input is a table with its
zeros replaced by
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
which has none and is used as it is.

propr (Quinn et al. 2017) replaces only the zeros, by 1, which is not
`pseudocount = 1` (that adds 1 to every value). To reproduce it, write
`x[x == 0] <- 1` and then call these functions on `x`.

Pseudocounts are not neutral: in sparse tables the many zeros, all given
the same value, can dominate the variances. Filtering out rare features
first helps.

## Size

The result has \\D(D - 1)/2\\ values: 20,000 features give 2e8 pairs,
1.6 GB. Filter out rare features first for large tables. The table
itself stays sparse.

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

Erb, I., & Notredame, C. (2016). How should we measure proportionality
on relative gene expression data? *Theory in Biosciences*, 135, 21–36.

Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., &
Bähler, J. (2015). Proportionality: a valid alternative to correlation
for relative data. *PLoS Computational Biology*, 11, e1004075.

Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
propr: an R-package for identifying proportionally abundant features
using compositional data analysis. *Scientific Reports*, 7, 16252.

## See also

[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
(`'clr'`),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md).

## Examples

``` r
    x   <- replace_zeros(ex_counts, 'mult')
    rho <- rho_p(x)
    rho
#> # A dist: 6 features (lower triangle)
#> # Proportionality (rho_p)
#>                   Streptococcus Bacteroides Corynebacterium Haemophilus
#> Bacteroides              -0.504
#> Corynebacterium           0.085      -0.646
#> Haemophilus               0.865      -0.215          -0.319
#> Propionibacterium        -0.935       0.307           0.203      -0.983
#> Staphylococcus           -0.447      -0.412           0.652      -0.536
#> # … 1 more column: Propionibacterium
#> # Attributes: method

    # Pairs above a cutoff
    m <- as.matrix(rho)
    which(m > 0.5 & upper.tri(m), arr.ind = TRUE)
#>                   row col
#> Streptococcus       1   4
#> Corynebacterium     3   6
#> Propionibacterium   5   6

    # Relative to one feature (alr): that feature is left out
    rho_p(x, denominator = 'Streptococcus')
#> # A dist: 5 features (lower triangle)
#> # Proportionality (rho_p)
#>                   Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Corynebacterium         0.328
#> Haemophilus             0.231          -0.225
#> Propionibacterium       0.762           0.786      0.0277
#> Staphylococcus          0.483           0.903      0.0318             0.887
#> # Attributes: method
```
