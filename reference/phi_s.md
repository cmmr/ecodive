# Symmetric Proportionality

The symmetric proportionality measure \\\phi_s\\ of Quinn et al. (2017),
propr's `"phs"`: the variance of the log ratio of two features relative
to that of their log product. `0` means perfectly proportional.

## Usage

``` r
phi_s(
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
  every distance, so it cannot be used with `pairs`.
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md)'s
  similarities are clustered as \\1 - \rho_p\\.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A `dist` object over the features, as
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) returns,
with `method` `"phi_s"`. It is a dissimilarity, so it can be given as
`dm` to [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) or
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

## Details

With \\z\\ the log-ratio-transformed abundances (clr or alr),
\$\$\phi\_{s,ij} = \frac{\mathrm{var}(z_i - z_j)}{\mathrm{var}(z_i +
z_j)} = \frac{1 - \rho\_{p,ij}}{1 + \rho\_{p,ij}}\$\$ where
\\\mathrm{var}(z_i - z_j)\\ is
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) and
\\\rho_p\\ is
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md). It lies
in \\\[0, \infty)\\: 0 for proportional features, 1 for unrelated ones,
and infinite for reciprocal ones (\\z_i = -z_j\\). It is computed as
\\\mathrm{vlr}\_{ij} / (2\mathrm{var}(z_i) + 2\mathrm{var}(z_j) -
\mathrm{vlr}\_{ij})\\, since \\\mathrm{var}(a + b) = 2\mathrm{var}(a) +
2\mathrm{var}(b) - \mathrm{var}(a - b)\\. Unlike
[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md), it depends
on the `denominator`. Near-reciprocal pairs (such as the two features of
a two-feature `denominator`, whose \\z\\ are exact opposites) have
\\\mathrm{var}(z_i + z_j)\\ near 0, where rounding dominates: their
value is `Inf` or merely very large.

Lovell et al.'s (2015) original \\\phi = \mathrm{vlr}\_{ij} /
\mathrm{var}(z_i)\\ is not symmetric (\\\phi\_{ij} \ne \phi\_{ji}\\), so
it can't be a `dist`; Quinn et al. (2017) replaced it with \\\phi_s\\.
This is the measure that propr computes with
`propr(counts, metric = "phs")` (see **Zeros** for its zero handling).

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

Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marguerat, S., &
Bähler, J. (2015). Proportionality: a valid alternative to correlation
for relative data. *PLoS Computational Biology*, 11, e1004075.

Quinn, T. P., Richardson, M. F., Lovell, D., & Crowley, T. M. (2017).
propr: an R-package for identifying proportionally abundant features
using compositional data analysis. *Scientific Reports*, 7, 16252.

## See also

[`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
[`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md).

## Examples

``` r
    x <- replace_zeros(ex_counts, 'mult')
    phi_s(x)
#> # A dist: 6 features (lower triangle)
#> # Proportionality (phi_s)
#>                   Streptococcus Bacteroides Corynebacterium Haemophilus
#> Bacteroides              3.03
#> Corynebacterium          0.843         4.65
#> Haemophilus              0.0725        1.55           1.94
#> Propionibacterium       29.7           0.53           0.662      119
#> Staphylococcus           2.62          2.4            0.211        3.31
#> # … 1 more column: Propionibacterium
#> # Attributes: method

    # The same as (1 - rho_p) / (1 + rho_p)
    rho <- rho_p(x)
    all.equal(c(phi_s(x)), c((1 - rho) / (1 + rho)))
#> [1] TRUE
```
