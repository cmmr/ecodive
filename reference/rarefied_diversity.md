# Diversity Averaged over Repeated Rarefactions

Rarefy a count table many times, compute alpha or beta diversity from
each rarefied table, and average the results. A single rarefaction
throws reads away at random, so its diversity values carry sampling
noise; averaging many keeps rarefaction's control for uneven sequencing
depth and removes most of that noise. This is mothur's
`dist.shared(subsample = T)`, vegan's `avgdist()`, and QIIME 2's
q2-boots.

## Usage

``` r
beta_div_rarefied(
  counts,
  metric,
  depth = NULL,
  times = 100L,
  average = "mean",
  sd = FALSE,
  margin = 1L,
  norm = "none",
  pseudocount = NULL,
  power = 1.5,
  alpha = 0.5,
  tree = NULL,
  pairs = NULL,
  clust = NULL,
  seed = 0L,
  cpus = n_cpus()
)

alpha_div_rarefied(
  counts,
  metric,
  depth = NULL,
  times = 100L,
  average = "mean",
  sd = FALSE,
  cutoff = 10L,
  digits = 3L,
  tree = NULL,
  margin = 1L,
  q = 2,
  bias_corrected = TRUE,
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A count table (samples \\\times\\ features) of whole numbers; see
  **Input Types**. A `phyloseq`, `rbiom` or `TreeSummarizedExperiment`
  object's tree is used for the phylogenetic metrics.

- metric:

  The name of a diversity metric: for `beta_div_rarefied()`, one of
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)'s
  metrics, such as `'bray'`; for `alpha_div_rarefied()`, one of
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)'s,
  such as `'shannon'`. Matched as those functions match it
  (case-insensitive, with prefixes and aliases). See
  [`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md).

- depth:

  The number of observations to rarefy every sample to: a single whole
  number of at least 1. Samples with fewer observations are dropped
  (with a warning), and are listed in the result's `dropped` attribute.
  `NULL` (the default) chooses the depth as
  [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) does
  (see "Choosing the depth").

- times:

  The number of rarefactions to average, a whole number of at least 1.
  Default: `100`, as `vegan::avgdist()`.

- average:

  How the rarefactions' values are combined: `"mean"` (the default),
  `"median"`, or for `beta_div_rarefied()` only, `"medoid"`.
  Case-insensitive, and a prefix is enough. See "Choosing how to
  average".

- sd:

  Also return the standard deviation of each value over the `times`
  rarefactions, as attribute `sd`. Default: `FALSE`.

- margin:

  The margin containing samples: `1` if samples are rows, `2` if samples
  are columns. Ignored for `phyloseq`, `rbiom` and
  `SummarizedExperiment` objects. Default: `1`.

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

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- power:

  Only used when `metric = 'minkowski'`. Scaling factor for the
  magnitude of differences between communities (\\p\\). Default: `1.5`

- alpha:

  Only used when `metric = 'generalized_unifrac'`. How much weight to
  give to relative abundances; a value between 0 and 1, inclusive.
  Setting `alpha=1` is equivalent to
  [`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md).

- tree:

  Only used by phylogeny-aware metrics. A `phylo`-class object
  representing the phylogenetic tree for the OTUs in `counts`. The OTU
  identifiers given by `colnames(counts)` must be present in `tree`. Can
  be omitted if a tree is embedded with the `counts` object or as
  `attr(counts, 'tree')`.

- pairs:

  Which pairs of samples to compute distances for, as in
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md):
  indexes into the distance matrix of the samples that are **kept** (not
  dropped for having fewer than `depth` observations), as with
  `beta_div(rarefy(counts, depth), metric, pairs = pairs)`. Default:
  `NULL` (all pairs).

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

- seed:

  An integer seed for the rarefactions: iteration \\t\\ uses
  `seed + t - 1`, as `rarefy(times = )` does. The same seed gives the
  same result, for any `cpus`. Default: `0`.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

- cutoff:

  The maximum number of observations to consider "rare". Default: `10`.

- digits:

  Precision of the returned values, in number of decimal places. E.g.
  the default `digits=3` could return `6.392`.

- q:

  The order of a Hill number, or of Rényi's or Tsallis's entropy: a
  single number of at least 0, where `Inf` is allowed. Only used by
  [`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
  [`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md).
  Default: `2`, as scikit-bio (which calls it `order`).

- bias_corrected:

  For Chao1 only: use the bias-corrected form (the default, as
  scikit-bio and QIIME 2), or the classic form when a sample has both
  singletons and doubletons (`FALSE`). See
  [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md).

## Value

`beta_div_rarefied()`: a distance matrix of class
`c("ecodive_dist", "dist")` over the kept samples, as
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
returns, with attribute `method` (the metric's ID). Its title (in its
`display` attribute; see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md))
names the metric and the rarefactions, such as
`"Bray-Curtis Dissimilarity (rarefied 100×)"`, so ordinations and tests
of it are titled after it. Further attributes:

- `depth`: the depth every sample was rarefied to.

- `dropped`: the names (or, without names, the positions) of the samples
  dropped for having fewer than `depth` observations; empty if none.

- `sd` (with `sd = TRUE`): the standard deviations, a distance matrix
  over the same samples, in the same order.

- `medoid` (with `average = "medoid"`): the iteration \\t\\ whose
  distance matrix was returned.

- `order` (with `clust`): see
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).

`alpha_div_rarefied()`: a named numeric vector over the kept samples, as
[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
returns, with attributes `depth`, `dropped`, and with `sd = TRUE`, `sd`
(a named numeric vector).

## Details

- `beta_div_rarefied()`: a distance matrix, as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).

- `alpha_div_rarefied()`: one value per sample, as from
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md).

Iteration \\t\\ (\\t = 1, \ldots,\\ `times`) rarefies `counts` with
[`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md), with
seed `seed + t - 1`, and computes the metric from it, exactly as
`beta_div(rarefy(counts, depth, seed = seed + t - 1), metric)` would. So
`times = 1` gives the same result as that call.

## When to use it

Rarefying controls for uneven sequencing depth: every sample is compared
at the same number of reads, so that differences in depth don't pass for
differences in diversity. McMurdie & Holmes (2014) called rarefying
"inadmissible" because a single rarefaction discards data and adds
noise; averaging many rarefactions answers the second objection. Schloss
(2024) found rarefaction better than the alternatives at controlling for
uneven sequencing effort in amplicon data, and Schloss (2026) found it
better than robust Aitchison PCA
([`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md)) and other
compositional methods. To compare both views on your own data, ordinate
both and test their agreement (see the examples).

## Choosing the depth

Every sample with fewer than `depth` observations is dropped, the same
samples from every iteration. A higher depth keeps more of each sample's
information; a lower one keeps more samples. With `depth = NULL`, the
depth is the smallest library size at which the observations kept (that
depth times the number of samples with at least that many) make up at
least 10% of all the observations in the table, as in
[`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md). Check
`attr(x, "dropped")`, and choose `depth` yourself when the samples
dropped matter.

## Choosing how to average

- `"mean"` (the default, as in `vegan::avgdist()` and mothur): the
  average of each value over the rarefactions. It needs one pass over
  them, and memory for a few distance matrices.

- `"median"`: the median of each value; less affected by an unusual
  rarefaction.

- `"medoid"` (beta only): the distance matrix of the one rarefaction
  closest to all the others (smallest sum of Euclidean distances between
  the iterations' distances). Unlike the mean and the median, the result
  is the distance matrix of one real rarefied table, so it keeps
  whatever properties the metric has (such as being Euclidean), which a
  cell-by-cell average can lose.

q2-boots' defaults are `average = "medoid"` for beta diversity and
`average = "median"` for alpha diversity; set them to get the same.

**Memory and time.** The mean rarefies `times` times. The median and the
medoid need every rarefaction's value of each pair of samples at once,
so for many samples they work through the pairs in blocks that fit in
256 MB (`8 * times` bytes per pair), rarefying the table `times` times
for each block, with the same seeds and so the same draws. With
`times = 100`, 500 samples (124,750 pairs) fit in one block, 1,000
samples need two, and 10,000 samples about 150. The medoid also needs a
`times` \\\times\\ `times` matrix, and rarefies once more at the end,
for the medoid's matrix.

## Is an ordination driven by depth?

[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) fits
library size onto any ordination, with a permutation p-value:
`envfit(ord, data.frame(log_depth = log(rowSums(counts))))`. A large
\\r^2\\ means that the ordination still follows sequencing depth.
`cor(ord, log(rowSums(counts)))` gives the correlation with each axis.
See the examples.

## Warnings

A warning raised by the metric (such as the Aitchison distance's choice
of pseudocount) is given once, after all the rarefactions, however many
of them raised it.

## Same names in other packages

None, but `vegan::avgdist()` computes the same as `beta_div_rarefied()`
with `average = "mean"` (or `"median"`), from R's random number
generator, so their results agree only up to Monte Carlo error.

## Differences from q2-boots

- No bootstrap mode:
  [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) draws
  only without replacement, so there is no equivalent of q2-boots'
  `replacement = True`.

- The default `average` is the mean, for alpha and beta diversity alike
  (q2-boots: the median for alpha, the medoid for beta).

- The rarefactions are seeded with ecodive's own PCG generator, so the
  draws differ from q2-boots' (and from vegan's) for the same seed.

## References

Raspet, I., Gehret, E., Herman, C., Meilander, J., Manley, A., Simard,
A., Bolyen, E., & Caporaso, J. G. (2025). Facilitating bootstrapped and
rarefaction-based microbiome diversity analysis with q2-boots.
*F1000Research*, 14, 87.
[doi:10.12688/f1000research.156295.1](https://doi.org/10.12688/f1000research.156295.1)

Schloss, P. D. (2024). Rarefaction is currently the best approach to
control for uneven sequencing effort in amplicon sequence analyses.
*mSphere*, 9, e00354-23.
[doi:10.1128/msphere.00354-23](https://doi.org/10.1128/msphere.00354-23)

Schloss, P. D. (2026). Rarefaction is better than robust Aitchison PCA
and other compositional data analysis methods at controlling for uneven
sequencing effort. *mSystems*, 11, e00930-26.
[doi:10.1128/msystems.00930-26](https://doi.org/10.1128/msystems.00930-26)

McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why rarefying
microbiome data is inadmissible. *PLoS Computational Biology*, 10,
e1003531.
[doi:10.1371/journal.pcbi.1003531](https://doi.org/10.1371/journal.pcbi.1003531)

## See also

[`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), and
[`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md)
to choose a depth.

## Examples

``` r
    # Bray-Curtis dissimilarity, averaged over 100 rarefactions to the
    # automatically chosen depth.
    dm <- beta_div_rarefied(ex_counts, 'bray')
    dm
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Dissimilarity (rarefied 100×)
#>       Saliva  Gums  Nose
#> Gums   0.426
#> Nose   0.972 0.973
#> Stool  0.991 0.991 0.992
#> # Attributes: method, depth, dropped
    attr(dm, 'depth')
#> [1] 345

    # With the spread over the rarefactions.
    dm <- beta_div_rarefied(ex_counts, 'bray', depth = 100, sd = TRUE)
    attr(dm, 'sd')
#> # A dist: 4 samples (lower triangle)
#> # SD of Bray-Curtis Dissimilarity (rarefied 100×)
#>        Saliva    Gums    Nose
#> Gums  0.0517
#> Nose  0.0163  0.0155
#> Stool 0.00795 0.00842 0.00757

    # q2-boots' defaults: the medoid for beta, the median for alpha.
    beta_div_rarefied(ex_counts, 'bray', depth = 100, average = 'medoid')
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Dissimilarity (rarefied 100×)
#>       Saliva Gums Nose
#> Gums    0.42
#> Nose    0.97 0.98
#> Stool   0.99 1    0.99
#> # Attributes: method, depth, dropped, medoid
    alpha_div_rarefied(ex_counts, 'shannon', depth = 100, average = 'median')
#>     Saliva       Gums       Nose      Stool 
#> 0.73365719 0.34836336 1.11804259 0.05600153 
#> attr(,"depth")
#> [1] 100
#> attr(,"dropped")
#> character(0)

    # Ordinations and tests are titled after the rarefied matrix.
    ord <- pcoa(dm, n_axes = 2)
    ord
#> # A matrix: 4 samples × 2 axes
#> # Principal Coordinates Analysis of Bray-Curtis Dissimilarity (rarefied 100×)
#> # % of variance: 53.5, 39.0
#>         PCoA1   PCoA2
#> Saliva -0.408 -0.0241
#> Gums   -0.41  -0.0244
#> Nose    0.369  0.519
#> Stool   0.449 -0.47
#> # Attributes: method, eig, prop, trace, correction, constant

    # Does the ordination still follow sequencing depth?
    kept <- attr(dm, 'Labels')
    envfit(ord, data.frame(log_depth = log(rowSums(ex_counts)[kept])))
#> # A data.frame: 1 × 7
#> # Fitted variables on Principal Coordinates Analysis of Bray-Curtis Dissimilarity (rarefied 100×)
#>           variable  level PCoA1 PCoA2    r2 p_value     z
#>           <chr>     <chr> <dbl> <dbl> <dbl>   <dbl> <dbl>
#> log_depth log_depth NA      0.6   0.8 0.358    0.83 -1.02
#> # Attributes: r2_perm, n_perm
```
