# PhILR Distance

Euclidean distances between samples' PhILR coordinates, with the weights
that Silverman et al. (2017) recommend: a compositional, phylogenetic
beta diversity.

## Usage

``` r
philr_distance(
  counts,
  tree = NULL,
  margin = 1L,
  pseudocount = NULL,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of counts or proportions (samples \\\times\\
  features, see `margin`), or an object of one of the types in **Input
  Types**. Its features must have names, matched to the tips of `tree`.
  The part weights assume counts.

- tree:

  A `phylo` object with branch lengths, whose tips include every feature
  of `counts`, or `NULL` to use the tree inside `counts` (a phyloseq,
  rbiom or TreeSummarizedExperiment object, or `attr(counts, 'tree')`).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Added to every value before taking logs. `NULL` (the default) adds
  nothing when `counts` has no zeros, and otherwise uses half the
  smallest non-zero value, with a warning. See **Pseudocount**.

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
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A distance matrix: a `dist` object of class `c("ecodive_dist", "dist")`,
as every metric function returns (see
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)).

## Details

The distance between samples \\x\\ and \\y\\ is \$\$\sqrt{\sum_i
\left(x^\*\_i - y^\*\_i\right)^2},\$\$ where \\x^\*\\ are the
[`philr()`](https://cmmr.github.io/ecodive/reference/philr.md)
coordinates of \\x\\ with part weights `'enorm.x.gm.counts'` and ILR
weights `'blw.sqrt'` (Silverman et al. 2017, Methods), so that

    philr_distance(x, tree) ==
      euclidean(philr(x, tree, part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt'))

For other weightings, use `euclidean(philr(x, tree, ...))`. With uniform
weights, PhILR distances are Aitchison distances
([`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)):
the weights are what make them phylogenetic. The part weights are
computed from all the samples, so adding or removing a sample changes
the distances between the others. `pairs` limits which distances are
computed, not the samples used.

The tree is handled as
[`philr()`](https://cmmr.github.io/ecodive/reference/philr.md) handles
it: tips that aren't features are dropped, and polytomies are resolved,
with a warning (see **The tree** in
[`philr()`](https://cmmr.github.io/ecodive/reference/philr.md)). It
needs branch lengths.

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

Logs of zero are undefined, so zeros need replacing first. The
pseudocount is added to every value (not only to zeros), as in
`normalize(norm = 'clr')`. With `pseudocount = NULL` (the default), no
pseudocount is used when `counts` has no zeros; otherwise half the
smallest non-zero value is, with a warning. To avoid a pseudocount,
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md):
its result has none.

The choice of pseudocount is not neutral: it can distort the results,
especially for sparse tables. See
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## References

Silverman, J. D., Washburne, A. D., Mukherjee, S., & David, L. A.
(2017). A phylogenetic transform enhances analysis of compositional
microbiota data. *eLife*, 6, e21887.
[doi:10.7554/eLife.21887](https://doi.org/10.7554/eLife.21887)

## See also

[`philr()`](https://cmmr.github.io/ecodive/reference/philr.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md),
[`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md),
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md),
[`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    philr_distance(ex_counts, tree = ex_tree, pseudocount = 1)
#> # A dist: 4 samples (lower triangle)
#> # PhILR Distance
#>       Saliva Gums Nose
#> Gums    45.4
#> Nose   175    173
#> Stool  281    263  183
#> # Attributes: method

    # The same, from beta_div()
    d <- beta_div(ex_counts, 'philr', tree = ex_tree, pseudocount = 1)
```
