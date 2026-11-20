# documentation

documentation

## Arguments

- alpha:

  How much weight to give to relative abundances; a value between 0 and
  1, inclusive. Setting `alpha=1` is equivalent to
  [`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md).

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

- cutoff:

  The maximum number of observations to consider "rare". Default: `10`.

- digits:

  Precision of the returned values, in number of decimal places. E.g.
  the default `digits=3` could return `6.392`.

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

  Default: `'none'`.

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

- pairwise:

  Also test each pair of groups: run the same test on just the two
  groups' samples, for every pair, and adjust the pairs' p-values for
  multiple testing with `p_adjust`. See **Pairwise comparisons**.
  Default: `FALSE`.

- p_adjust:

  With `pairwise = TRUE`, how to adjust the pairs' p-values for multiple
  testing: one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html)
  (`'holm'`, `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`,
  `'fdr'`, or `'none'`), matched as
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) matches
  them. Giving it without `pairwise = TRUE` is an error. Default:
  `'holm'`.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. See **Pseudocount** section.

- strata:

  Permute samples only within these groups (blocks), for designs such as
  repeated measures. With `data`, the name of one of its columns, or a
  vector with one value per row of `data`. Without `data`, a vector with
  one value per sample, matched to the samples by name when it has
  names, otherwise in the same order. Default: `NULL` (free
  permutation).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

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
