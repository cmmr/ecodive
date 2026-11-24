# Beta Diversity Wrapper Function

Beta Diversity Wrapper Function

## Usage

``` r
beta_div(
  counts,
  metric,
  margin = 1L,
  norm = "none",
  pseudocount = NULL,
  power = 1.5,
  alpha = 0.5,
  tree = NULL,
  pairs = NULL,
  clust = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted by some diversity metrics.

- metric:

  The name of a beta diversity metric. One of
  `c('aitchison', 'bhattacharyya', 'bray', 'bray_balanced', 'bray_gradient', 'canberra', 'chebyshev', 'chord', 'clark', 'divergence', 'euclidean', 'generalized_unifrac', 'gower', 'hamming', 'hellinger', 'horn', 'jaccard', 'jaccard_nestedness', 'jaccard_turnover', 'jensen', 'jsd', 'lorentzian', 'manhattan', 'matusita', 'minkowski', 'morisita', 'motyka', 'normalized_unifrac', 'ochiai', 'philr_distance', 'psym_chisq', 'robust_aitchison', 'soergel', 'soergel_balanced', 'soergel_gradient', 'sorensen', 'sorensen_nestedness', 'sorensen_turnover', 'squared_chisq', 'squared_chord', 'squared_euclidean', 'topsoe', 'unweighted_unifrac', 'variance_adjusted_unifrac', 'wave_hedges', 'weighted_unifrac')`.
  Flexible matching is supported (see below). Programmatic access via
  `list_metrics('beta')`.

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

## Value

A distance matrix: a `dist` object with one distance per pair of samples
(the packed lower triangle), labelled with the sample names, as
[`stats::dist()`](https://rdrr.io/r/stats/dist.html) returns. Its class
is `c("ecodive_dist", "dist")`, so that
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it as a
heatmap (see
[`plot.ecodive_dist()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_dist.md));
every function for `dist` objects, such as
[`as.matrix()`](https://rdrr.io/r/base/matrix.html), works on it as
before. Each metric's own function, such as
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), returns
the same.

Besides `dist`'s `Size`, `Labels`, `Diag` and `Upper`, it has attributes
`method`, the metric's ID (such as `"bray"`), as
[`stats::dist()`](https://rdrr.io/r/stats/dist.html) sets it; with
`clust`, `order` (see `clust`); and `display`, for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`print()`](https://rdrr.io/r/base/print.html): its title, the metric's
display name (such as `"Bray-Curtis Dissimilarity"`), which ordinations
and tests of the matrix name in theirs, and with `clust`, the clustering
method (see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)).

## Details

**List of Beta Diversity Metrics**

|  |  |
|----|----|
| Option / Function Name | Metric Name |
| `aitchison` | Aitchison distance |
| `bhattacharyya` | Bhattacharyya distance |
| `bray` | Bray-Curtis dissimilarity |
| `bray_balanced` | Bray-Curtis balanced variation |
| `bray_gradient` | Bray-Curtis abundance gradient |
| `canberra` | Canberra distance |
| `chebyshev` | Chebyshev distance |
| `chord` | Chord distance |
| `clark` | Clark's divergence distance |
| `divergence` | Divergence |
| `euclidean` | Euclidean distance |
| `generalized_unifrac` | Generalized UniFrac (GUniFrac) |
| `gower` | Gower distance |
| `hamming` | Hamming distance |
| `hellinger` | Hellinger distance |
| `horn` | Horn-Morisita dissimilarity |
| `jaccard` | Jaccard distance |
| `jaccard_nestedness` | Jaccard nestedness |
| `jaccard_turnover` | Jaccard turnover |
| `jensen` | Jensen-Shannon distance |
| `jsd` | Jesen-Shannon divergence (JSD) |
| `lorentzian` | Lorentzian distance |
| `manhattan` | Manhattan distance |
| `matusita` | Matusita distance |
| `minkowski` | Minkowski distance |
| `morisita` | Morisita dissimilarity |
| `motyka` | Motyka dissimilarity |
| `normalized_unifrac` | Normalized Weighted UniFrac |
| `ochiai` | Otsuka-Ochiai dissimilarity |
| `philr_distance` | PhILR distance |
| `psym_chisq` | Probabilistic Symmetric Chi-Squared distance |
| `soergel` | Soergel distance |
| `soergel_balanced` | Ruzicka (Soergel) balanced variation |
| `soergel_gradient` | Ruzicka (Soergel) abundance gradient |
| `sorensen` | Dice-Sorensen dissimilarity |
| `sorensen_nestedness` | Sorensen nestedness |
| `sorensen_turnover` | Simpson dissimilarity (Sorensen turnover) |
| `squared_chisq` | Squared Chi-Squared distance |
| `squared_chord` | Squared Chord distance |
| `squared_euclidean` | Squared Euclidean distance |
| `topsoe` | Topsoe distance |
| `unweighted_unifrac` | Unweighted UniFrac |
| `variance_adjusted_unifrac` | Variance-Adjusted Weighted UniFrac (VAW-UniFrac) |
| `wave_hedges` | Wave Hedges distance |
| `weighted_unifrac` | Weighted UniFrac |

**Flexible name matching**

Case insensitive and partial matching. Any runs of non-alpha characters
are converted to underscores. E.g. `metric = 'Weighted UniFrac` selects
`weighted_unifrac`.

UniFrac names can be shortened to the first letter plus "unifrac". E.g.
`uunifrac`, `w_unifrac`, or `V UniFrac`. These also support partial
matching.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the metrics list.

The turnover / nestedness and balanced / gradient components are
Baselga's partition of `sorensen`, `jaccard`, `bray` and `soergel`; see
[`sorensen_turnover()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md),
[`bray_balanced()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md)
and
[`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md).

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

## Examples

``` r
    # Example counts matrix
    ex_counts
#>        Streptococcus Bacteroides Corynebacterium Haemophilus Propionibacterium
#> Saliva           162           2               0         180                 1
#> Gums             793           4               0          87                 1
#> Nose              22           2             498           2               251
#> Stool              1         611               1           1                 0
#>        Staphylococcus
#> Saliva              0
#> Gums                1
#> Nose              236
#> Stool               1
    
    # Bray-Curtis distances
    dm <- beta_div(ex_counts, 'bray')
    dm
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Dissimilarity
#>       Saliva  Gums  Nose
#> Gums   0.591
#> Nose   0.96  0.97
#> Stool  0.992 0.991 0.993
#> # Attributes: method

    # As a heatmap.
    plot(dm)

    
    # Generalized UniFrac distances
    beta_div(ex_counts, 'GUniFrac', tree = ex_tree)
#> # A dist: 4 samples (lower triangle)
#> # Generalized UniFrac (GUniFrac)
#>       Saliva  Gums  Nose
#> Gums   0.447
#> Nose   0.822 0.761
#> Stool  0.973 0.978 0.973
#> # Attributes: method
    
```
