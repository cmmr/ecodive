# Beta Diversity Wrapper Function

Beta Diversity Wrapper Function

## Usage

``` r
beta_div(
  counts,
  metric,
  norm = "percent",
  power = 1.5,
  pseudocount = NULL,
  alpha = 0.5,
  tree = NULL,
  pairs = NULL,
  margin = 1L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- metric:

  The name of a beta diversity metric. One of
  `c('aitchison', 'bhattacharyya', 'bray', 'canberra', 'chebyshev', 'chord', 'clark', 'divergence', 'euclidean', 'generalized_unifrac', 'gower', 'hamming', 'hellinger', 'horn', 'jaccard', 'jensen', 'jsd', 'lorentzian', 'manhattan', 'matusita', 'minkowski', 'morisita', 'motyka', 'normalized_unifrac', 'ochiai', 'psym_chisq', 'soergel', 'sorensen', 'squared_chisq', 'squared_chord', 'squared_euclidean', 'topsoe', 'unweighted_unifrac', 'variance_adjusted_unifrac', 'wave_hedges', 'weighted_unifrac')`.
  Flexible matching is supported (see below). Programmatic access via
  `list_metrics('beta')`.

- norm:

  Normalize the incoming counts. Options are: \* `'percent'`: Relative
  abundance (sample abundances sum to 1). \* `'binary'`: Unweighted
  presence/absence (each count is either 0 or 1). \* `'clr'`: Centered
  log ratio. \* `'none'`: No transformation.

         Default: `'percent'`, which is the expected input for these formulas.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

- pseudocount:

  The value to add to all counts in `counts` to prevent taking `log(0)`
  for unobserved features. The default, `NULL`, selects the smallest
  non-zero value in `counts`.

- alpha:

  How much weight to give to relative abundances; a value between 0 and
  1, inclusive. Setting `alpha=1` is equivalent to
  [`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md).

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A numeric vector.

## Details

**List of Beta Diversity Metrics**

|                             |                                                  |
|-----------------------------|--------------------------------------------------|
| Option / Function Name      | Metric Name                                      |
| `aitchison`                 | Aitchison distance                               |
| `bhattacharyya`             | Bhattacharyya distance                           |
| `bray`                      | Bray-Curtis dissimilarity                        |
| `canberra`                  | Canberra distance                                |
| `chebyshev`                 | Chebyshev distance                               |
| `chord`                     | Chord distance                                   |
| `clark`                     | Clark's divergence distance                      |
| `divergence`                | Divergence                                       |
| `euclidean`                 | Euclidean distance                               |
| `generalized_unifrac`       | Generalized UniFrac (GUniFrac)                   |
| `gower`                     | Gower distance                                   |
| `hamming`                   | Hamming distance                                 |
| `hellinger`                 | Hellinger distance                               |
| `horn`                      | Horn-Morisita dissimilarity                      |
| `jaccard`                   | Jaccard distance                                 |
| `jensen`                    | Jensen-Shannon distance                          |
| `jsd`                       | Jesen-Shannon divergence (JSD)                   |
| `lorentzian`                | Lorentzian distance                              |
| `manhattan`                 | Manhattan distance                               |
| `matusita`                  | Matusita distance                                |
| `minkowski`                 | Minkowski distance                               |
| `morisita`                  | Morisita dissimilarity                           |
| `motyka`                    | Motyka dissimilarity                             |
| `normalized_unifrac`        | Normalized Weighted UniFrac                      |
| `ochiai`                    | Otsuka-Ochiai dissimilarity                      |
| `psym_chisq`                | Probabilistic Symmetric Chi-Squared distance     |
| `soergel`                   | Soergel distance                                 |
| `sorensen`                  | Dice-Sorensen dissimilarity                      |
| `squared_chisq`             | Squared Chi-Squared distance                     |
| `squared_chord`             | Squared Chord distance                           |
| `squared_euclidean`         | Squared Euclidean distance                       |
| `topsoe`                    | Topsoe distance                                  |
| `unweighted_unifrac`        | Unweighted UniFrac                               |
| `variance_adjusted_unifrac` | Variance-Adjusted Weighted UniFrac (VAW-UniFrac) |
| `wave_hedges`               | Wave Hedges distance                             |
| `weighted_unifrac`          | Weighted UniFrac                                 |

**Flexible name matching**

Case insensitive and partial matching. Any runs of non-alpha characters
are converted to underscores. E.g. `metric = 'Weighted UniFrac` selects
`weighted_unifrac`.

UniFrac names can be shortened to the first letter plus "unifrac". E.g.
`uunifrac`, `w_unifrac`, or `V UniFrac`. These also support partial
matching.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the metrics list.

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
    beta_div(ex_counts, 'bray')
#>          Saliva      Gums      Nose
#> Gums  0.4265973                    
#> Nose  0.9713843 0.9720256          
#> Stool 0.9909509 0.9911046 0.9915177
    
    # Generalized UniFrac distances
    beta_div(ex_counts, 'GUniFrac', tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4471644                    
#> Nose  0.8215129 0.7607876          
#> Stool 0.9727827 0.9784242 0.9730332
    
```
