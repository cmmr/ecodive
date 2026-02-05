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
  `c('aitchison', 'bhattacharyya', 'bray', 'canberra', 'chebyshev', 'chord', 'clark', 'divergence', 'euclidean', 'generalized_unifrac', 'gower', 'hamming', 'hellinger', 'horn', 'jaccard', 'jensen', 'jsd', 'lorentzian', 'manhattan', 'matusita', 'minkowski', 'morisita', 'motyka', 'normalized_unifrac', 'ochiai', 'psym_chisq', 'soergel', 'sorensen', 'squared_chisq', 'squared_chord', 'squared_euclidean', 'topsoe', 'unweighted_unifrac', 'variance_adjusted_unifrac', 'wave_hedges', 'weighted_unifrac')`.
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

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. See **Pseudocount** section.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

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

The `pseudocount` parameter is only relevant when `norm = 'clr'`.

Zeros are undefined in the centered log-ratio (CLR) transformation. If
`norm = 'clr'`, `pseudocount` is `NULL` (the default), and zeros are
detected, the function uses half the minimum non-zero value
(`min(x[x>0]) / 2`) and issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`).

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
    beta_div(ex_counts, 'bray')
#>          Saliva      Gums      Nose
#> Gums  0.5905768                    
#> Nose  0.9601770 0.9704797          
#> Stool 0.9916667 0.9906729 0.9926199
    
    # Generalized UniFrac distances
    beta_div(ex_counts, 'GUniFrac', tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4471644                    
#> Nose  0.8215129 0.7607876          
#> Stool 0.9727827 0.9784242 0.9730332
    
```
