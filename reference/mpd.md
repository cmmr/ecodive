# Mean Pairwise and Mean Nearest Taxon Distance

How closely related the taxa in each sample are, along a phylogeny (or
any dissimilarity among the taxa): the mean phylogenetic distance
between pairs of the sample's taxa (MPD), or from each taxon to its
nearest relative in the sample (MNTD) (Webb 2000; Webb et al. 2002). MPD
responds to structure deep in the tree, MNTD to structure near its tips.
The `weighted_*` forms weight each taxon by its abundance. To test them
against a null model (NRI and NTI), see
[`mpd_ses()`](https://cmmr.github.io/ecodive/reference/mpd_ses.md).

## Usage

``` r
mpd(counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus())

weighted_mpd(counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus())

mntd(counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus())

weighted_mntd(
  counts,
  tree = NULL,
  taxa_dm = NULL,
  margin = 1L,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\ taxa,
  see `margin`), or any of the other objects listed in "Input Types".
  Values need not be integers: `mpd()` and `mntd()` use only presence,
  and the weighted forms only relative abundances. The taxa must have
  names, to match them to the tips of `tree` or the labels of `taxa_dm`.

- tree:

  A `phylo` object with branch lengths (none negative or missing), whose
  tips include every taxon of `counts` (see
  [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)).
  Tips that aren't taxa of `counts` are ignored. Can be omitted if a
  tree is embedded with the `counts` object or as
  `attr(counts, 'tree')`. Distances are path lengths along the tree
  (patristic distances), as from
  [`ape::cophenetic.phylo()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html).
  The tree may be unrooted, or have polytomies or nodes with one child:
  none of these change path lengths.

- taxa_dm:

  Instead of `tree`: the dissimilarities among the taxa, as a `dist`
  object or symmetric matrix (picante's `dis`, which may be
  functional-trait distances, for instance). Labels are matched to the
  taxa of `counts`; without labels, rows are taken in the order of
  `counts`' taxa. An embedded tree is ignored when `taxa_dm` is given.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A numeric vector with one value per sample, named by sample. A sample
with fewer than two taxa is `NA`.

## Details

For a sample with taxa \\i = 1 \ldots S\\, abundances \\x_i\\ (total
\\W\\), and distances \\d\_{ij}\\:

- **MPD:** \\\sum\_{i\<j} d\_{ij} / (S(S-1)/2)\\, the mean over pairs of
  distinct taxa.

- **Weighted MPD:** \\\sum\_{i,j} x_i x_j d\_{ij} / W^2\\, the mean
  distance between two individuals drawn at random, with replacement, so
  that pairs within a taxon count too (at distance 0). This is picante's
  and phylocom's definition.

- **MNTD:** \\\sum_i m_i / S\\, with \\m_i = \min\_{j \ne i} d\_{ij}\\.

- **Weighted MNTD:** \\\sum_i x_i m_i / W\\.

The weighted forms are unchanged by scaling a sample, so counts and
proportions give the same values. "Weighted" here means weighting by
abundance; it is not the weighting of
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md).

These are also alpha diversity metrics:
`alpha_div(counts, 'mpd', tree = tree)`, and likewise `'weighted_mpd'`,
`'mntd'` and `'weighted_mntd'`.

**Speed and memory**

With a tree, no taxa \\\times\\ taxa matrix is formed. Each sample's
taxa are sorted into the tree's preorder, and the tree's lowest common
ancestors of neighbours (found in constant time from a table built once)
give the part of the tree that spans them, with at most \\2S - 1\\
nodes; MPD and MNTD are one or two passes over it. Each sample takes
time proportional to \\S \log S\\, however large the tree. With
`taxa_dm`, each sample reads its \\S^2\\ distances. Results are
identical for any number of `cpus`.

**Differences from picante**

Values agree with picante 1.8.4's `mpd()` and `mntd()` (with
`dis = cophenetic(tree)` and `abundance.weighted` set to match) to about
\\10^{-12}\\. picante computes `cophenetic(tree)` first, a dense matrix
of every pair of tips (3.2 GB for 20,000 tips); `mpd()` reads the tree
directly. phylocom's `comstruct` gives MNTD 0, not NA, for a sample with
one taxon.

**Same names in other packages**

`picante::mpd()` and `picante::mntd()` share these names (taking `samp`,
`dis` and `abundance.weighted`); whichever package was attached last
masks the other's. Write `ecodive::mpd()` to be sure.
`PhyloMeasures::mpd.query()` computes the unweighted MPD.

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

Webb, C. O. (2000). Exploring the phylogenetic structure of ecological
communities: an example for rain forest trees. *American Naturalist*,
156(2), 145-155. [doi:10.1086/303378](https://doi.org/10.1086/303378)

Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
Phylogenies and community ecology. *Annual Review of Ecology and
Systematics*, 33, 475-505.
[doi:10.1146/annurev.ecolsys.33.010802.150448](https://doi.org/10.1146/annurev.ecolsys.33.010802.150448)

Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H.,
Ackerly, D. D., Blomberg, S. P., & Webb, C. O. (2010). Picante: R tools
for integrating phylogenies and ecology. *Bioinformatics*, 26(11),
1463-1464.
[doi:10.1093/bioinformatics/btq166](https://doi.org/10.1093/bioinformatics/btq166)

## See also

[`mpd_ses()`](https://cmmr.github.io/ecodive/reference/mpd_ses.md) for
NRI and NTI;
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md);
[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md);
[`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md).

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md),
[`mpd_ses()`](https://cmmr.github.io/ecodive/reference/mpd_ses.md),
[`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md),
[`philr_distance()`](https://cmmr.github.io/ecodive/reference/philr_distance.md),
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md),
[`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    mpd(ex_counts, ex_tree)
#>   Saliva     Gums     Nose    Stool 
#> 92.16667 81.40000 79.06667 85.80000 
    weighted_mpd(ex_counts, ex_tree)
#>   Saliva     Gums     Nose    Stool 
#> 44.31146 16.49659 31.26699  1.43670 
    mntd(ex_counts, ex_tree)
#>   Saliva     Gums     Nose    Stool 
#> 74.50000 58.00000 52.66667 60.20000 
    weighted_mntd(ex_counts, ex_tree)
#>    Saliva      Gums      Nose     Stool 
#>  69.38261  34.36795  35.17606 103.64390 

    # The same through alpha_div().
    alpha_div(ex_counts, 'weighted_mntd', tree = ex_tree)
#>    Saliva      Gums      Nose     Stool 
#>  69.38261  34.36795  35.17606 103.64390 

    # Distances between (made-up) trait values instead of a tree.
    traits <- data.frame(
      row.names = colnames(ex_counts),
      gram_pos  = c(1, 0, 1, 0, 1, 1),
      size      = c(0.8, 1.5, 0.7, 1.2, 0.9, 1.0) )
    mpd(ex_counts, taxa_dm = dist(traits))
#>    Saliva      Gums      Nose     Stool 
#> 0.8179849 0.7345747 0.6896271 0.7734185 
```
