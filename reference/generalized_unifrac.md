# Generalized UniFrac (GUniFrac)

A unified UniFrac distance that balances the weight of abundant and rare
lineages.

## Usage

``` r
generalized_unifrac(
  counts,
  tree = NULL,
  alpha = 0.5,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- alpha:

  How much weight to give to relative abundances; a value between 0 and
  1, inclusive. Setting `alpha=1` is equivalent to
  [`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md).

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The Generalized UniFrac distance is defined as:
\$\$\frac{\sum\_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}\left\|\frac{P_i -
Q_i}{P_i + Q_i}\right\|}{\sum\_{i=1}^{n} L_i(P_i + Q_i)^{\alpha}}\$\$

Where:

- \\n\\ : The number of branches in the tree.

- \\L_i\\ : The length of the \\i\\-th branch.

- \\P_i\\, \\Q_i\\ : The proportion of the community descending from
  branch \\i\\ in sample P and Q.

- \\\alpha\\ : A scalable weighting factor.

**Parameter: alpha**

The `alpha` parameter controls the weight given to abundant lineages.
\\\alpha = 1\\ corresponds to Weighted UniFrac, while \\\alpha = 0\\
corresponds to Unweighted UniFrac.

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

Chen, J., Bittinger, K., Charlson, E. S., Hoffmann, C., Lewis, J., Wu,
G. D., ... & Li, H. (2012). Associating microbiome composition with
environmental covariates using generalized UniFrac distances.
*Bioinformatics*, 28(16), 2106-2113.
[doi:10.1093/bioinformatics/bts342](https://doi.org/10.1093/bioinformatics/bts342)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md),
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md),
[`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    generalized_unifrac(ex_counts, tree = ex_tree, alpha = 0.5)
#>          Saliva      Gums      Nose
#> Gums  0.4471644                    
#> Nose  0.8215129 0.7607876          
#> Stool 0.9727827 0.9784242 0.9730332
```
