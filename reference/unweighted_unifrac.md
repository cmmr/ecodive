# Unweighted UniFrac

A phylogenetic distance metric that accounts for the presence/absence of
lineages.

## Usage

``` r
unweighted_unifrac(
  counts,
  tree = NULL,
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

The Unweighted UniFrac distance is defined as:
\$\$\frac{1}{n}\sum\_{i=1}^{n} L_i\|A_i - B_i\|\$\$

Where:

- \\n\\ : The number of branches in the tree.

- \\L_i\\ : The length of the \\i\\-th branch.

- \\A_i\\, \\B_i\\ : Binary values (0 or 1) indicating if descendants of
  branch \\i\\ are present in sample A or B.

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

Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method
for comparing microbial communities. *Applied and Environmental
Microbiology*, 71(12), 8228-8235.
[doi:10.1128/AEM.71.12.8228-8235.2005](https://doi.org/10.1128/AEM.71.12.8228-8235.2005)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md),
[`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md),
[`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    unweighted_unifrac(ex_counts, tree = ex_tree)
#>           Saliva       Gums       Nose
#> Gums  0.05759162                      
#> Nose  0.16279070 0.11162791           
#> Stool 0.22325581 0.17209302 0.06046512
```
