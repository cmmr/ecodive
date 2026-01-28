# Normalized Weighted UniFrac

Weighted UniFrac normalized by the tree length to allow comparison
between trees.

## Usage

``` r
normalized_unifrac(
  counts,
  tree = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data where each column is a feature, and
  each row is a sample. Any object coercible with
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) can be given here,
  as well as `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. For optimal performance with very
  large datasets, see the guide in
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md).

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a `phyloseq`,
  `rbiom`, `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
  Default: `1L`

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

The Normalized Weighted UniFrac distance is defined as:
\$\$\frac{\sum\_{i=1}^{n} L_i\|P_i - Q_i\|}{\sum\_{i=1}^{n} L_i(P_i +
Q_i)}\$\$

Where:

- \\n\\ : The number of branches in the tree.

- \\L_i\\ : The length of the \\i\\-th branch.

- \\P_i\\, \\Q_i\\ : The proportion of the community descending from
  branch \\i\\ in sample P and Q.

## References

Lozupone, C. A., Hamady, M., Kelley, S. T., & Knight, R. (2007).
Quantitative and qualitative beta diversity measures lead to different
insights into factors that structure microbial communities. *Applied and
Environmental Microbiology*, 73(5), 1576-1585.
[doi:10.1128/AEM.01996-06](https://doi.org/10.1128/AEM.01996-06)

## See also

beta_div

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md),
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md),
[`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    normalized_unifrac(ex_counts, tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4343436                    
#> Nose  0.7891160 0.6714592          
#> Stool 0.9707293 0.9868636 0.9925679
```
