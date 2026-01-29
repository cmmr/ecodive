# Variance-Adjusted Weighted UniFrac

A weighted UniFrac that adjusts for the expected variance of the metric.

## Usage

``` r
variance_adjusted_unifrac(
  counts,
  tree = NULL,
  margin = 1L,
  pairs = NULL,
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

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

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

The Variance-Adjusted Weighted UniFrac distance is defined as:
\$\$\frac{\sum\_{i=1}^{n} L_i\frac{\|P_i - Q_i\|}{\sqrt{(P_i + Q_i)(2 -
P_i - Q_i)}} }{\sum\_{i=1}^{n} L_i\frac{P_i + Q_i}{\sqrt{(P_i + Q_i)(2 -
P_i - Q_i)}} }\$\$

Where:

- \\n\\ : The number of branches in the tree.

- \\L_i\\ : The length of the \\i\\-th branch.

- \\P_i\\, \\Q_i\\ : The proportion of the community descending from
  branch \\i\\ in sample P and Q.

## References

Chang, Q., Luan, Y., & Sun, F. (2011). Variance adjusted weighted
UniFrac: a powerful beta diversity measure for comparing communities
based on phylogeny. *BMC Bioinformatics*, 12, 118.
[doi:10.1186/1471-2105-12-118](https://doi.org/10.1186/1471-2105-12-118)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Phylogenetic metrics:
[`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
[`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md),
[`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md),
[`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md),
[`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)

## Examples

``` r
    variance_adjusted_unifrac(ex_counts, tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4242631                    
#> Nose  0.7753369 0.5565010          
#> Stool 0.9655749 0.9807634 0.9785147
```
