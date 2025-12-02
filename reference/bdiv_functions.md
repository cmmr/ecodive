# Beta Diversity Metrics

Beta Diversity Metrics

## Usage

``` r
aitchison(
  counts,
  pseudocount = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

bhattacharyya(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

bray(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

canberra(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

chebyshev(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

chord(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

clark(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

divergence(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

euclidean(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

gower(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

hellinger(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

horn(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

jensen(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

jsd(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

lorentzian(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

manhattan(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

matusita(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

minkowski(
  counts,
  norm = "percent",
  power = 1.5,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

morisita(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

motyka(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

psym_chisq(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

soergel(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

squared_chisq(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

squared_chord(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

squared_euclidean(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

topsoe(counts, norm = "percent", margin = 1L, pairs = NULL, cpus = n_cpus())

wave_hedges(
  counts,
  norm = "percent",
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

hamming(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

jaccard(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

ochiai(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

sorensen(counts, margin = 1L, pairs = NULL, cpus = n_cpus())

unweighted_unifrac(
  counts,
  tree = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

weighted_unifrac(
  counts,
  tree = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

normalized_unifrac(
  counts,
  tree = NULL,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

generalized_unifrac(
  counts,
  tree = NULL,
  alpha = 0.5,
  margin = 1L,
  pairs = NULL,
  cpus = n_cpus()
)

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

  A numeric matrix of count data where each column is a feature, and
  each row is a sample. Any object coercible with
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) can be given here,
  as well as `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. For optimal performance with very
  large datasets, see the guide in
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md).

- pseudocount:

  The value to add to all counts in `counts` to prevent taking `log(0)`
  for unobserved features. The default, `NULL`, selects the smallest
  non-zero value in `counts`.

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

- norm:

  Normalize the incoming counts. Options are:

  `norm = "percent"` -

  :   Relative abundance (sample abundances sum to 1).

  `norm = "binary"` -

  :   Unweighted presence/absence (each count is either 0 or 1).

  `norm = "clr"` -

  :   Centered log ratio.

  `norm = "none"` -

  :   No transformation.

  Default: `'percent'`, which is the expected input for these formulas.

- power:

  Scaling factor for the magnitude of differences between communities
  (\\p\\). Default: `1.5`

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- alpha:

  How much weight to give to relative abundances; a value between 0 and
  1, inclusive. Setting `alpha=1` is equivalent to
  `normalized_unifrac()`.

## Value

A `dist` object.

## Formulas

Given:

- \\n\\ : The number of features.

- \\X_i\\, \\Y_i\\ : Absolute counts for the \\i\\-th feature in samples
  \\X\\ and \\Y\\.

- \\X_T\\, \\Y_T\\ : Total counts in each sample. \\X_T =
  \sum\_{i=1}^{n} X_i\\

- \\P_i\\, \\Q_i\\ : Proportional abundances of \\X_i\\ and \\Y_i\\.
  \\P_i = X_i / X_T\\

- \\X_L\\, \\Y_L\\ : Mean log of abundances. \\X_L =
  \frac{1}{n}\sum\_{i=1}^{n} \ln{X_i}\\

- \\R_i\\ : The range of the \\i\\-th feature across all samples (max -
  min).

[TABLE]

### Presence / Absence

Given:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common.

[TABLE]

### Phylogenetic

Given \\n\\ branches with lengths \\L\\ and a pair of samples' binary
(\\A\\ and \\B\\) or proportional abundances (\\P\\ and \\Q\\) on each
of those branches.

[TABLE]

See `vignette('unifrac')` for detailed example UniFrac calculations.

## References

Levy, A., Shalom, B. R., & Chalamish, M. (2024). A guide to similarity
measures. *arXiv*.

Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures
between probability density functions. *International Journal of
Mathematical Models and Methods in Applied Sciences*, 1(4), 300–307.

## Examples

``` r
    # Example counts matrix
    t(ex_counts)
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  793   22     1
#> Bacteroides            2    4    2   611
#> Corynebacterium        0    0  498     1
#> Haemophilus          180   87    2     1
#> Propionibacterium      1    1  251     0
#> Staphylococcus         0    1  236     1
    
    bray(ex_counts)
#>          Saliva      Gums      Nose
#> Gums  0.4265973                    
#> Nose  0.9713843 0.9720256          
#> Stool 0.9909509 0.9911046 0.9915177
    
    jaccard(ex_counts)
#>          Saliva      Gums      Nose
#> Gums  0.2000000                    
#> Nose  0.3333333 0.1666667          
#> Stool 0.5000000 0.3333333 0.1666667
    
    generalized_unifrac(ex_counts, tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4471644                    
#> Nose  0.8215129 0.7607876          
#> Stool 0.9727827 0.9784242 0.9730332
    
    # Only calculate distances for Saliva vs all.
    bray(ex_counts, pairs = 1:3)
#>          Saliva      Gums      Nose
#> Gums  0.4265973                    
#> Nose  0.9713843        NA          
#> Stool 0.9909509        NA        NA
    
```
