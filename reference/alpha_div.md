# Alpha Diversity Wrapper Function

Alpha Diversity Wrapper Function

## Usage

``` r
alpha_div(
  counts,
  metric,
  norm = "percent",
  cutoff = 10L,
  digits = 3L,
  tree = NULL,
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

  The name of an alpha diversity metric. One of
  `c('ace', 'berger', 'brillouin', 'chao1', 'faith', 'fisher', 'inv_simpson', 'margalef', 'mcintosh', 'menhinick', 'observed', 'shannon', 'simpson', 'squares')`.
  Case-insensitive and partial name matching is supported. Programmatic
  access via `list_metrics('alpha')`.

- norm:

  Normalize the incoming counts. Options are:

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  - `'none'`: No transformation.

  Default: `'percent'`, which is the expected input for these formulas.

- cutoff:

  The maximum number of observations to consider "rare". Default: `10`.

- digits:

  Precision of the returned values, in number of decimal places. E.g.
  the default `digits=3` could return `6.392`.

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

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

### Integer Count Requirements

A frequent and critical error in alpha diversity analysis is providing
the wrong type of data to a metric's formula. Some indices are
mathematically defined based on counts of individuals and require raw,
integer abundance data. Others are based on proportional abundances and
can accept either integer counts (which are then converted to
proportions) or pre-normalized proportional data. Using proportional
data with a metric that requires integer counts will return an error
message.

#### Requires Integer Counts Only

- Chao1

- ACE

- Squares Richness Estimator

- Margalef's Index

- Menhinick's Index

- Fisher's Alpha

- Brillouin Index

#### Can Use Proportional Data

- Observed Features

- Shannon Index

- Gini-Simpson Index

- Inverse Simpson Index

- Berger-Parker Index

- McIntosh Index

- Faith's PD

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
    
    # Shannon diversity values
    alpha_div(ex_counts, 'Shannon')
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.36684449 1.14222899 0.04824952 
    
    # Chao1 diversity values
    alpha_div(ex_counts, 'c')
#> Saliva   Gums   Nose  Stool 
#>    4.5    Inf    6.0    Inf 
    
    # Faith PD values
    alpha_div(ex_counts, 'faith', tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    191    215    202 
    
    
```
