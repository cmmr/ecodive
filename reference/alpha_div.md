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
  q = 2,
  bias_corrected = TRUE,
  cpus = n_cpus()
)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- metric:

  The name of an alpha diversity metric. One of
  `c('ace', 'berger', 'brillouin', 'chao1', 'dominance', 'faith', 'fisher', 'gini', 'goods_coverage', 'heip', 'hill', 'inv_simpson', 'margalef', 'mcintosh', 'mcintosh_evenness', 'menhinick', 'observed', 'pielou', 'renyi', 'shannon', 'simpson', 'simpson_evenness', 'squares', 'strong', 'tsallis')`.
  Case-insensitive and partial name matching is supported, as are
  scikit-bio's names (such as `'pielou_e'`, `'simpson_d'`,
  `'observed_features'` and `'enspie'`). Programmatic access via
  `list_metrics('alpha')`.

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

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- q:

  The order of a Hill number, or of Rényi's or Tsallis's entropy: a
  single number of at least 0, where `Inf` is allowed. Only used by
  [`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
  [`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md).
  Default: `2`, as scikit-bio (which calls it `order`).

- bias_corrected:

  For Chao1 only: use the bias-corrected form (the default, as
  scikit-bio and QIIME 2), or the classic form when a sample has both
  singletons and doubletons (`FALSE`). See
  [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md).

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A numeric vector, one value per sample, named by the samples. To compute
several metrics, call `alpha_div()` once for each (see the examples).

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

- McIntosh Index

- McIntosh's Evenness

- Good's Coverage

#### Can Use Proportional Data

- Observed Features

- Shannon Index

- Gini-Simpson Index

- Inverse Simpson Index

- Hill Numbers, Rényi and Tsallis Entropies

- Pielou's, Heip's and Simpson's Evenness

- Berger-Parker Index

- Simpson's and Strong's Dominance Indices

- Gini Index

- Faith's PD

### Differences from scikit-bio

The metrics follow `skbio.diversity.alpha` (scikit-bio 0.7.4), which
QIIME 2's `diversity alpha` runs, with these differences:

- Logarithms are natural. scikit-bio's
  [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
  `pielou_e()` and
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) take a
  `base`; divide ecodive's
  [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md) and
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) by
  `log(base)` instead. Pielou's evenness doesn't depend on the base.

- No `finite` option for
  [`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md),
  [`dominance()`](https://cmmr.github.io/ecodive/reference/dominance.md)
  and
  [`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md)
  (scikit-bio's estimator without replacement).

- The [`gini()`](https://cmmr.github.io/ecodive/reference/gini.md) index
  is computed by rectangles only (scikit-bio's default
  `method = "rectangles"`, not `"trapezoids"`).

- A sample with no observations is `NA` for every metric. scikit-bio
  gives `nan` for most metrics, but `0` for
  [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md) and
  `gini_index()`.

- `chao1(bias_corrected = FALSE)` is the classic formula, falling back
  on the bias-corrected one when a sample has no singletons or no
  doubletons, exactly as scikit-bio does.

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
#>      4      6      6     11 
    
    # Faith PD values
    alpha_div(ex_counts, 'faith', tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    191    215    202 
    
    # Hill number of order 1 (the exponential of Shannon's entropy)
    alpha_div(ex_counts, 'hill', q = 1)
#>   Saliva     Gums     Nose    Stool 
#> 2.098450 1.443173 3.133746 1.049432 
    
    # Several metrics: one call each
    sapply(c('observed', 'shannon', 'pielou'), alpha_div, counts = ex_counts)
#>        observed    shannon     pielou
#> Saliva        4 0.74119910 0.53466213
#> Gums          5 0.36684449 0.22793330
#> Nose          6 1.14222899 0.63749014
#> Stool         5 0.04824952 0.02997911
    
    
```
