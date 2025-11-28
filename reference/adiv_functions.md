# Alpha Diversity Metrics

Alpha Diversity Metrics

## Usage

``` r
ace(counts, cutoff = 10L, margin = 1L, cpus = n_cpus())

berger(counts, norm = "percent", margin = 1L, cpus = n_cpus())

brillouin(counts, margin = 1L, cpus = n_cpus())

chao1(counts, margin = 1L, cpus = n_cpus())

faith(counts, tree = NULL, margin = 1L, cpus = n_cpus())

fisher(counts, digits = 3L, margin = 1L, cpus = n_cpus())

inv_simpson(counts, norm = "percent", margin = 1L, cpus = n_cpus())

margalef(counts, margin = 1L, cpus = n_cpus())

mcintosh(counts, margin = 1L, cpus = n_cpus())

menhinick(counts, margin = 1L, cpus = n_cpus())

observed(counts, margin = 1L, cpus = n_cpus())

shannon(counts, norm = "percent", margin = 1L, cpus = n_cpus())

simpson(counts, norm = "percent", margin = 1L, cpus = n_cpus())

squares(counts, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data where each column is a feature, and
  each row is a sample. Any object coercible with
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) can be given here,
  as well as `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects.

- cutoff:

  The maximum number of observations to consider "rare". Default: `10`.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a `phyloseq`,
  `rbiom`, `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
  Default: `1L`

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

- tree:

  A `phylo`-class object representing the phylogenetic tree for the OTUs
  in `counts`. The OTU identifiers given by `colnames(counts)` must be
  present in `tree`. Can be omitted if a tree is embedded with the
  `counts` object or as `attr(counts, 'tree')`.

- digits:

  Precision of the returned values, in number of decimal places. E.g.
  the default `digits=3` could return `6.392`.

## Value

A numeric vector.

## Formulas

Prerequisite: all counts are whole numbers.

Given:

- \\n\\ : The number of features (e.g. species, OTUs, ASVs, etc).

- \\X_i\\ : Integer count of the \\i\\-th feature.

- \\X_T\\ : Total of all counts (i.e. sequencing depth). \\X_T =
  \sum\_{i=1}^{n} X_i\\

- \\P_i\\ : Proportional abundance of the \\i\\-th feature. \\P_i = X_i
  / X_T\\

- \\F_1\\ : Number of features where \\X_i = 1\\ (i.e. singletons).

- \\F_2\\ : Number of features where \\X_i = 2\\ (i.e. doubletons).

[TABLE]

### Abundance-based Coverage Estimator (ACE)

Given:

- \\n\\ : The number of features (e.g. species, OTUs, ASVs, etc).

- \\r\\ : Rare cutoff. Features with \\\le r\\ counts are considered
  rare.

- \\X_i\\ : Integer count of the \\i\\-th feature.

- \\F_i\\ : Number of features with exactly \\i\\ counts.

- \\F_1\\ : Number of features where \\X_i = 1\\ (i.e. singletons).

- \\F\_{rare}\\ : Number of rare features where \\X_i \le r\\.

- \\F\_{abund}\\ : Number of abundant features where \\X_i \> r\\.

- \\X\_{rare}\\ : Total counts belonging to rare features.

- \\C\_{ace}\\ : The sample abundance coverage estimator, defined below.

- \\\gamma\_{ace}^2\\ : The estimated coefficient of variation, defined
  below.

- \\D\_{ace}\\ : Estimated number of features in the sample.

\\\displaystyle C\_{ace} = 1 - \frac{F_1}{X\_{rare}}\\

\\\displaystyle \gamma\_{ace}^2 = \max\left\[\frac{F\_{rare}
\sum\_{i=1}^{r}i(i-1)F_i}{C\_{ace}X\_{rare}(X\_{rare} - 1)} - 1,
0\right\]\\

\\\displaystyle D\_{ace} = F\_{abund} + \frac{F\_{rare}}{C\_{ace}} +
\frac{F_1}{C\_{ace}}\gamma\_{ace}^2 \\

### Faith's Phylogenetic Diversity (Faith's PD)

Given \\n\\ branches with lengths \\L\\ and a sample's abundances \\A\\
on each of those branches coded as 1 for present or 0 for absent:

\\\sum\_{i = 1}^{n} L_i A_i\\

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
    
    ace(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>    5.0    8.9    6.0    NaN 
    
    chao1(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>    4.5    Inf    6.0    Inf 
    
    squares(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#>  4.492762  8.243044  6.000000 20.793551 
    
```
