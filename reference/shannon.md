# Shannon Diversity Index

A commonly used diversity index accounting for both abundance and
evenness.

## Usage

``` r
shannon(counts, norm = "percent", margin = 1L, cpus = n_cpus())
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

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a `phyloseq`,
  `rbiom`, `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
  Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The Shannon index (entropy) is defined as: \$\$-\sum\_{i = 1}^{n} P_i
\times \ln(P_i)\$\$

Where:

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

**Base R Equivalent:**

    p <- x / sum(x)
    -sum(p * log(p))

## References

Shannon, C. E. (1948). A mathematical theory of communication. *Bell
System Technical Journal*, 27, 379-423.

## See also

alpha_div

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md)

## Examples

``` r
    shannon(ex_counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.36684449 1.14222899 0.04824952 
```
