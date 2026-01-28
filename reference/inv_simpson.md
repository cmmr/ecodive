# Inverse Simpson Index

A transformation of the Simpson index that represents the "effective
number of species".

## Usage

``` r
inv_simpson(counts, norm = "percent", margin = 1L, cpus = n_cpus())
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

The Inverse Simpson index is defined as: \$\$1 / \sum\_{i = 1}^{n}
P_i^2\$\$

Where:

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

**Base R Equivalent:**

    p <- x / sum(x)
    1 / sum(p ** 2)

## References

Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688.
[doi:10.1038/163688a0](https://doi.org/10.1038/163688a0)

## See also

alpha_div

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md)

## Examples

``` r
    inv_simpson(ex_counts)
#>   Saliva     Gums     Nose    Stool 
#> 2.029446 1.233425 2.783607 1.013125 
```
