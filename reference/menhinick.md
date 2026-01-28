# Menhinick's Richness Index

A richness metric that normalizes the number of species by the square
root of the total sample size.

## Usage

``` r
menhinick(counts, margin = 1L, cpus = n_cpus())
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

Menhinick's index is defined as: \$\$\frac{n}{\sqrt{X_T}}\$\$

Where:

- \\n\\ : The number of features.

- \\X_T\\ : Total of all counts.

**Base R Equivalent:**

    sum(x > 0) / sqrt(sum(x))

## References

Menhinick, E. F. (1964). A comparison of some species-individuals
diversity indices applied to samples of field insects. *Ecology*, 45(4),
859-861. [doi:10.2307/1934933](https://doi.org/10.2307/1934933)

## See also

alpha_div

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    menhinick(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.2153528 0.1679783 0.1887016 0.2016195 
```
