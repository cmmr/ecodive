# Shannon Diversity Index

A commonly used diversity index accounting for both abundance and
evenness.

## Usage

``` r
shannon(counts, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The Shannon index (entropy) is defined as: \$\$-\sum\_{i = 1}^{n} P_i
\times \ln(P_i)\$\$

Where:

- \\n\\ : The number of features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

**Base R Equivalent:**

    p <- x / sum(x)
    -sum(p * log(p))

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

Shannon, C. E. (1948). A mathematical theory of communication. *Bell
System Technical Journal*, 27, 379-423.

Shannon, C. E., & Weaver, W. (1949). *The Mathematical Theory of
Communication*. University of Illinois Press.

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

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
