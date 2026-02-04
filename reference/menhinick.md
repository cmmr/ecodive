# Menhinick's Richness Index

A richness metric that normalizes the number of species by the square
root of the total sample size.

## Usage

``` r
menhinick(counts, margin = 1L, cpus = n_cpus())
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

Menhinick's index is defined as: \$\$\frac{n}{\sqrt{X_T}}\$\$

Where:

- \\n\\ : The number of features.

- \\X_T\\ : Total of all counts.

**Base R Equivalent:**

    sum(x > 0) / sqrt(sum(x))

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

Menhinick, E. F. (1964). A comparison of some species-individuals
diversity indices applied to samples of field insects. *Ecology*, 45(4),
859-861. [doi:10.2307/1934933](https://doi.org/10.2307/1934933)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

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
