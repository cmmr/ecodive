# McIntosh Index

A dominance index based on the Euclidean distance from the origin.

## Usage

``` r
mcintosh(counts, margin = 1L, cpus = n_cpus())
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

The McIntosh index is defined as: \$\$\frac{X_T - \sqrt{\sum\_{i =
1}^{n} (X_i)^2}}{X_T - \sqrt{X_T}}\$\$

Where:

- \\n\\ : The number of features.

- \\X_i\\ : Integer count of the \\i\\-th feature.

- \\X_T\\ : Total of all counts.

**Base R Equivalent:**

    x <- ex_counts[1,]
    (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))

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

McIntosh, R. P. (1967). An index of diversity and the relation of
certain concepts to diversity. *Ecology*, 48(3), 392-404.
[doi:10.2307/1932674](https://doi.org/10.2307/1932674)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Dominance metrics:
[`berger()`](https://cmmr.github.io/ecodive/reference/berger.md)

## Examples

``` r
    mcintosh(ex_counts)
#>      Saliva        Gums        Nose       Stool 
#> 0.315000947 0.103044944 0.413637580 0.006771808 
```
