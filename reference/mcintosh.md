# McIntosh Index

A dominance index based on the Euclidean distance from the origin.

## Usage

``` r
mcintosh(counts, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a special object
  class (e.g. `phyloseq`). Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The McIntosh index is defined as: \$\$\frac{X_T - \sqrt{\sum\_{i =
1}^{n} (X_i)^2}}{X_T - \sqrt{X_T}}\$\$

Where:

- \\X_i\\ : Integer count of the \\i\\-th feature.

- \\X_T\\ : Total of all counts.

**Base R Equivalent:**

    (sum(x) - sqrt(sum(x^2))) / (sum(x) - sqrt(sum(x)))

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
