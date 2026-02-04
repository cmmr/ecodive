# Berger-Parker Index

A measure of the numerical importance of the most abundant species.

## Usage

``` r
berger(counts, margin = 1L, cpus = n_cpus())
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

The Berger-Parker index is defined as the proportional abundance of the
most dominant feature: \$\$\max(P_i)\$\$

Where:

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

**Base R Equivalent:**

    p <- x / sum(x)
    max(p)

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

Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic
foraminifera in deep-sea sediments. *Science*, 168(3937), 1345-1347.
[doi:10.1126/science.168.3937.1345](https://doi.org/10.1126/science.168.3937.1345)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Dominance metrics:
[`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md)

## Examples

``` r
    berger(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.5217391 0.8950339 0.4925816 0.9934959 
```
