# Berger-Parker Index

A measure of the numerical importance of the most abundant species.

## Usage

``` r
berger(counts, norm = "percent", margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features). Also
  supports `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. See
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
  for optimizing large datasets.

- norm:

  Normalize the incoming counts. Options are:

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'binary'`: Unweighted presence/absence (each count is either 0 or
    1).

  - `'clr'`: Centered log ratio.

  - `'none'`: No transformation.

  Default: `'percent'`, which is the expected input for these formulas.

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

    max(x / sum(x))

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
