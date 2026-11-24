# McIntosh's Evenness

An evenness index derived from the McIntosh index.

## Usage

``` r
mcintosh_evenness(counts, margin = 1L, cpus = n_cpus())
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
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Details

McIntosh's evenness is defined as: \$\$\frac{\sqrt{\sum\_{i = 1}^{n}
(X_i)^2}}{\sqrt{(X_T - n + 1)^2 + n - 1}}\$\$

Where:

- \\n\\ : The number of observed features.

- \\X_i\\ : Integer count of the \\i\\-th feature.

- \\X_T\\ : Total of all counts.

This follows scikit-bio's `mcintosh_e()` (Heip & Engels 1974).

**Base R Equivalent:**

    x <- ex_counts[1,]
    n <- sum(x > 0)
    sqrt(sum(x^2)) / sqrt((sum(x) - n + 1)^2 + n - 1)

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

Heip, C., & Engels, P. (1974). Comparing species diversity and evenness
indices. *Journal of the Marine Biological Association of the UK*, 54,
559-563.
[doi:10.1017/S0025315400022748](https://doi.org/10.1017/S0025315400022748)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Evenness metrics:
[`heip()`](https://cmmr.github.io/ecodive/reference/heip.md),
[`pielou()`](https://cmmr.github.io/ecodive/reference/pielou.md),
[`simpson_evenness()`](https://cmmr.github.io/ecodive/reference/simpson_evenness.md)

## Examples

``` r
    mcintosh_evenness(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.7081066 0.9044981 0.6023489 1.0000000 
```
