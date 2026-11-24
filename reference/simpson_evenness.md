# Simpson's Evenness

The inverse Simpson index relative to the number of observed features.

## Usage

``` r
simpson_evenness(counts, margin = 1L, cpus = n_cpus())
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

Simpson's evenness is defined as: \$\$\frac{1}{n \sum\_{i = 1}^{n}
P_i^2}\$\$

Where:

- \\n\\ : The number of observed features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

It is `inv_simpson() / observed()`.

**Base R Equivalent:**

    x <- ex_counts[1,]
    p <- x[x > 0] / sum(x)
    1 / (length(p) * sum(p ** 2))

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

Smith, B., & Wilson, J. B. (1996). A consumer's guide to evenness
indices. *Oikos*, 76, 70-82.
[doi:10.2307/3545749](https://doi.org/10.2307/3545749)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Evenness metrics:
[`heip()`](https://cmmr.github.io/ecodive/reference/heip.md),
[`mcintosh_evenness()`](https://cmmr.github.io/ecodive/reference/mcintosh_evenness.md),
[`pielou()`](https://cmmr.github.io/ecodive/reference/pielou.md)

## Examples

``` r
    simpson_evenness(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.5073616 0.2466850 0.4639345 0.2026251 
```
