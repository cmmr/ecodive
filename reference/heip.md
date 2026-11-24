# Heip's Evenness

An evenness index derived from the exponential of Shannon's entropy.

## Usage

``` r
heip(counts, margin = 1L, cpus = n_cpus())
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

Heip's evenness is defined as: \$\$\frac{e^H - 1}{n - 1}\$\$

Where:

- \\H\\ : Shannon's entropy,
  [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md).

- \\n\\ : The number of observed features.

A sample with one feature has evenness 1, as in scikit-bio.

**Base R Equivalent:**

    x <- ex_counts[1,]
    p <- x[x > 0] / sum(x)
    (exp(-sum(p * log(p))) - 1) / (length(p) - 1)

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

Heip, C. (1974). A new index measuring evenness. *Journal of the Marine
Biological Association of the UK*, 54, 555-557.
[doi:10.1017/S0025315400022736](https://doi.org/10.1017/S0025315400022736)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Evenness metrics:
[`mcintosh_evenness()`](https://cmmr.github.io/ecodive/reference/mcintosh_evenness.md),
[`pielou()`](https://cmmr.github.io/ecodive/reference/pielou.md),
[`simpson_evenness()`](https://cmmr.github.io/ecodive/reference/simpson_evenness.md)

## Examples

``` r
    heip(ex_counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.36615009 0.11079337 0.42674914 0.01235812 
```
