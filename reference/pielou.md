# Pielou's Evenness

Shannon's entropy relative to its largest possible value for the number
of observed features.

## Usage

``` r
pielou(counts, margin = 1L, cpus = n_cpus())
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

Pielou's evenness \\J\\ is defined as: \$\$\frac{H}{\ln n}\$\$

Where:

- \\H\\ : Shannon's entropy,
  [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md).

- \\n\\ : The number of observed features.

It is one of the alpha diversity outputs of QIIME 2's `core-metrics`
(`evenness_vector`). It doesn't depend on the base of the logarithm.

A sample with one feature has evenness 1, as in scikit-bio (where
`vegan::diversity(x) / log(vegan::specnumber(x))` gives `NaN`).

**Base R Equivalent:**

    x <- ex_counts[1,]
    p <- x[x > 0] / sum(x)
    -sum(p * log(p)) / log(length(p))

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

Pielou, E. C. (1966). The measurement of diversity in different types of
biological collections. *Journal of Theoretical Biology*, 13, 131-144.
[doi:10.1016/0022-5193(66)90013-0](https://doi.org/10.1016/0022-5193%2866%2990013-0)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Evenness metrics:
[`heip()`](https://cmmr.github.io/ecodive/reference/heip.md),
[`mcintosh_evenness()`](https://cmmr.github.io/ecodive/reference/mcintosh_evenness.md),
[`simpson_evenness()`](https://cmmr.github.io/ecodive/reference/simpson_evenness.md)

## Examples

``` r
    pielou(ex_counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.53466213 0.22793330 0.63749014 0.02997911 
```
