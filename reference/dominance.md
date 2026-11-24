# Simpson's Dominance Index

The probability that two entities taken at random (with replacement) are
of the same type.

## Usage

``` r
dominance(counts, margin = 1L, cpus = n_cpus())
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

Simpson's dominance index \\D\\ is defined as: \$\$\sum\_{i = 1}^{n}
P_i^2\$\$

Where:

- \\n\\ : The number of features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

It is `1 - simpson()` (the Gini-Simpson index) and `1 / inv_simpson()`.
scikit-bio calls it `dominance()` and `simpson_d()`; both names are
aliases in
[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md).

**Base R Equivalent:**

    x <- ex_counts[1,]
    p <- x / sum(x)
    sum(p ** 2)

## Same names in other packages

scikit-bio's `dominance()` is the same index.

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

Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688.
[doi:10.1038/163688a0](https://doi.org/10.1038/163688a0)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Dominance metrics:
[`berger()`](https://cmmr.github.io/ecodive/reference/berger.md),
[`gini()`](https://cmmr.github.io/ecodive/reference/gini.md),
[`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md),
[`strong()`](https://cmmr.github.io/ecodive/reference/strong.md)

## Examples

``` r
    dominance(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.4927452 0.8107506 0.3592461 0.9870447 
```
