# Brillouin Index

A diversity index derived from information theory, appropriate for fully
censused communities.

## Usage

``` r
brillouin(counts, margin = 1L, cpus = n_cpus())
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

The Brillouin index is defined as: \$\$\frac{\ln{\[(\sum\_{i = 1}^{n}
X_i)!\]} - \sum\_{i = 1}^{n} \ln{(X_i!)}}{\sum\_{i = 1}^{n} X_i}\$\$

Where:

- \\n\\ : The number of features.

- \\X_i\\ : Integer count of the \\i\\-th feature.

**Base R Equivalent:**

    x <- ex_counts[1,]
    # note: lgamma(x + 1) == log(x!)
    (lgamma(sum(x) + 1) - sum(lgamma(x + 1))) / sum(x)

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

Brillouin, L. (1956). Science and information theory. Academic Press.

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Diversity metrics:
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md)

## Examples

``` r
    brillouin(ex_counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.72541640 0.35924823 1.13029903 0.04175076 
```
