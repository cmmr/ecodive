# Brillouin Index

A diversity index derived from information theory, appropriate for fully
censused communities.

## Usage

``` r
brillouin(counts, margin = 1L, cpus = n_cpus())
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

The Brillouin index is defined as: \$\$\frac{\ln{\[(\sum\_{i = 1}^{n}
X_i)!\]} - \sum\_{i = 1}^{n} \ln{(X_i!)}}{\sum\_{i = 1}^{n} X_i}\$\$

Where:

- \\n\\ : The number of features.

- \\X_i\\ : Integer count of the \\i\\-th feature.

**Base R Equivalent:**

    # note: lgamma(x + 1) == log(x!)
    (lgamma(sum(x) + 1) - sum(lgamma(x + 1))) / sum(x)

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
