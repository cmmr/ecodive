# Chao1 Richness Estimator

A non-parametric estimator of the lower bound of species richness.

## Usage

``` r
chao1(counts, margin = 1L, cpus = n_cpus())
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

The Chao1 estimator uses the ratio of singletons to doubletons to
estimate the number of missing species: \$\$n + \frac{(F_1)^2}{2
F_2}\$\$

Where:

- \\n\\ : The number of observed features.

- \\F_1\\ : Number of features observed once (singletons).

- \\F_2\\ : Number of features observed twice (doubletons).

**Base R Equivalent:**

    sum(x>0) + (sum(x == 1) ** 2) / (2 * sum(x == 2))

## References

Chao, A. (1984). Nonparametric estimation of the number of classes in a
population. *Scandinavian Journal of Statistics*, 11, 265-270.

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    chao1(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>    4.5    Inf    6.0    Inf 
```
