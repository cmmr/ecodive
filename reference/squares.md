# Squares Richness Estimator

A richness estimator based on the concept of "squares" (counts of
species observed once or twice).

## Usage

``` r
squares(counts, margin = 1L, cpus = n_cpus())
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

The Squares estimator is defined as: \$\$n + \frac{(F_1)^2
\sum\_{i=1}^{n} (X_i)^2}{X_T^2 - nF_1}\$\$

Where:

- \\n\\ : The number of observed features.

- \\X_T\\ : Total of all counts.

- \\F_1\\ : Number of features observed once (singletons).

- \\X_i\\ : Integer count of the \\i\\-th feature.

**Base R Equivalent:**

    N  <- sum(x)      # sampling depth
    S  <- sum(x > 0)  # observed features
    F1 <- sum(x == 1) # singletons
    S + ((sum(x^2) * (F1^2)) / ((N^2) - F1 * S))

## References

Alroy, J. (2018). Limits to species richness estimates based on
subsampling. *Paleobiology*, 44(2), 177-194.
[doi:10.1017/pab.2017.38](https://doi.org/10.1017/pab.2017.38)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md)

## Examples

``` r
    squares(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#>  4.492762  8.243044  6.000000 20.793551 
```
