# Rényi Entropy

The generalized entropy of order \\q\\: the logarithm of the Hill number
of the same order.

## Usage

``` r
renyi(counts, q = 2, margin = 1L, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- q:

  The order of a Hill number, or of Rényi's or Tsallis's entropy: a
  single number of at least 0, where `Inf` is allowed. Only used by
  [`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
  `renyi()` and
  [`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md).
  Default: `2`, as scikit-bio (which calls it `order`).

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

Rényi's entropy of order \\q\\ is defined as: \$\${}^qH = \frac{1}{1 -
q} \ln \sum\_{i = 1}^{n} P_i^q\$\$

with its limits \\{}^1H = H\\ (Shannon's entropy) and \\{}^\infty H =
-\ln \max(P_i)\\; \\{}^0H\\ is the log of the number of observed
features.

Where:

- \\n\\ : The number of features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

It is `log(hill(counts, q))`. A sample with one feature has entropy 0.
Logarithms are natural; for another base \\b\\, divide by `log(b)`. The
default order, 2, is scikit-bio's.

**Base R Equivalent:**

    x <- ex_counts[1,]
    q <- 2
    p <- x[x > 0] / sum(x)
    log(sum(p ** q)) / (1 - q)

## Same names in other packages

scikit-bio's `renyi()` is the same (its `order` is `q`).
`vegan::renyi()` gives a profile over several orders (`scales`).

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

Rényi, A. (1961). On measures of entropy and information. *Proceedings
of the Fourth Berkeley Symposium on Mathematical Statistics and
Probability*, 1, 547-561.

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md),
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md)

## Examples

``` r
    renyi(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.7077630 0.2097948 1.0237476 0.0130399 
    renyi(ex_counts, q = 1) # Shannon's entropy
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.36684449 1.14222899 0.04824952 
```
