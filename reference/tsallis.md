# Tsallis Entropy

The non-additive generalized entropy of order \\q\\.

## Usage

``` r
tsallis(counts, q = 2, margin = 1L, cpus = n_cpus())
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
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
  `tsallis()`. Default: `2`, as scikit-bio (which calls it `order`).

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

Tsallis's entropy of order \\q\\ is defined as: \$\${}^qS = \frac{1}{q -
1} \left(1 - \sum\_{i = 1}^{n} P_i^q\right)\$\$

with its limits \\{}^1S = H\\ (Shannon's entropy) and \\{}^\infty S =
0\\.

Where:

- \\n\\ : The number of features.

- \\P_i\\ : Proportional abundance of the \\i\\-th feature.

At \\q = 0\\ it is the number of observed features minus 1, and at \\q =
2\\ the Gini-Simpson index,
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md). A
sample with one feature has entropy 0. The default order, 2, is
scikit-bio's.

**Base R Equivalent:**

    x <- ex_counts[1,]
    q <- 2
    p <- x[x > 0] / sum(x)
    (1 - sum(p ** q)) / (q - 1)

## Same names in other packages

scikit-bio's `tsallis()` is the same (its `order` is `q`).

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

Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs
statistics. *Journal of Statistical Physics*, 52, 479-487.
[doi:10.1007/BF01016429](https://doi.org/10.1007/BF01016429)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Diversity metrics:
[`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md),
[`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md),
[`hill()`](https://cmmr.github.io/ecodive/reference/hill.md),
[`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md),
[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md),
[`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md),
[`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md)

## Examples

``` r
    tsallis(ex_counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.50725478 0.18924937 0.64075388 0.01295525 
```
