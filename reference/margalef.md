# Margalef's Richness Index

A richness metric that normalizes the number of species by the log of
the total sample size.

## Usage

``` r
margalef(counts, margin = 1L, cpus = n_cpus())
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

Margalef's index is defined as: \$\$\frac{n - 1}{\ln{X_T}}\$\$

Where:

- \\n\\ : The number of features.

- \\X_T\\ : Total of all counts (sequencing depth).

**Base R Equivalent:**

    (sum(x > 0) - 1) / log(sum(x))

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

Margalef, R. (1958). Information theory in ecology. *General Systems*,
3, 36-71.

Gamito, S. (2010). Caution is needed when applying Margalef diversity
index. *Ecological Indicators*, 10(2), 550-551.
[doi:10.1016/j.ecolind.2009.07.006](https://doi.org/10.1016/j.ecolind.2009.07.006)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    margalef(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.5133870 0.5893866 0.7226796 0.6228956 
```
