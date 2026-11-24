# Good's Coverage

The estimated fraction of a community's individuals that belong to the
features observed in the sample: a check that sequencing was deep
enough.

## Usage

``` r
goods_coverage(counts, margin = 1L, cpus = n_cpus())
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

Good's coverage is defined as: \$\$1 - \frac{F_1}{X_T}\$\$

Where:

- \\F_1\\ : Number of features observed once (singletons).

- \\X_T\\ : Total of all counts (sequencing depth).

It is Good's (1953) estimate of the sample coverage: the probability
that the next individual drawn belongs to a feature already seen. A
sample of only singletons has coverage 0. Chao & Jost (2012) give a less
biased estimator, which uses doubletons too.

**Base R Equivalent:**

    x <- ex_counts[1,]
    1 - sum(x == 1) / sum(x)

## Same names in other packages

scikit-bio's `goods_coverage()` is the same estimate.

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

Good, I. J. (1953). The population frequencies of species and the
estimation of population parameters. *Biometrika*, 40, 237-264.
[doi:10.1093/biomet/40.3-4.237](https://doi.org/10.1093/biomet/40.3-4.237)

Chao, A., & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
*Ecology*, 93, 2533-2547.
[doi:10.1890/11-1952.1](https://doi.org/10.1890/11-1952.1)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

## Examples

``` r
    goods_coverage(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.9971014 0.9977427 1.0000000 0.9934959 
```
