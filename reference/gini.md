# Gini Index

The inequality of the features' abundances: 0 when all features of the
table are equally abundant, approaching 1 when one feature has them all.

## Usage

``` r
gini(counts, margin = 1L, cpus = n_cpus())
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

The Gini index is computed from the Lorenz curve of all \\K\\ features
of the table, by rectangles, as scikit-bio's `gini_index()`:
\$\$\max\left(0,\\ 1 - \frac{2}{K} \sum\_{j = 1}^{K} L_j\right)\$\$

Where:

- \\K\\ : The number of features in the table, including those absent
  from the sample.

- \\L_j\\ : The Lorenz curve: the share of the sample's total held by
  the \\j\\ least abundant features, \\\sum\_{i \le j} X\_{(i)} / X_T\\,
  with \\X\_{(1)} \le \ldots \le X\_{(K)}\\ the sorted abundances.

**The Gini index depends on how many features the table has,** including
features the sample lacks: each absent feature makes the sample look
less even. Filtering features out of a table (or adding all-zero ones)
changes it, so compare samples only from the same table. A sample of one
table with \\K\\ features, all equally abundant, has index 0; a sample
with all its counts in one feature has index \\1 - 2/K\\ (by rectangles,
the Lorenz curve's area is overestimated by \\1/(2K)\\).

`alpha_div(counts, 'gini')` is this index; the Gini-Simpson index is
`'simpson'`.

**Base R Equivalent:**

    x <- ex_counts[1,]
    L <- cumsum(sort(x)) / sum(x)
    max(0, 1 - 2 * sum(L) / length(x))

## Same names in other packages

scikit-bio's `gini_index()` is the same index (with its default
`method = "rectangles"`).

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

Gini, C. (1912). *Variabilità e mutabilità*. Cuppini, Bologna.

Wittebolle, L., Marzorati, M., Clement, L., Balloi, A., Daffonchio, D.,
Heylen, K., De Vos, P., Verstraete, W., & Boon, N. (2009). Initial
community evenness favours functionality under selective stress.
*Nature*, 458, 623-626.
[doi:10.1038/nature07840](https://doi.org/10.1038/nature07840)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Dominance metrics:
[`berger()`](https://cmmr.github.io/ecodive/reference/berger.md),
[`dominance()`](https://cmmr.github.io/ecodive/reference/dominance.md),
[`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md),
[`strong()`](https://cmmr.github.io/ecodive/reference/strong.md)

## Examples

``` r
    gini(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.5033816 0.6282919 0.4005935 0.6612466 
```
