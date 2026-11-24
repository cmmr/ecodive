# Strong's Dominance Index

The largest gap between the cumulative share of the most abundant
features and the share they would have if all were equally abundant.

## Usage

``` r
strong(counts, margin = 1L, cpus = n_cpus())
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

Strong's dominance index is defined as: \$\$\max\_{i} \left(
\frac{c_i}{X_T} - \frac{i}{n} \right)\$\$

Where:

- \\n\\ : The number of observed features.

- \\c_i\\ : The sum of the \\i\\ largest counts.

- \\X_T\\ : Total of all counts.

**Base R Equivalent:**

    x <- ex_counts[1,]
    x <- sort(x[x > 0], decreasing = TRUE)
    max(cumsum(x) / sum(x) - seq_along(x) / length(x))

## Same names in other packages

scikit-bio's `strong()` is the same index.

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

Strong, W. L. (2002). Assessing species abundance unevenness within and
between plant communities. *Community Ecology*, 3, 237-246.
[doi:10.1556/ComEc.3.2002.2.9](https://doi.org/10.1556/ComEc.3.2002.2.9)

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Dominance metrics:
[`berger()`](https://cmmr.github.io/ecodive/reference/berger.md),
[`dominance()`](https://cmmr.github.io/ecodive/reference/dominance.md),
[`gini()`](https://cmmr.github.io/ecodive/reference/gini.md),
[`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md)

## Examples

``` r
    strong(ex_counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.4913043 0.6950339 0.4742829 0.7934959 
```
