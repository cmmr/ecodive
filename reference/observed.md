# Observed Features

The count of unique features (richness) in a sample.

## Usage

``` r
observed(counts, margin = 1L, cpus = n_cpus())
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

Observed features is defined simply as the number of features with
non-zero abundance: \$\$n\$\$

**Base R Equivalent:**

    sum(x > 0)

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

## See also

[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
[`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

Other Richness metrics:
[`ace()`](https://cmmr.github.io/ecodive/reference/ace.md),
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md),
[`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md),
[`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md),
[`squares()`](https://cmmr.github.io/ecodive/reference/squares.md)

## Examples

``` r
    observed(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>      4      5      6      5 
```
