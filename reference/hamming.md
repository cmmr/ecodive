# Hamming distance

Measures the minimum number of substitutions required to change one
string into the other.

## Usage

``` r
hamming(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data where each column is a feature, and
  each row is a sample. Any object coercible with
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) can be given here,
  as well as `phyloseq`, `rbiom`, `SummarizedExperiment`, and
  `TreeSummarizedExperiment` objects. For optimal performance with very
  large datasets, see the guide in
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md).

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a `phyloseq`,
  `rbiom`, `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
  Default: `1L`

- pairs:

  Which combinations of samples should distances be calculated for? The
  default value (`NULL`) calculates all-vs-all. Provide a numeric or
  logical vector specifying positions in the distance matrix to
  calculate. See examples.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Details

The Hamming distance is defined as: \$\$(A + B) - 2J\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

Base R Equivalent:

    sum(xor(x, y))

## References

Hamming, R. W. (1950). Error detecting and error correcting codes. *Bell
System Technical Journal*, 29(2), 147-160.
[doi:10.1002/j.1538-7305.1950.tb00463.x](https://doi.org/10.1002/j.1538-7305.1950.tb00463.x)

## See also

beta_div

Other Presence/Absence metrics:
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
[`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md),
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)

## Examples

``` r
    hamming(ex_counts)
#>       Saliva Gums Nose
#> Gums       1          
#> Nose       2    1     
#> Stool      3    2    1
```
