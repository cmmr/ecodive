# Dice-Sorensen dissimilarity

A statistic used for comparing the similarity of two samples.

## Usage

``` r
sorensen(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Dice-Sorensen dissimilarity is defined as: \$\$\frac{2J}{(A +
B)}\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

Base R Equivalent:

    2 * sum(x & y) / sum(x>0, y>0)

## References

Sørensen, T. (1948). A method of establishing groups of equal amplitude
in plant sociology based on similarity of species content. *Kongelige
Danske Videnskabernes Selskab, Biologiske Skrifter*, 5, 1-34.

## See also

beta_div

Other Presence/Absence metrics:
[`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md),
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
[`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md)

## Examples

``` r
    sorensen(ex_counts)
#>           Saliva       Gums       Nose
#> Gums  0.11111111                      
#> Nose  0.20000000 0.09090909           
#> Stool 0.33333333 0.20000000 0.09090909
```
