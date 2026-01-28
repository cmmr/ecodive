# Otsuka-Ochiai dissimilarity

Also known as the cosine similarity for binary data.

## Usage

``` r
ochiai(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Otsuka-Ochiai dissimilarity is defined as: \$\$1 -
\frac{J}{\sqrt{AB}}\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

Base R Equivalent:

    1 - sum(x & y) / sqrt(sum(x>0) * sum(y>0))

## References

Ochiai, A. (1957). Zoogeographic studies on the soleoid fishes found in
Japan and its neighbouring regions. *Bulletin of the Japanese Society of
Scientific Fisheries*, 22, 526-530.

## See also

beta_div

Other Presence/Absence metrics:
[`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md),
[`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)

## Examples

``` r
    ochiai(ex_counts)
#>           Saliva       Gums       Nose
#> Gums  0.10557281                      
#> Nose  0.18350342 0.08712907           
#> Stool 0.32917961 0.20000000 0.08712907
```
