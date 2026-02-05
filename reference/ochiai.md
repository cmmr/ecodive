# Otsuka-Ochiai dissimilarity

Also known as the cosine similarity for binary data.

## Usage

``` r
ochiai(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

    x <- ex_counts[1,]
    y <- ex_counts[2,]
    1 - sum(x & y) / sqrt(sum(x>0) * sum(y>0))

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

Ochiai, A. (1957). Zoogeographic studies on the soleoid fishes found in
Japan and its neighbouring regions. *Bulletin of the Japanese Society of
Scientific Fisheries*, 22, 526-530.

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

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
