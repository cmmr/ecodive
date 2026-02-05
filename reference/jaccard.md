# Jaccard distance

Measures dissimilarity between sample sets.

## Usage

``` r
jaccard(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Jaccard distance is defined as: \$\$1 - \frac{J}{(A + B - J)}\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

Base R Equivalent:

    x <- ex_counts[1,]
    y <- ex_counts[2,]
    1 - sum(x & y) / sum(x | y)

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

Jaccard, P. (1912). The distribution of the flora in the alpine zone.
*New Phytologist*, 11(2), 37-50.
[doi:10.1111/j.1469-8137.1912.tb05611.x](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)

Jaccard, P. (1908). Nouvelles recherches sur la distribution florale.
*Bulletin de la Societe Vaudoise des Sciences Naturelles*, 44(163),
223-270.
[doi:10.5169/seals-268384](https://doi.org/10.5169/seals-268384)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

Other Presence/Absence metrics:
[`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md),
[`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md),
[`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md)

## Examples

``` r
    jaccard(ex_counts)
#>          Saliva      Gums      Nose
#> Gums  0.2000000                    
#> Nose  0.3333333 0.1666667          
#> Stool 0.5000000 0.3333333 0.1666667
```
