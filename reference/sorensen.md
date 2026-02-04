# Dice-Sorensen dissimilarity

A statistic used for comparing the similarity of two samples.

## Usage

``` r
sorensen(counts, margin = 1L, pairs = NULL, cpus = n_cpus())
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

The Dice-Sorensen dissimilarity is defined as: \$\$\frac{2J}{(A +
B)}\$\$

Where:

- \\A\\, \\B\\ : Number of features in each sample.

- \\J\\ : Number of features in common (intersection).

Base R Equivalent:

    2 * sum(x & y) / sum(x>0, y>0)

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

Sørensen, T. (1948). A method of establishing groups of equal amplitude
in plant sociology based on similarity of species content. *Kongelige
Danske Videnskabernes Selskab, Biologiske Skrifter*, 5, 1-34.

Dice, L. R. (1945). Measures of the amount of ecologic association
between species. *Ecology*, 26(3), 297–302.
[doi:10.2307/1932409](https://doi.org/10.2307/1932409)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md),
[`vignette('bdiv_guide')`](https://cmmr.github.io/ecodive/articles/bdiv_guide.md)

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
