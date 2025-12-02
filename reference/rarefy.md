# Rarefy OTU counts.

Sub-sample OTU observations such that all samples have an equal number.
If called on data with non-integer abundances, values will be re-scaled
to integers between 1 and `depth` such that they sum to `depth`.

## Usage

``` r
rarefy(
  counts,
  depth = 0.1,
  n_samples = NULL,
  seed = 0,
  times = NULL,
  drop = TRUE,
  margin = 1L,
  cpus = n_cpus()
)
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

- depth:

  How many observations to keep per sample. When `0 < depth < 1`, it is
  taken as the minimum percentage of the dataset's observations to keep.
  Ignored when `n_samples` is specified. Default: `0.1`

- n_samples:

  The number of samples to keep. When `0 < n_samples < 1`, it is taken
  as the percentage of samples to keep. If negative, that number of
  samples is dropped. If `0`, all samples are kept. If `NULL`, then
  `depth` is used instead. Default: `NULL`

- seed:

  An integer seed for randomizing which observations to keep or drop. If
  you need to create different random rarefactions of the same data, set
  the seed to a different number each time. Default: `0`

- times:

  How many independent rarefactions to perform. If set, `rarefy()` will
  return a list of matrices. The seeds for each matrix will be
  sequential, starting from `seed`. Default: `NULL`

- drop:

  Drop rows and columns with zero observations after rarefying. Default:
  `TRUE`

- margin:

  If your samples are in the matrix's rows, set to `1L`. If your samples
  are in columns, set to `2L`. Ignored when `counts` is a `phyloseq`,
  `rbiom`, `SummarizedExperiment`, or `TreeSummarizedExperiment` object.
  Default: `1L`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A rarefied matrix. `Matrix` and `slam` objects will be returned with the
same type; otherwise a base R `matrix` will be returned.

## Examples

``` r
    # A 4-sample x 5-OTU matrix with samples in rows.
    counts <- matrix(c(0,0,0,0,0,8,9,10,5,5,5,5,2,0,0,0,6,5,7,0), 4, 5,
      dimnames = list(LETTERS[1:4], paste0('OTU', 1:5)))
    counts
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> A    0    0    5    2    6
#> B    0    8    5    0    5
#> C    0    9    5    0    7
#> D    0   10    5    0    0
    rowSums(counts)
#>  A  B  C  D 
#> 13 18 21 15 
    
    # Rarefy all samples to a depth of 13.
    # Note that sample 'A' has 0 counts and is dropped.
    r_mtx <- rarefy(counts, depth = 13, seed = 1)
    r_mtx
#>   OTU2 OTU3 OTU4 OTU5
#> A    0    5    2    6
#> B    4    4    0    5
#> C    3    5    0    5
#> D   10    3    0    0
    rowSums(r_mtx)
#>  A  B  C  D 
#> 13 13 13 13 
    
    # Keep zero-sum rows and columns by setting `drop = FALSE`.
    rarefy(counts, depth = 13, drop = FALSE, seed = 1)
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> A    0    0    5    2    6
#> B    0    4    4    0    5
#> C    0    3    5    0    5
#> D    0   10    3    0    0
    
    # Rarefy to the depth of the 2nd most abundant sample (B, depth=22).
    rarefy(counts, n_samples = 2, seed = 1)
#>   OTU2 OTU3 OTU5
#> B    8    5    5
#> C    6    5    7
    
    # Perform 3 independent rarefactions.
    r_list <- rarefy(counts, depth = 13, times = 3, seed = 1)
    length(r_list)
#> [1] 3
    r_list[[1]]
#>   OTU2 OTU3 OTU4 OTU5
#> A    0    5    2    6
#> B    4    4    0    5
#> C    3    5    0    5
#> D   10    3    0    0
    
    # The class of the input matrix is preserved.
    if (requireNamespace('Matrix', quietly = TRUE)) {
      counts_dgC <- Matrix::Matrix(counts, sparse = TRUE)
      class(counts_dgC)
      r_dgC <- rarefy(counts_dgC, depth = 13, seed = 1)
      class(r_dgC)
    }
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```
