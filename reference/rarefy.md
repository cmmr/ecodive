# Rarefy Observation Counts

Sub-sample observations from a feature table such that all samples have
the same library size (depth). This is performed via random sampling
without replacement.

## Usage

``` r
rarefy(
  counts,
  depth = NULL,
  seed = 0,
  times = NULL,
  drop = TRUE,
  margin = 1L,
  cpus = n_cpus(),
  warn = interactive()
)
```

## Arguments

- counts:

  A numeric matrix or sparse matrix object (e.g., `dgCMatrix`). Counts
  must be integers (non-integer counts will be cast to integers).

- depth:

  The number of observations to keep per sample. If `NULL` (the
  default), a depth is auto-selected to maximize data retention.

- seed:

  An integer seed for the random number generator. Providing the same
  seed guarantees reproducible results. Default: `0`

- times:

  The number of independent rarefactions to perform. If set, returns a
  list of matrices. Seeds for subsequent iterations are sequential
  (`seed`, `seed + 1`, ...). Default: `NULL`

- drop:

  Logical. If `TRUE`, samples with fewer than `depth` observations are
  discarded. If `FALSE`, they are kept with their original counts.
  Default: `TRUE`

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

- warn:

  Logical. If `TRUE`, emits a warning when samples are dropped or
  returned unrarefied due to insufficient depth. Default:
  [`interactive()`](https://rdrr.io/r/base/interactive.html)

## Value

A rarefied matrix. The output class (`matrix`, `dgCMatrix`, etc.)
matches the input class.

## Details

**Auto-Depth Selection**  
If `depth` is `NULL`, the function defaults to the highest depth that
retains at least 10% of the total observations in the dataset.

**Dropping vs. Retaining Samples**  
If a sample has fewer observations than the specified `depth`:

- `drop = TRUE` (Default): The sample is removed from the output matrix.

- `drop = FALSE`: The sample is returned **unmodified** (with its
  original counts). It is *not* rarefied or zeroed out.

**Zero-Sum Features**  
Features (OTUs, ASVs, Genes) that lose all observations during
rarefaction are **always retained** as columns/rows of zeros. This
ensures the output matrix dimensions remain consistent with the input
(barring dropped samples).

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
    # Sample 'A' (0 counts) and 'D' (12 counts) will be dropped.
    r_mtx <- rarefy(counts, depth = 13, seed = 1)
    r_mtx
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> A    0    0    5    2    6
#> B    0    4    4    0    5
#> C    0    3    5    0    5
#> D    0   10    3    0    0
    rowSums(r_mtx)
#>  A  B  C  D 
#> 13 13 13 13 
    
    # Keep under-sampled samples by setting `drop = FALSE`.
    # Samples 'A' and 'D' are returned with their original counts.
    rarefy(counts, depth = 13, drop = FALSE, seed = 1)
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> A    0    0    5    2    6
#> B    0    4    4    0    5
#> C    0    3    5    0    5
#> D    0   10    3    0    0
    
    # Perform 3 independent rarefactions.
    r_list <- rarefy(counts, depth = 13, times = 3, seed = 1)
    length(r_list)
#> [1] 3
    
    # Sparse matrices are supported and their class is preserved.
    if (requireNamespace('Matrix', quietly = TRUE)) {
      counts_dgC <- Matrix::Matrix(counts, sparse = TRUE)
      class(rarefy(counts_dgC, depth = 13))
    }
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```
