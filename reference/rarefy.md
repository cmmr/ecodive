# Rarefy Observation Counts

Sub-sample observations from a feature table such that all samples have
the same library size (depth). This is performed via random sampling
without replacement.

## Usage

``` r
rarefy(
  counts,
  depth = NULL,
  seed = 0L,
  times = NULL,
  drop = TRUE,
  margin = 1L,
  cpus = n_cpus(),
  warn = interactive(),
  coverage = NULL
)
```

## Arguments

- counts:

  A numeric matrix or sparse matrix object (e.g., `dgCMatrix`). Counts
  must be integers.

- depth:

  The number of observations to keep per sample: one number for every
  sample, or one per sample (matched to the samples by name when both
  are named, else in the samples' order), as
  [`qmp_rarefy()`](https://cmmr.github.io/ecodive/reference/qmp.md)
  uses. If `NULL` (the default), a depth is auto-selected to maximize
  data retention.

- seed:

  An integer seed for the random number generator. Providing the same
  seed guarantees reproducible results. Default: `0`

- times:

  The number of independent rarefactions to perform. If set, returns a
  list of matrices. Seeds for subsequent iterations are sequential
  (`seed`, `seed + 1`, ...). Default: `NULL`

- drop:

  Logical. If `TRUE`, samples with fewer than `depth` observations
  (their own depth, when there is one per sample) are discarded. If
  `FALSE`, they are kept with their original counts. Default: `TRUE`

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

- warn:

  Logical. If `TRUE`, emits a warning when samples are dropped or
  returned unrarefied due to insufficient depth. Default:
  [`interactive()`](https://rdrr.io/r/base/interactive.html)

- coverage:

  Rarefy to equal coverage instead of equal depth: a sample coverage
  greater than 0 and at most 1, such as `0.95`. Each sample is rarefied
  to the smallest depth whose expected coverage reaches it (see
  "Coverage-Based Rarefaction"). Give `depth` or `coverage`, not both.
  Default: `NULL` (rarefy to `depth`).

## Value

A rarefied matrix. The output class (`matrix`, `dgCMatrix`, etc.)
matches the input class.

## Auto-Depth Selection

If `depth` is `NULL`, the depth is the smallest library size at which
the observations kept (that depth times the number of samples with at
least that many) make up at least 10% of all the observations in the
table.
[`beta_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md)
and
[`alpha_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md)
choose their default depth the same way.

## Coverage-Based Rarefaction

Chao & Jost (2012) compare samples at equal **coverage** (completeness:
the share of the community's individuals that belong to the features
seen) rather than at equal depth. With `coverage`, each sample is
rarefied to the smallest depth \\d\\ whose expected coverage \\C_d\\ (as
in
[`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md)'s
`coverage` column) is at least `coverage`. At the sample's own total,
\\C_d\\ is the estimated coverage of the whole sample (Chao et al. 2014,
Eq. 12). Samples whose whole-sample coverage is below `coverage` are
dropped (or, with `drop = FALSE`, returned unrarefied). The depths come
from the counts alone, so they are the same for every `seed`; each
sample is then rarefied as with one `depth` per sample.
[`rowSums()`](https://rdrr.io/r/base/colSums.html) of the result gives
the depths chosen.

Any metric can then be compared at equal coverage, not only richness or
Hill numbers.

## One Depth per Sample

With one `depth` per sample, each sample is rarefied to its own depth,
and gives exactly what rarefying the table to that sample's depth would
give it, with the same `seed`. `drop` and `warn` compare each sample
with its own depth.

## Dropping vs. Retaining Samples

If a sample has fewer observations than the specified `depth`:

- `drop = TRUE` (Default): The sample is removed from the output matrix.

- `drop = FALSE`: The sample is returned **unmodified** (with its
  original counts). It is *not* rarefied or zeroed out.

## Zero-Sum Features

Features (OTUs, ASVs, Genes) that lose all observations during
rarefaction are **always retained** as columns/rows of zeros. This
ensures the output matrix dimensions remain consistent with the input
(barring dropped samples).

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

Chao, A., & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
*Ecology*, 93, 2533–2547.
[doi:10.1890/11-1952.1](https://doi.org/10.1890/11-1952.1)

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H.,
Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation
with Hill numbers: a framework for sampling and estimation in species
diversity studies. *Ecological Monographs*, 84, 45–67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1)

## See also

[`qmp_rarefy()`](https://cmmr.github.io/ecodive/reference/qmp.md), which
rarefies each sample to a depth set by its measured microbial load;
[`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md),
to choose a depth or a coverage.

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
    
    # Rarefy all samples to a depth of 18.
    # Sample 'A' (13 counts) and 'D' (15 counts) will be dropped.
    r_mtx <- rarefy(counts, depth = 18)
    r_mtx
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> B    0    8    5    0    5
#> C    0    8    4    0    6
    rowSums(r_mtx)
#>  B  C 
#> 18 18 
    
    # Keep under-sampled samples by setting `drop = FALSE`.
    # Samples 'A' and 'D' are returned with their original counts.
    r_mtx <- rarefy(counts, depth = 18, drop = FALSE)
    r_mtx
#>   OTU1 OTU2 OTU3 OTU4 OTU5
#> A    0    0    5    2    6
#> B    0    8    5    0    5
#> C    0    8    4    0    6
#> D    0   10    5    0    0
    rowSums(r_mtx)
#>  A  B  C  D 
#> 13 18 18 15 
    
    # One depth per sample, matched by name.
    r_mtx <- rarefy(counts, depth = c(A = 10, B = 15, C = 20, D = 12))
    rowSums(r_mtx)
#>  A  B  C  D 
#> 10 15 20 12 
    
    # Rarefy each sample to an expected coverage of 80%.
    r_mtx <- rarefy(counts, coverage = 0.8)
    rowSums(r_mtx)
#> A B C D 
#> 4 5 5 3 
    
    # Perform 3 independent rarefactions.
    r_list <- rarefy(counts, times = 3)
    length(r_list)
#> [1] 3
    
    # Sparse matrices are supported and their class is preserved.
    if (requireNamespace('Matrix', quietly = TRUE)) {
      counts_dgC <- Matrix::Matrix(counts, sparse = TRUE)
      str(rarefy(counts_dgC))
    }
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:11] 1 2 3 0 1 2 3 0 0 1 ...
#>   ..@ p       : int [1:6] 0 0 3 7 8 11
#>   ..@ Dim     : int [1:2] 4 5
#>   ..@ Dimnames:List of 2
#>   .. ..$ : chr [1:4] "A" "B" "C" "D"
#>   .. ..$ : chr [1:5] "OTU1" "OTU2" "OTU3" "OTU4" ...
#>   ..@ x       : num [1:11] 8 6 9 5 2 3 4 2 6 3 ...
#>   ..@ factors : list()
```
