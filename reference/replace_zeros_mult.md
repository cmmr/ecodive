# Multiplicative Zero Replacement

Replaces every zero with a small proportion \\\delta\\, and scales each
sample's non-zero proportions down by the same factor so that the sample
still sums to 1 (Martín-Fernández et al. 2003). The ratios between
non-zero parts are unchanged.

## Usage

``` r
replace_zeros_mult(counts, margin = 1L, delta = NULL)
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts
  or proportions.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- delta:

  The proportion that every zero becomes: a single value, or one value
  per sample (matched to the samples by name when both are named, else
  in order), each between 0 and 1. `NULL` (the default) uses \\1 /
  D^2\\, with \\D\\ the number of features, as scikit-bio does. A
  sample's \\\delta\\ times its number of zeros must be below 1.

## Value

A dense numeric `matrix` of proportions with no zeros, in the
orientation of `counts` (samples \\\times\\ features for `margin = 1`,
features \\\times\\ samples for `margin = 2`, and for a single sample
given as a vector), with its names. Each sample sums to 1. Attributes:

- `method`: `'mult'`.

- `delta`: the \\\delta\\ used: one value, or one per sample (named
  after the samples).

- `n_replaced`: the number of zeros replaced.

## Details

For a sample with proportions \\x_j = c_j / \sum_k c_k\\ and \\Z\\
zeros, \$\$r_j = \delta \text{ if } x_j = 0, \qquad r_j = (1 -
Z\delta)\\ x_j \text{ otherwise.}\$\$

The \\\delta\\ is the same for every sample (unless given per sample),
whatever its depth: a zero in a sample of 1,000 reads and one in a
sample of 100,000 become the same proportion. Choose \\\delta\\ below
the smallest proportion that the samples could detect, e.g. a fraction
of `1 / max(depths)`, or give each sample its own, e.g.
`delta = 0.5 / rowSums(counts)` (half a read).

A feature that is zero in every sample is replaced like any other zero.
Samples with no non-zero values are an error: there is nothing to scale.

The result equals scikit-bio's `multi_replace()` (version 0.7.4) to the
last bit: each sample's total is summed in numpy's (pairwise) order.
Unlike scikit-bio, a \\\delta\\ that would leave a sample's non-zero
proportions at exactly 0 (\\Z\delta = 1\\) is an error too, so that the
result has no zeros.

The result feeds straight into
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) and any
function with `norm = 'clr'`, with no pseudocount.

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

Martín-Fernández, J. A., Barceló-Vidal, C., & Pawlowsky-Glahn, V.
(2003). Dealing with zeros and missing values in compositional data sets
using nonparametric imputation. *Mathematical Geology*, 35, 253–278.
[doi:10.1023/A:1023866030544](https://doi.org/10.1023/A%3A1023866030544)

## See also

[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
[`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)

## Examples

``` r
    x <- replace_zeros_mult(ex_counts)
    x[1:3, 1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva    0.44347826 0.005475040      0.02777778 0.492753623
#> Gums      0.87017181 0.004389265      0.02777778 0.095466516
#> Nose      0.02176063 0.001978239      0.49258160 0.001978239
    attr(x, 'delta')
#> [1] 0.02777778
    rowSums(x)
#> Saliva   Gums   Nose  Stool 
#>      1      1      1      1 

    # The ratios between non-zero parts are kept.
    nz <- ex_counts[1, ] > 0
    all.equal(x[1, nz] / sum(x[1, nz]), ex_counts[1, nz] / sum(ex_counts[1, nz]))
#> [1] TRUE

    # Half a read in each sample.
    x <- replace_zeros_mult(ex_counts, delta = 0.5 / rowSums(ex_counts))
    attr(x, 'delta')
#>       Saliva         Gums         Nose        Stool 
#> 0.0014492754 0.0005643341 0.0004945598 0.0008130081 
```
