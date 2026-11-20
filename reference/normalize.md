# Normalize a Count Table

Applies one of the transformations that
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md), and
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) offer
through their `norm` argument, and returns the result. The values are
exactly those that these functions compute internally.

## Usage

``` r
normalize(counts, norm, margin = 1L, pseudocount = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- norm:

  The transformation:

  - `'percent'`: Relative abundance: each sample's counts divided by
    their sum (so they sum to 1).

  - `'chord'`: Each sample's counts divided by their Euclidean norm (so
    their squares sum to 1).

  - `'hellinger'`: The square root of each count divided by the sample's
    sum (the Hellinger transformation). Euclidean distances between the
    results are Hellinger distances.

  - `'binary'`: Presence/absence: 1 for each non-zero count.

  - `'clr'`: Centered log ratio: \\\log(x + pseudocount)\\ minus the
    sample's mean of those logs, over every feature. Zeros get a value
    too, so the result is dense. See **Pseudocount** below.

  - `'rclr'`: Robust centered log ratio: \\\log(x)\\ minus the mean log
    of the sample's non-zero counts, for non-zero counts only; zeros
    stay 0.

  - `'log'`: Log relative abundance (log-TSS): \\\log\\ of each sample's
    proportions after adding `pseudocount` to every value, \\\log((x +
    pseudocount) / \sum (x + pseudocount))\\. It is the clr shifted by
    one constant per sample. Zeros get a value too, so the result is
    dense. See **Pseudocount** below.

  - `'none'`: No transformation.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

The normalized table, with the same orientation (samples in rows or
columns, as given by `margin`) and dimnames as `counts`, and in the same
form where zeros stay zero:

- A base `matrix` gives a numeric `matrix`.

- A Matrix `dgCMatrix`, `dgTMatrix`, or `dgeMatrix` gives the same
  class, and a slam `simple_triplet_matrix` gives a
  `simple_triplet_matrix`. Zeros are not stored.

- A `phyloseq`, `rbiom`, or `SummarizedExperiment` object gives its
  normalized count matrix (not the object), in that object's
  orientation.

`norm = 'clr'` (and `'log'`) always gives a dense base `matrix`, since
every zero becomes a (usually negative) value, which can take far more
memory than a sparse input.
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md), and
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) accept
`norm = 'clr'` themselves and never make that dense copy, so prefer them
where you can. `normalize()` warns when the result would have more than
\\10^8\\ values (800 MB).

## Details

Samples whose counts are already proportions (summing to 1, within
\\10^{-6}\\) are left as they are by `norm = 'percent'`, as in the
functions with a `norm` argument. Stored zeros in sparse input are
dropped, and repeated entries of a `dgTMatrix` are summed.

**Masking**: BiocGenerics, which is attached along with phyloseq or
SummarizedExperiment, also has a `normalize()`. Whichever package was
attached last masks the other's, so call `ecodive::normalize()` to be
sure of getting this one.

## Pseudocount

The `pseudocount` parameter is only relevant when `norm = 'clr'` or
`norm = 'log'`.

Zeros are undefined in the centered log-ratio (CLR) transformation and
in the log of proportions. If `norm` is `'clr'` or `'log'`,
`pseudocount` is `NULL` (the default), and zeros are detected, the
function uses half the minimum non-zero value (`min(x[x>0]) / 2`) and
issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`), or
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md),
whose result needs no pseudocount.

**Why this matters:** The choice of pseudocount is not neutral; it acts
as a weighting factor that can significantly distort downstream results,
especially for sparse datasets. See Gloor et al. (2017) and Kaul et al.
(2017) for open-access discussions on the mathematical implications, or
Costea et al. (2014) for the impact on community clustering.

See [`aitchison`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## Machine learning

Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
and compositional geometry; don't expect them to improve prediction.
With random forests predicting 65 metadata variables in four datasets,
Yerke et al. (2024) found that alr, clr and ilr (PhILR) "performed
generally slightly worse or only as well as compositionally naïve
transformations", and relative abundances did best, by a small margin.
Compare with `norm = 'percent'` (see `normalize()`) before choosing a
log-ratio predictor table.

## References

Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
normalizations outperform compositional data transformations in machine
learning applications. *Microbiome*, 12, 45.
[doi:10.1186/s40168-023-01747-z](https://doi.org/10.1186/s40168-023-01747-z)

## See also

[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md),
[`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md),
[`qmp_scale()`](https://cmmr.github.io/ecodive/reference/qmp.md)

## Examples

``` r
    # Relative abundances
    normalize(ex_counts, 'percent')[, 1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva   0.469565217 0.005797101     0.000000000 0.521739130
#> Gums     0.895033860 0.004514673     0.000000000 0.098194131
#> Nose     0.021760633 0.001978239     0.492581602 0.001978239
#> Stool    0.001626016 0.993495935     0.001626016 0.001626016

    # Centered log ratios, with a pseudocount of 1 for the zeros
    clr <- normalize(ex_counts, 'clr', pseudocount = 1)
    rowSums(clr)  # 0 for every sample
#>        Saliva          Gums          Nose         Stool 
#> -2.220446e-16 -6.661338e-16  3.552714e-15 -1.110223e-15 

    # The same values as pca()'s own normalization
    a <- pca(ex_counts, norm = 'clr', pseudocount = 1, n_axes = 2)
    b <- pca(clr, n_axes = 2)
    all.equal(a, b)
#> [1] TRUE
```
