# Isometric Log-Ratio Transform

Coordinates of each sample in an orthonormal basis of the Aitchison
simplex: Egozcue et al.'s isometric log ratio (ilr). Euclidean distances
between them are Aitchison distances, for any basis.

## Usage

``` r
ilr(counts, basis = NULL, margin = 1L, pseudocount = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts,
  proportions, and the output of a zero replacement are all accepted:
  each sample's scale cancels out of the log ratios.

- basis:

  The ilr basis, one of:

  - `NULL` (the default): the Helmert basis that scikit-bio's `ilr()`
    uses (built by Gram-Schmidt orthogonalization). Coordinate \\j\\ is
    the balance of the first \\j\\ features against feature \\j + 1\\:
    row \\j\\ has \\\frac{1}{j}\sqrt{j/(j+1)}\\ in columns \\1, \dots,
    j\\, \\-\sqrt{j/(j+1)}\\ in column \\j + 1\\, and zeros after.

  - A **sequential binary partition** (SBP): a \\(D - 1) \times D\\
    matrix of `1`, `-1` and `0`, one row per balance, `1` marking the
    features in its numerator and `-1` those in its denominator. Each
    row's \\r\\ entries of `1` become \\\sqrt{s / (r (r + s))}\\ and its
    \\s\\ entries of `-1` become \\-\sqrt{r / (s (r + s))}\\, as
    scikit-bio's `sbp_basis()`. Every row needs both signs, and the
    partition must be sequential (its balances orthogonal).

  - An **orthonormal basis** in clr coordinates: a numeric \\(D - 1)
    \times D\\ matrix with orthonormal rows that each sum to 0 (checked
    to \\10^{-6}\\, as scikit-bio checks orthonormality).

  A matrix has one row per coordinate and one column per feature
  (scikit-bio's orientation). Its columns are matched to the features by
  name when both have names, and by position otherwise. Its row names
  name the coordinates.

- margin:

  The margin containing samples. `1` if samples are rows, `2` if samples
  are columns. Ignored when `counts` is a special object class (e.g.
  `phyloseq`). Default: `1`

- pseudocount:

  Added to every value before taking logs. `NULL` (the default) adds
  nothing when `counts` has no zeros, and otherwise uses half the
  smallest non-zero value, with a warning. See **Pseudocount**.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A dense numeric `matrix` of the \\D - 1\\ coordinates, in the
orientation of `counts`: samples \\\times\\ coordinates for
`margin = 1`, coordinates \\\times\\ samples for `margin = 2` (and for a
single sample given as a vector). Sample names are kept. Coordinates are
named by `rownames(basis)`, or `ilr1`, `ilr2`, .... Attributes:

- `method`: `'ilr'`.

- `basis`: the \\(D - 1) \times D\\ basis used (an SBP converted to its
  basis), with its columns in the order of the features of `counts`, so
  that `coords %*% basis` gives back `normalize(counts, 'clr')` (for
  `margin = 1`). For the default basis of more than 1000 features, whose
  matrix would take more than 8 MB, it is the string `'helmert'`
  instead.

- `pseudocount`: the pseudocount used (`0` for none).

## Details

For a sample \\x\\ with clr coordinates \\\mathrm{clr}(x) = \log(x +
c) - \overline{\log(x + c)}\\ (natural logs, \\c\\ the pseudocount) and
a basis \\\Psi\\ with orthonormal rows, \$\$\mathrm{ilr}(x) =
\mathrm{clr}(x)\\\Psi^T.\$\$ The default basis is applied by a running
sum, in \\O(D)\\ per sample, and other bases through their non-zero
entries, so an SBP costs only its non-zeros.

Because the coordinates are orthonormal, `euclidean(ilr(x))` is
`aitchison(x)` (with the same pseudocount), whichever basis is used:
**for distances, use
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)**.
`ilr()` is for analyses of the coordinates themselves: regression,
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) or
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) without the
clr's singular covariance, and balances that can be interpreted, from an
SBP. The coordinates feed straight into
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
[`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md)
with `norm = 'none'`.

`ilr(x + c)` here equals scikit-bio's `ilr(x + c)` (and
`ilr(x + c, sbp_basis(sbp))`) to rounding.

**Same names in other packages**: `compositions::ilr()` and
`robCompositions::pivotCoord()` compute ilr coordinates with other
default bases (`compositions::ilr()`'s default is this one with the
opposite sign). Whichever package was attached last masks the other's
`ilr()`; write `ecodive::ilr()` to be sure.

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

## Pseudocount

Logs of zero are undefined, so zeros need replacing first. The
pseudocount is added to every value (not only to zeros), as in
`normalize(norm = 'clr')`. With `pseudocount = NULL` (the default), no
pseudocount is used when `counts` has no zeros; otherwise half the
smallest non-zero value is, with a warning. To avoid a pseudocount,
replace the zeros first with
[`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md):
its result has none.

The choice of pseudocount is not neutral: it can distort the results,
especially for sparse tables. See
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
for references.

## Machine learning

Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
and compositional geometry; don't expect them to improve prediction.
With random forests predicting 65 metadata variables in four datasets,
Yerke et al. (2024) found that alr, clr and ilr (PhILR) "performed
generally slightly worse or only as well as compositionally naïve
transformations", and relative abundances did best, by a small margin.
Compare with `norm = 'percent'` (see
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md))
before choosing a log-ratio predictor table.

## References

Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
Barceló-Vidal, C. (2003). Isometric logratio transformations for
compositional data analysis. *Mathematical Geology*, 35, 279–300.

Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
balances in compositional data analysis. *Mathematical Geology*, 37,
795–828.

Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
normalizations outperform compositional data transformations in machine
learning applications. *Microbiome*, 12, 45.
[doi:10.1186/s40168-023-01747-z](https://doi.org/10.1186/s40168-023-01747-z)

## See also

[`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
(`'clr'` and `'rclr'`),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)

## Examples

``` r
    x <- ilr(ex_counts, pseudocount = 1)
    dim(x)
#> [1] 4 5

    # Euclidean distances on ilr coordinates are Aitchison distances
    all.equal(
      as.vector(euclidean(x, norm = 'none')),
      as.vector(aitchison(ex_counts, pseudocount = 1)) )
#> [1] TRUE

    # Balances from a sequential binary partition of 4 features
    sbp <- rbind(
      c(1,  1, -1, -1),
      c(1, -1,  0,  0),
      c(0,  0,  1, -1) )
    ilr(ex_counts[, 1:4], basis = sbp, pseudocount = 1)[1:3, ]
#>              ilr1     ilr2      ilr3
#> Saliva  0.4969327 2.824989 -3.675893
#> Gums    1.9045923 3.583367 -3.165955
#> Nose   -1.5385559 1.440293  3.616140
```
