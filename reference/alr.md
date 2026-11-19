# Additive Log-Ratio Transform

The log of every feature's abundance relative to one chosen feature, the
`denominator`: Aitchison's additive log ratio (alr).

## Usage

``` r
alr(counts, denominator, margin = 1L, pseudocount = NULL, cpus = n_cpus())
```

## Arguments

- counts:

  A numeric matrix of non-negative abundances (samples \\\times\\
  features), or any of the objects listed under **Input Types**. Counts,
  proportions, and the output of a zero replacement are all accepted:
  each sample's scale cancels out of the log ratios.

- denominator:

  The feature to divide by: its name, or its (1-based) index among the
  features, or a logical vector with one `TRUE`. Required. Choose a
  feature that is present in every sample and abundant and stable across
  them, since every coordinate is relative to it.

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
single sample given as a vector). Sample names are kept. Each coordinate
is named after its numerator feature, in the order of `counts`; features
without names give `alr1`, `alr2`, ..., numbered by their position in
`counts` (so `denominator = 2` of 4 gives `alr1`, `alr3`, `alr4`).
Attributes:

- `method`: `'alr'`.

- `denominator`: the denominator feature's name, or its index when the
  features have no names.

- `pseudocount`: the pseudocount used (`0` for none).

## Details

For a sample \\x\\ of \\D\\ features and denominator \\x_d\\,
\$\$\mathrm{alr}(x) = \left\[\log\frac{x_j + c}{x_d + c}\right\]\_{j
\neq d}\$\$ with \\c\\ the pseudocount (natural logs). The coordinates
are easy to interpret (everything relative to one feature) and are what
multinomial regression models use, but the transform is not an isometry:
Euclidean distances between alr coordinates are not Aitchison distances.
For distances, use
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md);
for orthonormal coordinates,
[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md). The centered
log ratio is `normalize(counts, 'clr')` (or `norm = 'clr'` in any
function that takes `norm`).

`alr(x + c, d)` here equals scikit-bio's `alr(x + c, ref_idx = d - 1)`
exactly. scikit-bio defaults to the first feature and
`compositions::alr()` to the last: ecodive has no default, because the
choice shapes every coordinate.

The result feeds straight into
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
[`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md),
which take negative values with `norm = 'none'`.

**Same names in other packages**: `compositions::alr()` (the last part
as denominator) and `robCompositions::addLR()` compute the same
transform. Whichever package was attached last masks the other's
`alr()`; write `ecodive::alr()` to be sure.

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

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman & Hall.

## See also

[`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md),
[`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
(`'clr'` and `'rclr'`),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)

## Examples

``` r
    # Relative to the first feature
    x <- alr(ex_counts, denominator = 1, pseudocount = 1)
    x[1:3, 1:3]
#>        Bacteroides Corynebacterium Haemophilus
#> Saliva   -3.995138       -5.093750   0.1047468
#> Gums     -5.067646       -6.677083  -2.1997466
#> Nose     -2.036882        3.077112  -2.0368819
    attr(x, 'denominator')
#> [1] "Streptococcus"

    # The same as the log ratios themselves
    all.equal(
      unname(x[, 1]),
      unname(log((ex_counts[, 2] + 1) / (ex_counts[, 1] + 1))) )
#> [1] TRUE
```
