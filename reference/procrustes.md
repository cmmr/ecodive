# Procrustes Rotation of One Ordination onto Another

Rotates, reflects, and optionally dilates (uniformly scales) one
ordination, `ord`, to fit another of the same samples, `target`, as
closely as possible, and reports how close the fit is: the Procrustes
sum of squares, \\m^2\\. Use it to ask how similar two pictures of the
same samples are: the same samples under two beta diversity metrics,
say, or from 16S and from shotgun metagenomics. Gives the same results
as `vegan::procrustes()`.
[`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
tests whether the fit is better than chance (PROTEST).

## Usage

``` r
procrustes(
  ord,
  target,
  axes = NULL,
  symmetric = FALSE,
  dilate = TRUE,
  n_axes = 3L,
  cpus = n_cpus()
)
```

## Arguments

- ord:

  The ordination to rotate onto `target`: a numeric matrix (or data
  frame) with one row per sample and one column per axis, such as from
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md), or
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md); or a
  `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  which is first ordinated with
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).

- target:

  The ordination that `ord` is fitted to, in any of the forms `ord`
  takes. When both `ord` and `target` have sample (row) names,
  `target`'s rows are matched to `ord`'s by name; otherwise they must be
  in the same order.

- axes:

  Which axes (columns) of `ord` and `target` to use: whole numbers, as
  in [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md).
  When one of the two has fewer columns than asked for, it is padded
  with columns of zeros, as in `vegan::procrustes()`. Default: `NULL`,
  the first `min(ncol(ord), ncol(target))` axes.

- symmetric:

  Use the symmetric Procrustes statistic: scale both `ord` and `target`
  to unit sum of squares first, so that \\m^2\\ is between 0 and 1 and
  is the same whichever of the two is rotated. As
  `vegan::procrustes()`'s `symmetric`. Default: `FALSE`.

- dilate:

  Allow a uniform dilation (scaling) of `ord`, as well as the rotation
  and reflection. `vegan::procrustes()`'s `scale`, renamed because
  `scale` means "standardise each feature" elsewhere in ecodive.
  Default: `TRUE`.

- n_axes:

  When `ord` or `target` is a `dist` object, the number of PCoA axes to
  compute from it, as in
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md). Ignored
  otherwise. Default: `3`.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

The rotated `ord`: a numeric matrix of class `ecodive_ord`, with one row
per sample and one column per axis (named for `target`'s axes), as
`vegan::procrustes()`'s `Yrot`. It is centred, and in `target`'s frame:
`target`'s centred coordinates (scaled to unit sum of squares, if
`symmetric`) are attribute `target`. Its attributes:

- `method`: `'procrustes'`.

- `eig`, `prop`: `NA` for each axis, as Procrustes rotation has no
  eigenvalues.

- `ss`: The Procrustes sum of squares, \\m^2\\: the sum of squared
  distances between `target` and the rotated `ord`. vegan's `ss`.

- `rotation`: The \\k \times k\\ orthogonal matrix that rotates (and
  possibly reflects) the centred `ord` onto `target`. vegan's
  `rotation`.

- `translation`: The shift that takes the centroid of the rotated and
  dilated `ord` to that of `target`, one value per axis, so that
  `dilation * ord %*% rotation + translation` (adding `translation` to
  each row of the uncentred `ord`) gives the fitted coordinates. Zero
  when `symmetric`, as both are centred first. vegan's `translation`.

- `dilation`: The dilation factor applied to `ord`: 1 when
  `dilate = FALSE`. vegan's `scale`.

- `residuals`: Each sample's distance between `target` and the rotated
  `ord`, as [`residuals()`](https://rdrr.io/r/stats/residuals.html) of a
  `vegan::procrustes()` result.

- `symmetric`: Whether the symmetric statistic was used.

- `target`: `target`'s coordinates as fitted: centred, scaled to unit
  sum of squares if `symmetric`, and padded with zeros to the same axes.
  vegan's `X`.

## Details

**Method**

Procrustes analysis (Gower 1975; Gower & Dijksterhuis 2004) centres
`ord` (\\Y\\) and `target` (\\X\\), both \\n \times k\\, and finds the
orthogonal matrix \\A\\ (a rotation, possibly with a reflection) and the
dilation \\c\\ that minimise \$\$m^2 = \\ X - c Y A \\^2\$\$ (the sum of
squared differences). With \\X^\top Y = U \Sigma V^\top\\ (singular
value decomposition), \\A = V U^\top\\, as in SciPy's
`orthogonal_procrustes()`, and \\c = \mathrm{tr}(\Sigma) /
\mathrm{tr}(Y^\top Y)\\, or 1 when `dilate = FALSE`.

The ordinary (non-symmetric) \\m^2\\ depends on the scales of `ord` and
`target`, and on which one is rotated. With `symmetric = TRUE`, both are
first scaled to unit sum of squares, so that \\m^2 = 1 -
\mathrm{tr}(\Sigma)^2\\ lies between 0 and 1 and is the same either way
round: SciPy's `scipy.spatial.procrustes()` "disparity".
[`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
reports \\r = \sqrt{1 - m^2}\\, this symmetric Procrustes correlation.

Rotation can't make up for a different number of axes: compare the same
number from each ordination, as the default `axes` does.

**Differences from vegan::procrustes()**

- The argument order: `procrustes(ord, target)` rotates `ord` onto
  `target`, as `vegan::procrustes(target, ord)` (vegan's
  `procrustes(X, Y)` rotates `Y` onto `X`).

- `dilate` is vegan's `scale`.

- By default, the axes both have in common are used; vegan pads the one
  with fewer axes with zeros. Give `axes` to do the same.

- `ord` and `target` may be `dist` objects, and their samples are
  matched by name.

- The result is the rotated ordination itself, with the rest as
  attributes. There are no
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) methods.

**Same names in other packages**

`vegan::procrustes()` takes its arguments the other way round (see
above). When both packages are attached, the one attached last masks the
other: call `ecodive::procrustes()` to be sure.

## References

Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
40(1), 33-51.
[doi:10.1007/BF02291478](https://doi.org/10.1007/BF02291478)

Gower, J. C., & Dijksterhuis, G. B. (2004). *Procrustes Problems*.
Oxford University Press.

Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate
data sets match? The advantages of a Procrustean superimposition
approach over the Mantel test. *Oecologia*, 129(2), 169-178.
[doi:10.1007/s004420100720](https://doi.org/10.1007/s004420100720)

## See also

[`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md),
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)

## Examples

``` r
    # The same samples under two beta diversity metrics: 30 samples along a gradient that favours some genera at one end
    # and others at the other.
    set.seed(1)
    grad   <- seq(-2, 2, length.out = 30)
    counts <- sapply(1:20, function (j) rpois(30, 5 * exp(grad * (j - 10.5) / 5)))
    rownames(counts) <- paste0('S', 1:30)

    bc <- pcoa(bray(counts), n_axes = 2)
    jc <- pcoa(jaccard(counts), n_axes = 2)

    fit <- procrustes(jc, bc)
    attr(fit, 'ss')
#> [1] 2.520756
    attr(fit, 'residuals')
#>         S1         S2         S3         S4         S5         S6         S7 
#> 0.34932591 0.36099692 0.36411311 0.27373695 0.23874099 0.33038984 0.12135288 
#>         S8         S9        S10        S11        S12        S13        S14 
#> 0.08004099 0.30050059 0.30933607 0.22656374 0.28374082 0.29777517 0.28249482 
#>        S15        S16        S17        S18        S19        S20        S21 
#> 0.29095713 0.26651277 0.27342894 0.27743497 0.37654038 0.30411042 0.18633496 
#>        S22        S23        S24        S25        S26        S27        S28 
#> 0.24588617 0.36725760 0.23591621 0.18771922 0.23811780 0.60565469 0.23559225 
#>        S29        S30 
#> 0.18422628 0.10820729 

    # Symmetric: m^2 between 0 and 1.
    attr(procrustes(jc, bc, symmetric = TRUE), 'ss')
#> [1] 0.49545

    # Distances are ordinated with pcoa() first.
    attr(procrustes(jaccard(counts), bray(counts), n_axes = 2), 'ss')
#> [1] 2.520756
```
