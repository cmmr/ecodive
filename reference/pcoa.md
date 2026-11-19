# Principal Coordinates Analysis (PCoA)

Ordinates samples from their distances, for plotting beta diversity.
Also known as classical or metric multidimensional scaling. Gives the
same coordinates and eigenvalues as
[`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html),
`ape::pcoa()`, and `vegan::wcmdscale()`, but computes only the first few
axes, which takes far less time and memory for large datasets.

## Usage

``` r
pcoa(dm, n_axes = 3L, correction = "none", cpus = n_cpus())
```

## Arguments

- dm:

  A `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  or [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), or a
  symmetric numeric matrix of distances. Must not contain `NA`s or
  negative values.

- n_axes:

  The number of axes to compute. Default: `3`. Reduced, with a warning,
  to the number of axes with positive eigenvalues, which is at most
  \\n - 1\\ for \\n\\ samples.

- correction:

  How to deal with negative eigenvalues: `'none'` (the default),
  `'lingoes'`, or `'cailliez'`. See "Negative eigenvalue corrections"
  below.

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A numeric matrix with one row per sample and one column per axis, named
`"PCoA1"`, `"PCoA2"`, etc., as described in "Ordination result" below.
Its attributes are:

- `method`: `"pcoa"`.

- `eig`: The eigenvalue of each axis: of Gower's centered matrix, so on
  the scale of sums of squares.
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s `eig` are
  variances, these divided by \\n - 1\\;
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)'s are
  too when any distance is 4 or more.

- `prop`: The proportion of the total variance that each axis explains,
  `eig / trace`.

- `trace`: The total variance: the sum of all eigenvalues, including
  negative ones.

- `correction`: `"none"`, `"lingoes"`, or `"cailliez"`.

- `constant`: The constant that the correction added, \\c_1\\ or \\c_2\\
  below; `0` with `correction = "none"`, or when the distances needed no
  correction.

With a correction, `eig`, `prop`, and `trace` are those of the corrected
distances.

## Details

**Method**

With \\A\\ the matrix of squared distances and \\J = I -
\frac{1}{n}\mathbf{1}\mathbf{1}^T\\ the centering matrix, PCoA finds the
eigenvalues \\\lambda_i\\ and unit eigenvectors \\u_i\\ of \\B =
-\frac{1}{2}JAJ\\. Axis \\i\\'s coordinates are \\u_i\sqrt{\lambda_i}\\.

[`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html), `ape::pcoa()`,
and `vegan::wcmdscale()` build \\B\\ in full (\\n \times n\\) and
compute all \\n\\ of its eigenvalues, which takes time proportional to
\\n^3\\. `pcoa()` finds only the largest `n_axes` eigenvalues, by block
Lanczos iteration with thick restarts. That only needs products of \\B\\
with a few vectors at a time, which are computed in parallel straight
from `dm`, without ever forming \\B\\.

**Proportion of variance**

The total variance is \\trace(B) = \frac{1}{n}\sum\_{i\<j}d\_{ij}^2\\,
the sum of all \\n\\ eigenvalues. `prop` is each axis's eigenvalue
divided by it: the same as `ape::pcoa()`'s `Relative_eig` and `vegan`'s
"Proportion Explained", the values that `phyloseq::plot_ordination()`
shows on its axes.

Non-Euclidean distances, such as Bray-Curtis, can give negative
eigenvalues. They still count towards the total, so `prop` can sum to
more than 1 over the positive axes. Axes with eigenvalues that aren't
positive have no real coordinates and are dropped with a warning, as
[`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) does. So are axes
with eigenvalues at rounding-error level, below
`sqrt(.Machine$double.eps)` times the first one, as in
`vegan::wcmdscale()`.

**Negative eigenvalue corrections**

A correction adds a constant to the distances between different samples,
the smallest that makes them Euclidean, so that no eigenvalue is
negative (Gower & Legendre 1986; Legendre & Legendre 2012, sec. 9.3):

- `'lingoes'` (Lingoes 1971) adds \\2c_1\\ to the squared distances:
  \\d'\_{ij} = \sqrt{d\_{ij}^2 + 2c_1}\\, where \\c_1\\ is the magnitude
  of \\B\\'s most negative eigenvalue. That adds \\c_1\\ to every
  eigenvalue but the trivial 0 and leaves the eigenvectors unchanged, so
  the axes are those of `correction = 'none'`, with coordinates
  \\u_i\sqrt{\lambda_i + c_1}\\, and `trace` grows by \\(n-1)c_1\\. The
  axis whose eigenvalue was most negative gets eigenvalue 0.

- `'cailliez'` (Cailliez 1983) adds \\c_2\\ to the distances: \\d'\_{ij}
  = d\_{ij} + c_2\\, where \\c_2\\ is the largest real eigenvalue of the
  \\2n \times 2n\\ matrix \\\left\[\begin{smallmatrix}0 & 2B \\ -I &
  -4B_2\end{smallmatrix}\right\]\\, with \\B_2 = -\frac{1}{2}JDJ\\ the
  centered matrix of the unsquared distances. The axes change, and are
  computed from \\d'\\.

`eig`, `prop`, and `trace` are then those of the corrected distances,
and `prop` sums to at most 1. They equal `ape::pcoa()`'s `Corr_eig` and
`Rel_corr_eig` for the same correction, and
`vegan::wcmdscale(add = 'lingoes')` or `add = 'cailliez'`; \\c_2\\ is
`cmdscale(add = TRUE)`'s `ac`.

Negative eigenvalues above \\-\sqrt{\epsilon}\\ trace(B)\\ (with
\\\epsilon\\ = `.Machine$double.eps`, so about \\-1.5 \times 10^{-8}\\
of the total) count as rounding error: such distances, Euclidean ones
for instance, get no correction. `constant` is then 0, and the
coordinates, `eig`, `prop`, and `trace` are identical to those from
`correction = 'none'`. (`ape` and `vegan` correct those tiny negative
eigenvalues too, so they can add a tiny constant where `pcoa()` adds 0.)

With `n_axes` = \\n - 1\\, the axis whose eigenvalue the correction
makes 0 is dropped without a warning.

**When to use a correction**

Non-Euclidean distances, such as Bray-Curtis or Jaccard, have negative
eigenvalues, so the positive axes' `prop` can sum to more than 1. Most
plots simply show the first few positive axes, which is fine when the
negative eigenvalues are small next to the axes plotted. A correction
makes every `prop` a share of a total of 1, at a price: both corrections
inflate `trace`, and so shrink `prop`, more so for more samples. Which
inflates it more depends on the data. `'lingoes'` keeps the uncorrected
axes and only changes their eigenvalues; `'cailliez'` changes the axes
too. Another remedy is to transform the distances: the square root of
Bray-Curtis dissimilarities, for one, is Euclidean (Legendre & Legendre
2012).

**How the corrections are computed**

Neither the \\n \times n\\ nor the \\2n \times 2n\\ matrix is formed.
\\c_1\\ is the largest eigenvalue of \\-B\\, from the same Lanczos
solver. \\c_2\\ is the largest \\c\\ for which \\B + 2cB_2 +
\frac{c^2}{2}I\\ is singular (on vectors orthogonal to \\\mathbf{1}\\),
found by a Rayleigh-Ritz iteration that needs only products of \\B\\ and
\\B_2\\ with vectors, computed together from `dm`. That iteration can
settle on a root that isn't the largest (as with some exactly symmetric
designs), so \\c_2\\ is checked: the corrected matrix's smallest
eigenvalue, computed from a fresh start, must not be below
\\-\sqrt{\epsilon}\\ trace(B')\\. If it is, its eigenvector restarts the
iteration, for up to 10 rounds, after which there is a warning.

The smallest eigenvalues usually form a dense cluster, which Lanczos
resolves slowly, so corrections take longer than `correction = 'none'`.
For Bray-Curtis distances between 2,000 sparse samples, on one thread,
`'none'` took 0.4 seconds, `'lingoes'` 5 seconds, and `'cailliez'` 20
seconds; `cmdscale(add = TRUE)`, which builds the \\2n \times 2n\\
matrix, took 2 minutes. How much longer depends on how tightly the
smallest eigenvalues cluster. Each pass over `dm` takes time
proportional to \\n^2\\, while `cmdscale(add = TRUE)` and `ape::pcoa()`
take time proportional to \\n^3\\.

**Axis signs and reproducibility**

An eigenvector's sign is arbitrary, so different software can return
mirror images of the same axis. `pcoa()` flips each axis so that its
largest coordinate in absolute value is positive. Results are identical
from run to run and for any number of `cpus`.

`n_axes` larger than \\n - 1\\ (the most axes possible) is reduced to
\\n - 1\\.

## Ordination result

Every ordination function (`pcoa()`,
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capda()`](https://cmmr.github.io/ecodive/reference/capda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md))
returns the same kind of object: a numeric matrix with one row per
sample and one column per axis, so that it can go straight into a data
frame. Rows are named after the samples, when the input has names.
Columns are named after the method: `"PCoA1"`, `"PC1"`, `"RPCA1"`,
`"dbRDA1"` or `"CAP1"` (then `"MDS1"` for unconstrained axes), `"MDS1"`,
`"NMDS1"`, `"Isomap1"`, `"TSNE1"`, `"UMAP1"`, and so on.

Its attributes always begin with these three, in this order, each with
one value per column:

- `method`: The method's ID, e.g. `"pcoa"` or `"nmds"`.

- `eig`: Each axis's eigenvalue, on the method's own scale: `pcoa()` and
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) give
  the eigenvalues of Gower's centered matrix (sums of squares);
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) give
  variances, which are `pcoa()`'s divided by \\n - 1\\ for \\n\\
  samples; and
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) and
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  follow `vegan`, which also divides by \\n - 1\\ when any distance is 4
  or more. `NA` for
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), which
  have no eigenvalues.

- `prop`: The proportion of the total variance (or inertia) that each
  axis explains. Unlike `eig`, it is comparable between methods, so use
  it to label axes. `NA` where it isn't defined.

Attributes particular to the method follow them, and last, the records
of how it was made (see
[ecodive_steps](https://cmmr.github.io/ecodive/reference/ecodive_steps.md)):
`dist`, the distance matrix's, for the methods that take distances, and
`ord`, the ordination's own (its method, display name and settings). The
matrix has class `"ecodive_ord"` (as well as `"matrix"`), so that
[`print()`](https://rdrr.io/r/base/print.html) lists the other
attributes by name rather than in full, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
samples on the first two axes, labelled from the attributes, with any
loadings, constraints and centroids on top (see
[`plot.ecodive_ord()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md)).
`plot(ord, "Species", iris)` colours the samples by a column of
metadata. Subsetting with `[`, or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), gives a
plain matrix or data frame without the attributes.

## References

Cailliez, F. (1983). The analytical solution of the additive constant
problem. *Psychometrika*, 48(2), 305-308.
[doi:10.1007/BF02294026](https://doi.org/10.1007/BF02294026)

Gower, J. C. (1966). Some distance properties of latent root and vector
methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
[doi:10.1093/biomet/53.3-4.325](https://doi.org/10.1093/biomet/53.3-4.325)

Gower, J. C., & Legendre, P. (1986). Metric and Euclidean properties of
dissimilarity coefficients. *Journal of Classification*, 3(1), 5-48.
[doi:10.1007/BF01896809](https://doi.org/10.1007/BF01896809)

Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English
ed., sec. 9.3). Elsevier.

Lingoes, J. C. (1971). Some boundary conditions for a monotone analysis
of symmetric matrices. *Psychometrika*, 36(2), 195-203.
[doi:10.1007/BF02291398](https://doi.org/10.1007/BF02291398)

Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
Applications*, 22(2), 602-616.
[doi:10.1137/S0895479898334605](https://doi.org/10.1137/S0895479898334605)

## See also

[`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
[`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html)

## Examples

``` r
    dm  <- bray(iris[, 1:4])
    ord <- pcoa(dm)
    plot(ord, 'Species', iris)

    head(ord)
#>          PCoA1       PCoA2        PCoA3
#> [1,] 0.1692812 -0.02591990 -0.005053611
#> [2,] 0.1706354  0.02346821  0.030848760
#> [3,] 0.1870090 -0.00320247  0.032221661
#> [4,] 0.1801593  0.01132449  0.037211653
#> [5,] 0.1763516 -0.02998591 -0.006128668
#> [6,] 0.1363264 -0.03979963 -0.043074846
    str(ord)
#>  'ecodive_ord' num [1:150, 1:3] 0.169 0.171 0.187 0.18 0.176 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:3] "PCoA1" "PCoA2" "PCoA3"
#>  - attr(*, "method")= chr "pcoa"
#>  - attr(*, "eig")= num [1:3] 2.3473 0.2459 0.0766
#>  - attr(*, "prop")= num [1:3] 0.8901 0.0932 0.0291
#>  - attr(*, "trace")= num 2.64
#>  - attr(*, "correction")= chr "none"
#>  - attr(*, "constant")= num 0
#>  - attr(*, "dist")= Bray-Curtis Dissimilarity (metric = bray, norm = none, pseudocount = 0)
#>  - attr(*, "ord")= Principal Coordinates Analysis (method = pcoa, n_axes = 3, correction = none)
```
