# Principal Component Analysis (PCA)

Ordinates samples from their feature abundances, optionally after a
normalization such as the centered log ratio (CLR, giving an "Aitchison
PCA"). Gives the same scores, loadings, and variances as
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) on the
normalized table, but computes only the first few axes, straight from
sparse input, which takes far less time and memory for large datasets.

## Usage

``` r
pca(
  counts,
  margin = 1L,
  norm = "none",
  pseudocount = NULL,
  n_axes = 3L,
  scale = FALSE,
  cpus = n_cpus()
)
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

- norm:

  Normalize the incoming counts before the PCA. Options are:

  - `'none'`: No transformation.

  - `'percent'`: Relative abundance (sample abundances sum to 1).

  - `'chord'`: Each sample scaled to unit length (Euclidean norm 1).

  - `'hellinger'`: Square root of relative abundance (the Hellinger
    transformation).

  - `'binary'`: Presence/absence (each count is either 0 or 1).

  - `'clr'`: Centered log ratio, giving an Aitchison PCA.

  - `'rclr'`: Robust centered log ratio. Zeros stay 0.

  - `'log'`: Log relative abundance (log-TSS), with `pseudocount` added
    to every value first.

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'` or `'log'`.
  Ignored for other normalization methods. See **Pseudocount** section.

- n_axes:

  The number of axes to compute. Default: `3`

- scale:

  Scale each feature to unit variance (after `norm`, before the PCA), as
  in `prcomp(scale. = TRUE)`. This gives every feature equal weight, so
  that a PCA of the correlation matrix rather than the covariance matrix
  results. Default: `FALSE`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), is
  the CPU cores this R process may use, within any limits set on it.
  Results are the same for any `cpus`.

## Value

A numeric matrix with one row per sample and one column per axis, named
`"PC1"`, `"PC2"`, etc., as described in "Ordination result" below. Its
attributes are:

- `method`: `"pca"`.

- `eig`: The variance of each axis:
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)'s `sdev^2`.

- `prop`: The proportion of the total variance that each axis explains,
  `eig / trace`.

- `trace`: The total variance: the sum of all features' variances (after
  `norm` and `scale`), which is also the sum of all axes' `eig`.

- `loadings`: A matrix with one row per feature and one column per axis:
  the unit-length principal axes,
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)'s `rotation`.

- `norm`, `scale`: The transform the loadings are in (`norm`'s name and
  `scale`), which qualifies them:
  [plot()](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md)
  reads them to say whether the links of a compositional biplot are log
  ratios.

## Details

**Method**

With \\Z\\ the normalized table (samples in rows), \\J\\ the centering
matrix, and \\D\\ the diagonal matrix of the features' standard
deviations (or the identity matrix when `scale = FALSE`), PCA is the
singular value decomposition \\JZD^{-1} = USV^T\\. The scores are
\\US\\, the loadings \\V\\, and each axis's variance (`eig`) is \\s^2 /
(n - 1)\\, for \\n\\ samples.

[`prcomp()`](https://rdrr.io/r/stats/prcomp.html) builds the whole
centered table and decomposes it in full. `pca()` finds only the first
`n_axes` axes, by block Lanczos iteration, which only needs products of
the table with a few vectors at a time. Those are computed in parallel
from the non-zero counts alone: the zeros of each sample all share one
value after normalization (0, or \\\log(pseudocount)\\ minus the
sample's mean log for `norm = 'clr'`), and centering only subtracts each
feature's mean, so a sparse table is never made dense.

**Relation to
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)**

The scores are the same as those of
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) on the
Euclidean distances between the normalized samples: `pca(counts)` and
`pcoa(euclidean(counts))`, or `pca(counts, norm = 'clr')` and
`pcoa(aitchison(counts))`, give the same coordinates. PCoA's eigenvalues
are \\s^2\\, so `pca()`'s `eig` is
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)'s divided
by \\n - 1\\; `prop` is the same for both. What `pca()` adds is the
`loadings`: how much each feature contributes to each axis, e.g. for a
biplot.

**Biplots**

`plot(ord)` draws a form biplot: the samples at their scores, and the
loadings as arrows scaled to fit. `plot(ord, scaling = 'covariance')` of
a `norm = 'clr'` PCA is Aitchison and Greenacre's (2002) compositional
biplot, in which the rays' lengths and the `links` between their tips
read as standard deviations of clr values and log ratios;
`scaling = 'contribution'` marks the features that contribute most to
the plane. See "Compositional biplots" in
[`plot.ecodive_ord()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_ord.md).

**Robust Aitchison PCA**

`norm = 'rclr'` with `pca()` is a plain PCA of the rclr-transformed
table, with its zeros left as 0. DEICODE and gemelli's Robust Aitchison
PCA instead treat those zeros as missing values and fill them in by
low-rank matrix completion (OptSpace) before the decomposition, so their
results differ.

**Axis signs and reproducibility**

A principal axis's sign is arbitrary, so different software can return
mirror images of the same axis. As
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) does,
`pca()` flips each axis so that its largest score in absolute value is
positive, and flips its loadings to match. Results are identical from
run to run and for any number of `cpus`.

`n_axes` larger than the number of samples minus one, or the number of
features, is reduced to the smaller of those. Axes with variances at
rounding-error level, below `sqrt(.Machine$double.eps)` times the first
one, are dropped with a warning. With `scale = TRUE`, features with
constant values (such as features absent from every sample) are left
out, where [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) would stop
with an error.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md), `pca()`,
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

- `eig`: Each axis's eigenvalue, on the method's own scale:
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) and
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) give
  the eigenvalues of Gower's centered matrix (sums of squares); `pca()`
  and [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) give
  variances, which are
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)'s divided
  by \\n - 1\\ for \\n\\ samples; and
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

Attributes particular to the method follow them, and last `display`,
with the ordination's title for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`print()`](https://rdrr.io/r/base/print.html): the method's name, and
for the methods that take distances, the distance matrix's title, such
as "Principal Coordinates Analysis of Bray-Curtis Dissimilarity" (see
[ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)).
The matrix has class `"ecodive_ord"` (as well as `"matrix"`), so that
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

## References

Pearson, K. (1901). On lines and planes of closest fit to systems of
points in space. *The London, Edinburgh, and Dublin Philosophical
Magazine and Journal of Science*, 2(11), 559-572.
[doi:10.1080/14786440109462720](https://doi.org/10.1080/14786440109462720)

Aitchison, J. (1983). Principal component analysis of compositional
data. *Biometrika*, 70(1), 57-65.
[doi:10.1093/biomet/70.1.57](https://doi.org/10.1093/biomet/70.1.57)

## See also

[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)

## Examples

``` r
    ord <- pca(iris[, 1:4])
    plot(ord, 'Species', iris)

    ord
#> # A matrix: 150 samples × 3 axes
#> # Principal Component Analysis
#> # % of variance: 92.5, 5.3, 1.7
#>      PC1     PC2     PC3
#>  1 -2.68  0.319  -0.0279
#>  2 -2.71 -0.177  -0.21
#>  3 -2.89 -0.145   0.0179
#>  4 -2.75 -0.318   0.0316
#>  5 -2.73  0.327   0.0901
#>  6 -2.28  0.741   0.169
#>  7 -2.82 -0.0895  0.258
#>  8 -2.63  0.163  -0.0219
#>  9 -2.89 -0.578   0.0208
#> 10 -2.67 -0.114  -0.198
#> # … 140 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method, eig, prop, trace, loadings, norm, scale
    str(ord)
#>  'ecodive_ord' num [1:150, 1:3] -2.68 -2.71 -2.89 -2.75 -2.73 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:3] "PC1" "PC2" "PC3"
#>  - attr(*, "method")= chr "pca"
#>  - attr(*, "eig")= num [1:3] 4.2282 0.2427 0.0782
#>  - attr(*, "prop")= num [1:3] 0.9246 0.0531 0.0171
#>  - attr(*, "trace")= num 4.57
#>  - attr(*, "loadings")= num [1:4, 1:3] 0.3614 -0.0845 0.8567 0.3583 0.6566 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
#>   .. ..$ : chr [1:3] "PC1" "PC2" "PC3"
#>  - attr(*, "norm")= chr "none"
#>  - attr(*, "scale")= logi FALSE
#>  - attr(*, "display")=List of 1
#>   ..$ title: chr "Principal Component Analysis"
```
