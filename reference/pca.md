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
  k = 3L,
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

  - `'binary'`: Presence/absence (each count is either 0 or 1).

  - `'clr'`: Centered log ratio, giving an Aitchison PCA.

  - `'rclr'`: Robust centered log ratio. Zeros stay 0.

  Default: `'none'`.

- pseudocount:

  Value added to counts to handle zeros when `norm = 'clr'`. Ignored for
  other normalization methods. See **Pseudocount** section.

- k:

  The number of axes to compute. Default: `3`

- scale:

  Scale each feature to unit variance (after `norm`, before the PCA), as
  in `prcomp(scale. = TRUE)`. This gives every feature equal weight, so
  that a PCA of the correlation matrix rather than the covariance matrix
  results. Default: `FALSE`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

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
`k` axes, by block Lanczos iteration, which only needs products of the
table with a few vectors at a time. Those are computed in parallel from
the non-zero counts alone: the zeros of each sample all share one value
after normalization (0, or \\\log(pseudocount)\\ minus the sample's mean
log for `norm = 'clr'`), and centering only subtracts each feature's
mean, so a sparse table is never made dense.

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

`k` larger than the number of samples minus one, or the number of
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
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
[`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md))
returns the same kind of object: a numeric matrix with one row per
sample and one column per axis, so that it can go straight into
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
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

Attributes particular to the method follow them. The matrix has class
`"ecodive_ord"` (as well as `"matrix"`), only so that
[`print()`](https://rdrr.io/r/base/print.html) lists the other
attributes by name rather than in full. Subsetting with `[`, or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), gives a
plain matrix or data frame without them.

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

The `pseudocount` parameter is only relevant when `norm = 'clr'`.

Zeros are undefined in the centered log-ratio (CLR) transformation. If
`norm = 'clr'`, `pseudocount` is `NULL` (the default), and zeros are
detected, the function uses half the minimum non-zero value
(`min(x[x>0]) / 2`) and issues a warning.

To suppress the warning, provide an explicit value (e.g., `1`).

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
[`ordination()`](https://cmmr.github.io/ecodive/reference/ordination.md),
[`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)

## Examples

``` r
    ord <- pca(ex_counts, norm = 'clr', pseudocount = 1)
    ord
#>              PC1       PC2          PC3
#> Saliva -3.528246 -1.336813  0.879437784
#> Gums   -3.691532 -1.381718 -0.860746227
#> Nose    6.050248 -2.662859 -0.010805525
#> Stool   1.169530  5.381390 -0.007886032
#> 
#> pca: 4 samples, 3 axes; % of variance: 60.8, 37.7, 1.4
#> Attributes: method, eig, prop, trace, loadings

    # Percent of variance explained by each axis.
    round(100 * attr(ord, 'prop'), 1)
#> [1] 60.8 37.7  1.4

    # Features contributing most to the first axis.
    head(sort(abs(attr(ord, 'loadings')[, 1]), decreasing = TRUE))
#>       Haemophilus     Streptococcus   Corynebacterium    Staphylococcus 
#>        0.56399158        0.49259684        0.46351081        0.35322904 
#> Propionibacterium       Bacteroides 
#>        0.30818630        0.06833773 

    # Same scores as prcomp(), up to the sign of each axis.
    clr <- log(ex_counts + 1)
    clr <- clr - rowMeans(clr)
    head(prcomp(clr)$x[, 1:3])
#>              PC1       PC2          PC3
#> Saliva  3.528246  1.336813  0.879437784
#> Gums    3.691532  1.381718 -0.860746227
#> Nose   -6.050248  2.662859 -0.010805525
#> Stool  -1.169530 -5.381390 -0.007886032

    plot(ord[, 1:2], asp = 1)

```
