# Uniform Manifold Approximation and Projection (UMAP)

Ordinates samples from their distances, for plotting beta diversity, by
placing samples that are near neighbours in `dm` near each other on the
plot. Follows the same procedure as `uwot::umap(batch = TRUE)` given a
`dist` object, and gives the same layout from the same settings,
starting coordinates and random seeds, but runs in parallel C, without
ever copying `dm` into an \\n \times n\\ matrix.

## Usage

``` r
umap(
  dm,
  n_axes = 2L,
  n_neighbors = 15L,
  min_dist = 0.1,
  spread = 1,
  n_epochs = NULL,
  learning_rate = 1,
  negative_sample_rate = 5,
  repulsion_strength = 1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  init = "spectral",
  seed = 0L,
  cpus = n_cpus()
)
```

## Arguments

- dm:

  A `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  or [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), or a
  symmetric numeric matrix of distances. Must not contain `NA`s.

- n_axes:

  The number of dimensions. Default: `2`

- n_neighbors:

  How many nearest neighbours of each sample (counting the sample
  itself) to use. Larger values preserve more of the global structure;
  smaller ones more of the local detail. Values larger than the number
  of samples are reduced to it. Default: `15`

- min_dist:

  How close samples can be packed together on the plot, relative to
  `spread`. Default: `0.1`, as in umap-learn (uwot's is `0.01`).

- spread:

  The scale of the plot's distances, together with `min_dist`. Default:
  `1`

- n_epochs:

  The number of optimization epochs. `NULL`, the default, uses 500 for
  up to 10,000 samples and 200 for more, as `uwot::umap()` does.

- learning_rate:

  The initial learning rate (step size) of the optimization. It
  decreases linearly to 0 over the epochs. Default: `1`

- negative_sample_rate:

  How many non-neighbours (negative samples) to push away for each
  neighbour pulled closer. Default: `5`

- repulsion_strength:

  The weight of the repulsion from negative samples, relative to the
  attraction of neighbours. Default: `1`

- set_op_mix_ratio:

  How the neighbour graphs of the samples are combined, from `0` (fuzzy
  intersection) to `1` (fuzzy union). Default: `1`

- local_connectivity:

  How many nearest neighbours of each sample are assumed to be fully
  connected to it. Default: `1`

- bandwidth:

  Scales the effective number of neighbours, as in `uwot::umap()`.
  Default: `1`

- init:

  The starting coordinates: `"spectral"` (the default), `"pcoa"`,
  `"random"`, or a numeric matrix with one row per sample and `n_axes`
  columns. See "Starting coordinates" below.

- seed:

  An integer seed for the random parts of the method: the negative
  samples, random starting coordinates, and the small random jitter
  added to spectral ones. The same seed gives the same result. Default:
  `0`

- cpus:

  How many parallel processing threads should be used. The default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), will
  use all logical CPU cores.

## Value

A numeric matrix with one row per sample and one column per dimension,
named `"UMAP1"`, `"UMAP2"`, etc., as described in "Ordination result"
below. `method` is `"umap"`; `eig` and `prop` are `NA`, as UMAP has no
eigenvalues. Then come:

- `a`, `b`: The parameters of the curve \\1 / (1 + a d^{2b})\\ that
  turns plot distances \\d\\ into similarities, fitted from `min_dist`
  and `spread`.

- `init`: How the starting coordinates were chosen: `"spectral"`,
  `"pcoa"`, `"random"`, or `"matrix"` (given in `init`); `"none"` for
  fewer than 2 samples, which are placed at 0.

- `n_epochs`: The number of optimization epochs run.

- `n_neighbors`: The `n_neighbors` given.

## Details

**Method**

UMAP (McInnes et al. 2018) builds a weighted graph of each sample's
nearest neighbours, then finds plot coordinates whose own neighbour
graph is as similar to it as possible.

1.  **Neighbours.** Each sample's `n_neighbors` nearest neighbours,
    counting itself, are found exactly from `dm`. Ties go to the sample
    that comes first.

2.  **Graph.** A neighbour \\j\\ of sample \\i\\ at distance \\d\_{ij}\\
    gets weight \\\exp(-(d\_{ij} - \rho_i) / \sigma_i)\\, where
    \\\rho_i\\ is the distance to the nearest neighbour (other than
    duplicates) and \\\sigma_i\\ is chosen by binary search so that the
    weights sum to \\\log_2(\\`n_neighbors`\\)\\. The two weights
    \\w\_{ij}\\ and \\w\_{ji}\\ of each pair are combined by fuzzy
    union, \\w\_{ij} + w\_{ji} - w\_{ij} w\_{ji}\\.

3.  **Curve.** `a` and `b` are fitted by least squares
    (Levenberg-Marquardt, as umap-learn fits them) to a curve that is 1
    up to `min_dist` and decays exponentially, with scale `spread`,
    beyond it.

4.  **Layout.** Starting from `init`, stochastic gradient descent pulls
    the two samples of each edge together, sampling each edge in
    proportion to its weight, and pushes each sample away from
    `negative_sample_rate` randomly chosen samples per edge sampled. As
    in `uwot::umap()`'s batch mode, each epoch computes every sample's
    gradient from the coordinates at the start of the epoch, then moves
    all samples at once by the Adam optimizer (\\\beta_1 = 0.5\\,
    \\\beta_2 = 0.9\\, \\\epsilon = 10^{-7}\\, uwot's settings for batch
    mode). Finally, the coordinates are centered.

**Why batch mode**

`uwot::umap()`'s default optimization moves samples one edge at a time,
each move seeing the ones before it. That can only run on one thread if
it is to be reproducible, and it takes nearly all of the run time. In
batch mode, the samples' gradients within an epoch don't depend on each
other, so `umap()` computes them on `cpus` threads, and still gets the
same result for any number of threads. The layouts differ from those of
uwot's default mode only as different random draws would: neighbour
preservation, trustworthiness, and the final cross-entropy are about the
same.

**Starting coordinates**

- `"spectral"`: The eigenvectors of the graph's normalized Laplacian for
  its 2nd to \\(k+1)\\th smallest eigenvalues, scaled so that the
  largest absolute coordinate is 10, plus normal noise with standard
  deviation 0.0001. When the graph is not connected, the `"pcoa"` start
  is used instead, scaled to standard deviation 1, and with fewer than 3
  samples (or \\k \> n - 2\\), the `"random"` one, as in `uwot::umap()`.

- `"pcoa"`: The first `n_axes` axes of
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md), each
  centered and scaled to standard deviation 0.0001 (`uwot::umap()`'s
  `"spca"`).

- `"random"`: Coordinates uniform on \[-10, 10).

- A matrix: Used as is.

**Differences from `uwot::umap()`**

- `umap()` always optimizes in batch mode, so its layouts are those of
  `uwot::umap(batch = TRUE)`, not of uwot's default.

- `min_dist` defaults to `0.1`, the canonical value used by umap-learn
  and most other implementations, rather than uwot's `0.01`. Other
  defaults are the same.

- The random numbers come from `seed`, not from R's random number
  generator, so [`set.seed()`](https://rdrr.io/r/base/Random.html) has
  no effect, and R's random number stream is left untouched. In batch
  mode, `uwot::umap()` draws from R's generator two 32-bit seeds per
  sample per epoch, for each sample's stream of negative samples, and
  the jitter of spectral starts. `umap()` makes those from `seed`
  instead, with the same streams of negative samples (PCG) from the same
  seeds; given the same seeds and starting coordinates, the two give the
  same coordinates (to the last bit, for the same compiler and math
  library).

- The spectral start's eigenvectors are computed by block Lanczos
  iteration, the same solver as
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)'s, to
  relative residuals of \\10^{-14}\\ (at worst \\10^{-10}\\) rather than
  `RSpectra`'s tolerance of \\10^{-4}\\, and each one's sign is fixed
  (its largest absolute coordinate is positive).

- The `"pcoa"` start, and the fallback for unconnected graphs, come from
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) rather
  than [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html); they
  differ only by rounding and the signs of the axes.

- Finally, each axis is flipped so that its largest absolute coordinate
  is positive.

**Speed and memory**

`uwot::umap()` given a `dist` object copies it into an \\n \times n\\
matrix and sorts each of its columns in R. `umap()` reads the packed
triangle of `dm` directly, keeping only the `n_neighbors` closest of
each sample. Beyond `dm` itself, memory use is proportional to \\n\\
times `n_neighbors`. `cpus` threads find the neighbours, build the
graph, compute the spectral start, and compute each epoch's gradients.

**Reproducibility**

Results are identical from run to run and for any number of `cpus`.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capda()`](https://cmmr.github.io/ecodive/reference/capda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), `umap()`,
and so
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
  the eigenvalues of Gower's centered matrix (sums of squares);
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) give
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
  `umap()`, which have no eigenvalues.

- `prop`: The proportion of the total variance (or inertia) that each
  axis explains. Unlike `eig`, it is comparable between methods, so use
  it to label axes. `NA` where it isn't defined.

Attributes particular to the method follow them. The matrix has class
`"ecodive_ord"` (as well as `"matrix"`), so that
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

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold
approximation and projection for dimension reduction. *arXiv*,
1802.03426.
[doi:10.48550/arXiv.1802.03426](https://doi.org/10.48550/arXiv.1802.03426)

Melville, J. (2020). uwot: The Uniform Manifold Approximation and
Projection (UMAP) Method for Dimensionality Reduction. R package. Its
C++ header library (BSD 2-Clause License) is the basis of the graph and
layout code. <https://github.com/jlmelville/uwot>

Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
Applications*, 22(2), 602-616.
[doi:10.1137/S0895479898334605](https://doi.org/10.1137/S0895479898334605)

## See also

[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)

## Examples

``` r
    dm  <- bray(iris[, 1:4])
    ord <- umap(dm)
    plot(ord, 'Species', iris)

    head(ord)
#>         UMAP1       UMAP2
#> [1,] 5.254048 -1.77263912
#> [2,] 4.404295 -0.18979135
#> [3,] 5.162174  0.11891402
#> [4,] 4.950985 -0.01741511
#> [5,] 5.250620 -1.60047599
#> [6,] 4.974723 -2.73126598
    str(ord)
#>  'ecodive_ord' num [1:150, 1:2] 5.25 4.4 5.16 4.95 5.25 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "UMAP1" "UMAP2"
#>  - attr(*, "method")= chr "umap"
#>  - attr(*, "eig")= num [1:2] NA NA
#>  - attr(*, "prop")= num [1:2] NA NA
#>  - attr(*, "a")= num 1.58
#>  - attr(*, "b")= num 0.895
#>  - attr(*, "init")= chr "pcoa"
#>  - attr(*, "n_epochs")= int 500
#>  - attr(*, "n_neighbors")= int 15
```
