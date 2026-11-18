# Distance-based Ordination

Runs any of ecodive's distance-based ordination methods by name, as
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
does for beta diversity metrics and
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
does for statistical tests. Compute the distances once with
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md) (or
a metric's own function, such as
[`bray()`](https://cmmr.github.io/ecodive/reference/bray.md)), then pass
them to `dist_ord()` and
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md).

## Usage

``` r
dist_ord(dm, method, ...)
```

## Arguments

- dm:

  The distances to ordinate: a `dist` object such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  or a symmetric matrix of distances. Every method listed by
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  takes distances.

- method:

  The name of an ordination method: `"capscale"`, `"dbrda"`, `"isomap"`,
  `"mds"`, `"nmds"`, `"pcoa"`, `"tsne"`, or `"umap"`. Flexible matching
  is supported (see below). Programmatic access via
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md).

- ...:

  Further arguments for the method's function, such as `n_axes`, `seed`,
  `cpus`, or `pcoa`'s `correction`. Each method keeps its own defaults;
  see
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  and [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)
  (which also need `groups`, or `data` and `formula`, as in
  `dist_ord(dm, 'dbrda', 'treatment', meta)`),
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md). `seed`
  is ignored by the methods that have no random parts (`capscale`,
  `dbrda`, `isomap` and `pcoa`), so the same arguments can be given to
  every method.

## Value

The method's result, as described in "Ordination result" below.
`dist_ord(dm, 'pcoa', ...)` is the same as `pcoa(dm, ...)`.

## Details

**List of Ordination Methods**

|  |  |
|----|----|
| Option / Function Name | Method Name |
| `capscale` | Constrained Analysis of Principal Coordinates |
| `dbrda` | Distance-based Redundancy Analysis (db-RDA) |
| `isomap` | Isometric Feature Mapping (Isomap) |
| `mds` | Metric Multidimensional Scaling (MDS) |
| `nmds` | Non-metric Multidimensional Scaling (NMDS) |
| `pcoa` | Principal Coordinates Analysis (PCoA) |
| `tsne` | t-distributed Stochastic Neighbor Embedding (t-SNE) |
| `umap` | Uniform Manifold Approximation and Projection (UMAP) |

[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) ordinate a
count table rather than distances, so they are run by
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md),
or called directly. PCA is PCoA of Euclidean distances:
`dist_ord(euclidean(counts), 'pcoa')` gives the same scores as
`pca(counts)`, up to each axis's sign, but
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) also returns
the features' loadings and never forms the distance matrix.
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
and ecodive's other tests are run by
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md).

Default numbers of axes differ: `n_axes = 3` for the methods whose axes
have eigenvalues (`capscale`, `dbrda`, `isomap`, `pcoa`, and
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md)), and
`n_axes = 2` for those that fit a layout in a fixed number of dimensions
(`mds`, `nmds`, `tsne`, and `umap`).

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `method = 'PCoA'` or `'principal coord'` selects `pcoa`. Common
other names also work: `'PCO'`, `'cmdscale'`, or `'classical MDS'` for
`pcoa`; `'MDS'` or `'metric MDS'` for `mds`; `'non-metric MDS'`,
`'metaMDS'`, `'monoMDS'`, or `'isoMDS'` for `nmds`;
`'isometric mapping'` for `isomap`; `'Rtsne'`, `'bhtsne'`, or
`'Barnes-Hut t-SNE'` for `tsne`; and `'uwot'` for `umap`. `'dbRDA'` or
`'distance-based RDA'` select `dbrda`, and `'CAP'` selects `capscale`.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

**Same names in other packages**

`ape::pcoa()`, `vegan::pca()`, `vegan::dbrda()`, `vegan::capscale()`,
`vegan::isomap()` (whose `k` is the number of neighbours, not axes), and
`uwot::umap()` share names with ecodive's functions, and whichever
package was attached last masks the other's. Write
[`ecodive::pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)
and so on, or use `dist_ord()`, which always runs ecodive's, to be sure.

## Ordination result

Every ordination function
([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), and so
`dist_ord()`) returns the same kind of object: a numeric matrix with one
row per sample and one column per axis, so that it can go straight into
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

## See also

[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
[`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md),
[`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md);
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
for [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), which take
a count table;
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
for ecodive's tests.

## Examples

``` r
    dm  <- bray(iris[, 1:4])
    ord <- dist_ord(dm, 'pcoa')
    plot(ord, col = iris$Species)

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
```
