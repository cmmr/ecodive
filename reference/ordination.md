# Ordination Wrapper Function

Runs any of ecodive's ordination methods by name, as
[`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
does for beta diversity metrics.

## Usage

``` r
ordination(dm, method, ...)
```

## Arguments

- dm:

  The distances to ordinate: a `dist` object such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  or a symmetric matrix of distances. Every method listed by
  [`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md)
  takes distances.

- method:

  The name of an ordination method: `"capscale"`, `"dbrda"`, `"isomap"`,
  `"mds"`, `"nmds"`, `"pcoa"`, `"tsne"`, or `"umap"`. Flexible matching
  is supported (see below). Programmatic access via
  [`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md).

- ...:

  Further arguments for the method's function, such as `k`, `seed`,
  `cpus`, or `pcoa`'s `correction`. Each method keeps its own defaults;
  see
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  and [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)
  (which also need `groups`, or `data` and `formula`, as in
  `ordination(dm, 'dbrda', 'treatment', meta)`),
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
`ordination(dm, 'pcoa', ...)` is the same as `pcoa(dm, ...)`.

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

[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) ordinate a
count table rather than distances, so they are called directly, not
through `ordination()`.
[`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
and ecodive's other tests are run by
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md).

Default numbers of axes differ: `k = 3` for `capscale`, `dbrda`, `pcoa`
(and [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md)), and
`k = 2` for `isomap`, `mds`, `nmds`, `tsne`, and `umap`.

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
and so on, or use `ordination()`, which always runs ecodive's, to be
sure.

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
`ordination()`) returns the same kind of object: a numeric matrix with
one row per sample and one column per axis, so that it can go straight
into [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or a data
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
[`list_ordinations()`](https://cmmr.github.io/ecodive/reference/ordinations.md);
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), which are
called directly;
[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
for ecodive's tests.

## Examples

``` r
    # 30 samples in three groups.
    set.seed(1)
    counts <- matrix(rpois(30 * 20, 5), 30, 20)
    counts[11:20, 1:5] <- counts[11:20, 1:5] + 5
    counts[21:30, 6:10] <- counts[21:30, 6:10] + 5
    d <- bray(counts)

    ordination(d, 'pcoa')
#>               PCoA1        PCoA2         PCoA3
#>  [1,] -0.0360804114  0.091711440 -0.1105907710
#>  [2,] -0.0230948908  0.123523905 -0.0720542583
#>  [3,]  0.0523017762  0.023894547  0.0235011671
#>  [4,] -0.0009831621  0.069549555  0.0169203791
#>  [5,] -0.0019453437  0.054035906  0.0292109270
#>  [6,] -0.0371823183  0.093398842  0.1360915815
#>  [7,] -0.0495658911  0.022897755 -0.0384654522
#>  [8,]  0.0178610881  0.180607395 -0.0002660476
#>  [9,]  0.0484965214  0.067648469  0.0504743481
#> [10,]  0.0181776324  0.059231479 -0.0944284441
#> [11,] -0.1704793942 -0.017409877 -0.1143911894
#> [12,] -0.0815143578 -0.063620085  0.0198252902
#> [13,] -0.0969931110 -0.085572949  0.0926132764
#> [14,] -0.1065929852 -0.089396338 -0.0179647914
#> [15,] -0.1184732425 -0.116377213 -0.0303943691
#> [16,] -0.1665560361  0.027653615  0.0174461498
#> [17,] -0.2065533700  0.078005363  0.0479182367
#> [18,] -0.1541357844 -0.056715440  0.0690006947
#> [19,] -0.1668647980 -0.057014693 -0.0136648054
#> [20,] -0.1185943205 -0.071167788 -0.0275055595
#> [21,]  0.0607370055 -0.096812933  0.0213953915
#> [22,]  0.1124361421  0.009240374  0.0398499073
#> [23,]  0.1513725637 -0.010969491  0.0091507334
#> [24,]  0.1776768087 -0.036775431  0.0316047331
#> [25,]  0.1493110539 -0.066778910 -0.1082642886
#> [26,]  0.2149586662 -0.049343885 -0.0087576882
#> [27,]  0.1359001225  0.012315397  0.0378789660
#> [28,]  0.1551753035  0.002891378 -0.0037064635
#> [29,]  0.1120427050 -0.034752576 -0.0060420219
#> [30,]  0.1291620279 -0.063897809  0.0036143680
#> 
#> pcoa: 30 samples, 3 axes; % of variance: 39.1, 14.3, 9.0
#> Attributes: method, eig, prop, trace, correction, constant

    ordination(d, 'NMDS', k = 2)
#>              NMDS1        NMDS2
#>  [1,] -0.076781363  0.184513370
#>  [2,] -0.040414313  0.248387588
#>  [3,]  0.064948612 -0.002575346
#>  [4,] -0.003660619  0.096767066
#>  [5,]  0.003267326  0.043200706
#>  [6,] -0.091256926  0.118311810
#>  [7,] -0.068782527  0.016120041
#>  [8,]  0.093922221  0.254939142
#>  [9,]  0.073463986  0.040734910
#> [10,]  0.088556261  0.142786881
#> [11,] -0.240601041 -0.003062966
#> [12,] -0.101167871 -0.090308848
#> [13,] -0.119556208 -0.184893547
#> [14,] -0.142478350 -0.117483748
#> [15,] -0.155876934 -0.106334842
#> [16,] -0.188632895  0.019645493
#> [17,] -0.255225924  0.069545906
#> [18,] -0.203245632 -0.086022836
#> [19,] -0.266318367 -0.061424564
#> [20,] -0.153158074 -0.061986908
#> [21,]  0.063659131 -0.086858301
#> [22,]  0.163322207  0.012150440
#> [23,]  0.185608829 -0.043395410
#> [24,]  0.231877693 -0.045640625
#> [25,]  0.192147954 -0.100967879
#> [26,]  0.272631460 -0.090684559
#> [27,]  0.172236247 -0.009536038
#> [28,]  0.210499651  0.013269110
#> [29,]  0.132892496 -0.065269638
#> [30,]  0.158122970 -0.103926413
#> 
#> nmds: 30 samples, 2 axes
#> Attributes: method, eig, prop, stress, converged, n_iters, n_starts, best_start, n_best, ties

    # The same arguments for every method.
    ords <- lapply(c('pcoa', 'nmds', 'umap'), function (m)
      ordination(d, m, k = 2, seed = 1, cpus = 1) )
```
