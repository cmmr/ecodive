# Count Table Ordination

Runs any of ecodive's ordination methods that take a count table by
name, as
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
does for the methods that take distances.

## Usage

``` r
count_ord(counts, method, ...)
```

## Arguments

- counts:

  A numeric matrix of count data (samples \\\times\\ features).
  Typically contains absolute abundances (integer counts), though
  proportions are also accepted.

- method:

  The name of an ordination method: `"pca"`, `"rda"`, or `"rpca"`.
  Flexible matching is supported (see below). Programmatic access via
  [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md).

- ...:

  Further arguments for the method's function, such as `n_axes`,
  `margin`, `cpus`,
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s `norm`,
  `pseudocount`, and `scale`, or
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md)'s `groups`,
  `data`, and `formula`. Each method keeps its own defaults; see
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md), and
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md). `seed`
  is ignored by the methods that have no random parts (all of them, at
  present), so the same arguments can be given to every method.

## Value

The method's result, as described in "Ordination result" below.
`count_ord(counts, 'pca', ...)` is the same as `pca(counts, ...)`.

## Details

**List of Ordination Methods**

|                        |                                        |
|------------------------|----------------------------------------|
| Option / Function Name | Method Name                            |
| `pca`                  | Principal Component Analysis (PCA)     |
| `rda`                  | Redundancy Analysis (RDA), constrained |
| `rpca`                 | Robust Aitchison PCA                   |

`rda` is constrained by sample metadata, given as further arguments:
`count_ord(counts, 'rda', 'treatment', meta)` is
`rda(counts, 'treatment', meta)`. Without them, it is a PCA.

Methods that ordinate distances, such as
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) and
[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), are run by
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md).

**Flexible name matching**

Case insensitive and partial matching, ignoring anything but letters.
E.g. `method = 'PCA'` or `'principal comp'` selects `pca`. Common other
names also work: `'prcomp'` for `pca`; `'redundancy analysis'`,
`'tbRDA'`, or `'transformation-based RDA'` for `rda`; and
`'robust PCA'`, `'DEICODE'`, or `'gemelli'` for `rpca`.

Finished code should always use the full primary option name to avoid
ambiguity with future additions to the list of methods.

**Same names in other packages**

`vegan::pca()` and `vegan::rda()` share their names with ecodive's
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md), and
whichever package was attached last masks the other's. Write
[`ecodive::pca()`](https://cmmr.github.io/ecodive/reference/pca.md) or
[`ecodive::rda()`](https://cmmr.github.io/ecodive/reference/rda.md), or
use `count_ord()`, which always runs ecodive's, to be sure.

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
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md))
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

[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md);
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) for
the methods that take distances.

## Examples

``` r
    ord <- count_ord(ex_counts, 'pca', norm = 'clr', pseudocount = 1)
    head(ord)
#>              PC1       PC2          PC3
#> Saliva -3.528246 -1.336813  0.879437784
#> Gums   -3.691532 -1.381718 -0.860746227
#> Nose    6.050248 -2.662859 -0.010805525
#> Stool   1.169530  5.381390 -0.007886032
    str(ord)
#>  'ecodive_ord' num [1:4, 1:3] -3.53 -3.69 6.05 1.17 -1.34 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:4] "Saliva" "Gums" "Nose" "Stool"
#>   ..$ : chr [1:3] "PC1" "PC2" "PC3"
#>  - attr(*, "method")= chr "pca"
#>  - attr(*, "eig")= num [1:3] 21.35 13.249 0.505
#>  - attr(*, "prop")= num [1:3] 0.6082 0.3774 0.0144
#>  - attr(*, "trace")= num 35.1
#>  - attr(*, "loadings")= num [1:6, 1:3] -0.4926 -0.0683 0.4635 -0.564 0.3082 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:6] "Streptococcus" "Bacteroides" "Corynebacterium" "Haemophilus" ...
#>   .. ..$ : chr [1:3] "PC1" "PC2" "PC3"

    # The same arguments for every method.
    lapply(list_count_ords('id'), function (m) attr(count_ord(ex_counts, m, n_axes = 2), 'prop'))
#> [[1]]
#> [1] 0.5428410 0.3713532
#> 
#> [[2]]
#> [1] 0.5428410 0.3713532
#> 
#> [[3]]
#> [1] 0.5737685 0.4262315
#> 
```
