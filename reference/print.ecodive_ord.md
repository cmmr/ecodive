# Print an ordination result

Prints the coordinates, then one line on the axes and the names of the
other attributes, rather than every attribute in full (such as
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s loadings).

## Usage

``` r
# S3 method for class 'ecodive_ord'
print(x, ...)
```

## Arguments

- x:

  An ordination result, from
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md),
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md),
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md),
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), or
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md).

- ...:

  Passed on to [`print()`](https://rdrr.io/r/base/print.html) for the
  matrix of coordinates.

## Value

`x`, invisibly.
