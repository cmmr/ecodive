# Plot community types

Draws a result of
[`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)
or [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)
with base graphics.

## Usage

``` r
# S3 method for class 'ecodive_types'
plot(x, type = c("profiles", "selection", "scores"), n_features = 10, ...)
```

## Arguments

- x:

  A result of
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)
  or
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md).

- type:

  `"profiles"`, `"selection"` or `"scores"`.

- n_features:

  The number of features to show in `"profiles"`. Default `10`.

- ...:

  Passed on to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html),
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)
  or [`graphics::boxplot()`](https://rdrr.io/r/graphics/boxplot.html),
  such as `main`, `ylab` or `col`.

## Value

`x`, invisibly.

## Details

- `"profiles"` (the default): stacked bars of each type's composition
  (for
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
  its samples' mean proportions; for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
  the centroids), with the `n_features` features that reach the largest
  proportion in any type, and the rest pooled as "Other". The bars are
  labelled with the type and its number of samples.

- `"selection"`: the criterion that chose the number of types (for
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md),
  the Calinski-Harabasz index) against the number of types, with the one
  chosen circled. Not for
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md),
  whose centroids fix the types.

- `"scores"`
  ([`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)
  only): boxplots of the samples' scores (their similarity to their
  sub-CST's centroid), one per sub-CST with samples, in VALENCIA's
  colours, each beside a diamond and bar for the mean \\\pm\\ standard
  deviation of VALENCIA's reference samples, as VALENCIA's own plot
  shows them.

The types keep the order of `levels(x$type)` in every plot.

## Examples

``` r
    et <- enterotype(iris[, 1:4], n_types = 2:5)
    plot(et)

    plot(et, 'selection')


    # VALENCIA on thirteen samples made from its centroids.
    counts <- round(valencia_centroids * 10000)
    rownames(counts) <- paste0('S', 1:13)
    vl <- valencia(counts)
    plot(vl)

    plot(vl, 'scores')
```
