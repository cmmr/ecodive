# Plot an ordination

Draws the samples on two axes of an ordination from
[`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
[`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md)
or any other of ecodive's ordination methods, with the title, axis
labels and subtitle taken from the result's attributes. The overlays the
result holds are drawn on top of the samples: feature loadings, the
arrows of numeric constraints, factor or group centroids, and the target
of a Procrustes rotation. Nothing is recomputed: this is a quick look at
results already in hand.

## Usage

``` r
# S3 method for class 'ecodive_ord'
plot(
  x,
  groups = NULL,
  data = NULL,
  axes = c(1, 2),
  loadings = 10,
  biplot = TRUE,
  centroids = TRUE,
  target = TRUE,
  legend = "topright",
  col = NULL,
  pch = NULL,
  ...
)
```

## Arguments

- x:

  An ordination result: a matrix of class `ecodive_ord`, from
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md),
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md),
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md),
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md),
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md),
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md), or
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md).

- groups:

  Colours the samples by group, with a legend: the name of a column of
  `data`, or a vector with one value per sample (matched to the samples
  by name when both are named, as for
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)).
  A factor, character or logical vector gives each value its own colour,
  in the order of a factor's levels; a numeric vector gives a colour
  gradient (convert it with
  [`factor()`](https://rdrr.io/r/base/factor.html) for one colour per
  value). A factor's unused levels keep their colours, so that plots of
  subsets match, but are left out of the legend. By default, the groups
  that [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) and
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md) store
  in their `groups` attribute, if any.

- data:

  A data frame of sample metadata, with one row per sample (matched by
  row name, else in order), when `groups` names one of its columns.

- axes:

  Which two axes to plot, by number, as in
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md).
  Default `c(1, 2)`. A result with only one axis (such as
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) of two
  groups) is drawn as a strip chart, with one strip per group (in which
  case `axes` is not used).

- loadings:

  How many of the features in the `loadings` attribute to draw, those
  furthest from the origin on the plotted axes first: a number, `Inf`
  (or `TRUE`) for all, or `0` (or `FALSE`) for none. Default `10`. Drawn
  as arrows for
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) and
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md), and
  as labelled points for the species scores of
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) and
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md).

- biplot:

  Whether to draw the arrows of the constraints in the `biplot`
  attribute of
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) and
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md). Those of
  factor levels are left out when `centroids` are drawn. Default `TRUE`.

- centroids:

  Whether to label the points in the `centroids` attribute: the centroid
  of each factor level of a constrained ordination, or of each group of
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) and
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md).
  Default `TRUE`. Those of the plotted `groups` are labelled with the
  group and drawn in its colour.

- target:

  For
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md),
  whether to draw an arrow from each sample to its place in the `target`
  attribute, the ordination it was fitted to. Default `TRUE`.

- legend:

  Where to put the legend of `groups`: a keyword for
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html), such
  as `"topright"` (the default) or `"bottomleft"`, or `NULL` for none.

- col:

  The points' colours. Without `groups`, passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) as given, so
  that `col = iris$Species` works as for any matrix. With `groups`, one
  colour per group (in order, or named by group); for numeric `groups`,
  the colours of the gradient, from low to high. Defaults to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)'s
  `"Dark 3"` palette for groups and its `"Viridis"` palette for a
  gradient.

- pch:

  The points' symbols, as for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Default `1`
  (open circles), or `16` (filled circles) with `groups`.

- ...:

  Passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `sub`, `xlab`, `xlim`, `asp` or `cex`. They replace
  the defaults made from `x`. To leave out the axes, use `xaxt = "n"`
  and `yaxt = "n"`, since `axes` chooses which axes to plot.

## Value

`x`, invisibly.

## Details

The plot is made of what `x` holds:

- **Title**: the method's name, as
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  and
  [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  give it, and for an ordination of distances, the distances' metric,
  from the records of how `x` was made (see
  [ecodive_steps](https://cmmr.github.io/ecodive/reference/ecodive_steps.md)):
  such as "Principal Coordinates Analysis of Bray-Curtis Dissimilarity".

- **Axis labels**: the axes' names, with the percentage of the variance
  each explains (the `prop` attribute), for methods that have one.

- **Subtitle**: how well the method fit, for methods that measure it:
  the `stress` of
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) and
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), the
  `kl_divergence` of
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), the
  \\R^2\\ (`r_squared` and `adj_r_squared`) of a constrained ordination,
  `r2x` and `r2y` of
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md), and
  the `ss` of
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md);
  and for [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), their
  settings.

The axes have the same scale (`asp = 1`), so that distances on the plot
are true to the ordination; give `asp = NA` to let each axis fill the
plot.

Arrows (loadings and constraints) are scaled to fit among the samples,
as in any biplot: only their directions and relative lengths have
meaning. Loadings and constraints are scaled separately. Centroids, and
the species scores of
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) and
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) (`vegan`'s
scaling 2), are drawn where they are.

Base graphics draw on the current device, so further points, text or
arrows can be added on top. Use
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md)'s plot
method, `plot(fit, add = TRUE)`, to add fitted variables.

## See also

[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and its
plot method, `plot(fit, add = TRUE)`, to add fitted variables to an
ordination plot.

## Examples

``` r
    # Colour by a column of metadata, with a legend.
    ord <- pcoa(bray(iris[, 1:4]))
    plot(ord, 'Species', iris)


    # The third axis instead of the second.
    plot(ord, 'Species', iris, axes = c(1, 3))


    # A PCA biplot: the loadings as arrows.
    plot(pca(iris[, 1:4], scale = TRUE), iris$Species)


    # A constrained ordination: centroids of each species, and an arrow
    # for the numeric constraint.
    ord <- rda(iris[, 1:4], data = iris, formula = ~ Species + Sepal.Width)
    plot(ord, 'Species', iris)


    # Supervised methods colour by their own groups.
    plot(lda(iris[, 1:4], 'Species', iris))


    # Any argument of plot() still works.
    plot(nmds(bray(iris[, 1:4]), seed = 1), col = 'grey40', main = 'NMDS')
```
