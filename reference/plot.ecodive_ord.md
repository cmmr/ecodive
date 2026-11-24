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
  scaling = NULL,
  links = NULL,
  type = "samples",
  numerator = NULL,
  denominator = NULL,
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
  (or `TRUE`) for all, or `0` (or `FALSE`) for none. Default `10`. In a
  rank plot, how many features to name, those with the largest loadings
  (among `numerator` and `denominator`, when given). Drawn as arrows for
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) and
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md), and
  as labelled points for the species scores of
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) and
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md). With
  `scaling = 'contribution'`, those contributing most to the plane
  first.

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

- scaling:

  How the samples and features of a biplot are scaled against each other
  (see "Compositional biplots" and "Scalings of rda(), cca() and ca()"
  below):

  - For [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
    [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) and
    [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md):
    `'form'` (or `1`), the default; `'covariance'` (or `2`); or
    `'contribution'`.

  - For [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
    [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
    [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md): `vegan`'s
    scaling `1` (or `'sites'`), `2` (or `'species'`), the default, as
    the result is stored, or `3` (or `'symmetric'`).

  Default `NULL`: the default of the method, which draws the result as
  it is stored. Other ordinations take no `scaling`. Names may be
  abbreviated (`'cov'`), and `'form'` and `'covariance'` also stand for
  1 and 2 with
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md).

- links:

  Features to join with dashed links in a covariance biplot
  (`scaling = 'covariance'`) of
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) or
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md):
  a vector of at least two features (names or positions in the
  `loadings`), every pair of which is linked, or a two-column matrix of
  pairs of them. Each link is labelled with its length, which
  approximates the standard deviation of the log ratio of its two
  features (see "Compositional biplots"). The linked features are drawn
  whatever `loadings` is. Default `NULL`, for none.

- type:

  `"samples"` (the default): the samples on two axes, with the overlays
  below. `"rank"`: a rank plot of the feature loadings on axis `axes[1]`
  (see **Rank plots**), for ordinations of a count table with a
  `loadings` attribute. Any other value is passed on to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html),
  so `type = "n"` sets up the plot without drawing the samples, as for
  any matrix.

- numerator, denominator:

  For `type = "rank"`: the features to highlight as the numerator (red)
  and denominator (blue) of a log ratio, such as those that
  [`log_ratio_features()`](https://cmmr.github.io/ecodive/reference/log_ratio_features.md)
  picks: names (rows of `loadings`), indices, or a logical vector, as
  for
  [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md).
  They must not overlap.

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

- **Title**: the ordination's title (see
  [ecodive_display](https://cmmr.github.io/ecodive/reference/ecodive_display.md)):
  the method's name, as
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  and
  [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  give it, and for an ordination of distances, the distances' metric,
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
scaling 2), are drawn where they are. The covariance and contribution
biplots below draw the features as rays from the origin, at their own
coordinates.

**Compositional biplots**

A PCA of clr-transformed counts, `pca(counts, norm = 'clr')`, is
Aitchison and Greenacre's (2002) biplot of compositional data, which
Gloor et al. (2017) recommend for microbiome tables. With \\S\\ the
scores, \\V\\ the unit `loadings` and \\\lambda\\ each axis's variance
(`eig`), each `scaling` draws the samples (rows) and features (columns)
on the plotted axes as:

|  |  |  |  |
|----|----|----|----|
| `scaling` | Samples | Features | What reads off the plot |
| `'form'` | \\S\\ | \\V\\, scaled to fit | distances between samples approximate their Aitchison distances |
| `'covariance'` | \\S / \sqrt{\lambda}\\ (unit variance) | \\V \sqrt{\lambda}\\ | rays' lengths approximate the clr values' standard deviations; links, the log ratios' |
| `'contribution'` | \\S\\ | \\V\\, times a round number | features far from the origin contribute most to the plane |

In each, a sample's coordinates times a feature's (divided by the
features' multiplier) are what the PCA fits for that sample and feature,
\\S V^T\\. Aitchison and Greenacre's rules for the covariance biplot,
which hold for clr data (`norm = 'clr'`, `scale = FALSE`):

- A **link** between the tips of two rays approximates the standard
  deviation of the log ratio of those two features, \\\log(x_j / x_k)\\;
  the label on each link drawn by `links` is its length on the plane.
  Over all axes the lengths are exact, so the plane approximates them
  **from below**: a short link can hide a large standard deviation on
  the other axes, and a long one cannot. With
  [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  `sqrt(vlr(counts))` gives the exact values.

- Links are least-squares fits of the log ratios, as rays are of the clr
  values, and explain the same share of their variance. The subtitle
  gives that share (the plane's percentage of the variance). When it is
  under 50%, drawing `links` warns that their lengths may greatly
  underestimate the standard deviations.

- The cosine of the angle between two links approximates the correlation
  of their log ratios: perpendicular links are uncorrelated log ratios,
  and parallel ones correlated.

- Three features whose tips lie on a line, at distances \\\lambda\\ (A
  to B) and \\\mu\\ (B to C), suggest that the log contrast \\\mu \log
  A + \lambda \log C - (\lambda + \mu) \log B\\ is nearly constant
  across samples; four that form a parallelogram A, B, C, D, that \\\log
  A - \log B + \log C - \log D\\ is.

Links are log ratios only for clr data. Of a
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) with another
`norm`, or with `scale = TRUE`, they show standard deviations of
differences of the transformed values, and drawing them warns so.
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md)'s and
[`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md)'s
completed rclr tables are row-centred as clr data are, so their links
are log ratios, of the completed values; for
[`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md),
those of the tree's nodes (branch-length-weighted counts). Their `prop`,
and so the subtitle's share, is relative to the completed low-rank
table, which overstates the fit to the data.

The **contribution biplot** (Greenacre 2013) draws each feature at its
contribution coordinates, which for an unweighted PCA such as
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)'s are the
unit loadings \\V\\ themselves: the squared coordinate on an axis is the
feature's share of that axis's variance. It differs from `'form'` only
in using a fixed, round multiplier, given in the subtitle ("feature
coordinates \\\times\\ 20"), so that lengths can be compared between
plots, and in highlighting and labelling the features whose contribution
to either plotted axis is above the average, \\1/p\\ for \\p\\ features.
`loadings` then counts features by their contribution to the plane,
\\(\lambda_a v\_{ja}^2 + \lambda_b v\_{jb}^2) / (\lambda_a +
\lambda_b)\\, rather than by length.

**Scalings of
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md)**

[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) store `vegan`'s
scaling 2 (species scaled by their eigenvalues), drawn by default.
`scaling = 1` scales the sites by the eigenvalues instead (so that
distances between sites approximate their Euclidean or chi-square
distances), and `scaling = 3` splits the scaling between the two, as
`vegan::scores(x, scaling = 1)` and `3` do: each axis's site scores,
centroids and constraints are multiplied by \\f\\ or \\\sqrt{f}\\, and
its species scores divided by them, where \\f = \sqrt{\lambda /
\mathrm{trace}}\\ for
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) and
\\\sqrt{\lambda}\\ for
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md). The subtitle
names the scaling when one is given.

Base graphics draw on the current device, so further points, text or
arrows can be added on top. Use
[`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md)'s plot
method, `plot(fit, add = TRUE)`, to add fitted variables.

**Rank plots**

With `type = "rank"`, each feature is a bar from 0 to its loading on
axis `axes[1]`, sorted from the smallest at the left to the largest at
the right (ties in the order of `loadings`), as Qurro's rank plot of
feature loadings (Fedarko et al. 2020,
[doi:10.1093/nargab/lqaa023](https://doi.org/10.1093/nargab/lqaa023) ).
For a compositional ordination such as
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), the
features at the two ends pull the samples apart along the axis, and
their log ratio (Morton et al. 2019,
[doi:10.1038/s41467-019-10656-5](https://doi.org/10.1038/s41467-019-10656-5)
) follows the axis without depending on the samples' totals: pick them
with
[`log_ratio_features()`](https://cmmr.github.io/ecodive/reference/log_ratio_features.md)
and plot it with
[`plot.ecodive_log_ratio()`](https://cmmr.github.io/ecodive/reference/plot.ecodive_log_ratio.md).
The bars are grey; `numerator` and `denominator` colour theirs, with a
legend in the top-left corner (`legend` moves it, and `NULL` leaves it
out). The title is "Rank plot of" the ordination's title. The loadings
of
[`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md)
include the tree's internal nodes. `scaling` and `links` don't apply:
each scaling multiplies an axis's loadings by one positive number, which
leaves their ranking unchanged.

## References

Aitchison, J., & Greenacre, M. (2002). Biplots of compositional data.
*Journal of the Royal Statistical Society: Series C (Applied
Statistics)*, 51(4), 375-392.
[doi:10.1111/1467-9876.00275](https://doi.org/10.1111/1467-9876.00275)

Greenacre, M. (2013). Contribution biplots. *Journal of Computational
and Graphical Statistics*, 22(1), 107-122.
[doi:10.1080/10618600.2012.702494](https://doi.org/10.1080/10618600.2012.702494)

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
(2017). Microbiome datasets are compositional: and this is not optional.
*Frontiers in Microbiology*, 8, 2224.
[doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)

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


    # A compositional (covariance) biplot, with links whose lengths
    # approximate the standard deviations of log ratios.
    ord <- pca(ex_counts, norm = 'clr', pseudocount = 1)
    plot(ord, scaling = 'covariance',
         links = c('Streptococcus', 'Bacteroides', 'Corynebacterium'))


    # The features that contribute most to the plane.
    plot(ord, scaling = 'contribution')


    # vegan's scaling 1 of an RDA.
    plot(rda(iris[, 1:4]), iris$Species, scaling = 1)


    # A constrained ordination: centroids of each species, and an arrow
    # for the numeric constraint.
    ord <- rda(iris[, 1:4], data = iris, formula = ~ Species + Sepal.Width)
    plot(ord, 'Species', iris)


    # Supervised methods colour by their own groups.
    plot(lda(iris[, 1:4], 'Species', iris))


    # The features ranked by their loadings on the first axis.
    plot(pca(iris[, 1:4], scale = TRUE), type = 'rank')


    # Any argument of plot() still works.
    plot(nmds(bray(iris[, 1:4]), seed = 1), col = 'grey40', main = 'NMDS')
```
