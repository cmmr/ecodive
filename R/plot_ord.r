# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot an ordination
#'
#' Draws the samples on two axes of an ordination from [pcoa()], [pca()],
#' [dbrda()], [procrustes()] or any other of ecodive's ordination methods,
#' with the title, axis labels and subtitle taken from the result's
#' attributes. The overlays the result holds are drawn on top of the samples:
#' feature loadings, the arrows of numeric constraints, factor or group
#' centroids, and the target of a Procrustes rotation. Nothing is recomputed:
#' this is a quick look at results already in hand.
#'
#' @param x  An ordination result: a matrix of class `ecodive_ord`, from
#'        [dist_ord()], [count_ord()], [pcoa()], [mds()], [nmds()],
#'        [isomap()], [tsne()], [umap()], [dbrda()], [capscale()], [rda()],
#'        [cca()], [ca()], [pca()], [rpca()], [rpca_phylo()], [lda()],
#'        [plsda()], [oplsda()], or [procrustes()].
#'
#' @param groups  Colours the samples by group, with a legend: the name of a
#'        column of `data`, or a vector with one value per sample (matched to
#'        the samples by name when both are named, as for [permanova()]). A
#'        factor, character or logical vector gives each value its own
#'        colour, in the order of a factor's levels; a numeric vector gives a
#'        colour gradient (convert it with `factor()` for one colour per
#'        value). A factor's unused levels keep their colours, so that
#'        plots of subsets match, but are left out of the legend. By
#'        default, the groups that [lda()], [plsda()] and [oplsda()] store in
#'        their `groups` attribute, if any.
#'
#' @param data  A data frame of sample metadata, with one row per sample
#'        (matched by row name, else in order), when `groups` names one of
#'        its columns.
#'
#' @param axes  Which two axes to plot, by number, as in [envfit()]. Default
#'        `c(1, 2)`. A result with only one axis (such as [lda()] of two
#'        groups) is drawn as a strip chart, with one strip per group (in
#'        which case `axes` is not used).
#'
#' @param loadings  How many of the features in the `loadings` attribute to
#'        draw, those furthest from the origin on the plotted axes first: a
#'        number, `Inf` (or `TRUE`) for all, or `0` (or `FALSE`) for none.
#'        Default `10`. Drawn as arrows for [pca()], [rpca()],
#'        [rpca_phylo()], [rda()], [plsda()] and [oplsda()], and as labelled
#'        points for the species scores of [ca()] and [cca()]. With
#'        `scaling = 'contribution'`, those contributing most to the plane
#'        first.
#'
#' @param biplot  Whether to draw the arrows of the constraints in the
#'        `biplot` attribute of [dbrda()], [capscale()], [rda()] and [cca()].
#'        Those of factor levels are left out when `centroids` are drawn.
#'        Default `TRUE`.
#'
#' @param centroids  Whether to label the points in the `centroids`
#'        attribute: the centroid of each factor level of a constrained
#'        ordination, or of each group of [plsda()] and [oplsda()]. Default
#'        `TRUE`. Those of the plotted `groups` are labelled with the group
#'        and drawn in its colour.
#'
#' @param target  For [procrustes()], whether to draw an arrow from each
#'        sample to its place in the `target` attribute, the ordination it was
#'        fitted to. Default `TRUE`.
#'
#' @param legend  Where to put the legend of `groups`: a keyword for
#'        [graphics::legend()], such as `"topright"` (the default) or
#'        `"bottomleft"`, or `NULL` for none.
#'
#' @param col  The points' colours. Without `groups`, passed to `plot()` as
#'        given, so that `col = iris$Species` works as for any matrix. With
#'        `groups`, one colour per group (in order, or named by group); for
#'        numeric `groups`, the colours of the gradient, from low to high.
#'        Defaults to [grDevices::hcl.colors()]'s `"Dark 3"` palette for
#'        groups and its `"Viridis"` palette for a gradient.
#'
#' @param pch  The points' symbols, as for `plot()`. Default `1` (open
#'        circles), or `16` (filled circles) with `groups`.
#'
#' @param scaling  How the samples and features of a biplot are scaled
#'        against each other (see "Compositional biplots" and "Scalings of
#'        rda(), cca() and ca()" below):
#'   * For [pca()], [rpca()] and [rpca_phylo()]: `'form'` (or `1`), the
#'     default; `'covariance'` (or `2`); or `'contribution'`.
#'   * For [rda()], [cca()] and [ca()]: `vegan`'s scaling `1` (or
#'     `'sites'`), `2` (or `'species'`), the default, as the result is
#'     stored, or `3` (or `'symmetric'`).
#'
#'   Default `NULL`: the default of the method, which draws the result as it
#'   is stored. Other ordinations take no `scaling`. Names may be
#'   abbreviated (`'cov'`), and `'form'` and `'covariance'` also stand for 1
#'   and 2 with [rda()], [cca()] and [ca()].
#'
#' @param links  Features to join with dashed links in a covariance biplot
#'        (`scaling = 'covariance'`) of [pca()], [rpca()] or [rpca_phylo()]:
#'        a vector of at least two features (names or positions in the
#'        `loadings`), every pair of which is linked, or a two-column matrix
#'        of pairs of them. Each link is labelled with its length, which
#'        approximates the standard deviation of the log ratio of its two
#'        features (see "Compositional biplots"). The linked features are
#'        drawn whatever `loadings` is. Default `NULL`, for none.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim`, `asp` or `cex`. They replace the defaults
#'        made from `x`. To leave out the axes, use `xaxt = "n"` and `yaxt =
#'        "n"`, since `axes` chooses which axes to plot.
#'
#' @details
#' The plot is made of what `x` holds:
#'
#' * **Title**: the ordination's title (see [ecodive_display]): the
#'   method's name, as [list_dist_ords()] and [list_count_ords()] give it,
#'   and for an ordination of distances, the distances' metric, such as
#'   "Principal Coordinates Analysis of Bray-Curtis Dissimilarity".
#' * **Axis labels**: the axes' names, with the percentage of the variance
#'   each explains (the `prop` attribute), for methods that have one.
#' * **Subtitle**: how well the method fit, for methods that measure it: the
#'   `stress` of [mds()] and [nmds()], the `kl_divergence` of [tsne()], the
#'   \eqn{R^2} (`r_squared` and `adj_r_squared`) of a constrained
#'   ordination, `r2x` and `r2y` of [oplsda()], and the `ss` of
#'   [procrustes()]; and for [pcoa()], [isomap()] and [umap()], their
#'   settings.
#'
#' The axes have the same scale (`asp = 1`), so that distances on the plot
#' are true to the ordination; give `asp = NA` to let each axis fill the
#' plot.
#'
#' Arrows (loadings and constraints) are scaled to fit among the samples, as
#' in any biplot: only their directions and relative lengths have meaning.
#' Loadings and constraints are scaled separately. Centroids, and the
#' species scores of [ca()] and [cca()] (`vegan`'s scaling 2), are drawn
#' where they are. The covariance and contribution biplots below draw the
#' features as rays from the origin, at their own coordinates.
#'
#' **Compositional biplots**
#'
#' A PCA of clr-transformed counts, `pca(counts, norm = 'clr')`, is
#' Aitchison and Greenacre's (2002) biplot of compositional data, which
#' Gloor et al. (2017) recommend for microbiome tables. With \eqn{S} the
#' scores, \eqn{V} the unit `loadings` and \eqn{\lambda} each axis's
#' variance (`eig`), each `scaling` draws the samples (rows) and features
#' (columns) on the plotted axes as:
#'
#' | `scaling` | Samples | Features | What reads off the plot |
#' |:----------|:--------|:---------|:------------------------|
#' | `'form'` | \eqn{S} | \eqn{V}, scaled to fit | distances between samples approximate their Aitchison distances |
#' | `'covariance'` | \eqn{S / \sqrt{\lambda}} (unit variance) | \eqn{V \sqrt{\lambda}} | rays' lengths approximate the clr values' standard deviations; links, the log ratios' |
#' | `'contribution'` | \eqn{S} | \eqn{V}, times a round number | features far from the origin contribute most to the plane |
#'
#' In each, a sample's coordinates times a feature's (divided by the
#' features' multiplier) are what the PCA fits for that sample and feature,
#' \eqn{S V^T}. Aitchison and Greenacre's rules for the covariance biplot,
#' which hold for clr data (`norm = 'clr'`, `scale = FALSE`):
#'
#' * A **link** between the tips of two rays approximates the standard
#'   deviation of the log ratio of those two features, \eqn{\log(x_j /
#'   x_k)}; the label on each link drawn by `links` is its length on the
#'   plane. Over all axes the lengths are exact, so the plane approximates
#'   them **from below**: a short link can hide a large standard deviation
#'   on the other axes, and a long one cannot. With [vlr()],
#'   `sqrt(vlr(counts))` gives the exact values.
#' * Links are least-squares fits of the log ratios, as rays are of the clr
#'   values, and explain the same share of their variance. The subtitle
#'   gives that share (the plane's percentage of the variance). When it is
#'   under 50%, drawing `links` warns that their lengths may greatly
#'   underestimate the standard deviations.
#' * The cosine of the angle between two links approximates the
#'   correlation of their log ratios: perpendicular links are uncorrelated
#'   log ratios, and parallel ones correlated.
#' * Three features whose tips lie on a line, at distances \eqn{\lambda}
#'   (A to B) and \eqn{\mu} (B to C), suggest that the log contrast
#'   \eqn{\mu \log A + \lambda \log C - (\lambda + \mu) \log B} is
#'   nearly constant across samples; four that form a parallelogram A, B,
#'   C, D, that \eqn{\log A - \log B + \log C - \log D} is.
#'
#' Links are log ratios only for clr data. Of a `pca()` with another `norm`,
#' or with `scale = TRUE`, they show standard deviations of differences of
#' the transformed values, and drawing them warns so. [rpca()]'s and
#' [rpca_phylo()]'s completed rclr tables are row-centred as clr data are,
#' so their links are log ratios, of the completed values; for
#' [rpca_phylo()], those of the tree's nodes (branch-length-weighted
#' counts). Their `prop`, and so the subtitle's share, is relative to the
#' completed low-rank table, which overstates the fit to the data.
#'
#' The **contribution biplot** (Greenacre 2013) draws each feature at its
#' contribution coordinates, which for an unweighted PCA such as
#' `pca()`'s are the unit loadings \eqn{V} themselves: the squared
#' coordinate on an axis is the feature's share of that axis's variance. It
#' differs from `'form'` only in using a fixed, round multiplier, given in
#' the subtitle ("feature coordinates \eqn{\times} 20"), so that lengths
#' can be compared between plots, and in highlighting and labelling the
#' features whose contribution to either plotted axis is above the average,
#' \eqn{1/p} for \eqn{p} features. `loadings` then counts features by
#' their contribution to the plane, \eqn{(\lambda_a v_{ja}^2 + \lambda_b
#' v_{jb}^2) / (\lambda_a + \lambda_b)}, rather than by length.
#'
#' **Scalings of `rda()`, `cca()` and `ca()`**
#'
#' [rda()], [cca()] and [ca()] store `vegan`'s scaling 2 (species scaled by
#' their eigenvalues), drawn by default. `scaling = 1` scales the sites by
#' the eigenvalues instead (so that distances between sites approximate
#' their Euclidean or chi-square distances), and `scaling = 3` splits the
#' scaling between the two, as `vegan::scores(x, scaling = 1)` and `3` do:
#' each axis's site scores, centroids and constraints are multiplied by
#' \eqn{f} or \eqn{\sqrt{f}}, and its species scores divided by them,
#' where \eqn{f = \sqrt{\lambda / \mathrm{trace}}} for [rda()] and
#' \eqn{\sqrt{\lambda}} for [cca()] and [ca()]. The subtitle names the
#' scaling when one is given.
#'
#' Base graphics draw on the current device, so further points, text or
#' arrows can be added on top. Use [envfit()]'s plot method, `plot(fit, add =
#' TRUE)`, to add fitted variables.
#'
#' @return `x`, invisibly.
#'
#' @references
#' Aitchison, J., & Greenacre, M. (2002). Biplots of compositional data.
#' *Journal of the Royal Statistical Society: Series C (Applied
#' Statistics)*, 51(4), 375-392. \doi{10.1111/1467-9876.00275}
#'
#' Greenacre, M. (2013). Contribution biplots. *Journal of Computational
#' and Graphical Statistics*, 22(1), 107-122.
#' \doi{10.1080/10618600.2012.702494}
#'
#' Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
#' (2017). Microbiome datasets are compositional: and this is not optional.
#' *Frontiers in Microbiology*, 8, 2224. \doi{10.3389/fmicb.2017.02224}
#'
#' @seealso [envfit()] and its plot method, `plot(fit, add = TRUE)`, to add
#'   fitted variables to an ordination plot.
#'
#' @export
#' @examples
#'     # Colour by a column of metadata, with a legend.
#'     ord <- pcoa(bray(iris[, 1:4]))
#'     plot(ord, 'Species', iris)
#'
#'     # The third axis instead of the second.
#'     plot(ord, 'Species', iris, axes = c(1, 3))
#'
#'     # A PCA biplot: the loadings as arrows.
#'     plot(pca(iris[, 1:4], scale = TRUE), iris$Species)
#'
#'     # A compositional (covariance) biplot, with links whose lengths
#'     # approximate the standard deviations of log ratios.
#'     ord <- pca(ex_counts, norm = 'clr', pseudocount = 1)
#'     plot(ord, scaling = 'covariance',
#'          links = c('Streptococcus', 'Bacteroides', 'Corynebacterium'))
#'
#'     # The features that contribute most to the plane.
#'     plot(ord, scaling = 'contribution')
#'
#'     # vegan's scaling 1 of an RDA.
#'     plot(rda(iris[, 1:4]), iris$Species, scaling = 1)
#'
#'     # A constrained ordination: centroids of each species, and an arrow
#'     # for the numeric constraint.
#'     ord <- rda(iris[, 1:4], data = iris, formula = ~ Species + Sepal.Width)
#'     plot(ord, 'Species', iris)
#'
#'     # Supervised methods colour by their own groups.
#'     plot(lda(iris[, 1:4], 'Species', iris))
#'
#'     # Any argument of plot() still works.
#'     plot(nmds(bray(iris[, 1:4]), seed = 1), col = 'grey40', main = 'NMDS')
plot.ecodive_ord <- function (
    x, groups = NULL, data = NULL, axes = c(1, 2), loadings = 10,
    biplot = TRUE, centroids = TRUE, target = TRUE, legend = 'topright',
    col = NULL, pch = NULL, scaling = NULL, links = NULL, ... ) {

  args     <- list(...)
  att      <- attributes(x)
  pos      <- plot_legend_pos(legend)
  loadings <- plot_count(loadings, 'loadings')
  for (arg in c('biplot', 'centroids', 'target'))
    if (!isTRUE(get(arg)) && !isFALSE(get(arg)))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)

  if (ncol(x) > 1) axes <- plot_axes(axes, ncol(x))
  else if (ncol(x) == 0) stop('`x` has no axes to plot.', call. = FALSE)
  ord_scaling_arg(scaling)
  xy   <- unclass(x)
  attributes(xy) <- att[c('dim', 'dimnames')]

  g    <- ord_plot_groups(x, groups, data)
  look <- ord_plot_look(g, col, pch)

  main <- plot_title(x, ord_method_name(att$method))
  sub  <- ord_subtitle(x)
  labs <- ord_axis_labels(x)


  #--------------------------------------------------
  # One axis: a strip chart, one strip per group.
  #--------------------------------------------------
  if (ncol(x) == 1) {
    strip <- if (is.null(g$f)) rep(1L, nrow(x)) else match(as.integer(g$f), which(look$used))
    plot_do(
      graphics::plot.default, c(list(strip, xy[, 1]), args),
      list(
        col = look$col, pch = look$pch, xaxt = 'n',
        xlim = c(0.5, max(strip, 1L) + 0.5), main = main, sub = sub,
        xlab = '', ylab = labs[[1]] ))
    if (!is.null(g$f) && plot_x_axis(args))
      graphics::axis(1, at = seq_along(look$legend), labels = look$legend)
    if (!is.null(g$values))
      plot_legend(pos, look$legend, look$legend_col, look$legend_pch, g$title)
    return (invisible(x))
  }


  #--------------------------------------------------
  # Two axes: the samples, then the overlays.
  #--------------------------------------------------
  # The samples and features in the chosen scaling (see
  # ord_plot_scaling()): `load` holds the features to draw, in the plot's
  # coordinates.
  sc   <- ord_plot_scaling(x, axes, scaling, loadings, links)
  pts  <- sc$rows
  load <- sc$cols
  sub  <- paste(c(sub, sc$sub), collapse = '; ')
  if (!nzchar(sub)) sub <- NULL

  arr_con  <- if (biplot) ord_plot_biplot(att$biplot, axes, if (centroids) att$centroids)
  cen      <- if (centroids && !is.null(att$centroids) && nrow(att$centroids))
    sweep(att$centroids[, axes, drop = FALSE], 2, sc$row_mul, '*')
  tgt      <- if (target && !is.null(att$target))
    att$target[, axes, drop = FALSE]
  if (!is.null(arr_con)) arr_con <- sweep(arr_con, 2, sc$row_mul, '*')

  # The plot holds the samples, centroids, Procrustes targets, and the
  # features drawn where they are (species scores, and the rays of the
  # covariance and contribution biplots).
  all_xy <- rbind(
    pts, cen, tgt,
    if (sc$style == 'points') load,
    if (sc$style == 'rays') load * 1.1 )   # with room for the labels

  plot_do(
    graphics::plot.default, c(list(pts[, 1], pts[, 2]), args),
    list(
      col = look$col, pch = look$pch, asp = 1, main = main, sub = sub,
      xlab = labs[[axes[[1]]]], ylab = labs[[axes[[2]]]],
      xlim = range(all_xy[, 1], finite = TRUE),
      ylim = range(all_xy[, 2], finite = TRUE) ))

  if (!is.null(tgt))
    ord_plot_target(pts, tgt)

  ext <- c(max(abs(pts[, 1])), max(abs(pts[, 2])))

  if (!is.null(load) && nrow(load)) {
    if (sc$style == 'points') {
      graphics::points(load[, 1], load[, 2], pch = 4, col = PLOT_MUTED)
      plot_label(load[, 1], load[, 2], rownames(load), col = PLOT_MUTED, pos = 3, cex = 0.8)
    }
    else if (sc$style == 'rays') {
      ord_plot_rays(load, sc$col, sc$label)
    }
    else {
      plot_arrows(load, col = PLOT_MUTED)
    }
  }

  if (!is.null(sc$links))
    ord_plot_links(load, sc$links)

  if (!is.null(arr_con) && nrow(arr_con))
    plot_arrows(arr_con * plot_arrow_mul(arr_con, ext), col = PLOT_HIGHLIGHT)

  if (!is.null(cen))
    ord_plot_centroids(cen, g, look)

  if (!is.null(g$f) || !is.null(g$values))
    plot_legend(pos, look$legend, look$legend_col, look$legend_pch, g$title)

  invisible(x)
}



# Labels each centroid with a cross (see ord_centroid_labels()).
ord_plot_centroids <- function (cen, g, look) {
  lc <- ord_centroid_labels(rownames(cen), g, look)
  graphics::points(cen[, 1], cen[, 2], pch = 3, col = lc$col, cex = 1.2, lwd = 2)
  plot_label(cen[, 1], cen[, 2], lc$labels, col = lc$col, pos = 3, font = 2)
  invisible()
}


# The labels and colours of centroids named `nms`. The centroids of the
# plotted groups are labelled with the group and coloured as its samples:
# those of plsda() and oplsda() (named by group), and those of a
# constrained ordination's factor when `groups` names the same column of
# `data` (named by the column and the level pasted together, as in vegan).
# Others keep their names. Returns list(labels, col).
ord_centroid_labels <- function (nms, g, look) {

  lab <- nms
  col <- rep(PLOT_HIGHLIGHT, length(nms))

  if (!is.null(g$f)) {
    lev <- levels(g$f)
    i   <- match(lab, lev)
    if (!is.null(g$title) && all(is.na(i)))
      i <- match(lab, paste0(g$title, lev))
    hit <- !is.na(i)
    lab[hit] <- lev[i[hit]]
    col[hit] <- look$pal[i[hit]]
  }

  list(labels = lab, col = col)
}



# The method's name for the title: from the tables behind list_dist_ords()
# and list_count_ords().
ord_method_name <- function (method) {
  tab <- rbind(ORDINATIONS[c('id', 'name')], COUNT_ORDS[c('id', 'name')])
  i   <- match(method, tab$id)
  if (!is.na(i)) return (tab$name[[i]])
  if (identical(method, 'procrustes')) return ('Procrustes Rotation')
  if (is.null(method)) NULL else as.character(method)
}


# Axis labels: the column names, with each axis's percentage of the variance
# (`prop`) when known.
ord_axis_labels <- function (x) {
  labs <- colnames(x)
  if (is.null(labs)) labs <- paste('Axis', seq_len(ncol(x)))
  prop <- attr(x, 'prop', exact = TRUE)
  if (length(prop) == length(labs)) {
    ok <- !is.na(prop)
    labs[ok] <- sprintf('%s (%s%%)', labs[ok], format(round(100 * prop[ok], 1), nsmall = 1, trim = TRUE))
  }
  labs
}


# The subtitle: how well the method fit, or its settings, from attributes.
ord_subtitle <- function (x) {

  a <- function (nm) attr(x, nm, exact = TRUE)

  switch(
    as.character(a('method')),

    pcoa = if (!is.null(a('correction')) && !identical(a('correction'), 'none'))
      sprintf('%s correction, constant %s', plot_cap(a('correction')), plot_num(a('constant'))),

    mds  = ,
    nmds = sprintf('Stress = %s', plot_num(a('stress'))),

    isomap = if (!is.null(a('radius')) && !is.na(a('radius')))
      sprintf('Neighbours within a radius of %s', plot_num(a('radius')))
    else if (!is.null(a('n_neighbors')))
      sprintf('%d nearest neighbours', a('n_neighbors')),

    tsne = sprintf(
      'KL divergence = %s, perplexity = %s',
      plot_num(a('kl_divergence')), plot_num(a('perplexity')) ),

    umap = if (!is.null(a('n_neighbors')))
      sprintf('%d nearest neighbours', a('n_neighbors')),

    dbrda    = ,
    capscale = ,
    rda      = ,
    cca      = if (length(a('constrained_eig')))
      sprintf(
        'R-squared = %s, adj. R-squared = %s',
        plot_num(a('r_squared')), plot_num(a('adj_r_squared')) ),

    oplsda = sprintf('R2X = %s, R2Y = %s', plot_num(a('r2x')), plot_num(a('r2y'))),

    procrustes = sprintf(
      'Procrustes sum of squares = %s%s', plot_num(a('ss')),
      if (isTRUE(a('symmetric'))) ' (symmetric)' else '' ),

    NULL )
}



# The samples' groups, from `groups` and `data` (by per_sample()'s rule), or
# else from the `groups` attribute. Returns list(f = a factor (with NA as a
# level, if any), or NULL; values = numeric values, for a gradient, or NULL;
# numeric; title = the legend's title, or NULL).
ord_plot_groups <- function (x, groups, data) {
  plot_groups(nrow(x), rownames(x), groups, data, attr(x, 'groups', exact = TRUE))
}


# The groups of `n` samples named `labels` (or NULL), from `groups` and
# `data` (by per_sample()'s rule), or else `stored`. As ord_plot_groups().
plot_groups <- function (n, labels, groups, data, stored = NULL) {

  title <- NULL

  if (is.null(groups)) {
    v <- stored
    if (is.null(v)) return (list(f = NULL, values = NULL, numeric = FALSE, title = NULL))
  }
  else {
    idx <- NULL
    if (!is.null(data)) {
      if (!is.data.frame(data) && !is.matrix(data))
        stop('`data` must be a data.frame.', call. = FALSE)
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
      idx  <- match_rows(labels, n, dn, nrow(data), ref = '`x`')
      if (is.character(groups) && length(groups) == 1 && groups %in% names(data))
        title <- groups
    }
    v <- per_sample(groups, 'groups', data, idx, labels, n, ref = '`x`')
  }

  if (is.numeric(v) && !is.factor(v))
    return (list(f = NULL, values = as.double(v), numeric = TRUE, title = title))

  # A factor keeps its unused levels, for their colours.
  f <- if (is.factor(v)) v else factor(v)
  if (anyNA(f)) f <- addNA(f, ifany = TRUE)

  list(f = f, values = NULL, numeric = FALSE, title = title)
}



# Colours and symbols of the points, and of the legend. Returns list(col,
# pch (the points'), pal (one colour per level), used (which levels have
# samples), legend, legend_col, legend_pch). The legend has the used levels
# only.
ord_plot_look <- function (g, col, pch) {

  # No groups: the user's col and pch, as plot() takes them.
  if (is.null(g$f) && is.null(g$values))
    return (list(
      col = if (is.null(col)) graphics::par('fg') else col,
      pch = if (is.null(pch)) 1 else pch ))

  if (is.null(pch)) pch <- 16

  # A gradient for numeric groups, in 100 steps. Its legend shows a few
  # round values within the range.
  if (!is.null(g$values)) {
    v    <- g$values
    ramp <- if (is.null(col)) grDevices::hcl.colors(100, 'Viridis')
            else grDevices::colorRampPalette(col)(100)
    if (!any(is.finite(v)))
      return (list(col = rep(PLOT_MUTED, length(v)), pch = pch, legend = NULL))
    rng  <- range(v[is.finite(v)])
    wide <- rng[2] > rng[1]
    step <- if (wide) (v - rng[1]) / (rng[2] - rng[1]) else v * 0 + 0.5
    pcol <- ramp[1 + round(99 * step)]
    pcol[is.na(pcol)] <- PLOT_MUTED
    at   <- rng[1]
    if (wide) {
      at <- pretty(rng, 4)
      at <- at[at >= rng[1] & at <= rng[2]]
      if (length(at) < 2) at <- rng
    }
    lcol <- ramp[1 + round(99 * (if (wide) (at - rng[1]) / (rng[2] - rng[1]) else 0.5))]
    return (list(
      col = pcol, pch = pch,
      legend = format(at, trim = TRUE), legend_col = lcol, legend_pch = pch[[1]] ))
  }

  lev <- levels(g$f)
  k   <- length(lev)
  pal <- if (is.null(col)) plot_palette(k)
         else if (!is.null(names(col)) && all(lev[!is.na(lev)] %in% names(col)))
           col[ifelse(is.na(lev), NA, lev)]
         else rep_len(col, k)
  pal[is.na(pal)] <- PLOT_MUTED
  pal <- unname(pal)

  used <- tabulate(as.integer(g$f), k) > 0
  lpch <- if (length(pch) == k) pch[used] else pch[[1]]
  list(
    col = pal[as.integer(g$f)], pch = if (length(pch) == k) pch[as.integer(g$f)] else pch,
    pal = pal, used = used, legend = ifelse(is.na(lev), 'NA', lev)[used],
    legend_col = pal[used], legend_pch = lpch )
}



# The constraints' arrows on the plotted axes, leaving out those of factor
# levels when their centroids are drawn (vegan names both alike), or NULL.
ord_plot_biplot <- function (bp, axes, cen) {
  if (is.null(bp) || !is.matrix(bp) || nrow(bp) == 0 || ncol(bp) < max(axes))
    return (NULL)
  a <- bp[, axes, drop = FALSE]
  if (!is.null(cen)) a <- a[!rownames(a) %in% rownames(cen), , drop = FALSE]
  a
}


# procrustes(): an arrow from each sample to its place in the target (if
# far enough away to draw one).
ord_plot_target <- function (pts, tgt) {
  keep <- plot_arrow_keep(pts[, 1], pts[, 2], tgt[, 1], tgt[, 2])
  if (any(keep))
    graphics::arrows(
      pts[keep, 1], pts[keep, 2], tgt[keep, 1], tgt[keep, 2],
      length = 0.05, col = PLOT_MUTED )
  invisible()
}



#--------------------------------------------------
# Biplot scalings and links (see "Compositional biplots" in
# ?plot.ecodive_ord).
#--------------------------------------------------

# The methods whose results take each kind of `scaling`: unweighted PCAs
# (scores U S, unit loadings V, eig = s^2 / (n - 1)), and vegan's methods,
# stored in its scaling 2.
ORD_PCA_FAMILY   <- c('pca', 'rpca', 'rpca_phylo')
ORD_VEGAN_FAMILY <- c('rda', 'cca', 'ca')


# Checks `scaling`: NULL, a name, or vegan's number. Returns NULL, or the
# canonical name: 'form' (vegan's 1, "sites"), 'covariance' (2, "species"),
# 'symmetric' (3) or 'contribution'.
ord_scaling_arg <- function (scaling) {

  if (is.null(scaling)) return (NULL)

  canon <- c(
    form = 'form', covariance = 'covariance', contribution = 'contribution',
    symmetric = 'symmetric', sites = 'form', species = 'covariance' )

  if (is.numeric(scaling) && length(scaling) == 1 && isTRUE(scaling %in% 1:3))
    return (c('form', 'covariance', 'symmetric')[[scaling]])

  if (is.character(scaling) && length(scaling) == 1 && !is.na(scaling)) {
    i <- pmatch(tolower(scaling), names(canon))
    if (!is.na(i)) return (canon[[i]])
  }

  stop(
    "`scaling` must be 'form', 'covariance', 'contribution' or 'symmetric', ",
    "or vegan's scaling 1, 2 or 3 ('sites', 'species' or 'symmetric').",
    call. = FALSE )
}


# The samples and features of a biplot on `axes` in the chosen `scaling`,
# with the `n_load` features to draw and the pairs of `links`. Returns
# list(
#   rows     = the samples' coordinates,
#   cols     = the features to draw, in the plot's coordinates (times
#              `mult`), or NULL,
#   mult     = the multiplier of `cols` (1 when drawn where they are),
#   style    = 'arrows' (scaled to fit), 'rays' (from the origin, where they
#              are) or 'points' (species scores),
#   col      = each drawn feature's colour, label = whether it is labelled,
#   row_mul  = the per-axis multiplier of the rows (for centroids and
#              constraints), col_mul = the same for the features,
#   sub      = the subtitle's part, or NULL,
#   scaling  = the canonical scaling,
#   links    = NULL, or a data frame of the drawn links: `from` and `to`
#              (rows of `cols`) and `length` (planar) ).
#
# For pca(), rpca() and rpca_phylo(), with S the scores, V the unit
# loadings and l = sqrt(eig) on `axes`:
#   'form':         rows S,     cols V, times a factor that fits them in;
#   'covariance':   rows S / l, cols V l (true coordinates);
#   'contribution': rows S,     cols V, times a round factor (1, 2 or 5
#                   times a power of 10).
# so that rows %*% t(cols) / mult = S %*% t(V) in each. For rda(), cca()
# and ca(), stored in vegan's scaling 2, with f = sqrt(eig / trace) for rda
# and sqrt(eig) for cca and ca, vegan's scaling 1 multiplies the sites
# (and centroids and constraints) by f and the species by 1 / f, and its
# scaling 3 by sqrt(f) and 1 / sqrt(f).
ord_plot_scaling <- function (x, axes, scaling = NULL, n_load = 10, links = NULL) {

  att    <- attributes(x)
  method <- as.character(att$method)
  s      <- ord_scaling_arg(scaling)
  given  <- !is.null(s)
  pca_f  <- length(method) == 1 && method %in% ORD_PCA_FAMILY
  veg_f  <- length(method) == 1 && method %in% ORD_VEGAN_FAMILY

  if (given && !pca_f && !veg_f)
    stop(
      '`scaling` applies to the results of pca(), rpca(), rpca_phylo(), rda(), ',
      'cca() and ca(), not ', if (length(method) == 1) paste0(method, '()') else 'this one', '.',
      call. = FALSE )
  if (pca_f && identical(s, 'symmetric'))
    stop(
      "`scaling = 'symmetric'` (vegan's 3) is for rda(), cca() and ca(). ",
      "A ", method, "() biplot takes 'form', 'covariance' or 'contribution'.", call. = FALSE )
  if (veg_f && identical(s, 'contribution'))
    stop(
      "`scaling = 'contribution'` is for pca(), rpca() and rpca_phylo(). ",
      "An ", method, "() result takes vegan's scaling 1, 2 or 3.", call. = FALSE )
  if (!is.null(links)) {
    if (!pca_f)
      stop('`links` are drawn only for the results of pca(), rpca() and rpca_phylo().', call. = FALSE)
    if (!identical(s, 'covariance'))
      stop(
        "`links` needs `scaling = 'covariance'`: only there do link lengths ",
        "approximate the standard deviations of log ratios.", call. = FALSE )
  }

  # What the result holds is vegan's scaling 2, or else a form biplot.
  if (!given) s <- if (veg_f) 'covariance' else 'form'

  xy <- unclass(x)
  attributes(xy) <- att[c('dim', 'dimnames')]
  rows <- xy[, axes, drop = FALSE]
  eig  <- att$eig[axes]

  # Per-axis multipliers of the rows and of the features.
  rmul <- cmul <- c(1, 1)
  if (veg_f) {
    f    <- sqrt(if (method == 'rda') eig / att$trace else eig)
    rmul <- switch(s, form = f, covariance = c(1, 1), symmetric = sqrt(f))
    cmul <- 1 / rmul
  }
  else if (pca_f && s == 'covariance') {
    rmul <- 1 / sqrt(eig)
    cmul <- sqrt(eig)
  }
  rows <- sweep(rows, 2, rmul, '*')

  style <- if (veg_f && method != 'rda') 'points'
           else if (pca_f && s != 'form') 'rays'
           else 'arrows'

  # The features: those furthest from the origin (or, for 'contribution',
  # those contributing most to the plane), and those that `links` names.
  load <- att$loadings
  if (is.null(load) || !is.matrix(load) || nrow(load) == 0 || ncol(load) < max(axes))
    load <- NULL
  pairs <- NULL
  cols  <- NULL
  sel   <- integer(0)
  if (!is.null(load)) {
    v <- load[, axes, drop = FALSE]
    if (is.null(rownames(v))) rownames(v) <- as.character(seq_len(nrow(v)))
    if (!is.null(links)) pairs <- ord_link_pairs(links, rownames(v))
    v <- sweep(v, 2, cmul, '*')
    score <- if (s == 'contribution') as.vector(v ^ 2 %*% eig) / sum(eig) else rowSums(v ^ 2)
    if (n_load > 0) sel <- utils::head(order(score, decreasing = TRUE), n_load)
    sel  <- c(sel, setdiff(unique(as.vector(pairs)), sel))
    cols <- v[sel, , drop = FALSE]
  }

  # Arrows scaled to fit among the samples; contribution coordinates, by a
  # round factor that fits them.
  mult <- 1
  if (!is.null(cols) && nrow(cols) && (style == 'arrows' || s == 'contribution')) {
    ext  <- c(max(abs(rows[, 1])), max(abs(rows[, 2])))
    mult <- plot_arrow_mul(cols, ext)
    if (s == 'contribution') mult <- ord_round_mult(mult)
    cols <- cols * mult
  }

  # Colours and labels: in a contribution biplot, the features above the
  # average contribution (1/p) to either axis are highlighted and labelled.
  nc    <- if (is.null(cols)) 0L else nrow(cols)
  colv  <- rep(PLOT_MUTED, nc)
  label <- rep(TRUE, nc)
  if (s == 'contribution' && nc) {
    hi <- unname(rowSums(load[sel, axes, drop = FALSE] ^ 2 > 1 / nrow(load)) > 0)
    colv[hi] <- PLOT_HIGHLIGHT
    label    <- hi
  }

  # Links: their planar lengths, and the warnings about their meaning.
  lk <- NULL
  if (!is.null(pairs)) {
    from <- match(pairs[, 1], sel)
    to   <- match(pairs[, 2], sel)
    lk   <- data.frame(
      from   = from,
      to     = to,
      length = sqrt(rowSums((cols[from, , drop = FALSE] - cols[to, , drop = FALSE]) ^ 2)) )
    if (method == 'pca') {
      msg <- ord_links_norm_msg(att$norm, att$scale)
      if (!is.null(msg)) warning(msg, call. = FALSE)
    }
    plane <- sum(att$prop[axes])
    if (isTRUE(plane < 0.5))
      warning(
        sprintf('The two axes explain %s%% of the variance, ', format(round(100 * plane))),
        'so link lengths may greatly underestimate the log-ratio standard deviations.',
        call. = FALSE )
  }

  sub <- NULL
  if (given && pca_f) {
    plane <- sum(att$prop[axes])
    sub   <- switch(
      s,
      form         = 'form biplot',
      covariance   = paste0(
        'covariance biplot',
        if (!is.na(plane)) sprintf(', plane explains %s%% of the variance', format(round(100 * plane)))),
      contribution = sprintf(
        'contribution biplot, feature coordinates %s %s', print_sym()$times,
        format(mult, scientific = mult < 1e-4 || mult >= 1e5)) )
  }
  else if (given && veg_f) {
    sub <- sprintf('scaling %d', match(s, c('form', 'covariance', 'symmetric')))
  }

  list(
    rows = rows, cols = cols, mult = mult, style = style, col = colv, label = label,
    row_mul = rmul, col_mul = cmul, sub = sub, scaling = s, links = lk )
}


# `m` rounded down to 1, 2 or 5 times a power of 10.
ord_round_mult <- function (m) {
  if (!is.finite(m) || m <= 0) return (1)
  k <- floor(log10(m))
  d <- m / 10 ^ k
  d <- if (d >= 5 - 1e-9) 5 else if (d >= 2 - 1e-9) 2 else 1
  as.numeric(format(d * 10 ^ k, digits = 1))
}


# Why the links of a pca() result in the transform `norm` (scaled or not)
# are not log ratios, or NULL when they are.
ord_links_norm_msg <- function (norm, scale) {
  scaled <- isTRUE(scale)
  if (identical(norm, 'clr') && !scaled) return (NULL)
  if (!is.character(norm) || length(norm) != 1)
    return ('links of this pca() result may not show standard deviations of log ratios.')
  what <- if (scaled) paste('standardized', norm, 'values')
          else switch(
            norm,
            none    = 'the untransformed values',
            percent = 'relative abundances',
            rclr    = 'rclr values, whose zeros stay 0',
            paste(norm, 'values') )
  sprintf(
    "links of `pca(norm = '%s'%s)` show standard deviations of differences of %s, not of log ratios.",
    norm, if (scaled) ', scale = TRUE' else '', what )
}


# The pairs of features that `links` names, among `nms`: every pair of a
# vector of features (names or positions), or the rows of a two-column
# matrix of pairs. Returns a two-column integer matrix of positions.
ord_link_pairs <- function (links, nms) {

  bad <- function ()
    stop(
      '`links` must be a vector of at least two features (names or positions ',
      'in the loadings), or a two-column matrix of pairs of them.', call. = FALSE )

  is_mat <- is.matrix(links) || is.data.frame(links)
  if (is.data.frame(links)) links <- as.matrix(links)
  if (is_mat && (ncol(links) != 2 || nrow(links) == 0)) bad()
  v <- if (is.factor(links)) as.character(links) else as.vector(links)
  if (!(is.character(v) || is.numeric(v)) || anyNA(v)) bad()

  p   <- length(nms)
  idx <- if (is.character(v)) match(v, nms)
         else ifelse(v == round(v) & v >= 1 & v <= p, as.integer(v), NA_integer_)
  if (anyNA(idx))
    stop('`links` names features not in the loadings: ', name_list(unique(v[is.na(idx)])), '.', call. = FALSE)

  if (is_mat) {
    pairs <- matrix(idx, ncol = 2)
  }
  else {
    u <- unique(idx)
    if (length(u) < 2)
      stop('`links` needs at least two features: every pair among them is linked.', call. = FALSE)
    pairs <- t(utils::combn(u, 2))
  }

  if (any(pairs[, 1] == pairs[, 2]))
    stop('`links` pairs a feature with itself.', call. = FALSE)

  pairs
}


# Rays from the origin to `cols`, coloured by `col`, labelled where `label`.
ord_plot_rays <- function (cols, col, label) {
  keep <- plot_arrow_keep(0, 0, cols[, 1], cols[, 2])
  if (!any(keep)) return (invisible())
  graphics::arrows(0, 0, cols[keep, 1], cols[keep, 2], length = 0.06, col = col[keep])
  lab <- keep & label
  if (any(lab))
    plot_label(cols[lab, 1] * 1.08, cols[lab, 2] * 1.08, rownames(cols)[lab], col = col[lab], cex = 0.8)
  invisible()
}


# Dashed links between the tips of `cols`, each labelled at its midpoint
# with its length on the plane.
ord_plot_links <- function (cols, lk) {
  x0 <- cols[lk$from, 1]; y0 <- cols[lk$from, 2]
  x1 <- cols[lk$to,   1]; y1 <- cols[lk$to,   2]
  graphics::segments(x0, y0, x1, y1, lty = 2, col = PLOT_HIGHLIGHT)
  plot_label(
    (x0 + x1) / 2, (y0 + y1) / 2, vapply(lk$length, plot_num, '', digits = 2L),
    col = PLOT_HIGHLIGHT, cex = 0.75 )
  invisible()
}
