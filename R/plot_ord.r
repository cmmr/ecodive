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
#'        value). By default, the groups that [lda()], [plsda()] and
#'        [oplsda()] store in their `groups` attribute, if any.
#'
#' @param data  A data frame of sample metadata, with one row per sample
#'        (matched by row name, else in order), when `groups` names one of
#'        its columns.
#'
#' @param axes  Which two axes to plot, by number, as in [envfit()]. Default
#'        `c(1, 2)`. A result with only one axis (such as [lda()] of two
#'        groups) is drawn as a strip chart, with one strip per group.
#'
#' @param loadings  How many of the features in the `loadings` attribute to
#'        draw, those furthest from the origin on the plotted axes first: a
#'        number, `Inf` (or `TRUE`) for all, or `0` (or `FALSE`) for none.
#'        Default `10`. Drawn as arrows for [pca()], [rpca()],
#'        [rpca_phylo()], [rda()], [plsda()] and [oplsda()], and as labelled
#'        points for the species scores of [ca()] and [cca()].
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
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim`, `asp` or `cex`. They replace the defaults
#'        made from `x`. To leave out the axes, use `xaxt = "n"` and `yaxt =
#'        "n"`, since `axes` chooses which axes to plot.
#'
#' @details
#' The plot is made of what `x` holds:
#'
#' * **Title**: the method's name, as [list_dist_ords()] and
#'   [list_count_ords()] give it.
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
#' where they are.
#'
#' Base graphics draw on the current device, so further points, text or
#' arrows can be added on top. Use [envfit()]'s plot method, `plot(fit, add =
#' TRUE)`, to add fitted variables.
#'
#' @return `x`, invisibly.
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
    col = NULL, pch = NULL, ... ) {

  args     <- list(...)
  att      <- attributes(x)
  pos      <- plot_legend_pos(legend)
  loadings <- plot_count(loadings, 'loadings')
  for (arg in c('biplot', 'centroids', 'target'))
    if (!isTRUE(get(arg)) && !isFALSE(get(arg)))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)

  if (ncol(x) > 1) axes <- plot_axes(axes, ncol(x))
  else if (ncol(x) == 0) stop('`x` has no axes to plot.', call. = FALSE)
  xy   <- unclass(x)
  attributes(xy) <- att[c('dim', 'dimnames')]

  g    <- ord_plot_groups(x, groups, data)
  look <- ord_plot_look(g, col, pch)

  main <- ord_method_name(att$method)
  sub  <- ord_subtitle(x)
  labs <- ord_axis_labels(x)


  #--------------------------------------------------
  # One axis: a strip chart, one strip per group.
  #--------------------------------------------------
  if (ncol(x) == 1) {
    v <- xy[, 1]
    s <- if (is.null(g$f)) list(v) else split(v, g$f)
    plot_do(
      graphics::stripchart, c(list(s), args),
      list(
        method = 'overplot', vertical = TRUE, pch = look$pch,
        col = if (is.null(g$f)) look$col else look$pal,
        xlim = c(0.5, length(s) + 0.5), main = main, sub = sub, ylab = labs[[1]],
        group.names = if (is.null(g$f)) '' else levels(g$f) ))
    return (invisible(x))
  }


  #--------------------------------------------------
  # Two axes: the samples, then the overlays.
  #--------------------------------------------------
  pts <- xy[, axes, drop = FALSE]

  # ca() and cca() hold species scores, drawn where they are; the others,
  # loadings, drawn as arrows.
  sp_pts   <- isTRUE(att$method %in% c('ca', 'cca'))
  load     <- if (loadings > 0) ord_plot_loadings(att$loadings, axes, loadings)
  arr_con  <- if (biplot) ord_plot_biplot(att$biplot, axes, if (centroids) att$centroids)
  cen      <- if (centroids && !is.null(att$centroids) && nrow(att$centroids))
    att$centroids[, axes, drop = FALSE]
  tgt      <- if (target && !is.null(att$target))
    att$target[, axes, drop = FALSE]

  # The plot holds the samples, centroids, species and Procrustes targets.
  all_xy <- rbind(pts, cen, tgt, if (sp_pts) load)

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
    if (sp_pts) {
      graphics::points(load[, 1], load[, 2], pch = 4, col = PLOT_MUTED)
      plot_label(load[, 1], load[, 2], rownames(load), col = PLOT_MUTED, pos = 3, cex = 0.8)
    }
    else {
      plot_arrows(load * plot_arrow_mul(load, ext), col = PLOT_MUTED)
    }
  }

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

  n      <- nrow(x)
  labels <- rownames(x)
  title  <- NULL

  if (is.null(groups)) {
    v <- attr(x, 'groups', exact = TRUE)
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

  f <- if (is.factor(v)) droplevels(v) else factor(v)
  if (anyNA(f)) f <- addNA(f, ifany = TRUE)

  list(f = f, values = NULL, numeric = FALSE, title = title)
}



# Colours and symbols of the points, and of the legend. Returns list(col,
# pch (the points'), pal (one colour per level), legend, legend_col,
# legend_pch).
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
    rng  <- range(v, na.rm = TRUE)
    if (!all(is.finite(rng)))
      return (list(col = PLOT_MUTED, pch = pch, legend = NULL))
    step <- if (rng[2] > rng[1]) (v - rng[1]) / (rng[2] - rng[1]) else v * 0 + 0.5
    pcol <- ramp[1 + round(99 * step)]
    pcol[is.na(pcol)] <- PLOT_MUTED
    at   <- pretty(rng, 4)
    at   <- at[at >= rng[1] & at <= rng[2]]
    if (length(at) < 2) at <- rng
    lcol <- ramp[1 + round(99 * (if (rng[2] > rng[1]) (at - rng[1]) / (rng[2] - rng[1]) else 0.5))]
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

  lpch <- if (length(pch) == k) pch else pch[[1]]
  list(
    col = pal[as.integer(g$f)], pch = if (length(pch) == k) pch[as.integer(g$f)] else pch,
    pal = pal, legend = ifelse(is.na(lev), 'NA', lev), legend_col = pal, legend_pch = lpch )
}



# The `n` longest loadings on the plotted axes, as a two-column matrix, or
# NULL.
ord_plot_loadings <- function (load, axes, n) {
  if (is.null(load) || !is.matrix(load) || nrow(load) == 0 || ncol(load) < max(axes))
    return (NULL)
  a <- load[, axes, drop = FALSE]
  if (is.null(rownames(a))) rownames(a) <- as.character(seq_len(nrow(a)))
  len <- rowSums(a ^ 2)
  a[utils::head(order(len, decreasing = TRUE), n), , drop = FALSE]
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
