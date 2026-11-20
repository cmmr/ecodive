# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Plot an association network
#'
#' Draws a network from [assoc_network()], the two networks compared by
#' [assoc_network_compare()] side by side, or the significant differential
#' associations from [assoc_network_diff()], on a circular layout.
#'
#' @param x  A result of [assoc_network()], [assoc_network_compare()] or
#'        [assoc_network_diff()].
#'
#' @param layout  Where to put the nodes: `NULL` (the default) for the
#'        circular layout described below, or a matrix with two columns (x
#'        and y) and one row per feature, with the features' names as row
#'        names.
#'
#' @param isolated  Whether to draw the features that have no edges.
#'        Default `TRUE`.
#'
#' @param curved  Whether to draw edges as curves bowed toward the centre,
#'        rather than straight lines. Default `FALSE`.
#'
#' @param labels  Whether to label the nodes with the features' names.
#'        Default `TRUE`.
#'
#' @param legend  Where to put the legend of edge colours: a keyword such as
#'        `"bottomright"` (the default), or `NULL` for none.
#'
#' @param ...  Passed on to [graphics::plot()], such as `main` or `sub`.
#'
#' @details
#' **Layout.** The nodes go clockwise from 12 o'clock, grouped by cluster
#' (the largest clusters first) and, within a cluster, by decreasing
#' degree; features without edges come last. For [assoc_network_compare()]
#' both networks are drawn on `x`'s layout; for [assoc_network_diff()] the
#' features are in `x`'s input order.
#'
#' **Styling.** Nodes are coloured by cluster (grey for cluster 0) and sized by the first of the centralities that defined the hubs
#' (`hubs` in [assoc_network()]); hubs have a black border and a bold
#' label. Edges are red for positive and blue for negative associations
#' (from near the ends of [grDevices::hcl.colors()]'s "Blue-Red 3"
#' palette), with
#' widths proportional to \eqn{|r|}. For [assoc_network_diff()], an edge is
#' red when the association is larger in `x` and blue when it is larger in
#' `y`, and its width is proportional to \eqn{|z|}.
#'
#' The title is the result's (see [ecodive_display]); the subtitle gives
#' the cutoff and the number of edges. The two panels of a comparison are
#' titled with the `groups` it was given (or "x" and "y"), under the
#' comparison's title.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     basis <- matrix(rnorm(60 * 20), 60, 20)
#'     basis[, 2:5] <- basis[, 2:5] + basis[, 1]
#'     counts <- matrix(rpois(1200, exp(basis + 4)), 60, 20,
#'                      dimnames = list(NULL, paste0('OTU', 1:20)))
#'     net <- assoc_network(sparcc(counts, n_draws = 10))
#'     plot(net)
#'     plot(net, curved = TRUE, isolated = FALSE)
#'
plot.ecodive_network <- function (
    x, layout = NULL, isolated = TRUE, curved = FALSE, labels = TRUE,
    legend = 'bottomright', ... ) {

  args  <- list(...)
  nodes <- x$nodes
  edges <- x$edges
  legend <- plot_legend_pos(legend)

  pos <- network_plot_layout(
    layout, nodes$feature, nodes$cluster, nodes$degree,
    c(edges$from, edges$to), isolated )

  size_by <- attr(x, 'hubs', exact = TRUE)[1]
  if (is.null(size_by) || is.na(size_by)) size_by <- 'degree'
  sub <- sprintf(
    '|r| >= %s, %d edge%s', format(attr(x, 'threshold', exact = TRUE)),
    nrow(edges), if (nrow(edges) == 1) '' else 's' )

  network_plot_draw(
    pos, nodes, edges$from, edges$to,
    edge_col = network_sign_col(edges$sign),
    edge_lwd = 0.5 + 3 * abs(edges$association),
    size     = nodes[[size_by]],
    curved   = curved, labels = labels, args = args,
    main     = plot_title(x, 'Association Network'), sub = sub )

  if (!is.null(legend) && nrow(edges) > 0)
    plot_legend(
      legend, c('positive', 'negative'), network_sign_col(c(1L, -1L)),
      pch = NA, lty = 1 )

  invisible(x)
}



#' @rdname plot.ecodive_network
#' @export
plot.ecodive_network_compare <- function (
    x, layout = NULL, isolated = TRUE, curved = FALSE, labels = TRUE,
    legend = 'bottomright', ... ) {

  args   <- list(...)
  legend <- plot_legend_pos(legend)
  cl     <- x$clusters
  deg    <- x$nodes[x$nodes$centrality == 'degree', , drop = FALSE]
  deg_x  <- deg$x[match(cl$feature, deg$feature)]
  deg_y  <- deg$y[match(cl$feature, deg$feature)]

  pos <- network_plot_layout(
    layout, cl$feature, cl$cluster_x, deg_x,
    c(x$edges_x$from, x$edges_x$to, x$edges_y$from, x$edges_y$to), isolated )

  grp <- attr(x, 'groups', exact = TRUE)
  if (is.null(grp)) grp <- c('x', 'y')

  title <- plot_arg(args, 'main', plot_title(x, 'Comparison of Association Networks'))
  args  <- plot_drop(args, 'main')
  old <- graphics::par(mfrow = c(1L, 2L), oma = c(0, 0, 2.5, 0))
  on.exit(graphics::par(old))

  sides <- list(
    list(e = x$edges_x, cl = cl$cluster_x, hub = cl$hub_x, deg = deg_x),
    list(e = x$edges_y, cl = cl$cluster_y, hub = cl$hub_y, deg = deg_y) )
  for (k in 1:2) {
    s     <- sides[[k]]
    nodes <- data.frame(feature = cl$feature, cluster = s$cl, hub = s$hub, stringsAsFactors = FALSE)
    network_plot_draw(
      pos, nodes, s$e$from, s$e$to,
      edge_col = network_sign_col(s$e$sign),
      edge_lwd = 0.5 + 3 * abs(s$e$association),
      size     = s$deg,
      curved   = curved, labels = labels, args = args,
      main     = grp[[k]],
      sub      = sprintf('%d edge%s', nrow(s$e), if (nrow(s$e) == 1) '' else 's') )
    if (k == 2L && !is.null(legend))
      plot_legend(legend, c('positive', 'negative'), network_sign_col(c(1L, -1L)), pch = NA, lty = 1)
  }
  graphics::mtext(title, side = 3, outer = TRUE, line = 0.5, font = 2, cex = 1.2)

  invisible(x)
}



#' @rdname plot.ecodive_network
#' @export
plot.ecodive_network_diff <- function (
    x, layout = NULL, isolated = TRUE, curved = FALSE, labels = TRUE,
    legend = 'bottomright', ... ) {

  args   <- list(...)
  legend <- plot_legend_pos(legend)
  df     <- as.data.frame(unclass(x), stringsAsFactors = FALSE)

  feats <- unique(c(df$feature1, df$feature2))
  e     <- df[df$reject %in% TRUE, , drop = FALSE]

  nodes <- data.frame(feature = feats, cluster = 0L, hub = FALSE, stringsAsFactors = FALSE)
  pos   <- network_plot_layout(layout, feats, NULL, NULL, c(e$feature1, e$feature2), isolated)

  larger <- ifelse(e$association_x > e$association_y, 1L, -1L)
  sub    <- sprintf(
    '%d pair%s with %s-adjusted p <= %s', nrow(e), if (nrow(e) == 1) '' else 's',
    attr(x, 'p_adjust', exact = TRUE), format(attr(x, 'alpha', exact = TRUE)) )

  network_plot_draw(
    pos, nodes, e$feature1, e$feature2,
    edge_col = network_sign_col(larger),
    edge_lwd = 0.5 + 0.5 * pmin(abs(e$z), 8),
    size     = NULL,
    curved   = curved, labels = labels, args = args,
    main     = plot_title(x, 'Differential Associations'), sub = sub )

  if (!is.null(legend) && nrow(e) > 0)
    plot_legend(legend, c('larger in x', 'larger in y'), network_sign_col(c(1L, -1L)), pch = NA, lty = 1)

  invisible(x)
}



# Red for positive, blue for negative: from near the ends of "Blue-Red 3"
# (its very ends are almost black at the widths of thin edges).
network_sign_col <- function (sign) {
  pal <- grDevices::hcl.colors(11L, 'Blue-Red 3')
  ifelse(sign > 0, pal[[10]], pal[[2]])
}



# Node positions: a matrix (rows named by feature) for the features drawn.
# The circular layout orders clusters by size (largest first; ties by
# cluster id), nodes in a cluster by decreasing degree, then the features
# without edges (in input order).
network_plot_layout <- function (layout, features, cluster, degree, ends, isolated) {

  if (!isTRUE(isolated) && !isFALSE(isolated))
    stop('`isolated` must be TRUE or FALSE.', call. = FALSE)

  linked <- features %in% ends
  keep   <- if (isolated) rep(TRUE, length(features)) else linked

  if (!is.null(layout)) {
    if (!is.matrix(layout) || !is.numeric(layout) || ncol(layout) != 2L ||
        is.null(rownames(layout)) || !all(features[keep] %in% rownames(layout)))
      stop(
        '`layout` must be a numeric matrix with two columns and the ',
        'features\' names as row names.', call. = FALSE )
    pos <- layout[features[keep], , drop = FALSE]
    return (pos)
  }

  n <- length(features)
  if (is.null(cluster)) cluster <- integer(n)
  if (is.null(degree))  degree  <- numeric(n)
  cl   <- ifelse(linked, cluster, NA_integer_)
  size <- table(cl)
  csz  <- ifelse(is.na(cl), 0L, as.integer(size[as.character(cl)]))
  ord  <- order(!linked, -csz, cl, -degree, seq_len(n))
  ord  <- ord[keep[ord]]

  k   <- length(ord)
  ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / max(k, 1L)
  pos <- cbind(x = cos(ang), y = sin(ang))
  rownames(pos) <- features[ord]
  pos
}



# Draws nodes and edges. `nodes` has feature, cluster, hub; edges are given
# by their ends' names.
network_plot_draw <- function (
    pos, nodes, from, to, edge_col, edge_lwd, size, curved, labels, args, main, sub ) {

  plot_do(graphics::plot, plot_drop(args, c('col')), list(
    x = 0, y = 0, type = 'n', asp = 1, axes = FALSE, xlab = '', ylab = '',
    xlim = range(c(pos[, 1], -1.15, 1.15)), ylim = range(c(pos[, 2], -1.15, 1.15)),
    main = main, sub = sub ))

  drawn <- from %in% rownames(pos) & to %in% rownames(pos)
  from  <- from[drawn]; to <- to[drawn]
  edge_col <- edge_col[drawn]; edge_lwd <- edge_lwd[drawn]

  # Thinner edges first, so strong ones are on top.
  for (e in order(edge_lwd)) {
    p0 <- pos[from[[e]], ]
    p1 <- pos[to[[e]], ]
    if (isTRUE(curved)) {
      mid <- (p0 + p1) / 4
      graphics::xspline(
        c(p0[1], mid[1], p1[1]), c(p0[2], mid[2], p1[2]),
        shape = c(0, 1, 0), open = TRUE, border = edge_col[[e]], lwd = edge_lwd[[e]] )
    }
    else {
      graphics::segments(p0[1], p0[2], p1[1], p1[2], col = edge_col[[e]], lwd = edge_lwd[[e]])
    }
  }

  i   <- match(rownames(pos), nodes$feature)
  cl  <- nodes$cluster[i]
  k   <- max(c(0L, cl), na.rm = TRUE)
  pal <- plot_palette(k)
  bg  <- ifelse(cl > 0L, pal[pmax(cl, 1L)], 'grey80')
  hub <- nodes$hub[i] %in% TRUE

  cex <- if (is.null(size)) rep(1.6, length(i)) else {
    v <- size[i]
    r <- max(v, na.rm = TRUE)
    if (!is.finite(r) || r <= 0) rep(1.2, length(i)) else 1 + 1.6 * v / r
  }

  graphics::points(
    pos[, 1], pos[, 2], pch = 21, bg = bg, cex = cex,
    col = ifelse(hub, 'black', 'grey40'), lwd = ifelse(hub, 2.5, 1) )

  # Labels just outside each node, away from the centre of the layout.
  if (isTRUE(labels) && nrow(pos) > 0) {
    ctr <- colMeans(pos)
    ang <- atan2(pos[, 2] - ctr[2], pos[, 1] - ctr[1])
    off <- graphics::strheight('M', cex = 0.75) * (0.6 + cex / 2.5)
    for (v in seq_len(nrow(pos)))
      graphics::text(
        pos[v, 1] + off[v] * cos(ang[v]), pos[v, 2] + off[v] * sin(ang[v]),
        rownames(pos)[v], cex = 0.75, font = if (hub[v]) 2 else 1, xpd = NA,
        adj = c((1 - cos(ang[v])) / 2, (1 - sin(ang[v])) / 2) )
  }
  invisible()
}
