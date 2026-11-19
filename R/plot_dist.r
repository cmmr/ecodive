# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot a distance matrix
#'
#' A heatmap of a distance matrix from [beta_div()] or any of its metric
#' functions (such as [bray()] or [unweighted_unifrac()]), or of the
#' proportionality of features from [vlr()], [rho_p()] or [phi_s()]: one row
#' and one column per sample (or feature), each cell coloured by the
#' distance between its row's and its column's, with the samples ordered as
#' the matrix's `clust` attribute says.
#'
#' @param x  A distance matrix of class `ecodive_dist`, as ecodive's
#'        distance functions return. (`plot(as.dist(m))` of another matrix
#'        would draw base R's plot of the distances as numbers; give it the
#'        class with `class(d) <- c("ecodive_dist", "dist")`.)
#'
#' @param groups  Sorts the samples by group: the name of a column of
#'        `data`, or a vector with one value per sample (matched by name when
#'        both are named, as for [permanova()]). Factor, character or logical
#'        groups are drawn as blocks, in the order of a factor's levels,
#'        separated by white lines and labelled with the group; numeric ones
#'        sort the samples from smallest to largest.
#'
#' @param data  A data frame of sample metadata, with one row per sample, when
#'        `groups` names one of its columns.
#'
#' @param col  The colour scale, from the smallest value to the largest.
#'        Defaults to [grDevices::hcl.colors()]'s `"YlOrRd"` palette, light
#'        for near and dark for far; or, when some values are negative (as
#'        from [rho_p()]), its `"Blue-Red 3"` palette, symmetric about 0.
#'
#' @param key  Whether to draw the colour key on the right. Default `TRUE`.
#'
#' @param ...  Passed on to [graphics::image()], such as `main`, `sub`, or
#'        `useRaster`.
#'
#' @details
#' **Order.** The samples are ordered by the `clust` argument of the
#' function that made `x`, recorded in its `clust` attribute:
#'
#' * Not set (`clust = NULL`, the default): the plot clusters them itself,
#'   by average linkage (UPGMA, [stats::hclust()]), so that similar samples
#'   sit together. `x` is not changed.
#' * `"none"` (`clust = FALSE`): their order in `x`, which is the input's.
#' * A clustering method: `x` is already in that order, which is kept.
#'
#' To plot a matrix in its own order, set `attr(x, "clust") <- "none"`.
#'
#' **Title.** The metric's display name, from the `metric_name` attribute
#' (such as `"Bray-Curtis Dissimilarity"`), or else its `metric` or `method`
#' attribute; `"Distances"` without any of them.
#' With `groups`, the samples are sorted by group, keeping that order within
#' each group.
#'
#' The rows run from the first sample at the top to the last at the bottom,
#' and the columns from left to right in the same order, so the matrix is
#' symmetric about the diagonal from the top left; only the rows are
#' labelled. The diagonal, which a `dist` object doesn't store, is left
#' blank.
#'
#' The heatmap has one cell per pair of samples, so it needs memory for an
#' \eqn{n \times n}{n x n} matrix: subset or summarize thousands of samples
#' first. Large heatmaps are drawn as raster images, where the device can.
#'
#' The class `ecodive_dist` marks distances that ecodive computed. It comes
#' before `"dist"`, so every function for `dist` objects works on them
#' as before.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     dm <- bray(iris[, 1:4])
#'     plot(dm)
#'
#'     # Sorted by species: three blocks along the diagonal.
#'     plot(dm, 'Species', iris)
#'
#'     # In the input's order, not clustered.
#'     plot(bray(iris[, 1:4], clust = FALSE))
#'
#'     # Clustered by Ward's method when computed: the result is reordered.
#'     dm <- bray(iris[, 1:4], clust = 'ward.D2')
#'     head(attr(dm, 'order'))
#'     plot(dm, 'Species', iris)
#'
#'     # Proportionality, with negative values.
#'     plot(rho_p(ex_counts + 1))
plot.ecodive_dist <- function (x, groups = NULL, data = NULL, col = NULL, key = TRUE, ...) {

  args <- list(...)
  if (!isTRUE(key) && !isFALSE(key))
    stop('`key` must be TRUE or FALSE.', call. = FALSE)

  n <- attr(x, 'Size', exact = TRUE)
  if (length(n) != 1 || is.na(n) || n < 2)
    stop('`x` must have at least 2 samples to plot.', call. = FALSE)

  labels <- attr(x, 'Labels', exact = TRUE)
  names_ <- if (is.null(labels)) as.character(seq_len(n)) else as.character(labels)


  #--------------------------------------------------
  # The order of the samples: as `clust` says (clustered here when it isn't
  # set), then by group.
  #--------------------------------------------------
  cl   <- attr(x, 'clust', exact = TRUE)
  base <- if (is.null(cl) && !anyNA(x)) dist_clust_order(x, 'average') else seq_len(n)
  how  <- dist_plot_order_name(if (is.null(cl) && !anyNA(x)) 'average' else cl)

  g      <- plot_groups(n, labels, groups, data)
  ord    <- base
  blocks <- NULL
  if (!is.null(g$f)) {
    ord    <- base[order(as.integer(g$f)[base])]  # stable: keeps `base` within groups
    runs   <- rle(as.integer(g$f[ord]))
    lev    <- levels(g$f)
    blocks <- list(
      ends   = cumsum(runs$lengths),
      labels = ifelse(is.na(lev), 'NA', lev)[runs$values] )
    blocks$mids <- blocks$ends - (runs$lengths - 1) / 2
  }
  else if (!is.null(g$values)) {
    ord <- base[order(g$values[base])]
  }

  m <- as.matrix(x)[ord, ord, drop = FALSE]
  diag(m) <- NA
  names_ <- names_[ord]


  #--------------------------------------------------
  # The colour scale.
  #--------------------------------------------------
  v   <- as.vector(x)
  v   <- v[is.finite(v)]
  neg <- length(v) > 0 && any(v < 0)
  top <- if (length(v)) max(abs(v)) else 0
  if (top == 0) top <- 1
  lim <- if (neg) c(-top, top) else c(0, top)
  pal <- if (!is.null(col)) col
         else if (neg) grDevices::hcl.colors(65, 'Blue-Red 3')
         else grDevices::hcl.colors(64, 'YlOrRd', rev = TRUE)
  brk <- seq(lim[[1]], lim[[2]], length.out = length(pal) + 1L)


  #--------------------------------------------------
  # Margins: room for the row labels, and the key.
  #--------------------------------------------------
  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  ylab <- if (is.null(blocks)) names_ else blocks$labels
  mai  <- graphics::par('mai')
  mai[2L] <- max(mai[2L], max(graphics::strwidth(ylab, 'inch', cex = 0.7), 0, na.rm = TRUE) + 0.25)
  if (key) mai[4L] <- max(mai[4L], 1)
  graphics::par(mai = mai)

  # image() draws its first row at the bottom: rows reversed, to put the
  # first sample at the top.
  z    <- t(m[rev(seq_len(n)), , drop = FALSE])
  rast <- n > 100 && identical(grDevices::dev.capabilities('rasterImage')$rasterImage, 'yes')
  what <- dist_plot_name(x)

  plot_do(
    graphics::image, c(list(seq_len(n), seq_len(n), z), args),
    list(
      col = pal, breaks = brk, axes = FALSE, xlab = '', ylab = '',
      useRaster = rast, main = what$title,
      sub = paste0(
        n, ' ', what$unit,
        if (!is.null(g$title)) paste0(', by ', g$title)
        else if (!is.null(groups)) ', by group',
        if (!is.null(how)) paste0(', ', how) ) ))

  if (is.null(blocks)) {
    graphics::axis(2, at = n - seq_len(n) + 1, labels = names_, las = 1, tick = FALSE, cex.axis = 0.7)
  }
  else {
    graphics::axis(2, at = n - blocks$mids + 1, labels = blocks$labels, las = 1, tick = FALSE, cex.axis = 0.7)
    graphics::axis(1, at = blocks$mids, labels = blocks$labels, tick = FALSE, cex.axis = 0.7)
    cuts <- utils::head(blocks$ends, -1) + 0.5
    if (length(cuts))
      graphics::abline(v = cuts, h = n + 1 - cuts, col = 'white', lwd = 2)
  }
  graphics::box()

  if (key) dist_plot_key(pal, brk, what$key)

  invisible(x)
}



# How the samples were ordered, for the subtitle, from a `clust` value: NULL
# for their own order.
dist_plot_order_name <- function (clust) {
  if (is.null(clust) || identical(clust, 'none')) return (NULL)
  if (identical(clust, 'pcoa')) return ('ordered by PCoA axis 1')
  sprintf('clustered (%s linkage)', clust)
}


# The title, the unit (samples or features) and the colour key's title of a
# distance matrix: the title from its `metric_name` attribute (the metric's
# display name), the rest from `metric` (or `method`, as results from before
# ecodive 2.4.0 have it). Without either, "Distances".
dist_plot_name <- function (x) {
  m <- attr(x, 'metric', exact = TRUE)
  if (is.null(m)) m <- attr(x, 'method', exact = TRUE)
  m <- if (is.character(m) && length(m) == 1 && !is.na(m)) m else ''
  title <- attr(x, 'metric_name', exact = TRUE)
  if (!is.character(title) || length(title) != 1 || is.na(title))
    title <- if (nzchar(m)) m else 'Distances'
  features <- m %in% c('vlr', 'rho_p', 'phi_s')
  list(
    title = title,
    unit  = if (features) 'features' else 'samples',
    key   = if (features) m else 'Distance' )
}



# A vertical colour key in the right margin: one box per colour of `pal`
# between the `brk`s, labelled with round values.
dist_plot_key <- function (pal, brk, title) {

  usr <- graphics::par('usr')
  w   <- diff(usr[1:2])
  x0  <- usr[[2]] + 0.04 * w
  x1  <- usr[[2]] + 0.08 * w
  y   <- usr[[3]] + (brk - brk[[1]]) / diff(range(brk)) * diff(usr[3:4])
  k   <- length(pal)

  graphics::rect(x0, y[-(k + 1)], x1, y[-1], col = pal, border = NA, xpd = NA)
  graphics::rect(x0, y[[1]], x1, y[[k + 1]], xpd = NA)

  at <- pretty(range(brk), 5)
  at <- at[at >= min(brk) & at <= max(brk)]
  ya <- usr[[3]] + (at - brk[[1]]) / diff(range(brk)) * diff(usr[3:4])
  graphics::segments(x1, ya, x1 + 0.01 * w, ya, xpd = NA)
  graphics::text(x1 + 0.015 * w, ya, format(at, trim = TRUE), adj = c(0, 0.5), cex = 0.7, xpd = NA)
  graphics::text((x0 + x1) / 2, usr[[4]], title, pos = 3, cex = 0.75, xpd = NA)

  invisible()
}
