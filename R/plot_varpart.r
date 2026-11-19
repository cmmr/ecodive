# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot a variation partitioning
#'
#' A Venn diagram of [varpart()]'s result: a circle (or, for four sets, an
#' ellipse) for each set of variables, with the adjusted \eqn{R^2} of each
#' fraction in its region, and that of the residuals below.
#'
#' @param x  The result of [varpart()].
#'
#' @param digits  The number of decimal places of the fractions. Default
#'        `2`.
#'
#' @param col  The sets' colours, one per set. Defaults to
#'        [grDevices::hcl.colors()]'s `"Dark 3"` palette.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub` or `cex`.
#'
#' @details
#' The regions are all the same size, whatever their fractions: their values
#' are written in them, not drawn to scale. Fractions shared by several sets
#' can be negative, when the sets explain more together than the sum of what
#' they explain apart. The fractions are the `adj_r_squared` column of `x`,
#' and its row names (`[a]`, `[b]`, ...) follow the regions as in
#' `vegan::showvarparts()`: `[a]` is the first set's own fraction, `[b]` the
#' second's, and so on.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     vp <- varpart(
#'       iris[, 3:4], species = ~ Species, sepal = ~ Sepal.Length + Sepal.Width,
#'       data = iris, norm = 'hellinger' )
#'     plot(vp)
#'
#'     vp <- varpart(
#'       iris[, 3:4], species = ~ Species, length = ~ Sepal.Length,
#'       width = ~ Sepal.Width, data = iris, norm = 'hellinger' )
#'     plot(vp)
plot.ecodive_varpart <- function (x, digits = 2L, col = NULL, ...) {

  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0)
    stop('`digits` must be a number of at least 0.', call. = FALSE)

  labs <- attr(x, 'set_labels', exact = TRUE)
  k    <- length(labs)
  lay  <- varpart_layout(k)
  if (is.null(lay))
    stop('`x` must partition 2, 3 or 4 sets.', call. = FALSE)

  if (is.null(col)) col <- plot_palette(k)
  col <- rep_len(col, k)

  method <- attr(x, 'method', exact = TRUE)
  how    <- switch(
    method,
    rda   = sprintf('RDA%s', if (!identical(attr(x, 'norm'), 'none')) sprintf(' (%s)', attr(x, 'norm')) else ''),
    dbrda = 'db-RDA',
    cca   = 'CCA' )

  plot_do(
    graphics::plot.default, c(list(NA, NA), list(...)),
    list(
      type = 'n', asp = 1, xlim = lay$xlim, ylim = lay$ylim, axes = FALSE,
      xlab = '', ylab = '',
      main = sprintf('Variation partitioning by %s', how),
      sub  = as.expression(bquote(
        'Adjusted' ~ R^2 ~ 'of each fraction,' ~ .(attr(x, 'n', exact = TRUE)) ~ 'samples' )) ))

  # The sets' shapes, then their names.
  tt <- seq(0, 2 * pi, length.out = 181)
  for (i in seq_len(k)) {
    e <- lay$shapes[i, ]
    u <- e[['a']] * cos(tt)
    v <- e[['b']] * sin(tt)
    r <- e[['angle']] * pi / 180
    graphics::polygon(
      e[['x']] + u * cos(r) - v * sin(r), e[['y']] + u * sin(r) + v * cos(r),
      col = grDevices::adjustcolor(col[[i]], alpha.f = 0.15), border = col[[i]], lwd = 2 )
  }
  graphics::text(lay$names[, 1], lay$names[, 2], names(labs), adj = lay$names_adj,
                 col = col, font = 2, xpd = NA)

  # The fractions, each in its region.
  fmt <- function (v) {
    if (is.na(v)) return ('NA')
    v <- round(v, digits)
    formatC(if (v == 0) 0 else v, format = 'f', digits = digits)  # not "-0.00"
  }
  regions <- varpart_regions(k)
  for (i in seq_along(regions)) {
    at <- lay$labels[[paste(regions[[i]], collapse = ',')]]
    graphics::text(at[[1]], at[[2]], fmt(x$adj_r_squared[[i]]), cex = if (k == 4) 0.8 else 0.9)
  }

  res <- x$adj_r_squared[[nrow(x)]]
  graphics::text(
    lay$xlim[[2]], lay$ylim[[1]], paste('Residuals =', fmt(res)),
    adj = c(1, 0), cex = 0.9, xpd = NA )

  invisible(x)
}



# Where to draw the sets of a varpart() Venn diagram of k sets: their shapes
# (centre, semi-axes and angle of each circle or ellipse), the label point of
# each region (named by its sets, as "1,3"), and where to write the sets'
# names. Circles for two and three sets; for four, symmetric ellipses placed
# so that each of the 15 regions has room for its label.
varpart_layout <- function (k) {

  shapes <- function (...) {
    m <- rbind(...)
    colnames(m) <- c('x', 'y', 'a', 'b', 'angle')
    m
  }

  if (k == 2) return (list(
    shapes = shapes(c(-0.55, 0, 1, 1, 0), c(0.55, 0, 1, 1, 0)),
    labels = list(`1` = c(-1, 0), `2` = c(1, 0), `1,2` = c(0, 0)),
    names  = rbind(c(-0.55, 1.12), c(0.55, 1.12)),
    names_adj = c(0.5, 0),
    xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.35) ))

  if (k == 3) return (list(
    shapes = shapes(c(-0.55, 0.32, 1, 1, 0), c(0.55, 0.32, 1, 1, 0), c(0, -0.63, 1, 1, 0)),
    labels = list(
      `1`   = c(-0.96, 0.56), `2`   = c(0.96, 0.56),  `3`   = c(0, -1.1),
      `1,2` = c(0, 0.7),      `2,3` = c(0.6, -0.35),  `1,3` = c(-0.6, -0.35),
      `1,2,3` = c(0, 0) ),
    names  = rbind(c(-0.55, 1.44), c(0.55, 1.44), c(0, -1.75)),
    names_adj = c(0.5, 0.5),
    xlim = c(-1.6, 1.6), ylim = c(-1.95, 1.55) ))

  if (k == 4) return (list(
    shapes = shapes(
      c(-0.28, 0.03, 1, 0.41,  46), c(0, -0.2, 1, 0.41,  46),
      c( 0,   -0.2, 1, 0.41, -46), c(0.28, 0.03, 1, 0.41, -46) ),
    labels = list(
      `1`     = c(-0.79, -0.20), `2`     = c(-0.30, -0.77),
      `3`     = c( 0.30, -0.77), `4`     = c( 0.79, -0.20),
      `1,2`   = c(-0.55, -0.52), `2,3`   = c( 0,    -0.56),
      `1,3`   = c(-0.55,  0.12), `1,4`   = c( 0,     0.59),
      `2,4`   = c( 0.55,  0.12), `3,4`   = c( 0.55, -0.52),
      `1,2,4` = c( 0.29,  0.37), `1,2,3` = c(-0.33, -0.23),
      `2,3,4` = c( 0.33, -0.23), `1,3,4` = c(-0.29,  0.37),
      `1,2,3,4` = c(0, 0.06) ),
    names  = rbind(c(-1.05, 0.12), c(-0.65, -1.12), c(0.65, -1.12), c(1.05, 0.12)),
    names_adj = c(0.5, 0.5),
    xlim = c(-1.45, 1.45), ylim = c(-1.3, 0.95) ))

  NULL
}
