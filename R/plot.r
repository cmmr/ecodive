# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Helpers shared by ecodive's plot() methods (plot_ord.r, plot_tests.r and
# plot_tables.r). The methods only draw what a result already holds: their
# titles, labels and legends come from its columns and attributes, and they
# compute no new statistics. They use only the 'graphics' and 'grDevices'
# packages, and only functions that R 3.6.0 already had.



# Calls `fun` with `args` (a list, usually list(...) from the user), and with
# each of `defaults` that `args` doesn't name. So the user's main =, xlab =,
# col =, ... replace ecodive's defaults, as with any base plot().
plot_do <- function (fun, args, defaults) {
  given <- names(args)
  if (is.null(given)) given <- character(0)
  do.call(fun, c(args, defaults[setdiff(names(defaults), given)]))
}


# `args` without the elements named `drop`.
plot_drop <- function (args, drop) {
  args[setdiff(names(args), drop)]
}


# The user's value for argument `nm` from `args` (list(...)), else `default`.
plot_arg <- function (args, nm, default) {
  if (nm %in% names(args)) args[[nm]] else default
}


# A number for a title or label: `digits` significant digits, as
# print.ecodive_test() shows them.
plot_num <- function (v, digits = 3L) {
  if (length(v) == 0 || is.na(v[[1]])) return ('NA')
  trimws(formatC(v[[1]], digits = digits, format = 'g'))
}


# A p-value for a title or label, as print.ecodive_test() shows it.
plot_p <- function (p, digits = 3L) {
  if (length(p) == 0 || is.na(p[[1]])) return ('NA')
  trimws(format.pval(p[[1]], digits = digits))
}


# "Spearman" from "spearman".
plot_cap <- function (s) {
  paste0(toupper(substring(s, 1L, 1L)), substring(s, 2L))
}


# The colours of `k` groups: hcl.colors()'s "Dark 3" palette (R >= 3.6.0).
plot_palette <- function (k) {
  grDevices::hcl.colors(max(k, 1L), 'Dark 3')[seq_len(k)]
}


# Colour of the marks that stand out (observed statistics, rejected
# features, constraint arrows), and of those that don't.
PLOT_HIGHLIGHT <- 'firebrick'
PLOT_MUTED     <- 'grey50'


# Checks `legend`, a position keyword for graphics::legend(), or NULL, FALSE
# or NA for none. Returns the keyword, or NULL.
plot_legend_pos <- function (legend) {
  if (is.null(legend) || isFALSE(legend) || identical(legend, NA)) return (NULL)
  pos <- c('bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top',
           'topright', 'right', 'center')
  if (!is.character(legend) || length(legend) != 1 || !legend %in% pos)
    stop('`legend` must be one of ', paste0('"', pos, '"', collapse = ', '),
         ', or NULL for none.', call. = FALSE)
  legend
}


# Checks a count of things to draw or label (`loadings`, `labels`, `n_taxa`,
# ...): a number >= 0, Inf, TRUE (all) or FALSE (none). Returns a number.
plot_count <- function (v, arg) {
  if (isTRUE(v))  return (Inf)
  if (isFALSE(v)) return (0)
  if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0)
    stop('`', arg, '` must be a number of at least 0, Inf, TRUE or FALSE.', call. = FALSE)
  floor(v)
}


# Checks `axes`, which two of the `n` columns to plot. Returns them as
# integers.
plot_axes <- function (axes, n) {
  ok <- is.numeric(axes) && length(axes) == 2 && all(is.finite(axes)) &&
    all(axes == floor(axes)) && all(axes >= 1) && all(axes <= n) && axes[[1]] != axes[[2]]
  if (!ok)
    stop(sprintf('`axes` must be two different axis numbers from 1 to %d.', n), call. = FALSE)
  as.integer(axes)
}


# The factor that scales arrows `a` (a matrix with two columns: their tips)
# so that the longest reaches `frac` of the extent `ext` (two half-widths,
# measured from the origin) on each axis. As biplots do, only the arrows'
# directions and relative lengths are kept.
plot_arrow_mul <- function (a, ext, frac = 0.8) {
  top <- c(max(abs(a[, 1])), max(abs(a[, 2])))
  m   <- ext / top
  m   <- m[is.finite(m) & m > 0]
  if (length(m)) frac * min(m) else 1
}


# Draws arrows from the origin to the rows of `a` (already scaled), with
# their row names at their tips. Zero-length arrows are left out.
plot_arrows <- function (a, col, cex = 0.8) {
  keep <- rowSums(a ^ 2) > 0
  a    <- a[keep, , drop = FALSE]
  if (nrow(a) == 0) return (invisible())
  graphics::arrows(0, 0, a[, 1], a[, 2], length = 0.08, col = col)
  plot_label(a[, 1] * 1.08, a[, 2] * 1.08, rownames(a), col = col, cex = cex)
  invisible()
}


# Text with a thin halo in the background colour, so that labels stay
# readable over points of their own colour.
plot_label <- function (x, y, labels, col, pos = NULL, cex = 0.85, font = 1) {
  if (length(x) == 0) return (invisible())
  bg <- graphics::par('bg')
  if (identical(bg, 'transparent')) bg <- 'white'
  dx <- graphics::strwidth('x', cex = cex) * 0.15
  dy <- graphics::strheight('x', cex = cex) * 0.15
  for (a in seq(0, 2 * pi, length.out = 13)[-13])
    graphics::text(
      x + dx * cos(a), y + dy * sin(a), labels, col = bg, pos = pos,
      cex = cex, font = font, xpd = NA )
  graphics::text(x, y, labels, col = col, pos = pos, cex = cex, font = font, xpd = NA)
  invisible()
}


# A legend for colours `col` with names `labels`. `title` may be NULL.
plot_legend <- function (pos, labels, col, pch, title = NULL, lty = NULL) {
  if (is.null(pos) || length(labels) == 0) return (invisible())
  graphics::legend(
    pos, legend = labels, col = col, pch = pch, lty = lty, title = title,
    bg = 'white', cex = 0.8, inset = 0.01 )
  invisible()
}
