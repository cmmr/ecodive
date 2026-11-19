# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot a permutation test
#'
#' Plots the result of one of ecodive's permutation tests. By default,
#' [permanova()], [anosim()], [mantel()] and [mantel_partial()] results show
#' the permutation distribution of their statistic, with the observed value
#' marked; [mrpp()] shows each group's mean within-group distance;
#' [permdisp()] each group's dispersion; and [procrustes_test()] the
#' Procrustes fit. `type = "null"` shows the permutation distribution for any
#' of them. Nothing is recomputed: the plot is made from the stored results.
#'
#' @param x  A permutation test result, of class `ecodive_test`: from
#'        [permanova()], [anosim()], [mrpp()], [mantel()],
#'        [mantel_partial()], [permdisp()], [procrustes_test()], or
#'        [dist_test()].
#'
#' @param type  What to draw: `"null"` for the permutation distribution of
#'        the statistic; or for [mrpp()] and [permdisp()], `"groups"` (their
#'        default) for a bar for each group; or for [procrustes_test()],
#'        `"fit"` (its default) for the Procrustes rotation, drawn by
#'        [plot.ecodive_ord()].
#'
#' @param ...  Passed on to the plotting function: [graphics::hist()] for
#'        `"null"`, [graphics::barplot()] for `"groups"`, and
#'        [plot.ecodive_ord()] for `"fit"`. Such as `main`, `sub`, `xlab`,
#'        `col` or `breaks`; they replace the defaults made from `x`.
#'
#' @details
#' The title is the test's name, and the subtitle its statistic, p-value
#' and standardized effect size \eqn{z}, as printed.
#'
#' * `"null"`: a histogram of the statistic in each permutation (the
#'   `f_perm`, `r_perm` or `delta_perm` attribute), with a line at the
#'   observed statistic. The p-value is the proportion of permutations, with
#'   the observed one, at least as extreme as the observed statistic. It
#'   needs a test run with `permutations` above 0.
#' * `"groups"` for [mrpp()]: the mean within-group distance of each group
#'   (`class_delta`), with lines at the observed \eqn{\delta} (their
#'   weighted mean) and the \eqn{\delta} expected by chance
#'   (`expected_delta`, the mean of all distances).
#' * `"groups"` for [permdisp()]: each group's mean distance to its centroid
#'   or median (`dispersions`).
#' * `"fit"` for [procrustes_test()]: the symmetric Procrustes rotation it
#'   tested (its `procrustes` field), with an arrow from each sample to the
#'   target ordination.
#'
#' @return `x`, invisibly.
#'
#' @name plot.ecodive_test
#' @export
#' @examples
#'     dm <- bray(iris[, 1:4])
#'
#'     # The pseudo-F of each permutation, and the observed F.
#'     res <- permanova(dm, 'Species', iris, permutations = 199)
#'     plot(res)
#'
#'     # MRPP: each species' mean within-group distance.
#'     res <- mrpp(dm, 'Species', iris, permutations = 199)
#'     plot(res)
#'     plot(res, type = 'null')
#'
#'     # PERMDISP: the dispersion of each species.
#'     plot(permdisp(dm, 'Species', iris, permutations = 199))
#'
#'     # PROTEST: the Procrustes fit of two ordinations.
#'     res <- procrustes_test(pcoa(dm, n_axes = 2), pca(iris[, 1:4], n_axes = 2))
#'     plot(res)
plot.ecodive_test <- function (x, type = 'null', ...) {
  test_plot_type(type, 'null')
  test_plot_null(x, ...)
}


#' @rdname plot.ecodive_test
#' @export
plot.ecodive_mrpp <- function (x, type = c('groups', 'null'), ...) {

  type <- test_plot_type(type, c('groups', 'null'))
  if (type == 'null') return (test_plot_null(x, ...))

  d    <- x$class_delta
  exp_ <- x$expected_delta
  obs  <- x$statistic[[1]]
  top  <- max(d, exp_, obs, na.rm = TRUE)

  plot_do(
    graphics::barplot, c(list(d), list(...)),
    list(
      main = attr(x, 'title', exact = TRUE), sub = test_subtitle(x),
      ylab = 'Mean within-group distance', col = 'grey85',
      ylim = c(0, 1.5 * top), las = 1 ))

  graphics::abline(h = obs, col = PLOT_HIGHLIGHT, lwd = 2)
  graphics::abline(h = exp_, lty = 2, lwd = 2)
  plot_legend(
    'topright', c('Observed delta', 'Expected delta'),
    col = c(PLOT_HIGHLIGHT, 'black'), pch = NULL, lty = c(1, 2) )

  invisible(x)
}


#' @rdname plot.ecodive_test
#' @export
plot.ecodive_permdisp <- function (x, type = c('groups', 'null'), ...) {

  type <- test_plot_type(type, c('groups', 'null'))
  if (type == 'null') return (test_plot_null(x, ...))

  d <- x$dispersions

  plot_do(
    graphics::barplot, c(list(d), list(...)),
    list(
      main = attr(x, 'title', exact = TRUE), sub = test_subtitle(x),
      ylab = sprintf('Mean distance to group %s', if (identical(x$type, 'median')) 'median' else 'centroid'),
      col = 'grey85', ylim = c(0, 1.1 * max(d, na.rm = TRUE)), las = 1 ))

  invisible(x)
}


#' @rdname plot.ecodive_test
#' @export
plot.ecodive_procrustes_test <- function (x, type = c('fit', 'null'), ...) {

  type <- test_plot_type(type, c('fit', 'null'))
  if (type == 'null') return (test_plot_null(x, ...))

  plot_do(
    plot.ecodive_ord, c(list(x$procrustes), list(...)),
    list(main = attr(x, 'title', exact = TRUE), sub = test_subtitle(x)) )

  invisible(x)
}



# Checks `type` against `choices`; returns it.
test_plot_type <- function (type, choices) {
  if (identical(type, choices)) return (choices[[1]])
  if (!is.character(type) || length(type) != 1 || !type %in% choices)
    stop('`type` must be ', paste0('"', choices, '"', collapse = ' or '), '.', call. = FALSE)
  type
}


# "F = 533, p = 0.005, z = 551", as print.ecodive_test() gives it.
test_subtitle <- function (x) {
  s <- sprintf('%s = %s', names(x$statistic), plot_num(x$statistic))
  if (isTRUE(x$n_perm > 0))
    s <- sprintf('%s, p = %s, z = %s', s, plot_p(x$p_value), plot_num(x$z))
  s
}


# A histogram of the permuted statistics, with a line at the observed one.
test_plot_null <- function (x, ...) {

  stat <- x$statistic
  nm   <- names(stat)
  perm <- attr(x, paste0(tolower(nm), '_perm'), exact = TRUE)
  perm <- perm[is.finite(perm)]

  if (length(perm) == 0)
    stop(
      'There are no permutations to plot: run the test with `permutations` ',
      'above 0.', call. = FALSE )

  obs  <- stat[[1]]
  args <- list(...)

  # The x axis spans the histogram's breaks, and the observed statistic.
  brk  <- graphics::hist(perm, breaks = plot_arg(args, 'breaks', 'Sturges'), plot = FALSE)$breaks

  plot_do(
    graphics::hist, c(list(perm), args),
    list(
      main = attr(x, 'title', exact = TRUE), sub = test_subtitle(x),
      xlab = sprintf('%s (%d permutations, %s)', nm, x$n_perm, x$permutation),
      xlim = range(brk, if (is.finite(obs)) obs),
      col = 'grey85', border = 'grey45', las = 1 ))

  if (is.finite(obs)) {
    graphics::abline(v = obs, col = PLOT_HIGHLIGHT, lwd = 2)
    graphics::mtext('observed', side = 3, at = obs, line = 0.1, col = PLOT_HIGHLIGHT, cex = 0.8)
  }

  invisible(x)
}
