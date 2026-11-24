# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# The class of log_ratio() and balance() results, `ecodive_log_ratio`: a
# named numeric vector with `method`, `numerator`, `denominator` and
# `pseudocount` attributes. Subsetting, sort() and rev() already return a
# plain vector (base R drops the class); arithmetic does too, here, since
# `-x` or `x - mean(x)` is no longer log(N / D).



#' @export
print.ecodive_log_ratio <- function (x, ...) {
  cat(log_ratio_header(x), sep = '\n')
  print(log_ratio_plain(x), ...)
  invisible(x)
}


#' @export
Ops.ecodive_log_ratio <- function (e1, e2) {
  e1 <- log_ratio_plain(e1)
  if (missing(e2)) return (get(.Generic)(e1))
  get(.Generic)(e1, log_ratio_plain(e2))
}


#' @export
Math.ecodive_log_ratio <- function (x, ...) {
  get(.Generic)(log_ratio_plain(x), ...)
}


# A log-ratio result as a plain numeric vector, keeping its names; anything
# else unchanged.
log_ratio_plain <- function (x) {
  if (!inherits(x, 'ecodive_log_ratio')) return (x)
  nms <- names(x)
  attributes(x) <- NULL
  names(x) <- nms
  x
}


# print()'s header: "Log ratio of 2 over 3 features; pseudocount 0", the
# numerator, the denominator, and the samples with no value.
log_ratio_header <- function (x) {
  a   <- function (nm) attr(x, nm, exact = TRUE)
  num <- a('numerator')
  den <- a('denominator')
  pc  <- a('pseudocount')
  na  <- sum(is.na(x))
  c(
    sprintf(
      '%s of %d over %d %s%s', log_ratio_what(x), length(num), length(den),
      if (length(den) == 1L) 'feature' else 'features',
      if (is.numeric(pc) && length(pc) == 1L) sprintf('; pseudocount %s', plot_num(pc, 4L)) else '' ),
    sprintf('Numerator:   %s', name_list(num)),
    sprintf('Denominator: %s', name_list(den)),
    if (na > 0) sprintf(
      '%d %s NA (%s)', na, if (na == 1L) 'sample is' else 'samples are',
      if (identical(a('method'), 'balance')) 'a zero in a selected feature'
      else 'a zero sum' ) )
}


# "Log ratio" or "Balance", from the `method` attribute.
log_ratio_what <- function (x) {
  if (identical(attr(x, 'method', exact = TRUE), 'balance')) 'Balance' else 'Log ratio'
}



#' Plot Log Ratios
#'
#' Plots the per-sample values of [log_ratio()] or [balance()], as Qurro's
#' sample plot does (Fedarko et al. 2020): by group, as a box plot with the
#' samples drawn over it; against a numeric variable, as a scatter plot; or
#' with no `groups`, as one box of all the samples.
#'
#' @param x  The result of [log_ratio()] or [balance()].
#'
#' @param groups  The samples' groups or a numeric variable: the name of a
#'        column of `data`, or a vector with one value per sample (matched to
#'        the samples by name when both are named, as for [plot.ecodive_ord()]).
#'        A factor, character or logical vector gives a box per value (in the
#'        order of a factor's levels); a numeric vector puts the variable on
#'        the x axis. Default `NULL`: one box of all the samples.
#'
#' @param data  A data frame of sample metadata, with one row per sample
#'        (matched by row name, else in order), when `groups` names one of
#'        its columns.
#'
#' @param legend  For groups: where to put the legend of the samples'
#'        colours, a keyword for [graphics::legend()], such as `"topright"`
#'        (the default), or `NULL` for none.
#'
#' @param col  The points' colours: with groups, one per group (in order,
#'        or named by group), as for [plot.ecodive_ord()]; else for all the
#'        points.
#'
#' @param pch  The points' symbols. Default `16` (filled circles) with
#'        groups, else `1`.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `ylab`, `ylim` or `cex`.
#'
#' @details
#' The title is "Log ratio" or "Balance", from the result's `method`
#' attribute; the y axis label names the numerator's and the
#' denominator's features (up to three on each side, else their number),
#' as \eqn{\log(N / D)} for a log ratio; and the subtitle lists them (up to
#' five on each side).
#'
#' Samples whose value is undefined (`NA`: a numerator or denominator that
#' sums to 0 in [log_ratio()], a zero in a selected feature in [balance()]
#' with `pseudocount = 0`) are left out, as Qurro leaves them out, and
#' counted in the subtitle.
#'
#' Within a box, the samples are spread sideways so that ties stay visible:
#' the \eqn{i}th smallest of a group sits at \eqn{0.6\,((0.618034\,i) \bmod
#' 1 - 0.5)} from the box's centre. The spread is fixed, so it doesn't use
#' (or change) R's random number generator, and the same values are always
#' drawn in the same places.
#'
#' Nothing is fitted or tested: test the log ratio as any other per-sample
#' variable (see [log_ratio()]).
#'
#' @return `x`, invisibly.
#'
#' @references
#' Fedarko, M. W., Martino, C., Morton, J. T., et al. (2020). Visualizing
#' 'omic feature rankings and log-ratios using Qurro. *NAR Genomics and
#' Bioinformatics*, 2, lqaa023. \doi{10.1093/nargab/lqaa023}
#'
#' @seealso [log_ratio()], [balance()], and [log_ratio_features()] to choose
#'   their features from a ranking.
#'
#' @export
#' @examples
#'     lr <- log_ratio(iris[, 1:4], 'Petal.Width', c('Sepal.Length', 'Sepal.Width'))
#'     lr
#'
#'     # By group, against a numeric variable, and all together.
#'     plot(lr, 'Species', iris)
#'     plot(lr, 'Sepal.Length', iris)
#'     plot(lr)
#'
#'     # A balance.
#'     plot(balance(iris[, 1:4], 3:4, 1:2), iris$Species)
plot.ecodive_log_ratio <- function (
    x, groups = NULL, data = NULL, legend = 'topright', col = NULL, pch = NULL, ... ) {

  args <- list(...)
  pos  <- plot_legend_pos(legend)
  v    <- log_ratio_plain(x)
  nms  <- names(v)
  v    <- unname(as.double(v))
  n    <- length(v)
  if (n == 0) stop('`x` has no samples to plot.', call. = FALSE)

  g    <- plot_groups(n, nms, groups, data)
  look <- ord_plot_look(g, col, pch)
  ok   <- is.finite(v)
  main <- log_ratio_what(x)
  ylab <- log_ratio_ylab(x)
  sub  <- log_ratio_plot_sub(x)

  if (!any(ok)) stop('`x` has no defined values to plot.', call. = FALSE)


  # A numeric variable: a scatter plot.
  if (isTRUE(g$numeric)) {
    pc <- if (is.null(col)) graphics::par('fg') else col
    pp <- if (is.null(pch)) 1 else pch
    plot_do(
      graphics::plot.default, c(list(g$values[ok], v[ok]), args),
      list(
        col = if (length(pc) == n) pc[ok] else pc,
        pch = if (length(pp) == n) pp[ok] else pp, las = 1,
        main = main, sub = sub, ylab = ylab,
        xlab = if (is.null(g$title)) '' else g$title ))
    return (invisible(x))
  }


  # Groups (or none): a box per group, with the samples over it.
  if (is.null(g$f)) {
    box  <- rep(1L, n)
    labs <- 'All samples'
    pcol <- if (is.null(col)) PLOT_MUTED else col
    pcol <- rep_len(pcol, n)
    ppch <- rep_len(if (is.null(pch)) 1 else pch, n)
  }
  else {
    used <- which(tabulate(as.integer(g$f)[ok], nlevels(g$f)) > 0)
    box  <- match(as.integer(g$f), used)
    labs <- ifelse(is.na(levels(g$f)), 'NA', levels(g$f))[used]
    pcol <- rep_len(look$col, n)
    ppch <- rep_len(look$pch, n)
  }

  xy <- log_ratio_spread(box[ok], v[ok])

  plot_do(
    graphics::plot.default, c(list(NA_real_, NA_real_), args),
    list(
      xlim = c(0.5, length(labs) + 0.5), ylim = range(v[ok]), xaxt = 'n',
      las = 1, main = main, sub = sub, ylab = ylab,
      xlab = if (is.null(g$title)) '' else g$title ))

  graphics::boxplot(
    split(v[ok], factor(box[ok], levels = seq_along(labs))), add = TRUE,
    at = seq_along(labs), outline = FALSE, border = 'grey60', col = NA,
    names = labs, axes = FALSE, boxwex = 0.6 )

  if (plot_x_axis(args))
    graphics::axis(1, at = seq_along(labs), labels = labs)

  graphics::points(xy$x, xy$y, col = pcol[ok], pch = ppch[ok])

  if (!is.null(g$f) && !is.null(pos)) {
    lpch <- if (length(pch) == nlevels(g$f)) pch[used] else ppch[[1]]
    plot_legend(pos, labs, look$pal[used], lpch, g$title)
  }

  invisible(x)
}


# The points of a box plot: at x = `box` (a group number per sample), each
# spread sideways by its rank `i` within its group (by value; ties in
# sample order), 0.6 * ((0.618034 i) mod 1 - 0.5). Returns list(x, y).
log_ratio_spread <- function (box, v) {
  off <- numeric(length(v))
  for (b in unique(box)) {
    i <- which(box == b)
    i <- i[order(v[i])]
    off[i] <- 0.6 * ((seq_along(i) * 0.618034) %% 1 - 0.5)
  }
  list(x = box + off, y = v)
}


# The y axis label: "log((a + b) / c)" for a log ratio, "Balance of a + b
# vs c" for a balance; more than three features on a side are counted.
log_ratio_ylab <- function (x) {
  side <- function (f, paren) {
    if (length(f) > 3L) return (sprintf('%d features', length(f)))
    s <- paste(f, collapse = ' + ')
    if (paren && length(f) > 1L) sprintf('(%s)', s) else s
  }
  num <- attr(x, 'numerator',   exact = TRUE)
  den <- attr(x, 'denominator', exact = TRUE)
  if (identical(attr(x, 'method', exact = TRUE), 'balance'))
    sprintf('Balance of %s vs %s', side(num, FALSE), side(den, FALSE))
  else
    sprintf('log(%s / %s)', side(num, TRUE), side(den, TRUE))
}


# The subtitle: "Numerator: 'a', 'b'; denominator: 'c'" (up to five on
# each side), then the number of samples with no value.
log_ratio_plot_sub <- function (x) {
  na <- sum(!is.finite(log_ratio_plain(x)))
  s  <- sprintf(
    'Numerator: %s; denominator: %s',
    name_list(attr(x, 'numerator', exact = TRUE)),
    name_list(attr(x, 'denominator', exact = TRUE)) )
  if (na > 0)
    s <- sprintf(
      '%s; %d %s with an undefined %s not shown', s, na,
      if (na == 1L) 'sample' else 'samples', tolower(log_ratio_what(x)) )
  s
}
