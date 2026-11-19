# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot fitted variables
#'
#' Draws the variables that [envfit()] fitted onto an ordination: an arrow
#' for each numeric variable (or feature), pointing the way it increases,
#' with a length proportional to \eqn{\sqrt{r^2}}, and a label at the
#' centroid of each factor level. With `add = TRUE`, they go on top of an
#' ordination plot, such as `plot(ord)`.
#'
#' @param x  The result of [envfit()].
#'
#' @param add  Whether to add to the current plot, such as `plot(ord)` of the
#'        ordination that `x` was fitted to (`TRUE`), or start a new one
#'        (`FALSE`, the default).
#'
#' @param axes  Which two of `x`'s axes to draw, by number among them.
#'        Default `c(1, 2)`.
#'
#' @param p_max  Draw only the variables with a `p_value` at most `p_max`.
#'        Default `1`: all of them (including those with no p-value).
#'
#' @param col  The colour of the arrows and labels.
#'
#' @param ...  With `add = FALSE`, passed on to
#'        [graphics::plot.default()], such as `main`, `xlab`, `xlim` or
#'        `asp`.
#'
#' @details
#' The arrows are scaled to fit the plot: with `add = TRUE`, the longest
#' reaches 80% of the way from the origin to the nearest edge of the
#' plot; with `add = FALSE`, 80% of the way to the furthest centroid (or
#' they keep their own lengths, \eqn{\sqrt{r^2}} at most 1, when there are
#' no centroids). Only their directions and relative lengths have meaning,
#' as in any biplot. The centroids are in the ordination's units, so they are
#' drawn where they are. Factor levels are labelled with the level, or with
#' the row name (the variable and the level pasted together) when two
#' variables share a level's name.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     ord <- pcoa(bray(iris[, 1:4]))
#'     fit <- envfit(ord, iris, permutations = 99)
#'
#'     plot(ord, 'Species', iris)
#'     plot(fit, add = TRUE)
#'
#'     # On their own.
#'     plot(fit)
plot.ecodive_envfit <- function (
    x, add = FALSE, axes = c(1, 2), p_max = 1, col = 'blue3', ... ) {

  if (!isTRUE(add) && !isFALSE(add))
    stop('`add` must be TRUE or FALSE.', call. = FALSE)
  if (!is.numeric(p_max) || length(p_max) != 1 || is.na(p_max))
    stop('`p_max` must be a number.', call. = FALSE)

  cols <- setdiff(names(x), c('variable', 'level', 'r2', 'p_value', 'z'))
  if (length(cols) < 2)
    stop('`x` has only one axis: there is nothing to draw in two dimensions.', call. = FALSE)
  cols <- cols[plot_axes(axes, length(cols))]

  p    <- x$p_value
  keep <- if (p_max >= 1) rep(TRUE, nrow(x)) else !is.na(p) & p <= p_max
  df   <- x[keep, , drop = FALSE]

  is_vec <- is.na(df$level)

  arr <- as.matrix(df[is_vec, cols, drop = FALSE]) * sqrt(df$r2[is_vec])
  rownames(arr) <- rownames(df)[is_vec]
  attributes(arr) <- list(dim = dim(arr), dimnames = list(rownames(df)[is_vec], cols))

  cen <- as.matrix(df[!is_vec, cols, drop = FALSE])
  lev <- df$level[!is_vec]
  attributes(cen) <- list(dim = dim(cen), dimnames = list(
    if (anyDuplicated(x$level[!is.na(x$level)])) rownames(df)[!is_vec] else lev, cols ))

  if (!add) {
    # Arrows to 80% of the furthest centroid, or their own lengths.
    ext <- rep(if (nrow(cen)) max(abs(cen)) else 1.25, 2)
    mul <- if (nrow(arr) && nrow(cen)) plot_arrow_mul(arr, ext) else 1
    lim <- range(0, arr * mul * 1.15, cen, finite = TRUE)
    plot_do(
      graphics::plot.default, c(list(NA, NA), list(...)),
      list(
        type = 'n', asp = 1, xlim = lim, ylim = lim,
        xlab = cols[[1]], ylab = cols[[2]], main = 'Fitted variables',
        sub = sprintf('%d permutations', attr(x, 'n_perm', exact = TRUE)) ))
    graphics::abline(h = 0, v = 0, lty = 3, col = 'grey70')
  }
  else {
    # Arrows to 80% of the way to the plot's nearest edge.
    usr <- graphics::par('usr')
    ext <- c(
      if (usr[1] < 0 && usr[2] > 0) min(-usr[1], usr[2]) else diff(usr[1:2]) / 2,
      if (usr[3] < 0 && usr[4] > 0) min(-usr[3], usr[4]) else diff(usr[3:4]) / 2 )
    mul <- if (nrow(arr)) plot_arrow_mul(arr, ext) else 1
  }

  if (nrow(arr))
    plot_arrows(arr * mul, col = col)

  if (nrow(cen)) {
    graphics::points(cen[, 1], cen[, 2], pch = 3, col = col, cex = 1.2, lwd = 2)
    plot_label(cen[, 1], cen[, 2], rownames(cen), col = col, pos = 3, font = 2)
  }

  invisible(x)
}



#' Plot a SIMPER result
#'
#' A dot chart of the taxa that contribute most to the Bray-Curtis
#' dissimilarity between each pair of groups, from [simper()]: each taxon's
#' average contribution, largest first, in one panel per pair. Filled points
#' mark the taxa whose adjusted permutation p-value is at most `alpha`.
#'
#' @param x  The result of [simper()].
#'
#' @param n_taxa  How many taxa to show for each pair, those with the largest
#'        average contributions. Default `10`; `Inf` for all.
#'
#' @param pairs  Which pairs of groups to show, by name (such as
#'        `"setosa-virginica"`). Default: all of them.
#'
#' @param alpha  Taxa with a `p_adjusted` at most `alpha` are drawn with
#'        filled points. Default `0.05`.
#'
#' @param ...  Passed on to [graphics::dotchart()], such as `main`, `xlab`,
#'        `cex` or `color`; and `sub`, to [graphics::title()].
#'
#' @details
#' Each pair's label gives its mean between-group Bray-Curtis dissimilarity
#' (the `overall` attribute), which the contributions of all its taxa add up
#' to. Without permutations (`permutations = 0`), all points are open.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     res <- simper(round(iris[, 1:4] * 10), 'Species', iris, permutations = 99)
#'     plot(res)
#'     plot(res, pairs = 'setosa-virginica')
plot.ecodive_simper <- function (x, n_taxa = 10, pairs = NULL, alpha = 0.05, ...) {

  n_taxa <- plot_count(n_taxa, 'n_taxa')
  args   <- list(...)

  all_pairs <- unique(x$pair)
  if (is.null(pairs)) pairs <- all_pairs
  if (!is.character(pairs) || length(pairs) == 0 || !all(pairs %in% all_pairs))
    stop('`pairs` must be among ', name_list(all_pairs), '.', call. = FALSE)
  pairs <- all_pairs[all_pairs %in% pairs]

  overall <- attr(x, 'overall', exact = TRUE)

  # dotchart() draws each group's values from the bottom up and the groups
  # from the top down; so reverse each pair's rows to put its largest first.
  rows <- lapply(pairs, function (p) rev(utils::head(which(x$pair == p), n_taxa)))
  idx  <- unlist(rows)
  grp  <- factor(rep(pairs, lengths(rows)), levels = pairs)
  levels(grp) <- sprintf('%s (%s)', pairs, vapply(pairs, function (p) plot_num(overall[[p]]), ''))

  sig <- !is.na(x$p_adjusted[idx]) & x$p_adjusted[idx] <= alpha
  sub <- if (any(!is.na(x$p_adjusted)))
    sprintf('Filled: adjusted p <= %s', format(alpha))

  plot_do(
    graphics::dotchart, c(list(x$average[idx]), plot_drop(args, 'sub')),
    list(
      labels = x$taxon[idx], groups = grp, pch = ifelse(sig, 19, 1),
      xlim = c(0, max(x$average[idx], na.rm = TRUE)),
      main = 'SIMPER', xlab = 'Average contribution to Bray-Curtis dissimilarity' ))

  graphics::title(sub = plot_arg(args, 'sub', sub))

  invisible(x)
}



#' Plot a Dirichlet-multinomial t-test
#'
#' A volcano plot of [dirmult_ttest()]'s results: each feature's log2 fold
#' change against its adjusted p-value, with the features that differ (those
#' with `reject = TRUE`) highlighted, and the most significant of them named.
#'
#' @param x  The result of [dirmult_ttest()].
#'
#' @param labels  How many of the differentially abundant features to name,
#'        those with the smallest adjusted p-values. Default `10`; `0` for
#'        none.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `xlab`, `xlim` or `cex`.
#'
#' @details
#' The y axis is \eqn{-\log_{10}} of `p_adjusted` (an adjusted p-value of 0
#' is drawn at the smallest positive double). The dashed line marks `alpha`;
#' a feature also needs its confidence interval for the log2 fold change to
#' exclude 0 to count as differentially abundant, so some features above the
#' line may not be highlighted.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     res <- dirmult_ttest(round(iris[, 1:4] * 10), 'Species', iris)
#'     plot(res)
plot.ecodive_dirmult_ttest <- function (x, labels = 10, ...) {

  labels <- plot_count(labels, 'labels')

  ct    <- attr(x, 'contrast', exact = TRUE)
  n     <- attr(x, 'n', exact = TRUE)
  alpha <- attr(x, 'alpha', exact = TRUE)
  ref   <- if (is.na(ct[[2]])) 'all others' else ct[[2]]

  y   <- -log10(pmax(x$p_adjusted, .Machine$double.xmin))
  rej <- x$reject %in% TRUE

  plot_do(
    graphics::plot.default, c(list(x$log2_fc, y), list(...)),
    list(
      col  = ifelse(rej, PLOT_HIGHLIGHT, PLOT_MUTED), pch = 16, las = 1,
      ylim = range(0, y, finite = TRUE),
      main = 'Dirichlet-multinomial t-test',
      sub  = sprintf('%d of %d features differentially abundant', sum(rej), nrow(x)),
      xlab = sprintf('log2 fold change, %s (n = %d) / %s (n = %d)', ct[[1]], n[[1]], ref, n[[2]]),
      ylab = expression(-log[10] ~ (adjusted ~ italic(p))) ))

  graphics::abline(v = 0, lty = 3, col = 'grey70')
  if (length(alpha) && is.finite(alpha) && alpha > 0)
    graphics::abline(h = -log10(alpha), lty = 2, col = 'grey40')

  if (labels > 0 && any(rej)) {
    i <- which(rej)
    i <- utils::head(i[order(x$p_adjusted[i])], labels)
    graphics::text(x$log2_fc[i], y[i], x$feature[i], pos = 3, cex = 0.75, xpd = NA)
  }

  invisible(x)
}



#' Plot an ANCOM result
#'
#' A dot chart of [ancom()]'s \eqn{W} statistic for the features with the
#' largest, with the features declared differentially abundant
#' (`reject = TRUE`) highlighted, and a line at the cutoff that separates
#' them.
#'
#' @param x  The result of [ancom()].
#'
#' @param n_features  How many features to show, those with the largest
#'        \eqn{W}. Default `30`; `Inf` for all.
#'
#' @param ...  Passed on to [graphics::dotchart()], such as `main`, `xlab`,
#'        or `cex`; and `sub`, to [graphics::title()].
#'
#' @details
#' A feature's \eqn{W} is the number of other features whose log ratio with
#' it differs significantly between the groups, out of \eqn{D - 1} for
#' \eqn{D} features. ANCOM rejects the features with \eqn{W} above
#' `cutoff` \eqn{\times D} (the dashed line). There is no line when the
#' largest \eqn{W} was below `theta`, so that nothing was rejected.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     counts <- matrix(rpois(200, 20), nrow = 10)
#'     counts[1:5, 1:3] <- counts[1:5, 1:3] * 8
#'     colnames(counts) <- paste0('F', 1:20)
#'     res <- ancom(counts, rep(c('a', 'b'), each = 5))
#'     plot(res)
plot.ecodive_ancom <- function (x, n_features = 30, ...) {

  n_features <- plot_count(n_features, 'n_features')
  args       <- list(...)

  D   <- nrow(x)
  cut <- attr(x, 'cutoff', exact = TRUE)
  ord <- utils::head(order(x$W, decreasing = TRUE), n_features)
  idx <- rev(ord)  # dotchart() draws from the bottom up.
  rej <- x$reject[idx] %in% TRUE

  sub <- sprintf(
    '%d of %d features differentially abundant%s', sum(x$reject %in% TRUE), D,
    if (length(ord) < D) sprintf('; top %d shown', length(ord)) else '' )

  plot_do(
    graphics::dotchart, c(list(x$W[idx]), plot_drop(args, 'sub')),
    list(
      labels = x$feature[idx], pch = 19,
      color  = ifelse(rej, PLOT_HIGHLIGHT, PLOT_MUTED),
      xlim   = c(0, max(D - 1, 1)),
      main   = sprintf('ANCOM: %s', paste(attr(x, 'groups', exact = TRUE), collapse = ' vs ')),
      xlab   = sprintf('W (of %d log ratios)', D - 1) ))

  graphics::title(sub = plot_arg(args, 'sub', sub))
  if (length(cut) && !is.na(cut))
    graphics::abline(v = cut * D, lty = 2, col = 'grey40')

  invisible(x)
}



#' Plot a BIO-ENV result
#'
#' The correlation of the best subset of variables of each size with the
#' distances, from [bioenv()], with the best subset overall circled and named
#' in the subtitle.
#'
#' @param x  The result of [bioenv()].
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `ylab`, `ylim` or `col`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     res <- bioenv(bray(iris[, 1:4]), iris[, 1:4])
#'     plot(res)
plot.ecodive_bioenv <- function (x, ...) {

  best <- attr(x, 'best', exact = TRUE)
  i    <- match(length(best), x$size)

  plot_do(
    graphics::plot.default, c(list(x$size, x$correlation), list(...)),
    list(
      type = 'b', pch = 16, xaxt = 'n', las = 1,
      main = if (isTRUE(attr(x, 'partial', exact = TRUE))) 'BIO-ENV (partial)' else 'BIO-ENV',
      sub  = paste('Best:', paste(best, collapse = ', ')),
      xlab = 'Number of variables',
      ylab = sprintf('%s correlation with distances', plot_cap(attr(x, 'method', exact = TRUE))) ))

  graphics::axis(1, at = x$size)
  if (!is.na(i))
    graphics::points(x$size[[i]], x$correlation[[i]], cex = 2.5, lwd = 2, col = PLOT_HIGHLIGHT)

  invisible(x)
}



#' Plot a forward selection
#'
#' The cumulative adjusted \eqn{R^2} as [forward_select()] added each
#' variable, with a dashed line at that of the global model (all candidates
#' together), and the next best variable, which was not added, as an open
#' point.
#'
#' @param x  The result of [forward_select()].
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `ylab`, `ylim` or `col`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     fs <- forward_select(
#'       iris[, 3:4], ~ Species + Sepal.Length + Sepal.Width, iris,
#'       norm = 'hellinger', permutations = 99 )
#'     plot(fs)
plot.ecodive_forward_select <- function (x, ...) {

  g    <- attr(x, 'global', exact = TRUE)
  nb   <- attr(x, 'next_best', exact = TRUE)
  k    <- nrow(x)
  step <- 0:k
  y    <- c(0, x$adj_r_squared_cum)
  labs <- c('(none)', x$variable)

  nb_y <- if (!is.null(nb) && nrow(nb)) nb$adj_r_squared_cum[[1]]
  glob <- if (!is.null(g)) g$adj_r_squared[[1]]

  plot_do(
    graphics::plot.default, c(list(step, y), list(...)),
    list(
      type = 'b', pch = 16, xaxt = 'n', las = 1,
      xlim = c(-0.3, k + if (!is.null(nb_y)) 1.3 else 0.3),
      ylim = range(y, nb_y, glob, finite = TRUE) + c(-0.06, 0.06) * max(abs(c(y, nb_y, glob)), na.rm = TRUE),
      main = sprintf('Forward selection by %s', forward_select_how(x)),
      sub  = sprintf('Stopped: %s', attr(x, 'stop', exact = TRUE)),
      xlab = 'Variables selected', ylab = expression('Cumulative adjusted' ~ R^2) ))

  graphics::axis(1, at = c(step, if (!is.null(nb_y)) k + 1))
  graphics::text(step, y, labs, pos = 1, cex = 0.8, xpd = NA)

  if (!is.null(nb_y) && is.finite(nb_y)) {
    graphics::segments(k, y[[k + 1]], k + 1, nb_y, lty = 3)
    graphics::points(k + 1, nb_y, pch = 1)
    graphics::text(k + 1, nb_y, sprintf('%s\n(next best)', nb$variable[[1]]), pos = 1, cex = 0.8, col = PLOT_MUTED, xpd = NA)
  }

  if (length(glob) && is.finite(glob)) {
    graphics::abline(h = glob, lty = 2, col = 'grey40')
    graphics::text(graphics::par('usr')[1], glob, 'global model', adj = c(-0.05, -0.4), cex = 0.8, col = 'grey40')
  }

  invisible(x)
}


# "RDA (hellinger)" or "db-RDA", as print.ecodive_forward_select() says it.
forward_select_how <- function (x) {
  switch(
    attr(x, 'method', exact = TRUE),
    rda   = sprintf('RDA%s', if (!identical(attr(x, 'norm'), 'none')) sprintf(' (%s)', attr(x, 'norm')) else ''),
    dbrda = 'db-RDA' )
}
