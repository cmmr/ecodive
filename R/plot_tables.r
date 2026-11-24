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
        xlab = cols[[1]], ylab = cols[[2]],
        main = plot_title(x, 'Fitted variables'),
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
#'        `cex` or `col` (its `color`); and `sub`, to [graphics::title()].
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
  if (n_taxa < 1)
    stop('`n_taxa` must be at least 1.', call. = FALSE)

  all_pairs <- unique(x$pair)
  if (is.null(pairs)) pairs <- all_pairs
  if (!is.character(pairs) || length(pairs) == 0 || !all(pairs %in% all_pairs))
    stop('`pairs` must be among ', name_list(all_pairs), '.', call. = FALSE)
  pairs <- all_pairs[all_pairs %in% pairs]

  overall <- attr(x, 'overall', exact = TRUE)
  rows    <- simper_plot_rows(x$pair, pairs, n_taxa)
  idx     <- rows$idx
  grp     <- rows$groups
  levels(grp) <- sprintf('%s (%s)', pairs, vapply(pairs, function (p) plot_num(overall[[p]]), ''))

  sig <- !is.na(x$p_adjusted[idx]) & x$p_adjusted[idx] <= alpha
  sub <- if (any(!is.na(x$p_adjusted)))
    sprintf('Filled: adjusted p <= %s', format(alpha))

  plot_do(
    graphics::dotchart, c(list(x$average[idx]), plot_dot_args(args)),
    list(
      labels = x$taxon[idx], groups = grp, pch = ifelse(sig, 19, 1),
      xlim = c(0, max(x$average[idx], na.rm = TRUE)),
      main = 'SIMPER', xlab = 'Average contribution to Bray-Curtis dissimilarity' ))

  graphics::title(sub = plot_arg(args, 'sub', sub))

  invisible(x)
}



# The rows of a simper() result to plot (the first n of each pair, which has
# its taxa largest first) in the order dotchart() wants them, and their
# pairs, a factor with levels `pairs`. dotchart() draws each group's values
# from the bottom up, so each pair's rows are reversed to put its largest at
# the top. It also sorts the rows by group, last group first, and before R
# 3.6.2 didn't sort `pch` with them; so the pairs are given in that order
# already, and the sort leaves them in place. Returns list(idx, groups).
simper_plot_rows <- function (pair, pairs, n) {
  rows <- lapply(rev(pairs), function (p) rev(utils::head(which(pair == p), n)))
  list(
    idx    = unlist(rows),
    groups = factor(rep(rev(pairs), lengths(rows)), levels = pairs) )
}



#' Plot a Dirichlet-multinomial t-test
#'
#' Four views of [dirmult_ttest()]'s results, each feature a point (or a
#' row), with the features that differ (`reject = TRUE`) highlighted and
#' the most significant of them named:
#'
#' * `"volcano"` (the default): the log2 fold change against the adjusted
#'   p-value.
#' * `"effect"`: the log2 fold change against the pooled within-group
#'   standard deviation (Gloor et al. 2016; ALDEx2's `"MW"` plot).
#' * `"ma"`: the log2 fold change against the mean clr abundance (ALDEx2's
#'   `"MA"` plot).
#' * `"forest"`: the log2 fold changes and their intervals of the
#'   `n_features` most significant features, largest at the top.
#' * `"rank"`: every feature's log2 fold change, sorted (see **Rank
#'   plots**).
#'
#' @param x  The result of [dirmult_ttest()].
#'
#' @param type  The plot: `"volcano"` (the default), `"effect"`, `"ma"`,
#'        `"forest"`, or `"rank"`.
#'
#' @param labels  How many of the differentially abundant features to name
#'        (in the volcano, effect and MA plots), those with the smallest
#'        adjusted p-values. In the rank plot, the highlighted features with
#'        the largest effects. Default `10`; `0` for none.
#'
#' @param n_features  For `"forest"`: how many features to show, those with
#'        the smallest adjusted p-values (ties in the table's order).
#'        Default `30`; `Inf` for all.
#'
#' @inheritParams plot.ecodive_ancombc
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim` or `cex`.
#'
#' @inheritSection plot.ecodive_ancombc Rank plots
#'
#' @details
#' **Volcano plot.** The y axis is \eqn{-\log_{10}} of `p_adjusted` (an
#' adjusted p-value of 0 is drawn at the smallest positive double). The
#' dashed line marks `alpha`; a feature also needs its interval for the
#' log2 fold change to exclude 0 to count as differentially abundant, so
#' some features above the line may not be highlighted. Gloor et al.
#' (2016) argue that volcano plots mislead: with many samples, tiny
#' differences reach small p-values, so read them with the effect plot.
#'
#' **Effect plot.** Each feature's difference between the groups (y,
#' `log2_fc`) against its dispersion within them (x, `pooled_sd`), both in
#' log2 units. The dashed lines are \eqn{y = \pm x}: features beyond them
#' differ between the groups by more than they vary within them (an
#' average `effect` of about 1 or more). Features near the y axis vary
#' little, so a small difference can be significant; features far to the
#' right need a large difference. This is ALDEx2's `aldex.plot(type =
#' "MW")`, with ALDEx3's pooled standard deviation in place of ALDEx2's
#' `diff.win` (see [dirmult_ttest()]'s **Effect sizes**).
#'
#' **MA plot.** The log2 fold change against the feature's mean clr value
#' (`mean`): how the differences depend on abundance. Rare features (on
#' the left) vary most, from sampling alone.
#'
#' **Forest plot.** A point at each feature's `log2_fc` and a bar from
#' `ci_lower` to `ci_upper`, which are the envelope of the draws'
#' intervals (the lowest lower and highest upper bound), not an interval
#' of one draw: wider than a single confidence interval.
#'
#' Features whose value is `NaN` are left out of each plot. The subtitle
#' names the contrast and counts the features rejected.
#'
#' @return `x`, invisibly.
#'
#' @references
#' Gloor, G. B., Macklaim, J. M., & Fernandes, A. D. (2016). Displaying
#' variation in large datasets: plotting a visual summary of effect sizes.
#' *Journal of Computational and Graphical Statistics*, 25, 971-979.
#' \doi{10.1080/10618600.2015.1131161}
#'
#' @export
#' @examples
#'     res <- dirmult_ttest(round(iris[, 1:4] * 10), 'Species', iris)
#'     plot(res)
#'     plot(res, 'effect')
#'     plot(res, 'ma')
#'     plot(res, 'forest')
#'     plot(res, 'rank', numerator = 'Petal.Width', denominator = 'Sepal.Width')
plot.ecodive_dirmult_ttest <- function (
    x, type = c('volcano', 'effect', 'ma', 'forest', 'rank'), labels = 10,
    n_features = 30, numerator = NULL, denominator = NULL, ... ) {

  type <- test_plot_type(type, c('volcano', 'effect', 'ma', 'forest', 'rank'))
  da_plot_need(x, type, c(effect = 'pooled_sd', ma = 'mean'), column = TRUE)

  conf <- attr(x, 'conf_level', exact = TRUE)
  da_plot(
    x, type, labels, n_features, list(...),
    est = x$log2_fc, p = x$p_adjusted, reject = x$reject, feature = x$feature,
    sd = x$pooled_sd, mean = x$mean, lo = x$ci_lower, hi = x$ci_upper,
    what     = da_plot_contrast(x),
    est_lab  = lfc_label(x),
    sd_lab   = 'Pooled within-group standard deviation (log2)',
    mean_lab = 'Mean clr abundance (log2)',
    interval = sprintf("envelope of the draws' %s intervals", plot_pct(conf)),
    num = numerator, den = denominator )

  invisible(x)
}



# "Log2 fold change", "Log fold change" (natural logs) or "Log10 fold
# change", from a DA result's `log_base` attribute, for axis labels.
# Without one (or for another base), "Log fold change", with the base when
# it is known.
lfc_label <- function (x) {
  b <- attr(x, 'log_base', exact = TRUE)
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= 0 || b == 1)
    return('Log fold change')
  if (abs(b - exp(1)) < 1e-12) return('Log fold change')
  if (b %in% c(2, 10))         return(sprintf('Log%g fold change', b))
  sprintf('Log fold change (base %g)', b)
}


# "95%", from a confidence level.
plot_pct <- function (level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level)) return ('')
  paste0(format(100 * level, digits = 6), '%')
}


# The subtitle of a DA plot: the comparison `what` (or NULL), the number of
# features rejected out of all of them, and with a sensitivity analysis
# (`robust`) the number robust to the pseudocount (filled), then `extra`.
da_plot_sub <- function (what, reject, robust = NULL, extra = NULL) {
  s <- sprintf(
    '%d of %d features differentially abundant', sum(reject %in% TRUE),
    length(reject) )
  if (!is.null(robust))
    s <- sprintf('%s, %d robust to the pseudocount (filled)', s, sum(robust %in% TRUE))
  s <- paste(c(if (length(what)) paste0(what, ':'), s), collapse = ' ')
  paste(c(s, extra), collapse = '; ')
}


# The highlighting of DA plots: colours from `reject`; with a sensitivity
# analysis (`robust`), the robust decisions filled and the others that are
# significant hollow; and which features to name (robust ones only, then).
# Returns list(col, pch, named).
da_plot_marks <- function (reject, robust, pch = 16) {
  rej <- reject %in% TRUE
  col <- ifelse(rej, PLOT_HIGHLIGHT, PLOT_MUTED)
  if (is.null(robust)) return (list(col = col, pch = pch, named = rej))
  rob <- robust %in% TRUE
  list(col = col, pch = ifelse(rej & !rob, 1, pch), named = rob)
}


# Names the `labels` features with the smallest `p` among `named`, at
# (x, y).
da_plot_labels <- function (x, y, p, feature, named, labels) {
  if (labels > 0 && any(named)) {
    i <- which(named)
    i <- utils::head(i[order(p[i])], labels)
    graphics::text(x[i], y[i], feature[i], pos = 3, cex = 0.75, xpd = NA)
  }
}


# A volcano plot: log fold changes `lfc` against -log10 of the adjusted
# p-values `p`, with the features that `reject` highlighted, the `labels`
# most significant of them named, and dashed lines at `alpha` and 0. `args`
# are the user's arguments for graphics::plot.default(). `sub` is the
# subtitle (da_plot_sub()).
plot_volcano <- function (
    lfc, p, reject, feature, alpha, labels, args, main, xlab, sub,
    robust = NULL ) {

  y  <- -log10(pmax(p, .Machine$double.xmin))
  mk <- da_plot_marks(reject, robust)

  plot_do(
    graphics::plot.default, c(list(lfc, y), args),
    list(
      col  = mk$col, pch = mk$pch, las = 1,
      ylim = range(0, y, finite = TRUE), main = main, xlab = xlab, sub = sub,
      ylab = expression(-log[10] ~ (adjusted ~ italic(p))) ))

  graphics::abline(v = 0, lty = 3, col = 'grey70')
  if (length(alpha) && is.finite(alpha) && alpha > 0)
    graphics::abline(h = -log10(alpha), lty = 2, col = 'grey40')

  da_plot_labels(lfc, y, p, feature, mk$named, labels)

  invisible()
}


# The effect and MA plots: log fold changes `y` against `x` (a dispersion,
# from 0, with dashed lines at y = x and y = -x, for `guides = 'effect'`;
# a mean abundance for 'ma'), with a dotted line at y = 0, highlighting
# and names as plot_volcano().
plot_da_scatter <- function (
    x, y, p, reject, feature, labels, args, main, xlab, ylab, sub,
    guides = c('effect', 'ma'), robust = NULL ) {

  guides <- match.arg(guides)
  mk     <- da_plot_marks(reject, robust)
  xlim   <- if (guides == 'effect') range(0, x, finite = TRUE) else range(x, finite = TRUE)

  plot_do(
    graphics::plot.default, c(list(x, y), args),
    list(
      col = mk$col, pch = mk$pch, las = 1, xlim = xlim,
      main = main, xlab = xlab, ylab = ylab, sub = sub ))

  graphics::abline(h = 0, lty = 3, col = 'grey70')
  if (guides == 'effect') {
    graphics::abline(0,  1, lty = 2, col = 'grey40')
    graphics::abline(0, -1, lty = 2, col = 'grey40')
  }

  da_plot_labels(x, y, p, feature, mk$named, labels)

  invisible()
}


# A forest plot: for the `n` features with the smallest `p` (ties in their
# order), a point at `est` and a bar from `lo` to `hi`, one row each,
# sorted by `est` (largest at the top) and labelled with `feature`, with
# the features that `reject` highlighted (robust ones filled) and a dotted
# line at 0. Features whose `est` is NaN are left out; a bar whose ends are
# NaN isn't drawn. Like dotchart(), it widens the left margin to fit the
# names, and restores it after.
plot_forest <- function (
    est, lo, hi, p, reject, feature, n, args, main, xlab, sub, robust = NULL ) {

  n <- plot_count(n, 'n_features')
  if (n < 1)
    stop('`n_features` must be at least 1.', call. = FALSE)

  ok  <- which(is.finite(est))
  top <- utils::head(ok[order(p[ok])], n)
  idx <- top[order(est[top])]         # bottom to top
  k   <- length(idx)
  mk  <- da_plot_marks(reject[idx], if (!is.null(robust)) robust[idx])
  if (k < length(ok)) sub <- sprintf('%s; top %d shown', sub, k)

  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  linch <- max(graphics::strwidth(feature[idx], 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)

  xlim <- range(0, est[idx], lo[idx], hi[idx], finite = TRUE)
  y    <- seq_len(k)

  plot_do(
    graphics::plot.default, c(list(est[idx], y), args),
    list(
      type = 'n', yaxt = 'n', xlim = xlim, ylim = c(0.5, max(k, 1) + 0.5),
      main = main, xlab = xlab, ylab = '', sub = sub, las = 1 ))

  graphics::abline(v = 0, lty = 3, col = 'grey70')
  graphics::abline(h = y, lty = 3, col = 'grey90')
  bar <- is.finite(lo[idx]) & is.finite(hi[idx])
  if (any(bar))
    graphics::segments(lo[idx][bar], y[bar], hi[idx][bar], y[bar], col = mk$col[bar], lwd = 1.5)
  graphics::points(est[idx], y, pch = mk$pch, col = mk$col, bg = 'white')
  if (k) graphics::axis(2, at = y, labels = feature[idx], tick = FALSE, las = 1, cex.axis = 0.8)

  invisible()
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
#'        `cex` or `col` (its `color`, replacing the highlighting); and
#'        `sub`, to [graphics::title()].
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

  D   <- nrow(x)
  cut <- attr(x, 'cutoff', exact = TRUE)

  plot_ranked(
    x$W, x$feature, x$reject, n_features, list(...),
    main = sprintf('ANCOM: %s', paste(attr(x, 'groups', exact = TRUE), collapse = ' vs ')),
    xlab = sprintf('W (of %d log ratios)', D - 1),
    xlim = c(0, max(D - 1, 1)),
    line = if (length(cut) && !is.na(cut)) cut * D )

  invisible(x)
}



# A dot chart of the `n` features with the largest `stat` (largest at the
# top), labelled with `labels`, those that `reject` highlighted, and a
# dashed line at `line` (if not NULL). `args` are the user's arguments for
# graphics::dotchart(). The subtitle counts the features rejected.
plot_ranked <- function (
    stat, labels, reject, n, args, main, xlab, xlim = NULL, line = NULL,
    robust = NULL ) {

  n <- plot_count(n, 'n_features')
  if (n < 1)
    stop('`n_features` must be at least 1.', call. = FALSE)

  D   <- length(stat)
  ord <- utils::head(order(stat, decreasing = TRUE), n)
  idx <- rev(ord)  # dotchart() draws from the bottom up.
  rej <- reject[idx] %in% TRUE

  sub <- sprintf(
    '%d of %d features differentially abundant%s%s', sum(reject %in% TRUE), D,
    if (is.null(robust)) '' else sprintf(', %d robust to the pseudocount (filled)', sum(robust %in% TRUE)),
    if (length(ord) < D) sprintf('; top %d shown', length(ord)) else '' )

  # With a sensitivity analysis, significant but not robust: hollow.
  pch <- if (is.null(robust)) 19 else ifelse(rej & !(robust[idx] %in% TRUE), 1, 19)

  if (is.null(xlim)) xlim <- range(0, stat[idx], finite = TRUE)

  plot_do(
    graphics::dotchart, c(list(stat[idx]), plot_dot_args(args)),
    list(
      labels = labels[idx], pch = pch,
      color  = ifelse(rej, PLOT_HIGHLIGHT, PLOT_MUTED),
      xlim   = xlim, main = main, xlab = xlab ))

  graphics::title(sub = plot_arg(args, 'sub', sub))
  if (length(line) && is.finite(line))
    graphics::abline(v = line, lty = 2, col = 'grey40')

  invisible()
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

  args <- list(...)
  best <- attr(x, 'best', exact = TRUE)
  i    <- match(length(best), x$size)

  plot_do(
    graphics::plot.default, c(list(x$size, x$correlation), args),
    list(
      type = 'b', pch = 16, xaxt = 'n', las = 1,
      main = if (isTRUE(attr(x, 'partial', exact = TRUE))) 'BIO-ENV (partial)' else 'BIO-ENV',
      sub  = paste('Best:', paste(best, collapse = ', ')),
      xlab = 'Number of variables',
      ylab = sprintf('%s correlation with distances', plot_cap(attr(x, 'method', exact = TRUE))) ))

  if (plot_x_axis(args)) graphics::axis(1, at = x$size)
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

  args <- list(...)
  g    <- attr(x, 'global', exact = TRUE)
  nb   <- attr(x, 'next_best', exact = TRUE)
  k    <- nrow(x)
  step <- 0:k
  y    <- c(0, x$adj_r_squared_cum)
  labs <- c('(none)', x$variable)

  nb_y <- if (!is.null(nb) && nrow(nb)) nb$adj_r_squared_cum[[1]]
  glob <- if (!is.null(g)) g$adj_r_squared[[1]]

  plot_do(
    graphics::plot.default, c(list(step, y), args),
    list(
      type = 'b', pch = 16, xaxt = 'n', las = 1,
      xlim = c(-0.3, k + if (!is.null(nb_y)) 1.3 else 0.3),
      ylim = range(y, nb_y, glob, finite = TRUE) + c(-0.06, 0.06) * max(abs(c(y, nb_y, glob)), na.rm = TRUE),
      main = sprintf('Forward selection by %s', forward_select_how(x)),
      sub  = sprintf('Stopped: %s', attr(x, 'stop', exact = TRUE)),
      xlab = 'Variables selected', ylab = expression('Cumulative adjusted' ~ R^2) ))

  if (plot_x_axis(args)) graphics::axis(1, at = c(step, if (!is.null(nb_y)) k + 1))
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
