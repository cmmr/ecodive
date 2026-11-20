# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Plot an ANCOM-BC result
#'
#' Plots the results of [ancombc()] and [ancombc2()] and of their post-hoc
#' tests, from what each holds:
#'
#' * [ancombc()], [ancombc2()]: one covariate's bias-corrected log fold
#'   changes, as one of four plots (`type`): `"volcano"` (the default),
#'   against their adjusted p-values; `"effect"`, against each feature's
#'   residual standard deviation (Gloor et al. 2016; ALDEx2's `"MW"`
#'   plot), for a two-level comparison; `"ma"`, against each feature's mean
#'   clr abundance (ALDEx2's `"MA"` plot); or `"forest"`, with their
#'   confidence intervals, for the `n_features` most significant features.
#' * [ancombc_pairwise()], [ancombc_dunnett()]: a heatmap of the log fold
#'   change of each feature (rows) in each comparison (columns), with the
#'   significant ones in bold.
#' * [ancombc_global()], [ancombc_trend()]: a dot chart of each feature's
#'   test statistic, largest first.
#'
#' Features declared differentially abundant (`reject = TRUE`) are
#' highlighted in all of them. When the fit was made with a sensitivity
#' analysis (`sensitivity`, see [ancombc()]), the robust decisions
#' (`reject_robust`) are drawn filled (in bold in the heatmap) and named,
#' and the others that are significant hollow (in italics).
#'
#' @param x  The result of [ancombc()], [ancombc2()], [ancombc_global()],
#'        [ancombc_pairwise()], [ancombc_dunnett()], or [ancombc_trend()].
#'
#' @param covariate  For [ancombc()] and [ancombc2()]: which covariate's log
#'        fold changes to plot, as named in the result's `covariate` column
#'        (such as `"statusmoderate"` or `"age"`). Defaults to the first
#'        covariate that isn't the intercept.
#'
#' @param type  For [ancombc()] and [ancombc2()]: the plot, `"volcano"`
#'        (the default), `"effect"`, `"ma"`, or `"forest"`.
#'
#' @param labels  For [ancombc()] and [ancombc2()]: how many of the
#'        differentially abundant features to name (in the volcano, effect
#'        and MA plots), those with the smallest adjusted p-values. Default
#'        `10`; `0` for none.
#'
#' @param n_features  How many features to show: for the forest plot, those
#'        with the smallest adjusted p-value (ties in the table's order);
#'        for the heatmap, those with the smallest adjusted p-value in any
#'        comparison; for the dot chart, those with the largest statistic.
#'        Default `30`; `Inf` for all.
#'
#' @param conf_level  For the forest plot: the confidence level of the
#'        intervals. Default `0.95`.
#'
#' @param digits  For the heatmap: the number of decimal places of the log
#'        fold changes written in its cells. Default `2`.
#'
#' @param ...  Passed on to [graphics::plot.default()] for the volcano,
#'        effect, MA and forest plots, [graphics::image()] for the heatmap (where `col` is the colour
#'        scale, from most negative to most positive), and
#'        [graphics::dotchart()] for the dot chart (`col` as its `color`),
#'        such as `main`, `sub`, `xlab` or `cex`.
#'
#' @details
#' Log fold changes are natural logs, as ANCOM-BC estimates them; the axis
#' label and the heatmap's subtitle say so from the result's `log_base`
#' attribute. Coefficients that couldn't be estimated (`NaN`) are left out.
#'
#' * **Volcano plot.** The y axis is \eqn{-\log_{10}} of `p_adjusted`,
#'   with a dashed line at `alpha`.
#' * **Effect plot.** The x axis is the feature's residual standard
#'   deviation (the `resid_sd` attribute, natural logs); the dashed lines
#'   are \eqn{y = \pm x}, so features beyond them differ by more than they
#'   vary within the groups. Only for a covariate whose design column is 0
#'   or 1 (a factor's level against its reference): for a model of that
#'   factor alone, `resid_sd` is the pooled within-group standard
#'   deviation; with other covariates, it is the spread left after them.
#' * **MA plot.** The x axis is the feature's mean clr abundance (the `mean`
#'   attribute; the rclr with zeros as missing values).
#' * **Forest plot.** Bars over Wald intervals, \eqn{\mathrm{log\_fc} \pm
#'   q\,\mathrm{se}}, with \eqn{q} the \eqn{(1 + }`conf_level`\eqn{)/2}
#'   quantile of the standard normal (`ancombc()`) or of the t distribution
#'   with the feature's `dof` (`ancombc2()`), as the p-values; they aren't
#'   adjusted for testing many features, so a bar can miss 0 while the
#'   feature isn't rejected. The ANCOM-BC2 vignette's bar charts of log
#'   fold changes and standard errors show the same.
#'
#' In the heatmap, red cells are increases and blue cells
#' decreases, on a scale symmetric about 0, and the value of each cell is
#' written in it, in bold when that comparison is significant for that
#' feature. For [ancombc_trend()], each feature's label gives the pattern
#' with the largest statistic.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     x <- rbind(
#'       c(1, 10, 16, 0, 0, 12, 0, 7), c(0, 10, 45, 2, 1, 25, 1, 9),
#'       c(4,  8, 22, 1, 5, 14, 1, 5), c(1, 11, 22, 1, 4, 14, 1, 3),
#'       c(4,  7, 36, 0, 0, 19, 0, 9), c(2, 30, 26, 0, 0, 12, 1, 0),
#'       c(1, 36, 10, 1, 5, 15, 2, 9), c(1, 25, 11, 2, 2,  8, 0, 5),
#'       c(1, 36, 23, 0, 2, 11, 1, 2), c(3, 60, 33, 3, 4, 14, 1, 6),
#'       c(3, 57, 30, 3, 1,  3, 0, 4), c(0, 65, 23, 0, 2,  2, 0, 8),
#'       c(2, 46, 19, 0, 1,  2, 0, 1), c(0, 41, 26, 0, 2,  0, 0, 3),
#'       c(3, 46, 20, 0, 0,  3, 0, 0) )
#'     dimnames(x) <- list(paste0('S', 1:15), paste0('F', 1:8))
#'     meta <- data.frame(
#'       status = rep(c('mild', 'moderate', 'severe'), each = 5),
#'       age    = c(25, 48, 23, 35, 52, 51, 18, 29, 40, 29, 44, 39, 26, 37, 46),
#'       row.names = rownames(x) )
#'     fit <- ancombc2(x, ~ status + age, meta, groups = 'status')
#'
#'     # Moderate against mild, then age.
#'     plot(fit)
#'     plot(fit, 'age')
#'
#'     # Each feature's difference against its spread, its abundance, and
#'     # with its interval.
#'     plot(fit, type = 'effect')
#'     plot(fit, type = 'ma')
#'     plot(fit, 'age', type = 'forest')
#'
#'     # Every pair of groups.
#'     plot(ancombc_pairwise(fit))
#'
#'     # Any difference between the groups.
#'     plot(ancombc_global(fit))
plot.ecodive_ancombc <- function (
    x, covariate = NULL, type = c('volcano', 'effect', 'ma', 'forest'),
    labels = 10, n_features = 30, conf_level = 0.95, ... ) {

  type <- test_plot_type(type, c('volcano', 'effect', 'ma', 'forest'))
  da_plot_need(x, type, c(effect = 'resid_sd', ma = 'mean'), column = FALSE)

  covs <- unique(x$covariate)
  if (is.null(covariate)) {
    covariate <- setdiff(covs, '(Intercept)')
    covariate <- if (length(covariate)) covariate[[1]] else covs[[1]]
  }
  if (!is.character(covariate) || length(covariate) != 1 || !covariate %in% covs)
    stop('`covariate` must be one of ', name_list(covs), '.', call. = FALSE)
  if (type == 'effect') da_plot_binary(x, covariate)

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)

  # Wald intervals, with the p-values' distribution.
  r   <- x[x$covariate == covariate, , drop = FALSE]
  dof <- attr(x, 'dof', exact = TRUE)
  q   <- if (is.null(dof)) stats::qnorm((1 + conf_level) / 2)
         else suppressWarnings(stats::qt((1 + conf_level) / 2, dof))
  v2  <- !is.null(dof)

  mean_lab <- sprintf(
    'Mean abundance (%s, natural log)',
    if (v2 && isTRUE(attr(x, 'pseudocount', exact = TRUE) == 0)) 'rclr' else 'clr' )

  da_plot(
    x, type, labels, n_features, list(...),
    est = r$log_fc, p = r$p_adjusted, reject = r$reject, feature = r$feature,
    sd = attr(x, 'resid_sd', exact = TRUE), mean = attr(x, 'mean', exact = TRUE),
    lo = r$log_fc - q * r$se, hi = r$log_fc + q * r$se, robust = r$reject_robust,
    what     = covariate,
    est_lab  = lfc_label(x),
    sd_lab   = 'Residual standard deviation (natural log)',
    mean_lab = mean_lab,
    interval = sprintf('unadjusted %s Wald intervals', plot_pct(conf_level)) )

  invisible(x)
}


#' @rdname plot.ecodive_ancombc
#' @export
plot.ecodive_ancombc_pairwise <- function (x, n_features = 30, digits = 2L, ...) {
  ancombc_plot_heatmap(x, n_features, digits, list(...), 'pairwise comparisons')
  invisible(x)
}


#' @rdname plot.ecodive_ancombc
#' @export
plot.ecodive_ancombc_dunnett <- function (x, n_features = 30, digits = 2L, ...) {
  ancombc_plot_heatmap(x, n_features, digits, list(...), "Dunnett's comparisons")
  invisible(x)
}


#' @rdname plot.ecodive_ancombc
#' @export
plot.ecodive_ancombc_global <- function (x, n_features = 30, ...) {
  plot_ranked(
    x$w, x$feature, x$reject, n_features, list(...), robust = x$reject_robust,
    main = sprintf('%s global test', ancombc_plot_name(x)),
    xlab = sprintf('W, any difference between the levels of %s', ancombc_plot_groups(x)) )
  invisible(x)
}


#' @rdname plot.ecodive_ancombc
#' @export
plot.ecodive_ancombc_trend <- function (x, n_features = 30, ...) {
  plot_ranked(
    x$w, sprintf('%s (%s)', x$feature, x$pattern), x$reject, n_features, list(...),
    robust = x$reject_robust,
    main = sprintf('%s trend test', ancombc_plot_name(x)),
    xlab = sprintf('W, a trend over the levels of %s', ancombc_plot_groups(x)) )
  invisible(x)
}



# "ANCOM-BC" or "ANCOM-BC2", from the result's `method` attribute.
ancombc_plot_name <- function (x) {
  m <- attr(x, 'method', exact = TRUE)
  switch(paste0('', m), ancombc = 'ANCOM-BC', ancombc2 = 'ANCOM-BC2', toupper(paste0('', m)))
}


# The grouping variable's name, from the `groups` attribute, or "the groups".
ancombc_plot_groups <- function (x) {
  g <- attr(x, 'groups', exact = TRUE)
  if (is.character(g) && length(g) == 1 && !is.na(g)) g else 'the groups'
}



# A heatmap of the log fold changes of ancombc_pairwise() and
# ancombc_dunnett() results: a row per feature (the `n` with the smallest
# adjusted p-value in any comparison, the smallest at the top), a column per
# comparison, each cell coloured by its value on a scale symmetric about 0,
# and its value written in it (bold when rejected). Like dotchart(), it
# widens the left margin to fit the features' names, and restores it after.
ancombc_plot_heatmap <- function (x, n, digits, args, what) {

  n <- plot_count(n, 'n_features')
  if (n < 1)
    stop('`n_features` must be at least 1.', call. = FALSE)
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) || digits < 0)
    stop('`digits` must be a number of at least 0.', call. = FALSE)

  feats <- unique(x$feature)
  comps <- unique(x$comparison)

  # Features by their smallest adjusted p-value, in their order for ties.
  best  <- vapply(feats, function (f) {
    p <- x$p_adjusted[x$feature == f]
    if (any(!is.na(p))) min(p, na.rm = TRUE) else Inf
  }, 0)
  keep  <- utils::head(feats[order(best)], n)

  cell <- function (col) {
    m <- matrix(NA, length(keep), length(comps), dimnames = list(keep, comps))
    i <- x$feature %in% keep
    m[cbind(match(x$feature[i], keep), match(x$comparison[i], comps))] <- col[i]
    m
  }
  lfc <- cell(x$log_fc)
  rej <- cell(x$reject %in% TRUE)
  rej[is.na(rej)] <- FALSE
  # With a sensitivity analysis, significant but not robust: italic.
  sens <- !is.null(x$reject_robust)
  rob  <- if (sens) cell(x$reject_robust %in% TRUE) else rej
  rob[is.na(rob)] <- FALSE

  # Rows from the top: image() draws its first row at the bottom.
  rows <- rev(seq_along(keep))
  fin  <- is.finite(lfc)
  top  <- if (any(fin)) max(abs(lfc[fin])) else 0
  if (top == 0) top <- 1
  pal  <- plot_arg(args, 'col', grDevices::hcl.colors(21, 'Blue-Red 3'))

  # Room for the features' names on the left, as dotchart() makes it.
  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  linch <- max(graphics::strwidth(keep, 'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  graphics::par(mai = mai)

  n_rej <- length(unique(x$feature[x$reject %in% TRUE]))
  sub   <- sprintf(
    '%ss; bold: differentially abundant (%d of %d features)%s%s',
    lfc_label(x), n_rej, length(feats),
    if (sens) '; italic: not robust to the pseudocount' else '',
    if (length(keep) < length(feats)) sprintf('; top %d shown', length(keep)) else '' )

  plot_do(
    graphics::image, c(list(seq_along(comps), seq_along(keep), t(lfc[rows, , drop = FALSE])), args),
    list(
      col = pal, breaks = seq(-top, top, length.out = length(pal) + 1L),
      axes = FALSE, xlab = '', ylab = '', sub = sub,
      main = sprintf('%s %s of %s', ancombc_plot_name(x), what, ancombc_plot_groups(x)) ))

  graphics::axis(1, at = seq_along(comps), labels = comps, tick = FALSE, cex.axis = 0.8)
  graphics::axis(2, at = seq_along(keep), labels = keep[rows], tick = FALSE, las = 1, cex.axis = 0.8)
  graphics::box()

  # Each cell's value, bold when rejected, and white on the darkest cells.
  xy  <- expand.grid(r = seq_along(keep), c = seq_along(comps))
  val <- lfc[cbind(rows[xy$r], xy$c)]
  hit <- rej[cbind(rows[xy$r], xy$c)]
  rbs <- rob[cbind(rows[xy$r], xy$c)]
  ok  <- is.finite(val)
  if (any(ok)) {
    v   <- round(val[ok], digits)
    v[v == 0] <- 0  # not "-0.00"
    graphics::text(
      xy$c[ok], xy$r[ok], formatC(v, format = 'f', digits = digits),
      font = ifelse(hit[ok], ifelse(rbs[ok], 2, 3), 1), cex = 0.75,
      col  = ifelse(abs(val[ok]) > 0.55 * top, 'white', ifelse(hit[ok], 'black', 'grey30')) )
  }

  invisible()
}
