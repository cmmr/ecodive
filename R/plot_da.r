# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Plot a Per-feature Test Result
#'
#' Views of the results of [da_wilcoxon()], [da_ttest()] or [da_lm()], each
#' feature a point (or a row), with the features declared differentially
#' abundant (`reject = TRUE`) highlighted and the most significant of them
#' named:
#'
#' * `"volcano"` (the default): the effect against the adjusted p-value.
#' * `"effect"` ([da_ttest()], [da_lm()]): the effect against the
#'   within-group standard deviation (Gloor et al. 2016; ALDEx2's `"MW"`
#'   plot).
#' * `"ma"`: the effect against the mean abundance (ALDEx2's `"MA"` plot).
#' * `"forest"` ([da_ttest()], [da_lm()]): the effects and their confidence
#'   intervals of the `n_features` most significant features, largest at
#'   the top.
#'
#' @param x  The result of [da_wilcoxon()], [da_ttest()], or [da_lm()].
#'
#' @param type  The plot: `"volcano"` (the default), `"effect"`, `"ma"`, or
#'        `"forest"` (for [da_wilcoxon()], `"volcano"` or `"ma"`).
#'
#' @param covariate  For [da_lm()]: which covariate's coefficients to plot,
#'        as named in the result's `covariate` column (such as
#'        `"statusmoderate"` or `"age"`). Defaults to the first column of the
#'        term of interest (`groups`), or else to the first covariate that
#'        isn't the intercept.
#'
#' @param labels  How many of the differentially abundant features to name
#'        (in the volcano, effect and MA plots), those with the smallest
#'        adjusted p-values. Default `10`; `0` for none.
#'
#' @param n_features  For `"forest"`: how many features to show, those with
#'        the smallest adjusted p-values (ties in the table's order).
#'        Default `30`; `Inf` for all.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim` or `cex`.
#'
#' @details
#' The effect is the difference in medians ([da_wilcoxon()]), in means
#' ([da_ttest()]), or the coefficient ([da_lm()]), in the units of the
#' result's `norm`: a log fold change for `'clr'`, `'rclr'` and `'log'`
#' (natural logs, from the `log_base` attribute), a difference in
#' proportions for `'percent'`.
#'
#' * **Volcano plot.** The y axis is \eqn{-\log_{10}} of `p_adjusted` (an
#'   adjusted p-value of 0 is drawn at the smallest positive double), with
#'   a dashed line at `alpha`.
#' * **Effect plot.** The x axis is the within-group dispersion:
#'   [da_ttest()]'s `pooled_sd`, or [da_lm()]'s `resid_sd` attribute (the
#'   residual standard deviation, which for a model of one two-level factor
#'   is the pooled within-group standard deviation). The dashed lines are
#'   \eqn{y = \pm x}: features beyond them differ between the groups by
#'   more than they vary within them. For [da_lm()], only a covariate whose
#'   design column is 0 or 1 (a factor's level against its reference) has
#'   one; use `"ma"` or `"forest"` for others.
#' * **MA plot.** The x axis is the feature's mean transformed abundance:
#'   the `mean` column ([da_wilcoxon()], [da_ttest()]) or attribute
#'   ([da_lm()]).
#' * **Forest plot.** Points at the effects, bars over their confidence
#'   intervals (`ci_lower`, `ci_upper`, at the result's `conf_level`, not
#'   adjusted for testing many features).
#'
#' Features whose values are `NaN` are left out. The subtitle names the
#' comparison and counts the features rejected.
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
#'     x   <- ex_counts + 1
#'     grp <- c('oral', 'oral', 'other', 'other')
#'     res <- da_ttest(x, grp)
#'     plot(res)
#'     plot(res, 'effect')
#'     plot(res, 'forest', n_features = 10)
#'     plot(da_lm(x, groups = grp), type = 'ma')
plot.ecodive_da_ttest <- function (
    x, type = c('volcano', 'effect', 'ma', 'forest'), labels = 10,
    n_features = 30, ... ) {

  type <- test_plot_type(type, c('volcano', 'effect', 'ma', 'forest'))
  da_plot_need(x, type, c(effect = 'pooled_sd', ma = 'mean'), column = TRUE)

  units <- da_norm_units(x)
  da_plot(
    x, type, labels, n_features, list(...),
    est = x$diff, p = x$p_adjusted, reject = x$reject, feature = x$feature,
    sd = x$pooled_sd, mean = x$mean, lo = x$ci_lower, hi = x$ci_upper,
    what     = da_plot_contrast(x),
    est_lab  = da_effect_label(x, 'Difference in means'),
    sd_lab   = sprintf('Pooled within-group standard deviation (%s)', units),
    mean_lab = sprintf('Mean abundance (%s)', units),
    interval = sprintf('%s confidence intervals', plot_pct(attr(x, 'conf_level', exact = TRUE))) )

  invisible(x)
}


#' @rdname plot.ecodive_da_ttest
#' @export
plot.ecodive_da_wilcoxon <- function (x, type = c('volcano', 'ma'), labels = 10, ...) {

  type <- test_plot_type(type, c('volcano', 'ma'))
  da_plot_need(x, type, c(ma = 'mean'), column = TRUE)

  da_plot(
    x, type, labels, 30, list(...),
    est = x$median_diff, p = x$p_adjusted, reject = x$reject, feature = x$feature,
    mean     = x$mean,
    what     = da_plot_contrast(x),
    est_lab  = da_effect_label(x, 'Difference in medians'),
    mean_lab = sprintf('Mean abundance (%s)', da_norm_units(x)) )

  invisible(x)
}


#' @rdname plot.ecodive_da_ttest
#' @export
plot.ecodive_da_lm <- function (
    x, covariate = NULL, type = c('volcano', 'effect', 'ma', 'forest'),
    labels = 10, n_features = 30, ... ) {

  type <- test_plot_type(type, c('volcano', 'effect', 'ma', 'forest'))
  da_plot_need(x, type, c(effect = 'resid_sd', ma = 'mean'), column = FALSE)

  covs <- unique(x$covariate)
  if (is.null(covariate)) {
    g    <- attr(x, 'groups', exact = TRUE)
    asg  <- attr(attr(x, 'design', exact = TRUE), 'assign')
    term <- attr(attr(x, 'design', exact = TRUE), 'term.labels')
    covariate <- NULL
    if (is.character(g) && length(g) == 1L && !is.null(asg) && g %in% term)
      covariate <- utils::head(covs[asg == match(g, term)], 1L)
    if (!length(covariate)) {
      covariate <- setdiff(covs, '(Intercept)')
      covariate <- if (length(covariate)) covariate[[1]] else covs[[1]]
    }
  }
  if (!is.character(covariate) || length(covariate) != 1 || !covariate %in% covs)
    stop('`covariate` must be one of ', name_list(covs), '.', call. = FALSE)
  if (type == 'effect') da_plot_binary(x, covariate)

  r     <- x[x$covariate == covariate, , drop = FALSE]
  units <- da_norm_units(x)
  da_plot(
    x, type, labels, n_features, list(...),
    est = r$estimate, p = r$p_adjusted, reject = r$reject, feature = r$feature,
    sd = attr(x, 'resid_sd', exact = TRUE), mean = attr(x, 'mean', exact = TRUE),
    lo = r$ci_lower, hi = r$ci_upper,
    what     = covariate,
    est_lab  = da_effect_label(x, 'Coefficient'),
    sd_lab   = sprintf('Residual standard deviation (%s)', units),
    mean_lab = sprintf('Mean abundance (%s)', units),
    interval = sprintf('%s confidence intervals', plot_pct(attr(x, 'conf_level', exact = TRUE))) )

  invisible(x)
}



# Draws a DA result's plot of `type` ('volcano', 'effect', 'ma' or
# 'forest') from per-feature vectors: the effects `est`, adjusted p-values
# `p`, decisions `reject` (and `robust`, with a sensitivity analysis),
# names `feature`, dispersions `sd`, mean abundances `mean`, and interval
# bounds `lo`, `hi`. `what` is the comparison (for the subtitle), the
# `*_lab` are axis labels and `interval` describes the forest plot's bars
# (in its x axis label).
# The title is da_plot_title()'s.
da_plot <- function (
    x, type, labels, n_features, args, est, p, reject, feature,
    sd = NULL, mean = NULL, lo = NULL, hi = NULL, robust = NULL, what = NULL,
    est_lab = 'Log fold change', sd_lab = NULL, mean_lab = NULL,
    interval = NULL ) {

  labels <- plot_count(labels, 'labels')
  main   <- da_plot_title(x, type)
  sub    <- da_plot_sub(what, reject, robust)
  rob    <- function (i) if (!is.null(robust)) robust[i]

  if (type == 'volcano') {
    i <- which(is.finite(est))
    plot_volcano(
      est[i], p[i], reject[i], feature[i], attr(x, 'alpha', exact = TRUE),
      labels, args, main = main, xlab = est_lab, sub = sub, robust = rob(i) )
  }
  else if (type == 'forest') {
    plot_forest(
      est, lo, hi, p, reject, feature, n_features, args, main = main,
      xlab = sprintf('%s; bars: %s', est_lab, interval), sub = sub,
      robust = robust )
  }
  else {
    xv <- if (type == 'effect') sd else mean
    i  <- which(is.finite(est) & is.finite(xv))
    plot_da_scatter(
      xv[i], est[i], p[i], reject[i], feature[i], labels, args, main = main,
      xlab = if (type == 'effect') sd_lab else mean_lab, ylab = est_lab,
      sub = sub, guides = type, robust = rob(i) )
  }

  invisible()
}


# "Volcano plot of Welch t-test", from the plot's type and the result's
# `method`.
da_plot_title <- function (x, type = 'volcano') {
  lab <- c(volcano = 'Volcano plot', effect = 'Effect plot', ma = 'MA plot', forest = 'Forest plot')
  paste(lab[[type]], 'of', da_method_name(attr(x, 'method', exact = TRUE)))
}


# "treatment (n = 3) vs all others (n = 4)", from the result's `contrast`
# and `n`.
da_plot_contrast <- function (x) {
  ct  <- attr(x, 'contrast', exact = TRUE)
  n   <- attr(x, 'n', exact = TRUE)
  ref <- if (is.na(ct[[2]])) 'all others' else ct[[2]]
  sprintf('%s (n = %d) vs %s (n = %d)', ct[[1]], n[[1]], ref, n[[2]])
}


# The units of a da_*() result's values: "clr, natural log" for a log norm
# (from `log_base`), else the norm's name ("percent").
da_norm_units <- function (x) {
  norm <- attr(x, 'norm', exact = TRUE)
  if (!is.character(norm) || length(norm) != 1L) norm <- 'none'
  b <- attr(x, 'log_base', exact = TRUE)
  if (is.numeric(b) && length(b) == 1L && isTRUE(abs(b - exp(1)) < 1e-12))
    sprintf('%s, natural log', norm)
  else norm
}


# The effect's axis label for a da_*() result: "Log fold change (difference
# in means, clr)" for a log norm, else "Difference in means (percent)".
da_effect_label <- function (x, what) {
  norm <- attr(x, 'norm', exact = TRUE)
  if (!is.character(norm) || length(norm) != 1L) norm <- 'none'
  b <- attr(x, 'log_base', exact = TRUE)
  if (is.numeric(b) && length(b) == 1L && is.finite(b))
    sprintf('%s (%s, %s)', lfc_label(x), tolower(what), norm)
  else if (norm == 'none') what
  else sprintf('%s (%s)', what, norm)
}


# Stops unless the result has what the plot of `type` reads: `needs` names
# a column (column = TRUE) or attribute per type. Results made before
# ecodive 2.4.0 have neither.
da_plot_need <- function (x, type, needs, column) {
  if (!type %in% names(needs)) return (invisible())
  nm  <- needs[[type]]
  has <- if (column) nm %in% names(x) else !is.null(attr(x, nm, exact = TRUE))
  if (!has)
    stop(
      'This result has no `', nm, '` ', if (column) 'column' else 'attribute',
      ', which `type = "', type, '"` plots: fit it again with this version ',
      'of ecodive.', call. = FALSE )
  invisible()
}


# Stops unless `covariate`'s design column is 0 or 1 (a factor's level
# against its reference), which the effect plot needs.
da_plot_binary <- function (x, covariate) {
  X   <- attr(x, 'design', exact = TRUE)
  col <- if (!is.null(X) && covariate %in% colnames(X)) X[, covariate]
  if (is.null(col) || !all(col %in% c(0, 1)) || length(unique(col)) != 2L)
    stop(
      '`type = "effect"` needs a two-level comparison (a design column of ',
      '0s and 1s), and `', covariate, '` isn\'t one; use "ma" or "forest".',
      call. = FALSE )
  invisible()
}
