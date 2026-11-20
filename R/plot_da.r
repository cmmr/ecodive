# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Plot a Per-feature Test Result
#'
#' A volcano plot of the results of [da_wilcoxon()], [da_ttest()] or
#' [da_lm()]: each feature's effect against its adjusted p-value, with the
#' features declared differentially abundant (`reject = TRUE`) highlighted
#' and the most significant of them named.
#'
#' @param x  The result of [da_wilcoxon()], [da_ttest()], or [da_lm()].
#'
#' @param covariate  For [da_lm()]: which covariate's coefficients to plot,
#'        as named in the result's `covariate` column (such as
#'        `"statusmoderate"` or `"age"`). Defaults to the first column of the
#'        term of interest (`groups`), or else to the first covariate that
#'        isn't the intercept.
#'
#' @param labels  How many of the differentially abundant features to name,
#'        those with the smallest adjusted p-values. Default `10`; `0` for
#'        none.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim` or `cex`.
#'
#' @details
#' The x axis is the effect: the difference in medians ([da_wilcoxon()]),
#' in means ([da_ttest()]), or the coefficient ([da_lm()]), in the units of
#' the result's `norm` (natural logs for `'clr'`, `'rclr'` and `'log'`).
#' The y axis is \eqn{-\log_{10}} of `p_adjusted` (an adjusted p-value of
#' 0 is drawn at the smallest positive double), with a dashed line at
#' `alpha`. Features whose effect or p-value is `NaN` are left out.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     x   <- ex_counts + 1
#'     grp <- c('oral', 'oral', 'other', 'other')
#'     plot(da_ttest(x, grp))
#'     plot(da_lm(x, groups = grp))
plot.ecodive_da_ttest <- function (x, labels = 10, ...) {
  da_plot_two_sample(x, x$diff, 'means', labels, list(...))
  invisible(x)
}


#' @rdname plot.ecodive_da_ttest
#' @export
plot.ecodive_da_wilcoxon <- function (x, labels = 10, ...) {
  da_plot_two_sample(x, x$median_diff, 'medians', labels, list(...))
  invisible(x)
}


#' @rdname plot.ecodive_da_ttest
#' @export
plot.ecodive_da_lm <- function (x, covariate = NULL, labels = 10, ...) {

  labels <- plot_count(labels, 'labels')

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

  r <- x[x$covariate == covariate & is.finite(x$estimate), , drop = FALSE]

  plot_volcano(
    r$estimate, r$p_adjusted, r$reject, r$feature, attr(x, 'alpha', exact = TRUE),
    labels, list(...),
    main = da_plot_title(x),
    xlab = sprintf('Coefficient%s, %s', da_plot_units(x), covariate) )

  invisible(x)
}



# The two-sample tests' volcano plot.
da_plot_two_sample <- function (x, effect, what, labels, args) {

  labels <- plot_count(labels, 'labels')

  ct  <- attr(x, 'contrast', exact = TRUE)
  n   <- attr(x, 'n', exact = TRUE)
  ref <- if (is.na(ct[[2]])) 'all others' else ct[[2]]
  ok  <- is.finite(effect)

  plot_volcano(
    effect[ok], x$p_adjusted[ok], x$reject[ok], x$feature[ok],
    attr(x, 'alpha', exact = TRUE), labels, args,
    main = da_plot_title(x),
    xlab = sprintf(
      'Difference in %s%s, %s (n = %d) - %s (n = %d)',
      what, da_plot_units(x), ct[[1]], n[[1]], ref, n[[2]] ) )
}


# "Volcano plot of Welch t-test", from the result's `method`.
da_plot_title <- function (x) {
  paste('Volcano plot of', da_method_name(attr(x, 'method', exact = TRUE)))
}


# " (clr)", the result's norm for the axis label.
da_plot_units <- function (x) {
  norm <- attr(x, 'norm', exact = TRUE)
  if (is.character(norm) && length(norm) == 1L && norm != 'none') sprintf(' (%s)', norm) else ''
}
