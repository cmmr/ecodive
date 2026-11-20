# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# ancombc_global(), ancombc_pairwise(), ancombc_dunnett(), ancombc_trend()
# and structural_zeros() are an R/C port of ANCOMBCResult.global_test(),
# .pairwise_test(), .dunnett_test() and .trend_test() (_global_test(),
# _global_stats(), _pairwise_test(), _mdfdr_pairwise(), _dunnett_test(),
# _mdfdr_dunnett(), _dunn_global(), _trend_test(), _constrain_est(),
# _constrain_est_identity(), _safe_inverse_spd(), _l_infty()) and
# struc_zero() from scikit-bio 0.7.4 (skbio/stats/composition/
# _ancombc.py), Copyright (c) 2013--, scikit-bio development team, under the
# BSD-3 license given in src/ancombc.c and inst/COPYRIGHTS. scikit-bio's
# code is itself derived from the R package ANCOMBC, Copyright (c) Huang
# Lin, under the Artistic License 2.0 (full text in inst/COPYRIGHTS). This
# file is a modified version of that code (Artistic License 2.0, section
# 4): the bootstraps draw from ecodive's PCG samplers (src/rng_dist.c), not
# NumPy's or R's; the trend test's constrained estimates come from a QP
# solver written for ecodive (src/ancombc.c) in place of SLSQP or
# quadprog; comparisons are labelled with level names; a fit's grouping can
# be replaced; and the results are ecodive's per-feature tables. The
# user-visible differences are listed under "Differences from scikit-bio
# and ANCOMBC" in the documentation below, and the rest in
# plans/ancombc.md.



#' Post-hoc tests for ANCOM-BC and ANCOM-BC2
#'
#' Multigroup tests on a fit from [ancombc()] or [ancombc2()], for a
#' categorical variable with three or more levels (Lin & Peddada 2024), as
#' scikit-bio's `ANCOMBCResult` methods:
#'
#' * `ancombc_global()`: does the feature differ between any of the
#'   levels? A Wald-type test of all the variable's coefficients at once.
#' * `ancombc_pairwise()`: every level against the reference, and every
#'   pair of the other levels, with the mixed directional false discovery
#'   rate (mdFDR) screened by the global test.
#' * `ancombc_dunnett()`: every level against the reference only, with
#'   mdFDR screened by a bootstrapped global test (Dunnett's type).
#' * `ancombc_trend()`: does the feature follow an ordered pattern across
#'   the levels (by default increasing or decreasing)? Constrained
#'   estimates, with a bootstrapped p-value.
#'
#' @param fit  A result of [ancombc()] or [ancombc2()].
#'
#' @param groups  The categorical variable to compare: the name of a
#'        main-effect term of the fit's formula with at least three
#'        levels. Default: the fit's own `groups`. (Any such term can be
#'        given, since the fit keeps every coefficient's covariances.)
#'
#' @param alpha  The significance level of the adjusted p-values. Default:
#'        the fit's.
#'
#' @param p_adjust  The adjustment of the p-values (any method of
#'        [stats::p.adjust()]): across the features for `ancombc_global()`
#'        and `ancombc_trend()`; within each feature, as the family-wise
#'        error control of mdFDR, for `ancombc_pairwise()` and
#'        `ancombc_dunnett()`. Default: the fit's.
#'
#' @param n_bootstraps  The number of bootstrap draws of the null
#'        statistics, per feature. Default: `100`, as scikit-bio and
#'        ANCOMBC. Each p-value of the trend test is a proportion of these
#'        draws, so it needs many more (about 10,000) to be stable.
#'
#' @param patterns  `ancombc_trend()`: the ordered patterns to test, as a
#'        named list of constraint matrices, each with one column per
#'        coefficient of `groups` (in the design's order) and one row per
#'        constraint `C %*% beta >= 0`. Default: `increasing` and
#'        `decreasing`, the coefficients (against the reference level, 0)
#'        rising, or falling, from level to level.
#'
#' @param nodes  `ancombc_trend()`: for each pattern, the coefficient (an
#'        index into the `groups` coefficients) whose contrast with the
#'        reference and with the last level gives the test statistic.
#'        Default: the last, for every pattern.
#'
#' @inheritParams dirmult_ttest
#'
#' @return A data frame of class `c('ecodive_ancombc_<test>', 'data.frame')`.
#'
#'   `ancombc_global()` and `ancombc_trend()`: one row per feature, with
#'   `feature`, (`ancombc_trend()`: `pattern`, the pattern with the largest
#'   statistic,) `w` (the statistic), `p_value`, `p_adjusted` and `reject`.
#'
#'   `ancombc_pairwise()` and `ancombc_dunnett()`: one row per feature and
#'   comparison, with `feature`, `comparison` (`'b - a'`: level b against
#'   level a; with no intercept, a level alone is its own coefficient),
#'   `log_fc`, `se`, `w`, `p_value`, `p_adjusted` and `reject`. As in
#'   scikit-bio and ANCOMBC, `p_value` is set to 1 for the features that
#'   the screening global test didn't find (and wherever it would be 0 or
#'   `NaN`), and `p_adjusted` is the mdFDR-adjusted value.
#'
#'   When the fit was made with `sensitivity` (see [ancombc()]), each test
#'   is repeated on the fit at every pseudocount, and the table gains
#'   `passed_ss` (the decision is the same in every fit; for the trend test,
#'   the global test's decision, as in ANCOMBC) and `reject_robust`
#'   (`reject & passed_ss`), the decisions to use.
#'
#'   Attributes: `method` (the fit's), `groups`, `alpha`, `p_adjust`, for
#'   `ancombc_pairwise()` and `ancombc_dunnett()` `log_base` (`exp(1)`:
#'   `log_fc` is a natural log), for the bootstrapped tests `n_bootstraps`
#'   and `seed`, and with the sensitivity analysis `sensitivity` (its
#'   pseudocounts).
#'
#' @details
#' With \eqn{\hat\beta_f} the feature's \eqn{k} coefficients of `groups`
#' and \eqn{V_f} their covariance matrix (from the fit):
#'
#' * **Global test:** \eqn{W = \hat\beta' V^+ \hat\beta}, referred to
#'   \eqn{\chi^2_k} (`ancombc()`) or \eqn{F_{k, \nu}} (`ancombc2()`, with
#'   the feature's residual degrees of freedom), with the two-sided p-value
#'   \eqn{2 \min(P, 1 - P)}, as in ANCOMBC and scikit-bio.
#' * **Pairwise and Dunnett's tests:** each comparison's \eqn{W =
#'   \hat\delta / \mathrm{se}}, against the normal or t distribution as in
#'   the fit. mdFDR (Guo, Sarkar & Peddada 2010): a global screen at
#'   `alpha`, with Benjamini-Hochberg across features (the global test for
#'   pairwise; for Dunnett's, the largest \eqn{|W|} against `n_bootstraps`
#'   draws of the largest of \eqn{k} null t or normal values); then each
#'   screened feature's comparisons are adjusted with `p_adjust` for a
#'   family of \eqn{m D / R} tests (\eqn{m} comparisons, \eqn{D} features,
#'   \eqn{R} passing the screen).
#' * **Trend test:** for each pattern, the estimate closest to
#'   \eqn{\hat\beta} in the metric \eqn{V^{-1}} that satisfies its
#'   constraints (a small quadratic programme); its statistic is
#'   \eqn{\max(|\beta_{node}|, |\beta_{node} - \beta_k|)}, and \eqn{W} the
#'   largest over the patterns. The p-value is the proportion of
#'   `n_bootstraps` null draws (standard normal coefficients, projected the
#'   same way and scaled by the standard errors) whose statistic exceeds
#'   \eqn{W}; then `p_adjust` across features.
#'
#' Features whose coefficients can't be estimated get `NaN` statistics and a
#' p-value of 1.
#'
#' **Random numbers.** The bootstraps use ecodive's own generator (PCG),
#' with `seed`: one stream per feature, so results are identical for any
#' `cpus` and a feature's results don't depend on the others.
#' `ancombc_global()` and `ancombc_pairwise()` draw nothing.
#'
#' @section Differences from scikit-bio and ANCOMBC:
#'   The statistics and deterministic p-values are scikit-bio's (0.7.4), to
#'   about `1e-10` given the same fit; the bootstrapped p-values agree with
#'   scikit-bio's in distribution (they agree exactly given the same null
#'   draws), not draw for draw, since the random numbers differ. The trend
#'   test's constrained estimates are the exact solutions of its quadratic
#'   programmes, which scikit-bio's SLSQP and ANCOMBC's `quadprog`
#'   approximate. Comparisons are labelled with level names (scikit-bio:
#'   `status[T.severe]_status[T.moderate]`), and `groups` can differ from
#'   the fit's.
#'
#' @references
#' Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
#' microbiomes with covariate adjustments and repeated measures. *Nature
#' Methods*, 21, 83–91.
#'
#' Guo, W., Sarkar, S. K., & Peddada, S. D. (2010). Controlling false
#' discoveries in multidimensional directional decisions, with applications
#' to gene expression data on ordered categories. *Biometrics*, 66, 485–492.
#'
#' @seealso [ancombc()], [structural_zeros()].
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
#'     res <- ancombc_global(fit)
#'     res
#'     plot(res)
#'
#'     res <- ancombc_pairwise(fit)
#'     res[res$reject, ]
#'     plot(res)
#'
#'     ancombc_dunnett(fit)
#'     ancombc_trend(fit, n_bootstraps = 1000)
#'
ancombc_global <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), cpus = n_cpus() ) {

  ph  <- posthoc_setup(environment())
  st  <- posthoc_global_stats(ph)
  res <- posthoc_adjust_global(st$p, ph$p_adjust, ph$alpha)

  out <- data.frame(
    feature    = ph$features,
    w          = st$W,
    p_value    = st$p,
    p_adjusted = res$q,
    reject     = res$reject,
    stringsAsFactors = FALSE )
  out <- posthoc_result(out, ph, 'global')
  posthoc_sensitivity(out, ph, function (f)
    ancombc_global(f, ph$groups, ph$alpha, ph$p_adjust, ph$cpus) )
}



#' @rdname ancombc_global
#' @export
ancombc_pairwise <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), cpus = n_cpus() ) {

  ph <- posthoc_setup(environment())
  k  <- ph$k
  D  <- ph$D
  b  <- ph$beta
  V  <- ph$vcov

  # Each level against the reference, then each pair (j, i), i < j, as
  # itertools.combinations(): b_j - b_i.
  pairs <- if (k > 1L) utils::combn(k, 2L) else matrix(0L, 2L, 0L)
  m     <- k + ncol(pairs)
  beta  <- matrix(NA_real_, D, m)
  var   <- matrix(NA_real_, D, m)
  beta[, seq_len(k)] <- b
  var[, seq_len(k)]  <- ph$var
  for (q in seq_len(ncol(pairs))) {
    i <- pairs[1L, q]
    j <- pairs[2L, q]
    beta[, k + q] <- b[, j] - b[, i]
    var[, k + q]  <- V[, j, j] + V[, i, i] - 2 * V[, j, i]
  }
  labels <- c(ph$labels, if (ncol(pairs))
    paste(ph$levels[pairs[2L, ]], '-', ph$levels[pairs[1L, ]]))

  if (!is.null(ph$valid)) {
    beta[!ph$valid, ] <- NaN
    var[!ph$valid, ]  <- NaN
  }

  # The screen: the global test, with Benjamini-Hochberg.
  st     <- posthoc_global_stats(ph)
  screen <- posthoc_adjust_global(st$p, 'BH', ph$alpha)$reject

  out <- posthoc_mdfdr(ph, beta, var, labels, screen, 'pairwise')
  posthoc_sensitivity(out, ph, function (f)
    ancombc_pairwise(f, ph$groups, ph$alpha, ph$p_adjust, ph$cpus) )
}



#' @rdname ancombc_global
#' @export
ancombc_dunnett <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), n_bootstraps = 100L, seed = 0L,
    cpus = n_cpus() ) {

  ancombc_dunnett_fit(environment(), null_draws = NULL)
}



#' @rdname ancombc_global
#' @export
ancombc_trend <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), patterns = NULL, nodes = NULL,
    n_bootstraps = 100L, seed = 0L, cpus = n_cpus() ) {

  ancombc_trend_fit(environment(), null_draws = NULL)
}



# The bootstrapped tests with given null draws (tests): `null_draws` is a
# features x k x n_bootstraps array (Dunnett: t or normal values; trend:
# standard normal coefficients), used in place of the sampler.
ancombc_dunnett_details <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), n_bootstraps = 100L, seed = 0L,
    cpus = n_cpus(), null_draws = NULL ) {

  ancombc_dunnett_fit(environment(), null_draws)
}

ancombc_trend_details <- function (
    fit, groups = attr(fit, 'groups'), alpha = attr(fit, 'alpha'),
    p_adjust = attr(fit, 'p_adjust'), patterns = NULL, nodes = NULL,
    n_bootstraps = 100L, seed = 0L, cpus = n_cpus(), null_draws = NULL ) {

  ancombc_trend_fit(environment(), null_draws)
}



ancombc_dunnett_fit <- function (env, null_draws) {

  ph <- posthoc_setup(env, boot = TRUE)
  D  <- ph$D
  k  <- ph$k

  beta <- ph$beta
  var  <- ph$var
  if (!is.null(ph$valid)) {
    beta[!ph$valid, ] <- NaN
    var[!ph$valid, ]  <- NaN
  }
  W <- beta / sqrt(pmax(var, 0))

  draws <- posthoc_draws(null_draws, D, k, ph$n_bootstraps)
  exceed <- .Call(
    C_ancombc_dunnett, W, if (!is.null(ph$dof)) as.double(ph$dof),
    ph$n_bootstraps, ph$seed, draws, ph$cpus )

  # The screen (_dunn_global()).
  p_global <- exceed / ph$n_bootstraps
  if (!is.null(ph$valid)) p_global[!ph$valid] <- 1
  screen <- posthoc_adjust_global(p_global, 'BH', ph$alpha)$reject

  out <- posthoc_mdfdr(ph, beta, var, ph$labels, screen, 'dunnett')
  posthoc_sensitivity(out, ph, function (f)
    ancombc_dunnett(f, ph$groups, ph$alpha, ph$p_adjust, ph$n_bootstraps, ph$seed, ph$cpus) )
}



ancombc_trend_fit <- function (env, null_draws) {

  ph <- posthoc_setup(env, boot = TRUE)
  D  <- ph$D
  k  <- ph$k


  #--------------------------------------------------------
  # Patterns and nodes.
  #--------------------------------------------------------
  patterns <- env$patterns
  nodes    <- env$nodes

  if (is.null(patterns)) {
    inc <- diag(k)
    if (k > 1L) inc[cbind(2:k, 1:(k - 1L))] <- -1
    patterns <- list(increasing = inc, decreasing = -inc)
  }
  if (is.matrix(patterns)) patterns <- list(patterns)
  if (!is.list(patterns) || length(patterns) == 0L)
    stop('`patterns` must be a list of constraint matrices.', call. = FALSE)
  for (i in seq_along(patterns)) {
    p <- patterns[[i]]
    if (!is.numeric(p) || !is.matrix(p) || ncol(p) != k || nrow(p) < 1L || !all(is.finite(p)))
      stop(
        'Each of `patterns` must be a numeric matrix with one column per ',
        'coefficient of `groups` (', k, ').', call. = FALSE )
    storage.mode(p) <- 'double'
    patterns[[i]] <- unname(p)
  }
  nms <- names(patterns)
  if (is.null(nms)) nms <- rep('', length(patterns))
  nms[!nzchar(nms)] <- paste0('pattern', which(!nzchar(nms)))
  names(patterns) <- nms

  if (is.null(nodes)) nodes <- rep(k, length(patterns))
  if (!is.null(names(nodes)) && all(names(nodes) %in% nms)) nodes <- nodes[nms]
  if (!is.numeric(nodes) || length(nodes) != length(patterns) || anyNA(nodes) ||
      any(nodes %% 1 != 0) || any(nodes < 1 | nodes > k))
    stop('`nodes` must give each pattern a coefficient between 1 and ', k, '.', call. = FALSE)


  #--------------------------------------------------------
  # Statistics and bootstrap (C).
  #--------------------------------------------------------
  draws <- posthoc_draws(null_draws, D, k, ph$n_bootstraps)
  res   <- .Call(
    C_ancombc_trend, ph$beta, ph$var, ph$vcov,
    if (!is.null(ph$valid)) ph$valid, patterns, as.integer(nodes) - 1L,
    ph$n_bootstraps, ph$seed, draws, ph$cpus )

  pval <- res$counts[seq_len(D)] / ph$n_bootstraps
  if (!is.null(ph$valid)) pval[!ph$valid] <- 1
  adj  <- posthoc_adjust_global(pval, ph$p_adjust, ph$alpha)

  out <- data.frame(
    feature    = ph$features,
    pattern    = nms[res$counts[D + seq_len(D)]],
    w          = res$W,
    p_value    = pval,
    p_adjusted = adj$q,
    reject     = adj$reject,
    stringsAsFactors = FALSE )

  out <- posthoc_result(out, ph, 'trend')

  # The sensitivity analysis uses the global test's decisions, as ANCOMBC's
  # (its trend test is random).
  glob <- function (f) ancombc_global(f, ph$groups, ph$alpha, ph$p_adjust, ph$cpus)
  if (!is.null(ph$sensitivity))
    out <- posthoc_sensitivity(out, ph, glob, main_reject = glob(ph$fit)$reject)
  attr(out, 'patterns') <- patterns
  attr(out, 'nodes')    <- stats::setNames(as.integer(nodes), nms)
  beta <- res$beta
  dimnames(beta) <- list(ph$features, ph$labels)
  attr(out, 'coef') <- beta
  out
}



#------------------------------------------------------------------------------
# Shared pieces.
#------------------------------------------------------------------------------

# Checks the arguments of a post-hoc test (in env) and extracts the fit's
# grouping coefficients. Returns a list: D, k, features, beta (D x k), var
# (D x k), vcov (D x k x k), dof (NULL or D), valid (NULL or logical D),
# levels (the non-reference levels, or every level without an intercept),
# labels (each coefficient's comparison), groups, alpha, p_adjust, method,
# cpus, and with boot = TRUE n_bootstraps and seed.
posthoc_setup <- function (env, boot = FALSE) {

  fit <- env$fit
  if (!inherits(fit, 'ecodive_ancombc') || is.null(attr(fit, 'coef')) ||
      is.null(attr(attr(fit, 'design'), 'xlevels')))
    stop('`fit` must be a result of ancombc() or ancombc2().', call. = FALSE)

  validate_cpus(env)

  groups <- env$groups
  design <- attr(fit, 'design')
  xlev   <- attr(design, 'xlevels')
  terms  <- attr(design, 'term.labels')

  if (is.null(groups))
    stop(
      '`groups` is required: the fit was made without one, so name a ',
      'categorical variable of its formula.', call. = FALSE )
  if (!is.character(groups) || length(groups) != 1L || is.na(groups))
    stop('`groups` must be the name of a variable in the fit\'s formula.', call. = FALSE)
  term <- match(groups, terms)
  if (is.na(term))
    stop('`groups` must be a term in the fit\'s formula: ', name_list(groups), ' is not.', call. = FALSE)
  if (!groups %in% names(xlev))
    stop('`groups` (', name_list(groups), ') must be categorical.', call. = FALSE)
  lev <- xlev[[groups]]
  if (length(lev) < 3L)
    stop('`groups` must have at least three observed levels.', call. = FALSE)
  gi <- which(attr(design, 'assign') == term)
  if (length(gi) < 2L)
    stop(
      '`groups` (', name_list(groups), ') gives fewer than two design columns: ',
      'it must be categorical.', call. = FALSE )

  # Treatment coding (with an intercept): the levels after the first, each
  # against the first. Full coding: every level, as itself.
  if (length(gi) == length(lev) - 1L) {
    levels <- lev[-1L]
    labels <- paste(levels, '-', lev[[1L]])
  }
  else {
    levels <- lev
    labels <- lev
  }

  alpha <- env$alpha
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !(alpha > 0 && alpha < 1))
    stop('`alpha` must be a single number strictly between 0 and 1.', call. = FALSE)

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)

  coef <- attr(fit, 'coef')
  est  <- attr(fit, 'estimable')
  valid <- NULL
  if (!is.null(est)) {
    valid <- apply(est[, gi, drop = FALSE], 1L, all)
    if (all(valid)) valid <- NULL
  }

  out <- list(
    D        = nrow(coef),
    k        = length(gi),
    features = rownames(coef),
    beta     = unname(coef[, gi, drop = FALSE]),
    var      = unname(attr(fit, 'var')[, gi, drop = FALSE]),
    vcov     = unname(attr(fit, 'vcov')[, gi, gi, drop = FALSE]),
    dof      = if (!is.null(attr(fit, 'dof'))) unname(attr(fit, 'dof')),
    valid    = unname(valid),
    levels   = levels,
    labels   = labels,
    groups   = groups,
    alpha    = as.double(alpha),
    p_adjust = p_adjust,
    method   = attr(fit, 'method'),
    cpus     = env$cpus,
    fit      = fit,
    sensitivity = attr(fit, 'sensitivity') )

  if (boot) {
    n_boot <- env$n_bootstraps
    if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
        n_boot < 1 || n_boot %% 1 != 0 || n_boot > .Machine$integer.max)
      stop('`n_bootstraps` must be a single positive integer.', call. = FALSE)
    validate_seed(env)
    out$n_bootstraps <- as.integer(n_boot)
    out$seed         <- env$seed
  }

  out
}


# The global test's statistics and p-values (_global_stats()).
posthoc_global_stats <- function (ph) {
  W <- .Call(C_ancombc_global, ph$beta, ph$vcov, if (!is.null(ph$valid)) ph$valid, ph$cpus)
  k <- ph$k
  if (is.null(ph$dof)) {
    lo <- stats::pchisq(W, k)
    hi <- stats::pchisq(W, k, lower.tail = FALSE)
  }
  else {
    lo <- suppressWarnings(stats::pf(W, k, ph$dof))
    hi <- suppressWarnings(stats::pf(W, k, ph$dof, lower.tail = FALSE))
  }
  p <- 2 * pmin(lo, hi)
  p[is.na(p)] <- NaN
  if (!is.null(ph$valid)) p[!ph$valid] <- 1
  list(W = W, p = p)
}


# Adjustment across features (p.adjust() over the p-values that aren't NaN),
# NaN as 1, and the decisions.
posthoc_adjust_global <- function (p, p_adjust, alpha) {
  ok <- !is.na(p)
  q  <- p
  q[ok] <- stats::p.adjust(p[ok], p_adjust)
  q[!ok] <- 1
  list(q = q, reject = q <= alpha)
}


# The mdFDR step shared by the pairwise and Dunnett's tests
# (_mdfdr_pairwise(), _mdfdr_dunnett()), and their table.
posthoc_mdfdr <- function (ph, beta, var, labels, screen, test) {

  D  <- ph$D
  m  <- ncol(beta)
  se <- sqrt(pmax(var, 0))
  W  <- beta / se

  p <- if (is.null(ph$dof)) 2 * stats::pnorm(abs(W), lower.tail = FALSE)
       else suppressWarnings(2 * stats::pt(abs(W), ph$dof, lower.tail = FALSE))
  p <- p * screen
  p[!is.na(p) & p == 0] <- 1
  p[is.na(p)] <- 1

  R <- sum(screen)
  q <- matrix(1, D, m)
  if (R > 0) {
    n <- m * D / R
    for (f in seq_len(D)) q[f, ] <- stats::p.adjust(p[f, ], ph$p_adjust, n = n)
  }

  byrow <- function (x) as.vector(t(x))
  out <- data.frame(
    feature    = rep(ph$features, each = m),
    comparison = rep(labels, times = D),
    log_fc     = byrow(beta),
    se         = byrow(se),
    w          = byrow(W),
    p_value    = byrow(p),
    p_adjusted = byrow(q),
    reject     = byrow(q <= ph$alpha),
    stringsAsFactors = FALSE )

  posthoc_result(out, ph, test)
}


# The sensitivity analysis of a post-hoc result, when the fit has one: the
# same test (`rerun`, a function of a fit) on the fit at each pseudocount;
# `passed_ss` is TRUE where its decision (or, for the trend test, the global
# test's, `main_reject`) is the same in every fit, and `reject_robust` is
# reject & passed_ss.
posthoc_sensitivity <- function (out, ph, rerun, main_reject = out$reject) {
  sens <- ph$sensitivity
  if (is.null(sens)) return (out)
  passed <- rep(TRUE, nrow(out))
  for (f in sens$fits) passed <- passed & (rerun(f)$reject == main_reject)
  out$passed_ss     <- passed
  out$reject_robust <- out$reject & passed
  attr(out, 'sensitivity') <- sens$pseudocounts
  out
}


posthoc_result <- function (out, ph, test) {
  attr(out, 'method')   <- ph$method
  attr(out, 'groups')   <- ph$groups
  attr(out, 'alpha')    <- ph$alpha
  attr(out, 'p_adjust') <- ph$p_adjust
  if ('log_fc' %in% names(out)) attr(out, 'log_base') <- exp(1)
  if (!is.null(ph$n_bootstraps)) {
    attr(out, 'n_bootstraps') <- ph$n_bootstraps
    attr(out, 'seed')         <- ph$seed
  }
  class(out) <- c(paste0('ecodive_ancombc_', test), 'data.frame')
  out
}


posthoc_draws <- function (draws, D, k, B) {
  if (is.null(draws)) return (NULL)
  stopifnot(is.numeric(draws), identical(as.integer(dim(draws)), as.integer(c(D, k, B))))
  storage.mode(draws) <- 'double'
  draws
}



#' Structural zeros
#'
#' Finds the features that are absent, or nearly absent, from a group of
#' samples (Kaul et al. 2017), as scikit-bio's `struc_zero()` and ANCOMBC's
#' `struc_zero` option. A feature that is a structural zero in a group is
#' differentially (less) abundant there, whatever a test of its abundance
#' says, and should be left out of analyses of that group; one that is a
#' structural zero in every group should be removed.
#'
#' @inheritParams ancom
#'
#' @param counts  A numeric matrix of abundances (samples \eqn{\times}
#'        features), or any of the objects listed under **Input Types**:
#'        the table **before** any pseudocount or zero replacement.
#'
#' @param neg_lb  `FALSE` (the default): a feature is a structural zero in
#'        a group when it is zero in every sample of that group. `TRUE`:
#'        when the lower bound of the 95% interval of its proportion of
#'        non-zero samples, \eqn{\hat p - 1.96 \sqrt{\hat p (1 - \hat p) /
#'        n}}, is at most 0; recommended for large groups.
#'
#' @return A logical matrix, features \eqn{\times} groups, `TRUE` where the
#'   feature is a structural zero in the group. Groups are in sorted order
#'   (`sort(method = 'radix')`, as scikit-bio's), or in a factor's level
#'   order.
#'
#' @references
#' Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
#' microbiome data in the presence of excess zeros. *Frontiers in
#' Microbiology*, 8, 2114.
#'
#' @inheritSection documentation Input Types
#'
#' @seealso [ancombc()], [ancombc_global()].
#'
#' @export
#' @examples
#'     x <- rbind(
#'       c( 7, 1,  0, 11,  3,  1), c( 1, 1, 0, 13, 13,  0),
#'       c(11, 5,  0,  1,  4,  1), c( 2, 2, 0, 16,  4,  0),
#'       c( 0, 1,  0,  0,  6,  0), c(14, 8, 7,  9,  0,  5),
#'       c( 0, 7,  4,  1,  0, 26), c( 8, 1, 4, 28,  0, 10),
#'       c( 2, 2,  2,  4,  0,  5), c( 6, 4, 10,  1,  0,  9) )
#'     colnames(x) <- paste0('f', 0:5)
#'     grp <- rep(c('treatment', 'placebo'), each = 5)
#'     structural_zeros(x, grp)
#'     structural_zeros(x, grp, neg_lb = TRUE)
#'
structural_zeros <- function (counts, groups, data = NULL, margin = 1L, neg_lb = FALSE) {

  if (missing(groups)) groups <- NULL
  validate_counts()
  validate_margin()

  if (!isTRUE(neg_lb) && !isFALSE(neg_lb))
    stop('`neg_lb` must be TRUE or FALSE.', call. = FALSE)

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  pkg  <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]]) else as.character(seq_len(D))
  samples  <- if (length(dn[[margin]])) as.character(dn[[margin]])

  des <- test_group_values(groups, data, n, samples, ref = '`counts`')
  val <- des$groups
  lev <- if (is.factor(val)) levels(val) else sort(unique(val), method = 'radix')
  gi  <- if (is.factor(val)) as.integer(val) else match(val, lev)

  nz   <- .Call(C_ancombc_nonzero, counts, margin, gi - 1L, length(lev))
  size <- rep(tabulate(gi, length(lev)), each = D)
  p    <- nz / size
  if (neg_lb) p <- p - 1.96 * (p * (1 - p) / size)^0.5

  out <- matrix(p <= 0, D, length(lev), dimnames = list(features, as.character(lev)))
  out
}
