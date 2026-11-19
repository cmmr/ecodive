# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# ancom() follows ancom(), _log_compare(), _adjust_pvalues() and
# _check_grouping() from scikit-bio 0.7.4 (skbio/stats/composition/_ancom.py
# and _utils.py), Copyright (c) 2013--, scikit-bio development team, under
# the BSD-3 license given in src/ancom.c and inst/COPYRIGHTS.



#' Analysis of Composition of Microbiomes (ANCOM)
#'
#' Finds features whose abundance differs between groups without assuming
#' that the samples' totals mean anything (Mandal et al. 2015). Every log
#' ratio between two features is tested for a difference between the
#' groups; a feature whose abundance really changed makes most of its log
#' ratios differ. Each feature's \eqn{W} is the number of its log ratios that
#' differ significantly, and the features with a large \eqn{W} are declared
#' differentially abundant.
#'
#' @inheritParams vlr
#' @inheritParams simper
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of [replace_zeros()]
#'        are all accepted: each sample's scale cancels out of the log
#'        ratios.
#'
#' @param alpha  The significance level of each log ratio's (adjusted)
#'        p-value. Strictly between 0 and 1. Default: `0.05`.
#'
#' @param tau  How much the proportion of features above a cutoff may change
#'        between neighbouring cutoffs for the cutoff to count as stable
#'        (see **Details**); small values are conservative. Strictly between
#'        0 and 1. Default: `0.02`.
#'
#' @param theta  The smallest largest \eqn{W}, as a proportion of the number
#'        of features, for anything to be declared differentially abundant.
#'        Strictly between 0 and 1. Default: `0.1`.
#'
#' @param p_adjust  How to adjust each feature's p-values for testing its
#'        log ratios with every other feature: `'holm'` (the default),
#'        `'bonferroni'`, `'BH'` (or `'fdr'`), `'BY'`, or `'none'`, as in
#'        [stats::p.adjust()]. Each feature's p-values are one family (see
#'        **Details**).
#'
#' @param sig_test  The test of each log ratio across the groups:
#'        `'anova'` (the default; scikit-bio's `"f_oneway"`), a one-way
#'        ANOVA, which is Student's t-test for two groups; or `'kruskal'`,
#'        the Kruskal-Wallis rank sum test, which the original ANCOM paper
#'        used.
#'
#' @return A data frame of class `c('ecodive_ancom', 'data.frame')`, one row
#'   per feature, in the order of `counts`, with columns:
#'
#'   * `feature`: the feature's name (its index, as text, when `counts` has
#'     no feature names).
#'   * `W`: the number of other features whose log ratio with this one
#'     differs significantly between the groups (integer, `0` to \eqn{D -
#'     1}).
#'   * `reject`: `TRUE` for the features declared differentially abundant
#'     (scikit-bio's `Signif`).
#'
#'   Attributes: `groups` (the groups, in order), `sig_test`, `p_adjust`,
#'   `alpha`, `tau`, `theta`, `cutoff` (the chosen \eqn{\nu}, as a proportion
#'   of the number of features, or `NA` when nothing could be rejected
#'   because the largest \eqn{W} was below `theta`), and `pseudocount` (the
#'   pseudocount used; `0` for none).
#'
#'   ANCOM gives **no effect size and no p-value per feature**, only \eqn{W}
#'   and a decision. For effect sizes and their uncertainty, see
#'   `dirmult_ttest()`; for bias-corrected estimates of each feature's log
#'   fold change, with p-values, see [ancombc()] and [ancombc2()].
#'
#' @details
#' With \eqn{L = \log(x + c)} (natural logs, \eqn{c} the pseudocount) and
#' \eqn{D} features:
#'
#' 1. For every pair of features \eqn{i < j}, the log ratio \eqn{L_i - L_j}
#'    of each sample is tested across the groups with `sig_test`, giving
#'    \eqn{p_{ij} = p_{ji}}.
#' 2. Each feature's \eqn{D - 1} p-values are adjusted for multiple testing
#'    as one family (`p_adjust`).
#' 3. \eqn{W_i} is the number of features \eqn{j} with adjusted
#'    \eqn{p_{ij} < \alpha}.
#' 4. If the largest \eqn{W}, as a proportion of \eqn{D}, is below `theta`,
#'    nothing is rejected. Otherwise, with \eqn{c_0 = \max(W)/D}, five
#'    cutoffs \eqn{c_0 - 0.05, \ldots, c_0 - 0.25} are tried: the
#'    proportions of features with \eqn{W > D \cdot \mathrm{cutoff}}, and the
#'    changes between neighbouring proportions, decide which cutoff
#'    \eqn{\nu} is stable (the first change below `tau` picks \eqn{c_0 -
#'    0.10}, and so on, as in scikit-bio). Features with \eqn{W \ge \nu D}
#'    are rejected.
#'
#' **Matching scikit-bio.** The results (\eqn{W} and the decisions) are
#' those of scikit-bio's `ancom()` (0.7.4), whose rules are followed
#' exactly, including two of its quirks:
#'
#' * Each feature's family of p-values for the adjustment includes a 0 for
#'   the feature against itself. For Holm this changes nothing; for
#'   Bonferroni, BH and BY it makes the adjustment slightly more stringent
#'   (\eqn{D} tests instead of \eqn{D - 1}). In R terms, feature \eqn{i}'s
#'   adjusted p-values are `p.adjust(c(0, p[i, -i]), p_adjust)[-1]`.
#' * Its cutoffs are \eqn{c_0 - } `numpy.linspace(0.05, 0.25, 5)`, whose
#'   third value is one rounding step above `0.15`.
#'
#' A log ratio that is constant within every group but differs between
#' groups has \eqn{F = \infty} and \eqn{p = 0} with `sig_test = 'anova'`, as
#' with SciPy's `f_oneway()`; one that is the same in every sample has no
#' p-value (`NaN`, as in SciPy), and is left out of the feature's family and
#' never counted in \eqn{W}. The Kruskal-Wallis test uses average ranks for
#' ties and SciPy's tie correction, with ties by exact equality, as SciPy's
#' `kruskal()`; a log ratio with every value tied has no p-value (SciPy
#' sometimes returns an infinite statistic there, by rounding, and so a
#' p-value of 0). The F and \eqn{\chi^2} p-values are R's [stats::pf()] and
#' [stats::pchisq()], which agree with SciPy's to about `1e-13`. The F
#' statistics come from two-pass sums of squares; SciPy's one-pass ones
#' differ from them by up to about `1e-12` (relative), from rounding.
#'
#' Groups are sorted as scikit-bio sorts them (`sort(method = 'radix')`),
#' unless `groups` is a factor, whose levels' order is kept. The order
#' changes nothing but the `groups` attribute.
#'
#' scikit-bio's `percentiles` table (each group's percentiles of each
#' feature's raw abundance) is not returned: use [stats::quantile()].
#' ANCOM-II's covariate adjustment and structural-zero detection (Kaul et
#' al. 2017) are not offered.
#'
#' **Computation.** The \eqn{D(D - 1)/2} tests run in parallel C, each
#' straight from its log ratios (no \eqn{n \times D} log table is made for a
#' sparse `counts`), and the p-values take \eqn{4 D^2} bytes: 5,000
#' features need 100 MB. Results are identical for any number of `cpus`.
#'
#' @section Zeros:
#'   Logs of zero are undefined, so zeros need replacing first. With
#'   `pseudocount = NULL` (the default), no pseudocount is used when
#'   `counts` has no zeros; otherwise half the smallest non-zero value is
#'   added to every value, with a warning. The recommended input is a table
#'   with its zeros replaced by [replace_zeros()], which has none and is
#'   used as it is. (scikit-bio's `ancom()` refuses zeros; the original
#'   ANCOM paper added a pseudocount of 0.001.)
#'
#' @section Same names in other packages:
#'   scikit-bio (Python) has `ancom()`, which this matches. R's ANCOMBC
#'   package has `ancom()` too, an implementation of ANCOM-II with its own
#'   zero handling and defaults; it is not what this function reproduces.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., &
#' Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel
#' method for studying microbial composition. *Microbial Ecology in Health
#' and Disease*, 26, 27663.
#'
#' Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
#' microbiome data in the presence of excess zeros. *Frontiers in
#' Microbiology*, 8, 2114.
#'
#' @seealso [replace_zeros()], [vlr()], [simper()].
#'
#' @export
#' @examples
#'     # scikit-bio's example: b2 is more abundant under placebo.
#'     x <- rbind(
#'       c(12, 11, 10, 10, 10, 10, 10),
#'       c( 9, 11, 12, 10, 10, 10, 10),
#'       c( 1, 11, 10, 11, 10,  5,  9),
#'       c(22, 21,  9, 10, 10, 10, 10),
#'       c(20, 22, 10, 10, 13, 10, 10),
#'       c(23, 21, 14, 10, 10, 10, 10) )
#'     colnames(x) <- paste0('b', 1:7)
#'     grp <- rep(c('treatment', 'placebo'), each = 3)
#'     res <- ancom(x, grp)
#'     res
#'     plot(res)
#'
#'     # A table with zeros: replace them first.
#'     grp <- c('oral', 'oral', 'other', 'other')
#'     ancom(replace_zeros(ex_counts, 'mult'), grp, sig_test = 'kruskal')
#'
ancom <- function (
    counts, groups, data = NULL, margin = 1L, pseudocount = NULL,
    alpha = 0.05, tau = 0.02, theta = 0.1, p_adjust = 'holm',
    sig_test = 'anova', cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  ancom_fit(environment(), details = FALSE)
}



# ancom()'s statistics and p-values, for the tests: the same arguments, and
# returns list(result, statistic (packed triangle over the features, row by
# row as a `dist`), p_value (D x D, symmetric, NaN on the diagonal),
# p_adjusted (D x D, row i is feature i's family; 1 on the diagonal)).
ancom_details <- function (
    counts, groups, data = NULL, margin = 1L, pseudocount = NULL,
    alpha = 0.05, tau = 0.02, theta = 0.1, p_adjust = 'holm',
    sig_test = 'anova', cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  ancom_fit(environment(), details = TRUE)
}



ancom_fit <- function (env, details) {

  info <- log_ratio_input('ANCOM', env = env, dense_out = FALSE)
  counts      <- env$counts
  margin      <- env$margin
  pseudocount <- env$pseudocount
  cpus        <- env$cpus

  for (arg in c('alpha', 'tau', 'theta')) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !(v > 0 && v < 1))
      stop('`', arg, '` must be a single number strictly between 0 and 1.', call. = FALSE)
  }
  alpha <- as.double(env$alpha)
  tau   <- as.double(env$tau)
  theta <- as.double(env$theta)

  p_adjust <- ancom_p_adjust(env$p_adjust)
  sig_test <- ancom_sig_test(env$sig_test)


  #--------------------------------------------------------
  # Each sample's group, in counts' order.
  #--------------------------------------------------------
  des <- test_group_values(env$groups, env$data, info$n, info$samples, ref = '`counts`')
  val <- des$groups
  lev <- if (is.factor(val)) levels(val) else sort(unique(val), method = 'radix')
  gi  <- if (is.factor(val)) as.integer(val) else match(val, lev)

  if (length(lev) < 2L)
    stop('`groups` must have at least two different values.', call. = FALSE)
  if (length(lev) == info$n)
    stop(
      'Every sample is in a group of its own: ANCOM needs replicates ',
      'within at least one group.', call. = FALSE )


  D   <- info$D
  res <- .Call(
    C_ancom, counts, margin, gi - 1L, length(lev),
    switch(sig_test, anova = 1L, kruskal = 2L),
    switch(p_adjust, none = 0L, bonferroni = 1L, holm = 2L, BH = 3L, BY = 4L),
    alpha, pseudocount, details, cpus )

  W <- if (details) res$W else res


  cut    <- ancom_cutoff(W, D, tau, theta)
  reject <- cut$reject
  nu     <- cut$nu

  features <- if (is.null(info$features)) as.character(seq_len(D)) else info$features

  out <- data.frame(
    feature = features, W = W, reject = reject,
    stringsAsFactors = FALSE )

  attr(out, 'groups')      <- if (is.factor(val)) lev else as.character(lev)
  attr(out, 'sig_test')    <- sig_test
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'alpha')       <- alpha
  attr(out, 'tau')         <- tau
  attr(out, 'theta')       <- theta
  attr(out, 'cutoff')      <- nu
  attr(out, 'pseudocount') <- pseudocount
  class(out) <- c('ecodive_ancom', 'data.frame')

  if (!details) return (out)

  p <- matrix(NaN, D, D, dimnames = list(features, features))
  p[lower.tri(p)] <- res$p_value      # as a `dist` stores it
  p[upper.tri(p)] <- t(p)[upper.tri(p)]

  padj <- t(res$p_adjusted)
  dimnames(padj) <- list(features, features)

  list(result = out, statistic = res$statistic, p_value = p, p_adjusted = padj)
}



# scikit-bio's choice of cutoff from the W of D features. Returns
# list(reject, nu); nu is NA when the largest W is below theta * D.
ancom_cutoff <- function (W, D, tau, theta) {

  c_start <- max(W) / D
  if (c_start < theta)
    return (list(reject = logical(length(W)), nu = NA_real_))

  # numpy.linspace(0.05, 0.25, 5): its third value is not 0.15.
  steps    <- (0:4) * ((0.25 - 0.05) / 4) + 0.05
  steps[5] <- 0.25
  cutoff   <- c_start - steps
  prop_cut <- vapply(cutoff, function (cut) sum(W > D * cut) / D, numeric(1))
  dels     <- abs(prop_cut - c(prop_cut[-1], prop_cut[1]))
  dels[5]  <- 0

  nu <- if      (dels[1] <  tau && dels[2] < tau && dels[3] < tau) cutoff[2]
        else if (dels[1] >= tau && dels[2] < tau && dels[3] < tau) cutoff[3]
        else if (dels[2] >= tau && dels[3] < tau && dels[4] < tau) cutoff[4]
        else                                                      cutoff[5]

  list(reject = W >= nu * D, nu = nu)
}



# ancom()'s `p_adjust`: the methods with ANCOM's own adjustment in C.
ancom_p_adjust <- function (p_adjust) {
  ok <- c('holm', 'bonferroni', 'BH', 'BY', 'none')
  m  <- tryCatch(pairwise_method(TRUE, p_adjust, TRUE), error = function (e) NA)
  if (identical(m, 'fdr')) m <- 'BH'
  if (!isTRUE(m %in% ok))
    stop(
      '`p_adjust` must be one of ', paste0("'", ok, "'", collapse = ', '),
      " (or 'fdr', the same as 'BH').", call. = FALSE )
  m
}



# ancom()'s `sig_test`: 'anova' (or scikit-bio's 'f_oneway') or 'kruskal'.
ancom_sig_test <- function (sig_test) {
  msg <- "`sig_test` must be 'anova' or 'kruskal'."
  if (!is.character(sig_test) || length(sig_test) != 1L || is.na(sig_test))
    stop(msg, call. = FALSE)
  s <- tolower(sig_test)
  if (s == 'f_oneway') return ('anova')
  m <- pmatch(s, c('anova', 'kruskal'))
  if (is.na(m)) stop(msg, call. = FALSE)
  c('anova', 'kruskal')[[m]]
}
