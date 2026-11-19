# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# dirmult_ttest() follows dirmult_ttest() from scikit-bio 0.7.4
# (skbio/stats/composition/_dirmult.py) and _check_trt_ref_groups() and
# _adjust_pvalues() (_utils.py), Copyright (c) 2013--, scikit-bio
# development team, under the BSD-3 license given in src/dirmult.c and
# inst/COPYRIGHTS.



#' Dirichlet-multinomial Welch t-tests (ALDEx2-style)
#'
#' Finds features whose abundance differs between two sets of samples,
#' carrying the uncertainty of low counts and zeros through the test
#' (Fernandes et al. 2013, 2014). Each sample's proportions are drawn many
#' times from their Dirichlet posterior given its counts; every draw is
#' clr-transformed and each feature is compared with Welch's t-test; and the
#' results are averaged over the draws. This is the ALDEx2 approach, as
#' scikit-bio implements it.
#'
#' @inheritParams ancom
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Raw counts are recommended: the Dirichlet posterior's
#'        spread depends on the counts' size, so proportions or rescaled
#'        values give wrong uncertainties.
#'
#' @param contrast  Which samples to compare, as `c(treatment, reference)`,
#'        two values of `groups`. A positive `log2_fc` means more abundant in
#'        the treatment group. Use `NA` as the reference to compare the
#'        treatment group with all other samples pooled; with two named
#'        groups, samples in other groups are left out. With `NULL` (the
#'        default), the treatment is the **first** group (in sorted order,
#'        or a factor's first level) and the reference is all other samples:
#'        with two groups, the first is compared with the second, the
#'        reverse of R's usual "second level versus first".
#'
#' @param pseudocount  The Dirichlet prior: a positive number added to every
#'        value, zero or not, before drawing. Default: `0.5`, as scikit-bio
#'        and ALDEx2. This is not the zero patch that other functions'
#'        `pseudocount` is (no `NULL` rule): every value gets it.
#'
#' @param n_draws  The number of Monte Carlo draws from the Dirichlet
#'        posterior (scikit-bio's `draws`). Default: `128`.
#'
#' @param alpha  The significance level for `reject`, applied to the
#'        adjusted p-values. Strictly between 0 and 1. Default: `0.05`, as
#'        scikit-bio (which fixes it).
#'
#' @param conf_level  The confidence level of each draw's interval for the
#'        difference. Strictly between 0 and 1. Default: `0.95`, as
#'        scikit-bio (which fixes it).
#'
#' @param p_adjust  How to adjust the averaged p-values for testing every
#'        feature: one of [stats::p.adjust.methods] (`'holm'`, the default,
#'        `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`, `'fdr'`,
#'        or `'none'`).
#'
#' @param seed  An integer seed for the random draws. The same `seed` gives
#'        the same results, for any number of `cpus`.
#'
#' @return A data frame of class `c('ecodive_dirmult_ttest', 'data.frame')`,
#'   one row per feature, in the order of `counts`, with columns:
#'
#'   * `feature`: the feature's name (its index, as text, when `counts` has
#'     no feature names).
#'   * `t`: Welch's t statistic, averaged over the draws (scikit-bio's
#'     `T-statistic`).
#'   * `log2_fc`: the difference between the groups' mean clr values,
#'     treatment minus reference, in log2 units, averaged over the draws
#'     (`Log2(FC)`).
#'   * `ci_lower`, `ci_upper`: the lowest lower bound and the highest upper
#'     bound of the draws' intervals for that difference, in log2 units
#'     (`CI(2.5)`, `CI(97.5)` for the default `conf_level`).
#'   * `p_value`: Welch's two-sided p-value, averaged over the draws
#'     (`pvalue`).
#'   * `p_adjusted`: `p_value` adjusted across the features with `p_adjust`
#'     (`qvalue`).
#'   * `reject`: `TRUE` when `p_adjusted <= alpha` and the interval
#'     (`ci_lower`, `ci_upper`) excludes 0 (`Signif`).
#'
#'   Attributes: `contrast` (the treatment and reference groups; the
#'   reference is `NA` when all other samples were pooled), `n` (the number
#'   of samples on each side), `n_draws`, `pseudocount`, `alpha`,
#'   `conf_level`, `p_adjust` and `seed`.
#'
#' @details
#' With \eqn{x_{sf}} the abundance of feature \eqn{f} in sample \eqn{s} and
#' \eqn{c} the `pseudocount`, each of the `n_draws` draws:
#'
#' 1. Draws \eqn{G_{sf} \sim \mathrm{Gamma}(x_{sf} + c, 1)}, independently.
#'    Each sample's \eqn{G_{s\cdot}}, divided by its sum, is a draw from the
#'    Dirichlet posterior of its proportions; the clr makes the division
#'    unnecessary.
#' 2. Takes the clr: \eqn{z_{sf} = \log G_{sf} - \frac{1}{D}\sum_g \log
#'    G_{sg}}.
#' 3. For each feature, compares the treatment's \eqn{z} with the
#'    reference's by Welch's t-test: the difference of means \eqn{\Delta},
#'    its standard error \eqn{se = \sqrt{s_1^2/n_1 + s_2^2/n_2}}, the
#'    Welch-Satterthwaite degrees of freedom \eqn{\nu}, \eqn{t = \Delta /
#'    se}, the two-sided p-value, and the interval \eqn{\Delta \pm
#'    t_{\nu}\, se}, with \eqn{t_\nu} the \eqn{(1 + }`conf_level`\eqn{)/2}
#'    quantile of Student's t distribution.
#'
#' `t`, `log2_fc` and `p_value` are the averages over the draws (`log2_fc`
#' divided by \eqn{\log 2}). The interval is **not** a quantile of the
#' draws: it runs from the lowest lower bound to the highest upper bound of
#' any draw, as scikit-bio's, so a feature is only rejected when no draw's
#' interval includes 0. This conservative envelope guards against false
#' positives among rare features, whose draws vary most.
#'
#' **The reference frame is the geometric mean.** A clr difference is a log
#' fold change only if the groups' geometric means are the same, i.e. if the
#' typical feature doesn't change; when the total load differs between the
#' groups, every `log2_fc` is shifted by the same amount (the order of the
#' features stays meaningful). [ancom()] tests log ratios without that
#' assumption, but gives no effect sizes.
#'
#' A side with a single sample has no variance: its features' `t`,
#' `p_value`, intervals and `p_adjusted` are `NaN`, with a warning, and none
#' are rejected (as scikit-bio's, which gives them silently).
#'
#' **Matching scikit-bio.** Given the same Gamma draws, the results are
#' those of scikit-bio's `dirmult_ttest()` (0.7.4) to about `1e-13`: the
#' tests feed both the same draws. The draws themselves come from ecodive's
#' own Gamma sampler (Marsaglia & Tsang 2000, on PCG), not NumPy's, so for
#' the same data the results differ from scikit-bio's by Monte Carlo error
#' only. The t distribution is R's [stats::pt()] and [stats::qt()] (SciPy's
#' `stdtr()` and `stdtrit()` in scikit-bio), and `p_adjust` is
#' [stats::p.adjust()], with `NaN` p-values left out of the family, as
#' scikit-bio. Differences from scikit-bio: `alpha` and `conf_level` are
#' arguments (scikit-bio fixes them at 0.05 and 0.95); `pseudocount` must be
#' positive; and `seed` defaults to `0`, not to a random seed.
#'
#' **ALDEx2** (Bioconductor) follows the same idea, but differs in its
#' details: it reports both Welch's and Wilcoxon's tests, averages
#' Benjamini-Hochberg adjusted p-values over the draws, and gives its own
#' effect size. Its results agree with this function's in direction and
#' rough size, not exactly.
#'
#' **Computation.** The draws run in parallel C, one draw at a time, with
#' memory for one draw (8 bytes per value of the samples in the contrast).
#' Every sample of every draw has its own random stream, so results are
#' identical for any number of `cpus`. 500 samples \eqn{\times} 5,000
#' features, 128 draws, take a few seconds on 4 cores.
#'
#' @section Same names in other packages:
#'   scikit-bio (Python) has `dirmult_ttest()`, which this follows. ALDEx2's
#'   `aldex()` with `test = "t"` is the original method (see above).
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Fernandes, A. D., Macklaim, J. M., Linn, T. G., Reid, G., & Gloor, G. B.
#' (2013). ANOVA-like differential expression (ALDEx) analysis for mixed
#' population RNA-Seq. *PLoS ONE*, 8, e67019.
#'
#' Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A.,
#' Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of
#' high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA
#' gene sequencing and selective growth experiments by compositional data
#' analysis. *Microbiome*, 2, 15.
#'
#' Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating
#' gamma variables. *ACM Transactions on Mathematical Software*, 26,
#' 363-372.
#'
#' @seealso [ancom()], [replace_zeros()].
#'
#' @export
#' @examples
#'     # scikit-bio's example: b1 and b2 are less abundant under treatment.
#'     x <- rbind(
#'       c( 20, 110, 100, 101, 100, 103, 104),
#'       c( 33, 110, 120, 100, 101, 100, 102),
#'       c( 12, 110, 100, 110, 100,  50,  90),
#'       c(202, 201,   9,  10,  10,  11,  11),
#'       c(200, 202,  10,  10,  13,  10,  10),
#'       c(203, 201,  14,  10,  10,  13,  12) )
#'     colnames(x) <- paste0('b', 1:7)
#'     grp <- rep(c('treatment', 'placebo'), each = 3)
#'     dirmult_ttest(x, grp, contrast = c('treatment', 'placebo'))
#'
#'     # The default contrast: the first group ('placebo', in sorted order)
#'     # against the rest.
#'     dirmult_ttest(x, grp, n_draws = 32)$log2_fc
#'
dirmult_ttest <- function (
    counts, groups, data = NULL, contrast = NULL, margin = 1L,
    pseudocount = 0.5, n_draws = 128L, alpha = 0.05, conf_level = 0.95,
    p_adjust = 'holm', seed = 0L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  dirmult_ttest_fit(environment())
}



# dirmult_ttest() for the tests: the same arguments, plus `gamma_draws` (an
# n x D x n_draws array of Gamma draws to use instead of sampling, samples
# and features as in `counts`' default orientation); returns list(result,
# draws), `draws` holding each draw's diff, se, dof, t, p_tail (one-sided),
# lower and upper (n_draws x D matrices, natural logs).
dirmult_ttest_details <- function (
    counts, groups, data = NULL, contrast = NULL, margin = 1L,
    pseudocount = 0.5, n_draws = 128L, alpha = 0.05, conf_level = 0.95,
    p_adjust = 'holm', seed = 0L, cpus = n_cpus(), gamma_draws = NULL ) {

  if (missing(groups)) groups <- NULL
  dirmult_ttest_fit(environment(), details = TRUE, gamma_draws = gamma_draws)
}



dirmult_ttest_fit <- function (env, details = FALSE, gamma_draws = NULL) {

  info <- dirmult_input(env)
  counts <- env$counts
  margin <- env$margin
  cpus   <- env$cpus
  validate_seed(env)
  seed <- env$seed

  pc <- env$pseudocount
  if (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc <= 0)
    stop('`pseudocount` must be a single positive number.', call. = FALSE)
  pc <- as.double(pc)

  n_draws <- env$n_draws
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1 || n_draws %% 1 != 0 || n_draws > .Machine$integer.max)
    stop('`n_draws` must be a single positive integer.', call. = FALSE)
  n_draws <- as.integer(n_draws)

  for (arg in c('alpha', 'conf_level')) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !(v > 0 && v < 1))
      stop('`', arg, '` must be a single number strictly between 0 and 1.', call. = FALSE)
  }
  alpha      <- as.double(env$alpha)
  conf_level <- as.double(env$conf_level)

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)


  #--------------------------------------------------------
  # The two sides of the contrast.
  #--------------------------------------------------------
  des <- test_group_values(env$groups, env$data, info$n, info$samples, ref = '`counts`')
  ct  <- dirmult_contrast(des$groups, env$contrast)
  n1  <- sum(ct$side == 1L)
  n2  <- sum(ct$side == 2L)

  if (n1 < 2L || n2 < 2L)
    warning(
      'The ', if (n1 < 2L) 'treatment' else 'reference', ' side of `contrast` ',
      'has a single sample, so no variance: `t`, the p-values and the ',
      'intervals are NaN.', call. = FALSE )


  #--------------------------------------------------------
  # Draws given by the tests.
  #--------------------------------------------------------
  if (!is.null(gamma_draws)) {
    stopifnot(
      is.numeric(gamma_draws),
      identical(as.integer(dim(gamma_draws)), c(info$n, info$D, n_draws)) )
    storage.mode(gamma_draws) <- 'double'
  }

  res <- .Call(
    C_dirmult_ttest, counts, margin, ct$side, pc, n_draws,
    (1 + conf_level) / 2, seed, gamma_draws, details, cpus )

  p_adjusted <- stats::p.adjust(res$p_value, p_adjust)

  lo <- res$ci_lower
  hi <- res$ci_upper
  reject <- p_adjusted <= alpha & ((lo > 0 & hi > 0) | (lo < 0 & hi < 0))
  reject[is.na(reject)] <- FALSE

  D <- info$D
  features <- if (is.null(info$features)) as.character(seq_len(D)) else info$features

  out <- data.frame(
    feature    = features,
    t          = res$t,
    log2_fc    = res$log2_fc,
    ci_lower   = lo,
    ci_upper   = hi,
    p_value    = res$p_value,
    p_adjusted = p_adjusted,
    reject     = reject,
    stringsAsFactors = FALSE )

  attr(out, 'contrast')    <- ct$contrast
  attr(out, 'n')           <- c(treatment = n1, reference = n2)
  attr(out, 'n_draws')     <- n_draws
  attr(out, 'pseudocount') <- pc
  attr(out, 'alpha')       <- alpha
  attr(out, 'conf_level')  <- conf_level
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'seed')        <- seed
  class(out) <- c('ecodive_dirmult_ttest', 'data.frame')

  if (!details) return (out)

  draws <- lapply(res$draws, function (m) { colnames(m) <- features; m })
  list(result = out, draws = draws)
}



# The table's checks: any input type, numeric, finite and non-negative, at
# least two features. Returns list(D, n, features, samples).
dirmult_input <- function (env) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- env$counts
  margin <- env$margin

  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  pkg  <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  if (D < 2L)
    stop('`counts` must have at least 2 features for a Dirichlet-multinomial t-test.', call. = FALSE)

  list(
    D        = D,
    n        = n,
    features = if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]]),
    samples  = if (length(dn[[margin]]))      as.character(dn[[margin]]) )
}



# Each sample's side of the contrast (1 treatment, 2 reference, 0 left
# out), as scikit-bio's _check_trt_ref_groups(), and the contrast's labels.
# `val` is each sample's group (unused factor levels dropped).
dirmult_contrast <- function (val, contrast) {

  lev <- if (is.factor(val)) levels(val) else sort(unique(val), method = 'radix')
  lab <- as.character(lev)
  gi  <- if (is.factor(val)) as.integer(val) else match(val, lev)

  if (length(lev) < 2L)
    stop('`groups` must have at least two different values.', call. = FALSE)

  if (is.null(contrast)) {
    trt <- 1L
    ref <- NA_integer_
  }
  else {
    if (is.factor(contrast)) contrast <- as.character(contrast)
    if (!is.atomic(contrast) || length(contrast) != 2L || is.na(contrast[[1]]))
      stop(
        '`contrast` must be two groups, c(treatment, reference), or ',
        'c(treatment, NA) to compare with all other samples.', call. = FALSE )

    find <- function (x) {
      i <- match(as.character(x), lab)
      if (is.na(i))
        stop(
          '`contrast` group "', x, '" is not one of the groups in `groups` ',
          '(with samples): ', paste0('"', lab, '"', collapse = ', '), '.',
          call. = FALSE )
      i
    }
    trt <- find(contrast[[1]])
    ref <- if (is.na(contrast[[2]])) NA_integer_ else find(contrast[[2]])
    if (identical(trt, ref))
      stop('The two groups in `contrast` must be different.', call. = FALSE)
  }

  side <- integer(length(gi))
  side[gi == trt] <- 1L
  if (is.na(ref)) side[gi != trt] <- 2L
  else            side[gi == ref] <- 2L

  list(
    side     = side,
    contrast = c(treatment = lab[[trt]], reference = if (is.na(ref)) NA_character_ else lab[[ref]]) )
}
