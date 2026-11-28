# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# fastANCOM, written from Zhou et al. (2022)'s abstract, the documentation
# of fastANCOM 0.0.4 (its Rd page, README and manual) and its outputs
# (fastANCOM is GPL-3; its code was not read). How each step was settled
# from outputs is recorded in plans/fastancom.md.



#' ANCOM from Marginal Models (fastANCOM)
#'
#' Tests each feature for differential abundance as ANCOM does, by testing
#' its log ratio with every other feature, but from one linear model per
#' feature instead of one per pair (Zhou et al. 2022). Each feature also
#' gets a log2 fold change against a reference set of features that seem
#' not to change, with a standard error and a p-value. Covariates are
#' allowed.
#'
#' @inheritParams linda
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of [replace_zeros()]
#'        are all accepted: each sample's scale cancels out of the log
#'        ratios.
#'
#' @param formula  The model, a one-sided formula over the columns of
#'        `data`, such as `~ status + age`, or the same as a string. Every
#'        factor (and character or logical column) is coded with treatment
#'        contrasts, its first level (in sorted order, or a factor's own
#'        order) the reference. The model must have an intercept. `NULL`
#'        (the default) is `~ groups`.
#'
#' @param groups  The term of interest, whose effect is tested: with
#'        `formula`, the name of one of its terms (required when `formula`
#'        has more than one); alone, the name of a column of `data` or a
#'        vector with one value per sample, and the model is `~ groups`.
#'        fastANCOM's `x`; the other terms are its covariates `Z`.
#'
#' @param pseudocount  The value added to every abundance before taking
#'        logs. `NULL` (the default) follows the shared rule: nothing is
#'        added when `counts` has no zeros, else half its smallest non-zero
#'        value (0.5 for counts), with a warning. fastANCOM adds 0.5 always.
#'        `0` is allowed for a table without zeros.
#'
#' @param alpha  The significance level of each log ratio's test, for
#'        counting \eqn{W}. fastANCOM's `sig`. Default `0.05`.
#'
#' @param detect_rate  A feature with \eqn{W} at least `detect_rate`
#'        \eqn{\times D} is declared differentially abundant (and more may
#'        be, by ANCOM's own cutoff; see **Details**). Strictly between 0
#'        and 1. fastANCOM's `detect.rate`. Default `0.7`.
#'
#' @param ref_rate  The proportion of the features, those with the smallest
#'        \eqn{W}, that make the reference set for the log2 fold changes.
#'        Above 0 and at most 0.5. fastANCOM's `ref.rate`. Default `0.05`.
#'
#' @param p_adjust  How to adjust the log2 fold changes' p-values for
#'        testing every feature, separately for each comparison: any method
#'        of [stats::p.adjust()] (`'holm'`, the default, as ecodive's other
#'        DA methods; `'BH'`, as fastANCOM's `log2FC.qval`; `'hochberg'`,
#'        `'hommel'`, `'bonferroni'`, `'BY'`, `'fdr'`, `'none'`). `NaN`
#'        p-values are left out of the family. The decisions (`reject`)
#'        don't use them.
#'
#' @return A data frame of class `c('ecodive_fastancom', 'data.frame')` with
#'   a row per feature and comparison, in the order of `counts`:
#'
#'   * `feature`: the feature's name, or its index as text.
#'   * `covariate`: the comparison. For a numeric term or a factor of two
#'     levels, its design column, as [stats::model.matrix()] names it
#'     (`'groupb'`, `'age'`). For a factor of more levels, first the global
#'     test of all its levels (the term's name, `'status'`), then each level
#'     against the reference level (`'statusb'`, `'statusc'`), then each
#'     other pair of levels (`'statusc - statusb'`).
#'   * `log2_fc`: the log2 fold change against the reference set (`NaN`
#'     for a global test).
#'   * `se`: its standard error (`NaN` for a global test).
#'   * `statistic`: `log2_fc / se`, a t statistic; for a global test, the F
#'     statistic of the feature's coefficients against the reference set's.
#'   * `w`: the number of other features whose log ratio with this one
#'     differs significantly (integer, `0` to \eqn{D - 1}).
#'   * `p_value`: the t-test of `log2_fc`, or the global F-test, on \eqn{n
#'     - p} degrees of freedom.
#'   * `p_adjusted`: adjusted with `p_adjust`, per comparison.
#'   * `reject`: `TRUE` for the features declared differentially abundant,
#'     from \eqn{W} (fastANCOM's `REJECT`; see **Details**), not from
#'     `p_adjusted`. For decisions from the p-values instead, use
#'     `p_adjusted <= alpha`.
#'
#'   Attributes: `method` (`'fastancom'`), `formula`, `groups` (the term of
#'   interest), `design` (the model matrix), `alpha`, `detect_rate`,
#'   `ref_rate`, `p_adjust`, `pseudocount` (the value added; `0` for none),
#'   `log_base` (`2`), `dof` (\eqn{n - p}), `cutoff` (per comparison, the
#'   smallest \eqn{W} that is rejected), `reference` (per comparison, the
#'   features of the reference set), `coef` (features \eqn{\times} design
#'   columns: each feature's own log2 coefficients, fastANCOM's marginal
#'   `beta`) and `mean` (each feature's mean log2 abundance). [plot()] draws
#'   the \eqn{W} of the features, a volcano plot, an MA plot, a forest plot
#'   or a rank plot ([plot.ecodive_fastancom()]).
#'
#' @details
#' **Log ratios from marginal models.** With \eqn{L_k = \log_2(x_k + c)}
#' the logs of feature \eqn{k}'s abundances (\eqn{c} the pseudocount) and
#' \eqn{X} the design (\eqn{n} samples, \eqn{p} columns), each feature gets
#' its own least squares fit, \eqn{L_k = X b_k + e_k}. Because every
#' feature has the same design, the fit of the log ratio \eqn{L_k - L_j} is
#' the difference of the two fits: coefficients \eqn{b_k - b_j}, residuals
#' \eqn{e_k - e_j}. So the test of the term of interest on any log ratio,
#' which ANCOM would fit separately for each of the \eqn{D(D - 1)/2} pairs,
#' follows from the \eqn{D} fits: a t-test on \eqn{n - p} degrees of
#' freedom, with variance \eqn{c_{xx} \lVert e_k - e_j \rVert^2 / (n - p)}
#' (\eqn{c_{xx}} the term's entry of \eqn{(X^\top X)^{-1}}), or for a
#' factor of more than two levels the F-test of all its coefficients.
#' \eqn{W_k} is the number of features \eqn{j} with \eqn{p_{kj} <}
#' `alpha`, with no adjustment. These are exactly ANCOM's tests: without
#' covariates and for a categorical term, \eqn{W} is [ancom()]'s with
#' `sig_test = 'anova'` and `p_adjust = 'none'`.
#'
#' **Decisions.** A feature is rejected when \eqn{W_k \ge} `detect_rate`
#' \eqn{\times D}, or when ANCOM's own empirical cutoff (Mandal et al.
#' 2015) detects it and \eqn{W_k \ge D/2}. That cutoff takes \eqn{c_0 =
#' \max(W)/D} and the five cutoffs \eqn{c_0 - 0.05, \ldots, c_0 - 0.25};
#' with \eqn{\pi_i} the proportion of features with \eqn{W \ge D c_i} and
#' \eqn{\delta_i = |\pi_i - \pi_{i+1}|}, it chooses \eqn{c_1} if
#' \eqn{\delta_1, \delta_2, \delta_3 < 0.02}, else \eqn{c_2} if only
#' \eqn{\delta_1 \ge 0.02} of them, else \eqn{c_3} if \eqn{\delta_2 \ge
#' 0.02} and \eqn{\delta_3, \delta_4 < 0.02}, else \eqn{c_4}; and detects
#' the features with \eqn{W} at least the smallest \eqn{W \ge D \nu}. It
#' detects nothing for fewer than 10 features or when \eqn{c_0 < 0.26}. The
#' result's `cutoff` attribute gives the resulting smallest rejected
#' \eqn{W}.
#'
#' **Log fold changes.** The reference set \eqn{R} is the `ref_rate`
#' \eqn{\times D} features (rounded down, at least one) with the smallest
#' \eqn{W}, ties in the order of `counts`: those whose abundance most likely
#' didn't change. Each feature's log2 fold change is its
#' coefficient less the reference set's mean, \eqn{b_k - \bar b_R}: the
#' coefficient of the log ratio of the feature to the geometric mean of the
#' reference set. Its standard error comes from the residuals \eqn{e_k -
#' \bar e_R} in the same way, and its p-value from the t distribution on
#' \eqn{n - p} degrees of freedom. A feature alone in the reference set has
#' a log2 fold change of 0 and no p-value (`NaN`).
#'
#' **Differences from fastANCOM.** The results are those of fastANCOM
#' 0.0.4's `fastANCOM()` (to about \eqn{10^{-12}}; \eqn{W} and `reject`
#' identical) with `zero_cut = 1` and `pseudo = 0.5` given, and `p_adjust =
#' 'BH'`, for a numeric term or a factor of two levels, except that:
#'
#' * fastANCOM's standard errors take the reference set's mean residuals
#'   as their sum divided by `ref_rate` \eqn{\times D}, not by the size of
#'   the set, while its log2 fold changes divide by the size; they agree
#'   when `ref_rate` \eqn{\times D} is a whole number (or below 1).
#'   `fastancom()` divides by the size for both, so that `se` is the
#'   standard error of `log2_fc`.
#' * With covariates, fastANCOM tests on \eqn{n - 2} degrees of freedom
#'   (as if there were none) while estimating the variances on \eqn{n - p};
#'   `fastancom()` uses \eqn{n - p} for both, as a linear model does. The
#'   log fold changes and standard errors are the same; the p-values,
#'   \eqn{W} and decisions can differ.
#' * fastANCOM's filters (`zero_cut`, `lib_cut`), structural zeros
#'   (`struc_zero`), outlier removal (`outlier`) and random effects
#'   (`rand`) are not carried over: filter `counts` first.
#' * `pseudocount = NULL` adds nothing to a table without zeros (fastANCOM
#'   adds `pseudo = 0.5` always).
#' * Covariates may be any columns of `data` (fastANCOM: a numeric matrix
#'   `Z`), coded as in [linda()].
#' * For a factor of more than two levels, fastANCOM 0.0.4's `global` and
#'   pairwise comparisons fail with an error on every input tried, so these
#'   are ecodive's own: the F-test of all the term's coefficients for the
#'   global test (\eqn{W} as [ancom()]'s with `sig_test = 'anova'` without
#'   covariates), and each pair of levels as a contrast of the same fits
#'   (the variance from all the samples, as in a linear model).
#' * `p_adjust` defaults to `'holm'` (fastANCOM: `'BH'`).
#'
#' **Computation.** The \eqn{D} fits use one QR decomposition of the
#' design; the \eqn{D(D - 1)/2} tests then take \eqn{O(n)} each, in
#' parallel C, without forming any log ratio or the \eqn{D \times D}
#' covariance matrix. Results are identical for any number of `cpus`.
#'
#' In benchmarks on real data, fastANCOM was the only microbiome-specific
#' method that controlled false discoveries at a relatively high
#' sensitivity (Wirbel et al. 2024), though with a lower sensitivity than
#' others in Pelto et al. (2025).
#'
#' @section Same names in other packages:
#'   fastANCOM (an R package on GitHub) has `fastANCOM()`, the original.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Zhou, C., Wang, H., Zhao, H., & Wang, T. (2022). fastANCOM: a fast
#' method for analysis of compositions of microbiomes. *Bioinformatics*,
#' 38(7), 2039-2041. \doi{10.1093/bioinformatics/btac060}
#'
#' Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., &
#' Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel
#' method for studying microbial composition. *Microbial Ecology in Health
#' and Disease*, 26, 27663.
#'
#' Wirbel, J., Essex, M., Forslund, S. K., & Zeller, G. (2024). A realistic
#' benchmark for differential abundance testing and confounder adjustment in
#' human microbiome studies. *Genome Biology*, 25, 247.
#' \doi{10.1186/s13059-024-03390-9}
#'
#' Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary
#' methods provide more replicable results in microbial differential
#' abundance analysis. *Briefings in Bioinformatics*, 26, bbaf130.
#' \doi{10.1093/bib/bbaf130}
#'
#' @seealso [da_test()], [ancom()], [ancombc()], [linda()],
#'   [replace_zeros()].
#'
#' @export
#' @examples
#'     # fastANCOM's example: 60 samples, 100 features; the first 10
#'     # are more abundant in group 1.
#'     set.seed(123)
#'     x <- matrix(rpois(100 * 60, 10), 60)
#'     b <- rnorm(10, 4)
#'     x[31:60, 1:10] <- t(t(x[31:60, 1:10]) * 2^b)
#'     colnames(x) <- paste0('f', 1:100)
#'     grp <- rep(c(0, 1), each = 30)
#'
#'     res <- fastancom(x, groups = grp, pseudocount = 0.5)
#'     res
#'     plot(res)
#'     plot(res, type = 'volcano')
#'
#'     # A covariate, and three groups.
#'     meta <- data.frame(
#'       grp = rep(c('a', 'b', 'c'), each = 20),
#'       age = round(runif(60, 20, 60)) )
#'     fastancom(x, ~ grp + age, meta, groups = 'grp', pseudocount = 0.5)
#'
fastancom <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    pseudocount = NULL, alpha = 0.05, detect_rate = 0.7, ref_rate = 0.05,
    p_adjust = 'holm', cpus = n_cpus() ) {

  fastancom_fit(environment())
}



fastancom_fit <- function (env) {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  pc <- env$pseudocount
  if (!is.null(pc) && (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc < 0))
    stop('`pseudocount` must be NULL or a single non-negative number.', call. = FALSE)

  info <- log_ratio_input('fastANCOM', env = env, zero_ok = TRUE, dense_out = FALSE)
  counts <- env$counts
  pc     <- as.double(env$pseudocount)
  if (pc == 0 && fastancom_has_zeros(counts))
    stop('`counts` has zeros, and `pseudocount = 0` leaves their logs undefined.', call. = FALSE)

  for (arg in c('alpha', 'detect_rate')) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !(v > 0 && v < 1))
      stop('`', arg, '` must be a single number strictly between 0 and 1.', call. = FALSE)
  }
  alpha       <- as.double(env$alpha)
  detect_rate <- as.double(env$detect_rate)

  ref_rate <- env$ref_rate
  if (!is.numeric(ref_rate) || length(ref_rate) != 1L || is.na(ref_rate) ||
      !(ref_rate > 0 && ref_rate <= 0.5))
    stop('`ref_rate` must be a single number above 0 and at most 0.5.', call. = FALSE)
  ref_rate <- as.double(ref_rate)

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)


  #--------------------------------------------------------
  # The design and the term of interest.
  #--------------------------------------------------------
  rank_def <- FALSE
  des <- withCallingHandlers(
    da_design(env$formula, env$data, env$groups, info$samples, info$n, 'fastANCOM'),
    warning = function (w) {
      if (grepl('rank deficient', conditionMessage(w), fixed = TRUE)) {
        rank_def <<- TRUE
        invokeRestart('muffleWarning')
      }
    })
  if (rank_def)
    stop(
      'The design is rank deficient: some of its columns are combinations ',
      'of others. Remove the redundant terms from `formula`.', call. = FALSE )

  X <- des$X
  n <- info$n
  D <- info$D
  P <- ncol(X)

  if (des$intercept < 0L)
    stop(
      'fastANCOM needs a model with an intercept: `formula` must not ',
      'remove it (`- 1` or `+ 0`).', call. = FALSE )

  labels <- attr(X, 'term.labels')
  term   <- des$groups
  if (is.null(term)) {
    if (length(labels) != 1L)
      stop(
        '`groups` is required when `formula` has more than one term: give ',
        'the term of interest, one of ', name_list(labels), '.', call. = FALSE )
    term <- labels
  }
  tcols <- which(attr(X, 'assign') == match(term, labels))

  cmp <- fastancom_comparisons(X, term, tcols)


  #--------------------------------------------------------
  # Each comparison's inverse covariance factor and F
  # critical value.
  #--------------------------------------------------------
  dof  <- n - P
  dX   <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'
  qx   <- qr(dX)
  xtxi <- matrix(0, P, P)
  xtxi[qx$pivot, qx$pivot] <- chol2inv(qx$qr[seq_len(P), seq_len(P), drop = FALSE])

  cmp_c <- lapply(cmp, function (cm) {
    L <- cm$L
    V <- crossprod(L, xtxi %*% L)
    list(L = L, A = solve(V), crit = stats::qf(1 - alpha, ncol(L), dof))
  })

  # As fastANCOM: the count of 1:(ref_rate * D), whose `:` adds a fuzz.
  nref <- length(1:(ref_rate * D))

  res <- .Call(
    C_fastancom, counts, env$margin, pc, dX, des$intercept, cmp_c, nref,
    env$cpus )


  #--------------------------------------------------------
  # p-values and decisions.
  #--------------------------------------------------------
  features <- info$features
  if (is.null(features)) features <- as.character(seq_len(D))
  M     <- length(cmp)
  nms   <- vapply(cmp, `[[`, '', 'name')
  q     <- vapply(cmp, function (cm) ncol(cm$L), integer(1))

  pval <- padj <- res$stat
  rej  <- matrix(FALSE, D, M)
  cut  <- stats::setNames(numeric(M), nms)
  ref  <- stats::setNames(vector('list', M), nms)

  for (c in seq_len(M)) {
    st <- res$stat[, c]
    pval[, c] <- if (q[[c]] == 1L) 2 * stats::pt(-abs(st), dof)
                 else stats::pf(st, q[[c]], dof, lower.tail = FALSE)
    padj[, c] <- da_adjust(pval[, c], p_adjust, alpha)$p_adjusted
    dec       <- fastancom_reject(res$W[, c], D, detect_rate)
    rej[, c]  <- dec$reject
    cut[[c]]  <- dec$cutoff
    ref[[c]]  <- features[res$ref[[c]]]
  }

  byrow <- function (m) { v <- as.vector(t(m)); v[is.na(v)] <- NaN; v }

  out <- data.frame(
    feature    = rep(features, each = M),
    covariate  = rep(nms, times = D),
    log2_fc    = byrow(res$est),
    se         = byrow(res$se),
    statistic  = byrow(res$stat),
    w          = as.vector(t(res$W)),
    p_value    = byrow(pval),
    p_adjusted = byrow(padj),
    reject     = as.vector(t(rej)),
    stringsAsFactors = FALSE )

  coef <- res$coef
  dimnames(coef) <- list(features, colnames(X))

  attr(out, 'method')      <- 'fastancom'
  attr(out, 'formula')     <- des$formula
  attr(out, 'groups')      <- term
  attr(out, 'design')      <- X
  attr(out, 'alpha')       <- alpha
  attr(out, 'detect_rate') <- detect_rate
  attr(out, 'ref_rate')    <- ref_rate
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'pseudocount') <- pc
  attr(out, 'log_base')    <- 2
  attr(out, 'dof')         <- dof
  attr(out, 'cutoff')      <- cut
  attr(out, 'reference')   <- ref
  attr(out, 'coef')        <- coef
  attr(out, 'mean')        <- stats::setNames(res$mean, features)
  class(out) <- c('ecodive_fastancom', 'data.frame')

  out
}



# Whether a table has zeros (stored or, for a sparse table, unstored).
fastancom_has_zeros <- function (counts) {
  vals <- stored_values(counts)
  if (any(vals == 0)) return (TRUE)
  switch(
    get_matrix_package(counts),
    'base'   = FALSE,
    'slam'   = length(vals) < as.double(counts$nrow) * counts$ncol,
    'Matrix' = length(vals) < prod(as.double(dim(counts))) )
}



# The comparisons for the term of interest (columns `tcols` of the design
# X): a list of list(name, L), L the P x q contrasts. One column: itself. A
# factor of L > 2 levels (treatment contrasts): the global test of its
# columns, each column (a level against the reference), then each other
# pair of levels. Other terms of several columns: the global test only.
fastancom_comparisons <- function (X, term, tcols) {

  P   <- ncol(X)
  cn  <- colnames(X)
  e   <- function (j) { v <- numeric(P); v[j] <- 1; matrix(v, P, 1L) }

  if (length(tcols) == 1L)
    return (list(list(name = cn[[tcols]], L = e(tcols))))

  glob <- matrix(0, P, length(tcols))
  glob[cbind(tcols, seq_along(tcols))] <- 1
  out  <- list(list(name = term, L = glob))

  lev <- attr(X, 'xlevels')[[term]]
  if (is.null(lev) || length(lev) != length(tcols) + 1L) return (out)

  for (j in tcols) out[[length(out) + 1L]] <- list(name = cn[[j]], L = e(j))
  for (a in seq_along(tcols)) for (b in seq_along(tcols)) {
    if (b <= a) next
    out[[length(out) + 1L]] <- list(
      name = paste(cn[[tcols[[b]]]], '-', cn[[tcols[[a]]]]),
      L    = e(tcols[[b]]) - e(tcols[[a]]) )
  }
  out
}



# The original ANCOM's empirical cutoff (Mandal et al. 2015), as
# fastANCOM applies it: the smallest W detected, or NA when nothing is
# (fewer than 10 features, or max(W) / K below 0.26).
fastancom_detect <- function (W, K) {

  if (K < 10L) return (NA_real_)
  c0 <- max(W) / K
  if (!(c0 >= 0.26)) return (NA_real_)

  cut <- c0 - c(0.05, 0.10, 0.15, 0.20, 0.25)
  pc  <- vapply(cut, function (cc) mean(W >= K * cc), numeric(1))
  del <- abs(diff(pc))
  tau <- 0.02

  nu <- if      (del[1] <  tau && del[2] < tau && del[3] < tau) cut[1]
        else if (del[1] >= tau && del[2] < tau && del[3] < tau) cut[2]
        else if (del[2] >= tau && del[3] < tau && del[4] < tau) cut[3]
        else                                                   cut[4]

  as.double(min(W[W >= nu * K]))
}



# fastANCOM's REJECT from the W of K features: W >= detect_rate * K, or
# detected by ANCOM's cutoff with W >= K / 2. Returns list(reject, cutoff),
# the cutoff being the smallest W that is rejected.
fastancom_reject <- function (W, K, detect_rate) {
  up  <- fastancom_detect(W, K)
  thr <- detect_rate * K
  if (is.na(up)) return (list(reject = W >= thr, cutoff = thr))
  list(
    reject = W >= thr | (W >= up & W >= 0.5 * K),
    cutoff = min(thr, max(up, 0.5 * K)) )
}



#' Plot a fastANCOM Result
#'
#' Views of a [fastancom()] result for one comparison, with the features
#' declared differentially abundant (`reject = TRUE`) highlighted:
#'
#' * `"w"` (the default): a dot chart of the features' \eqn{W}, those with
#'   the largest at the top, with a dashed line at the smallest \eqn{W}
#'   that is rejected (the `cutoff` attribute).
#' * `"volcano"`: the log2 fold change against the adjusted p-value.
#' * `"ma"`: the log2 fold change against the mean log2 abundance.
#' * `"forest"`: the log2 fold changes and their confidence intervals of
#'   the `n_features` most significant features.
#' * `"rank"`: every feature's log2 fold change, sorted.
#'
#' @param x  The result of [fastancom()].
#'
#' @param covariate  Which comparison to plot, as named in the result's
#'        `covariate` column. Defaults to the global test of the term of
#'        interest when there is one, else to its first column.
#'
#' @param type  The plot: `"w"` (the default), `"volcano"`, `"ma"`,
#'        `"forest"`, or `"rank"`. A global test has \eqn{W} but no log2
#'        fold change, so only `"w"`.
#'
#' @param n_features  For `"w"`: how many features to show, those with the
#'        largest \eqn{W}; for `"forest"`, those with the smallest adjusted
#'        p-values. Default `30`; `Inf` for all.
#'
#' @inheritParams plot.ecodive_linda
#'
#' @param ...  Passed on to [graphics::dotchart()] (`"w"`) or
#'        [graphics::plot.default()], such as `main`, `sub`, `xlab` or
#'        `cex`.
#'
#' @inheritSection plot.ecodive_ancombc Rank plots
#'
#' @details
#' Log fold changes are in log2, against the reference set (see
#' [fastancom()]). In the volcano plot, the y axis is \eqn{-\log_{10}} of
#' `p_adjusted`, with a dashed line at `alpha`; fastANCOM's decisions come
#' from \eqn{W}, so the highlighted features needn't be those above the
#' line. The forest plot's bars are \eqn{\mathrm{log2\_fc} \pm
#' q\,\mathrm{se}}, with \eqn{q} the \eqn{(1 + }`conf_level`\eqn{)/2}
#' quantile of the t distribution with the result's `dof`, not adjusted for
#' testing many features.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     x <- matrix(rpois(20 * 40, 40), 20, 40)
#'     x[11:20, 1:4] <- x[11:20, 1:4] * 5
#'     res <- fastancom(x, groups = rep(c('a', 'b'), each = 10))
#'     plot(res)
#'     plot(res, type = 'volcano')
#'     plot(res, type = 'forest', n_features = 10)
plot.ecodive_fastancom <- function (
    x, covariate = NULL, type = c('w', 'volcano', 'ma', 'forest', 'rank'),
    labels = 10, n_features = 30, conf_level = 0.95, numerator = NULL,
    denominator = NULL, ... ) {

  type <- test_plot_type(type, c('w', 'volcano', 'ma', 'forest', 'rank'))

  if (is.null(covariate)) {
    g <- attr(x, 'groups', exact = TRUE)
    covariate <- if (is.character(g) && length(g) == 1L && g %in% x$covariate) g
                 else da_rank_covariate(x, NULL)
  }
  covariate <- da_rank_covariate(x, covariate)

  r <- x[x$covariate == covariate, , drop = FALSE]
  D <- nrow(r)

  if (type == 'w') {
    rank_only(type, numerator, denominator)
    cut <- attr(x, 'cutoff', exact = TRUE)
    cut <- if (covariate %in% names(cut)) cut[[covariate]]
    plot_ranked(
      r$w, r$feature, r$reject, n_features, list(...),
      main = paste('W of', da_method_name(attr(x, 'method', exact = TRUE))),
      xlab = sprintf('W (of %d log ratios): %s', D - 1, covariate),
      xlim = c(0, max(D - 1, 1)),
      line = cut )
    return (invisible(x))
  }

  if (!any(is.finite(r$log2_fc)))
    stop(
      '`covariate = "', covariate, '"` is a global test, with no log2 fold ',
      'changes: use `type = "w"`, or a pairwise comparison.', call. = FALSE )

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)

  q <- stats::qt((1 + conf_level) / 2, attr(x, 'dof', exact = TRUE))

  da_plot(
    x, type, labels, n_features, list(...),
    est = r$log2_fc, p = r$p_adjusted, reject = r$reject, feature = r$feature,
    mean = attr(x, 'mean', exact = TRUE),
    lo = r$log2_fc - q * r$se, hi = r$log2_fc + q * r$se,
    what     = covariate,
    est_lab  = lfc_label(x),
    mean_lab = 'Mean abundance (log2)',
    interval = sprintf('unadjusted %s confidence intervals', plot_pct(conf_level)),
    num = numerator, den = denominator )

  invisible(x)
}
