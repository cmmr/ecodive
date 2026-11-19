# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# ancombc() and ancombc2() are an R/C port of ancombc(), ancombc2() and
# _ancombc_core() from scikit-bio 0.7.4 (skbio/stats/composition/
# _ancombc.py, with _check_metadata(), _build_dmatrix() and
# _adjust_pvalues() from _utils.py), Copyright (c) 2013--, scikit-bio
# development team, under the BSD-3 license given in src/ancombc.c and
# inst/COPYRIGHTS. scikit-bio's code is itself derived from the R package
# ANCOMBC, Copyright (c) Huang Lin, under the Artistic License 2.0
# (full text in inst/COPYRIGHTS). This file is a modified version of that
# code (Artistic License 2.0, section 4): it builds the design with R's
# model.matrix() rather than patsy, computes the per-covariate statistics
# here and the fits in src/ancombc.c, and returns ecodive's per-feature
# table; the user-visible differences are listed under "Differences from
# scikit-bio and ANCOMBC" in the documentation below, and the rest in
# plans/ancombc.md.



#' Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC)
#'
#' Estimates each feature's log fold change for every covariate of a linear
#' model, corrected for the bias that differing sampling fractions cause,
#' with standard errors and p-values (Lin & Peddada 2020, 2024).
#' `ancombc()` is ANCOM-BC; `ancombc2()` is ANCOM-BC2, which also corrects
#' each sample's logs for its estimated sampling fraction before refitting,
#' adds a variance floor for rare features, uses t tests, and by default
#' treats zeros as missing values instead of adding a pseudocount. Both
#' follow scikit-bio's `ancombc()` and `ancombc2()`.
#'
#' @inheritParams ancom
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of [replace_zeros()]
#'        are all accepted: each sample's scale is absorbed by its sampling
#'        fraction.
#'
#' @param formula  The model, a one-sided formula over the columns of `data`,
#'        such as `~ status + age`, or the same as a string (`'status +
#'        age'`, with or without the `~`). Anything [stats::model.matrix()]
#'        accepts may be used, including interactions (`status * age`) and
#'        a model without an intercept (`~ 0 + status`). Every factor is
#'        coded with treatment contrasts (see **Covariates**).
#'
#' @param data  A data frame of sample metadata, with the variables of
#'        `formula`. Its rows are matched to the samples of `counts` by
#'        name when both have names (rows for other samples are left out),
#'        else by position. A single vector or factor is taken as one
#'        variable named `group`.
#'
#' @param groups  `NULL` (the default), or the name of a categorical
#'        variable of `data` that is a term of `formula`, with at least
#'        three levels: the grouping the post-hoc tests (planned:
#'        `ancombc_global()`, `ancombc_pairwise()`, `ancombc_dunnett()`,
#'        `ancombc_trend()`) compare. It doesn't change the results
#'        (scikit-bio's `grouping`).
#'
#' @param pseudocount  The value added to every abundance before taking
#'        logs. For `ancombc()`, `NULL` (the default) follows the shared
#'        rule: nothing is added when `counts` has no zeros, else half the
#'        smallest non-zero value, with a warning. For `ancombc2()`, the
#'        default `0` adds nothing and treats zeros as missing values, as
#'        scikit-bio and ANCOMBC do; a positive value is added to every
#'        value instead.
#'
#' @param var_quantile  `ancombc2()` only: the quantile of each covariate's
#'        coefficient variances (over the features) added to every variance
#'        of that covariate, which keeps rare features with tiny standard
#'        errors from looking significant (ANCOMBC's `s0_perc`). Between 0
#'        and 1; `0` turns it off. Default: `0.05`.
#'
#' @param max_iters  The most iterations of the bias-estimating EM
#'        algorithm, per covariate. Default: `100`.
#'
#' @param tol  The EM's convergence tolerance, on the Euclidean norm of the
#'        change in its eight parameters. Default: `1e-5`.
#'
#' @param alpha  The significance level of the adjusted p-values: `reject`
#'        is `p_adjusted <= alpha`. Strictly between 0 and 1. Default:
#'        `0.05`.
#'
#' @param p_adjust  How to adjust the p-values for multiple testing, over
#'        the features, separately for each covariate: any method of
#'        [stats::p.adjust()] (`'holm'`, the default, `'hochberg'`,
#'        `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`, `'fdr'`, `'none'`).
#'
#' @return A data frame of class `c('ecodive_ancombc', 'data.frame')`, one
#'   row per feature and covariate (the covariates of the first feature,
#'   then those of the second, and so on), with columns:
#'
#'   * `feature`: the feature's name (its index, as text, when `counts` has
#'     no feature names).
#'   * `covariate`: the column of the design, as [stats::model.matrix()]
#'     names it: `(Intercept)`, `statusmoderate`, `age`, `statusmoderate:age`,
#'     ... (scikit-bio's `Intercept`, `status[T.moderate]`, ...).
#'   * `log_fc`: the bias-corrected coefficient, a natural-log fold change
#'     (scikit-bio's `Log(FC)`).
#'   * `se`: its standard error.
#'   * `w`: the test statistic, `log_fc / se`.
#'   * `p_value`: two-sided, from the standard normal (`ancombc()`) or a t
#'     distribution (`ancombc2()`).
#'   * `p_adjusted`: `p_value` adjusted with `p_adjust`, within each
#'     covariate.
#'   * `reject`: `TRUE` when `p_adjusted <= alpha`.
#'
#'   A coefficient that can't be estimated (with zeros as missing values, a
#'   feature whose observed samples leave a factor with fewer than two
#'   levels, or no samples at all) has `NaN` for `log_fc`, `se` and `w`, a
#'   `p_value` of 1 (which counts in the adjustment) and `reject = FALSE`,
#'   as in scikit-bio and ANCOMBC.
#'
#'   Attributes:
#'
#'   * `method`: `'ancombc'` or `'ancombc2'`.
#'   * `formula`, `groups`, `alpha`, `p_adjust`, `max_iters`, `tol`,
#'     `pseudocount` (the value used), and for `ancombc2()` `var_quantile`.
#'   * The fit, for the post-hoc tests: `coef` (features \eqn{\times}
#'     covariates: the bias-corrected coefficients, including those that
#'     can't be estimated), `var` (their variances; for `ancombc2()`
#'     adjusted for the bias's uncertainty and `var_quantile`), `vcov`
#'     (features \eqn{\times} covariates \eqn{\times} covariates:
#'     each feature's covariance matrix, with `var` on its diagonal), `dof`
#'     (`NULL` for `ancombc()`; each feature's residual degrees of freedom
#'     for `ancombc2()`), `estimable` (`NULL` when every coefficient could
#'     be estimated, else a logical features \eqn{\times} covariates
#'     matrix), and `design` (the samples \eqn{\times} covariates model
#'     matrix, with [stats::model.matrix()]'s `assign` and `contrasts`
#'     attributes).
#'   * `bias`: covariates \eqn{\times} 3, each covariate's estimated bias
#'     (`delta_em`, subtracted from the coefficients), its weighted least
#'     squares estimate (`delta_wls`) and that estimate's variance
#'     (`var_delta`); `n_iters`: the EM iterations run per covariate.
#'   * `sample_fractions` (`ancombc2()`): each sample's estimated log
#'     sampling fraction.
#'
#' @details
#' **The sampling fraction.** A sample's sequencing reads are an unknown
#' fraction of its community's absolute abundance, and that fraction differs
#' from sample to sample. In the log scale it shifts every feature of a
#' sample equally, so a naive log fold change between groups is off by the
#' groups' difference in mean log sampling fraction, the same for every
#' feature: ANCOM-BC estimates that shift (the bias) from the distribution
#' of all features' coefficients and removes it.
#'
#' With \eqn{L} the \eqn{n \times D} table of logs (natural logs of the
#' abundance plus `pseudocount`) and \eqn{X} the design:
#'
#' 1. Each feature's coefficients are fitted by least squares,
#'    \eqn{\hat\beta_f = X^+ L_f}, and each sample's mean residual over the
#'    features estimates its log sampling fraction. Each coefficient's
#'    variance is a sandwich estimate from the residuals around those.
#' 2. For each covariate, an EM algorithm fits a three-part normal mixture to
#'    the features' coefficients: null features, centred on the bias
#'    \eqn{\delta}, and features that decrease or increase. The bias is
#'    subtracted from that covariate's coefficients.
#' 3. `ancombc()` tests \eqn{W = \hat\beta / \mathrm{se}} against the
#'    standard normal.
#'
#' `ancombc2()` centres each feature's logs first; after step 2 it estimates
#' every sample's sampling fraction, subtracts it from the sample's logs,
#' and fits again (now without sample effects); it then widens each variance
#' to \eqn{(\mathrm{se} + \sqrt{\mathrm{var}_\delta})^2}, for the bias's
#' uncertainty, adds the `var_quantile` quantile of the covariate's
#' variances, and tests against a t distribution with each feature's
#' residual degrees of freedom (observed samples minus the rank of its
#' design). With zeros as missing values, each feature is fitted on its
#' observed samples, the sample effects and coefficients being refined in
#' turn (up to 20 rounds, to a tolerance of 0.01), as in ANCOMBC. Because
#' that fit stops early, rescaling some samples' abundances moves the results
#' a little (by up to about 0.1 in `log_fc` in a test with 30 samples), as it
#' does in scikit-bio and ANCOMBC; without missing values, any per-sample
#' scaling is absorbed exactly.
#'
#' The EM stops after `max_iters` iterations even if it hasn't converged
#' (`attr(, 'n_iters')` shows how many ran). Its estimate then depends on the
#' exact path of its minimisations, so it can differ from ANCOMBC's by much
#' more than when it converges (by 0.02 in one test).
#'
#' @section Covariates:
#'   Character and logical variables become factors; every factor is coded
#'   with treatment contrasts, even an ordered one, and its unused levels are
#'   dropped. The **reference level** is a factor's first level, or for
#'   character variables the first value in sorted order (`sort(method =
#'   'radix')`, C-locale order, as scikit-bio's), so that the results match
#'   scikit-bio's. To choose the reference, make the variable a factor with
#'   that level first, e.g. with [stats::relevel()]. A numeric variable is
#'   one covariate, its slope.
#'
#'   With zeros as missing values (`ancombc2()`), a feature that is missing
#'   from every sample of a factor's reference level is fitted against its
#'   first observed level instead, and the coefficients re-expressed, as
#'   ANCOMBC does (not for factors in an interaction). A feature whose
#'   observed samples leave a factor with fewer than two levels can't be
#'   estimated (see **Value**).
#'
#' @section Differences from scikit-bio and ANCOMBC:
#'   The results are those of scikit-bio's `ancombc()` and `ancombc2()`
#'   (0.7.4), whose algorithm is followed step by step, including its
#'   reproduction of ANCOMBC's handling of zeros, failed fits and missing
#'   reference levels. The EM's one-dimensional minimisations use SciPy's
#'   Nelder-Mead rules, as scikit-bio does; ANCOMBC uses `nloptr`'s, so its
#'   estimates differ slightly, by about `1e-5`, when the EM converges.
#'   Also, unlike ANCOMBC:
#'
#'   * `var_quantile = 0` adds nothing, as in scikit-bio; ANCOMBC's
#'     `s0_perc = 0` adds the smallest variance.
#'   * The pseudoinverse is of the design itself, with numpy's cutoff;
#'     ANCOMBC inverts \eqn{X'X} with `MASS::ginv()`, whose tolerance drops
#'     columns of a badly scaled design (e.g. a covariate in the tens of
#'     thousands), and then gives wrong standard errors.
#'   * With zeros as missing values, a factor in an interaction is not
#'     re-expressed against a feature-local reference level (scikit-bio's
#'     choice), so such features' coefficients differ from ANCOMBC's.
#'
#'   And unlike both:
#'
#'   * The covariates are named as [stats::model.matrix()] names them, and
#'     ordered as in R (patsy puts categorical terms before numeric ones).
#'   * Zeros in `ancombc()`: scikit-bio refuses them; ecodive applies its
#'     shared pseudocount rule (see `pseudocount`).
#'   * Not carried over: random effects (ANCOMBC's `rand_formula`),
#'     scikit-bio's `aggregator`, ANCOMBC's data preprocessing (`prv_cut`,
#'     `lib_cut`, `tax_level`; see the examples), its structural-zero
#'     detection and its pseudocount sensitivity analysis (repeat the call
#'     with other `pseudocount` values).
#'   * A design with as many columns as samples, or more, is an error (the
#'     residuals would all be zero); a rank-deficient design gives a
#'     warning.
#'   * The p-values of each covariate are adjusted with [stats::p.adjust()]
#'     over the features whose p-value isn't `NaN`, as scikit-bio does.
#'
#' **Computation.** The per-feature fits and variances run in parallel C,
#' over blocks of features, and the EM in parallel over covariates. Results
#' are identical for any number of `cpus`, and for dense and sparse
#' `counts`. Memory: the \eqn{n \times D} table of logs, plus \eqn{P^2}
#' doubles per feature for the covariances (and, with zeros as missing,
#' each feature's masked-design decomposition).
#'
#' @section Same names in other packages:
#'   ANCOMBC (Bioconductor) has `ancombc()` and `ancombc2()`, the original
#'   implementations, which take a phyloseq or TreeSummarizedExperiment
#'   object and preprocess it; scikit-bio (Python) has the functions this
#'   follows.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Lin, H., & Peddada, S. D. (2020). Analysis of compositions of
#' microbiomes with bias correction. *Nature Communications*, 11, 3514.
#'
#' Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
#' microbiomes with covariate adjustments and repeated measures. *Nature
#' Methods*, 21, 83–91.
#'
#' @seealso [ancom()], [dirmult_ttest()], [replace_zeros()].
#'
#' @export
#' @examples
#'     # scikit-bio's example: three disease states and age.
#'     x <- rbind(
#'       c(2, 0, 7, 0, 0, 2, 3, 2), c(0, 2, 0, 1, 1, 2, 0, 0),
#'       c(3, 1, 0, 9, 0, 1, 1, 0), c(2, 0, 1, 8, 1, 0, 11, 46),
#'       c(1, 0, 1, 1, 0, 0, 2, 2), c(0, 3, 22, 1, 1, 1, 3, 0),
#'       c(1, 7, 16, 1, 0, 0, 2, 2), c(0, 5, 6, 1, 2, 1, 0, 1),
#'       c(1, 7, 0, 2, 1, 0, 3, 2), c(0, 6, 4, 2, 0, 2, 2, 1),
#'       c(3, 13, 7, 0, 0, 0, 3, 9), c(1, 8, 5, 1, 0, 0, 0, 0),
#'       c(0, 5, 14, 1, 0, 1, 0, 1), c(5, 26, 3, 2, 0, 3, 1, 3),
#'       c(0, 18, 7, 0, 0, 3, 1, 0) )
#'     dimnames(x) <- list(paste0('S', 1:15), paste0('F', 1:8))
#'     meta <- data.frame(
#'       status = rep(c('mild', 'moderate', 'severe'), each = 5),
#'       age    = c(39, 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36),
#'       row.names = rownames(x) )
#'
#'     res <- ancombc(x + 1, ~ status + age, meta)
#'     res[res$reject & res$covariate != '(Intercept)', ]
#'
#'     # ANCOM-BC2, with zeros as missing values.
#'     res2 <- ancombc2(x, ~ status + age, meta, groups = 'status')
#'     head(res2)
#'
#'     # ANCOMBC's preprocessing, done first: keep features present in at
#'     # least 10% of samples (prv_cut = 0.1), and samples with at least 10
#'     # reads (lib_cut = 10).
#'     x <- x[, colMeans(x > 0) >= 0.1, drop = FALSE]
#'     x <- x[rowSums(x) >= 10, , drop = FALSE]
#'
ancombc <- function (
    counts, formula, data, groups = NULL, margin = 1L, pseudocount = NULL,
    max_iters = 100L, tol = 1e-5, alpha = 0.05, p_adjust = 'holm',
    cpus = n_cpus() ) {

  ancombc_fit(environment(), v2 = FALSE)
}



#' @rdname ancombc
#' @export
ancombc2 <- function (
    counts, formula, data, groups = NULL, margin = 1L, pseudocount = 0,
    var_quantile = 0.05, max_iters = 100L, tol = 1e-5, alpha = 0.05,
    p_adjust = 'holm', cpus = n_cpus() ) {

  ancombc_fit(environment(), v2 = TRUE)
}



ancombc_fit <- function (env, v2) {

  what <- if (v2) 'ANCOM-BC2' else 'ANCOM-BC'

  if (missing_arg(env, 'formula'))
    stop('`formula` is required, such as `~ status + age`.', call. = FALSE)
  if (missing_arg(env, 'data'))
    stop('`data` is required: a data.frame of sample metadata.', call. = FALSE)

  if (v2) {
    pc <- env$pseudocount
    if (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc < 0)
      stop('`pseudocount` must be a single non-negative number.', call. = FALSE)
  }

  info <- log_ratio_input(what, env = env, zero_ok = TRUE, dense_out = FALSE)
  counts      <- env$counts
  margin      <- env$margin
  pseudocount <- env$pseudocount
  cpus        <- env$cpus

  max_iters <- env$max_iters
  if (!is.numeric(max_iters) || length(max_iters) != 1L || !is.finite(max_iters) ||
      max_iters < 1 || max_iters %% 1 != 0 || max_iters > .Machine$integer.max)
    stop('`max_iters` must be a single positive integer.', call. = FALSE)
  max_iters <- as.integer(max_iters)

  tol <- env$tol
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop('`tol` must be a single positive number.', call. = FALSE)
  tol <- as.double(tol)

  alpha <- env$alpha
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !(alpha > 0 && alpha < 1))
    stop('`alpha` must be a single number strictly between 0 and 1.', call. = FALSE)
  alpha <- as.double(alpha)

  var_quantile <- NULL
  if (v2) {
    var_quantile <- env$var_quantile
    if (!is.numeric(var_quantile) || length(var_quantile) != 1L ||
        is.na(var_quantile) || var_quantile < 0 || var_quantile > 1)
      stop('`var_quantile` must be a single number between 0 and 1.', call. = FALSE)
    var_quantile <- as.double(var_quantile)
  }

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)


  #--------------------------------------------------------
  # The design.
  #--------------------------------------------------------
  des <- ancombc_design(env$formula, env$data, info$samples, info$n, env$groups)
  X   <- des$X
  n   <- info$n
  D   <- info$D
  P   <- ncol(X)
  covariates <- colnames(X)

  # The global pseudoinverse and the inverse Gram matrix's row sums, from
  # one SVD (_lstsq_dense()), with numpy.linalg.pinv()'s cutoff.
  sv   <- svd(X)
  dinv <- ifelse(sv$d > 1e-15 * max(sv$d), 1 / sv$d, 0)
  Xinv <- sv$v %*% (dinv * t(sv$u))
  gsum <- as.vector(sv$v %*% (dinv^2 * colSums(sv$v)))

  dX <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'
  storage.mode(Xinv) <- 'double'

  res <- .Call(
    C_ancombc, counts, margin, as.double(pseudocount), v2, dX, Xinv, gsum,
    des$factors, des$intercept, max_iters, tol, cpus )

  bad <- which(res$n_used == 0L)
  if (length(bad))
    stop(
      'No estimable features remain for covariate ', name_list(covariates[bad]),
      '. Consider a simpler model or a positive `pseudocount`.', call. = FALSE )


  #--------------------------------------------------------
  # Bias correction, variances and degrees of freedom.
  #--------------------------------------------------------
  bias <- res$bias
  dimnames(bias) <- list(covariates, c('delta_em', 'delta_wls', 'var_delta'))

  if (!v2) {
    coef <- res$beta - rep(bias[, 'delta_em'], each = D)
    var  <- res$var
    dof  <- NULL
  }
  else {
    coef <- res$beta
    var  <- (sqrt(res$var) + rep(sqrt(bias[, 'var_delta']), each = D))^2
    if (var_quantile > 0)
      var <- var + rep(apply(var, 2L, np_quantile, var_quantile), each = D)

    if (!is.null(res$rank)) {
      n_obs <- n - res$n_missing
      dof   <- ifelse(n_obs > res$rank, as.double(n_obs - res$rank), NaN)
    }
    else {
      dof <- rep(if (n > P) as.double(n - P) else NaN, D)
    }
  }

  vcov <- res$vcov
  for (p in seq_len(P)) vcov[, p, p] <- var[, p]

  estimable <- NULL
  if (!is.null(res$valid))
    estimable <- matrix(res$valid, D, P)


  #--------------------------------------------------------
  # Statistics (_calc_statistics()).
  #--------------------------------------------------------
  lfc <- coef
  se  <- sqrt(var)
  if (!is.null(estimable)) {
    lfc[!estimable] <- NaN
    se[!estimable]  <- NaN
  }
  w <- lfc / se

  pval <- if (v2) suppressWarnings(2 * stats::pt(abs(w), dof, lower.tail = FALSE))
          else    2 * stats::pnorm(abs(w), lower.tail = FALSE)
  pval[is.na(pval)] <- NaN
  if (!is.null(estimable)) pval[!estimable] <- 1

  padj <- pval
  for (p in seq_len(P)) {
    ok <- !is.na(pval[, p])
    padj[ok, p] <- stats::p.adjust(pval[ok, p], p_adjust)
  }
  reject <- !is.na(padj) & padj <= alpha


  #--------------------------------------------------------
  # The table, feature by feature.
  #--------------------------------------------------------
  features <- if (is.null(info$features)) as.character(seq_len(D)) else info$features
  byrow    <- function (m) as.vector(t(m))

  out <- data.frame(
    feature    = rep(features, each = P),
    covariate  = rep(covariates, times = D),
    log_fc     = byrow(lfc),
    se         = byrow(se),
    w          = byrow(w),
    p_value    = byrow(pval),
    p_adjusted = byrow(padj),
    reject     = byrow(reject),
    stringsAsFactors = FALSE )

  dn <- list(features, covariates)
  dimnames(coef) <- dn
  dimnames(var)  <- dn
  dimnames(vcov) <- list(features, covariates, covariates)
  if (!is.null(estimable)) dimnames(estimable) <- dn
  if (!is.null(dof)) names(dof) <- features

  attr(out, 'method')      <- if (v2) 'ancombc2' else 'ancombc'
  attr(out, 'formula')     <- des$formula
  attr(out, 'groups')      <- des$groups
  attr(out, 'alpha')       <- alpha
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'max_iters')   <- max_iters
  attr(out, 'tol')         <- tol
  if (v2) attr(out, 'var_quantile') <- var_quantile
  attr(out, 'pseudocount') <- pseudocount
  attr(out, 'coef')        <- coef
  attr(out, 'var')         <- var
  attr(out, 'vcov')        <- vcov
  attr(out, 'dof')         <- dof
  attr(out, 'estimable')   <- estimable
  attr(out, 'design')      <- X
  attr(out, 'bias')        <- bias
  attr(out, 'n_iters')     <- stats::setNames(res$n_iters, covariates)
  if (v2) {
    sf <- res$theta
    if (!is.null(info$samples)) names(sf) <- info$samples
    attr(out, 'sample_fractions') <- sf
  }
  class(out) <- c('ecodive_ancombc', 'data.frame')

  out
}



# Whether a function's argument was left out (it has no default).
missing_arg <- function (env, arg) {
  eval(call('missing', as.name(arg)), envir = env)
}



# numpy.nanquantile(x, q) with numpy's default ('linear') method, including
# its form of the interpolation, which differs from quantile(type = 7) by a
# rounding step.
np_quantile <- function (x, q) {
  x <- sort(x[!is.na(x)], method = 'radix')
  m <- length(x)
  if (m == 0L) return (NaN)
  vi <- (m - 1) * q
  if (vi >= m - 1) return (x[[m]])
  lo <- floor(vi)
  a  <- x[[lo + 1]]
  b  <- x[[lo + 2]]
  t  <- vi - lo
  d  <- b - a
  if (t >= 0.5) b - d * (1 - t) else a + d * t
}



#------------------------------------------------------------------------------
# The design of ancombc() and ancombc2(): `data` aligned with the samples,
# `formula` turned into a model matrix with an intercept (unless removed)
# and treatment contrasts for every factor, as patsy's in scikit-bio.
#
# Returns a list:
#   X          n x P model matrix (assign and contrasts attributes kept).
#   formula    The formula used.
#   groups     The post-hoc term, or NULL.
#   intercept  0-based column of the intercept, or -1.
#   factors    One integer vector per factor with a main-effect term, for
#              C_ancombc: c(standard, rebase, L, level_col, codes) (see
#              src/ancombc.c).
#------------------------------------------------------------------------------
ancombc_design <- function (formula, data, labels, n, groups) {

  if (is.null(data))
    stop('`data` is required: a data.frame of sample metadata.', call. = FALSE)

  if (is.null(dim(data)) && is.atomic(data)) {
    nms  <- names(data)
    data <- data.frame(group = data, stringsAsFactors = FALSE)
    if (!is.null(nms) && !anyDuplicated(nms)) rownames(data) <- nms
  }
  if (!is.data.frame(data)) {
    if (!is.matrix(data))
      stop('`data` must be a data.frame, matrix, vector, or factor.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
  idx  <- match_rows(labels, n, dn, nrow(data), ref = '`counts`')
  data <- data[idx, , drop = FALSE]


  #--------------------------------------------------------
  # The formula.
  #--------------------------------------------------------
  env <- globalenv()
  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    txt <- if (grepl('~', formula, fixed = TRUE)) formula else paste('~', formula)
    f   <- tryCatch(str2lang(txt), error = function (e) NULL)
    if (!is.call(f) || !identical(f[[1]], as.name('~')))
      stop('`formula` must be a formula, such as `~ status + age`.', call. = FALSE)
    formula <- eval(f)
    environment(formula) <- env
  }
  if (!inherits(formula, 'formula'))
    stop('`formula` must be a formula, such as `~ status + age`.', call. = FALSE)
  if (length(formula) == 3L) {
    env     <- environment(formula)
    formula <- formula[-2L]
    environment(formula) <- if (is.null(env)) globalenv() else env
  }

  tt <- tryCatch(
    stats::terms(formula, data = data),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  vars    <- all.vars(tt)
  missing <- vars[!vars %in% names(data)]
  if (length(missing))
    stop('Not found in `data`: ', name_list(missing), call. = FALSE)

  na <- vars[vapply(data[vars], anyNA, logical(1))]
  if (length(na))
    stop('Missing values (NA) in `data` column', if (length(na) > 1) 's', ': ',
         name_list(na), call. = FALSE)

  inf <- vars[vapply(data[vars], function (v) is.numeric(v) && !all(is.finite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in `data` column', if (length(inf) > 1) 's', ': ',
         name_list(inf), call. = FALSE)


  #--------------------------------------------------------
  # Categorical variables: factors with treatment contrasts,
  # levels in their own order (factors) or sorted.
  #--------------------------------------------------------
  as_levels <- function (x) {
    if (is.factor(x))    return (factor(as.character(x), levels = levels(x)[levels(x) %in% as.character(x)]))
    if (is.logical(x))   return (factor(x, levels = c(FALSE, TRUE)[c(FALSE, TRUE) %in% x]))
    if (is.character(x)) return (factor(x, levels = sort(unique(x), method = 'radix')))
    x
  }
  for (v in vars) data[[v]] <- as_levels(data[[v]])

  mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
  for (i in seq_along(mf)) mf[[i]] <- as_levels(mf[[i]])

  is_fac <- vapply(mf, is.factor, logical(1))
  one    <- names(mf)[is_fac][vapply(mf[is_fac], nlevels, integer(1)) < 2L]
  if (length(one))
    stop('A categorical variable in `formula` has a single value: ', name_list(one), call. = FALSE)

  ctr <- stats::setNames(rep(list('contr.treatment'), sum(is_fac)), names(mf)[is_fac])
  X   <- tryCatch(
    stats::model.matrix(tt, mf, contrasts.arg = if (length(ctr)) ctr),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  P <- ncol(X)
  if (P == 0L)
    stop('`formula` has no terms.', call. = FALSE)
  if (n <= P)
    stop(
      'ANCOM-BC needs more samples than design columns: ', n, ' sample',
      if (n != 1L) 's', ', ', P, ' column', if (P != 1L) 's', '.', call. = FALSE )
  if (qr(X)$rank < P)
    warning(
      'The design is rank deficient: some of its columns are combinations of ',
      'others, so their coefficients are not unique.', call. = FALSE )


  #--------------------------------------------------------
  # Factors with a main-effect term (_r_fit_info()).
  #--------------------------------------------------------
  asg       <- attr(X, 'assign')
  fac       <- attr(tt, 'factors')
  ord       <- attr(tt, 'order')
  intercept <- if (attr(tt, 'intercept') == 1L) 0L else -1L
  factors   <- list()

  for (v in names(mf)[is_fac]) {
    if (!v %in% rownames(fac)) next
    main <- which(fac[v, ] > 0 & colSums(fac > 0) == 1L)
    if (length(main) != 1L) next
    cols <- which(asg == main) - 1L
    L    <- nlevels(mf[[v]])
    if      (length(cols) == L - 1L) { standard <- 1L; lcol <- c(-1L, cols) }
    else if (length(cols) == L)      { standard <- 0L; lcol <- cols }
    else next
    in_ia  <- any(fac[v, ] > 0 & ord > 1L)
    rebase <- as.integer(standard == 1L && !in_ia && intercept >= 0L)
    factors[[length(factors) + 1L]] <- as.integer(c(
      standard, rebase, L, lcol, as.integer(mf[[v]]) - 1L ))
  }


  #--------------------------------------------------------
  # The post-hoc grouping (_validate_grouping()).
  #--------------------------------------------------------
  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) != 1L || is.na(groups))
      stop('`groups` must be the name of a variable in `data`, or NULL.', call. = FALSE)
    if (!groups %in% names(data))
      stop('`groups` must be a column of `data`: ', name_list(groups), ' is not.', call. = FALSE)
    term <- match(groups, attr(tt, 'term.labels'))
    if (is.na(term))
      stop('`groups` must be a term in `formula`: ', name_list(groups), ' is not.', call. = FALSE)
    if (length(unique(data[[groups]])) < 3L)
      stop('`groups` must have at least three observed levels.', call. = FALSE)
    if (sum(asg == term) < 2L)
      stop(
        '`groups` (', name_list(groups), ') gives fewer than two design columns: ',
        'it must be categorical.', call. = FALSE )
  }

  list(
    X         = X,
    formula   = formula,
    groups    = groups,
    intercept = intercept,
    factors   = factors )
}
