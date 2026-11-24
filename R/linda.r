# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# LinDA, written from Zhou et al. (2022) and from the documentation and
# outputs of MicrobiomeStat's linda() (GPL-3; its code was not read). How
# each step was settled from outputs is recorded in plans/linda.md.



#' Linear Models for Differential Abundance (LinDA)
#'
#' Tests each feature for differential abundance with LinDA (Zhou et al.
#' 2022): an ordinary linear model of the feature's centred log2 ratios
#' (clr, in log2), whose coefficients are then corrected for the bias that
#' the composition puts on all of them, estimated as their mode across
#' features. Fast, and robust to strong compositional effects.
#'
#' @inheritParams da_wilcoxon
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Whole numbers are counts; anything else (such as the
#'        proportions of [replace_zeros()]) is taken as proportions. See
#'        **Zeros**.
#'
#' @param formula  The model, a one-sided formula over the columns of
#'        `data`, such as `~ status + age`, or the same as a string.
#'        Anything [stats::model.matrix()] accepts may be used. Every factor
#'        (and character or logical column) is coded with treatment
#'        contrasts, its first level (in sorted order, or a factor's own
#'        order) the reference. One random intercept may be added, as in
#'        lme4: `~ status + age + (1 | subject)`, for samples grouped by a
#'        column of `data` (see **Repeated measures**). `NULL` (the
#'        default) is `~ groups`.
#'
#' @param groups  The term of interest (a fixed effect): with `formula`, the name of one of
#'        its terms (the one [plot()] shows by default); alone, the name of a
#'        column of `data` or a vector with one value per sample, and the
#'        model is `~ groups`.
#'
#' @param zero_handling  How zeros of counts are replaced before taking
#'        logs: `'pseudocount'` (add `pseudocount` to every value),
#'        `'imputation'` (a zero becomes the sample's total divided by the
#'        largest total among the samples where the feature is zero), or
#'        `'adaptive'` (the default): imputation if the samples' totals are
#'        associated with the design (see **Zeros**), else the pseudocount.
#'        Ignored for proportions.
#'
#' @param pseudocount  For `zero_handling = 'pseudocount'` (or `'adaptive'`
#'        choosing it): the value added to every count. `NULL` (the default)
#'        follows the shared rule: nothing is added when the (winsorized)
#'        table has no zeros, else half its smallest non-zero value (0.5 for
#'        counts), with a warning. LinDA adds 0.5 always.
#'
#' @param corr_cut  For `zero_handling = 'adaptive'`: imputation is chosen
#'        when a design column's p-value for its association with the
#'        samples' log totals is at most `corr_cut`. Default `0.1`.
#'
#' @param outlier_pct  The proportion of each feature's values that is
#'        expected to be outliers: values above its `1 - outlier_pct`
#'        quantile are winsorized (brought down to it). `0` for none.
#'        Default `0.03`.
#'
#' @param standardize  `TRUE` (the default, as LinDA) to centre and scale
#'        every numeric variable of `formula` to unit standard deviation
#'        before fitting, so that their coefficients are per standard
#'        deviation; `FALSE` for coefficients per unit.
#'
#' @param p_adjust  How to adjust the p-values for testing every feature,
#'        separately for each covariate: any method of [stats::p.adjust()]
#'        (`'holm'`, the default, as ecodive's other DA methods; `'BH'`, as
#'        LinDA; `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BY'`, `'fdr'`,
#'        `'none'`). `NaN` p-values are left out of the family.
#'
#' @return A data frame of class `c('ecodive_linda', 'data.frame')` with a
#'   row per feature and covariate (every design column but the intercept),
#'   in the order of `counts`:
#'
#'   * `feature`: the feature's name, or its index as text.
#'   * `covariate`: the design column, as [stats::model.matrix()] names it
#'     (`'statusmoderate'`, `'age'`).
#'   * `log2_fc`: the bias-corrected coefficient, a log2 fold change in
#'     abundance (per standard deviation of a numeric covariate when
#'     `standardize = TRUE`).
#'   * `se`: its standard error.
#'   * `t`: `log2_fc / se`.
#'   * `dof`: the degrees of freedom, the number of samples minus the rank
#'     of the design; with a random intercept, Satterthwaite's, per feature
#'     and covariate.
#'   * `p_value`: \eqn{2 F_t(-|t|)}.
#'   * `p_adjusted`: adjusted with `p_adjust`, per covariate.
#'   * `reject`: `p_adjusted <= alpha` (`FALSE` where `NaN`).
#'
#'   Coefficients that can't be estimated (an aliased design column) are
#'   `NaN`. Attributes: `method` (`'linda'`), `formula`, `groups` (the term
#'   of interest, or `NULL`), `design` (the model matrix, standardized),
#'   `standardize`, `zero_handling` (what was done: `'pseudocount'`,
#'   `'imputation'`, or `'half_min'` for proportions), `pseudocount` (the
#'   value added; 0 when none was), `outlier_pct`, `alpha`, `p_adjust`,
#'   `log_base` (`2`), `bias` (the estimated bias of each covariate's
#'   coefficients: their mode, which was subtracted), `vcov` (a feature
#'   \eqn{\times} design column \eqn{\times} design column array, each
#'   feature's covariance matrix of the coefficients), and per feature
#'   (named by feature): `mean` (the mean of its centred log2 values, its
#'   mean clr abundance) and `resid_sd` (the residual standard deviation of
#'   its fit). With a random intercept, also `random` (its grouping
#'   variable) and per feature `random_sd` (the random intercepts' standard
#'   deviation; `vcov` and `resid_sd` are then the mixed model's). [plot()]
#'   draws volcano, effect, MA and forest plots of them
#'   ([plot.ecodive_linda()]).
#'
#' @details
#' **The model.** For each feature \eqn{i}, the centred log2 ratios
#' \eqn{W_{is} = \log_2 Y_{is} - \frac{1}{m} \sum_j \log_2 Y_{js}} of its
#' (winsorized, zero-free) abundances \eqn{Y} are fitted to the design by
#' least squares: \eqn{W_i = X \tilde\beta_i + \varepsilon_i}. If the
#' absolute abundances follow a linear model in \eqn{X}, each
#' \eqn{\tilde\beta_{ik}} is the feature's true log2 fold change plus one
#' bias, the same for every feature, that comes from the composition (the
#' clr's geometric mean moves with \eqn{X} too). When most features don't
#' change, the coefficients cluster around that bias, so LinDA estimates it
#' as the **mode** of \eqn{\{\tilde\beta_{ik}\}_i}, and reports
#' \eqn{\hat\beta_{ik} = \tilde\beta_{ik} - \text{mode}_k}, with the
#' standard error of the least squares fit and a t-test on \eqn{n -
#' \text{rank}(X)} degrees of freedom. The correction assumes that most
#' features are not differentially abundant (or that the changes are
#' balanced); see Zhou et al. (2022) for its asymptotic justification.
#'
#' The mode is that of a kernel density estimate: the mean shift
#' algorithm, with a Gaussian kernel and [stats::bw.nrd0()]'s bandwidth,
#' started at the shorth (the mean of the shortest half of the
#' coefficients), until a step is smaller than \eqn{\sqrt{\epsilon}} times
#' the estimate (at most 1000 steps). This is what MicrobiomeStat's
#' `linda()` does (it calls modeest's `meanshift()`).
#'
#' **Zeros.** Counts (whole numbers) are first winsorized: each feature's
#' proportions (count / sample total) above their \eqn{1 - }`outlier_pct`
#' quantile are brought down to it and turned back into counts (rounded).
#' Then zeros are replaced by `zero_handling`. With `'adaptive'`, the log
#' totals of the winsorized samples are regressed on the design's columns
#' (`lm(log(totals) ~ X)`); if any column's t-test p-value is at most
#' `corr_cut`, sequencing depth is associated with the design and a
#' pseudocount would bias the fold changes, so zeros are imputed instead
#' (Zhou et al. 2022, Methods). Values that aren't whole numbers are taken
#' as proportions: they are winsorized as they are, and each feature's
#' zeros become half its smallest non-zero value (a feature with no
#' non-zero value is an error).
#'
#' **Differences from LinDA.** The results are those of MicrobiomeStat
#' 1.4's `linda()` (to about \eqn{10^{-12}}) with `zero_handling` and
#' `pseudocount = 0.5` given, and `p_adjust = 'BH'`, except that:
#'
#' * MicrobiomeStat 1.4 decides `'adaptive'` as above, and reports its
#'   choice, but then adds the pseudocount whatever it chose; `linda()`
#'   does what it reports.
#' * LinDA always adds `pseudo.cnt = 0.5`; `linda()`'s `pseudocount = NULL`
#'   adds nothing to a table without zeros, and half the smallest non-zero
#'   value (0.5 for counts with ones) to one with zeros.
#' * LinDA's feature filters (`prev.filter`, `mean.abund.filter`,
#'   `max.abund.filter`) are not carried over: filter `counts` first, for
#'   example `counts[, colMeans(counts > 0) >= 0.1]` for a prevalence of
#'   10\%.
#' * Whether `counts` are counts or proportions is decided from the values
#'   (LinDA's `feature.dat.type`).
#' * `p_adjust` defaults to `'holm'` (LinDA: `'BH'`).
#' * Factor levels are sorted in the C locale (LinDA: in the session's).
#' * Of LinDA's mixed models, one random intercept, `(1 | id)`; random
#'   slopes, several random terms and nested terms are not supported.
#' * Mixed models are fitted at the exact REML optimum, and their degrees
#'   of freedom from exact derivatives (see **Repeated measures**); lme4's
#'   optimizer stops a little short of the optimum, and lmerTest takes
#'   derivatives numerically, so these differ from LinDA's by about
#'   \eqn{10^{-6}} (\eqn{10^{-4}} in places where the likelihood is
#'   very flat).
#' * The result is one data frame (LinDA: a list per covariate), without
#'   `baseMean` (the `mean` attribute gives each feature's mean clr
#'   abundance instead).
#'
#' **Repeated measures.** With `(1 | id)` in `formula`, each feature's
#' clr values get a linear mixed model, with a random intercept for each
#' value of `id` (a subject, say), fitted by restricted maximum likelihood
#' (REML) as lme4's `lmer()` does (Bates et al. 2015), and the same mode
#' correction (Zhou et al. 2022, Remark 1). Coefficients are tested with
#' Satterthwaite's degrees of freedom, as by lmerTest (Kuznetsova et al.
#' 2017), from the variance parameters' REML covariance. The optimum is
#' found exactly (the root of the REML criterion's derivative in the ratio
#' of the two standard deviations, by Brent's method then bracketing), and
#' every derivative is analytic. A random intercept whose variance is
#' estimated as 0 (lme4's "boundary (singular) fit") gives the fixed-effects
#' fit, with \eqn{n - p} degrees of freedom. With `zero_handling =
#' 'adaptive'`, the samples' log totals get the same mixed model.
#'
#' **Computation.** The table is read by feature from its sparse form,
#' winsorized and transformed per feature in parallel C, and every feature
#' is fitted to one QR decomposition of the design; the per-sample means of
#' the clr are summed in a fixed order, so results are identical for any
#' number of `cpus`.
#'
#' @section Same names in other packages:
#' MicrobiomeStat and LinDA have `linda()`, the same method.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Zhou, H., He, K., Chen, J., & Zhang, X. (2022). LinDA: linear models for
#' differential abundance analysis of microbiome compositional data.
#' *Genome Biology*, 23, 95. \doi{10.1186/s13059-022-02655-5}
#'
#' Yang, L., & Chen, J. (2023). Benchmarking differential abundance analysis
#' methods for correlated microbiome sequencing data. *Briefings in
#' Bioinformatics*, 24, bbac607. \doi{10.1093/bib/bbac607}
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting
#' linear mixed-effects models using lme4. *Journal of Statistical
#' Software*, 67(1), 1-48. \doi{10.18637/jss.v067.i01}
#'
#' Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017).
#' lmerTest package: tests in linear mixed effects models. *Journal of
#' Statistical Software*, 82(13), 1-26. \doi{10.18637/jss.v082.i13}
#'
#' Fukunaga, K., & Hostetler, L. (1975). The estimation of the gradient of
#' a density function, with applications in pattern recognition. *IEEE
#' Transactions on Information Theory*, 21, 32-40.
#' \doi{10.1109/TIT.1975.1055330}
#'
#' @seealso [da_test()], [ancombc()], [da_lm()], [replace_zeros()].
#'
#' @export
#' @examples
#'     # 30 samples, 40 features; f1 to f4 are more abundant in group b.
#'     set.seed(1)
#'     x <- matrix(rpois(30 * 40, 50), 30, 40)
#'     colnames(x) <- paste0('f', 1:40)
#'     meta <- data.frame(
#'       grp = rep(c('a', 'b'), each = 15),
#'       age = round(runif(30, 20, 60)) )
#'     x[16:30, 1:4] <- x[16:30, 1:4] * 4
#'
#'     res <- linda(x, ~ grp + age, meta, groups = 'grp')
#'     res
#'     plot(res)
#'
#'     # The same, as LinDA's defaults would give it.
#'     linda(x, ~ grp + age, meta, zero_handling = 'pseudocount',
#'           pseudocount = 0.5, p_adjust = 'BH')
#'
#'     # Two samples from each of 15 subjects: a random intercept per subject.
#'     meta$subject <- rep(sprintf('s%02d', 1:15), 2)
#'     linda(x, ~ grp + age + (1 | subject), meta, groups = 'grp')
#'
linda <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    zero_handling = 'adaptive', pseudocount = NULL, corr_cut = 0.1,
    outlier_pct = 0.03, standardize = TRUE, alpha = 0.05, p_adjust = 'holm',
    cpus = n_cpus() ) {

  linda_fit(environment())
}



linda_fit <- function (env) {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  pc <- env$pseudocount
  env$norm <- 'none'
  info <- da_input(env)
  env$pseudocount <- pc

  counts <- env$counts
  vals   <- stored_values(counts)
  if (any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)
  if (info$D < 2L)
    stop('`counts` must have at least two features.', call. = FALSE)
  is_count <- all(vals == round(vals))
  remove('vals')

  zh <- env$zero_handling
  if (!is.character(zh) || length(zh) != 1L || is.na(zh))
    stop("`zero_handling` must be 'adaptive', 'pseudocount', or 'imputation'.", call. = FALSE)
  zh <- match.arg(tolower(zh), c('adaptive', 'pseudocount', 'imputation'))

  for (arg in c('corr_cut', 'outlier_pct')) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v >= 1)
      stop('`', arg, '` must be a single number in [0, 1).', call. = FALSE)
  }
  corr_cut    <- as.double(env$corr_cut)
  outlier_pct <- as.double(env$outlier_pct)

  standardize <- env$standardize
  if (!isTRUE(standardize) && !isFALSE(standardize))
    stop('`standardize` must be TRUE or FALSE.', call. = FALSE)

  if (!is.null(pc) && (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc < 0))
    stop('`pseudocount` must be NULL or a single non-negative number.', call. = FALSE)


  #--------------------------------------------------------
  # The design.
  #--------------------------------------------------------
  rnd <- linda_random(env$formula)
  des <- da_design(
    rnd$fixed, env$data, env$groups, info$samples, info$n, 'LinDA',
    standardize = standardize )
  X <- des$X
  P <- ncol(X)
  keep <- setdiff(seq_len(P), des$intercept + 1L)
  if (!length(keep))
    stop('`formula` must have at least one term besides the intercept.', call. = FALSE)

  dX <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'

  # A random intercept: each sample's group.
  mixed <- !is.null(rnd$random)
  codes <- NULL
  if (mixed) {
    v <- rnd$random
    if (!v %in% names(des$data))
      stop('Not found in `data`: ', name_list(v), call. = FALSE)
    idv <- des$data[[v]]
    if (anyNA(idv))
      stop('Missing values (NA) in `data` column ', name_list(v), call. = FALSE)
    codes <- as.integer(factor(as.character(idv)))
    J <- max(codes)
    if (J < 2L)
      stop('The random intercept\'s `', v, '` must have at least two values.', call. = FALSE)
    if (J >= info$n)
      stop(
        'The random intercept\'s `', v, '` must have fewer values than there ',
        'are samples (', J, ' values, ', info$n, ' samples).', call. = FALSE )
    if (qr(dX)$rank < P)
      stop('A mixed model needs a design of full rank.', call. = FALSE)
  }


  #--------------------------------------------------------
  # Winsorization, and how zeros are replaced.
  #--------------------------------------------------------
  prep <- .Call(
    C_linda_prep, counts, env$margin, is_count, outlier_pct, env$cpus )

  pc_used <- 0
  if (is_count) {

    if (any(prep$depth <= 0))
      stop('Every sample must have a positive total count.', call. = FALSE)

    if (zh == 'adaptive') {
      if (any(prep$depth_w <= 0))
        stop(
          'A sample has no counts left after winsorization, so the ',
          'adaptive zero handling can\'t test its total: give ',
          '`zero_handling`.', call. = FALSE )
      Xk <- cbind(1, dX[, keep, drop = FALSE])
      if (mixed) {
        fit <- .Call(C_feature_lmm, Xk, matrix(log(prep$depth_w)), codes, NULL)
        p   <- (2 * stats::pt(-abs(fit$coef / fit$se), fit$df))[-1L]
      }
      else {
        p <- linda_lm_p(stats::lm.fit(Xk, log(prep$depth_w)))
      }
      zh  <- if (any(p <= corr_cut, na.rm = TRUE)) 'imputation' else 'pseudocount'
    }

    if (zh == 'imputation') {
      zero_mode <- 2L
      if (any(prep$depth_w <= 0))
        stop(
          'A sample has no counts left after winsorization, so its zeros ',
          'can\'t be imputed.', call. = FALSE )
    }
    else {
      zero_mode <- 1L
      has_zero  <- any(prep$n_zero > 0L)
      if (is.null(pc)) {
        pc_used <- if (has_zero)
          log_pseudocount(matrix(c(0, min(prep$min_nz))), NULL, 'LinDA', constant_ok = TRUE)
          else 0
      }
      else {
        if (pc == 0 && has_zero)
          stop('`counts` has zeros, and `pseudocount = 0` leaves their logs undefined.', call. = FALSE)
        pc_used <- as.double(pc)
      }
    }
  }

  else {
    zh <- 'half_min'
    zero_mode <- 3L
    none <- which(!is.finite(prep$min_nz))
    if (length(none))
      stop(
        'Feature', if (length(none) > 1) 's', ' with no non-zero value ',
        '(after winsorization): ', name_list(info$features[none]),
        '. Their zeros can\'t be replaced; leave them out.', call. = FALSE )
  }


  #--------------------------------------------------------
  # The fits, and the mode of each coefficient.
  #--------------------------------------------------------
  res <- .Call(
    C_linda_fit, counts, env$margin, is_count, outlier_pct, prep$depth,
    prep$depth_w, zero_mode, as.double(pc_used), dX, 1e-7, codes, env$cpus )

  D   <- info$D
  dof <- if (mixed) res$df[, keep, drop = FALSE] else info$n - res$rank
  B   <- res$coef[, keep, drop = FALSE]

  start <- h <- rep(NaN, length(keep))
  for (j in seq_along(keep)) {
    b <- B[, j]
    b <- b[!is.na(b)]
    if (length(b) < 2L) next
    start[[j]] <- linda_shorth(b)
    h[[j]]     <- stats::bw.nrd0(b)
  }
  mode <- .Call(C_linda_mode, B, start, h, sqrt(.Machine$double.eps), 1000L)$mode
  covariates <- colnames(X)[keep]
  names(mode) <- covariates

  est <- sweep(B, 2L, mode)
  se  <- res$se[, keep, drop = FALSE]
  if (!mixed && !(dof > 0)) se[] <- NaN
  tstat <- est / se
  pval  <- 2 * stats::pt(-abs(tstat), dof)

  K <- length(keep)
  padj   <- pval
  reject <- matrix(FALSE, D, K)
  for (k in seq_len(K)) {
    adj <- da_adjust(pval[, k], env$p_adjust, env$alpha)
    padj[, k]   <- adj$p_adjusted
    reject[, k] <- adj$reject
  }

  byrow    <- function (m) { v <- as.vector(t(m)); v[is.na(v)] <- NaN; v }
  features <- info$features

  out <- data.frame(
    feature    = rep(features, each = K),
    covariate  = rep(covariates, times = D),
    log2_fc    = byrow(est),
    se         = byrow(se),
    t          = byrow(tstat),
    dof        = if (mixed) byrow(dof) else if (dof > 0) as.double(dof) else NaN,
    p_value    = byrow(pval),
    p_adjusted = byrow(padj),
    reject     = as.vector(t(reject)),
    stringsAsFactors = FALSE )


  #--------------------------------------------------------
  # Each feature's covariance matrix: s^2 (X'X)^-1, or the
  # mixed model's s^2 (X'H^-1 X)^-1.
  #--------------------------------------------------------
  if (mixed) {
    vcov     <- res$vcov
    dimnames(vcov) <- list(features, colnames(X), colnames(X))
    resid_sd <- res$rss
  }
  else {
    s2   <- if (dof > 0) res$rss / dof else rep(NaN, D)
    qx   <- qr(dX, tol = 1e-7, LAPACK = FALSE)
    r    <- qx$rank
    xtxi <- matrix(NaN, P, P)
    if (r > 0) {
      piv <- qx$pivot[seq_len(r)]
      xtxi[piv, piv] <- chol2inv(qx$qr[seq_len(r), seq_len(r), drop = FALSE])
    }
    vcov <- array(rep(s2, times = P * P) * rep(as.vector(xtxi), each = D),
                  dim = c(D, P, P), dimnames = list(features, colnames(X), colnames(X)))
    resid_sd <- sqrt(s2)
  }
  resid_sd[is.na(resid_sd)] <- NaN

  attr(out, 'method')        <- 'linda'
  attr(out, 'formula')       <- if (mixed) linda_formula(des$formula, rnd$random) else des$formula
  attr(out, 'groups')        <- des$groups
  attr(out, 'design')        <- X
  attr(out, 'standardize')   <- standardize
  attr(out, 'zero_handling') <- zh
  attr(out, 'pseudocount')   <- pc_used
  attr(out, 'outlier_pct')   <- outlier_pct
  attr(out, 'alpha')         <- env$alpha
  attr(out, 'p_adjust')      <- env$p_adjust
  attr(out, 'log_base')      <- 2
  attr(out, 'bias')          <- mode
  attr(out, 'vcov')          <- vcov
  attr(out, 'mean')          <- stats::setNames(res$mean, features)
  attr(out, 'resid_sd')      <- stats::setNames(resid_sd, features)
  if (mixed) {
    random_sd <- res$theta * resid_sd
    random_sd[is.na(random_sd)] <- NaN
    attr(out, 'random')    <- rnd$random
    attr(out, 'random_sd') <- stats::setNames(random_sd, features)
  }
  class(out) <- c('ecodive_linda', 'data.frame')

  out
}



# The shorth: the mean of the shortest interval holding ceiling(n / 2) of
# the sorted values; among equally short ones, the one starting at the
# (rounded down) mean of their positions. MicrobiomeStat's start for the
# mean shift (modeest's shorth(), whose results this matches).
linda_shorth <- function (x) {
  y <- sort(x)
  n <- length(y)
  k <- ceiling(n / 2) - 1
  w <- y[(k + 1):n] - y[1:(n - k)]
  i <- which(w == min(w))
  j <- floor(mean(i))
  mean(y[j:(j + k)])
}


# Splits a random intercept, `(1 | id)`, from `formula`: list(fixed (the
# formula without it, or `formula` as given), random (the grouping
# variable's name, or NULL)).
linda_random <- function (formula) {

  f <- formula
  if (is.character(f) && length(f) == 1L && !is.na(f)) {
    txt <- if (grepl('~', f, fixed = TRUE)) f else paste('~', f)
    f   <- tryCatch(eval(str2lang(txt)), error = function (e) NULL)
    if (inherits(f, 'formula')) environment(f) <- globalenv()
  }
  if (!inherits(f, 'formula')) return (list(fixed = formula, random = NULL))

  rhs   <- f[[length(f)]]
  split <- function (e)
    if (is.call(e) && identical(e[[1L]], as.name('+')) && length(e) == 3L)
      c(split(e[[2L]]), split(e[[3L]]))
    else list(e)
  parts <- split(rhs)
  bar   <- vapply(parts, function (e) any(all.names(e) %in% c('|', '||')), logical(1))
  if (!any(bar)) return (list(fixed = formula, random = NULL))

  ok <- function (e) {
    if (!is.call(e) || !identical(e[[1L]], as.name('('))) return (FALSE)
    e <- e[[2L]]
    is.call(e) && identical(e[[1L]], as.name('|')) && length(e) == 3L &&
      identical(e[[2L]], 1) && is.name(e[[3L]])
  }
  if (!all(vapply(parts[bar], ok, logical(1))) || sum(bar) > 1L)
    stop(
      'Of random effects, `formula` may have one random intercept, such as ',
      '`(1 | subject)`, of a column of `data`.', call. = FALSE )
  if (all(bar))
    stop('`formula` must have at least one fixed effect besides `(1 | ...)`.', call. = FALSE)

  fixed <- Reduce(function (a, b) call('+', a, b), parts[!bar])
  out   <- stats::as.formula(call('~', fixed), env = environment(f))
  list(fixed = out, random = as.character(parts[bar][[1L]][[2L]][[3L]]))
}


# The fixed formula with the random intercept added back: ~ g + (1 | id).
linda_formula <- function (fixed, random) {
  rhs <- call('+', fixed[[length(fixed)]], call('(', call('|', 1, as.name(random))))
  stats::as.formula(call('~', rhs), env = environment(fixed))
}


# The t-test p-values of the non-intercept coefficients of lm.fit()'s fit
# (summary.lm()'s), NaN where aliased or without residual degrees of freedom.
linda_lm_p <- function (fit) {
  r   <- fit$rank
  dfr <- length(fit$residuals) - r
  p   <- length(fit$coefficients)
  out <- rep(NaN, p - 1L)
  if (dfr <= 0 || r < 1) return (out)
  piv <- fit$qr$pivot[seq_len(r)]
  rinv <- chol2inv(fit$qr$qr[seq_len(r), seq_len(r), drop = FALSE])
  s2  <- sum(fit$residuals^2) / dfr
  se  <- sqrt(diag(rinv) * s2)
  tt  <- fit$coefficients[piv] / se
  pv  <- rep(NaN, p)
  pv[piv] <- 2 * stats::pt(-abs(tt), dfr)
  pv[-1L]
}



#' Plot a LinDA Result
#'
#' Views of a [linda()] result for one covariate, each feature a point (or
#' a row), with the features declared differentially abundant (`reject =
#' TRUE`) highlighted and the most significant of them named:
#'
#' * `"volcano"` (the default): the log2 fold change against the adjusted
#'   p-value.
#' * `"effect"`: the log2 fold change against the residual standard
#'   deviation (for a two-level comparison).
#' * `"ma"`: the log2 fold change against the mean clr abundance.
#' * `"forest"`: the log2 fold changes and their confidence intervals of
#'   the `n_features` most significant features.
#' * `"rank"`: every feature's log2 fold change, sorted (see **Rank
#'   plots**).
#'
#' @param x  The result of [linda()].
#'
#' @param covariate  Which covariate's log2 fold changes to plot, as named
#'        in the result's `covariate` column (such as `"statusmoderate"` or
#'        `"age"`). Defaults to the first column of the term of interest
#'        (`groups`), or else to the first covariate.
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
#' @param conf_level  For `"forest"`: the confidence level of the
#'        intervals. Default `0.95`.
#'
#' @inheritParams plot.ecodive_ancombc
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `sub`, `xlab`, `xlim` or `cex`.
#'
#' @inheritSection plot.ecodive_ancombc Rank plots
#'
#' @details
#' Log fold changes are in log2, as LinDA estimates them (the result's
#' `log_base`), bias-corrected.
#'
#' * **Volcano plot.** The y axis is \eqn{-\log_{10}} of `p_adjusted` (an
#'   adjusted p-value of 0 is drawn at the smallest positive double), with
#'   a dashed line at `alpha`.
#' * **Effect plot.** The x axis is the feature's residual standard
#'   deviation (the `resid_sd` attribute, log2); the dashed lines are
#'   \eqn{y = \pm x}. Only for a covariate whose design column is 0 or 1 (a
#'   factor's level against its reference).
#' * **MA plot.** The x axis is the feature's mean centred log2 ratio (the
#'   `mean` attribute).
#' * **Forest plot.** Bars over \eqn{\mathrm{log2\_fc} \pm q\,\mathrm{se}},
#'   with \eqn{q} the \eqn{(1 + }`conf_level`\eqn{)/2} quantile of the t
#'   distribution with `dof` degrees of freedom, as the p-values; they
#'   aren't adjusted for testing many features.
#'
#' Features whose values are `NaN` are left out. The subtitle names the
#' covariate and counts the features rejected.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     x <- matrix(rpois(20 * 30, 40), 20, 30)
#'     x[11:20, 1:3] <- x[11:20, 1:3] * 5
#'     res <- linda(x, groups = rep(c('a', 'b'), each = 10))
#'     plot(res)
#'     plot(res, type = 'effect')
#'     plot(res, type = 'forest', n_features = 10)
#'     plot(res, type = 'rank', numerator = 1:3)
plot.ecodive_linda <- function (
    x, covariate = NULL, type = c('volcano', 'effect', 'ma', 'forest', 'rank'),
    labels = 10, n_features = 30, conf_level = 0.95, numerator = NULL,
    denominator = NULL, ... ) {

  type <- test_plot_type(type, c('volcano', 'effect', 'ma', 'forest', 'rank'))

  covariate <- da_rank_covariate(x, covariate)
  if (type == 'effect') da_plot_binary(x, covariate)

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)

  r <- x[x$covariate == covariate, , drop = FALSE]
  q <- suppressWarnings(stats::qt((1 + conf_level) / 2, r$dof))

  da_plot(
    x, type, labels, n_features, list(...),
    est = r$log2_fc, p = r$p_adjusted, reject = r$reject, feature = r$feature,
    sd = attr(x, 'resid_sd', exact = TRUE), mean = attr(x, 'mean', exact = TRUE),
    lo = r$log2_fc - q * r$se, hi = r$log2_fc + q * r$se,
    what     = covariate,
    est_lab  = lfc_label(x),
    sd_lab   = 'Residual standard deviation (log2)',
    mean_lab = 'Mean abundance (clr, log2)',
    interval = sprintf('unadjusted %s confidence intervals', plot_pct(conf_level)),
    num = numerator, den = denominator )

  invisible(x)
}
