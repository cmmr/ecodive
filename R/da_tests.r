# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Per-feature Wilcoxon, t and Linear Model Tests
#'
#' Tests each feature of a normalized count table for differential
#' abundance with a classical test: Wilcoxon's rank-sum test
#' (`da_wilcoxon()`), Welch's t-test (`da_ttest()`), or a linear model
#' (`da_lm()`), with the p-values adjusted for testing every feature. Each
#' gives, feature by feature, the results of base R's [stats::wilcox.test()],
#' [stats::t.test()] or [stats::lm()], computed in parallel C from a sparse
#' table without a dense copy.
#'
#' @inheritParams ancom
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Negative values are allowed with `norm = 'none'`, for a
#'        table that is already transformed.
#'
#' @param contrast  Which samples to compare, as `c(treatment, reference)`,
#'        two values of `groups`. A positive effect means more abundant in
#'        the treatment group. Use `NA` as the reference to compare the
#'        treatment group with all other samples pooled; with two named
#'        groups, samples in other groups are left out. With `NULL` (the
#'        default), the treatment is the **first** group (in sorted order,
#'        or a factor's first level) and the reference is all other samples,
#'        as for [dirmult_ttest()].
#'
#' @param formula  `da_lm()`: the model, a one-sided formula over the
#'        columns of `data`, such as `~ status + age`, or the same as a
#'        string. Anything [stats::model.matrix()] accepts may be used.
#'        Every factor (and character or logical column) is coded with
#'        treatment contrasts, its first level (in sorted order, or a
#'        factor's own order) the reference. `NULL` (the default) is
#'        `~ groups`.
#'
#' @param data  A data frame of sample metadata, with the variables of
#'        `groups` or `formula`. Its rows are matched to the samples of
#'        `counts` by name when both have names, else by position.
#'
#' @param groups  The grouping to compare: the name of a column of `data`,
#'        or a vector with one value per sample. `da_lm()`: with `formula`,
#'        the name of its term of interest (the one [plot()] shows by
#'        default); alone, the model `~ groups`.
#'
#' @param norm  How to transform each sample before testing: `'clr'` (the
#'        default), `'log'` (log relative abundance, "log-TSS"), `'rclr'`,
#'        `'percent'` (relative abundance, "TSS"), or any other choice of
#'        [normalize()]. With `'rclr'`, zeros are missing values, left out
#'        of their feature's test. See **Normalization**.
#'
#' @param pseudocount  For `norm = 'clr'` or `'log'`: the value added to
#'        every abundance before taking logs. `NULL` (the default) follows
#'        the shared rule: nothing is added when `counts` has no zeros, else
#'        half the smallest non-zero value, with a warning (see
#'        **Pseudocount** below). Ignored for other norms.
#'
#' @param var_equal  `da_ttest()`: `FALSE` (the default) for Welch's t-test,
#'        `TRUE` for Student's, with a pooled variance (`t.test()`'s
#'        `var.equal`).
#'
#' @param conf_level  The confidence level of the intervals (`ci_lower`,
#'        `ci_upper`). Strictly between 0 and 1. Default: `0.95`.
#'
#' @param alpha  The significance level of the adjusted p-values: `reject`
#'        is `p_adjusted <= alpha`. Strictly between 0 and 1. Default:
#'        `0.05`.
#'
#' @param p_adjust  How to adjust the p-values for testing every feature
#'        (separately for each covariate, in `da_lm()`): any method of
#'        [stats::p.adjust()] (`'holm'`, the default, as ecodive's other DA
#'        methods; `'BH'`, as the benchmarks below; `'hochberg'`,
#'        `'hommel'`, `'bonferroni'`, `'BY'`, `'fdr'`, `'none'`). `NaN`
#'        p-values are left out of the family.
#'
#' @return A data frame with a row per feature (`da_lm()`: per feature and
#'   covariate), in the order of `counts`, of class
#'   `c('ecodive_da_wilcoxon', 'data.frame')`,
#'   `c('ecodive_da_ttest', 'data.frame')`, or
#'   `c('ecodive_da_lm', 'data.frame')`. Every result starts with `feature`
#'   (the feature's name, or its index as text) and ends with `p_value`,
#'   `p_adjusted` (adjusted with `p_adjust`), and `reject` (`p_adjusted <=
#'   alpha`; `FALSE` where `NaN`). In between:
#'
#'   * `da_wilcoxon()`: `w`, Wilcoxon's \eqn{W} (`wilcox.test()`'s
#'     statistic: the treatment's rank sum minus its smallest possible
#'     value); `auc`, \eqn{W / (n_1 n_2)}, the probability that a
#'     treatment sample's value is above a reference sample's (ties count
#'     half); and `median_diff`, the treatment's median minus the
#'     reference's.
#'   * `da_ttest()`: `diff`, the treatment's mean minus the reference's;
#'     `ci_lower`, `ci_upper`, its confidence interval; `t`; `dof`, the
#'     degrees of freedom (Welch-Satterthwaite's, or \eqn{n_1 + n_2 - 2});
#'     `pooled_sd`, the pooled within-group standard deviation; and `mean`,
#'     the mean of both groups' values.
#'   * `da_lm()`: `covariate`, the design column (`'(Intercept)'` first,
#'     then as [stats::model.matrix()] names them); `estimate`, its
#'     coefficient; `se`; `t`; and `ci_lower`, `ci_upper`. Coefficients
#'     that can't be estimated (an aliased column, as `lm()`'s `NA`) are
#'     `NaN`.
#'
#'   Effects are differences in the units of `norm`: of natural logs for
#'   `'clr'`, `'rclr'` and `'log'` (log fold changes; `'clr'` and `'rclr'`
#'   relative to each sample's geometric mean), of proportions for
#'   `'percent'`.
#'
#'   Attributes: `method` (the function's name), `norm`, `pseudocount`,
#'   `alpha`, `p_adjust`, `log_base` (`exp(1)` for the log norms, `NA`
#'   otherwise), and for `da_wilcoxon()` and `da_ttest()`: `contrast` (the
#'   treatment and reference groups; the reference is `NA` when all other
#'   samples were pooled) and `n` (the number of samples on each side);
#'   `da_ttest()` and `da_lm()`: `conf_level`; `da_ttest()`: `var_equal`;
#'   `da_lm()`: `formula`, `groups` (the term of interest, or `NULL`),
#'   `design` (the model matrix), and per feature (named by feature):
#'   `dof` (the residual degrees of freedom), `mean` (the mean of its
#'   transformed values) and `resid_sd` (the residual standard deviation,
#'   `sigma()` of its `lm()`).
#'
#' @details
#' **The tests.** For each feature, the transformed values of the treatment
#' and reference samples are compared:
#'
#' * `da_wilcoxon()`: `wilcox.test(trt, ref)`, two-sided, with its defaults
#'   as of R 4.6.0: when both sides have fewer than 50 values, the exact
#'   p-value ([stats::pwilcox()] without ties; with ties, the exact
#'   distribution of the rank sum given the ties, as the Streitberg-Rohmel
#'   shift algorithm computes it, twice the smaller tail); otherwise the
#'   normal approximation, with the variance corrected for ties and a
#'   continuity correction. Before R 4.6.0, `wilcox.test()` used the normal
#'   approximation whenever there were ties, so its p-values differ there.
#'   Ties are equal values (with `norm = 'clr'`, all of a sample's zeros
#'   are one value, and zeros tie across samples only with `'percent'` or
#'   `'none'`). A feature whose values are all equal gets a p-value of 1.
#' * `da_ttest()`: `t.test(trt, ref, var.equal = var_equal, conf.level =
#'   conf_level)`.
#' * `da_lm()`: `lm(y ~ ...)` of the feature's values `y` on the design of
#'   `formula`, every coefficient tested against 0 with its t statistic
#'   (`summary.lm()`), with intervals as [stats::confint()]. A factor with
#'   several levels gets a test per level against the reference (as
#'   ANCOM-BC), not an F-test of the factor as a whole. With `norm =
#'   'rclr'`, each feature is fitted to the samples it was observed in, so
#'   its degrees of freedom (and, if a level has none of them, its aliased
#'   coefficients) are its own.
#'
#' Where base R stops with an error, these give `NaN` statistics and
#' p-values (never rejected) for that feature: a side with no values (with
#' `'rclr'`) or a single value (Welch's t-test, with a warning when a whole
#' side of the contrast has one sample), or, for the t-test, a feature
#' whose values are constant ("data are essentially constant").
#'
#' **Normalization.** The tests compare relative abundances, not absolute
#' ones: with `'clr'`, relative to each sample's geometric mean, so that a
#' change in a few abundant features shifts every other feature a little
#' the other way; with `'percent'` or `'log'`, relative to each sample's
#' total. [ancombc()], [ancom()] and [dirmult_ttest()] deal with that
#' differently. In the benchmarks of Wirbel et al. (2024) on realistic
#' simulations, the Wilcoxon test, t-test and linear models were among the
#' few methods to control false discoveries, and adjusting the linear model
#' for confounders (`formula = ~ group + confounder`) mitigated them; they
#' found the clr *lowered* the Wilcoxon test's sensitivity, while log
#' relative abundances (`norm = 'log'`) raised the t-test's and the linear
#' model's. Pelto et al. (2025) found the Wilcoxon test on relative
#' abundances (`norm = 'percent'`) and the t-test among the most
#' replicable methods across 69 gut datasets. The default (`'clr'`)
#' follows the compositional reading of the other DA methods; try the
#' others with [da_test()].
#'
#' **Computation.** The table is transformed per sample in C, without a
#' dense copy, and each feature is tested in parallel; each feature's
#' result depends only on its own values, so results are identical for any
#' number of `cpus`. The p-values, from the statistics, are those of base R
#' to rounding (about \eqn{10^{-12}}); `da_lm()` decomposes the design once
#' (a pivoted Householder QR, as `lm()`), and again per feature only for
#' features with missing values.
#'
#' @inheritSection pseudocount_section Pseudocount
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Streitberg, B., & Rohmel, J. (1986). Exact distributions for
#' permutation and rank tests: an introduction to some recently published
#' algorithms. *Statistical Software Newsletter*, 12, 10-17.
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
#' @seealso [da_test()], [ancombc()], [dirmult_ttest()], [ancom()],
#'   [normalize()], [replace_zeros()].
#'
#' @export
#' @examples
#'     # b1 and b2 are less abundant under treatment.
#'     x <- rbind(
#'       c( 20, 110, 100, 101, 100, 103, 104),
#'       c( 33, 110, 120, 100, 101, 100, 102),
#'       c( 12, 110, 100, 110, 100,  50,  90),
#'       c(202, 201,   9,  10,  10,  11,  11),
#'       c(200, 202,  10,  10,  13,  10,  10),
#'       c(203, 201,  14,  10,  10,  13,  12) )
#'     colnames(x) <- paste0('b', 1:7)
#'     grp <- rep(c('treatment', 'placebo'), each = 3)
#'
#'     da_ttest(x, grp, contrast = c('treatment', 'placebo'))
#'
#'     # With relative abundances (TSS), as the benchmarks.
#'     da_wilcoxon(x, grp, norm = 'percent', p_adjust = 'BH')
#'
#'     # A linear model, adjusted for a covariate.
#'     meta <- data.frame(trt = grp, age = c(30, 41, 52, 33, 45, 50))
#'     res <- da_lm(x, ~ trt + age, meta, groups = 'trt')
#'     res
#'     plot(res)
#'
da_wilcoxon <- function (
    counts, groups, data = NULL, contrast = NULL, margin = 1L, norm = 'clr',
    pseudocount = NULL, alpha = 0.05, p_adjust = 'holm', cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  da_two_sample(environment(), 'da_wilcoxon')
}


#' @rdname da_wilcoxon
#' @export
da_ttest <- function (
    counts, groups, data = NULL, contrast = NULL, margin = 1L, norm = 'clr',
    pseudocount = NULL, var_equal = FALSE, conf_level = 0.95, alpha = 0.05,
    p_adjust = 'holm', cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL
  da_two_sample(environment(), 'da_ttest')
}


#' @rdname da_wilcoxon
#' @export
da_lm <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    norm = 'clr', pseudocount = NULL, conf_level = 0.95, alpha = 0.05,
    p_adjust = 'holm', cpus = n_cpus() ) {

  da_lm_fit(environment())
}



#------------------------------------------------------------------------------
# Checks of the arguments the three tests share. Sets env's norm (to its
# code), pseudocount, margin, cpus, alpha, p_adjust (and conf_level when
# env has one). Returns list(D, n, features, samples, norm (the name),
# log_base).
#------------------------------------------------------------------------------
da_input <- function (env) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- env$counts
  margin <- env$margin

  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  norm <- env$norm
  if (!is.character(norm) || length(norm) != 1L || is.na(norm))
    stop(
      "`norm` must be one of 'clr', 'log', 'rclr', 'percent', 'chord', ",
      "'hellinger', 'binary', or 'none'.", call. = FALSE )
  norm_name <- match.arg(
    tolower(norm), c('none', 'percent', 'chord', 'hellinger', 'binary', 'clr', 'rclr', 'log') )

  pkg  <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (norm_name != 'none' && any(vals < 0))
    stop("`counts` must not contain negative values (except with norm = 'none').", call. = FALSE)
  remove('vals')

  validate_norm(env)
  validate_pseudocount(env)

  for (arg in c('alpha', 'conf_level')) {
    if (!exists(arg, envir = env, inherits = FALSE)) next
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !(v > 0 && v < 1))
      stop('`', arg, '` must be a single number strictly between 0 and 1.', call. = FALSE)
    env[[arg]] <- as.double(v)
  }
  env$p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  if (D < 1L) stop('`counts` must have at least one feature.', call. = FALSE)

  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])
  list(
    D        = D,
    n        = n,
    features = if (is.null(features)) as.character(seq_len(D)) else features,
    samples  = if (length(dn[[margin]])) as.character(dn[[margin]]),
    norm     = norm_name,
    log_base = if (norm_name %in% c('clr', 'rclr', 'log')) exp(1) else NA_real_ )
}



# stats::pwilcox(q, m, n, lower.tail), called once per (m, n): R 4.6.2's
# pwilcox() gives wrong values when vectorised over different sizes (as
# rclr's missing values make them), and right ones for scalar m and n.
da_pwilcox <- function (q, m, n, lower.tail) {
  out <- numeric(length(q))
  key <- paste(m, n)
  for (k in unique(key)) {
    i <- which(key == k)
    out[i] <- stats::pwilcox(q[i], m[[i[[1]]]], n[[i[[1]]]], lower.tail = lower.tail)
  }
  out
}


# p.adjust() over the non-NaN p-values; reject where p_adjusted <= alpha.
da_adjust <- function (p, p_adjust, alpha) {
  padj <- p
  ok   <- !is.na(p)
  padj[ok] <- stats::p.adjust(p[ok], p_adjust)
  list(p_adjusted = padj, reject = !is.na(padj) & padj <= alpha)
}



#------------------------------------------------------------------------------
# da_wilcoxon() and da_ttest().
#------------------------------------------------------------------------------
da_two_sample <- function (env, id) {

  info <- da_input(env)
  ttest <- id == 'da_ttest'

  if (ttest) {
    ve <- env$var_equal
    if (!isTRUE(ve) && !isFALSE(ve))
      stop('`var_equal` must be TRUE or FALSE.', call. = FALSE)
  }

  des <- test_group_values(env$groups, env$data, info$n, info$samples, ref = '`counts`')
  ct  <- dirmult_contrast(des$groups, env$contrast)
  n1  <- sum(ct$side == 1L)
  n2  <- sum(ct$side == 2L)

  if (ttest && !env$var_equal && (n1 < 2L || n2 < 2L))
    warning(
      'The ', if (n1 < 2L) 'treatment' else 'reference', ' side of `contrast` ',
      'has a single sample, so no variance: `t`, the p-values and the ',
      'intervals are NaN.', call. = FALSE )

  res <- .Call(
    C_da_tests, env$counts, env$margin, env$norm, as.double(env$pseudocount),
    ct$side, if (ttest) 2L else 1L, env$cpus )

  m1 <- res[, 1]
  m2 <- res[, 2]


  if (!ttest) {

    #--------------------------------------------------------
    # wilcox.test(x, y): exact without ties when both sides
    # are under 50, else normal with continuity correction.
    #--------------------------------------------------------
    W    <- res[, 3]
    ties <- res[, 4]
    none <- m1 == 0 | m2 == 0
    W[none] <- NaN

    p     <- rep(NaN, info$D)
    exact <- !none & m1 < 50 & m2 < 50 & ties == 0
    if (any(exact)) {
      w  <- W[exact]; a <- m1[exact]; b <- m2[exact]
      up <- w > a * b / 2
      pe <- ifelse(up, da_pwilcox(w - 1, a, b, FALSE), da_pwilcox(w, a, b, TRUE))
      p[exact] <- pmin(2 * pe, 1)
    }
    # With ties: exact too, from the conditional distribution given the
    # ties (computed in C), as wilcox.test() from R 4.6.0 on.
    tied <- !none & !exact & m1 < 50 & m2 < 50 & !is.na(res[, 7])
    if (any(tied))
      p[tied] <- pmin(1, 2 * pmin(res[tied, 7], res[tied, 8]))
    approx <- !none & !exact & !tied
    if (any(approx)) {
      w  <- W[approx]; a <- m1[approx]; b <- m2[approx]; tt <- ties[approx]
      z  <- w - a * b / 2
      sg <- sqrt((a * b / 12) * ((a + b + 1) - tt / ((a + b) * (a + b - 1))))
      z  <- (z - sign(z) * 0.5) / sg
      p[approx] <- 2 * pmin(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
    }

    adj <- da_adjust(p, env$p_adjust, env$alpha)

    out <- data.frame(
      feature     = info$features,
      w           = W,
      auc         = W / (m1 * m2),
      median_diff = res[, 5] - res[, 6],
      p_value     = p,
      p_adjusted  = adj$p_adjusted,
      reject      = adj$reject,
      stringsAsFactors = FALSE )
  }

  else {

    #--------------------------------------------------------
    # t.test(x, y, var.equal, conf.level).
    #--------------------------------------------------------
    mx <- res[, 3]; my <- res[, 4]; vx <- res[, 5]; vy <- res[, 6]
    q1 <- ifelse(m1 > 1, (m1 - 1) * vx, 0)
    q2 <- ifelse(m2 > 1, (m2 - 1) * vy, 0)

    if (env$var_equal) {
      df     <- m1 + m2 - 2
      stderr <- sqrt((q1 + q2) / df * (1 / m1 + 1 / m2))
    }
    else {
      sx     <- sqrt(vx / m1)
      sy     <- sqrt(vy / m2)
      stderr <- sqrt(sx^2 + sy^2)
      df     <- stderr^4 / (sx^4 / (m1 - 1) + sy^4 / (m2 - 1))
    }

    # t.test() stops with "data are essentially constant".
    const <- !is.na(stderr) & stderr < 10 * .Machine$double.eps * pmax(abs(mx), abs(my))
    bad   <- const | is.na(stderr) | stderr == 0 | m1 < 1 | m2 < 1 | !(df > 0)
    stderr[bad] <- NaN
    df[bad]     <- NaN

    tstat <- (mx - my) / stderr
    p     <- 2 * stats::pt(-abs(tstat), df)
    cint  <- stats::qt(1 - (1 - env$conf_level) / 2, df)

    adj <- da_adjust(p, env$p_adjust, env$alpha)

    pooled <- sqrt((q1 + q2) / (m1 + m2 - 2))
    pooled[!(m1 + m2 > 2)] <- NaN

    out <- data.frame(
      feature    = info$features,
      diff       = mx - my,
      ci_lower   = (tstat - cint) * stderr,
      ci_upper   = (tstat + cint) * stderr,
      t          = tstat,
      dof        = df,
      pooled_sd  = pooled,
      mean       = res[, 7],
      p_value    = p,
      p_adjusted = adj$p_adjusted,
      reject     = adj$reject,
      stringsAsFactors = FALSE )
  }

  for (j in seq_along(out)) if (is.double(out[[j]])) out[[j]][is.na(out[[j]])] <- NaN

  attr(out, 'method')      <- id
  attr(out, 'contrast')    <- ct$contrast
  attr(out, 'n')           <- c(treatment = n1, reference = n2)
  attr(out, 'norm')        <- info$norm
  attr(out, 'pseudocount') <- env$pseudocount
  if (ttest) {
    attr(out, 'var_equal')  <- env$var_equal
    attr(out, 'conf_level') <- env$conf_level
  }
  attr(out, 'alpha')       <- env$alpha
  attr(out, 'p_adjust')    <- env$p_adjust
  attr(out, 'log_base')    <- info$log_base
  class(out) <- c(paste0('ecodive_', id), 'data.frame')

  out
}



#------------------------------------------------------------------------------
# da_lm().
#------------------------------------------------------------------------------
da_lm_fit <- function (env) {

  info <- da_input(env)

  des <- da_design(env$formula, env$data, env$groups, info$samples, info$n, 'A linear model')
  X   <- des$X
  P   <- ncol(X)
  D   <- info$D
  covariates <- colnames(X)

  dX <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'

  # lm()'s tolerance for aliased columns; classical standard errors.
  res <- .Call(
    C_da_lm, env$counts, env$margin, env$norm, as.double(env$pseudocount),
    dX, 1e-7, 0L, env$cpus )

  est <- res$coef
  se  <- res$se
  dof <- res$dof
  dof[!(dof > 0)] <- NaN

  tstat <- est / se
  pval  <- 2 * stats::pt(abs(tstat), dof, lower.tail = FALSE)
  a     <- (1 - env$conf_level) / 2
  lo    <- est + se * stats::qt(a, dof)
  hi    <- est + se * stats::qt(1 - a, dof)

  padj   <- pval
  reject <- matrix(FALSE, D, P)
  for (p in seq_len(P)) {
    adj <- da_adjust(pval[, p], env$p_adjust, env$alpha)
    padj[, p]   <- adj$p_adjusted
    reject[, p] <- adj$reject
  }

  byrow <- function (m) { v <- as.vector(t(m)); v[is.na(v)] <- NaN; v }
  features <- info$features

  out <- data.frame(
    feature    = rep(features, each = P),
    covariate  = rep(covariates, times = D),
    estimate   = byrow(est),
    se         = byrow(se),
    t          = byrow(tstat),
    ci_lower   = byrow(lo),
    ci_upper   = byrow(hi),
    p_value    = byrow(pval),
    p_adjusted = byrow(padj),
    reject     = as.vector(t(reject)),
    stringsAsFactors = FALSE )

  resid_sd <- sqrt(res$rss / dof)
  resid_sd[is.na(resid_sd)] <- NaN

  attr(out, 'method')      <- 'da_lm'
  attr(out, 'formula')     <- des$formula
  attr(out, 'groups')      <- des$groups
  attr(out, 'norm')        <- info$norm
  attr(out, 'pseudocount') <- env$pseudocount
  attr(out, 'conf_level')  <- env$conf_level
  attr(out, 'alpha')       <- env$alpha
  attr(out, 'p_adjust')    <- env$p_adjust
  attr(out, 'log_base')    <- info$log_base
  attr(out, 'design')      <- X
  attr(out, 'dof')         <- stats::setNames(dof, features)
  attr(out, 'mean')        <- stats::setNames(res$mean, features)
  attr(out, 'resid_sd')    <- stats::setNames(resid_sd, features)
  class(out) <- c('ecodive_da_lm', 'data.frame')

  out
}
