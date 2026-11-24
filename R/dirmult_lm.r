# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# dirmult_lm() is a port of the linear-model path of aldex() from ALDEx3
# 1.4.0 (aldex(), aldex.sampler(), aldex.pvals(), fflm(), clr.sm(),
# tss.sm(), sample.sm(), rLogDirichlet() and summary.aldex()), Copyright (c)
# 2026 ALDEx3 authors, under the MIT license given in src/dirmult_lm.c and
# inst/COPYRIGHTS; restructured for ecodive (see plans/dirmult_lm.md).


#' Dirichlet-multinomial linear models with scale models (ALDEx3)
#'
#' Fits a linear model to every feature's log abundance, carrying the
#' uncertainty of low counts and zeros (Dirichlet posterior draws) and of
#' the samples' unmeasured total loads (a **scale model**) through to the
#' estimates and p-values: the linear models of ALDEx3 (McGovern &
#' Silverman 2026; Nixon et al. 2025). It is the regression form of
#' [dirmult_ttest()]: any design, with covariates, numeric variables and an
#' optional random intercept, and a scale model that can use measured loads.
#'
#' @inheritParams da_wilcoxon
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Raw counts are recommended: the Dirichlet posterior's
#'        spread depends on the counts' size, so proportions or rescaled
#'        values give wrong uncertainties.
#'
#' @param formula  The model, a one-sided formula over the columns of
#'        `data`, such as `~ status + age`, or the same as a string.
#'        Anything [stats::model.matrix()] accepts may be used. Every factor
#'        (and character or logical column) is coded with treatment
#'        contrasts, its first level (in sorted order, or a factor's own
#'        order) the reference. One random intercept may be added, as in
#'        lme4: `~ status + age + (1 | subject)` (see **Repeated
#'        measures**). `NULL` (the default) is `~ groups`.
#'
#' @param groups  The term of interest: with `formula`, the name of one of
#'        its terms (the one [plot()] shows by default); alone, the name of a
#'        column of `data` or a vector with one value per sample, and the
#'        model is `~ groups`.
#'
#' @param pseudocount  The Dirichlet prior: a positive number added to every
#'        value, zero or not, before drawing. Default: `0.5`, as ALDEx3 and
#'        [dirmult_ttest()].
#'
#' @param scale_model  What is assumed about the samples' total loads,
#'        without measurements: `'clr'` (the default; ALDEx3's `clr.sm()`)
#'        centres them where the clr puts them, so that the typical feature
#'        doesn't change; `'tss'` (`tss.sm()`) assumes no systematic
#'        difference in total load, so that relative abundances
#'        (proportions) are compared. Either way `scale_sd` says how
#'        uncertain that is. See **Scale models**.
#'
#' @param scale_sd  The uncertainty of the scale model: the standard
#'        deviation, in log2 units, of a random shift in log2 load drawn for
#'        each design column in each draw (ALDEx3's `gamma`). Default:
#'        `0.5`, as ALDEx3's `clr.sm()`. `0` takes the assumption as exact.
#'
#' @param load  Measured (or predicted) total loads, one per sample (flow
#'        cytometry, qPCR, spike-ins, or the load predicted from the profile
#'        as in Nishijima et al. 2025): numbers greater than 0, in any unit,
#'        matched to the samples by name when both are named; or the name of
#'        a column of `data`. Replaces `scale_model` and `scale_sd` (ALDEx3's
#'        `sample.sm()`). Default: `NULL` (none).
#'
#' @param load_sd  With `load`: the standard deviation of each sample's log2
#'        load's error, one number for all samples or one per sample (or a
#'        column of `data`); `0` takes the loads as exact. Required with
#'        `load`.
#'
#' @param n_draws  The number of Monte Carlo draws from the Dirichlet
#'        posterior and the scale model (ALDEx3's `nsample`). Default:
#'        `1000`, the least ALDEx3 recommends with a scale model (its
#'        default is 2,000).
#'
#' @param se_type  The standard errors of each draw's coefficients:
#'        `'hc3'` (the default, as ALDEx3) or `'hc0'`, heteroskedasticity-
#'        consistent ones (Long & Ervin 2000), or `'classical'` (ordinary
#'        least squares). With a random intercept the standard errors are
#'        the mixed model's, and only `'classical'` may be given.
#'
#' @param alpha  The significance level for `reject`, applied to
#'        `p_adjusted`. Strictly between 0 and 1. Default: `0.05`.
#'
#' @param p_adjust  How to adjust each draw's p-values for testing every
#'        feature, separately for each covariate: any method of
#'        [stats::p.adjust()] (`'holm'`, the default, as ecodive's other DA
#'        methods; `'BH'`, as ALDEx3; `'hochberg'`, `'hommel'`,
#'        `'bonferroni'`, `'BY'`, `'fdr'`, `'none'`).
#'
#' @param seed  An integer seed for the random draws. The same `seed` gives
#'        the same results, for any number of `cpus`.
#'
#' @return A data frame of class `c('ecodive_dirmult_lm', 'data.frame')`
#'   with a row per feature and design column (the intercept included), in
#'   the order of `counts`:
#'
#'   * `feature`: the feature's name, or its index as text.
#'   * `covariate`: the design column, as [stats::model.matrix()] names it
#'     (`"(Intercept)"`, `"statusmoderate"`, `"age"`, ...).
#'   * `estimate`: the coefficient, in log2 units (a log2 fold change for a
#'     factor's level against its reference, per unit for a numeric
#'     variable), averaged over the draws.
#'   * `se`: its standard error, averaged over the draws.
#'   * `dof`: the degrees of freedom of each draw's t-tests: \eqn{n - P}
#'     (samples less design columns); with a random intercept,
#'     Satterthwaite's, averaged over the draws.
#'   * `p_value`: the two-sided p-value, from the draws' one-sided ones (see
#'     **Details**).
#'   * `p_adjusted`: the same from the draws' adjusted p-values.
#'   * `reject`: `TRUE` when `p_adjusted <= alpha`.
#'
#'   Attributes: `method` (`"dirmult_lm"`), `formula`, `groups`, `design`
#'   (the model matrix), `random` (the random intercept's variable, if
#'   any), `scale_model` (`"clr"`, `"tss"` or `"load"`), `scale_sd` (`NA`
#'   with `load`), `load_sd` (per sample, with `load`), `se_type`,
#'   `n_draws`, `pseudocount`, `alpha`, `p_adjust`, `seed`, `log_base` (`2`),
#'   and per feature `mean` (its mean clr value over the samples, log2,
#'   averaged over the draws: the x axis of an MA plot) and `resid_sd` (its
#'   residual standard deviation, log2, averaged over the draws; for a
#'   model of one two-level factor, the pooled within-group standard
#'   deviation).
#'
#' @details
#' With \eqn{x_{sf}} the abundance of feature \eqn{f} in sample \eqn{s},
#' each of the `n_draws` draws:
#'
#' 1. Draws each sample's composition from its Dirichlet posterior,
#'    \eqn{G_{sf} \sim \mathrm{Gamma}(x_{sf} +} `pseudocount`\eqn{, 1)}, and
#'    takes \eqn{\log_2} of the proportions \eqn{G_{sf} / \sum_g G_{sg}}.
#' 2. Draws each sample's log2 scale (total load) from the scale model, and
#'    adds it: \eqn{Y_{sf} = \log_2 \mathrm{composition}_{sf} + \log_2
#'    \mathrm{scale}_s}, an estimate of the feature's log2 absolute
#'    abundance.
#' 3. Fits the linear model \eqn{Y_{\cdot f} = X\beta_f + e} to every
#'    feature by least squares (or, with a random intercept, by REML), with
#'    the `se_type` standard errors and t statistics.
#'
#' `estimate` and `se` average the draws' coefficients and standard errors.
#' The p-values follow ALDEx3's `aldex.pvals()`, which allows for draws in
#' which a coefficient changes sign (Nixon et al. 2025): each draw's
#' one-sided p-values for a negative and for a positive coefficient,
#' doubled and capped at 1, are averaged over the draws separately, and
#' `p_value` is the smaller average. For `p_adjusted`, each draw's doubled
#' one-sided p-values are first adjusted across the features (with
#' `p_adjust`, separately for each design column), then averaged and the
#' smaller taken. So `p_adjusted` is an average of adjusted p-values, as in
#' ALDEx2 and ALDEx3, where [dirmult_ttest()] adjusts its averaged
#' p-values, as scikit-bio does. A feature is rejected when `p_adjusted <=
#' alpha`; there is no interval criterion (unlike [dirmult_ttest()]), so
#' decisions settle as `n_draws` grows.
#'
#' **Scale models.** Sequencing gives each sample's composition, not its
#' total load (its scale), and a log fold change of absolute abundance is
#' a fold change of composition plus a fold change of load. Every
#' compositional method assumes something about the second term; a scale
#' model says what, and how sure that is (Nixon et al. 2025):
#'
#' | Information available | Scale model |
#' | :-------------------- | :---------- |
#' | The clr's centre (the geometric mean of the features) is plausible: most features don't change | `scale_model = 'clr'` (default), `scale_sd = 0.5` |
#' | Total load plausibly doesn't differ systematically between groups | `scale_model = 'tss'`, `scale_sd = 0.5` |
#' | Measured or predicted total load per sample (flow cytometry, qPCR, spike-ins, a load predictor) | `load`, with `load_sd` its error (log2) |
#' | Certainty in either assumption (a normalization, as in most methods) | `scale_sd = 0` |
#'
#' With `'clr'`, each sample's log2 scale is minus its mean log2
#' proportion (the clr), plus \eqn{x_s^\top\Lambda}; with `'tss'`, it is
#' \eqn{x_s^\top\Lambda} alone; in both, \eqn{\Lambda \sim N(0,}
#' `scale_sd`\eqn{^2 I)} is drawn afresh in each draw, one value per design
#' column, so that samples with the same design share a shift and every
#' coefficient gets its own uncertainty, added to its estimate's spread over
#' the draws. With `load`, each sample's log2 scale is \eqn{N(\log_2}
#' `load`\eqn{_s,} `load_sd`\eqn{_s^2)}, independently per sample. A larger
#' `scale_sd` (or `load_sd`) widens the averaged p-values, so only changes
#' clear of plausible changes in load are rejected; it controls false
#' discoveries when the assumed scale is wrong, at a cost in power (Nixon et
#' al. 2025; Dos Santos et al. 2026). Check how the results depend on it:
#' [da_sensitivity()] with `scale_sds = TRUE` reruns the model at
#' `scale_sd = 0`, `0.25`, `0.5` and `1`. With a scale model, use many
#' draws: the p-values average over the rare large shifts.
#'
#' **Units.** Every estimate is in log2 units of absolute abundance, as the
#' scale model defines it (the `log_base` attribute). With `'clr'` and
#' `scale_sd = 0` the coefficients are those of a regression of each
#' feature's clr values (log2); the intercept is then a feature's mean clr
#' value in the reference samples. With `load`, they are log2 changes of
#' absolute abundance; the unit of `load` shifts only the intercept.
#'
#' **Relation to [dirmult_ttest()].** For the same `seed` the two functions
#' draw the same Gamma variables for a sample (one random stream per draw
#' and sample of the table), so with `scale_model = 'clr'`, `scale_sd = 0`
#' and `~ group` for two groups, `estimate` for the group's column is
#' `dirmult_ttest()`'s `log2_fc` for `contrast = c(<second group>, <first
#' group>)` (to rounding), and `resid_sd` and `mean` are its `pooled_sd`
#' and `mean`. The tests differ: an OLS t-test with HC3 standard errors here
#' (`se_type = 'classical'` gives the pooled-variance t-test), Welch's
#' t-test there; and their scale models draw different shifts (one per
#' design column here, one for the contrast there).
#'
#' **Repeated measures.** With a random intercept, `~ status + (1 |
#' subject)`, each draw's model is a linear mixed model, fitted by REML,
#' with Satterthwaite's degrees of freedom: ALDEx3's `method = "lme4"`,
#' whose fits use lme4 and lmerTest (McGovern & Silverman 2026), with
#' ecodive's REML fit instead (as [linda()]'s): one random intercept, of a
#' column of `data`, with at least two values and fewer values than
#' samples; random slopes and several random terms are not supported.
#' ecodive finds the exact REML optimum, where lme4's optimizer stops close
#' to it, so results agree with ALDEx3's to about `1e-5`. A mixed model
#' takes longer: about 0.1 ms of CPU time per feature and draw for 100
#' samples in 25 groups, where a linear model takes a few microseconds.
#'
#' **Matching ALDEx3.** Given the same Dirichlet and scale draws, the
#' results are those of ALDEx3 1.4.0's `aldex()` (`method = "lm"`, `test =
#' "t.HC3"`, `"t.HC0"` or `"t"`) and `summary()` to about `1e-12`: the
#' tests feed both the same draws. The draws themselves come from ecodive's
#' own samplers (on PCG, not R's `rgamma()` and `rnorm()`), so for the same
#' data the results differ from ALDEx3's by Monte Carlo error only.
#' Differences from ALDEx3: least squares by a QR decomposition, not the
#' normal equations (more accurate for ill-conditioned designs); the
#' defaults `p_adjust = 'holm'` (ALDEx3: `'BH'`) and `n_draws = 1000`
#' (2,000); ALDEx3 needs a scale model to be named, here `'clr'` with
#' `scale_sd = 0.5` is the default; the intercept's rows are kept (ALDEx3's
#' `summary()` drops them); a rank-deficient design is an error (as in
#' ALDEx3, where `solve()` fails); with HC3 standard errors, a sample of
#' leverage 1 (the only sample of a group) makes the standard errors and
#' p-values `NaN`, with a warning, where ALDEx3 gives meaningless numbers;
#' ALDEx3's `coefficient.sm()` (a prior on each coefficient's shift),
#' custom scale functions, correlated load errors (`s.cov`), one-sided
#' tests, the `nlme` and `blmm` mixed-model engines and the effect-size
#' diagnostics of `aldex.effect()` are not offered; and there is no
#' streaming, as no draws are kept.
#'
#' **Computation.** The draws run in parallel C, one at a time, with memory
#' for one draw (8 bytes per value of the table); every sample of every
#' draw has its own random stream, so results are identical for any number
#' of `cpus`. The p-values' adjustment calls [stats::p.adjust()] for each
#' draw and covariate.
#'
#' @section Same names in other packages:
#'   ALDEx3's `aldex()` is the original method (see **Matching ALDEx3**).
#'   ALDEx2's `aldex.glm()` (Bioconductor) also fits a linear model to each
#'   clr draw, without a scale model by default: with `scale_sd = 0` and
#'   `se_type = 'classical'`, its estimates and standard errors agree with
#'   this function's to Monte Carlo error (its p-values are averages of the
#'   draws' p-values, not ALDEx3's aggregation of their two tails).
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' McGovern, K. C., & Silverman, J. D. (2026). Scale reliant mixed effects
#' models enhance microbiome data analysis. *Microbiome*, 14, 117.
#' \doi{10.1186/s40168-026-02377-x}
#'
#' Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating
#' scale uncertainty in microbiome and gene expression analysis as an
#' extension of normalization. *Genome Biology*, 26, 139.
#' \doi{10.1186/s13059-025-03609-3}
#'
#' Dos Santos, S. J., Murariu, A. C., Silverman, J. D., & Gloor, G. B.
#' (2026). Not every gene is special: modelling scale controls the false
#' discovery rate when analysing high-throughput sequencing data. *PLOS
#' Computational Biology*, 22, e1014728.
#' \doi{10.1371/journal.pcbi.1014728}
#'
#' Nishijima, S., Stankevic, E., Aasmets, O., et al. (2025). Fecal microbial
#' load is a major determinant of gut microbiome variation and a confounder
#' for disease associations. *Cell*, 188, 222-236.
#' \doi{10.1016/j.cell.2024.10.022}
#'
#' Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
#' standard errors in the linear regression model. *The American
#' Statistician*, 54, 217-224.
#'
#' @seealso [dirmult_ttest()], [da_test()], [da_sensitivity()], [linda()],
#'   [qmp_scale()] for loads.
#'
#' @export
#' @examples
#'     # Two groups and a numeric covariate.
#'     set.seed(1)
#'     x <- matrix(rpois(16 * 20, 50), 16, 20, dimnames = list(NULL, paste0('f', 1:20)))
#'     x[9:16, 1:3] <- x[9:16, 1:3] * 4
#'     meta <- data.frame(trt = rep(c('control', 'treated'), each = 8), age = 21:36)
#'
#'     res <- dirmult_lm(x, ~ trt + age, meta, n_draws = 200)
#'     res
#'     plot(res)
#'
#'     # Measured loads, with a log2 error of 0.25.
#'     meta$load <- c(rep(1e9, 8), rep(2e9, 8))
#'     dirmult_lm(x, ~ trt, meta, load = 'load', load_sd = 0.25, n_draws = 200)
#'
#'     # Which decisions hold over a range of scale uncertainties.
#'     sens <- da_sensitivity(x, 'dirmult_lm', 'trt', meta, scale_sds = TRUE, n_draws = 100)
#'     plot(sens, 'traces')
#'
dirmult_lm <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    pseudocount = 0.5, scale_model = c('clr', 'tss'), scale_sd = 0.5,
    load = NULL, load_sd = NULL, n_draws = 1000L,
    se_type = c('hc3', 'hc0', 'classical'), alpha = 0.05, p_adjust = 'holm',
    seed = 0L, cpus = n_cpus() ) {

  dirmult_lm_fit(environment(), given = c(
    scale_model = !missing(scale_model), scale_sd = !missing(scale_sd),
    se_type = !missing(se_type) ))
}



# dirmult_lm() for the tests: the same arguments, plus `gamma_draws` (an
# n x D x n_draws array of Gamma draws to use instead of sampling, samples
# and features as in `counts`' default orientation) and `scale_draws` (an
# n x n_draws matrix of log2 scales to add to the log2 compositions instead
# of a scale model, as aldex(scale = <matrix>)). Returns list(result,
# draws): `draws` holds each draw's estimate, se, df, p_lower and p_upper
# (P x D x n_draws arrays, ALDEx3's p.lower and p.upper), scale (n x n_draws,
# log2: the scale added to each sample's log2 composition) and shift (the
# scale model's normal draws times scale_sd, P x n_draws; or, with `load`,
# times load_sd, n x n_draws).
dirmult_lm_details <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    pseudocount = 0.5, scale_model = c('clr', 'tss'), scale_sd = 0.5,
    load = NULL, load_sd = NULL, n_draws = 1000L,
    se_type = c('hc3', 'hc0', 'classical'), alpha = 0.05, p_adjust = 'holm',
    seed = 0L, cpus = n_cpus(), gamma_draws = NULL, scale_draws = NULL ) {

  dirmult_lm_fit(
    environment(), details = TRUE, gamma_draws = gamma_draws,
    scale_draws = scale_draws, given = c(
      scale_model = !missing(scale_model), scale_sd = !missing(scale_sd),
      se_type = !missing(se_type) ))
}



dirmult_lm_fit <- function (
    env, given, details = FALSE, gamma_draws = NULL, scale_draws = NULL ) {

  info <- dirmult_input(env, 'a Dirichlet-multinomial linear model')
  validate_seed(env)
  seed <- env$seed
  n    <- info$n

  pc <- env$pseudocount
  if (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc <= 0)
    stop('`pseudocount` must be a single positive number.', call. = FALSE)
  pc <- as.double(pc)

  n_draws <- env$n_draws
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1 || n_draws %% 1 != 0 || n_draws > .Machine$integer.max)
    stop('`n_draws` must be a single positive integer.', call. = FALSE)
  n_draws <- as.integer(n_draws)

  alpha <- env$alpha
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !(alpha > 0 && alpha < 1))
    stop('`alpha` must be a single number strictly between 0 and 1.', call. = FALSE)
  alpha <- as.double(alpha)

  p_adjust <- pairwise_method(TRUE, env$p_adjust, TRUE)

  se_type <- env$se_type
  if (!is.character(se_type) || length(se_type) < 1L || anyNA(se_type))
    stop("`se_type` must be 'hc3', 'hc0', or 'classical'.", call. = FALSE)
  se_type <- tryCatch(
    match.arg(tolower(se_type[[1L]]), c('hc3', 'hc0', 'classical')),
    error = function (e) stop("`se_type` must be 'hc3', 'hc0', or 'classical'.", call. = FALSE) )


  #--------------------------------------------------------
  # The design, with an optional random intercept.
  #--------------------------------------------------------
  rnd <- linda_random(env$formula)
  des <- withCallingHandlers(
    da_design(rnd$fixed, env$data, env$groups, info$samples, n, 'dirmult_lm()'),
    warning = function (w)
      if (startsWith(conditionMessage(w), 'The design is rank deficient'))
        invokeRestart('muffleWarning') )
  X <- des$X
  P <- ncol(X)
  dX <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'
  if (qr(dX, tol = 1e-7)$rank < P)
    stop(
      'The design is rank deficient: some of its columns are combinations of ',
      'others, so their coefficients are not unique. Remove the redundant ',
      'terms from `formula`.', call. = FALSE )

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
    if (J >= n)
      stop(
        'The random intercept\'s `', v, '` must have fewer values than there ',
        'are samples (', J, ' values, ', n, ' samples).', call. = FALSE )
    if (given[['se_type']] && se_type != 'classical')
      stop(
        '`se_type = "', se_type, '"`: with a random intercept, the standard ',
        'errors are the mixed model\'s (REML), not heteroskedasticity-',
        'consistent ones.', call. = FALSE )
    se_type <- 'classical'
  }


  #--------------------------------------------------------
  # The scale model.
  #--------------------------------------------------------
  load <- env$load
  mu <- lsd <- NULL
  scale_sd <- env$scale_sd

  if (!is.null(load)) {
    if (given[['scale_model']] || given[['scale_sd']])
      stop(
        'With `load`, the scale model is the measured load and its error ',
        '(`load_sd`): give `scale_model` and `scale_sd` only without ',
        '`load`.', call. = FALSE )
    if (is.null(env$load_sd))
      stop(
        '`load_sd` is required with `load`: the standard deviation of each ',
        'log2 load\'s error, one value or one per sample.', call. = FALSE )
    load <- qmp_load(load, env$data, info$samples, n)
    mu   <- log2(load)
    lsd  <- env$load_sd
    if (is.character(lsd) || length(lsd) != 1L) {
      idx <- NULL
      dat <- env$data
      if (!is.null(dat)) {
        dat <- as.data.frame(dat, stringsAsFactors = FALSE)
        dn  <- if (.row_names_info(dat) > 0) rownames(dat) else NULL
        idx <- match_rows(info$samples, n, dn, nrow(dat), ref = '`counts`')
      }
      lsd <- per_sample(lsd, 'load_sd', dat, idx, info$samples, n, ref = '`counts`')
    }
    if (!is.numeric(lsd) || is.factor(lsd) || length(lsd) < 1L || !all(is.finite(lsd)) || any(lsd < 0))
      stop('`load_sd` must be non-negative numbers: one, or one per sample.', call. = FALSE)
    lsd <- rep_len(as.double(lsd), n)
    mode <- 2L
    scale_model <- 'load'
    scale_sd    <- NA_real_
  }
  else {
    if (!is.null(env$load_sd))
      stop('`load_sd` is only for `load`, which is not given.', call. = FALSE)
    sm <- env$scale_model
    if (!is.character(sm) || length(sm) < 1L || anyNA(sm))
      stop("`scale_model` must be 'clr' or 'tss'.", call. = FALSE)
    scale_model <- tryCatch(
      match.arg(tolower(sm[[1L]]), c('clr', 'tss')),
      error = function (e) stop("`scale_model` must be 'clr' or 'tss'.", call. = FALSE) )
    if (!is.numeric(scale_sd) || length(scale_sd) != 1L || !is.finite(scale_sd) || scale_sd < 0)
      stop('`scale_sd` must be a single non-negative number.', call. = FALSE)
    scale_sd <- as.double(scale_sd)
    mode <- if (scale_model == 'clr') 0L else 1L
  }


  #--------------------------------------------------------
  # Draws given by the tests.
  #--------------------------------------------------------
  if (!is.null(gamma_draws)) {
    stopifnot(
      is.numeric(gamma_draws),
      identical(as.integer(dim(gamma_draws)), c(n, info$D, n_draws)) )
    storage.mode(gamma_draws) <- 'double'
  }
  if (!is.null(scale_draws)) {
    stopifnot(
      is.numeric(scale_draws),
      identical(as.integer(dim(scale_draws)), c(n, n_draws)) )
    storage.mode(scale_draws) <- 'double'
    mode <- 3L
  }

  hc <- match(se_type, c('classical', 'hc0', 'hc3')) - 1L
  if (hc == 2L) {
    h <- stats::hat(dX, intercept = FALSE)
    if (any(1 - h <= 1e-10))
      warning(
        'A sample has leverage 1 (such as the only sample of a group), so ',
        'the HC3 standard errors are undefined: they and the p-values are ',
        'NaN. Use `se_type = "classical"`.', call. = FALSE )
  }

  padj <- if (p_adjust == 'none') NULL else {
    method <- p_adjust
    function (p) stats::p.adjust(p, method)
  }

  res <- .Call(
    C_dirmult_lm, env$counts, env$margin, pc, n_draws, dX, hc, mode,
    if (mode <= 1L) scale_sd else 0, mu, lsd, codes, padj, seed,
    gamma_draws, scale_draws, details, env$cpus )


  #--------------------------------------------------------
  # The table: a row per feature and design column.
  #--------------------------------------------------------
  D <- info$D
  features   <- if (is.null(info$features)) as.character(seq_len(D)) else info$features
  covariates <- colnames(X)
  vec <- function (m) { v <- as.vector(m); v[is.na(v)] <- NaN; v }

  p_adjusted <- vec(res$p_adjusted)
  out <- data.frame(
    feature    = rep(features, each = P),
    covariate  = rep(covariates, times = D),
    estimate   = vec(res$estimate),
    se         = vec(res$se),
    dof        = vec(res$df),
    p_value    = vec(res$p_value),
    p_adjusted = p_adjusted,
    reject     = !is.na(p_adjusted) & p_adjusted <= alpha,
    stringsAsFactors = FALSE )

  resid_sd <- res$resid_sd
  resid_sd[is.na(resid_sd)] <- NaN

  attr(out, 'method')      <- 'dirmult_lm'
  attr(out, 'formula')     <- if (mixed) linda_formula(des$formula, rnd$random) else des$formula
  attr(out, 'groups')      <- des$groups
  attr(out, 'design')      <- X
  if (mixed) attr(out, 'random') <- rnd$random
  attr(out, 'scale_model') <- scale_model
  attr(out, 'scale_sd')    <- scale_sd
  if (!is.null(lsd)) attr(out, 'load_sd') <- lsd
  attr(out, 'se_type')     <- se_type
  attr(out, 'n_draws')     <- n_draws
  attr(out, 'pseudocount') <- pc
  attr(out, 'alpha')       <- alpha
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'seed')        <- seed
  attr(out, 'log_base')    <- 2
  attr(out, 'mean')        <- stats::setNames(res$mean, features)
  attr(out, 'resid_sd')    <- stats::setNames(resid_sd, features)
  class(out) <- c('ecodive_dirmult_lm', 'data.frame')

  if (!details) return (out)

  dr <- res$details
  for (k in c('estimate', 'se', 'df', 'p_lower', 'p_upper'))
    dimnames(dr[[k]]) <- list(covariates, features, NULL)
  list(result = out, draws = dr)
}


#' Plot a Dirichlet-multinomial Linear Model
#'
#' Views of a [dirmult_lm()] result for one covariate, each feature a point
#' (or a row), with the features declared differentially abundant (`reject
#' = TRUE`) highlighted and the most significant of them named:
#'
#' * `"volcano"` (the default): the estimate against the adjusted p-value.
#' * `"effect"`: the estimate against the residual standard deviation (for
#'   a two-level comparison).
#' * `"ma"`: the estimate against the mean clr abundance.
#' * `"forest"`: the estimates and approximate confidence intervals of the
#'   `n_features` most significant features.
#' * `"rank"`: every feature's estimate, sorted (see **Rank plots**).
#'
#' @param x  The result of [dirmult_lm()].
#'
#' @param covariate  Which covariate's estimates to plot, as named in the
#'        result's `covariate` column (such as `"statusmoderate"` or
#'        `"age"`). Defaults to the first column of the term of interest
#'        (`groups`), or else to the first covariate that isn't the
#'        intercept.
#'
#' @param conf_level  For `"forest"`: the confidence level of the
#'        intervals. Default `0.95`.
#'
#' @inheritParams plot.ecodive_linda
#'
#' @inheritSection plot.ecodive_ancombc Rank plots
#'
#' @details
#' Estimates are in log2 units (the result's `log_base`). The subtitle names
#' the covariate, counts the features rejected, and gives the scale model.
#'
#' * **Volcano plot.** The y axis is \eqn{-\log_{10}} of `p_adjusted` (an
#'   adjusted p-value of 0 is drawn at the smallest positive double), with
#'   a dashed line at `alpha`.
#' * **Effect plot.** The x axis is the feature's residual standard
#'   deviation (the `resid_sd` attribute, log2), which for a model of one
#'   two-level factor is the pooled within-group standard deviation; the
#'   dashed lines are \eqn{y = \pm x}. Only for a covariate whose design
#'   column is 0 or 1.
#' * **MA plot.** The x axis is the feature's mean clr value (the `mean`
#'   attribute, log2).
#' * **Forest plot.** Bars over \eqn{\mathrm{estimate} \pm q\,\mathrm{se}},
#'   with \eqn{q} the \eqn{(1 + }`conf_level`\eqn{)/2} quantile of the t
#'   distribution with `dof` degrees of freedom: approximate, as the
#'   estimate and standard error are averages over the draws, and the
#'   p-values come from the draws instead. They aren't adjusted for testing
#'   many features.
#'
#' Features whose values are `NaN` are left out.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     x <- matrix(rpois(16 * 20, 50), 16, 20)
#'     x[9:16, 1:3] <- x[9:16, 1:3] * 4
#'     res <- dirmult_lm(x, groups = rep(c('a', 'b'), each = 8), n_draws = 200)
#'     plot(res)
#'     plot(res, type = 'effect')
#'     plot(res, type = 'forest', n_features = 10)
plot.ecodive_dirmult_lm <- function (
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
    est = r$estimate, p = r$p_adjusted, reject = r$reject, feature = r$feature,
    sd = attr(x, 'resid_sd', exact = TRUE), mean = attr(x, 'mean', exact = TRUE),
    lo = r$estimate - q * r$se, hi = r$estimate + q * r$se,
    what     = covariate,
    est_lab  = lfc_label(x),
    sd_lab   = 'Residual standard deviation (log2)',
    mean_lab = 'Mean abundance (clr, log2)',
    interval = sprintf('approximate %s confidence intervals', plot_pct(conf_level)),
    num = numerator, den = denominator, extra = dirmult_lm_scale_label(x) )

  invisible(x)
}


# The scale model, for plots' subtitles: "scale model: clr, scale_sd = 0.5",
# or "scale model: measured load"; NULL if the result has none.
dirmult_lm_scale_label <- function (x) {
  sm <- attr(x, 'scale_model', exact = TRUE)
  if (!is.character(sm) || length(sm) != 1L) return (NULL)
  if (sm == 'load') return ('scale model: measured load')
  sd <- attr(x, 'scale_sd', exact = TRUE)
  sprintf('scale model: %s, scale_sd = %s', sm, format(sd))
}
