# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# maaslin3() is a port of MaAsLin 3 (maaslin3 1.4.0, Bioconductor; MIT,
# Copyright (c) 2025 Harvard School of Public Health; the notice is in
# src/maaslin3.c and inst/COPYRIGHTS): R/maaslin3.R maaslin_fit(),
# R/fit.R add_qvals(), create_combined_pval(),
# flag_abundance_turned_prevalence(), run_median_comparison(),
# restructured for ecodive (see plans/maaslin3.md). The fits are in
# src/maaslin3.c.

#' Multivariable Association with Linear Models (MaAsLin 3)
#'
#' Tests each feature for association with the covariates of a model, twice
#' (Nickols et al. 2026): its **abundance** where it is present (a linear
#' model of its log2 relative abundance, zeros left out) and its
#' **prevalence** (a logistic model of whether it is present). Abundance
#' coefficients are tested against the median coefficient over all features
#' (the *median comparison*), which turns a change in relative abundance
#' into a statement about absolute abundance when the typical feature
#' doesn't change. A joint p-value combines the two models. A port of
#' maaslin3 1.4.0.
#'
#' @inheritParams da_wilcoxon
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), counts or proportions, or any of the
#'        objects listed under **Input Types**. A value of 0 means the
#'        feature is absent.
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
#' @param groups  The term of interest (a fixed effect): with `formula`,
#'        the name of one of its terms (the one [plot()] shows by default);
#'        alone, the name of a column of `data` or a vector with one value
#'        per sample, and the model is `~ groups`.
#'
#' @param norm  How abundances are normalized for the abundance model:
#'        `'percent'` (the default; each sample's values divided by its
#'        total, maaslin3's `'TSS'`), `'rclr'` (log ratios to the geometric
#'        mean of the sample's non-zero values, natural log; maaslin3's
#'        `'CLR'`), or `'none'`. Zeros are left out of the abundance model
#'        in every case.
#'
#' @param transform  `'log'` (the default): \eqn{\log_2} of the normalized
#'        values; or `'none'`. `'log'` can't be used with `norm = 'rclr'`,
#'        whose values are logs already.
#'
#' @param standardize  `TRUE` (the default, as maaslin3) to centre and scale
#'        every numeric variable of `formula` to unit standard deviation
#'        before fitting, so that their coefficients are per standard
#'        deviation; `FALSE` for coefficients per unit.
#'
#' @param augment  `TRUE` (the default, as maaslin3) to fit the prevalence
#'        model to augmented data, which keeps its estimates finite when a
#'        covariate separates the samples where a feature is present from
#'        those where it is absent (see **The prevalence model**).
#'
#' @param median_comparison_abundance,median_comparison_prevalence  Whether
#'        to test each model's coefficients against their median over the
#'        features (see **The median comparison**), rather than against 0.
#'        Defaults `TRUE` for abundance and `FALSE` for prevalence, as
#'        maaslin3.
#'
#' @param median_comparison_abundance_threshold,median_comparison_prevalence_threshold
#'        A coefficient within this distance of the median gets a p-value of
#'        1 in the median comparison (divided by the standard deviation of
#'        the design column when it isn't 0/1). Default `0`.
#'
#' @param subtract_median  `TRUE` to report the coefficients of a model
#'        with a median comparison minus that median; `FALSE` (the default)
#'        for the coefficients as fitted.
#'
#' @param warn_prevalence  `TRUE` (the default, as maaslin3) to flag
#'        prevalence associations that a stronger abundance association in
#'        the same direction may have caused (see **Flags**). Needs
#'        `transform = 'log'` and a `norm` other than `'rclr'`.
#'
#' @param small_random_effects  With a random intercept, `TRUE` to fit the
#'        prevalence model with its grouping as a fixed effect instead (as
#'        maaslin3's option of that name), which is better when the groups
#'        have few samples each (fewer than about 4). Default `FALSE`.
#'
#' @param n_sims  The number of Monte Carlo draws for the median
#'        comparison's covariances. Default `10000`, as maaslin3.
#'
#' @param alpha  The significance level for `reject`, `reject_abundance`
#'        and `reject_prevalence` (maaslin3's `max_significance`). Default
#'        `0.05` (maaslin3: `0.1`).
#'
#' @param p_adjust  How to adjust the p-values for testing many features
#'        and covariates: any method of [stats::p.adjust()] (`'holm'`, the
#'        default, as ecodive's other DA methods; `'BH'`, as maaslin3;
#'        `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BY'`, `'fdr'`,
#'        `'none'`). The abundance and prevalence p-values are adjusted
#'        together, as one family; the joint p-values separately.
#'
#' @param seed  An integer seed for the median comparison's Monte Carlo
#'        draws. Default `0`.
#'
#' @return A data frame of class `c('ecodive_maaslin3', 'data.frame')` with
#'   a row per feature and covariate (every design column but the
#'   intercept), in the order of `counts`:
#'
#'   * `feature`, `covariate`: the feature's name (or its index as text),
#'     and the design column, as [stats::model.matrix()] names it
#'     (`'statusmoderate'`, `'age'`).
#'   * `coef_abundance`, `se_abundance`: the abundance model's coefficient
#'     (log2 fold change in relative abundance where present, with the
#'     default `norm` and `transform`) and its standard error.
#'   * `dof_abundance`: its degrees of freedom: the residual's, or with a
#'     random intercept Satterthwaite's.
#'   * `null_abundance`: the value it was tested against: the median, with
#'     a median comparison, else 0.
#'   * `p_abundance`, `q_abundance`: its p-value and adjusted p-value.
#'   * `coef_prevalence`, `se_prevalence`, `null_prevalence`,
#'     `p_prevalence`, `q_prevalence`: the same for the prevalence model
#'     (log odds ratios; Wald z-tests).
#'   * `n_nonzero`: the number of samples where the feature is present.
#'   * `flag_abundance`, `flag_prevalence`: maaslin3's error or warning for
#'     the fit, or `NA` (see **Flags**).
#'   * `p_value`, `p_adjusted`: the joint p-value of the two models and its
#'     adjustment.
#'   * `reject`: `TRUE` when the joint adjusted p-value is at most `alpha`,
#'     or either model's is (and that model has no flag): maaslin3's
#'     significant associations.
#'   * `reject_abundance`, `reject_prevalence`: each model's own decision
#'     (adjusted p-value at most `alpha`, no flag).
#'
#'   Values that can't be estimated are `NaN`. Attributes: `method`
#'   (`'maaslin3'`), `formula`, `groups` (the term of interest, or `NULL`),
#'   `design` (the model matrix, standardized), `random` (with a random
#'   intercept: its grouping variable), `norm`, `transform`, `alpha`,
#'   `p_adjust`, `log_base` (of the abundance coefficients: `2` with
#'   `transform = 'log'`, `exp(1)` with `norm = 'rclr'`, else `NA`), and
#'   per feature (named by feature) `mean` (its mean transformed abundance
#'   where present), `prevalence` (the proportion of samples where it is
#'   present), `resid_sd` (the abundance model's residual standard
#'   deviation) and, with a random intercept, `random_sd` (the abundance
#'   model's random intercepts' standard deviation). [plot()] draws volcano,
#'   MA, forest and rank plots of either model
#'   ([plot.ecodive_maaslin3()]).
#'
#' @details
#' **The abundance model.** For each feature, its values are normalized
#' (`norm`), its zeros set aside, and the rest transformed (`transform`);
#' then a linear model of `formula` is fitted by least squares to the
#' samples where the feature is present, and each coefficient tested with
#' a t-test. A feature whose non-zero values are all the same (or that has
#' none) is flagged "All linear values are the same"; one present in no
#' sample at a factor's reference level, "No data points have the baseline
#' factor level"; one present at a single level of a factor, "contrasts can
#' be applied only to factors with 2 or more levels". Their abundance
#' results are `NaN`. A design column that is constant among a feature's
#' samples can't be estimated, and gets `NaN` with the flag "Fitting error
#' (NA p-value returned from fitting procedure)".
#'
#' **The prevalence model.** Each feature's presence (1) or absence (0) is
#' fitted to `formula` by logistic regression over all samples
#' (iteratively reweighted least squares, as [stats::glm()]), and tested
#' with Wald z-tests. With `augment = TRUE` the data are augmented first
#' (Nickols et al. 2026): every sample twice more, once present and once
#' absent, each copy with weight \eqn{p / (2n)} (\eqn{p} the number of
#' design columns besides the intercept, \eqn{n} the samples), which keeps
#' the estimates finite under separation. A feature present in every sample
#' or in none is flagged "All logistic values are the same". glm()'s
#' warnings become flags: "glm.fit: algorithm did not converge" (25
#' iterations), "glm.fit: fitted probabilities numerically 0 or 1
#' occurred".
#'
#' **The median comparison.** If most features don't change in absolute
#' abundance, the coefficients of relative abundance cluster around the
#' shift that the composition puts on all of them, so each coefficient
#' \eqn{\hat\beta_j} is tested against their median \eqn{m} over the
#' features whose own p-value is below 0.95: \eqn{t_j = (\hat\beta_j -
#' m)/\sqrt{se_j^2 + v_m - 2 c_j}}, with \eqn{v_m = \frac{\pi}{2}\,
#' \mathrm{var}(\hat\beta)/n} the median's asymptotic variance and \eqn{c_j}
#' its covariance with \eqn{\hat\beta_j}, estimated from `n_sims` Monte
#' Carlo draws of every coefficient from \eqn{N(\hat\beta_j, se_j^2)}. The
#' t-test has the model's degrees of freedom (the normal distribution for
#' prevalence). It is on by default for abundance only. Turn it off
#' (`median_comparison_abundance = FALSE`) for absolute abundances, such as
#' those of [qmp_scale()], or when many features are expected to change in
#' the same direction.
#'
#' **The joint p-value.** For each feature and covariate, \eqn{p =
#' 1 - (1 - \min(p_{abundance}, p_{prevalence}))^2} (the distribution of the
#' smaller of two independent uniform p-values, `pbeta(p, 1, 2)`), or the
#' p-value of the model that has one without a flag; adjusted with
#' `p_adjust` over all features and covariates.
#'
#' **Flags.** Besides the fitting problems above, as maaslin3:
#'
#' * A prevalence coefficient larger than 15 in absolute value, or a
#'   p-value below \eqn{10^{-10}}, of a feature present in under 5\% of
#'   the samples (and fewer than 50), or absent from under 5\% (and fewer
#'   than 50), is flagged as a likely misfit.
#' * With `warn_prevalence`: a prevalence association in the same direction
#'   as a stronger abundance association (of log2 proportions, without a
#'   median comparison) whose adjusted p-value is below `alpha` is flagged
#'   "Prevalence association possibly induced by stronger abundance
#'   association": a feature whose abundance drops can fall below detection
#'   in more samples. The flag leaves the prevalence p-value out of the joint
#'   p-value and out of `reject_prevalence`.
#' * With a random intercept whose groups have fewer than 4 samples on
#'   average, the prevalence rows are flagged "<4 average observations per
#'   random effect group often inflates coefficients and deflates p-values:
#'   consider setting small_random_effects=TRUE and see tutorial". This flag
#'   alone doesn't stop a rejection.
#'
#' A flag other than these two leaves the model's p-value out of the
#' adjustments.
#'
#' **Repeated measures.** With `(1 | id)` in `formula`, both models get a
#' random intercept for each value of `id` (a subject, say). The abundance
#' model is a linear mixed model fitted by restricted maximum likelihood
#' (REML), as lme4's `lmer()`, with Satterthwaite's degrees of freedom, as
#' lmerTest; ecodive finds the REML optimum exactly (see [linda()]). The
#' prevalence model is a logistic mixed model fitted by the Laplace
#' approximation, as lme4's `glmer()`, minimised over the random
#' intercept's standard deviation and the coefficients by Newton's method;
#' its standard errors are those of glmer (from the criterion's Hessian, or
#' for a singular fit, whose standard deviation is below \eqn{10^{-4}},
#' from the penalized least squares). With few samples per group, consider
#' `small_random_effects = TRUE`.
#'
#' **Differences from maaslin3.** The results are those of maaslin3 1.4.0
#' with its filters off (to about \eqn{10^{-10}}; the median comparison's
#' p-values to about \eqn{10^{-8}}, the precision of the numerical t
#' integral maaslin3 uses), when its median comparison's draws are given
#' (see **Computation**), except that:
#'
#' * Every feature is reported. maaslin3's filters (`min_abundance`,
#'   `min_prevalence`, `max_prevalence`, `min_variance`) aren't ported:
#'   filter `counts` first. Even with them off, maaslin3 leaves out
#'   features with fewer than two non-zero values or all non-zero values
#'   equal (its `min_variance = 0`), and with `normalization = 'CLR'` those
#'   whose CLR values are all at most 0; here they are kept, and enter the
#'   adjustments and the median comparison through the model they have.
#' * The median comparison's draws come from PCG streams (`seed`), not R's
#'   random numbers, so its p-values agree with maaslin3's to Monte Carlo
#'   error (about 1\% relative), and are the same for any `cpus`.
#' * `p_adjust` defaults to `'holm'` and `alpha` to `0.05` (maaslin3:
#'   `'BH'` and `0.1`).
#' * Factor levels: a factor's own order, else sorted in the C locale;
#'   maaslin3 requires `reference` for a variable of more than two levels
#'   that isn't a factor.
#' * Missing values in `data` are an error (maaslin3 leaves those samples
#'   out of each fit).
#' * Of maaslin3's model terms, one random intercept, `(1 | id)`; not
#'   `group()`, `ordered()` or `strata()` terms, random slopes or several
#'   random effects. Not ported either: `feature_specific_covariate`,
#'   `unscaled_abundance`, contrast tests, `evaluate_only`, the plots and
#'   files, `zero_threshold` (zero is absent), and the `PLOG` transform
#'   (with zeros set aside it equals `LOG`).
#' * Mixed models are fitted at their exact optimum; lme4's optimizers stop
#'   a little short (maaslin3 keeps the first optimizer's fit that gives no
#'   warning), so results differ from maaslin3's by about \eqn{10^{-6}} for
#'   abundance and \eqn{10^{-4}} to \eqn{10^{-3}} for prevalence, more
#'   where glmer's first optimizer stops early on a flat criterion. lme4's
#'   convergence warnings, which maaslin3 turns into flags, aren't
#'   reproduced; a prevalence fit that doesn't converge here is flagged
#'   "The mixed model did not converge".
#' * The result is one data frame, a row per feature and covariate with
#'   both models side by side (maaslin3: a table per model, sorted by
#'   q-value).
#'
#' **Computation.** The table is read by feature from its sparse form, and
#' the features are fitted in parallel C, each from its own samples, so
#' results are identical for any number of `cpus`. The median comparison
#' draws each Monte Carlo draw from its own PCG stream (keyed by `seed`,
#' the model, the covariate and the draw), and sums the covariances in
#' blocks of draws in a fixed order.
#'
#' @section Same names in other packages:
#' maaslin3 (Bioconductor) has `maaslin3()`, the original.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Nickols, W. A., Kuntz, T., Shen, J., et al. (2026). MaAsLin 3: refining
#' and extending generalized multivariable linear models for meta-omic
#' association discovery. *Nature Methods*, 23(3), 554-564.
#' \doi{10.1038/s41592-025-02923-9}
#'
#' Mallick, H., Rahnavard, A., McIver, L. J., et al. (2021). Multivariable
#' association discovery in population-scale meta-omics studies. *PLoS
#' Computational Biology*, 17(11), e1009442.
#' \doi{10.1371/journal.pcbi.1009442}
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting
#' linear mixed-effects models using lme4. *Journal of Statistical
#' Software*, 67(1), 1-48. \doi{10.18637/jss.v067.i01}
#'
#' Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017).
#' lmerTest package: tests in linear mixed effects models. *Journal of
#' Statistical Software*, 82(13), 1-26. \doi{10.18637/jss.v082.i13}
#'
#' @seealso [da_test()], [linda()], [da_lm()], [ancombc2()].
#'
#' @export
#' @examples
#'     # 40 samples, 30 features; f1 to f4 are more abundant in group b,
#'     # and f5 is more often absent there.
#'     set.seed(1)
#'     x <- matrix(rpois(40 * 30, 8) * rbinom(40 * 30, 1, 0.7), 40, 30)
#'     colnames(x) <- paste0('f', 1:30)
#'     meta <- data.frame(
#'       grp = rep(c('a', 'b'), each = 20),
#'       age = round(runif(40, 20, 60)) )
#'     x[21:40, 1:4] <- x[21:40, 1:4] * 5
#'     x[21:40, 5]   <- x[21:40, 5] * rbinom(20, 1, 0.2)
#'
#'     res <- maaslin3(x, ~ grp + age, meta, groups = 'grp')
#'     res
#'     plot(res)
#'     plot(res, model = 'prevalence')
#'
#'     # As maaslin3's defaults would give it.
#'     maaslin3(x, ~ grp + age, meta, p_adjust = 'BH', alpha = 0.1)
#'
#'     # Two samples from each of 20 subjects: a random intercept per subject.
#'     meta$subject <- rep(sprintf('s%02d', 1:20), 2)
#'     maaslin3(x, ~ grp + age + (1 | subject), meta, groups = 'grp')
#'
maaslin3 <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    norm = 'percent', transform = 'log', standardize = TRUE, augment = TRUE,
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
    median_comparison_abundance_threshold = 0,
    median_comparison_prevalence_threshold = 0, subtract_median = FALSE,
    warn_prevalence = TRUE, small_random_effects = FALSE, n_sims = 10000L,
    alpha = 0.05, p_adjust = 'holm', seed = 0L, cpus = n_cpus() ) {

  maaslin3_fit(environment())$result
}



# maaslin3() for the tests: the same arguments, plus `sims`, a list with
# elements `abundance` and `prevalence`, each NULL or a list named by
# covariate of matrices of standard normal draws (one row per feature
# entering that covariate's median comparison, in table order; one column
# per draw) to use instead of PCG's. Returns list(result, abundance,
# prevalence): per model, the fits (D x K matrices coef, se, df, p, error)
# before the median comparison, and the median comparison's `median`,
# `sd_median`, `cov` (per covariate).
maaslin3_details <- function (
    counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
    norm = 'percent', transform = 'log', standardize = TRUE, augment = TRUE,
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
    median_comparison_abundance_threshold = 0,
    median_comparison_prevalence_threshold = 0, subtract_median = FALSE,
    warn_prevalence = TRUE, small_random_effects = FALSE, n_sims = 10000L,
    alpha = 0.05, p_adjust = 'holm', seed = 0L, cpus = n_cpus(), sims = NULL ) {

  maaslin3_fit(environment(), sims = sims)
}



# maaslin3's error texts, by the status src/maaslin3.c returns.
M3_TEXT <- c(
  same_abundance  = 'All linear values are the same',
  same_prevalence = 'All logistic values are the same',
  baseline        = 'No data points have the baseline factor level',
  one_level       = 'contrasts can be applied only to factors with 2 or more levels',
  one_group       = 'grouping factors must have > 1 sampled level',
  many_groups     = 'number of levels of each grouping factor must be < number of observations',
  failed          = 'The model could not be fitted',
  no_conv         = 'The mixed model did not converge',
  glm_conv        = 'glm.fit: algorithm did not converge',
  glm_01          = 'glm.fit: fitted probabilities numerically 0 or 1 occurred',
  na_p            = 'Fitting error (NA p-value returned from fitting procedure)',
  median_na       = 'P-value became NA in median comparison, try rerunning without the median comparison',
  rare            = 'A large coefficient (>15 in absolute value) or small p-value (< 10^-10) was obtained from a feature present in <5% of samples. Check this is intended.',
  common          = 'A large coefficient (>15 in absolute value) or small p-value (< 10^-10) was obtained from a feature present in >95% of samples. Check this is intended.',
  induced         = 'Prevalence association possibly induced by stronger abundance association',
  small_groups    = '<4 average observations per random effect group often inflates coefficients and deflates p-values: consider setting small_random_effects=TRUE and see tutorial' )



maaslin3_fit <- function (env, sims = NULL) {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  norm <- env$norm
  env$norm <- 'none'
  info <- da_input(env)
  counts <- env$counts
  n <- info$n
  D <- info$D

  if (any(stored_values(counts) < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)

  norm <- m3_choice(norm, 'norm', c('percent', 'rclr', 'none'))
  transform <- m3_choice(env$transform, 'transform', c('log', 'none'))
  if (norm == 'rclr' && transform == 'log')
    stop(
      "`transform = 'log'` needs positive values, and rclr values are logs ",
      "already (with negative values): use `transform = 'none'`.", call. = FALSE )

  if (isTRUE(env$warn_prevalence) && (norm == 'rclr' || transform != 'log'))
    stop(
      '`warn_prevalence = TRUE` compares with the abundance model of log2 ',
      "proportions, so it needs `transform = 'log'` (and not `norm = 'rclr'`), ",
      'as in maaslin3: set `warn_prevalence = FALSE`.', call. = FALSE )

  flags <- c(
    'standardize', 'augment', 'median_comparison_abundance',
    'median_comparison_prevalence', 'subtract_median', 'warn_prevalence',
    'small_random_effects' )
  for (a in flags)
    if (!isTRUE(env[[a]]) && !isFALSE(env[[a]]))
      stop('`', a, '` must be TRUE or FALSE.', call. = FALSE)
  for (a in c('median_comparison_abundance_threshold', 'median_comparison_prevalence_threshold')) {
    v <- env[[a]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0)
      stop('`', a, '` must be a single non-negative number.', call. = FALSE)
  }

  n_sims <- env$n_sims
  if (!is.numeric(n_sims) || length(n_sims) != 1L || !is.finite(n_sims) ||
      n_sims < 2 || n_sims %% 1 != 0 || n_sims > .Machine$integer.max)
    stop('`n_sims` must be a single integer, at least 2.', call. = FALSE)
  n_sims <- as.integer(n_sims)

  validate_seed(env)
  alpha    <- env$alpha
  p_adjust <- env$p_adjust
  cpus     <- env$cpus


  #--------------------------------------------------------
  # The design, with an optional random intercept.
  #--------------------------------------------------------
  rnd <- linda_random(env$formula)
  des <- da_design(
    rnd$fixed, env$data, env$groups, info$samples, n, 'maaslin3()',
    standardize = env$standardize )
  X  <- des$X
  P  <- ncol(X)
  dX <- X
  attributes(dX) <- list(dim = dim(X))
  storage.mode(dX) <- 'double'

  keep <- setdiff(seq_len(P), des$intercept + 1L)
  if (!length(keep))
    stop('`formula` must have at least one term besides the intercept.', call. = FALSE)
  covariates <- colnames(X)[keep]
  K <- length(keep)

  mixed <- !is.null(rnd$random)
  codes <- NULL
  small <- FALSE
  if (mixed) {
    v <- rnd$random
    if (!v %in% names(des$data))
      stop('Not found in `data`: ', name_list(v), call. = FALSE)
    idv <- des$data[[v]]
    if (anyNA(idv))
      stop('Missing values (NA) in `data` column ', name_list(v), call. = FALSE)
    idf   <- factor(as.character(idv), levels = sort(unique(as.character(idv)), method = 'radix'))
    codes <- as.integer(idf)
    J <- max(codes)
    if (J < 2L)
      stop('The random intercept\'s `', v, '` must have at least two values.', call. = FALSE)
    if (J >= n)
      stop(
        'The random intercept\'s `', v, '` must have fewer values than there ',
        'are samples (', J, ' values, ', n, ' samples).', call. = FALSE )
    small <- env$small_random_effects
  }
  else if (env$small_random_effects) {
    stop('`small_random_effects = TRUE` needs a random intercept, `(1 | id)`, in `formula`.', call. = FALSE)
  }

  # The prevalence model's design: with small_random_effects, the grouping
  # as a fixed effect (treatment contrasts), as maaslin3's.
  Xp <- dX
  if (small) {
    Z  <- stats::model.matrix(~ idf)[, -1L, drop = FALSE]
    attributes(Z) <- list(dim = dim(Z))
    Xp <- cbind(dX, Z)
  }

  # Each factor's levels (1: the reference), for maaslin3's check that a
  # feature is present in some sample at the reference level, and lm()'s
  # that it is present at two levels.
  base <- NULL
  for (v in all.vars(stats::terms(des$formula))) {
    x <- des$data[[v]]
    if (is.factor(x)) base <- cbind(base, as.integer(x))
  }
  if (!is.null(base)) storage.mode(base) <- 'integer'

  ncode  <- match(norm, c('none', 'percent', 'rclr')) - 1L
  is_log <- as.integer(transform == 'log')
  a_aug  <- if (env$augment) (ncol(Xp) - 1) / (2 * n) else 0


  #--------------------------------------------------------
  # The fits.
  #--------------------------------------------------------
  ab <- .Call(
    C_maaslin3_fit, counts, env$margin, ncode, is_log, dX, 0L, codes, base, 0, cpus )
  pr <- .Call(
    C_maaslin3_fit, counts, env$margin, ncode, is_log, Xp, 1L,
    if (mixed && !small) codes, NULL, as.double(a_aug), cpus )

  abt <- m3_table(ab, keep, 'abundance', mixed, covariates)
  prt <- m3_table(pr, keep, 'prevalence', mixed && !small, covariates)
  ab_raw <- abt

  # The median comparisons.
  mca <- mcp <- NULL
  if (env$median_comparison_abundance) {
    mca <- m3_median(
      abt, ab$xsd[, keep, drop = FALSE], env$median_comparison_abundance_threshold,
      env$subtract_median, n_sims, env$seed, 0L, sims$abundance, covariates, cpus )
    abt <- mca$table
  }
  if (env$median_comparison_prevalence) {
    mcp <- m3_median(
      prt, pr$xsd[, keep, drop = FALSE], env$median_comparison_prevalence_threshold,
      env$subtract_median, n_sims, env$seed, 1L, sims$prevalence, covariates, cpus )
    prt <- mcp$table
  }

  # Prevalence fits that look like misfits: a feature in under 5% (or over
  # 95%) of the samples, and fewer than 50 of them, with |coef| > 15 or
  # p < 1e-10.
  nz   <- matrix(ab$n_obs, D, K)
  big  <- (!is.na(prt$coef) & abs(prt$coef) > 15) | (!is.na(prt$p) & prt$p < 1e-10)
  rare <- nz < 50 & nz / n < 0.05 & big
  prt$error[rare & is.na(prt$error)] <- M3_TEXT[['rare']]
  comm <- (n - nz) < 50 & nz / n > 0.95 & big
  prt$error[comm & is.na(prt$error)] <- M3_TEXT[['common']]

  # q-values: abundance and prevalence together.
  q <- m3_adjust(abt, prt, p_adjust)
  abt$q <- q$a
  prt$q <- q$b

  # Prevalence associations that a stronger abundance association may have
  # induced: the abundance model as percent + log, without a median
  # comparison.
  if (env$warn_prevalence) {
    nab <- if (norm == 'percent' && transform == 'log') ab_raw
           else m3_table(
             .Call(C_maaslin3_fit, counts, env$margin, 1L, 1L, dX, 0L, codes, base, 0, cpus),
             keep, 'abundance', mixed, covariates )
    nq  <- m3_adjust(nab, prt, p_adjust)$a
    tst <- !is.na(prt$coef) & !is.na(nab$coef) & !is.na(prt$p) & !is.na(nab$p)
    inner <- is.na(prt$error) & is.na(nab$error) & !is.na(nq) & nq < alpha &
      sign(prt$coef) == sign(nab$coef) & abs(nab$coef) > abs(prt$coef)
    prt$error[tst & !is.na(inner) & inner] <- M3_TEXT[['induced']]
  }

  # The joint p-value: the smaller of the two, as the minimum of two
  # independent uniforms (pbeta(p, 1, 2)), or the one without an error.
  ab_bad <- is.na(abt$p) | !is.na(abt$error)
  pr_bad <- is.na(prt$p) | !is.na(prt$error)
  pj <- stats::pbeta(pmin(abt$p, prt$p), 1, 2)
  pj[ab_bad] <- prt$p[ab_bad]
  pj[pr_bad] <- abt$p[pr_bad]
  pj[ab_bad & pr_bad] <- NA
  qj <- pj
  ok <- !is.na(pj)
  qj[ok] <- stats::p.adjust(pj[ok], p_adjust)

  # Too few samples per random intercept group for the logistic model.
  if (mixed && !small && mean(table(codes)) < 4)
    prt$error[is.na(prt$error)] <- M3_TEXT[['small_groups']]

  rej_ab <- !is.na(abt$q) & abt$q <= alpha & is.na(abt$error)
  rej_pr <- !is.na(prt$q) & prt$q <= alpha &
    (is.na(prt$error) | prt$error == M3_TEXT[['small_groups']])
  rej    <- (!is.na(qj) & qj <= alpha) | rej_ab | rej_pr


  #--------------------------------------------------------
  # The result, a row per feature and covariate.
  #--------------------------------------------------------
  byrow <- function (m) { v <- as.vector(t(m)); if (is.double(v)) v[is.na(v)] <- NaN; v }
  bychr <- function (m) as.vector(t(m))
  features <- info$features

  out <- data.frame(
    feature           = rep(features, each = K),
    covariate         = rep(covariates, times = D),
    coef_abundance    = byrow(abt$coef),
    se_abundance      = byrow(abt$se),
    dof_abundance     = byrow(abt$df),
    null_abundance    = byrow(abt$null),
    p_abundance       = byrow(abt$p),
    q_abundance       = byrow(abt$q),
    coef_prevalence   = byrow(prt$coef),
    se_prevalence     = byrow(prt$se),
    null_prevalence   = byrow(prt$null),
    p_prevalence      = byrow(prt$p),
    q_prevalence      = byrow(prt$q),
    n_nonzero         = rep(as.integer(ab$n_obs), each = K),
    flag_abundance    = bychr(abt$error),
    flag_prevalence   = bychr(prt$error),
    p_value           = byrow(matrix(pj, D, K)),
    p_adjusted        = byrow(matrix(qj, D, K)),
    reject            = as.vector(t(rej)),
    reject_abundance  = as.vector(t(rej_ab)),
    reject_prevalence = as.vector(t(rej_pr)),
    stringsAsFactors  = FALSE )

  nan <- function (v) { v[is.na(v)] <- NaN; stats::setNames(v, features) }

  attr(out, 'method')      <- 'maaslin3'
  attr(out, 'formula')     <- if (mixed) linda_formula(des$formula, rnd$random) else des$formula
  attr(out, 'groups')      <- des$groups
  attr(out, 'design')      <- X
  if (mixed) attr(out, 'random') <- rnd$random
  attr(out, 'norm')        <- norm
  attr(out, 'transform')   <- transform
  attr(out, 'alpha')       <- alpha
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'log_base')    <- if (transform == 'log') 2 else if (norm == 'rclr') exp(1) else NA_real_
  attr(out, 'mean')        <- nan(ab$mean)
  attr(out, 'prevalence')  <- nan(ab$n_obs / n)
  attr(out, 'resid_sd')    <- nan(ab$resid_sd)
  if (mixed) attr(out, 'random_sd') <- nan(ab$theta * ab$resid_sd)
  class(out) <- c('ecodive_maaslin3', 'data.frame')

  list(
    result     = out,
    abundance  = c(ab_raw, list(median = mca$median, sd_median = mca$sd_median, cov = mca$cov)),
    prevalence = c(m3_table(pr, keep, 'prevalence', mixed && !small, covariates),
                   list(median = mcp$median, sd_median = mcp$sd_median, cov = mcp$cov)) )
}



# A single choice among `opts`, by name.
m3_choice <- function (x, what, opts) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !tolower(x) %in% opts)
    stop('`', what, '` must be ', paste0("'", opts[-length(opts)], "'", collapse = ', '),
         ", or '", opts[[length(opts)]], "'.", call. = FALSE)
  tolower(x)
}


# One model's fits as D x K matrices (coef, se, df, p, error, null), for
# the columns `keep` of the design.
m3_table <- function (fit, keep, model, mixed, names) {

  D <- nrow(fit$coef)
  K <- length(keep)
  coef <- fit$coef[, keep, drop = FALSE]
  se   <- fit$se[, keep, drop = FALSE]
  df   <- fit$df[, keep, drop = FALSE]
  p    <- 2 * stats::pt(-abs(coef / se), df)

  st  <- fit$status
  txt <- rep(NA_character_, D)
  txt[st == 1L] <- if (model == 'abundance') M3_TEXT[['same_abundance']] else M3_TEXT[['same_prevalence']]
  txt[st == 2L] <- M3_TEXT[['baseline']]
  txt[st == 3L] <- M3_TEXT[['one_group']]
  txt[st == 4L] <- M3_TEXT[['many_groups']]
  txt[st == 5L] <- M3_TEXT[['failed']]
  txt[st == 6L] <- M3_TEXT[['no_conv']]
  txt[st == 7L] <- M3_TEXT[['one_level']]
  txt[fit$warn == 1L] <- M3_TEXT[['glm_conv']]
  txt[bitwAnd(fit$warn, 2L) > 0L] <- M3_TEXT[['glm_01']]

  error <- matrix(txt, D, K)
  error[is.na(error) & is.na(p)] <- M3_TEXT[['na_p']]

  # A failed fit has no estimates.
  fail <- st %in% c(1L, 2L, 3L, 4L, 5L, 7L)
  coef[fail, ] <- se[fail, ] <- df[fail, ] <- p[fail, ] <- NA

  out <- list(coef = coef, se = se, df = df, p = p, error = error, null = matrix(0, D, K))
  lapply(out, function (m) { colnames(m) <- names; m })
}


# Adjusts the p-values of two models together, leaving out those with an
# error: list(a, b) of q-values.
m3_adjust <- function (ta, tb, p_adjust) {
  pa <- ta$p; pa[!is.na(ta$error)] <- NA
  pb <- tb$p; pb[!is.na(tb$error)] <- NA
  p  <- c(as.vector(pa), as.vector(pb))
  q  <- p
  ok <- !is.na(p)
  q[ok] <- stats::p.adjust(p[ok], p_adjust)
  list(a = matrix(q[seq_along(pa)], nrow(pa)), b = matrix(q[-seq_along(pa)], nrow(pb)))
}


# The median comparison of one model's coefficients (maaslin3's
# run_median_comparison()): each covariate's coefficients are tested
# against their median over the features whose p-value is below 0.95,
# allowing for the median's own variance and its covariance with each.
m3_median <- function (tab, xsd, threshold, subtract, n_sims, seed, model, sims, covariates, cpus) {

  K <- ncol(tab$coef)
  med <- sdm <- stats::setNames(rep(NA_real_, K), covariates)
  cvl <- stats::setNames(vector('list', K), covariates)
  tab$null[] <- NA

  for (k in seq_len(K)) {

    ok <- which(is.na(tab$error[, k]))
    if (!length(ok)) next
    b  <- tab$coef[ok, k]
    s  <- tab$se[ok, k]
    p0 <- tab$p[ok, k]
    use <- !is.na(p0) & p0 < 0.95
    m   <- stats::median(b[use])
    med[[k]] <- m

    if (is.na(m)) {
      pn <- rep(NA_real_, length(ok))
    }
    else {
      sd_median <- sqrt(0.25 * 2 * base::pi * stats::var(b[use]) / sum(use))
      sdm[[k]] <- sd_median
      z <- if (!is.null(sims)) sims[[covariates[[k]]]]
      if (!is.null(z)) {
        z <- as.matrix(z)
        storage.mode(z) <- 'double'
        if (nrow(z) != length(ok) || ncol(z) != n_sims)
          stop('`sims` for ', covariates[[k]], ' must be ', length(ok), ' x ', n_sims, '.', call. = FALSE)
      }
      mc <- .Call(
        C_maaslin3_median, as.double(b), as.double(s), use, n_sims, seed,
        as.integer(model * 65536L + k), z, as.double(m), cpus )
      cvl[[k]] <- mc$cov

      offsets <- abs(m - b) * sqrt(s^2 / (s^2 + sd_median^2 - 2 * mc$cov)) + b
      thr <- ifelse(is.na(xsd[ok, k]), threshold, threshold / xsd[ok, k])
      tt  <- (b - offsets) / s
      pn  <- 2 * stats::pt(-abs(tt), tab$df[ok, k])
      pn[!is.na(b) & abs(b - m) < thr] <- 1
      pn[is.na(b)] <- NA
      if (subtract) tab$coef[ok, k] <- b - m
    }

    became <- is.na(pn) & !is.na(p0)
    tab$error[ok[became], k] <- M3_TEXT[['median_na']]
    tab$p[ok, k]    <- pn
    tab$null[ok, k] <- m
  }

  list(table = tab, median = med, sd_median = sdm, cov = cvl)
}



#' Plot a MaAsLin 3 Result
#'
#' Views of a [maaslin3()] result for one covariate and one of its two
#' models, each feature a point (or a row), with the features that the
#' model declares associated (`reject_abundance` or `reject_prevalence`)
#' highlighted and the most significant of them named:
#'
#' * `"volcano"` (the default): the coefficient against its adjusted
#'   p-value.
#' * `"ma"`: the coefficient against the feature's mean transformed
#'   abundance where present.
#' * `"forest"`: the coefficients and their confidence intervals of the
#'   `n_features` most significant features.
#' * `"rank"`: every feature's coefficient, sorted (see **Rank plots**).
#'
#' @param x  The result of [maaslin3()].
#'
#' @param covariate  Which covariate to plot, as named in the result's
#'        `covariate` column (such as `"statusmoderate"` or `"age"`).
#'        Defaults to the first column of the term of interest (`groups`),
#'        or else to the first covariate.
#'
#' @param model  `"abundance"` (the default) or `"prevalence"`.
#'
#' @param type  The plot: `"volcano"` (the default), `"ma"`, `"forest"`,
#'        or `"rank"`.
#'
#' @inheritParams plot.ecodive_linda
#'
#' @inheritSection plot.ecodive_ancombc Rank plots
#'
#' @details
#' Abundance coefficients are in the units of the result's `norm` and
#' `transform` (log2 by default, its `log_base`); prevalence coefficients
#' are log odds ratios.
#'
#' * **Volcano plot.** The y axis is \eqn{-\log_{10}} of the model's
#'   adjusted p-value (`q_abundance` or `q_prevalence`; 0 is drawn at the
#'   smallest positive double), with a dashed line at `alpha`.
#' * **MA plot.** The x axis is the feature's mean transformed abundance
#'   where present (the `mean` attribute).
#' * **Forest plot.** Bars over \eqn{\mathrm{coef} \pm q\,\mathrm{se}},
#'   with \eqn{q} the \eqn{(1 + }`conf_level`\eqn{)/2} quantile of the t
#'   distribution with `dof_abundance` degrees of freedom (abundance), or
#'   of the normal distribution (prevalence); they aren't adjusted for
#'   testing many features. With a median comparison, the coefficients were
#'   tested against the median (named in the subtitle), not against 0.
#'
#' Features whose values are `NaN` are left out. The subtitle names the
#' covariate and the model, and counts the features rejected.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     x <- matrix(rpois(30 * 25, 10) * rbinom(30 * 25, 1, 0.7), 30, 25)
#'     x[16:30, 1:3] <- x[16:30, 1:3] * 6
#'     res <- maaslin3(x, groups = rep(c('a', 'b'), each = 15))
#'     plot(res)
#'     plot(res, model = 'prevalence')
#'     plot(res, type = 'forest', n_features = 10)
#'     plot(res, type = 'rank', numerator = 1:3)
plot.ecodive_maaslin3 <- function (
    x, covariate = NULL, model = c('abundance', 'prevalence'),
    type = c('volcano', 'ma', 'forest', 'rank'), labels = 10, n_features = 30,
    conf_level = 0.95, numerator = NULL, denominator = NULL, ... ) {

  type  <- test_plot_type(type, c('volcano', 'ma', 'forest', 'rank'))
  model <- test_plot_type(model, c('abundance', 'prevalence'))

  covariate <- da_rank_covariate(x, covariate)

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)

  r   <- x[x$covariate == covariate, , drop = FALSE]
  ab  <- model == 'abundance'
  est <- if (ab) r$coef_abundance else r$coef_prevalence
  se  <- if (ab) r$se_abundance else r$se_prevalence
  q   <- if (ab) suppressWarnings(stats::qt((1 + conf_level) / 2, r$dof_abundance))
         else stats::qnorm((1 + conf_level) / 2)
  nul <- if (ab) r$null_abundance else r$null_prevalence
  nul <- unique(nul[!is.na(nul)])

  est_lab <- if (!ab) 'Log odds ratio (prevalence)'
    else {
      b <- attr(x, 'log_base', exact = TRUE)
      if (isTRUE(b == 2)) 'Coefficient (abundance, log2)'
      else if (isTRUE(abs(b - exp(1)) < 1e-12)) 'Coefficient (abundance, natural log)'
      else 'Coefficient (abundance)'
    }

  da_plot(
    x, type, labels, n_features, list(...),
    est = est, p = if (ab) r$q_abundance else r$q_prevalence,
    reject = if (ab) r$reject_abundance else r$reject_prevalence,
    feature = r$feature,
    mean = attr(x, 'mean', exact = TRUE),
    lo = est - q * se, hi = est + q * se,
    what     = paste0(covariate, ', ', model),
    est_lab  = est_lab,
    mean_lab = 'Mean abundance where present (transformed)',
    interval = sprintf('unadjusted %s confidence intervals', plot_pct(conf_level)),
    num = numerator, den = denominator,
    extra = if (length(nul) == 1L && nul != 0) sprintf('tested against the median, %.3g', nul) )

  invisible(x)
}
