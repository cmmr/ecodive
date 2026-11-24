# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# neutral_model(): Sloan et al.'s (2006) neutral community model, fitted by
# least squares as Burns et al. (2016) do, with src/neutral_model.c for the
# per-feature sums. The model function (the probability that a feature of
# metacommunity abundance p is detected above 1/N), the definitions of the
# mean abundance and the occurrence frequency, the Wilson band around the
# fit and the above / below rule follow neufit.py (Michael Sieber; notice
# below), restructured for ecodive, see plans/neutral_model.md: m is fitted
# by Brent's method over log m and polished by solving for a zero of the
# gradient, without neufit's bound m <= 1; neufit's rarefaction step is
# left to the user (rarefy()). Burns et al.'s code (no licence stated) was
# used only through its outputs (iCAMP::snm(), tyRa's fit_sncm()), never
# read.
#
# neufit's notice:
#
#   neufit: Fit a neutral community model to species abundances, e.g. from
#   an OTU table
#
#   Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
#   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' Sloan's Neutral Community Model
#'
#' Fits Sloan et al.'s (2006) neutral community model to a count table, as
#' Burns et al. (2016) do: each feature's **occurrence frequency** (the
#' share of samples it is detected in) against its **mean relative
#' abundance**, with one free parameter, the immigration probability
#' \eqn{m}. Features detected more often than the fit predicts (above its
#' confidence band) may be favoured by the habitat or disperse well; those
#' detected less often (below it) may be dispersal-limited or selected
#' against.
#'
#' @inheritParams documentation
#'
#' @param counts  A count table (samples \eqn{\times} features) of whole
#'        numbers, or any of the objects listed under **Input Types**,
#'        **rarefied to equal depth** (see [rarefy()]). The depth \eqn{N}
#'        sets the detection limit \eqn{1/N}, so proportions are not
#'        accepted. With unequal depths, \eqn{N} is their mean, with a
#'        warning (as Burns et al.'s code does). Samples with no counts are
#'        left out, with a warning.
#'
#' @param groups  Optional: fit the model separately within each group of
#'        samples (a habitat, a host stage, ...). One value per sample,
#'        matched to the samples by name when both are named, else in the
#'        samples' order; or, with `data`, the name of one of its columns.
#'        Samples whose group is `NA` are left out. Each group needs at
#'        least 2 samples. Default: `NULL`, one fit to every sample.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `groups` as a column. Rows are matched to the
#'        samples in `counts` by row name; without row names, they must be in
#'        the same order as the samples. Default: `NULL`.
#'
#' @param pool  Optional count table (samples \eqn{\times} features, the
#'        same `margin` as `counts`) of the **source metacommunity**: each
#'        feature's mean relative abundance \eqn{p} is its share of the
#'        pool's counts instead of the share of the counts in `counts`.
#'        Matched to the features of `counts` by name (or by position when
#'        either has no feature names); every feature present in `counts`
#'        must be in `pool`. Default: `NULL`, the samples themselves.
#'
#' @param conf_level  The confidence level of the interval for \eqn{m} and
#'        of the Wilson band around the fit: a number strictly between 0
#'        and 1. Default: `0.95`, as Burns et al. and neufit.
#'
#' @return A data frame of class `c("ecodive_neutral_model", "data.frame")`,
#'   one row per feature fitted (with `groups`, per group and feature),
#'   ordered by mean abundance (within each group), with columns:
#'   \describe{
#'     \item{`group`}{The group (a factor; only with `groups`).}
#'     \item{`feature`}{The feature's name (its index if it has none).}
#'     \item{`mean_abundance`}{\eqn{p}, its mean relative abundance in the
#'       metacommunity: its share of all the counts (of the group, or of
#'       `pool`).}
#'     \item{`occurrence`}{The share of the samples it is present in.}
#'     \item{`predicted`}{The neutral model's occurrence at the fitted
#'       \eqn{m}.}
#'     \item{`lower`, `upper`}{The Wilson interval (Wilson 1927) of a
#'       proportion `predicted` among `n_samples` samples, at `conf_level`.}
#'     \item{`partition`}{`"above"` (`occurrence > upper`), `"below"`
#'       (`occurrence < lower`) or `"neutral"`: a factor. A difference
#'       must exceed \eqn{10^{-12}}, so that rounding doesn't decide at the
#'       band's edges (at occurrence 1, the upper end is 1 to within
#'       rounding).}
#'     \item{`predicted_binomial`, `predicted_poisson`}{The occurrence that
#'       random sampling alone predicts: \eqn{1 - (1 - p)^N} (binomial) and
#'       \eqn{1 - e^{-Np}} (Poisson), with no parameter to fit.}
#'   }
#'   Features that are absent from every sample (of the group), or whose
#'   \eqn{p} is 0, can't be fitted and are left out (as Burns et al. do);
#'   their names are in attribute `dropped`.
#'
#'   Attributes, each a number (with `groups`, a named vector, one value per
#'   group):
#'   \describe{
#'     \item{`method`}{`"neutral_model"`.}
#'     \item{`m`, `m_lower`, `m_upper`}{The least-squares immigration
#'       probability and its profile confidence interval at `conf_level`
#'       (0 or `Inf` when the residual sum of squares stays within the
#'       cutoff on that side).}
#'     \item{`m_se`}{The standard error of `m` by the delta method, as
#'       neufit (lmfit) reports it.}
#'     \item{`nm`}{\eqn{Nm}, the number of immigrants per generation.}
#'     \item{`n_individuals`}{\eqn{N}, the reads per sample (their mean).}
#'     \item{`detection`}{The detection limit \eqn{1/N}.}
#'     \item{`n_samples`}{The number of samples fitted.}
#'     \item{`r_squared`, `rmse`, `aic`, `bic`}{The fit's \eqn{R^2}
#'       (\eqn{1 - SS_{res} / SS_{tot}}; it can be negative, and is `NA`
#'       when every feature has the same occurrence), the root mean squared
#'       error \eqn{\sqrt{SS_{res} / (K - 1)}} for \eqn{K} features, and
#'       the AIC and BIC of a Gaussian likelihood of the residuals (2
#'       parameters: \eqn{m} and the variance).}
#'     \item{`r_squared_binomial`, `rmse_binomial`, `aic_binomial`,
#'       `bic_binomial`, and the same with `_poisson`}{The same for the
#'       binomial and Poisson models, whose likelihoods have a fitted mean
#'       and variance (2 parameters), as Burns et al.'s. A lower AIC than
#'       the neutral model's means random sampling alone fits better.}
#'     \item{`conf_level`}{As given.}
#'     \item{`dropped`}{The features left out (with `groups`, a list by
#'       group).}
#'   }
#'
#' @details
#' **The model.** In a local community of \eqn{N} individuals that loses
#' one individual at random each step and replaces it, with probability
#' \eqn{m}, by an immigrant from a metacommunity where the feature has
#' relative abundance \eqn{p} (and otherwise by local reproduction), the
#' feature's local relative abundance \eqn{x} is, at equilibrium,
#' Beta-distributed with parameters \eqn{Nmp} and \eqn{Nm(1 - p)} (Sloan
#' et al. 2006). A feature is detected when \eqn{x > 1/N}, so its expected
#' occurrence is
#' \deqn{F(p) = 1 - I_{1/N}(Nmp,\; Nm(1 - p)),}
#' with \eqn{I} the regularized incomplete beta function
#' ([stats::pbeta()]). Small \eqn{m} (rare immigration) lets drift spread
#' the feature's abundance among communities, so it is often missing; large
#' \eqn{m} keeps every community close to the metacommunity, and the curve
#' approaches random sampling.
#'
#' **The fit.** \eqn{N} is fixed at the reads per sample, and \eqn{m} is
#' the value that minimizes the residual sum of squares \eqn{\sum_i
#' (f_i - F(p_i))^2} over the \eqn{K} features: Brent's method
#' ([stats::optimize()]) over \eqn{\log m} in \eqn{[10^{-12}, 10^4]},
#' started from the best of a grid of 61 values (a warning when two minima
#' are equally good), then polished by solving for the zero of the
#' gradient ([stats::uniroot()]). The profile interval holds the \eqn{m}
#' where \eqn{(SS(m) - SS(\hat m)) / s^2} is at most the `conf_level`
#' quantile of \eqn{F_{1, K-1}}, with \eqn{s^2 = SS(\hat m) / (K - 1)}:
#' what R's `confint()` gives for an `nls()` fit (Burns et al.'s
#' interval), with its ends solved for exactly instead of interpolated.
#' It assumes residuals of equal variance; occurrences are binomial, so in
#' simulations of the model its coverage is somewhat below `conf_level`
#' (about 87\% for a 95\% interval).
#'
#' **What it means.** The fit describes how far occurrence follows
#' abundance as drift and dispersal alone would make it; it is not a test of
#' neutrality, and a high \eqn{R^2} doesn't show that a community is
#' neutral. The model assumes equal sampling depths, a single
#' metacommunity, and communities at equilibrium. Features outside the band
#' are candidates for selection or dispersal effects, not proof of them.
#' \eqn{m} is a probability, but least squares can give \eqn{m > 1} on
#' near-saturated data (most features in most samples); it is then
#' reported, with a warning.
#'
#' @section Differences from Burns et al.'s code and neufit:
#' * \eqn{m} is the exact least-squares value (it agrees with a 40-digit
#'   computation to about \eqn{10^{-13}}). Burns et al.'s
#'   Levenberg-Marquardt `nlsLM()` and neufit's lmfit stop earlier: their
#'   \eqn{m} agree with it to about \eqn{10^{-5}} and \eqn{10^{-6}}.
#' * The interval for \eqn{m} has exact ends (R's `confint()` interpolates a
#'   sampled profile, so Burns et al.'s ends match to about \eqn{10^{-5}}).
#'   Burns et al. report the lower half-width only.
#' * \eqn{m} is not bounded above (neufit bounds it to \eqn{[0, 1]}), as in
#'   Burns et al.'s code, but a warning says when it exceeds 1.
#' * The AICs and BICs are computed in closed form (Burns et al. maximize
#'   the likelihoods numerically, with `stats4::mle()`, which agrees to
#'   about \eqn{10^{-6}}).
#' * The binomial prediction is \eqn{1 - (1 - p)^N} even when \eqn{N}, a
#'   mean depth, isn't a whole number.
#' * With `pool`, \eqn{p} is the feature's share of all the pool's counts.
#'   Burns et al. divide its mean count in the pool by the \eqn{N} of
#'   `counts`, which is the same when the pool is rarefied to the same
#'   depth.
#' * neufit rarefies the table itself (at random); `neutral_model()` warns
#'   about unequal depths and leaves rarefying to [rarefy()].
#' * A feature is above or below the band only by more than
#'   \eqn{10^{-12}}; the references' strict comparisons can flip on rounding
#'   at occurrence 1.
#' * `groups` fits each group on its own, as subsetting `counts` would.
#' * iCAMP's `snm()`, which builds on Burns et al.'s code, rounds a mean
#'   depth \eqn{N} to a whole number, and keeps features whose \eqn{p} in
#'   the pool is 0; Burns et al.'s `fit_sncm()` (in tyRa) and
#'   `neutral_model()` do neither.
#'
#' MicEco's `neutral.fit()` fits \eqn{m} by maximum likelihood (of a
#' binomial model of the occurrences) instead of least squares, so its
#' \eqn{m} differs a little from these.
#'
#' @section Same names in other packages:
#' The model is also known as the neutral community model (NCM) or Sloan's
#' neutral community model (SNCM): MicEco's `neutral.fit()`, tyRa's
#' `fit_sncm()` and iCAMP's `snm()` fit it too.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S., & Curtis,
#' T. P. (2006). Quantifying the roles of immigration and chance in shaping
#' prokaryote community structure. *Environmental Microbiology*, 8,
#' 732–740. \doi{10.1111/j.1462-2920.2005.00956.x}
#'
#' Burns, A. R., Stephens, W. Z., Stagaman, K., Wong, S., Rawls, J. F.,
#' Guillemin, K., & Bohannan, B. J. M. (2016). Contribution of neutral
#' processes to the assembly of gut microbial communities in the zebrafish
#' over host development. *The ISME Journal*, 10, 655–664.
#' \doi{10.1038/ismej.2015.142}
#'
#' Sieber, M., Pita, L., Weiland-Bräuer, N., et al. (2019). Neutrality in
#' the metaorganism. *PLoS Biology*, 17, e3000298.
#' \doi{10.1371/journal.pbio.3000298}
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference. *Journal of the American Statistical
#' Association*, 22, 209–212. \doi{10.1080/01621459.1927.10502953}
#'
#' @seealso [plot.ecodive_neutral_model()] for the fitted curve, and
#'   [rarefy()].
#'
#' @export
#' @examples
#'     # A neutral community: each sample draws its composition around the
#'     # metacommunity's (p), with immigration m = 0.05, then 1000 reads.
#'     set.seed(1)
#'     p <- rgamma(60, 0.3); p <- p / sum(p)
#'     counts <- t(replicate(30, {
#'       x <- rgamma(60, 1000 * 0.05 * p)
#'       as.vector(rmultinom(1, 1000, x / sum(x)))
#'     }))
#'     dimnames(counts) <- list(paste0('S', 1:30), paste0('F', 1:60))
#'
#'     res <- neutral_model(counts)
#'     res
#'     attr(res, 'm')
#'     plot(res)
#'     plot(res, 'binomial')
#'
#'     # One fit per group of samples
#'     res <- neutral_model(counts, groups = rep(c('A', 'B'), each = 15))
#'     attr(res, 'm')
#'     plot(res)
#'
neutral_model <- function (
    counts, groups = NULL, data = NULL, pool = NULL, margin = 1L,
    conf_level = 0.95, cpus = n_cpus() ) {

  info <- neutral_input()

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !(conf_level > 0 && conf_level < 1))
    stop('`conf_level` must be a single number strictly between 0 and 1.', call. = FALSE)
  conf_level <- as.double(conf_level)

  n  <- info$n
  D  <- info$D
  features <- if (is.null(info$features)) as.character(seq_len(D)) else info$features

  # Groups of samples (one group without `groups`); NA groups are left out.
  grouped <- !is.null(groups)
  if (grouped) {
    src <- if (!is.null(info$data) && is.character(groups) && length(groups) == 1 &&
               !is.na(groups) && groups %in% names(info$data)) groups
    g <- per_sample(groups, 'groups', info$data, info$idx, info$labels, n, ref = '`counts`')
    f <- if (is.factor(g)) droplevels(g) else factor(g)
    if (nlevels(f) == 0)
      stop(
        if (is.null(src)) '`groups`' else paste0('`data` column ', name_list(src)),
        ' has no groups: every value is NA.', call. = FALSE )
    lev  <- levels(f)
    code <- ifelse(is.na(f), -1L, as.integer(f) - 1L)
  }
  else {
    lev  <- NULL
    code <- rep(0L, n)
  }

  # Samples with no counts are left out.
  depth <- info$depth
  empty <- depth == 0 & code >= 0L
  if (all(depth[code >= 0L] == 0))
    stop('Every sample in `counts` is empty (no counts).', call. = FALSE)
  if (any(empty)) {
    warning(
      'Left out ', sum(empty), ' sample', if (sum(empty) == 1) '' else 's',
      ' with no counts: ', name_list(if (is.null(info$labels)) which(empty) else info$labels[empty]),
      '.', call. = FALSE )
    code[empty] <- -1L
  }

  G  <- max(1L, length(lev))
  st <- .Call(C_neutral_stats, info$counts, info$margin, as.integer(code), G, info$cpus)

  # The pool's relative abundances, in the order of `counts`' features.
  p_pool <- if (!is.null(pool)) neutral_pool(pool, info, features, colSums(st$occ) > 0)

  fits  <- vector('list', G)
  rows  <- vector('list', G)
  uneq  <- character(0)
  for (k in seq_len(G)) {

    gl   <- if (grouped) lev[[k]]
    in_g <- code == (k - 1L)
    what <- if (grouped) paste0(' in group ', name_list(gl)) else ''
    n_g  <- sum(in_g)
    if (n_g < 2L)
      stop(
        'The neutral model needs at least 2 samples', what, ', but ',
        if (n_g == 1L) 'there is 1.' else 'there are none.', call. = FALSE )

    d_g <- depth[in_g]
    N   <- mean(d_g)
    if (any(d_g != d_g[[1]])) uneq <- c(uneq, if (grouped) gl else '')

    occ <- st$occ[k, ] / n_g
    p   <- if (is.null(p_pool)) st$total[k, ] / sum(d_g) else p_pool
    keep <- p > 0 & occ > 0
    if (sum(keep) < 2L)
      stop(
        'The neutral model needs at least 2 features with counts', what,
        ', but there ', if (sum(keep) == 1L) 'is 1.' else 'are none.', call. = FALSE )

    ord <- which(keep)[order(p[keep], method = 'radix')]
    fit <- neutral_fit(p[ord], occ[ord], N, n_g, conf_level, what)
    fit$dropped <- features[!keep]

    rows[[k]] <- data.frame(
      feature            = features[ord],
      mean_abundance     = p[ord],
      occurrence         = occ[ord],
      predicted          = fit$predicted,
      lower              = fit$lower,
      upper              = fit$upper,
      partition          = fit$partition,
      predicted_binomial = fit$predicted_binomial,
      predicted_poisson  = fit$predicted_poisson,
      stringsAsFactors = FALSE )
    if (grouped)
      rows[[k]] <- cbind(
        data.frame(group = factor(rep(gl, length(ord)), levels = lev)), rows[[k]] )
    fits[[k]] <- fit
  }

  if (length(uneq))
    warning(
      'neutral_model() assumes equal depths (see rarefy()), but the samples',
      if (grouped) paste0(' in group', if (length(uneq) > 1) 's', ' ', name_list(uneq)),
      ' have different totals; N is their mean.', call. = FALSE )

  res <- do.call(rbind, rows)
  rownames(res) <- NULL

  one <- function (nm) {
    v <- vapply(fits, function (f) as.double(f[[nm]]), double(1))
    if (grouped) names(v) <- lev
    v
  }

  attr(res, 'method')        <- 'neutral_model'
  for (nm in c('m', 'm_lower', 'm_upper', 'm_se', 'nm', 'n_individuals', 'detection'))
    attr(res, nm) <- one(nm)
  attr(res, 'n_samples')     <- as.integer(one('n_samples'))
  if (grouped) names(attr(res, 'n_samples')) <- lev
  for (nm in c(
      'r_squared', 'rmse', 'aic', 'bic',
      'r_squared_binomial', 'rmse_binomial', 'aic_binomial', 'bic_binomial',
      'r_squared_poisson', 'rmse_poisson', 'aic_poisson', 'bic_poisson' ))
    attr(res, nm) <- one(nm)
  attr(res, 'conf_level')    <- conf_level
  attr(res, 'dropped')       <- if (grouped) stats::setNames(lapply(fits, `[[`, 'dropped'), lev)
                                else fits[[1]]$dropped
  class(res) <- c('ecodive_neutral_model', 'data.frame')
  res
}



#------------------------------------------------------------------------------
# Checks neutral_model()'s `counts`, `margin`, `cpus` and `data` (in the
# caller's environment, as validate_*() do). Returns list(counts, margin,
# cpus, n, D, labels, features, data, idx (data's row for each sample),
# depth (each sample's total)).
#------------------------------------------------------------------------------
neutral_input <- function (env = parent.frame()) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- neutral_check_counts(env$counts, '`counts`')
  margin <- env$margin

  smp  <- count_samples(counts, margin)
  dn   <- switch(
    get_matrix_package(counts),
    base = dimnames(counts), slam = counts$dimnames, Matrix = counts@Dimnames )
  dims <- switch(
    get_matrix_package(counts),
    base = dim(counts), slam = c(counts$nrow, counts$ncol), Matrix = counts@Dim )
  D <- as.integer(dims[[3L - margin]])
  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])
  if (D == 0L)
    stop('`counts` has no features.', call. = FALSE)

  data <- env$data
  idx  <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(smp$labels, smp$n, rn, nrow(data), ref = '`counts`')
  }

  list(
    counts = counts, margin = margin, cpus = env$cpus,
    n = smp$n, D = D, labels = smp$labels, features = features,
    data = data, idx = idx, depth = sample_totals(counts, margin) )
}


# A validated count table, checked to hold whole numbers >= 0 (logical
# matrices become numbers). `what` names it in messages.
neutral_check_counts <- function (counts, what) {

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop(what, ' must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop(what, ' must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(vals < 0))
    stop(what, ' must not contain negative values.', call. = FALSE)
  if (any(vals != round(vals)))
    stop(
      what, ' must be counts (whole numbers), not proportions: the reads per ',
      'sample, N, set the detection limit 1/N.', call. = FALSE )
  counts
}


# `pool`'s relative abundances (its features' shares of all its counts), for
# the features of `counts` (`features`; `present`: whether each occurs in
# `counts`).
neutral_pool <- function (pool, info, features, present) {

  env <- new.env()
  env$counts <- pool
  env$margin <- info$margin
  tryCatch(
    { validate_counts(env); validate_margin(env) },
    error = function (e) stop('`pool`: ', conditionMessage(e), call. = FALSE) )
  pool   <- neutral_check_counts(env$counts, '`pool`')
  margin <- env$margin

  dn <- switch(
    get_matrix_package(pool),
    base = dimnames(pool), slam = pool$dimnames, Matrix = pool@Dimnames )
  pf <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])

  n_pool <- count_samples(pool, margin)$n
  st     <- .Call(C_neutral_stats, pool, margin, rep(0L, n_pool), 1L, info$cpus)
  total  <- st$total[1, ]
  if (!(sum(total) > 0))
    stop('`pool` is empty (no counts).', call. = FALSE)
  rel <- total / sum(total)

  if (!is.null(info$features) && !is.null(pf)) {
    if (anyDuplicated(pf))
      stop('`pool` has duplicate feature names.', call. = FALSE)
    idx  <- match(features, pf)
    miss <- is.na(idx) & present
    if (any(miss))
      stop(
        sum(miss), ' feature', if (sum(miss) == 1) '' else 's', ' in `counts` ',
        if (sum(miss) == 1) 'is' else 'are', ' not in `pool`: ',
        name_list(features[miss]), '.', call. = FALSE )
    out <- rel[idx]
    out[is.na(out)] <- 0
    return (out)
  }

  if (length(rel) != length(features))
    stop(
      sprintf(
        '`pool` has %d features, but `counts` has %d. Name the features of both to match them by name.',
        length(rel), length(features) ),
      call. = FALSE )
  rel
}



#------------------------------------------------------------------------------
# The fit, for K features with mean abundances p (> 0) and occurrences occ
# (> 0), N reads per sample and n samples. `what` (" in group 'A'", or "")
# goes into warnings. Returns a list of the per-feature columns and the
# statistics.
#------------------------------------------------------------------------------
NEUTRAL_M_RANGE <- c(1e-12, 1e4)

neutral_fit <- function (p, occ, N, n, conf_level, what = '') {

  K   <- length(p)
  fit <- neutral_fit_m(p, occ, N, what)
  m   <- fit$m

  pred <- neutral_curve(p, N, m)
  res  <- occ - pred
  rss  <- sum(res ^ 2)
  sst  <- sum((occ - mean(occ)) ^ 2)

  if (!(rss > 0))
    warning(
      'The model fits every occurrence', what, ' exactly (every feature ',
      'that is in every sample is predicted to be), so m is not identifiable: ',
      'any larger m fits as well. Its interval is NA.', call. = FALSE )
  if (!(sst > 0))
    warning(
      'Every feature', what, ' has the same occurrence (', format(occ[[1]]),
      '), so R-squared is undefined (NA).', call. = FALSE )
  if (m > 1)
    warning(
      'The fitted immigration probability m', what, ' is ', format(signif(m, 4)),
      ', above 1: the occurrences are closer to random sampling than the ',
      'model allows (near-saturated data?).', call. = FALSE )

  ci <- neutral_m_interval(p, occ, N, m, rss, conf_level)
  band <- neutral_wilson(pred, n, conf_level)

  partition <- neutral_partition(occ, band)

  pred_b <- -expm1(N * log1p(-p))
  pred_p <- -expm1(-N * p)

  c(
    list(
      predicted          = pred,
      lower              = band$lower,
      upper              = band$upper,
      partition          = factor(partition, levels = c('above', 'neutral', 'below')),
      predicted_binomial = pred_b,
      predicted_poisson  = pred_p,
      m                  = m,
      m_lower            = ci$lower,
      m_upper            = ci$upper,
      m_se               = ci$se,
      nm                 = N * m,
      n_individuals      = N,
      detection          = 1 / N,
      n_samples          = n ),
    neutral_gof(res, sst, K, mu = FALSE, ''),
    neutral_gof(occ - pred_b, sst, K, mu = TRUE, '_binomial'),
    neutral_gof(occ - pred_p, sst, K, mu = TRUE, '_poisson') )
}


# Above or below the band by more than NEUTRAL_EDGE, else neutral: at
# occurrence 1, the band's upper end is 1 to within rounding, and the strict
# comparison of Burns et al. and neufit would decide by the last bit.
NEUTRAL_EDGE <- 1e-12

neutral_partition <- function (occ, band) {
  ifelse(
    occ - band$upper > NEUTRAL_EDGE, 'above',
    ifelse(band$lower - occ > NEUTRAL_EDGE, 'below', 'neutral') )
}


# The neutral model's occurrence: the probability that a feature with
# metacommunity abundance p exceeds the detection limit 1/N (neufit's
# beta_cdf(), Sloan et al. 2006).
neutral_curve <- function (p, N, m) {
  stats::pbeta(1 / N, N * m * p, N * m * (1 - p), lower.tail = FALSE)
}


# R-squared, RMSE, AIC and BIC of residuals `r` (a Gaussian likelihood:
# with mean 0 for the neutral model, whose parameters are m and the
# variance; with a fitted mean for the binomial and Poisson models, whose
# parameters are the mean and the variance; 2 in each case).
neutral_gof <- function (r, sst, K, mu, suffix) {
  ss <- sum(r ^ 2)
  s2 <- if (mu) sum((r - mean(r)) ^ 2) / K else ss / K
  ll <- if (s2 > 0) -K / 2 * (log(2 * pi * s2) + 1) else NA_real_
  out <- list(
    r_squared = if (sst > 0) 1 - ss / sst else NA_real_,
    rmse      = if (K > 1) sqrt(ss / (K - 1)) else NA_real_,
    aic       = -2 * ll + 2 * 2,
    bic       = -2 * ll + 2 * log(K) )
  names(out) <- paste0(names(out), suffix)
  out
}


# The Wilson interval of a proportion `pred` among n trials (Wilson 1927),
# as statsmodels' proportion_confint(method = 'wilson') and Hmisc's
# binconf() give it. The lower end of x successes, (x + z^2/2 - z sqrt(S))
# / (n + z^2) with S = x (n - x) / n + z^2 / 4, is computed as x^2 / (n (x
# + z^2/2 + z sqrt(S))), the same without cancellation; the upper end is 1
# minus the lower end of n - x. So the ends are exactly 0 and 1 at x = 0
# and x = n.
neutral_wilson <- function (pred, n, conf_level) {
  z  <- stats::qnorm(1 - (1 - conf_level) / 2)
  z2 <- z * z
  low <- function (x) x * x / (n * (x + z2 / 2 + z * sqrt(pmax(0, x * (n - x) / n) + z2 / 4)))
  x  <- pred * n
  list(lower = low(x), upper = 1 - low(n - x))
}


# The least-squares m: the best of a grid over log m, refined by Brent's
# method around each of the grid's local minima, then polished by solving
# for the zero of the gradient. Returns list(m).
neutral_fit_m <- function (p, occ, N, what) {

  rss <- function (u) sum((occ - neutral_curve(p, N, exp(u))) ^ 2)

  lim  <- log(NEUTRAL_M_RANGE)
  grid <- seq(lim[[1]], lim[[2]], length.out = 61L)
  vals <- vapply(grid, rss, double(1))
  L    <- length(grid)

  # The grid's local minima (plateaus count once, at their first point).
  is_min <- vapply(seq_len(L), function (i) {
    v <- vals[[i]]
    (i == 1L || v < vals[[i - 1L]]) && (i == L || v <= vals[[i + 1L]])
  }, logical(1))
  cand <- which(is_min)

  opt <- lapply(cand, function (i) {
    br <- grid[c(max(1L, i - 1L), min(L, i + 1L))]
    o  <- stats::optimize(rss, br, tol = 1e-12)
    if (vals[[i]] < o$objective) list(minimum = grid[[i]], objective = vals[[i]]) else o
  })
  obj  <- vapply(opt, `[[`, double(1), 'objective')
  best <- which.min(obj)
  u    <- opt[[best]]$minimum
  r0   <- obj[[best]]

  # Another minimum, elsewhere, as good.
  far <- abs(vapply(opt, `[[`, double(1), 'minimum') - u) > 1e-3
  if (any(far & obj <= r0 + 1e-8 * max(r0, .Machine$double.xmin)))
    warning(
      'The residual sum of squares', what, ' has more than one minimum in m ',
      'that fit equally well (within 1e-8); m is the best of them.', call. = FALSE )

  if (u <= lim[[1]] + 1e-6 || u >= lim[[2]] - 1e-6) {
    warning(
      'The least-squares m', what, ' is at the edge of the search interval [',
      format(NEUTRAL_M_RANGE[[1]]), ', ', format(NEUTRAL_M_RANGE[[2]]), ']: ',
      'the fit keeps improving as m ', if (u < 0) 'falls' else 'grows',
      ', so m is not identifiable.', call. = FALSE )
    return (list(m = exp(u)))
  }

  # The gradient, d SS / d log m (up to a factor -2), with each feature's
  # dF / d log m by central differences with Richardson extrapolation.
  grad <- function (u) {
    h  <- 1e-3
    d1 <- (neutral_curve(p, N, exp(u + h))     - neutral_curve(p, N, exp(u - h)))     / (2 * h)
    d2 <- (neutral_curve(p, N, exp(u + h / 2)) - neutral_curve(p, N, exp(u - h / 2))) / h
    sum((occ - neutral_curve(p, N, exp(u))) * (4 * d2 - d1) / 3)
  }
  br <- u + c(-1e-3, 1e-3)
  gb <- c(grad(br[[1]]), grad(br[[2]]))
  if (all(is.finite(gb)) && gb[[1]] > 0 && gb[[2]] < 0) {
    root <- stats::uniroot(grad, br, f.lower = gb[[1]], f.upper = gb[[2]], tol = 1e-15)$root
    if (rss(root) <= r0 * (1 + 1e-12)) u <- root
  }

  list(m = exp(u))
}


# The profile interval for m: where (SS(m) - SS(m_hat)) / s2 reaches the
# conf_level quantile of F(1, K - 1); and the delta method's standard
# error, sqrt(s2 / sum((dF / dm)^2)). NA with fewer than 2 features or a
# perfect fit.
neutral_m_interval <- function (p, occ, N, m, rss, conf_level) {

  K  <- length(p)
  na <- list(lower = NA_real_, upper = NA_real_, se = NA_real_)
  if (K < 2L || !(rss > 0)) return (na)

  s2  <- rss / (K - 1)
  cut <- stats::qf(conf_level, 1, K - 1)
  tau <- function (u) sum((occ - neutral_curve(p, N, exp(u))) ^ 2) / s2 - (K - 1) - cut

  u   <- log(m)
  lim <- log(NEUTRAL_M_RANGE)
  end <- function (a, b, out) {
    fa <- tau(a); fb <- tau(b)
    if (!(fa * fb < 0)) return (out)
    exp(stats::uniroot(tau, c(a, b), f.lower = fa, f.upper = fb, tol = 1e-13)$root)
  }
  lower <- if (u > lim[[1]]) end(lim[[1]], u, 0)   else 0
  upper <- if (u < lim[[2]]) end(u, lim[[2]], Inf) else Inf

  h  <- 1e-4 * m
  J  <- (neutral_curve(p, N, m + h) - neutral_curve(p, N, m - h)) / (2 * h)
  se <- sqrt(s2 / sum(J ^ 2))

  list(lower = lower, upper = upper, se = se)
}



#------------------------------------------------------------------------------
# print(): the header lines.
#------------------------------------------------------------------------------
neutral_model_label <- function (x) {

  m <- attr(x, 'm', exact = TRUE)
  if (is.null(m)) return ('Neutral Community Model')

  nm  <- attr(x, 'nm', exact = TRUE)
  r2  <- attr(x, 'r_squared', exact = TRUE)
  lo  <- attr(x, 'm_lower', exact = TRUE)
  hi  <- attr(x, 'm_upper', exact = TRUE)
  cl  <- attr(x, 'conf_level', exact = TRUE)
  grp <- names(m)

  part <- if ('group' %in% names(x)) split(x$partition, x$group) else list(x$partition)

  lines <- vapply(seq_along(m), function (k) {
    tab <- table(part[[k]])
    sprintf(
      '%sm = %s (%s%% CI %s to %s), Nm = %s, R2 = %s; %d above, %d neutral, %d below',
      if (is.null(grp)) '' else paste0(grp[[k]], ': '),
      plot_num(m[[k]], 4L), format(100 * cl), plot_num(lo[[k]], 4L), plot_num(hi[[k]], 4L),
      plot_num(nm[[k]], 4L), plot_num(r2[[k]], 3L),
      tab[['above']], tab[['neutral']], tab[['below']] )
  }, character(1))

  c('Neutral Community Model', lines)
}



#' Plot a Neutral Community Model Fit
#'
#' The occurrence frequency of each feature against its mean relative
#' abundance (log scale), with the neutral model's fitted curve and its
#' Wilson band, as Burns et al. (2016) and neufit draw it. Points are
#' coloured by partition: above the band, within it, or below it.
#'
#' The curve and the band are drawn from the result's attributes
#' (`n_individuals`, `m`, `n_samples` and `conf_level`) over a grid of
#' abundances from a tenth of the smallest to 1, as neufit does.
#'
#' @param x  A [neutral_model()] result.
#'
#' @param type  `'neutral'` (the default): the neutral model's curve and
#'        band. `'binomial'`: the binomial (random sampling) model's curve
#'        too, \eqn{1 - (1 - p)^N}, dashed.
#'
#' @param group  With `groups`, which groups to draw (names or positions),
#'        one panel each. Default: `NULL`, all of them.
#'
#' @param legend  Where the legend goes, as a [graphics::legend()] keyword
#'        (`"topleft"`, ...), or `NULL` for none.
#'
#' @param ...  Passed on to [graphics::plot()], such as `main`, `pch`,
#'        `cex` or `xlim`.
#'
#' @return `x`, invisibly.
#'
#' @seealso [neutral_model()]
#'
#' @export
#' @examples
#'     set.seed(1)
#'     p <- rgamma(60, 0.3); p <- p / sum(p)
#'     counts <- t(replicate(30, {
#'       x <- rgamma(60, 1000 * 0.05 * p)
#'       as.vector(rmultinom(1, 1000, x / sum(x)))
#'     }))
#'
#'     res <- neutral_model(counts)
#'     plot(res)
#'     plot(res, 'binomial', legend = 'topleft')
#'
plot.ecodive_neutral_model <- function (
    x, type = c('neutral', 'binomial'), group = NULL, legend = 'bottomright', ...) {

  type <- test_plot_type(type, c('neutral', 'binomial'))
  pos  <- plot_legend_pos(legend)

  m <- attr(x, 'm', exact = TRUE)
  if (is.null(m) || is.null(attr(x, 'n_individuals', exact = TRUE)))
    stop(
      '`x` has lost the attributes of its fit (rows taken with `[` drop ',
      'them): plot the whole neutral_model() result.', call. = FALSE )

  grouped <- 'group' %in% names(x)
  lev     <- if (grouped) names(m) else NULL
  ks      <- if (grouped) which(lev %in% as.character(x$group)) else 1L
  if (!is.null(group)) {
    if (!grouped)
      stop('`group` is only for a result made with `groups`.', call. = FALSE)
    ks <- if (is.numeric(group)) group else match(group, lev)
    if (length(ks) == 0 || anyNA(ks) || any(ks < 1 | ks > length(lev) | ks != round(ks)))
      stop('`group` must name groups of `x`: ', name_list(lev), '.', call. = FALSE)
  }
  if (length(ks) == 0 || nrow(x) == 0)
    stop('`x` has no rows to plot.', call. = FALSE)

  if (length(ks) > 1L) {
    old <- graphics::par(mfrow = rev(grDevices::n2mfrow(length(ks))))
    on.exit(graphics::par(old), add = TRUE)
  }

  for (k in ks)
    neutral_plot_one(
      x = if (grouped) x[x$group == lev[[k]], , drop = FALSE] else x,
      at = function (nm) attr(x, nm, exact = TRUE)[[k]],
      cl = attr(x, 'conf_level', exact = TRUE),
      title = if (grouped) paste0('Neutral Community Model: ', lev[[k]]) else 'Neutral Community Model',
      type = type, pos = pos, ... )

  invisible(x)
}


neutral_plot_one <- function (x, at, cl, title, type, pos, ...) {

  args <- list(...)
  N <- at('n_individuals'); m <- at('m'); n <- at('n_samples')

  cols <- c(above = plot_palette(2L)[[1]], neutral = PLOT_MUTED, below = plot_palette(2L)[[2]])
  pmin_ <- min(x$mean_abundance)
  grid  <- 10 ^ seq(log10(pmin_ / 10), 0, length.out = 400L)
  curve <- neutral_curve(grid, N, m)
  band  <- neutral_wilson(curve, n, cl)

  sub <- sprintf(
    'm = %s, Nm = %s, R2 = %s', plot_num(m, 3L), plot_num(at('nm'), 3L),
    plot_num(at('r_squared'), 3L) )

  plot_do(
    graphics::plot, c(list(x$mean_abundance, x$occurrence), args),
    list(
      main = title, sub = sub, log = 'x', type = 'n', las = 1,
      xlim = range(grid), ylim = c(0, 1),
      xlab = 'Mean relative abundance', ylab = 'Occurrence frequency' ))
  graphics::polygon(c(grid, rev(grid)), c(band$lower, rev(band$upper)), col = 'grey90', border = NA)
  graphics::lines(grid, band$lower, lty = 2, col = 'grey40')
  graphics::lines(grid, band$upper, lty = 2, col = 'grey40')
  graphics::points(
    x$mean_abundance, x$occurrence, pch = plot_arg(args, 'pch', 16),
    cex = plot_arg(args, 'cex', 1), col = cols[as.character(x$partition)] )
  graphics::lines(grid, curve, lwd = 2, col = 'black')
  if (type == 'binomial')
    graphics::lines(grid, -expm1(N * log1p(-grid)), lwd = 2, lty = 3, col = PLOT_HIGHLIGHT)

  if (!is.null(pos)) {
    tab <- table(factor(x$partition, levels = names(cols)))
    lab <- sprintf('%s (%d)', c('Above', 'Neutral', 'Below'), as.integer(tab))
    graphics::legend(
      pos,
      legend = c(lab, 'Neutral model', if (type == 'binomial') 'Binomial model'),
      col    = c(cols, 'black', if (type == 'binomial') PLOT_HIGHLIGHT),
      pch    = c(16, 16, 16, NA, if (type == 'binomial') NA),
      lty    = c(NA, NA, NA, 1, if (type == 'binomial') 3),
      lwd    = c(NA, NA, NA, 2, if (type == 'binomial') 2),
      bg = 'white', cex = 0.8, inset = 0.01 )
  }
  invisible(NULL)
}
