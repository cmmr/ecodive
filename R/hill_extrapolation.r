# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Interpolated and extrapolated Hill numbers (iNEXT-style), written from the
# papers in a clean room: Chao et al. (2014), Chao & Jost (2012), Chao, Wang
# & Jost (2013) and Hsieh, Ma & Chao (2016). iNEXT (GPL) was run only for
# its outputs. The computations are in src/hill_extrapolation.c; the
# bootstrap's multinomial draws in src/rng_binomial.c (ported from NumPy).
# See plans/hill_extrapolation.md.



#' Interpolated and extrapolated Hill numbers
#'
#' Diversity grows with sequencing depth, so samples of unequal depth can't
#' be compared directly. `hill_extrapolation()` gives each sample's Hill
#' numbers of order \eqn{q = 0, 1, 2} (richness, Shannon diversity and
#' Simpson diversity) at a range of depths, by **rarefaction**
#' (interpolation) below the sample's size and **extrapolation** above it,
#' with the **sample coverage** at each depth and the **asymptotic**
#' estimates (Chao et al. 2014), as `iNEXT::iNEXT()` does.
#' `hill_standardize()` gives the Hill numbers of every sample at one
#' common depth or, by default, one common coverage (as
#' `iNEXT::estimateD()`): values that can be compared, and tested with
#' `alpha_test()`.
#'
#' @section Rarefaction, extrapolation and coverage:
#' For a sample of \eqn{n} observations with \eqn{f_k} features seen
#' \eqn{k} times:
#'
#' * **Rarefaction** to \eqn{m < n}: the expected Hill number of a random
#'   subsample of \eqn{m} of the observations (Chao et al. 2014, Table 1):
#'   richness \eqn{S_{obs} - \sum_i \binom{n - X_i}{m} / \binom{n}{m}} (as
#'   [rarefaction_curve()]), Shannon diversity from the expected frequency
#'   counts of the subsample, and Simpson diversity \eqn{1 / [1/m + (m-1)/m
#'   \sum_i X_i (X_i - 1) / (n (n - 1))]}.
#' * **Extrapolation** to \eqn{n + m^*}: richness \eqn{S_{obs} + \hat f_0
#'   [1 - (1 - f_1 / (n \hat f_0 + f_1))^{m^*}]}, with Chao1's \eqn{\hat
#'   f_0}; Shannon diversity \eqn{\exp[\frac{n}{n + m^*} H_{obs} +
#'   \frac{m^*}{n + m^*} \hat H]} (Chao et al. 2014, Eq. 10c), with Chao,
#'   Wang & Jost's (2013) entropy \eqn{\hat H}; Simpson diversity as for
#'   rarefaction, with \eqn{n + m^*} for \eqn{m}.
#' * **Sample coverage**, the share of the community's individuals that
#'   belong to features seen in a sample of that depth (Chao & Jost 2012),
#'   by its rarefied and extrapolated estimates, with Chao et al.'s (2014,
#'   Eq. 12) estimate at \eqn{n}.
#' * **Asymptotes**, the estimates for a complete sample: Chao1 (bias
#'   corrected), \eqn{\exp(\hat H)}, and \eqn{n (n - 1) / \sum_i X_i (X_i -
#'   1)}.
#'
#' Between two whole depths (when a coverage level falls between them),
#' every value is the straight line between its values at those depths, as
#' in iNEXT.
#'
#' **Reliability:** richness can be extrapolated reliably to about twice the
#' sample's size; Shannon and Simpson diversity, if the data aren't too
#' sparse, to their asymptotes (Chao et al. 2014). `max_depth` beyond twice
#' a sample's size is allowed, but read richness there with care.
#'
#' @section Equal depth or equal coverage:
#' Comparing samples at an equal **depth** says how diverse equally large
#' samples would be; comparing them at an equal **coverage** compares
#' equally *complete* samples, which Chao & Jost (2012) recommend: a more
#' diverse community needs more observations to reach the same coverage,
#' and comparing at equal depth understates its difference from a less
#' diverse one. `hill_standardize()` uses coverage unless you give
#' `depth`; its default coverage, the smallest among the samples
#' extrapolated to twice their size, is as far as every sample can be
#' extrapolated reliably.
#'
#' @section Confidence intervals:
#' With `n_boot > 0`, each estimate gets a bootstrap standard error and the
#' interval estimate \eqn{\pm z_{(1 + \text{conf\_level})/2}} SE (Chao et
#' al. 2014, Appendix G): `n_boot` times, a sample of the same size is drawn
#' from the estimated community (the observed features, with their
#' relative abundances adjusted for coverage, plus Chao1's undetected
#' ones), and every estimate is recomputed from it. At a fixed coverage,
#' each bootstrap sample's own depth for that coverage is used. Lower
#' limits stop at 0 (and the asymptotic richness's at the observed
#' richness), and coverage limits at 0 and 1, as in iNEXT.
#'
#' Replicate \eqn{b} of sample \eqn{i} has its own random stream, from
#' `seed`, so the results are the same for any `cpus`. The cost grows with
#' `n_boot + 1` and with the number of depths (`n_depths`): set `n_boot = 0`
#' for estimates alone.
#'
#' @section Differences from iNEXT:
#' * **Shannon diversity beyond the sample's size** uses Chao et al.'s
#'   (2014) Eq. 10c, the formula of the iNEXT paper. iNEXT 3 uses a newer
#'   extrapolation for \eqn{q = 1} that these papers don't give; the two
#'   agree at the sample's size and at the asymptote, and differ in
#'   between (by under 1% at twice the size, for iNEXT's spider data).
#'   Everything else equals iNEXT 3.0.2's estimates to rounding.
#' * **The bootstrap** draws from ecodive's own random streams, so standard
#'   errors and intervals agree with iNEXT's only within Monte Carlo error.
#'   The asymptotes' intervals use the same `n_boot` replicates as the
#'   curves; iNEXT's always use 100 (whatever `nboot`).
#' * **Coverage-based depths** are found exactly; iNEXT finds them by a
#'   numerical search, which leaves its depths (and its estimates there)
#'   off by about \eqn{10^{-5}}. Where coverage is flat over a range of
#'   depths, ecodive takes the smallest, and gives one row; iNEXT's search
#'   lands anywhere in the range, and repeats the row.
#' * With `depths`, each sample's own size is added (iNEXT also adds
#'   \eqn{n - 1} and \eqn{n + 1}). With `n_depths` below 4, iNEXT's grid
#'   differs.
#' * Abundance data only, \eqn{q \in \{0, 1, 2\}}, and the results' shapes
#'   (below). `hill_standardize()` holds coverage equal by default;
#'   `estimateD()` holds depth equal (`base = "size"`).
#'
#' @section Same names in other packages:
#' None. The same estimates: `iNEXT::iNEXT()`, `iNEXT::estimateD()`,
#' `iNEXT::DataInfo()`, `iNEXT::ChaoRichness()`, `iNEXT::ChaoShannon()` and
#' `iNEXT::ChaoSimpson()`.
#'
#' @inheritParams documentation
#' @inheritParams rarefaction_curve
#'
#' @param counts  A count table (samples \eqn{\times} features) of whole
#'        numbers; see **Input Types**. Each sample is one assemblage.
#'
#' @param q  The orders of the Hill numbers: any of `0` (richness), `1`
#'        (Shannon diversity) and `2` (Simpson diversity). Default:
#'        `c(0, 1, 2)`. Shannon diversity costs the most.
#'
#' @param depths  The depths to compute the curves at, for every sample:
#'        whole numbers of at least 1. Each sample's own size is added.
#'        Default: `NULL`, `n_depths` depths from 1 to `max_depth` (see
#'        "The depths").
#'
#' @param max_depth  The largest depth: `NULL` (the default) for twice each
#'        sample's size, or a whole number for the same for all samples.
#'        With `depths`, larger depths are dropped.
#'
#' @param n_depths  The number of depths in the default grid, a whole number
#'        of at least 2. Default: `40`.
#'
#' @param n_boot  The number of bootstrap replicates for the standard errors
#'        and intervals, a whole number; `0` for none. Default: `50`.
#'
#' @param conf_level  The confidence level of the intervals, between 0 and 1.
#'        Default: `0.95`.
#'
#' @param seed  An integer seed for the bootstrap. Default: `0`.
#'
#' @return
#' `hill_extrapolation()`: a list of class `ecodive_hill_extrapolation`, of
#' four data frames, with attributes `conf_level` and `n_boot`:
#' * `size`: the curves by depth (iNEXT's `size_based`), with columns
#'   `sample`, `q`, `depth`, `method` (`"rarefaction"`, `"observed"` or
#'   `"extrapolation"`), `estimate`, `lower`, `upper`, `coverage`,
#'   `coverage_lower` and `coverage_upper`. Ordered by sample, `q` and
#'   depth.
#' * `coverage`: the same curves by coverage (iNEXT's `coverage_based`):
#'   columns `sample`, `q`, `coverage`, `depth`, `method`, `estimate`,
#'   `lower` and `upper`. Each row has a `size` row's coverage, at the
#'   smallest depth that reaches it (coverage can be flat over a range of
#'   depths, which then gives one row); its interval is at that fixed
#'   coverage.
#' * `asymptotic`: columns `sample`, `q`, `observed`, `estimate`, `se`,
#'   `lower` and `upper` (iNEXT's `AsyEst`).
#' * `info`: columns `sample`, `n` (the size), `s_obs` (the richness),
#'   `coverage` (at the sample's size), and `f1` to `f10`, the numbers of
#'   features seen once to ten times (iNEXT's `DataInfo`).
#'
#' Without bootstrap (`n_boot` of 0 or 1), the interval and `se` columns are
#' `NA`. An empty sample has no rows in `size` and `coverage`, and `NA`s in
#' `asymptotic` and `info` (with a warning).
#'
#' @section The depths:
#' With `depths = NULL`, each sample of size \eqn{n} with endpoint \eqn{M}
#' (`max_depth`, else \eqn{2n}) gets iNEXT's grid of \eqn{K} = `n_depths`
#' knots: if \eqn{M > n}, \eqn{\lfloor K/2 \rfloor - 1} depths evenly from 1
#' to \eqn{n - 1}, then \eqn{n}, then \eqn{\lfloor K/2 \rfloor} from \eqn{n +
#' 1} to \eqn{M}; else \eqn{K} from 1 to \eqn{M}. Each is rounded down, and
#' duplicates dropped.
#'
#' @references
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
#' R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
#' numbers: a framework for sampling and estimation in species diversity
#' studies. *Ecological Monographs*, 84, 45--67. \doi{10.1890/13-0133.1}
#'
#' Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. *Ecology*, 93,
#' 2533--2547. \doi{10.1890/11-1952.1}
#'
#' Chao, A., Wang, Y. T., & Jost, L. (2013). Entropy and the species
#' accumulation curve: a novel entropy estimator via discovery rates of new
#' species. *Methods in Ecology and Evolution*, 4, 1091--1100.
#' \doi{10.1111/2041-210X.12108}
#'
#' Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for
#' rarefaction and extrapolation of species diversity (Hill numbers).
#' *Methods in Ecology and Evolution*, 7, 1451--1456.
#' \doi{10.1111/2041-210X.12613}
#'
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a
#' population. *Scandinavian Journal of Statistics*, 11, 265--270.
#'
#' @seealso [plot.ecodive_hill_extrapolation()], [rarefaction_curve()],
#'   [hill()], [rarefy()]
#'
#' @export
#' @examples
#'     # Richness, Shannon and Simpson diversity of each sample, by depth.
#'     he <- hill_extrapolation(ex_counts, n_boot = 10)
#'     he
#'     head(he$asymptotic)
#'     plot(he, q = 1)
#'     plot(he, type = 'coverage')
#'
#'     # Every sample at an equal coverage, then at an equal depth.
#'     hill_standardize(ex_counts, n_boot = 10)
#'     hill_standardize(ex_counts, depth = 1000, n_boot = 0)
#'
hill_extrapolation <- function (
    counts,
    q          = c(0, 1, 2),
    depths     = NULL,
    max_depth  = NULL,
    n_depths   = 40L,
    n_boot     = 50L,
    conf_level = 0.95,
    margin     = 1L,
    seed       = 0L,
    cpus       = n_cpus() ) {

  counts     <- rarefied_counts(environment())
  q          <- hx_check_q(q)
  depths     <- rc_check_depths(depths)
  if (!is.null(max_depth))
    max_depth <- rc_check_whole(max_depth, 'max_depth', 1, 'NULL or a whole number of at least 1')
  n_depths   <- rc_check_whole(n_depths, 'n_depths', 2, 'a whole number of at least 2')
  n_boot     <- hx_check_n_boot(n_boot)
  conf_level <- hx_check_conf(conf_level)

  info   <- count_samples(counts, margin)
  n      <- info$n
  labels <- if (is.null(info$labels)) as.character(seq_len(n)) else info$labels
  totals <- sample_sums(counts, margin)
  hx_warn_empty(totals, labels)

  per_sam <- lapply(totals, function (N) hx_grid(N, depths, max_depth, n_depths))

  res <- hx_run(
    counts, margin, size = per_sam,
    cov_level = lapply(per_sam, function (d) rep(NA_real_, length(d))),
    cov_given = per_sam,
    q = q, n_boot = n_boot, seed = seed, cpus = cpus )

  z   <- stats::qnorm((1 + conf_level) / 2)
  sam <- rep(seq_len(n), lengths(per_sam))
  dep <- unlist(per_sam, use.names = FALSE)
  if (is.null(dep)) dep <- numeric(0)
  N   <- totals[sam]
  met <- hx_method(dep, N)
  qi  <- match(q, 0:2)


  # By depth: one block of rows per q, then ordered by sample, q, depth.
  size <- do.call(rbind, lapply(seq_along(q), function (j) data.frame(
    sample         = labels[sam],
    q              = rep(q[[j]], length(dep)),
    depth          = dep,
    method         = met,
    estimate       = res$size[, qi[[j]]],
    lower          = hx_lower(res$size[, qi[[j]]], res$size_sd[, qi[[j]]], z),
    upper          = res$size[, qi[[j]]] + z * res$size_sd[, qi[[j]]],
    coverage       = res$size[, 4],
    coverage_lower = pmax(0, res$size[, 4] - z * res$size_sd[, 4]),
    coverage_upper = pmin(1, res$size[, 4] + z * res$size_sd[, 4]),
    stringsAsFactors = FALSE )))
  size <- hx_order(size, sam, length(q), dep)


  # By coverage: the same rows, with intervals at each fixed coverage.
  coverage <- do.call(rbind, lapply(seq_along(q), function (j) data.frame(
    sample   = labels[sam],
    q        = rep(q[[j]], length(dep)),
    coverage = res$cov_level,
    depth    = res$cov_depth,
    method   = hx_method(res$cov_depth, N),
    estimate = res$cov[, qi[[j]]],
    lower    = hx_lower(res$cov[, qi[[j]]], res$cov_sd[, qi[[j]]], z),
    upper    = res$cov[, qi[[j]]] + z * res$cov_sd[, qi[[j]]],
    stringsAsFactors = FALSE )))
  coverage <- hx_order(coverage, sam, length(q), dep)
  # Where coverage is flat over a range of depths, one row (the smallest).
  coverage <- coverage[!duplicated(coverage[c('sample', 'q', 'coverage')]), , drop = FALSE]
  rownames(coverage) <- NULL


  # Asymptotes.
  obs <- hx_observed(res, per_sam, n)
  asymptotic <- do.call(rbind, lapply(seq_along(q), function (j) {
    est <- res$asy[, qi[[j]]]
    se  <- res$asy_sd[, qi[[j]]]
    lo  <- est - z * se
    lo  <- if (q[[j]] == 0) pmax(obs[, 1], lo) else pmax(0, lo)
    data.frame(
      sample   = labels,
      q        = rep(q[[j]], n),
      observed = obs[, qi[[j]]],
      estimate = est,
      se       = se,
      lower    = lo,
      upper    = est + z * se,
      stringsAsFactors = FALSE )
  }))
  asymptotic <- asymptotic[order(rep(seq_len(n), length(q))), , drop = FALSE]
  rownames(asymptotic) <- NULL


  out <- list(
    size       = size,
    coverage   = coverage,
    asymptotic = asymptotic,
    info       = hx_info(res$info, labels) )
  attr(out, 'conf_level') <- conf_level
  attr(out, 'n_boot')     <- n_boot
  class(out) <- 'ecodive_hill_extrapolation'
  out
}



#' @rdname hill_extrapolation
#'
#' @param coverage  The sample coverage to compare the samples at, a number
#'        in (0, 1]. At most one of `coverage` and `depth`; with neither,
#'        the smallest coverage among the samples extrapolated to twice
#'        their size (iNEXT's default for `base = "coverage"`). A coverage of
#'        1 gives the asymptotes.
#'
#' @param depth  The depth (number of observations) to compare the samples
#'        at instead, a whole number of at least 1.
#'
#' @return
#' `hill_standardize()`: a numeric matrix with one row per sample and one
#' column per order (`q0`, `q1`, `q2`), of class `ecodive_hill_standardize`
#' (it prints a preview; `unclass()` gives the plain matrix), with
#' attributes:
#' * `lower`, `upper`: matrices of the interval limits (only with
#'   `n_boot > 1`).
#' * `depth`, `coverage`: each sample's depth and coverage there (one of
#'   them is the common level; the depth can be fractional at a fixed
#'   coverage, and `Inf` at coverage 1).
#' * `method`: each sample's `"rarefaction"`, `"observed"` or
#'   `"extrapolation"`.
#' * `level` and `base` (`"coverage"` or `"depth"`): what was held equal.
#' * `conf_level`.
#'
#' @export
hill_standardize <- function (
    counts,
    q          = c(0, 1, 2),
    coverage   = NULL,
    depth      = NULL,
    n_boot     = 50L,
    conf_level = 0.95,
    margin     = 1L,
    seed       = 0L,
    cpus       = n_cpus() ) {

  counts     <- rarefied_counts(environment())
  q          <- hx_check_q(q)
  n_boot     <- hx_check_n_boot(n_boot)
  conf_level <- hx_check_conf(conf_level)
  if (!is.null(coverage) && !is.null(depth))
    stop('Give `coverage` or `depth`, not both.', call. = FALSE)
  if (!is.null(coverage) &&
      (!is.numeric(coverage) || length(coverage) != 1 || is.na(coverage) ||
       coverage <= 0 || coverage > 1))
    stop('`coverage` must be NULL or a single number in (0, 1].', call. = FALSE)
  if (!is.null(depth))
    depth <- rc_check_whole(depth, 'depth', 1, 'NULL or a whole number of at least 1')

  info   <- count_samples(counts, margin)
  n      <- info$n
  labels <- if (is.null(info$labels)) as.character(seq_len(n)) else info$labels
  totals <- sample_sums(counts, margin)
  hx_warn_empty(totals, labels)
  some   <- totals >= 1

  base <- if (is.null(depth)) 'coverage' else 'depth'
  if (base == 'coverage' && is.null(coverage)) {
    # The smallest coverage of the samples extrapolated to twice their size.
    two <- hx_run(
      counts, margin, size = lapply(totals, function (N) if (N >= 1) 2 * N else numeric(0)),
      cov_level = lapply(totals, function (N) numeric(0)),
      cov_given = lapply(totals, function (N) numeric(0)),
      q = numeric(0), n_boot = 0L, seed = seed, cpus = cpus )
    coverage <- if (any(some)) min(two$size[, 4]) else 1
  }
  level <- if (base == 'coverage') as.double(coverage) else as.double(depth)

  empty <- lapply(totals, function (N) numeric(0))
  if (base == 'depth') {
    res <- hx_run(
      counts, margin, size = lapply(totals, function (N) if (N >= 1) level else numeric(0)),
      cov_level = empty, cov_given = empty,
      q = q, n_boot = n_boot, seed = seed, cpus = cpus )
    est <- res$size[, 1:3, drop = FALSE]
    sd  <- res$size_sd[, 1:3, drop = FALSE]
    dep <- rep(level, sum(some))
    cvg <- res$size[, 4]
  }
  else {
    res <- hx_run(
      counts, margin, size = empty,
      cov_level = lapply(totals, function (N) if (N >= 1) level else numeric(0)),
      cov_given = lapply(totals, function (N) if (N >= 1) NA_real_ else numeric(0)),
      q = q, n_boot = n_boot, seed = seed, cpus = cpus )
    est <- res$cov
    sd  <- res$cov_sd
    dep <- res$cov_depth
    cvg <- res$cov_reach
  }

  qi  <- match(q, 0:2)
  col <- paste0('q', q)
  fill <- function (v) {
    m <- matrix(NA_real_, n, length(q), dimnames = list(labels, col))
    m[some, ] <- v[, qi, drop = FALSE]
    m
  }
  out <- fill(est)

  z  <- stats::qnorm((1 + conf_level) / 2)
  if (n_boot > 1) {
    lo <- est - z * sd
    lo[] <- pmax(0, lo)
    attr(out, 'lower') <- fill(lo)
    attr(out, 'upper') <- fill(est + z * sd)
  }
  dv <- cv <- rep(NA_real_, n)
  dv[some] <- dep
  cv[some] <- cvg
  names(dv) <- names(cv) <- labels
  mt <- rep(NA_character_, n)
  mt[some] <- hx_method(dep, totals[some])
  names(mt) <- labels
  attr(out, 'depth')      <- dv
  attr(out, 'coverage')   <- cv
  attr(out, 'method')     <- mt
  attr(out, 'level')      <- level
  attr(out, 'base')       <- base
  attr(out, 'conf_level') <- conf_level
  class(out) <- c('ecodive_hill_standardize', 'matrix', 'array')
  out
}




#------------------------------------------------------------------------------
# Internals.
#------------------------------------------------------------------------------

# Runs C_hill_ext(). `size`, `cov_level`, `cov_given`: lists with one
# numeric vector per sample (depths; coverage levels, NA for the coverage at
# the given depth; the depth of the point estimate, NA to solve). `boot`:
# NULL, or a list of one numeric matrix per sample of bootstrap counts (for
# tests). Returns the estimates and their bootstrap standard deviations,
# one row per size or coverage row.
hx_run <- function (counts, margin, size, cov_level, cov_given, q, n_boot,
                    seed, cpus, boot = NULL) {

  size_pos <- as.integer(c(0L, cumsum(lengths(size))))
  cov_pos  <- as.integer(c(0L, cumsum(lengths(cov_level))))
  sd_      <- as.double(unlist(size, use.names = FALSE))
  cl       <- as.double(unlist(cov_level, use.names = FALSE))
  cg       <- as.double(unlist(cov_given, use.names = FALSE))
  need     <- as.integer(0:2 %in% q)
  if (!is.null(boot)) boot <- lapply(boot, function (b) { storage.mode(b) <- 'double'; b })

  r <- .Call(
    C_hill_ext, counts, margin, size_pos, sd_, cov_pos, cl, cg, need,
    as.integer(n_boot), as.integer(seed), boot, as.integer(cpus) )

  R  <- length(sd_)
  C  <- length(cl)
  n  <- length(size)
  out <- list(
    size      = t(r$size),
    cov_depth = r$cov_depth,
    cov       = t(r$cov)[, 1:3, drop = FALSE],
    cov_reach = t(r$cov)[, 4],
    cov_level = r$cov_level,
    asy       = t(r$asy),
    info      = t(r$info) )
  dim(out$size) <- c(R, 4)
  dim(out$cov)  <- c(C, 3)
  dim(out$asy)  <- c(n, 3)

  # Each sample's replicate values: 4 per size row, 3 per coverage row,
  # then 3 asymptotes.
  size_sd <- matrix(NA_real_, R, 4)
  cov_sd  <- matrix(NA_real_, C, 3)
  asy_sd  <- matrix(NA_real_, n, 3)
  if (n_boot > 0) {
    vp <- r$val_pos
    for (s in seq_len(n)) {
      v  <- r$sd[(vp[[s]] + 1L):vp[[s + 1L]]]
      nr <- size_pos[[s + 1L]] - size_pos[[s]]
      nc <- cov_pos[[s + 1L]]  - cov_pos[[s]]
      if (nr) size_sd[size_pos[[s]] + seq_len(nr), ] <- matrix(v[seq_len(4 * nr)], nr, 4, byrow = TRUE)
      if (nc) cov_sd[cov_pos[[s]] + seq_len(nc), ]   <- matrix(v[4 * nr + seq_len(3 * nc)], nc, 3, byrow = TRUE)
      asy_sd[s, ] <- v[4 * nr + 3 * nc + 1:3]
    }
  }
  out$size_sd <- size_sd
  out$cov_sd  <- cov_sd
  out$asy_sd  <- asy_sd
  out
}


# Each sample's depths: iNEXT's knots (see "The depths"), or `depths` with
# the sample's size added, at most `max_depth`.
hx_grid <- function (N, depths, max_depth, n_depths) {
  if (N < 1) return (numeric(0))
  if (!is.null(depths)) {
    d <- sort(unique(c(depths, N)))
    if (!is.null(max_depth)) d <- d[d <= max(max_depth, N)]
    return (as.double(d))
  }
  M <- if (is.null(max_depth)) 2 * N else max_depth
  K <- n_depths
  if (M <= N) {
    d <- floor(seq(1, M, length.out = K))
  }
  else {
    h  <- K %/% 2
    lo <- if (N > 1 && h > 1) floor(seq(1, N - 1, length.out = h - 1)) else numeric(0)
    d  <- c(lo, N, floor(seq(N + 1, M, length.out = h)))
  }
  as.double(unique(d))
}


hx_method <- function (depth, N) {
  ifelse(depth < N, 'rarefaction', ifelse(depth == N, 'observed', 'extrapolation'))
}


# Lower limits, which stop at 0.
hx_lower <- function (est, se, z) pmax(0, est - z * se)


# Rows in `df` (one block of rows per q, each in `sam`, `dep` order), ordered
# by sample, q, then depth.
hx_order <- function (df, sam, k, dep) {
  R   <- length(sam)
  ord <- order(rep(sam, k), rep(seq_len(k), each = R), rep(dep, k))
  df  <- df[ord, , drop = FALSE]
  rownames(df) <- NULL
  df
}


# The observed Hill numbers of each sample (from its row at its own size).
hx_observed <- function (res, per_sam, n) {
  obs <- matrix(NA_real_, n, 3)
  pos <- c(0L, cumsum(lengths(per_sam)))
  for (s in seq_len(n)) {
    d <- per_sam[[s]]
    i <- which(d == res$info[s, 1])
    if (length(i)) obs[s, ] <- res$size[pos[[s]] + i[[1]], 1:3]
  }
  obs
}


hx_info <- function (info, labels) {
  df <- data.frame(
    sample   = labels,
    n        = info[, 1],
    s_obs    = info[, 2],
    coverage = info[, 3],
    stringsAsFactors = FALSE )
  for (k in 1:10) df[[paste0('f', k)]] <- info[, 3 + k]
  df
}


hx_warn_empty <- function (totals, labels) {
  e <- which(totals < 1)
  if (length(e))
    warning(
      sprintf('%d empty sample%s (%s): no estimates.', length(e),
              if (length(e) == 1) '' else 's',
              paste(utils::head(labels[e], 5), collapse = ', ')),
      call. = FALSE )
}


hx_check_q <- function (q) {
  if (!is.numeric(q) || length(q) == 0 || anyNA(q) || !all(q %in% c(0, 1, 2)))
    stop('`q` must be one or more of 0, 1 and 2.', call. = FALSE)
  unique(as.double(q))
}


hx_check_n_boot <- function (n_boot) {
  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) ||
      n_boot < 0 || n_boot %% 1 != 0 || n_boot > .Machine$integer.max)
    stop('`n_boot` must be a whole number of at least 0.', call. = FALSE)
  as.integer(n_boot)
}


hx_check_conf <- function (conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1 || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop('`conf_level` must be a single number between 0 and 1.', call. = FALSE)
  as.double(conf_level)
}




#------------------------------------------------------------------------------
# Printing.
#------------------------------------------------------------------------------

#' @export
print.ecodive_hill_extrapolation <- function (x, ...) {
  info <- x$info
  ns   <- nrow(info)
  q    <- unique(x$asymptotic$q)
  nb   <- attr(x, 'n_boot', exact = TRUE)
  cl   <- attr(x, 'conf_level', exact = TRUE)
  cat('Interpolated and extrapolated Hill numbers\n')
  cat(sprintf(
    '%d sample%s; q = %s; %d depth%s by size; %s\n\n',
    ns, if (ns == 1) '' else 's', paste(q, collapse = ', '),
    nrow(x$size), if (nrow(x$size) == 1) '' else 's',
    if (isTRUE(nb > 1)) sprintf('%g%% intervals from %d bootstrap replicates', 100 * cl, nb)
    else 'no intervals' ))
  a <- x$asymptotic
  k <- min(nrow(a), 3L * length(q))
  if (k) {
    cat('Asymptotic estimates:\n')
    df <- a[seq_len(k), c('sample', 'q', 'observed', 'estimate', 'lower', 'upper')]
    print(format(df, digits = 4), row.names = FALSE)
    if (nrow(a) > k) cat(sprintf('# ... %d more rows\n', nrow(a) - k))
  }
  cat('\nTables: $size, $coverage, $asymptotic, $info\n')
  invisible(x)
}


#' @rdname print_ecodive
#' @export
print.ecodive_hill_standardize <- function (x, n = NULL, width = NULL, digits = 3L, ...) {
  att <- attributes(x)
  sym <- print_sym()
  lv  <- att$level
  head <- c(
    sprintf('A matrix: %d sample%s %s %d order%s', nrow(x), if (nrow(x) == 1) '' else 's',
            sym$times, ncol(x), if (ncol(x) == 1) '' else 's'),
    if (identical(att$base, 'coverage'))
      sprintf('Hill numbers at coverage %s', format(signif(lv, 4)))
    else sprintf('Hill numbers at depth %s', format(lv, scientific = FALSE)) )
  m <- unclass(x)
  attributes(m) <- att[c('dim', 'dimnames')]
  df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
  cat(
    print_table_lines(
      df, head, n = n, width = width, digits = digits, types = FALSE,
      rows = 'samples', attrs = print_attrs(x, c('dim', 'dimnames', 'class')) ),
    sep = '\n' )
  invisible(x)
}




#------------------------------------------------------------------------------
# Plotting.
#------------------------------------------------------------------------------

#' Plot interpolated and extrapolated Hill numbers
#'
#' Draws the curves of [hill_extrapolation()]: one per sample, rarefaction
#' as a solid line, extrapolation dashed, each sample's own size as a point,
#' with interval bands when `x` has them.
#'
#' @param x  The result of [hill_extrapolation()].
#'
#' @param type  What to draw: `"depth"` (the default; `"size"` too), the
#'        Hill number against depth (iNEXT's type 1); `"completeness"`,
#'        sample coverage against depth (type 2); or `"coverage"`, the Hill
#'        number against coverage (type 3). A prefix is enough.
#'
#' @param q  Which order to draw, one of `x`'s. Default: the first.
#'
#' @param samples  Which samples to draw, by name or position. Default: all.
#'
#' @param legend  Where to put the legend of samples (drawn for up to 10): a
#'        keyword for [graphics::legend()], or `NULL` for none. Default:
#'        `"bottomright"`.
#'
#' @param col  The samples' colours, one per sample drawn. Default:
#'        [grDevices::hcl.colors()]'s `"Dark 3"` palette.
#'
#' @param ...  Passed on to [graphics::plot.default()], such as `main`,
#'        `xlab`, `ylab`, `xlim`, `ylim` or `log`.
#'
#' @details
#' The titles are "Sample-size-based Rarefaction and Extrapolation", "Sample
#' Completeness Curves" and "Coverage-based Rarefaction and Extrapolation";
#' the y axis is "Species richness" (\eqn{q = 0}), "Shannon diversity"
#' (\eqn{q = 1}) or "Simpson diversity" (\eqn{q = 2}), as in Chao et al.
#' (2014). The subtitle names the bands' confidence level.
#'
#' @return `x`, invisibly.
#'
#' @references
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
#' R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
#' numbers: a framework for sampling and estimation in species diversity
#' studies. *Ecological Monographs*, 84, 45--67. \doi{10.1890/13-0133.1}
#'
#' @export
#' @examples
#'     he <- hill_extrapolation(ex_counts, n_boot = 10)
#'     plot(he)
#'     plot(he, q = 2, type = 'coverage')
#'     plot(he, type = 'completeness', samples = c('Saliva', 'Stool'))
plot.ecodive_hill_extrapolation <- function (
    x, type = 'depth', q = NULL, samples = NULL, legend = 'bottomright',
    col = NULL, ... ) {

  pos  <- plot_legend_pos(legend)
  type <- hx_plot_type(type)

  have <- unique(x$size$q)
  if (length(have) == 0)
    stop('`x` has no rows: there is nothing to draw.', call. = FALSE)
  if (is.null(q)) q <- have[[1]]
  if (!is.numeric(q) || length(q) != 1 || !q %in% have)
    stop('`q` must be one of ', paste(have, collapse = ', '), ', the orders in `x`.', call. = FALSE)

  labs <- x$info$sample
  if (is.null(samples)) {
    keep <- labs
  }
  else {
    keep <- if (is.numeric(samples)) labs[samples] else as.character(samples)
    if (anyNA(keep) || !all(keep %in% labs))
      stop('`samples` must name samples of `x`, or give their positions.', call. = FALSE)
    keep <- unique(keep)
  }
  keep <- keep[keep %in% x$size$sample]
  if (length(keep) == 0)
    stop('None of the samples has any rows to draw.', call. = FALSE)
  k    <- length(keep)
  cols <- if (is.null(col)) plot_palette(k) else rep_len(col, k)

  tab <- if (type == 'coverage') x$coverage else x$size
  tab <- tab[tab$q == q & tab$sample %in% keep, , drop = FALSE]
  tab <- tab[is.finite(tab$depth), , drop = FALSE]

  if (type == 'completeness') {
    xv <- tab$depth; yv <- tab$coverage; lo <- tab$coverage_lower; hi <- tab$coverage_upper
  }
  else if (type == 'coverage') {
    xv <- tab$coverage; yv <- tab$estimate; lo <- tab$lower; hi <- tab$upper
  }
  else {
    xv <- tab$depth; yv <- tab$estimate; lo <- tab$lower; hi <- tab$upper
  }
  bands <- any(is.finite(lo)) && any(is.finite(hi))

  lab  <- hx_plot_labels(type, q, attr(x, 'conf_level', exact = TRUE), bands)
  args <- list(...)
  yr   <- range(c(yv, if (bands) c(lo, hi)), finite = TRUE)
  plot_do(
    graphics::plot.default, c(list(NA, NA), args),
    list(
      type = 'n',
      xlim = range(xv, finite = TRUE),
      ylim = yr,
      main = lab$main, sub = lab$sub, xlab = lab$xlab, ylab = lab$ylab ))

  lwd <- plot_arg(args, 'lwd', 2)
  for (i in seq_len(k)) {
    j <- which(tab$sample == keep[[i]])
    j <- j[order(xv[j])]
    if (length(j) == 0) next
    if (bands) {
      b <- j[is.finite(lo[j]) & is.finite(hi[j])]
      if (length(b) > 1)
        graphics::polygon(
          c(xv[b], rev(xv[b])), c(lo[b], rev(hi[b])),
          col = grDevices::adjustcolor(cols[[i]], alpha.f = 0.2), border = NA )
    }
    r <- j[tab$method[j] != 'extrapolation']
    e <- j[tab$method[j] != 'rarefaction']
    if (length(r) > 1) graphics::lines(xv[r], yv[r], col = cols[[i]], lwd = lwd, lty = 1)
    if (length(e) > 1) graphics::lines(xv[e], yv[e], col = cols[[i]], lwd = lwd, lty = 2)
    o <- j[tab$method[j] == 'observed']
    if (length(o)) graphics::points(xv[o], yv[o], col = cols[[i]], pch = 19)
  }

  if (k <= 10) plot_legend(pos, keep, cols, NA, lty = 1)

  invisible(x)
}


# The default title, subtitle and axis labels of plot(x, type, q).
hx_plot_labels <- function (type, q, conf, bands) {
  list(
    main = switch(
      type,
      depth        = 'Sample-size-based Rarefaction and Extrapolation',
      completeness = 'Sample Completeness Curves',
      coverage     = 'Coverage-based Rarefaction and Extrapolation' ),
    sub  = if (bands) sprintf('Bands: %g%% confidence intervals', 100 * conf),
    xlab = if (type == 'coverage') 'Sample Coverage' else 'Depth',
    ylab = if (type == 'completeness') 'Sample Coverage'
           else c('Species richness', 'Shannon diversity', 'Simpson diversity')[q + 1] )
}


# `type`, matched (exactly, else by a unique prefix, ignoring case); "size"
# is "depth".
hx_plot_type <- function (type) {
  choices <- c('depth', 'size', 'completeness', 'coverage')
  if (is.character(type) && length(type) == 1 && !is.na(type) && nzchar(type)) {
    needle <- tolower(type)
    hit    <- if (needle %in% choices) needle else choices[startsWith(choices, needle)]
    if (length(hit) == 1) return (if (hit == 'size') 'depth' else hit)
  }
  stop('`type` must be one of "depth" (or "size"), "completeness", "coverage".', call. = FALSE)
}
