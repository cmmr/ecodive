# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# SparCC (Friedman & Alm 2012), written from the paper alone; see
# src/sparcc.c and plans/sparcc.md.



#' SparCC Correlation
#'
#' Estimates the correlations between the log **absolute** abundances of
#' features from their relative abundances alone, with SparCC (Friedman &
#' Alm 2012): the "co-occurrence network" correlations most often reported
#' for microbiome data. They can be estimated only approximately, by
#' assuming that most pairs of features are uncorrelated and that there are
#' many features.
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of raw counts (samples \eqn{\times}
#'        features), or any of the objects listed under **Input Types**.
#'        Correlations are taken between the **features** (the columns when
#'        `margin = 1`), and need at least 4 features and 3 samples. Use the
#'        counts as sequenced, not proportions, rarefied or zero-replaced
#'        values: the uncertainty that SparCC allows for depends on their
#'        size (see **Details**). Values that are not whole numbers are used
#'        as they are, with a warning.
#'
#' @param n_draws  How many draws of the samples' fractions to take from
#'        their Dirichlet posterior. The result is each pair's median over
#'        the draws. Default `100L`, as in the paper.
#'
#' @param max_exclusions  At most how many strongly correlated pairs to
#'        exclude, one at a time, when estimating the features' variances.
#'        `0L` excludes none (the "basic" SparCC). Default `20L`, as in the
#'        paper.
#'
#' @param exclusion_threshold  A pair is strongly correlated, and may be
#'        excluded, when the absolute value of its correlation is above
#'        this. Default `0.1`, as in the paper.
#'
#' @param permutations  How many null datasets to make for p-values: `0L`
#'        (the default) for none. Each null dataset resamples every
#'        feature's counts across the samples on its own, and is run
#'        through all of SparCC with the same settings, so this multiplies
#'        the run time by `permutations + 1`. See **P-values**.
#'
#' @param seed  An integer seed for the random draws and resampling. The same
#'        `seed` gives the same results, for any number of `cpus`.
#'
#' @param clust  Whether and how to cluster the features (rows and columns)
#'        of the result, as for [vlr()] (see [beta_div()]'s `clust`). The
#'        correlations are clustered as the dissimilarities \eqn{1 - r}.
#'
#' @return A `dist` object **over the features** (the packed lower
#'   triangle, as [rho_p()] returns), labelled with the feature names, of
#'   class `c("ecodive_dist", "dist")`, so that `plot()` draws it as a
#'   heatmap ([plot.ecodive_dist()]). Besides `dist`'s `Size`, `Labels`,
#'   `Diag` and `Upper`, its attributes are:
#'
#'   * `method`: `"sparcc"`.
#'   * `n_samples`: the number of samples.
#'   * `variances`: the estimated variance of each feature's log basis
#'     abundance (\eqn{\omega^2}), its median over the draws, named by
#'     feature and in the order of the result's labels.
#'   * `p_value`: with `permutations > 0`, one two-sided permutation
#'     p-value per pair, **parallel to the triangle** (same length and
#'     order: `p_value[k]` belongs to `x[k]`). Not adjusted for testing
#'     every pair: see [stats::p.adjust()].
#'   * `display`: its title, `"SparCC Correlation"` (see
#'     [ecodive_display]), and `clust`; and `order` when clustered.
#'
#'   The correlations are a **similarity** in a `dist` container:
#'   `as.matrix()` puts 0 on the diagonal, where they are 1, and they must
#'   not be given as `dm` to an ordination or test. Threshold them to find
#'   strongly associated pairs (`which(abs(r) > 0.3)`, or `abs(as.matrix(r))
#'   > 0.3` for an adjacency matrix), or use `1 - r` for a dissimilarity.
#'   With p-values, `r[p.adjust(attr(r, "p_value"), "BH") <= 0.05]` keeps
#'   the significant pairs.
#'
#' @details
#' **What is estimated.** Let \eqn{y_i} be the log of feature \eqn{i}'s
#' unobserved absolute abundance (the "basis"), with variance
#' \eqn{\omega_i^2}. The variance of a log ratio is the same for relative
#' and absolute abundances:
#' \deqn{t_{ij} = \mathrm{var}\left(\log\frac{x_i}{x_j}\right) =
#' \omega_i^2 + \omega_j^2 - 2\rho_{ij}\omega_i\omega_j}
#' (the variation matrix, as [vlr()]). SparCC solves this for the
#' correlations \eqn{\rho_{ij}} of the \eqn{y}. With \eqn{D} features there
#' are more unknowns than equations, so it assumes that most correlations
#' are about 0: then \eqn{t_i = \sum_{j \ne i} t_{ij} \approx (D - 2)
#' \omega_i^2 + \sum_j \omega_j^2}, which gives the variances, and
#' \deqn{\rho_{ij} = \frac{\omega_i^2 + \omega_j^2 - t_{ij}}{2 \omega_i
#' \omega_j}.}
#' The pair with the largest \eqn{|\rho|} above `exclusion_threshold` is
#' then left out of these sums, the variances are estimated again, and so
#' on, at most `max_exclusions` times.
#'
#' The approximation is better with more features, and when most pairs are
#' truly uncorrelated. As a rule of thumb, Friedman & Alm (2012) ask for an
#' effective number of features (a measure of how evenly abundant they
#' are) of at least about 50 before ordinary correlations of relative
#' abundances can be trusted, and show SparCC to be more robust to fewer.
#' Values beyond \eqn{\pm 1} are possible when too many pairs are excluded
#' for the variances to be estimated well, and are kept with a warning.
#'
#' **Zeros and counting noise.** The fractions of each sample are not
#' taken as observed: they are drawn `n_draws` times from their posterior
#' given the counts, \eqn{\mathrm{Dirichlet}(N_s + 1)} (a uniform prior,
#' so every count gets 1 added), and SparCC is run on each draw's log
#' fractions. The result is each pair's median over the draws. This
#' handles zeros without a pseudocount, and lets rare features, whose
#' fractions are uncertain, show weaker correlations. It is why raw counts
#' are needed.
#'
#' **Why not [rho_p()]?** Both start from the variation matrix \eqn{t_{ij}}
#' ([vlr()]), but they scale it differently and answer different questions.
#' \eqn{\rho_p = 1 - t_{ij} / (v_i + v_j)}, with \eqn{v_i} the variance of
#' the clr coordinate, needs no assumption, but is 1 only when the two
#' features keep a constant ratio (proportionality): it asks "do these two
#' keep a constant ratio?". SparCC replaces \eqn{v_i} with the estimated
#' basis variances, and asks "are their absolute abundances correlated?",
#' which can be answered only under its sparsity assumption. Pairs related
#' by a power law with an exponent other than 1, or with unequal
#' variances, have a high SparCC correlation but a low \eqn{\rho_p}.
#'
#' @section P-values:
#'   With `permutations = B`, each null dataset resamples the counts of
#'   every feature across the samples with replacement, independently of
#'   the other features, as FastSpar's `fastspar_bootstrap` does: all
#'   features are then uncorrelated, and each keeps the distribution of its
#'   counts. SparCC is run on each
#'   one with the same settings (and draws of its own), and a pair's
#'   p-value is the fraction of the \eqn{B + 1} datasets, the observed one
#'   included, whose correlation is at least as large in absolute value:
#'   \eqn{(1 + k) / (B + 1)} when \eqn{k} of the null correlations have
#'   \eqn{|r_b| \ge |r|}. It is never 0 (Phipson & Smyth 2010), as FastSpar
#'   computes it; the paper used the fraction of the null datasets alone.
#'
#'   This runs SparCC \eqn{B + 1} times, each with `n_draws` draws. When
#'   that is too slow, a smaller `n_draws` (such as `20L`) for the whole
#'   call is better than fewer permutations: the null datasets then have the
#'   same Monte Carlo noise as the observed one, so the test stays valid.
#'   The smallest possible p-value is \eqn{1 / (B + 1)}, so many
#'   permutations are needed to survive an adjustment for testing every
#'   pair (for 1,000 pairs and Bonferroni at 0.05, \eqn{B \ge 20000}).
#'
#' @section Other implementations:
#'   FastSpar's `--iterations`, `--exclude_iterations` and `--threshold`
#'   (Watts et al. 2019; defaults 50, 10 and 0.1), and SpiecEasi's
#'   `sparcc(iter, inner_iter, th)` (defaults 20, 10 and 0.1), are
#'   `n_draws`, `max_exclusions` and `exclusion_threshold`. With the same
#'   values, results agree with both to within the noise of the random
#'   draws, which no two implementations share (about 0.002 for counts in
#'   the thousands and 100 draws), with one exception: FastSpar makes at
#'   most \eqn{D - 3} exclusions, so with fewer than `max_exclusions + 3`
#'   features, give `max_exclusions = D - 3` to reproduce it. SpiecEasi has
#'   no such limit, like ecodive. Otherwise:
#'
#'   * ecodive draws with its own generator (PCG), so results are the same
#'     on every platform and for any number of `cpus`.
#'   * Each pair's result is the median of its values over the draws (the
#'     mean of the middle two for an even number of draws).
#'   * Within a draw, a feature whose estimated variance is not positive
#'     has no correlations; they are left out of the medians, with a
#'     warning giving how often it happened.
#'   * With `exclusion_threshold = 0`, pairs whose correlation is 0 but for
#'     rounding error can be excluded, until a feature has no pair left
#'     (it is then left out of that draw) or the pairs still included no
#'     longer determine the variances. The loop then stops, when fewer than
#'     4 features are left or before that last exclusion. With any
#'     positive threshold this can't happen: such pairs' correlations are
#'     exactly 0.
#'   * The p-values are FastSpar's: its resampling, and its estimate
#'     \eqn{(1 + k) / (B + 1)}. SpiecEasi's `sparccboot()` permutes each
#'     feature's counts instead.
#'
#' @section Same names in other packages:
#'   `SpiecEasi::sparcc()` computes the same estimate, with other defaults
#'   (see **Other implementations**).
#'
#' @section Size:
#'   The result has \eqn{D(D - 1)/2} values: 2,000 features give 2e6 pairs,
#'   16 MB. Every value of the table is drawn, so each thread works on a
#'   dense copy (\eqn{n \times D} values) and the \eqn{D(D - 1)/2} log-ratio
#'   variances of one draw at a time. The medians need every draw's value
#'   of a pair: up to 256 MB of them are kept at once, and beyond that the
#'   pairs are taken in blocks, each one drawing again (with the same
#'   numbers), which takes longer. Filter out rare features first for large
#'   tables.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from
#' genomic survey data. *PLoS Computational Biology*, 8, e1002687.
#' \doi{10.1371/journal.pcbi.1002687}
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology*, 9, 39.
#' \doi{10.2202/1544-6115.1585}
#'
#' Watts, S. C., Ritchie, S. C., Inouye, M., & Holt, K. E. (2019). FastSpar:
#' rapid and scalable correlation estimation for compositional data.
#' *Bioinformatics*, 35, 1064–1066. \doi{10.1093/bioinformatics/bty734}
#'
#' @seealso [rho_p()], proportionality of features; [vlr()], the variation
#'   matrix; [plot.ecodive_dist()].
#'
#' @export
#' @examples
#'     # 60 samples of 20 taxa, whose first two are correlated
#'     set.seed(1)
#'     basis <- matrix(rnorm(60 * 20), 60, 20)
#'     basis[, 2] <- basis[, 1] + rnorm(60, sd = 0.3)
#'     counts <- matrix(
#'       rpois(1200, exp(basis + 4)), 60, 20,
#'       dimnames = list(NULL, paste0('OTU', 1:20)) )
#'
#'     r <- sparcc(counts, n_draws = 20)
#'     r
#'
#'     # Each taxon's estimated log-abundance variance
#'     attr(r, 'variances')
#'
#'     # Pairs above a cutoff
#'     m <- as.matrix(r)
#'     which(abs(m) > 0.5 & upper.tri(m), arr.ind = TRUE)
#'
#'     # With p-values from 99 null datasets
#'     r <- sparcc(counts, n_draws = 10, permutations = 99)
#'     head(attr(r, 'p_value'))
#'
#'     plot(r)
#'
sparcc <- function (
    counts, margin = 1L, n_draws = 100L, max_exclusions = 20L,
    exclusion_threshold = 0.1, permutations = 0L, clust = NULL,
    seed = 0L, cpus = n_cpus() ) {

  sparcc_fit(environment())
}



# sparcc() with the details the tests need. `log_fractions`: a list of
# n x D matrices (samples x features), one per draw, of log abundances to
# use in place of the Dirichlet draws. `block_pairs`: pairs per stored
# block (0: the default for the memory). `details = FALSE` returns what
# sparcc() would.
sparcc_details <- function (
    counts, margin = 1L, n_draws = 100L, max_exclusions = 20L,
    exclusion_threshold = 0.1, permutations = 0L, seed = 0L,
    cpus = n_cpus(), log_fractions = NULL, block_pairs = 0, details = TRUE ) {

  clust <- NULL
  sparcc_fit(environment(), details, log_fractions, block_pairs)
}



sparcc_fit <- function (env, details = FALSE, log_fractions = NULL, block_pairs = 0) {

  validate_clust(env)
  info   <- sparcc_input(env)
  counts <- env$counts
  margin <- env$margin
  cpus   <- env$cpus
  validate_seed(env)
  seed   <- env$seed
  n      <- info$n
  D      <- info$D

  count_arg <- function (arg, min) {
    v <- env[[arg]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < min ||
        v %% 1 != 0 || v > .Machine$integer.max)
      stop(
        '`', arg, '` must be a single ',
        if (min > 0) 'positive' else 'non-negative', ' integer.', call. = FALSE )
    as.integer(v)
  }
  n_draws        <- count_arg('n_draws', 1)
  max_exclusions <- count_arg('max_exclusions', 0)
  permutations   <- count_arg('permutations', 0)

  thr <- env$exclusion_threshold
  if (!is.numeric(thr) || length(thr) != 1L || is.na(thr) || thr < 0)
    stop('`exclusion_threshold` must be a single non-negative number.', call. = FALSE)
  thr <- as.double(thr)

  given <- NULL
  if (!is.null(log_fractions)) {
    stopifnot(
      is.list(log_fractions), length(log_fractions) >= 1L,
      all(vapply(log_fractions, function (m) identical(as.integer(dim(m)), c(n, D)), logical(1))) )
    given   <- array(as.double(unlist(log_fractions)), c(n, D, length(log_fractions)))
    stopifnot(all(is.finite(given)))
    n_draws <- length(log_fractions)
  }

  res <- .Call(
    C_sparcc, counts, margin, n_draws, max_exclusions, thr, permutations,
    seed, given, as.double(block_pairs), isTRUE(details), cpus )

  r <- res$rho
  attributes(r) <- list(
    Size   = D,
    Labels = info$features,
    Diag   = FALSE,
    Upper  = FALSE,
    class  = c('ecodive_dist', 'dist') )
  attr(r, 'method')    <- 'sparcc'
  attr(r, 'n_samples') <- n
  attr(r, 'variances') <- stats::setNames(res$variances, info$features)
  if (permutations > 0L) attr(r, 'p_value') <- res$p_value
  r <- set_display(r, title = FEATURE_ASSOCS$name[match('sparcc', FEATURE_ASSOCS$id)])


  #--------------------------------------------------------
  # Warnings.
  #--------------------------------------------------------
  if (res$n_nonpos > 0L)
    warning(sprintf(
      paste0(
        'In %d of the %d draws, some features had an estimated variance ',
        'that was not positive (%d in all); their correlations in those ',
        'draws were left out of the medians.' ),
      res$n_draws_nonpos, n_draws, res$n_nonpos), call. = FALSE )

  n_na <- sum(is.na(r))
  if (n_na > 0L)
    warning(sprintf(
      '%d of the %d pairs have no correlation in any draw, and are NA.',
      n_na, length(r) ), call. = FALSE )

  n_big <- sum(abs(r) > 1, na.rm = TRUE)
  if (n_big > 0L)
    warning(sprintf(
      paste0(
        '%d of the %d correlations are beyond -1 or 1: the variances are ',
        'poorly estimated (too few features, or too many strongly ',
        'correlated pairs excluded).' ),
      n_big, length(r) ), call. = FALSE )

  r <- dist_clust(r, env$clust, cpus)
  ord <- attr(r, 'order', exact = TRUE)
  if (!is.null(ord)) {
    v <- attr(r, 'variances')[ord]
    names(v) <- attr(r, 'Labels')     # positions, for unnamed features
    attr(r, 'variances') <- v
  }

  if (!isTRUE(details)) return (r)

  list(
    rho            = r,
    draw_rho       = res$draw_rho,
    draw_variances = t(res$draw_variances),
    n_exclusions   = res$n_exclusions,
    exclusions     = res$exclusions,
    stop           = c('max_exclusions', 'threshold', 'three', 'singular')[res$stop],
    null_rho       = res$null_rho )
}



# Checks the count table of sparcc() in `env`: numeric, finite,
# non-negative, at least 3 samples and 4 features; warns for values that
# are not whole numbers.
sparcc_input <- function (env) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- env$counts

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
  if (any(vals %% 1 != 0))
    warning(
      '`counts` are not all whole numbers. SparCC expects raw counts: the ',
      'spread of its Dirichlet draws, and so its results, depend on their ',
      'size.', call. = FALSE )
  remove('vals')

  margin <- env$margin
  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  if (D < 4L)
    stop('`counts` must have at least 4 features for SparCC.', call. = FALSE)
  if (n < 3L)
    stop('`counts` must have at least 3 samples for SparCC.', call. = FALSE)

  list(
    D        = D,
    n        = n,
    features = if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]]) )
}
