# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Homogeneity of Multivariate Dispersions (PERMDISP)
#'
#' Tests whether groups of samples differ in how spread out they are: their
#' multivariate dispersion, or beta diversity. Gives the same results as
#' `vegan::betadisper()` followed by `vegan::permutest()`: each sample's
#' distance to its group's centre, each group's mean distance, and a
#' permutation test of whether those distances differ between groups. Note
#' that the default centre differs: see `type`.
#'
#' A significant PERMANOVA ([permanova()]) can come from groups that differ
#' in location, in dispersion, or both. A non-significant `permdisp()` result
#' supports reading a PERMANOVA difference as a difference in location.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values.
#'
#' @param groups  Each sample's group: a vector or factor with one value per
#'        sample, matched to the samples in `dm` by name when both are named,
#'        else in the same order as `dm`. Or, when `data` is given, the name of
#'        one of its columns (or a vector with one value per row of `data`).
#'        Values are treated as categories, even if they are numbers.
#'        Groups with no samples left are dropped.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `groups` (and `strata`) as columns. Rows are matched
#'        to the samples in `dm` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as `dm`.
#'        Default: `NULL`.
#'
#' @param type  The group centre to measure distances from:
#'
#'   * `'centroid'`: The centroid (multivariate mean), as in Anderson's
#'     (2006) original PERMDISP. The default: exact, and much faster for
#'     large data sets.
#'   * `'median'`: The spatial median, the point with the smallest sum of
#'     distances to the group's samples. More robust to skewed data and
#'     outlying samples. The default of `vegan::betadisper()` and
#'     scikit-bio's `permdisp`, so pass `type = 'median'` to compare with
#'     them.
#'
#' @param bias_adjust  Correct the distances for small group sizes, by
#'        multiplying each by \eqn{\sqrt{n_g / (n_g - 1)}}, for a group of
#'        \eqn{n_g} samples (Stier et al. 2013), as `betadisper(bias.adjust =
#'        TRUE)`. Default: `FALSE`.
#'
#' @param pairwise  Also test each pair of groups with a \eqn{t} statistic,
#'        as `permutest(pairwise = TRUE)`. Default: `FALSE`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A list of class `c("ecodive_permdisp", "ecodive_test")`, with
#'   these fields:
#'
#'   * `statistic`: The \eqn{F} statistic, named `F`.
#'   * `p_value`: The permutation p-value (`NA` without permutations).
#'   * `z`: How far \eqn{F} lies from its permuted values, in standard
#'     deviations: \eqn{(F - mean(F_{perm})) / sd(F_{perm})}, as
#'     `summary(vegan::permustats())`. `NA` with fewer than 2 permutations.
#'   * `n_perm`: The number of permutations.
#'   * `table`: The one-way ANOVA of the distances on the groups, as
#'     `permutest()$tab`: a data frame (of class `"anova"`) with rows
#'     `"Groups"` and `"Residuals"`, and columns `Df`, `Sum Sq`, `Mean Sq`,
#'     `F`, `N.Perm`, and `Pr(>F)`.
#'   * `distances`: Each sample's distance to its group's centre (named by
#'     sample, when `dm` has names), as `betadisper()$distances`.
#'   * `dispersions`: Each group's mean distance, as
#'     `betadisper()$group.distances`.
#'   * `type`: `'median'` or `'centroid'`.
#'   * `bias_adjust`: `TRUE` or `FALSE`.
#'   * `pairwise`: With `pairwise = TRUE`, a data frame with one row per pair
#'     of groups, named like `"a-b"`, and columns `t` (the \eqn{t} statistic
#'     of group `a` minus group `b`), `Df`, `p_value` (from the \eqn{t}
#'     distribution), and `p_value_perm` (from the permutations). Otherwise
#'     `NULL`.
#'   * `permutation`: How samples were permuted: `'free'`, `'within
#'     strata'`, or `'supplied matrix'`.
#'
#'   Attribute `f_perm` holds the \eqn{F} statistic of each permutation, as
#'   `permutest()$perm`, and, with `pairwise = TRUE`, attribute `t_perm` is a
#'   matrix of the \eqn{t} statistics of each permutation (rows) for each
#'   pair (columns).
#'
#' @details
#'
#' **Method**
#'
#' PERMDISP (Anderson 2006) places the samples in the principal coordinates
#' of `dm` ([pcoa()]), finds each group's centre, and measures each sample's
#' distance to it. For distances that are not Euclidean, such as
#' Bray-Curtis, some principal coordinates are imaginary (their eigenvalues
#' are negative); following Anderson (2006), a sample's squared distance to
#' its centre is the squared distance on the real axes minus that on the
#' imaginary ones, and is set to zero, with a warning, when that is negative.
#'
#' The test is the one-way ANOVA \eqn{F} of these distances. As in
#' `permutest()`, the p-value comes from permuting the ANOVA's residuals (each
#' distance minus its group's mean) among the samples: the proportion of
#' permutations, counting the observed one, with an \eqn{F} at least as
#' large, \eqn{(1 + \#\{F_{perm} \ge F\}) / (1 + permutations)}. Permuted
#' statistics within `sqrt(.Machine$double.eps)` of the observed one count as
#' ties. Pairwise \eqn{t} statistics permute the distances themselves, and
#' count \eqn{|t_{perm}| \ge |t|}.
#'
#'
#' **Choosing `type`**
#'
#' * `'centroid'`, the default, measures dispersion around the same centres
#'   that [permanova()] does: the squared distances to the centroids sum to
#'   PERMANOVA's within-group sum of squares. So it is the direct check of
#'   whether a significant PERMANOVA could come from differences in
#'   dispersion (Anderson 2006; Anderson et al. 2006). It is also exact, and
#'   needs only the distances within each group, so it is much faster for
#'   many samples (see below).
#' * `'median'` is the robust alternative for skewed or outlier-prone data,
#'   as the median-based (Brown-Forsythe) version of Levene's test is
#'   (Anderson 2006). It is the default of `vegan::betadisper()` and of
#'   scikit-bio's (and so QIIME 2's) `permdisp`: to compare with them, pass
#'   `type = 'median'`.
#'
#' `betadisper()` stops its search for the spatial median early, so its
#' distances to medians are accurate only to about \eqn{10^{-5}};
#' `permdisp()` solves for the median exactly.
#'
#'
#' **Speed and memory**
#'
#' A sample's distance to its group's centroid needs only the distances
#' within its group: \deqn{z_i^2 = \frac{1}{n_g}\sum_{j \in g} d_{ij}^2 -
#' \frac{1}{2n_g^2}\sum_{j,k \in g} d_{jk}^2,} which already accounts for
#' the imaginary axes. So `type = 'centroid'` needs no principal coordinates
#' at all. The spatial median does need them: a full eigendecomposition of
#' an \eqn{n \times n} matrix, which takes about as long as in
#' `betadisper()` (a few seconds for 1000 samples). Each permutation then
#' costs only one pass over the \eqn{n} distances, and permutations run in
#' parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. Results are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `permutest()`'s, so
#' p-values agree with it only up to sampling error, unless the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from betadisper() and permutest()**
#'
#' * `type = 'centroid'` gives the same distances and \eqn{F} as
#'   `betadisper()`, to rounding error.
#' * `type = 'median'` finds each spatial median to full precision.
#'   `betadisper()` stops its optimizer (`optim()`'s BFGS) early, so its
#'   distances to the median are off by up to about \eqn{10^{-5}} (relative),
#'   and its \eqn{F} agrees with `permdisp()`'s to about 4 or 5 significant
#'   digits.
#' * Samples with an `NA` group are an error, where `betadisper()` drops
#'   them.
#' * With `bias_adjust = TRUE`, a group of one sample keeps a distance of 0,
#'   where `betadisper()` gives `NaN`.
#' * Pairwise permutation p-values count permuted \eqn{|t|} within
#'   `sqrt(.Machine$double.eps)` of the observed one as ties, as for
#'   \eqn{F}.
#' * `strata` permutes freely within strata. For other restricted designs,
#'   pass a `permutations` matrix from `permute::shuffleSet()`.
#' * `betadisper()`'s `sqrt.dist` and `add` options transform the distances
#'   first: pass `sqrt(dm)` instead, for instance.
#'
#' @references
#' Anderson, M. J. (2006). Distance-based tests for homogeneity of
#' multivariate dispersions. *Biometrics*, 62(1), 245-253.
#' \doi{10.1111/j.1541-0420.2005.00440.x}
#'
#' Anderson, M. J., Ellingsen, K. E., & McArdle, B. H. (2006). Multivariate
#' dispersion as a measure of beta diversity. *Ecology Letters*, 9(6),
#' 683-693. \doi{10.1111/j.1461-0248.2006.00926.x}
#'
#' Stier, A. C., Geange, S. W., Hanson, K. M., & Bolker, B. M. (2013).
#' Predator density and timing of arrival affect reef fish community
#' assembly. *Ecology*, 94(5), 1057-1068. \doi{10.1890/11-1983.1}
#'
#' Vardi, Y., & Zhang, C.-H. (2000). The multivariate L1-median and
#' associated data depth. *Proceedings of the National Academy of Sciences*,
#' 97(4), 1423-1426. \doi{10.1073/pnas.97.4.1423}
#'
#' @seealso [permanova()], [beta_div()], [pcoa()]
#'
#' @export
#' @examples
#'     # 30 samples in three groups. The 'high' group is more variable.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       treatment = rep(c('ctrl', 'low', 'high'), each = 10) )
#'     high <- meta$treatment == 'high'
#'     counts[high, ] <- counts[high, ] * rpois(10 * 20, 1)
#'
#'     dm <- bray(counts)
#'
#'     # Do the treatment groups differ in dispersion?
#'     res <- permdisp(dm, 'treatment', meta)
#'     res
#'
#'     # Each group's mean distance to its centroid.
#'     res$dispersions
#'
#'     # Distances to spatial medians, as vegan::betadisper()'s default,
#'     # with each pair of groups tested.
#'     res <- permdisp(dm, meta$treatment, type = 'median', pairwise = TRUE)
#'     res$pairwise
#'
permdisp <- function (
    dm, groups, data = NULL, type = 'centroid', bias_adjust = FALSE,
    pairwise = FALSE, permutations = 999L, strata = NULL, seed = 0L,
    cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  type <- match.arg(type, c('centroid', 'median'))

  if (!isTRUE(bias_adjust) && !isFALSE(bias_adjust))
    stop('`bias_adjust` must be TRUE or FALSE.', call. = FALSE)
  if (!isTRUE(pairwise) && !isFALSE(pairwise))
    stop('`pairwise` must be TRUE or FALSE.', call. = FALSE)

  dv <- dist_vector(dm)
  n  <- dv$n

  if (n < 2)
    stop('`dm` must have at least 2 samples.', call. = FALSE)


  #--------------------------------------------------------
  # Each sample's group, aligned with the samples of `dm`.
  #--------------------------------------------------------
  if (missing(groups)) groups <- NULL

  des   <- test_group_values(groups, data, n, dv$labels, strata)
  grp_f <- des$groups
  if (!is.factor(grp_f)) grp_f <- factor(grp_f, exclude = NULL)
  grp_f <- droplevels(grp_f)
  ng    <- nlevels(grp_f)
  grp   <- as.integer(grp_f) - 1L
  cnt   <- tabulate(grp + 1L, ng)

  pm <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  #--------------------------------------------------------
  # Distances to the group centres.
  #--------------------------------------------------------
  res <- .Call(C_permdisp, dv$dm, grp, as.integer(c(ng, type == 'median', cpus)))
  z2  <- res$z2

  if (!res$converged)
    warning('A spatial median did not converge in 1000 iterations.', call. = FALSE)

  if (any(z2 < -1e-12 * max(abs(z2))))
    warning('Some squared distances to group centres are negative, and were set to 0.', call. = FALSE)

  z <- sqrt(pmax(z2, 0))

  if (bias_adjust) {
    ni <- cnt[grp + 1L]
    z  <- ifelse(ni > 1, z * sqrt(ni / (ni - 1)), z)
  }

  names(z) <- dv$labels

  disp <- vapply(seq_len(ng), function (g) mean(z[grp == g - 1L]), numeric(1))
  names(disp) <- levels(grp_f)


  #--------------------------------------------------------
  # One-way ANOVA of the distances.
  #--------------------------------------------------------
  mu     <- disp[grp + 1L]
  df_grp <- ng - 1L
  df_res <- n - ng
  ss_grp <- sum(cnt * (disp - mean(z))^2)
  ss_res <- sum((z - mu)^2)
  f_obs  <- (ss_grp / df_grp) / (ss_res / df_res)

  # Degenerate fits. Sums of squares within ~1e-12 of the total are rounding
  # noise: with none left within the groups, F is infinite (and cannot be
  # tested), whichever `type` made the distances.
  ss_tot <- ss_grp + ss_res
  if (ng < 2) {
    warning('`groups` has only one group, so there is nothing to test.', call. = FALSE)
    f_obs <- NA_real_
  }
  else if (df_res < 1) {
    warning('Every group has only one sample, so there is nothing to test.', call. = FALSE)
    f_obs <- NA_real_
  }
  else if (!isTRUE(ss_tot > 0)) {
    warning('Every sample is as far from its group centre as every other, so there is nothing to test.', call. = FALSE)
    f_obs <- NA_real_
  }
  else if (ss_res <= 1e-12 * ss_tot) {
    warning('The samples of each group are all equally far from its centre (no residual variation), so F is infinite and cannot be tested.', call. = FALSE)
    f_obs <- Inf
  }


  #--------------------------------------------------------
  # Permutations.
  #--------------------------------------------------------
  pr <- .Call(
    C_permdisp_perm, unname(z - mu), unname(z), grp, pm$perms,
    if (is.null(des$strata)) NULL else as.integer(des$strata),
    as.integer(c(ng, pm$n_perm, seed, cpus, pairwise)) )

  f_perm <- pr$f
  if (!is.finite(f_obs)) f_perm[] <- NA_real_
  p_obs  <- if (is.finite(f_obs)) perm_p(f_obs, f_perm) else NA_real_


  #--------------------------------------------------------
  # Pairwise t tests, as permutest(pairwise = TRUE).
  #--------------------------------------------------------
  pw <- t_perm <- NULL

  if (pairwise && ng >= 2) {

    combos <- utils::combn(ng, 2)
    labs   <- paste(levels(grp_f)[combos[1, ]], levels(grp_f)[combos[2, ]], sep = '-')

    t_obs <- apply(combos, 2, function (uv) {
      x <- z[grp == uv[[1]] - 1L]
      y <- z[grp == uv[[2]] - 1L]
      m <- length(x)
      k <- length(y)
      pooled <- sqrt(((m - 1) * stats::var(x) + (k - 1) * stats::var(y)) / (m + k - 2))
      # No variance within the two groups beyond rounding: t is undefined.
      if (!isTRUE(pooled > 1e-12 * max(abs(c(x, y))))) return (NA_real_)
      (mean(x) - mean(y)) / (pooled * sqrt(1 / m + 1 / k))
    })
    t_obs[is.nan(t_obs)] <- NA_real_
    t_df <- cnt[combos[1, ]] + cnt[combos[2, ]] - 2

    t_perm <- t(matrix(pr$t, nrow = ncol(combos)))
    t_perm[is.nan(t_perm)] <- NA_real_
    # As var(), a group of one sample has no variance.
    t_perm[, cnt[combos[1, ]] < 2 | cnt[combos[2, ]] < 2] <- NA_real_
    dimnames(t_perm) <- list(NULL, labs)

    p_t <- perm_p(t_obs, t_perm, 'two.sided')

    pw <- data.frame(
      t            = t_obs,
      Df           = as.double(t_df),
      p_value      = 2 * stats::pt(-abs(t_obs), t_df),
      p_value_perm = unname(p_t),
      row.names    = labs )
  }


  #--------------------------------------------------------
  # The table.
  #--------------------------------------------------------
  tbl <- data.frame(
    Df     = as.double(c(df_grp, df_res)),
    SS     = c(ss_grp, ss_res),
    MS     = c(ss_grp / df_grp, ss_res / df_res),
    F      = c(f_obs, NA),
    N.Perm = c(as.double(pm$n_perm), NA),
    P      = c(p_obs, NA),
    row.names   = c('Groups', 'Residuals'),
    check.names = FALSE )
  names(tbl) <- c('Df', 'Sum Sq', 'Mean Sq', 'F', 'N.Perm', 'Pr(>F)')

  res <- new_ecodive_test(
    class       = 'ecodive_permdisp',
    title       = 'Homogeneity of multivariate dispersions (PERMDISP)',
    statistic   = c(F = f_obs),
    p_value     = p_obs,
    z           = perm_z(f_obs, f_perm),
    n_perm      = pm$n_perm,
    table       = structure(tbl, class = c('anova', 'data.frame')),
    distances   = z,
    dispersions = disp,
    type        = type,
    bias_adjust = bias_adjust,
    pairwise    = pw,
    perm        = f_perm,
    permutation = if (!is.null(pm$perms)) 'supplied matrix'
                  else if (!is.null(des$strata)) 'within strata'
                  else 'free' )

  attr(res, 't_perm') <- t_perm
  res
}



# Lines added to print.ecodive_test()'s output.
#' @export
test_details.ecodive_permdisp <- function (x, digits) {

  out <- function (v) utils::capture.output(print(v, digits = digits))

  c(
    paste0(
      'Distances to group ',
      if (x$type == 'median') 'spatial medians' else 'centroids',
      if (x$bias_adjust) ', adjusted for small groups', ':' ),
    out(x$dispersions),
    '',
    out(x$table),
    if (!is.null(x$pairwise)) c('', 'Pairwise t tests:', out(x$pairwise)) )
}
