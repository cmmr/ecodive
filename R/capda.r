# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Discriminant Canonical Analysis of Principal Coordinates (CAP)
#'
#' Supervised ordination of samples from any distances: linear discriminant
#' analysis ([lda()]) of their first \eqn{m} principal coordinates
#' ([pcoa()]), with \eqn{m} chosen by leave-one-out classification success,
#' and a permutation test of the separation of the groups (Anderson &
#' Willis 2003; "CAP" in PRIMER-PERMANOVA+). Gives the same number of axes,
#' squared canonical correlations and classification success as
#' `BiodiversityR::CAPdiscrim()`. It is not `vegan::capscale()`, which is
#' distance-based redundancy analysis.
#'
#' @inheritParams anosim
#' @inheritParams documentation
#'
#' @param groups  Each sample's group, which the axes separate: a vector or
#'        factor with one value per sample, matched to the samples in `dm` by
#'        name when both are named, else in the same order as `dm`. Or, when
#'        `data` is given, the name of one of its columns (or a vector with
#'        one value per row of `data`). Values are treated as categories,
#'        even if they are numbers.
#'
#' @param n_pcoa_axes  The number of principal coordinates (PCoA axes),
#'        \eqn{m}, to run the discriminant analysis on. Default: `NULL`,
#'        which chooses the \eqn{m} with the highest leave-one-out
#'        classification success, as Anderson & Willis (2003) do (see
#'        "Choosing m"). At most \eqn{n - g} for \eqn{n} samples in
#'        \eqn{g} groups, and the number of PCoA axes with positive
#'        eigenvalues.
#'
#' @param max_pcoa_axes  With `n_pcoa_axes = NULL`, the largest \eqn{m} to
#'        consider. Default: `NULL`, which means 50. Reduced, silently, to
#'        \eqn{n - g - 1} (the most that leave-one-out classification can
#'        use) and to the number of PCoA axes with positive eigenvalues.
#'        `BiodiversityR::CAPdiscrim()`'s `mmax` (whose default is 10).
#'        Ignored when `n_pcoa_axes` is given.
#'
#' @param n_axes  The number of discriminant (canonical) axes to return.
#'        Default: `NULL`, which means as many as there are, up to 3. There
#'        are at most \eqn{g - 1}, and no more than \eqn{m}; a larger
#'        `n_axes` is reduced to that.
#'
#' @param permutations  The number of permutations for the test of the
#'        trace statistic. Default: `999`. Alternatively, a matrix with one
#'        permutation of `1:n` per row, as from `permute::shuffleSet()`, to
#'        use instead of random ones. `0` skips the test.
#'
#' @return A numeric matrix with one row per sample and one column per
#'   canonical axis, named `"CAP1"`, `"CAP2"`, etc., as described in
#'   "Ordination result" below. The values are the discriminant scores of
#'   [lda()] run on the first \eqn{m} PCoA axes. Attributes:
#'
#'   * `method`: `"capda"`.
#'   * `eig`: Each axis's squared canonical correlation \eqn{\delta^2}
#'     between the \eqn{m} PCoA axes and the groups: the proportion of the
#'     variance of the axis's scores that lies between the groups (between
#'     0 and 1).
#'   * `prop`: Each axis's proportion of `trace_stat`, the sum of
#'     \eqn{\delta^2} over all the canonical axes.
#'   * `n_pcoa_axes`: \eqn{m}, the number of PCoA axes used.
#'   * `success`: Leave-one-out classification success: a data frame with
#'     one row per \eqn{m} considered (every \eqn{m} from 1 to the maximum
#'     when it was chosen, else just the given one), and columns
#'     `n_pcoa_axes` (\eqn{m}), `overall` (the proportion of all samples
#'     classified correctly) and one per group (the proportion of that
#'     group's samples classified correctly), named by group.
#'   * `loo_table`: The leave-one-out classification table at the chosen
#'     \eqn{m}: the number of samples of each group (rows) given to each
#'     group (columns), with a column `NA` for samples that couldn't be
#'     classified, if any.
#'   * `loo_class`: Each sample's leave-one-out class at the chosen
#'     \eqn{m}: a factor with the levels of `groups`, `NA` where it couldn't
#'     be classified.
#'   * `trace_stat`: The trace statistic: the sum of \eqn{\delta^2} over
#'     all the canonical axes (Pillai's trace, for a MANOVA of the
#'     \eqn{m} PCoA axes).
#'   * `p_value`: The permutation p-value of `trace_stat`: the proportion
#'     of permutations, counting the observed one, with a trace at least as
#'     large. `NA` without permutations.
#'   * `n_perm`: The number of permutations.
#'   * `trace_perm`: The trace of each permutation.
#'   * `scaling`: The coefficients of the canonical axes on the PCoA axes:
#'     one row per PCoA axis and one column per canonical axis, as
#'     [lda()]'s `scaling`.
#'   * `groups`: Each sample's group, as a factor.
#'
#' @details
#'
#' **Method**
#'
#' Canonical analysis of principal coordinates (Anderson & Willis 2003)
#' finds the axes that best separate known groups of samples in the space
#' of any distance measure. It computes the principal coordinates of the
#' distances ([pcoa()]), and runs linear discriminant analysis ([lda()]) on
#' the first \eqn{m} of them. `capda(dm, groups, n_pcoa_axes = m)` gives
#' the same scores as `lda(pcoa(dm, n_axes = m), groups)`.
#'
#' The canonical correlations \eqn{\delta} are the correlations between
#' each axis's scores and the groups; \eqn{\delta^2}{delta^2} is the
#' proportion of the variance of the scores that lies between the groups.
#' They relate to [lda()]'s `eig` (\eqn{F}, the ratio of the between- to
#' the within-group mean squares) by
#' \eqn{\delta^2 = F (g - 1) / (F (g - 1) + n - g)}{delta^2 = F (g - 1) / (F (g - 1) + n - g)}.
#'
#'
#' **Choosing m**
#'
#' Few PCoA axes may miss the differences between the groups; too many fit
#' noise, and the discriminant analysis needs fewer axes than samples. With
#' `n_pcoa_axes = NULL`, `capda()` classifies each sample by the
#' discriminant analysis of all the other samples (leave-one-out), for each
#' \eqn{m} from 1 to `max_pcoa_axes`, and uses the \eqn{m} that classifies
#' the most samples correctly; on ties, the smallest. A sample that can't
#' be classified counts as misclassified: one alone in its group (which,
#' left out, leaves its group without a mean), and every sample once the
#' within-group covariance of the other samples is singular. The success
#' for every \eqn{m} is in the `success` attribute, so you can check that
#' the chosen \eqn{m} is not an isolated peak. The classification uses
#' each group's proportion of the samples as its prior probability, as
#' `MASS::lda(CV = TRUE)` does by default.
#'
#'
#' **Permutation test**
#'
#' The trace statistic, \eqn{\sum \delta^2}{sum(delta^2)}, measures how
#' well the \eqn{m} PCoA axes separate the groups. The test permutes the
#' samples' groups (freely, or within `strata`) and recomputes the trace on
#' the same \eqn{m} axes, as Anderson & Willis (2003) do. The choice of
#' \eqn{m} is not repeated for each permutation. `BiodiversityR::CAPdiscrim()`
#' tests the classification success instead.
#'
#'
#' **Validate the separation**
#'
#' Like [lda()] and [plsda()], `capda()` is supervised: its plots separate
#' the groups, if not as well, even when the groups don't differ. The
#' leave-one-out success and the permutation test are the evidence; the
#' plot is not.
#'
#'
#' **Matching BiodiversityR**
#'
#' `BiodiversityR::CAPdiscrim(dm ~ groups, data, m = 0, mmax = 10)` gives
#' the same \eqn{m}, classification success (its `percent` and
#' `percent.level`, as percentages), leave-one-out classes (`CV`) and
#' discriminant scores (`x`, up to each axis's sign) as `capda(dm, groups,
#' data, max_pcoa_axes = 10)`, and its `F` is [lda()]'s `eig`. Its
#' `manova` holds `trace_stat`, as Pillai's trace. It stops with an error
#' when a group has a single sample; `capda()` counts that sample as
#' misclassified. PRIMER-PERMANOVA+'s CAP scales the axes differently, and
#' does not use prior probabilities.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::capscale()` ("constrained analysis of principal coordinates")
#' is distance-based redundancy analysis, ecodive's [capscale()] and
#' [dbrda()], not this method. In [dist_ord()], `'CAP'` selects
#' `capscale`; use `'capda'` for this one. BiodiversityR's version is
#' `CAPdiscrim()`.
#'
#'
#' **Axis signs and reproducibility**
#'
#' An axis's sign is arbitrary, so different software can return mirror
#' images of the same axis. `capda()` flips each axis so that its largest
#' score in absolute value is positive, and flips its column of `scaling`
#' with it. Results are identical from run to run and for any number of
#' `cpus`; random permutations depend only on `seed`.
#'
#' @references
#' Anderson, M. J., & Willis, T. J. (2003). Canonical analysis of principal
#' coordinates: a useful method of constrained ordination for ecology.
#' *Ecology*, 84(2), 511-525. \doi{10.1890/0012-9658(2003)084[0511:CAOPCA]2.0.CO;2}
#'
#' Anderson, M. J., Gorley, R. N., & Clarke, K. R. (2008). *PERMANOVA+ for
#' PRIMER: Guide to Software and Statistical Methods*. PRIMER-E, Plymouth.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [lda()], [pcoa()], [plsda()], [permanova()], [dist_ord()]
#'
#' @export
#' @examples
#'     # CAP of iris, from Bray-Curtis distances.
#'     dm  <- bray(iris[, 1:4])
#'     ord <- capda(dm, 'Species', iris, permutations = 99)
#'     plot(ord, 'Species', iris)
#'
#'     # The number of PCoA axes chosen, and the classification success.
#'     attr(ord, 'n_pcoa_axes')
#'     head(attr(ord, 'success'))
#'     attr(ord, 'loo_table')
#'
#'     # The squared canonical correlations, and the permutation test.
#'     attr(ord, 'eig')
#'     attr(ord, 'p_value')
#'
capda <- function (
    dm, groups, data = NULL, n_pcoa_axes = NULL, max_pcoa_axes = NULL,
    n_axes = NULL, permutations = 999L, strata = NULL, seed = 0L,
    cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL

  validate_seed()
  validate_cpus()
  if (!is.null(n_axes))        n_axes        <- ord_k(n_axes)
  if (!is.null(n_pcoa_axes))   n_pcoa_axes   <- ord_arg(n_pcoa_axes,
    '`n_pcoa_axes` must be a whole number of at least 1, or NULL.', whole = TRUE, lo = 1)
  if (!is.null(max_pcoa_axes)) max_pcoa_axes <- ord_arg(max_pcoa_axes,
    '`max_pcoa_axes` must be a whole number of at least 1, or NULL.', whole = TRUE, lo = 1)

  dv     <- dist_vector(dm)
  n      <- dv$n
  labels <- dv$labels


  # The groups, as a factor, and the strata.
  des <- test_group_values(groups, data, n, labels, strata)
  grp <- des$groups
  grp <- if (is.factor(grp)) droplevels(grp) else factor(grp)
  g   <- nlevels(grp)

  if (g < 2)
    stop('`groups` must have at least two different values.', call. = FALSE)

  pm <- perm_matrix(permutations, n)
  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  # How many PCoA axes to compute.
  if (!is.null(n_pcoa_axes)) {
    if (n_pcoa_axes > n - g)
      stop(sprintf(paste0(
        '`n_pcoa_axes` must be at most the number of samples minus the ',
        'number of groups: %d - %d = %d.'), n, g, n - g), call. = FALSE )
    k <- n_pcoa_axes
  }
  else {
    if (n - g - 1 < 1)
      stop(sprintf(paste0(
        'capda() needs at least two more samples than groups to choose ',
        '`n_pcoa_axes` (%d samples, %d groups); or give `n_pcoa_axes`.'), n, g),
        call. = FALSE )
    k <- min(if (is.null(max_pcoa_axes)) 50L else max_pcoa_axes, n - g - 1L)
  }


  # The PCoA axes with positive eigenvalues, as pcoa().
  res <- .Call(C_pcoa, dv$dm, n, as.integer(k), cpus)
  if (!res$converged) # nocov start
    warning('PCoA did not converge; coordinates may be inexact.', call. = FALSE)
  # nocov end
  n_pos  <- sum(res$eig > max(0, sqrt(.Machine$double.eps) * res$eig[1]), na.rm = TRUE)
  points <- res$points[, seq_len(n_pos), drop = FALSE]

  if (!is.null(n_pcoa_axes) && n_pos < n_pcoa_axes)
    stop(sprintf(
      'Only %d PCoA axes have positive eigenvalues: `n_pcoa_axes` must be at most %d.',
      n_pos, n_pos ), call. = FALSE )
  if (n_pos < 1)
    stop('The distances have no PCoA axes with positive eigenvalues.', call. = FALSE)


  # Leave-one-out classification success for each m, NA (misclassified)
  # from the first m at which the groups barely vary within.
  prior   <- lda_prior(NULL, grp)
  loo     <- lda_loo(points, grp, prior, cpus)$class
  m_sing  <- capda_singular_from(points, grp)
  if (m_sing <= n_pos) loo[, m_sing:n_pos] <- NA_integer_
  ms      <- if (is.null(n_pcoa_axes)) seq_len(n_pos) else n_pcoa_axes
  code    <- as.integer(grp)
  correct <- !is.na(loo[, ms, drop = FALSE]) & loo[, ms, drop = FALSE] == code
  hits    <- colSums(correct)
  m       <- ms[[which.max(hits)]]

  success <- data.frame(n_pcoa_axes = ms, overall = hits / n)
  by_grp  <- rowsum(correct * 1, grp, reorder = TRUE) / as.vector(table(grp))
  for (lev in levels(grp)) success[[lev]] <- unname(by_grp[lev, ])
  names(success) <- make.unique(names(success))

  cls       <- factor(levels(grp)[loo[, m]], levels = levels(grp))
  loo_table <- table(group = grp, loo_class = cls, useNA = 'ifany')
  if (!is.null(labels)) names(cls) <- labels


  # The discriminant analysis of the first m axes.
  x   <- points[, seq_len(m), drop = FALSE]
  fit <- if (m < m_sing) lda_fit(x, grp, prior, cpus)
  if (is.null(fit))
    stop(sprintf(
      'The pooled within-group covariance of the first %d PCoA axes is singular%s',
      m, if (m > 1) ': give a smaller `n_pcoa_axes`.' else '.' ), call. = FALSE )

  r <- length(fit$sv)
  if (r == 0)
    warning('the groups have the same means: there are no discriminant axes', call. = FALSE)

  # Squared canonical correlations, from the canonical F.
  f      <- fit$sv^2 * (g - 1)
  delta2 <- f / (f + (n - g))

  n_ax <- min(if (is.null(n_axes)) 3L else n_axes, r)
  axes <- if (n_ax > 0) paste0('CAP', seq_len(n_ax)) else NULL
  flip <- lda_flip(
    fit$scores[,  seq_len(n_ax), drop = FALSE],
    fit$scaling[, seq_len(n_ax), drop = FALSE] )
  dimnames(flip$scaling) <- list(paste0('PCoA', seq_len(m)), axes)


  # The trace statistic and its permutation test, on the m axes scaled to
  # unit length.
  q <- sweep(x, 2L, sqrt(colSums(x^2)), '/')
  tr <- .Call(
    C_capda_trace, q, code - 1L, g, pm$perms, des$strata,
    as.integer(c(pm$n_perm, seed, cpus)) )
  trace_perm <- tr[-1L]

  if (!is.null(labels)) names(grp) <- labels

  ord_result(
    flip$scores, 'capda', labels, prefix = 'CAP', from = dv$title,
    eig         = delta2[seq_len(n_ax)],
    prop        = if (r > 0) delta2[seq_len(n_ax)] / sum(delta2) else numeric(0),
    n_pcoa_axes = m,
    success     = success,
    loo_table   = loo_table,
    loo_class   = cls,
    trace_stat  = tr[[1L]],
    p_value     = perm_p(tr[[1L]], trace_perm),
    n_perm      = pm$n_perm,
    trace_perm  = trace_perm,
    scaling     = flip$scaling,
    groups      = grp )
}



# The first m at which the groups barely vary within the first m PCoA axes
# (`points`, n x M), or M + 1. lda_fit() and lda_loo() scale each axis by
# its within-group standard deviation, so they can't tell an axis along
# which the groups don't vary within, but for rounding error, from any
# other. With the axes scaled to unit length, W = I - C is the within-group
# sum of squares, and the pivots of its Cholesky (LDL') factorization are,
# in turn, the proportion of each axis's variation that lies within the
# groups and not along the earlier axes. The fit is singular from the first
# pivot below 1e-8: a within-group standard deviation below 1e-4 times the
# total, as lda()'s tolerance.
capda_singular_from <- function (points, grp) {

  q <- sweep(points, 2L, sqrt(colSums(points^2)), '/')
  s <- rowsum(q, grp) / sqrt(as.vector(table(grp)))
  w <- crossprod(q) - crossprod(s)

  M <- ncol(w)
  L <- matrix(0, M, M)
  d <- numeric(M)

  for (j in seq_len(M)) {
    h    <- seq_len(j - 1L)
    lj   <- L[j, h]
    d[j] <- w[j, j] - sum(lj^2 * d[h])
    if (!(d[j] > 1e-8)) return (j)
    if (j < M) {
      i <- (j + 1L):M
      L[i, j] <- (w[i, j] - L[i, h, drop = FALSE] %*% (lj * d[h])) / d[j]
    }
  }

  M + 1L
}
