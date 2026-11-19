# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Linear Discriminant Analysis (LDA)
#'
#' Supervised ordination of samples from a few features: Fisher's linear
#' discriminants, the axes along which known groups are furthest apart
#' relative to the variation within them. Gives the same scores,
#' coefficients and proportion of trace as `MASS::lda()`, and, up to a
#' constant factor, the same scores as scikit-learn's
#' `LinearDiscriminantAnalysis`.
#'
#' @inherit pseudocount_section
#' @inheritParams plsda
#'
#' @param counts  A numeric matrix of feature values (samples \eqn{\times}
#'        features), or any of the table types listed under "Input Types".
#'        It needs more samples than features: at least the number of
#'        features plus the number of groups (see "Singular tables").
#'
#' @param prior   The prior probability of each group: a numeric vector
#'        with one positive value per group, named by group or else in the
#'        order of the groups (a factor's levels, otherwise sorted), adding
#'        up to 1. It weights each group's mean in the centre of the scores
#'        and in the between-group variance, as in `MASS::lda()`. Default:
#'        `NULL`, the proportion of samples in each group.
#'
#' @param n_axes   The number of discriminant axes to return. Default:
#'        `NULL`, which means as many as there are, up to 3. There are at
#'        most \eqn{g - 1} for \eqn{g} groups, and no more than the
#'        number of features; a larger `n_axes` is reduced to that.
#'
#' @return A numeric matrix with one row per sample and one column per
#'   discriminant axis, named `"LD1"`, `"LD2"`, etc., as described in
#'   "Ordination result" below. The values are the discriminant scores, as
#'   `predict(MASS::lda(...))$x` gives them. Attributes:
#'
#'   * `method`: `"lda"`.
#'   * `eig`: The ratio of the between-group to the within-group variance
#'     (mean squares) of the scores on each axis: `MASS::lda()`'s `svd`,
#'     squared (the "canonical F statistics").
#'   * `prop`: Each axis's proportion of the sum of `eig` over all the
#'     discriminant axes: `MASS::lda()`'s "proportion of trace", and
#'     scikit-learn's `explained_variance_ratio_`.
#'   * `scaling`: The coefficients of the discriminants: one row per
#'     feature and one column per axis (`MASS::lda()`'s `scaling`). The
#'     scores are the table, minus `prior %*% means`, times `scaling`.
#'   * `means`: The mean of each feature in each group: one row per group
#'     and one column per feature.
#'   * `prior`: The prior probabilities, named by group.
#'   * `groups`: Each sample's group, as a factor whose levels are the
#'     groups, in the order of the rows of `means`.
#'
#' @details
#'
#' **Method**
#'
#' Let \eqn{X} be the normalized table, \eqn{m_k} the mean of group \eqn{k}
#' (of \eqn{n_k} samples, with prior \eqn{\pi_k}), and
#' \eqn{\bar{x} = \sum_k \pi_k m_k}{xbar = sum(pi_k m_k)}. With \eqn{W}
#' the pooled within-group covariance (each sample minus its group's mean,
#' divisor \eqn{n - g}) and \eqn{B} the between-group covariance
#' \eqn{\sum_k n \pi_k (m_k - \bar{x})(m_k - \bar{x})^T / (g - 1)}{sum(n pi_k (m_k - xbar)(m_k - xbar)') / (g - 1)}
#' (with the default priors, \eqn{n \pi_k = n_k}), the discriminant axes
#' \eqn{a} maximize \eqn{a^T B a / a^T W a}, the ratio of the
#' between-group to the within-group variance of the scores (Fisher 1936).
#' They are the eigenvectors of \eqn{W^{-1} B}, and the ratios its
#' eigenvalues, `eig`. Each axis is scaled so that the pooled within-group
#' variance of its scores is 1, and the scores are
#' \eqn{(x - \bar{x})^T a}{(x - xbar)' a}, for every sample \eqn{x}.
#'
#' As in scikit-learn's `"svd"` solver, \eqn{W} and \eqn{B} are never
#' formed: the axes come from the SVD of the within-group residuals (each
#' feature scaled by its within-group standard deviation), which spheres
#' \eqn{W}, then the SVD of the weighted group means in that sphered
#' space. Axes whose singular value is below \eqn{10^{-4}} times the first
#' are dropped, as both references do.
#'
#'
#' **Singular tables**
#'
#' LDA needs \eqn{W} to be invertible: more samples than features plus
#' groups, and no feature that is constant within every group or a linear
#' combination of others within groups (as relative abundances that add up
#' to 1 are). Otherwise `lda()` stops with an error (scikit-learn and
#' `MASS::lda()` drop the null directions with a warning, but the axes
#' they then report depend on rounding). Microbiome tables have far more
#' taxa than samples, so use [capda()] instead: canonical analysis of
#' principal coordinates (Anderson & Willis 2003) runs LDA on the first
#' few axes of [pcoa()], choosing how many by how well they classify the
#' samples, and tests the separation of the groups by permutation.
#' `capda(dm, groups, n_pcoa_axes = 5)` gives the same scores as
#' `lda(pcoa(dm, n_axes = 5), groups)`. A few axes of [pca()] also work.
#'
#'
#' **Validate the separation**
#'
#' LDA is supervised: it looks for the directions that separate the
#' groups, so its plots separate them, if not as readily as [plsda()]'s,
#' even when the groups don't differ. Don't read the plot as evidence of a
#' difference. Test it with a permutation test of the distances, such as
#' [permanova()], or estimate the error rate of classifying new samples by
#' cross-validation.
#'
#'
#' **Matching MASS and scikit-learn**
#'
#' `lda()` gives the same scores, `scaling`, group means and proportion of
#' trace as `MASS::lda()` and its `predict()`, up to the sign of each axis,
#' and `eig` is the square of its `svd`. scikit-learn's
#' `LinearDiscriminantAnalysis(solver = "svd")` divides the within-group
#' covariance by \eqn{n} rather than \eqn{n - g}, so its `transform()`
#' scores and `scalings_` are larger by \eqn{\sqrt{n / (n - g)}}{sqrt(n / (n - g))};
#' its `explained_variance_ratio_` is `prop`. Its priors that don't add up
#' to 1 are rescaled with a warning; here they are an error.
#'
#'
#' **Same names in other packages**
#'
#' `MASS::lda()` shares its name with ecodive's `lda()`, and MASS is often
#' attached: whichever package was attached last masks the other's. Write
#' `ecodive::lda()`, or use [count_ord()], which always runs ecodive's, to
#' be sure. Elsewhere, "LDA" also means latent Dirichlet allocation (topic
#' models), and LEfSe's "LDA score" is an effect size from this method.
#'
#'
#' **Axis signs and reproducibility**
#'
#' An axis's sign is arbitrary, so different software can return mirror
#' images of the same axis. `lda()` flips each axis so that its largest
#' score in absolute value is positive, and flips its column of `scaling`
#' with it. Results are identical from run to run and for any number of
#' `cpus`.
#'
#' @references
#' Fisher, R. A. (1936). The use of multiple measurements in taxonomic
#' problems. *Annals of Eugenics*, 7(2), 179-188.
#' \doi{10.1111/j.1469-1809.1936.tb02137.x}
#'
#' Anderson, M. J., & Willis, T. J. (2003). Canonical analysis of principal
#' coordinates: a useful method of constrained ordination for ecology.
#' *Ecology*, 84(2), 511-525. \doi{10.1890/0012-9658(2003)084[0511:CAOPCA]2.0.CO;2}
#'
#' Venables, W. N., & Ripley, B. D. (2002). *Modern Applied Statistics with
#' S* (4th ed.). Springer. \doi{10.1007/978-0-387-21706-2}
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830.
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [capda()], [plsda()], [pca()], [pcoa()], [permanova()], [count_ord()]
#'
#' @export
#' @examples
#'     # LDA of iris: three groups, so two discriminant axes.
#'     ord <- lda(iris[, 1:4], 'Species', iris)
#'     plot(ord)
#'     attr(ord, 'prop')
#'     attr(ord, 'scaling')
#'
#'     # Microbiome tables have more taxa than samples: use a few PCoA axes,
#'     # or capda(), which chooses how many.
#'     grp <- rep(c('a', 'b'), length.out = nrow(ex_counts))
#'     pc  <- pcoa(bray(ex_counts), n_axes = 2)
#'     lda(pc, grp)
#'
#'     # Supervised: always validate the separation, e.g. with PERMANOVA.
#'     permanova(bray(ex_counts), grp, permutations = 99)
#'
lda <- function (
    counts, groups, data = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, prior = NULL, n_axes = NULL, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL

  validate_counts()
  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()
  if (!is.null(n_axes)) n_axes <- ord_k(n_axes)

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')


  # Dimensions and names, with samples first.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  n <- as.integer(dims[[1]])
  p <- as.integer(dims[[2]])


  # The groups, as a factor.
  tg  <- test_groups(groups, data, n, dn[[1]], ref = '`counts`')
  grp <- tg$data[[tg$col]]
  grp <- if (is.factor(grp)) droplevels(grp) else factor(grp)
  g   <- nlevels(grp)

  if (g < 2)
    stop('`groups` must have at least two different values.', call. = FALSE)
  if (n <= g)
    stop(sprintf(
      'lda() needs more samples than groups: %d samples, %d groups.', n, g ),
      call. = FALSE )

  prior <- lda_prior(prior, grp)

  if (p > n - g) lda_singular(n, p, g)


  # The normalized table, dense, samples x features.
  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  x <- .Call(C_normalize, counts, margin, norm, as.double(pseudocount), TRUE, cpus)
  if (margin == 2L) x <- t(x)
  storage.mode(x) <- 'double'
  if (!all(is.finite(x)))
    stop('The normalized `counts` contain NA, NaN, or infinite values.', call. = FALSE) # nocov

  fit <- lda_fit(x, grp, prior, cpus)
  if (is.null(fit)) lda_singular(n, p, g)

  r <- length(fit$sv)
  if (r == 0)
    warning('the groups have the same means: there are no discriminant axes', call. = FALSE)
  k <- min(if (is.null(n_axes)) 3L else n_axes, r)
  axes <- if (k > 0) paste0('LD', seq_len(k)) else NULL

  flip    <- lda_flip(
    fit$scores[,  seq_len(k), drop = FALSE],
    fit$scaling[, seq_len(k), drop = FALSE] )
  scores  <- flip$scores
  scaling <- flip$scaling
  eig     <- fit$sv[seq_len(k)]^2

  dimnames(scaling)   <- list(dn[[2]], axes)
  dimnames(fit$means) <- list(levels(grp), dn[[2]])
  if (!is.null(dn[[1]])) names(grp) <- dn[[1]]

  ord_result(
    scores, 'lda', dn[[1]], prefix = 'LD',
    eig     = eig,
    prop    = eig / fit$sv_sum2,
    scaling = scaling,
    means   = fit$means,
    prior   = prior,
    groups  = grp )
}



# lda()'s `prior`: NULL for the group proportions, or one positive value
# per group (named by group, or in level order) adding up to 1. Returns the
# priors, named, in level order, rescaled to add up to exactly 1.
lda_prior <- function (prior, grp) {

  lev <- levels(grp)

  if (is.null(prior)) {
    prior <- as.vector(table(grp)) / length(grp)
  }
  else {
    if (!is.numeric(prior) || !is.null(dim(prior)) || length(prior) != length(lev) ||
        !all(is.finite(prior)) || any(prior <= 0))
      stop(sprintf(
        '`prior` must be %d positive numbers, one per group.', length(lev) ),
        call. = FALSE )
    if (!is.null(names(prior))) {
      if (!setequal(names(prior), lev) || anyDuplicated(names(prior)))
        stop(
          '`prior` names must be the groups: ',
          paste0('"', lev, '"', collapse = ', '), '.', call. = FALSE )
      prior <- prior[lev]
    }
    if (abs(sum(prior) - 1) > 1e-5)
      stop('`prior` must add up to 1.', call. = FALSE)
    prior <- as.vector(prior) / sum(prior)
  }

  names(prior) <- lev
  prior
}



# Flips each axis so that its largest score (in absolute value) is
# positive, with its column of `scaling`: lda() and capda().
lda_flip <- function (scores, scaling) {
  for (a in seq_len(ncol(scores))) {
    if (!identical(ord_flip(scores[, a, drop = FALSE]), scores[, a, drop = FALSE])) {
      scores[,  a] <- -scores[,  a]
      scaling[, a] <- -scaling[, a]
    }
  }
  list(scores = scores, scaling = scaling)
}



lda_singular <- function (n, p, g) {
  stop(
    'The pooled within-group covariance is singular ',
    sprintf('(%d samples, %d groups, %d features): ', n, g, p),
    'LDA needs more samples than features plus groups, and no feature that ',
    'is constant within groups or collinear with others. Use fewer ',
    'features, or capda(), which runs LDA on the first few PCoA axes of ',
    'the distances between samples.', call. = FALSE )
}



#------------------------------------------------------------------------------
# The dense core, for lda() and capda(). `x` is a dense n x p double matrix
# (samples x features), `grp` a factor without unused levels (n > nlevels),
# `prior` one positive value per level, adding up to 1.
#
# lda_fit(): list(scores = n x r, scaling = p x r, means = g x p, sv = the r
# between-group singular values (eig = sv^2), sv_sum2 = their sum of squares
# over all axes), with r the number of discriminant axes and signs not yet
# fixed; or NULL if the within-group covariance is singular.
#
# lda_loo(): leave-one-out classification, as MASS::lda(CV = TRUE): sample i
# classified by the fit to the other samples (with the same priors), for the
# first m columns of x, for every m = 1..p at once. Returns list(class = an
# n x p integer matrix of level codes, column m for the first m columns, NA
# where that fit is singular; posterior = n x g, for all p columns).
#------------------------------------------------------------------------------
lda_fit <- function (x, grp, prior, cpus = 1L) {
  res <- .Call(C_lda, x, as.integer(grp), nlevels(grp), as.double(prior), as.integer(cpus))
  if (res$rank < ncol(x)) return (NULL)
  res[c('scores', 'scaling', 'means', 'sv', 'sv_sum2')]
}

lda_loo <- function (x, grp, prior, cpus = 1L) {
  .Call(C_lda_loo, x, as.integer(grp), nlevels(grp), as.double(prior), as.integer(cpus))
}
