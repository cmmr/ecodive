# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' Supervised ordination of samples from two groups: one predictive
#' component that separates the groups, and orthogonal components that
#' describe the rest of the table's variation, unrelated to the groups.
#' OPLS-DA (Trygg & Wold 2002) gives the same predictions as a one-component
#' [plsda()], rotated so that the plot of the first two components shows the
#' group difference along the first axis alone. Gives the same results as
#' pyopls and `ropls::opls()`, computed straight from sparse input without
#' forming any dense table or \eqn{n \times n}{n x n} matrix.
#'
#' @inherit pseudocount_section
#' @inheritParams plsda
#'
#' @param groups  The known group of each sample: a vector or factor with
#'        one value per sample, matched to the samples in `counts` by name
#'        when both are named, else in the same order. Or, when `data` is
#'        given, the name of one of its columns (or a vector with one value
#'        per row of `data`). There must be exactly two groups (see "Two
#'        groups only"). A factor's unused levels are dropped.
#'
#' @param scale   Scale each feature to unit variance (after `norm`, before
#'        the analysis), as `ropls::opls()` (`scaleC = 'standard'`) and
#'        pyopls (`scale = True`) do by default, so that every feature has
#'        equal weight. Features with constant values are left out: their
#'        weights, loadings and VIP are 0. Default: `TRUE`
#'
#' @param n_orthogonal   The number of orthogonal components. Default: `1`.
#'        `0` gives a one-component [plsda()]. Reduced to \eqn{n - 2} for
#'        \eqn{n} samples, or to the number of features less one, if fewer;
#'        and, with a warning, to the number that the table can fit (see
#'        "Method").
#'
#' @return A numeric matrix with one row per sample and one column per
#'   component: the predictive component `"p1"`, then the orthogonal ones
#'   `"o1"`, `"o2"`, etc., as described in "Ordination result" below. The
#'   values are the scores (`ropls`' `scoreMN` and `orthoScoreMN`). So
#'   `plot(ord)` draws the usual OPLS-DA score plot, `p1` against `o1`.
#'   Attributes:
#'
#'   * `method`: `"oplsda"`.
#'   * `eig`: The variance of the table that each component explains,
#'     \eqn{|t|^2 |p|^2 / (n - 1)} for scores \eqn{t} and loadings \eqn{p}.
#'     (As [plsda()]'s: OPLS components are not eigenvectors.)
#'   * `prop`: The proportion of the table's total variance that each
#'     component explains, `eig / trace`: \eqn{R^2X} of each component, as
#'     `ropls` reports it (`modelDF$R2X`).
#'   * `trace`: The total variance of the table: the sum of all features'
#'     variances (after `norm` and `scale`), as [pca()]'s `trace`.
#'   * `r2x`: \eqn{R^2X}, the proportion of the table's variance that all the
#'     components explain together: `sum(prop)` (`ropls`' `R2X(cum)`).
#'   * `r2y`: \eqn{R^2Y}, the proportion of the variance of the groups'
#'     indicator that the model explains (`ropls`' `R2Y(cum)`; the \eqn{R^2}
#'     of the fitted values, as in pyopls). Only the predictive component
#'     contributes.
#'   * `loadings`: The loadings \eqn{p}: one row per feature and one column
#'     per component (`ropls`' `loadingMN` and `orthoLoadingMN`).
#'   * `weights`: The weights \eqn{w}, unit-length columns (`ropls`'
#'     `weightMN` and `orthoWeightMN`, pyopls's `W_ortho_`).
#'   * `y_loading`: The Y loading \eqn{q} of the predictive component, for
#'     the indicator of the second group as centred (and scaled).
#'   * `vip`: Variable importance in projection (VIP) of the predictive
#'     component: one value per feature (`ropls`' `vipVn`). See "Method".
#'   * `groups`: Each sample's group, as a factor with two levels, in the
#'     order of the rows of `centroids`.
#'   * `centroids`: The mean score of each group: two rows and one column
#'     per component.
#'
#' @details
#'
#' **Method**
#'
#' Let \eqn{X} be the normalized table, its features centred (and, with
#' `scale = TRUE`, scaled to unit variance), and \eqn{y} the indicator of the
#' second group (1 for its samples, 0 for the first group's), centred (and
#' scaled) the same way. The predictive weights are
#' \eqn{w = X^T y / |X^T y|}{w = X'y / |X'y|}, as for a one-component PLS.
#' O-PLS (Trygg & Wold 2002) then removes, one at a time, the variation of
#' \eqn{X} that is orthogonal to \eqn{y}. With \eqn{X_a} what is left of
#' \eqn{X} before orthogonal component \eqn{a}:
#'
#' * \eqn{t = X_a w} and \eqn{p = X_a^T t / |t|^2}{p = X_a't / |t|^2};
#' * the orthogonal weights \eqn{w_o} are the part of \eqn{p} orthogonal to
#'   \eqn{w}, \eqn{p - (w^T p) w}{p - (w'p) w}, scaled to unit length;
#' * the orthogonal scores are \eqn{t_o = X_a w_o} and the loadings
#'   \eqn{p_o = X_a^T t_o / |t_o|^2}{p_o = X_a't_o / |t_o|^2};
#' * \eqn{X_{a+1} = X_a - t_o p_o^T}{X_(a+1) = X_a - t_o p_o'}.
#'
#' The predictive component is fitted last, to the filtered table
#' \eqn{X_k}: scores \eqn{t_p = X_k w}, loadings
#' \eqn{p_p = X_k^T t_p / |t_p|^2}{p_p = X_k't_p / |t_p|^2}, and Y loading
#' \eqn{q = y^T t_p / |t_p|^2}{q = y't_p / |t_p|^2}. Its scores are those of
#' the one-component PLS of \eqn{y} on \eqn{X_k}, and \eqn{t_p q} are the
#' fitted values of \eqn{y}, the same as a one-component [plsda()]'s up to
#' scaling: \eqn{R^2Y = |t_p|^2 q^2 / |y|^2}. All the scores are
#' orthogonal, so the components' `eig`s add up to what they explain
#' together.
#'
#' Orthogonal components stop early, with a warning, when \eqn{|w_o|} (before
#' scaling) falls to \eqn{\sqrt{\epsilon} \approx 1.5 \times
#' 10^{-8}}{sqrt(eps) = 1.5e-8} times \eqn{|p|}: then what is left of the
#' table has no variation orthogonal to the groups, beyond rounding error.
#'
#' VIP (Galindo-Prieto et al. 2014) is, for one predictive component,
#' \eqn{VIP_j = \sqrt{P} |p_{pj}| / |p_p|}{VIP_j = sqrt(P) |p_pj| / |p_p|}
#' for \eqn{P} features (constant ones included, when `scale` leaves them
#' out): each feature's share of the predictive loadings. The mean of the
#' squared VIPs is 1, so features with a VIP above 1 are often called
#' important.
#'
#' The table is never centred in memory: products with \eqn{X_a} are
#' computed from the non-zero values alone, in parallel, as in [pca()], so
#' a sparse table is never made dense.
#'
#'
#' **Two groups only**
#'
#' OPLS-DA, as pyopls, `ropls` and most tools implement it, models one
#' response: here, membership of the second group. With more than two groups,
#' use [plsda()].
#'
#'
#' **Choosing `norm`**
#'
#' For microbiome counts, `norm = 'clr'` (with a `pseudocount`) is the usual
#' choice: OPLS-DA then works on log-ratios, which suit compositional data.
#'
#'
#' **Validate the separation**
#'
#' OPLS-DA is supervised: its predictive component separates the groups even
#' when they don't differ (with more features than samples, almost always).
#' Don't read the plot, or a high \eqn{R^2Y}, as evidence of a difference.
#' Test it with a permutation test of the distances, such as [permanova()],
#' or estimate the error rate of classifying new samples by cross-validation
#' (\eqn{Q^2}, which `oplsda()` doesn't compute).
#'
#'
#' **Matching pyopls and ropls**
#'
#' `oplsda()` gives the same orthogonal scores, weights and loadings as
#' pyopls's `OPLS(n_orthogonal, scale = True)` (`T_ortho_`, `W_ortho_` and
#' `P_ortho_`), and the same predictive scores and loadings as scikit-learn's
#' `PLSRegression(1, scale = False)` fitted to pyopls's filtered table, up to
#' the signs of the components. It also gives the same scores, loadings,
#' \eqn{R^2X}, \eqn{R^2Y} and VIP as `ropls::opls(predI = 1, orthoI =
#' n_orthogonal)`. pyopls's own \eqn{R^2X} (its `score()`) is the proportion
#' of variance left in the filtered table, which is a different measure:
#' `1 - sum(prop[-1])`.
#'
#'
#' **Same names in other packages**
#'
#' `ropls::opls()` does OPLS-DA under another name; no other common package
#' has an `oplsda()`.
#'
#'
#' **Component signs and reproducibility**
#'
#' A component's sign is arbitrary, so different software can return mirror
#' images of the same component. `oplsda()` flips each component so that its
#' largest score in absolute value is positive, and flips its weights and
#' loadings (and the Y loading) with it. Results are identical from run to
#' run and for any number of `cpus`.
#'
#' @references
#' Trygg, J., & Wold, S. (2002). Orthogonal projections to latent structures
#' (O-PLS). *Journal of Chemometrics*, 16(3), 119-128. \doi{10.1002/cem.695}
#'
#' Bylesjö, M., Rantalainen, M., Cloarec, O., Nicholson, J. K., Holmes, E.,
#' & Trygg, J. (2006). OPLS discriminant analysis: combining the strengths
#' of PLS-DA and SIMCA classification. *Journal of Chemometrics*, 20,
#' 341-351. \doi{10.1002/cem.1006}
#'
#' Galindo-Prieto, B., Eriksson, L., & Trygg, J. (2014). Variable influence
#' on projection (VIP) for orthogonal projections to latent structures
#' (OPLS). *Journal of Chemometrics*, 28(8), 623-632. \doi{10.1002/cem.2627}
#'
#' Thévenot, E. A., Roux, A., Xu, Y., Ezan, E., & Junot, C. (2015). Analysis
#' of the human adult urinary metabolome variations with age, body mass
#' index, and gender by implementing a comprehensive workflow for univariate
#' and OPLS statistical analyses. *Journal of Proteome Research*, 14(8),
#' 3322-3335. \doi{10.1021/acs.jproteome.5b00354}
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [plsda()], [pca()], [permanova()], [count_ord()]
#'
#' @export
#' @examples
#'     # OPLS-DA of two iris species, as if their measurements were counts.
#'     two    <- droplevels(subset(iris, Species != 'setosa'))
#'     counts <- round(two[, 1:4] * 10)
#'     ord    <- oplsda(counts, 'Species', two)
#'     plot(ord)
#'     attr(ord, 'r2x')
#'     attr(ord, 'r2y')
#'     attr(ord, 'vip')
#'
#'     # CLR-transformed microbiome counts, with two orthogonal components.
#'     grp <- rep(c('a', 'b'), length.out = nrow(ex_counts))
#'     ord <- oplsda(ex_counts, grp, norm = 'clr', pseudocount = 1, n_orthogonal = 2)
#'     attr(ord, 'centroids')
#'
#'     # Supervised: always validate the separation, e.g. with PERMANOVA.
#'     permanova(aitchison(ex_counts, pseudocount = 1), grp, permutations = 99)
#'
oplsda <- function (
    counts, groups, data = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, scale = TRUE, n_orthogonal = 1L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL

  validate_counts()
  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()
  k <- ord_arg(
    n_orthogonal, '`n_orthogonal` must be a whole number of at least 0.',
    whole = TRUE, lo = 0 )

  if (!isTRUE(scale) && !isFALSE(scale))
    stop('`scale` must be TRUE or FALSE.', call. = FALSE)


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


  # The groups, as a factor, and the indicator of the second: centred and,
  # with `scale`, scaled to unit variance (pyopls's _center_scale_xy()).
  tg  <- test_groups(groups, data, n, dn[[1]], ref = '`counts`')
  grp <- tg$data[[tg$col]]
  grp <- if (is.factor(grp)) droplevels(grp) else factor(grp)
  lev <- levels(grp)

  if (length(lev) < 2)
    stop('`groups` must have two different values.', call. = FALSE)
  if (length(lev) > 2)
    stop(
      '`groups` has ', length(lev), ' different values, but oplsda() ',
      'takes exactly two; use plsda() for more.', call. = FALSE )

  y <- (as.integer(grp) == 2L) + 0
  y <- y - mean(y)
  if (scale) y <- y / sqrt(sum(y^2) / (n - 1))

  k_max <- max(0L, min(k, n - 2L, p - 1L))

  res <- .Call(C_oplsda, counts, margin, norm, pseudocount, scale, y, k_max, cpus)

  if (is.null(res))
    stop(
      'the table has no covariance with `groups` ',
      '(for instance, every feature is constant).', call. = FALSE )

  kk <- ncol(res$scores) - 1L
  if (kk < k_max)
    warning(sprintf(
      paste(
        'only %d of the first %d orthogonal components could be fitted:',
        'the rest of the table has no variation orthogonal to `groups`' ),
      kk, k_max ), call. = FALSE )

  axes <- c('p1', if (kk > 0) paste0('o', seq_len(kk)))


  # VIP of the predictive component (Galindo-Prieto et al. 2014).
  pp  <- res$loadings[, 1]
  vip <- sqrt(p) * abs(pp) / sqrt(sum(pp^2))
  names(vip) <- dn[[2]]

  # R2Y: the share of y's sum of squares that t_p q explains.
  r2y <- sum(res$scores[, 1]^2) * res$y_loading^2 / sum(y^2)

  # The mean score of each group.
  centroids <- rowsum(res$scores, grp, reorder = TRUE) / as.vector(table(grp))

  feat <- list(dn[[2]], axes)
  dimnames(res$loadings) <- feat
  dimnames(res$weights)  <- feat
  dimnames(centroids)    <- list(lev, axes)
  if (!is.null(dn[[1]])) names(grp) <- dn[[1]]

  trace <- res$total / (n - 1)
  prop  <- if (trace > 0) res$eig / trace else res$eig

  ord_result(
    res$scores, 'oplsda', dn[[1]], cols = axes,
    eig       = res$eig,
    prop      = prop,
    trace     = trace,
    r2x       = sum(prop),
    r2y       = r2y,
    loadings  = res$loadings,
    weights   = res$weights,
    y_loading = res$y_loading,
    vip       = vip,
    groups    = grp,
    centroids = centroids )
}
