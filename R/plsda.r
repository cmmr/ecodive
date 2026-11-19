# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Partial Least Squares Discriminant Analysis (PLS-DA)
#'
#' Supervised ordination of samples from their feature abundances: the
#' components of a (normalized) count table that best separate known groups.
#' PLS-DA is partial least squares (PLS) regression of a dummy matrix of the
#' groups on the table. Gives the same results as scikit-learn's
#' `PLSRegression` and `mixOmics::plsda()`, computed straight from sparse
#' input without forming any dense table or \eqn{n \times n}{n x n} matrix.
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of feature abundances (samples
#'        \eqn{\times} features), or any of the table types listed under
#'        "Input Types".
#'
#' @param groups  The known group of each sample: a vector or factor with
#'        one value per sample, matched to the samples in `counts` by name
#'        when both are named, else in the same order. Or, when `data` is
#'        given, the name of one of its columns (or a vector with one value
#'        per row of `data`). Every distinct value is a group, numbers
#'        included; there must be at least two. A factor's unused levels are
#'        dropped.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the column named in `groups`. Rows are matched to the
#'        samples in `counts` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as the
#'        samples. Default: `NULL`.
#'
#' @param scale   Scale each feature to unit variance (after `norm`, before
#'        the analysis), and each column of the dummy matrix too, as
#'        scikit-learn's `PLSRegression(scale = True)` and
#'        `mixOmics::plsda(scale = TRUE)` do, so that every feature has equal
#'        weight. Features with constant values are left out: their
#'        weights, loadings and VIP are 0. Default: `TRUE`
#'
#' @param n_axes   The number of components to return. Default: `2`.
#'        Reduced to \eqn{n - 1} for \eqn{n} samples, or to the number of
#'        features, if fewer; and, with a warning, to the number of
#'        components that the table can fit (see "Method").
#'
#' @return A numeric matrix with one row per sample and one column per
#'   component, named `"PLSDA1"`, `"PLSDA2"`, etc., as described in
#'   "Ordination result" below. The values are the X scores \eqn{t}
#'   (scikit-learn's `x_scores_`, `mixOmics`' `variates$X`). Attributes:
#'
#'   * `method`: `"plsda"`.
#'   * `eig`: The variance of the table that each component explains,
#'     \eqn{|t|^2 |p|^2 / (n - 1)} for scores \eqn{t} and loadings \eqn{p}.
#'     (PLS components are not eigenvectors; this is the analogue of
#'     [pca()]'s eigenvalues, and it adds up over components.)
#'   * `prop`: The proportion of the table's total variance that each
#'     component explains, `eig / trace`: `mixOmics`' `prop_expl_var$X`.
#'   * `trace`: The total variance of the table: the sum of all features'
#'     variances (after `norm` and `scale`), as [pca()]'s `trace`.
#'   * `loadings`: The X loadings \eqn{p}: one row per feature and one
#'     column per component (scikit-learn's `x_loadings_`). Each column
#'     regresses the (deflated) table on the component's scores.
#'   * `weights`: The X weights \eqn{w}, unit-length columns (scikit-learn's
#'     `x_weights_`; what `mixOmics` calls `loadings$X`).
#'   * `y_loadings`: The Y loadings \eqn{q}: one row per group and one
#'     column per component (scikit-learn's `y_loadings_`), for the dummy
#'     matrix as centred (and scaled).
#'   * `vip`: Variable importance in projection (VIP): one row per feature, and
#'     one column per number of components, the column for \eqn{a}
#'     components counting the first \eqn{a} (as `mixOmics::vip()`). See
#'     "Method".
#'   * `groups`: Each sample's group, as a factor whose levels are the
#'     groups, in the order of the rows of `y_loadings` and `centroids`.
#'   * `centroids`: The mean score of each group: one row per group and one
#'     column per component.
#'
#' @details
#'
#' **Method**
#'
#' Let \eqn{X} be the normalized table, its features centred (and, with
#' `scale = TRUE`, scaled to unit variance), and \eqn{Y} the dummy matrix of
#' the groups (one column per group, 1 for that group's samples and 0 for
#' the others), centred (and scaled) the same way. PLS regression (Wold et
#' al. 2001) finds one component at a time. With \eqn{X_a} and \eqn{Y_a}
#' what is left of \eqn{X} and \eqn{Y} before component \eqn{a}:
#'
#' * the weights \eqn{w} are the top left singular vector of
#'   \eqn{X_a^T Y_a}, the direction in which \eqn{X_a} covaries most with
#'   \eqn{Y_a};
#' * the scores are \eqn{t = X_a w};
#' * the loadings are \eqn{p = X_a^T t / |t|^2} and the Y loadings
#'   \eqn{q = Y_a^T t / |t|^2};
#' * both tables are deflated by the scores: \eqn{X_{a+1} = X_a - t p^T}
#'   and \eqn{Y_{a+1} = Y_a - t q^T} (regression mode).
#'
#' Barker & Rayens (2003) give the link to discriminant analysis. Fitting
#' stops early, with a warning, when \eqn{X_a^T Y_a} falls to
#' \eqn{\sqrt{\epsilon} \approx 1.5 \times 10^{-8}}{sqrt(eps) = 1.5e-8} times
#' its value for the first component: then what is left of the table has
#' no covariance with the groups, beyond rounding error (for instance, once
#' the components span the whole table).
#'
#' VIP (Wold et al. 2001; Mehmood et al. 2012) sums each feature's squared
#' weights over the first \eqn{a} components, each weighted by the
#' component's share of the explained sum of squares of \eqn{Y},
#' \eqn{SSY_b = |t_b|^2 |q_b|^2}:
#' \deqn{VIP_j = \sqrt{P \sum_{b \le a} SSY_b w_{jb}^2 / \sum_{b \le a} SSY_b}}{VIP_j = sqrt(P * sum(SSY_b * w_jb^2) / sum(SSY_b))}
#' for \eqn{P} features (constant ones included, when `scale` leaves them
#' out). The mean of the squared VIPs is 1, so features with a VIP above 1
#' are often called important.
#'
#' The table is never centred in memory: products with \eqn{X_a} are
#' computed from the non-zero values alone, in parallel, as in [pca()], so
#' a sparse table is never made dense.
#'
#'
#' **Choosing `norm`**
#'
#' For microbiome counts, `norm = 'clr'` (with a `pseudocount`) is the usual
#' choice, as in `mixOmics`' `logratio = 'CLR'`: PLS-DA then works on
#' log-ratios, which suit compositional data.
#'
#'
#' **Validate the separation**
#'
#' PLS-DA is supervised: it looks for the directions that separate the
#' groups, so its plots separate them even when the groups don't differ
#' (with more features than samples, almost always). Don't read the plot as
#' evidence of a difference. Test it with a permutation test of the
#' distances, such as [permanova()], or estimate the error rate of
#' classifying new samples by cross-validation.
#'
#'
#' **Matching scikit-learn and mixOmics**
#'
#' `plsda()` gives the same scores, weights, loadings and Y loadings as
#' scikit-learn's `PLSRegression(n_components, scale = True)` fitted to the
#' normalized table and the dummy matrix, up to the signs of the
#' components, and the same scores, weights (`loadings$X`), explained
#' variance (`prop_expl_var$X`) and VIP (`mixOmics::vip()`) as
#' `mixOmics::plsda()`. (With `scale = FALSE`, `mixOmics::vip()` weights each
#' component by its squared correlations with the dummy columns rather than
#' by \eqn{SSY}, which is the same only when they are scaled; `plsda()`
#' keeps the textbook formula.) Both fit each component by NIPALS, a power iteration that converges to
#' the weights above; `plsda()` computes them directly instead, as
#' scikit-learn's `algorithm = "svd"` does, so it matches them when they
#' are run with a small tolerance (`tol`). With their default tolerances,
#' their results with three groups or more can differ from the converged
#' ones by about 0.1%. With two groups, NIPALS converges at once.
#'
#'
#' **Same names in other packages**
#'
#' `mixOmics::plsda()` and `caret::plsda()` share their names with
#' ecodive's `plsda()`, and whichever package was attached last masks the
#' other's. Write `ecodive::plsda()`, or use [count_ord()], which always runs
#' ecodive's, to be sure.
#'
#'
#' **Component signs and reproducibility**
#'
#' A component's sign is arbitrary, so different software can return mirror
#' images of the same component. `plsda()` flips each component so that its
#' largest score in absolute value is positive, and flips its weights,
#' loadings and Y loadings with it. Results are identical from run to run
#' and for any number of `cpus`.
#'
#' @references
#' Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic
#' tool of chemometrics. *Chemometrics and Intelligent Laboratory Systems*,
#' 58(2), 109-130. \doi{10.1016/S0169-7439(01)00155-1}
#'
#' Barker, M., & Rayens, W. (2003). Partial least squares for
#' discrimination. *Journal of Chemometrics*, 17(3), 166-173.
#' \doi{10.1002/cem.785}
#'
#' Mehmood, T., Liland, K. H., Snipen, L., & Sæbø, S. (2012). A review of
#' variable selection methods in partial least squares regression.
#' *Chemometrics and Intelligent Laboratory Systems*, 118, 62-69.
#' \doi{10.1016/j.chemolab.2012.07.010}
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830.
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [pca()], [rda()], [permanova()], [count_ord()]
#'
#' @export
#' @examples
#'     # PLS-DA of iris, as if its measurements were counts.
#'     counts <- round(iris[, 1:4] * 10)
#'     ord    <- plsda(counts, 'Species', iris)
#'     plot(ord, col = iris$Species)
#'     attr(ord, 'prop')
#'     attr(ord, 'vip')
#'
#'     # CLR-transformed microbiome counts.
#'     grp <- rep(c('a', 'b'), length.out = nrow(ex_counts))
#'     ord <- plsda(ex_counts, grp, norm = 'clr', pseudocount = 1)
#'     attr(ord, 'centroids')
#'
#'     # Supervised: always validate the separation, e.g. with PERMANOVA.
#'     permanova(aitchison(ex_counts, pseudocount = 1), grp, permutations = 99)
#'
plsda <- function (
    counts, groups, data = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, scale = TRUE, n_axes = 2L, cpus = n_cpus() ) {

  if (missing(groups)) groups <- NULL

  validate_counts()
  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()
  k <- ord_k(n_axes)

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


  # The groups, as a factor, and their dummy matrix: centred and, with
  # `scale`, each column scaled to unit variance (scikit-learn's
  # _center_scale_xy()).
  tg  <- test_groups(groups, data, n, dn[[1]], ref = '`counts`')
  grp <- tg$data[[tg$col]]
  grp <- if (is.factor(grp)) droplevels(grp) else factor(grp)
  lev <- levels(grp)
  ng  <- length(lev)

  if (ng < 2)
    stop('`groups` must have at least two different values.', call. = FALSE)

  y <- outer(as.integer(grp), seq_len(ng), '==') + 0
  y <- sweep(y, 2L, colMeans(y))
  if (scale) y <- sweep(y, 2L, sqrt(colSums(y^2) / (n - 1)), '/')

  k_max <- min(k, n - 1L, p)

  res <- .Call(C_plsda, counts, margin, norm, pseudocount, scale, y, k_max, cpus)

  kk <- ncol(res$scores)
  if (kk < k_max)
    warning(sprintf(
      paste(
        'only %d of the first %d components could be fitted:',
        'the rest of the table has no covariance with `groups`' ),
      kk, k_max ), call. = FALSE )

  axes <- if (kk > 0) paste0('PLSDA', seq_len(kk)) else NULL


  # VIP for the first a components, a = 1..kk: each component weighted by
  # the sum of squares of Y it explains, |t|^2 |q|^2.
  w   <- res$weights
  ssy <- colSums(res$scores^2) * colSums(res$y_loadings^2)
  num <- w^2 * rep(ssy, each = p)
  for (a in seq_len(kk)[-1]) num[, a] <- num[, a - 1L] + num[, a]
  vip <- sqrt(p * sweep(num, 2L, cumsum(ssy), '/'))

  # The mean score of each group.
  centroids <- rowsum(res$scores, grp, reorder = TRUE) / as.vector(table(grp))

  feat <- list(dn[[2]], axes)
  dimnames(res$loadings)   <- feat
  dimnames(w)              <- feat
  dimnames(vip)            <- feat
  dimnames(res$y_loadings) <- list(lev, axes)
  dimnames(centroids)      <- list(lev, axes)
  if (!is.null(dn[[1]])) names(grp) <- dn[[1]]

  trace <- res$total / (n - 1)

  ord_result(
    res$scores, 'plsda', dn[[1]], prefix = 'PLSDA',
    eig        = res$eig,
    prop       = if (trace > 0) res$eig / trace else res$eig,
    trace      = trace,
    loadings   = res$loadings,
    weights    = w,
    y_loadings = res$y_loadings,
    vip        = vip,
    groups     = grp,
    centroids  = centroids )
}
