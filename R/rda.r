# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Redundancy Analysis (RDA)
#'
#' Constrained ordination of samples from their feature abundances: the
#' axes of a (normalized) count table that sample metadata explain, followed
#' by principal components of what they leave unexplained. Gives the same
#' results as `vegan::rda()`, computed straight from sparse input without
#' forming any \eqn{n \times n}{n x n} or dense table, so it takes far less
#' time and memory for large datasets. With `norm = 'hellinger'` or
#' `norm = 'chord'` it is transformation-based RDA (tb-RDA); with
#' `Condition()` terms, partial RDA.
#'
#' @inherit pseudocount_section
#' @inheritParams pca
#' @inheritParams dbrda
#'
#' @param counts  A numeric matrix of feature abundances (samples
#'        \eqn{\times} features), or any of the table types listed under
#'        "Input Types". Alternatively, as in `vegan::rda()`, a formula with
#'        the table on its left-hand side and the constraints on its right,
#'        such as `counts ~ treatment`, in place of `counts` and `formula`;
#'        the metadata may then be given second, as in
#'        `rda(counts ~ treatment, meta)`.
#'
#' @param groups  The one variable that constrains the ordination: a vector
#'        or factor with one value per sample, matched to the samples in
#'        `counts` by name when both are named, else in the same order. Or,
#'        when `data` is given, the name of one of its columns (or a vector
#'        with one value per row of `data`). Categories (factor, character,
#'        or logical) or a numeric covariate. The same as
#'        `formula = ~ groups`. Give `groups` or `formula`, or neither for an
#'        unconstrained RDA, which is a PCA.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the columns named in `groups` or `formula`. Rows are
#'        matched to the samples in `counts` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as the samples. Default: `NULL`.
#'
#' @param formula  Instead of `groups`, for models of several variables: a
#'        one-sided formula naming the variables in `data` that constrain the
#'        ordination, such as `~ treatment + ph`. Factors, numeric
#'        covariates, and interactions are all allowed. Wrap variables in
#'        `Condition()` to remove their effect first (partial RDA), as in
#'        `vegan::rda()`: `~ treatment + Condition(site)`. With only
#'        `Condition()` terms, the result is the PCA of the table after
#'        removing their effect. A left-hand side, if given, is ignored.
#'        Default: `NULL`.
#'
#' @param scale   Scale each feature to unit variance (after `norm`, before
#'        the analysis), as in `vegan::rda(scale = TRUE)` and [pca()]'s
#'        `scale`, so that every feature has equal weight. Features with
#'        constant values are left out. Default: `FALSE`
#'
#' @param n_axes   The number of axes to return. Default: `3`. Reduced, with a
#'        warning, to the number of axes with positive eigenvalues, which is
#'        at most \eqn{n - 1} for \eqn{n} samples.
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   as described in "Ordination result" below. Constrained axes come first,
#'   named `"RDA1"`, `"RDA2"`, etc. If the constraints give fewer than
#'   `n_axes` of them, unconstrained axes follow, named `"PC1"`, `"PC2"`,
#'   etc. The values are `vegan`'s site scores ("wa" scores, scaling 2):
#'   what `vegan::scores(x, display = "sites")` returns and what
#'   `phyloseq::plot_ordination()` plots. Its attributes are those of
#'   [dbrda()], with the same meanings, plus `loadings`:
#'
#'   * `method`: `"rda"`.
#'   * `eig`: The variance of each axis.
#'   * `prop`: The proportion of the variance that each axis explains, not
#'     counting the conditional variance:
#'     `eig / (trace - inertia[["conditional"]])`. `vegan`'s "Proportion
#'     Explained".
#'   * `trace`: The total variance: the sum of all features' variances
#'     (after `norm` and `scale`), as [pca()]'s `trace`.
#'   * `inertia`: The total variance split into three parts, which add up to
#'     `trace`: a named vector of
#'     - `conditional`: explained by the `Condition()` variables (`0` if
#'       there are none);
#'     - `constrained`: explained by the constraints, after the conditions
#'       (the sum of all constrained eigenvalues);
#'     - `unconstrained`: the rest.
#'   * `constrained_eig`: Every non-zero constrained eigenvalue.
#'   * `r_squared`: The proportion of the total variance that the
#'     constraints explain, `inertia[["constrained"]] / trace`, as
#'     `vegan::RsquareAdj(x)$r.squared`.
#'   * `adj_r_squared`: `r_squared` adjusted for the number of constraint
#'     columns, as `vegan::RsquareAdj(x)$adj.r.squared`; see [dbrda()] for
#'     the formula, with and without conditions. `NA` when the model leaves
#'     no residual degrees of freedom.
#'   * `biplot`: Biplot scores of the constraints: one row per column of the
#'     model matrix and one column per axis (`0` for unconstrained axes).
#'     `vegan::scores(x, display = "bp")`.
#'   * `centroids`: The centroid of each factor level on each axis, as
#'     `vegan::scores(x, display = "cn")`, or `NULL` if no variable is a
#'     factor.
#'   * `loadings`: The species scores: one row per feature and one column
#'     per axis, `vegan::scores(x, display = "species")` (scaling 2). Unlike
#'     [pca()]'s `loadings`, these are not unit-length: each axis's are
#'     scaled by the square root of its share of the variance, so that the
#'     sites and species form a correlation biplot, as in `vegan`.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{Z} the normalized table (samples in rows), \eqn{J} the
#' centering matrix, and \eqn{D} the diagonal matrix of the features'
#' standard deviations (or the identity when `scale = FALSE`), RDA (Rao
#' 1964; Legendre & Legendre 2012, §11.1) analyzes the centered table
#' \eqn{Y = JZD^{-1} / \sqrt{n - 1}}. With \eqn{H} the projection onto the
#' (centered) constraints, after removing any conditions (Borcard et al.
#' 1992), the constrained axes are the principal components of the fitted
#' values \eqn{HY}, and the unconstrained axes those of the residuals
#' \eqn{(I - H)Y}. Each axis's eigenvalue is a variance, so they add up to
#' the total variance of \eqn{Y}.
#'
#' `rda()` shares its solver with [dbrda()]: since \eqn{YY^T} is Gower's
#' centered matrix of the Euclidean distances between the rows of \eqn{Y},
#' the constrained axes come from a small eigenproblem with one row per
#' constraint column, and the unconstrained axes, when needed, from block
#' Lanczos iteration. Both need only products of the table with a few
#' vectors at a time, computed in parallel from the non-zero values alone,
#' as in [pca()], so a sparse table is never made dense.
#'
#'
#' **Transformation-based RDA**
#'
#' RDA of raw counts measures differences between samples by Euclidean
#' distances, which suit species abundance data poorly (Legendre & Gallagher
#' 2001). The Hellinger and chord transformations fix that:
#' `rda(counts, ..., norm = 'hellinger')` is an RDA that preserves Hellinger
#' distances between samples, and `norm = 'chord'` chord distances. `norm =
#' 'clr'` (with a `pseudocount`) preserves Aitchison distances.
#'
#'
#' **Relation to dbrda(), pca() and permanova()**
#'
#' * `rda(counts, ..., norm = 'hellinger')` has the same axes, `prop`,
#'   \eqn{R^2} and biplot scores as
#'   `dbrda(euclidean(counts, norm = 'hellinger'), ...)`: db-RDA of
#'   Euclidean distances is RDA. `dbrda()` analyzes \eqn{YY^T} without
#'   dividing by \eqn{n - 1}, so its `trace`, `eig`, `inertia`, and
#'   `constrained_eig` are `rda()`'s times \eqn{n - 1}, and its scores
#'   `rda()`'s times \eqn{(n - 1)^{1/4}}. When any distance is 4 or more
#'   (see [dbrda()]'s "Scale"), as Euclidean distances between raw counts
#'   often are, `dbrda()` divides by \eqn{n - 1} too, and all of these are
#'   the same as `rda()`'s.
#' * With no constraints, `rda(counts)` is a PCA: its `eig`, `prop`, and
#'   `trace` are those of [pca()] with the same `norm` and `scale`. Its
#'   scores point the same way as `pca()`'s, but are in `vegan`'s scaling
#'   2: each axis's scores are `pca()`'s divided by
#'   \eqn{\sqrt{(n - 1) \cdot eig}} and multiplied by
#'   \eqn{\sqrt[4]{(n - 1) \cdot trace}}.
#' * The permutation test of an RDA (`vegan::anova.cca()`) is PERMANOVA of
#'   the Euclidean distances between the normalized samples: e.g.
#'   `permanova(euclidean(counts, norm = 'hellinger'), 'treatment', meta)`.
#'   Its \eqn{R^2} is `rda()`'s `r_squared`, and its pseudo-F and
#'   p-value those of `anova(vegan::rda(...))` for the same permutations.
#'
#' `varpart()`, `forward_select()` and `principal_response_curves()`, which
#' are planned for ecodive, will build on `rda()`.
#'
#'
#' **Matching vegan**
#'
#' `rda()` gives the same scores, eigenvalues, inertia, \eqn{R^2}, biplot
#' scores, centroids and species scores as `vegan::rda()` on the normalized
#' table, up to the signs of the axes. Columns of the model matrix that are
#' linear combinations of earlier ones are dropped, with a message, as in
#' `vegan`. With `scale = TRUE`, features with constant values are left out
#' (their loadings are 0), where `vegan` would stop with an error.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::rda()` shares its name with ecodive's `rda()`, and whichever
#' package was attached last masks the other's. Write `ecodive::rda()`, or
#' use [count_ord()], which always runs ecodive's, to be sure.
#'
#'
#' **Axis signs and reproducibility**
#'
#' An eigenvector's sign is arbitrary, so different software can return
#' mirror images of the same axis. `rda()` flips each axis so that its
#' largest score in absolute value is positive, and flips the biplot scores,
#' centroids and loadings with it. Results are identical from run to run
#' and for any number of `cpus`.
#'
#' @references
#' Rao, C. R. (1964). The use and interpretation of principal component
#' analysis in applied research. *Sankhya A*, 26(4), 329-358.
#'
#' Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
#' transformations for ordination of species data. *Oecologia*, 129(2),
#' 271-280. \doi{10.1007/s004420100716}
#'
#' Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the
#' spatial component of ecological variation. *Ecology*, 73(3), 1045-1055.
#' \doi{10.2307/1940179}
#'
#' Legendre, P., Oksanen, J., & ter Braak, C. J. F. (2011). Testing the
#' significance of canonical axes in redundancy analysis. *Methods in
#' Ecology and Evolution*, 2(3), 269-277.
#' \doi{10.1111/j.2041-210X.2010.00078.x}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.).
#' Elsevier.
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [dbrda()], [capscale()], [pca()], [permanova()], [count_ord()]
#'
#' @export
#' @examples
#'     # tb-RDA: RDA of Hellinger-transformed abundances.
#'     counts <- round(iris[, 1:4] * 10)
#'     ord    <- rda(counts, 'Species', iris, norm = 'hellinger')
#'     plot(ord, 'Species', iris)
#'     attr(ord, 'r_squared')
#'     attr(ord, 'loadings')
#'
#'     # Chord-transformed, in vegan's style.
#'     ord <- rda(counts ~ Species, iris, norm = 'chord', n_axes = 2)
#'     str(ord)
#'
#'     # Its permutation test: PERMANOVA of the Hellinger distances.
#'     permanova(euclidean(counts, norm = 'hellinger'), 'Species', iris, permutations = 99)
#'
#'     # Without constraints, a PCA.
#'     rda(counts)
#'
rda <- function (
    counts, groups, data = NULL, formula = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, scale = FALSE, n_axes = 3L, cpus = n_cpus() ) {

  env       <- parent.frame()
  no_groups <- missing(groups)
  if (no_groups) groups <- NULL

  # vegan's style: the table on the formula's left-hand side.
  if (inherits(counts, 'formula') && length(counts) == 3) {
    if (!no_groups && is.null(data) && is.data.frame(groups)) {
      data      <- groups
      groups    <- NULL
      no_groups <- TRUE
    }
    ord_formula_check('rda', 'counts', groups, data)
    if (!is.null(formula))
      stop(
        'Give either a formula with the table on its left-hand side, ',
        'or `counts` and `formula`, not both.', call. = FALSE )
    fd      <- formula_dm(counts, NULL)
    counts  <- fd$dm
    formula <- fd$formula
  }
  else {
    ord_formula_check('rda', 'counts', groups, data)
  }

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

  model <- ord_model(groups, no_groups, data, formula, n, dn[[1]], env, optional = TRUE)
  d     <- cord_design(model$x, model$z, n)


  if (n < 2) {
    res          <- cord_empty(n, d)
    res$loadings <- matrix(numeric(0), p, 0)
  }
  else {
    res <- .Call(
      C_rda, counts, margin, norm, pseudocount, scale,
      d$design, ncol(d$z), k, cpus )

    if (!res$converged) # nocov start
      warning('RDA did not converge; coordinates may be inexact.', call. = FALSE)
    # nocov end
  }

  # vegan's species scores, scaling 2: v * sqrt(eig / trace) * const, with
  # v = Y' u / sqrt(eig) and const = ((n - 1) * trace)^(1/4), where Y is the
  # table divided by sqrt(n - 1). That is A' u * const / sqrt((n - 1) *
  # trace), for C_rda()'s loadings A' u.
  loadings <- res$loadings
  if (ncol(loadings) > 0)
    loadings <- loadings * (sqrt(sqrt((n - 1) * res$total)) / sqrt((n - 1) * res$total))

  ord <- cord_result(
    res, d, model$mf, k, n, dn[[1]],
    method = 'rda', prefix = 'RDA', prefix_u = 'PC',
    loadings = loadings )

  dimnames(attr(ord, 'loadings')) <- list(dn[[2]], colnames(ord))

  ord
}
