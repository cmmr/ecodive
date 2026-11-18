# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Canonical Correspondence Analysis (CCA)
#'
#' Constrained ordination of samples from their feature abundances in the
#' chi-square metric: the axes of correspondence analysis that sample
#' metadata explain, followed by correspondence analysis axes of what they
#' leave unexplained (ter Braak 1986). Gives the same results as
#' `vegan::cca()`, computed straight from sparse input without forming any
#' \eqn{n \times n}{n x n} or dense table, so it takes far less time and
#' memory for large datasets. With `Condition()` terms, partial CCA.
#'
#' This is canonical **correspondence** analysis, a unimodal method for
#' species abundances. It is not canonical *correlation* analysis, which
#' scikit-learn, MultivariateStats.jl and linfa also call "CCA".
#'
#' @inheritParams rda
#' @inheritParams permanova
#'
#' @param counts  A numeric matrix of non-negative feature abundances
#'        (samples \eqn{\times} features), or any of the table types listed
#'        under "Input Types". Every sample must have at least one non-zero
#'        abundance. Alternatively, as in `vegan::cca()`, a formula with the
#'        table on its left-hand side and the constraints on its right, such
#'        as `counts ~ treatment`, in place of `counts` and `formula`; the
#'        metadata may then be given second, as in
#'        `cca(counts ~ treatment, meta)`.
#'
#' @param groups  The one variable that constrains the ordination: a vector
#'        or factor with one value per sample, matched to the samples in
#'        `counts` by name when both are named, else in the same order. Or,
#'        when `data` is given, the name of one of its columns (or a vector
#'        with one value per row of `data`). Categories (factor, character,
#'        or logical) or a numeric covariate. The same as
#'        `formula = ~ groups`. Give `groups` or `formula`, or neither for an
#'        unconstrained correspondence analysis.
#'
#' @param formula  Instead of `groups`, for models of several variables: a
#'        one-sided formula naming the variables in `data` that constrain the
#'        ordination, such as `~ treatment + ph`. Factors, numeric
#'        covariates, and interactions are all allowed. Wrap variables in
#'        `Condition()` to remove their effect first (partial CCA), as in
#'        `vegan::cca()`: `~ treatment + Condition(site)`. With only
#'        `Condition()` terms, the result is the correspondence analysis of
#'        the table after removing their effect. A left-hand side, if given,
#'        is ignored. Default: `NULL`.
#'
#' @param permutations  The number of permutations for `adj_r_squared`.
#'        Default: `1000`, as `vegan::RsquareAdj()`. Alternatively, a matrix
#'        with one permutation of `1:n` per row, as from
#'        `permute::shuffleSet()`, to use instead of random ones. `0` gives
#'        an `adj_r_squared` of `NA`. Unused for partial CCA and
#'        unconstrained CA, whose `adj_r_squared` is always `NA`.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same `adj_r_squared`. Default: `0`
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   as described in "Ordination result" below. Constrained axes come first,
#'   named `"CCA1"`, `"CCA2"`, etc. If the constraints give fewer than
#'   `n_axes` of them, unconstrained axes follow, named `"CA1"`, `"CA2"`,
#'   etc. The values are `vegan`'s site scores ("wa" scores, scaling 2):
#'   what `vegan::scores(x, display = "sites")` returns and what
#'   `phyloseq::plot_ordination()` plots. Its attributes are those of
#'   [rda()], with the same meanings, but for inertia in place of variance:
#'
#'   * `method`: `"cca"`.
#'   * `eig`: The eigenvalue (inertia) of each axis.
#'   * `prop`: The proportion of the inertia that each axis explains, not
#'     counting the conditional inertia:
#'     `eig / (trace - inertia[["conditional"]])`. `vegan`'s "Proportion
#'     Explained".
#'   * `trace`: The total inertia: Pearson's chi-square statistic of the
#'     table divided by its grand total (`vegan`'s `tot.chi`).
#'   * `inertia`: The total inertia split into three parts, which add up to
#'     `trace`: a named vector of
#'     - `conditional`: explained by the `Condition()` variables (`0` if
#'       there are none);
#'     - `constrained`: explained by the constraints, after the conditions
#'       (the sum of all constrained eigenvalues);
#'     - `unconstrained`: the rest.
#'   * `constrained_eig`: Every non-zero constrained eigenvalue.
#'   * `r_squared`: The proportion of the total inertia that the
#'     constraints explain, `inertia[["constrained"]] / trace`, as
#'     `vegan::RsquareAdj(x)$r.squared`.
#'   * `adj_r_squared`: `r_squared` adjusted by permutation, as
#'     `vegan::RsquareAdj(x)$adj.r.squared` (see "Adjusted R-squared"
#'     below). `NA` for partial CCA, unconstrained CA, or
#'     `permutations = 0`.
#'   * `biplot`: Biplot scores of the constraints: their weighted
#'     correlations with the constrained axes' "lc" scores, one row per
#'     column of the model matrix and one column per axis (`0` for
#'     unconstrained axes). `vegan::scores(x, display = "bp")`.
#'   * `centroids`: The weighted centroid of each factor level on each
#'     axis, as `vegan::scores(x, display = "cn")`, or `NULL` if no
#'     variable is a factor.
#'   * `loadings`: The species scores: one row per feature and one column
#'     per axis, `vegan::scores(x, display = "species")` (scaling 2).
#'     Features with no observations get `0`.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{P} the table divided by its grand total, and \eqn{r} and
#' \eqn{c} its row and column sums (the sample and feature weights), CCA
#' (ter Braak 1986; Legendre & Legendre 2012, §11.2) analyzes the matrix of
#' chi-square residuals \eqn{\bar{Q} = D_r^{-1/2} (P - rc^T) D_c^{-1/2}},
#' as correspondence analysis does (§9.2). The constraints are centered and
#' regressed with the weights \eqn{r}. With \eqn{H} the projection onto the
#' weighted constraints, after removing any conditions, the constrained
#' axes are the principal components of \eqn{H\bar{Q}}, and the
#' unconstrained axes those of \eqn{(I - H)\bar{Q}}. Site scores are
#' divided by \eqn{\sqrt{r}}, and species scores by \eqn{\sqrt{c}}.
#'
#' `cca()` shares its solver with [rda()] and [dbrda()]: the constrained
#' axes come from a small eigenproblem with one row per constraint column,
#' and the unconstrained axes, when needed, from block Lanczos iteration.
#' Both need only products of \eqn{\bar{Q}} with a few vectors at a time.
#' \eqn{\bar{Q}} is a sparse matrix less a rank-one term, so these are
#' computed in parallel from the table's non-zero values alone, and a
#' sparse table is never made dense.
#'
#'
#' **CCA or RDA?**
#'
#' CCA and correspondence analysis assume that each species responds
#' unimodally to the gradients, peaking at some optimum; RDA and PCA assume
#' linear responses. Over long gradients, where many species are absent
#' from much of the range, unimodal methods suit better; over short ones,
#' linear methods on transformed abundances (tb-RDA:
#' `rda(counts, ..., norm = 'hellinger')`) are often preferred (Legendre &
#' Gallagher 2001). The chi-square metric gives rare species a large
#' weight, so CCA can be sensitive to them.
#'
#'
#' **Adjusted R-squared**
#'
#' Ezekiel's adjustment, which [rda()] uses, is biased for CCA (Peres-Neto
#' et al. 2006). So, as `vegan::RsquareAdj()` does, `cca()` estimates it by
#' permutation: \eqn{R^2_{adj} = 1 - (1 - R^2) / (1 - \bar{R}^2_{perm})},
#' where \eqn{\bar{R}^2_{perm}} is the mean \eqn{R^2} over `permutations`
#' random permutations of the samples (their rows of the table, against the
#' constraints). It is random: the same `seed` gives the same value, and
#' a `permutations` matrix the same value as
#' `vegan::RsquareAdj(x, permutations = matrix)`. As in `vegan`, it is `NA`
#' for partial CCA.
#'
#'
#' **Matching vegan**
#'
#' `cca()` gives the same scores, eigenvalues, inertia, \eqn{R^2}, biplot
#' scores, centroids and species scores as `vegan::cca()`, up to the signs
#' of the axes. Columns of the model matrix that are linear combinations of
#' earlier ones are dropped, with a message, as in `vegan`. Samples with no
#' observations are an error, as in `vegan`; features with no observations
#' have no weight, and get species scores of 0.
#'
#' With no constraints, `cca(counts)` is a correspondence analysis, as
#' `vegan::cca(counts)`.
#'
#'
#' **Permutation tests**
#'
#' `cca()` has no permutation test of its own; `vegan::anova.cca()` gives
#' one. (Unlike RDA's, CCA's test is not a PERMANOVA of any distance.)
#'
#'
#' **Same names in other packages**
#'
#' `vegan::cca()` shares its name with ecodive's `cca()`, and whichever
#' package was attached last masks the other's. Write `ecodive::cca()`, or
#' use [count_ord()], which always runs ecodive's, to be sure.
#'
#'
#' **Axis signs and reproducibility**
#'
#' An eigenvector's sign is arbitrary, so different software can return
#' mirror images of the same axis. `cca()` flips each axis so that its
#' largest score in absolute value is positive, and flips the biplot scores,
#' centroids and loadings with it. Results are identical from run to run
#' and for any number of `cpus`.
#'
#' @references
#' ter Braak, C. J. F. (1986). Canonical correspondence analysis: a new
#' eigenvector technique for multivariate direct gradient analysis.
#' *Ecology*, 67(5), 1167-1179. \doi{10.2307/1938672}
#'
#' Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
#' transformations for ordination of species data. *Oecologia*, 129(2),
#' 271-280. \doi{10.1007/s004420100716}
#'
#' Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006). Variation
#' partitioning of species data matrices: estimation and comparison of
#' fractions. *Ecology*, 87(10), 2614-2625.
#' \doi{10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.).
#' Elsevier.
#'
#' Greenacre, M. (2007). *Correspondence Analysis in Practice* (2nd ed.).
#' Chapman & Hall/CRC.
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [rda()], [dbrda()], [pca()], [count_ord()]
#'
#' @export
#' @examples
#'     counts <- round(iris[, 1:4] * 10)
#'     ord    <- cca(counts, 'Species', iris)
#'     plot(ord, col = iris$Species)
#'     attr(ord, 'r_squared')
#'     attr(ord, 'loadings')
#'
#'     # vegan's style, with a numeric constraint.
#'     ord <- cca(counts ~ Petal.Width, iris, n_axes = 2)
#'     str(ord)
#'
#'     # Without constraints, a correspondence analysis.
#'     head(cca(counts))
#'
cca <- function (
    counts, groups, data = NULL, formula = NULL, margin = 1L, n_axes = 3L,
    permutations = 1000L, seed = 0L, cpus = n_cpus() ) {

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
    ord_formula_check('cca', 'counts', groups, data)
    if (!is.null(formula))
      stop(
        'Give either a formula with the table on its left-hand side, ',
        'or `counts` and `formula`, not both.', call. = FALSE )
    fd      <- formula_dm(counts, NULL)
    counts  <- fd$dm
    formula <- fd$formula
  }
  else {
    ord_formula_check('cca', 'counts', groups, data)
  }

  validate_counts()
  validate_margin()
  validate_seed()
  validate_cpus()
  k <- ord_k(n_axes)


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
  pm    <- perm_matrix(permutations, n)

  # vegan gives no adjusted R-squared for partial models.
  n_perm <- if (ncol(d$z) > 0 || ncol(d$x) == 0) 0L else pm$n_perm


  if (n < 2) {
    if (n == 1 && p > 0) cca_empty_check(counts, margin, dn[[1]])
    res          <- cord_empty(n, d)
    res$loadings <- matrix(numeric(0), p, 0)
    res$w        <- rep(1, n)
    res$perm_inertia <- numeric(0)
  }
  else {
    res <- .Call(
      C_cca, counts, margin, d$design, ncol(d$z), k,
      if (n_perm > 0) pm$perms, n_perm, seed, cpus )

    if (!is.null(res$empty))
      stop(cca_empty_msg(res$empty, dn[[1]]), call. = FALSE)

    if (!res$converged) # nocov start
      warning('CCA did not converge; coordinates may be inexact.', call. = FALSE)
    # nocov end
  }

  r2 <- cca_r_squared(res, n_perm, ncol(d$x) > 0)

  ord <- cord_result(
    res, d, model$mf, k, n, dn[[1]],
    method = 'cca', prefix = 'CCA', prefix_u = 'CA',
    scale = 1, r2 = r2, weights = res$w,
    loadings = res$loadings )

  dimnames(attr(ord, 'loadings')) <- list(dn[[2]], colnames(ord))

  ord
}



# R-squared, and vegan's permutation-based adjusted R-squared for cca
# (RsquareAdj.cca): 1 - (1 - R2) / (1 - mean permuted R2). Without
# constraints, both are 0, as rda()'s; NA without permutations (partial
# models, or permutations = 0).
cca_r_squared <- function (res, n_perm, constrained) {

  total <- res$total
  if (!(total > 0)) return (c(NA_real_, NA_real_))

  r2 <- res$constrained / total
  if (!constrained) return (c(r2, 0))
  if (n_perm == 0)  return (c(r2, NA_real_))

  # Saturated models (no residual inertia) have no adjustment.
  mp <- mean(res$perm_inertia) / total
  if (!(1 - mp > sqrt(.Machine$double.eps))) return (c(r2, NA_real_))

  c(r2, 1 - (1 - r2) / (1 - mp))
}



cca_empty_msg <- function (idx, labels) {
  which_ <- if (is.null(labels)) idx else labels[idx]
  if (length(which_) > 5) which_ <- c(which_[1:5], '...')
  paste0(
    '`counts` has ', length(idx), ' sample', if (length(idx) > 1) 's',
    ' with no observations (all zero), which CCA cannot place: ',
    paste(which_, collapse = ', '), '. Remove ',
    if (length(idx) > 1) 'them' else 'it', ' first.' )
}


# For one sample (n = 1, so no C call): is it empty?
cca_empty_check <- function (counts, margin, labels) {
  tot <- switch(
    EXPR = get_matrix_package(counts),
    'base'   = sum(counts),
    'slam'   = sum(counts$v),
    'Matrix' = sum(counts) )
  if (!isTRUE(tot > 0)) stop(cca_empty_msg(1L, labels), call. = FALSE)
}
