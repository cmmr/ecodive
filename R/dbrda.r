# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Distance-based Redundancy Analysis (db-RDA)
#'
#' Constrained ordination of samples from their distances: the axes of beta
#' diversity that sample metadata explain. Gives the same results as
#' `vegan::dbrda()`, without ever forming an \eqn{n \times n}{n x n} matrix,
#' so it takes far less time and memory for large datasets. [capscale()] is
#' the same analysis of the distances' positive eigenvalues only, as in
#' `vegan::capscale()`.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values. Alternatively, as in `vegan::dbrda()`, a formula
#'        with the distances on its left-hand side and the constraints on its
#'        right, such as `dm ~ treatment`, in place of `dm` and `formula`;
#'        the metadata may then be given second, as in
#'        `dbrda(dm ~ treatment, meta)`.
#'
#' @param groups  The one variable that constrains the ordination: a vector
#'        or factor with one value per sample, matched to the samples in `dm`
#'        by name when both are named, else in the same order as `dm`. Or,
#'        when `data` is given, the name of one of its columns (or a vector
#'        with one value per row of `data`). Categories (factor, character,
#'        or logical) or a numeric covariate. The same as
#'        `formula = ~ groups`. Give `groups` or `formula`, not both.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the columns named in `groups` or `formula`. Rows are
#'        matched to the samples in `dm` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as `dm`. Default: `NULL`.
#'
#' @param formula  Instead of `groups`, for models of several variables: a
#'        one-sided formula naming the variables in `data` that constrain the
#'        ordination, such as `~ treatment + ph`. Factors, numeric
#'        covariates, and interactions are all allowed. Wrap variables in
#'        `Condition()` to remove their effect first (a partial analysis), as
#'        in `vegan::dbrda()`: `~ treatment + Condition(site)`. With only
#'        `Condition()` terms, the result is the unconstrained ordination
#'        after removing their effect. A left-hand side, if given, is
#'        ignored. Give `groups` or `formula`: there is no default model.
#'        Default: `NULL`.
#'
#' @param n_axes   The number of axes to return. Default: `3`. Reduced, with a
#'        warning, to the number of axes with positive eigenvalues, which is
#'        at most \eqn{n - 1} for \eqn{n} samples.
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   as described in "Ordination result" below. Constrained axes come first,
#'   named `"dbRDA1"`, `"dbRDA2"`, etc. If the constraints give fewer than `n_axes`
#'   of them, unconstrained axes follow, named `"MDS1"`, `"MDS2"`, etc. The
#'   values are `vegan`'s site scores ("wa" scores, scaling 2): what
#'   `vegan::scores(x, display = "sites")` returns and what
#'   `phyloseq::plot_ordination()` plots. Its attributes are:
#'
#'   * `method`: `"dbrda"`.
#'   * `eig`: The eigenvalue of each axis.
#'   * `prop`: The proportion of the inertia that each axis explains, not
#'     counting the conditional inertia:
#'     `eig / (trace - inertia[["conditional"]])`. `vegan`'s "Proportion
#'     Explained".
#'   * `trace`: The total inertia: the sum of all eigenvalues of the
#'     unconstrained ordination, including negative ones. On the same scale
#'     as `eig`; see "Scale" below for how it relates to [pcoa()]'s.
#'   * `inertia`: The total inertia split into three parts, which add up to
#'     `trace`: a named vector of
#'     - `conditional`: explained by the `Condition()` variables (`0` if
#'       there are none);
#'     - `constrained`: explained by the constraints, after the conditions
#'       (the sum of all constrained eigenvalues);
#'     - `unconstrained`: the rest.
#'   * `constrained_eig`: Every non-zero constrained eigenvalue, including
#'     negative ones.
#'   * `r_squared`: The proportion of the total inertia that the
#'     constraints explain, `inertia[["constrained"]] / trace`: PERMANOVA's
#'     \eqn{R^2} for the same model ([permanova()]), and
#'     `vegan::RsquareAdj(x)$r.squared`.
#'   * `adj_r_squared`: `r_squared` adjusted for the number of constraint
#'     columns \eqn{p} (after dropping redundant ones), as
#'     `vegan::RsquareAdj(x)$adj.r.squared`: \eqn{1 - (1 - R^2)(n - 1) /
#'     (n - p - 1)}. With conditions, of \eqn{q} columns explaining a
#'     proportion \eqn{R^2_c} of `trace`, it is the adjusted \eqn{R^2} of
#'     constraints and conditions together less that of the conditions
#'     alone (Peres-Neto et al. 2006). `NA` when the model leaves no
#'     residual degrees of freedom (\eqn{n - p - q - 1 = 0}).
#'   * `biplot`: Biplot scores of the constraints: one row per column of the
#'     model matrix and one column per axis (`0` for unconstrained axes).
#'     `vegan::scores(x, display = "bp")`.
#'   * `centroids`: The centroid of each factor level on each axis, as
#'     `vegan::scores(x, display = "cn")`, or `NULL` if no variable is a
#'     factor.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{A} the matrix of squared distances and \eqn{J} the centering
#' matrix, db-RDA (McArdle & Anderson 2001) analyzes Gower's centered matrix
#' \eqn{G = -\frac{1}{2}JAJ}. With \eqn{H} the projection onto the (centered)
#' constraints, after removing any conditions, the constrained axes are the
#' eigenvectors of \eqn{HGH}, and the unconstrained axes are those of
#' \eqn{(I-H)G(I-H)}. Non-Euclidean distances, such as Bray-Curtis, can give
#' negative eigenvalues. They are kept in the inertia, as `vegan::dbrda()`
#' does, but have no axes.
#'
#' `vegan::dbrda()` forms \eqn{G} and \eqn{HGH} in full and computes all
#' \eqn{n} of their eigenvalues, which takes time proportional to \eqn{n^3}.
#' `dbrda()` instead projects \eqn{G} onto an orthonormal basis of the
#' constraints, which only needs a few products of \eqn{G} with vectors,
#' computed in parallel straight from `dm`, and a small eigenproblem with one
#' row per constraint column. Unconstrained axes are only computed when
#' needed, by block Lanczos iteration as in [pcoa()].
#'
#'
#' **dbrda() or capscale()?**
#'
#' [capscale()] analyzes \eqn{G_+}, the part of \eqn{G} with positive
#' eigenvalues, where `dbrda()` analyzes all of \eqn{G}. For Euclidean
#' distances (such as Euclidean, Hellinger, or Aitchison distances) \eqn{G}
#' has no negative eigenvalues, and the two give the same results, apart
#' from the axis names. For others, such as Bray-Curtis:
#'
#' * `dbrda()` analyzes the distances as they are. Its `r_squared` is
#'   PERMANOVA's, and it is fast for any number of samples. But some
#'   constrained eigenvalues can be negative, and these have no axes.
#' * `capscale()` analyzes the closest Euclidean configuration, the principal
#'   coordinates with positive eigenvalues, so every eigenvalue is a
#'   variance. But it ignores the negative eigenvalues' part of the
#'   distances, and it needs all of the eigenvalues of \eqn{G}, which takes
#'   time proportional to \eqn{n^3}.
#'
#' `phyloseq::ordinate(method = "CAP")` calls `vegan::capscale()`, so its
#' results match [capscale()]'s, not `dbrda()`'s.
#'
#'
#' **Matching vegan**
#'
#' `dbrda()` gives the same scores, eigenvalues, inertia, \eqn{R^2}, biplot
#' scores and centroids as `vegan::dbrda()`, up to the signs of the axes.
#' Columns of the model matrix that are linear combinations of earlier ones
#' are dropped, with a message, as in `vegan`.
#'
#'
#' **Scale**
#'
#' Like `vegan::dbrda()`, `dbrda()` divides the distances by
#' \eqn{\sqrt{n-1}} when any of them is 4 or more ("mean squared
#' distances"). That divides `trace`, `eig`, `inertia` and `constrained_eig`
#' by \eqn{n-1}, but leaves `prop` and \eqn{R^2} unchanged. So for the same
#' `dm`, `attr(dbrda(dm, ...), 'trace')` equals `attr(pcoa(dm), 'trace')`
#' when every distance is below 4 (e.g. Bray-Curtis or Jaccard), and equals
#' it divided by \eqn{n-1} otherwise (e.g. Euclidean distances between
#' counts).
#'
#'
#' **Axis signs and reproducibility**
#'
#' An eigenvector's sign is arbitrary, so different software can return
#' mirror images of the same axis. `dbrda()` flips each axis so that its
#' largest score in absolute value is positive, and flips the biplot scores
#' and centroids with it. Results are identical from run to run and for any
#' number of `cpus`.
#'
#' @references
#' McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
#' community data: a comment on distance-based redundancy analysis. *Ecology*,
#' 82(1), 290-297. \doi{10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2}
#'
#' Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy
#' analysis: testing multispecies responses in multifactorial ecological
#' experiments. *Ecological Monographs*, 69(1), 1-24.
#' \doi{10.1890/0012-9615(1999)069[0001:DBRATM]2.0.CO;2}
#'
#' Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006). Variation
#' partitioning of species data matrices: estimation and comparison of
#' fractions. *Ecology*, 87(10), 2614-2625.
#' \doi{10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [capscale()], [rda()], [pcoa()], [pca()], [permanova()],
#'   [dist_ord()], [beta_div()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- dbrda(dm, 'Species', iris)
#'     plot(ord, 'Species', iris)
#'     head(ord)
#'     str(ord)
#'
dbrda <- function (
    dm, groups, data = NULL, formula = NULL, n_axes = 3L, cpus = n_cpus() ) {

  validate_cpus()
  k <- ord_k(n_axes)

  no_groups <- missing(groups)
  constrained_ord(
    'dbrda', dm, if (no_groups) NULL else groups, no_groups,
    data, formula, k, cpus, env = parent.frame() )
}



# The front shared by dbrda() and capscale(), as permanova()'s: `groups`,
# or `data` and `formula`, or vegan's `dm ~ terms` (with the metadata second,
# where `groups` is). `no_groups`: whether the caller's `groups` was missing,
# so that an explicit NULL (such as a misspelled column, meta$tretment) is an
# error rather than no model. Then the design, and dbrda_core().
constrained_ord <- function (
    fn, dm, groups, no_groups, data, formula, k, cpus, env ) {

  lhs <- inherits(dm, 'formula') && length(dm) == 3
  if (lhs && !no_groups && is.null(data) && is.data.frame(groups)) {
    data      <- groups
    groups    <- NULL
    no_groups <- TRUE
  }

  ord_formula_check(fn, 'dm', groups, data)

  fd <- formula_dm(dm, formula)
  dv <- dist_vector(fd$dm)

  model <- ord_model(groups, no_groups, data, fd$formula, dv$n, dv$labels, env)

  dbrda_core(dv, model$x, model$z, model$mf, k, cpus, positive = (fn == 'capscale'))
}



# A formula given in place of `data` or `groups` is an error.
ord_formula_check <- function (fn, first, groups, data) {
  if (inherits(data, 'formula') || inherits(groups, 'formula'))
    stop(
      'Name the formula argument: ', fn, '(', first, ', data = meta, formula = ~ a + b), ',
      'or give the groups: ', fn, '(', first, ', groups, data).', call. = FALSE )
}



# The constraints (x), conditions (z) and the constraints' model frame (mf,
# for factor-level centroids) from `groups`, or `data` and `formula`, for n
# samples named `labels`. With `optional`, no model at all (no `groups`, no
# `formula`) gives list(x = NULL, z = NULL, mf = NULL), an unconstrained
# ordination, rather than an error.
ord_model <- function (
    groups, no_groups, data, formula, n, labels, env, optional = FALSE ) {

  if (!is.null(groups)) {

    if (!is.null(formula))
      stop('Give `groups` or `formula`, not both.', call. = FALSE)

    g       <- test_groups(groups, data, n, labels)
    data    <- g$data
    formula <- stats::as.formula(call('~', as.name(g$col)), env = env)
  }
  else if (is.null(formula)) {
    if (optional && no_groups) return (list(x = NULL, z = NULL, mf = NULL))
    stop(
      if (no_groups) '`groups` is required: one value per sample, or a column name of `data`. '
      else '`groups` is NULL (a misspelled column, such as meta$tretment?). ',
      'Or give `data` and a `formula`, such as `formula = ~ a + b`.', call. = FALSE )
  }
  else if (is.null(data)) {
    stop('`data` is required with `formula`: the sample metadata.', call. = FALSE)
  }

  design <- model_design(
    data, formula, labels = labels, n = n, allow_empty = TRUE, env = env )

  # The constraints' variables, for factor-level centroids (vegan leaves out
  # the conditions' variables).
  mf <- NULL
  if (!is.null(design$terms)) {
    mf <- stats::model.frame(design$terms, design$data, na.action = stats::na.fail)
    attr(mf, 'terms') <- NULL
  }

  list(x = design$x, z = design$z, mf = mf)
}



# The numeric design: constraints `x` and conditions `z` as matrices with n
# rows (zero columns if NULL), and both together, conditions first.
cord_design <- function (x, z, n) {

  if (is.null(x)) x <- matrix(0, n, 0)
  if (is.null(z)) z <- matrix(0, n, 0)
  x <- as.matrix(x)
  z <- as.matrix(z)
  if (is.null(colnames(x)) && ncol(x)) colnames(x) <- paste0('X', seq_len(ncol(x)))

  design <- cbind(z, x)
  storage.mode(design) <- 'double'

  list(x = x, z = z, design = design)
}



# The empty C result for n < 2 samples.
cord_empty <- function (n, d) {
  list(
    scores = matrix(numeric(0), n, 0), eig = numeric(0), n_con = 0L,
    con_eig = numeric(0), total = 0, conditional = 0, constrained = 0,
    kept = rep(TRUE, ncol(d$design)), biplot = matrix(numeric(0), ncol(d$x), 0),
    u = matrix(numeric(0), n, 0), real = 0, imaginary = 0 )
}



# db-RDA (or, with `positive`, capscale) from a numeric design. `dv` is from
# dist_vector(); `x` (constraints) and `z` (conditions) are numeric matrices
# with one row per sample, in the order of dv$labels, or NULL. `mf` is a
# data.frame of the constraints' variables, for factor-level centroids, or
# NULL.
dbrda_core <- function (
    dv, x, z = NULL, mf = NULL, k = 3L, cpus = n_cpus(), positive = FALSE ) {

  n <- dv$n
  d <- cord_design(x, z, n)


  # vegan's dbrda() and capscale() divide distances by sqrt(n - 1) when any
  # is 4 or more ("mean squared distances"), which divides G by n - 1.
  adj <- 1
  if (n >= 2 && isTRUE(max(dv$dm) >= 4 + .Machine$double.eps))
    adj <- 1 / (n - 1)


  if (n < 2) {
    res <- cord_empty(n, d)
  }
  else {
    res <- .Call(C_dbrda, dv$dm, n, d$design, ncol(d$z), k, adj, positive, cpus)

    if (!res$converged) # nocov start
      warning(
        if (positive) 'capscale' else 'db-RDA',
        ' did not converge; coordinates may be inexact.', call. = FALSE )
    # nocov end
  }

  cord_result(
    res, d, mf, k, n, dv$labels,
    method   = if (positive) 'capscale' else 'dbrda',
    prefix   = if (positive) 'CAP' else 'dbRDA',
    prefix_u = 'MDS',
    positive = positive )
}



# The ecodive_ord result of a constrained ordination, from C_dbrda()'s,
# C_rda()'s or C_cca()'s list `res`, the design `d` (cord_design()), and
# `mf`. Columns are named `prefix` (constrained) and `prefix_u`
# (unconstrained); `...` adds attributes at the end. `scale` is vegan's
# scores() constant (NULL: rda()'s and dbrda()'s), `r2` replaces
# dbrda_r_squared()'s pair, and `weights` are cca()'s sqrt(r), for weighted
# centroids.
cord_result <- function (
    res, d, mf, k, n, labels, method, prefix, prefix_u, positive = FALSE,
    scale = NULL, r2 = NULL, weights = NULL, ... ) {

  x <- d$x
  z <- d$z
  design <- d$design

  if (any(aliased <- !res$kept))
    message(
      'Dropped model columns that are linear combinations of the ones ',
      'before them (or of the conditions): ',
      paste(sQuote(colnames(design)[aliased]), collapse = ', ') )


  total <- res$total     # vegan's tot.chi (includes imaginary inertia)
  real  <- res$real      # trace(G+) for capscale; total for dbrda
  kc    <- res$n_con
  m     <- ncol(res$scores)
  if (is.null(scale))
    scale <- sqrt(sqrt((n - 1) * total)) # vegan's scores() constant

  # Axes with eigenvalues that aren't positive have no coordinates.
  if (m < min(k, n - 1L))
    warning(sprintf('only %d of the first %d eigenvalues are > 0', m, min(k, n - 1L)), call. = FALSE)

  cols <- c(sprintf('%s%d', prefix, seq_len(kc)), sprintf('%s%d', prefix_u, seq_len(m - kc)))

  kept_x <- res$kept[ncol(z) + seq_len(ncol(x))]
  biplot <- matrix(0, nrow(res$biplot), m, dimnames = list(
    colnames(x)[kept_x],
    if (m) cols ))
  biplot[, seq_len(kc)] <- res$biplot

  # Rounding can leave a tiny negative remainder.
  unconstrained <- real - res$conditional - res$constrained
  if (abs(unconstrained) < sqrt(.Machine$double.eps) * abs(total)) unconstrained <- 0

  inertia <- c(
    conditional   = res$conditional,
    constrained   = res$constrained,
    unconstrained = unconstrained,
    imaginary     = if (positive) res$imaginary )

  centroids <- dbrda_centroids(res$u, mf, kc, m, scale, cols, weights)

  if (is.null(r2)) r2 <- dbrda_r_squared(
    con = res$constrained, cond = res$conditional, real = real, total = total,
    n = n, p = sum(kept_x), q = sum(res$kept[seq_len(ncol(z))]) )

  ord_result(
    res$scores * scale, method, labels, cols = cols,
    eig             = res$eig,
    prop            = res$eig / (real - res$conditional),
    trace           = total,
    inertia         = inertia,
    constrained_eig = res$con_eig,
    r_squared       = r2[[1]],
    adj_r_squared   = r2[[2]],
    biplot          = biplot,
    centroids       = centroids,
    ... )
}



# R-squared and adjusted R-squared, as vegan::RsquareAdj() gives them for
# dbrda() and capscale() results: the constrained inertia over the real
# total (for capscale, trace(G+); for dbrda, the total), adjusted for p
# constraint columns. With q > 0 condition columns, explaining cond / total,
# the adjusted R-squared of both together less that of the conditions
# alone. NA without residual degrees of freedom.
dbrda_r_squared <- function (con, cond, real, total, n, p, q) {

  if (!(real > 0) || !(total > 0)) return (c(NA_real_, NA_real_))

  r2  <- con / real
  adj <- function (r, m) 1 - (1 - r) * (n - 1) / (n - m - 1)

  if (n - p - q - 1 <= 0) return (c(r2, NA_real_))

  if (q == 0) return (c(r2, adj(r2, p)))

  rc <- cond / total
  c(r2, adj(r2 + rc, p + q) - adj(rc, q))
}



# vegan's factor-level centroids (centroids.cca() and getCentroids()): the
# mean unit constrained eigenvector of each level of each factor or character
# variable in the model frame, dropping near-zero rows, then scaled. With
# `weights` (cca()'s sqrt(r)), the weighted mean (weights r) of the lc
# scores u / sqrt(r).
dbrda_centroids <- function (u, mf, kc, m, scale, cols, weights = NULL) {

  if (is.null(mf) || ncol(u) == 0) return (NULL)

  facts <- vapply(mf, function (v) is.factor(v) || is.character(v), logical(1))
  if (!any(facts)) return (NULL)

  out <- NULL
  for (nm in names(mf)[facts]) {
    f   <- droplevels(as.factor(mf[[nm]]))
    if (is.null(weights)) {
      cnt <- tabulate(f, nlevels(f))
      cn  <- rowsum(u, f, reorder = TRUE) / cnt
    }
    else {
      cnt <- as.vector(rowsum(weights ^ 2, f, reorder = TRUE))
      cn  <- rowsum(u * weights, f, reorder = TRUE) / cnt
    }
    rownames(cn) <- paste0(nm, levels(f))
    out <- rbind(out, cn)
  }

  out <- unique(out)
  out <- out[rowSums(out ^ 2) > 1e-4, , drop = FALSE]
  if (nrow(out) == 0) return (NULL)

  cn <- matrix(0, nrow(out), m, dimnames = list(rownames(out), cols))
  cn[, seq_len(kc)] <- out[, seq_len(kc)] * scale
  cn
}
