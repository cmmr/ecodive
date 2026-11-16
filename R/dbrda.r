# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Distance-based Redundancy Analysis (db-RDA)
#'
#' Constrained ordination of samples from their distances: the axes of beta
#' diversity that sample metadata explain. Gives the same results as
#' `vegan::dbrda()` (not `vegan::capscale()`; see "Matching vegan" below),
#' without ever forming an \eqn{n \times n}{n x n} matrix, so it takes far less
#' time and memory for large datasets.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param data  The sample metadata: a data frame with one row per sample.
#'        Rows are matched to the samples in `dm` by row name; without row
#'        names, they must be in the same order as `dm`. A single vector or
#'        factor is also accepted, as one variable named `group`.
#'
#' @param formula  A one-sided formula naming the variables in `data` that
#'        constrain the ordination, such as `~ treatment + ph`. Factors,
#'        numeric covariates, and interactions are all allowed. Wrap
#'        variables in `Condition()` to remove their effect first (a partial
#'        db-RDA), as in `vegan::dbrda()`: `~ treatment + Condition(site)`.
#'        A left-hand side, if any, is ignored. Default: `NULL`, for every
#'        column of `data`.
#'
#' @param k   The number of axes to return. Default: `3`
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   as described in "Ordination result" below. Constrained axes come first,
#'   named `"dbRDA1"`, `"dbRDA2"`, etc. If the constraints give fewer than `k`
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
#' **Matching vegan**
#'
#' `dbrda()` matches `vegan::dbrda()`. It does not match `vegan::capscale()`
#' or `phyloseq::ordinate(method = "CAP")` (which calls `capscale()`):
#' `capscale()` first removes the negative eigenvalues of non-Euclidean
#' distances, such as Bray-Curtis, so its eigenvalues and scores differ from
#' `dbrda()`'s. For Euclidean distances, the two agree.
#'
#' Columns of the model matrix that are linear combinations of earlier ones
#' (aliased) are dropped with a message, as in `vegan`.
#'
#'
#' **Scale**
#'
#' Like `vegan::dbrda()`, `dbrda()` divides the distances by
#' \eqn{\sqrt{n-1}} when any of them is 4 or more ("mean squared
#' distances"). That divides `trace`, `eig`, `inertia` and `constrained_eig`
#' by \eqn{n-1}, but leaves `prop` unchanged. So for the same `dm`,
#' `attr(dbrda(dm, ...), 'trace')` equals `attr(pcoa(dm), 'trace')` when every
#' distance is below 4 (e.g. Bray-Curtis or Jaccard), and equals it divided
#' by \eqn{n-1} otherwise (e.g. Euclidean distances between counts).
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
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [permanova()], [beta_div()]
#'
#' @export
#' @examples
#'     d    <- bray(ex_counts)
#'     meta <- data.frame(
#'       row.names = rownames(ex_counts),
#'       treated   = c(TRUE, TRUE, FALSE, FALSE) )
#'
#'     ord <- dbrda(d, meta, ~ treated)
#'     ord
#'
#'     # Percent of the inertia explained by each axis.
#'     round(100 * attr(ord, 'prop'), 1)
#'
#'     plot(ord[, 1:2])
#'
dbrda <- function (dm, data, formula = NULL, k = 3L, cpus = n_cpus()) {

  validate_cpus()

  tryCatch(
    expr = {
      stopifnot(is.numeric(k), length(k) == 1, !is.na(k), k >= 1, k %% 1 == 0)
      k <- as.integer(min(k, .Machine$integer.max))
    },
    error = function (e)
      stop(e$message, '\n`k` must be a whole number of at least 1.', call. = FALSE)
  )

  dv     <- dist_vector(dm)
  design <- model_design(data, formula, labels = dv$labels, n = dv$n)

  # The constraints' variables, for factor-level centroids (vegan leaves out
  # the conditions' variables).
  mf <- stats::model.frame(design$terms, design$data, na.action = stats::na.fail)
  attr(mf, 'terms') <- NULL

  dbrda_core(dv, design$x, design$z, mf, k, cpus)
}



# db-RDA from a numeric design. `dv` is from dist_vector(); `x` (constraints)
# and `z` (conditions) are numeric matrices with one row per sample, in the
# order of dv$labels, or NULL. `mf` is a data.frame of the constraints'
# variables, for factor-level centroids, or NULL.
dbrda_core <- function (dv, x, z = NULL, mf = NULL, k = 3L, cpus = n_cpus()) {

  n      <- dv$n
  labels <- dv$labels

  if (is.null(x)) x <- matrix(0, n, 0)
  if (is.null(z)) z <- matrix(0, n, 0)
  x <- as.matrix(x)
  z <- as.matrix(z)
  if (is.null(colnames(x)) && ncol(x)) colnames(x) <- paste0('X', seq_len(ncol(x)))

  design <- cbind(z, x)
  storage.mode(design) <- 'double'


  # vegan's dbrda() divides distances by sqrt(n - 1) when any is 4 or more
  # ("mean squared distances"), which divides G by n - 1.
  adj <- 1
  if (n >= 2 && isTRUE(max(dv$dm) >= 4 + .Machine$double.eps))
    adj <- 1 / (n - 1)


  if (n < 2) {
    res <- list(
      scores = matrix(numeric(0), n, 0), eig = numeric(0), n_con = 0L,
      con_eig = numeric(0), total = 0, conditional = 0, constrained = 0,
      kept = rep(TRUE, ncol(design)), biplot = matrix(numeric(0), ncol(x), 0),
      u = matrix(numeric(0), n, 0) )
  }
  else {
    res <- .Call(C_dbrda, dv$dm, n, design, ncol(z), k, adj, cpus)

    if (!res$converged) # nocov start
      warning('db-RDA did not converge; coordinates may be inexact.', call. = FALSE)
    # nocov end
  }

  if (any(aliased <- !res$kept))
    message(
      'Some constraints or conditions were aliased because they were ',
      'redundant: ', paste(sQuote(colnames(design)[aliased]), collapse = ', ') )


  total <- res$total
  kc    <- res$n_con
  m     <- ncol(res$scores)
  scale <- sqrt(sqrt((n - 1) * total)) # vegan's scores() constant

  points <- res$scores * scale
  dimnames(points) <- list(labels, c(
    sprintf('dbRDA%d', seq_len(kc)),
    sprintf('MDS%d',   seq_len(m - kc)) ))

  biplot <- matrix(0, nrow(res$biplot), m, dimnames = list(
    colnames(x)[res$kept[ncol(z) + seq_len(ncol(x))]],
    colnames(points) ))
  biplot[, seq_len(kc)] <- res$biplot

  inertia <- c(
    conditional   = res$conditional,
    constrained   = res$constrained,
    unconstrained = total - res$conditional - res$constrained )

  structure(
    .Data           = points,
    method          = 'dbrda',
    eig             = res$eig,
    prop            = res$eig / (total - res$conditional),
    trace           = total,
    inertia         = inertia,
    constrained_eig = res$con_eig,
    biplot          = biplot,
    centroids       = dbrda_centroids(res$u, mf, kc, m, scale, colnames(points)) )
}



# vegan's factor-level centroids (centroids.cca() and getCentroids()): the
# mean unit constrained eigenvector of each level of each factor or character
# variable in the model frame, dropping near-zero rows, then scaled.
dbrda_centroids <- function (u, mf, kc, m, scale, cols) {

  if (is.null(mf) || ncol(u) == 0) return (NULL)

  facts <- vapply(mf, function (v) is.factor(v) || is.character(v), logical(1))
  if (!any(facts)) return (NULL)

  out <- NULL
  for (nm in names(mf)[facts]) {
    f   <- droplevels(as.factor(mf[[nm]]))
    cnt <- tabulate(f, nlevels(f))
    cn  <- rowsum(u, f, reorder = TRUE) / cnt
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
