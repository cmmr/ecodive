# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Principal Coordinates Analysis (PCoA)
#'
#' Ordinates samples from their distances, for plotting beta diversity. Also
#' known as classical or metric multidimensional scaling. Gives the same
#' coordinates and eigenvalues as `stats::cmdscale()`, `ape::pcoa()`, and
#' `vegan::wcmdscale()`, but computes only the first few axes, which takes far
#' less time and memory for large datasets.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values.
#'
#' @param k   The number of axes to compute. Default: `3`. Reduced, with a
#'        warning, to the number of axes with positive eigenvalues, which is
#'        at most \eqn{n - 1} for \eqn{n} samples.
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   named `"PCoA1"`, `"PCoA2"`, etc., as described in "Ordination result"
#'   below. Its attributes are:
#'
#'   * `method`: `"pcoa"`.
#'   * `eig`: The eigenvalue of each axis: of Gower's centered matrix, so on
#'     the scale of sums of squares. [pca()]'s `eig` are variances, these
#'     divided by \eqn{n - 1}; [dbrda()]'s are too when any distance is 4 or
#'     more.
#'   * `prop`: The proportion of the total variance that each axis explains,
#'     `eig / trace`.
#'   * `trace`: The total variance: the sum of all eigenvalues, including
#'     negative ones.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{A} the matrix of squared distances and \eqn{J = I - \frac{1}{n}\mathbf{1}\mathbf{1}^T}
#' the centering matrix, PCoA finds the eigenvalues \eqn{\lambda_i} and unit
#' eigenvectors \eqn{u_i} of \eqn{B = -\frac{1}{2}JAJ}. Axis \eqn{i}'s
#' coordinates are \eqn{u_i\sqrt{\lambda_i}}.
#'
#' `cmdscale()`, `ape::pcoa()`, and `vegan::wcmdscale()` build \eqn{B} in full
#' (\eqn{n \times n}) and compute all \eqn{n} of its eigenvalues, which takes
#' time proportional to \eqn{n^3}. `pcoa()` finds only the largest `k`
#' eigenvalues, by block Lanczos iteration with thick restarts. That only needs
#' products of \eqn{B} with a few vectors at a time, which are computed in
#' parallel straight from `dm`, without ever forming \eqn{B}.
#'
#'
#' **Proportion of variance**
#'
#' The total variance is \eqn{trace(B) = \frac{1}{n}\sum_{i<j}d_{ij}^2}, the
#' sum of all \eqn{n} eigenvalues. `prop` is each axis's eigenvalue divided by
#' it: the same as `ape::pcoa()`'s `Relative_eig` and `vegan`'s "Proportion
#' Explained", the values that `phyloseq::plot_ordination()` shows on its axes.
#'
#' Non-Euclidean distances, such as Bray-Curtis, can give negative eigenvalues.
#' They still count towards the total, so `prop` can sum to more than 1 over
#' the positive axes. Axes with eigenvalues that aren't positive have no real
#' coordinates and are dropped with a warning, as `cmdscale()` does. So are
#' axes with eigenvalues at rounding-error level, below
#' `sqrt(.Machine$double.eps)` times the first one, as in `vegan::wcmdscale()`.
#'
#'
#' **Axis signs and reproducibility**
#'
#' An eigenvector's sign is arbitrary, so different software can return
#' mirror images of the same axis. `pcoa()` flips each axis so that its
#' largest coordinate in absolute value is positive. Results are identical
#' from run to run and for any number of `cpus`.
#'
#' `k` larger than \eqn{n - 1} (the most axes possible) is reduced to
#' \eqn{n - 1}.
#'
#' @references
#' Gower, J. C. (1966). Some distance properties of latent root and vector
#' methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
#' \doi{10.1093/biomet/53.3-4.325}
#'
#' Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
#' symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
#' Applications*, 22(2), 602-616. \doi{10.1137/S0895479898334605}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pca()], [dbrda()], [nmds()], [ordination()], [beta_div()],
#'   `stats::cmdscale()`
#'
#' @export
#' @examples
#'     ord <- pcoa(bray(ex_counts))
#'     ord
#'
#'     # Percent of variance explained by each axis.
#'     round(100 * attr(ord, 'prop'), 1)
#'
#'     # Same coordinates as cmdscale(), up to the sign of each axis.
#'     cmdscale(bray(ex_counts), k = 3)
#'
#'     plot(ord[, 1:2], asp = 1)
#'
pcoa <- function (dm, k = 3L, cpus = n_cpus()) {

  validate_cpus()
  k <- ord_k(k)

  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  k <- min(k, n - 1L)

  if (k < 1) {

    # One sample or none: no axes.
    points <- matrix(numeric(0), nrow = n, ncol = 0)
    eig    <- numeric(0)
    trace  <- 0

  } else {

    res <- .Call(C_pcoa, dm, n, k, cpus)

    if (!res$converged) # nocov start
      warning('PCoA did not converge; coordinates may be inexact.', call. = FALSE)
    # nocov end

    keep   <- ord_keep(res$eig, k)
    points <- res$points[, keep, drop = FALSE]
    eig    <- res$eig[keep]
    trace  <- res$trace
  }

  ord_result(
    points, 'pcoa', labels, prefix = 'PCoA',
    eig   = eig,
    prop  = if (trace > 0) eig / trace else eig,
    trace = trace )
}
