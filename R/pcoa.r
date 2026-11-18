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
#' @param correction   How to deal with negative eigenvalues: `'none'` (the
#'        default), `'lingoes'`, or `'cailliez'`. See "Negative eigenvalue
#'        corrections" below.
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
#'   * `correction`: `"none"`, `"lingoes"`, or `"cailliez"`.
#'   * `constant`: The constant that the correction added, \eqn{c_1} or
#'     \eqn{c_2} below; `0` with `correction = "none"`, or when the distances
#'     needed no correction.
#'
#'   With a correction, `eig`, `prop`, and `trace` are those of the corrected
#'   distances.
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
#' **Negative eigenvalue corrections**
#'
#' A correction adds a constant to the distances between different samples,
#' the smallest that makes them Euclidean, so that no eigenvalue is negative
#' (Gower & Legendre 1986; Legendre & Legendre 2012, sec. 9.3):
#'
#' * `'lingoes'` (Lingoes 1971) adds \eqn{2c_1} to the squared distances:
#'   \eqn{d'_{ij} = \sqrt{d_{ij}^2 + 2c_1}}, where \eqn{c_1} is the magnitude
#'   of \eqn{B}'s most negative eigenvalue. That adds \eqn{c_1} to every
#'   eigenvalue but the trivial 0 and leaves the eigenvectors unchanged, so
#'   the axes are those of `correction = 'none'`, with coordinates
#'   \eqn{u_i\sqrt{\lambda_i + c_1}}, and `trace` grows by \eqn{(n-1)c_1}.
#'   The axis whose eigenvalue was most negative gets eigenvalue 0.
#' * `'cailliez'` (Cailliez 1983) adds \eqn{c_2} to the distances:
#'   \eqn{d'_{ij} = d_{ij} + c_2}, where \eqn{c_2} is the largest real
#'   eigenvalue of the \eqn{2n \times 2n} matrix
#'   \eqn{\left[\begin{smallmatrix}0 & 2B \\ -I & -4B_2\end{smallmatrix}\right]}{[[0, 2B], [-I, -4 B2]]},
#'   with \eqn{B_2 = -\frac{1}{2}JDJ} the centered matrix of the unsquared
#'   distances. The axes change, and are computed from \eqn{d'}.
#'
#' `eig`, `prop`, and `trace` are then those of the corrected distances, and
#' `prop` sums to at most 1. They equal `ape::pcoa()`'s `Corr_eig` and
#' `Rel_corr_eig` for the same correction, and `vegan::wcmdscale(add =
#' 'lingoes')` or `add = 'cailliez'`; \eqn{c_2} is `cmdscale(add = TRUE)`'s
#' `ac`.
#'
#' Negative eigenvalues above \eqn{-\sqrt{\epsilon}\, trace(B)} (with
#' \eqn{\epsilon} = `.Machine$double.eps`, so about \eqn{-1.5 \times 10^{-8}}
#' of the total) count as rounding error: such distances, Euclidean ones for
#' instance, get no correction. `constant` is then 0, and the coordinates,
#' `eig`, `prop`, and `trace` are identical to those from
#' `correction = 'none'`. (`ape` and `vegan` correct those tiny negative
#' eigenvalues too, so they can add a tiny constant where `pcoa()` adds 0.)
#'
#' With `k` = \eqn{n - 1}, the axis whose eigenvalue the correction makes 0
#' is dropped without a warning.
#'
#'
#' **When to use a correction**
#'
#' Non-Euclidean distances, such as Bray-Curtis or Jaccard, have negative
#' eigenvalues, so the positive axes' `prop` can sum to more than 1. Most
#' plots simply show the first few positive axes, which is fine when the
#' negative eigenvalues are small next to the axes plotted. A correction
#' makes every `prop` a share of a total of 1, at a price: both corrections
#' inflate `trace`, and so shrink `prop`, more so for more samples. Which
#' inflates it more depends on the data. `'lingoes'` keeps the uncorrected
#' axes and only changes their eigenvalues; `'cailliez'` changes the axes
#' too. Another remedy is to transform the distances: the square root of
#' Bray-Curtis dissimilarities, for one, is Euclidean (Legendre & Legendre
#' 2012).
#'
#'
#' **How the corrections are computed**
#'
#' Neither the \eqn{n \times n} nor the \eqn{2n \times 2n} matrix is formed.
#' \eqn{c_1} is the largest eigenvalue of \eqn{-B}, from the same Lanczos
#' solver. \eqn{c_2} is the largest \eqn{c} for which
#' \eqn{B + 2cB_2 + \frac{c^2}{2}I} is singular (on vectors orthogonal to
#' \eqn{\mathbf{1}}), found by a Rayleigh-Ritz iteration that needs only
#' products of \eqn{B} and \eqn{B_2} with vectors, computed together from
#' `dm`. That iteration can settle on a root that isn't the largest (as with
#' some exactly symmetric designs), so \eqn{c_2} is checked: the corrected
#' matrix's smallest eigenvalue, computed from a fresh start, must not be
#' below \eqn{-\sqrt{\epsilon}\, trace(B')}. If it is, its eigenvector
#' restarts the iteration, for up to 10 rounds, after which there is a
#' warning.
#'
#' The smallest eigenvalues usually form a dense cluster, which Lanczos
#' resolves slowly, so corrections take longer than `correction = 'none'`.
#' For Bray-Curtis distances between 2,000 sparse samples, on one thread,
#' `'none'` took 0.4 seconds, `'lingoes'` 5 seconds, and `'cailliez'` 20
#' seconds; `cmdscale(add = TRUE)`, which builds the \eqn{2n \times 2n}
#' matrix, took 2 minutes. How much longer depends on how tightly the
#' smallest eigenvalues cluster. Each pass over `dm` takes time proportional
#' to \eqn{n^2}, while `cmdscale(add = TRUE)` and `ape::pcoa()` take time
#' proportional to \eqn{n^3}.
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
#' Cailliez, F. (1983). The analytical solution of the additive constant
#' problem. *Psychometrika*, 48(2), 305-308. \doi{10.1007/BF02294026}
#'
#' Gower, J. C. (1966). Some distance properties of latent root and vector
#' methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
#' \doi{10.1093/biomet/53.3-4.325}
#'
#' Gower, J. C., & Legendre, P. (1986). Metric and Euclidean properties of
#' dissimilarity coefficients. *Journal of Classification*, 3(1), 5-48.
#' \doi{10.1007/BF01896809}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English ed.,
#' sec. 9.3). Elsevier.
#'
#' Lingoes, J. C. (1971). Some boundary conditions for a monotone analysis of
#' symmetric matrices. *Psychometrika*, 36(2), 195-203.
#' \doi{10.1007/BF02291398}
#'
#' Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
#' symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
#' Applications*, 22(2), 602-616. \doi{10.1137/S0895479898334605}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa_biplot()], [pca()], [dbrda()], [nmds()], [dist_ord()],
#'   [beta_div()], `stats::cmdscale()`
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- pcoa(dm)
#'     plot(ord, col = iris$Species)
#'     head(ord)
#'     str(ord)
#'
pcoa <- function (dm, k = 3L, correction = 'none', cpus = n_cpus()) {

  validate_cpus()
  k <- ord_k(k)

  if (!is.character(correction) || length(correction) != 1 || is.na(correction) ||
      is.na(correction <- pmatch(correction, c('none', 'lingoes', 'cailliez'))))
    stop('`correction` must be "none", "lingoes", or "cailliez".', call. = FALSE)
  id         <- correction - 1L
  correction <- c('none', 'lingoes', 'cailliez')[[correction]]

  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  k <- min(k, n - 1L)

  constant <- 0

  if (k < 1) {

    # One sample or none: no axes.
    points <- matrix(numeric(0), nrow = n, ncol = 0)
    eig    <- numeric(0)
    trace  <- 0

  } else {

    res <- if (id == 0L) .Call(C_pcoa, dm, n, k, cpus)
           else          .Call(C_pcoa_correct, dm, n, k, id, cpus)

    if (!res$converged) # nocov start
      warning('PCoA did not converge; coordinates may be inexact.', call. = FALSE)
    if (isFALSE(res$checked))
      warning(
        'The Cailliez constant could not be confirmed as the largest root; ',
        'some corrected eigenvalues may be negative.', call. = FALSE )
    # nocov end

    if (id > 0L) constant <- res$constant

    # A correction makes one eigenvalue 0 by design, the last of all n - 1:
    # with k = n - 1, it is dropped without a warning.
    keep <- if (constant > 0 && k == n - 1L) {
      c(ord_keep(res$eig[-k], k - 1L), FALSE)
    } else {
      ord_keep(res$eig, k)
    }

    points <- res$points[, keep, drop = FALSE]
    eig    <- res$eig[keep]
    trace  <- res$trace
  }

  ord_result(
    points, 'pcoa', labels, prefix = 'PCoA',
    eig        = eig,
    prop       = if (trace > 0) eig / trace else eig,
    trace      = trace,
    correction = correction,
    constant   = constant )
}
