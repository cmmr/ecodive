# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Constrained Analysis of Principal Coordinates (capscale)
#'
#' Constrained ordination of samples from their distances, like [dbrda()],
#' but of the principal coordinates with positive eigenvalues only. Gives
#' the same results as `vegan::capscale()` (and so as
#' `phyloseq::ordinate(method = "CAP")`). For Euclidean distances it gives
#' the same results as [dbrda()].
#'
#' @inheritParams dbrda
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   as described in "Ordination result" below. Constrained axes come first,
#'   named `"CAP1"`, `"CAP2"`, etc., then, if the constraints give fewer than
#'   `n_axes` of them, unconstrained axes named `"MDS1"`, `"MDS2"`, etc. The values
#'   are `vegan`'s site scores ("wa" scores, scaling 2), as for [dbrda()].
#'   Its attributes are:
#'
#'   * `method`: `"capscale"`.
#'   * `eig`: The eigenvalue of each axis.
#'   * `prop`: The proportion of the real inertia (see `inertia`) that each
#'     axis explains, not counting the conditional inertia:
#'     `eig / (trace - inertia[["imaginary"]] - inertia[["conditional"]])`.
#'     `vegan`'s "Proportion Explained".
#'   * `trace`: The total inertia, `vegan`'s: the sum of all eigenvalues of
#'     \eqn{G}, including negative ones, as [dbrda()]'s.
#'   * `inertia`: The total inertia split into four parts, which add up to
#'     `trace`: a named vector of
#'     - `conditional`: explained by the `Condition()` variables (`0` if
#'       there are none);
#'     - `constrained`: explained by the constraints, after the conditions;
#'     - `unconstrained`: the rest of the real inertia;
#'     - `imaginary`: the sum of the negative eigenvalues of \eqn{G}, which
#'       `capscale()` leaves out (`0` for Euclidean distances).
#'
#'     The first three add up to the real inertia, `trace - imaginary`, the
#'     sum of the positive eigenvalues of \eqn{G}.
#'   * `constrained_eig`: Every non-zero constrained eigenvalue.
#'   * `r_squared`: The proportion of the real inertia that the constraints
#'     explain, `inertia[["constrained"]] / (trace - inertia[["imaginary"]])`,
#'     as `vegan::RsquareAdj(x)$r.squared`.
#'   * `adj_r_squared`: `r_squared` adjusted for the number of constraint
#'     columns, as `vegan::RsquareAdj(x)$adj.r.squared` and as for
#'     [dbrda()]. With conditions, `vegan` takes their proportion of `trace`,
#'     not of the real inertia, and so does `capscale()`.
#'   * `biplot`, `centroids`: As for [dbrda()].
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{G = -\frac{1}{2}JAJ} as in [pcoa()] and [dbrda()], and its
#' eigen-decomposition \eqn{G = V\Lambda V^T}, the principal coordinates with
#' positive eigenvalues are \eqn{Y = V_+\Lambda_+^{1/2}}, and
#' \eqn{G_+ = YY^T} is the positive part of \eqn{G}. `capscale()` is the
#' redundancy analysis of \eqn{Y} on the constraints, which is [dbrda()] with
#' \eqn{G_+} in place of \eqn{G}. As in `vegan`, eigenvalues within
#' `max(e, e * lambda1)` of 0, for \eqn{e} = `sqrt(.Machine$double.eps)` and
#' \eqn{\lambda_1} the largest, count as 0.
#'
#' For Euclidean distances (such as Euclidean, Hellinger, or Aitchison
#' distances), \eqn{G} has no negative eigenvalues, so \eqn{G_+ = G}.
#' `capscale()` first looks for a negative eigenvalue, by block Lanczos
#' iteration on \eqn{-G} with products straight from `dm` (as [pcoa()]
#' finds the largest eigenvalues). It stops as soon as one is certain, and
#' otherwise runs until the smallest eigenvalue of \eqn{G} has converged.
#' Having found none, the rest is exactly [dbrda()]'s, in as little memory.
#' That is quick when \eqn{G} has many eigenvalues at 0, such as for
#' Euclidean distances between samples with fewer features than samples.
#' When the smallest eigenvalues are close together, convergence can take
#' many steps; after about \eqn{n/10} of them (each one pass over `dm`), when
#' going on would cost about as much as the full eigen-decomposition,
#' `capscale()` stops looking and decomposes \eqn{G}, as below, which gives
#' the same results either way.
#'
#' Otherwise, \eqn{G_+} needs every eigenvalue of \eqn{G}, as in
#' `vegan::capscale()`: `capscale()` forms \eqn{G} (\eqn{8n^2} bytes, such as
#' 200 MB for 5000 samples), reduces it to tridiagonal form, finds all of its
#' eigenvalues, and then the eigenvectors of the positive ones, or of the
#' others if they are fewer, by inverse iteration. That takes time
#' proportional to \eqn{n^3}, though several times less than a full
#' eigen-decomposition: for 2000 Bray-Curtis distances, about 2 seconds on
#' one core, against 12 for `vegan::capscale()`. The reduction and the
#' eigenvectors are computed in parallel on `cpus` threads.
#'
#'
#' **capscale() or dbrda()?**
#'
#' See [dbrda()]. In short: `dbrda()` analyzes the distances as they are,
#' and quickly for any number of samples; `capscale()` analyzes the closest
#' Euclidean configuration, so its eigenvalues are never negative, but it
#' ignores the rest and takes much longer for large non-Euclidean data.
#' Another way to avoid negative eigenvalues is to make the distances
#' Euclidean first, such as by the square roots of Bray-Curtis distances
#' (Legendre & Legendre 2012, Table 7.2), or with [pcoa()]'s `correction`,
#' and then to use `dbrda()`.
#'
#'
#' **Matching vegan**
#'
#' `capscale()` gives the same scores, eigenvalues, inertia (including the
#' imaginary inertia), \eqn{R^2}, biplot scores and centroids as
#' `vegan::capscale()` given distances, up to the signs of the axes. As
#' there, distances are divided by \eqn{\sqrt{n-1}} when any is 4 or more
#' (see "Scale" in [dbrda()]). `vegan::capscale()` can also take a community
#' table, add a constant to the distances (`add`), or take their square roots
#' (`sqrt.dist`), and gives species scores; `capscale()` does not.
#'
#' `capscale()` is named after `vegan::capscale()`. Anderson & Willis's (2003)
#' "canonical analysis of principal coordinates" (CAP) is a related but
#' different method.
#'
#'
#' **Axis signs and reproducibility**
#'
#' As for [dbrda()]: each axis's largest score in absolute value is
#' positive, and results are identical from run to run and for any number
#' of `cpus`.
#'
#' @references
#' Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy
#' analysis: testing multispecies responses in multifactorial ecological
#' experiments. *Ecological Monographs*, 69(1), 1-24.
#' \doi{10.1890/0012-9615(1999)069[0001:DBRATM]2.0.CO;2}
#'
#' Anderson, M. J., & Willis, T. J. (2003). Canonical analysis of principal
#' coordinates: a useful method of constrained ordination for ecology.
#' *Ecology*, 84(2), 511-525.
#' \doi{10.1890/0012-9658(2003)084[0511:CAOPCA]2.0.CO;2}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd English
#' ed.). Elsevier.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [dbrda()], [pcoa()], [permanova()], [dist_ord()], [beta_div()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- capscale(dm, 'Species', iris)
#'     plot(ord, col = iris$Species)
#'     head(ord)
#'     str(ord)
#'
capscale <- function (
    dm, groups, data = NULL, formula = NULL, n_axes = 3L, cpus = n_cpus() ) {

  validate_cpus()
  k <- ord_k(n_axes)

  no_groups <- missing(groups)
  constrained_ord(
    'capscale', dm, if (no_groups) NULL else groups, no_groups,
    data, formula, k, cpus, env = parent.frame() )
}
