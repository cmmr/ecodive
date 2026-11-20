# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The algorithm of mds() follows scikit-learn's metric SMACOF (BSD-3-Clause,
# Copyright (c) 2007-2026 The scikit-learn developers); its C port is in
# src/nmds.c, which carries scikit-learn's license in full.



#' Metric Multidimensional Scaling (MDS)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing them so that their distances on the plot match the distances in
#' `dm` as closely as possible, in least squares. Uses the SMACOF algorithm
#' of scikit-learn's metric `MDS()`, ported to parallel C that works from `dm`
#' directly, with far less time and memory.
#'
#' **Which ordination?** [pcoa()], `mds()`, and [nmds()] all fit a
#' configuration to `dm` in `n_axes` dimensions. [pcoa()] (classical scaling) is exact and fast, has
#' eigenvalues and "variance explained", but fits the squared distances after
#' double centering, so large distances dominate, and non-Euclidean
#' distances such as Bray-Curtis lose part of their structure to negative
#' eigenvalues. `mds()` fits the distances themselves, in their own units,
#' so the plot's distances can be read as the data's, with `stress` saying
#' how well. [nmds()] fits only their rank order, which suits distances
#' whose values matter less than their order, such as the saturated
#' Bray-Curtis values of very different samples; it is the usual choice for
#' community data. `mds()` and `nmds()` iterate from several starts and take
#' longer than `pcoa()`.
#'
#' @inheritParams nmds
#'
#' @param n_starts The number of random starts to run after the start from the
#'        PCoA solution, from 0 to 1000. Default: `4`
#'
#' @param tol   Each start stops once an iteration lowers the raw stress by
#'        no more than `tol` times the sum of squared distances on the plot.
#'        Default: `1e-7`
#'
#' @return A numeric matrix with one row per sample and one column per
#'   dimension, named `"MDS1"`, `"MDS2"`, etc., as described in "Ordination
#'   result" below. `method` is `"mds"`; `eig` and `prop` are `NA`, as metric
#'   MDS has no eigenvalues. Then come:
#'
#'   * `stress`: The stress of the solution, from 0 (a perfect fit) to 1 (no
#'     better than putting every point in one place). See "Stress" below.
#'   * `converged`: `TRUE` if the best start stopped because a plain step met
#'     `tol`, `FALSE` if it stopped at `max_iters`.
#'   * `n_iters`: The number of SMACOF iterations the best start ran.
#'   * `n_starts`: The number of random starts run.
#'   * `best_start`: Which start gave the solution: `0` for the PCoA start,
#'     else the random start's number.
#'   * `n_best`: How many starts (of `n_starts + 1`) ended within 0.1% of the
#'     lowest stress.
#'
#'   With fewer than two samples there is nothing to fit: the result has no
#'   columns, `stress` is 0, and no starts are run.
#'
#' @details
#'
#' **Method**
#'
#' Metric MDS minimizes the raw stress
#' \eqn{\sigma = \sum_{i<j}(d_{ij} - \delta_{ij})^2}, where \eqn{d_{ij}} are
#' the distances between points on the plot and \eqn{\delta_{ij}} the
#' distances in `dm`. It does so by SMACOF (de Leeuw, 1977; Borg & Groenen,
#' 2005), as scikit-learn's metric `smacof()` does: each iteration moves the
#' points by the Guttman transform, \eqn{X \leftarrow \Gamma(X) =
#' \frac{1}{n} B(X) X}, which never increases \eqn{\sigma}. A start stops
#' when an iteration lowers \eqn{\sigma} by no more than `tol` times
#' \eqn{\sum d_{ij}^2}, or after `max_iters` iterations.
#'
#' As in [nmds()] (see there), iterations take relaxed steps,
#' \eqn{X \leftarrow 2\Gamma(X) - X} (de Leeuw & Heiser, 1980), when `k > 1`,
#' until one makes no more than `tol` progress, and then plain steps until a
#' plain step does; `n_axes = 1` uses plain steps only, and so does the last
#' iteration allowed by `max_iters`. So the returned configuration is a fixed
#' point of \eqn{\Gamma} to within `tol`, at the right scale
#' (\eqn{\sum d_{ij}\delta_{ij} = \sum d_{ij}^2}), typically in about half the
#' iterations that plain steps take. With plain steps only (an internal
#' setting used by the tests), a start follows scikit-learn's iterations from
#' the same starting configuration to within rounding error (about
#' \eqn{10^{-13}}, relatively, in every iteration's raw stress).
#'
#' SMACOF only finds a local minimum of stress, so, as in [nmds()]:
#'
#' 1. The first start is the PCoA solution, from [pcoa()]. For distances
#'    that are Euclidean in `n_axes` dimensions, PCoA recovers the configuration,
#'    and SMACOF stops at once. (An axis without a positive eigenvalue gets
#'    small random coordinates instead.)
#' 2. Then come `n_starts` random starts, with coordinates uniform on \[0, 1).
#'    On beta diversity distances, a random start sometimes ends with a
#'    stress up to a few percent lower than the PCoA start's.
#' 3. The start with the lowest stress is kept; if several tie exactly, the
#'    earliest. Its configuration is centered, rotated to its principal
#'    axes, and each axis is flipped so that its coordinate with the largest
#'    absolute value is positive (the first such sample, if several are
#'    within a relative \eqn{10^{-8}} of it). It is not rescaled: distances
#'    on the plot are in the units of `dm`.
#'
#'
#' **Stress**
#'
#' The `stress` reported is
#' \eqn{\sqrt{\sigma / \sum \delta_{ij}^2}}: the raw stress, normalized by
#' the distances in `dm`. It only depends on \eqn{\sigma}, so the lowest
#' stress is the best fit, and it doesn't depend on the units of `dm`. It is
#' 1 when all points coincide, and at most 1 after any SMACOF iteration
#' (`max_iters = 0` returns the starting configuration, whose stress can be
#' higher). Pairs with a distance of 0 in `dm`, such as duplicate samples,
#' are pulled towards the same point.
#'
#' scikit-learn reports \eqn{\sigma} itself by default for metric MDS, or,
#' with `normalized_stress = True`, \eqn{\sqrt{\sigma / \sum d_{ij}^2}}
#' (normalized by the distances on the plot instead). For a result `x`, that
#' is `sqrt(sum((dist(x) - dm)^2) / sum(dist(x)^2))`.
#'
#'
#' **Differences from scikit-learn**
#'
#' * scikit-learn takes plain Guttman steps; `mds()` takes relaxed ones
#'   (see above).
#' * scikit-learn runs `n_init` random starts (`init = "random"`, its
#'   current default) or one start from classical MDS (`init =
#'   "classical_mds"`); `mds()` runs both kinds, and ranks them by the stress
#'   described above.
#' * scikit-learn's defaults are `max_iter = 300` and `eps = 1e-6`; `mds()`'s
#'   are `max_iters = 1000` and `tol = 1e-7`.
#' * scikit-learn returns the configuration as SMACOF leaves it; `mds()`
#'   centers, rotates, and flips it, which doesn't change the stress.
#' * scikit-learn stops when stress falls by *less* than `eps`; `mds()` stops
#'   when it falls by no more than `tol`, so that `tol = 0` stops once stress
#'   no longer falls.
#' * scikit-learn computes distances on the plot with the expansion
#'   \eqn{|x|^2 - 2 x \cdot y + |y|^2}, and replaces distances of 0 by
#'   \eqn{10^{-5}} in the Guttman transform; `mds()` computes them directly,
#'   and leaves out pairs at distance 0, whose contribution
#'   \eqn{r_{ij}(x_i - x_j)} is 0 either way.
#' * The random starts come from `seed`, not from R's or NumPy's random
#'   number generators, so `set.seed()` has no effect, and R's random number
#'   stream is left untouched.
#'
#'
#' **Speed and memory**
#'
#' `mds()` shares its C engine with [nmds()], without the sorting and the
#' isotonic regression: `dm` is read in its own order, once converted to a
#' plain vector of doubles. Each start running at the same time needs
#' \eqn{256 n k} bytes; nothing of size \eqn{n \times n} is formed, and
#' there is no limit on the number of samples. Each SMACOF iteration makes one pass over the pairs, which takes
#' time proportional to \eqn{n^2 k}. Several starts run at once, one per
#' thread, while one iteration takes under about \eqn{2^{26}} pairs of work
#' (\eqn{n} up to about 11,500); otherwise starts run one at a time, with
#' each iteration split across `cpus` threads (for at least \eqn{2^{18}}
#' pairs), in parts between checks for user interrupts.
#'
#'
#' **Reproducibility**
#'
#' The result, with all its attributes, is bit-identical from run to run, in
#' new R sessions, and for any value of `cpus`, given the same `dm`, `n_axes`,
#' `n_starts`, `max_iters`, `tol`, and `seed`, for the reasons given in [nmds()].
#' `dm * 2^p` gives exactly the result for `dm`, times \eqn{2^p}. (As for
#' every ordination, distances above \eqn{10^{150}} are rejected.) Results can
#' still differ in the last bits between computers or compilers, as for any
#' floating point code.
#'
#' `n_axes` larger than \eqn{n - 1} is reduced to \eqn{n - 1}.
#'
#' @references
#' de Leeuw, J. (1977). Applications of convex analysis to multidimensional
#' scaling. In J. R. Barra, F. Brodeau, G. Romier, & B. van Cutsem (Eds.),
#' *Recent Developments in Statistics* (pp. 133-145). North Holland.
#'
#' de Leeuw, J., & Heiser, W. J. (1980). Multidimensional scaling with
#' restrictions on the configurations. In P. R. Krishnaiah (Ed.),
#' *Multivariate Analysis V* (pp. 501-522). North Holland.
#'
#' Borg, I., & Groenen, P. J. F. (2005). *Modern Multidimensional Scaling:
#' Theory and Applications* (2nd ed.). Springer.
#' \doi{10.1007/0-387-28981-X}
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830. The C code of
#' `mds()` is ported from scikit-learn 1.9.1's `sklearn.manifold.smacof()`,
#' under scikit-learn's BSD-3-Clause license (Copyright (c) 2007-2026 The
#' scikit-learn developers).
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [nmds()], [beta_div()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- mds(dm)
#'     plot(ord, 'Species', iris)
#'     ord
#'     str(ord)
#'
mds <- function (
    dm, n_axes = 2L, n_starts = 4L, max_iters = 1000L, tol = 1e-7,
    seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  smacof_fit(dm, n_axes, n_starts, max_iters, tol, seed, cpus, 'mds')
}
