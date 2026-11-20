# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing samples that are near neighbours in `dm` near each other on the
#' plot. Follows the same procedure as `uwot::umap(batch = TRUE)` given a
#' `dist` object, and gives the same layout from the same settings, starting
#' coordinates and random seeds, but runs in parallel C, without ever copying
#' `dm` into an \eqn{n \times n} matrix.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param n_axes   The number of dimensions. Default: `2`
#'
#' @param n_neighbors   How many nearest neighbours of each sample (counting
#'        the sample itself) to use. Larger values preserve more of the
#'        global structure; smaller ones more of the local detail. Values
#'        larger than the number of samples are reduced to it. Default: `15`
#'
#' @param min_dist   How close samples can be packed together on the plot,
#'        relative to `spread`. Default: `0.1`, as in umap-learn (uwot's is
#'        `0.01`).
#'
#' @param spread   The scale of the plot's distances, together with
#'        `min_dist`. Default: `1`
#'
#' @param n_epochs   The number of optimization epochs. `NULL`, the
#'        default, uses 500 for up to 10,000 samples and 200 for more, as
#'        `uwot::umap()` does.
#'
#' @param learning_rate   The initial learning rate (step size) of the
#'        optimization. It decreases linearly to 0 over the epochs.
#'        Default: `1`
#'
#' @param negative_sample_rate   How many non-neighbours (negative samples)
#'        to push away for each neighbour pulled closer. Default: `5`
#'
#' @param repulsion_strength   The weight of the repulsion from negative
#'        samples, relative to the attraction of neighbours. Default: `1`
#'
#' @param set_op_mix_ratio   How the neighbour graphs of the samples are
#'        combined, from `0` (fuzzy intersection) to `1` (fuzzy union).
#'        Default: `1`
#'
#' @param local_connectivity   How many nearest neighbours of each sample
#'        are assumed to be fully connected to it. Default: `1`
#'
#' @param bandwidth   Scales the effective number of neighbours, as in
#'        `uwot::umap()`. Default: `1`
#'
#' @param init   The starting coordinates: `"spectral"` (the default),
#'        `"pcoa"`, `"random"`, or a numeric matrix with one row per sample
#'        and `n_axes` columns. See "Starting coordinates" below.
#'
#' @param seed   An integer seed for the random parts of the method: the
#'        negative samples, random starting coordinates, and the small
#'        random jitter added to spectral ones. The same seed gives the same
#'        result. Default: `0`
#'
#' @return A numeric matrix with one row per sample and one column per
#'   dimension, named `"UMAP1"`, `"UMAP2"`, etc., as described in "Ordination
#'   result" below. `method` is `"umap"`; `eig` and `prop` are `NA`, as UMAP
#'   has no eigenvalues. Then come:
#'
#'   * `a`, `b`: The parameters of the curve
#'     \eqn{1 / (1 + a d^{2b})} that turns plot distances \eqn{d} into
#'     similarities, fitted from `min_dist` and `spread`.
#'   * `init`: How the starting coordinates were chosen: `"spectral"`,
#'     `"pcoa"`, `"random"`, or `"matrix"` (given in `init`); `"none"` for
#'     fewer than 2 samples, which are placed at 0.
#'   * `n_epochs`: The number of optimization epochs run.
#'   * `n_neighbors`: The `n_neighbors` given.
#'
#' @details
#'
#' **Method**
#'
#' UMAP (McInnes et al. 2018) builds a weighted graph of each sample's
#' nearest neighbours, then finds plot coordinates whose own neighbour graph
#' is as similar to it as possible.
#'
#' 1. **Neighbours.** Each sample's `n_neighbors` nearest neighbours,
#'    counting itself, are found exactly from `dm`. Ties go to the sample
#'    that comes first.
#' 2. **Graph.** A neighbour \eqn{j} of sample \eqn{i} at distance
#'    \eqn{d_{ij}} gets weight \eqn{\exp(-(d_{ij} - \rho_i) / \sigma_i)},
#'    where \eqn{\rho_i} is the distance to the nearest neighbour (other than
#'    duplicates) and \eqn{\sigma_i} is chosen by binary search so that the
#'    weights sum to \eqn{\log_2(}`n_neighbors`\eqn{)}. The two weights
#'    \eqn{w_{ij}} and \eqn{w_{ji}} of each pair are combined by fuzzy union,
#'    \eqn{w_{ij} + w_{ji} - w_{ij} w_{ji}}.
#' 3. **Curve.** `a` and `b` are fitted by least squares (Levenberg-Marquardt,
#'    as umap-learn fits them) to a curve that is 1 up to `min_dist` and decays
#'    exponentially, with scale `spread`, beyond it.
#' 4. **Layout.** Starting from `init`, stochastic gradient descent pulls
#'    the two samples of each edge together, sampling each edge in proportion
#'    to its weight, and pushes each sample away from `negative_sample_rate`
#'    randomly chosen samples per edge sampled. As in `uwot::umap()`'s batch
#'    mode, each epoch computes every sample's gradient from the coordinates
#'    at the start of the epoch, then moves all samples at once by the Adam
#'    optimizer (\eqn{\beta_1 = 0.5}, \eqn{\beta_2 = 0.9},
#'    \eqn{\epsilon = 10^{-7}}, uwot's settings for batch mode). Finally, the
#'    coordinates are centered.
#'
#'
#' **Why batch mode**
#'
#' `uwot::umap()`'s default optimization moves samples one edge at a time,
#' each move seeing the ones before it. That can only run on one thread if
#' it is to be reproducible, and it takes nearly all of the run time. In
#' batch mode, the samples' gradients within an epoch don't depend on each
#' other, so `umap()` computes them on `cpus` threads, and still gets the
#' same result for any number of threads. The layouts differ from those of
#' uwot's default mode only as different random draws would: neighbour
#' preservation, trustworthiness, and the final cross-entropy are about the
#' same.
#'
#'
#' **Starting coordinates**
#'
#' * `"spectral"`: The eigenvectors of the graph's normalized Laplacian for
#'   its 2nd to \eqn{(k+1)}th smallest eigenvalues, scaled so that the
#'   largest absolute coordinate is 10, plus normal noise with standard
#'   deviation 0.0001. When the graph is not connected, the `"pcoa"` start is
#'   used instead, scaled to standard deviation 1, and with fewer than 3
#'   samples (or \eqn{k > n - 2}), the `"random"` one, as in `uwot::umap()`.
#' * `"pcoa"`: The first `n_axes` axes of [pcoa()], each centered and scaled to
#'   standard deviation 0.0001 (`uwot::umap()`'s `"spca"`).
#' * `"random"`: Coordinates uniform on \[-10, 10).
#' * A matrix: Used as is.
#'
#'
#' **Differences from `uwot::umap()`**
#'
#' * `umap()` always optimizes in batch mode, so its layouts are those of
#'   `uwot::umap(batch = TRUE)`, not of uwot's default.
#' * `min_dist` defaults to `0.1`, the canonical value used by umap-learn and
#'   most other implementations, rather than uwot's `0.01`. Other defaults
#'   are the same.
#' * The random numbers come from `seed`, not from R's random number
#'   generator, so `set.seed()` has no effect, and R's random number stream
#'   is left untouched. In batch mode, `uwot::umap()` draws from R's
#'   generator two 32-bit seeds per sample per epoch, for each sample's
#'   stream of negative samples, and the jitter of spectral starts. `umap()`
#'   makes those from `seed` instead, with the same streams of negative
#'   samples (PCG) from the same seeds; given the same seeds and starting
#'   coordinates, the two give the same coordinates (to the last bit, for
#'   the same compiler and math library).
#' * The spectral start's eigenvectors are computed by block Lanczos
#'   iteration, the same solver as [pcoa()]'s, to relative residuals of
#'   \eqn{10^{-14}} (at worst \eqn{10^{-10}}) rather than `RSpectra`'s
#'   tolerance of \eqn{10^{-4}}, and each one's sign is fixed (its largest
#'   absolute coordinate is positive).
#' * The `"pcoa"` start, and the fallback for unconnected graphs, come from
#'   [pcoa()] rather than `cmdscale()`; they differ only by rounding and the
#'   signs of the axes.
#' * Finally, each axis is flipped so that its largest absolute coordinate
#'   is positive.
#'
#'
#' **Speed and memory**
#'
#' `uwot::umap()` given a `dist` object copies it into an \eqn{n \times n}
#' matrix and sorts each of its columns in R. `umap()` reads the packed
#' triangle of `dm` directly, keeping only the `n_neighbors` closest of each
#' sample. Beyond `dm` itself, memory use is proportional to \eqn{n} times
#' `n_neighbors`. `cpus` threads find the neighbours, build the graph,
#' compute the spectral start, and compute each epoch's gradients.
#'
#'
#' **Reproducibility**
#'
#' Results are identical from run to run and for any number of `cpus`.
#'
#' @references
#' McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold
#' approximation and projection for dimension reduction. *arXiv*,
#' 1802.03426. \doi{10.48550/arXiv.1802.03426}
#'
#' Melville, J. (2020). uwot: The Uniform Manifold Approximation and Projection
#' (UMAP) Method for Dimensionality Reduction. R package. Its C++ header
#' library (BSD 2-Clause License) is the basis of the graph and layout code.
#' \url{https://github.com/jlmelville/uwot}
#'
#' Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
#' symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
#' Applications*, 22(2), 602-616. \doi{10.1137/S0895479898334605}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [tsne()], [pcoa()], [nmds()], [dist_ord()], [beta_div()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- umap(dm)
#'     plot(ord, 'Species', iris)
#'     ord
#'     str(ord)
#'
umap <- function (
    dm, n_axes = 2L, n_neighbors = 15L, min_dist = 0.1, spread = 1,
    n_epochs = NULL, learning_rate = 1, negative_sample_rate = 5,
    repulsion_strength = 1, set_op_mix_ratio = 1, local_connectivity = 1,
    bandwidth = 1, init = 'spectral', seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  k           <- ord_k(n_axes)
  n_neighbors <- ord_arg(n_neighbors, '`n_neighbors` must be a whole number of at least 2.', whole = TRUE, lo = 2)
  min_dist    <- ord_arg(min_dist, '`min_dist` must be a number of at least 0.')
  spread      <- ord_arg(spread, '`spread` must be a number greater than 0.', lo = .Machine$double.xmin)
  learning_rate <- ord_arg(learning_rate, '`learning_rate` must be a number greater than 0.', lo = .Machine$double.xmin)
  negative_sample_rate <- ord_arg(negative_sample_rate, '`negative_sample_rate` must be a number greater than 0.', lo = .Machine$double.xmin)
  repulsion_strength   <- ord_arg(repulsion_strength, '`repulsion_strength` must be a number of at least 0.')
  set_op_mix_ratio     <- ord_arg(set_op_mix_ratio, '`set_op_mix_ratio` must be a number from 0 to 1.', hi = 1)
  local_connectivity   <- ord_arg(local_connectivity, '`local_connectivity` must be a number of at least 1.', lo = 1)
  bandwidth            <- ord_arg(bandwidth, '`bandwidth` must be a number greater than 0.', lo = .Machine$double.xmin)

  if (!is.null(n_epochs))
    n_epochs <- ord_arg(n_epochs, '`n_epochs` must be NULL or a whole number of at least 0.', whole = TRUE)

  ia       <- ord_init_arg(init, c('spectral', 'pcoa', 'random'))
  init     <- ia$init
  init_mtx <- ia$mtx


  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  # A constant column would stay constant: its gradients are all zero.
  if (!is.null(init_mtx)) {
    init_mtx <- ord_init_matrix(init_mtx, n, k)
    if (n >= 2 && any(apply(init_mtx, 2L, function (v) all(v == v[[1]]))))
      stop('`init` must not have a column whose values are all the same.', call. = FALSE)
  }

  # Coordinates (n x k) and the graph's edges (up to 2n(n_neighbors - 1))
  # are indexed with ints in C.
  if (as.double(n) * k > .Machine$integer.max)
    stop('`n_axes` is too large: ', n, ' samples times `n_axes` dimensions must be ',
         'at most ', .Machine$integer.max, '.', call. = FALSE)
  if (2 * as.double(n) * (min(n_neighbors, n) - 1) > .Machine$integer.max)
    stop('`n_neighbors` is too large for ', n, ' samples.', call. = FALSE)

  if (is.null(n_epochs))
    n_epochs <- if (n <= 10000) 500L else 200L

  ab <- umap_ab(spread, min_dist)

  result <- function (x, init) {
    x <- unname(x)
    if (!all(is.finite(x)))
      stop(
        'UMAP diverged: coordinates became infinite or NaN. ',
        'Try a smaller `learning_rate` or `repulsion_strength`.', call. = FALSE )
    if (n > 0) x <- ord_flip(scale(x, center = TRUE, scale = FALSE))
    ord_result(
      matrix(x, n, k), 'umap', labels, prefix = 'UMAP', from = dv$title,
      a        = ab[['a']],
      b        = ab[['b']],
      init        = init,
      n_epochs    = n_epochs,
      n_neighbors = n_neighbors )
  }

  if (n < 2) {
    # One sample or none: nothing to arrange.
    x <- if (is.null(init_mtx)) matrix(0, n, k) else init_mtx
    return (result(x, if (init == 'matrix') init else 'none'))
  }

  n_neighbors <- min(n_neighbors, n)


  #--------------------------------------------------
  # Neighbours and fuzzy graph.
  #--------------------------------------------------
  graph <- .Call(
    C_umap_graph, dm, n, n_neighbors,
    c(local_connectivity, bandwidth, set_op_mix_ratio), cpus )


  #--------------------------------------------------
  # Starting coordinates.
  #--------------------------------------------------
  if (init == 'spectral') {
    if (n < 3 || k > n - 2) {
      init <- 'random'
    }
    else if (graph$n_components > 1) {
      x    <- umap_pcoa_init(dm, n, k, 1, cpus)
      init <- if (is.null(x)) 'random' else 'pcoa'
    }
    else {
      eig <- .Call(C_umap_spectral, graph, k, cpus)
      if (!eig$converged) # nocov start
        warning('Spectral initialization did not converge; coordinates may be inexact.', call. = FALSE)
      # nocov end
      x <- eig$vectors
      x <- x * (10 / max(abs(x))) + 1e-4 * .Call(C_umap_rand, n * k, 1L, seed)
    }
  }
  else if (init == 'pcoa') {
    x <- umap_pcoa_init(dm, n, k, 1e-4, cpus)
    if (is.null(x)) init <- 'random'
  }
  else if (init == 'matrix') {
    x <- init_mtx
  }

  if (init == 'random')
    x <- matrix(.Call(C_umap_rand, n * k, 0L, seed), n, k)


  #--------------------------------------------------
  # Layout: uwot's batch mode, with Adam's defaults
  # there (beta1 = 0.5, beta2 = 0.9, eps = 1e-7).
  #--------------------------------------------------
  if (n_epochs > 0)
    x <- .Call(
      C_umap_sgd, graph, x, n_epochs,
      c(ab[['a']], ab[['b']], repulsion_strength, learning_rate,
        negative_sample_rate, 0.5, 0.9, 1e-7),
      seed, cpus )

  result(x, init)
}



# The a and b parameters of the curve 1 / (1 + a * d^(2b)): a least-squares
# fit, in C, to 1 below min_dist and exp(-(d - min_dist) / spread) beyond,
# the target defined by McInnes et al. (2018).
umap_ab <- function (spread, min_dist) {
  ab <- .Call(C_umap_ab, as.double(spread), as.double(min_dist))
  c(a = ab[[1]], b = ab[[2]])
}



# PCoA starting coordinates, each axis centered and scaled to standard
# deviation `sdev`, as uwot's "spca". NULL if there are fewer than k axes.
umap_pcoa_init <- function (dm, n, k, sdev, cpus) {
  x <- pcoa_start(dm, n, k, cpus)
  if (is.null(x)) return (NULL)
  scale(x, scale = apply(x, 2, stats::sd) / sdev)
}
