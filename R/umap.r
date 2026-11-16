# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing samples that are near neighbours in `dm` near each other on the
#' plot. Follows the same procedure, with the same defaults, as
#' `uwot::umap()` given a `dist` object, and gives the same layout from the
#' same starting coordinates and random seeds, but finds the neighbours in
#' parallel C, without ever copying `dm` into an \eqn{n \times n} matrix.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param k   The number of dimensions. Default: `2`
#'
#' @param n_neighbors   How many nearest neighbours of each sample (counting
#'        the sample itself) to use. Larger values preserve more of the
#'        global structure; smaller ones more of the local detail. Values
#'        larger than the number of samples are reduced to it. Default: `15`
#'
#' @param min_dist   How close samples can be packed together on the plot,
#'        relative to `spread`. Default: `0.01`
#'
#' @param spread   The scale of the plot's distances, together with
#'        `min_dist`. Default: `1`
#'
#' @param n_epochs   The number of optimization epochs. `NULL`, the
#'        default, uses 500 for up to 10,000 samples and 200 for more, as
#'        `uwot::umap()` does.
#'
#' @param learning_rate   The initial learning rate of the optimization. It
#'        decreases linearly to 0 over the epochs. Default: `1`
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
#'        and `k` columns. See "Starting coordinates" below.
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
#'     `"pcoa"`, `"random"`, or `"matrix"` (given in `init`).
#'   * `n_epochs`: The number of optimization epochs run.
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
#' 3. **Curve.** `a` and `b` are fitted by nonlinear least squares
#'    ([stats::nls()]) to a curve that is 1 up to `min_dist` and decays
#'    exponentially, with scale `spread`, beyond it.
#' 4. **Layout.** Starting from `init`, stochastic gradient descent pulls
#'    the two samples of each edge together, sampling each edge in proportion
#'    to its weight, and pushes each sample away from `negative_sample_rate`
#'    randomly chosen samples per edge sampled. Finally, the coordinates are
#'    centered.
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
#' * `"pcoa"`: The first `k` axes of [pcoa()], each centered and scaled to
#'   standard deviation 0.0001 (`uwot::umap()`'s `"spca"`).
#' * `"random"`: Coordinates uniform on \[-10, 10).
#' * A matrix: Used as is.
#'
#'
#' **Differences from `uwot::umap()`**
#'
#' * The random numbers come from `seed`, not from R's random number
#'   generator, so `set.seed()` has no effect, and R's random number stream
#'   is left untouched. `uwot::umap()` draws from R's generator a 32-bit seed
#'   for each epoch's stream of negative samples, and the jitter of spectral
#'   starts. `umap()` makes those from `seed` instead, with the same stream
#'   of negative samples (PCG) from the same per-epoch seeds; given the same
#'   seeds and starting coordinates, the two give the same coordinates (to
#'   the last bit, for the same compiler and math library).
#' * The spectral start's eigenvectors are computed to a tolerance of about
#'   \eqn{10^{-10}} by Lanczos iteration, rather than `RSpectra`'s
#'   \eqn{10^{-4}}, and each one's sign is fixed (its largest absolute
#'   coordinate is positive).
#' * The `"pcoa"` start, and the fallback for unconnected graphs, come from
#'   [pcoa()] rather than `cmdscale()`; they differ only by rounding and the
#'   signs of the axes.
#' * The optimization runs on one thread, as `uwot::umap()`'s does by
#'   default (`n_sgd_threads = 0`), so that the result is the same for any
#'   number of `cpus`. `cpus` threads find the neighbours, the graph, and
#'   the spectral start.
#' * Finally, each axis is flipped so that its largest absolute coordinate
#'   is positive.
#'
#'
#' **Speed and memory**
#'
#' `uwot::umap()` given a `dist` object copies it into an \eqn{n \times n}
#' matrix and sorts each of its columns in R. `umap()` reads the packed
#' triangle of `dm` directly, keeping only the `n_neighbors` closest of each
#' sample, in parallel. Beyond `dm` itself, memory use is proportional to
#' \eqn{n} times `n_neighbors`.
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
#' Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large
#' symmetric eigenvalue problems. *SIAM Journal on Matrix Analysis and
#' Applications*, 22(2), 602-616. \doi{10.1137/S0895479898334605}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [nmds()], [beta_div()]
#'
#' @export
#' @examples
#'     # 60 samples from three kinds of community.
#'     set.seed(1)
#'     rates  <- rbind(c(50, 20, 5, 1), c(5, 50, 20, 1), c(1, 5, 20, 50))
#'     counts <- matrix(rpois(240, rates[rep(1:3, each = 20), ]), 60, 4)
#'     rownames(counts) <- paste0('S', 1:60)
#'
#'     ord <- umap(bray(counts), n_neighbors = 10)
#'     head(ord)
#'
#'     plot(ord, col = rep(1:3, each = 20), asp = 1)
#'
umap <- function (
    dm, k = 2L, n_neighbors = 15L, min_dist = 0.01, spread = 1,
    n_epochs = NULL, learning_rate = 1, negative_sample_rate = 5,
    repulsion_strength = 1, set_op_mix_ratio = 1, local_connectivity = 1,
    bandwidth = 1, init = 'spectral', seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  k           <- nmds_arg(k, '`k` must be a whole number of at least 1.', whole = TRUE, lo = 1)
  n_neighbors <- nmds_arg(n_neighbors, '`n_neighbors` must be a whole number of at least 2.', whole = TRUE, lo = 2)
  min_dist    <- nmds_arg(min_dist, '`min_dist` must be a number of at least 0.')
  spread      <- nmds_arg(spread, '`spread` must be a number greater than 0.', lo = .Machine$double.xmin)
  learning_rate <- nmds_arg(learning_rate, '`learning_rate` must be a number greater than 0.', lo = .Machine$double.xmin)
  negative_sample_rate <- nmds_arg(negative_sample_rate, '`negative_sample_rate` must be a number greater than 0.', lo = .Machine$double.xmin)
  repulsion_strength   <- nmds_arg(repulsion_strength, '`repulsion_strength` must be a number of at least 0.')
  set_op_mix_ratio     <- nmds_arg(set_op_mix_ratio, '`set_op_mix_ratio` must be a number from 0 to 1.', hi = 1)
  local_connectivity   <- nmds_arg(local_connectivity, '`local_connectivity` must be a number of at least 1.', lo = 1)
  bandwidth            <- nmds_arg(bandwidth, '`bandwidth` must be a number greater than 0.', lo = .Machine$double.xmin)

  if (!is.null(n_epochs))
    n_epochs <- nmds_arg(n_epochs, '`n_epochs` must be NULL or a whole number of at least 0.', whole = TRUE)

  init_mtx <- NULL
  if (is.matrix(init) || is.data.frame(init)) {
    init_mtx <- init
    init     <- 'matrix'
  }
  else if (!is.character(init) || length(init) != 1 || !init %in% c('spectral', 'pcoa', 'random')) {
    stop('`init` must be "spectral", "pcoa", "random", or a numeric matrix.', call. = FALSE)
  }


  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  if (!is.null(init_mtx)) {
    init_mtx <- unname(as.matrix(init_mtx))
    if (!is.numeric(init_mtx) || !identical(dim(init_mtx), c(n, k)) || !all(is.finite(init_mtx)))
      stop('`init` must be a numeric matrix with ', n, ' rows and ', k,
           ' columns, without NA or infinite values.', call. = FALSE)
    storage.mode(init_mtx) <- 'double'
  }

  if (!all(is.finite(dm)))
    stop('`dm` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(dm < 0))
    stop('`dm` must not contain negative distances.', call. = FALSE)

  if (is.null(n_epochs))
    n_epochs <- if (n <= 10000) 500L else 200L

  ab <- umap_ab(spread, min_dist)

  result <- function (x, init) {
    x <- unname(x)
    if (n > 0) x <- scale(x, center = TRUE, scale = FALSE)
    for (j in seq_len(k)) {
      v <- x[, j]
      if (n > 0 && v[[which.max(abs(v))]] < 0) x[, j] <- -v
    }
    x <- matrix(as.double(x), n, k, dimnames = list(labels, paste0('UMAP', seq_len(k))))
    structure(
      .Data    = x,
      method   = 'umap',
      eig      = rep(NA_real_, k),
      prop     = rep(NA_real_, k),
      a        = ab[['a']],
      b        = ab[['b']],
      init     = init,
      n_epochs = n_epochs )
  }

  if (n < 2) {
    # One sample or none: nothing to arrange.
    x <- if (is.null(init_mtx)) matrix(0, n, k) else init_mtx
    return (result(x, init))
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
  # Layout.
  #--------------------------------------------------
  if (n_epochs > 0)
    x <- .Call(
      C_umap_sgd, graph, x, n_epochs,
      c(ab[['a']], ab[['b']], repulsion_strength, learning_rate, negative_sample_rate),
      seed )

  result(x, init)
}



# The a and b parameters of the curve 1 / (1 + a * d^(2b)), fitted to
# 1 for d < min_dist and exp(-(d - min_dist) / spread) beyond, at 300
# points from 0 to 3 * spread, as in uwot and umap-learn.
umap_ab <- function (spread, min_dist) {

  xv <- seq(from = 0, to = spread * 3, length.out = 300)
  yv <- rep(0, length(xv))
  yv[xv <  min_dist] <- 1
  yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)

  res <- tryCatch(
    stats::nls(yv ~ 1 / (1 + a * xv^(2 * b)), start = list(a = 1, b = 1))$m$getPars(),
    error = function (e) NULL )

  if (is.null(res))
    stop('Cannot fit a and b for spread = ', spread, ' and min_dist = ', min_dist, '.', call. = FALSE)

  res
}



# PCoA starting coordinates, each axis centered and scaled to standard
# deviation `sdev`, as uwot's "spca". NULL if there are fewer than k axes.
umap_pcoa_init <- function (dm, n, k, sdev, cpus) {

  if (k > n - 1) return (NULL)

  res <- .Call(C_pcoa, dm, n, k, cpus)
  eig <- res$eig
  if (length(eig) < k || !all(eig > max(0, sqrt(.Machine$double.eps) * eig[[1]])))
    return (NULL)

  x <- res$points
  scale(x, scale = apply(x, 2, stats::sd) / sdev)
}
