# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' t-distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing samples that are near each other in `dm` near each other on the
#' plot. Follows the same procedure as `Rtsne::Rtsne()` given a distance
#' matrix, and from the same starting coordinates gives the same Barnes-Hut
#' embedding, but runs in parallel C with far less memory.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param k   The number of dimensions. Must be 1, 2, or 3 unless
#'        `theta = 0`. Default: `2`
#'
#' @param perplexity   Roughly, how many close neighbors each sample has.
#'        Must be at most \eqn{(n - 1) / 3} for \eqn{n} samples. Default:
#'        `30`, or \eqn{(n - 1) / 3} if that is smaller.
#'
#' @param theta   Speed/accuracy trade-off for the Barnes-Hut
#'        approximation, from `0` to `1`. `0` computes t-SNE exactly, which
#'        takes time proportional to \eqn{n^2} per iteration. Larger values
#'        approximate it faster. Default: `0.5`
#'
#' @param maxit   The number of iterations. Default: `1000`
#'
#' @param init   The starting coordinates: `"pcoa"` (the default),
#'        `"random"`, or a numeric matrix with one row per sample, in the same
#'        order as `dm`, and `k` columns. See "Starting coordinates" below.
#'
#' @param learning_rate   The step size of the gradient descent (`eta` in
#'        `Rtsne()`). `"auto"`, the default, uses
#'        `max(n / exaggeration_factor, 200)` for `n` samples, as
#'        scikit-learn and FIt-SNE do: `Rtsne()`'s fixed `200` for up to 2400
#'        samples, and larger steps for more, which need them to converge in
#'        `maxit` iterations.
#'
#' @param momentum,final_momentum   The momentum of the gradient descent
#'        before and after iteration `mom_switch_iter`. Defaults: `0.5` and
#'        `0.8`
#'
#' @param mom_switch_iter   The iteration from which `final_momentum` is
#'        used. `NULL`, the default, means `250`, or `0` when `init` is a
#'        matrix.
#'
#' @param exaggeration_factor   How much the input similarities are
#'        multiplied by during the early iterations ("early exaggeration"),
#'        which helps clusters separate. Default: `12`
#'
#' @param stop_lying_iter   The iteration at which early exaggeration stops.
#'        `NULL`, the default, means `250`, or `0` (no exaggeration) when
#'        `init` is a matrix.
#'
#' @param seed   An integer seed for the random starting coordinates. The
#'        same seed gives the same result. Not used when `init` is `"pcoa"`
#'        or a matrix. Default: `0`
#'
#' @return A numeric matrix with one row per sample and one column per
#'   dimension, named `"TSNE1"`, `"TSNE2"`, etc., as described in "Ordination
#'   result" below. `method` is `"tsne"`; `eig` and `prop` are `NA`, as t-SNE
#'   has no eigenvalues. Then come:
#'
#'   * `kl_divergence`: The Kullback-Leibler divergence of the embedding's
#'     similarities from the input similarities at the end: the cost that
#'     t-SNE minimizes. As in `Rtsne()`, it is approximate when `theta > 0`.
#'   * `perplexity`: The perplexity used.
#'   * `init`: How the starting coordinates were chosen: `"random"`,
#'     `"pcoa"`, or `"matrix"` (given in `init`).
#'
#' @details
#'
#' **Method**
#'
#' t-SNE turns the distances \eqn{d_{ij}} into similarities: each sample
#' \eqn{i} gets a Gaussian kernel \eqn{p_{j|i} \propto \exp(-\beta_i
#' d_{ij}^2)}, whose \eqn{\beta_i} is found by binary search so that the
#' kernel's perplexity (2 to the power of its entropy in bits) equals
#' `perplexity`. Then \eqn{p_{ij} = (p_{j|i} + p_{i|j}) / 2n}. Samples are
#' placed in `k` dimensions so that the heavy-tailed similarities
#' \eqn{q_{ij} \propto (1 + |y_i - y_j|^2)^{-1}} between their coordinates
#' match \eqn{p_{ij}} as well as possible: the coordinates minimize the
#' Kullback-Leibler divergence \eqn{\sum p_{ij} \log(p_{ij} / q_{ij})} by
#' gradient descent with momentum and per-coordinate adaptive gains.
#'
#' With `theta > 0`, two Barnes-Hut approximations make this fast for large
#' \eqn{n} (van der Maaten 2014): \eqn{p_{j|i}} is computed only for each
#' sample's \eqn{\lfloor 3 \times} `perplexity`\eqn{\rfloor} nearest
#' neighbors (others are 0), and the repulsion between distant groups of
#' samples is computed from their centers of mass in a quadtree (an octree
#' for `k = 3`). A group is used whole when its cell's width divided by its
#' distance is below `theta`.
#'
#' t-SNE preserves neighborhoods, not distances: the sizes of clusters and
#' the distances between them on the plot are not meaningful.
#'
#'
#' **Starting coordinates**
#'
#' * `"pcoa"`: The first `k` axes of [pcoa()], scaled so that the first has
#'   standard deviation \eqn{10^{-4}}, as scikit-learn, openTSNE and FIt-SNE
#'   scale their PCA starts by default. This preserves more of the
#'   large-scale layout than a random start, and makes the result depend
#'   only on `dm`. If fewer than `k` axes have positive eigenvalues,
#'   `"random"` is used instead.
#' * `"random"`: Normal, with mean 0 and standard deviation \eqn{10^{-4}}, as
#'   in `Rtsne()`, drawn from `seed`.
#' * A matrix: used as given. As in `Rtsne()`, the defaults then skip early
#'   exaggeration and the low-momentum phase, as suits continuing from an
#'   earlier result.
#'
#'
#' **Speed and memory**
#'
#' `Rtsne()` copies a dist object into an \eqn{n \times n} matrix, and for
#' `theta = 0` keeps three more \eqn{n \times n} matrices. `tsne()` works
#' straight from `dm`. For `theta = 0` it keeps one array of \eqn{p_{ij}} the
#' size of `dm`; for `theta > 0`, only the nearest neighbors'. The neighbor
#' search, the calibration of each \eqn{\beta_i}, and each iteration's
#' gradient run on `cpus` threads.
#'
#'
#' **Differences from `Rtsne()`**
#'
#' * The defaults `init = "pcoa"` and `learning_rate = "auto"` follow
#'   scikit-learn (since version 1.2), openTSNE and FIt-SNE, which have
#'   moved on from `Rtsne()`'s random start and fixed learning rate. Use
#'   `init = "random", learning_rate = 200` for `Rtsne()`'s procedure.
#' * The random starting coordinates come from `seed`, not from R's random
#'   number generator, so `set.seed()` has no effect, and R's random number
#'   stream is left untouched. They have the same distribution as
#'   `Rtsne()`'s: normal with standard deviation \eqn{10^{-4}}.
#' * Nearest neighbors are found exactly, by scanning `dm`. `Rtsne()` finds
#'   them with a vantage-point tree, which can miss some when the distances
#'   aren't metric (e.g. Bray-Curtis), so its results can then differ. Ties
#'   between neighbors' distances are broken by sample order.
#' * `perplexity` too large for the number of samples is reduced (with a
#'   warning, if it was given), not an error.
#'
#' Otherwise every step is computed as `Rtsne()` computes it, in the same
#' order, with one exception: for `theta = 0`, the sum that normalizes
#' \eqn{q_{ij}} is added up in parallel, which changes its rounding. So from
#' the same starting coordinates (`init`; `Y_init` in `Rtsne()`), the
#' Barnes-Hut (`theta > 0`) coordinates and cost are exactly `Rtsne()`'s.
#' With `theta = 0` they agree closely at first, but t-SNE magnifies tiny
#' rounding differences over its iterations, so the two can settle in
#' different, equally good, local minima. The same can happen with
#' `theta > 0` on other hardware (e.g. with fused multiply-adds).
#'
#'
#' **Reproducibility**
#'
#' Results are identical from run to run and for any number of `cpus`. The
#' coordinates are centered on 0, but not otherwise rotated or flipped.
#'
#' @references
#' van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE.
#' *Journal of Machine Learning Research*, 9, 2579-2605.
#'
#' van der Maaten, L. (2014). Accelerating t-SNE using tree-based
#' algorithms. *Journal of Machine Learning Research*, 15, 3221-3245.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [nmds()], [beta_div()]
#'
#' @export
#' @examples
#'     # Four samples only allow a perplexity of 1.
#'     tsne(bray(ex_counts))
#'
#'     # A larger example: iris flowers, by Euclidean distance.
#'     flowers <- unique(iris)
#'     d       <- dist(flowers[, 1:4])
#'     ord     <- tsne(d, perplexity = 20)
#'     attr(ord, 'kl_divergence')
#'
#'     plot(ord, asp = 1, col = flowers$Species)
#'
#'     # Exact gradients, from a random start.
#'     ord <- tsne(d, perplexity = 20, init = 'random', theta = 0)
#'
tsne <- function (
    dm, k = 2L, perplexity = 30, theta = 0.5, maxit = 1000L, init = 'pcoa',
    learning_rate = 'auto', momentum = 0.5, final_momentum = 0.8,
    mom_switch_iter = NULL,
    exaggeration_factor = 12,
    stop_lying_iter = NULL,
    seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  init_mtx <- NULL
  if (is.matrix(init) || is.data.frame(init)) {
    init_mtx <- init
    init     <- 'matrix'
  }
  else if (!is.character(init) || length(init) != 1 || !init %in% c('random', 'pcoa')) {
    stop('`init` must be "random", "pcoa", or a numeric matrix.', call. = FALSE)
  }

  # As in Rtsne, a given start skips early exaggeration and low momentum.
  if (is.null(mom_switch_iter)) mom_switch_iter <- if (init == 'matrix') 0L else 250L
  if (is.null(stop_lying_iter)) stop_lying_iter <- if (init == 'matrix') 0L else 250L

  k               <- nmds_arg(k,               '`k` must be a whole number of at least 1.', whole = TRUE, lo = 1)
  perp_given      <- !missing(perplexity)
  perplexity      <- nmds_arg(perplexity,      '`perplexity` must be a number greater than 0.')
  theta           <- nmds_arg(theta,           '`theta` must be a number from 0 to 1.', hi = 1)
  maxit           <- nmds_arg(maxit,           '`maxit` must be a whole number of at least 0.', whole = TRUE)
  if (!identical(learning_rate, 'auto'))
    learning_rate <- nmds_arg(learning_rate, '`learning_rate` must be "auto" or a number greater than 0.')
  momentum        <- nmds_arg(momentum,        '`momentum` must be a number of at least 0.')
  final_momentum  <- nmds_arg(final_momentum,  '`final_momentum` must be a number of at least 0.')
  mom_switch_iter <- nmds_arg(mom_switch_iter, '`mom_switch_iter` must be a whole number of at least 0.', whole = TRUE)
  exaggeration_factor <- nmds_arg(exaggeration_factor, '`exaggeration_factor` must be a number greater than 0.')
  stop_lying_iter <- nmds_arg(stop_lying_iter, '`stop_lying_iter` must be a whole number of at least 0.', whole = TRUE)

  if (perplexity <= 0)
    stop('`perplexity` must be a number greater than 0.', call. = FALSE)
  if (is.numeric(learning_rate) && learning_rate <= 0)
    stop('`learning_rate` must be "auto" or a number greater than 0.', call. = FALSE)
  if (exaggeration_factor <= 0)
    stop('`exaggeration_factor` must be a number greater than 0.', call. = FALSE)
  if (theta > 0 && k > 3)
    stop('`k` must be 1, 2, or 3 when `theta` > 0; use `theta = 0` for more.', call. = FALSE)


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

  result <- function (x, kl_divergence) {
    dimnames(x) <- list(labels, paste0('TSNE', seq_len(k)))
    structure(
      .Data         = x,
      method        = 'tsne',
      eig           = rep(NA_real_, k),
      prop          = rep(NA_real_, k),
      kl_divergence = kl_divergence,
      perplexity    = perplexity,
      init          = init )
  }

  if (n < 2) {
    # One sample or none: nothing to embed.
    return (result(matrix(0, n, k), 0))
  }

  if (!all(is.finite(dm)))
    stop('`dm` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (any(dm < 0))
    stop('`dm` must not contain negative distances.', call. = FALSE)


  # Each point needs floor(3 * perplexity) neighbors.
  if (3 * perplexity > n - 1) {
    if (perp_given)
      warning(
        '`perplexity` reduced from ', perplexity, ' to ', signif((n - 1) / 3, 4),
        ', the most that ', n, ' samples allow.', call. = FALSE )
    perplexity <- (n - 1) / 3
  }


  if (identical(learning_rate, 'auto'))
    learning_rate <- max(n / exaggeration_factor, 200)

  if (init == 'pcoa') {
    init_mtx <- tsne_pcoa_init(dm, n, k, cpus)
    if (is.null(init_mtx)) init <- 'random'
  }

  res <- .Call(
    C_tsne, dm, n, k, init_mtx,
    as.integer(c(maxit, stop_lying_iter, mom_switch_iter, cpus, seed)),
    c(perplexity, theta, learning_rate, momentum, final_momentum, exaggeration_factor) )

  result(res$points, res$kl_divergence)
}



# PCoA starting coordinates, scaled so that the first axis has standard
# deviation 1e-4, as openTSNE and scikit-learn scale their PCA starts. NULL
# if there are fewer than k axes with positive eigenvalues.
tsne_pcoa_init <- function (dm, n, k, cpus) {

  if (k > n - 1) return (NULL)

  res <- .Call(C_pcoa, dm, n, k, cpus)
  eig <- res$eig
  if (length(eig) < k || !all(eig > max(0, sqrt(.Machine$double.eps) * eig[[1]])))
    return (NULL)

  x <- res$points
  x / stats::sd(x[, 1]) * 1e-4
}
