# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Ported in part from scikit-learn (BSD-3-Clause License, Copyright (c)
# 2007-2026 The scikit-learn developers); see src/tsne.c and
# inst/COPYRIGHTS for the license text.



#' t-distributed Stochastic Neighbor Embedding (t-SNE)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing samples that are near neighbors in `dm` near each other on the
#' plot. t-SNE is good at showing clusters, but the distances between
#' clusters, and their sizes, mean little. A port of scikit-learn's `TSNE`
#' (with `metric = "precomputed"`) to parallel C, with exact and Barnes-Hut
#' gradients, that never copies `dm` into an \eqn{n \times n} matrix.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values. Its scale doesn't matter (see "Scale").
#'
#' @param n_axes   The number of dimensions. Default: `2`
#'
#' @param perplexity   Roughly, how many neighbors each sample's
#'        neighborhood counts. Must be at least 1; values above
#'        \eqn{(n - 1) / 3} are reduced to that (see "Perplexity and small
#'        datasets" below). Typical values are 5 to 50. Default: `30`
#'
#' @param theta   Accuracy of the Barnes-Hut gradient, from 0 (exact, and
#'        slow) to 1 (fast, and rough): a group of samples is treated as one
#'        point when its cell's width, divided by its distance, is below
#'        `theta`. scikit-learn calls it `angle`. Not used by the exact
#'        gradient. Default: `0.5`
#'
#' @param gradient   How to compute the gradient: `"exact"`, `"barnes_hut"`,
#'        or `"auto"` (the default), which uses `"exact"` for up to 500
#'        samples, or when `k > 3`, and `"barnes_hut"` otherwise.
#'        (scikit-learn calls it `method`.)
#'
#' @param init   The starting layout: `"pcoa"` (the default), `"random"`, or
#'        a numeric matrix with one row per sample, in the same order as
#'        `dm`, and `n_axes` columns, none of them constant (a constant column
#'        could never move). See "Starting layout" below.
#'
#' @param max_iters   The most iterations to run, including the first 250 with
#'        early exaggeration. The run stops sooner when it stops improving
#'        (see "Method"). Default: `1000`
#'
#' @param early_exaggeration   How much to multiply the affinities by during
#'        the first 250 iterations, which pulls clusters apart early on.
#'        Must be at least 1. Default: `12`
#'
#' @param learning_rate   The step size: a positive number, or `"auto"` (the
#'        default) for \eqn{\max(n / early\_exaggeration / 4, 50)}, as in
#'        scikit-learn.
#'
#' @param seed   An integer seed for the random starting layout. The same
#'        seed gives the same result. Not used when `init` is a matrix, nor
#'        by a PCoA start with `n_axes` axes. Default: `0`
#'
#' @return A numeric matrix with one row per sample and one column per
#'   dimension, named `"TSNE1"`, `"TSNE2"`, etc., as described in "Ordination
#'   result" below. `method` is `"tsne"`; `eig` and `prop` are `NA`, as t-SNE
#'   has no eigenvalues. Then come:
#'
#'   * `kl_divergence`: The Kullback-Leibler divergence of the result's
#'     affinities from the input affinities: the quantity t-SNE minimizes.
#'     With the Barnes-Hut gradient, it is the Barnes-Hut estimate.
#'   * `perplexity`: The perplexity used, after any reduction.
#'   * `init`: How the starting layout was chosen: `"pcoa"`, `"random"`, or
#'     `"matrix"` (given in `init`); `"none"` for fewer than 2 samples, which
#'     are all placed at 0.
#'   * `n_iters`: The number of iterations run.
#'   * `learning_rate`: The learning rate used.
#'   * `gradient`: `"exact"` or `"barnes_hut"`.
#'
#' @details
#'
#' **Method**
#'
#' Each sample \eqn{i} gets a Gaussian distribution over the other samples,
#' \eqn{p_{j|i} \propto \exp(-\beta_i d_{ij}^2)}, with \eqn{\beta_i} chosen by
#' binary search so that its perplexity (2 to the power of its entropy, in
#' bits) equals `perplexity`. The input affinities are
#' \eqn{p_{ij} = (p_{j|i} + p_{i|j}) / 2n}. On the plot, the affinities are
#' \eqn{q_{ij} \propto (1 + |y_i - y_j|^2)^{-1}}, a Student t distribution
#' with one degree of freedom, whose heavy tail lets dissimilar samples sit
#' far apart. t-SNE moves the points \eqn{y_i} to minimize the
#' Kullback-Leibler divergence \eqn{\sum p_{ij} \log(p_{ij} / q_{ij})} by
#' gradient descent, with momentum and per-coordinate gains: momentum 0.5
#' for the first 250 iterations, during which \eqn{p_{ij}} is multiplied by
#' `early_exaggeration`, then 0.8. Every 50 iterations, it stops if the
#' gradient's norm is below \eqn{10^{-7}}, or if the divergence hasn't
#' improved for 300 iterations (after the first 250). These are
#' scikit-learn's defaults.
#'
#'
#' **Exact and Barnes-Hut gradients**
#'
#' The exact gradient sums over every pair of samples, which takes time
#' proportional to \eqn{n^2} per iteration. The Barnes-Hut approximation
#' (van der Maaten 2014) takes time proportional to \eqn{n \log n}:
#'
#' * Each sample's \eqn{p_{j|i}} is computed over its
#'   \eqn{\lfloor 3 \times perplexity + 1 \rfloor} nearest neighbors only
#'   (found exactly from `dm`), and is 0 for the rest.
#' * The repulsion between points on the plot is computed from a tree of
#'   cells: a quadtree for `n_axes = 2`, an octree for `n_axes = 3`, and a binary tree
#'   for `n_axes = 1`. A cell far enough away (see `theta`) acts as one point at
#'   its center of mass.
#'
#' The tree only exists for up to 3 dimensions, so larger `n_axes` always uses
#' the exact gradient. `gradient = "auto"` uses the exact gradient for up to
#' 500 samples, where it is fast enough, and Barnes-Hut for more.
#'
#'
#' **Perplexity and small datasets**
#'
#' Each sample needs at least \eqn{3 \times perplexity} neighbors, so
#' `perplexity` is reduced to \eqn{(n - 1) / 3} when it is larger (but not
#' below 1), with a warning if you set `perplexity` yourself. For example,
#' 40 samples allow a perplexity of up to 13. `n_axes` is used as given, even
#' when it is \eqn{n} or more.
#'
#'
#' **Scale**
#'
#' t-SNE's affinities don't depend on the scale of `dm`, but finding them
#' does: each sample's binary search (as in scikit-learn) starts at a
#' precision of 1 and doubles or halves it at most 100 times, so it can't
#' reach the precision that distances far from 1 need, and their affinities
#' would come out uniform. So `dm` is first multiplied by a power of 2
#' (which is exact) that brings its largest distance to between 1 and 2;
#' multiplying `dm` by any power of 2 gives exactly the same result, and by
#' any other number, the same up to the search's tolerance. A sample whose
#' search still doesn't converge in those 100 steps (e.g. one whose nearest
#' neighbors are \eqn{10^{-20}} times closer than the farthest sample) is
#' searched again, starting from a better guess and with no limit on the
#' range.
#'
#'
#' **Starting layout**
#'
#' t-SNE only finds a local minimum, and which one depends on the starting
#' layout. The default, `"pcoa"`, starts from the first `n_axes` axes of
#' [pcoa()], scaled so that the first has a standard deviation of
#' \eqn{10^{-4}}, as scikit-learn scales its PCA start: this keeps more of
#' the data's global structure than a random start, and doesn't depend on a
#' seed (Kobak & Linderman 2021). Without `n_axes` axes with positive
#' eigenvalues, it starts at random instead. `"random"` draws each
#' coordinate from a normal distribution with a standard deviation of
#' \eqn{10^{-4}}, from `seed`.
#'
#'
#' **Speed and memory**
#'
#' Finding the affinities reads each distance twice (once for each of its
#' two samples), in parallel. Each sample's binary search then makes up to
#' 100 passes over its \eqn{3 \times perplexity + 1} neighbors (Barnes-Hut)
#' or all \eqn{n - 1} other samples (exact). Each iteration takes time
#' proportional to \eqn{k n^2} (exact) or roughly
#' \eqn{k n (\log n + 3 \times perplexity)} (Barnes-Hut), divided among
#' `cpus` threads. The affinities are found in batches, between which the
#' computation can be interrupted.
#'
#' Besides `dm`, the exact gradient keeps \eqn{p_{ij}} for every pair, 8
#' bytes per pair, the same as `dm` itself. Barnes-Hut keeps
#' \eqn{n (3 \times perplexity + 1)} to twice that many affinities, 12 bytes
#' each (and, while building them, two more lists of
#' \eqn{n (3 \times perplexity + 1)}), and a tree of at most \eqn{2n}
#' cells. While reading distances for the exact gradient, each thread also
#' holds 32 rows of them (\eqn{256 n} bytes).
#'
#'
#' **Differences from scikit-learn**
#'
#' This is a port of scikit-learn 1.9.1's `TSNE` (BSD-3-Clause License,
#' Copyright (c) 2007-2026 The scikit-learn developers). From the same
#' affinities and starting layout, the exact gradient gives the same
#' iterations as scikit-learn's to within rounding; and, as in
#' scikit-learn, `dm` is squared before the affinities are found. These are
#' the differences:
#'
#' * Everything is computed in double precision. scikit-learn rounds the
#'   squared distances to single precision before finding the affinities,
#'   and runs Barnes-Hut in single precision, which changes affinities and
#'   gradients by about 1 part in \eqn{10^7}. Since t-SNE is chaotic, two
#'   runs that differ by that much end at different (but equally good)
#'   layouts.
#' * \eqn{q_{ij}} always has one degree of freedom. scikit-learn uses
#'   \eqn{\max(k - 1, 1)}, i.e. 2 for `n_axes = 3`, with a gradient that isn't
#'   exactly the divergence's for more than one degree of freedom. For
#'   `n_axes = 1` and `n_axes = 2`, the two agree.
#' * The default start is PCoA. scikit-learn doesn't allow its PCA start with
#'   precomputed distances, so starts at random there. Random starts come
#'   from `seed` (see "Reproducibility"), so they differ from scikit-learn's.
#' * `gradient = "auto"` uses the exact gradient for up to 500 samples;
#'   scikit-learn defaults to Barnes-Hut.
#' * `dm` is rescaled by a power of 2, and a binary search that doesn't
#'   converge is redone (see "Scale"). Where scikit-learn's search does
#'   converge, the affinities are the same to within its tolerance
#'   (identical if the largest distance is already between 1 and 2).
#' * A `perplexity` too large for \eqn{n} is reduced instead of being an
#'   error; it must be at least 1. Any `max_iters` of at least 1 is allowed
#'   (scikit-learn requires at least 250). A layout that diverges (e.g.
#'   from a far too large `learning_rate`) is an error.
#' * The Barnes-Hut tree handles identical points exactly. scikit-learn's
#'   merges points within \eqn{10^{-6}} of each other on every axis and
#'   leaves out their repulsion from each other. A cell that holds the point
#'   being moved is never summarized (scikit-learn's can be, when `theta` is
#'   above \eqn{1/\sqrt{2}} for `n_axes = 2`, or \eqn{1/\sqrt{3}} for `n_axes = 3`).
#' * Nearest neighbors with equal distances are chosen by sample order.
#' * The gradient and the divergence leave out scikit-learn's lower limits on
#'   \eqn{q_{ij}} (\eqn{2.2 \times 10^{-16}} for the exact gradient, and
#'   single precision's smallest number for Barnes-Hut). They only matter
#'   for points about \eqn{10^8} times farther apart than usual.
#' * The result is centered, and each axis is flipped so that its
#'   coordinate with the largest absolute value is positive, as for the
#'   other ordination methods. `kl_divergence` is that of the result, rather
#'   than of the layout before the last step, and `n_iters` counts the
#'   iterations run (scikit-learn's `n_iter_` is one less).
#'
#'
#' **Compared with `Rtsne()`**
#'
#' `Rtsne::Rtsne()` (van der Maaten's bhtsne) runs the same kind of
#' optimization with different defaults: a random start, a learning rate of
#' 200, Barnes-Hut whenever `theta > 0`, \eqn{\lfloor 3 \times perplexity
#' \rfloor} neighbors found with a vantage-point tree, and always exactly
#' `max_iters` iterations. `init = "random", learning_rate = 200` comes close
#' to its procedure. The two don't give the same coordinates, but layouts of
#' the same quality: the tests check that the KL divergence and
#' trustworthiness are as good as `Rtsne()`'s on the same data.
#'
#'
#' **Reproducibility**
#'
#' The result, including its attributes, is bit-for-bit identical from run
#' to run, in the same or a new R session, and for any number of `cpus`:
#' each quantity is computed by one thread in a fixed order, sums over
#' samples are taken in sample order, and all of it is in double precision
#' (not R's long double, whose precision varies between platforms). Random
#' starts come from a PCG random number generator seeded with `seed`, not
#' from R's, so `set.seed()` has no effect and R's random number stream is
#' left untouched. Different compilers, processors, or math libraries can
#' change the last bits of a result, which t-SNE amplifies, so other
#' platforms can give a different (but equally good) layout.
#'
#' @references
#' van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE.
#' *Journal of Machine Learning Research*, 9(86), 2579-2605.
#' <https://jmlr.org/papers/v9/vandermaaten08a.html>
#'
#' van der Maaten, L. (2014). Accelerating t-SNE using tree-based
#' algorithms. *Journal of Machine Learning Research*, 15(93), 3221-3245.
#' <https://jmlr.org/papers/v15/vandermaaten14a.html>
#'
#' Belkina, A. C., Ciccolella, C. O., Anno, R., Halpert, R., Spidlen, J., &
#' Snyder-Cappione, J. E. (2019). Automated optimized parameters for
#' T-distributed stochastic neighbor embedding improve visualization and
#' analysis of large datasets. *Nature Communications*, 10, 5415.
#' \doi{10.1038/s41467-019-13055-y}
#'
#' Kobak, D., & Linderman, G. C. (2021). Initialization is critical for
#' preserving global data structure in both t-SNE and UMAP. *Nature
#' Biotechnology*, 39, 156-157. \doi{10.1038/s41587-020-00809-z}
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830.
#' <https://jmlr.org/papers/v12/pedregosa11a.html>
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [umap()], [pcoa()], [nmds()], [dist_ord()], [beta_div()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- tsne(dm)
#'     plot(ord, 'Species', iris)
#'     head(ord)
#'     str(ord)
#'
tsne <- function (
    dm, n_axes = 2L, perplexity = 30, theta = 0.5, gradient = 'auto',
    init = 'pcoa', max_iters = 1000L, early_exaggeration = 12,
    learning_rate = 'auto', seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  k          <- ord_k(n_axes)
  perp_given <- !missing(perplexity)
  perplexity <- ord_arg(perplexity, '`perplexity` must be a number of at least 1.', lo = 1)
  theta      <- ord_arg(theta,      '`theta` must be a number from 0 to 1.', hi = 1)
  max_iters  <- ord_arg(max_iters,  '`max_iters` must be a whole number of at least 1.', whole = TRUE, lo = 1)
  early_exaggeration <- ord_arg(
    early_exaggeration, '`early_exaggeration` must be a number of at least 1.', lo = 1 )
  if (!identical(learning_rate, 'auto'))
    learning_rate <- ord_arg(
      learning_rate, '`learning_rate` must be "auto" or a number greater than 0.',
      lo = .Machine$double.xmin )

  choices <- c('auto', 'exact', 'barnes_hut')
  if (!is.character(gradient) || length(gradient) != 1 || !gradient %in% choices)
    stop('`gradient` must be "auto", "exact", or "barnes_hut".', call. = FALSE)
  if (gradient == 'barnes_hut' && k > 3)
    stop('`gradient = "barnes_hut"` supports at most `n_axes = 3` dimensions.', call. = FALSE)

  ia       <- ord_init_arg(init, c('pcoa', 'random'))
  init     <- ia$init
  init_mtx <- ia$mtx


  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  # Coordinates (n x k) are indexed with ints in C.
  if (as.double(n) * k > .Machine$integer.max)
    stop('`n_axes` is too large: ', n, ' samples times `n_axes` dimensions must be ',
         'at most ', .Machine$integer.max, '.', call. = FALSE)

  if (!is.null(init_mtx)) {
    init_mtx <- ord_init_matrix(init_mtx, n, k)
    if (n >= 2 && any(apply(init_mtx, 2L, function (v) all(v == v[[1]]))))
      stop('`init` must not have a column whose values are all the same.', call. = FALSE)
  }


  # Each sample needs at least 3 * perplexity neighbors.
  perp_max <- max(1, (n - 1) / 3)
  if (perplexity > perp_max) {
    if (perp_given)
      warning(sprintf(
        '`perplexity` reduced from %g to %g, the most that %d samples allow.',
        perplexity, perp_max, n ), call. = FALSE)
    perplexity <- perp_max
  }

  exact <- switch(
    EXPR       = gradient,
    exact      = TRUE,
    barnes_hut = FALSE,
    auto       = k > 3 || n <= TSNE_EXACT_MAX )

  if (identical(learning_rate, 'auto'))
    learning_rate <- max(n / early_exaggeration / 4, 50)


  result <- function (x, kl_divergence, n_iter) {
    ord_result(
      ord_flip(x), 'tsne', labels, prefix = 'TSNE', from = dv$title,
      kl_divergence = kl_divergence,
      perplexity    = perplexity,
      init          = init,
      n_iters       = n_iter,
      learning_rate = learning_rate,
      gradient      = if (exact) 'exact' else 'barnes_hut' )
  }

  if (n < 2) {
    # One sample or none: nothing to embed.
    if (init != 'matrix') init <- 'none'
    return (result(matrix(0, n, k), 0, 0L))
  }


  # t-SNE doesn't depend on the scale of `dm`, but its arithmetic does: the
  # affinities' binary search starts at a precision of 1. So distances are
  # multiplied by a power of 2 (exactly) to make the largest one from 1 to 2.
  top   <- max(dm)
  e     <- if (top > 0) floor(log2(top)) else 0
  scale <- 2^-min(max(e, -1000), 1000)


  # Starting layout: 0 = given, 1 = PCoA (scaled in C), 2 = random.
  start <- 0L
  if (init == 'pcoa') {
    # Far from 1, PCoA's squared distances could underflow; the start is
    # rescaled anyway, so scale a copy first.
    dp       <- if (scale < 2^-64 || scale > 2^64) dm * scale else dm
    init_mtx <- pcoa_start(dp, n, k, cpus)
    rm(dp)
    if (is.null(init_mtx)) init <- 'random'
    else                   start <- 1L
  }
  if (init == 'random') start <- 2L


  # Input affinities, then the optimization.
  K <- if (exact) 0L else as.integer(min(n - 1, floor(3 * perplexity + 1)))
  P <- .Call(C_tsne_p, dm, n, perplexity, K, FALSE, scale, cpus)
  rm(dm)

  iopt <- as.integer(c(n, k, max_iters, 250L, cpus, 50L, 300L, 1L, start, seed))
  dopt <- c(theta, early_exaggeration, learning_rate, 1e-7, 0.5, 0.8)
  res  <- .Call(C_tsne_run, P, init_mtx, iopt, dopt)

  if (res$diverged)
    stop(
      't-SNE diverged: coordinates became infinite or NaN. ',
      'Try a smaller `learning_rate` or `early_exaggeration`.', call. = FALSE )

  result(res$points, res$kl, res$n_iter)
}



# Largest number of samples that gradient = "auto" runs with the exact
# gradient.
TSNE_EXACT_MAX <- 500L
