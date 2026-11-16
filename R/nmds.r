# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Non-metric Multidimensional Scaling (NMDS)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' placing them so that their distances on the plot are in the same rank order
#' as the distances in `dm`, as nearly as possible. Follows the same procedure
#' as `vegan::metaMDS()` with its default `monoMDS()` engine, and gives the
#' same solutions from the same starting configurations, but runs the random
#' starts in parallel C with far less memory.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param k   The number of dimensions. Default: `2`
#'
#' @param try,trymax   The minimum and maximum number of random starts.
#'        Random starts continue past `try` until the best solution is found
#'        again, or `trymax` is reached. Default: `20` and `20`
#'
#' @param maxit   The maximum number of iterations for each start.
#'        Default: `200`
#'
#' @param weakties   How to treat tied distances in `dm`. `TRUE`, the default,
#'        lets them have different distances on the plot (Kruskal's primary
#'        approach); `FALSE` requires their plot distances to be equal too
#'        (secondary approach).
#'
#' @param stress   Kruskal's stress formula to minimize: `1` (the default) or
#'        `2`.
#'
#' @param smin,sfgrmin,sratmax   Each start stops when its stress falls
#'        below `smin`, when the gradient's scale factor falls below
#'        `sfgrmin`, or when the stress ratio between iterations rises above
#'        `sratmax` (i.e., stress stops improving). Defaults: `1e-4`, `1e-7`
#'        and `0.999999`, as in `vegan::monoMDS()`.
#'
#' @param seed   An integer seed for the random starting configurations.
#'        The same seed gives the same result. Default: `0`
#'
#' @return A numeric matrix with one row per sample and one column per
#'   dimension, named `"NMDS1"`, `"NMDS2"`, etc. Its attributes are:
#'
#'   * `stress`: The stress of the solution (formula `stress`), from 0
#'     (perfect) to 1.
#'   * `converged`: `TRUE` if a random start found the best solution again.
#'   * `tries`: The number of random starts run.
#'   * `best_try`: Which start gave the solution: `0` for the PCoA start,
#'     else the random start's number.
#'
#' @details
#'
#' **Method**
#'
#' NMDS minimizes Kruskal's stress: how far the distances between points
#' \eqn{d_{ij}} are from the closest sequence \eqn{\hat{d}_{ij}} that has the
#' same rank order as the distances in `dm` (found by monotone regression).
#' Formula 1 is
#' \eqn{\sqrt{\sum(d_{ij} - \hat{d}_{ij})^2 / \sum d_{ij}^2}}; formula 2
#' divides by \eqn{\sum(d_{ij} - \bar{d})^2} instead.
#'
#' Each start minimizes stress by steepest descent, with Kruskal's step-size
#' rules, as in Peter Minchin's `monoMDS` (the default engine of
#' `vegan::metaMDS()`). This is a line-by-line port of it, so from the same
#' starting configuration, a start gives exactly the same stress and
#' coordinates as `vegan::monoMDS()`.
#'
#' Steepest descent only finds a local minimum, so, as in `metaMDS()`:
#'
#' 1. The first start is the PCoA solution, from [pcoa()].
#' 2. Then come random starts, with coordinates uniform on \[0, 1). A start
#'    whose stress is within 0.0005 of the best so far is compared with it by
#'    a Procrustes rotation; the best is "found again" when the root mean
#'    squared difference is below 0.005 and no point differs by more than
#'    0.01. At least `try` random starts are run, and more, up to `trymax`,
#'    until the best is found again.
#' 3. The best solution is centered, rotated to its principal components, and
#'    scaled so that its largest distance equals the largest distance in
#'    `dm`. When `dm` has a `maxdist` attribute (as from `vegan::vegdist()`),
#'    it is instead scaled to half-change units, as in `vegan::postMDS()`.
#'
#'
#' **Speed and memory**
#'
#' `metaMDS()` and `monoMDS()` copy `dm` into an \eqn{n \times n} matrix, and
#' keep several arrays with an element for every pair of samples, for every
#' start. `nmds()` sorts the pairs once, into one shared array of 4 bytes per
#' pair, and each start that runs at the same time needs about 16 more bytes
#' per pair. `cpus` starts run at once, each on its own thread.
#'
#'
#' **Reproducibility**
#'
#' The random starting configurations come from `seed`, not from R's random
#' number generator, so `set.seed()` has no effect, and R's random number
#' stream is left untouched. They differ from the ones `metaMDS()` draws, so
#' the two can settle on different, but equally good, local minima.
#' Results are identical from run to run and for any number of `cpus`.
#'
#' Axes are flipped so that each one's coordinate with the largest absolute
#' value is positive.
#'
#' `k` larger than \eqn{n - 1} is reduced to \eqn{n - 1}. A stress below
#' 0.001 gives a warning, as in `metaMDS()`: there may be too few samples for
#' `k` dimensions.
#'
#' @references
#' Kruskal, J. B. (1964). Nonmetric multidimensional scaling: A numerical
#' method. *Psychometrika*, 29(2), 115-129. \doi{10.1007/BF02289694}
#'
#' Minchin, P. R. (1987). An evaluation of the relative robustness of
#' techniques for ecological ordination. *Vegetatio*, 69, 89-107.
#' \doi{10.1007/BF00038690}
#'
#' @seealso [pcoa()], [beta_div()]
#'
#' @export
#' @examples
#'     ord <- nmds(bray(ex_counts))
#'     ord
#'
#'     attr(ord, 'stress')
#'
#'     plot(ord, asp = 1)
#'
nmds <- function (
    dm, k = 2L, try = 20L, trymax = 20L, maxit = 200L, weakties = TRUE,
    stress = 1L, smin = 1e-4, sfgrmin = 1e-7, sratmax = 0.999999,
    seed = 0L, cpus = n_cpus() ) {

  validate_cpus()
  validate_seed()

  k       <- nmds_arg(k,       '`k` must be a whole number of at least 1.',  whole = TRUE, lo = 1)
  try     <- nmds_arg(try,     '`try` must be a whole number of at least 0.', whole = TRUE)
  trymax  <- nmds_arg(trymax,  '`trymax` must be a whole number of at least 0.', whole = TRUE)
  maxit   <- nmds_arg(maxit,   '`maxit` must be a whole number of at least 0.', whole = TRUE)
  smin    <- nmds_arg(smin,    '`smin` must be a number of at least 0.')
  sfgrmin <- nmds_arg(sfgrmin, '`sfgrmin` must be a number of at least 0.')
  sratmax <- nmds_arg(sratmax, '`sratmax` must be a number from 0 to 1.', hi = 1)

  if (!isTRUE(weakties) && !isFALSE(weakties))
    stop('`weakties` must be TRUE or FALSE.', call. = FALSE)

  if (!is.numeric(stress) || length(stress) != 1 || !stress %in% 1:2)
    stop('`stress` must be 1 or 2.', call. = FALSE)

  try <- min(try, trymax)


  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  k <- min(k, n - 1L)

  if (k < 1) {
    # One sample or none: no dimensions.
    return(structure(
      .Data     = matrix(numeric(0), nrow = n, ncol = 0, dimnames = list(labels, NULL)),
      stress    = 0,
      converged = FALSE,
      tries     = 0L,
      best_try  = 0L ))
  }

  if (n > 65536)
    stop('nmds() supports at most 65,536 samples.', call. = FALSE) # nocov

  rng <- range(dm)
  if (!all(is.finite(rng)))
    stop('`dm` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (rng[[2]] <= 0)
    stop('`dm` must have at least one distance greater than 0.', call. = FALSE)
  if (rng[[1]] < -sqrt(.Machine$double.eps))
    warning('`dm` has negative distances.', call. = FALSE)


  # Pairs sorted by distance, shared by all starts.
  prep <- .Call(C_nmds_prep, dm, n)

  opts <- function (first, m)
    as.integer(c(maxit, weakties, stress, cpus, seed, first, m))
  tols <- c(smin, sfgrmin, sratmax)

  points <- function (res, i)
    matrix(res$points[seq_len(n * k) + (i - 1) * n * k], n, k)


  # Run 0, from the PCoA solution. Axes with eigenvalues <= 0 are all 0.
  init <- .Call(C_pcoa, dm, n, k, cpus)$points
  res  <- .Call(C_nmds_run, prep, n, k, init, opts(0L, 1L), tols)

  best      <- list(points = points(res, 1), stress = res$stress)
  best_try  <- 0L
  tries     <- 0L
  converged <- 0L


  # Random starts, `cpus` at a time. Each batch is evaluated in order, so the
  # outcome is what running them one by one would give.
  more <- function () tries < try || (tries < trymax && converged == 0)

  while (more()) {

    first <- tries + 1L
    m     <- if (tries < try) min(cpus, try - tries) else min(cpus, trymax - tries)
    res   <- .Call(C_nmds_run, prep, n, k, NULL, opts(first, m), tols)

    for (i in seq_len(m)) {

      if (!more()) break
      tries <- tries + 1L

      if (best$stress - res$stress[[i]] > -0.0005) {

        x   <- points(res, i)
        pro <- procrustes_fit(best$points, x)

        if (res$stress[[i]] < best$stress) {
          best      <- list(points = x, stress = res$stress[[i]])
          best_try  <- tries
          converged <- 0L
        }

        if (pro$rmse < 0.005 && max(pro$resid) < 0.01)
          converged <- converged + 1L
      }
    }
  }

  if (best$stress < 0.001)
    warning('stress is (nearly) zero: there may be too few samples', call. = FALSE)


  #--------------------------------------------------
  # Center, rotate to principal components, and scale.
  #--------------------------------------------------
  x <- scale(best$points, scale = FALSE)
  x <- unname(stats::prcomp(x, center = FALSE)$x)

  maxdist    <- dv$maxdist
  halfchange <- is.numeric(maxdist) && length(maxdist) == 1 && !is.na(maxdist)
  threshold  <- if (halfchange) maxdist * 0.8 else -Inf

  if (halfchange && threshold < sqrt(.Machine$double.eps))
    warning('dubious maxdist ', maxdist, ': similarities instead of dissimilarities?', call. = FALSE)

  sc <- .Call(C_nmds_scale, dm, x, threshold, cpus)

  if (halfchange && sc[[2]] < 10) {
    warning('skipping half-change scaling: too few points below threshold', call. = FALSE)
    halfchange <- FALSE
  }

  hc <- if (halfchange) (maxdist - sc[[3]]) / 2 / sc[[4]] else sc[[1]] / rng[[2]]
  x  <- x / hc

  # Each axis's largest absolute coordinate is positive.
  for (j in seq_len(k)) {
    v <- x[, j]
    if (v[[which.max(abs(v))]] < 0) x[, j] <- -v
  }

  dimnames(x) <- list(labels, paste0('NMDS', seq_len(k)))

  structure(
    .Data     = x,
    stress    = best$stress,
    converged = converged > 0,
    tries     = tries,
    best_try  = best_try )
}



# Validate a numeric scalar argument; return it as an integer if `whole`.
nmds_arg <- function (x, msg, whole = FALSE, lo = 0, hi = Inf) {

  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x >= lo && x <= hi &&
    (!whole || x %% 1 == 0)

  if (!ok) stop(msg, call. = FALSE)

  if (whole) as.integer(min(x, .Machine$integer.max)) else as.double(x)
}



# Symmetric Procrustes comparison of configuration `y` with `x`, as
# `summary(vegan::procrustes(x, y, symmetric = TRUE))`: the root mean
# squared error and each point's residual.
procrustes_fit <- function (x, y) {

  x <- scale(x, scale = FALSE)
  y <- scale(y, scale = FALSE)
  x <- x / sqrt(sum(x^2))
  y <- y / sqrt(sum(y^2))

  sol  <- svd(crossprod(x, y))
  rot  <- sol$v %*% t(sol$u)
  c    <- sum(sol$d) / sum(y^2)
  yrot <- c * y %*% rot
  ss   <- sum(x^2) + c * c * sum(y^2) - 2 * c * sum(sol$d)

  list(
    rmse  = sqrt(max(ss, 0) / nrow(x)),
    resid = sqrt(rowSums((yrot - x)^2)) )
}
