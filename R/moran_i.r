# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# moran_i() follows esda 2.10.0's Moran (esda/moran.py): its statistic,
# expected value, variances under normality and randomisation, and analytic
# p-values; and libpysal 4.15.1's weight transformations and the sums s0, s1
# and s2 (libpysal/weights/weights.py), which src/moran_i.c computes. Those
# sources are under this license:
#
#   esda: Copyright 2018 PySAL Developers.
#   libpysal: Copyright (c) 2007-, PySAL Developers.
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   3. Neither the name of the copyright holder nor the names of its
#      contributors may be used to endorse or promote products derived from
#      this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
#   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' Moran's I Test of Spatial Autocorrelation
#'
#' Tests whether a numeric variable, such as a taxon's abundance or a
#' sample's alpha diversity, is more similar (or less similar) at nearby
#' sites than at distant ones. Gives the same statistic and analytic tests
#' as PySAL's `esda.Moran`, with a permutation test computed in parallel C.
#'
#' @inheritParams documentation
#'
#' @param x  The variable: a numeric vector with one value per sample,
#'        matched to the samples of `space` by name when both have names,
#'        otherwise in the same order. Or, with `data`, the name of one of
#'        its numeric columns, or a vector with one value per row of `data`.
#'
#' @param space  The samples' positions, in one of three forms:
#'   * **Spatial weights:** a square numeric matrix, samples \eqn{\times}
#'     samples, where `space[i, j]` is the weight of sample `j` as a
#'     neighbour of sample `i`: larger for closer samples, and 0 for samples
#'     that aren't neighbours. Weights must be finite and not negative, and
#'     needn't be symmetric. The diagonal is ignored, so that, for instance,
#'     `1 / as.matrix(d)` can be given as is.
#'   * **Geographic distances:** a `dist` object. The weights are the inverse
#'     distances, \eqn{w_{ij} = 1 / d_{ij}}.
#'   * **Coordinates:** a data frame, or a matrix that isn't square, with one
#'     row per sample and one numeric column per coordinate; or, with
#'     `data`, the names of its coordinate columns. The weights are the
#'     inverse Euclidean distances. Project longitudes and latitudes first,
#'     or give great-circle distances as a `dist`.
#'
#'   Samples are matched to `x` and `data` by the names in `space` (its row
#'   names, or a `dist`'s labels) when it has them.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, whose columns `x` and `space` can name. Rows are matched to
#'        the samples of `space` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order.
#'        Default: `NULL`.
#'
#' @param style  How the weights are transformed before use, as esda's
#'        `transformation` (and spdep's `style`):
#'   * `'row'`: row-standardised, so that each sample's weights sum to 1
#'     (esda's and spdep's default);
#'   * `'binary'`: every non-zero weight becomes 1;
#'   * `'variance'`: variance-stabilising (Tiefelsdorf et al. 1999): each
#'     sample's weights are divided by the square root of their sum of
#'     squares, then all are scaled to sum to the number of samples;
#'   * `'none'`: as given.
#'
#'   Default: `'row'`.
#'
#' @param alternative  Which values of I count as more extreme than the
#'        observed one, for both the permutation and the analytic p-values:
#'        `'greater'` (more positive: neighbours are more alike than by
#'        chance), `'less'` (more negative: neighbours are less alike), or
#'        `'two.sided'`. Default: `'greater'`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones: each
#'        permutation `p` gives sample `i` the value of sample `p[i]`.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A list of class `c('ecodive_moran_i', 'ecodive_test')`, with:
#'
#'   * `statistic`: Moran's I, `c(I = ...)`.
#'   * `p_value`: The permutation p-value.
#'   * `z`: The standardized effect size, `(I - mean(i_perm)) / sd(i_perm)`.
#'   * `n_perm`: The number of permutations.
#'   * `expected`: The expected value of I under no spatial
#'     autocorrelation, \eqn{-1 / (n - 1)}.
#'   * `sd_norm`, `p_value_norm`: The standard deviation of I, and the
#'     analytic p-value of \eqn{(I - E[I]) / sd}, under the normality
#'     assumption (esda's `seI_norm` and `p_norm`).
#'   * `sd_rand`, `p_value_rand`: The same under the randomisation
#'     assumption (esda's `seI_rand` and `p_rand`; spdep's default; and
#'     `ape::Moran.I()`'s `sd` and `p.value`).
#'   * `style`: How the weights were transformed.
#'   * `alternative`: `'greater'`, `'less'`, or `'two.sided'`.
#'   * `n`: The number of samples.
#'   * `space_type`: What `space` was, which says how the weights were made:
#'     `'weights'`, `'distances'` or `'coordinates'`.
#'   * `permutation`: `'free'` or `'supplied matrix'`.
#'
#'   Attribute `i_perm` is a vector of the statistics of the permutations,
#'   as esda's `sim`. `as.data.frame()` gives the main results as one row.
#'
#' @details
#'
#' **Method**
#'
#' Moran's I (Moran 1950; Cliff & Ord 1981) of values \eqn{x} with spatial
#' weights \eqn{w_{ij}} is
#' \deqn{I = \frac{n}{S_0} \frac{\sum_i \sum_j w_{ij} z_i z_j}{\sum_i z_i^2}}
#' where \eqn{z_i = x_i - \bar{x}} and \eqn{S_0 = \sum_i \sum_j w_{ij}}.
#' It is near its expected value \eqn{-1 / (n - 1)} when values are
#' arranged at random, positive when neighbours have similar values, and
#' negative when they have dissimilar ones. With row-standardised weights,
#' it is the slope of the regression of each sample's neighbours' mean
#' value (the spatial lag) on its own.
#'
#' The analytic tests compare \eqn{(I - E[I]) / sd} with the standard normal
#' distribution, with the variance of I under one of two assumptions (Cliff
#' & Ord 1981): normality, where the values are a sample from a normal
#' distribution; or randomisation, where they are fixed and arranged at
#' random, and the variance depends on their kurtosis. The variance under
#' randomisation needs at least 4 samples. When the weights make I the same
#' for every arrangement of the values (such as equal weights between every
#' pair of samples), both standard deviations are 0, and their p-values are
#' `NA`.
#'
#' The permutation test shuffles the values among the samples, holding the
#' weights fixed, and recomputes I. The p-value is the proportion of
#' permutations, counting the observed one, with an I at least as extreme:
#' \eqn{(1 + \#\{I_{perm} \ge I\}) / (1 + permutations)} for
#' `alternative = 'greater'`. Permuted statistics within
#' `sqrt(.Machine$double.eps)` of the observed one count as ties.
#'
#' A sample whose weights are all 0 (an "island", with no neighbours) is
#' kept, with a warning: it counts in \eqn{n} and in \eqn{\bar{x}}, but
#' contributes nothing to the sums, as in esda. When `x` is the same for
#' every sample, I is undefined: the statistic and p-values are `NA`, with a
#' warning.
#'
#'
#' **Speed, memory, and reproducibility**
#'
#' Only the packed lower triangle of the weights' symmetric part,
#' \eqn{(W + W^T) / 2}, is kept, which is all that I needs, and each
#' permutation is one pass over it: \eqn{O(n^2)} time per permutation.
#' Permutations run in parallel. They come from `seed`, not from R's random
#' number generator, so `set.seed()` has no effect, and p-values are
#' identical from run to run and for any number of `cpus`. They agree with
#' esda's (or spdep's) exactly when the same permutations are used.
#'
#'
#' **Differences from esda.Moran**
#'
#' * `alternative` sets the tail of both the analytic and the permutation
#'   p-values. esda's analytic p-values are two-sided by default
#'   (`two_tailed`), and its permutation p-value is one-sided in the
#'   direction of the observed I (the smaller tail).
#' * `z` is \eqn{(I - mean) / sd} of the permuted values with the sample
#'   standard deviation, as for ecodive's other tests; esda's `z_sim` uses
#'   the population standard deviation.
#' * `space` may also be distances or coordinates, made into inverse-distance
#'   weights.
#' * Weights on the diagonal (self-weights) are ignored.
#'
#' `ape::Moran.I()` gives the same I as `style = 'row'` (it row-standardises
#' the weights), with the standard deviation under randomisation (`sd_rand`)
#' and a two-sided p-value by default. `spdep::moran.test()` gives the same
#' I as `style = 'row'` for `mat2listw()`'s default style `'W'`, with the
#' variance under randomisation by default (`randomisation = TRUE`); its
#' `adjust.n` (the default) leaves islands out of \eqn{n}, which `moran_i()`
#' doesn't.
#'
#' @references
#' Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models and
#' Applications*. Pion.
#'
#' Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
#' *Biometrika*, 37(1/2), 17-23. \doi{10.2307/2332142}
#'
#' Tiefelsdorf, M., Griffith, D. A., & Boots, B. (1999). A
#' variance-stabilizing coding scheme for spatial link matrices.
#' *Environment and Planning A*, 31(1), 165-180. \doi{10.1068/a310165}
#'
#' @seealso [mantel()], which tests distances against distances.
#'
#' @export
#' @examples
#'     # 40 sites, with a taxon more abundant towards the east.
#'     set.seed(1)
#'     sites <- data.frame(
#'       row.names = paste0('S', 1:40),
#'       east  = runif(40, 0, 10),
#'       north = runif(40, 0, 10) )
#'     sites$abund <- rpois(40, 2 + sites$east)
#'
#'     # Inverse-distance weights from the coordinates.
#'     res <- moran_i('abund', c('east', 'north'), sites)
#'     res
#'     plot(res)
#'
#'     # The same from a dist object of the distances.
#'     geo <- dist(sites[, c('east', 'north')])
#'     moran_i(sites$abund, geo)
#'
#'     # Binary weights to each site's 5 nearest neighbours.
#'     d <- as.matrix(geo)
#'     w <- t(apply(d, 1, function (r) rank(r, ties.method = 'first') %in% 2:6))
#'     dimnames(w) <- dimnames(d)
#'     moran_i(sites$abund, w * 1)
#'
moran_i <- function (
    x, space, data = NULL, style = 'row', alternative = 'greater',
    permutations = 999L, seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  style       <- match.arg(style, c('row', 'binary', 'variance', 'none'))
  alternative <- match.arg(alternative, c('greater', 'less', 'two.sided'))

  if (missing(space) || is.null(space))
    stop('`space` is required: spatial weights, distances, or coordinates.', call. = FALSE)
  if (missing(x) || is.null(x))
    stop('`x` is required: one value per sample, or a column name of `data`.', call. = FALSE)

  sp <- moran_space(space, data)
  n  <- sp$n

  if (n < 3)
    stop('`space` must have at least 3 samples.', call. = FALSE)

  x  <- moran_values(x, data, sp$labels, n)
  pm <- perm_matrix(permutations, n)
  wt <- moran_weights(sp, style)

  i_all <- moran_stat(x, wt, pm, seed, cpus)
  mom   <- moran_moments(x, wt, n)

  i_obs  <- i_all[[1]]
  i_perm <- i_all[-1]

  p_norm <- function (sd) {
    if (is.na(sd) || is.na(i_obs)) return (NA_real_)
    z <- (i_obs - mom$expected) / sd
    switch(alternative,
      greater   = stats::pnorm(z, lower.tail = FALSE),
      less      = stats::pnorm(z),
      two.sided = 2 * stats::pnorm(-abs(z)) )
  }

  new_ecodive_test(
    class        = 'ecodive_moran_i',
    title        = "Moran's I",
    statistic    = c(I = i_obs),
    p_value      = perm_p(i_obs, i_perm, alternative),
    z            = perm_z(i_obs, i_perm),
    n_perm       = pm$n_perm,
    expected     = mom$expected,
    sd_norm      = mom$sd_norm,
    p_value_norm = p_norm(mom$sd_norm),
    sd_rand      = mom$sd_rand,
    p_value_rand = p_norm(mom$sd_rand),
    style        = style,
    alternative  = alternative,
    n            = n,
    space_type   = sp$type,
    perm         = i_perm,
    permutation  = if (is.null(pm$perms)) 'free' else 'supplied matrix' )
}



#------------------------------------------------------------------------------
# Lines added to print.ecodive_test()'s printout.
#------------------------------------------------------------------------------
#' @export
test_details.ecodive_moran_i <- function (x, digits) {

  fmt <- function (v) if (is.na(v)) 'NA' else formatC(v, digits = digits, format = 'g')
  pv  <- function (v) if (is.na(v)) 'NA' else format.pval(v, digits = digits)

  c(
    paste0('Expected I: ', fmt(x$expected)),
    paste0('Under normality: sd = ', fmt(x$sd_norm), ', p = ', pv(x$p_value_norm)),
    paste0('Under randomisation: sd = ', fmt(x$sd_rand), ', p = ', pv(x$p_value_rand)),
    paste0('Weights: ', switch(x$style,
      row      = 'row-standardised',
      binary   = 'binary',
      variance = 'variance-stabilising',
      none     = 'as given' ),
      switch(x$space_type,
        weights     = '',
        distances   = ' inverse distances',
        coordinates = ' inverse distances between coordinates' )),
    paste0('Alternative: ', switch(x$alternative,
      greater   = 'greater (neighbours are more alike than by chance)',
      less      = 'less (neighbours are less alike than by chance)',
      two.sided = 'two-sided' )),
    paste0('Samples: ', x$n) )
}



#------------------------------------------------------------------------------
# `space` as list(type, n, labels, w, inverse): the raw weights `w`, either
# an n x n matrix (inverse = FALSE) or the packed distances whose inverses
# are the weights (inverse = TRUE).
#------------------------------------------------------------------------------
moran_space <- function (space, data) {

  coords <- NULL

  if (is.character(space) && is.null(dim(space))) {

    if (is.null(data))
      stop('`space` names columns, but `data` is not given.', call. = FALSE)
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    miss <- setdiff(space, names(data))
    if (length(miss) || anyNA(space))
      stop('`space` names columns that are not in `data`: ', name_list(miss), call. = FALSE)
    coords <- data[space]
  }

  else if (inherits(space, 'dist')) {

    v <- tryCatch(
      dist_vector(space),
      error = function (e)
        stop(gsub('`dm`', '`space`', conditionMessage(e), fixed = TRUE), call. = FALSE) )

    return (moran_inverse(v$dm, v$n, v$labels, 'distances'))
  }

  else if (is.data.frame(space) || (is.matrix(space) && nrow(space) != ncol(space))) {
    coords <- space
  }

  else if (is.matrix(space)) {

    if (!is.numeric(space) && !is.logical(space))
      stop('A spatial weights matrix (`space`) must be numeric.', call. = FALSE)

    rn <- rownames(space)
    cn <- colnames(space)
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn))
      stop('`space` has different row and column names.', call. = FALSE)

    w <- space
    storage.mode(w) <- 'double'
    attributes(w) <- list(dim = dim(w))
    off <- w[row(w) != col(w)]
    if (!all(is.finite(off)) || any(off < 0))
      stop('Spatial weights (`space`) must be finite and not negative.', call. = FALSE)

    return (list(
      type = 'weights', n = nrow(w), labels = if (is.null(rn)) cn else rn,
      w = w, inverse = FALSE ))
  }

  else {
    stop(
      '`space` must be a matrix of spatial weights, a dist object of ',
      'distances, or coordinates (a data frame, or the names of columns of ',
      '`data`).', call. = FALSE )
  }


  #--------------------------------------------------------
  # Coordinates: inverse Euclidean distances.
  #--------------------------------------------------------
  labels <- if (is.data.frame(coords)) {
    if (.row_names_info(coords) > 0) rownames(coords) else NULL
  } else rownames(coords)

  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  if (ncol(coords) == 0 || !all(vapply(coords, is.numeric, logical(1))))
    stop('Coordinates (`space`) must be numeric columns.', call. = FALSE)
  coords <- as.matrix(coords)
  storage.mode(coords) <- 'double'
  if (!all(is.finite(coords)))
    stop('Coordinates (`space`) must not contain NA, NaN, or infinite values.', call. = FALSE)

  d <- as.vector(stats::dist(coords))
  if (!all(is.finite(d)))
    stop('The distances between coordinates (`space`) are too large.', call. = FALSE)

  moran_inverse(d, nrow(coords), labels, 'coordinates')
}


# Packed distances, whose inverses are the weights: none may be 0.
moran_inverse <- function (d, n, labels, type) {

  zero <- which(d == 0)
  if (length(zero)) {
    j  <- rep.int(seq_len(n - 1L), (n - 1L):1)
    i  <- j + sequence((n - 1L):1)
    ab <- c(j[zero[[1]]], i[zero[[1]]])
    nm <- if (is.null(labels)) ab else labels[ab]
    stop(
      sprintf(
        '`space` puts samples %s and %s at the same place, where the inverse-distance weight is infinite. ',
        name_list(nm[[1]]), name_list(nm[[2]]) ),
      'Merge them, or give spatial weights as a matrix.', call. = FALSE )
  }

  list(type = type, n = as.integer(n), labels = labels, w = as.double(d), inverse = TRUE)
}



#------------------------------------------------------------------------------
# `x`, one finite value per sample of `space`, in its order.
#------------------------------------------------------------------------------
moran_values <- function (x, data, labels, n) {

  if (!is.atomic(x) || !is.null(dim(x)))
    stop(
      '`x` must be a numeric vector with one value per sample',
      if (!is.null(data)) ', or the name of a column of `data`', '.', call. = FALSE )

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, dn, nrow(data), ref = '`space`')
  }

  src <- is.character(x) && length(x) == 1 && !is.null(data) && x %in% names(data)
  v   <- per_sample(x, 'x', data, idx, labels, n, ref = '`space`')

  if (!is.numeric(v) || is.factor(v))
    stop(
      if (src) paste0('`data` column ', name_list(x), ' is not numeric.')
      else '`x` must be numeric.', call. = FALSE )

  if (!all(is.finite(v)))
    stop(
      if (src) paste0('`data` column ', name_list(x), ' must not contain NA, NaN, or infinite values.')
      else '`x` must not contain NA, NaN, or infinite values.', call. = FALSE )

  unname(as.double(v))
}



#------------------------------------------------------------------------------
# The transformed weights, as list(S, s0, s1, s2): S, the packed lower
# triangle of their symmetric part, and libpysal's sums. Warns of islands.
#------------------------------------------------------------------------------
moran_weights <- function (sp, style) {

  code <- match(style, c('row', 'binary', 'variance', 'none')) - 1L
  res  <- .Call(C_moran_weights, sp$w, sp$n, sp$inverse, code)

  if (!(res$s[[1]] > 0))
    stop('Every spatial weight (`space`) is 0.', call. = FALSE)

  if (res$islands > 0)
    warning(
      res$islands, ' sample', if (res$islands == 1) ' has' else 's have',
      ' no neighbours (all their weights are 0), so ',
      if (res$islands == 1) 'it adds' else 'they add',
      ' nothing to the sums of Moran\'s I.', call. = FALSE )

  list(S = res$S, s0 = res$s[[1]], s1 = res$s[[2]], s2 = res$s[[3]], n = sp$n)
}



#------------------------------------------------------------------------------
# Moran's I of `x` and of each permutation (from perm_matrix()), or NA (with
# a warning) when x is constant. Shared with dbmem(): with pm =
# perm_matrix(0L, n), the observed I alone.
#------------------------------------------------------------------------------
moran_stat <- function (x, wt, pm, seed = 0L, cpus = 1L) {

  n  <- length(x)
  z  <- x - mean(x)
  ss <- sum(z * z)

  if (all(x == x[[1]]) || !(ss > 0)) {
    warning(
      '`x` has the same value for every sample, so Moran\'s I is undefined.',
      call. = FALSE )
    return (rep(NA_real_, pm$n_perm + 1L))
  }

  q <- .Call(
    C_moran_perm, wt$S, z, pm$perms,
    as.integer(c(pm$n_perm, seed, cpus)) )

  n / wt$s0 * q / ss
}



#------------------------------------------------------------------------------
# The expected value of I, and its standard deviations under normality and
# randomisation, as esda's Moran.__moments(). A variance that is 0 to
# rounding error (I is the same for every arrangement) gives NA.
#------------------------------------------------------------------------------
moran_moments <- function (x, wt, n) {

  s0 <- wt$s0; s1 <- wt$s1; s2 <- wt$s2
  ei  <- -1 / (n - 1)
  n2  <- n * n
  s02 <- s0 * s0

  sd_of <- function (lead) {
    v <- lead - ei * ei
    if (!is.finite(v) || !(v > 1e-10 * abs(lead))) NA_real_ else sqrt(v)
  }

  sd_norm <- sd_of((n2 * s1 - n * s2 + 3 * s02) / ((n - 1) * (n + 1) * s02))

  z  <- x - mean(x)
  m2 <- sum(z^2) / n
  sd_rand <- if (n < 4 || !(m2 > 0)) NA_real_ else {
    k <- (sum(z^4) / n) / (m2 * m2)
    a <- n * ((n2 - 3 * n + 3) * s1 - n * s2 + 3 * s02)
    b <- k * ((n2 - n) * s1 - 2 * n * s2 + 6 * s02)
    sd_of((a - b) / ((n - 1) * (n - 2) * (n - 3) * s02))
  }

  list(expected = ei, sd_norm = sd_norm, sd_rand = sd_rand)
}
