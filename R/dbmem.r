# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# dbmem() is written from Borcard & Legendre (2002) and Dray, Legendre &
# Peres-Neto (2006); adespatial's dbmem() (GPL) was run only to check its
# outputs. The computation is in src/dbmem.c.



#' Distance-Based Moran's Eigenvector Maps
#'
#' Turns the samples' positions into a set of uncorrelated spatial
#' variables, from broad-scale patterns to fine-scale ones: distance-based
#' Moran's eigenvector maps (dbMEM, formerly PCNM, principal coordinates of
#' neighbour matrices). Used as explanatory variables, in [rda()],
#' [dbrda()], [varpart()] or [forward_select()], they model the spatial
#' structure of a community. Gives the same MEMs as `adespatial::dbmem()`,
#' computed in parallel C.
#'
#' @inheritParams documentation
#'
#' @param space  The samples' positions, as for [moran_i()] but without
#'        spatial weights:
#'   * **Coordinates:** a data frame, or a matrix that isn't square, with
#'     one row per sample and one numeric column per coordinate; or, with
#'     `data`, the names of its coordinate columns. Distances between them
#'     are Euclidean. Project longitudes and latitudes first, or give
#'     great-circle distances as a `dist`.
#'   * **Geographic distances:** a `dist` object.
#'
#'   The MEMs' row names are `space`'s row names (or a `dist`'s labels).
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, whose coordinate columns `space` names. Default: `NULL`.
#'
#' @param threshold  The truncation distance \eqn{t}: samples at most
#'        \eqn{t} apart are neighbours. `NULL` (the default) uses the
#'        longest edge of the minimum spanning tree of the distances, the
#'        smallest distance that leaves no sample, or group of samples,
#'        without a neighbour outside it. A smaller value is allowed, with
#'        a warning, but leaves groups of samples unconnected to the rest.
#'
#' @param autocor  Which MEMs to return, by the sign of their spatial
#'        autocorrelation (Moran's I), as adespatial's `MEM.autocor`:
#'   * `'positive'`: positive autocorrelation (the default), the MEMs
#'     usually used to model spatial structure;
#'   * `'negative'`: negative autocorrelation;
#'   * `'non-null'`: either;
#'   * `'all'`: all \eqn{n - 1}, including any whose eigenvalue is 0.
#'
#' @return A numeric matrix of class `ecodive_dbmem`, with one row per
#'   sample and one column per MEM, named `MEM1`, `MEM2`, ..., in order of
#'   decreasing eigenvalue (from broad-scale positive autocorrelation to
#'   fine-scale negative). Each column has mean 0 and sum of squares
#'   \eqn{n}, and the columns are orthogonal. Its attributes are:
#'
#'   * `values`: The eigenvalues, as adespatial gives them: those of the
#'     doubly centred spatial weights \eqn{W}, divided by \eqn{n}.
#'   * `moran_i`: Each MEM's Moran's I, with the weights \eqn{W} as
#'     they are (`moran_i(x, style = 'none')`): \eqn{n^2 \lambda / S_0},
#'     for eigenvalue \eqn{\lambda} (`values`) and \eqn{S_0} the sum of the
#'     weights.
#'   * `threshold`: The truncation distance \eqn{t}.
#'   * `method`: `'dbmem'`.
#'
#'   The matrix can be given as is wherever a table of explanatory variables
#'   is taken (or converted with `as.data.frame()`). Columns taken with `[`
#'   are a plain matrix, without the attributes.
#'
#' @details
#'
#' **Method**
#'
#' dbMEM (Borcard & Legendre 2002; Dray et al. 2006) describes the spatial
#' arrangement of the samples at every scale the sampling can resolve:
#'
#' 1. Samples at most \eqn{t} apart are neighbours. Each pair of neighbours
#'    has spatial weight \eqn{w_{ij} = 1 - (d_{ij} / 4t)^2}, and every other
#'    pair 0.
#' 2. The weights matrix \eqn{W} is doubly centred (its rows and columns
#'    made to sum to 0), and its eigenvectors are the MEMs: orthogonal
#'    variables, each with Moran's I proportional to its eigenvalue. Those
#'    with large positive eigenvalues vary smoothly across the study area
#'    (neighbours have similar values); those with large negative ones
#'    alternate between neighbours.
#'
#' This is principal coordinates analysis ([pcoa()]) of the truncated
#' distances (distances larger than \eqn{t} replaced by \eqn{4t}), the
#' original PCNM method, whose eigenvectors are the same: the PCoA
#' eigenvalues are \eqn{8 t^2 (n \lambda + 1)} for each eigenvalue
#' \eqn{\lambda} in `values`. Keeping the MEMs with positive eigenvalues
#' keeps those with positive Moran's I, which PCNM's positive PCoA
#' eigenvalues did not (Dray et al. 2006).
#'
#' Eigenvalues within \eqn{\sqrt{\epsilon}} (about \eqn{1.5 \times 10^{-8}})
#' of the largest one in absolute value are taken to be 0: such MEMs are
#' neither positive nor negative, and are returned only with `autocor =
#' 'all'`, with `values` of exactly 0. Pairs whose distance exceeds \eqn{t}
#' by less than one part in a million (rounding in computed distances) are
#' neighbours too.
#'
#' The MEMs are found only up to sign: each column is made positive where
#' its largest absolute value is. MEMs whose eigenvalues are equal (as on a
#' regular grid) are any orthogonal basis of their shared space, as in any
#' eigendecomposition.
#'
#' Samples at the same place (distance 0) are allowed. A sample with no
#' neighbours (possible only with a `threshold` below the default) gives a
#' warning; it still has values in every MEM.
#'
#'
#' **Speed, memory, and reproducibility**
#'
#' The eigenvectors of an \eqn{n \times n} matrix take \eqn{O(n^3)} time and
#' \eqn{O(n^2)} memory, about \eqn{16 n^2} bytes (64 MB at \eqn{n = 2000}),
#' plus 16 bytes for each pair of neighbours.
#' The work is shared among `cpus` threads, by a Householder reduction to
#' tridiagonal form and the QL algorithm, in C, without the BLAS. Results
#' are bit-identical for any `cpus`.
#'
#'
#' **Differences from adespatial**
#'
#' `adespatial::dbmem()` returns a data frame (of class `orthobasisSp`)
#' with the same columns, up to sign, and the same `values`. It stores the
#' spatial weights (`listw`) instead of each MEM's Moran's I, which is here
#' attribute `moran_i`. Its `MEM.autocor` is `autocor`, `thresh` is
#' `threshold`, and `xyORdist` is `space`; it has no `data`.
#'
#'
#' **Same names in other packages**
#'
#' `dbmem()` masks `adespatial::dbmem()` when both packages are attached:
#' call either as `ecodive::dbmem()` or `adespatial::dbmem()`.
#'
#' @references
#' Borcard, D., & Legendre, P. (2002). All-scale spatial analysis of
#' ecological data by means of principal coordinates of neighbour matrices.
#' *Ecological Modelling*, 153(1-2), 51-68.
#' \doi{10.1016/S0304-3800(01)00501-4}
#'
#' Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbour
#' matrices (PCNM). *Ecological Modelling*, 196(3-4), 483-493.
#' \doi{10.1016/j.ecolmodel.2006.02.015}
#'
#' @seealso [moran_i()], which tests a variable for spatial
#'   autocorrelation; [forward_select()] and [varpart()], which use MEMs as
#'   explanatory variables; [plot.ecodive_dbmem()].
#'
#' @export
#' @examples
#'     # 40 sites, with taxa that follow the east-west gradient and pH.
#'     set.seed(1)
#'     sites <- data.frame(
#'       row.names = paste0('S', 1:40),
#'       east  = runif(40, 0, 10),
#'       north = runif(40, 0, 10),
#'       ph    = runif(40, 5, 8) )
#'     counts <- matrix(rpois(40 * 12, 3), 40, 12)
#'     counts[, 1:4] <- counts[, 1:4] + round(2 * sites$east)
#'     counts[, 5:8] <- counts[, 5:8] + round(3 * (sites$ph - 5))
#'
#'     mem <- dbmem(c('east', 'north'), sites)
#'     mem
#'     plot(mem)
#'     plot(mem, sites[, c('east', 'north')])
#'
#'     # The MEMs that explain the community (a formula of them), then the
#'     # variation explained by pH and by space.
#'     fs <- forward_select(counts, data = mem, norm = 'hellinger')
#'     fs
#'     varpart(
#'       counts, ph = ~ ph, space = attr(fs, 'selected'),
#'       data = cbind(sites, mem), norm = 'hellinger' )
#'
dbmem <- function (
    space, data = NULL, threshold = NULL, autocor = 'positive', cpus = n_cpus() ) {

  validate_cpus()

  if (!is.character(autocor) || length(autocor) != 1 || is.na(autocor) ||
      is.na(code <- pmatch(autocor, c('positive', 'negative', 'non-null', 'all'))))
    stop('`autocor` must be "positive", "negative", "non-null", or "all".', call. = FALSE)

  if (missing(space) || is.null(space))
    stop('`space` is required: coordinates, or distances.', call. = FALSE)

  if (!is.null(threshold) && (
      !is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
      !is.finite(threshold) || !(threshold > 0) ))
    stop('`threshold` must be a positive number, or NULL.', call. = FALSE)

  sp <- moran_space(space, data, distances = TRUE)
  n  <- sp$n

  if (n < 3)
    stop('`space` must have at least 3 samples.', call. = FALSE)

  res <- .Call(
    C_dbmem, sp$w, n, if (is.null(threshold)) NA_real_ else as.double(threshold),
    code - 1L, cpus )

  if (res$n_pairs == 0) {
    if (res$mst == 0)
      stop('Every sample of `space` is at the same place.', call. = FALSE)
    stop(
      sprintf(
        'No two samples are within `threshold` (%s) of each other. ',
        format(res$threshold, digits = 4) ),
      'The default, the longest edge of the minimum spanning tree, is ',
      format(res$mst, digits = 4), '.', call. = FALSE )
  }

  # With the same allowance for rounding as the neighbours (src/dbmem.c).
  if (res$threshold * (1 + 1e-6) < res$mst)
    warning(
      sprintf(
        '`threshold` (%s) is below the longest edge of the minimum spanning tree (%s), ',
        format(res$threshold, digits = 4), format(res$mst, digits = 4) ),
      'so some samples aren\'t connected to the others through neighbours',
      if (res$islands > 0) sprintf(
        ', and %d %s none', res$islands,
        if (res$islands == 1) 'sample has' else 'samples have' ),
      '.', call. = FALSE )

  if (!res$converged) # nocov start
    warning('The eigendecomposition did not converge; MEMs may be inexact.', call. = FALSE)
  # nocov end

  k  <- ncol(res$vectors)
  nm <- if (k > 0) paste0('MEM', seq_len(k)) else character(0)

  if (k == 0)
    warning(
      'No MEM has ', switch(code, 'positive', 'negative', 'non-zero', ''),
      ' spatial autocorrelation.', call. = FALSE )

  out <- res$vectors
  dimnames(out) <- list(sp$labels, nm)

  attr(out, 'values')    <- stats::setNames(res$values, nm)
  attr(out, 'moran_i')   <- stats::setNames(res$moran, nm)
  attr(out, 'threshold') <- res$threshold
  attr(out, 'method')    <- 'dbmem'
  class(out) <- c('ecodive_dbmem', 'matrix', 'array')
  out
}



#------------------------------------------------------------------------------
# Printing and plotting.
#------------------------------------------------------------------------------

#' @rdname print_ecodive
#' @export
print.ecodive_dbmem <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  att <- attributes(x)
  sym <- print_sym()
  thr <- att$threshold

  head <- c(
    sprintf('A matrix: %d samples %s %d MEMs', nrow(x), sym$times, ncol(x)),
    "Moran's Eigenvector Maps",
    if (!is.null(thr)) sprintf('Threshold: %s', format(thr, digits = digits)) )

  m <- unclass(x)
  attributes(m) <- att[c('dim', 'dimnames')]
  df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
  if (is.null(colnames(m))) names(df) <- sprintf('[,%d]', seq_len(ncol(m)))

  cat(
    print_table_lines(
      df, head, n = n, width = width, digits = digits, types = FALSE,
      rows = 'samples', attrs = print_attrs(x, c('dim', 'dimnames', 'class')) ),
    sep = '\n' )
  invisible(x)
}


#' Plot Moran's eigenvector maps
#'
#' Each MEM's Moran's I, or maps of the MEMs over the samples' coordinates.
#'
#' @param x  The result of [dbmem()].
#'
#' @param coords  `NULL` (the default) for a bar chart of each MEM's Moran's
#'        I (attribute `moran_i`). Or the samples' coordinates, in the order
#'        of `x`'s rows: a matrix or data frame whose first two numeric
#'        columns are drawn as x and y, to map MEMs `mems`. Each sample is
#'        a square, filled for positive values and open for negative ones,
#'        sized by its absolute value.
#'
#' @param mems  Which MEMs to map with `coords`: column numbers or names of
#'        `x`. Default: the first four (or fewer).
#'
#' @param ...  Passed on to [graphics::barplot()] (for the bar chart) or
#'        [graphics::plot()] (for each map), such as `main`, `col` or `cex`.
#'
#' @return `x`, invisibly.
#'
#' @details
#' The bar chart reads only attribute `moran_i`, so columns taken with `[`
#' can't be drawn as one. Maps of several MEMs share the device, in a grid
#' of panels (restored afterwards).
#'
#' @seealso [dbmem()]
#'
#' @export
#' @examples
#'     g   <- expand.grid(x = 1:10, y = 1:6)
#'     mem <- dbmem(g)
#'     plot(mem)
#'     plot(mem, g, mems = 1:6)
plot.ecodive_dbmem <- function (x, coords = NULL, mems = NULL, ...) {

  args <- list(...)
  main <- "Moran's Eigenvector Maps"
  k    <- ncol(x)
  nms  <- colnames(x)
  if (is.null(nms) && k > 0) nms <- paste0('MEM', seq_len(k))

  if (is.null(coords)) {

    if (!is.null(mems))
      stop('`mems` needs `coords`: give the samples\' coordinates to map MEMs.', call. = FALSE)

    mi <- attr(x, 'moran_i', exact = TRUE)
    if (is.null(mi))
      stop(
        '`x` has lost its Moran\'s I (parts taken with `[` drop it): ',
        'plot the whole result.', call. = FALSE )
    if (k == 0)
      stop('`x` has no MEMs to plot.', call. = FALSE)

    plot_do(
      graphics::barplot, c(list(unname(mi)), args),
      list(
        main = main, names.arg = nms, las = 2,
        ylab = "Moran's I", col = ifelse(mi >= 0, 'grey40', 'grey80'),
        border = NA ))
    graphics::abline(h = 0, col = PLOT_MUTED)
    return (invisible(x))
  }


  #--------------------------------------------------------
  # Maps over the coordinates.
  #--------------------------------------------------------
  if (k == 0)
    stop('`x` has no MEMs to plot.', call. = FALSE)

  if (!is.data.frame(coords) && !is.matrix(coords))
    stop('`coords` must be a matrix or data frame of the samples\' coordinates.', call. = FALSE)
  cd  <- as.data.frame(coords, stringsAsFactors = FALSE)
  num <- which(vapply(cd, is.numeric, logical(1)))
  if (length(num) < 2 || nrow(cd) != nrow(x))
    stop(
      '`coords` must have one row per sample of `x` (', nrow(x), ') and ',
      'at least two numeric columns.', call. = FALSE )
  xy <- as.matrix(cd[, num[1:2]])
  if (!all(is.finite(xy)))
    stop('`coords` must not contain NA, NaN, or infinite values.', call. = FALSE)
  lab <- colnames(xy)

  if (is.null(mems)) mems <- seq_len(min(4L, k))
  if (is.character(mems)) {
    idx <- match(mems, nms)
    if (anyNA(idx))
      stop('`mems` names MEMs that are not in `x`: ', name_list(mems[is.na(idx)]), call. = FALSE)
    mems <- idx
  }
  if (!is.numeric(mems) || !length(mems) || anyNA(mems) || any(mems %% 1 != 0) ||
      any(mems < 1) || any(mems > k))
    stop('`mems` must be column numbers or names of `x`, from 1 to ', k, '.', call. = FALSE)
  mems <- as.integer(mems)

  if (length(mems) > 1) {
    nr  <- floor(sqrt(length(mems)))
    nc  <- ceiling(length(mems) / nr)
    old <- graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.5, 1), mgp = c(1.5, 0.5, 0))
    on.exit(graphics::par(old), add = TRUE)
  }

  m   <- unclass(x)
  col <- plot_arg(args, 'col', 'black')
  cex <- plot_arg(args, 'cex', 2)

  for (j in mems) {
    v  <- m[, j]
    sz <- cex * sqrt(abs(v) / max(abs(v), .Machine$double.xmin))
    plot_do(
      graphics::plot, c(list(xy[, 1], xy[, 2]), plot_drop(args, c('col', 'cex', 'pch'))),
      list(
        main = if (length(mems) == 1) paste0(main, ': ', nms[[j]]) else nms[[j]],
        type = 'n', asp = 1, las = 1, xlab = lab[[1]], ylab = lab[[2]] ))
    graphics::points(
      xy[, 1], xy[, 2], pch = ifelse(v >= 0, 15, 0), cex = sz, col = col )
  }

  invisible(x)
}
