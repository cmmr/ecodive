# Writes tests/tinytests/expected/dbmem.rds: inputs and reference results
# for test-dbmem.r, so that the tests need no reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/dbmem.r
#
# References, all run as black boxes (outputs only; their code is not read):
#   * adespatial 0.3.30's dbmem() (GPL-2 | GPL-3), with each MEM.autocor:
#     its MEMs (columns of the data frame it returns) and attribute
#     `values`. The MEMs of "positive", "negative" and "non-null" are checked
#     here to be the matching columns of those of "all" (up to sign, and up
#     to rotation within eigenvalues that are equal), so only those of "all"
#     are stored.
#   * spdep 1.3.1's moran(), with the spatial weights that adespatial
#     stores with its result (attribute `listw`): each MEM's Moran's I.
#   * vegan 2.6.4's spantree() (GPL-2): the longest edge of the minimum
#     spanning tree, adespatial's default threshold.
#   * vegan 2.6.4's pcnm() (GPL-2), with the same threshold, as a check of
#     the truncation: its `values` (the PCoA eigenvalues of the truncated
#     distances), `threshold`, and `vectors` (unit length; for "mite" only).
#     Not for mite_t2 and mite_dist_t08, where distances fall just above the
#     threshold by rounding (see one_case()).
#
# Each case holds space (coordinates, or a dist), threshold (the argument:
# NULL or a number), mst (vegan's longest MST edge), t (the threshold used),
# values (a list of adespatial's `values` for each MEM.autocor), vectors
# (adespatial's MEMs for "all", a plain matrix, or for "positive" in cases
# mite_t2, mite_dist_t08 and mite_dup), moran (their Moran's I),
# and pcnm (list(values, threshold, vectors), or NULL).

library(adespatial)
library(spdep)
library(vegan)

stopifnot(
  packageVersion('adespatial') == '0.3.30',
  packageVersion('spdep')      == '1.3.1',
  packageVersion('vegan')      == '2.6.4' )

data(mite.xy, package = 'vegan')

autocors <- c('positive', 'negative', 'non-null', 'all')

# The largest difference between the column spaces of a and b over runs of
# equal eigenvalues (`vals`), with b's columns those of a.
same_spaces <- function (a, b, vals) {
  cl <- cumsum(c(TRUE, abs(diff(vals)) > 1e-9 * max(abs(vals))))
  max(vapply(unique(cl), function (g) {
    j <- which(cl == g)
    max(abs(tcrossprod(a[, j, drop = FALSE]) - tcrossprod(b[, j, drop = FALSE])))
  }, 0))
}

one_case <- function (space, threshold = NULL, all_vectors = TRUE, with_pcnm = TRUE) {

  d   <- if (inherits(space, 'dist')) space else dist(space)
  mst <- max(spantree(d)$dist)
  t   <- if (is.null(threshold)) mst else threshold
  n   <- attr(d, 'Size')

  run <- function (a) {
    if (is.null(threshold)) dbmem(space, MEM.autocor = a, silent = TRUE)
    else                    dbmem(space, thresh = threshold, MEM.autocor = a, silent = TRUE)
  }

  all <- run('all')
  va  <- attr(all, 'values')
  A   <- unname(as.matrix(all))
  lw  <- attr(all, 'listw')
  mi  <- vapply(seq_len(ncol(A)), function (k) moran(A[, k], lw, n, Szero(lw))$I, 0)

  values <- list()
  tol    <- sqrt(.Machine$double.eps) * max(abs(va))
  for (a in autocors) {
    r <- if (a == 'all') all else run(a)
    v <- attr(r, 'values')
    j <- switch(a,
      positive   = which(va >  tol),
      negative   = which(va < -tol),
      `non-null` = which(abs(va) > tol),
      all        = seq_along(va) )
    stopifnot(
      length(j) == ncol(r),
      max(abs(v - va[j])) < 1e-14,
      same_spaces(as.matrix(r), A[, j, drop = FALSE], v) < 1e-12 )
    values[[a]] <- v
  }

  # pcnm() makes no allowance for rounding at the threshold: where computed
  # distances fall just above it (2.0000000000000004 for 2), it truncates
  # pairs that adespatial keeps, so those cases have no pcnm() results.
  if (with_pcnm) {
    stopifnot(!any(d > t & d <= t * (1 + 1e-6)))
    p <- pcnm(d, threshold = t)
  } else {
    stopifnot(any(d > t & d <= t * (1 + 1e-6)))
    p <- NULL
  }

  # To keep the file small, some cases keep only the MEMs (and Moran's I)
  # of "positive", the first columns of those of "all".
  if (!all_vectors) {
    j  <- seq_along(values$positive)
    A  <- A[, j, drop = FALSE]
    mi <- mi[j]
  }

  list(
    space = space, threshold = threshold, mst = mst, t = t,
    values = values, vectors = A, moran = mi,
    pcnm = if (!is.null(p)) list(values = p$values, threshold = p$threshold, vectors = NULL) )
}


cases <- list()

xy <- as.matrix(mite.xy)
rownames(xy) <- paste0('m', seq_len(nrow(xy)))
cases$mite <- one_case(xy)
cases$mite$pcnm$vectors <- unname(pcnm(dist(xy))$vectors)

cases$mite_t2 <- one_case(xy, 2, all_vectors = FALSE, with_pcnm = FALSE)

# Below the default threshold: a disconnected neighbour graph.
cases$mite_dist_t08 <- one_case(dist(xy), 0.8, all_vectors = FALSE, with_pcnm = FALSE)

# Duplicate positions.
cases$mite_dup <- one_case(
  rbind(xy, d1 = xy[1, ], d2 = xy[2, ], d3 = xy[40, ]), all_vectors = FALSE )

# A regular grid: equal eigenvalues, and eigenvalues of 0.
grid <- as.data.frame(expand.grid(x = 1:8, y = 1:5))
cases$grid    <- one_case(grid)
cases$grid_15 <- one_case(grid, 1.5)

# An irregular transect (collinear points).
cases$transect <- one_case(cbind(x = c(0, 1, 2, 3.5, 4, 6, 7, 7.5, 9, 10), y = 0))

saveRDS(cases, file.path('tests', 'tinytests', 'expected', 'dbmem.rds'), version = 2)
