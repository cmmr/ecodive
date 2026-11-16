# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Helpers shared by the ordination methods (pcoa, pca, dbrda, nmds, tsne,
# umap) and by permanova().
#
#   dist_vector()    a dist object or symmetric matrix -> its packed lower
#                    triangle, checked for NA, negative and huge values.
#   formula_dm()     splits vegan-style `dm ~ terms` calls.
#   ord_arg()        checks a numeric scalar argument.
#   ord_k()          checks `k`.
#   ord_keep()       which eigenvalues are positive, with a warning if not all.
#   ord_result()     the shared result: a matrix of class "ecodive_ord".



# The packed lower triangle of a dist object or symmetric matrix, as a
# double vector `dm`, with the number of samples `n`, their `labels`, and
# the dist object's `maxdist` attribute (from vegan::vegdist()), if any.
#
# Every method rejects NA, NaN, infinite, and negative distances. Negative
# values within rounding error of 0 (down to -1e-12 times the largest
# distance) are set to 0. Distances above 1e150 are rejected, since their
# squares would overflow.
dist_vector <- function (dm) {

  maxdist <- NULL

  if (inherits(dm, 'formula'))
    stop(
      '`dm` must be a dist object or a symmetric matrix, not a formula. ',
      'Give the formula separately, e.g. permanova(dm, data, ~ group).', call. = FALSE )

  if (inherits(dm, 'dist')) {

    n       <- as.integer(attr(dm, 'Size', exact = TRUE))
    labels  <- attr(dm, 'Labels', exact = TRUE)
    maxdist <- attr(dm, 'maxdist', exact = TRUE)

    if (length(n) != 1 || is.na(n) || n < 0 || length(dm) != n * (n - 1) / 2)
      stop('`dm` is not a valid dist object.', call. = FALSE)
  }
  else {

    dm <- tryCatch(as.matrix(dm), error = function (e) NULL)
    n  <- nrow(dm)

    if (!is.numeric(dm) || is.null(n) || n != ncol(dm) || !isSymmetric(unname(dm)))
      stop('`dm` must be a dist object or a symmetric numeric matrix.', call. = FALSE)

    labels <- if (is.null(rownames(dm))) colnames(dm) else rownames(dm)
    dm     <- dm[lower.tri(dm)]
  }

  if (!is.numeric(dm))
    stop('`dm` must be numeric.', call. = FALSE)
  if (!is.double(dm)) storage.mode(dm) <- 'double'
  attributes(dm) <- NULL

  if (length(dm)) {

    if (!all(is.finite(dm)))
      stop('`dm` must not contain NA, NaN, or infinite values.', call. = FALSE)

    top <- max(dm)
    if (top > 1e150)
      stop('`dm` has distances too large to square (over 1e150); rescale them.', call. = FALSE)

    low <- min(dm)
    if (low < 0) {
      if (low < -1e-12 * max(top, 0))
        stop('`dm` must not contain negative distances.', call. = FALSE)
      dm[dm < 0] <- 0  # rounding error
    }
  }

  list(dm = dm, n = n, labels = labels, maxdist = maxdist)
}



# vegan-style calls give the distances as a formula's left-hand side:
# permanova(dm ~ group, data). Returns list(dm, formula): the left-hand side
# evaluated in the formula's environment, and a one-sided formula. Anything
# else is returned unchanged.
formula_dm <- function (dm, formula) {

  if (!inherits(dm, 'formula') || length(dm) != 3)
    return (list(dm = dm, formula = formula))

  if (!is.null(formula))
    stop(
      'Give either a formula with the distances on its left-hand side (`dm`), ',
      'or `dm` and `formula`, not both.', call. = FALSE )

  env <- environment(dm)
  lhs <- eval(dm[[2]], envir = if (is.null(env)) globalenv() else env)
  rhs <- dm[-2]
  environment(rhs) <- if (is.null(env)) globalenv() else env

  list(dm = lhs, formula = rhs)
}



# Checks a numeric scalar argument against [lo, hi]; returns it as an
# integer if `whole`, else as a double. `msg` is the error message.
ord_arg <- function (x, msg, whole = FALSE, lo = 0, hi = Inf) {

  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && x >= lo && x <= hi &&
    (!whole || x %% 1 == 0)

  if (!isTRUE(ok)) stop(msg, call. = FALSE)

  if (whole) as.integer(min(x, .Machine$integer.max)) else as.double(x)
}


# Checks `k`, the number of axes.
ord_k <- function (k) {
  ord_arg(k, '`k` must be a whole number of at least 1.', whole = TRUE, lo = 1)
}



# Which of `eig` (the first `k` eigenvalues, largest first) to keep: the
# positive ones, ignoring rounding-error ones below sqrt(eps) times the first.
# Warns if any are dropped, as cmdscale() does.
ord_keep <- function (eig, k = length(eig)) {

  keep <- eig > max(0, sqrt(.Machine$double.eps) * eig[1])
  keep[is.na(keep)] <- FALSE

  if (sum(keep) < k)
    warning(sprintf('only %d of the first %d eigenvalues are > 0', sum(keep), k), call. = FALSE)

  keep
}



# tsne() and umap()'s `init`: one of `choices`, or a numeric matrix (or data
# frame). Returns list(init = the choice, or 'matrix'; mtx = the matrix, or
# NULL).
ord_init_arg <- function (init, choices) {

  if (is.matrix(init) || is.data.frame(init))
    return (list(init = 'matrix', mtx = init))

  if (!is.character(init) || length(init) != 1 || !init %in% choices)
    stop('`init` must be ', paste0('"', choices, '"', collapse = ', '),
         ', or a numeric matrix.', call. = FALSE)

  list(init = init, mtx = NULL)
}


# A given matrix of starting coordinates, checked to be finite and n x k.
# Returns a plain double matrix.
ord_init_matrix <- function (mtx, n, k) {

  mtx <- unname(as.matrix(mtx))
  if (!is.numeric(mtx) || !identical(dim(mtx), c(n, k)) || !all(is.finite(mtx)))
    stop('`init` must be a numeric matrix with ', n, ' rows and ', k,
         ' columns, without NA or infinite values.', call. = FALSE)

  storage.mode(mtx) <- 'double'
  mtx
}


# The first k PCoA axes of `dm` (packed, n samples), as starting coordinates
# for tsne() and umap(), which scale them their own ways. NULL if there are
# fewer than k axes with positive eigenvalues.
pcoa_start <- function (dm, n, k, cpus) {

  if (k > n - 1) return (NULL)

  res <- .Call(C_pcoa, dm, n, k, cpus)
  eig <- res$eig
  if (length(eig) < k || !all(eig > max(0, sqrt(.Machine$double.eps) * eig[[1]])))
    return (NULL)

  res$points
}



# Flips each column of `x` so that its first element whose magnitude is
# within a relative 1e-8 of the largest is positive, as lz_sign_index() does
# in C. Allowing for rounding keeps last-bit differences between platforms
# from breaking a tie the other way.
ord_flip <- function (x) {
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (length(v) == 0) next
    a <- abs(v)
    i <- which(a >= max(a) * (1 - 1e-8))[[1]]
    if (v[[i]] < 0) x[, j] <- -v
  }
  x
}



# The shared ordination result: `points`, an n x k numeric matrix, with one
# row per sample (named `labels`) and one column per axis, named `prefix`
# followed by 1, 2, ..., or else `cols`. Of class "ecodive_ord". Attributes
# `method`, `eig` and `prop` come first (NA_real_ for each axis by default),
# then those in `...`.
ord_result <- function (
    points, method, labels, prefix = NULL, cols = NULL,
    eig = NULL, prop = NULL, ... ) {

  points <- matrix(as.double(points), nrow = NROW(points), ncol = NCOL(points))
  k      <- ncol(points)

  if (!is.null(prefix)) cols <- paste0(prefix, seq_len(k))
  if (k == 0) cols <- NULL

  dimnames(points) <- list(labels, cols)

  if (is.null(eig))  eig  <- rep(NA_real_, k)
  if (is.null(prop)) prop <- rep(NA_real_, k)

  structure(
    .Data  = points,
    method = method,
    eig    = as.double(eig),
    prop   = as.double(prop),
    ...,
    class  = c('ecodive_ord', 'matrix', 'array') )
}



#' Print an ordination result
#'
#' Prints the coordinates, then one line on the axes and the names of the
#' other attributes, rather than every attribute in full (such as `pca()`'s
#' loadings).
#'
#' @param x  An ordination result, from [ordination()], [pcoa()], [pca()],
#'        [dbrda()], [nmds()], [tsne()], or [umap()].
#'
#' @param ...  Passed on to `print()` for the matrix of coordinates.
#'
#' @return `x`, invisibly.
#'
#' @keywords internal
#' @export
print.ecodive_ord <- function (x, ...) {

  att <- attributes(x)
  m   <- unclass(x)
  attributes(m) <- att[c('dim', 'dimnames')]
  print(m, ...)

  k    <- ncol(x)
  prop <- att$prop
  line <- sprintf('\n%s: %d sample%s, %d ax%s', att$method,
    nrow(x), if (nrow(x) == 1) '' else 's', k, if (k == 1) 'is' else 'es')
  if (k > 0 && !anyNA(prop))
    line <- paste0(line, sprintf('; %% of variance: %s',
      paste(format(round(100 * prop, 1), nsmall = 1, trim = TRUE), collapse = ', ')))
  cat(line, '\n', sep = '')

  more <- setdiff(names(att), c('dim', 'dimnames', 'class'))
  cat('Attributes: ', paste(more, collapse = ', '), '\n', sep = '')

  invisible(x)
}
