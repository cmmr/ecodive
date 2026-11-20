# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# alr(), ilr() and the ILR basis helpers (helmert_basis(), sbp_basis(),
# ilr_basis()) follow alr(), ilr(), _gram_schmidt_basis(), sbp_basis() and
# _check_basis() from scikit-bio 0.7.4 (skbio/stats/composition/_base.py),
# Copyright (c) 2013--, scikit-bio development team, under the BSD-3 license
# given in src/log_ratio_transforms.c and inst/COPYRIGHTS.



#' Additive Log-Ratio Transform
#'
#' The log of every feature's abundance relative to one chosen feature, the
#' `denominator`: Aitchison's additive log ratio (alr).
#'
#' @inheritParams documentation
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**. Counts, proportions, and the output of a zero replacement
#'        are all accepted: each sample's scale cancels out of the log
#'        ratios.
#'
#' @param denominator  The feature to divide by: its name, or its (1-based)
#'        index among the features, or a logical vector with one `TRUE`.
#'        Required. Choose a feature that is present in every sample and
#'        abundant and stable across them, since every coordinate is relative
#'        to it.
#'
#' @param pseudocount  Added to every value before taking logs. `NULL` (the
#'        default) adds nothing when `counts` has no zeros, and otherwise
#'        uses half the smallest non-zero value, with a warning. See
#'        **Pseudocount**.
#'
#' @return A dense numeric `matrix` of the \eqn{D - 1} coordinates, in the
#'   orientation of `counts`: samples \eqn{\times} coordinates for `margin =
#'   1`, coordinates \eqn{\times} samples for `margin = 2` (and for a single
#'   sample given as a vector). Sample names are kept. Each coordinate is
#'   named after its numerator feature, in the order of `counts`; features
#'   without names give `alr1`, `alr2`, ..., numbered by their position in
#'   `counts` (so `denominator = 2` of 4 gives `alr1`, `alr3`, `alr4`).
#'   Attributes:
#'
#'   * `method`: `'alr'`.
#'   * `denominator`: the denominator feature's name, or its index when the
#'     features have no names.
#'   * `pseudocount`: the pseudocount used (`0` for none).
#'
#' @details
#' For a sample \eqn{x} of \eqn{D} features and denominator \eqn{x_d},
#' \deqn{\mathrm{alr}(x) = \left[\log\frac{x_j + c}{x_d + c}\right]_{j \neq d}}
#' with \eqn{c} the pseudocount (natural logs). The coordinates are easy to
#' interpret (everything relative to one feature) and are what multinomial
#' regression models use, but the transform is not an isometry: Euclidean
#' distances between alr coordinates are not Aitchison distances. For
#' distances, use [aitchison()]; for orthonormal coordinates, [ilr()]. The
#' centered log ratio is `normalize(counts, 'clr')` (or `norm = 'clr'` in any
#' function that takes `norm`).
#'
#' `alr(x + c, d)` here equals scikit-bio's `alr(x + c, ref_idx = d - 1)`
#' exactly. scikit-bio defaults to the first feature and `compositions::alr()`
#' to the last: ecodive has no default, because the choice shapes every
#' coordinate.
#'
#' The result feeds straight into [pca()], [rda()], [lda()], [envfit()] and
#' [euclidean()], which take negative values with `norm = 'none'`.
#'
#' **Same names in other packages**: `compositions::alr()` (the last part as
#' denominator) and `robCompositions::addLR()` compute the same transform.
#' Whichever package was attached last masks the other's `alr()`; write
#' `ecodive::alr()` to be sure.
#'
#' @inheritSection documentation Input Types
#' @inheritSection log_ratio_pseudocount_section Pseudocount
#' @inheritSection ml_section Machine learning
#'
#' @references
#' Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
#' Chapman & Hall.
#'
#' Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
#' normalizations outperform compositional data transformations in machine
#' learning applications. *Microbiome*, 12, 45.
#' \doi{10.1186/s40168-023-01747-z}
#'
#' @seealso [ilr()], [normalize()] (`'clr'` and `'rclr'`), [aitchison()]
#'
#' @export
#' @examples
#'     # Relative to the first feature
#'     x <- alr(ex_counts, denominator = 1, pseudocount = 1)
#'     x[1:3, 1:3]
#'     attr(x, 'denominator')
#'
#'     # The same as the log ratios themselves
#'     all.equal(
#'       unname(x[, 1]),
#'       unname(log((ex_counts[, 2] + 1) / (ex_counts[, 1] + 1))) )
#'
alr <- function (
    counts, denominator, margin = 1L, pseudocount = NULL, cpus = n_cpus() ) {

  if (missing(denominator))
    stop(
      '`denominator` is required: the name or index of the feature that the ',
      'others are divided by. There is no default, because every coordinate is ',
      'relative to it: choose a feature that is present in every sample, abundant ',
      'and stable.', call. = FALSE )

  info <- log_ratio_input('the additive log-ratio transform')
  den  <- feature_index(denominator, info$features, info$D, 'denominator')

  res <- .Call(
    C_log_ratio_coords, counts, margin, den - 1L, NULL, pseudocount, cpus )

  coords <- if (is.null(info$features)) paste0('alr', seq_len(info$D))[-den]
            else info$features[-den]
  dimnames(res) <- if (margin == 1L) list(info$samples, coords)
                   else              list(coords, info$samples)

  attr(res, 'method')      <- 'alr'
  attr(res, 'denominator') <- if (is.null(info$features)) den else info$features[[den]]
  attr(res, 'pseudocount') <- pseudocount
  res
}



#' Isometric Log-Ratio Transform
#'
#' Coordinates of each sample in an orthonormal basis of the Aitchison
#' simplex: Egozcue et al.'s isometric log ratio (ilr). Euclidean distances
#' between them are Aitchison distances, for any basis.
#'
#' @inheritParams alr
#'
#' @param basis  The ilr basis, one of:
#'
#'   * `NULL` (the default): the Helmert basis that scikit-bio's `ilr()`
#'     uses (built by Gram-Schmidt orthogonalization). Coordinate \eqn{j}
#'     is the balance of the first \eqn{j} features against feature
#'     \eqn{j + 1}: row \eqn{j} has \eqn{\frac{1}{j}\sqrt{j/(j+1)}} in
#'     columns \eqn{1, \dots, j}, \eqn{-\sqrt{j/(j+1)}} in column \eqn{j +
#'     1}, and zeros after.
#'   * A **sequential binary partition** (SBP): a \eqn{(D - 1) \times D}
#'     matrix of `1`, `-1` and `0`, one row per balance, `1` marking the
#'     features in its numerator and `-1` those in its denominator. Each
#'     row's \eqn{r} entries of `1` become \eqn{\sqrt{s / (r (r + s))}} and
#'     its \eqn{s} entries of `-1` become \eqn{-\sqrt{r / (s (r + s))}}, as
#'     scikit-bio's `sbp_basis()`. Every row needs both signs, and the
#'     partition must be sequential (its balances orthogonal).
#'   * An **orthonormal basis** in clr coordinates: a numeric \eqn{(D - 1)
#'     \times D} matrix with orthonormal rows that each sum to 0 (checked
#'     to \eqn{10^{-6}}, as scikit-bio checks orthonormality).
#'
#'   A matrix has one row per coordinate and one column per feature
#'   (scikit-bio's orientation). Its columns are matched to the features by
#'   name when both have names, and by position otherwise. Its row names
#'   name the coordinates.
#'
#' @return A dense numeric `matrix` of the \eqn{D - 1} coordinates, in the
#'   orientation of `counts`: samples \eqn{\times} coordinates for `margin =
#'   1`, coordinates \eqn{\times} samples for `margin = 2` (and for a single
#'   sample given as a vector). Sample names are kept. Coordinates are named
#'   by `rownames(basis)`, or `ilr1`, `ilr2`, .... Attributes:
#'
#'   * `method`: `'ilr'`.
#'   * `basis`: the \eqn{(D - 1) \times D} basis used (an SBP converted to
#'     its basis), with its columns in the order of the features of
#'     `counts`, so that `coords %*% basis` gives back
#'     `normalize(counts, 'clr')` (for `margin = 1`). For the default basis
#'     of more than 1000 features, whose matrix would take more than 8 MB,
#'     it is the string `'helmert'` instead.
#'   * `pseudocount`: the pseudocount used (`0` for none).
#'
#' @details
#' For a sample \eqn{x} with clr coordinates \eqn{\mathrm{clr}(x) = \log(x +
#' c) - \overline{\log(x + c)}} (natural logs, \eqn{c} the pseudocount) and a
#' basis \eqn{\Psi} with orthonormal rows,
#' \deqn{\mathrm{ilr}(x) = \mathrm{clr}(x)\,\Psi^T.}
#' The default basis is applied by a running sum, in \eqn{O(D)} per sample,
#' and other bases through their non-zero entries, so an SBP costs only its
#' non-zeros.
#'
#' Because the coordinates are orthonormal, `euclidean(ilr(x))` is
#' `aitchison(x)` (with the same pseudocount), whichever basis is used: **for
#' distances, use [aitchison()]**. `ilr()` is for analyses of the coordinates
#' themselves: regression, [lda()] or [rda()] without the clr's singular
#' covariance, and balances that can be interpreted, from an SBP. The coordinates feed straight into [pca()],
#' [rda()], [lda()], [envfit()] and [euclidean()] with `norm = 'none'`.
#'
#' `ilr(x + c)` here equals scikit-bio's `ilr(x + c)` (and
#' `ilr(x + c, sbp_basis(sbp))`) to rounding.
#'
#' **Same names in other packages**: `compositions::ilr()` and
#' `robCompositions::pivotCoord()` compute ilr coordinates with other default
#' bases (`compositions::ilr()`'s default is this one with the opposite
#' sign). Whichever package was attached
#' last masks the other's `ilr()`; write `ecodive::ilr()` to be sure.
#'
#' @inheritSection documentation Input Types
#' @inheritSection log_ratio_pseudocount_section Pseudocount
#' @inheritSection ml_section Machine learning
#'
#' @references
#' Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal,
#' C. (2003). Isometric logratio transformations for compositional data
#' analysis. *Mathematical Geology*, 35, 279–300.
#'
#' Egozcue, J. J., & Pawlowsky-Glahn, V. (2005). Groups of parts and their
#' balances in compositional data analysis. *Mathematical Geology*, 37,
#' 795–828.
#'
#' Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based
#' normalizations outperform compositional data transformations in machine
#' learning applications. *Microbiome*, 12, 45.
#' \doi{10.1186/s40168-023-01747-z}
#'
#' @seealso [alr()], [normalize()] (`'clr'` and `'rclr'`), [aitchison()]
#'
#' @export
#' @examples
#'     x <- ilr(ex_counts, pseudocount = 1)
#'     dim(x)
#'
#'     # Euclidean distances on ilr coordinates are Aitchison distances
#'     all.equal(
#'       as.vector(euclidean(x, norm = 'none')),
#'       as.vector(aitchison(ex_counts, pseudocount = 1)) )
#'
#'     # Balances from a sequential binary partition of 4 features
#'     sbp <- rbind(
#'       c(1,  1, -1, -1),
#'       c(1, -1,  0,  0),
#'       c(0,  0,  1, -1) )
#'     ilr(ex_counts[, 1:4], basis = sbp, pseudocount = 1)[1:3, ]
#'
ilr <- function (
    counts, basis = NULL, margin = 1L, pseudocount = NULL, cpus = n_cpus() ) {

  info <- log_ratio_input('the isometric log-ratio transform')
  b    <- ilr_basis(basis, info$D, info$features)

  res <- .Call(
    C_log_ratio_coords, counts, margin, -1L, b$csr, pseudocount, cpus )

  dimnames(res) <- if (margin == 1L) list(info$samples, b$names)
                   else              list(b$names, info$samples)

  attr(res, 'method')      <- 'ilr'
  attr(res, 'basis')       <- b$basis
  attr(res, 'pseudocount') <- pseudocount
  res
}



#' documentation
#'
#' @name log_ratio_pseudocount_section
#' @keywords internal
#'
#' @section Pseudocount:
#'   Logs of zero are undefined, so zeros need replacing first. The
#'   pseudocount is added to every value (not only to zeros), as in
#'   `normalize(norm = 'clr')`. With `pseudocount = NULL` (the default), no
#'   pseudocount is used when `counts` has no zeros; otherwise half the
#'   smallest non-zero value is, with a warning. To avoid a pseudocount,
#'   replace the zeros first with [replace_zeros()]: its result has none.
#'
#'   The choice of pseudocount is not neutral: it can distort the results,
#'   especially for sparse tables. See [aitchison()] for references.
NULL



# Checks the count table of a log-ratio transform in the calling function,
# and sets its `counts` (a matrix the C code reads), `margin`, `cpus` (if any) and
# `pseudocount`. Returns the table's sizes and names. `what` names the
# transform in error messages. `zero_ok` accepts `pseudocount = 0` (see
# log_pseudocount()); `dense_out = FALSE` (one value per sample) skips the
# warning about a large dense result. `check`, if given, is called with the
# feature names (or NULL) before the pseudocount is resolved, so that its
# errors come before the warning about zeros; its value is returned as
# element `check`.
log_ratio_input <- function (
    what, env = parent.frame(), zero_ok = FALSE, dense_out = TRUE, check = NULL ) {

  validate_counts(env)
  validate_margin(env)
  if (exists('cpus', envir = env, inherits = FALSE)) validate_cpus(env)

  counts <- env$counts
  margin <- env$margin

  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  pkg  <- get_matrix_package(counts)
  vals <- switch(pkg, base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  D <- as.integer(dims[[3L - margin]])
  n <- as.integer(dims[[margin]])
  if (D < 2L)
    stop('`counts` must have at least 2 features for ', what, '.', call. = FALSE)

  features <- if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]])
  checked  <- if (is.function(check)) check(features)

  env$pseudocount <- log_pseudocount(
    counts, env$pseudocount, what, constant_ok = TRUE, zero_ok = zero_ok )

  if (dense_out && as.double(n) * (D - 1) > 1e8)
    warning(
      sprintf('The result is a dense %s x %s matrix (%.1f GB).',
        n, D - 1, as.double(n) * (D - 1) * 8 / 1e9), call. = FALSE )

  list(
    D        = D,
    n        = n,
    features = features,
    samples  = if (length(dn[[margin]])) as.character(dn[[margin]]),
    check    = checked )
}



# The 1-based index of one feature, given by name, index or a logical
# vector with one TRUE.
feature_index <- function (x, features, D, arg) {

  msg <- sprintf('`%s` must be one feature: its name, its 1-based index, or a logical vector with one TRUE.', arg)

  if (is.factor(x)) x <- as.character(x)
  if (length(x) == 0L || anyNA(x)) stop(msg, call. = FALSE)

  if (is.character(x)) {
    if (length(x) != 1L) stop(msg, call. = FALSE)
    if (is.null(features))
      stop(sprintf('`counts` has no feature names; give `%s` as an index.', arg), call. = FALSE)
    i <- match(x, features)
    if (is.na(i))
      stop(sprintf('`%s` "%s" is not a feature of `counts`.', arg, x), call. = FALSE)
    return (i)
  }

  if (is.logical(x)) {
    if (length(x) != D || sum(x) != 1L)
      stop(msg, ' A logical `', arg, '` needs one value per feature.', call. = FALSE)
    return (which(x))
  }

  if (!is.numeric(x) || length(x) != 1L || x %% 1 != 0)
    stop(msg, call. = FALSE)
  if (x < 1 || x > D)
    stop(sprintf('`%s` = %s is not a feature index: `counts` has %d features.', arg, x, D), call. = FALSE)
  as.integer(x)
}



# The default ILR basis for D features, (D - 1) x D: scikit-bio's
# _gram_schmidt_basis(D), to the last bit.
helmert_basis <- function (D) {
  .Call(C_helmert_basis, as.integer(D))
}



# An orthonormal basis from a sequential binary partition: scikit-bio's
# sbp_basis(), in its arithmetic. Rows without both signs give NaN or Inf;
# ilr_basis() checks them first.
sbp_basis <- function (sbp) {
  psi <- matrix(0, nrow(sbp), ncol(sbp), dimnames = dimnames(sbp))
  for (i in seq_len(nrow(sbp))) {
    row   <- as.double(sbp[i, ])
    n_pos <- sum(row ==  1)
    n_neg <- sum(row == -1)
    psi[i, ] <- row * sqrt((n_neg / n_pos) ^ row / sum(abs(row)))
  }
  psi
}



# Checks and resolves an ILR `basis` for D features (named `features`, or
# NULL). `basis` is NULL (the Helmert basis), a sequential binary partition,
# or an orthonormal basis in clr coordinates, (D - 1) x D, columns matched to
# the features by name when both have names. Returns
#
#   basis  the (D - 1) x D basis, columns in the features' order, with
#          dimnames; for the default with D > 1000, the string 'helmert'.
#   names  the coordinates' names.
#   csr    for C_log_ratio_coords(): NULL (Helmert), or list(p, i, x), the
#          basis's non-zero entries by row, 0-based.
#   sbp    TRUE when `basis` was an SBP.
ilr_basis <- function (basis, D, features = NULL) {

  nc <- D - 1L

  if (is.null(basis)) {
    nm <- paste0('ilr', seq_len(nc))
    b  <- if (as.double(D) * nc <= 1e6) {
            structure(helmert_basis(D), dimnames = list(nm, features))
          } else 'helmert'
    return (list(basis = b, names = nm, csr = NULL, sbp = FALSE))
  }

  if (inherits(basis, 'data.frame') || !is.null(attr(basis, 'Dim')) && !is.matrix(basis))
    basis <- as.matrix(basis)

  if (!is.matrix(basis) || !(is.numeric(basis) || is.logical(basis)) || !all(is.finite(basis)))
    stop('`basis` must be NULL or a numeric matrix with finite values.', call. = FALSE)
  storage.mode(basis) <- 'double'

  if (!identical(dim(basis), c(nc, D)))
    stop(sprintf(
      '`basis` must be a %d x %d matrix (one row per coordinate, one column per feature), not %d x %d.%s',
      nc, D, nrow(basis), ncol(basis),
      if (identical(dim(basis), c(D, nc))) ' Transposed? Try t(basis).' else '' ),
      call. = FALSE )


  # Columns to features, by name when both have them.
  if (!is.null(colnames(basis)) && !is.null(features)) {
    if (anyDuplicated(colnames(basis)))
      stop('`basis` has duplicated column names.', call. = FALSE)
    if (!setequal(colnames(basis), features)) {
      miss <- setdiff(features, colnames(basis))
      extra <- setdiff(colnames(basis), features)
      stop(
        "`basis`'s column names don't match the features of `counts`.",
        if (length(miss))  paste0(' Not in `basis`: ',  paste(utils::head(miss,  5), collapse = ', '), if (length(miss) > 5)  ', ...', '.'),
        if (length(extra)) paste0(' Not in `counts`: ', paste(utils::head(extra, 5), collapse = ', '), if (length(extra) > 5) ', ...', '.'),
        call. = FALSE )
    }
    basis <- basis[, features, drop = FALSE]
  }
  else if (!is.null(features)) {
    colnames(basis) <- features
  }

  nm <- rownames(basis)
  if (is.null(nm)) nm <- paste0('ilr', seq_len(nc))
  rownames(basis) <- nm


  # A sequential binary partition?
  sbp <- all(basis == 1 | basis == -1 | basis == 0)
  if (sbp) {
    no_pos <- which(rowSums(basis ==  1) == 0)
    no_neg <- which(rowSums(basis == -1) == 0)
    bad    <- sort(union(no_pos, no_neg))
    if (length(bad))
      stop(sprintf(
        '`basis` is a sequential binary partition (entries 1, -1 and 0), but row %d has no %s. Each row needs both: 1 for the numerator, -1 for the denominator.',
        bad[[1]], if (bad[[1]] %in% no_pos) '1' else '-1' ), call. = FALSE )
    basis <- sbp_basis(basis)
  }


  # Orthonormal, as scikit-bio's _check_basis(), and in clr coordinates.
  eye <- diag(nc)
  if (!all(abs(tcrossprod(basis) - eye) < 1e-4 * eye + 1e-6))
    stop(
      if (sbp) '`basis` is not a valid sequential binary partition: its balances are not orthogonal.'
      else     '`basis` is not orthonormal.', call. = FALSE )
  if (!all(abs(rowSums(basis)) < 1e-6 * sqrt(D)))
    stop('Each row of `basis` must sum to 0: a basis of clr coordinates.', call. = FALSE)


  # Non-zero entries by row.
  bt  <- t(basis)
  nz  <- which(bt != 0)
  row <- (nz - 1L) %/% D
  csr <- list(
    p = c(0L, cumsum(tabulate(row + 1L, nbins = nc))),
    i = as.integer((nz - 1L) %% D),
    x = bt[nz] )

  list(basis = basis, names = nm, csr = csr, sbp = sbp)
}
