# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Normalize a Count Table
#'
#' Applies one of the transformations that [beta_div()], [pca()], and
#' [envfit()] offer through their `norm` argument, and returns the result.
#' The values are exactly those that these functions compute internally.
#'
#' @inheritParams documentation
#'
#' @param norm  The transformation:
#'
#'   * `'percent'`: Relative abundance: each sample's counts divided by
#'     their sum (so they sum to 1).
#'   * `'chord'`: Each sample's counts divided by their Euclidean norm
#'     (so their squares sum to 1).
#'   * `'binary'`: Presence/absence: 1 for each non-zero count.
#'   * `'clr'`: Centered log ratio: \eqn{\log(x + pseudocount)} minus the
#'     sample's mean of those logs, over every feature. Zeros get a value
#'     too, so the result is dense. See **Pseudocount** below.
#'   * `'rclr'`: Robust centered log ratio: \eqn{\log(x)} minus the mean
#'     log of the sample's non-zero counts, for non-zero counts only; zeros
#'     stay 0.
#'   * `'none'`: No transformation.
#'
#' @return The normalized table, with the same orientation (samples in rows
#'   or columns, as given by `margin`) and dimnames as `counts`, and in the
#'   same form where zeros stay zero:
#'
#'   * A base `matrix` gives a numeric `matrix`.
#'   * A \pkg{Matrix} `dgCMatrix`, `dgTMatrix`, or `dgeMatrix` gives the
#'     same class, and a \pkg{slam} `simple_triplet_matrix` gives a
#'     `simple_triplet_matrix`. Zeros are not stored.
#'   * A `phyloseq`, `rbiom`, or `SummarizedExperiment` object gives its
#'     normalized count matrix (not the object), in that object's
#'     orientation.
#'
#'   `norm = 'clr'` always gives a dense base `matrix`, since every zero
#'   becomes a (usually negative) value, which can take far more memory than
#'   a sparse input. [beta_div()], [pca()], and [envfit()] accept `norm =
#'   'clr'` themselves and never make that dense copy, so prefer them where
#'   you can. `normalize()` warns when the result would have more than
#'   \eqn{10^8} values (800 MB).
#'
#' @inheritSection pseudocount_section Pseudocount
#'
#' @details
#' Samples whose counts are already proportions (summing to 1, within
#' \eqn{10^{-6}}) are left as they are by `norm = 'percent'`, as in the
#' functions with a `norm` argument. Stored zeros in sparse input are
#' dropped, and repeated entries of a `dgTMatrix` are summed.
#'
#' **Masking**: \pkg{BiocGenerics}, which is attached along with
#' \pkg{phyloseq} or \pkg{SummarizedExperiment}, also has a `normalize()`.
#' Whichever package was attached last masks the other's, so call
#' `ecodive::normalize()` to be sure of getting this one.
#'
#' @seealso [beta_div()], [pca()], [envfit()], [rarefy()]
#'
#' @export
#' @examples
#'     # Relative abundances
#'     normalize(ex_counts, 'percent')[, 1:4]
#'
#'     # Centered log ratios, with a pseudocount of 1 for the zeros
#'     clr <- normalize(ex_counts, 'clr', pseudocount = 1)
#'     rowSums(clr)  # 0 for every sample
#'
#'     # The same values as pca()'s own normalization
#'     a <- pca(ex_counts, norm = 'clr', pseudocount = 1, k = 2)
#'     b <- pca(clr, k = 2)
#'     all.equal(a, b)
#'
normalize <- function (
    counts, norm, margin = 1L, pseudocount = NULL, cpus = n_cpus() ) {

  if (missing(norm))
    stop("`norm` is required: 'percent', 'chord', 'binary', 'clr', 'rclr', or 'none'.",
         call. = FALSE)
  if (!is.character(norm) || length(norm) != 1 || is.na(norm))
    stop("`norm` must be one of 'percent', 'chord', 'binary', 'clr', 'rclr', or 'none'.",
         call. = FALSE)

  validate_counts()

  # As pca() and envfit(): NA would pass through, or become a presence.
  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()

  pkg <- get_matrix_package(counts)
  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  clr   <- identical(norm, NORM_CLR)
  dense <- clr || pkg == 'base' || inherits(counts, 'dgeMatrix')

  if (clr && prod(as.double(dims)) > 1e8)
    warning(
      sprintf('normalize(norm = "clr") gives a dense %s x %s matrix (%.1f GB). ',
        dims[[1]], dims[[2]], prod(as.double(dims)) * 8 / 1e9),
      'beta_div(), pca() and envfit() accept norm = "clr" without it.', call. = FALSE )

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'

  res <- .Call(
    C_normalize, counts, margin, norm, as.double(pseudocount), dense, cpus )

  if (dense) {
    dimnames(res) <- if (is.null(dn) || all(vapply(dn, is.null, logical(1)))) NULL else dn
    if (!clr && inherits(counts, 'dgeMatrix')) {
      counts@x <- as.vector(res)
      counts@factors <- list()
      return (counts)
    }
    return (res)
  }

  if (pkg == 'slam')
    return (structure(
      list(
        i = res$i + 1L, j = res$j + 1L, v = res$x,
        nrow = as.integer(dims[[1]]), ncol = as.integer(dims[[2]]),
        dimnames = counts$dimnames ),
      class = 'simple_triplet_matrix' ))

  # dgCMatrix or dgTMatrix: same object, new entries.
  if (inherits(counts, 'dgTMatrix')) {
    counts@i <- res$i
    counts@j <- res$j
  }
  else {
    counts@i <- res$i
    counts@p <- c(0L, cumsum(tabulate(res$j + 1L, nbins = dims[[2]])))
  }
  counts@x       <- res$x
  counts@factors <- list()
  counts
}
