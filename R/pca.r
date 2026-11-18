# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Principal Component Analysis (PCA)
#'
#' Ordinates samples from their feature abundances, optionally after a
#' normalization such as the centered log ratio (CLR, giving an "Aitchison
#' PCA"). Gives the same scores, loadings, and variances as `stats::prcomp()`
#' on the normalized table, but computes only the first few axes, straight
#' from sparse input, which takes far less time and memory for large
#' datasets.
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#'
#' @param norm   Normalize the incoming counts before the PCA. Options are:
#'
#'   * `'none'`: No transformation.
#'   * `'percent'`: Relative abundance (sample abundances sum to 1).
#'   * `'chord'`: Each sample scaled to unit length (Euclidean norm 1).
#'   * `'hellinger'`: Square root of relative abundance (the Hellinger
#'     transformation).
#'   * `'binary'`: Presence/absence (each count is either 0 or 1).
#'   * `'clr'`: Centered log ratio, giving an Aitchison PCA.
#'   * `'rclr'`: Robust centered log ratio. Zeros stay 0.
#'
#'   Default: `'none'`.
#'
#' @param k   The number of axes to compute. Default: `3`
#'
#' @param scale   Scale each feature to unit variance (after `norm`, before
#'        the PCA), as in `prcomp(scale. = TRUE)`. This gives every feature
#'        equal weight, so that a PCA of the correlation matrix rather than
#'        the covariance matrix results. Default: `FALSE`
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   named `"PC1"`, `"PC2"`, etc., as described in "Ordination result"
#'   below. Its attributes are:
#'
#'   * `method`: `"pca"`.
#'   * `eig`: The variance of each axis: `prcomp()`'s `sdev^2`.
#'   * `prop`: The proportion of the total variance that each axis explains,
#'     `eig / trace`.
#'   * `trace`: The total variance: the sum of all features' variances (after
#'     `norm` and `scale`), which is also the sum of all axes' `eig`.
#'   * `loadings`: A matrix with one row per feature and one column per axis:
#'     the unit-length principal axes, `prcomp()`'s `rotation`.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{Z} the normalized table (samples in rows), \eqn{J} the centering
#' matrix, and \eqn{D} the diagonal matrix of the features' standard
#' deviations (or the identity matrix when `scale = FALSE`), PCA is the
#' singular value decomposition \eqn{JZD^{-1} = USV^T}. The scores are
#' \eqn{US}, the loadings \eqn{V}, and each axis's variance (`eig`) is
#' \eqn{s^2 / (n - 1)}, for \eqn{n} samples.
#'
#' `prcomp()` builds the whole centered table and decomposes it in full.
#' `pca()` finds only the first `k` axes, by block Lanczos iteration, which
#' only needs products of the table with a few vectors at a time. Those are
#' computed in parallel from the non-zero counts alone: the zeros of each
#' sample all share one value after normalization (0, or
#' \eqn{\log(pseudocount)} minus the sample's mean log for `norm = 'clr'`),
#' and centering only subtracts each feature's mean, so a sparse table is
#' never made dense.
#'
#'
#' **Relation to `pcoa()`**
#'
#' The scores are the same as those of [pcoa()] on the Euclidean distances
#' between the normalized samples: `pca(counts)` and
#' `pcoa(euclidean(counts))`, or `pca(counts, norm = 'clr')` and
#' `pcoa(aitchison(counts))`, give the same coordinates. PCoA's eigenvalues
#' are \eqn{s^2}, so `pca()`'s `eig` is `pcoa()`'s divided by \eqn{n - 1};
#' `prop` is the same for both. What `pca()` adds is the `loadings`: how much
#' each feature contributes to each axis, e.g. for a biplot.
#'
#'
#' **Robust Aitchison PCA**
#'
#' `norm = 'rclr'` with `pca()` is a plain PCA of the rclr-transformed table,
#' with its zeros left as 0. DEICODE and gemelli's Robust Aitchison PCA
#' instead treat those zeros as missing values and fill them in by low-rank
#' matrix completion (OptSpace) before the decomposition, so their results
#' differ.
#'
#'
#' **Axis signs and reproducibility**
#'
#' A principal axis's sign is arbitrary, so different software can return
#' mirror images of the same axis. As [pcoa()] does, `pca()` flips each axis
#' so that its largest score in absolute value is positive, and flips its
#' loadings to match. Results are identical from run to run and for any
#' number of `cpus`.
#'
#' `k` larger than the number of samples minus one, or the number of
#' features, is reduced to the smaller of those. Axes with variances at
#' rounding-error level, below `sqrt(.Machine$double.eps)` times the first
#' one, are dropped with a warning. With `scale = TRUE`, features with
#' constant values (such as features absent from every sample) are left out,
#' where `prcomp()` would stop with an error.
#'
#' @references
#' Pearson, K. (1901). On lines and planes of closest fit to systems of
#' points in space. *The London, Edinburgh, and Dublin Philosophical Magazine
#' and Journal of Science*, 2(11), 559-572. \doi{10.1080/14786440109462720}
#'
#' Aitchison, J. (1983). Principal component analysis of compositional data.
#' *Biometrika*, 70(1), 57-65. \doi{10.1093/biomet/70.1.57}
#'
#' @inheritSection ordination_section Ordination result
#' @inheritSection documentation Input Types
#'
#' @seealso [pcoa()], [dbrda()], [dist_ord()], [aitchison()], `stats::prcomp()`
#'
#' @export
#' @examples
#'     ord <- pca(iris[, 1:4])
#'     plot(ord, col = iris$Species)
#'     head(ord)
#'     str(ord)
#'
pca <- function (
    counts,
    margin      = 1L,
    norm        = 'none',
    pseudocount = NULL,
    k           = 3L,
    scale       = FALSE,
    cpus        = n_cpus() ) {

  validate_counts()
  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()
  k <- ord_k(k)

  if (!isTRUE(scale) && !isFALSE(scale))
    stop('`scale` must be TRUE or FALSE.', call. = FALSE)


  # Dimensions and names, with samples first.
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                     dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol);     dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                      dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  n <- as.integer(dims[[1]])
  p <- as.integer(dims[[2]])

  k <- min(k, n - 1L, p)

  if (k < 1) {

    # One sample: no axes.
    scores   <- matrix(numeric(0), nrow = n, ncol = 0)
    loadings <- matrix(numeric(0), nrow = p, ncol = 0)
    eig      <- numeric(0)
    total    <- 0

  } else {

    res <- .Call(C_pca, counts, margin, norm, pseudocount, k, scale, cpus)

    if (!res$converged) # nocov start
      warning('PCA did not converge; results may be inexact.', call. = FALSE)
    # nocov end

    keep     <- ord_keep(res$eig, k)
    scores   <- res$scores[, keep, drop = FALSE]
    loadings <- res$loadings[, keep, drop = FALSE]
    eig      <- res$eig[keep]
    total    <- res$total
  }

  axes <- if (ncol(scores) > 0) paste0('PC', seq_len(ncol(scores))) else NULL
  dimnames(loadings) <- list(dn[[2]], axes)

  denom <- max(1L, n - 1L)

  ord_result(
    scores, 'pca', dn[[1]], prefix = 'PC',
    eig      = eig / denom,
    prop     = if (total > 0) eig / total else eig,
    trace    = total / denom,
    loadings = loadings )
}
