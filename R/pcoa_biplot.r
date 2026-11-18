# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Follows pcoa_biplot() from scikit-bio, copyright (c) 2013--,
# scikit-bio development team, under the BSD 3-Clause License; see
# src/pcoa_biplot.c and inst/COPYRIGHTS for its terms.



#' PCoA Biplot Scores
#'
#' Projects the features of a table, such as taxa, onto the axes of a
#' [pcoa()] ordination, so that the features driving each axis can be drawn
#' as arrows over the samples. Computes the same scores as `pcoa_biplot()`
#' in the Python package scikit-bio.
#'
#' @inherit pseudocount_section
#' @inheritParams documentation
#' @inheritSection documentation Input Types
#'
#' @param ord   The result of [pcoa()], with any `correction`.
#'
#' @param counts   A table of the samples in `ord` (samples \eqn{\times}
#'        features when `margin = 1`): usually the feature table that `ord`'s
#'        distances came from, or environmental variables. Any of the forms
#'        described in "Input Types" below. Must not contain `NA`s.
#'
#' @param norm   How to transform each sample's values before projecting
#'   them:
#'
#'   * `'none'`: As given, as in scikit-bio.
#'   * `'percent'`: Relative abundance (sample abundances sum to 1).
#'   * `'chord'`: Each sample scaled to unit length (Euclidean norm 1).
#'   * `'binary'`: Presence/absence (each count is either 0 or 1).
#'   * `'clr'`: Centered log ratio.
#'   * `'rclr'`: Robust centered log ratio. Zeros stay 0.
#'
#'   Default: `'none'`.
#'
#' @return A numeric matrix with one row per feature, named after the
#'   features, and one column per axis of `ord`, named as `ord`'s columns
#'   (`"PCoA1"`, `"PCoA2"`, ...). Each row is the tip of that feature's
#'   arrow. Its attribute `method` is `"pcoa_biplot"`.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{Y} the feature table (\eqn{n} samples by \eqn{p} features) with
#' each feature centered, \eqn{Z} the PCoA coordinates with each axis
#' standardized (centered and divided by its standard deviation), and
#' \eqn{\lambda_j} axis \eqn{j}'s eigenvalue (`attr(ord, 'eig')`), the scores
#' are Legendre and Legendre's (1998) equations 9.44 and 9.55, as in
#' scikit-bio:
#' \deqn{S = \frac{1}{n-1} Y^T Z, \qquad
#'       U_{proj} = \sqrt{n-1}\, S\, \mathrm{diag}(\lambda_j^{-1/2})}{S = Y' Z / (n - 1),   U_proj = sqrt(n - 1) S diag(lambda_j^(-1/2))}
#'
#' A PCoA axis \eqn{x_j} has \eqn{\|x_j\|^2 = \lambda_j}, so this is
#' \eqn{Y^T x_j / \lambda_j}: the covariance of each feature with the axis,
#' divided by the axis's variance. That is the least-squares slope of the
#' feature on the axis: how much the feature changes (in the units of
#' `counts`, after `norm`) per unit along the axis. The axes are
#' uncorrelated, so it is also the feature's coefficient for that axis in a
#' multiple regression on all of `ord`'s axes. An arrow points in the
#' direction in which the feature increases, and is longer the faster it
#' does.
#'
#' The work is \eqn{O(\mathrm{nnz} \times k)}{O(nnz * k)} for a table with
#' \eqn{\mathrm{nnz}}{nnz} non-zero values and \eqn{k} axes: sparse tables stay
#' sparse, and features are processed in parallel. Results are identical
#' from run to run, for any number of `cpus`, and for dense and sparse
#' versions of the same table.
#'
#'
#' **Relation to `envfit()`**
#'
#' [envfit()] fits each feature to the ordination too. Its arrows point in
#' the same direction as these scores, but are scaled to unit length (and
#' then often by \eqn{\sqrt{r^2}}), and it tests each fit by permutation.
#' `pcoa_biplot()` gives the unscaled slopes, with no permutations, so it is
#' far faster for large tables.
#'
#'
#' **Drawing a biplot**
#'
#' The scores are in units of `counts` per ordination unit, so arrows
#' usually need a common scale factor to fit the plot. Pick the features
#' with the longest arrows on the axes drawn, as in the example below.
#'
#'
#' **Which table to use**
#'
#' Pass the table in the form the distances saw: the rarefied table if the
#' distances came from one (see [rarefy()]), and the same `norm` as the
#' distances used. Raw counts give long arrows to abundant taxa, while
#' `norm = 'percent'` or `'clr'` put taxa on a more even footing. With a
#' phylogenetic metric, such as UniFrac, the axes reflect the tree as well as
#' the counts; the scores still show how each feature varies along the axes,
#' but not how much it contributes to the distances. Environmental variables
#' in different units should be standardized first, e.g. with `scale()`.
#'
#'
#' **Samples**
#'
#' Samples are matched by name, as in [envfit()]: `counts` must have every
#' sample of `ord` (`rownames(ord)`), in any order; any others are ignored.
#' If either has no sample names, they are matched by position, and must be
#' equal in number.
#'
#'
#' **Differences from scikit-bio**
#'
#' * scikit-bio adds the scores to its ordination object; `pcoa_biplot()`
#'   returns them as a matrix.
#' * scikit-bio multiplies the uncentered table. The two are the same in
#'   exact arithmetic (the columns of \eqn{Z} sum to 0); centering keeps
#'   rounding error proportional to each feature's spread rather than to its
#'   mean, and a feature with the same value in every sample scores exactly
#'   0.
#' * `norm` optionally transforms the table first; scikit-bio uses it as
#'   given.
#' * scikit-bio replaces `NaN` scores (from `NaN`s in the table) with 0;
#'   here a table with `NA`, `NaN`, or infinite values is an error.
#' * Axes that `pcoa()` dropped (eigenvalues that aren't positive) are
#'   absent here; scikit-bio returns columns of zeros for them.
#'
#' As in scikit-bio, an axis whose coordinates have zero standard deviation
#' is divided by 1.
#'
#' [nmds()], [tsne()], and [umap()] have no eigenvalues, so this projection
#' doesn't apply to them; use [envfit()].
#'
#' @references
#' Legendre, P., & Legendre, L. (1998). *Numerical Ecology* (2nd English
#' ed., chap. 9). Elsevier.
#'
#' scikit-bio development team. scikit-bio: A Bioinformatics Library for
#' Data Scientists, Students, and Developers. \url{https://scikit.bio}
#'
#' @seealso [pcoa()], [envfit()], [pca()]
#'
#' @export
#' @examples
#'     dm     <- bray(iris[, 1:4])
#'     ord    <- pcoa(dm)
#'     scores <- pcoa_biplot(ord, iris[, 1:4])
#'     scores
#'
#'     # Scale the arrows to the samples' spread, and draw them.
#'     xy <- scores[, 1:2] * 0.8 * max(abs(ord[, 1:2])) / max(abs(scores[, 1:2]))
#'     plot(ord, col = iris$Species)
#'     arrows(0, 0, xy[, 1], xy[, 2], length = 0.1)
#'     text(xy, labels = rownames(xy), pos = 3)
#'
pcoa_biplot <- function (
    ord, counts, margin = 1L, norm = 'none', pseudocount = NULL,
    cpus = n_cpus() ) {

  validate_cpus()


  #--------------------------------------------------
  # The ordination.
  #--------------------------------------------------
  eig <- attr(ord, 'eig', exact = TRUE)

  ok <- is.matrix(ord) && is.numeric(ord) &&
    identical(attr(ord, 'method', exact = TRUE), 'pcoa') &&
    is.numeric(eig) && length(eig) == ncol(ord) && nrow(ord) >= 1 &&
    ncol(ord) <= nrow(ord) - 1 && all(is.finite(eig)) && all(eig > 0) &&
    all(is.finite(ord))

  if (!isTRUE(ok))
    stop('`ord` must be a result from pcoa().', call. = FALSE)

  n      <- nrow(ord)
  k      <- ncol(ord)
  labels <- rownames(ord)


  #--------------------------------------------------
  # The table, its samples matched to ord's rows.
  #--------------------------------------------------
  validate_counts()
  validate_margin()
  validate_norm()
  validate_pseudocount()

  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  ns  <- as.integer(dims[[1]])
  idx <- match_rows(labels, n, dn[[1]], ns, ref = '`ord`', what = '`counts`', unit = 'sample')

  if (is.matrix(counts) && !is.double(counts)) {
    if (!is.numeric(counts) && !is.logical(counts))
      stop('`counts` must be numeric.', call. = FALSE)
    storage.mode(counts) <- 'double'
  }

  map <- rep(-1L, ns)
  map[idx] <- seq_along(idx) - 1L

  feat_names <- dn[[2]]


  #--------------------------------------------------
  # Scores: Yc' P, with P = Z diag(eig^-1/2) / sqrt(n - 1).
  # A PCoA axis's standard deviation is sqrt(eig / (n - 1)),
  # but as in scikit-bio it is taken from the coordinates.
  #--------------------------------------------------
  if (k == 0) {
    scores <- matrix(0, dims[[2]], 0)
  }
  else {
    x  <- matrix(as.double(ord), n, k)
    xc <- x - rep(colMeans(x), each = n)
    nm <- sqrt(colSums(xc^2))
    nm[nm == 0] <- sqrt(n - 1)  # sd 0 -> 1
    P  <- xc / rep(nm * sqrt(eig), each = n)

    scores <- .Call(C_pcoa_biplot, counts, margin, norm, pseudocount, P, map, cpus)

    if (!all(is.finite(scores)))
      stop(
        '`counts` gives scores beyond the range of double precision numbers; ',
        'rescale it.', call. = FALSE )
  }

  dimnames(scores) <- list(feat_names, colnames(ord))
  attr(scores, 'method') <- 'pcoa_biplot'
  scores
}
