# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The rotation follows scipy.linalg.orthogonal_procrustes(), and the
# symmetric statistic scipy.spatial.procrustes(), from SciPy 1.18.1
# (BSD-3-Clause; see src/procrustes.c and inst/COPYRIGHTS).



#' Procrustes Rotation of One Ordination onto Another
#'
#' Rotates, reflects, and optionally dilates (uniformly scales) one
#' ordination, `ord`, to fit another of the same samples, `target`, as
#' closely as possible, and reports how close the fit is: the Procrustes sum
#' of squares, \eqn{m^2}. Use it to ask how similar two pictures of the same
#' samples are: the same samples under two beta diversity metrics, say, or
#' from 16S and from shotgun metagenomics. Gives the same results as
#' `vegan::procrustes()`. [procrustes_test()] tests whether the fit is better
#' than chance (PROTEST).
#'
#' @inheritParams documentation
#'
#' @param ord  The ordination to rotate onto `target`: a numeric matrix (or
#'        data frame) with one row per sample and one column per axis, such
#'        as from [pcoa()], [nmds()], [pca()], or [umap()]; or a `dist`
#'        object, such as from [beta_div()], which is first ordinated with
#'        [pcoa()].
#'
#' @param target  The ordination that `ord` is fitted to, in any of the
#'        forms `ord` takes. When both `ord` and `target` have sample (row)
#'        names, `target`'s rows are matched to `ord`'s by name; otherwise
#'        they must be in the same order.
#'
#' @param axes  Which axes (columns) of `ord` and `target` to use: whole
#'        numbers, as in [envfit()]. When one of the two has fewer columns
#'        than asked for, it is padded with columns of zeros, as in
#'        `vegan::procrustes()`. Default: `NULL`, the first
#'        `min(ncol(ord), ncol(target))` axes.
#'
#' @param n_axes  When `ord` or `target` is a `dist` object, the number of
#'        PCoA axes to compute from it, as in [pcoa()]. Ignored otherwise.
#'        Default: `3`.
#'
#' @param symmetric  Use the symmetric Procrustes statistic: scale both `ord`
#'        and `target` to unit sum of squares first, so that \eqn{m^2} is
#'        between 0 and 1 and is the same whichever of the two is rotated.
#'        As `vegan::procrustes()`'s `symmetric`. Default: `FALSE`.
#'
#' @param dilate  Allow a uniform dilation (scaling) of `ord`, as well as the
#'        rotation and reflection. `vegan::procrustes()`'s `scale`, renamed
#'        because `scale` means "standardise each feature" elsewhere in
#'        ecodive. Default: `TRUE`.
#'
#' @return The rotated `ord`: a numeric matrix of class `ecodive_ord`, with
#'   one row per sample and one column per axis (named for `target`'s axes),
#'   as `vegan::procrustes()`'s `Yrot`. It is centred, and in `target`'s
#'   frame: `target`'s centred coordinates (scaled to unit sum of squares, if
#'   `symmetric`) are attribute `target`. Its attributes:
#'
#'   * `method`: `'procrustes'`.
#'   * `eig`, `prop`: `NA` for each axis, as Procrustes rotation has no
#'     eigenvalues.
#'   * `ss`: The Procrustes sum of squares, \eqn{m^2}: the sum of squared
#'     distances between `target` and the rotated `ord`. vegan's `ss`.
#'   * `rotation`: The \eqn{k \times k} orthogonal matrix that rotates (and
#'     possibly reflects) the centred `ord` onto `target`. vegan's
#'     `rotation`.
#'   * `translation`: The shift that takes the centroid of the rotated and
#'     dilated `ord` to that of `target`, one value per axis, so that
#'     `dilation * ord %*% rotation + translation` (adding `translation` to
#'     each row of the uncentred `ord`) gives the fitted coordinates. Zero
#'     when `symmetric`, as both are centred first. vegan's `translation`.
#'   * `dilation`: The dilation factor applied to `ord`: 1 when `dilate =
#'     FALSE`. vegan's `scale`.
#'   * `residuals`: Each sample's distance between `target` and the rotated
#'     `ord`, as `residuals()` of a `vegan::procrustes()` result.
#'   * `symmetric`: Whether the symmetric statistic was used.
#'   * `target`: `target`'s coordinates as fitted: centred, scaled to unit
#'     sum of squares if `symmetric`, and padded with zeros to the same
#'     axes. vegan's `X`.
#'
#' @details
#'
#' **Method**
#'
#' Procrustes analysis (Gower 1975; Gower & Dijksterhuis 2004) centres `ord`
#' (\eqn{Y}) and `target` (\eqn{X}), both \eqn{n \times k}, and finds the
#' orthogonal matrix \eqn{A} (a rotation, possibly with a reflection) and
#' the dilation \eqn{c} that minimise
#' \deqn{m^2 = \| X - c Y A \|^2}
#' (the sum of squared differences). With \eqn{X^\top Y = U \Sigma V^\top}
#' (singular value decomposition), \eqn{A = V U^\top}, as in SciPy's
#' `orthogonal_procrustes()`, and
#' \eqn{c = \mathrm{tr}(\Sigma) / \mathrm{tr}(Y^\top Y)}, or 1 when
#' `dilate = FALSE`.
#'
#' The ordinary (non-symmetric) \eqn{m^2} depends on the scales of `ord`
#' and `target`, and on which one is rotated. With `symmetric = TRUE`, both
#' are first scaled to unit sum of squares, so that
#' \eqn{m^2 = 1 - \mathrm{tr}(\Sigma)^2} lies between 0 and 1 and is the
#' same either way round: SciPy's `scipy.spatial.procrustes()` "disparity".
#' [procrustes_test()] reports \eqn{r = \sqrt{1 - m^2}}, this symmetric
#' Procrustes correlation.
#'
#' Rotation can't make up for a different number of axes: compare the same
#' number from each ordination, as the default `axes` does.
#'
#'
#' **Differences from vegan::procrustes()**
#'
#' * The argument order: `procrustes(ord, target)` rotates `ord` onto
#'   `target`, as `vegan::procrustes(target, ord)` (vegan's `procrustes(X,
#'   Y)` rotates `Y` onto `X`).
#' * `dilate` is vegan's `scale`.
#' * By default, the axes both have in common are used; vegan pads the one
#'   with fewer axes with zeros. Give `axes` to do the same.
#' * `ord` and `target` may be `dist` objects, and their samples are
#'   matched by name.
#' * The result is the rotated ordination itself, with the rest as
#'   attributes. There are no `plot()` or `fitted()` methods.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::procrustes()` takes its arguments the other way round (see
#' above). When both packages are attached, the one attached last masks the
#' other: call `ecodive::procrustes()` to be sure.
#'
#' @references
#' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
#' 40(1), 33-51. \doi{10.1007/BF02291478}
#'
#' Gower, J. C., & Dijksterhuis, G. B. (2004). *Procrustes Problems*. Oxford
#' University Press.
#'
#' Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate data
#' sets match? The advantages of a Procrustean superimposition approach over
#' the Mantel test. *Oecologia*, 129(2), 169-178.
#' \doi{10.1007/s004420100720}
#'
#' @seealso [procrustes_test()], [pcoa()], [nmds()], [mantel()]
#'
#' @export
#' @examples
#'     # The same samples under two beta diversity metrics: 30 samples along a gradient that favours some genera at one end
#'     # and others at the other.
#'     set.seed(1)
#'     grad   <- seq(-2, 2, length.out = 30)
#'     counts <- sapply(1:20, function (j) rpois(30, 5 * exp(grad * (j - 10.5) / 5)))
#'     rownames(counts) <- paste0('S', 1:30)
#'
#'     bc <- pcoa(bray(counts), n_axes = 2)
#'     jc <- pcoa(jaccard(counts), n_axes = 2)
#'
#'     fit <- procrustes(jc, bc)
#'     attr(fit, 'ss')
#'     attr(fit, 'residuals')
#'
#'     # Symmetric: m^2 between 0 and 1.
#'     attr(procrustes(jc, bc, symmetric = TRUE), 'ss')
#'
#'     # Distances are ordinated with pcoa() first.
#'     attr(procrustes(jaccard(counts), bray(counts), n_axes = 2), 'ss')
#'
procrustes <- function (
    ord, target, axes = NULL, symmetric = FALSE, dilate = TRUE,
    n_axes = 3L, cpus = n_cpus() ) {

  validate_cpus()

  if (!is.logical(symmetric) || length(symmetric) != 1 || is.na(symmetric))
    stop('`symmetric` must be TRUE or FALSE.', call. = FALSE)
  if (!is.logical(dilate) || length(dilate) != 1 || is.na(dilate))
    stop('`dilate` must be TRUE or FALSE.', call. = FALSE)

  cf <- procrustes_inputs(ord, target, axes, n_axes, cpus)

  procrustes_fit(cf, symmetric, dilate)
}



#' Procrustes Randomization Test (PROTEST)
#'
#' Tests whether two ordinations of the same samples agree better than
#' chance: PROTEST (Jackson 1995; Peres-Neto & Jackson 2001). The statistic
#' is the symmetric Procrustes correlation, \eqn{r = \sqrt{1 - m^2}}, from
#' [procrustes()] with `symmetric = TRUE`, and its p-value comes from
#' permutations of `ord`'s samples. Gives the same statistic as
#' `vegan::protest()`, and the same p-values for the same permutations,
#' computed in parallel C.
#'
#' @inheritParams procrustes
#' @inheritParams mantel
#'
#' @param strata  Permute samples only within these groups (blocks), for
#'        designs such as repeated measures: a vector with one value per
#'        sample, matched to `ord`'s samples by name when it has names,
#'        otherwise in the same order. Default: `NULL` (free permutation).
#'
#' @return A list of class `c('ecodive_procrustes_test', 'ecodive_test')`,
#'   with:
#'
#'   * `statistic`: The symmetric Procrustes correlation, `c(r = ...)`:
#'     \eqn{\sqrt{1 - m^2}}. `vegan::protest()`'s `t0`.
#'   * `p_value`: The permutation p-value: the proportion of permutations,
#'     counting the observed one, with \eqn{r} at least as large.
#'   * `z`: The standardized effect size, `(r - mean(r_perm)) / sd(r_perm)`.
#'   * `n_perm`: The number of permutations.
#'   * `ss`: The symmetric Procrustes sum of squares, \eqn{m^2}.
#'     `vegan::protest()`'s `ss`.
#'   * `n`: The number of samples.
#'   * `n_axes`: The number of axes compared (after any zero-padding).
#'   * `procrustes`: The fit: [procrustes()]'s result, with `symmetric =
#'     TRUE`.
#'   * `permutation`: How samples were permuted: `'free'`, `'within strata'`,
#'     or `'supplied matrix'`.
#'
#'   Attribute `r_perm` is a vector of the statistics of the permutations, as
#'   `vegan::protest()`'s `t`. `as.data.frame()` gives the main results as
#'   one row.
#'
#' @details
#'
#' **Method**
#'
#' Both configurations are centred and scaled to unit sum of squares. The
#' rows of `ord` are then permuted, and the symmetric Procrustes correlation
#' recomputed: the sum of the singular values of the \eqn{k \times k}
#' cross-product of the two. The p-value is
#' \eqn{(1 + \#\{r_{perm} \ge r\}) / (1 + permutations)}, where, as in
#' vegan, permuted statistics within `sqrt(.Machine$double.eps)` of the
#' observed one count as ties.
#'
#' The test is the same whichever of `ord` and `target` is permuted, and,
#' because the statistic is symmetric, whichever is rotated.
#'
#'
#' **Speed and reproducibility**
#'
#' Each permutation is one pass over the samples, \eqn{O(n k^2)}, and a
#' \eqn{k \times k} singular value decomposition; permutations run in
#' parallel. They come from `seed`, not from R's random number generator, so
#' p-values are identical from run to run and for any number of `cpus`, and
#' agree with `vegan::protest()`'s exactly when the same `permutations`
#' matrix is given to both.
#'
#'
#' **Differences from vegan::protest()**
#'
#' * The name: `procrustes_test()`, not `protest()`, so that it sorts beside
#'   [procrustes()].
#' * The argument order, as for [procrustes()]: `procrustes_test(ord,
#'   target)` is `vegan::protest(target, ord)`. The statistic and p-value
#'   are the same either way round.
#' * The default `axes`, as for [procrustes()].
#'
#' @references
#' Jackson, D. A. (1995). PROTEST: a PROcrustean Randomization TEST of
#' community environment concordance. *Ecoscience*, 2(3), 297-303.
#' \doi{10.1080/11956860.1995.11682297}
#'
#' Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate data
#' sets match? The advantages of a Procrustean superimposition approach over
#' the Mantel test. *Oecologia*, 129(2), 169-178.
#' \doi{10.1007/s004420100720}
#'
#' @seealso [procrustes()], [mantel()]
#'
#' @export
#' @examples
#'     # 30 samples along a gradient that favours some genera at one end
#'     # and others at the other.
#'     set.seed(1)
#'     grad   <- seq(-2, 2, length.out = 30)
#'     counts <- sapply(1:20, function (j) rpois(30, 5 * exp(grad * (j - 10.5) / 5)))
#'     rownames(counts) <- paste0('S', 1:30)
#'
#'     bc <- pcoa(bray(counts), n_axes = 2)
#'     jc <- pcoa(jaccard(counts), n_axes = 2)
#'
#'     # Do the Bray-Curtis and Jaccard PCoAs agree?
#'     procrustes_test(jc, bc)
#'
procrustes_test <- function (
    ord, target, axes = NULL, n_axes = 3L, permutations = 999L,
    strata = NULL, seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  cf <- procrustes_inputs(ord, target, axes, n_axes, cpus)
  n  <- nrow(cf$y)

  strata <- strata_codes(per_sample(strata, 'strata', NULL, NULL, cf$labels, n, ref = '`ord`'))
  pm     <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  fit <- procrustes_fit(cf, symmetric = TRUE, dilate = TRUE)
  ss  <- attr(fit, 'ss')
  r   <- sqrt(max(0, 1 - ss))

  # Both centred and scaled to unit sum of squares, one sample per column.
  std <- function (m) {
    m <- sweep(m, 2L, colMeans(m))
    t(m / sqrt(sum(m * m)))
  }

  r_all  <- .Call(
    C_protest, std(cf$x), std(cf$y), pm$perms, strata,
    as.integer(c(pm$n_perm, seed, cpus)) )
  r_perm <- r_all[-1]

  new_ecodive_test(
    class       = 'ecodive_procrustes_test',
    title       = 'Procrustes randomization test (PROTEST)',
    statistic   = c(r = r),
    p_value     = perm_p(r, r_perm, 'greater'),
    z           = perm_z(r, r_perm),
    n_perm      = pm$n_perm,
    ss          = ss,
    n           = as.integer(n),
    n_axes      = ncol(cf$y),
    procrustes  = fit,
    perm        = r_perm,
    permutation = if (!is.null(pm$perms)) 'supplied matrix'
                  else if (!is.null(strata)) 'within strata'
                  else 'free' )
}



#------------------------------------------------------------------------------
# Lines added to print.ecodive_test()'s printout.
#------------------------------------------------------------------------------
#' @export
test_details.ecodive_procrustes_test <- function (x, digits) {
  c(
    paste0('Symmetric Procrustes sum of squares (m^2): ',
           formatC(x$ss, digits = digits, format = 'g')),
    paste0('Samples: ', x$n, ', axes: ', x$n_axes) )
}



#------------------------------------------------------------------------------
# `ord` and `target` as n x k double matrices, with `ord`'s samples' order:
# list(y = ord, x = target, labels, cols (the axes' names, from target where
# it has them)).
#------------------------------------------------------------------------------
procrustes_inputs <- function (ord, target, axes, n_axes, cpus) {

  y <- procrustes_coords(ord,    'ord',    n_axes, cpus)
  x <- procrustes_coords(target, 'target', n_axes, cpus)


  #--------------------------------------------------------
  # Samples: matched by name when both have names.
  #--------------------------------------------------------
  ry <- rownames(y)
  rx <- rownames(x)

  if (!is.null(ry) && !is.null(rx)) {

    if (anyDuplicated(ry)) stop('`ord` has duplicate sample names.', call. = FALSE)
    if (anyDuplicated(rx)) stop('`target` has duplicate sample names.', call. = FALSE)

    miss_x <- setdiff(ry, rx)
    miss_y <- setdiff(rx, ry)
    if (length(miss_x) || length(miss_y)) {
      which_arg <- if (length(miss_x)) c('target', 'ord') else c('ord', 'target')
      miss      <- if (length(miss_x)) miss_x else miss_y
      stop('`', which_arg[1], '` has no samples named ',
           paste0('"', utils::head(miss, 5), '"', collapse = ', '),
           if (length(miss) > 5) ', ...', ', which `', which_arg[2], '` has. ',
           'Both must have the same samples.', call. = FALSE)
    }

    x <- x[match(ry, rx), , drop = FALSE]
  }
  else if (nrow(x) != nrow(y)) {
    stop('`ord` has ', nrow(y), ' samples, but `target` has ', nrow(x), '.', call. = FALSE)
  }

  labels <- if (!is.null(ry)) ry else rx
  n      <- nrow(y)

  if (n < 2)
    stop('`ord` and `target` must have at least 2 samples.', call. = FALSE)


  #--------------------------------------------------------
  # Axes: the first min(ky, kx), or those asked for, with
  # zero columns where one has too few.
  #--------------------------------------------------------
  ky <- ncol(y)
  kx <- ncol(x)

  if (is.null(axes)) {
    axes <- seq_len(min(ky, kx))
  }
  else {
    if (!is.numeric(axes) || length(axes) < 1 || anyNA(axes) || any(axes %% 1 != 0) ||
        any(axes < 1 | axes > max(ky, kx)) || anyDuplicated(axes))
      stop('`axes` must be distinct whole numbers between 1 and ', max(ky, kx),
           ' (the number of columns of `ord` or `target`, whichever is more).',
           call. = FALSE)
    axes <- as.integer(axes)
  }

  pick <- function (m, k) {
    out <- matrix(0, n, length(axes))
    ok  <- axes <= k
    out[, ok] <- m[, axes[ok]]
    out
  }

  cols <- ifelse(axes <= kx, colnames(x)[pmin(axes, kx)], colnames(y)[pmin(axes, ky)])

  list(y = pick(y, ky), x = pick(x, kx), labels = labels, cols = cols)
}



#------------------------------------------------------------------------------
# An ordination, numeric matrix or data frame, or a dist object (ordinated
# with pcoa()), as a double matrix with column names (Dim1, ... if it had
# none) and row names (or NULL).
#------------------------------------------------------------------------------
procrustes_coords <- function (v, arg, n_axes, cpus) {

  if (inherits(v, 'dist')) {
    v <- tryCatch(
      pcoa(v, n_axes = n_axes, cpus = cpus),
      error = function (e) stop(
        gsub('`dm`', paste0('`', arg, '`'), conditionMessage(e), fixed = TRUE),
        call. = FALSE ))
  }

  if (is.data.frame(v)) {
    if (!all(vapply(v, is.numeric, logical(1))))
      stop('`', arg, '` must be a numeric matrix or data frame of coordinates, ',
           'or a dist object.', call. = FALSE)
    v <- as.matrix(v)
  }

  if (!is.matrix(v) || !is.numeric(v) || nrow(v) < 1 || ncol(v) < 1)
    stop('`', arg, '` must be a numeric matrix of coordinates, one row per ',
         'sample, or a dist object.', call. = FALSE)

  dn <- dimnames(v)
  m  <- unclass(v)
  attributes(m) <- list(dim = dim(v))
  storage.mode(m) <- 'double'

  if (!all(is.finite(m)))
    stop('`', arg, '` must not contain NA, NaN, or infinite values.', call. = FALSE)

  cn <- dn[[2]]
  if (is.null(cn)) cn <- rep('', ncol(m))
  bad <- is.na(cn) | cn == ''
  cn[bad] <- paste0('Dim', seq_len(ncol(m)))[bad]

  dimnames(m) <- list(dn[[1]], cn)
  m
}



#------------------------------------------------------------------------------
# The Procrustes fit of cf$y (ord) onto cf$x (target), as procrustes()
# returns it.
#------------------------------------------------------------------------------
procrustes_fit <- function (cf, symmetric, dilate) {

  x <- cf$x
  y <- cf$y

  xmean <- colMeans(x)
  ymean <- colMeans(y)
  x <- sweep(x, 2L, xmean)
  y <- sweep(y, 2L, ymean)

  ssx <- sum(x * x)
  ssy <- sum(y * y)

  if (!(ssx > 0)) stop('`target` has no variation: every sample is at the same point.', call. = FALSE)
  if (!(ssy > 0)) stop('`ord` has no variation: every sample is at the same point.', call. = FALSE)

  if (symmetric) {
    x <- x / sqrt(ssx)
    y <- y / sqrt(ssy)
    ssy <- sum(y * y)
  }

  # SciPy's orthogonal_procrustes(): X'Y = U S V', rotation V U'.
  s   <- svd(crossprod(x, y))
  rot <- s$v %*% t(s$u)
  dil <- if (dilate) sum(s$d) / ssy else 1

  yrot <- dil * y %*% rot
  res  <- x - yrot
  ss   <- sum(res * res)

  trans <- if (symmetric) rep(0, ncol(x)) else drop(xmean - dil * ymean %*% rot)

  cols <- cf$cols
  dimnames(rot) <- list(cols, cols)
  names(trans)  <- cols
  dimnames(x)   <- list(cf$labels, cols)

  resid <- sqrt(rowSums(res * res))
  names(resid) <- cf$labels

  ord_result(
    yrot, 'procrustes', cf$labels, cols = cols,
    ss          = ss,
    rotation    = rot,
    translation = trans,
    dilation    = dil,
    residuals   = resid,
    symmetric   = symmetric,
    target      = x )
}
