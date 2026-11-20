# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Fit Variables or Taxa onto an Ordination
#'
#' Finds the direction in an ordination along which each numeric variable
#' (or taxon) increases fastest, and how strongly it correlates with the
#' ordination (\eqn{r^2}); and, for each factor, the centroid of each of its
#' levels and how well the factor separates the samples. Permutation tests
#' give p-values. Gives the same arrows, centroids, and \eqn{r^2} as
#' `vegan::envfit()`, and the same p-values for the same permutations, but
#' runs in parallel C, so that fitting thousands of taxa with 999
#' permutations each takes seconds. Use it for the arrows and labels of a
#' biplot.
#'
#' @inheritParams documentation
#'
#' @param ord  An ordination: a numeric matrix with one row per sample and
#'        one column per axis, such as from [pcoa()], [pca()], [nmds()],
#'        [dbrda()], [tsne()], or [umap()]. A data frame of numeric columns
#'        also works.
#'
#' @param data  Sample metadata to fit: a data frame with one row per
#'        sample. Numeric columns are fitted as vectors (arrows), and factor,
#'        character, and logical columns as factors (centroids). Rows are
#'        matched to the samples in `ord` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as `ord`. A matrix, or a single vector or factor (named
#'        `group`), is also accepted.
#'
#' @param counts  A feature table to fit, such as taxa abundances: each
#'        feature is fitted as a vector. Accepts the same objects as
#'        [beta_div()]: a numeric matrix (samples in rows, or in columns with
#'        `margin = 2`), a sparse matrix from \pkg{Matrix} or \pkg{slam}, or
#'        a `phyloseq`, `rbiom`, or `SummarizedExperiment` object. Sparse
#'        tables are never made dense. Samples are matched by name to `ord`'s
#'        row names, and samples not in `ord` are ignored; without names,
#'        they must be in the same order as `ord`. Give `data`, `counts`, or
#'        both.
#'
#' @param axes  Which columns of `ord` to fit onto: indices or column names.
#'        The default, `NULL`, uses every column. `vegan::envfit()` uses the
#'        first two by default (its `choices = c(1, 2)`), so give `axes =
#'        1:2` to compare with it.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param strata  Permute samples only within these groups (blocks). With
#'        `data`, the name of one of its columns (which is then not fitted),
#'        or a vector with one value per row of `data`. Without `data`, a
#'        vector with one value per row of `ord`, matched to them by name when
#'        it has names, otherwise in the same order. Default: `NULL` (free
#'        permutation).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @param norm  Normalize `counts` before fitting them, as in
#'        [beta_div()]: `'none'`, `'percent'`, `'chord'`, `'hellinger'`,
#'        `'binary'`, `'clr'`, `'rclr'`, or `'log'`. Default: `'none'`. `margin`, `norm`, and
#'        `pseudocount` apply only to `counts`.
#'
#' @param na.rm  Leave out samples with missing values (`NA`) in `ord` or in
#'        any column of `data`, as `vegan::envfit()` does. `counts` must not
#'        have any.
#'        Default: `FALSE`, which stops with an error if there are any.
#'
#' @return A data frame with one row per numeric variable and per level of
#'   each factor of `data`, in the order of its columns, followed by one row
#'   per feature of `counts`. Its columns are:
#'
#'   * `variable`: The variable's name.
#'   * `level`: The factor level, or `NA` for numeric variables.
#'   * One column per axis of `ord` fitted, named as in `ord` (such as
#'     `PCoA1`, `PCoA2`). For a numeric variable, the direction of its
#'     arrow, as unit-length coordinates (multiply by `sqrt(r2)` for
#'     `vegan::scores(fit, "vectors")`'s arrows). For a factor level, the
#'     mean coordinates of its samples (the centroid).
#'   * `r2`: The squared correlation of a numeric variable with its fitted
#'     values, or for a factor, the proportion of the ordination's variance
#'     that is between its levels (repeated on each level's row).
#'   * `p_value`: The permutation p-value of `r2` (also repeated for
#'     factors).
#'   * `z`: The standardized effect size, `(r2 - mean(r2_perm)) /
#'     sd(r2_perm)`, for that variable's column of `r2_perm`.
#'
#'   Row names are the variables' names, and for factor levels, the
#'   variable's name pasted to the level's, as in `vegan::envfit()`. The
#'   data frame has class `"ecodive_envfit"`, only to mark what made it.
#'
#'   Attributes: `r2_perm`, a matrix of the \eqn{r^2} of each permutation
#'   (rows) for each variable (columns); and `n_perm`, the number of
#'   permutations.
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{X} the \eqn{n \times k} centered coordinates of the axes and
#' \eqn{v} a variable, centered, the multiple regression \eqn{v \approx X
#' \beta} (Dargie 1984; Oksanen et al.) gives the variable's direction:
#' \eqn{\beta}, scaled to unit length. Its strength is
#' \eqn{r^2 = \|\hat{v}\|^2 / \|v\|^2}, the squared correlation between the
#' variable and its fitted values. For a factor, \eqn{r^2} is the sum of
#' squares between its levels' centroids divided by the total sum of squares
#' of \eqn{X}, as in a one-way MANOVA.
#'
#' The p-value is the proportion of permutations of the samples, counting
#' the observed data, with an \eqn{r^2} at least as large: \eqn{(1 +
#' \#\{r^2_{perm} \ge r^2\}) / (1 + permutations)}. As in `vegan::envfit()`,
#' permuted values within `sqrt(.Machine$double.eps)` of the observed one
#' count as ties.
#'
#'
#' **Speed and memory**
#'
#' With \eqn{Q} an orthonormal basis of \eqn{X}'s columns, \eqn{\|\hat{v}\|^2
#' = \|Q^T v\|^2}. A permutation of the samples permutes the rows of
#' \eqn{Q} instead of every variable, so each permutation costs \eqn{k}
#' multiply-adds per non-zero value of the table, with no refitting, and
#' zeros of a sparse table cost nothing. A factor's level sums cost \eqn{k}
#' additions per sample. Permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own: the same as [permanova()]'s for the same `seed` and
#' number of samples. P-values are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `vegan`'s, so
#' p-values agree with it only up to sampling error, unless the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from vegan::envfit()**
#'
#' * All of `ord`'s columns are fitted by default, not the first two.
#' * The result is one data frame, not a list of vectors and factors.
#' * A variable that has the same value in every sample (such as a taxon
#'   found in none of them) has \eqn{r^2 = 0}, an arrow of zeros, and a
#'   p-value and `z` of `NA`, where `vegan` gives a p-value of 1. So does a
#'   factor with one level.
#' * When the axes of `ord` are linearly dependent, the aliased axes get 0
#'   in the arrows, where `vegan` gives `NA` for every arrow.
#' * There are no weights (`w`): no ecodive ordination has them. `vegan`
#'   uses them only for its CA and CCA results.
#' * `strata` permutes freely within strata. For other restricted designs,
#'   pass a `permutations` matrix from `permute::shuffleSet()`.
#' * When there are fewer possible permutations than requested, `vegan`
#'   (via `permute`) runs each of them once instead. `envfit()` always draws
#'   random permutations.
#'
#' @references
#' Dargie, T. C. D. (1984). On the integrated interpretation of indirect
#' site ordinations: a case study using semi-arid vegetation in
#' southeastern Spain. *Vegetatio*, 55, 37-55. \doi{10.1007/BF00039980}
#'
#' Oksanen, J., et al. vegan: Community Ecology Package. `envfit()`.
#' \url{https://CRAN.R-project.org/package=vegan}
#'
#' @seealso [pcoa()], [pca()], [nmds()], [permanova()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- pcoa(dm, n_axes = 2)
#'
#'     # An arrow for each measurement, a centroid for each species.
#'     fit <- envfit(ord, iris)
#'     fit
#'
#'     # A biplot: the fit drawn on the ordination.
#'     plot(ord, 'Species', iris)
#'     plot(fit, add = TRUE)
#'
envfit <- function (
    ord, data = NULL, counts = NULL, margin = 1L, norm = 'none',
    pseudocount = NULL, axes = NULL, na.rm = FALSE, permutations = 999L,
    strata = NULL, seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  if (!isTRUE(na.rm) && !isFALSE(na.rm))
    stop('`na.rm` must be TRUE or FALSE.', call. = FALSE)

  if (is.null(data) && is.null(counts))
    stop('Give `data` (sample metadata), `counts` (a feature table), or both.', call. = FALSE)

  if (is.null(counts) && (!identical(norm, 'none') || !identical(as.numeric(margin), 1) ||
                          !is.null(pseudocount)))
    stop('`margin`, `norm`, and `pseudocount` only apply to `counts`.', call. = FALSE)


  #--------------------------------------------------------
  # The ordination's coordinates, and the metadata aligned
  # with its rows.
  #--------------------------------------------------------
  x    <- envfit_ord(ord, axes)
  n0   <- nrow(x)
  lab0 <- rownames(x)
  keep <- stats::complete.cases(x)

  # `strata`, as every per-sample argument: with `data`, a column name or
  # one value per row of `data`; without, one value per row of `ord`,
  # matched by name when named.
  env <- NULL
  if (!is.null(data)) {
    env    <- envfit_env(data, x, strata, has_counts = !is.null(counts))
    strata <- env$strata
    if (ncol(env$data) > 0) keep <- keep & stats::complete.cases(env$data)
    else env <- NULL
  }
  else {
    strata <- per_sample(strata, 'strata', NULL, NULL, lab0, n0, ref = '`ord`')
  }

  # With na.rm, a sample without a stratum is left out like any other.
  if (!is.null(strata) && na.rm) keep <- keep & !is.na(strata)

  if (!all(keep)) {
    if (!na.rm)
      stop('Missing values (NA) in `ord`', if (!is.null(data)) ' or `data`',
           '; use na.rm = TRUE to leave out those samples.', call. = FALSE)
    x <- x[keep, , drop = FALSE]
  }

  n <- nrow(x)
  if (n < 2)
    stop('At least 2 samples are needed.', call. = FALSE)


  # Strata codes for the samples kept.
  if (!is.null(strata)) strata <- strata_codes(strata[keep])

  pm <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  #--------------------------------------------------------
  # Orthonormal basis of the centered coordinates, with
  # qr()'s rank decision, as vegan's vectorfit().
  #--------------------------------------------------------
  xc <- x - rep(colMeans(x), each = n)
  attributes(xc) <- list(dim = dim(x))
  qx <- qr(xc, tol = 1e-7)
  rk <- qx$rank
  if (rk == 0)
    stop('All samples have the same coordinates in `ord`.', call. = FALSE)
  q  <- qr.Q(qx)[, seq_len(rk), drop = FALSE]

  fit <- list(
    x = x, xc = xc, q = q, qx = qx, pm = pm, strata = strata,
    seed = seed, cpus = cpus )


  #--------------------------------------------------------
  # Metadata rows, then feature rows. Both use the same
  # permutations (same seed or matrix, same samples).
  #--------------------------------------------------------
  parts <- list()

  if (!is.null(env))
    parts[[1]] <- envfit_fit(envfit_vars(env$data[keep, , drop = FALSE]), fit)

  if (!is.null(counts))
    parts[[length(parts) + 1]] <- envfit_fit(
      envfit_table(counts, lab0, n0, keep, margin, norm, pseudocount), fit )

  out     <- do.call(rbind, lapply(parts, `[[`, 'out'))
  rn      <- unlist(lapply(parts, `[[`, 'rn'))
  r2_perm <- do.call(cbind, lapply(parts, `[[`, 'r2_perm'))

  rownames(out) <- make.unique(as.character(rn))
  attr(out, 'r2_perm') <- r2_perm
  attr(out, 'n_perm')  <- pm$n_perm
  class(out) <- c('ecodive_envfit', 'data.frame')

  # The title names the ordination (?ecodive_display).
  on_ <- display_title(ord)
  set_display(out, title = if (is.null(on_)) 'Fitted variables' else paste('Fitted variables on', on_))
}



#------------------------------------------------------------------------------
# Fits one set of variables (from envfit_vars() or envfit_table()). `fit`
# holds the coordinates, their basis, and the permutations. Returns
# list(out = data frame, rn = its row names, r2_perm).
#------------------------------------------------------------------------------
envfit_fit <- function (vars, fit) {

  x  <- fit$x
  n  <- nrow(x)
  qx <- fit$qx
  rk <- qx$rank
  pm <- fit$pm

  facs <- vars$factors
  grp  <- NULL
  nlev <- integer(0)
  if (length(facs)) {
    # Codes by first appearance, not by level, so that sums over levels run
    # in an order that the level order (or character vs factor) cannot change.
    grp  <- vapply(facs, function (f) { i <- as.integer(f); match(i, unique(i)) - 1L }, integer(n))
    grp  <- matrix(grp, nrow = n)
    nlev <- vapply(facs, nlevels, integer(1))
  }


  #--------------------------------------------------------
  # r2 for the identity and each permutation.
  #--------------------------------------------------------
  has_vec <- !is.null(vars$counts)

  res <- .Call(
    C_envfit,
    if (has_vec) fit$q,
    vars$counts,
    as.integer(vars$margin),
    as.integer(vars$norm),
    as.double(vars$pseudocount),
    as.integer(vars$map),
    fit$xc, grp, as.integer(nlev),
    pm$perms, fit$strata,
    as.integer(c(pm$n_perm, fit$seed, fit$cpus)) )

  pvals <- function (r2) perm_p(r2[1, ], r2[-1, , drop = FALSE])
  zvals <- function (r2) perm_z(r2[1, ], r2[-1, , drop = FALSE])

  axis_names <- colnames(x)


  #--------------------------------------------------------
  # Vectors: unit arrows, from the regression coefficients.
  #--------------------------------------------------------
  vec <- NULL
  if (has_vec) {

    m     <- length(vars$vec_names)
    const <- res$const
    r2    <- res$r2v[1, ]
    p     <- pvals(res$r2v)
    z     <- zvals(res$r2v)
    r2[const] <- 0
    p[const]  <- NA_real_
    z[const]  <- NA_real_

    coef <- matrix(0, ncol(x), m)
    coef[qx$pivot[seq_len(rk)], ] <- backsolve(
      qr.R(qx)[seq_len(rk), seq_len(rk), drop = FALSE], res$qtv )
    len  <- pmax(sqrt(colSums(coef^2)), .Machine$double.eps)
    arr  <- t(coef) / len
    arr[const, ] <- 0
    colnames(arr) <- axis_names

    vec <- list(
      names = vars$vec_names,
      rows  = data.frame(
        variable = vars$vec_names, level = NA_character_, arr,
        r2 = r2, p_value = p, z = z, check.names = FALSE, stringsAsFactors = FALSE ),
      perm  = res$r2v[-1, , drop = FALSE] )
    rownames(vec$rows) <- NULL
  }


  #--------------------------------------------------------
  # Factors: centroids of the (uncentered) coordinates.
  #--------------------------------------------------------
  fac <- NULL
  if (length(facs)) {

    r2 <- res$r2f[1, ]
    p  <- pvals(res$r2f)
    z  <- zvals(res$r2f)
    p[nlev < 2] <- NA_real_
    z[nlev < 2] <- NA_real_
    r2[nlev < 2] <- 0

    rows <- lapply(seq_along(facs), function (i) {
      f   <- facs[[i]]
      cen <- rowsum(x, f, reorder = TRUE) / as.vector(table(f))
      cen <- matrix(cen, ncol = ncol(x), dimnames = list(NULL, axis_names))
      data.frame(
        variable = names(facs)[[i]], level = levels(f), cen,
        r2 = r2[[i]], p_value = p[[i]], z = z[[i]], check.names = FALSE, stringsAsFactors = FALSE )
    })

    fac <- list(
      names = names(facs),
      rows  = rows,
      perm  = res$r2f[-1, , drop = FALSE] )
  }


  #--------------------------------------------------------
  # One data frame, in the order of data's columns.
  #--------------------------------------------------------
  type <- vapply(vars$order, `[[`, character(1), 'type')
  idx  <- vapply(vars$order, `[[`, integer(1),   'index')

  if (all(type == 'vector')) {
    out     <- vec$rows
    rn      <- vec$names
    r2_perm <- vec$perm
  }
  else {
    frows <- do.call(rbind, fac$rows)
    nl    <- vapply(fac$rows, nrow, integer(1))
    fbeg  <- cumsum(c(0L, nl))
    nv    <- if (is.null(vec)) 0L else nrow(vec$rows)

    # Each variable's rows in rbind(vector rows, factor level rows).
    rows <- unlist(lapply(seq_along(type), function (i) {
      if (type[[i]] == 'vector') idx[[i]]
      else nv + fbeg[[idx[[i]]]] + seq_len(nl[[idx[[i]]]]) }))

    out <- rbind(if (nv) vec$rows, frows)[rows, , drop = FALSE]
    rn  <- c(if (nv) vec$names, paste0(frows$variable, frows$level))[rows]

    cols    <- ifelse(type == 'vector', idx, nv + idx)
    r2_perm <- cbind(vec$perm, fac$perm)[, cols, drop = FALSE]
  }

  colnames(r2_perm) <- vapply(vars$order, `[[`, character(1), 'name')

  list(out = out, rn = rn, r2_perm = r2_perm)
}



#------------------------------------------------------------------------------
# `ord` -> a numeric n x k matrix of the chosen axes, with row names (or
# NULL) and column names (Dim1, Dim2, ... if it had none).
#------------------------------------------------------------------------------
envfit_ord <- function (ord, axes) {

  if (is.data.frame(ord)) {
    if (!all(vapply(ord, is.numeric, logical(1))))
      stop('`ord` must be a numeric matrix or data frame of coordinates.', call. = FALSE)
    ord <- as.matrix(ord)
  }

  if (!is.matrix(ord) || !is.numeric(ord) || nrow(ord) < 1 || ncol(ord) < 1)
    stop('`ord` must be a numeric matrix of coordinates, one row per sample.', call. = FALSE)

  dn <- dimnames(ord)
  x  <- unclass(ord)
  attributes(x) <- list(dim = dim(ord))
  storage.mode(x) <- 'double'

  rn <- dn[[1]]
  cn <- dn[[2]]
  if (is.null(cn)) cn <- rep('', ncol(x))
  cn[is.na(cn) | cn == ''] <- paste0('Dim', seq_len(ncol(x)))[is.na(cn) | cn == '']

  if (!is.null(axes)) {
    if (is.character(axes)) {
      idx <- match(axes, cn)
      if (anyNA(idx))
        stop('`axes` not in `ord`: ', name_list(axes[is.na(idx)]), call. = FALSE)
      axes <- idx
    }
    if (!is.numeric(axes) || length(axes) < 1 || anyNA(axes) || any(axes %% 1 != 0) ||
        any(axes < 1 | axes > ncol(x)) || anyDuplicated(axes))
      stop('`axes` must be distinct column numbers or names of `ord`, ',
           'between 1 and ', ncol(x), '.', call. = FALSE)
    axes <- as.integer(axes)
  }
  else {
    axes <- seq_len(ncol(x))
  }

  x <- x[, axes, drop = FALSE]
  dimnames(x) <- list(rn, cn[axes])

  if (any(is.infinite(x)))
    stop('`ord` must not contain infinite values.', call. = FALSE)

  x
}



#------------------------------------------------------------------------------
# Metadata (a data frame, matrix, vector, or factor), with its rows aligned to
# ord's rows. Returns list(data, strata): `strata` holds each sample's stratum
# (per_sample()'s rule: a column of `data`, which is then dropped from it, or
# one value per row of `data`), or NULL.
#------------------------------------------------------------------------------
envfit_env <- function (data, x, strata, has_counts = FALSE) {

  n <- nrow(x)

  if (is.null(dim(data))) {
    if (!is.atomic(data))
      stop('`data` must be a data frame, matrix, vector, or factor.', call. = FALSE)
    nms  <- names(data)
    data <- data.frame(group = data, stringsAsFactors = FALSE)
    if (!is.null(nms) && !anyDuplicated(nms)) rownames(data) <- nms
  }
  if (!is.data.frame(data)) {
    if (!is.matrix(data))
      stop('`data` must be a data frame, matrix, vector, or factor.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  dn  <- if (.row_names_info(data) > 0) rownames(data) else NULL
  idx <- match_rows(rownames(x), n, dn, nrow(data), ref = '`ord`')

  # A strata column is a blocking factor, not a variable to fit.
  col <- if (is.character(strata) && length(strata) == 1 && strata %in% names(data)) strata
  strata <- per_sample(strata, 'strata', data, idx, rownames(x), n, ref = '`ord`')
  if (!is.null(col)) data[[col]] <- NULL

  if (ncol(data) == 0 && !has_counts)
    stop('`data` has no variables to fit.', call. = FALSE)

  data <- data[idx, , drop = FALSE]

  ok <- vapply(data, function (v) is.null(dim(v)) &&
    (is.numeric(v) || is.factor(v) || is.character(v) || is.logical(v)), logical(1))
  if (!all(ok))
    stop('Columns of `data` must be numeric, factor, character, or logical: ',
         name_list(names(data)[!ok]), call. = FALSE)

  inf <- names(data)[vapply(data, function (v) is.numeric(v) && any(is.infinite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in `data` column', if (length(inf) > 1) 's', ': ',
         name_list(inf), call. = FALSE)

  list(data = data, strata = strata)
}



#------------------------------------------------------------------------------
# The variables of aligned, complete metadata, for envfit_fit(): numeric
# columns as a dense table of vectors, the others as factors.
#------------------------------------------------------------------------------
envfit_vars <- function (data) {

  num  <- vapply(data, is.numeric, logical(1))
  vnms <- names(data)
  if (is.null(vnms) || any(is.na(vnms) | vnms == '')) vnms <- paste0('V', seq_along(data))

  order <- vector('list', length(data))
  iv <- 0L; jf <- 0L
  for (i in seq_along(data)) {
    if (num[[i]]) { iv <- iv + 1L; order[[i]] <- list(type = 'vector', index = iv, name = vnms[[i]]) }
    else          { jf <- jf + 1L; order[[i]] <- list(type = 'factor', index = jf, name = vnms[[i]]) }
  }

  vars <- list(order = order, factors = NULL, counts = NULL)

  if (any(num)) {
    vars$counts      <- matrix(as.double(unlist(data[num], use.names = FALSE)), nrow = nrow(data))
    vars$margin      <- 1L
    vars$norm        <- 0L
    vars$pseudocount <- 0
    vars$map         <- seq_len(nrow(data)) - 1L
    vars$vec_names   <- vnms[num]
  }

  if (any(!num)) {
    facs <- lapply(data[!num], function (v) droplevels(factor(v, exclude = NULL)))
    names(facs) <- vnms[!num]
    vars$factors <- facs
  }

  vars
}



#------------------------------------------------------------------------------
# A feature table: every feature is a vector. `labels` and `n` are ord's row
# names and rows; `keep` which rows are fitted. Returns the vars for
# envfit_fit().
#------------------------------------------------------------------------------
envfit_table <- function (counts, labels, n, keep, margin, norm, pseudocount) {

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
  idx <- match_rows(labels, n, dn[[1]], ns, ref = '`ord`', what = '`counts`', unit = 'sample')[keep]

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'

  map <- rep(-1L, ns)
  map[idx] <- seq_along(idx) - 1L

  vnms <- dn[[2]]
  if (is.null(vnms)) vnms <- paste0('V', seq_len(dims[[2]]))

  list(
    order       = lapply(seq_along(vnms), function (i) list(type = 'vector', index = i, name = vnms[[i]])),
    factors     = NULL,
    counts      = counts,
    margin      = margin,
    norm        = norm,
    pseudocount = pseudocount,
    map         = map,
    vec_names   = vnms )
}
