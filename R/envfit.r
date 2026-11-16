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
#' @param data  The variables to fit. Either:
#'
#'   * A data frame of sample metadata, with one row per sample. Numeric
#'     columns are fitted as vectors (arrows), and factor, character, and
#'     logical columns as factors (centroids). Rows are matched to the
#'     samples in `ord` by row name, and rows for other samples are ignored;
#'     without row names, they must be in the same order as `ord`. A single
#'     vector or factor is also accepted.
#'   * A table of counts or abundances, with samples in rows (`margin = 1`)
#'     or columns (`margin = 2`), to fit each feature (such as each taxon)
#'     as a vector. Accepts the same objects as [beta_div()]: a numeric
#'     matrix, a sparse matrix from \pkg{Matrix} or \pkg{slam}, or a
#'     `phyloseq`, `rbiom`, or `SummarizedExperiment` object. Sparse tables
#'     are never made dense. Samples are matched by name to `ord`'s row
#'     names, and samples not in `ord` are ignored; without names, they
#'     must be in the same order as `ord`.
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
#' @param strata  Permute samples only within these groups (blocks): the
#'        name of a column of `data` (when `data` is a data frame, which is
#'        then not fitted), or a vector with one value per row of `ord`, in
#'        the same order. Default: `NULL` (free permutation).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @param norm  Normalize a count table in `data` before fitting it, as in
#'        [beta_div()]: `'none'`, `'percent'`, `'chord'`, `'binary'`,
#'        `'clr'`, or `'rclr'`. Default: `'none'`. Only for count tables.
#'
#' @param na.rm  Leave out samples with missing values (`NA`) in `ord` or in
#'        any variable of a data frame `data`, as `vegan::envfit()` does.
#'        Default: `FALSE`, which stops with an error if there are any.
#'
#' @return A data frame with one row per numeric variable (or feature) and
#'   one row per level of each factor, in the order of `data`'s columns. Its
#'   columns are:
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
#'     # 30 samples in three treatment groups. The 'high' group has more
#'     # of the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     dimnames(counts) <- list(paste0('S', 1:30), paste0('G', 1:20))
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       treatment = rep(c('ctrl', 'low', 'high'), each = 10),
#'       ph        = round(rnorm(30, 7, 0.3), 1) )
#'     high <- meta$treatment == 'high'
#'     counts[high, 1:5] <- counts[high, 1:5] + 5
#'
#'     ord <- pcoa(bray(counts), k = 2)
#'
#'     # Metadata: an arrow for pH, a centroid for each treatment.
#'     envfit(ord, meta)
#'
#'     # Each genus's relative abundance as an arrow.
#'     fit <- envfit(ord, counts, norm = 'percent')
#'     head(fit[order(fit$p_value), ])
#'
#'     # A biplot of the genera with p <= 0.05.
#'     top <- fit[fit$p_value <= 0.05, ]
#'     plot(ord, asp = 1)
#'     arrows(0, 0, top$PCoA1 * sqrt(top$r2) * 0.3, top$PCoA2 * sqrt(top$r2) * 0.3)
#'
envfit <- function (
    ord, data, axes = NULL, permutations = 999L, strata = NULL, seed = 0L,
    margin = 1L, norm = 'none', pseudocount = NULL, na.rm = FALSE,
    cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  if (!isTRUE(na.rm) && !isFALSE(na.rm))
    stop('`na.rm` must be TRUE or FALSE.', call. = FALSE)

  if (missing(data) || is.null(data))
    stop('`data` is required: a data frame of variables, or a table of counts.', call. = FALSE)


  #--------------------------------------------------------
  # The ordination's coordinates.
  #--------------------------------------------------------
  x <- envfit_ord(ord, axes)
  n <- nrow(x)


  #--------------------------------------------------------
  # Count table or data frame?
  #--------------------------------------------------------
  is_tbl <- !is.data.frame(data) && !(is.null(dim(data)) && is.atomic(data)) &&
    !(is.matrix(data) && !is.numeric(data) && !is.logical(data))

  if (!is_tbl && (!identical(norm, 'none') || !identical(as.numeric(margin), 1)))
    stop('`norm` and `margin` only apply when `data` is a table of counts.', call. = FALSE)

  res  <- if (is_tbl) envfit_table(data, x, margin, norm, pseudocount, na.rm)
          else        envfit_env(data, x, strata, na.rm)
  vars <- res$vars
  x    <- res$x
  keep <- res$keep

  if (!is.null(res$strata_col)) strata <- res$strata_col
  else if (is.character(strata) && length(strata) == 1 && length(keep) > 1)
    stop('`strata` names no column of `data`.', call. = FALSE)

  n <- nrow(x)
  if (n < 2)
    stop('At least 2 samples are needed.', call. = FALSE)


  #--------------------------------------------------------
  # Strata, parallel to ord's rows (before dropping NAs).
  #--------------------------------------------------------
  if (!is.null(strata)) {
    if (length(strata) != length(keep) || !is.null(dim(strata)))
      stop('`strata` must have one value per row of `ord` (', length(keep),
           ')', if (!is_tbl) ', or name a column of `data`', '.', call. = FALSE)
    strata <- strata[keep]
    if (anyNA(strata)) stop('`strata` must not contain NAs.', call. = FALSE)
    strata <- as.integer(factor(strata)) - 1L
  }

  pm <- perm_matrix(permutations, n) # TODO: use shared helper after Batch C

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


  #--------------------------------------------------------
  # Factors: 0-based level codes.
  #--------------------------------------------------------
  facs <- vars$factors
  grp  <- NULL
  nlev <- integer(0)
  if (length(facs)) {
    grp  <- vapply(facs, function (f) as.integer(f) - 1L, integer(n))
    grp  <- matrix(grp, nrow = n)
    nlev <- vapply(facs, nlevels, integer(1))
  }


  #--------------------------------------------------------
  # r2 for the identity and each permutation.
  #--------------------------------------------------------
  has_vec <- !is.null(vars$counts)

  res <- .Call(
    C_envfit,
    if (has_vec) q,
    vars$counts,
    as.integer(vars$margin),
    as.integer(vars$norm),
    as.double(vars$pseudocount),
    as.integer(vars$map),
    xc, grp, as.integer(nlev),
    pm$perms, strata,
    as.integer(c(pm$n_perm, seed, cpus)) )

  eps <- sqrt(.Machine$double.eps)

  pvals <- function (r2) {
    if (pm$n_perm == 0) return(rep(NA_real_, ncol(r2)))
    obs <- r2[1, ]
    (colSums(r2[-1, , drop = FALSE] >= rep(obs - eps, each = pm$n_perm)) + 1) / (pm$n_perm + 1)
  }

  # (r2 - mean(r2_perm)) / sd(r2_perm), as the other tests' z.
  zvals <- function (r2) {
    if (pm$n_perm < 2) return(rep(NA_real_, ncol(r2)))
    prm <- r2[-1, , drop = FALSE]
    mu  <- colMeans(prm)
    sdv <- sqrt(colSums((prm - rep(mu, each = pm$n_perm))^2) / (pm$n_perm - 1))
    z   <- (r2[1, ] - mu) / sdv
    z[!is.finite(z)] <- NA_real_
    z
  }

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

  rownames(out) <- make.unique(as.character(rn))
  dimnames(r2_perm) <- list(NULL, vapply(vars$order, `[[`, character(1), 'name'))

  attr(out, 'r2_perm') <- r2_perm
  attr(out, 'n_perm')  <- pm$n_perm
  class(out) <- c('ecodive_envfit', 'data.frame')
  out
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
# For each of ord's rows, the index of its row in `data` (by name, or by
# position). Rows of `data` for other samples are left out.
# TODO: use shared metadata matching after Batch C.
#------------------------------------------------------------------------------
envfit_match <- function (labels, n, data_names, n_data, what = '`data`') {

  named <- !is.null(labels) && !is.null(data_names)
  if (named && all(grepl("^[0-9]+$", data_names)))
    named <- all(labels %in% data_names)

  if (named) {

    if (anyDuplicated(labels))
      stop('`ord` has duplicate row names, so ', what, ' cannot be matched to them.', call. = FALSE)
    if (anyDuplicated(data_names))
      stop(what, ' has duplicate sample names.', call. = FALSE)

    idx <- match(labels, data_names)

    if (anyNA(idx)) {
      miss <- labels[is.na(idx)]
      stop(
        sprintf('%d sample%s in `ord` %s no match in %s: %s',
          length(miss), if (length(miss) == 1) '' else 's',
          if (length(miss) == 1) 'has' else 'have', what, name_list(miss) ),
        call. = FALSE )
    }
    return (idx)
  }

  if (n_data != n)
    stop(
      sprintf('%s has %d samples, but `ord` has %d rows.', what, n_data, n),
      if (!is.null(labels)) '\nGive them names to match them to the rows of `ord`.',
      call. = FALSE )

  seq_len(n)
}



#------------------------------------------------------------------------------
# A data frame (or vector) of variables. Returns list(vars, x, keep,
# strata_col): `vars` as for C_envfit, `x` the coordinates of the samples
# kept, `keep` which of ord's rows are kept (na.rm).
#------------------------------------------------------------------------------
envfit_env <- function (data, x, strata, na.rm) {

  n <- nrow(x)

  if (is.null(dim(data))) {
    nms  <- names(data)
    data <- data.frame(group = data, stringsAsFactors = FALSE)
    if (!is.null(nms) && !anyDuplicated(nms)) rownames(data) <- nms
  }
  if (!is.data.frame(data))
    data <- as.data.frame(data, stringsAsFactors = FALSE)

  # `strata` as a column name: its values, and no longer a variable.
  strata_col <- NULL
  if (is.character(strata) && length(strata) == 1 && strata %in% names(data)) {
    strata_col <- data[[strata]]
    data[[strata]] <- NULL
  }

  dn  <- if (.row_names_info(data) > 0) rownames(data) else NULL
  idx <- envfit_match(rownames(x), n, dn, nrow(data))
  if (!is.null(strata_col)) strata_col <- strata_col[idx]

  if (ncol(data) == 0)
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

  keep <- stats::complete.cases(x) & stats::complete.cases(data)
  if (!all(keep)) {
    if (!na.rm)
      stop('Missing values (NA) in `ord` or `data`; use na.rm = TRUE to leave out those samples.',
           call. = FALSE)
    x    <- x[keep, , drop = FALSE]
    data <- data[keep, , drop = FALSE]
  }

  num  <- vapply(data, is.numeric, logical(1))
  vnms <- names(data)
  if (is.null(vnms) || any(vnms == '')) vnms <- paste0('V', seq_along(data))

  order <- vector('list', length(data))
  iv <- 0L; jf <- 0L
  for (i in seq_along(data)) {
    if (num[[i]]) { iv <- iv + 1L; order[[i]] <- list(type = 'vector', index = iv, name = vnms[[i]]) }
    else          { jf <- jf + 1L; order[[i]] <- list(type = 'factor', index = jf, name = vnms[[i]]) }
  }

  vars <- list(order = order, factors = NULL, counts = NULL)

  if (any(num)) {
    mtx <- matrix(as.double(unlist(data[num], use.names = FALSE)), nrow = nrow(data))
    vars$counts    <- mtx
    vars$margin    <- 1L
    vars$norm      <- 0L
    vars$pseudocount <- 0
    vars$map       <- seq_len(nrow(data)) - 1L
    vars$vec_names <- vnms[num]
  }

  if (any(!num)) {
    facs <- lapply(data[!num], function (v) droplevels(factor(v)))
    names(facs) <- vnms[!num]
    vars$factors <- facs
  }

  list(vars = vars, x = x, keep = keep, strata_col = strata_col)
}



#------------------------------------------------------------------------------
# A table of counts: every feature is a vector. Returns list(vars, x, keep).
#------------------------------------------------------------------------------
envfit_table <- function (counts, x, margin, norm, pseudocount, na.rm) {

  n <- nrow(x)
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
  idx <- envfit_match(rownames(x), n, dn[[1]], ns, '`data`')

  keep <- stats::complete.cases(x)
  if (!all(keep)) {
    if (!na.rm)
      stop('Missing values (NA) in `ord`; use na.rm = TRUE to leave out those samples.',
           call. = FALSE)
    x   <- x[keep, , drop = FALSE]
    idx <- idx[keep]
  }

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'

  map <- rep(-1L, ns)
  map[idx] <- seq_along(idx) - 1L

  vnms <- dn[[2]]
  if (is.null(vnms)) vnms <- paste0('V', seq_len(dims[[2]]))

  order <- lapply(seq_along(vnms), function (i) list(type = 'vector', index = i, name = vnms[[i]]))

  vars <- list(
    order       = order,
    factors     = NULL,
    counts      = counts,
    margin      = margin,
    norm        = norm,
    pseudocount = pseudocount,
    map         = map,
    vec_names   = vnms )

  list(vars = vars, x = x, keep = keep)
}
