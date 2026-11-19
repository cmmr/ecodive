# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Follows bioenv() from scikit-bio, copyright (c) 2013--, scikit-bio
# development team, under the BSD 3-Clause License; see src/bioenv.c and
# inst/COPYRIGHTS for its terms.



#' Best Subset of Environmental Variables (BIOENV)
#'
#' Finds the subset of environmental variables whose distances among the
#' samples best match the community distances: for each subset size, the
#' subset with the highest (rank) correlation between its distances and
#' `dm`. This is BIOENV (Clarke & Ainsworth 1993), as `vegan::bioenv()` and
#' scikit-bio's `bioenv()` (QIIME 2's `diversity bioenv`) compute it, with
#' the same results. Every subset is searched, in parallel C.
#'
#' @inheritParams documentation
#'
#' @param dm  The community distances: a `dist` object, such as from
#'        [beta_div()] or [bray()], or a symmetric numeric matrix of
#'        distances.
#'
#' @param data  The environmental variables: a data frame (or matrix) with
#'        one row per sample. Rows are matched to the samples in `dm` by row
#'        name, and rows for other samples are ignored; without row names,
#'        they must be in the same order as `dm`.
#'
#' @param formula  Which variables to search: a one-sided formula over the
#'        columns of `data`, such as `~ N + P + K + pH` or
#'        `~ log(N) + P`, each of whose terms is one variable. They must be
#'        numeric. The default, `NULL`, searches every numeric column of
#'        `data`, skipping any others (and the column `given` names, if
#'        any).
#'
#' @param method  The correlation coefficient: `'spearman'`, `'pearson'`,
#'        or `'kendall'` (tau-b), as in [stats::cor()]. Default:
#'        `'spearman'`.
#'
#' @param metric  The distances among the samples on a subset of variables:
#'   * `'euclidean'`: Euclidean distances on the variables scaled to unit
#'     variance (and centred).
#'   * `'mahalanobis'`: Mahalanobis distances: Euclidean distances on the
#'     scaled variables, after also making them uncorrelated.
#'   * `'manhattan'`: Manhattan distances on the variables scaled to unit
#'     range (from 0 to 1).
#'   * `'gower'`: Gower distances: the Manhattan distances divided by the
#'     number of variables. Their correlations equal `'manhattan'`'s.
#'
#'   Default: `'euclidean'`, as in scikit-bio and vegan.
#'
#' @param max_vars  The largest subset size to search, as `vegan::bioenv()`'s
#'        `upto`. Default: `NULL`, every size up to the number of
#'        variables.
#'
#' @param given  Distances whose effect is removed from each correlation,
#'        as `vegan::bioenv()`'s `partial`: the partial correlation of `dm`
#'        and each subset's distances, given these. A `dist` object or a
#'        symmetric matrix, matched to `dm`'s samples by name when both
#'        have names; numeric values, one per row of `data`, whose absolute
#'        differences are the distances; or the name of a numeric column of
#'        `data`. Default: `NULL`.
#'
#' @return A data frame of class `c('ecodive_bioenv', 'data.frame')`, with
#'   one row per subset size, as `summary()` of `vegan::bioenv()` and
#'   scikit-bio's `bioenv()`:
#'
#'   * `size`: The number of variables.
#'   * `variables`: The best subset of that size: the variables' names,
#'     joined by `", "`.
#'   * `correlation`: Its correlation with `dm`.
#'
#'   Attributes:
#'
#'   * `best`: The best subset overall (the smallest, of equally good
#'     ones): a character vector of the variables' names.
#'   * `method`, `metric`: As given.
#'   * `partial`: Whether `given` was given.
#'   * `n_subsets`: The number of subsets searched.
#'
#' @details
#'
#' **Method**
#'
#' Each variable is standardized as `metric` says. For each subset of the
#' variables, the distances among the samples on just those variables are
#' correlated with the community distances `dm`, over the \eqn{n(n-1)/2}
#' pairs of samples, as in [mantel()]. For each subset size, the subset
#' with the highest correlation is kept. Subsets of a size are compared in
#' the order of `utils::combn()` (the order of the variables), and the first
#' of equally good ones is kept, as in scikit-bio and vegan.
#'
#' With `given`, each correlation is the partial correlation of `dm` (\eqn{x})
#' and the subset's distances (\eqn{y}) given those distances (\eqn{z}),
#' \deqn{r_{xy.z} = \frac{r_{xy} - r_{xz} r_{yz}}{\sqrt{1 - r_{xz}^2}
#'   \sqrt{1 - r_{yz}^2}}}
#' as in [mantel_partial()]. As vegan says, this is strictly correct only
#' for Pearson.
#'
#' A subset's correlation is undefined (and never best) when its distances
#' are all equal, or, for `'mahalanobis'`, when its variables are linearly
#' dependent. A size with no defined correlation has `NA` values.
#'
#' The best subset is chosen to fit `dm`, so its correlation is biased
#' upwards: a Mantel test of it against `dm` doesn't give a valid p-value.
#'
#'
#' **Speed**
#'
#' There are \eqn{2^p - 1} subsets of \eqn{p} variables (with `max_vars =
#' NULL`), so the search doubles in time with each variable: about a million
#' subsets for 20 variables. `max_vars` limits it. Each subset takes one pass
#' over the pairs of samples (and a sort, for Spearman and Kendall), with no
#' \eqn{n \times n} matrix, and subsets run in parallel. Results are identical
#' for any number of `cpus`.
#'
#'
#' **Differences from vegan::bioenv()**
#'
#' * `dm` must be distances: there's no community table and `index`.
#'   Compute them first, as with [bray()] (which equals `vegdist()`'s
#'   default).
#' * Variables come from `data` through `formula`, or are its numeric
#'   columns, rather than being a separate `env` table. Factors aren't
#'   searched, even with `'gower'`.
#' * `given` may also be one value per sample, or a column name of `data`.
#' * Samples are matched by name between `dm`, `data`, and `given`.
#' * The result is the table that `summary()` of vegan's result gives.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::bioenv()` has the same name. When both packages are attached,
#' the one attached last masks the other; call `ecodive::bioenv()` to be
#' sure.
#'
#' @references
#' Clarke, K. R., & Ainsworth, M. (1993). A method of linking multivariate
#' community structure to environmental variables. *Marine Ecology
#' Progress Series*, 92, 205-219. \doi{10.3354/meps092205}
#'
#' @seealso [mantel()], [mantel_partial()], [envfit()]
#'
#' @export
#' @examples
#'     # 30 samples, with five environmental variables; pH and N shift
#'     # the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     env <- data.frame(
#'       row.names = rownames(counts),
#'       ph    = runif(30, 6, 8),
#'       n     = runif(30, 0, 10),
#'       p     = runif(30, 0, 10),
#'       k     = runif(30, 0, 10),
#'       depth = runif(30, 0, 5) )
#'     counts[, 1:5] <- counts[, 1:5] + round(4 * (env$ph - 6) + env$n / 2)
#'
#'     dm <- bray(counts)
#'
#'     # The best subset of each size.
#'     bioenv(dm, env)
#'
#'     # Only some variables, transformed, with subsets of up to two.
#'     bioenv(dm, env, ~ ph + log(n) + p, max_vars = 2)
#'
#'     # Removing the effect of depth.
#'     bioenv(dm, env, ~ ph + n + p + k, given = 'depth')
#'
bioenv <- function (
    dm, data, formula = NULL, method = 'spearman', metric = 'euclidean',
    max_vars = NULL, given = NULL, cpus = n_cpus() ) {

  validate_cpus()

  method <- match.arg(method, c('spearman', 'pearson', 'kendall'))
  metric <- match.arg(metric, c('euclidean', 'mahalanobis', 'manhattan', 'gower'))

  x <- mantel_dist(dm, 'dm')
  n <- x$n

  if (n < 3)
    stop('`dm` must have at least 3 samples.', call. = FALSE)


  #--------------------------------------------------------
  # `data`'s rows, aligned with dm's samples.
  #--------------------------------------------------------
  if (missing(data) || is.null(data))
    stop('`data` is required: a data.frame of environmental variables.', call. = FALSE)
  if (!is.data.frame(data) && !is.matrix(data))
    stop('`data` must be a data.frame.', call. = FALSE)

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
  idx  <- match_rows(x$labels, n, dn, nrow(data))


  #--------------------------------------------------------
  # `given`: distances, values per row of `data`, or a
  # column name of `data`.
  #--------------------------------------------------------
  given_col <- NULL
  if (!is.null(given)) {
    if (is.character(given) && length(given) == 1) {
      if (is.na(given) || !given %in% names(data))
        stop('`given` is not a column of `data`: ', name_list(given), call. = FALSE)
      given_col <- given
      given     <- data[[given]][idx]
      if (!is.numeric(given))
        stop('`data` column ', name_list(given_col), ' is not numeric.', call. = FALSE)
      if (anyNA(given))
        stop('Missing values (NA) in `data` column ', name_list(given_col), '.', call. = FALSE)
      given <- unname(as.double(given))
    }
    else if (is.atomic(given) && is.null(dim(given)) && !inherits(given, 'dist')) {
      if (!is.numeric(given))
        stop('`given` must be distances, numeric values, or a column name of `data`.', call. = FALSE)
      given <- unname(per_sample(given, 'given', data, idx, x$labels, n))
    }
    given <- mantel_align(mantel_dist(given, 'given'), 'given', x$labels, n)
  }


  #--------------------------------------------------------
  # The candidate variables.
  #--------------------------------------------------------
  env <- bioenv_vars(data[idx, , drop = FALSE], formula, given_col)
  p   <- ncol(env)

  if (is.null(max_vars)) {
    max_vars <- p
  }
  else {
    if (!is.numeric(max_vars) || length(max_vars) != 1 || is.na(max_vars) ||
        max_vars < 1 || max_vars %% 1 != 0)
      stop('`max_vars` must be a whole number, 1 or more.', call. = FALSE)
    max_vars <- as.integer(min(max_vars, p))
  }


  #--------------------------------------------------------
  # Standardized variables: scaled to unit variance, or to
  # unit range.
  #--------------------------------------------------------
  z <- as.matrix(env)
  if (metric %in% c('euclidean', 'mahalanobis')) {
    z <- scale(z)
  }
  else {
    lo <- apply(z, 2L, min)
    z  <- sweep(z, 2L, lo)
    z  <- sweep(z, 2L, apply(z, 2L, max), '/')
  }
  z <- unname(z)
  attributes(z) <- list(dim = dim(z))
  storage.mode(z) <- 'double'

  if (!all(is.finite(z)))
    stop('The variables are too large or small to standardize.', call. = FALSE)


  #--------------------------------------------------------
  # The targets: dm and the given distances.
  #--------------------------------------------------------
  for (v in list(list(x$dm, 'dm'), list(given, 'given')))
    if (!is.null(v[[1]]) && all(v[[1]] == v[[1]][[1]]))
      stop('`', v[[2]], '` has the same distance between every pair of ',
           'samples, so it has no correlation with anything.', call. = FALSE)

  targets <- if (is.null(given)) list(x$dm) else list(x$dm, given)
  ties    <- NULL
  rxz     <- NULL

  if (!is.null(given)) {
    rxz <- mantel_cor(x$dm, list(given), method, perm_matrix(0L, n), NULL, 0L, 1L)[[1]]
    if (abs(rxz) >= 1 - 1e-12)
      stop('`dm` and `given` are perfectly correlated, so there is no ',
           'partial correlation.', call. = FALSE)
  }

  if (method == 'kendall') {
    targets <- lapply(targets, mantel_codes)
    ties    <- vapply(targets, mantel_tied, 0)
  }
  else {
    targets <- lapply(targets, mantel_std, method = method)
  }

  res <- .Call(
    C_bioenv, z, targets, ties, rxz,
    c(match(method, c('pearson', 'spearman', 'kendall')) - 1L,
      match(metric, c('euclidean', 'manhattan', 'gower', 'mahalanobis')) - 1L,
      max_vars, cpus) )


  #--------------------------------------------------------
  # The best subset of each size, and overall.
  #--------------------------------------------------------
  vars <- names(env)
  best <- lapply(seq_len(max_vars), function (k) {
    i <- res[[2]][seq_len(k), k]
    if (anyNA(i)) NULL else vars[i]
  })

  out <- data.frame(
    size        = seq_len(max_vars),
    variables   = vapply(best, function (b) if (is.null(b)) NA_character_ else paste(b, collapse = ', '), ''),
    correlation = res[[1]],
    stringsAsFactors = FALSE )

  top <- if (all(is.na(out$correlation))) NA_integer_ else which.max(out$correlation)

  attr(out, 'best')      <- if (is.na(top)) character(0) else best[[top]]
  attr(out, 'method')    <- method
  attr(out, 'metric')    <- metric
  attr(out, 'partial')   <- !is.null(given)
  attr(out, 'n_subsets') <- sum(choose(p, seq_len(max_vars)))
  class(out) <- c('ecodive_bioenv', 'data.frame')
  out
}



#------------------------------------------------------------------------------
# The candidate variables, from `data` (aligned with the samples): a data
# frame of numeric columns, named. With no formula, every numeric column but
# `given`'s; otherwise one per term of the one-sided formula, each of which
# must be numeric.
#------------------------------------------------------------------------------
bioenv_vars <- function (data, formula, given_col) {

  if (is.null(formula)) {

    keep <- vapply(data, function (v) is.numeric(v) && is.null(dim(v)), logical(1))
    keep[names(data) %in% given_col] <- FALSE
    if (!any(keep))
      stop('`data` has no numeric columns to search.', call. = FALSE)
    env <- data[keep]
  }
  else {

    if (is.character(formula)) {
      f <- tryCatch(str2lang(formula), error = function (e) NULL)
      if (!is.call(f) || !identical(f[[1]], as.name('~')))
        stop('`formula` must be a formula, such as `~ N + P + pH`.', call. = FALSE)
      formula <- eval(f)
      environment(formula) <- parent.frame(2)
    }
    if (!inherits(formula, 'formula'))
      stop('`formula` must be a formula, such as `~ N + P + pH`.', call. = FALSE)

    env <- environment(formula)
    if (length(formula) == 3) formula <- formula[-2]
    environment(formula) <- if (is.null(env)) globalenv() else env

    tt <- tryCatch(
      stats::terms(formula, data = data, specials = 'Condition'),
      error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

    if (length(attr(tt, 'specials')$Condition))
      stop('`formula` cannot have Condition() terms: give the distances ',
           'to partial out as `given`.', call. = FALSE)

    labs <- attr(tt, 'term.labels')
    if (!length(labs))
      stop('`formula` has no variables.', call. = FALSE)
    if (any(attr(tt, 'order') > 1))
      stop('`formula` cannot have interactions: ',
           name_list(labs[attr(tt, 'order') > 1]), call. = FALSE)

    vars    <- all.vars(tt)
    missing <- vars[!vars %in% names(data)]
    missing <- missing[!vapply(missing, function (v) {
      obj <- get0(v, envir = environment(formula), inherits = TRUE)
      !is.null(obj) && !is.function(obj) && length(obj) == 1 }, logical(1))]
    if (length(missing))
      stop('Not found in `data`: ', name_list(missing), call. = FALSE)

    mf <- tryCatch(
      stats::model.frame(tt, data, na.action = stats::na.pass),
      error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

    bad <- labs[!vapply(mf[labs], function (v) is.numeric(v) && is.null(dim(v)), logical(1))]
    if (length(bad))
      stop('Not numeric: ', name_list(bad), '. BIOENV standardizes ',
           'numeric variables; leave factors out of `formula`.', call. = FALSE)

    env <- as.data.frame(lapply(mf[labs], as.double), optional = TRUE)
    names(env) <- labs
  }

  na <- names(env)[vapply(env, anyNA, logical(1))]
  if (length(na))
    stop('Missing values (NA) in ', name_list(na), '.', call. = FALSE)

  inf <- names(env)[vapply(env, function (v) any(is.infinite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in ', name_list(inf), '.', call. = FALSE)

  const <- names(env)[vapply(env, function (v) all(v == v[[1]]), logical(1))]
  if (length(const))
    stop('The same value for every sample, so it cannot be standardized: ',
         name_list(const), '.', call. = FALSE)

  env
}
