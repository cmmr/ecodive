# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Mantel Test
#'
#' Tests whether two distance matrices are correlated: for instance, whether
#' samples with more different communities are also further apart, or more
#' different in pH. Gives the same statistic as `vegan::mantel()` (Pearson,
#' Spearman, or Kendall correlation), with permutation p-values. Computes
#' them in parallel C, straight from the distances, without the
#' \eqn{n \times n} matrices that `mantel()` builds for every permutation.
#' [mantel_partial()] removes the effect of a third distance matrix.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values. Its samples are the ones permuted.
#'
#' @param y  The distances to correlate with `dm`: a `dist` object or a
#'        symmetric numeric matrix, like `dm`, matched to `dm`'s samples by
#'        name when both have names, otherwise in the same order. Or numeric
#'        values, such as age or pH, whose absolute differences (as
#'        `stats::dist(y)`) are the distances: with `data`, the name of one
#'        of its numeric columns or a vector with one value per row of
#'        `data`; without `data`, a vector with one value per sample, matched
#'        to the samples by name when it has names, otherwise in the same
#'        order.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, whose columns `y` and `strata` can name. Rows are
#'        matched to the samples in `dm` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as `dm`. Default: `NULL`.
#'
#' @param method  The correlation coefficient: `'pearson'`, `'spearman'`,
#'        or `'kendall'` (tau-b), as in [stats::cor()]. Default:
#'        `'pearson'`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param alternative  Which correlations count as more extreme than the
#'        observed one: `'greater'` (more positive, as in `vegan::mantel()`),
#'        `'less'` (more negative), or `'two.sided'` (larger in absolute
#'        value). Default: `'greater'`.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A list of class `c('ecodive_mantel', 'ecodive_test')`, with:
#'
#'   * `statistic`: The Mantel statistic, `c(r = ...)`: the correlation of
#'     `dm` and `y`. `vegan::mantel()`'s `statistic`.
#'   * `p_value`: The permutation p-value.
#'   * `z`: The standardized effect size, `(r - mean(r_perm)) / sd(r_perm)`.
#'   * `n_perm`: The number of permutations.
#'   * `method`: The correlation coefficient: `'pearson'`, `'spearman'`, or
#'     `'kendall'`.
#'   * `alternative`: `'greater'`, `'less'`, or `'two.sided'`.
#'   * `partial`: `FALSE` (it is `TRUE` from [mantel_partial()]).
#'   * `n`: The number of samples.
#'   * `permutation`: How samples were permuted: `'free'`, `'within strata'`,
#'     or `'supplied matrix'`.
#'
#'   Attribute `r_perm` is a vector of the statistics of the permutations, as
#'   `vegan::mantel()`'s `perm`. `as.data.frame()` gives the main results as
#'   one row.
#'
#' @details
#'
#' **Method**
#'
#' The Mantel test (Mantel 1967; Legendre & Legendre 2012) correlates the
#' \eqn{n(n-1)/2} distances of `dm` with the matching distances of `y`. The
#' distances aren't independent of one another, so the p-value comes from
#' permutations of the samples: the rows and columns of `dm` are shuffled
#' together, `y` is held fixed, and the correlation is recomputed. The
#' p-value is the proportion of permutations, counting the observed one,
#' with a correlation at least as extreme:
#' \eqn{(1 + \#\{r_{perm} \ge r\}) / (1 + permutations)} for
#' `alternative = 'greater'`. As in `vegan::mantel()`, permuted statistics
#' within `sqrt(.Machine$double.eps)` of the observed one count as ties.
#'
#'
#' **Speed and memory**
#'
#' `dm` and `y` are centred and scaled once (after ranking them, for
#' Spearman), so each permuted correlation is one pass over the distances,
#' reading `dm` through the permutation: \eqn{O(n^2)} per permutation, with
#' no \eqn{n \times n} matrix. Kendall's tau-b is counted with Knight's (1966)
#' merge sort, in \eqn{O(n^2 \log n)} per permutation, where
#' `stats::cor()` takes \eqn{O(n^4)}. Permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. P-values are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `vegan::mantel()`'s,
#' so p-values agree with it only up to sampling error, unless the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from vegan::mantel()**
#'
#' * `y` may also be a numeric vector, one value per sample, or the name of
#'   a column of `data`.
#' * Samples are matched by name between `dm`, `y`, `data`, and `strata`.
#' * `alternative` adds lower-tail and two-sided tests.
#' * There's no `na.rm`: distances must not be `NA`.
#' * When there are fewer possible permutations than requested,
#'   `vegan::mantel()` (via `permute`) runs each of them once instead.
#'   `mantel()` always draws random permutations.
#'
#' @references
#' Knight, W. R. (1966). A computer method for calculating Kendall's tau
#' with ungrouped data. *Journal of the American Statistical Association*,
#' 61(314), 436-439. \doi{10.1080/01621459.1966.10480879}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
#' section 10.5. Elsevier.
#'
#' Mantel, N. (1967). The detection of disease clustering and a generalized
#' regression approach. *Cancer Research*, 27(2), 209-220.
#'
#' @seealso [mantel_partial()], [permanova()], [dist_test()], [beta_div()]
#'
#' @export
#' @examples
#'     # 30 samples along a pH gradient, which shifts the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     ph <- setNames(round(runif(30, 6, 8), 1), rownames(counts))
#'     counts[, 1:5] <- counts[, 1:5] + round(4 * (ph - 6))
#'
#'     dm <- bray(counts)
#'
#'     # Are samples that differ more in pH also more different in
#'     # community composition?
#'     res <- mantel(dm, ph)
#'     res
#'     plot(res)
#'
#'     # With rank correlation.
#'     mantel(dm, ph, method = 'spearman')
#'
#'     # The same as a column of sample metadata.
#'     meta <- data.frame(row.names = rownames(counts), ph = ph)
#'     mantel(dm, 'ph', data = meta, method = 'spearman')
#'
#'     # Two distance matrices.
#'     mantel(dm, dist(counts[, 6:20]))
#'
mantel <- function (
    dm, y, data = NULL, method = 'pearson',
    permutations = 999L, strata = NULL, alternative = 'greater',
    seed = 0L, cpus = n_cpus() ) {

  mantel_core(
    dm, y, data = data, given = NULL, method = method,
    permutations = permutations, strata = strata, alternative = alternative,
    seed = seed, cpus = cpus )
}



#' Partial Mantel Test
#'
#' Tests whether two distance matrices are correlated after removing the
#' effect of a third: for instance, whether samples with more different
#' communities are also more different in pH, beyond what geographic
#' distance between them explains. Gives the same statistic as
#' `vegan::mantel.partial()` (Pearson, Spearman, or Kendall partial
#' correlation), with permutation p-values, computed in parallel C straight
#' from the distances. [mantel()] is the simple Mantel test.
#'
#' @inheritParams mantel
#'
#' @param given  The distances whose effect is removed from the correlation
#'        of `dm` and `y`, as `vegan::mantel.partial()`'s `zdis`: such as
#'        geographic distances between the samples. Given in any of the
#'        forms `y` takes, including as a column name of `data`. Required.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, whose columns `y`, `given`, and `strata` can name. Rows are
#'        matched to the samples in `dm` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the same
#'        order as `dm`. Default: `NULL`.
#'
#' @return A list of class `c('ecodive_mantel', 'ecodive_test')`, as from
#'   [mantel()], with:
#'
#'   * `statistic`: The partial Mantel statistic, `c(r = ...)`: the partial
#'     correlation of `dm` and `y` given `given`.
#'     `vegan::mantel.partial()`'s `statistic`.
#'   * `p_value`: The permutation p-value.
#'   * `z`: The standardized effect size, `(r - mean(r_perm)) / sd(r_perm)`.
#'   * `n_perm`: The number of permutations.
#'   * `method`: The correlation coefficient: `'pearson'`, `'spearman'`, or
#'     `'kendall'`.
#'   * `alternative`: `'greater'`, `'less'`, or `'two.sided'`.
#'   * `partial`: `TRUE`.
#'   * `n`: The number of samples.
#'   * `permutation`: How samples were permuted: `'free'`, `'within strata'`,
#'     or `'supplied matrix'`.
#'
#'   Attribute `r_perm` is a vector of the statistics of the permutations, as
#'   `vegan::mantel.partial()`'s `perm`. `as.data.frame()` gives the main
#'   results as one row, with `test = 'mantel_partial'`.
#'
#' @details
#'
#' **Method**
#'
#' The partial Mantel test (Smouse et al. 1986; Legendre & Legendre 2012)
#' correlates the \eqn{n(n-1)/2} distances of `dm` (\eqn{x}) and `y` after
#' removing the effect of `given` (\eqn{z}), as the partial correlation
#' \deqn{r_{xy.z} = \frac{r_{xy} - r_{xz} r_{yz}}{\sqrt{1 - r_{xz}^2}
#'   \sqrt{1 - r_{yz}^2}}}
#' where each \eqn{r} is the Pearson, Spearman, or Kendall correlation of
#' two sets of distances. As in `vegan::mantel.partial()`, the p-value comes
#' from permutations of the samples of `dm`, with `y` and `given` held fixed:
#' \eqn{r_{xy}} and \eqn{r_{xz}} are recomputed for each permutation, and
#' \eqn{r_{yz}} is computed once. The p-value is the proportion of
#' permutations, counting the observed one, with a partial correlation at
#' least as extreme, as in [mantel()].
#'
#' Its p-values can be too small when the distances are spatially
#' autocorrelated (Legendre & Fortin 2010; Guillot & Rousset 2013).
#'
#' When `y` and `given`, or `dm` and `given`, are perfectly correlated, there
#' is no partial correlation: the statistic and p-value are `NA`, with a
#' warning.
#'
#'
#' **Speed, memory, and reproducibility**
#'
#' As for [mantel()]: each permutation is one pass over the distances (two
#' correlations, \eqn{r_{xy}} and \eqn{r_{xz}}, in the same pass), with no
#' \eqn{n \times n} matrix, and permutations run in parallel. They come from
#' `seed`, so p-values are identical from run to run and for any number of
#' `cpus`, and agree with `vegan::mantel.partial()`'s exactly when the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from vegan::mantel.partial()**
#'
#' * `y` and `given` may also be numeric vectors, one value per sample, or
#'   names of columns of `data`.
#' * Samples are matched by name between `dm`, `y`, `given`, `data`, and
#'   `strata`.
#' * `alternative` adds lower-tail and two-sided tests.
#' * There's no `na.rm`: distances must not be `NA`.
#'
#' @references
#' Guillot, G., & Rousset, F. (2013). Dismantling the Mantel tests.
#' *Methods in Ecology and Evolution*, 4(4), 336-344.
#' \doi{10.1111/2041-210x.12018}
#'
#' Legendre, P., & Fortin, M.-J. (2010). Comparison of the Mantel test and
#' alternative approaches for detecting complex multivariate relationships
#' in the spatial analysis of genetic data. *Molecular Ecology Resources*,
#' 10(5), 831-844. \doi{10.1111/j.1755-0998.2010.02866.x}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
#' section 10.5. Elsevier.
#'
#' Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression
#' and correlation extensions of the Mantel test of matrix correspondence.
#' *Systematic Zoology*, 35(4), 627-632. \doi{10.2307/2413122}
#'
#' @seealso [mantel()], [dist_test()], [beta_div()]
#'
#' @export
#' @examples
#'     # 30 samples along a pH gradient, which shifts the first five genera,
#'     # at random locations.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       ph  = round(runif(30, 6, 8), 1),
#'       lat = runif(30, 30, 35),
#'       lon = runif(30, -100, -95) )
#'     counts[, 1:5] <- counts[, 1:5] + round(4 * (meta$ph - 6))
#'
#'     dm  <- bray(counts)
#'     geo <- dist(meta[, c('lat', 'lon')])
#'
#'     # Are samples that differ more in pH also more different in
#'     # community composition, beyond the effect of geographic distance?
#'     mantel_partial(dm, 'ph', given = geo, data = meta)
#'
#'     # With rank correlation.
#'     mantel_partial(dm, 'ph', given = geo, data = meta, method = 'spearman')
#'
#'     # The same, through dist_test().
#'     dist_test(dm, 'mantel_partial', 'ph', meta, given = geo, method = 'spearman')
#'
mantel_partial <- function (
    dm, y, given, data = NULL, method = 'pearson',
    permutations = 999L, strata = NULL, alternative = 'greater',
    seed = 0L, cpus = n_cpus() ) {

  if (missing(given) || is.null(given))
    stop(
      '`given` is required: the distances whose effect is removed, as a dist ',
      'object, a matrix, one value per sample, or a column name of `data`. ',
      'For a simple Mantel test, use mantel().', call. = FALSE )

  mantel_core(
    dm, y, data = data, given = given, method = method,
    permutations = permutations, strata = strata, alternative = alternative,
    seed = seed, cpus = cpus )
}



#------------------------------------------------------------------------------
# mantel() (given = NULL) and mantel_partial().
#------------------------------------------------------------------------------
mantel_core <- function (
    dm, y, data, given, method, permutations, strata, alternative, seed,
    cpus ) {

  validate_seed()
  validate_cpus()

  method      <- match.arg(method, c('pearson', 'spearman', 'kendall'))
  alternative <- match.arg(alternative, c('greater', 'less', 'two.sided'))

  x <- mantel_dist(dm, 'dm')
  n <- x$n

  if (n < 3)
    stop('`dm` must have at least 3 samples.', call. = FALSE)

  if (missing(y) || is.null(y))
    stop('`y` is required: distances, one value per sample, or a column name of `data`.', call. = FALSE)


  #--------------------------------------------------------
  # With `data`, `y`, `given` and `strata` may name its
  # columns; its rows are matched to dm's samples.
  #--------------------------------------------------------
  if (is.null(data)) {
    for (arg in c('y', 'given')) {
      v <- get(arg)
      if (is.character(v) && length(v) == 1)
        stop('`', arg, '` names a column, but `data` is not given.', call. = FALSE)
    }
    strata <- strata_codes(per_sample(strata, 'strata', NULL, NULL, x$labels, n))
  }
  else {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    des <- model_design(
      data, formula = ~ 1, labels = x$labels, n = n, strata = strata,
      allow_empty = TRUE )
    strata <- des$strata

    # A vector has one value per row of `data` (as `groups` and `strata`).
    df <- as.data.frame(data, stringsAsFactors = FALSE)
    for (arg in c('y', 'given')) {
      v <- get(arg)
      if (is.atomic(v) && is.null(dim(v)) && !inherits(v, 'dist') &&
          !(is.character(v) && length(v) == 1))
        assign(arg, unname(per_sample(v, arg, df, des$idx, x$labels, n)))
    }

    y      <- mantel_column(y,     'y',     des$data)
    given  <- mantel_column(given, 'given', des$data)
  }

  y <- mantel_align(mantel_dist(y, 'y'), 'y', x$labels, n)
  g <- if (is.null(given)) NULL else mantel_align(mantel_dist(given, 'given'), 'given', x$labels, n)

  pm <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  x <- x$dm

  result <- function (r_all) mantel_result(
    r_all, method, alternative, partial = !is.null(g), n = n,
    permutation = if (!is.null(pm$perms)) 'supplied matrix'
                  else if (!is.null(strata)) 'within strata'
                  else 'free' )

  undefined <- rep(NA_real_, pm$n_perm + 1L)


  #--------------------------------------------------------
  # Constant distances have no correlation.
  #--------------------------------------------------------
  for (v in list(list(x, 'dm'), list(y, 'y'), list(g, 'given')))
    if (!is.null(v[[1]]) && all(v[[1]] == v[[1]][[1]])) {
      warning('`', v[[2]], '` has the same distance between every pair of ',
              'samples, so it has no correlation with anything.', call. = FALSE)
      return (result(undefined))
    }


  #--------------------------------------------------------
  # Correlations of permuted x with y (and the given
  # distances, g), and of y with g.
  #--------------------------------------------------------
  targets <- if (is.null(g)) list(y) else list(y, g)
  r_all   <- mantel_cor(x, targets, method, pm, strata, seed, cpus)

  if (is.null(g)) {
    r_all <- r_all[1, ]
  }
  else {
    ryg <- mantel_cor(y, list(g), method, perm_matrix(0L, n), NULL, seed, 1L)[[1]]
    if (abs(ryg) >= 1 - 1e-12) {
      warning('`y` and `given` are perfectly correlated, so there is no ',
              'partial correlation.', call. = FALSE)
      return (result(undefined))
    }
    rxy   <- r_all[1, ]
    rxg   <- r_all[2, ]
    if (abs(rxg[[1]]) >= 1 - 1e-12) {
      warning('`dm` and `given` are perfectly correlated, so there is no ',
              'partial correlation.', call. = FALSE)
      return (result(undefined))
    }
    r_all <- (rxy - rxg * ryg) / sqrt(1 - rxg * rxg) / sqrt(1 - ryg * ryg)  # as vegan
    r_all[!is.finite(r_all) | abs(rxg) >= 1 - 1e-12] <- NA_real_
  }

  result(r_all)
}



#------------------------------------------------------------------------------
# Lines added to print.ecodive_test()'s printout.
#------------------------------------------------------------------------------
#' @export
test_details.ecodive_mantel <- function (x, digits) {
  c(
    paste0('Correlation: ', switch(x$method,
      pearson  = "Pearson's product-moment correlation",
      spearman = "Spearman's rank correlation rho",
      kendall  = "Kendall's rank correlation tau-b" )),
    paste0('Alternative: ', switch(x$alternative,
      greater   = 'greater (r is more positive than by chance)',
      less      = 'less (r is more negative than by chance)',
      two.sided = 'two-sided (|r| is larger than by chance)' )),
    paste0('Samples: ', x$n) )
}



#------------------------------------------------------------------------------
# One row, as for every test, named for the function that ran it (and its
# dist_test() id): 'mantel' or 'mantel_partial'.
#------------------------------------------------------------------------------
#' @export
as.data.frame.ecodive_mantel <- function (x, row.names = NULL, optional = FALSE, ...) {
  df <- NextMethod()
  if (isTRUE(x$partial)) df$test <- 'mantel_partial'
  df
}



#------------------------------------------------------------------------------
# A dist object, symmetric matrix, or (except for dm) a numeric vector of
# values per sample, as dist_vector()'s list(dm, n, labels).
#------------------------------------------------------------------------------
mantel_dist <- function (v, arg) {

  if (arg != 'dm' && is.numeric(v) && is.null(dim(v)) && !inherits(v, 'dist')) {

    if (!all(is.finite(v)))
      stop('`', arg, '` must not contain NA, NaN, or infinite values.', call. = FALSE)

    # |v[i] - v[j]|, in stats::dist(v)'s order (which vegan users pass):
    # the same values, without its square and root (which overflow or
    # underflow for very large or small values).
    nms <- names(v)
    v   <- unname(as.double(v))
    m <- length(v)
    j <- rep.int(seq_len(m - 1L), (m - 1L):1)
    d <- abs(v[j + sequence((m - 1L):1)] - v[j])
    if (!all(is.finite(d)))
      stop('The differences between values of `', arg, '` are too large.', call. = FALSE)
    return (list(dm = d, n = m, labels = nms))
  }

  # dist_vector()'s errors name `dm`.
  tryCatch(
    dist_vector(v),
    error = function (e) {
      msg <- gsub('`dm`', paste0('`', arg, '`'), conditionMessage(e), fixed = TRUE)
      if (arg != 'dm')
        msg <- sub(
          'a symmetric numeric matrix.',
          'a symmetric numeric matrix, or a numeric vector with one value per sample.',
          msg, fixed = TRUE )
      stop(msg, call. = FALSE)
    })
}



#------------------------------------------------------------------------------
# v's packed distances in the order of dm's samples, matched by label when
# both have labels.
#------------------------------------------------------------------------------
mantel_align <- function (v, arg, labels, n) {

  if (v$n != n)
    stop('`', arg, '` has ', v$n, ' samples, but `dm` has ', n, '.', call. = FALSE)

  if (is.null(labels) || is.null(v$labels) || identical(labels, v$labels))
    return (v$dm)

  if (anyDuplicated(labels))
    stop('`dm` has duplicate sample names.', call. = FALSE)
  if (anyDuplicated(v$labels))
    stop('`', arg, '` has duplicate sample names.', call. = FALSE)

  idx <- match(labels, v$labels)
  if (anyNA(idx))
    stop('`', arg, '` has no samples named ',
         paste0('"', utils::head(labels[is.na(idx)], 5), '"', collapse = ', '),
         if (sum(is.na(idx)) > 5) ', ...', '.', call. = FALSE)

  # Pair (i, j), i > j, of dm is pair (idx[i], idx[j]) of v.
  j  <- rep.int(seq_len(n - 1L), (n - 1L):1)
  i  <- j + sequence((n - 1L):1)
  a  <- idx[i]
  b  <- idx[j]
  hi <- as.double(pmax(a, b))
  lo <- as.double(pmin(a, b))

  v$dm[(lo - 1) * n - (lo - 1) * lo / 2 + hi - lo]
}



#------------------------------------------------------------------------------
# `v` (`y` or `given`) as given, or, if it names a column of `data` (already
# aligned with dm's samples), that column's values in dm's sample order.
#------------------------------------------------------------------------------
mantel_column <- function (v, arg, data) {

  if (!is.character(v) || length(v) != 1) return (v)

  if (is.na(v) || !v %in% colnames(data))
    stop('`', arg, '` is not a column of `data`: ', name_list(v), call. = FALSE)

  col <- data[[v]]

  if (!is.numeric(col))
    stop(
      '`data` column ', name_list(v), ' is not numeric. A Mantel test correlates ',
      'distances; to test a categorical variable, use permanova(), mrpp(), or ',
      'anosim().', call. = FALSE )

  if (anyNA(col))
    stop('Missing values (NA) in `data` column ', name_list(v), '.', call. = FALSE)

  unname(as.double(col))
}



#------------------------------------------------------------------------------
# Correlations of x, permuted, with each of `targets` (packed distances, in
# the same sample order): an m x (n_perm + 1) matrix, column 1 observed.
#------------------------------------------------------------------------------
mantel_cor <- function (x, targets, method, pm, strata, seed, cpus) {

  n <- (1 + sqrt(1 + 8 * length(x))) / 2

  if (method == 'kendall') {

    # x's pairs in order, with its tie groups; the targets' codes.
    xc  <- mantel_codes(x)
    cnt <- tabulate(xc + 1L)
    tcs <- lapply(targets, mantel_codes)

    args <- list(
      order(xc, method = 'radix') - 1L,
      tcs,
      as.integer(c(0, cumsum(cnt))),
      c(mantel_tied(xc), vapply(tcs, mantel_tied, 0)) )
    meth <- 1L
  }
  else {

    args <- list(
      mantel_std(x, method),
      lapply(targets, mantel_std, method = method), NULL, NULL )
    meth <- 0L
  }

  .Call(
    C_mantel, args[[1]], args[[2]], args[[3]], args[[4]], pm$perms, strata,
    as.integer(c(n, meth, pm$n_perm, seed, cpus)) )
}



#------------------------------------------------------------------------------
# Distances (or any values) centred and scaled to unit length, after ranking
# for Spearman: their Pearson correlation is then a dot product. Also
# bioenv()'s.
#------------------------------------------------------------------------------
mantel_std <- function (v, method) {
  if (method == 'spearman') v <- rank(v)
  v  <- v - mean(v)
  ss <- sum(v * v)
  # Very small or large values: scale first, so that the sum of squares
  # neither underflows nor overflows.
  if (!(is.finite(ss) && ss >= 1e-290) && max(abs(v)) > 0) {
    v  <- v / max(abs(v))
    ss <- sum(v * v)
  }
  v / sqrt(ss)
}


# Kendall: integer codes 0 .. K - 1 of the distinct values, and the pairs of
# pairs tied in codes. Also bioenv()'s.
mantel_codes <- function (v) match(v, sort.int(unique(v))) - 1L
mantel_tied  <- function (codes) {
  cnt <- as.double(tabulate(codes + 1L))
  sum(cnt * (cnt - 1) / 2)
}



#------------------------------------------------------------------------------
# The result, from the observed and permuted statistics (r_all[1] and
# r_all[-1]).
#------------------------------------------------------------------------------
mantel_result <- function (r_all, method, alternative, partial, n, permutation) {

  r      <- r_all[[1]]
  r_perm <- r_all[-1]
  n_perm <- length(r_perm)
  p      <- perm_p(r, r_perm, alternative)

  new_ecodive_test(
    class       = 'ecodive_mantel',
    title       = paste0(
      if (partial) 'Partial Mantel test (' else 'Mantel test (',
      switch(method, pearson = 'Pearson', spearman = 'Spearman', kendall = 'Kendall'), ')' ),
    statistic   = c(r = r),
    p_value     = p,
    z           = perm_z(r, r_perm),
    n_perm      = n_perm,
    method      = method,
    alternative = alternative,
    partial     = partial,
    n           = as.integer(n),
    perm        = r_perm,
    permutation = permutation )
}
