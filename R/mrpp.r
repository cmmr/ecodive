# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Multi-Response Permutation Procedure (MRPP)
#'
#' Tests whether samples are more similar to the other samples in their group
#' than to samples in other groups. Gives the same statistics as
#' `vegan::mrpp()`: the observed and expected weighted mean within-group
#' distance (\eqn{\delta}), the chance-corrected within-group agreement
#' \eqn{A}, the mean distance within each group, and a permutation p-value.
#' Also gives the standardized effect size \eqn{z} that
#' `summary(vegan::permustats())` reports. Computes them in parallel C,
#' straight from the distances.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values.
#'
#' @param groups  Each sample's group: a vector or factor with one value per
#'        sample, matched to the samples in `dm` by name when both are named,
#'        else in the same order as `dm`. Or, when `data` is given, the name of
#'        one of its columns. Values are treated as categories, even if they
#'        are numbers.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `groups` (and `strata`) as columns. Rows are matched
#'        to the samples in `dm` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as `dm`.
#'        Default: `NULL`.
#'
#' @param weights  How to weight each group's mean within-group distance in
#'        \eqn{\delta}, by its number of samples \eqn{n}:
#'
#'   * `'n'`: \eqn{n}, as recommended by Mielke & Berry (2001). The default,
#'     and `vegan::mrpp()`'s `weight.type = 1`.
#'   * `'n-1'`: \eqn{n - 1}. `weight.type = 2`.
#'   * `'n(n-1)'`: \eqn{n(n - 1)}, the number of pairs of samples, so that
#'     \eqn{\delta} is the mean of all within-group distances.
#'     `weight.type = 3`.
#'
#'   `1`, `2`, or `3` are also accepted, as in `vegan::mrpp()`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param strata  Permute samples only within these groups (blocks), for
#'        designs such as repeated measures. The name of a column of `data`,
#'        or a vector with one value per sample (per row of `data`, when
#'        given; else in the same order as `groups`). Default: `NULL` (free
#'        permutation).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A list of class `c("ecodive_mrpp", "ecodive_test")`, which prints
#'   as a summary. `as.data.frame()` gives its first four fields as one row,
#'   with the same columns as ecodive's other permutation tests. Its fields:
#'
#'   * `statistic`: The observed \eqn{\delta}, named `delta`: the weighted
#'     mean of the groups' mean within-group distances.
#'   * `p_value`: The permutation p-value: the proportion of permutations,
#'     counting the observed one, with a \eqn{\delta} as small as or smaller
#'     than the observed one.
#'   * `z`: The standardized effect size,
#'     `(delta - mean(delta_perm)) / sd(delta_perm)`. Negative when samples
#'     are closer within groups than expected by chance.
#'   * `n_perm`: The number of permutations.
#'   * `expected_delta`: The expected \eqn{\delta} under the null hypothesis:
#'     the mean of all the distances.
#'   * `A`: Chance-corrected within-group agreement,
#'     `1 - delta / expected_delta`.
#'   * `class_delta`: Each group's mean within-group distance (`NA` for a
#'     group of one sample).
#'   * `class_n`: Each group's number of samples.
#'   * `weights`: The `weights` used.
#'   * `permutation`: How samples were permuted: `'free'`, `'within strata'`,
#'     or `'supplied matrix'`.
#'
#'   Attribute `delta_perm` holds the \eqn{\delta} of each permutation, as
#'   `vegan::mrpp()`'s `boot.deltas`.
#'
#' @details
#'
#' **Method**
#'
#' MRPP (Mielke et al. 1976; Mielke & Berry 2001) asks whether samples in the
#' same group are closer to each other than samples drawn at random. With
#' \eqn{\bar{d}_g} the mean distance between the \eqn{n_g (n_g - 1) / 2}
#' pairs of samples in group \eqn{g}, the test statistic is their weighted
#' mean
#' \deqn{\delta = \sum_g C_g \bar{d}_g, \quad C_g = w_g / \sum_h w_h}
#' with weights \eqn{w_g} from `weights`. Groups of one sample have no pairs,
#' so they are left out, as in `vegan::mrpp()`. The expected value of
#' \eqn{\delta} over all permutations of the samples is the mean of all the
#' distances, and the chance-corrected within-group agreement is
#' \eqn{A = 1 - \delta / E(\delta)} (McCune & Grace 2002): 1 when samples are
#' identical within groups, 0 when groups are no more alike than chance, and
#' negative when they are less alike.
#'
#' For each permutation, the group labels are shuffled among the samples and
#' \eqn{\delta} is recomputed. The p-value is
#' \eqn{(1 + \#\{\delta_{perm} \le \delta\}) / (1 + permutations)}. As in
#' `vegan::mrpp()`, permuted values within `sqrt(.Machine$double.eps)` of the
#' observed one count as ties. The standardized effect size \eqn{z} compares
#' the observed \eqn{\delta} with the permuted ones, using their mean and
#' (sample) standard deviation, as `summary(vegan::permustats())` does. It is
#' `NA` with fewer than two permutations, or when every permutation gives the
#' same \eqn{\delta}.
#'
#'
#' **Speed and memory**
#'
#' With \eqn{S_g} the sum of the distances within group \eqn{g},
#' \eqn{\delta = \sum_g c_g S_g} for constants \eqn{c_g}. So each
#' permutation costs one comparison per pair of samples, read straight from
#' `dm`, with no \eqn{n \times n} matrix. Many permutations share each pass
#' over the distances, and permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. Results are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `vegan::mrpp()`'s,
#' so p-values and \eqn{z} agree with it only up to sampling error, unless
#' the same `permutations` matrix is given to both.
#'
#'
#' **Differences from vegan::mrpp()**
#'
#' * `mrpp()` takes distances, not a community table: compute them first
#'   with [beta_div()] or a metric function such as [bray()].
#' * `A` is reported even with `permutations = 0`; `vegan::mrpp()` gives
#'   `NA` then.
#' * A group of one sample has a `class_delta` of `NA`, where
#'   `vegan::mrpp()` has `NaN`.
#' * `mrpp()` does not compute the classification strength (`CS`), which
#'   `vegan::mrpp()` also leaves `NA`.
#' * `strata` permutes freely within strata. For other restricted designs,
#'   pass a `permutations` matrix from `permute::shuffleSet()`.
#'
#' @references
#' McCune, B., & Grace, J. B. (2002). *Analysis of Ecological Communities*.
#' MjM Software Design, Gleneden Beach, Oregon. Chapter 24.
#'
#' Mielke, P. W., Berry, K. J., & Johnson, E. S. (1976). Multi-response
#' permutation procedures for a priori classifications. *Communications in
#' Statistics - Theory and Methods*, 5(14), 1409-1424.
#' \doi{10.1080/03610927608827451}
#'
#' Mielke, P. W., & Berry, K. J. (2001). *Permutation Methods: A Distance
#' Function Approach*. Springer, New York. \doi{10.1007/978-1-4757-3449-2}
#'
#' @seealso [permanova()], [beta_div()]
#'
#' @export
#' @examples
#'     # 30 samples in three treatment groups. The 'high' group has more
#'     # of the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       treatment = rep(c('ctrl', 'low', 'high'), each = 10),
#'       site      = rep(c('A', 'B'), 15) )
#'     high <- meta$treatment == 'high'
#'     counts[high, 1:5] <- counts[high, 1:5] + 5
#'
#'     dm <- bray(counts)
#'
#'     # Are samples more alike within treatment groups?
#'     mrpp(dm, 'treatment', meta)
#'
#'     # The same, with the groups as a vector.
#'     res <- mrpp(dm, meta$treatment)
#'     as.data.frame(res)
#'
#'     # Permuting only within sites.
#'     mrpp(dm, 'treatment', meta, strata = 'site')
#'
mrpp <- function (
    dm, groups, data = NULL, weights = 'n', permutations = 999L,
    strata = NULL, seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  wt  <- mrpp_weights(weights)
  dv  <- dist_vector(dm)
  n   <- dv$n

  if (n < 2)
    stop('`dm` must have at least 2 samples.', call. = FALSE)

  if (missing(groups) || is.null(groups))
    stop('`groups` is required: one group per sample, or a column name of `data`.', call. = FALSE)

  pm <- perm_matrix(permutations, n)


  #--------------------------------------------------------
  # Each sample's group, in the order of dm's samples.
  #--------------------------------------------------------
  if (is.null(data)) {

    if (!is.atomic(groups) || !is.null(dim(groups)))
      stop('`groups` must be a vector or factor with one value per sample.', call. = FALSE)
    if (is.character(groups) && length(groups) == 1)
      stop('`groups` names a column, but `data` is not given.', call. = FALSE)
    if (anyNA(groups))
      stop('`groups` must not contain NAs.', call. = FALSE)
    if (length(groups) != n && (is.null(names(groups)) || is.null(dv$labels)))
      stop(sprintf('`groups` has %d values, but `dm` has %d samples.', length(groups), n), call. = FALSE)

    if (!is.null(strata) && !(is.atomic(strata) && length(strata) == length(groups)))
      stop('Without `data`, `strata` must have one value per value of `groups`.', call. = FALSE)

    des <- model_design(
      groups, formula = ~ 1, labels = dv$labels, n = n, strata = strata,
      allow_empty = TRUE )
    grp <- des$data[[1]]
  }
  else {

    if (!is.character(groups) || length(groups) != 1 || is.na(groups))
      stop('With `data`, `groups` must be the name of one of its columns.', call. = FALSE)
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    if (!groups %in% colnames(data))
      stop('`groups` is not a column of `data`: ', name_list(groups), call. = FALSE)

    des <- model_design(
      data, formula = ~ 1, labels = dv$labels, n = n, strata = strata,
      allow_empty = TRUE )
    grp <- des$data[[groups]]

    if (anyNA(grp))
      stop('Missing values (NA) in `data` column ', name_list(groups), '.', call. = FALSE)
  }

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  grp   <- factor(grp)
  lvls  <- levels(grp)
  gi    <- as.integer(grp)
  ng    <- tabulate(gi, length(lvls))
  pairs <- ng * (ng - 1) / 2


  #--------------------------------------------------------
  # delta = sum_g c_g S_g, over groups with 2+ samples.
  #--------------------------------------------------------
  w     <- switch(wt, n = ng, `n-1` = ng - 1, `n(n-1)` = pairs)
  keep  <- ng >= 2
  coef  <- ifelse(keep, w / (sum(w[keep]) * pairs), 0)

  if (!any(keep))
    warning('No group has two or more samples, so delta is undefined.', call. = FALSE)
  else if (length(lvls) == 1)
    warning('All samples are in one group, so there is nothing to test.', call. = FALSE)

  res <- .Call(
    C_mrpp, dv$dm, as.integer(n), gi - 1L, as.double(coef), pm$perms,
    des$strata, as.integer(c(pm$n_perm, seed, cpus)) )

  delta      <- res$delta[[1]]
  delta_perm <- res$delta[-1]
  e_delta    <- mean(dv$dm)

  if (!any(keep)) {
    delta      <- NA_real_
    delta_perm <- rep(NA_real_, pm$n_perm)
  }

  class_delta <- ifelse(keep, res$class_sum / pairs, NA_real_)
  names(class_delta) <- names(ng) <- lvls

  A <- if (isTRUE(e_delta > 0)) 1 - delta / e_delta else NA_real_
  if (isTRUE(e_delta == 0))
    warning('All distances are 0, so there is nothing to test.', call. = FALSE)


  #--------------------------------------------------------
  # p-value and z-score from the permuted deltas.
  #--------------------------------------------------------
  eps <- sqrt(.Machine$double.eps)

  p <- if (pm$n_perm > 0 && !is.na(delta))
    (1 + sum(delta_perm <= delta + eps)) / (pm$n_perm + 1)
  else NA_real_

  new_ecodive_test(
    class          = 'ecodive_mrpp',
    title          = 'Multi-response permutation procedure (MRPP)',
    statistic      = c(delta = delta),
    p_value        = p,
    z              = perm_z(delta, delta_perm),
    n_perm         = pm$n_perm,
    expected_delta = e_delta,
    A              = A,
    class_delta    = class_delta,
    class_n        = ng,
    weights        = wt,
    perm           = delta_perm,
    permutation    = if (!is.null(pm$perms)) 'supplied matrix'
                     else if (!is.null(des$strata)) 'within strata'
                     else 'free' )
}



# `weights` as one of 'n', 'n-1', 'n(n-1)'. Also takes vegan's weight.type
# (1, 2, 3) and 'n(n-1)/2'.
mrpp_weights <- function (weights) {

  opts <- c('n', 'n-1', 'n(n-1)')

  ok <- length(weights) == 1 && !is.na(weights)
  if (ok && is.numeric(weights) && weights %in% 1:3)
    return (opts[[weights]])

  if (ok && is.character(weights)) {
    w <- gsub(' ', '', weights, fixed = TRUE)
    if (w == 'n(n-1)/2') w <- 'n(n-1)'
    if (w %in% opts) return (w)
  }

  stop("`weights` must be 'n', 'n-1', or 'n(n-1)' (or 1, 2, or 3).", call. = FALSE)
}



# Printed after the core results: each group's mean distance and size, and
# the expected delta and A.
#' @export
test_details.ecodive_mrpp <- function (x, digits) {

  fmt <- function (v) ifelse(is.na(v), 'NA', formatC(v, digits = digits, format = 'g'))

  tbl <- rbind(delta = fmt(x$class_delta), n = as.character(x$class_n))

  c(
    paste0('Mean within-group distances (weights: ', x$weights, '):'),
    utils::capture.output(print(noquote(tbl), right = TRUE)),
    '',
    paste0('Expected delta: ', fmt(x$expected_delta)),
    paste0('Chance-corrected within-group agreement A: ', fmt(x$A)) )
}
