# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Strongly Associated Pairs of Features
#'
#' Lists the pairs of features whose correlation is at least `threshold` in
#' absolute value, and optionally significant, from a [sparcc()] or
#' [rho_p()] result: one row per pair, strongest first.
#'
#' @param x  Correlations between features: a result of [sparcc()] or
#'        [rho_p()], or a symmetric matrix of correlations with 1 on the
#'        diagonal (such as from [stats::cor()]).
#'
#' @param threshold  Keep the pairs whose correlation is at least this in
#'        absolute value (\eqn{|r| \ge} `threshold`), both positive and
#'        negative ones. Default `0.3`, the SparCC paper's cutoff for a
#'        network. `0` keeps every pair.
#'
#' @param alpha  Also keep only the pairs whose adjusted p-value is at most
#'        this. `NULL` (the default) for no p-value filter. Needs p-values:
#'        a [sparcc()] result made with `permutations`.
#'
#' @param p_adjust  How to adjust the p-values for testing every pair: one
#'        of [stats::p.adjust.methods]. Default `'BH'` (Benjamini &
#'        Hochberg's false discovery rate). The adjustment is over all the
#'        pairs of `x`, not just those kept. Used only when `x` has
#'        p-values.
#'
#' @return A data frame of class `c("ecodive_assoc_pairs", "data.frame")`,
#'   one row per pair kept, ordered from the largest \eqn{|r|} to the
#'   smallest, with columns:
#'
#'   * `from`, `to`: the two features' names (their positions, as text,
#'     when the features have no names), `from` being the one that comes
#'     first in the input table.
#'   * `association`: their correlation, with its sign.
#'   * `p_value`, `p_adjusted`: when `x` has p-values, the pair's p-value
#'     and its adjustment by `p_adjust` over all pairs.
#'
#'   Attributes: `method` (the method that made `x`, such as `"sparcc"`;
#'   absent for a matrix), `threshold`, and, when used, `alpha` and
#'   `p_adjust`.
#'
#'   With no pair kept, the data frame has no rows.
#'
#' @details
#' Pairs whose correlation is `NA` or `NaN` are never kept. `threshold` and
#' `alpha` both apply when both are given; use `threshold = 0` to filter on
#' p-values alone. The pairs kept are the edges of a correlation network
#' with this cutoff.
#'
#' **Enough permutations.** A permutation p-value is at least \eqn{1 / (B +
#' 1)} for \eqn{B} null datasets, and adjusting it for testing every pair
#' multiplies it by up to the number of pairs \eqn{m}. So with too few
#' permutations no pair can pass `alpha`, however strong. With
#' Benjamini-Hochberg, \eqn{k} pairs at the smallest p-value pass only
#' when \eqn{m / (k (B + 1)) \le} `alpha`: 20 features (190 pairs) and
#' one strong pair need about 3,800.
#'
#' Only correlations are accepted. [vlr()] and [phi_s()] are dissimilarities
#' (smaller means more proportional), so "above a cutoff" doesn't apply to
#' them: use [rho_p()], their correlation-scale counterpart.
#'
#' @seealso [sparcc()], [rho_p()], [stats::p.adjust()].
#'
#' @export
#' @examples
#'     # 60 samples of 20 taxa, with five correlated pairs
#'     set.seed(1)
#'     basis <- matrix(rnorm(60 * 20), 60, 20)
#'     for (i in c(1, 3, 5, 7, 9))
#'       basis[, i + 1] <- basis[, i] + rnorm(60, sd = 0.5)
#'     counts <- matrix(
#'       rpois(1200, exp(basis + 4)), 60, 20,
#'       dimnames = list(NULL, paste0('OTU', 1:20)) )
#'
#'     r <- sparcc(counts, n_draws = 20)
#'     assoc_pairs(r)
#'
#'     # Only the strongest
#'     assoc_pairs(r, threshold = 0.5)
#'
#'     # Significant pairs, from 999 null datasets
#'     r <- sparcc(counts, n_draws = 10, permutations = 999)
#'     assoc_pairs(r, alpha = 0.05)
#'
#'     # Correlations from elsewhere
#'     assoc_pairs(cor(log(counts + 1)), threshold = 0.5)
#'
assoc_pairs <- function (x, threshold = 0.3, alpha = NULL, p_adjust = 'BH') {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0)
    stop('`threshold` must be a single non-negative number.', call. = FALSE)
  threshold <- as.double(threshold)

  if (!is.null(alpha) && (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha > 1))
    stop('`alpha` must be NULL or a single number above 0 and at most 1.', call. = FALSE)

  if (!is.character(p_adjust) || length(p_adjust) != 1L || !p_adjust %in% stats::p.adjust.methods)
    stop(
      '`p_adjust` must be one of ',
      paste0('"', stats::p.adjust.methods, '"', collapse = ', '), '.', call. = FALSE )

  a <- assoc_pairs_input(x)

  if (!is.null(alpha) && is.null(a$p_value))
    stop(
      '`alpha` needs p-values, and `x` has none. ',
      'Run sparcc() with `permutations` to get them.', call. = FALSE )


  sel <- assoc_select(a, threshold, alpha, p_adjust)
  res <- sel$pairs
  res$pos_from <- res$pos_to <- NULL

  class(res) <- c('ecodive_assoc_pairs', 'data.frame')
  if (!is.null(a$method))   attr(res, 'method')    <- a$method
  attr(res, 'threshold') <- threshold
  if (!is.null(alpha))      attr(res, 'alpha')     <- as.double(alpha)
  if (!is.null(a$p_value))  attr(res, 'p_adjust')  <- p_adjust
  res
}



# The pairs of `a` (from assoc_pairs_input()) that pass the filters:
# |r| >= threshold (|r| > threshold when `strict`), and, when `alpha` is
# given, a p-value adjusted by `p_adjust` over all pairs of at most alpha.
# Returns list(pairs, p_adjusted): `pairs` is a data frame with `from`, `to`
# (labels; `from` is the one first in the input table), `association`,
# `p_value` and `p_adjusted` (when `a` has p-values), and `pos_from`,
# `pos_to` (the two features' input positions), strongest first, ties in
# packed order; `p_adjusted` holds every pair's adjusted p-value (packed
# order), or NULL. Shared by assoc_pairs() and assoc_network().
assoc_select <- function (a, threshold, alpha = NULL, p_adjust = 'BH', strict = FALSE) {

  v    <- a$values
  keep <- !is.na(v) & (if (strict) abs(v) > threshold else abs(v) >= threshold)
  p_adj <- NULL
  if (!is.null(a$p_value)) {
    p_adj <- stats::p.adjust(a$p_value, p_adjust)
    if (!is.null(alpha)) keep <- keep & !is.na(p_adj) & p_adj <= alpha
  }
  k <- which(keep)
  k <- k[order(-abs(v[k]), k)]

  # Packed index k (1-based) -> positions i < j in x's order.
  n      <- a$n
  starts <- c(0, cumsum(seq.int(n - 1L, 1L, by = -1L)))
  i      <- findInterval(k - 1, starts)
  j      <- k - starts[i] + i

  # `from` is the one that comes first in the input table.
  pi    <- a$position[i]
  pj    <- a$position[j]
  first <- pi < pj

  res <- data.frame(
    from        = as.character(ifelse(first, a$labels[i], a$labels[j])),
    to          = as.character(ifelse(first, a$labels[j], a$labels[i])),
    association = v[k],
    stringsAsFactors = FALSE )
  if (!is.null(a$p_value)) {
    res$p_value    <- a$p_value[k]
    res$p_adjusted <- p_adj[k]
  }
  res$pos_from <- as.integer(pmin(pi, pj))
  res$pos_to   <- as.integer(pmax(pi, pj))
  rownames(res) <- NULL

  list(pairs = res, p_adjusted = p_adj)
}



# The correlations of `x` for assoc_pairs(): the packed values (as a dist
# stores them), the number of features `n`, their labels and input
# positions (in x's order), the method, and the parallel p-values (or NULL).
assoc_pairs_input <- function (x) {

  msg <- paste0(
    '`x` must be correlations between features: a result of sparcc() or ',
    'rho_p(), or a symmetric matrix of correlations with 1 on the diagonal.' )

  if (inherits(x, 'dist')) {

    m    <- dist_metric_of(x)
    info <- feature_assoc_info(m)
    if (is.null(info)) stop(msg, call. = FALSE)
    if (!identical(info$kind, 'similarity'))
      stop(
        '`x` is from ', m, '(), a dissimilarity, not a correlation. ',
        'Use rho_p() for proportionality on a correlation scale.', call. = FALSE )

    n   <- as.integer(attr(x, 'Size', exact = TRUE))
    lab <- attr(x, 'Labels', exact = TRUE)
    ord <- attr(x, 'order', exact = TRUE)
    pos <- if (is.numeric(ord) && length(ord) == n) as.integer(ord) else seq_len(n)
    if (is.null(lab)) lab <- as.character(pos)
    p   <- attr(x, 'p_value', exact = TRUE)
    if (!is.numeric(p) || length(p) != length(x)) p <- NULL

    return (list(
      values   = as.vector(unclass(x)),
      n        = n,
      labels   = as.character(lab),
      position = pos,
      method   = m,
      p_value  = p ))
  }

  if (is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x) && nrow(x) >= 2L) {

    d <- diag(x)
    if (!isTRUE(all.equal(unname(x), unname(t(x)))) || anyNA(d) || any(abs(d - 1) > 1e-8))
      stop(msg, call. = FALSE)

    n   <- nrow(x)
    lab <- colnames(x)
    if (is.null(lab)) lab <- rownames(x)
    if (is.null(lab)) lab <- as.character(seq_len(n))

    return (list(
      values   = x[lower.tri(x)],
      n        = n,
      labels   = as.character(lab),
      position = seq_len(n),
      method   = NULL,
      p_value  = NULL ))
  }

  stop(msg, call. = FALSE)
}



# The header of print(): e.g. "SparCC Correlation: 12 pairs with |r| >= 0.3".
assoc_pairs_label <- function (x) {
  thr <- attr(x, 'threshold', exact = TRUE)
  if (is.null(thr)) return (NULL)   # rows subset with `[` drop attributes
  m    <- attr(x, 'method', exact = TRUE)
  name <- if (is.null(m)) 'Correlations' else FEATURE_ASSOCS$name[match(m, FEATURE_ASSOCS$id)]
  out  <- sprintf('%s: pairs with |r| >= %s', name, format(thr))
  a    <- attr(x, 'alpha', exact = TRUE)
  if (!is.null(a))
    out <- sprintf('%s and %s-adjusted p <= %s', out, attr(x, 'p_adjust', exact = TRUE), format(a))
  out
}
