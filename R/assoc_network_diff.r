# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Differential associations between two groups' correlations, by Fisher's
# z. Written from the textbook formula; see plans/assoc_network.md.



#' Differential Associations
#'
#' Tests, for every pair of features, whether their correlation differs
#' between two groups of samples, with Fisher's \eqn{z} test, as NetCoMi's
#' `diffnet(diffMethod = "fisherTest")` does.
#'
#' @param x,y  The two groups' correlations between the same features: two
#'        results of the same method, one whose values are Pearson-type
#'        correlations (today, [sparcc()]), or two symmetric matrices of
#'        Pearson correlations with 1 on the diagonal (with `n_samples`).
#'        Features are matched by name.
#'
#' @param alpha  The adjusted p-value at most which a difference is called
#'        significant (the `reject` column). Default `0.05`.
#'
#' @param p_adjust  How to adjust the p-values for testing every pair: one
#'        of [stats::p.adjust.methods]. Default `'BH'`.
#'
#' @param n_samples  The numbers of samples that `x` and `y` were estimated
#'        from: one number for both, or two. Needed for matrices; `NULL`
#'        (the default) takes them from the results' `n_samples`
#'        attributes.
#'
#' @return A data frame of class `c("ecodive_network_diff", "data.frame")`,
#'   one row per pair of features, with columns:
#'
#'   * `feature1`, `feature2`: the pair, `feature1` first in `x`'s input
#'     table.
#'   * `association_x`, `association_y`: their correlation in each group.
#'   * `z`: Fisher's \eqn{z} statistic.
#'   * `p_value`, `p_adjusted`: its two-sided p-value, and that adjusted
#'     by `p_adjust` over all pairs.
#'   * `reject`: whether `p_adjusted` is at most `alpha`.
#'
#'   Attributes: `method` (the method that made `x` and `y`, when they are
#'   results), `alpha`, `p_adjust`, and `display` (the title; see
#'   [ecodive_display]). `plot()` draws the pairs with `reject` as a
#'   network ([plot.ecodive_network_diff()]).
#'
#' @details
#' Each correlation is transformed with Fisher's \eqn{z = \mathrm{atanh}(r)},
#' and the difference is divided by its standard error:
#'
#' \deqn{z = \frac{\mathrm{atanh}\,r_x - \mathrm{atanh}\,r_y}{\sqrt{1/(n_x
#' - 3) + 1/(n_y - 3)}}}
#'
#' which is standard normal when the two population correlations are equal.
#' Pairs with \eqn{|r| \ge 1} in either group, or a missing correlation,
#' get `NA` (with a warning for \eqn{|r| \ge 1}).
#'
#' For [sparcc()] this is an approximation: its correlations are medians
#' over Dirichlet draws of an estimate, not sample Pearson correlations, so
#' their sampling variance isn't exactly \eqn{1/(n - 3)}. [rho_p()] isn't
#' accepted (a concordance, not a correlation); neither are partial
#' correlations, whose test needs other degrees of freedom.
#'
#' **Names in NetCoMi.** `p_adjust` is `adjust`; `alpha` is `alpha`.
#' NetCoMi's default adjustment is the local false discovery rate; here
#' only [stats::p.adjust()]'s methods are offered. NetCoMi's
#' resampling-based tests (`"permute"`, `"discordant"`) aren't offered.
#'
#' @seealso [assoc_network()], [assoc_network_compare()], [sparcc()].
#'
#' @export
#' @examples
#'     # Two groups of 60 samples; taxa 1 and 2 correlated only in the first
#'     set.seed(1)
#'     make <- function (linked) {
#'       basis <- matrix(rnorm(60 * 10), 60, 10)
#'       if (linked) basis[, 2] <- basis[, 1] + rnorm(60, sd = 0.3)
#'       matrix(rpois(600, exp(basis + 4)), 60, 10,
#'              dimnames = list(NULL, paste0('OTU', 1:10)))
#'     }
#'     rx <- sparcc(make(TRUE),  n_draws = 20)
#'     ry <- sparcc(make(FALSE), n_draws = 20)
#'
#'     d <- assoc_network_diff(rx, ry)
#'     d
#'     d[d$reject, ]
#'
assoc_network_diff <- function (x, y, alpha = 0.05, p_adjust = 'BH', n_samples = NULL) {

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha > 1)
    stop('`alpha` must be a single number above 0 and at most 1.', call. = FALSE)
  p_adjust <- network_p_adjust(p_adjust)

  if (!is.null(n_samples) &&
      (!is.numeric(n_samples) || !length(n_samples) %in% 1:2 || anyNA(n_samples)))
    stop('`n_samples` must be NULL, or one or two numbers.', call. = FALSE)

  ax <- diff_input(x, 'x')
  ay <- diff_input(y, 'y')

  if (!identical(ax$method, ay$method))
    stop(
      '`x` and `y` must come from the same method; they are ',
      diff_method_name(ax$method), ' and ', diff_method_name(ay$method), '.', call. = FALSE )

  # Features, matched by name; x's input order.
  lx <- character(ax$n); lx[ax$position] <- ax$labels
  ly <- character(ay$n); ly[ay$position] <- ay$labels
  if (ax$n != ay$n || !setequal(lx, ly) || anyDuplicated(lx))
    stop('`x` and `y` must be over the same features (matched by name).', call. = FALSE)

  n  <- ax$n
  vx <- packed_in_order(ax, ax$position)
  vy <- packed_in_order(ay, match(ay$labels, lx))

  # Sample sizes: per pair (n_pairs) or per result.
  nx <- diff_n(ax, vx, n_samples, 1L, 'x')
  ny <- diff_n(ay, vy, n_samples, 2L, 'y')

  rx <- vx$values
  ry <- vy$values
  bad <- (!is.na(rx) & abs(rx) >= 1) | (!is.na(ry) & abs(ry) >= 1)
  if (any(bad))
    warning(
      sum(bad), ' pair', if (sum(bad) > 1) 's have' else ' has',
      ' a correlation of 1 or more in absolute value, so no Fisher z (NA).', call. = FALSE )

  z  <- (atanh(rx) - atanh(ry)) / sqrt(1 / (nx - 3) + 1 / (ny - 3))
  z[bad | !is.finite(z)] <- NA_real_
  p  <- 2 * stats::pnorm(-abs(z))
  pa <- stats::p.adjust(p, p_adjust)

  # Pair k of the packed order -> features i < j (input positions).
  starts <- c(0, cumsum(seq.int(n - 1L, 1L, by = -1L)))
  k <- seq_len(length(rx))
  i <- findInterval(k - 1, starts)
  j <- k - starts[i] + i

  res <- data.frame(
    feature1      = lx[i],
    feature2      = lx[j],
    association_x = rx,
    association_y = ry,
    z             = z,
    p_value       = p,
    p_adjusted    = pa,
    reject        = !is.na(pa) & pa <= alpha,
    stringsAsFactors = FALSE )

  class(res) <- c('ecodive_network_diff', 'data.frame')
  if (!is.null(ax$method)) attr(res, 'method') <- ax$method
  attr(res, 'alpha')    <- as.double(alpha)
  attr(res, 'p_adjust') <- p_adjust
  from <- if (inherits(x, 'dist')) display_title(x)
  set_display(res, title = title_of('Differential Associations', from))
}



# The correlations of `x` for assoc_network_diff(), refusing methods whose
# values aren't Pearson-type correlations.
diff_input <- function (x, nm) {

  if (inherits(x, 'dist')) {
    m    <- dist_metric_of(x)
    info <- feature_assoc_info(m)
    if (!is.null(info) && !isTRUE(info$fisher_z))
      stop(
        '`', nm, '` is from ', m, '(), whose values are not Pearson-type ',
        'correlations, so Fisher\'s z does not apply.', call. = FALSE )
  }

  a <- assoc_pairs_input(x)
  a$n_samples <- attr(x, 'n_samples', exact = TRUE)
  np <- attr(x, 'n_pairs', exact = TRUE)
  a$n_pairs <- if (is.numeric(np) && length(np) == length(a$values)) np
  a
}

diff_method_name <- function (m) if (is.null(m)) 'a matrix' else paste0(m, '()')


# The packed values (and n_pairs) of `a` in the order of features
# `target` (each of a's features' new positions, in a's order).
packed_in_order <- function (a, target) {
  n      <- a$n
  starts <- c(0, cumsum(seq.int(n - 1L, 1L, by = -1L)))
  k      <- seq_along(a$values)
  i      <- findInterval(k - 1, starts)
  j      <- k - starts[i] + i
  ti     <- target[i]
  tj     <- target[j]
  lo     <- pmin(ti, tj)
  hi     <- pmax(ti, tj)
  dest   <- starts[lo] + (hi - lo)
  vals   <- numeric(length(k)); vals[dest] <- a$values
  np     <- if (!is.null(a$n_pairs)) { v <- numeric(length(k)); v[dest] <- a$n_pairs; v }
  list(values = vals, n_pairs = np)
}


# The sample size(s) of one input to assoc_network_diff().
diff_n <- function (a, v, n_samples, which, nm) {
  n <- if (!is.null(v$n_pairs)) v$n_pairs
       else if (!is.null(n_samples)) n_samples[[min(which, length(n_samples))]]
       else a$n_samples
  if (is.null(n))
    stop(
      '`', nm, '` has no `n_samples`: give the numbers of samples with ',
      '`n_samples`.', call. = FALSE )
  if (any(n <= 3, na.rm = TRUE))
    stop('Fisher\'s z needs more than 3 samples in each group.', call. = FALSE)
  n
}
