# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Mantel Correlogram
#'
#' How does community similarity change with distance? Splits the distances
#' between samples (geographic, or any other) into classes, and for each
#' class runs a Mantel test of `dm` against a matrix marking the pairs of
#' samples in that class. Gives the same classes, statistics, and p-values
#' as `vegan::mantel.correlog()`, computing every class's test in one pass
#' of parallel C per permutation.
#'
#' @inheritParams mantel
#'
#' @param y  The distances to split into classes: a `dist` object or a
#'        symmetric numeric matrix, matched to `dm`'s samples by name when
#'        both have names, otherwise in the same order; or the samples'
#'        coordinates, a numeric matrix or data frame with one row per sample
#'        (matched by row name, else in the same order), whose Euclidean
#'        distances are used, as `vegan::mantel.correlog()`'s `XY` (a square,
#'        symmetric matrix with a zero diagonal is taken as distances); or
#'        numeric values, one per sample, whose absolute differences are
#'        used. With `data`, the name of one of its numeric columns (one value
#'        per sample), or the names of several (coordinates).
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, whose columns `y` and `strata` can name. Rows are matched to
#'        the samples in `dm` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as `dm`.
#'        Default: `NULL`.
#'
#' @param breaks  The boundaries of the distance classes: an increasing
#'        numeric vector, `n_classes + 1` values from the start of the first
#'        class to the end of the last. Distances outside them are in no
#'        class. Default: `NULL`, for `n_classes` classes of equal width from
#'        the smallest distance in `y` to the largest.
#'
#' @param n_classes  The number of classes of equal width, when `breaks`
#'        isn't given. Default: `NULL`, for Sturges' rule,
#'        \eqn{\lceil \log_2 N + 1 \rceil} for \eqn{N} pairs of samples, as
#'        in `vegan::mantel.correlog()`.
#'
#' @param cutoff  `TRUE` to test, in the second half of the classes, only the
#'        classes whose pairs include every sample, as
#'        `vegan::mantel.correlog()` does; `FALSE` to test every class.
#'        Default: `TRUE`.
#'
#' @param alternative  Which correlations count as more extreme than each
#'        class's observed one: `'sign'`, in the direction of its sign (more
#'        positive when it is positive, more negative otherwise), as in
#'        `vegan::mantel.correlog()`; `'greater'` (more positive); `'less'`
#'        (more negative); or `'two.sided'` (larger in absolute value).
#'        Default: `'sign'`.
#'
#' @param p_adjust  How to adjust the classes' p-values for multiple testing:
#'        any method of [stats::p.adjust()] (`'none'` for no adjustment), or
#'        `'sidak'`. Default: `'holm'`.
#'
#' @param progressive  `TRUE` for a progressive adjustment, as
#'        `vegan::mantel.correlog()` does by default (Legendre & Legendre
#'        2012): the first tested class's p-value is not adjusted, and the
#'        k-th's is adjusted together with those of the k - 1 tested classes
#'        before it, at shorter distances. `FALSE` adjusts all the tested
#'        classes' p-values together. Default: `TRUE`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones. Each
#'        permutation shuffles the samples of `y` (the distance classes),
#'        with `dm` held fixed, as in `vegan::mantel.correlog()`.
#'
#' @return A data frame of class
#'   `c('ecodive_mantel_correlogram', 'data.frame')`, with one row per
#'   distance class, in order of distance:
#'
#'   * `class_from`, `class_to`: The class's boundaries. A class holds the
#'     pairs whose distance in `y` is greater than `class_from` and at most
#'     `class_to`; the first class also holds those equal to its
#'     `class_from`.
#'   * `class_mid`: The class's midpoint, `vegan::mantel.correlog()`'s
#'     `class.index`.
#'   * `n_pairs`: The number of pairs of samples in the class.
#'     `vegan::mantel.correlog()`'s `n.dist` is twice this: it counts each
#'     pair twice, as in a full distance matrix.
#'   * `r`: The Mantel statistic, the correlation of `dm` with a matrix that
#'     is 0 for the pairs in the class and 1 for the others. It is positive
#'     when samples in the class are more similar than other pairs (positive
#'     autocorrelation). `NA` for a class that isn't tested.
#'   * `p_value`: Its permutation p-value.
#'   * `p_adjusted`: `p_value` adjusted for multiple testing, with
#'     `p_adjust` and `progressive`.
#'
#'   A class is tested when it has at least one pair of samples, but not
#'   every pair; with `cutoff = TRUE`, the classes after the first half
#'   (after class \eqn{\lfloor K/2 \rfloor} of \eqn{K}) are tested only when
#'   their pairs include every sample.
#'
#'   Attributes: `breaks` (the classes' boundaries), `method`,
#'   `alternative`, `p_adjust`, `progressive`, `n_perm`, `permutation`
#'   (`'free'`, `'within strata'`, or `'supplied matrix'`), and `r_perm`, a
#'   matrix of the permutations' statistics, one row per class (`NA` for
#'   classes that aren't tested). `plot()` draws the correlogram.
#'
#' @details
#'
#' **Method**
#'
#' A Mantel correlogram (Sokal 1986; Oden & Sokal 1986; Legendre & Legendre
#' 2012) shows how the similarity of samples changes with the distance
#' between them. The pairs of samples are grouped into classes by their
#' distance in `y`, and for each class \eqn{k}, the Mantel statistic is the
#' correlation of the distances `dm` with a model matrix that is 0 for the
#' pairs in class \eqn{k} and 1 for the others. So a positive statistic
#' means that samples in that class are more similar than other pairs, and a
#' negative one that they are less similar. Each statistic is tested by
#' permuting the samples, as in [mantel()]; all classes share the same
#' permutations. Borcard & Legendre (2012) found that the test has correct
#' type I error, unlike simple and partial Mantel tests of spatial
#' structure.
#'
#' By default, the p-values are one-tailed in the direction of each
#' statistic's sign, \eqn{(1 + \#\{r_{perm} \ge r\}) / (1 + permutations)}
#' when \eqn{r > 0}, and \eqn{(1 + \#\{r_{perm} \le r\}) / (1 + permutations)}
#' otherwise (a statistic within `1e-12` of 0 counts as 0), as in
#' `vegan::mantel.correlog()`. As in `vegan::mantel()`,
#' permuted statistics within `sqrt(.Machine$double.eps)` of the observed one
#' count as ties. More permutations are needed than for a single test, for
#' the adjusted p-values of the later classes to reach significance.
#'
#'
#' **Speed and memory**
#'
#' The correlation of `dm` with a class's 0/1 matrix depends on the
#' permutation only through the sum of `dm` over the pairs in the class:
#' Pearson's r on the standardized distances (Spearman's on their ranks),
#' and Kendall's tau-b on their midranks, as a Mann-Whitney statistic. So
#' each permutation is one pass over the \eqn{n(n-1)/2} pairs, which sums
#' every class at once, with no \eqn{n \times n} matrix, where
#' `vegan::mantel.correlog()` runs a separate Mantel test per class (and
#' Kendall's tau in \eqn{O(n^4)} each). Permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and p-values are identical from run to run
#' and for any number of `cpus`. They differ from `vegan::mantel.correlog()`'s
#' random permutations, which differ between classes, so p-values agree with
#' it only up to sampling error.
#'
#'
#' **Differences from vegan::mantel.correlog()**
#'
#' * `y` may also be numeric values per sample, or columns of `data`.
#'   Samples are matched by name between `dm`, `y`, `data`, and `strata`.
#' * `n_pairs` counts each pair of samples once; `n.dist` counts it twice.
#' * `alternative` adds upper-tail, lower-tail, and two-sided tests.
#' * vegan lowers the first boundary by `.Machine$double.eps` to put the
#'   smallest distance in the first class, which fails for distances of
#'   about 2 or more: their closest pairs then fall in no class. Here the
#'   first class includes its lower boundary.
#' * When a class is left untested between tested classes (it has no pairs,
#'   or, with `cutoff`, lacks a sample), vegan adjusts only the p-values
#'   before it, and gives the later classes no adjusted p-value. Here every
#'   tested class's p-value is adjusted, in order of distance.
#' * Permutations can be restricted to `strata`, or given as a matrix,
#'   which every class shares.
#'
#' @references
#' Borcard, D., & Legendre, P. (2012). Is the Mantel correlogram powerful
#' enough to be useful in ecological analysis? A simulation study.
#' *Ecology*, 93(6), 1473-1481. \doi{10.1890/11-1737.1}
#'
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.),
#' section 13.1. Elsevier.
#'
#' Oden, N. L., & Sokal, R. R. (1986). Directional autocorrelation: an
#' extension of spatial correlograms to two dimensions. *Systematic
#' Zoology*, 35(4), 608-617. \doi{10.2307/2413124}
#'
#' Sokal, R. R. (1986). Spatial data analysis and historical processes. In
#' E. Diday et al. (Eds.), *Data analysis and informatics IV* (pp. 29-43).
#' North-Holland.
#'
#' Sturges, H. A. (1926). The choice of a class interval. *Journal of the
#' American Statistical Association*, 21(153), 65-66.
#' \doi{10.1080/01621459.1926.10502161}
#'
#' @seealso [mantel()], [mantel_partial()], [plot.ecodive_mantel_correlogram()]
#'
#' @export
#' @examples
#'     # 40 samples on a 2-D landscape, whose communities change gradually
#'     # across it.
#'     set.seed(1)
#'     meta <- data.frame(
#'       row.names = paste0('S', 1:40),
#'       x = runif(40, 0, 10),
#'       y = runif(40, 0, 10) )
#'     counts <- matrix(rpois(40 * 20, 5), 40, 20)
#'     counts[, 1:5] <- counts[, 1:5] + round(2 * meta$x)
#'     rownames(counts) <- rownames(meta)
#'
#'     dm <- bray(counts)
#'
#'     # Are nearby samples more similar, and distant ones less?
#'     res <- mantel_correlogram(dm, c('x', 'y'), meta)
#'     res
#'     plot(res)
#'
#'     # Five classes, tested with Spearman's rank correlation.
#'     mantel_correlogram(dm, dist(meta), n_classes = 5, method = 'spearman')
#'
mantel_correlogram <- function (
    dm, y, data = NULL, breaks = NULL, n_classes = NULL, cutoff = TRUE,
    method = 'pearson', permutations = 999L, p_adjust = 'holm',
    progressive = TRUE, strata = NULL, alternative = 'sign', seed = 0L,
    cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  method      <- match.arg(method, c('pearson', 'spearman', 'kendall'))
  alternative <- match.arg(alternative, c('sign', 'greater', 'less', 'two.sided'))
  p_adjust    <- correlogram_p_adjust(p_adjust)

  for (arg in c('cutoff', 'progressive')) {
    v <- get(arg)
    if (!isTRUE(v) && !isFALSE(v))
      stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)
  }

  if (!is.null(breaks) && !is.null(n_classes))
    stop('Give `breaks` or `n_classes`, not both.', call. = FALSE)

  if (missing(y)) y <- NULL
  inp    <- mantel_inputs(dm, y, data, NULL, strata, coords = TRUE)
  n      <- inp$n
  x      <- inp$x
  yv     <- inp$y
  strata <- inp$strata
  N      <- length(x)

  pm <- perm_matrix(permutations, n)
  if (!is.null(pm$perms) && !is.null(strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  #--------------------------------------------------------
  # The classes: breaks, each pair's class, and which
  # classes are tested.
  #--------------------------------------------------------
  breaks <- correlogram_breaks(breaks, n_classes, yv)
  K      <- length(breaks) - 1L

  cls <- findInterval(yv, breaks, left.open = TRUE)
  cls[yv == breaks[[1]]] <- 1L
  cls[cls < 1L | cls > K] <- NA_integer_

  n_pairs <- tabulate(cls, K)

  tested <- n_pairs > 0 & n_pairs < N
  if (cutoff) {
    late <- which(tested & seq_len(K) > K %/% 2)
    if (length(late)) {
      pj <- rep.int(seq_len(n - 1L), (n - 1L):1)
      pi <- pj + sequence((n - 1L):1)
      for (k in late) {
        in_k <- which(cls == k)
        tested[[k]] <- all(tabulate(c(pi[in_k], pj[in_k]), n) > 0)
      }
    }
  }


  #--------------------------------------------------------
  # The classes' statistics, observed and permuted.
  #--------------------------------------------------------
  r_all <- matrix(NA_real_, K, pm$n_perm + 1L)
  idx   <- which(tested)

  if (all(x == x[[1]])) {
    warning('`dm` has the same distance between every pair of samples, so ',
            'it has no correlation with anything.', call. = FALSE)
    idx <- integer(0)
  }

  if (length(idx)) {

    code <- match(cls, idx) - 1L
    code[is.na(code)] <- -1L

    nk <- as.double(n_pairs[idx])

    if (method == 'kendall') {
      # Twice each distance's midrank, less N + 1: their sum over a class
      # is the class's Kendall score, concordant less discordant pairs.
      v  <- as.integer(rank(x) * 2 - (N + 1))
      n0 <- N * (N - 1) / 2
      n1 <- mantel_tied(mantel_codes(x))
      n2 <- nk * (nk - 1) / 2 + (N - nk) * (N - nk - 1) / 2
      den <- sqrt(n0 - n1) * sqrt(n0 - n2)
      cen <- 0
    }
    else {
      # Standardized (ranked for Spearman): the correlation with a class's
      # indicator is its centred sum over the class, scaled.
      v   <- mantel_std(x, method)
      den <- sqrt(nk * (N - nk) / N)
      cen <- nk / N * sum(v)
    }

    s <- .Call(
      C_mantel_classes, v, code, pm$perms, strata,
      as.integer(c(n, length(idx), pm$n_perm, seed, cpus)) )

    # vegan's sign: the model matrix is 0 in the class, 1 elsewhere.
    r_all[idx, ] <- -(s - cen) / den
  }

  r      <- r_all[, 1]
  r_perm <- r_all[, -1, drop = FALSE]


  #--------------------------------------------------------
  # P-values, and their adjustment over the tested classes
  # in order of distance.
  #--------------------------------------------------------
  p <- rep(NA_real_, K)
  for (k in idx) {
    # The direction of r's sign; r within rounding of 0 counts as 0.
    alt  <- if (alternative != 'sign') alternative else if (r[[k]] > 1e-12) 'greater' else 'less'
    p[k] <- perm_p(r[[k]], r_perm[k, ], alt)
  }

  p_adj <- rep(NA_real_, K)
  if (length(idx)) {
    adj <- function (pv) {
      if (p_adjust == 'sidak') 1 - (1 - pv) ^ length(pv)
      else stats::p.adjust(pv, p_adjust)
    }
    if (progressive) {
      # The first class's p-value is unadjusted, exactly.
      p_adj[idx[[1]]] <- p[idx[[1]]]
      for (j in seq_along(idx)[-1])
        p_adj[idx[[j]]] <- utils::tail(adj(p[idx[seq_len(j)]]), 1)
    }
    else {
      p_adj[idx] <- adj(p[idx])
    }
  }


  out <- data.frame(
    class_from = breaks[-(K + 1L)],
    class_to   = breaks[-1],
    class_mid  = (breaks[-(K + 1L)] + breaks[-1]) / 2,
    n_pairs    = n_pairs,
    r          = r,
    p_value    = p,
    p_adjusted = p_adj )

  attr(out, 'breaks')      <- breaks
  attr(out, 'method')      <- method
  attr(out, 'alternative') <- alternative
  attr(out, 'p_adjust')    <- p_adjust
  attr(out, 'progressive') <- progressive
  attr(out, 'n_perm')      <- pm$n_perm
  attr(out, 'permutation') <- if (!is.null(pm$perms)) 'supplied matrix'
                              else if (!is.null(strata)) 'within strata'
                              else 'free'
  attr(out, 'r_perm')      <- r_perm
  class(out) <- c('ecodive_mantel_correlogram', 'data.frame')

  set_display(out, title = title_of('Mantel correlogram', inp$from))
}



#------------------------------------------------------------------------------
# `p_adjust`: a method of stats::p.adjust(), or vegan's 'sidak'.
#------------------------------------------------------------------------------
correlogram_p_adjust <- function (p_adjust) {
  opts <- c(stats::p.adjust.methods, 'sidak')
  m    <- if (is.character(p_adjust) && length(p_adjust) == 1 && !is.na(p_adjust))
            pmatch(p_adjust, opts)
          else NA_integer_
  if (is.na(m))
    stop(
      '`p_adjust` must be one of ', paste0("'", opts, "'", collapse = ', '),
      ' (see stats::p.adjust()).', call. = FALSE )
  opts[[m]]
}



#------------------------------------------------------------------------------
# The classes' boundaries: `breaks` as given, or `n_classes` (by default,
# Sturges' rule) classes of equal width over the range of y's distances, as
# vegan::mantel.correlog().
#------------------------------------------------------------------------------
correlogram_breaks <- function (breaks, n_classes, yv) {

  if (!is.null(breaks)) {
    if (!is.numeric(breaks) || length(breaks) < 2 || !all(is.finite(breaks)) ||
        any(diff(breaks) <= 0))
      stop('`breaks` must be at least two finite numbers, in increasing order.', call. = FALSE)
    return (as.double(breaks))
  }

  if (is.null(n_classes)) {
    n_classes <- ceiling(log2(length(yv)) + 1)
  }
  else if (!is.numeric(n_classes) || length(n_classes) != 1 || !is.finite(n_classes) ||
           n_classes < 1 || n_classes != round(n_classes)) {
    stop('`n_classes` must be a whole number of at least 1.', call. = FALSE)
  }

  lo <- min(yv)
  hi <- max(yv)
  if (!(hi > lo))
    stop('`y` has the same distance between every pair of samples, so it ',
         'has no distance classes.', call. = FALSE)

  seq(lo, hi, length.out = n_classes + 1)
}



#------------------------------------------------------------------------------
#' Plot a Mantel correlogram
#'
#' Draws the Mantel statistic of each tested distance class against the
#' class's midpoint, as `plot()` of `vegan::mantel.correlog()` does: filled
#' squares for the classes whose adjusted p-value is at most `alpha`, and
#' open squares for the others.
#'
#' @param x  A [mantel_correlogram()] result.
#'
#' @param alpha  The significance level of the filled squares, compared with
#'        `p_adjusted` (with `p_value` when there are no adjusted p-values).
#'        Default: `0.05`.
#'
#' @param ...  Passed on to [graphics::plot()]: `main`, `xlab`, `ylab`,
#'        `col`, `xlim`, `ylim`, and so on.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     xy <- matrix(runif(60, 0, 10), 30, 2)
#'     counts <- matrix(rpois(30 * 10, 5), 30, 10)
#'     counts[, 1:3] <- counts[, 1:3] + round(2 * xy[, 1])
#'
#'     res <- mantel_correlogram(bray(counts), xy)
#'     plot(res)
#'     plot(res, alpha = 0.01, main = 'Distance decay')
#'
plot.ecodive_mantel_correlogram <- function (x, alpha = 0.05, ...) {

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1)
    stop('`alpha` must be a single number from 0 to 1.', call. = FALSE)

  if (!all(c('class_mid', 'r', 'p_value', 'p_adjusted') %in% names(x)))
    stop('`x` must be a mantel_correlogram() result.', call. = FALSE)

  ok <- !is.na(x$r)
  if (!any(ok))
    stop('`x` has no tested classes to plot.', call. = FALSE)

  mid <- x$class_mid[ok]
  r   <- x$r[ok]
  p   <- if (any(!is.na(x$p_adjusted[ok]))) x$p_adjusted[ok] else x$p_value[ok]
  sig <- !is.na(p) & p <= alpha

  args <- list(...)
  col  <- plot_arg(args, 'col', 'black')
  meth <- attr(x, 'method', exact = TRUE)

  plot_do(graphics::plot, plot_drop(args, 'col'), list(
    x    = mid,
    y    = r,
    type = 'n',
    main = plot_title(x, 'Mantel correlogram'),
    xlab = 'Distance class',
    ylab = if (is.character(meth) && length(meth) == 1)
             sprintf('Mantel r (%s)', plot_cap(meth)) else 'Mantel r' ))

  graphics::abline(h = 0, col = PLOT_MUTED, lty = 2)
  graphics::lines(mid, r, col = col)
  graphics::points(mid, r, pch = ifelse(sig, 15, 0), col = col, bg = col)

  invisible(x)
}
