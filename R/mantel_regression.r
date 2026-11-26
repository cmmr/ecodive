# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Multiple Regression on Distance Matrices
#'
#' Regresses one distance matrix on several others: for instance, how much
#' of the differences in community composition between samples is explained
#' by their geographic distances and by their differences in pH and
#' temperature, and which of these matter once the others are accounted
#' for. This is multiple regression on distance matrices (MRM), the
#' multiple-regression extension of the [mantel()] test. It gives the same
#' coefficients, \eqn{R^2}, \eqn{F}, and permutation p-values as
#' `ecodist::MRM()`, computed in parallel C straight from the distances.
#'
#' @inheritParams documentation
#'
#' @param dm  The response: a `dist` object, such as from [beta_div()] or
#'        [bray()], or a symmetric numeric matrix of distances. Must not
#'        contain `NA`s. Its samples are the ones permuted. Alternatively,
#'        a formula with the distances on its left-hand side and the
#'        explanatory variables on its right, such as `dm ~ ph + temp`,
#'        in place of `formula`.
#'
#' @param ...  Explanatory distances, each one term of the model, named as
#'        the term is to be called in the results (unnamed ones are called
#'        by what was written): such as `geo = dist(xy)`. Each is a `dist`
#'        object or a symmetric numeric matrix of distances, matched to `dm`'s
#'        samples by name when both have names, otherwise in the same order;
#'        or one value per sample (a vector, or, with `data`, the name of a
#'        column), whose differences are the distances, as for `formula`'s
#'        terms; or a numeric matrix or data frame of coordinates, one row per
#'        sample, whose Euclidean distances are used.
#'
#' @param data  Sample metadata: a data frame with one row per sample, whose
#'        columns `formula`, `...`, and `strata` can name. Rows are matched to
#'        the samples in `dm` by row name, and rows for other samples are
#'        ignored; without row names, they must be in the same order as `dm`.
#'        Default: `NULL`.
#'
#' @param formula  A one-sided formula naming variables in `data`, each one
#'        term of the model, such as `~ ph + temp + cbind(lat, lon)`. A
#'        numeric variable's distances are its absolute differences between
#'        samples, \eqn{|x_i - x_j|}, as for [mantel()]'s `y`; a numeric
#'        matrix's (such as `cbind(lat, lon)`) are its Euclidean distances;
#'        and a factor's, character's, or logical's are 0 for two samples with
#'        the same value and 1 for two with different values. Terms may be
#'        expressions of the variables, such as `log(depth)`, but not
#'        interactions. `formula`'s terms come first in the results, then
#'        those of `...`. Default: `NULL`.
#'
#' @param method  `'pearson'`, for ordinary least squares regression of the
#'        distances, or `'spearman'`, for the same regression of their ranks:
#'        `ecodist::MRM()`'s `mrank = FALSE` and `mrank = TRUE`. Default:
#'        `'pearson'`.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones. Give
#'        `0` for no p-values.
#'
#' @param backward  Whether to simplify the model by backward elimination
#'        (Lichstein 2007): remove the least significant term, refit, and
#'        repeat until every term left has a p-value of at most `alpha`. See
#'        **Backward elimination**. Default: `FALSE`.
#'
#' @param alpha  With `backward = TRUE`, the significance level a term's
#'        p-value must reach to stay in the model. Default: `0.05`.
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A list of class `ecodive_mantel_regression`, with:
#'
#'   * `coefficients`: A data frame with one row for the intercept and one
#'     per term, and columns:
#'     * `term`: `'(Intercept)'`, then the terms' names.
#'     * `estimate`: The regression coefficient, in units of `dm`'s
#'       distances per unit of the term's distances (of their ranks, for
#'       `'spearman'`). `ecodist::MRM()`'s `coef`.
#'     * `std_estimate`: The standardized coefficient,
#'       `estimate * sd(term) / sd(dm)` over the pairs, which can be
#'       compared between terms in different units. `NA` for the intercept.
#'     * `t`: The pseudo-t statistic, `estimate / se(estimate)`.
#'     * `p_value`: The permutation p-value of `t`, two-sided.
#'       `ecodist::MRM()`'s `pval`.
#'   * `r_squared`: \eqn{R^2}, the proportion of the variance of `dm`'s
#'     distances that the terms explain.
#'   * `F`: The \eqn{F} statistic of the whole model.
#'   * `p_value`: The permutation p-value of the whole model (of \eqn{F}, and
#'     equally of \eqn{R^2}).
#'   * `df`: The degrees of freedom of \eqn{F}: the number of terms, and the
#'     number of pairs less the number of terms less 1. They are not those
#'     of the test, which is by permutation: the pairs aren't independent.
#'   * `n_perm`: The number of permutations.
#'   * `method`: `'pearson'` or `'spearman'`.
#'   * `n`: The number of samples.
#'   * `permutation`: How samples were permuted: `'free'`, `'within strata'`,
#'     or `'supplied matrix'`.
#'
#'   With `backward = TRUE`, these describe the final model, and there are
#'   also:
#'
#'   * `alpha`: The significance level.
#'   * `elimination`: A data frame with one row per term removed, in order
#'     (no rows if none was): `step`, `term`, its `t` and `p_value` in the
#'     model it was removed from, and that model's \eqn{R^2}
#'     (`r_squared_before`) and the next one's (`r_squared_after`).
#'
#'   Attribute `f_perm` is a vector of the \eqn{F} statistics of the
#'   permutations, and attribute `t_perm` a matrix of their pseudo-t
#'   statistics, one row per permutation and one column per coefficient.
#'   `as.data.frame()` gives `coefficients`, and [plot()][plot.ecodive_mantel_regression()]
#'   draws the standardized coefficients.
#'
#' @details
#'
#' **Method**
#'
#' MRM (Smouse et al. 1986; Legendre et al. 1994; Lichstein 2007) unfolds
#' each matrix into the vector of its \eqn{N = n(n-1)/2} distances between
#' pairs of samples, and fits the linear regression of `dm`'s distances on
#' the terms', with an intercept, by least squares. The pairs aren't
#' independent of one another, so the p-values come from permutations of
#' the samples: the rows and columns of `dm` are shuffled together, the
#' explanatory distances are held fixed, and the regression is refitted
#' (Legendre et al. 1994). As in `ecodist::MRM()`:
#'
#' * Each coefficient is tested with its pseudo-t statistic (Legendre et al.
#'   1994), the coefficient divided by its standard error, rather than with
#'   the coefficient itself, two-sided: its p-value is the proportion of
#'   permutations, counting the observed one, with \eqn{|t|} at least as
#'   large.
#' * The whole model is tested with \eqn{F} (equivalently, with \eqn{R^2}):
#'   the proportion of permutations, counting the observed one, with an
#'   \eqn{F} at least as large.
#'
#' With `method = 'spearman'`, each set of distances is replaced by its
#' ranks (average ranks for ties) before the regression: `ecodist::MRM()`'s
#' `mrank = TRUE`.
#'
#' With a single term, \eqn{R^2} is the square of [mantel()]'s statistic,
#' and the test of \eqn{F} is the two-sided Mantel test.
#'
#' Terms whose distances are collinear (one a linear combination of the
#' others, or the same for every pair) are an error.
#'
#'
#' **Backward elimination**
#'
#' With `backward = TRUE`, the model is simplified by backward elimination,
#' as in Lichstein (2007): the term with the largest p-value is removed if that p-value
#' is above `alpha`, the model is refitted without it, and so on, until
#' every term left has a p-value of at most `alpha`. Among terms with equal
#' p-values (permutation p-values often tie), the one with the smallest
#' \eqn{|t|} goes first. Every refit uses the same permutations, so the
#' steps, like the p-values, are the same from run to run and for any number
#' of `cpus`. If every term is removed, only the intercept is left: then
#' \eqn{R^2} is 0, and \eqn{F} and its p-value are `NA`. The intercept is
#' never removed.
#'
#' The final model's p-values don't account for the selection, so they are
#' too small, as after any stepwise selection: report them as descriptive.
#' `ecodist::MRM()` has no backward elimination.
#'
#'
#' **Speed and memory**
#'
#' The explanatory distances are centred and orthonormalized once (by QR).
#' Every statistic of a permutation then follows from their \eqn{m} products
#' with `dm`'s permuted distances, computed in one pass over the pairs,
#' reading `dm` through the permutation: \eqn{O(n^2 m)} per permutation, for
#' \eqn{m} terms, with no \eqn{n \times n} matrix and no refit. Permutations
#' run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. P-values are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `ecodist::MRM()`'s,
#' so p-values agree with it only up to sampling error, unless the same
#' permutations are given to both.
#'
#'
#' **Differences from ecodist::MRM()**
#'
#' * The explanatory distances can be given as a formula of sample
#'   metadata, as one value per sample, or as coordinates, as well as
#'   distance matrices, and samples are matched by name between them all.
#' * `strata` restricts the permutations.
#' * `permutations` counts the permutations besides the observed samples;
#'   `ecodist::MRM()`'s `nperm` counts the observed samples as one of them.
#'   So `permutations = 999` matches `nperm = 1000`.
#' * There's no `method = 'logistic'`, for logistic regression of a binary
#'   response.
#' * `backward` adds backward elimination of terms.
#'
#' @references
#' Legendre, P., Lapointe, F.-J., & Casgrain, P. (1994). Modeling brain
#' evolution from behavior: a permutational regression approach.
#' *Evolution*, 48(5), 1487-1499. \doi{10.1111/j.1558-5646.1994.tb02191.x}
#'
#' Lichstein, J. W. (2007). Multiple regression on distance matrices: a
#' multivariate spatial analysis tool. *Plant Ecology*, 188(2), 117-131.
#' \doi{10.1007/s11258-006-9126-3}
#'
#' Smouse, P. E., Long, J. C., & Sokal, R. R. (1986). Multiple regression
#' and correlation extensions of the Mantel test of matrix correspondence.
#' *Systematic Zoology*, 35(4), 627-632. \doi{10.2307/2413122}
#'
#' @seealso [mantel()], [mantel_partial()], [mantel_correlogram()],
#'   [varpart()], [dbrda()]
#'
#' @export
#' @examples
#'     # 30 samples at random locations, along a pH gradient that shifts
#'     # the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       ph   = round(runif(30, 6, 8), 1),
#'       temp = round(runif(30, 10, 20)),
#'       soil = sample(c('clay', 'sand'), 30, replace = TRUE),
#'       lat  = runif(30, 30, 35),
#'       lon  = runif(30, -100, -95) )
#'     counts[, 1:5] <- counts[, 1:5] + round(4 * (meta$ph - 6))
#'
#'     dm <- bray(counts)
#'
#'     # Community differences explained by differences in pH and
#'     # temperature, soil type (same or different), and geographic distance.
#'     res <- mantel_regression(dm ~ ph + temp + soil + cbind(lat, lon), data = meta)
#'     res
#'     plot(res)
#'
#'     # The same with distance matrices, and with ranks.
#'     geo <- dist(meta[, c('lat', 'lon')])
#'     mantel_regression(dm, ph = dist(meta$ph), geo = geo, method = 'spearman')
#'
#'     # Keep only the terms that matter, by backward elimination.
#'     mantel_regression(
#'       dm ~ ph + temp + soil + cbind(lat, lon), data = meta, backward = TRUE )
#'
mantel_regression <- function (
    dm, ..., data = NULL, formula = NULL, method = 'pearson',
    permutations = 999L, backward = FALSE, alpha = 0.05, strata = NULL,
    seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  if (!isTRUE(backward) && !isFALSE(backward))
    stop('`backward` must be TRUE or FALSE.', call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1)
    stop('`alpha` must be a single number from 0 to 1.', call. = FALSE)

  method <- match.arg(method, c('pearson', 'spearman'))

  dots  <- list(...)
  exprs <- as.list(substitute(list(...)))[-1]

  if (any(vapply(dots, function (v) inherits(v, 'formula'), logical(1))) ||
      inherits(data, 'formula'))
    stop(
      'Name the formula argument: mantel_regression(dm, data = meta, ',
      'formula = ~ a + b), or give mantel_regression(dm ~ a + b, data = meta).',
      call. = FALSE )

  fd      <- formula_dm(dm, formula)
  dm      <- fd$dm
  formula <- fd$formula

  x <- mantel_dist(dm, 'dm')
  n <- x$n
  if (n < 3)
    stop('`dm` must have at least 3 samples.', call. = FALSE)
  labels <- x$labels


  #--------------------------------------------------------
  # The terms: formula's, then those in `...`, each as
  # packed distances in dm's sample order.
  #--------------------------------------------------------
  term_labels <- character(0)
  if (!is.null(formula)) {
    if (is.null(data))
      stop('`data` is required with `formula`: the sample metadata.', call. = FALSE)
    if (is.character(formula)) {
      f <- tryCatch(str2lang(formula), error = function (e) NULL)
      if (!is.call(f) || !identical(f[[1]], as.name('~')))
        stop('`formula` must be a formula, such as `~ ph + temp`.', call. = FALSE)
      formula <- eval(f)
      environment(formula) <- parent.frame()
    }
    term_labels <- mrm_term_labels(formula)
  }

  if (is.null(data)) {
    des    <- NULL
    strata <- strata_codes(per_sample(strata, 'strata', NULL, NULL, labels, n))
  }
  else {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    des <- model_design(
      data, formula = if (is.null(formula)) ~ 1 else formula, labels = labels,
      n = n, strata = strata, allow_empty = TRUE, env = parent.frame() )
    strata <- des$strata
  }

  terms <- list()

  for (lab in term_labels) {
    v <- tryCatch(
      eval(str2lang(lab), des$data, environment(formula)),
      error = function (e) stop('Cannot evaluate term `', lab, '`: ', conditionMessage(e), call. = FALSE) )
    terms[[length(terms) + 1L]] <- mrm_values_dist(v, lab, n)
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) dot_names <- rep('', length(dots))
  for (i in seq_along(dots)) {
    v  <- dots[[i]]
    nm <- dot_names[[i]]
    if (!nzchar(nm))
      nm <- if (is.character(v) && length(v) == 1 && !is.na(v)) v
            else paste(deparse(exprs[[i]], width.cutoff = 500L), collapse = ' ')
    terms[[length(terms) + 1L]] <- mrm_dot_dist(v, nm, des, data, labels, n)
    term_labels <- c(term_labels, nm)
  }

  p <- length(terms)
  if (p == 0)
    stop(
      'Give at least one explanatory term: distances in `...`, or a ',
      '`formula` with `data`.', call. = FALSE )

  if (anyDuplicated(term_labels))
    stop('Two terms are named ', name_list(term_labels[duplicated(term_labels)][[1]]),
         '. Name each term in `...` differently.', call. = FALSE)

  N  <- length(x$dm)
  df <- N - p - 1
  if (df < 1)
    stop(
      'There are too few pairs of samples (', N, ') for ', p, ' term',
      if (p > 1) 's', ': at least ', p + 2, ' are needed.', call. = FALSE )

  pm <- perm_matrix(permutations, n)
  if (!is.null(pm$perms) && !is.null(strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)
  if (backward && pm$n_perm == 0)
    stop('Backward elimination needs permutations: its decisions use the ',
         'permutation p-values.', call. = FALSE)


  #--------------------------------------------------------
  # The response and the terms' distances (their ranks, for
  # Spearman); the terms centred and orthonormalized.
  #--------------------------------------------------------
  y <- x$dm
  X <- matrix(unlist(terms, use.names = FALSE), N, p)
  if (method == 'spearman') {
    y <- rank(y)
    for (j in seq_len(p)) X[, j] <- rank(X[, j])
  }

  if (all(y == y[[1]]))
    stop('`dm` has the same distance between every pair of samples, so ',
         'there is nothing to explain.', call. = FALSE)

  qrx <- qr(X - rep(colMeans(X), each = N))
  if (qrx$rank < p) {
    bad <- term_labels[qrx$pivot[seq.int(qrx$rank + 1L, p)]]
    stop(
      'Collinear terms: ', name_list(bad), ' ', if (length(bad) > 1) 'are' else 'is',
      ' a linear combination of the other terms (or the same for every ',
      'pair of samples).', call. = FALSE )
  }


  #--------------------------------------------------------
  # The full model; then, for backward elimination, the
  # least significant term removed and the rest refitted
  # (with the same permutations), until every term left
  # has p <= alpha.
  #--------------------------------------------------------
  fit <- function (keep) mrm_fit(
    y, X[, keep, drop = FALSE], term_labels[keep], pm, strata, n, seed, cpus )

  keep <- seq_len(p)
  res  <- fit(keep)
  elim <- NULL

  if (backward) {
    repeat {
      cf <- res$coefficients[-1, , drop = FALSE]
      if (nrow(cf) == 0 || all(cf$p_value <= alpha)) break

      # The largest p-value; among ties, the smallest |t|, then the last.
      worst <- order(-cf$p_value, abs(cf$t), -seq_len(nrow(cf)))[[1]]
      r2    <- res$stats$r_squared
      keep  <- keep[-worst]
      res   <- fit(keep)

      elim <- rbind(elim, data.frame(
        step             = NROW(elim) + 1L,
        term             = cf$term[[worst]],
        t                = cf$t[[worst]],
        p_value          = cf$p_value[[worst]],
        r_squared_before = r2,
        r_squared_after  = res$stats$r_squared,
        stringsAsFactors = FALSE ))
    }
    if (is.null(elim))
      elim <- data.frame(
        step = integer(0), term = character(0), t = numeric(0),
        p_value = numeric(0), r_squared_before = numeric(0),
        r_squared_after = numeric(0), stringsAsFactors = FALSE )
  }

  out <- c(
    res$stats,
    list(
      n_perm      = pm$n_perm,
      method      = method,
      n           = as.integer(n),
      permutation = if (!is.null(pm$perms)) 'supplied matrix'
                    else if (!is.null(strata)) 'within strata'
                    else 'free' ),
    if (backward) list(alpha = alpha, elimination = elim) )
  out <- c(list(coefficients = res$coefficients), out)

  attr(out, 'f_perm') <- res$f_perm
  attr(out, 't_perm') <- res$t_perm
  class(out) <- 'ecodive_mantel_regression'

  set_display(out, title = title_of('Multiple Regression on Distance Matrices', x$title))
}



#------------------------------------------------------------------------------
# The regression of y on the columns of X (N pairs x p terms, of full rank
# after centring; p may be 0), its statistics for the observed samples, and
# its permutation p-values.
#
# Returns list(coefficients, stats = list(r_squared, F, p_value, df),
# f_perm, t_perm).
#------------------------------------------------------------------------------
mrm_fit <- function (y, X, term_labels, pm, strata, n, seed, cpus) {

  N  <- length(y)
  p  <- ncol(X)
  df <- N - p - 1

  ybar <- mean(y)
  yc   <- y - ybar
  tss  <- sum(yc * yc)
  sdy  <- sqrt(tss / (N - 1))

  # No terms (all eliminated): the intercept is the mean, the same for
  # every permutation, so its pseudo-t is too.
  if (p == 0) {
    t0     <- ybar / sqrt(tss / df / N)
    t_perm <- matrix(t0, pm$n_perm, 1L, dimnames = list(NULL, '(Intercept)'))
    coefs  <- data.frame(
      term = '(Intercept)', estimate = ybar, std_estimate = NA_real_, t = t0,
      p_value = perm_p(t0, t_perm, 'two.sided'), stringsAsFactors = FALSE )
    return (list(
      coefficients = coefs,
      stats  = list(r_squared = 0, F = NA_real_, p_value = NA_real_, df = c(0, df)),
      f_perm = rep(NA_real_, pm$n_perm),
      t_perm = t_perm ))
  }

  xbar <- colMeans(X)
  Xc   <- X - rep(xbar, each = N)
  qrx  <- qr(Xc)
  Q    <- qr.Q(qrx)
  Rinv <- backsolve(qr.R(qrx), diag(p))  # (Xc'Xc)^-1 = Rinv Rinv'


  #--------------------------------------------------------
  # Q'y for the observed samples and each permutation;
  # every statistic follows from it.
  #--------------------------------------------------------
  U <- .Call(
    C_mantel_regression, yc, Q, pm$perms, strata,
    as.integer(c(n, p, pm$n_perm, seed, cpus)) )

  ssr <- colSums(U * U)
  rss <- pmax(tss - ssr, 0)
  # The observed fit's residuals, exactly.
  rss[[1]] <- sum((yc - Q %*% U[, 1])^2)

  B  <- Rinv %*% U                     # Slopes, p x (n_perm + 1).
  b0 <- ybar - colSums(B * xbar)       # Intercepts.
  s2 <- rss / df
  vb  <- rowSums(Rinv * Rinv)          # diag((Xc'Xc)^-1)
  vb0 <- 1 / N + sum(crossprod(Rinv, xbar)^2)

  t_all <- rbind(b0 / sqrt(s2 * vb0), B / sqrt(outer(vb, s2)))
  f_all <- (ssr / p) / s2

  t_obs  <- t_all[, 1]
  f_obs  <- f_all[[1]]
  t_perm <- t(t_all[, -1, drop = FALSE])
  f_perm <- f_all[-1]
  colnames(t_perm) <- c('(Intercept)', term_labels)

  sdx <- sqrt(colSums(Xc * Xc) / (N - 1))

  coefs <- data.frame(
    term         = c('(Intercept)', term_labels),
    estimate     = c(b0[[1]], B[, 1]),
    std_estimate = c(NA_real_, B[, 1] * sdx / sdy),
    t            = t_obs,
    p_value      = perm_p(t_obs, t_perm, 'two.sided'),
    stringsAsFactors = FALSE )
  rownames(coefs) <- NULL

  list(
    coefficients = coefs,
    stats = list(
      r_squared = ssr[[1]] / tss,
      F         = f_obs,
      p_value   = perm_p(f_obs, f_perm, 'greater'),
      df        = c(p, df) ),
    f_perm = f_perm,
    t_perm = t_perm )
}



#------------------------------------------------------------------------------
# The term labels of `formula`: no Condition() terms, and no interactions.
#------------------------------------------------------------------------------
mrm_term_labels <- function (formula) {

  if (!inherits(formula, 'formula'))
    stop('`formula` must be a formula, such as `~ ph + temp`.', call. = FALSE)

  tt <- tryCatch(
    stats::terms(formula, specials = 'Condition'),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  if (length(attr(tt, 'specials')$Condition))
    stop('`formula` cannot have Condition() terms: every term of a ',
         'regression is fitted together.', call. = FALSE)

  if (any(attr(tt, 'order') > 1))
    stop('`formula` cannot have interactions.', call. = FALSE)

  attr(tt, 'term.labels')
}



#------------------------------------------------------------------------------
# Packed distances (in dm's sample order) of `v`, one value per sample,
# already in dm's sample order: |v_i - v_j| for numbers, Euclidean distances
# for a numeric matrix, and 0 (same) or 1 (different) for a factor,
# character, or logical vector. `arg` names it in messages.
#------------------------------------------------------------------------------
mrm_values_dist <- function (v, arg, n) {

  if (is.data.frame(v)) {
    if (!all(vapply(v, is.numeric, logical(1))))
      stop('Term ', name_list(arg), ' must be numeric to be coordinates.', call. = FALSE)
    v <- as.matrix(v)
  }

  if (is.matrix(v) && ncol(v) == 1) v <- v[, 1]

  if (is.matrix(v)) {
    if (!is.numeric(v))
      stop('Term ', name_list(arg), ' must be numeric to be coordinates.', call. = FALSE)
    if (nrow(v) != n)
      stop('Term ', name_list(arg), ' has ', nrow(v), ' rows, but `dm` has ', n,
           ' samples.', call. = FALSE)
    if (!all(is.finite(v)))
      stop('Term ', name_list(arg), ' must not contain NA, NaN, or infinite values.', call. = FALSE)
    storage.mode(v) <- 'double'
    d <- as.vector(stats::dist(unname(v)))
    if (!all(is.finite(d)))
      stop('The distances of term ', name_list(arg), ' are too large.', call. = FALSE)
    return (d)
  }

  if (!is.atomic(v) || length(v) != n)
    stop('Term ', name_list(arg), ' must have one value per sample (', n, ').', call. = FALSE)

  j <- rep.int(seq_len(n - 1L), (n - 1L):1)
  i <- j + sequence((n - 1L):1)

  if (is.factor(v) || is.character(v) || is.logical(v)) {
    if (anyNA(v))
      stop('Term ', name_list(arg), ' must not contain NAs.', call. = FALSE)
    v <- as.character(v)
    return (as.double(v[i] != v[j]))
  }

  if (!is.numeric(v))
    stop('Term ', name_list(arg), ' must be numeric, a factor, character, ',
         'or logical.', call. = FALSE)

  v <- unname(as.double(v))
  if (!all(is.finite(v)))
    stop('Term ', name_list(arg), ' must not contain NA, NaN, or infinite values.', call. = FALSE)

  d <- abs(v[i] - v[j])
  if (!all(is.finite(d)))
    stop('The differences between values of term ', name_list(arg), ' are too large.', call. = FALSE)
  d
}



#------------------------------------------------------------------------------
# Packed distances (in dm's sample order) of an argument in `...`: a dist
# object or distance matrix; coordinates (a data frame or a matrix that
# isn't one of distances); a column name of `data`; or one value per sample.
#------------------------------------------------------------------------------
mrm_dot_dist <- function (v, arg, des, data, labels, n) {

  # Distances.
  if (inherits(v, 'dist') || (is.matrix(v) && mantel_is_dist(v)))
    return (mantel_align(mantel_dist(v, arg), arg, labels, n))

  # Coordinates, one row per sample: matched by row name.
  if (is.data.frame(v) || is.matrix(v)) {
    rn  <- if (.row_names_info(v) > 0 || (is.matrix(v) && !is.null(rownames(v)))) rownames(v) else NULL
    idx <- match_rows(labels, n, rn, nrow(v), what = paste0('`', arg, '`'))
    return (mrm_values_dist(v[idx, , drop = FALSE], arg, n))
  }

  if (!is.atomic(v) || is.null(v))
    stop('`', arg, '` must be distances, coordinates, one value per sample, ',
         'or a column name of `data`.', call. = FALSE)

  df <- if (is.null(data)) NULL else as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.null(df) && is.character(v) && length(v) == 1 && !is.na(v) && !v %in% names(df))
    stop('`', arg, '` is not a column of `data`: ', name_list(v), call. = FALSE)

  vals <- per_sample(v, arg, df, des$idx, labels, n)
  mrm_values_dist(vals, arg, n)
}



#------------------------------------------------------------------------------
#' @export
print.ecodive_mantel_regression <- function (x, digits = max(3L, getOption('digits') - 3L), ...) {

  fmt <- function (v) if (is.na(v)) 'NA' else trimws(formatC(v, digits = digits, format = 'g'))

  title <- display_title(x)
  if (is.null(title)) title <- 'Multiple Regression on Distance Matrices'
  cat(title, '\n\n', sep = '')

  stat <- paste0('R2 = ', fmt(x$r_squared), ', F = ', fmt(x$F))
  if (x$n_perm > 0)
    stat <- paste0(stat, ', p = ', if (is.na(x$p_value)) 'NA' else format.pval(x$p_value, digits = digits))
  cat(stat, '\n', sep = '')
  cat('Permutations: ', x$n_perm, if (!is.null(x$permutation)) paste0(' (', x$permutation, ')'), '\n', sep = '')
  cat('Regression: ', switch(x$method,
        pearson  = 'least squares on the distances',
        spearman = 'least squares on the ranks of the distances' ), '\n', sep = '')
  cat('Samples: ', x$n, ' (', x$n * (x$n - 1) / 2, ' pairs)\n\n', sep = '')

  cf <- x$coefficients
  for (col in c('estimate', 'std_estimate', 't'))
    cf[[col]] <- ifelse(is.na(cf[[col]]), '', formatC(cf[[col]], digits = digits, format = 'g'))
  cf$p_value <- if (x$n_perm > 0) format.pval(cf$p_value, digits = digits) else NA
  print(cf, row.names = FALSE, right = TRUE)

  el <- x$elimination
  if (is.data.frame(el)) {
    cat('\nBackward elimination (alpha = ', fmt(x$alpha), '): ', sep = '')
    if (nrow(el) == 0) {
      cat('no term removed\n')
    }
    else {
      cat(nrow(el), if (nrow(el) == 1) ' term' else ' terms', ' removed\n', sep = '')
      el$t                <- formatC(el$t, digits = digits, format = 'g')
      el$p_value          <- format.pval(el$p_value, digits = digits)
      el$r_squared_before <- formatC(el$r_squared_before, digits = digits, format = 'g')
      el$r_squared_after  <- formatC(el$r_squared_after, digits = digits, format = 'g')
      print(el, row.names = FALSE, right = TRUE)
    }
  }

  invisible(x)
}



#------------------------------------------------------------------------------
#' @export
as.data.frame.ecodive_mantel_regression <- function (x, row.names = NULL, optional = FALSE, ...) {
  df <- x$coefficients
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}



#------------------------------------------------------------------------------
#' Plot a multiple regression on distance matrices
#'
#' Draws each term's standardized coefficient (`std_estimate`), which can be
#' compared between terms in different units: filled for the terms whose
#' p-value is at most `alpha`, open for the others, with a dashed line at 0.
#' The subtitle gives the model's \eqn{R^2} and p-value.
#'
#' @param x  A [mantel_regression()] result.
#'
#' @param alpha  The significance level of the filled points, compared with
#'        each term's `p_value`. Default: `0.05`.
#'
#' @param ...  Passed on to [graphics::dotchart()], such as `main`, `xlab`,
#'        `col`, `xlim`, or `sub`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     meta <- data.frame(ph = runif(20, 6, 8), temp = runif(20, 10, 20))
#'     counts <- matrix(rpois(20 * 10, 5), 20, 10)
#'     counts[, 1:3] <- counts[, 1:3] + round(4 * (meta$ph - 6))
#'
#'     res <- mantel_regression(bray(counts), data = meta, formula = ~ ph + temp)
#'     plot(res)
#'
plot.ecodive_mantel_regression <- function (x, alpha = 0.05, ...) {

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1)
    stop('`alpha` must be a single number from 0 to 1.', call. = FALSE)

  cf <- x$coefficients
  if (!is.data.frame(cf) || !all(c('term', 'std_estimate', 'p_value') %in% names(cf)))
    stop('`x` must be a mantel_regression() result.', call. = FALSE)

  cf  <- cf[cf$term != '(Intercept)', , drop = FALSE]
  if (nrow(cf) == 0)
    stop('`x` has no terms to plot: backward elimination removed them all.', call. = FALSE)
  idx <- rev(seq_len(nrow(cf)))  # dotchart() draws from the bottom up.
  est <- cf$std_estimate[idx]
  sig <- !is.na(cf$p_value[idx]) & cf$p_value[idx] <= alpha

  args <- list(...)
  col  <- plot_arg(args, 'col', 'black')
  sub  <- sprintf(
    'R2 = %s%s', plot_num(x$r_squared),
    if (isTRUE(x$n_perm > 0)) sprintf(', p = %s', plot_p(x$p_value)) else '' )

  plot_do(
    graphics::dotchart, c(list(est), plot_dot_args(plot_drop(args, 'col'))),
    list(
      labels = cf$term[idx], pch = ifelse(sig, 19, 1), color = col,
      xlim   = range(0, est, finite = TRUE),
      main   = plot_title(x, 'Multiple Regression on Distance Matrices'),
      xlab   = if (identical(x$method, 'spearman'))
                 'Standardized coefficient (ranks)' else 'Standardized coefficient' ))

  graphics::title(sub = plot_arg(args, 'sub', sub))
  graphics::abline(v = 0, lty = 2, col = PLOT_MUTED)

  invisible(x)
}
