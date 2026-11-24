# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# The p-values follow scipy.stats.pearsonr(), spearmanr() and
# kendalltau(method = 'asymptotic') (SciPy 1.18, BSD-3; see
# src/cross_assoc.c and inst/COPYRIGHTS).


# Display names of cross_assoc()'s methods, keyed by its `method` attribute.
CROSS_ASSOCS <- c(
  spearman = 'Spearman Correlation',
  pearson  = 'Pearson Correlation',
  kendall  = "Kendall's tau-b" )

CA_METHOD_CODES <- c(pearson = 0L, spearman = 1L, kendall = 2L)



#' Associations Between Two Tables
#'
#' Correlates every feature of one table (taxa) with every variable of
#' another measured on the same samples (environmental factors, metabolites,
#' clinical variables): the correlation, its p-value and a p-value adjusted
#' for testing every pair. `plot()` draws the familiar heatmap with stars on
#' the significant cells.
#'
#' @inheritParams documentation
#'
#' @param x  A count table (samples \eqn{\times} features), or any of the
#'        objects listed under **Input Types**, or any numeric table
#'        (metabolites, host gene expression) with `norm = 'none'`.
#'
#' @param y  The second table, with **samples in rows**: a numeric matrix,
#'        a data frame of numeric or logical columns, a numeric vector (one
#'        variable), or a count table. Or, with `data`, the names of columns
#'        of `data`. Logical columns become 0/1; factor and character columns
#'        are an error. May contain `NA`s (see **Missing values**). Samples
#'        are matched to `x`'s by name when both have names, otherwise taken
#'        in the same order.
#'
#' @param data  A data frame of sample metadata, with one row per sample,
#'        for `y`, `given` and `strata` given as column names (or `given` as
#'        a formula). Rows are matched to `x`'s samples by name when it has
#'        row names, otherwise taken in the same order. Default: `NULL`.
#'
#' @param method  The correlation coefficient: `'spearman'` (the default),
#'        `'pearson'`, or `'kendall'` (tau-b).
#'
#' @param norm  How `x` is transformed first: any of [normalize()]'s norms
#'        (`'percent'`, `'clr'`, `'log'`, `'rclr'`, `'hellinger'`, ...,
#'        `'none'`). Default: `'percent'` (relative abundances).
#'
#' @param y_norm  How `y` is transformed first, as `norm`. Default: `'none'`.
#'        Any other choice needs `y` to be a table of non-negative values
#'        without `NA`s, such as a second count table.
#'
#' @param pseudocount  For `norm` or `y_norm` `'clr'` or `'log'`: the value
#'        added to every count first. See **Pseudocount**.
#'
#' @param given  Covariates to partial out, for partial correlations: a
#'        one-sided formula over the columns of `data` (such as `~ age +
#'        sex`), the names of columns of `data`, or a data frame or matrix
#'        with one row per sample. Not available for Kendall. Samples with
#'        `NA` in `given` are left out of everything, with a warning.
#'        Default: `NULL`.
#'
#' @param permutations  `0` (the default) for p-values from the
#'        coefficient's distribution (see **P-values**); a number of
#'        permutations for permutation p-values; or a matrix with one
#'        permutation of `1:n` per row (`n` the samples), to use instead of
#'        random ones.
#'
#' @param seed  An integer seed for the random permutations. Default: `0`.
#'
#' @param alpha  The significance level of the adjusted p-values: `reject`
#'        is `p_adjusted <= alpha`. Default: `0.05`.
#'
#' @param p_adjust  How to adjust the p-values for testing every pair: any
#'        method of [stats::p.adjust()], over all the pairs at once. Default:
#'        `'BH'` (Benjamini-Hochberg false discovery rate).
#'
#' @param min_samples  Pairs with fewer samples (after leaving out `NA`s)
#'        get `NA`. At least 3 (the default), more with `given`.
#'
#' @return A data frame of class `c("ecodive_cross_assoc", "data.frame")`
#'   with one row per pair, in `x`-major order (all of `y` for the first
#'   feature of `x`, then the next), and columns:
#'
#'   * `x`, `y`: the feature and the variable (their input positions, as
#'     text, when unnamed);
#'   * `association`: the correlation;
#'   * `statistic`: Student's t (Pearson, Spearman) or the normal z (Kendall)
#'     behind the p-value;
#'   * `p_value`, `p_adjusted` and `reject` (`p_adjusted <= alpha`);
#'   * `n`: the samples used, when this differs between pairs (with `NA`s in
#'     `y`).
#'
#'   Attributes: `method` (`"spearman"`, `"pearson"` or `"kendall"`),
#'   `alpha`, `p_adjust`, `n_samples`, `given` (the covariates' names, with
#'   `given`), `n_perm` (with permutations), and `x_features`, `y_features`
#'   (each table's names in input order).
#'
#'   `as.matrix(res)` gives the \eqn{p \times q} matrix of correlations
#'   (features of `x` in rows), and `as.matrix(res, 'p_adjusted')` that of
#'   any other column. `plot(res)` draws it as a heatmap: see
#'   [plot.ecodive_cross_assoc()].
#'
#' @details
#' Pearson's \eqn{r} correlates the (transformed) values; Spearman's
#' \eqn{\rho} is Pearson's \eqn{r} of their ranks (ties given their average
#' rank); Kendall's \eqn{\tau_b} counts concordant and discordant pairs of
#' samples, corrected for ties. All three are the values of
#' [stats::cor()].
#'
#' **Which transformation.** Spearman and Kendall depend only on the order
#' of each feature's values across samples, so transformations of one
#' feature at a time (logs, ranks) don't change them. Transformations of one
#' sample at a time do: relative abundances (`norm = 'percent'`, the
#' default, and what most papers correlate) rank the samples differently
#' from raw counts. Relative abundances are compositional, so a taxon can
#' appear to correlate with a variable only because other taxa do (closure;
#' see [vlr()]); `norm = 'clr'` (with `pseudocount`) is the safer choice,
#' with any method. With many zeros, Spearman's \eqn{\rho} on relative
#' abundances is dominated by presence and absence: every zero gets the same
#' rank.
#'
#' **Partial correlations** (`given`): both tables' values (their ranks,
#' for Spearman, and the covariates' ranks too) are residualised on the
#' covariates (with an intercept) by least squares, and the result is
#' Pearson's \eqn{r} of the residuals, with \eqn{n - 2 - k} degrees of
#' freedom for \eqn{k} covariates. This is `ppcor::pcor.test()`'s
#' definition.
#'
#' **Speed.** The correlations are computed in C, in parallel over the
#' features of `x`, without forming the joint \eqn{(p + q) \times (p + q)}
#' matrix that [stats::cor()] of the two tables bound together computes.
#' Sparse tables stay sparse: the sums run over each feature's non-zero
#' samples only.
#'
#' @section P-values:
#'   With `permutations = 0`, Pearson's and Spearman's p-values are
#'   two-sided, from Student's \eqn{t = r \sqrt{(n - 2 - k) / (1 - r^2)}}
#'   with \eqn{n - 2 - k} degrees of freedom (\eqn{k = 0} without `given`),
#'   as `scipy.stats.pearsonr()` and `spearmanr()`, `Hmisc::rcorr()` and
#'   `psych::corr.test()` compute them. \eqn{|r| = 1} gives \eqn{p = 0}.
#'   Kendall's is the normal approximation with the tie-corrected variance
#'   (Kendall 1970), as `scipy.stats.kendalltau(method = 'asymptotic')`.
#'
#'   With `permutations`, the rows of `y` are permuted, the same permutation
#'   for every pair (which keeps the correlations within each table), and
#'   \eqn{p = (1 + \#\{|r_b| \ge |r|\}) / (B + 1)} for \eqn{B}
#'   permutations. With `given`, the residuals of `y` are permuted (ter
#'   Braak 1992). With `NA`s in `y`, the variables with the same samples
#'   missing share their permutations, of their own samples; a
#'   `permutations` matrix then can't be used.
#'
#'   The p-values are then adjusted over **all** \eqn{p \times q} pairs at
#'   once (Benjamini & Hochberg 1995 by default). None of these p-values
#'   allows for the correlations between taxa, or between variables, which
#'   make the pairs' tests far from independent.
#'
#' @section Missing values:
#'   `x` can't contain `NA`s. Each pair uses the samples where its `y`
#'   variable is not `NA` (pairwise complete observations, as
#'   `Hmisc::rcorr()` and `cor(use = 'pairwise.complete.obs')`), and the
#'   ranks are taken over those samples. Pairs with fewer than
#'   `min_samples` samples, and pairs with a feature or variable that is
#'   constant over their samples, get `NA` (the latter with a warning, as
#'   [stats::cor()]).
#'
#' @section Differences from cor.test(), Hmisc, psych and microeco:
#'   * [stats::cor.test()] gives Spearman's p-value from the exact
#'     distribution of \eqn{\rho} (AS 89) for fewer than 1290 samples without
#'     ties, and Kendall's from its exact distribution for fewer than 50
#'     samples without ties; `cor.test(exact = FALSE)` gives the
#'     approximations used here (the same values, ties included).
#'   * `Hmisc::rcorr()` correlates the two tables bound together, the whole
#'     \eqn{(p + q) \times (p + q)} matrix, with the same p-values; it
#'     computes them as one minus a probability, so those below about
#'     \eqn{10^{-8}} are accurate only to about \eqn{10^{-16}} absolutely.
#'     It doesn't adjust them.
#'   * `psych::corr.test(x, y, adjust = 'BH')` gives the same values and
#'     p-values, adjusted over the \eqn{p \times q} pairs, as here.
#'   * microeco's `trans_env$cal_cor(method = 'spearman')` uses
#'     [stats::cor.test()], so its p-values are exact (above) for pairs
#'     without ties; by default (`p_adjust_type = 'All'`) it adjusts over
#'     all pairs too, and it can instead adjust within each taxon or each
#'     variable.
#'   * `ppcor::pcor.test()` gives the same partial correlations and
#'     p-values, one pair at a time (ranking the covariates too, for
#'     Spearman).
#'
#' @inheritSection pseudocount_section Pseudocount
#' @inheritSection documentation Input Types
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. *Journal of
#' the Royal Statistical Society B*, 57, 289–300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Kendall, M. G. (1970). *Rank Correlation Methods* (4th ed.). Griffin.
#'
#' Knight, W. R. (1966). A computer method for calculating Kendall's tau
#' with ungrouped data. *Journal of the American Statistical Association*,
#' 61, 436–439. \doi{10.1080/01621459.1966.10480879}
#'
#' Spearman, C. (1904). The proof and measurement of association between
#' two things. *American Journal of Psychology*, 15, 72–101.
#' \doi{10.2307/1412159}
#'
#' ter Braak, C. J. F. (1992). Permutation versus bootstrap significance
#' tests in multiple regression and ANOVA. In K.-H. Jöckel, G. Rothe, &
#' W. Sendler (Eds.), *Bootstrapping and Related Techniques* (pp. 79–86).
#' Springer.
#'
#' @seealso [plot.ecodive_cross_assoc()]; [feature_assoc()], [rho_p()] and
#'   [sparcc()] for associations between the features of one table;
#'   [envfit()] to fit variables onto an ordination.
#'
#' @export
#' @examples
#'     # Taxa against three simulated environmental variables
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 8, 20), 30, 8,
#'       dimnames = list(paste0('S', 1:30), paste0('Taxon', 1:8)))
#'     env <- data.frame(
#'       pH   = counts[, 1] / 10 + rnorm(30),
#'       temp = rnorm(30, 15),
#'       salt = c(rnorm(28), NA, NA),
#'       row.names = rownames(counts) )
#'
#'     res <- cross_assoc(counts, env)
#'     res
#'
#'     # The p x q matrix of correlations
#'     round(as.matrix(res), 2)
#'
#'     plot(res)
#'
cross_assoc <- function (
    x, y, data = NULL, method = 'spearman', margin = 1L,
    norm = 'percent', y_norm = 'none', pseudocount = NULL,
    given = NULL, permutations = 0L, strata = NULL,
    alpha = 0.05, p_adjust = 'BH', min_samples = 3L,
    seed = 0L, cpus = n_cpus() ) {

  if (!is.character(method) || length(method) != 1 || is.na(method) ||
      !tolower(method) %in% names(CA_METHOD_CODES))
    stop("`method` must be 'spearman', 'pearson', or 'kendall'.", call. = FALSE)
  method <- tolower(method)

  validate_alpha()
  validate_seed()
  validate_cpus()
  p_adjust <- pairwise_method(TRUE, p_adjust, TRUE)

  if (!is.numeric(min_samples) || length(min_samples) != 1 || is.na(min_samples) ||
      min_samples < 3 || min_samples %% 1 != 0)
    stop('`min_samples` must be a whole number of at least 3.', call. = FALSE)

  if (!is.null(data) && !is.data.frame(data))
    stop('`data` must be a data frame.', call. = FALSE)

  if (!is.null(given) && method == 'kendall')
    stop('`given` (partial correlation) is not available for Kendall; use ',
         "method = 'spearman' or 'pearson'.", call. = FALSE)


  #--------------------------------------------------
  # The tables, on x's samples.
  #--------------------------------------------------
  xt     <- cross_x(x, margin, norm, pseudocount, cpus)
  n      <- xt$n
  labels <- xt$samples

  idx <- if (is.null(data)) NULL else
    match_rows(labels, n, cross_row_names(data), nrow(data), ref = '`x`', what = '`data`')

  ym <- cross_y(y, data, idx, labels, n)
  if (!identical(tolower(y_norm), 'none')) {
    if (anyNA(ym))
      stop('`y` must not contain NAs when `y_norm` is not "none".', call. = FALSE)
    ym <- normalize(ym, y_norm, margin = 1L, pseudocount = pseudocount, cpus = cpus)
  }
  q <- ncol(ym)

  G <- cross_given(given, data, idx, labels, n)

  keep <- rep(TRUE, n)
  if (!is.null(G)) {
    keep <- stats::complete.cases(G)
    if (!all(keep))
      warning(sum(!keep), ' sample', if (sum(!keep) == 1) ' has' else 's have',
              ' NA in `given`, and ', if (sum(!keep) == 1) 'is' else 'are',
              ' left out.', call. = FALSE)
    G <- G[keep, , drop = FALSE]
  }
  n_used <- sum(keep)

  strata <- strata_codes(per_sample(strata, 'strata', data, idx, labels, n, ref = '`x`'))
  if (!is.null(strata)) strata <- strata[keep]

  pm <- perm_matrix(permutations, n_used)


  #--------------------------------------------------
  # y's columns, grouped by their missing samples.
  #--------------------------------------------------
  used <- which(keep)
  pats <- cross_assoc_patterns(
    ym[used, , drop = FALSE], used, G, method, min_samples, pm, strata, seed )

  res <- .Call(
    C_cross_assoc, xt$x, pats$c_pats, CA_METHOD_CODES[[method]], q,
    pm$n_perm, cpus )


  #--------------------------------------------------
  # The table.
  #--------------------------------------------------
  r <- res[[1]]
  r[is.nan(r)] <- NA_real_
  p <- xt$p

  # Constant y columns: r is 0 from C.
  col_of <- rep(seq_len(q), times = p)
  r[col_of %in% pats$const_y] <- NA_real_
  # Pairs in patterns with too few samples.
  r[col_of %in% pats$few] <- NA_real_

  ok_cols <- setdiff(seq_len(q), c(pats$const_y, pats$few))
  const_x <- if (length(ok_cols) == 0) integer(0) else
    which(colSums(!is.na(matrix(r, nrow = q)[ok_cols, , drop = FALSE])) == 0)
  n_const <- length(pats$const_y)
  if (length(const_x) || n_const)
    warning(
      'The correlations of ',
      paste(c(
        if (length(const_x)) sprintf('%d feature%s of `x`', length(const_x), if (length(const_x) == 1) '' else 's'),
        if (n_const) sprintf('%d column%s of `y`', n_const, if (n_const == 1) '' else 's') ),
        collapse = ' and ' ),
      ' are NA: constant over their samples (the standard deviation is zero).',
      call. = FALSE )

  n_pair <- pats$n[col_of]
  df     <- pats$df[col_of]

  if (method == 'kendall') {
    stat <- res[[2]]
    stat[is.na(r)] <- NA_real_
    pval <- 2 * stats::pnorm(-abs(stat))
  }
  else {
    stat <- r * sqrt(df / ((1 + r) * (1 - r)))
    stat[!is.na(r) & abs(r) == 1] <- sign(r[!is.na(r) & abs(r) == 1]) * Inf
    pval <- 2 * stats::pt(-abs(stat), df)
  }

  if (pm$n_perm > 0) {
    pval <- (res[[3]] + 1) / (pm$n_perm + 1)
    pval[is.na(r)] <- NA_real_
  }

  padj <- stats::p.adjust(pval, p_adjust)

  out <- data.frame(
    x           = rep(xt$features, each = q),
    y           = rep(colnames(ym), times = p),
    association = r,
    statistic   = stat,
    p_value     = pval,
    p_adjusted  = padj,
    reject      = padj <= alpha,
    stringsAsFactors = FALSE )
  if (length(unique(pats$n)) > 1) out$n <- n_pair

  attr(out, 'method')     <- method
  attr(out, 'alpha')      <- alpha
  attr(out, 'p_adjust')   <- p_adjust
  attr(out, 'n_samples')  <- n_used
  if (!is.null(G))      attr(out, 'given')  <- colnames(G)
  if (pm$n_perm > 0)    attr(out, 'n_perm') <- pm$n_perm
  attr(out, 'x_features') <- xt$features
  attr(out, 'y_features') <- colnames(ym)
  class(out) <- c('ecodive_cross_assoc', 'data.frame')
  out
}



#------------------------------------------------------------------------------
# Row names of a data frame or matrix, or NULL for automatic ones.
#------------------------------------------------------------------------------
cross_row_names <- function (d) {
  if (is.data.frame(d) && .row_names_info(d) < 0) return (NULL)
  rownames(d)
}



#------------------------------------------------------------------------------
# `x`, normalized, for C_cross_assoc(): list(x (the C list), n, p, samples,
# features).
#------------------------------------------------------------------------------
cross_x <- function (x, margin, norm, pseudocount, cpus) {

  counts <- x
  remove('x')
  if (is.data.frame(counts)) {
    ok <- vapply(counts, function (v) is.numeric(v) || is.logical(v), logical(1))
    if (!all(ok))
      stop('`x` column ', name_list(names(counts)[!ok]), ' is not numeric.', call. = FALSE)
  }

  tryCatch(validate_counts(), error = function (e)
    stop(sub('`counts`', '`x`', conditionMessage(e), fixed = TRUE), call. = FALSE) )

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  if (is.matrix(counts) && !is.numeric(counts))
    stop('`x` must be numeric.', call. = FALSE)

  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!all(is.finite(vals)))
    stop('`x` must not contain NA, NaN, or infinite values.', call. = FALSE)
  remove('vals')

  validate_margin()
  validate_norm()
  tryCatch(validate_pseudocount(), error = function (e)
    stop(sub('`counts`', '`x`', conditionMessage(e), fixed = TRUE), call. = FALSE) )

  pkg <- get_matrix_package(counts)
  switch(
    EXPR = pkg,
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  if (is.null(dn)) dn <- list(NULL, NULL)

  clr   <- identical(norm, NORM_CLR) || identical(norm, NORM_LOG)
  dense <- clr || pkg == 'base' || inherits(counts, 'dgeMatrix')

  res <- .Call(C_normalize, counts, margin, norm, as.double(pseudocount), dense, cpus)

  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }
  n <- as.integer(dims[[1]])
  p <- as.integer(dims[[2]])

  if (dense) {
    m <- res
    if (margin == 2L) m <- t(m)
    storage.mode(m) <- 'double'
    attributes(m) <- list(dim = c(n, p))
    cx <- list(n, p, m, NULL, NULL, NULL)
  }
  else {
    # Non-zero entries, column-major in counts' orientation.
    s <- if (margin == 1L) res[[1]] else res[[2]]
    f <- if (margin == 1L) res[[2]] else res[[1]]
    v <- res[[3]]
    keep <- v != 0
    s <- s[keep]; f <- f[keep]; v <- v[keep]
    if (margin == 2L) { o <- order(f, s); s <- s[o]; f <- f[o]; v <- v[o] }
    starts <- c(0L, cumsum(tabulate(f + 1L, nbins = p)))
    cx <- list(n, p, NULL, as.integer(starts), as.integer(s), as.double(v))
  }

  features <- dn[[2]]
  if (is.null(features)) features <- as.character(seq_len(p))

  list(x = cx, n = n, p = p, samples = dn[[1]], features = as.character(features))
}



#------------------------------------------------------------------------------
# `y` as a numeric n x q matrix, in x's sample order, with column names.
#------------------------------------------------------------------------------
cross_y <- function (y, data, idx, labels, n) {

  if (is.character(y) && is.null(dim(y))) {
    if (is.null(data))
      stop('`y` names columns, but `data` is not given.', call. = FALSE)
    miss <- setdiff(y, names(data))
    if (length(miss))
      stop('`y` names ', name_list(miss), ', not column', if (length(miss) > 1) 's',
           ' of `data`.', call. = FALSE)
    ym <- cross_numeric(data[idx, y, drop = FALSE], '`y`')
    return (ym)
  }

  if (is.null(y)) stop('`y` is required.', call. = FALSE)

  if (inherits(y, c('dgCMatrix', 'dgTMatrix', 'dgeMatrix')))
    y <- as.matrix(y)
  else if (inherits(y, 'simple_triplet_matrix'))
    y <- as.matrix(y)

  if (is.atomic(y) && is.null(dim(y))) {
    nm <- names(y)
    if (!is.numeric(y) && !is.logical(y))
      stop('`y` is not numeric; use da_test() for groups.', call. = FALSE)
    y <- matrix(as.double(y), ncol = 1, dimnames = list(nm, 'y'))
  }

  if (!is.data.frame(y) && !is.matrix(y))
    stop('`y` must be a numeric matrix or data frame, a numeric vector, ',
         'or the names of columns of `data`.', call. = FALSE)

  rn <- cross_row_names(y)
  if (!is.null(rn) && (anyNA(rn) || anyDuplicated(rn))) rn <- NULL
  ym <- cross_numeric(y, '`y`')
  rows <- match_rows(labels, n, rn, nrow(ym), ref = '`x`', what = '`y`')
  ym[rows, , drop = FALSE]
}


# A data frame or matrix as a double matrix with column names; factor and
# character columns are an error.
cross_numeric <- function (d, what) {

  if (ncol(d) == 0) stop(what, ' has no columns.', call. = FALSE)
  nms <- colnames(d)
  if (is.null(nms)) nms <- as.character(seq_len(ncol(d)))

  if (is.data.frame(d)) {
    for (k in seq_along(d)) {
      v <- d[[k]]
      if (!(is.numeric(v) || is.logical(v)) || !is.null(dim(v)))
        stop(what, " column '", nms[[k]], "' is not numeric; use da_test() for groups.",
             call. = FALSE)
    }
    m <- matrix(as.double(unlist(lapply(d, as.double), use.names = FALSE)), nrow(d), ncol(d))
  }
  else {
    if (!(is.numeric(d) || is.logical(d)))
      stop(what, ' is not numeric; use da_test() for groups.', call. = FALSE)
    m <- matrix(as.double(d), nrow(d), ncol(d))
  }

  m[is.nan(m)] <- NA_real_
  if (any(is.infinite(m)))
    stop(what, ' must not contain infinite values.', call. = FALSE)
  dimnames(m) <- list(NULL, as.character(nms))
  m
}



#------------------------------------------------------------------------------
# `given` as an n x k numeric matrix (model.matrix() columns, without the
# intercept), in x's sample order, with NAs kept; or NULL.
#------------------------------------------------------------------------------
cross_given <- function (given, data, idx, labels, n) {

  if (is.null(given)) return (NULL)

  if (inherits(given, 'formula')) {
    if (is.null(data))
      stop('`given` is a formula, but `data` is not given.', call. = FALSE)
    if (length(given) != 2)
      stop('`given` must be a one-sided formula, such as ~ age + sex.', call. = FALSE)
    df <- data[idx, , drop = FALSE]
  }
  else if (is.character(given) && is.null(dim(given))) {
    if (is.null(data))
      stop('`given` names columns, but `data` is not given.', call. = FALSE)
    miss <- setdiff(given, names(data))
    if (length(miss))
      stop('`given` names ', name_list(miss), ', not column', if (length(miss) > 1) 's',
           ' of `data`.', call. = FALSE)
    df    <- data[idx, given, drop = FALSE]
    given <- NULL
  }
  else {
    if (is.atomic(given) && is.null(dim(given))) {
      nm    <- names(given)
      given <- data.frame(given = given)
      if (!is.null(nm) && !anyNA(nm) && !anyDuplicated(nm)) rownames(given) <- nm
    }
    if (is.matrix(given)) given <- as.data.frame(given)
    if (!is.data.frame(given))
      stop('`given` must be a formula, column names of `data`, or a data frame or ',
           'matrix with one row per sample.', call. = FALSE)
    rn <- cross_row_names(given)
    if (!is.null(rn) && (anyNA(rn) || anyDuplicated(rn))) rn <- NULL
    df    <- given[match_rows(labels, n, rn, nrow(given), ref = '`x`', what = '`given`'), , drop = FALSE]
    given <- NULL
  }

  tt <- if (is.null(given)) stats::terms(~ ., data = df) else stats::terms(given, data = df)
  mf <- stats::model.frame(tt, df, na.action = stats::na.pass)
  # model.matrix() drops rows with NA; keep them, as NA.
  ok <- stats::complete.cases(mf)
  mm <- matrix(NA_real_, n, 0)
  if (any(ok)) {
    mok <- stats::model.matrix(tt, mf[ok, , drop = FALSE])
    mok <- mok[, colnames(mok) != '(Intercept)', drop = FALSE]
    mm  <- matrix(NA_real_, n, ncol(mok), dimnames = list(NULL, colnames(mok)))
    mm[ok, ] <- mok
  }
  if (ncol(mm) == 0)
    stop('`given` has no covariates.', call. = FALSE)
  mm
}



#------------------------------------------------------------------------------
# y's columns, grouped by their pattern of NAs, ready for C_cross_assoc().
#
#   ym      The used samples' y (n_used x q).
#   used    Their 1-based positions among x's samples.
#   G       Their covariates (n_used x k), or NULL.
#
# Returns list(c_pats (for C), n (samples of each column), df (degrees of
# freedom of each column), const_y (constant columns), few (columns with too
# few samples)).
#------------------------------------------------------------------------------
cross_assoc_patterns <- function (ym, used, G, method, min_samples, pm, strata, seed) {

  q    <- ncol(ym)
  miss <- is.na(ym)
  key  <- apply(miss, 2L, function (v) paste(which(v), collapse = ','))
  pat  <- match(key, unique(key))
  n_of <- integer(q)
  df   <- rep(NA_real_, q)
  const_y <- integer(0)
  few  <- integer(0)
  pats <- list()

  if (!is.null(pm$perms) && any(miss))
    stop('A `permutations` matrix needs `y` without NAs: with NAs, each set of ',
         'samples gets its own permutations.', call. = FALSE)

  for (g in seq_len(max(pat))) {

    cols <- which(pat == g)
    rows <- which(!miss[, cols[[1]]])
    m    <- length(rows)
    n_of[cols] <- m

    k  <- 0L
    Q  <- NULL
    if (!is.null(G) && m > 0) {
      Gp <- G[rows, , drop = FALSE]
      if (method == 'spearman') Gp <- matrix(apply(Gp, 2L, rank), nrow = m)
      qz <- qr(cbind(1, Gp), tol = 1e-7)
      k  <- qz$rank - 1L
      Q  <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]
    }
    df[cols] <- m - 2 - k

    if (m < min_samples || m - 2 - k < 1) { few <- c(few, cols); next }

    if (method == 'kendall') {
      codes <- matrix(0L, m, length(cols))
      ystat <- matrix(0, 3, length(cols))
      for (t in seq_along(cols)) {
        v   <- ym[rows, cols[[t]]]
        u   <- sort.int(unique(v))
        codes[, t] <- match(v, u) - 1L
        cnt <- tabulate(codes[, t] + 1L, nbins = length(u))
        cnt <- as.double(cnt[cnt > 1])
        ystat[, t] <- c(sum(cnt * (cnt - 1) / 2), sum(cnt * (cnt - 1) * (cnt - 2)),
                        sum(cnt * (cnt - 1) * (2 * cnt + 5)))
        if (length(u) == 1) const_y <- c(const_y, cols[[t]])
      }
      Y <- NULL
      codes <- t(codes)
    }
    else {
      Y <- matrix(0, m, length(cols))
      for (t in seq_along(cols)) {
        v <- ym[rows, cols[[t]]]
        if (all(v == v[[1]])) { const_y <- c(const_y, cols[[t]]); next }
        if (method == 'spearman') v <- rank(v)
        v  <- v - mean(v)
        s0 <- sqrt(sum(v * v))
        if (!is.null(Q)) v <- v - drop(Q %*% crossprod(Q, v))
        s  <- sqrt(sum(v * v))
        if (!(s > 1e-10 * s0)) { const_y <- c(const_y, cols[[t]]); next }
        Y[, t] <- v / s
      }
      Y     <- t(Y)
      codes <- NULL
      ystat <- NULL
    }

    st <- if (is.null(strata)) NULL else as.integer(factor(strata[rows])) - 1L
    sd <- if (g == 1L) seed else pair_seed(seed, g - 1L)

    pats[[length(pats) + 1L]] <- list(
      as.integer(used[rows] - 1L), as.integer(cols - 1L),
      if (is.null(Y)) NULL else as.double(Y),
      if (is.null(codes)) NULL else as.integer(codes),
      if (is.null(ystat)) NULL else as.double(ystat),
      Q, pm$perms, st, as.integer(sd) )
  }

  list(c_pats = pats, n = n_of, df = df, const_y = const_y, few = few)
}



#------------------------------------------------------------------------------
# Methods.
#------------------------------------------------------------------------------

# Header of print(): e.g. "Spearman Correlation: 120 x 12 features, 37 pairs
# with BH-adjusted p <= 0.05".
cross_assoc_label <- function (x) {
  m <- attr(x, 'method', exact = TRUE)
  if (is.null(m)) return (NULL)   # rows subset with `[` drop attributes
  a  <- attr(x, 'alpha', exact = TRUE)
  pa <- attr(x, 'p_adjust', exact = TRUE)
  sprintf(
    '%s: %d x %d features, %d pair%s with %s-adjusted p <= %s',
    cross_assoc_name(x),
    length(attr(x, 'x_features', exact = TRUE)),
    length(attr(x, 'y_features', exact = TRUE)),
    sum(x$reject, na.rm = TRUE), if (sum(x$reject, na.rm = TRUE) == 1) '' else 's',
    pa, format(a) )
}


# The display name of a result's method.
cross_assoc_name <- function (x) {
  m  <- attr(x, 'method', exact = TRUE)
  nm <- if (is.null(m) || !m %in% names(CROSS_ASSOCS)) 'Correlation' else CROSS_ASSOCS[[m]]
  if (!is.null(attr(x, 'given', exact = TRUE))) nm <- paste('Partial', nm)
  nm
}



ca_or <- function (a, b) if (is.null(a)) b else a


#' @export
as.matrix.ecodive_cross_assoc <- function (x, what = 'association', ...) {

  if (!is.character(what) || length(what) != 1 || !what %in% names(x))
    stop('`what` must be the name of a column of `x`, such as "association" ',
         'or "p_adjusted".', call. = FALSE)

  xf <- attr(x, 'x_features', exact = TRUE)
  yf <- attr(x, 'y_features', exact = TRUE)
  if (is.null(xf)) xf <- unique(x$x)
  if (is.null(yf)) yf <- unique(x$y)

  v <- x[[what]]
  m <- matrix(if (is.logical(v)) NA else NA_real_, length(xf), length(yf),
              dimnames = list(xf, yf))
  m[cbind(match(x$x, xf), match(x$y, yf))] <- v
  m
}



#' Plot Associations Between Two Tables
#'
#' The heatmap of a [cross_assoc()] result: one row per feature of `x`, one
#' column per variable of `y`, coloured by the correlation, with stars on
#' the cells whose adjusted p-values are small.
#'
#' @param x  A [cross_assoc()] result.
#'
#' @param n_x,n_y  The most features of `x` (rows) and variables of `y`
#'        (columns) to show: those with the smallest adjusted p-value (ties
#'        broken by the largest \eqn{|r|}). Default: `30` each.
#'
#' @param clust  `'average'` (the default) orders the rows and the columns by
#'        average-linkage clustering ([stats::hclust()]) on \eqn{1 - r}
#'        between their profiles of correlations in the shown block;
#'        any other method of [stats::hclust()] also works. `'none'` keeps
#'        the input order.
#'
#' @param stars  Thresholds of the adjusted p-values for `*`, `**`, `***`
#'        (any number of them, largest first). Default: `c(0.05, 0.01,
#'        0.001)`. `NULL` marks the cells with `reject` by a dot instead.
#'
#' @param col  Colours for the scale from -1 to 1. Default: the diverging
#'        `hcl.colors(65, "Blue-Red 3")`.
#'
#' @param key  Whether to draw the colour key. Default: `TRUE`.
#'
#' @param ...  Passed on to [graphics::image()], such as `main`, `xlab` and
#'        `ylab`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 8, 20), 30, 8,
#'       dimnames = list(NULL, paste0('Taxon', 1:8)))
#'     env <- data.frame(pH = counts[, 1] / 10 + rnorm(30), temp = rnorm(30))
#'     plot(cross_assoc(counts, env))
#'
plot.ecodive_cross_assoc <- function (
    x, n_x = 30, n_y = 30, clust = 'average',
    stars = c(0.05, 0.01, 0.001), col = NULL, key = TRUE, ... ) {

  args <- list(...)
  if (!isTRUE(key) && !isFALSE(key))
    stop('`key` must be TRUE or FALSE.', call. = FALSE)
  for (a in c('n_x', 'n_y')) {
    v <- get(a)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 1)
      stop('`', a, '` must be a number of at least 1.', call. = FALSE)
  }
  if (!is.character(clust) || length(clust) != 1 || is.na(clust))
    stop("`clust` must be 'average', another method of hclust(), or 'none'.", call. = FALSE)
  if (!is.null(stars) && (!is.numeric(stars) || anyNA(stars) || any(stars <= 0 | stars > 1)))
    stop('`stars` must be NULL or p-value thresholds between 0 and 1.', call. = FALSE)

  r  <- as.matrix(x, 'association')
  pa <- as.matrix(x, 'p_adjusted')
  rj <- as.matrix(x, 'reject')
  P  <- nrow(r); Q <- ncol(r)
  if (P == 0 || Q == 0) stop('`x` has no pairs to plot.', call. = FALSE)


  #--------------------------------------------------
  # The rows and columns shown.
  #--------------------------------------------------
  pick <- function (pmin, rmax, k) {
    o <- order(pmin, -rmax, na.last = TRUE)
    sort.int(o[seq_len(min(k, length(o)))])
  }
  na_inf <- function (v) ifelse(is.na(v), Inf, v)
  rmax_of <- function (m, dim) apply(abs(m), dim, function (v) if (all(is.na(v))) -Inf else max(v, na.rm = TRUE))
  pmin_of <- function (m, dim) apply(m, dim, function (v) if (all(is.na(v))) Inf else min(v, na.rm = TRUE))
  rows <- pick(pmin_of(pa, 1L), rmax_of(r, 1L), n_x)
  cols <- pick(pmin_of(pa[rows, , drop = FALSE], 2L), rmax_of(r[rows, , drop = FALSE], 2L), n_y)
  r  <- r[rows, cols, drop = FALSE]
  pa <- pa[rows, cols, drop = FALSE]
  rj <- rj[rows, cols, drop = FALSE]

  if (clust != 'none') {
    clust_order <- function (m) {
      if (nrow(m) < 3) return (seq_len(nrow(m)))
      s <- suppressWarnings(stats::cor(t(m), use = 'pairwise.complete.obs'))
      d <- 1 - s
      d[!is.finite(d)] <- 1
      diag(d) <- 0
      stats::hclust(stats::as.dist(d), method = clust)$order
    }
    ro <- clust_order(r)
    co <- clust_order(t(r))
    r  <- r[ro, co, drop = FALSE]
    pa <- pa[ro, co, drop = FALSE]
    rj <- rj[ro, co, drop = FALSE]
  }
  p <- nrow(r); q <- ncol(r)


  #--------------------------------------------------
  # Colours and margins.
  #--------------------------------------------------
  pal <- if (!is.null(col)) col else grDevices::hcl.colors(65, 'Blue-Red 3')
  brk <- seq(-1, 1, length.out = length(pal) + 1L)

  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  mai <- graphics::par('mai')
  mai[2L] <- max(mai[2L], max(graphics::strwidth(rownames(r), 'inch', cex = 0.7), 0, na.rm = TRUE) + 0.25)
  mai[1L] <- max(mai[1L], max(graphics::strwidth(colnames(r), 'inch', cex = 0.7), 0, na.rm = TRUE) + 0.6)
  if (key) mai[4L] <- max(mai[4L], 1)
  graphics::par(mai = mai)

  sub <- c(
    if (P > p) sprintf('%d of %d features shown', p, P),
    if (Q > q) sprintf('%d of %d variables shown', q, Q),
    if (!is.null(stars))
      sprintf('stars: %s-adjusted p < %s', ca_or(attr(x, 'p_adjust', exact = TRUE), 'adjusted'),
              paste(vapply(sort(stars, decreasing = TRUE), format, ''), collapse = ', '))
    else
      sprintf('dots: %s-adjusted p <= %s', ca_or(attr(x, 'p_adjust', exact = TRUE), 'adjusted'),
              format(ca_or(attr(x, 'alpha', exact = TRUE), 0.05))) )

  # image() draws its first row at the bottom: rows reversed, to put the
  # first feature at the top.
  z <- t(r[rev(seq_len(p)), , drop = FALSE])
  z[!is.na(z) & z > 1] <- 1
  z[!is.na(z) & z < -1] <- -1

  plot_do(
    graphics::image, c(list(seq_len(q), seq_len(p), z), args),
    list(
      col = pal, breaks = brk, axes = FALSE, xlab = '', ylab = '',
      main = cross_assoc_name(x), sub = paste(sub, collapse = '; ') ))

  graphics::axis(2, at = p - seq_len(p) + 1, labels = rownames(r), las = 1, tick = FALSE, cex.axis = 0.7)
  graphics::axis(1, at = seq_len(q), labels = colnames(r), las = 2, tick = FALSE, cex.axis = 0.7)
  graphics::box()


  #--------------------------------------------------
  # Marks.
  #--------------------------------------------------
  yy <- p - row(r) + 1
  xx <- col(r)
  if (!is.null(stars)) {
    s   <- sort(stars, decreasing = TRUE)
    lev <- matrix(0L, p, q)
    for (k in seq_along(s)) lev <- lev + (!is.na(pa) & pa < s[[k]])
    lab <- strrep('*', lev)
    on  <- lev > 0
    if (any(on))
      graphics::text(xx[on], yy[on], lab[on], cex = 0.8, xpd = NA)
  }
  else {
    on <- !is.na(rj) & rj
    if (any(on)) graphics::points(xx[on], yy[on], pch = 20, cex = 0.8)
  }

  if (key) dist_plot_key(pal, brk, 'r')

  invisible(x)
}
