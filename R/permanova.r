# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Permutational Multivariate Analysis of Variance (PERMANOVA)
#'
#' Tests whether sample metadata, such as treatment groups, explain the
#' differences between samples in a distance matrix. Gives the same table as
#' `vegan::adonis2()`: degrees of freedom, sums of squares, \eqn{R^2}, and
#' pseudo-\eqn{F} statistics, with permutation p-values. Computes them in
#' parallel C, straight from the distances, without the \eqn{n \times n}
#' matrices that `adonis2()` builds for every permutation.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s.
#'
#' @param data  The sample metadata: a data frame with one row per sample.
#'        Rows are matched to the samples in `dm` by row name, and rows for
#'        other samples are ignored; without row names, they must be in the
#'        same order as `dm`. A single vector or factor is also accepted, as
#'        one variable named `group`.
#'
#' @param formula  A one-sided formula naming the variables in `data` to
#'        test, such as `~ treatment` or `~ treatment * site + ph`. Factors,
#'        numeric covariates, and interactions are all allowed. Wrap
#'        variables in `Condition()` to remove their effect before the other
#'        terms are tested, as in `vegan::dbrda()`. A left-hand side, if
#'        given, is ignored. The default, `NULL`, uses every column of
#'        `data` (`~ .`).
#'
#' @param by  How to test the terms of `formula`:
#'
#'   * `'overall'`: One test of all the terms together. The default, as in
#'     `adonis2()` since vegan 2.6-8; earlier versions of `adonis2()`
#'     defaulted to `by = "terms"`, so compare with those using
#'     `by = 'terms'`.
#'   * `'terms'`: Each term in turn, after the ones before it (sequential, or
#'     type I, sums of squares). Their order in `formula` matters.
#'   * `'margin'`: Each term after all the others (marginal sums of
#'     squares), for the terms that can be dropped without breaking
#'     marginality: in `~ a * b`, only `a:b`.
#'   * `'onedf'`: Each column of the model matrix in turn, such as each
#'     contrast of a factor.
#'
#' @param permutations  The number of permutations. Default: `999`.
#'        Alternatively, a matrix with one permutation of `1:n` per row, as
#'        from `permute::shuffleSet()`, to use instead of random ones.
#'
#' @param strata  Permute samples only within these groups (blocks), for
#'        designs such as repeated measures. The name of a column of `data`,
#'        or a vector with one value per row of `data`. Default: `NULL`
#'        (free permutation).
#'
#' @param seed  An integer seed for the random permutations. The same seed
#'        gives the same p-values. Default: `0`
#'
#' @return A data frame (of class `"anova"`, so it prints like
#'   `stats::anova()` tables) with one row per term tested (or one `"Model"`
#'   row when `by = 'overall'`), then `"Residual"` and `"Total"` rows. When
#'   `formula` has `Condition()` terms, a `"Conditional"` row comes first. Its
#'   columns are:
#'
#'   * `Df`: Degrees of freedom.
#'   * `SumOfSqs`: Sum of squares.
#'   * `R2`: Proportion of the total sum of squares, `SumOfSqs / Total`.
#'   * `F`: Pseudo-\eqn{F} statistic.
#'   * `Pr(>F)`: Permutation p-value.
#'
#'   Attribute `f_perm` is a matrix of the \eqn{F} statistics of each
#'   permutation (rows) for each test (columns), as `adonis2()`'s `F.perm`.
#'
#' @details
#'
#' **Method**
#'
#' PERMANOVA (Anderson 2001) partitions the variation in the distances among
#' the terms of a linear model, following McArdle & Anderson (2001): with
#' \eqn{A} the squared distances, \eqn{G = -\frac{1}{2}JAJ} Gower's centred
#' matrix, and \eqn{H} the hat matrix of the (centred) model, a model's sum
#' of squares is \eqn{trace(HGH)}, out of a total \eqn{trace(G) =
#' \frac{1}{n}\sum_{i<j}d_{ij}^2}. The pseudo-\eqn{F} of a term with
#' \eqn{q} degrees of freedom is
#' \deqn{F = \frac{SS_{term} / q}{SS_{residual} / (n - 1 - p)}}
#' where \eqn{p} counts the degrees of freedom of the whole model. For each
#' permutation, the samples' distances are shuffled relative to the metadata
#' and \eqn{F} is recomputed. The p-value is the proportion of permutations,
#' counting the observed one, with an \eqn{F} at least as large:
#' \eqn{(1 + \#\{F_{perm} \ge F\}) / (1 + permutations)}. As in `adonis2()`,
#' permuted statistics within `sqrt(.Machine$double.eps)` of the observed one
#' count as ties.
#'
#' Partial models, with `Condition()` terms, follow `vegan::dbrda()` and
#' `vegan::anova.cca()`'s default "reduced" model: the conditions' effect is
#' removed from the distances, and the residuals are permuted.
#'
#' Every term is tested with the same permutations of the (conditioned)
#' distances, as in `adonis2()`. PERMANOVA+ for PRIMER instead permutes the
#' residuals of a separate reduced model for each term (Freedman & Lane
#' 1983), which can give slightly different p-values for terms after the
#' first.
#'
#'
#' **Speed and memory**
#'
#' With \eqn{Q} an orthonormal basis of the model matrix, \eqn{trace(HGH) =
#' \sum_c q_c^T G q_c}, and each \eqn{q^T G q = -\sum_{i<j} q_i q_j
#' d_{ij}^2} is read straight from `dm`. A permutation of the samples
#' permutes the rows of \eqn{Q}, so each permutation costs one pass over the
#' distances, shared by many vectors at once, with no \eqn{n \times n}
#' matrix. For a single factor, the sum of squares within groups is summed
#' over the pairs of samples in the same group instead (Anderson 2001), which
#' is faster still. Permutations run in parallel.
#'
#'
#' **Reproducibility**
#'
#' The permutations come from `seed`, not from R's random number generator,
#' so `set.seed()` has no effect, and each permutation has a random number
#' stream of its own. P-values are identical from run to run and for any
#' number of `cpus`. The random permutations differ from `adonis2()`'s, so
#' p-values agree with it only up to sampling error, unless the same
#' `permutations` matrix is given to both.
#'
#'
#' **Differences from adonis2()**
#'
#' * `permanova()` takes a `dist` object, not a formula with one on the
#'   left-hand side; use [beta_div()] or a metric function to compute it.
#' * `strata` permutes freely within strata. For other restricted designs,
#'   pass a `permutations` matrix from `permute::shuffleSet()`.
#' * When there are fewer possible permutations than requested, `adonis2()`
#'   (via `permute`) runs each of them once instead. `permanova()` always
#'   draws random permutations.
#'
#' @references
#' Anderson, M. J. (2001). A new method for non-parametric multivariate
#' analysis of variance. *Austral Ecology*, 26(1), 32-46.
#' \doi{10.1111/j.1442-9993.2001.01070.pp.x}
#'
#' Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of
#' reported significance levels. *Journal of Business & Economic
#' Statistics*, 1(4), 292-298. \doi{10.1080/07350015.1983.10509354}
#'
#' McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
#' community data: a comment on distance-based redundancy analysis.
#' *Ecology*, 82(1), 290-297.
#' \doi{10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2}
#'
#' @seealso [beta_div()], [pcoa()]
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
#'       ph        = round(rnorm(30, 7, 0.3), 1) )
#'     high <- meta$treatment == 'high'
#'     counts[high, 1:5] <- counts[high, 1:5] + 5
#'
#'     dm <- bray(counts)
#'
#'     # Do treatment and pH explain differences in community composition?
#'     permanova(dm, meta, ~ treatment + ph)
#'
#'     # Each term in turn, after the ones before it.
#'     permanova(dm, meta, ~ treatment + ph, by = 'terms')
#'
#'     # pH, after removing the effect of treatment.
#'     permanova(dm, meta, ~ ph + Condition(treatment))
#'
permanova <- function (
    dm, data, formula = NULL, by = 'overall', permutations = 999L,
    strata = NULL, seed = 0L, cpus = n_cpus() ) {

  validate_seed()
  validate_cpus()

  by <- match.arg(by, c('overall', 'terms', 'margin', 'onedf'))

  dv  <- dist_vector(dm)
  n   <- dv$n
  des <- model_design(data, formula, dv$labels, n, strata)
  pm  <- perm_matrix(permutations, n)

  if (n < 2)
    stop('`dm` must have at least 2 samples.', call. = FALSE)

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)


  #--------------------------------------------------------
  # Bases for the full model and, for by = 'margin', each
  # model without one of its droppable terms.
  #--------------------------------------------------------
  qz   <- condition_basis(des$z)
  rz   <- if (is.null(qz)) 0L else ncol(qz)
  full <- design_basis(des$x, qz)

  ass <- des$assign[full$keep]  # each basis column's term

  sets  <- list(full)
  scope <- character(0)

  if (by == 'margin' && full$rank > 0) {

    scope <- intersect(stats::drop.scope(des$terms), des$term_labels[unique(ass)])

    for (lab in scope) {
      i <- match(lab, des$term_labels)
      sets[[length(sets) + 1]] <- design_basis(des$x[, des$assign != i, drop = FALSE], qz)
    }
  }


  #--------------------------------------------------------
  # Sums of squares for the identity and each permutation.
  # One factor: the fast within-group path.
  #--------------------------------------------------------
  grp <- NULL
  if (is.null(qz) && by != 'onedf' && length(des$term_labels) == 1 && full$rank > 0) {
    v <- des$data[[des$term_labels]]
    if (is.factor(v) || is.character(v) || is.logical(v))
      grp <- as.integer(factor(v)) - 1L
  }

  if (full$rank == 0) {
    res <- perm_ss(dv$dm, n, cols = matrix(0, n, 0), pm = perm_matrix(0L, n), cpus = cpus)
    res$ss <- matrix(0, pm$n_perm + 1L, 0)
  }
  else if (!is.null(grp)) {
    res <- perm_ss(dv$dm, n, grp = grp, pm = pm, strata = des$strata, seed = seed, cpus = cpus)
    sets[[1]]$set_ss <- res$ss[, 1]
    for (s in seq_along(sets)[-1]) sets[[s]]$set_ss <- numeric(pm$n_perm + 1L)
  }
  else {
    widths <- vapply(sets, `[[`, integer(1), 'rank')
    cols   <- do.call(cbind, lapply(sets, `[[`, 'cols'))
    res    <- perm_ss(
      dv$dm, n, cols = cols, sets = rep(seq_along(sets), widths), qz = qz,
      pm = pm, strata = des$strata, seed = seed, cpus = cpus )
    ends <- cumsum(widths)
    for (s in seq_along(sets)) {
      j <- seq_len(widths[[s]]) + ends[[s]] - widths[[s]]
      sets[[s]]$col_ss <- res$ss[, j, drop = FALSE]
      sets[[s]]$set_ss <- rowSums(res$ss[, j, drop = FALSE])
    }
  }

  total <- res$total

  # The conditions' sum of squares (observed only).
  ss_z <- 0
  if (!is.null(qz))
    ss_z <- sum(perm_ss(dv$dm, n, cols = qz, pm = perm_matrix(0L, n), cpus = cpus)$ss)

  ss_full <- if (full$rank > 0) sets[[1]]$set_ss else numeric(pm$n_perm + 1L)
  df_res  <- n - 1L - rz - full$rank
  ss_res  <- total - ss_z - ss_full  # per permutation


  #--------------------------------------------------------
  # Numerators: one column per test, one row per permutation.
  #--------------------------------------------------------
  if (full$rank == 0) {
    num <- matrix(0, pm$n_perm + 1L, 0)
    dfs <- integer(0)
  }
  else if (by == 'overall') {
    num <- matrix(ss_full, ncol = 1, dimnames = list(NULL, 'Model'))
    dfs <- full$rank
  }
  else if (by == 'terms') {
    if (is.null(grp)) {
      terms <- unique(ass)
      num   <- vapply(terms, function (t) rowSums(sets[[1]]$col_ss[, ass == t, drop = FALSE]),
                      numeric(pm$n_perm + 1L))
      num   <- matrix(num, ncol = length(terms))
      dfs   <- vapply(terms, function (t) sum(ass == t), integer(1))
    }
    else {
      terms <- 1L
      num   <- matrix(ss_full, ncol = 1)
      dfs   <- full$rank
    }
    colnames(num) <- des$term_labels[terms]
  }
  else if (by == 'onedf') {
    num <- sets[[1]]$col_ss
    dfs <- rep(1L, full$rank)
    colnames(num) <- colnames(full$cols)
  }
  else { # margin
    num <- vapply(sets[-1], function (s) ss_full - s$set_ss, numeric(pm$n_perm + 1L))
    num <- matrix(num, ncol = length(scope), dimnames = list(NULL, scope))
    dfs <- vapply(sets[-1], function (s) full$rank - s$rank, integer(1))
  }


  #--------------------------------------------------------
  # F statistics and p-values.
  #--------------------------------------------------------
  f_all <- sweep(num, 2, dfs, '/') / (ss_res / df_res)
  f_all <- matrix(f_all, nrow = pm$n_perm + 1L, dimnames = list(NULL, colnames(num)))
  f_obs <- f_all[1, ]
  f_prm <- f_all[-1, , drop = FALSE]

  if (df_res <= 0) {
    warning('The model has no residual degrees of freedom, so it cannot be tested.', call. = FALSE)
    f_obs[] <- NA_real_
  }

  eps <- sqrt(.Machine$double.eps)
  p   <- if (pm$n_perm > 0 && df_res > 0)
    (colSums(sweep(f_prm, 2, f_obs - eps, '>=')) + 1) / (pm$n_perm + 1)
  else rep(NA_real_, ncol(num))


  #--------------------------------------------------------
  # The table.
  #--------------------------------------------------------
  rows <- c(
    if (rz > 0) 'Conditional',
    colnames(num), 'Residual', 'Total' )

  tbl_df <- c(if (rz > 0) rz, dfs,      df_res,    n - 1L)
  tbl_ss <- c(if (rz > 0) ss_z, num[1, ], ss_res[[1]], total)
  tbl_f  <- c(if (rz > 0) NA, f_obs,    NA,        NA)
  tbl_p  <- c(if (rz > 0) NA, p,        NA,        NA)

  tbl <- data.frame(
    Df       = as.double(tbl_df),
    SumOfSqs = tbl_ss,
    R2       = tbl_ss / total,
    F        = unname(tbl_f),
    P        = unname(tbl_p),
    row.names   = rows,
    check.names = FALSE )
  names(tbl)[5] <- 'Pr(>F)'

  rownames(f_prm) <- NULL

  heading <- paste0(
    'Permutation test for PERMANOVA\n',
    switch(by,
      overall = 'Model tested as a whole\n',
      terms   = 'Terms added sequentially (first to last)\n',
      margin  = 'Marginal effects of terms\n',
      onedf   = 'Sequential test for contrasts\n' ),
    if (!is.null(pm$perms)) 'Permutation: supplied matrix\n'
    else if (!is.null(des$strata)) 'Permutation: free, within strata\n'
    else 'Permutation: free\n',
    'Number of permutations: ', pm$n_perm, '\n' )

  structure(
    .Data   = tbl,
    heading = heading,
    f_perm  = f_prm,
    class   = c('anova', 'data.frame') )
}
