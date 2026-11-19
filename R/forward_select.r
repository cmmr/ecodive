# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Forward Selection of Explanatory Variables
#'
#' Picks, one at a time, the explanatory variables that best explain a
#' community, by redundancy analysis of a count table ([rda()]) or
#' distance-based redundancy analysis of distances ([dbrda()]). At each step
#' the variable that raises the adjusted \eqn{R^2} most is added, if a
#' permutation test of it (given the variables already chosen) is
#' significant, and if the adjusted \eqn{R^2} stays within that of all the
#' candidates together: the double stopping rule of Blanchet et al. (2008),
#' as in `vegan::ordiR2step()`, with the same results. Usually run before
#' [varpart()].
#'
#' @inheritParams documentation
#' @inheritParams varpart
#' @inheritParams permanova
#' @inherit pseudocount_section
#' @inheritSection documentation Input Types
#'
#' @param x  The community: a `dist` object of distances among the samples,
#'        such as from [beta_div()] or [bray()], for db-RDA; or a count
#'        table (samples \eqn{\times} features, or any of the table types
#'        listed under "Input Types"), for RDA. A plain matrix is always a
#'        table; give distances as a `dist` object (`as.dist(m)`).
#'
#' @param formula  A one-sided formula naming the candidate variables in
#'        `data`, such as `~ ph + temp + depth` (vegan's `scope`). Each
#'        term is one candidate: a factor enters whole, with all its
#'        contrasts, and an interaction becomes a candidate only once the
#'        terms it contains are chosen. Variables wrapped in `Condition()`
#'        are not candidates: every model, test, and adjusted \eqn{R^2} is
#'        conditioned on them (partial RDA), as when `ordiR2step()` starts
#'        from a partial model. Default: `NULL`, every column of `data`
#'        (except `strata`) as a candidate.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the variables in `formula` (and `strata`). Rows are
#'        matched to the samples of `x` by row name, and rows for other
#'        samples are ignored; without row names, they must be in the
#'        samples' order.
#'
#' @param alpha  The significance level: a variable is added only if its
#'        permutation p-value is at most `alpha` (vegan's `Pin`). Default:
#'        `0.05`.
#'
#' @param max_vars  The most terms to select, or `NULL` for no limit (as
#'        adespatial's `forward.sel(K)`). Default: `NULL`.
#'
#' @param global_test  Test the global model (all the candidates together)
#'        first, and select nothing if it is not significant (p-value above
#'        `alpha`), as Blanchet et al. (2008) recommend. `FALSE` skips the
#'        test, as `ordiR2step()` does. Default: `TRUE`.
#'
#' @param permutations  The number of permutations for each test. Default:
#'        `999`. Alternatively, a matrix with one permutation of `1:n` per
#'        row, as from `permute::shuffleSet()`, to use instead of random
#'        ones; the same permutations then serve every test, as when the
#'        same matrix is given to `ordiR2step()`. `0` runs no tests:
#'        variables are then chosen by adjusted \eqn{R^2} alone (rules 1, 2
#'        and 4 below), `alpha` is not used, and p-values are `NA`.
#'
#' @param strata  Permute samples only within these groups (blocks): the
#'        name of a column of `data`, or a vector with one value per row of
#'        `data` (or, without `data` row names, per sample). For other
#'        restricted designs, give a `permutations` matrix instead. Default:
#'        `NULL`.
#'
#' @return A data frame of class `c('ecodive_forward_select',
#'   'data.frame')`, with one row per variable selected, in the order
#'   selected (no rows if none was), and columns:
#'
#'   * `variable`: The term added.
#'   * `df`: The number of model columns it adds.
#'   * `r_squared`: The \eqn{R^2} it adds.
#'   * `r_squared_cum`: The \eqn{R^2} of the model with it and every
#'     variable before it.
#'   * `adj_r_squared`: The adjusted \eqn{R^2} it adds.
#'   * `adj_r_squared_cum`: The adjusted \eqn{R^2} of the model with it and
#'     every variable before it (`ordiR2step()`'s `R2.adj`).
#'   * `F`: The pseudo-\eqn{F} of the variable, given the variables before
#'     it (and any `Condition()` terms).
#'   * `p_value`: Its permutation p-value.
#'
#'   Attributes:
#'
#'   * `selected`: A one-sided formula of the variables selected, in order,
#'     ready for [varpart()] or [rda()] (`~ 1` if none was). `Condition()`
#'     terms are not included.
#'   * `global`: A one-row data frame of the global model, all candidates
#'     together (given any `Condition()` terms): `df`, `r_squared`,
#'     `adj_r_squared`, `F`, and `p_value` (`NA` without `global_test`).
#'   * `stop`: Why selection stopped: `'global test not significant'`,
#'     `'no candidates left'`, `'adjusted R-squared does not increase'`,
#'     `'adjusted R-squared would exceed the global model\'s'`, `'not
#'     significant'`, or `'max_vars reached'`.
#'   * `next_best`: The candidate that would have been added next, as a
#'     one-row data frame with the same columns as the result (`p_value`
#'     is `NA` if it was not tested), or `NULL`.
#'   * `method`: `'rda'` or `'dbrda'`.
#'   * `norm`: The normalization, for RDA.
#'   * `n`: The number of samples.
#'   * `n_perm`: The number of permutations per test.
#'   * `permutation`: How samples were permuted: `'free'`, `'within
#'     strata'`, or `'supplied matrix'`.
#'
#' @details
#'
#' **Method**
#'
#' With the variables chosen so far, \eqn{S} (and any `Condition()` terms,
#' \eqn{Z}), each remaining candidate \eqn{t} gives a model \eqn{S + t},
#' whose adjusted \eqn{R^2} is Ezekiel's, \eqn{1 - (1 - R^2)(n - 1)/(n - m
#' - 1)} for \eqn{n} samples and \eqn{m} model columns (with conditions, the
#' adjusted \eqn{R^2} of \eqn{Z + S + t} less that of \eqn{Z}, as
#' `vegan::RsquareAdj()` gives it). The best candidate is added if:
#'
#' 1. its adjusted \eqn{R^2} is higher than that of \eqn{S};
#' 2. its adjusted \eqn{R^2} is no higher than that of the global model,
#'    with every candidate (skipped if the global model has no residual
#'    degrees of freedom);
#' 3. its permutation test has a p-value of at most `alpha`; and
#' 4. fewer than `max_vars` variables have been chosen.
#'
#' Otherwise selection stops. The test of rule 3 is that of vegan's
#' `anova()` of the two nested models, \eqn{Z + S} and \eqn{Z + S + t}:
#' the pseudo-\eqn{F} of the sum of squares that \eqn{t} adds, over the
#' residual of the larger model, where for each permutation each model's
#' sum of squares is recomputed with the samples permuted (after removing
#' \eqn{Z}'s effect, as in [permanova()]). Without `Condition()` terms this
#' is a plain permutation of the samples. (It differs slightly from a test
#' of `~ t + Condition(S + Z)`, which permutes the residuals of \eqn{S}
#' too.) The global model is tested as [permanova()] tests it. Every test
#' uses the same permutations.
#'
#' Rules 1 to 3 are `ordiR2step()`'s; rule 2
#' guards against the inflated type I error and overestimated explained
#' variation of forward selection (Blanchet et al. 2008).
#'
#' **Speed**
#'
#' One pass over the table's non-zeros (or the distances) gives the cross
#' products of every candidate's columns, from which every candidate's
#' adjusted \eqn{R^2}, at every step, follows without another pass (as in
#' [varpart()]). Only the chosen candidate of each step is tested, with the
#' permutation kernel of [permanova()], in parallel, reading the distances
#' as a packed triangle; for RDA, the Euclidean distances among the
#' (normalized) samples, which give the same test. Results are identical
#' for any number of `cpus`.
#'
#'
#' **Differences from vegan::ordiR2step()**
#'
#' * The result is a data frame of the steps, not the final model; the
#'   steps' table is close to `ordiR2step()`'s `anova` component.
#' * The global model is tested first, by default (`global_test`).
#' * `max_vars` limits the number of variables.
#' * There is no CCA (`ordiR2step()` on `cca()` models); use RDA of a
#'   Hellinger-normalized table instead.
#' * vegan's `ordiR2step()` chooses the scope from a formula with the
#'   table on its left; here the community is `x`, and `formula` lists the
#'   candidates.
#'
#' @references
#' Blanchet, F. G., Legendre, P., & Borcard, D. (2008). Forward selection
#' of explanatory variables. *Ecology*, 89, 2623-2632.
#' \doi{10.1890/07-0986.1}
#'
#' @seealso [varpart()], [rda()], [dbrda()], [permanova()], [bioenv()]
#'
#' @export
#' @examples
#'     # 40 samples: four environmental variables, two of which matter.
#'     set.seed(1)
#'     meta <- data.frame(
#'       ph    = runif(40, 5, 8),
#'       temp  = runif(40, 10, 20),
#'       depth = runif(40),
#'       site  = factor(rep(c('a', 'b', 'c', 'd'), 10)) )
#'     counts <- matrix(rpois(40 * 15, 4), 40, 15)
#'     counts[, 1:5] <- counts[, 1:5] + round(10 * (meta$ph - 5))
#'     counts[, 6:9] <- counts[, 6:9] + round(2 * (meta$temp - 10))
#'
#'     # tb-RDA on the Hellinger-normalized table.
#'     fs <- forward_select(counts, ~ ph + temp + depth + site, meta, norm = 'hellinger')
#'     fs
#'
#'     # db-RDA on Bray-Curtis distances, then partition the variation.
#'     fs <- forward_select(bray(counts), ~ ph + temp + depth + site, meta)
#'     attr(fs, 'selected')
#'     varpart(bray(counts), env = attr(fs, 'selected'), site = ~ site, data = meta)
#'
#'     # Conditioned on site.
#'     forward_select(bray(counts), ~ ph + temp + depth + Condition(site), meta)
#'
forward_select <- function (
    x, formula = NULL, data, norm = 'none', pseudocount = NULL, margin = 1L,
    alpha = 0.05, max_vars = NULL, global_test = TRUE, permutations = 999L,
    strata = NULL, seed = 0L, cpus = n_cpus() ) {

  env <- parent.frame()

  validate_seed()
  validate_cpus()

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha > 1)
    stop('`alpha` must be a number greater than 0 and at most 1.', call. = FALSE)

  if (!isTRUE(global_test) && !isFALSE(global_test))
    stop('`global_test` must be TRUE or FALSE.', call. = FALSE)

  if (!is.null(max_vars) && (
      !is.numeric(max_vars) || length(max_vars) != 1 || is.na(max_vars) ||
      max_vars < 1 || max_vars %% 1 != 0 ))
    stop('`max_vars` must be a whole number, 1 or more, or NULL.', call. = FALSE)


  #--------------------------------------------------------
  # The response and the candidates.
  #--------------------------------------------------------
  y <- varpart_response(x, norm, pseudocount, margin, FALSE)
  n <- y$n

  if (n < 3)
    stop('`x` must have at least 3 samples.', call. = FALSE)

  if (missing(data) || is.null(data))
    stop('`data` is required: a data.frame of the candidate variables.', call. = FALSE)
  if (inherits(data, 'formula'))
    stop('Give `formula` before `data`: forward_select(x, ~ a + b, data).', call. = FALSE)
  if (!is.data.frame(data) && !is.matrix(data))
    stop('`data` must be a data.frame.', call. = FALSE)

  des <- model_design(data, formula, y$labels, n, strata, env = env)
  pm  <- perm_matrix(permutations, n)

  if (!is.null(pm$perms) && !is.null(des$strata))
    stop('Give `strata` or a `permutations` matrix, not both.', call. = FALSE)

  terms_ <- des$term_labels
  n_t    <- length(terms_)
  xz     <- cbind(des$x, des$z)
  storage.mode(xz) <- 'double'
  zcols  <- ncol(des$x) + seq_len(if (is.null(des$z)) 0L else ncol(des$z))
  tcols  <- lapply(seq_len(n_t), function (t) which(des$assign == t))

  model <- response_cross(y, xz, cpus)
  total <- model$total
  if (!(total > 0))
    stop('`x` has no variation: every sample is the same.', call. = FALSE)


  #--------------------------------------------------------
  # R-squared, adjusted R-squared, and rank of a model with
  # terms `ts`, conditioned on Z (as vegan::RsquareAdj()).
  #--------------------------------------------------------
  fit_z <- subset_inertia(model, zcols)
  rz    <- fit_z$rank
  r2z   <- fit_z$inertia / total
  adjz  <- if (rz > 0) adj_r_squared(r2z, n, rz) else 0

  fit_terms <- function (ts) {
    cols <- c(zcols, unlist(tcols[ts]))
    f    <- subset_inertia(model, cols)
    rk   <- f$rank - rz
    r2   <- f$inertia / total - r2z
    adj  <- adj_r_squared(f$inertia / total, n, f$rank) - adjz
    list(r2 = r2, adj = adj, rank = rk, ss = f$inertia)
  }

  # The pseudo-F of model `big` against the smaller `small` (fit_terms()).
  obs_f <- function (big, small) {
    q    <- big$rank - small$rank
    df_r <- n - 1L - rz - big$rank
    if (q <= 0 || df_r <= 0) return (NA_real_)
    ((big$ss - small$ss) / q) / ((total - big$ss) / df_r)
  }
  base <- list(r2 = 0, adj = 0, rank = 0L, ss = fit_z$inertia)


  #--------------------------------------------------------
  # Permutation tests: the p-value of the F of `ts` given
  # `given` (and Z), with the same permutations every time.
  # As vegan's anova() of two nested models: each model's
  # sum of squares is permuted as in permanova() (under the
  # reduced model, with Z), and F is their difference over
  # the larger model's residual.
  #--------------------------------------------------------
  testing <- pm$n_perm > 0
  dm      <- NULL
  if (testing) {
    dm <- if (y$kind == 'dbrda') y$dm
          else dist_vector(euclidean(
            y$counts, margin = y$margin, norm = y$norm_name,
            pseudocount = y$pseudocount, cpus = cpus ))$dm
  }

  # The conditions' basis and sum of squares, fixed for every test.
  qz   <- if (testing) condition_basis(if (length(zcols)) xz[, zcols, drop = FALSE])
  ss_z <- if (!is.null(qz)) sum(perm_ss(dm, n, cols = qz, pm = perm_matrix(0L, n), cpus = cpus)$ss) else 0

  step_test <- function (ts, given) {

    if (!testing) return (NA_real_)

    # The columns of the smaller model (`given`) first, then the added
    # terms'. Gram-Schmidt within the one set makes the first r_s
    # (permuted, Z-orthogonal) columns span the smaller model.
    gx   <- unlist(tcols[given])
    xs   <- xz[, c(gx, unlist(tcols[ts])), drop = FALSE]
    b    <- design_basis(xs, qz)
    r_s  <- sum(b$keep <= length(gx))
    q    <- b$rank - r_s
    df_r <- n - 1L - rz - b$rank

    if (q == 0 || df_r <= 0) return (NA_real_)

    res  <- perm_ss(dm, n, cols = b$cols, qz = qz, pm = pm,
                    strata = des$strata, seed = seed, cpus = cpus)
    ss_f <- rowSums(res$ss)
    ss_s <- rowSums(res$ss[, seq_len(r_s), drop = FALSE])
    f    <- ((ss_f - ss_s) / q) / ((res$total - ss_z - ss_f) / df_r)
    perm_p(f[[1]], f[-1])
  }


  #--------------------------------------------------------
  # The global model.
  #--------------------------------------------------------
  all_fit <- fit_terms(seq_len(n_t))
  if (all_fit$rank == 0)
    stop('No candidate varies among the samples (or all are aliased ',
         'with the `Condition()` terms).', call. = FALSE)

  global <- list(F = obs_f(all_fit, base), p = NA_real_)
  if (global_test && testing) {
    global$p <- step_test(seq_len(n_t), integer(0))
    if (is.na(global$p))
      warning('The global model has no residual degrees of freedom, so it ',
              'cannot be tested; selecting without the global test.', call. = FALSE)
  }
  if (is.na(all_fit$adj))
    warning('The global model has no residual degrees of freedom, so its ',
            'adjusted R-squared is undefined and is not used as a stopping ',
            'rule.', call. = FALSE)


  #--------------------------------------------------------
  # Forward selection.
  #--------------------------------------------------------
  row_of <- function (t, fit, prev, F, p) data.frame(
    variable          = terms_[[t]],
    df                = as.integer(fit$rank - prev$rank),
    r_squared         = fit$r2 - prev$r2,
    r_squared_cum     = fit$r2,
    adj_r_squared     = fit$adj - prev$adj,
    adj_r_squared_cum = fit$adj,
    F                 = F,
    p_value           = p,
    stringsAsFactors  = FALSE )

  chosen  <- integer(0)
  prev    <- base
  rows    <- list()
  nxt     <- NULL
  reason  <- NULL
  max_v   <- if (is.null(max_vars)) Inf else max_vars
  tt      <- des$terms

  if (global_test && testing && !is.na(global$p) && global$p > alpha) {
    reason <- 'global test not significant'
    warning(sprintf(
      'The global model is not significant (p = %s > alpha = %s), so no variable was selected.',
      format(global$p, digits = 4), format(alpha) ), call. = FALSE)
  }

  while (is.null(reason)) {

    # Candidates that respect marginality (stats::add.scope()).
    cur  <- if (length(chosen)) stats::reformulate(terms_[chosen]) else ~ 1
    cand <- match(stats::add.scope(cur, tt), terms_)
    cand <- setdiff(cand[!is.na(cand)], chosen)

    if (!length(cand)) { reason <- 'no candidates left'; break }

    fits <- lapply(cand, function (t) fit_terms(c(chosen, t)))
    adj  <- vapply(fits, `[[`, 0, 'adj')
    if (all(is.na(adj))) { reason <- 'adjusted R-squared does not increase'; break }
    best <- which.max(adj)
    t    <- cand[[best]]
    fit  <- fits[[best]]

    if (length(chosen) >= max_v) {
      nxt    <- row_of(t, fit, prev, obs_f(fit, prev), NA_real_)
      reason <- 'max_vars reached'
      break
    }
    if (!(fit$adj > prev$adj)) {
      nxt    <- row_of(t, fit, prev, obs_f(fit, prev), NA_real_)
      reason <- 'adjusted R-squared does not increase'
      break
    }
    if (!is.na(all_fit$adj) && fit$adj > all_fit$adj) {
      nxt    <- row_of(t, fit, prev, obs_f(fit, prev), NA_real_)
      reason <- 'adjusted R-squared would exceed the global model\'s'
      break
    }

    f_t <- obs_f(fit, prev)
    p_t <- step_test(t, chosen)
    if (testing && !(!is.na(p_t) && p_t <= alpha)) {
      nxt    <- row_of(t, fit, prev, f_t, p_t)
      reason <- 'not significant'
      break
    }

    rows[[length(rows) + 1L]] <- row_of(t, fit, prev, f_t, p_t)
    chosen <- c(chosen, t)
    prev   <- fit
  }

  out <- if (length(rows)) do.call(rbind, rows)
         else row_of(1L, prev, prev, 0, 0)[0, , drop = FALSE]
  rownames(out) <- NULL

  permutation <- if (!is.null(pm$perms)) 'supplied matrix'
                 else if (!is.null(des$strata)) 'within strata'
                 else 'free'

  selected <- if (length(chosen)) stats::reformulate(terms_[chosen]) else ~ 1
  environment(selected) <- env
  attr(out, 'selected') <- selected
  attr(out, 'global') <- data.frame(
    df = as.integer(all_fit$rank), r_squared = all_fit$r2,
    adj_r_squared = all_fit$adj, F = global$F, p_value = global$p )
  attr(out, 'stop')        <- reason
  attr(out, 'next_best')   <- nxt
  attr(out, 'method')      <- y$kind
  attr(out, 'norm')        <- if (y$kind == 'rda') y$norm_name
  attr(out, 'n')           <- n
  attr(out, 'n_perm')      <- pm$n_perm
  attr(out, 'permutation') <- permutation
  class(out) <- c('ecodive_forward_select', 'data.frame')
  out
}



#' @export
print.ecodive_forward_select <- function (x, digits = 4L, ...) {

  method <- attr(x, 'method')
  how    <- switch(
    method,
    rda   = sprintf('RDA%s', if (!identical(attr(x, 'norm'), 'none')) sprintf(' (%s)', attr(x, 'norm')) else ''),
    dbrda = 'db-RDA' )
  g      <- attr(x, 'global')
  fmt    <- function (v) if (is.na(v)) 'NA' else formatC(v, format = 'f', digits = digits)
  n_perm <- attr(x, 'n_perm')

  cat(sprintf('Forward selection by %s: %d samples\n', how, attr(x, 'n')))
  cat(sprintf(
    'Permutations: %s\n',
    if (n_perm == 0) 'none (selected by adjusted R-squared alone)'
    else sprintf('%d (%s)', n_perm, if (attr(x, 'permutation') == 'within strata') 'free, within strata' else attr(x, 'permutation')) ))
  cat(sprintf(
    'Global model: %d df, adj. R-squared %s%s\n\n', g$df, fmt(g$adj_r_squared),
    if (is.na(g$p_value)) '' else sprintf(', F = %s, p = %s', format(signif(g$F, digits)), format(signif(g$p_value, digits))) ))

  if (nrow(x) == 0) {
    cat('No variable selected.\n')
  }
  else {
    df <- as.data.frame(unclass(x)[names(x)], stringsAsFactors = FALSE)
    for (j in c('r_squared', 'r_squared_cum', 'adj_r_squared', 'adj_r_squared_cum'))
      df[[j]] <- round(df[[j]], digits)
    df$F <- signif(df$F, digits)
    print(df)
  }

  cat(sprintf('\nStopped: %s.\n', attr(x, 'stop')))
  nb <- attr(x, 'next_best')
  if (!is.null(nb))
    cat(sprintf(
      'Next best: %s (adj. R-squared %s%s)\n', nb$variable, fmt(nb$adj_r_squared_cum),
      if (is.na(nb$p_value)) '' else sprintf(', F = %s, p = %s', format(signif(nb$F, digits)), format(signif(nb$p_value, digits))) ))
  cat(sprintf('Selected: %s\n', paste(deparse(attr(x, 'selected'), width.cutoff = 500L), collapse = ' ')))

  invisible(x)
}
