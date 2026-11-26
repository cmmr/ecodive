# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Principal Response Curves (PRC)
#'
#' How the communities of treated samples diverge from those of a control
#' over time, summarized as one curve per treatment (van den Brink & ter
#' Braak 1999). The method is a partial redundancy analysis,
#' `rda(counts ~ treatment:time + Condition(time))`: the effect of time is
#' removed, and what remains is the treatments' effect at each time point.
#' Each treatment's curve is its canonical coefficient at each time point,
#' measured from the control (whose curve is 0), and each feature's weight
#' says how closely its abundance follows the curves. Gives the same results
#' as `vegan::prc()`, with permutation tests of the first axis and of all
#' axes as `anova(vegan::prc(...), first = TRUE)` and `anova(...)`.
#'
#' @inherit pseudocount_section
#' @inheritParams rda
#'
#' @param counts  A numeric matrix of feature abundances (samples
#'        \eqn{\times} features), or any of the table types listed under
#'        "Input Types".
#'
#' @param treatment  Each sample's treatment: a vector or factor with one
#'        value per sample (matched to the samples by name when both are
#'        named, else in the same order), or, with `data`, the name of one of
#'        its columns. One of its values is the control (see `control`).
#'
#' @param time  Each sample's time point, given as `treatment` is. A factor's
#'        levels set the order of the time points; numbers are sorted, and
#'        other values are taken in sorted order (give a factor to set their
#'        order). Every time point is a level: time is not treated as a
#'        number.
#'
#' @param data  Sample metadata: a data frame with one row per sample,
#'        holding the columns named in `treatment`, `time`, `strata` or
#'        `plots`. Rows are matched to the samples in `counts` by row name,
#'        and rows for other samples are ignored; without row names, they
#'        must be in the same order as the samples. Default: `NULL`.
#'
#' @param control  The control treatment: one of the values of `treatment`.
#'        Default: `NULL`, the first level of `treatment` (its first value in
#'        sorted order, if not a factor), as in `vegan::prc()`.
#'
#' @param n_axes  The number of principal response curves (axes) to return.
#'        Default: `1`, the method's convention: one curve per treatment.
#'        This is the exception to ecodive's default of 3 axes for
#'        ordinations. Reduced, with a warning, to the number of constrained
#'        axes.
#'
#' @param scaling  How the coefficients and weights are scaled, as in
#'        `vegan`'s `summary()` and `plot()` of a `prc()` result: `vegan`'s
#'        scaling `1` (`'sites'`), `2` (`'species'`) or `3` (`'symmetric'`).
#'        `'form'` and `'covariance'` are accepted for 1 and 2. The product
#'        of a coefficient and a weight is the same in every scaling.
#'        Default: `'symmetric'`, `vegan`'s default for `plot()` (its
#'        `summary()` defaults to `'sites'`).
#'
#' @param permutations  The number of random permutations for the tests, or a
#'        matrix with one permutation of the samples per row (such as
#'        `permute::shuffleSet()` returns). `0` runs no tests. Default: `999`.
#'
#' @param strata  Permute samples only within these groups (blocks): the
#'        name of a column of `data`, or a vector with one value per sample.
#'        `strata = time`, for instance, permutes the samples of each time
#'        point among themselves. Give at most one of `strata`, `plots` and a
#'        `permutations` matrix. Default: `NULL`.
#'
#' @param plots  Permute whole plots (the experimental units, such as the
#'        ponds or ditches sampled at each time point) instead of samples:
#'        the name of a column of `data`, or a vector with one value per
#'        sample. Each permutation shuffles the plots, keeping each plot's
#'        series of samples together, as `permute::how(plots =
#'        Plots(strata = plots, type = "free"), within = Within(type =
#'        "none"))` does in `vegan`. Every plot must have one sample at each
#'        time point. Default: `NULL`.
#'
#' @param seed  The random seed for the permutations: an integer. Default:
#'        `0`.
#'
#' @return A list of class `ecodive_prc`:
#'
#'   * `coefficients`: A data frame with one row per treatment and time
#'     point, control first and in time order within each treatment:
#'     `treatment`, `time` (factors), and one column of canonical
#'     coefficients per axis, `PRC1`, `PRC2`, etc. The control's are 0.
#'     `NA` for a treatment with no samples at that time point, and for
#'     every treatment at a time point with no control samples.
#'     `vegan::summary(prc)$coefficients`, in long form.
#'   * `weights`: The features' weights, one row per feature and one column
#'     per axis: `vegan::summary(prc)$sp`.
#'   * `eig`: Each axis's eigenvalue (a variance).
#'   * `prop`: The proportion of the variance that each axis explains, not
#'     counting the variance that time explains, as in [rda()].
#'   * `test`: The permutation tests: a data frame with rows `PRC1` (the
#'     first axis) and `all` (all constrained axes), and columns `df`,
#'     `variance`, `F`, `z` and `p_value`. `NA` when `permutations = 0`, and
#'     `F` too when the model leaves no residual degrees of freedom.
#'   * `df_residual`, `residual`: The residual degrees of freedom and
#'     variance, the denominator of both tests.
#'   * `F_perm`: The permuted F statistics: one row per permutation and one
#'     column per test.
#'   * `n_perm`: The number of permutations.
#'   * `permutation`: How samples were permuted: `'free'`, `'within
#'     strata'`, `'whole plots'` or `'supplied matrix'`.
#'   * `scaling`: The scaling of `coefficients` and `weights`: `'sites'`,
#'     `'species'` or `'symmetric'`.
#'   * `control`: The control treatment.
#'   * `rda`: The underlying partial RDA, an [rda()] result with `n_axes`
#'     axes, whose signs follow the curves' (see "Signs").
#'
#' @details
#'
#' **Method**
#'
#' With \eqn{Z} the time points' indicators (the conditions) and \eqn{X} the
#' indicators of each treatment other than the control at each time point
#' (the constraints), PRC is the partial RDA of the normalized table on
#' \eqn{X} given \eqn{Z} (van den Brink & ter Braak 1999): `rda(counts,
#' data = d, formula = ~ treatment:time + Condition(time))`, which spans the
#' same model. Its fitted values are each cell's (treatment \eqn{\times}
#' time point's) mean, less the mean of its time point, so the site scores
#' are constant within a cell. A treatment's canonical coefficient at a time
#' point is its cell's site score less the control's at that time point
#' (vegan's \eqn{c_{dt}}), here averaged over the cell's samples (the "wa"
#' scores of [rda()], whose cell means are the fitted scores). The features'
#' weights \eqn{b_k} are the species scores, so that
#' \eqn{b_k c_{dt}} approximates the change in feature \eqn{k}'s
#' (normalized) abundance under treatment \eqn{d} at time \eqn{t}, relative
#' to the control.
#'
#' Use `norm` to analyze transformed abundances, such as `norm = 'clr'` (with
#' a `pseudocount`) or `norm = 'hellinger'`. van den Brink & ter Braak
#' analyzed \eqn{\log(10 x + 1)} counts; transform such values yourself and
#' give them with `norm = 'none'`, the default.
#'
#'
#' **Permutation tests**
#'
#' Both tests are those of `vegan::anova.cca()` under the reduced model:
#'
#' * all axes: \eqn{F = (\mathrm{SS}_X / q) / (\mathrm{RSS} / \nu)},
#'   for \eqn{q} constraint columns and \eqn{\nu = n - 1 - \mathrm{rank}(Z) -
#'   q} residual degrees of freedom. This is PERMANOVA of the Euclidean
#'   distances between the normalized samples ([permanova()] of
#'   [euclidean()]), with the indicators of the treatments other than the
#'   control at each time point as constraints and time as the condition:
#'   the columns of `model.matrix(~ time + time:treatment)` after time's.
#'   (`formula = ~ treatment:time + Condition(time)` spans the same model,
#'   so it has the same F, but it keeps other columns of the aliased
#'   model matrix, whose permutations give other permuted F values.)
#' * the first axis (`first = TRUE` in `vegan`): \eqn{F = \lambda_1 /
#'   (\mathrm{RSS} / \nu)}, with \eqn{\lambda_1} the first constrained
#'   eigenvalue and RSS the residual of the full model.
#'
#' Each permutation shuffles the rows of the constraints (equivalently, the
#' samples), removes the conditions from them, and refits the model, so that
#' the p-values equal `vegan`'s for the same permutations. Samples sampled
#' repeatedly are not independent, so permute whole plots (`plots`) or
#' within time points (`strata`) rather than freely; `vegan`'s own
#' restricted designs (such as `Within(type = "series")`) can be given as a
#' `permutations` matrix from `permute::shuffleSet()`.
#'
#'
#' **Signs**
#'
#' Each axis is flipped so that its coefficient with the largest absolute
#' value is positive (the first such, if tied), and the features' weights
#' and the `rda` element's scores, loadings, biplot scores and centroids are
#' flipped with it. So the curves mostly rise from the control, and their
#' signs may be the reverse of `vegan`'s.
#'
#'
#' **Matching vegan**
#'
#' `coefficients` and `weights` are `vegan`'s `summary(prc(...), scaling =
#' scaling)` (up to each axis's sign), and `eig` its eigenvalues. Where a
#' time point has no control samples, `vegan::prc()` drops some other
#' treatment's column and reports coefficients relative to that treatment;
#' here they are `NA`, with a warning. The test statistics and p-values are
#' those of `anova(prc, first = TRUE)` and `anova(prc)` for the same
#' permutations.
#'
#' @references
#' van den Brink, P. J., & ter Braak, C. J. F. (1999). Principal response
#' curves: analysis of time-dependent multivariate responses of biological
#' community to stress. *Environmental Toxicology and Chemistry*, 18(2),
#' 138-148. \doi{10.1002/etc.5620180207}
#'
#' Legendre, P., Oksanen, J., & ter Braak, C. J. F. (2011). Testing the
#' significance of canonical axes in redundancy analysis. *Methods in
#' Ecology and Evolution*, 2(3), 269-277.
#' \doi{10.1111/j.2041-210X.2010.00078.x}
#'
#' @inheritSection documentation Input Types
#'
#' @seealso [rda()], [permanova()], [plot.ecodive_prc()]
#'
#' @export
#' @examples
#'     # Four treatments and a control, each in two plots, at five times.
#'     set.seed(1)
#'     meta <- expand.grid(plot = 1:10, time = c(-1, 1, 2, 4, 8))
#'     meta$dose <- rep(c(0, 0.1, 1, 10, 100), 2)[meta$plot]
#'     effect <- outer(log1p(meta$dose) * (meta$time > 0) / (1 + meta$time / 4), 1:12)
#'     counts <- matrix(rpois(50 * 12, exp(2 + 0.05 * outer(meta$time, 1:12 - 6) - 0.05 * effect)), 50)
#'     colnames(counts) <- paste0('Sp', 1:12)
#'     rownames(counts) <- rownames(meta) <- paste0('S', 1:50)
#'
#'     prc <- principal_response_curves(
#'       log1p(counts), 'dose', 'time', meta, plots = 'plot', permutations = 99 )
#'     prc
#'     plot(prc)
#'
#'     # The coefficients, ready for ggplot2.
#'     head(prc$coefficients)
#'
principal_response_curves <- function (
    counts, treatment, time, data = NULL, control = NULL, margin = 1L,
    norm = 'none', pseudocount = NULL, n_axes = 1L, scaling = 'symmetric',
    permutations = 999L, strata = NULL, plots = NULL, seed = 0L, cpus = n_cpus() ) {

  if (missing(treatment)) stop('`treatment` is required: one value per sample, or a column name of `data`.', call. = FALSE)
  if (missing(time))      stop('`time` is required: one value per sample, or a column name of `data`.', call. = FALSE)

  norm_arg <- norm

  validate_counts()
  validate_norm()
  validate_pseudocount()
  validate_margin()
  validate_cpus()
  validate_seed()
  k <- ord_k(n_axes)

  scaling <- ord_scaling_arg(scaling)
  if (identical(scaling, 'contribution'))
    stop(
      "`scaling` must be vegan's scaling 1, 2 or 3 ('sites', 'species' or ",
      "'symmetric').", call. = FALSE )
  scaling <- c(form = 'sites', covariance = 'species', symmetric = 'symmetric')[[scaling]]


  #--------------------------------------------------------
  # Samples, and each one's treatment, time, stratum, plot.
  #--------------------------------------------------------
  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )

  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  n      <- as.integer(dims[[1]])
  labels <- dn[[1]]

  if (n < 3)
    stop('`counts` must have at least 3 samples.', call. = FALSE)

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, rn, nrow(data), ref = '`counts`')
  }

  value <- function (v, arg) {
    v <- per_sample(v, arg, data, idx, labels, n, ref = '`counts`')
    if (is.null(v))
      stop('`', arg, '` is NULL (a misspelled column, such as meta$tretment?).', call. = FALSE)
    if (anyNA(v))
      stop('`', arg, '` must not contain NAs.', call. = FALSE)
    v
  }

  trt <- droplevels(as.factor(value(treatment, 'treatment')))
  tim <- droplevels(as.factor(value(time, 'time')))

  if (is.null(control)) {
    control <- levels(trt)[[1]]
  }
  else if (!(is.atomic(control) && length(control) == 1 && !is.na(control) &&
             as.character(control) %in% levels(trt))) {
    stop('`control` must be one of the treatments: ', name_list(levels(trt)), '.', call. = FALSE)
  }
  control <- as.character(control)

  if (nlevels(trt) < 2)
    stop('`treatment` must have at least two values: the control and a treatment.', call. = FALSE)

  st <- per_sample(strata, 'strata', data, idx, labels, n, ref = '`counts`')
  pl <- per_sample(plots,  'plots',  data, idx, labels, n, ref = '`counts`')
  if (!is.null(pl) && anyNA(pl)) stop('`plots` must not contain NAs.', call. = FALSE)

  pm <- perm_matrix(permutations, n)
  if (sum(!is.null(st), !is.null(pl), !is.null(pm$perms)) > 1)
    stop('Give only one of `strata`, `plots` and a `permutations` matrix.', call. = FALSE)


  #--------------------------------------------------------
  # The partial RDA. Its formula's model matrix has more
  # columns than its rank; the aliased ones are dropped
  # quietly.
  #--------------------------------------------------------
  d <- data.frame(treatment = stats::relevel(trt, control), time = tim)
  if (!is.null(labels) && !anyNA(labels) && !anyDuplicated(labels)) rownames(d) <- labels

  ord <- withCallingHandlers(
    rda(
      counts, data = d,
      formula = if (nlevels(tim) > 1) ~ treatment:time + Condition(time) else ~ treatment,
      margin = margin, norm = norm_arg, pseudocount = pseudocount,
      n_axes = k, cpus = cpus ),
    message = function (m) {
      if (startsWith(conditionMessage(m), 'Dropped model columns')) invokeRestart('muffleMessage')
    },
    warning = function (w) {
      if (grepl('eigenvalues are > 0', conditionMessage(w), fixed = TRUE)) invokeRestart('muffleWarning')
    })

  n_con <- length(attr(ord, 'constrained_eig'))
  if (n_con == 0)
    stop('The treatments explain none of the variation left by time: there are no response curves.', call. = FALSE)
  if (k > n_con) {
    warning(sprintf('`n_axes` reduced to %d, the number of constrained axes.', n_con), call. = FALSE)
    k <- n_con
  }
  if (ncol(ord) > k) ord <- prc_keep_axes(ord, k)

  axes <- sprintf('PRC%d', seq_len(k))


  #--------------------------------------------------------
  # Coefficients: each cell's mean site score less the
  # control's at that time point.
  #--------------------------------------------------------
  cell   <- interaction(trt, tim, drop = FALSE, lex.order = TRUE, sep = '\r')
  tlev   <- levels(trt)
  wlev   <- levels(tim)
  rows   <- expand.grid(time = wlev, treatment = c(control, setdiff(tlev, control)),
                        stringsAsFactors = FALSE)[, 2:1]
  g      <- as.integer(cell)
  cnt    <- tabulate(g, nlevels(cell))
  sums   <- rowsum(unclass(ord)[, seq_len(k), drop = FALSE], g, reorder = TRUE)
  have   <- as.integer(rownames(sums))
  means  <- matrix(NA_real_, nlevels(cell), k, dimnames = list(levels(cell), NULL))
  means[have, ] <- sums / cnt[have]

  key    <- function (t, w) paste(t, w, sep = '\r')
  cf     <- means[key(rows$treatment, rows$time), , drop = FALSE] -
            means[key(control, rows$time), , drop = FALSE]
  cf[rows$treatment == control & !is.na(cf[, 1]), ] <- 0

  no_ctrl <- wlev[cnt[match(key(control, wlev), levels(cell))] == 0]
  if (length(no_ctrl))
    warning(
      'No control samples at time ', name_list(no_ctrl),
      ': its coefficients are NA.', call. = FALSE )


  #--------------------------------------------------------
  # Signs: each axis's largest coefficient is positive.
  #--------------------------------------------------------
  flip <- vapply(seq_len(k), function (a) {
    v <- cf[, a]
    i <- which.max(abs(replace(v, is.na(v), 0)))
    if (length(i) && isTRUE(v[[i]] < 0)) -1 else 1
  }, numeric(1))
  cf  <- cf * rep(flip, each = nrow(cf))
  ord <- prc_flip(ord, flip)


  #--------------------------------------------------------
  # vegan's scalings, from the stored scaling 2.
  #--------------------------------------------------------
  eig   <- attr(ord, 'eig')[seq_len(k)]
  f     <- sqrt(eig / attr(ord, 'trace'))
  mul   <- switch(scaling, sites = f, species = rep(1, k), symmetric = sqrt(f))
  cf    <- cf * rep(mul, each = nrow(cf))
  wts   <- attr(ord, 'loadings')[, seq_len(k), drop = FALSE] / rep(mul, each = nrow(attr(ord, 'loadings')))
  colnames(wts) <- axes

  coefficients <- data.frame(
    treatment = factor(rows$treatment, levels = tlev),
    time      = factor(rows$time, levels = wlev),
    stringsAsFactors = FALSE )
  for (a in seq_len(k)) coefficients[[axes[[a]]]] <- unname(cf[, a])


  #--------------------------------------------------------
  # Permutation tests.
  #--------------------------------------------------------
  tst <- prc_test(
    counts, margin, norm, pseudocount, trt, tim, control, st, pl, pm, seed, cpus )

  out <- list(
    coefficients = coefficients,
    weights      = wts,
    eig          = stats::setNames(eig, axes),
    prop         = stats::setNames(attr(ord, 'prop')[seq_len(k)], axes),
    test         = tst$test,
    df_residual  = tst$df_residual,
    residual     = tst$residual,
    F_perm       = tst$F_perm,
    n_perm       = tst$n_perm,
    permutation  = tst$permutation,
    scaling      = scaling,
    control      = control,
    rda          = ord )

  class(out) <- 'ecodive_prc'
  out
}



# The first k axes of an rda() result, with its per-axis attributes.
prc_keep_axes <- function (ord, k) {
  at  <- attributes(ord)
  out <- unclass(ord)[, seq_len(k), drop = FALSE]
  for (nm in setdiff(names(at), c('dim', 'dimnames')))
    attr(out, nm) <- at[[nm]]
  attr(out, 'eig')  <- at$eig[seq_len(k)]
  attr(out, 'prop') <- at$prop[seq_len(k)]
  for (nm in c('biplot', 'centroids', 'loadings'))
    if (!is.null(at[[nm]])) attr(out, nm) <- at[[nm]][, seq_len(k), drop = FALSE]
  class(out) <- class(ord)
  out
}


# Flips the axes of an rda() result whose `flip` is -1, with its biplot
# scores, centroids and loadings.
prc_flip <- function (ord, flip) {
  k <- length(flip)
  if (all(flip == 1)) return (ord)
  at <- attributes(ord)
  for (a in which(flip < 0)) ord[, a] <- -ord[, a]
  for (nm in c('biplot', 'centroids', 'loadings')) {
    v <- at[[nm]]
    if (is.null(v)) next
    v[, seq_len(k)] <- v[, seq_len(k), drop = FALSE] * rep(flip, each = nrow(v))
    attr(ord, nm) <- v
  }
  ord
}



# The tests of the first axis and of all axes, with C_prc_perm().
prc_test <- function (
    counts, margin, norm, pseudocount, trt, tim, control, strata, plots, pm, seed, cpus ) {

  n <- length(trt)

  # Conditions: time; constraints: each non-control cell.
  z <- if (nlevels(tim) > 1) stats::model.matrix(~ tim)[, -1, drop = FALSE] else NULL
  if (!is.null(z)) z <- z - rep(colMeans(z), each = n)

  cells <- interaction(trt, tim, drop = TRUE, lex.order = TRUE, sep = '\r')
  lev   <- levels(cells)
  lev   <- lev[!startsWith(lev, paste0(control, '\r'))]
  x     <- vapply(lev, function (l) as.numeric(cells == l), numeric(n))
  x     <- matrix(x, n)
  x     <- x - rep(colMeans(x), each = n)

  qz   <- condition_basis(z)
  rz   <- if (is.null(qz)) 0L else ncol(qz)
  full <- design_basis(x, qz)
  q    <- full$rank

  if (!is.null(qz)) {
    qz <- unname(qz[, , drop = FALSE])
    attributes(qz) <- list(dim = dim(qz))
  }
  cols <- unname(full$cols)
  attributes(cols) <- list(dim = dim(cols))
  storage.mode(cols) <- 'double'

  permutation <- if (!is.null(pm$perms)) 'supplied matrix'
                 else if (!is.null(strata)) 'within strata'
                 else if (!is.null(plots)) 'whole plots'
                 else 'free'

  if (!is.null(plots)) pm <- prc_plot_perms(plots, tim, pm$n_perm, seed)

  res <- .Call(
    C_prc_perm, counts, margin, norm, pseudocount, cols, qz, pm$perms,
    strata_codes(strata), as.integer(c(pm$n_perm, seed, cpus)) )

  stats  <- res$stats
  df_res <- n - 1L - rz - q
  rss    <- res$total - res$ss_z - stats[2, ]
  den    <- if (df_res > 0) rss / df_res else rep(NA_real_, ncol(stats))
  f1     <- stats[1, ] / den
  fall   <- (stats[2, ] / q) / den

  F_perm <- cbind(PRC1 = f1[-1], all = fall[-1])
  stat   <- c(f1[[1]], fall[[1]])
  test   <- data.frame(
    df       = c(1L, as.integer(q)),
    variance = c(stats[1, 1], stats[2, 1]) / (n - 1),
    F        = stat,
    z        = if (pm$n_perm > 0) perm_z(stat, F_perm) else NA_real_,
    p_value  = if (pm$n_perm > 0) perm_p(stat, F_perm) else NA_real_,
    row.names = c('PRC1', 'all') )

  list(
    test        = test,
    df_residual = as.integer(df_res),
    residual    = rss[[1]] / (n - 1),
    F_perm      = F_perm,
    n_perm      = pm$n_perm,
    permutation = permutation )
}



# Permutations of whole plots: each plot's series moves to another plot,
# keeping its time points. Returns perm_matrix()'s list.
prc_plot_perms <- function (plots, tim, n_perm, seed) {

  pf <- as.factor(plots)
  np <- nlevels(pf)
  nt <- nlevels(tim)
  n  <- length(pf)

  pos <- matrix(NA_integer_, np, nt)
  ok  <- n == np * nt && !anyDuplicated(cbind(as.integer(pf), as.integer(tim)))
  if (ok) pos[cbind(as.integer(pf), as.integer(tim))] <- seq_len(n)
  if (!ok || anyNA(pos))
    stop(
      '`plots` must give each plot one sample at each time point, so that ',
      'whole plots can be permuted.', call. = FALSE )
  if (np < 2)
    stop('`plots` must have at least two plots to permute.', call. = FALSE)

  if (n_perm == 0) return (list(n_perm = 0L, perms = NULL))

  pp    <- .Call(C_perm_indices, np, n_perm, NULL, seed) # np x n_perm, 1-based
  perms <- matrix(0L, n, n_perm)
  for (k in seq_len(n_perm))
    perms[pos, k] <- pos[pp[, k], ] - 1L

  list(n_perm = as.integer(n_perm), perms = perms)
}



#' @export
print.ecodive_prc <- function (x, digits = 4L, ...) {

  ord  <- x$rda
  cf   <- x$coefficients
  k    <- ncol(x$weights)
  iner <- attr(ord, 'inertia')
  tr   <- attr(ord, 'trace')
  sym  <- print_sym()
  fmt  <- function (v) format(signif(v, digits))
  pct  <- function (v) sprintf('%.1f%%', 100 * v)

  trts <- levels(cf$treatment)
  cat(sprintf(
    'Principal Response Curves: %d samples, %d treatments (control %s), %d time points\n',
    nrow(ord), length(trts) - 1L, sQuote(x$control, FALSE), nlevels(cf$time) ))
  cat(sprintf(
    'Variance: time %s (%s), treatment %s time %s (%s), residual %s (%s)\n',
    fmt(iner[['conditional']]), pct(iner[['conditional']] / tr), sym$times,
    fmt(iner[['constrained']]), pct(iner[['constrained']] / tr),
    fmt(iner[['unconstrained']]), pct(iner[['unconstrained']] / tr) ))
  cat(sprintf(
    '%s: eigenvalue %s, %s of the variance within time points\n',
    names(x$eig), vapply(x$eig, fmt, ''), pct(x$prop) ), sep = '')
  cat(sprintf('Scaling: %s\n', x$scaling))

  cat('\n')
  if (x$n_perm > 0) {
    cat(sprintf('Permutation tests (%d permutations, %s):\n', x$n_perm, x$permutation))
    num <- function (v) vapply(v, function (a) if (is.na(a)) '' else format(signif(a, digits)), '')
    tab <- data.frame(
      df       = c(x$test$df, x$df_residual),
      variance = num(c(x$test$variance, x$residual)),
      F        = num(c(x$test$F, NA)),
      p_value  = num(c(x$test$p_value, NA)),
      row.names = c('PRC1', 'All axes', 'Residual'),
      stringsAsFactors = FALSE )
    print(tab, right = TRUE)
  }
  else {
    cat('No permutation tests (permutations = 0).\n')
  }

  for (a in seq_len(min(k, 2L))) {
    nm  <- colnames(x$weights)[[a]]
    sub <- cf[cf$treatment != x$control, , drop = FALSE]
    m   <- matrix(
      sub[[nm]], ncol = nlevels(cf$time), byrow = TRUE,
      dimnames = list(as.character(unique(sub$treatment)), levels(cf$time)) )
    cat(sprintf('\nCoefficients on %s (treatment %s time):\n', nm, sym$times))
    print(round(m, digits))

    w   <- x$weights[, a]
    top <- order(-abs(w))[seq_len(min(8L, length(w)))]
    lab <- if (is.null(names(w))) as.character(top) else names(w)[top]
    cat(sprintf(
      'Largest weights on %s: %s%s\n', nm,
      paste(sprintf('%s %s', lab, formatC(w[top], digits = 3, format = 'f')), collapse = ', '),
      if (length(w) > length(top)) paste0(', ', sym$more) else '' ))
  }
  if (k > 2) cat(sprintf('(and %d more axes)\n', k - 2L))

  invisible(x)
}



#' Plot principal response curves
#'
#' The PRC diagram: each treatment's canonical coefficients over time, with
#' the control's as the flat line at 0, and the features' weights as ticks
#' on a scale at the right, labelled for those with the largest weights.
#' A feature's weight times a treatment's coefficient approximates the
#' change in its (normalized) abundance relative to the control.
#'
#' @param x  A result of [principal_response_curves()].
#'
#' @param axis  Which curve (axis) to draw. Default: `1`.
#'
#' @param species  How many features' weights to label: the ones with the
#'        largest absolute weights. Every feature's weight is drawn as a tick
#'        unless `species = FALSE`. `TRUE` labels them all. Default: `10`.
#'
#' @param legend  Where to put the legend: a position keyword for
#'        [graphics::legend()], or `NULL` for none. Default: `'topleft'`.
#'
#' @param col  The treatments' colours, in the order of their levels, the
#'        control's included. Default: `NULL`, grey for the control and a
#'        palette for the others.
#'
#' @param ...  Passed to [graphics::plot.default()], such as `main`, `xlab`,
#'        `ylab`, `sub`, `lwd` or `lty`.
#'
#' @details
#' Time points are placed at their values when every level of `time` is a
#' number, and evenly otherwise. The title is "Principal Response Curves";
#' the subtitle gives the axis's share of the variance within time points
#' and, when the tests were run, the first axis's p-value.
#'
#' @return `x`, invisibly.
#'
#' @seealso [principal_response_curves()]
#'
#' @export
#' @examples
#'     set.seed(1)
#'     meta <- expand.grid(plot = 1:10, time = c(-1, 1, 2, 4, 8))
#'     meta$dose <- rep(c(0, 0.1, 1, 10, 100), 2)[meta$plot]
#'     effect <- outer(log1p(meta$dose) * (meta$time > 0) / (1 + meta$time / 4), 1:12)
#'     counts <- matrix(rpois(50 * 12, exp(2 + 0.05 * outer(meta$time, 1:12 - 6) - 0.05 * effect)), 50)
#'     colnames(counts) <- paste0('Sp', 1:12)
#'
#'     prc <- principal_response_curves(log1p(counts), meta$dose, meta$time, permutations = 0)
#'     plot(prc, species = 5)
#'
plot.ecodive_prc <- function (
    x, axis = 1L, species = 10L, legend = 'topleft', col = NULL, ... ) {

  k <- ncol(x$weights)
  if (!is.numeric(axis) || length(axis) != 1 || !is.finite(axis) ||
      axis != round(axis) || axis < 1 || axis > k)
    stop(sprintf('`axis` must be an axis number from 1 to %d.', k), call. = FALSE)
  axis <- as.integer(axis)
  pos  <- plot_legend_pos(legend)
  n_sp <- plot_count(species, 'species')
  ticks <- !isFALSE(species)

  cf   <- x$coefficients
  nm   <- colnames(x$weights)[[axis]]
  tlev <- levels(cf$treatment)
  wlev <- levels(cf$time)

  xt <- suppressWarnings(as.numeric(wlev))
  if (anyNA(xt)) xt <- seq_along(wlev)
  xv <- xt[as.integer(cf$time)]
  yv <- cf[[nm]]

  if (is.null(col)) {
    col <- rep(PLOT_MUTED, length(tlev))
    oth <- which(tlev != x$control)
    col[oth] <- plot_palette(length(oth))
  }
  col <- rep_len(col, length(tlev))

  w    <- x$weights[, axis]
  xr   <- range(xt)
  span <- if (diff(xr) > 0) diff(xr) else 1
  xlim <- c(xr[[1]], xr[[2]] + if (ticks) 0.3 * span else 0)
  yr   <- range(c(0, yv, if (ticks) w), finite = TRUE)

  sub <- sprintf('%s explains %.1f%% of the variance within time points', nm, 100 * x$prop[[axis]])
  if (x$n_perm > 0 && axis == 1L)
    sub <- sprintf('%s; p = %s', sub, plot_p(x$test['PRC1', 'p_value']))

  args <- list(...)
  plot_do(
    graphics::plot.default, c(list(NA, NA), plot_drop(args, c('lwd', 'lty'))),
    list(
      type = 'n', xlim = xlim, ylim = yr, xaxt = 'n',
      xlab = 'Time', ylab = sprintf('Canonical coefficient (%s)', nm),
      main = 'Principal Response Curves', sub = sub, las = 1 ))

  if (plot_x_axis(args)) graphics::axis(1, at = xt, labels = wlev)
  graphics::abline(h = 0, col = 'grey85')

  lwd <- plot_arg(args, 'lwd', 1.5)
  lty <- plot_arg(args, 'lty', 1)
  for (i in seq_along(tlev)) {
    j <- which(cf$treatment == tlev[[i]])
    j <- j[order(xv[j])]
    graphics::lines(xv[j], yv[j], col = col[[i]], lwd = lwd, lty = lty)
    graphics::points(xv[j], yv[j], col = col[[i]], pch = 16, cex = 0.7)
  }

  # The weights: ticks on a scale at the right, the largest labelled.
  if (ticks && length(w)) {
    x0 <- xr[[2]] + 0.06 * span
    graphics::segments(x0, w, x0 + 0.03 * span, w, col = PLOT_MUTED)
    graphics::segments(x0, yr[[1]], x0, yr[[2]], col = 'grey85')
    graphics::text(x0, graphics::par('usr')[[4]], 'Weights', pos = 1, cex = 0.7, col = 'grey40')
    top <- order(-abs(w))[seq_len(min(n_sp, length(w)))]
    if (length(top)) {
      lab <- if (is.null(names(w))) as.character(top) else names(w)[top]
      ly  <- prc_spread(w[top], 1.1 * graphics::strheight('M', cex = 0.7))
      graphics::segments(x0 + 0.03 * span, w[top], x0 + 0.06 * span, ly, col = 'grey70')
      plot_label(x0 + 0.06 * span, ly, lab, col = 'grey20', pos = 4, cex = 0.7)
    }
  }

  plot_legend(pos, tlev, col, 16, 'Treatment', lty = 1)

  invisible(x)
}


# Label positions near y, at least `gap` apart, in the same order.
prc_spread <- function (y, gap) {
  o  <- order(y)
  v  <- y[o]
  if (length(v) > 1)
    for (i in 2:length(v)) if (v[[i]] - v[[i - 1]] < gap) v[[i]] <- v[[i - 1]] + gap
  # Centre the shifted block on the original values.
  v <- v - (mean(v) - mean(y))
  out <- numeric(length(y))
  out[o] <- v
  out
}
