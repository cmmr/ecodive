# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Sensitivity of a Differential Abundance Analysis
#'
#' Reruns a differential-abundance (DA) method of [da_test()] over the
#' settings that have no right value (the pseudocount, the treatment of
#' zeros, the assumed variation in total load) and reports, feature by
#' feature, whether its decision holds in every run.
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**.
#'
#' @param method   The name of a DA method, as for [da_test()]: `"ancom"`,
#'        `"ancombc"`, `"ancombc2"`, `"da_lm"`, `"da_ttest"`,
#'        `"da_wilcoxon"`, `"dirmult_ttest"`, or `"linda"`. Flexible
#'        matching is supported.
#'
#' @param ...   The base run: everything [da_test()] takes after `method`
#'        (`groups`, `data`, `formula`) and the method's own arguments
#'        (`alpha`, `p_adjust`, `norm`, `pseudocount`, `scale_sd`, `seed`,
#'        `cpus`, ...). Every run is this call with one setting changed.
#'
#' @param pseudocounts  Pseudocounts to try: `TRUE` for `c(0.1, 0.5, 1)`
#'        (ANCOMBC's), or a vector of positive numbers. Each run sets the
#'        method's `pseudocount`. Default: `NULL` (none).
#'
#' @param zeros  Zero replacements to try: `TRUE` for every method of
#'        [list_zero_replacements()], or their names. Each run gives the
#'        method `replace_zeros(counts, z)` instead of `counts`, with the
#'        method's default `pseudocount` (which adds nothing to a table
#'        without zeros). Not for the methods that model counts
#'        (`zero_free_ok = FALSE` in [list_da_tests()]: [dirmult_ttest()]).
#'        Default: `NULL` (none).
#'
#' @param scale_sds  Scale uncertainties to try: `TRUE` for
#'        `c(0, 0.25, 0.5, 1)`, or a vector of non-negative numbers. Each
#'        run sets the method's `scale_sd`, so only for methods that have
#'        one ([dirmult_ttest()]). Default: `NULL` (none).
#'
#' @param crossed  `FALSE` (the default) to vary one setting at a time
#'        around the base run; `TRUE` for every combination of the
#'        settings. The pseudocounts and the zero replacements are
#'        alternative treatments of the zeros, so they are not crossed with
#'        each other: each of them is crossed with each `scale_sds`.
#'
#' @return A data frame of class `c('ecodive_da_sensitivity',
#'   'data.frame')`, one row per feature (and covariate, for the model
#'   methods) and run, with columns:
#'
#'   * `feature`: the feature's name, as in the method's result.
#'   * `covariate`: for the model methods ([ancombc()], [ancombc2()],
#'     [da_lm()]), the column of the design; every covariate but the
#'     intercept is kept. `NA` for the group methods.
#'   * `run`: the run's number, `0` for the base run.
#'   * `parameter`: what the run changed: `'base'`, `'pseudocount'`,
#'     `'zeros'`, `'scale_sd'`, or, with `crossed = TRUE`, `'crossed'` for
#'     runs that change more than one.
#'   * `setting`: the run's setting in words, e.g. `"pseudocount = 0.5"`,
#'     `"zeros = mult"`, `"pseudocount = 0.5, scale_sd = 1"`, or `"base"`.
#'   * `effect`: the method's effect size, the column of its result named
#'     by `effect` in [list_da_tests()] (`log_fc`, `log2_fc`,
#'     `median_diff`, `diff`, `estimate`); `NA` for [ancom()], which has
#'     none.
#'   * `p_adjusted`: the method's adjusted p-value (`NA` for [ancom()]).
#'   * `reject`: the method's own decision, its `reject` column (for
#'     ANCOM-BC, not `reject_robust`).
#'
#'   Attributes:
#'
#'   * `method`: the method's id, e.g. `'ancombc2'`.
#'   * `effect_name`: the name of the method's effect column (`NA` for
#'     ANCOM).
#'   * `alpha`: the base run's significance level.
#'   * `log_base`: the base of the effects' logs (`exp(1)`, `2`, or `NA`
#'     when the effects are not logs).
#'   * `covariate`: for the model methods, the covariate of the term of
#'     interest (the first column of `groups`'s term, else the first
#'     covariate that isn't the intercept), which `plot()` shows by default;
#'     `NULL` for the group methods.
#'   * `runs`: a data frame of the runs, with columns `run`, `parameter`,
#'     `pseudocount`, `zeros`, `scale_sd` (`NA` where the run keeps the base
#'     run's) and `setting`.
#'   * `fits`: the method's result of every run, a list in run order (the
#'     base run first). This can take much memory: one full result per
#'     run.
#'   * `summary`: a data frame with a row per feature (and covariate):
#'     `feature`, `covariate`, `n_runs` (the runs with a result for it),
#'     `n_reject` (those that reject it), `prop_reject` (their proportion),
#'     `passed` (`TRUE` when every run makes the same decision, rejected or
#'     not), `robust` (`passed` and rejected in the base run: the decisions
#'     to report), and `sign_stable` (`TRUE` when the effect has the same
#'     sign in every run where it is finite; `NA` without effects).
#'
#' @details
#'
#' **Why.** A DA result can hinge on a setting with no right value. The
#' ANCOM-BC2 authors strongly recommend their pseudocount sensitivity
#' analysis (Lin & Peddada 2024), and benchmarks find ANCOM-BC2's results
#' among the least replicable without it (Pelto et al. 2025). Assumptions
#' about the samples' total microbial load (the "scale") can break the
#' control of false discoveries in common tools (McGovern et al. 2023; Dos
#' Santos et al. 2026), and Nixon et al. (2025) recommend testing over a
#' range of scale uncertainties. A finding that holds over these settings is
#' firmer than one that holds at a single setting.
#'
#' **The runs.** The base run is `da_test(counts, method, ...)`. Each other
#' run is the same call with one setting changed (or, with `crossed = TRUE`,
#' one combination of settings). Every run uses the base run's `seed`, so
#' that the runs differ by their settings and not by random draws, and its
#' `cpus`; the results are identical for any `cpus`. Runs other than the
#' base run of [ancombc()] and [ancombc2()] are made with `sensitivity =
#' FALSE`. Warnings from the runs are given once each, and not when the base
#' run already gave them.
#'
#' **What each setting means.**
#'
#' * `pseudocounts`: what `pseudocount` means differs by method. For
#'   [ancom()], [da_wilcoxon()], [da_ttest()] and [da_lm()] it is added to
#'   every value before the logs, to patch the zeros (for the last three,
#'   only with `norm = 'clr'` or `'log'`, the norms that take logs). For
#'   [ancombc()]
#'   and [ancombc2()] it is added to every value, as in ANCOMBC's own
#'   sensitivity analysis. For [dirmult_ttest()] it is the Dirichlet
#'   prior, added to the counts before the posterior draws. For [linda()]
#'   it is added to every count when the zeros get a pseudocount
#'   (`zero_handling = 'pseudocount'`, or `'adaptive'` choosing it); when
#'   they are imputed, every run is the same.
#' * `zeros`: replaces the zeros with [replace_zeros()] instead of a
#'   pseudocount, for the methods that take any table of positive values.
#' * `scale_sds`: [dirmult_ttest()]'s scale model, the standard deviation of
#'   the samples' log total load (ALDEx2's gamma): `0` assumes the loads are
#'   the same, larger values allow more variation and give wider intervals.
#'
#' **Relation to ANCOMBC's `passed_ss`.** `summary$passed` is ANCOMBC's
#' `passed_ss` generalised to any method and setting: with `method =
#' 'ancombc2'` and `pseudocounts = TRUE`, it is exactly [ancombc2()]'s
#' `passed_ss` (and `robust` its `reject_robust`). These pseudocount runs
#' are the fits behind `passed_ss`: when the base run already made them
#' (`sensitivity = TRUE`, [ancombc2()]'s default), they are taken from it
#' rather than fitted again. The ANCOMBC2 vignette's "sensitivity score"
#' (the proportion of the pseudocount runs whose adjusted p-value is above
#' `alpha`) is `1 - prop_reject` over those runs.
#'
#' **Plots.** `plot(x)` draws one of two displays (see
#' [plot.ecodive_da_sensitivity()]): `"tiles"` (the default) shows each
#' feature's decision (rows) in every run (columns), coloured by the sign of
#' the effect where the run rejects it, with the features whose decision
#' never changes in bold; `"traces"` shows each feature's effect against one
#' setting's values, with filled points where the run rejects it (for
#' `scale_sds`, the display of ALDEx2's `plotGamma()`).
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
#' microbiomes with covariate adjustment and repeated measures. *Nature
#' Methods*, 21, 83-91. \doi{10.1038/s41592-023-02092-7}
#'
#' Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating
#' scale uncertainty in microbiome and gene expression analysis as an
#' extension of normalization. *Genome Biology*, 26, 139.
#' \doi{10.1186/s13059-025-03609-3}
#'
#' Dos Santos, S. J., Murariu, A. C., Silverman, J. D., & Gloor, G. B.
#' (2026). Not every gene is special: modelling scale controls the false
#' discovery rate. *PLoS Computational Biology*, 22(9), e1014728.
#' \doi{10.1371/journal.pcbi.1014728}
#'
#' McGovern, K. C., Nixon, M. P., & Silverman, J. D. (2023). Addressing
#' erroneous scale assumptions in microbe and gene set enrichment analysis.
#' *PLoS Computational Biology*, 19(11), e1011659.
#' \doi{10.1371/journal.pcbi.1011659}
#'
#' Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary
#' methods provide more replicable results in microbial differential
#' abundance analysis. *Briefings in Bioinformatics*, 26, bbaf130.
#' \doi{10.1093/bib/bbaf130}
#'
#' @seealso [da_test()], [replace_zeros()], [plot.ecodive_da_sensitivity()];
#'   [ancombc2()]'s `sensitivity`.
#'
#' @export
#' @examples
#'     x <- rbind(
#'       c( 20, 110, 100, 101, 100, 103, 104),
#'       c( 33, 110, 120, 100, 101, 100, 102),
#'       c( 12, 110, 100, 110, 100,  50,  90),
#'       c(202, 201,   9,  10,  10,  11,  11),
#'       c(200, 202,  10,  10,  13,  10,  10),
#'       c(203, 201,  14,  10,  10,  13,  12) )
#'     colnames(x) <- paste0('b', 1:7)
#'     grp <- rep(c('treatment', 'placebo'), each = 3)
#'
#'     # Does the t-test's answer depend on the pseudocount?
#'     res <- da_sensitivity(x, 'da_ttest', grp, pseudocounts = TRUE, cpus = 1)
#'     res
#'     attr(res, 'summary')
#'     plot(res)
#'
#'     # The Dirichlet-multinomial t-test over scale uncertainties (gamma).
#'     res <- da_sensitivity(
#'       x, 'dirmult_ttest', grp, scale_sds = TRUE, n_draws = 32, cpus = 1 )
#'     plot(res, 'traces')
#'
da_sensitivity <- function (
    counts, method, ..., pseudocounts = NULL, zeros = NULL,
    scale_sds = NULL, crossed = FALSE ) {

  method <- match_da_test(method)
  id     <- method$id
  args   <- da_sens_args(list(...), method)

  if (!(isTRUE(crossed) || isFALSE(crossed)))
    stop('`crossed` must be TRUE or FALSE.', call. = FALSE)

  grids <- da_sens_grids(method, args, pseudocounts, zeros, scale_sds)
  runs  <- da_sens_runs(grids, crossed)
  v2    <- id %in% c('ancombc', 'ancombc2')

  # The default pseudocounts are on the scale of counts.
  input <- NULL
  if (length(grids$zeros) || isTRUE(pseudocounts))
    input <- da_sens_input(counts, if (is.null(args$margin)) 1L else args$margin)
  if (isTRUE(pseudocounts)) {
    vals <- stored_values(input$counts)
    if (any(vals != round(vals)))
      warning(
        '`counts` are not whole numbers, but the default `pseudocounts` ',
        '(0.1, 0.5, 1) are on the scale of counts: give `pseudocounts` on ',
        'the scale of your data.', call. = FALSE )
  }


  #--------------------------------------------------------
  # The base run; its warnings are given as usual.
  #--------------------------------------------------------
  seen <- character(0)
  base <- withCallingHandlers(
    da_sens_call(counts, id, args),
    warning = function (w) seen <<- c(seen, conditionMessage(w)) )

  # ANCOM-BC's own sensitivity analysis: fits already made by the base run.
  own <- if (v2) attr(base, 'sensitivity', exact = TRUE) else NULL


  #--------------------------------------------------------
  # The other runs; their warnings are collected and given once.
  #--------------------------------------------------------
  msgs <- character(0)
  fits <- vector('list', nrow(runs))
  fits[[1L]] <- base

  for (i in seq_len(nrow(runs))[-1L]) {

    r  <- runs[i, , drop = FALSE]
    pc <- r$pseudocount
    z  <- r$zeros
    s  <- r$scale_sd

    # A pseudocount the base run's sensitivity analysis already fitted.
    if (!is.null(own) && !is.na(pc) && is.na(z) && is.na(s)) {
      k <- match(pc, own$pseudocounts)
      if (!is.na(k)) { fits[[i]] <- own$fits[[k]]; next }
    }

    a <- args
    if (v2)        a[['sensitivity']] <- FALSE
    if (!is.na(pc)) a[['pseudocount']] <- pc
    if (!is.na(s))  a[['scale_sd']]    <- s
    cnt <- counts
    if (!is.na(z)) {
      cnt <- replace_zeros(input$counts, z, margin = input$margin)
      a[['pseudocount']] <- NULL
      a[['margin']]      <- input$margin
    }

    fits[[i]] <- withCallingHandlers(
      tryCatch(
        da_sens_call(cnt, id, a),
        error = function (e) stop(
          'The run with ', r$setting, ' failed: ', conditionMessage(e),
          call. = FALSE )),
      warning = function (w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart('muffleWarning')
      })
  }

  for (m in setdiff(unique(msgs), seen)) warning(m, call. = FALSE)


  #--------------------------------------------------------
  # The table: a row per feature [x covariate] x run.
  #--------------------------------------------------------
  eff <- method$effect
  out <- do.call(rbind, lapply(seq_along(fits), function (i) {
    f    <- fits[[i]]
    cov  <- if (is.null(f$covariate)) rep(NA_character_, nrow(f)) else as.character(f$covariate)
    keep <- is.na(cov) | cov != '(Intercept)'
    data.frame(
      feature    = as.character(f$feature)[keep],
      covariate  = cov[keep],
      run        = runs$run[[i]],
      parameter  = runs$parameter[[i]],
      setting    = runs$setting[[i]],
      effect     = if (is.na(eff)) rep(NA_real_, sum(keep)) else as.double(f[[eff]][keep]),
      p_adjusted = if (is.null(f$p_adjusted)) rep(NA_real_, sum(keep)) else as.double(f$p_adjusted[keep]),
      reject     = as.logical(f$reject[keep]),
      stringsAsFactors = FALSE )
  }))
  rownames(out) <- NULL

  log_base <- attr(base, 'log_base', exact = TRUE)
  if (is.null(log_base))
    log_base <- switch(paste0('', method$base), e = exp(1), '2' = 2, NA_real_)

  attr(out, 'method')      <- id
  attr(out, 'effect_name') <- eff
  attr(out, 'alpha')       <- attr(base, 'alpha', exact = TRUE)
  attr(out, 'log_base')    <- log_base
  attr(out, 'covariate')   <- if (method$design == 'formula') da_term_covariate(base) else NULL
  attr(out, 'runs')        <- runs
  attr(out, 'fits')        <- fits
  attr(out, 'summary')     <- da_sens_summary(out)
  class(out) <- c('ecodive_da_sensitivity', 'data.frame')

  out
}



# `...` with each name completed to the method's argument it partially
# matches (`pseudo = 1` becomes `pseudocount = 1`), so that a run's setting
# replaces it instead of being given twice.
da_sens_args <- function (args, method) {
  nms <- names(args)
  if (is.null(nms)) return (args)
  full <- unique(c('groups', 'data', 'formula', setdiff(method$params, '...')))
  for (i in which(nzchar(nms) & !nms %in% full)) {
    hit <- pmatch(nms[[i]], full)
    if (!is.na(hit)) nms[[i]] <- full[[hit]]
  }
  names(args) <- nms
  args
}


# One run: da_test() with these arguments. quote(counts), so that a large
# table isn't copied into the call.
da_sens_call <- function (counts, id, args) {
  do.call(da_test, c(list(quote(counts), id), args))
}


# The table as the methods read it (any input type as a matrix) and its
# margin, for replace_zeros() and the whole-number check.
da_sens_input <- function (counts, margin) {
  validate_counts()
  validate_margin()
  list(counts = counts, margin = margin)
}


# The grids, checked against the method: a list of `pseudocount`, `zeros`
# and `scale_sd`, each a vector (empty when not given).
da_sens_grids <- function (method, args, pseudocounts, zeros, scale_sds) {

  id   <- method$id
  none <- function (x) is.null(x) || isFALSE(x)

  if (none(pseudocounts) && none(zeros) && none(scale_sds))
    stop(
      'Give at least one of `pseudocounts`, `zeros` or `scale_sds`: the ',
      'settings to vary.', call. = FALSE )

  having <- function (arg) {
    ok <- vapply(list_da_tests('func'), function (f) arg %in% names(formals(f)), logical(1))
    paste(names(ok)[ok], collapse = ', ')
  }

  grid_num <- function (x, arg, default, min_ok) {
    if (none(x)) return (numeric(0))
    if (isTRUE(x)) return (default)
    if (!is.numeric(x) || length(x) < 1L || any(!is.finite(x)) || any(x < 0) ||
        (!min_ok && any(x == 0)))
      stop(
        '`', arg, '` must be TRUE or a vector of ',
        if (min_ok) 'non-negative' else 'positive', ' numbers.', call. = FALSE )
    if (anyDuplicated(x))
      stop('`', arg, '` must not repeat a value.', call. = FALSE)
    as.double(x)
  }

  pcs <- grid_num(pseudocounts, 'pseudocounts', c(0.1, 0.5, 1),       FALSE)
  sds <- grid_num(scale_sds,    'scale_sds',    c(0, 0.25, 0.5, 1),   TRUE)

  if (length(pcs)) {
    if (!'pseudocount' %in% method$params)
      stop(
        '`pseudocounts`: `', id, '` has no `pseudocount` argument. Methods ',
        'with one: ', having('pseudocount'), '.', call. = FALSE )

    # A norm without logs adds no pseudocount: every run would be the same.
    if ('norm' %in% method$params) {
      norm <- if ('norm' %in% names(args)) args$norm else eval(formals(method$func)$norm)
      norm <- tryCatch(
        match.arg(tolower(norm), c('none', 'percent', 'chord', 'hellinger', 'binary', 'clr', 'rclr', 'log')),
        error = function (e) NULL )
      if (!is.null(norm) && !norm %in% c('clr', 'log'))
        stop(
          '`pseudocounts`: with `norm = "', norm, '"`, `', id, '` adds no ',
          'pseudocount, so every run would be the same. Use `norm = "clr"` ',
          'or `"log"`.', call. = FALSE )
    }

    # linda() with imputed zeros adds no pseudocount either.
    if ('zero_handling' %in% method$params && 'zero_handling' %in% names(args)) {
      zh <- args$zero_handling
      if (is.character(zh) && length(zh) == 1L && !is.na(zh) &&
          startsWith('imputation', tolower(zh)) && nzchar(zh))
        stop(
          '`pseudocounts`: with `zero_handling = "imputation"`, `', id, '` adds ',
          'no pseudocount, so every run would be the same. Use ',
          '`zero_handling = "pseudocount"`.', call. = FALSE )
    }
  }

  if (length(sds) && !'scale_sd' %in% method$params)
    stop(
      '`scale_sds`: `', id, '` has no `scale_sd` argument. Methods with one: ',
      having('scale_sd'), '.', call. = FALSE )

  zs <- character(0)
  if (!none(zeros)) {
    if (!isTRUE(method$zero_free_ok)) {
      ok <- list_da_tests()
      stop(
        '`zeros`: `', id, '` models counts, so it can\'t take a table with ',
        'its zeros replaced. Methods that can: ',
        paste(ok$id[ok$zero_free_ok], collapse = ', '), '.', call. = FALSE )
    }
    if (isTRUE(zeros)) {
      zs <- list_zero_replacements('id')
    }
    else {
      if (!is.character(zeros) || length(zeros) < 1L || anyNA(zeros))
        stop(
          '`zeros` must be TRUE or the names of zero replacement methods ',
          '(list_zero_replacements()).', call. = FALSE )
      zs <- vapply(zeros, function (z) match_zero_replacement(z)$id, '', USE.NAMES = FALSE)
      if (anyDuplicated(zs))
        stop('`zeros` must not name a method twice.', call. = FALSE)
    }
  }

  list(pseudocount = pcs, zeros = zs, scale_sd = sds)
}


# The runs: a data frame of `run` (0 = base), `parameter`, the settings
# (`NA`: the base run's) and `setting`, the run in words.
da_sens_runs <- function (grids, crossed) {

  one <- function (pc = NA_real_, z = NA_character_, s = NA_real_)
    data.frame(pseudocount = pc, zeros = z, scale_sd = s, stringsAsFactors = FALSE)

  pcs <- grids$pseudocount
  zs  <- grids$zeros
  sds <- grids$scale_sd
  set <- list(one())

  if (!crossed) {
    if (length(pcs)) set <- c(set, list(one(pc = pcs)))
    if (length(zs))  set <- c(set, list(one(z  = zs)))
    if (length(sds)) set <- c(set, list(one(s  = sds)))
  }
  else {
    # The zero treatments (pseudocounts or replacements: alternatives),
    # each crossed with every scale_sd.
    treat <- list()
    if (length(pcs)) treat <- c(treat, list(one(pc = pcs)))
    if (length(zs))  treat <- c(treat, list(one(z  = zs)))
    treat <- if (length(treat)) do.call(rbind, treat) else one()
    if (!length(sds)) sds <- NA_real_
    grid  <- treat[rep(seq_len(nrow(treat)), each = length(sds)), , drop = FALSE]
    grid$scale_sd <- rep(sds, times = nrow(treat))
    set <- c(set, list(grid))
  }

  runs <- do.call(rbind, set)
  rownames(runs) <- NULL

  given <- cbind(!is.na(runs$pseudocount), !is.na(runs$zeros), !is.na(runs$scale_sd))
  n_set <- rowSums(given)
  param <- c('pseudocount', 'zeros', 'scale_sd')
  runs$parameter <- ifelse(
    n_set == 0L, 'base',
    ifelse(n_set == 1L, param[max.col(given, 'first')], 'crossed') )

  runs$setting <- vapply(seq_len(nrow(runs)), function (i) {
    s <- c(
      if (given[i, 1L]) paste('pseudocount =', as.character(runs$pseudocount[[i]])),
      if (given[i, 2L]) paste('zeros =',       runs$zeros[[i]]),
      if (given[i, 3L]) paste('scale_sd =',    as.character(runs$scale_sd[[i]])) )
    if (length(s)) paste(s, collapse = ', ') else 'base'
  }, '')

  runs$run <- seq_len(nrow(runs)) - 1L
  runs[, c('run', 'parameter', 'pseudocount', 'zeros', 'scale_sd', 'setting')]
}


# The per-feature [x covariate] summary of a da_sensitivity() table.
da_sens_summary <- function (x) {

  key <- paste(x$feature, x$covariate, sep = '\r')
  lev <- unique(key)
  grp <- factor(key, lev)
  first <- match(lev, key)

  rej    <- x$reject %in% TRUE
  n_runs <- tabulate(grp, length(lev))
  n_rej  <- tabulate(grp[rej], length(lev))
  passed <- n_rej == 0L | n_rej == n_runs

  base_rej <- rep(FALSE, length(lev))
  b <- x$run == 0L
  base_rej[match(key[b], lev)] <- rej[b]

  sgn <- split(x$effect, grp)
  sign_stable <- vapply(sgn, function (e) {
    e <- e[is.finite(e)]
    if (length(e)) length(unique(sign(e))) == 1L else NA
  }, NA, USE.NAMES = FALSE)

  data.frame(
    feature     = x$feature[first],
    covariate   = x$covariate[first],
    n_runs      = n_runs,
    n_reject    = n_rej,
    prop_reject = n_rej / n_runs,
    passed      = passed,
    robust      = passed & base_rej,
    sign_stable = sign_stable,
    stringsAsFactors = FALSE )
}


# A model method's covariate of interest: the first design column of the
# term `groups` names, else the first column that isn't the intercept.
da_term_covariate <- function (fit) {
  covs <- unique(as.character(fit$covariate))
  g    <- attr(fit, 'groups', exact = TRUE)
  des  <- attr(fit, 'design', exact = TRUE)
  asg  <- attr(des, 'assign')
  term <- attr(des, 'term.labels')
  if (is.null(term)) term <- attr(stats::terms(attr(fit, 'formula', exact = TRUE)), 'term.labels')
  if (is.character(g) && length(g) == 1L && !is.null(asg) && g %in% term) {
    hit <- colnames(des)[asg == match(g, term)]
    hit <- hit[hit %in% covs]
    if (length(hit)) return (hit[[1L]])
  }
  rest <- setdiff(covs, '(Intercept)')
  if (length(rest)) rest[[1L]] else covs[[1L]]
}



#' Plot a Differential Abundance Sensitivity Analysis
#'
#' Two displays of a [da_sensitivity()] result: which features each run
#' rejects (`"tiles"`), and how each feature's effect changes with one
#' setting (`"traces"`).
#'
#' @param x  A [da_sensitivity()] result.
#'
#' @param type  `"tiles"` (the default) or `"traces"`.
#'
#' @param covariate  For the model methods, which covariate to show: one of
#'        `x$covariate`. Default: the covariate of the term of interest
#'        (`attr(x, "covariate")`).
#'
#' @param parameter  Which runs to show. For `"tiles"`, the base run and the
#'        runs of one `parameter` (`"pseudocount"`, `"zeros"`, `"scale_sd"`
#'        or `"crossed"`); default: every run. For `"traces"`, the setting
#'        to draw the effects against (`"pseudocount"`, `"zeros"` or
#'        `"scale_sd"`); default: the first numeric one varied. With
#'        `crossed = TRUE`, the traces hold the other setting at its first
#'        value, as the subtitle says.
#'
#' @param n_features  The most features to show. Default: `30`.
#'
#' @param ...  Further graphical parameters for [graphics::image()]
#'        (`"tiles"`) or [graphics::plot()] (`"traces"`), such as `main`.
#'
#' @return `x`, invisibly.
#'
#' @details
#'
#' **Tiles.** A row per feature (those rejected in at least one run: the
#' most often rejected first, then by the size of their effect in the base
#' run), a column per run (the base run first, then the runs of each
#' setting, separated by vertical lines, with the setting named above).
#' A cell is red when the run rejects the feature with a positive effect,
#' blue with a negative one, dark grey when the method gives no effect
#' (ANCOM), light grey when the run doesn't reject it, and white when the
#' run has no result for it. A feature's name is bold when every run makes
#' the same decision (`summary$passed`).
#'
#' **Traces.** For one setting, a line per feature (those rejected in the
#' base run, else in any run, the largest effects first) through its effect
#' at each of the setting's values, with filled points where the run
#' rejects it and hollow ones where it doesn't, and a dotted line at 0.
#' Names are bold as in the tiles. Against `scale_sd`, this is the display of
#' ALDEx2's `plotGamma()`. ANCOM has no effects to draw.
#'
#' @seealso [da_sensitivity()]
#'
#' @export
#' @examples
#'     x <- rbind(
#'       c( 20, 110, 100, 101, 100, 103, 104),
#'       c( 33, 110, 120, 100, 101, 100, 102),
#'       c( 12, 110, 100, 110, 100,  50,  90),
#'       c(202, 201,   9,  10,  10,  11,  11),
#'       c(200, 202,  10,  10,  13,  10,  10),
#'       c(203, 201,  14,  10,  10,  13,  12) )
#'     colnames(x) <- paste0('b', 1:7)
#'     grp <- rep(c('treatment', 'placebo'), each = 3)
#'
#'     res <- da_sensitivity(
#'       x, 'da_ttest', grp, pseudocounts = c(0.01, 0.1, 0.5, 1, 5),
#'       zeros = TRUE, cpus = 1 )
#'     plot(res)
#'     plot(res, 'traces')
#'
plot.ecodive_da_sensitivity <- function (
    x, type = c('tiles', 'traces'), covariate = NULL, parameter = NULL,
    n_features = 30, ...) {

  type <- match.arg(type)
  n    <- plot_count(n_features, 'n_features')
  if (n < 1)
    stop('`n_features` must be at least 1.', call. = FALSE)

  runs <- attr(x, 'runs', exact = TRUE)
  summ <- attr(x, 'summary', exact = TRUE)
  if (is.null(runs) || is.null(summ))
    stop('`x` has lost its `runs` or `summary` attribute.', call. = FALSE)

  # The covariate (model methods).
  covs <- unique(x$covariate)
  if (all(is.na(covs))) {
    covariate <- NA_character_
  }
  else {
    if (is.null(covariate)) covariate <- attr(x, 'covariate', exact = TRUE)
    if (is.null(covariate)) covariate <- covs[[1L]]
    if (!is.character(covariate) || length(covariate) != 1 || !covariate %in% covs)
      stop('`covariate` must be one of ', name_list(covs), '.', call. = FALSE)
  }
  in_cov <- function (v) if (is.na(covariate)) is.na(v) else v %in% covariate
  d    <- data.frame(unclass(x)[c('feature', 'run', 'effect', 'reject')], stringsAsFactors = FALSE)
  d    <- d[in_cov(x$covariate), , drop = FALSE]
  summ <- summ[in_cov(summ$covariate), , drop = FALSE]

  info <- list(
    alpha       = attr(x, 'alpha',       exact = TRUE),
    effect_name = attr(x, 'effect_name', exact = TRUE),
    log_base    = attr(x, 'log_base',    exact = TRUE),
    method      = attr(x, 'method',      exact = TRUE) )
  if (is.null(info$effect_name)) info$effect_name <- NA_character_

  args <- list(...)
  main <- da_sens_plot_title(x)

  if (type == 'tiles') da_sens_plot_tiles(d, runs, summ, covariate, parameter, n, args, main, info)
  else                 da_sens_plot_traces(d, runs, summ, covariate, parameter, n, args, main, info)

  invisible(x)
}


# "Sensitivity of ANCOM-BC2": the method's display name from DA_TESTS.
da_sens_plot_title <- function (x) {
  paste('Sensitivity of', da_method_name(attr(x, 'method', exact = TRUE)))
}


# Each feature's effect in the base run, by the summary's rows.
da_sens_base_effect <- function (x, summ) {
  b <- x[x$run == 0L, , drop = FALSE]
  b$effect[match(summ$feature, b$feature)]
}


da_sens_plot_tiles <- function (x, runs, summ, covariate, parameter, n, args, main, info) {

  params <- setdiff(unique(runs$parameter), 'base')
  if (!is.null(parameter)) {
    if (!is.character(parameter) || length(parameter) != 1 || !parameter %in% params)
      stop('`parameter` must be one of ', name_list(params), '.', call. = FALSE)
    runs <- runs[runs$parameter %in% c('base', parameter), , drop = FALSE]
  }

  # Features rejected in any run: most often first, then by base effect.
  eff0 <- da_sens_base_effect(x, summ)
  ord  <- order(-summ$prop_reject, -abs(eff0), seq_len(nrow(summ)), na.last = TRUE)
  ord  <- ord[summ$n_reject[ord] > 0L]
  none <- !length(ord)
  if (none) ord <- order(-abs(eff0), seq_len(nrow(summ)), na.last = TRUE)
  keep <- utils::head(ord, n)
  feats <- summ$feature[keep]

  # Cell codes: 1 rejected, negative; 2 not rejected; 3 rejected, positive;
  # 4 rejected, no effect; NA no result.
  m <- matrix(NA_integer_, length(feats), nrow(runs))
  i <- x$feature %in% feats & x$run %in% runs$run
  e <- x$effect[i]
  code <- ifelse(
    !(x$reject[i] %in% TRUE), 2L,
    ifelse(is.finite(e) & e < 0, 1L, ifelse(is.finite(e) & e > 0, 3L, 4L)) )
  m[cbind(match(x$feature[i], feats), match(x$run[i], runs$run))] <- code

  pal  <- grDevices::hcl.colors(11, 'Blue-Red 3')
  cols <- c(pal[[1L]], 'grey90', pal[[11L]], 'grey35')

  # Column labels: the setting's value (the full setting for crossed runs).
  labs <- ifelse(
    runs$parameter == 'pseudocount', as.character(runs$pseudocount),
    ifelse(runs$parameter == 'zeros', runs$zeros,
    ifelse(runs$parameter == 'scale_sd', as.character(runs$scale_sd), runs$setting)) )

  # Room for the names on the left and the labels below, as dotchart() makes it.
  opar <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  csi   <- graphics::par('csi')
  linch <- max(graphics::strwidth(feats, 'inch', cex = 0.8, font = 2), 0, na.rm = TRUE)
  binch <- max(graphics::strwidth(labs,  'inch', cex = 0.8), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[2L] <- max(mai[2L], linch + 0.3)
  mai[1L] <- max(mai[1L], binch + 2.4 * csi)
  mai[3L] <- max(mai[3L], 3.2 * csi)
  graphics::par(mai = mai)

  rows <- rev(seq_along(feats))  # image() draws its first row at the bottom
  plot_do(
    graphics::image,
    c(list(seq_len(nrow(runs)), seq_along(feats), t(m[rows, , drop = FALSE])), args),
    list(
      col = cols, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
      axes = FALSE, xlab = '', ylab = '', main = main ))

  # Thin white lines between the cells.
  if (nrow(runs) > 1L)  graphics::abline(v = seq_len(nrow(runs) - 1L) + 0.5, col = 'white')
  if (length(feats) > 1L) graphics::abline(h = seq_len(length(feats) - 1L) + 0.5, col = 'white')

  graphics::axis(1, at = seq_len(nrow(runs)), labels = labs, tick = FALSE, las = 2, cex.axis = 0.8)
  bold <- summ$passed[keep][rows] %in% TRUE
  at   <- seq_along(feats)
  if (any(bold))
    graphics::axis(2, at = at[bold], labels = feats[rows][bold], tick = FALSE, las = 1, cex.axis = 0.8, font = 2)
  if (any(!bold))
    graphics::axis(2, at = at[!bold], labels = feats[rows][!bold], tick = FALSE, las = 1, cex.axis = 0.8)
  graphics::box()

  # The settings' groups, separated and named above.
  grp   <- runs$parameter
  edges <- which(grp[-1L] != grp[-length(grp)])
  if (length(edges)) graphics::abline(v = edges + 0.5, lwd = 1.5)
  starts <- c(1L, edges + 1L)
  ends   <- c(edges, length(grp))
  graphics::mtext(grp[starts], side = 3, line = 0.2, at = (starts + ends) / 2, cex = 0.8)

  n_feat <- nrow(summ)
  n_pass <- sum(summ$passed & summ$n_reject > 0L)
  n_any  <- sum(summ$n_reject > 0L)
  alpha  <- info$alpha
  key <- sprintf(
    '%s; bold: the same decision in every run%s%s%s',
    if (is.na(info$effect_name)) 'Dark: rejected'
    else 'Red / blue: rejected, positive / negative effect',
    if (none) '' else sprintf(' (%d of the %d rejected in any run)', n_pass, n_any),
    if (!is.na(covariate)) sprintf('; covariate %s', covariate) else '',
    if (none) '; no feature rejected in any run'
    else if (length(keep) < n_any) sprintf('; top %d of %d shown', length(keep), n_any)
    else '' )
  if (!is.null(alpha) && !is.na(info$effect_name))
    key <- sprintf('%s; alpha = %s', key, plot_num(alpha))
  if (!'sub' %in% names(args))
    graphics::mtext(
      key, side = 1, cex = 0.8,
      line = binch / csi + 1.2 )

  invisible(n_feat)
}


da_sens_plot_traces <- function (x, runs, summ, covariate, parameter, n, args, main, info) {

  eff_name <- info$effect_name
  if (is.null(eff_name) || is.na(eff_name))
    stop(
      'The "traces" plot draws effect sizes, and `',
      info$method, '` gives none: use type = "tiles".',
      call. = FALSE )

  cols   <- c('pseudocount', 'zeros', 'scale_sd')
  varied <- cols[vapply(cols, function (p) any(!is.na(runs[[p]])), logical(1))]
  if (is.null(parameter)) {
    num <- intersect(c('pseudocount', 'scale_sd'), varied)
    parameter <- if (length(num)) num[[1L]] else varied[[1L]]
  }
  if (!is.character(parameter) || length(parameter) != 1 || !parameter %in% varied)
    stop('`parameter` must be one of ', name_list(varied), '.', call. = FALSE)

  # The runs of this setting; with crossed runs, the other setting held at
  # its first value.
  sel <- runs$run != 0L & !is.na(runs[[parameter]])
  held <- NULL
  if (parameter == 'scale_sd') {
    tr  <- paste(runs$pseudocount, runs$zeros)
    sel <- sel & tr == tr[sel][[1L]]
    if (any(!is.na(c(runs$pseudocount[sel], runs$zeros[sel]))))
      held <- sub(', scale_sd = .*$', '', runs$setting[sel][[1L]])
  }
  else {
    s   <- runs$scale_sd
    sel <- sel & s %in% s[sel][[1L]]
    if (!is.na(s[sel][[1L]])) held <- paste('scale_sd =', s[sel][[1L]])
  }
  r  <- runs[sel, , drop = FALSE]
  xv <- if (parameter == 'zeros') seq_len(nrow(r)) else r[[parameter]]
  o  <- order(xv)
  r  <- r[o, , drop = FALSE]
  xv <- xv[o]

  # Features rejected in the base run, else in any run, else all; the
  # largest base effects first.
  eff0 <- da_sens_base_effect(x, summ)
  base_rej <- vapply(summ$feature, function (f) any(x$reject[x$run == 0L & x$feature == f] %in% TRUE), NA)
  pick <- if (any(base_rej)) base_rej else if (any(summ$n_reject > 0L)) summ$n_reject > 0L else rep(TRUE, nrow(summ))
  ord  <- order(-abs(eff0), seq_len(nrow(summ)), na.last = TRUE)
  ord  <- ord[pick[ord]]
  keep <- utils::head(ord, n)
  feats <- summ$feature[keep]

  y   <- matrix(NA_real_, length(feats), nrow(r))
  rej <- matrix(FALSE,    length(feats), nrow(r))
  i   <- x$feature %in% feats & x$run %in% r$run
  ij  <- cbind(match(x$feature[i], feats), match(x$run[i], r$run))
  y[ij]   <- x$effect[i]
  rej[ij] <- x$reject[i] %in% TRUE

  lb    <- info$log_base
  units <- if (is.null(lb) || is.na(lb)) '' else if (isTRUE(all.equal(lb, exp(1)))) ' (natural log)'
           else sprintf(' (log%s)', plot_num(lb))
  ylab  <- sprintf('%s%s%s', eff_name, units, if (is.na(covariate)) '' else paste(',', covariate))
  xlab  <- switch(parameter, pseudocount = 'Pseudocount', zeros = 'Zero replacement', scale_sd = 'scale_sd (gamma)')
  sub   <- sprintf(
    'Filled: rejected; bold: the same decision in every run%s%s',
    if (is.null(held)) '' else sprintf('; with %s', held),
    if (length(keep) < sum(pick)) sprintf('; top %d of %d shown', length(keep), sum(pick)) else '' )

  colr <- grDevices::hcl.colors(max(length(feats), 1L), 'Dark 3')
  yr   <- range(c(0, y[is.finite(y)]))

  # Room for the names on the right.
  opar  <- graphics::par('mai', 'mar')
  on.exit(graphics::par(opar))
  rinch <- max(graphics::strwidth(feats, 'inch', cex = 0.7, font = 2), 0, na.rm = TRUE)
  mai   <- graphics::par('mai')
  mai[4L] <- max(mai[4L], rinch + 0.2)
  graphics::par(mai = mai)

  plot_do(
    graphics::plot, c(list(range(xv), yr), args),
    list(
      type = 'n', xlab = xlab, ylab = ylab, main = main, sub = sub,
      log  = if (parameter == 'pseudocount') 'x' else '',
      xaxt = if (parameter == 'zeros') 'n' else 's' ))
  if (parameter == 'zeros') graphics::axis(1, at = xv, labels = r$zeros)
  graphics::abline(h = 0, lty = 3)

  bold <- summ$passed[keep] %in% TRUE
  for (k in seq_along(feats)) {
    graphics::lines(xv, y[k, ], col = colr[[k]])
    graphics::points(xv, y[k, ], col = colr[[k]], pch = ifelse(rej[k, ], 19, 1))
  }
  last <- y[, ncol(y)]
  ok   <- is.finite(last)
  if (any(ok))
    graphics::text(
      if (graphics::par('xlog')) 10^graphics::par('usr')[[2L]] else graphics::par('usr')[[2L]],
      last[ok], feats[ok], pos = 4, xpd = TRUE,
      cex = 0.7, col = colr[ok], font = ifelse(bold[ok], 2, 1) )

  invisible()
}
