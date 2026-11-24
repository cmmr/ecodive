# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Batch correction: what batch_mmuphin() and batch_combat() (and, later,
# batch_combat_seq() and batch_conqur()) share.
#
#   BATCH_CORRECTIONS   Facts about each method, keyed by `method` (the id):
#                       its display name (plot titles), and what the "means"
#                       plot shows. The batch_correct() wrapper, with
#                       list_batch_corrections() and match_batch_correction(),
#                       is planned to read the same table (with an alt_ids
#                       column for its aliases).
#   batch_input()       Checks `counts`, `margin` and `cpus`.
#   batch_design()      `batch`, `formula` and `data`: the batch codes and
#                       the design (batch indicators, then covariates).
#   batch_result()      The corrected table, with its attributes.
#   print() and plot() methods for class `ecodive_batch`.
BATCH_CORRECTIONS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                       id,       means
      MMUPHin Batch Correction,   mmuphin,  relative abundance
      ComBat,                     combat,   value
  " )

  return (df)
})


# The display name of a batch correction method, or NULL.
batch_method_name <- function (method) {
  if (!is.character(method) || length(method) != 1) return (NULL)
  i <- match(method, BATCH_CORRECTIONS$id)
  if (is.na(i)) NULL else BATCH_CORRECTIONS$name[[i]]
}



#------------------------------------------------------------------------------
# Checks the caller's `counts`, `margin` and `cpus` (as validate_*() do, in
# the caller's environment). `nonneg`: refuse negative values. Returns
# list(counts, margin, n, D, labels (samples' names or NULL), features
# (names or NULL), dn (counts' dimnames), vals (the stored values)).
#------------------------------------------------------------------------------
batch_input <- function (env = parent.frame(), nonneg = TRUE) {

  validate_counts(env)
  validate_margin(env)
  validate_cpus(env)

  counts <- env$counts
  if (is.logical(counts) && is.matrix(counts)) {
    storage.mode(counts) <- 'double'
    env$counts <- counts
  }

  vals <- stored_values(counts)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  if (!all(is.finite(vals)))
    stop('`counts` must not contain NA, NaN, or infinite values.', call. = FALSE)
  if (nonneg && any(vals < 0))
    stop('`counts` must not contain negative values.', call. = FALSE)

  margin <- env$margin
  pkg    <- get_matrix_package(counts)
  dims   <- switch(pkg, base = dim(counts), slam = c(counts$nrow, counts$ncol), Matrix = counts@Dim)
  dn     <- switch(pkg, base = dimnames(counts), slam = counts$dimnames, Matrix = counts@Dimnames)
  if (is.null(dn)) dn <- list(NULL, NULL)

  list(
    counts   = counts,
    margin   = margin,
    n        = as.integer(dims[[margin]]),
    D        = as.integer(dims[[3L - margin]]),
    labels   = if (length(dn[[margin]])) as.character(dn[[margin]]),
    features = if (length(dn[[3L - margin]])) as.character(dn[[3L - margin]]),
    dn       = dn,
    vals     = vals,
    sparse   = pkg != 'base' && length(vals) < prod(as.double(dims)) )
}



#------------------------------------------------------------------------------
# `batch`, `formula` and `data`, by ecodive's per-sample rules.
#
#   batch     One value per sample, or the name of a column of `data`.
#   formula   NULL, or a one-sided formula over `data` naming the covariates
#             whose effects are kept.
#   data      NULL, or a data.frame of sample metadata (rows matched to the
#             samples by name, else by position).
#   info      batch_input()'s result.
#   min_size  The fewest samples a batch may have.
#
# Returns list(batch (a factor, one value per sample), levels, codes
# (0-based), sizes, mod (n x q covariates, no intercept), X (n x (B + q):
# the batch indicators, then mod)).
#------------------------------------------------------------------------------
batch_design <- function (batch, formula, data, info, min_size = 2L) {

  n      <- info$n
  labels <- info$labels

  if (is.null(batch))
    stop('`batch` is required: the batch of each sample, or the name of a column of `data`.', call. = FALSE)

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame of sample metadata.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    rn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, rn, nrow(data), ref = '`counts`')
  }


  #--------------------------------------------------------
  # The batches: a factor's levels in order, other values
  # sorted.
  #--------------------------------------------------------
  b <- per_sample(batch, 'batch', data, idx, labels, n, ref = '`counts`')
  if (anyNA(b))
    stop('`batch` must not contain missing values (NA).', call. = FALSE)
  fac <- if (is.factor(b)) droplevels(b)
         else factor(as.character(b), levels = as.character(sort(unique(b), method = 'radix')))

  lev <- levels(fac)
  B   <- length(lev)
  if (B < 2)
    stop('`batch` must have at least two batches; it has one: ', name_list(lev), '.', call. = FALSE)

  sizes <- tabulate(fac, B)
  small <- lev[sizes < min_size]
  if (length(small))
    stop(
      'Each batch needs at least ', min_size, ' samples; ',
      if (length(small) == 1) 'this one has' else 'these have', ' fewer: ',
      name_list(small), '.', call. = FALSE )

  codes <- as.integer(fac) - 1L
  bmod  <- matrix(0, n, B, dimnames = list(NULL, lev))
  bmod[cbind(seq_len(n), codes + 1L)] <- 1


  #--------------------------------------------------------
  # The covariates, and confounding.
  #--------------------------------------------------------
  mod <- batch_covariates(formula, data, idx, n)

  if (ncol(mod) && qr(mod)$rank < ncol(mod))
    stop(
      'The covariates in `formula` are confounded: some of their columns are ',
      'combinations of the others. Remove one or more of them.', call. = FALSE )

  if (ncol(mod) && qr(cbind(bmod, mod))$rank < B + ncol(mod))
    stop(
      'Covariates and batch are confounded: a covariate in `formula` is a ',
      'combination of the batch indicators (for example, constant within each ',
      'batch), so its effect can\'t be told apart from the batches\'. Remove it.',
      call. = FALSE )

  list(
    batch  = fac,
    levels = lev,
    codes  = codes,
    sizes  = sizes,
    mod    = mod,
    X      = cbind(bmod, mod) )
}


# The covariates of `formula`: an n x q model matrix without the intercept
# (treatment contrasts; characters and logicals as factors with sorted
# levels).
batch_covariates <- function (formula, data, idx, n) {

  if (is.null(formula)) return (matrix(0, n, 0))

  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    txt <- if (grepl('~', formula, fixed = TRUE)) formula else paste('~', formula)
    f   <- tryCatch(str2lang(txt), error = function (e) NULL)
    if (!is.call(f) || !identical(f[[1]], as.name('~')))
      stop('`formula` must be a one-sided formula, such as `~ disease + age`.', call. = FALSE)
    formula <- eval(f)
    environment(formula) <- globalenv()
  }
  if (!inherits(formula, 'formula'))
    stop('`formula` must be a one-sided formula, such as `~ disease + age`.', call. = FALSE)
  if (length(formula) == 3L) {
    env     <- environment(formula)
    formula <- formula[-2L]
    environment(formula) <- if (is.null(env)) globalenv() else env
  }

  if (is.null(data))
    stop(
      '`formula` needs `data`: a data.frame of sample metadata with the ',
      'variables it names.', call. = FALSE )
  data <- data[idx, , drop = FALSE]

  tt <- tryCatch(
    stats::terms(formula, data = data),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )

  vars    <- all.vars(tt)
  missing <- vars[!vars %in% names(data)]
  if (length(missing))
    stop('Not found in `data`: ', name_list(missing), call. = FALSE)

  na <- vars[vapply(data[vars], anyNA, logical(1))]
  if (length(na))
    stop('Missing values (NA) in `data` column', if (length(na) > 1) 's', ': ',
         name_list(na), call. = FALSE)

  inf <- vars[vapply(data[vars], function (v) is.numeric(v) && !all(is.finite(v)), logical(1))]
  if (length(inf))
    stop('Infinite values in `data` column', if (length(inf) > 1) 's', ': ',
         name_list(inf), call. = FALSE)

  if (length(attr(tt, 'term.labels')) == 0L) return (matrix(0, n, 0))
  attr(tt, 'intercept') <- 1L

  for (v in vars) {
    x <- data[[v]]
    if (is.factor(x))         data[[v]] <- droplevels(x)
    else if (is.logical(x))   data[[v]] <- factor(x, levels = c(FALSE, TRUE)[c(FALSE, TRUE) %in% x])
    else if (is.character(x)) data[[v]] <- factor(x, levels = sort(unique(x), method = 'radix'))
  }

  mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
  is_fac <- vapply(mf, is.factor, logical(1))
  one    <- names(mf)[is_fac][vapply(mf[is_fac], nlevels, integer(1)) < 2L]
  if (length(one))
    stop(
      'A categorical variable in `formula` has a single value, so its effect ',
      'can\'t be told apart from the batches\': ', name_list(one), call. = FALSE )

  X <- tryCatch(
    stats::model.matrix(tt, mf),
    error = function (e) stop('Invalid `formula`: ', conditionMessage(e), call. = FALSE) )
  X <- X[, colnames(X) != '(Intercept)', drop = FALSE]
  attr(X, 'assign') <- attr(X, 'contrasts') <- NULL
  rownames(X) <- NULL
  storage.mode(X) <- 'double'
  X
}



#------------------------------------------------------------------------------
# The corrected table, from the C code's list `res`: dense, in the
# orientation of `counts` with its dimnames, of class `ecodive_batch`, with
# its attributes (see ?batch_mmuphin, "Value").
#
#   res     adjusted (n x D), gamma_hat, gamma_star, delta_hat, delta_star
#           (D x B), means_before, means_after (D x B), overall (D).
#   info    batch_input()'s result.
#   design  batch_design()'s result.
#   method  The method's id.
#   adjusted  D logicals: the features adjusted.
#   attrs   Further attributes (named list), after batch_effects and
#           adjusted, and before the means.
#------------------------------------------------------------------------------
batch_result <- function (res, info, design, method, adjusted, attrs = list()) {

  out <- res$adjusted
  if (info$margin == 2L) out <- t(out)
  dn <- info$dn
  dimnames(out) <- if (all(vapply(dn, is.null, logical(1)))) NULL else dn

  D     <- info$D
  lev   <- design$levels
  B     <- length(lev)
  feats <- if (is.null(info$features)) seq_len(D) else info$features

  effects <- data.frame(
    feature    = rep(feats, B),
    batch      = factor(rep(lev, each = D), levels = lev),
    gamma_hat  = as.vector(res$gamma_hat),
    gamma_star = as.vector(res$gamma_star),
    delta_hat  = as.vector(res$delta_hat),
    delta_star = as.vector(res$delta_star),
    stringsAsFactors = FALSE )
  for (k in 3:6) effects[[k]][is.nan(effects[[k]])] <- NA_real_

  names(adjusted) <- info$features

  dn2 <- list(info$features, lev)
  m0  <- matrix(res$means_before, D, B, dimnames = dn2)
  m1  <- matrix(res$means_after,  D, B, dimnames = dn2)
  ov  <- res$overall
  names(ov) <- info$features

  attr(out, 'method')        <- method
  attr(out, 'batch_effects') <- effects
  attr(out, 'adjusted')      <- adjusted
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  attr(out, 'batch_means')          <- m0
  attr(out, 'batch_means_adjusted') <- m1
  attr(out, 'overall_means')        <- ov
  class(out) <- c('ecodive_batch', 'matrix', 'array')
  out
}


# `tol` and `max_iters`, checked.
batch_iter_args <- function (tol, max_iters) {
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) || !(tol > 0) || !is.finite(tol))
    stop('`tol` must be a single positive number.', call. = FALSE)
  if (!is.numeric(max_iters) || length(max_iters) != 1 || is.na(max_iters) ||
      max_iters < 1 || max_iters %% 1 != 0 || max_iters > .Machine$integer.max)
    stop('`max_iters` must be a single positive integer.', call. = FALSE)
  list(tol = as.double(tol), max_iters = as.integer(max_iters))
}


# TRUE or FALSE, for `arg`.
batch_flag <- function (v, arg) {
  if (!is.logical(v) || length(v) != 1 || is.na(v))
    stop('`', arg, '` must be TRUE or FALSE.', call. = FALSE)
  v
}



#------------------------------------------------------------------------------
# Printing and plotting.
#------------------------------------------------------------------------------

#' @rdname print_ecodive
#' @export
print.ecodive_batch <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  att <- attributes(x)
  sym <- print_sym()
  adj <- att$adjusted

  # Samples are rows unless the features' names (or number) say otherwise.
  rows_are_features <- if (!is.null(names(adj)) && !is.null(rownames(x)) &&
                           identical(rownames(x), names(adj)) &&
                           !identical(colnames(x), names(adj))) TRUE
                       else if (is.null(names(adj))) length(adj) == nrow(x) && length(adj) != ncol(x)
                       else FALSE
  rows <- if (rows_are_features) 'features' else 'samples'
  cols <- if (rows_are_features) 'samples'  else 'features'
  nm   <- batch_method_name(att$method)

  head <- c(
    sprintf('A matrix: %d %s %s %d %s', nrow(x), rows, sym$times, ncol(x), cols),
    if (!is.null(nm)) sprintf('Batch-corrected (%s)', nm),
    if (!is.null(adj)) sprintf('%d of %d features adjusted', sum(adj), length(adj)) )

  m <- unclass(x)
  attributes(m) <- att[c('dim', 'dimnames')]
  df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
  if (is.null(colnames(m))) names(df) <- sprintf('[,%d]', seq_len(ncol(m)))

  cat(
    print_table_lines(
      df, head, n = n, width = width, digits = digits, types = FALSE,
      rows = rows, attrs = print_attrs(x, c('dim', 'dimnames', 'class')) ),
    sep = '\n' )
  invisible(x)
}


#' Plot a batch correction
#'
#' MMUPHin's diagnostic figure for a corrected table, in base graphics:
#' how far the empirical Bayes step shrank each batch effect, and each
#' feature's batch means before and after the correction.
#'
#' @param x  The result of [batch_mmuphin()] or [batch_combat()].
#'
#' @param type  `"shrinkage"` (the default): each feature's location
#'        effect in each batch, \eqn{\hat\gamma} (x), against its shrunken
#'        value \eqn{\gamma^*} (y), coloured by batch, with the 1:1 line.
#'        Points below the line in the upper right (and above it in the
#'        lower left) were pulled towards the batch's average effect.
#'        `"means"`: for each feature (one vertical line of points per
#'        feature), its mean in each batch (y) against its mean over all
#'        samples (x), before the correction (black) and after (red, just to
#'        the right), with the 1:1 line. A batch effect spreads a feature's
#'        batch means away from the line; the correction draws them
#'        together. [batch_mmuphin()]'s means are of relative abundances;
#'        [batch_combat()]'s are of the values it corrected.
#'
#' @param legend  Where to put the legend: a position keyword for
#'        [graphics::legend()] (such as `"topleft"`), or `NULL` for none.
#'        Default: `"topleft"`.
#'
#' @param ...  Passed on to [graphics::plot()], such as `main`, `xlab`,
#'        `ylab`, `log`, `pch`, `cex`, `col` (`"shrinkage"`: one colour per
#'        batch; `"means"`: before and after).
#'
#' @return `x`, invisibly.
#'
#' @details
#' The plot reads only the attributes stored with the result, so a table
#' taken apart with `[` can't be plotted. The real check of a correction is
#' an ordination coloured by batch, and the PERMANOVA \eqn{R^2} of batch,
#' before and after (see the examples of [batch_mmuphin()]).
#'
#' @seealso [batch_mmuphin()], [batch_combat()]
#'
#' @export
#' @examples
#'     set.seed(1)
#'     counts <- matrix(rpois(20 * 30, 20), 20, 30)
#'     batch  <- rep(c('run1', 'run2'), each = 10)
#'     counts[11:20, 1:10] <- counts[11:20, 1:10] * 3
#'
#'     x <- batch_mmuphin(counts, batch)
#'     plot(x)
#'     plot(x, type = 'means')
plot.ecodive_batch <- function (x, type = c('shrinkage', 'means'), legend = 'topleft', ...) {

  type <- test_plot_type(type, c('shrinkage', 'means'))
  pos  <- plot_legend_pos(legend)
  args <- list(...)

  eff <- attr(x, 'batch_effects', exact = TRUE)
  m0  <- attr(x, 'batch_means', exact = TRUE)
  m1  <- attr(x, 'batch_means_adjusted', exact = TRUE)
  ov  <- attr(x, 'overall_means', exact = TRUE)
  if (is.null(eff) || is.null(m0) || is.null(m1) || is.null(ov))
    stop(
      '`x` has lost the attributes of its correction (parts taken with `[` ',
      'drop them): plot the whole result.', call. = FALSE )

  method <- attr(x, 'method', exact = TRUE)
  main   <- batch_method_name(method)
  if (is.null(main)) main <- 'Batch Correction'

  if (type == 'shrinkage') {

    ok  <- !is.na(eff$gamma_hat) & !is.na(eff$gamma_star)
    lev <- levels(eff$batch)
    col <- plot_arg(args, 'col', plot_palette(length(lev)))
    col <- rep_len(col, length(lev))
    sub <- if (isTRUE(attr(x, 'mean_only', exact = TRUE))) 'Locations only (mean_only)' else NULL
    rng <- range(c(eff$gamma_hat[ok], eff$gamma_star[ok]), finite = TRUE)
    if (!all(is.finite(rng))) rng <- c(-1, 1)

    plot_do(
      graphics::plot, c(list(eff$gamma_hat[ok], eff$gamma_star[ok]), plot_drop(args, 'col')),
      list(
        main = main, sub = sub, type = 'p', las = 1, pch = 16,
        col  = col[as.integer(eff$batch[ok])],
        xlim = rng, ylim = rng,
        xlab = 'Batch location effect (gamma hat)',
        ylab = 'Shrunken (gamma star)' ))
    graphics::abline(0, 1, col = PLOT_MUTED, lty = 2)
    plot_legend(pos, lev, col, pch = 16, title = 'Batch')
  }

  else {

    what <- BATCH_CORRECTIONS$means[match(method, BATCH_CORRECTIONS$id)]
    if (is.na(what)) what <- 'value'
    col  <- plot_arg(args, 'col', c('black', PLOT_HIGHLIGHT))
    col  <- rep_len(col, 2L)
    D    <- length(ov)
    B    <- ncol(m0)
    logx <- grepl('x', plot_arg(args, 'log', ''), fixed = TRUE)
    logy <- grepl('y', plot_arg(args, 'log', ''), fixed = TRUE)

    # MMUPHin draws the adjusted means 1% of the largest overall mean to the
    # right, so that the two sets don't overlap.
    x0 <- ov
    x1 <- if (logx) ov * 1.02 else ov + 0.01 * max(abs(ov[is.finite(ov)]), 0)
    xx <- c(rep(x0, B), rep(x1, B))
    yy <- c(as.vector(m0), as.vector(m1))
    ok <- is.finite(xx) & is.finite(yy)
    if (logx) ok <- ok & xx > 0
    if (logy) ok <- ok & yy > 0

    plot_do(
      graphics::plot, c(list(xx[ok], yy[ok]), plot_drop(args, 'col')),
      list(
        main = main, type = 'n', las = 1,
        xlab = paste('Overall mean', what),
        ylab = paste('Batch mean', what) ))
    graphics::abline(0, 1, col = PLOT_MUTED, lty = 2)

    for (f in seq_len(D)) {
      for (k in 1:2) {
        xv <- if (k == 1) x0[[f]] else x1[[f]]
        yv <- if (k == 1) m0[f, ] else m1[f, ]
        keep <- is.finite(yv) & is.finite(xv)
        if (logx) keep <- keep & xv > 0
        if (logy) keep <- keep & yv > 0
        if (sum(keep) > 1) graphics::lines(c(xv, xv), range(yv[keep]), col = col[[k]])
      }
    }
    graphics::points(
      xx[ok], yy[ok], col = rep(rep(col, each = D * B))[ok],
      pch = plot_arg(args, 'pch', 16), cex = plot_arg(args, 'cex', 0.7) )
    plot_legend(pos, c('Original', 'Adjusted'), col, pch = 16)
  }

  invisible(x)
}
