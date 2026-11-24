# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Rank plots (`type = "rank"` of the DA and ordination plot methods) and the
# rankings behind them, which log_ratio_features() also reads. Qurro's
# display rules (Fedarko et al. 2020; Qurro 0.9.0, BSD-3), not its code:
# features sorted ascending by their differential or loading, one bar each,
# with the numerator and denominator of a log ratio highlighted.



# The ranking that a result defines: one effect per feature. Returns
# list(effect, feature, reject (or NULL), what (the comparison or axis, for
# a subtitle), kind ('da' or 'ord'), axis (an ordination's axis name)).
#
#   * DA results: the effect column that DA_TESTS names for the result's
#     `method`, if it is `signed`, for one covariate of a model result (`covariate`, by
#     da_rank_covariate()'s default); `reject` from the result.
#   * Ordinations: column `axes` of the `loadings` attribute.
#   * A numeric vector: itself.
#
# `effect` keeps every feature, non-finite ones included; `feature` is the
# features' names, or their indices when unnamed.
rank_source <- function (x, covariate = NULL, axes = 1L, has_cov = !is.null(covariate), has_axes = FALSE) {

  # Ordinations: a column of the loadings.
  if (inherits(x, 'ecodive_ord')) {
    if (has_cov)
      stop('`covariate` is for differential abundance results; give `axes` for an ordination.', call. = FALSE)
    load <- attr(x, 'loadings', exact = TRUE)
    if (is.null(load) || !is.matrix(load) || nrow(load) == 0L)
      stop(
        'This ordination has no feature loadings to rank: ranks need an ',
        'ordination of a count table with a `loadings` attribute, such as ',
        'pca(), rpca(), rpca_phylo() or rda().', call. = FALSE )
    axis <- rank_axis(axes, ncol(load))
    nm   <- colnames(load)[axis]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) nm <- colnames(x)[axis]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) nm <- paste('Axis', axis)
    feat <- rownames(load)
    if (is.null(feat)) feat <- seq_len(nrow(load))
    return (list(
      effect = unname(as.double(load[, axis])), feature = feat, reject = NULL,
      what = nm, kind = 'ord', axis = nm ))
  }

  if (has_axes)
    stop('`axes` is for ordinations; for a differential abundance result, give `covariate`.', call. = FALSE)

  # DA results: the registry's effect column.
  if (is.data.frame(x)) {
    method <- attr(x, 'method', exact = TRUE)
    i      <- match(paste0('', method), DA_TESTS$id)
    if (is.na(i))
      stop(
        '`x` must be the result of a differential abundance method (see ',
        'list_da_tests()), an ordination with loadings, or a named numeric ',
        'vector.', call. = FALSE )
    col <- DA_TESTS$effect[[i]]
    if (is.na(col) || !isTRUE(DA_TESTS$signed[[i]]) || !col %in% names(x))
      stop(
        method, '() results have no signed effect, so they rank no features ',
        'into a reference frame.', call. = FALSE )
    r   <- x
    cov <- NULL
    if ('covariate' %in% names(x)) {
      cov <- da_rank_covariate(x, covariate)
      r   <- x[x$covariate == cov, , drop = FALSE]
    }
    else if (has_cov)
      stop(
        '`covariate` is for model results (with a `covariate` column); ',
        method, '() results have none.', call. = FALSE )
    return (list(
      effect = as.double(r[[col]]), feature = as.character(r$feature),
      reject = if ('reject' %in% names(r)) r$reject %in% TRUE,
      what = if (is.null(cov)) da_rank_contrast(x) else cov, kind = 'da',
      axis = NULL ))
  }

  # A vector of effects.
  if (is.numeric(x) && is.null(dim(x)) && !is.object(x)) {
    if (has_cov)
      stop('`covariate` is for differential abundance results.', call. = FALSE)
    feat <- names(x)
    if (is.null(feat)) feat <- seq_along(x)
    return (list(
      effect = unname(as.double(x)), feature = feat, reject = NULL,
      what = NULL, kind = 'vector', axis = NULL ))
  }

  stop(
    '`x` must be the result of a differential abundance method (see ',
    'list_da_tests()), an ordination with loadings, or a named numeric ',
    'vector.', call. = FALSE )
}


# A DA result's comparison for a subtitle ("a (n = 6) vs b (n = 6)"), when
# it has one, else NULL.
da_rank_contrast <- function (x) {
  ct <- attr(x, 'contrast', exact = TRUE)
  n  <- attr(x, 'n', exact = TRUE)
  if (length(ct) == 2L && length(n) == 2L) da_plot_contrast(x) else NULL
}


# Checks `axes`, one axis number from 1 to `n`, for a ranking.
rank_axis <- function (axes, n) {
  ok <- is.numeric(axes) && length(axes) >= 1L && is.finite(axes[[1]]) &&
    axes[[1]] == floor(axes[[1]]) && axes[[1]] >= 1 && axes[[1]] <= n
  if (!ok)
    stop(sprintf('`axes` must be an axis number from 1 to %d.', n), call. = FALSE)
  as.integer(axes[[1]])
}


# The covariate of a model DA result (with a `covariate` column) to plot or
# rank: `covariate` if given, checked; else ANCOM-BC's first covariate that
# isn't the intercept, or for the others the first column of the term of
# interest (da_term_covariate()).
da_rank_covariate <- function (x, covariate) {
  covs <- unique(as.character(x$covariate))
  if (is.null(covariate)) {
    if (isTRUE(attr(x, 'method', exact = TRUE) %in% c('ancombc', 'ancombc2'))) {
      covariate <- setdiff(covs, '(Intercept)')
      covariate <- if (length(covariate)) covariate[[1]] else covs[[1]]
    }
    else covariate <- da_term_covariate(x)
  }
  if (!is.character(covariate) || length(covariate) != 1 || !covariate %in% covs)
    stop('`covariate` must be one of ', name_list(covs), '.', call. = FALSE)
  covariate
}



# The layout of a rank plot: the features with a finite effect, sorted
# ascending (ties in the features' order, as Qurro's stable sort leaves
# them), at x = 1, 2, ..., m. `num` and `den` are the 1-based indices of
# the selected features (or NULL); `reject` the DA decisions (or NULL); and
# `labels` how many features to name, the selected ones (else the rejected
# ones, else, for `name_any`, all) with the largest |effect|.
#
# Returns a data frame with a row per bar: x, effect, feature, class
# ('Numerator', 'Denominator', 'Rejected' or 'None'), col, and labelled.
# Its `omitted` attribute counts the features left out.
rank_layout <- function (effect, feature, reject = NULL, num = NULL, den = NULL,
                         labels = 10, name_any = FALSE) {

  ok  <- which(is.finite(effect))
  o   <- ok[order(effect[ok])]
  cls <- rep('None', length(effect))
  if (!is.null(reject)) cls[reject %in% TRUE] <- 'Rejected'
  sel <- !is.null(num) || !is.null(den)
  if (sel) {
    cls[]    <- 'None'
    cls[num] <- 'Numerator'
    cls[den] <- 'Denominator'
  }

  # Unselected bars are light; with neither a selection nor decisions to
  # show (ordinations), all are mid-grey.
  pal <- c(
    None = if (is.null(reject) && !sel) PLOT_MUTED else 'grey85',
    Rejected = 'grey30', Numerator = PLOT_HIGHLIGHT, Denominator = 'royalblue3' )

  out <- data.frame(
    x = seq_along(o), effect = effect[o], feature = feature[o],
    class = cls[o], stringsAsFactors = FALSE )
  out$col <- unname(pal[out$class])

  # The features to name: the largest |effect| among the candidates.
  cand <- if (sel) out$class %in% c('Numerator', 'Denominator')
          else if (!is.null(reject)) out$class == 'Rejected'
          else rep(name_any, nrow(out))
  i <- which(cand)
  i <- utils::head(i[order(-abs(out$effect[i]))], labels)
  out$labelled <- seq_len(nrow(out)) %in% i

  attr(out, 'omitted') <- length(effect) - length(ok)
  out
}


# Draws a rank plot from rank_layout()'s `lay`: a bar per feature from 0
# to its effect, a line at 0, the labelled features named at their bars'
# ends, and with a selection, a legend of the numerator and denominator at
# `legend`. `args` are the user's arguments for graphics::plot.default().
plot_rank_bars <- function (lay, args, main, ylab, sub = NULL, legend = 'topleft') {

  m <- nrow(lay)
  if (m == 0)
    stop('No feature has a finite effect to rank.', call. = FALSE)

  # Room for the names at the bars' ends.
  ylim <- range(0, lay$effect)
  if (any(lay$labelled)) ylim <- ylim + c(-1, 1) * 0.05 * diff(ylim)

  plot_do(
    graphics::plot.default, c(list(NA_real_, NA_real_), args),
    list(
      xlim = c(0.5, m + 0.5), ylim = ylim, las = 1,
      xaxs = 'i', main = main, sub = sub, ylab = ylab,
      xlab = sprintf('Feature rank (of %d)', m) ))

  graphics::rect(lay$x - 0.5, 0, lay$x + 0.5, lay$effect, col = lay$col, border = NA)
  graphics::abline(h = 0, col = 'grey40')

  i <- which(lay$labelled)
  if (length(i))
    plot_label(
      lay$x[i], lay$effect[i], lay$feature[i], col = 'black', cex = 0.7,
      pos = ifelse(lay$effect[i] >= 0, 3, 1) )

  k <- c(sum(lay$class == 'Numerator'), sum(lay$class == 'Denominator'))
  if (any(k > 0))
    plot_legend(
      legend, sprintf(c('Numerator (%d)', 'Denominator (%d)'), k),
      c(PLOT_HIGHLIGHT, 'royalblue3'), 15 )

  invisible()
}


# The 1-based indices among `feature` of a rank plot's `numerator` and
# `denominator` (either may be NULL): any form feature_set() takes, matched
# to the result's features, and not overlapping. Returns list(num, den).
rank_selection <- function (numerator, denominator, feature) {
  D    <- length(feature)
  feat <- if (is.character(feature)) feature else NULL
  num  <- if (!is.null(numerator))   feature_set(numerator,   feat, D, 'numerator',   'the result')
  den  <- if (!is.null(denominator)) feature_set(denominator, feat, D, 'denominator', 'the result')
  shared <- intersect(num, den)
  if (length(shared))
    stop(sprintf(
      '`numerator` and `denominator` must not share features, but they share %s: %s.',
      if (length(shared) == 1L) '1 feature' else paste(length(shared), 'features'),
      name_list(feature[shared]) ), call. = FALSE )
  list(num = num, den = den)
}


# The subtitle's note of features left out of a rank plot, or NULL.
rank_omitted <- function (lay, what = 'no estimate') {
  k <- attr(lay, 'omitted', exact = TRUE)
  if (!isTRUE(k > 0)) return (NULL)
  sprintf('%d %s with %s not shown', k, if (k == 1) 'feature' else 'features', what)
}


# A DA result's rank plot, from da_plot()'s per-feature vectors.
da_plot_rank <- function (x, est, reject, feature, labels, args, num, den, what, est_lab) {
  sel <- rank_selection(num, den, feature)
  lay <- rank_layout(est, feature, reject, sel$num, sel$den, labels)
  sub <- da_plot_sub(what, reject[is.finite(est)], extra = rank_omitted(lay))
  plot_rank_bars(lay, args, main = da_plot_title(x, 'rank'), ylab = est_lab, sub = sub)
  invisible()
}


# An ordination's rank plot: its loadings on axis `axes[[1]]`.
ord_plot_rank <- function (x, axes, labels, args, num, den, legend) {
  src <- rank_source(x, axes = axes[[1]])
  sel <- rank_selection(num, den, src$feature)
  lay <- rank_layout(src$effect, src$feature, NULL, sel$num, sel$den, labels, name_any = TRUE)
  plot_rank_bars(
    lay, args,
    main   = paste('Rank plot of', plot_title(x, ord_method_name(attr(x, 'method', exact = TRUE)))),
    ylab   = sprintf('Loading on %s', src$axis),
    sub    = rank_omitted(lay, 'a non-finite loading'),
    legend = legend )
  invisible()
}


# `numerator` and `denominator` are only for rank plots.
rank_only <- function (type, num, den) {
  if (!identical(type, 'rank') && (!is.null(num) || !is.null(den)))
    stop('`numerator` and `denominator` are for `type = "rank"`.', call. = FALSE)
  invisible()
}
