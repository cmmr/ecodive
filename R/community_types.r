# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Community types: the result that enterotype() and valencia() (and the
# planned dmm()) return, with its print() and plot() methods. See
# ?ecodive_types.



#' Community types
#'
#' [enterotype()] and [valencia()] (and the planned `dmm()`) assign each
#' sample to one of a few community types, and return a data frame of class
#' `c("ecodive_<method>", "ecodive_types", "data.frame")`, with one row per
#' sample, in the input's order, and the sample names as row names.
#'
#' @section Columns:
#' * `sample`: the sample's name (or its position, for unnamed samples).
#' * `type`: the type assigned to it, a factor whose levels are the types in
#'   order: `"1"`, `"2"`, ... for [enterotype()]; the sub-CSTs `"I-A"` to
#'   `"V"` for [valencia()].
#' * Columns of the method: for [enterotype()], `silhouette`, the sample's
#'   silhouette width; for [valencia()], `cst` (the higher-order community
#'   state type), `score` (the similarity to the assigned centroid) and
#'   `matched`.
#'
#' @section Attributes:
#' * `types`: a data frame, one row per type: `type`, `n_samples`,
#'   `top_feature` (the feature of largest mean proportion, Arumugam et
#'   al.'s "driver"; `NA` without abundances), and columns of the method,
#'   such as [enterotype()]'s `medoid` and `silhouette`, or [valencia()]'s
#'   `cst` and reference scores.
#' * `profiles`: a types \eqn{\times} features matrix of each type's mean
#'   proportions (`NULL` without abundances); for [valencia()], the
#'   centroids.
#' * `selection`: a data frame, one row per number of types tried, with the
#'   criterion that chose among them, and `selected`, `TRUE` for the one
#'   chosen. Absent (`NULL`) for [valencia()], whose centroids fix the
#'   types.
#' * Attributes of the method, such as [valencia()]'s `scores` (every
#'   sample's similarity to every centroid).
#' * `method`: the function that made the result, such as `"enterotype"`.
#'
#' @section Using the types:
#' `res$type` is a plain factor, so it can be given as `groups` to any test,
#' such as [permanova()], or to an ordination's `plot()`, to colour the
#' samples by type.
#'
#' @name ecodive_types
#' @seealso [enterotype()], [valencia()], [plot.ecodive_types()]
#' @examples
#'     et <- enterotype(iris[, 1:4], n_types = 2:5)
#'     et
#'     table(et$type)
#'     attr(et, 'types')
#'     attr(et, 'selection')
#'
#'     # The types as groups.
#'     permanova(jensen(iris[, 1:4]), et$type)
NULL



# Facts about each method that makes community types (README: "Facts about a
# method live in a registry keyed by `method`"): its display name, and the
# `selection` column that chose the number of types, with its label and
# whether its largest ('max') or smallest ('min') value wins.
# valencia() has no criterion: its centroids fix the number of types, and
# its `selection` is NULL.
COMMUNITY_TYPES <- data.frame(
  stringsAsFactors = FALSE,
  id              = c('enterotype',              'valencia'),
  name            = c('Enterotypes',             'VALENCIA Community State Types'),
  criterion       = c('ch',                      NA),
  criterion_label = c('Calinski-Harabasz index', NA),
  best            = c('max',                     NA) )


community_types_row <- function (x) {
  i <- match(attr(x, 'method', exact = TRUE), COMMUNITY_TYPES$id)
  if (is.na(i)) NULL else COMMUNITY_TYPES[i, , drop = FALSE]
}



# The result: `df` (sample, type, method columns) with its attributes, in
# the order listed in ?ecodive_types, and its classes.
new_ecodive_types <- function (df, method, types, profiles, selection, ...) {
  rownames(df) <- df$sample
  extra <- list(...)
  attr(df, 'types')     <- types
  attr(df, 'profiles')  <- profiles
  attr(df, 'selection') <- selection
  for (nm in names(extra)) attr(df, nm) <- extra[[nm]]
  attr(df, 'method')    <- method
  class(df) <- c(paste0('ecodive_', method), 'ecodive_types', 'data.frame')
  df
}



# Each type's mean proportions, from the samples' proportions `p` (samples x
# features) and their types `type` (a factor). A types x features matrix.
types_profiles <- function (p, type) {
  k   <- nlevels(type)
  out <- matrix(0, k, ncol(p), dimnames = list(levels(type), colnames(p)))
  for (t in seq_len(k)) {
    rows <- which(as.integer(type) == t)
    if (length(rows)) out[t, ] <- colMeans(p[rows, , drop = FALSE])
  }
  out
}


# The feature of largest mean proportion in each type: its name, or its
# position when the features are unnamed. NA without profiles.
types_top_feature <- function (profiles, k) {
  if (is.null(profiles)) return (rep(NA_character_, k))
  nms <- colnames(profiles)
  if (is.null(nms)) nms <- as.character(seq_len(ncol(profiles)))
  nms[apply(profiles, 1L, which.max)]
}



#------------------------------------------------------------------------------
# print()
#------------------------------------------------------------------------------

#' @rdname print_ecodive
#' @export
print.ecodive_types <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  sym   <- print_sym()
  types <- attr(x, 'types', exact = TRUE)
  sel   <- attr(x, 'selection', exact = TRUE)
  reg   <- community_types_row(x)
  title <- display_title(x)
  if (is.null(title) && !is.null(reg)) title <- reg$name

  k    <- if (is.data.frame(types)) nrow(types) else nlevels(x$type)
  used <- if (is.data.frame(types)) sum(types$n_samples > 0) else length(unique(stats::na.omit(x$type)))
  out  <- character(0)

  how <- ''
  if (is.data.frame(sel) && nrow(sel) > 1 && !is.null(reg) && reg$criterion %in% names(sel))
    how <- sprintf(
      ', chosen by the %s from %s',
      reg$criterion_label, print_range(sel$n_types) )
  if (used < k) how <- sprintf(' (%d with samples)%s', used, how)
  n_na <- sum(is.na(x$type))
  if (n_na) how <- sprintf('%s; %d without a type', how, n_na)
  out <- c(out, sprintf('# %s', title), sprintf('# %d samples in %d types%s', nrow(x), k, how))

  if (is.data.frame(types))
    out <- c(out, '#', print_table_lines(
      types, 'Types:', n = Inf, width = width, digits = digits, types = FALSE ))

  if (is.data.frame(sel) && nrow(sel) > 1)
    out <- c(out, '#', print_table_lines(
      sel, 'Selection:', n = Inf, width = width, digits = digits, types = FALSE ))

  df <- x
  class(df) <- 'data.frame'
  out <- c(out, '#', print_table_lines(
    df, sprintf('Samples: a data.frame: %d %s %d', nrow(df), sym$times, length(df)),
    n = n, width = width, digits = digits,
    attrs = print_attrs(x, c('names', 'row.names', 'class')) ))

  cat(out, sep = '\n')
  invisible(x)
}


# "2 to 10", "2, 4 and 7", or "3".
print_range <- function (v) {
  v <- sort(unique(v))
  if (length(v) == 1) return (as.character(v))
  if (all(diff(v) == 1)) return (sprintf('%d to %d', v[[1]], v[[length(v)]]))
  paste(paste(v[-length(v)], collapse = ', '), 'and', v[[length(v)]])
}



#------------------------------------------------------------------------------
# plot()
#------------------------------------------------------------------------------

#' Plot community types
#'
#' Draws a result of [enterotype()] or [valencia()] with base graphics.
#'
#' * `"profiles"` (the default): stacked bars of each type's composition
#'   (for [enterotype()], its samples' mean proportions; for [valencia()],
#'   the centroids), with the `n_features` features that reach the largest
#'   proportion in any type, and the rest pooled as "Other". The bars are
#'   labelled with the type and its number of samples.
#' * `"selection"`: the criterion that chose the number of types (for
#'   [enterotype()], the Calinski-Harabasz index) against the number of
#'   types, with the one chosen circled. Not for [valencia()], whose
#'   centroids fix the types.
#' * `"scores"` ([valencia()] only): boxplots of the samples' scores (their
#'   similarity to their sub-CST's centroid), one per sub-CST with samples,
#'   in VALENCIA's colours, each beside a diamond and bar for the mean
#'   \eqn{\pm} standard deviation of VALENCIA's reference samples, as
#'   VALENCIA's own plot shows them.
#'
#' The types keep the order of `levels(x$type)` in every plot.
#'
#' @param x  A result of [enterotype()] or [valencia()].
#'
#' @param type  `"profiles"`, `"selection"` or `"scores"`.
#'
#' @param n_features  The number of features to show in `"profiles"`.
#'        Default `10`.
#'
#' @param ...  Passed on to [graphics::barplot()],
#'        [graphics::plot.default()] or [graphics::boxplot()], such as
#'        `main`, `ylab` or `col`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     et <- enterotype(iris[, 1:4], n_types = 2:5)
#'     plot(et)
#'     plot(et, 'selection')
#'
#'     # VALENCIA on thirteen samples made from its centroids.
#'     counts <- round(valencia_centroids * 10000)
#'     rownames(counts) <- paste0('S', 1:13)
#'     vl <- valencia(counts)
#'     plot(vl)
#'     plot(vl, 'scores')
plot.ecodive_types <- function (x, type = c('profiles', 'selection', 'scores'), n_features = 10, ...) {

  type <- match.arg(type)
  args <- list(...)
  reg  <- community_types_row(x)
  main <- plot_title(x, if (!is.null(reg)) reg$name else 'Community types')

  if (type == 'scores') return (invisible(plot_types_scores(x, args, main)))

  if (type == 'selection') {

    sel <- attr(x, 'selection', exact = TRUE)
    if (is.null(reg) || !is.data.frame(sel) || !reg$criterion %in% names(sel))
      stop('This result has no criterion to plot.', call. = FALSE)

    v <- sel[[reg$criterion]]
    i <- which(sel$selected)

    plot_do(
      graphics::plot.default, c(list(sel$n_types, v), args),
      list(
        type = 'b', pch = 16, xaxt = 'n', las = 1, main = main,
        sub  = sprintf('Chosen: %d types', sel$n_types[i]),
        xlab = 'Number of types',
        ylab = reg$criterion_label ))

    if (plot_x_axis(args)) graphics::axis(1, at = sel$n_types)
    if (length(i))
      graphics::points(sel$n_types[i], v[i], cex = 2.5, lwd = 2, col = PLOT_HIGHLIGHT)

    return (invisible(x))
  }


  # Profiles: stacked bars.
  prof <- attr(x, 'profiles', exact = TRUE)
  if (is.null(prof))
    stop(
      'This result has no profiles: give enterotype() a count table, ',
      'or `counts` along with distances.', call. = FALSE )

  n_features <- plot_count(n_features, 'n_features')
  if (is.null(colnames(prof))) colnames(prof) <- as.character(seq_len(ncol(prof)))

  top  <- apply(prof, 2L, max)
  keep <- utils::head(order(-top, seq_along(top)), min(n_features, ncol(prof)))
  bars <- prof[, keep, drop = FALSE]
  rest <- 1 - rowSums(bars)
  if (any(rest > 1e-12)) bars <- cbind(bars, Other = pmax(rest, 0))

  types <- attr(x, 'types', exact = TRUE)
  sizes <- if (is.data.frame(types)) types$n_samples else as.vector(table(x$type))
  lab   <- sprintf('%s\n(n = %d)', rownames(prof), sizes)
  cex_n <- if (nrow(prof) > 8) 0.7 else 1

  col <- plot_palette(length(keep))
  if (ncol(bars) > length(keep)) col <- c(col, 'grey85')

  # Room for the legend in the right margin.
  leg  <- rev(colnames(bars))
  lcol <- rev(col)
  mar  <- graphics::par('mar')
  w    <- max(graphics::strwidth(leg, units = 'inches', cex = 0.8)) + 0.45
  lines_w <- w / graphics::par('csi')
  op <- graphics::par(mar = c(mar[1:3], max(mar[4], lines_w + 1)))
  on.exit(graphics::par(op), add = TRUE)

  plot_do(
    graphics::barplot, c(list(t(bars)), args),
    list(
      col = col, names.arg = lab, las = 1, border = NA, cex.names = cex_n,
      main = main, xlab = 'Type', ylab = 'Mean proportion',
      ylim = c(0, 1) ))

  graphics::legend(
    x = graphics::par('usr')[2], y = graphics::par('usr')[4],
    legend = leg, fill = lcol, border = NA, bty = 'n', cex = 0.8, xpd = NA )

  invisible(x)
}



# plot(type = 'scores'), for valencia(): each sub-CST's scores (boxes, in
# VALENCIA's colours), beside VALENCIA's reference mean +/- sd (diamond and
# bar), for the sub-CSTs with samples, as Valencia.py's plot draws them.
plot_types_scores <- function (x, args, main) {

  types <- attr(x, 'types', exact = TRUE)
  if (is.null(attr(x, 'scores', exact = TRUE)) || !is.data.frame(types) ||
      !'score' %in% names(x))
    stop('Only valencia() results have scores to plot.', call. = FALSE)

  lev <- as.character(types$type[types$n_samples > 0])
  if (!length(lev))
    stop('No sample was assigned a type, so there are no scores to plot.', call. = FALSE)

  vals <- lapply(lev, function (l) x$score[which(x$type == l)])
  at   <- seq_along(lev)
  i    <- match(lev, as.character(types$type))
  m    <- types$reference_score_mean[i]
  s    <- types$reference_score_sd[i]
  col  <- valencia_colours(lev)

  plot_do(
    graphics::boxplot, c(list(vals), args),
    list(
      at = at - 0.15, boxwex = 0.4, col = col, names = lev, xaxt = 'n',
      xlim = c(0.5, length(lev) + 0.5), ylim = c(0, 1), las = 1,
      main = main, xlab = 'Sub-CST', ylab = 'Similarity to assigned sub-CST' ))

  if (plot_x_axis(args)) graphics::axis(1, at = at, labels = lev)

  ref <- which(!is.na(m))
  if (length(ref)) {
    graphics::segments(at[ref] + 0.25, m[ref] - s[ref], at[ref] + 0.25, m[ref] + s[ref], lwd = 1.5)
    graphics::points(at[ref] + 0.25, m[ref], pch = 23, bg = col[ref], cex = 1.2)
    graphics::legend(
      'bottomleft', legend = c('These samples', 'VALENCIA reference: mean +/- sd'),
      pch = c(22, 23), pt.bg = 'grey85', pt.cex = c(2, 1.2), lty = c(NA, 1),
      bty = 'n', cex = 0.8 )
  }

  x
}
