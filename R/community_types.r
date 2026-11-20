# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Community types: the result that enterotype() (and the planned valencia()
# and dmm()) return, with its print() and plot() methods. See
# ?ecodive_types.



#' Community types
#'
#' [enterotype()] (and the planned `valencia()` and `dmm()`) assign each
#' sample to one of a few community types, and return a data frame of class
#' `c("ecodive_<method>", "ecodive_types", "data.frame")`, with one row per
#' sample, in the input's order, and the sample names as row names.
#'
#' @section Columns:
#' * `sample`: the sample's name (or its position, for unnamed samples).
#' * `type`: the type assigned to it, a factor whose levels (`"1"`, `"2"`,
#'   ...) are the types in order.
#' * Columns of the method: for [enterotype()], `silhouette`, the sample's
#'   silhouette width.
#'
#' @section Attributes:
#' * `types`: a data frame, one row per type: `type`, `n_samples`,
#'   `top_feature` (the feature of largest mean proportion, Arumugam et
#'   al.'s "driver"; `NA` without abundances), and columns of the method,
#'   such as [enterotype()]'s `medoid` and `silhouette`.
#' * `profiles`: a types \eqn{\times} features matrix of each type's mean
#'   proportions (`NULL` without abundances).
#' * `selection`: a data frame, one row per number of types tried, with the
#'   criterion that chose among them, and `selected`, `TRUE` for the one
#'   chosen.
#' * `method`: the function that made the result, such as `"enterotype"`.
#'
#' @section Using the types:
#' `res$type` is a plain factor, so it can be given as `groups` to any test,
#' such as [permanova()], or to an ordination's `plot()`, to colour the
#' samples by type.
#'
#' @name ecodive_types
#' @seealso [enterotype()], [plot.ecodive_types()]
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
COMMUNITY_TYPES <- data.frame(
  stringsAsFactors = FALSE,
  id              = c('enterotype'),
  name            = c('Enterotypes'),
  criterion       = c('ch'),
  criterion_label = c('Calinski-Harabasz index'),
  best            = c('max') )


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

  k   <- if (is.data.frame(types)) nrow(types) else nlevels(x$type)
  out <- character(0)

  how <- ''
  if (is.data.frame(sel) && nrow(sel) > 1 && !is.null(reg) && reg$criterion %in% names(sel))
    how <- sprintf(
      ', chosen by the %s from %s',
      reg$criterion_label, print_range(sel$n_types) )
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
#' Draws a result of [enterotype()] with base graphics.
#'
#' * `"profiles"` (the default): stacked bars of each type's mean
#'   composition, with the `n_features` features that reach the largest
#'   mean proportion in any type, and the rest pooled as "Other". The bars
#'   are labelled with the type and its number of samples.
#' * `"selection"`: the criterion that chose the number of types (for
#'   [enterotype()], the Calinski-Harabasz index) against the number of
#'   types, with the one chosen circled.
#'
#' The types keep the order and colours of `levels(x$type)` in every plot.
#'
#' @param x  A result of [enterotype()].
#'
#' @param type  `"profiles"` or `"selection"`.
#'
#' @param n_features  The number of features to show in `"profiles"`.
#'        Default `10`.
#'
#' @param ...  Passed on to [graphics::barplot()] or
#'        [graphics::plot.default()], such as `main`, `ylab` or `col`.
#'
#' @return `x`, invisibly.
#'
#' @export
#' @examples
#'     et <- enterotype(iris[, 1:4], n_types = 2:5)
#'     plot(et)
#'     plot(et, 'selection')
plot.ecodive_types <- function (x, type = c('profiles', 'selection'), n_features = 10, ...) {

  type <- match.arg(type)
  args <- list(...)
  reg  <- community_types_row(x)
  main <- plot_title(x, if (!is.null(reg)) reg$name else 'Community types')

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
      col = col, names.arg = lab, las = 1, border = NA,
      main = main, xlab = 'Type', ylab = 'Mean proportion',
      ylim = c(0, 1) ))

  graphics::legend(
    x = graphics::par('usr')[2], y = graphics::par('usr')[4],
    legend = leg, fill = lcol, border = NA, bty = 'n', cex = 0.8, xpd = NA )

  invisible(x)
}
