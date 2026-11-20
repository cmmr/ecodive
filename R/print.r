# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Print ecodive's results
#'
#' ecodive's tables, ordinations and distance matrices print a preview, in
#' the manner of the tibble package, rather than every value: a header
#' saying what the object is (a data frame, matrix or `dist`) and how big,
#' and its title (see [ecodive_display]); the first rows, with each
#' column's type under its name; as many columns as fit the console; and a
#' footer counting the rows and columns left out and naming the object's
#' other attributes. The object itself is unchanged: `as.data.frame(x)`,
#' `as.matrix(x)` or `unclass(x)` print it in full.
#'
#' @param x  A result of class `ecodive_ord` (an ordination), `ecodive_dist`
#'        (a distance matrix), or one of ecodive's data frames: from
#'        [ancom()], [ancombc()], [ancombc2()], [ancombc_global()],
#'        [ancombc_pairwise()], [ancombc_dunnett()], [ancombc_trend()],
#'        [dirmult_ttest()], [da_sensitivity()], [simper()], [bioenv()],
#'        [envfit()], [assoc_pairs()] or [enterotype()], or a source-tracking
#'        matrix of class `ecodive_source` (from [sourcetracker()]). A community-types result (see
#'        [ecodive_types]) prints its types and how their number was chosen
#'        before its first samples.
#'
#' @param n  The number of rows to print (samples, for a distance matrix).
#'        `NULL`, the default, prints every row of up to 20, and the first
#'        10 of more; `Inf` prints them all.
#'
#' @param width  The width of the output, in characters: columns that don't
#'        fit are left out and named in the footer. Default
#'        `getOption("width")`.
#'
#' @param digits  The significant digits of each number. Default `3`.
#'
#' @param ...  Not used.
#'
#' @return `x`, invisibly.
#'
#' @name print_ecodive
#' @examples
#'     # An ordination: the first rows of its coordinates.
#'     ord <- pcoa(bray(iris[, 1:4]))
#'     ord
#'
#'     # A distance matrix: the top-left corner of its lower triangle.
#'     bray(iris[, 1:4])
#'
#'     # The values as a plain matrix.
#'     head(unclass(ord))
NULL



# The rows of `total` to print, by `n` (see ?print_ecodive).
print_n_rows <- function (total, n) {
  if (is.null(n)) return (if (total <= 20) total else min(total, 10))
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0)
    stop('`n` must be a non-negative number, or Inf.', call. = FALSE)
  as.integer(min(total, n))
}


# The output width, from `width`.
print_width <- function (width) {
  if (is.null(width)) width <- getOption('width', 80L)
  if (!is.numeric(width) || length(width) != 1 || is.na(width) || width < 10)
    stop('`width` must be a number of at least 10.', call. = FALSE)
  as.integer(width)
}


# The multiplication sign and ellipsis, where the console can show them.
print_sym <- function () {
  if (isTRUE(l10n_info()[['UTF-8']])) list(times = '\u00d7', more = '\u2026')
  else list(times = 'x', more = '...')
}


# A column's type, abbreviated as tibble does: "dbl", "int", "chr", ...
print_type <- function (v) {
  if (is.ordered(v))   return ('ord')
  if (is.factor(v))    return ('fct')
  if (inherits(v, 'Date')) return ('date')
  if (is.matrix(v) || is.data.frame(v)) return ('mat')
  switch(
    typeof(v),
    double    = 'dbl',
    integer   = 'int',
    character = 'chr',
    logical   = 'lgl',
    complex   = 'cpl',
    list      = 'list',
    class(v)[[1]] )
}


# Strings no wider than `w`, cut with an ellipsis (one character wide in
# UTF-8, "..." otherwise).
print_cut <- function (s, w = 40L) {
  long <- !is.na(s) & nchar(s, type = 'width') > w
  if (any(long)) {
    more <- print_sym()$more
    s[long] <- paste0(substr(s[long], 1L, w - nchar(more, type = 'width')), more)
  }
  s
}


# The cells of column `v`, as strings, and whether to right-align them.
# Numbers have `digits` significant digits each, lined up on the decimal
# point unless any is in scientific notation. Cells where `blank` is TRUE
# are empty, and don't count towards the alignment.
print_cells <- function (v, digits, blank = NULL) {

  if (any(blank)) {
    cl <- print_cells(v[!blank], digits)
    s  <- rep('', length(v))
    s[!blank] <- cl$s
    return (list(s = s, right = cl$right))
  }

  if (is.matrix(v) || is.data.frame(v)) {
    k <- NCOL(v)
    return (list(s = rep(sprintf('<%s [%d]>', print_type(v), k), NROW(v)), right = FALSE))
  }

  if (is.list(v)) {
    s <- vapply(v, function (z) sprintf('<%s [%d]>', class(z)[[1]], length(z)), '')
    return (list(s = s, right = FALSE))
  }

  if (is.numeric(v) && !is.factor(v)) {
    s <- if (is.double(v)) as.character(signif(v, digits)) else as.character(v)
    s[is.na(s)] <- 'NA'
    if (is.double(v) && !any(grepl('e', s, fixed = TRUE))) {
      dot <- regexpr('.', s, fixed = TRUE)
      int <- ifelse(dot > 0, substr(s, 1L, dot - 1L), s)
      dec <- ifelse(dot > 0, substring(s, dot), '')
      s   <- paste0(
        formatC(int, width = max(0L, nchar(int))),
        formatC(dec, width = -max(0L, nchar(dec))) )
    }
    return (list(s = s, right = TRUE))
  }

  if (is.logical(v)) {
    s <- as.character(v)
    s[is.na(s)] <- 'NA'
    return (list(s = s, right = TRUE))
  }

  s <- as.character(v)
  s[is.na(s)] <- 'NA'
  list(s = print_cut(s), right = FALSE)
}


# Pads strings to width `w`, on the left (`right`) or the right.
print_pad <- function (s, w, right) {
  gap <- pmax(0L, w - nchar(s, type = 'width'))
  pad <- strrep(' ', gap)
  if (right) paste0(pad, s) else paste0(s, pad)
}


# The lines that print the data frame `df` (a plain one): `head`, lines
# that go after "# " at the top; a column-name row, a type row (`types`),
# and the first rows (by `n`), with as many columns as fit `width`; then
# the rows and columns left out, and `attrs`, the names of the object's
# other attributes. `rows` names the rows ("rows", "samples", ...);
# `extra_rows` more were left out before `df` was made. `blank`, a logical
# matrix the shape of `df`, marks cells to leave empty.
print_table_lines <- function (
    df, head, n = NULL, width = NULL, digits = 3L, attrs = NULL,
    types = TRUE, rows = 'rows', extra_rows = 0L, blank = NULL ) {

  sym   <- print_sym()
  width <- print_width(width)
  nr    <- nrow(df)
  nc    <- length(df)
  k     <- print_n_rows(nr, n)
  sub   <- df[seq_len(k), , drop = FALSE]

  # Row labels: the row names, or else the row numbers.
  named <- .row_names_info(df) > 0L
  rl    <- if (named) print_cut(rownames(sub), 30L) else as.character(seq_len(k))
  rw    <- max(0L, nchar(rl, type = 'width'))

  # The columns that fit, each as wide as its name, type and cells.
  cols <- list()
  used <- rw
  for (j in seq_len(nc)) {
    cl  <- print_cells(sub[[j]], digits, if (!is.null(blank)) blank[seq_len(k), j])
    nm  <- names(df)[[j]]
    ty  <- if (types) paste0('<', print_type(df[[j]]), '>') else ''
    w   <- max(nchar(c(nm, ty, cl$s), type = 'width'))
    if (length(cols) && used + 1L + w > width) break
    cols[[length(cols) + 1L]] <- list(name = nm, type = ty, s = cl$s, right = cl$right, w = w)
    used <- used + 1L + w
  }

  head <- head[!is.na(head) & nzchar(head)]
  out  <- paste('#', head)
  if (length(cols)) {
    row_line <- function (lab, cells) paste0(
      print_pad(lab, rw, !named), ' ',
      paste(cells, collapse = ' ') )
    out <- c(out, row_line('', vapply(cols, function (cc) print_pad(cc$name, cc$w, cc$right), '')))
    if (types)
      out <- c(out, row_line('', vapply(cols, function (cc) print_pad(cc$type, cc$w, cc$right), '')))
    for (i in seq_len(k))
      out <- c(out, row_line(rl[[i]], vapply(cols, function (cc) print_pad(cc$s[[i]], cc$w, cc$right), '')))
  }
  out <- sub(' +$', '', out)

  # What was left out.
  more <- nr - k + extra_rows
  if (more > 0)
    out <- c(out, sprintf('# %s %d more %s', sym$more, more, if (more == 1) sub('s$', '', rows) else rows))
  left <- setdiff(seq_len(nc), seq_along(cols))
  if (length(left)) {
    what <- vapply(utils::head(left, 20L), function (j) {
      nm <- names(df)[[j]]
      if (types) sprintf('%s <%s>', nm, print_type(df[[j]])) else nm
    }, '')
    if (length(left) > 20L) what <- c(what, sprintf('and\001%d\001more', length(left) - 20L))
    out <- c(out, print_wrap(sprintf(
      '%s %d more column%s: %s', sym$more, length(left),
      if (length(left) == 1) '' else 's', paste(what, collapse = ', ') ), width))
  }
  if (more > 0)
    out <- c(out, sprintf('# Use print(x, n = ...) to see more %s', rows))
  if (length(attrs))
    out <- c(out, print_wrap(paste('Attributes:', paste(attrs, collapse = ', ')), width))
  out
}


# A footer line, wrapped to `width` after "# ", without parting a column's
# name from its type.
print_wrap <- function (s, width) {
  s <- gsub(' <', '\001<', s, fixed = TRUE)
  gsub('\001', ' ', strwrap(s, width = width - 1L, exdent = 2L, prefix = '# '), fixed = TRUE)
}


# The names of the attributes of `x` besides `skip` (and `display`).
print_attrs <- function (x, skip) {
  setdiff(names(attributes(x)), c(skip, 'display'))
}



#------------------------------------------------------------------------------
# Data frames.
#------------------------------------------------------------------------------

# What a data frame result is, for its header: its title (?ecodive_display),
# or else the name its plot() gives it.
print_df_label <- function (x) {

  title <- display_title(x)
  if (!is.null(title)) return (title)

  cls <- class(x)[[1]]
  if (startsWith(cls, 'ecodive_ancombc')) {
    # Subsetting rows with `[` drops the `method` attribute.
    if (is.null(attr(x, 'method', exact = TRUE))) return (NULL)
    nm   <- ancombc_plot_name(x)
    test <- sub('^ecodive_ancombc_?', '', cls)
    return (switch(
      test,
      global   = paste(nm, 'global test'),
      pairwise = paste(nm, 'pairwise tests'),
      dunnett  = paste(nm, "Dunnett's tests"),
      trend    = paste(nm, 'trend test'),
      nm ))
  }

  switch(
    cls,
    ecodive_ancom = {
      g <- attr(x, 'groups', exact = TRUE)
      if (is.character(g) && length(g)) paste('ANCOM:', paste(g, collapse = ' vs ')) else 'ANCOM'
    },
    ecodive_dirmult_ttest = 'Dirichlet-multinomial t-test',
    ecodive_da_wilcoxon   = ,
    ecodive_da_ttest      = ,
    ecodive_da_lm         = ,
    ecodive_linda         = {
      # Subsetting rows with `[` drops the `method` attribute.
      m <- attr(x, 'method', exact = TRUE)
      if (is.null(m)) NULL else da_method_name(m)
    },
    ecodive_da_sensitivity = {
      m <- attr(x, 'method', exact = TRUE)
      if (is.null(m)) NULL else paste('Sensitivity of', da_method_name(m))
    },
    ecodive_assoc_pairs   = assoc_pairs_label(x),
    ecodive_simper        = 'SIMPER',
    ecodive_bioenv        = if (isTRUE(attr(x, 'partial', exact = TRUE))) 'BIO-ENV (partial)' else 'BIO-ENV',
    NULL )
}


# The print() method of every data frame result.
print_table_result <- function (x, n = NULL, width = NULL, digits = 3L, ...) {
  df <- x
  class(df) <- 'data.frame'
  sym  <- print_sym()
  head <- c(
    sprintf('A data.frame: %d %s %d', nrow(df), sym$times, length(df)),
    print_df_label(x) )
  cat(
    print_table_lines(
      df, head, n = n, width = width, digits = digits,
      attrs = print_attrs(x, c('names', 'row.names', 'class')) ),
    sep = '\n' )
  invisible(x)
}

#' @export
print.ecodive_ancom            <- print_table_result
#' @export
print.ecodive_assoc_pairs      <- print_table_result
#' @export
print.ecodive_ancombc          <- print_table_result
#' @export
print.ecodive_ancombc_global   <- print_table_result
#' @export
print.ecodive_ancombc_pairwise <- print_table_result
#' @export
print.ecodive_ancombc_dunnett  <- print_table_result
#' @export
print.ecodive_ancombc_trend    <- print_table_result
#' @export
print.ecodive_dirmult_ttest    <- print_table_result
#' @export
print.ecodive_da_wilcoxon      <- print_table_result
#' @export
print.ecodive_da_ttest         <- print_table_result
#' @export
print.ecodive_da_lm            <- print_table_result
#' @export
print.ecodive_linda            <- print_table_result
#' @export
print.ecodive_da_sensitivity   <- print_table_result
#' @export
print.ecodive_simper           <- print_table_result
#' @export
print.ecodive_bioenv           <- print_table_result
#' @export
print.ecodive_envfit           <- print_table_result



#------------------------------------------------------------------------------
# Distance matrices: the top-left corner of the lower triangle.
#------------------------------------------------------------------------------

#' @rdname print_ecodive
#' @export
print.ecodive_dist <- function (x, n = NULL, width = NULL, digits = 3L, ...) {

  size <- attr(x, 'Size', exact = TRUE)
  unit <- dist_plot_name(x)$unit
  how  <- dist_plot_order_name(dist_clust_of(x))
  head <- c(
    sprintf('A dist: %d %s (lower triangle)', size, if (size == 1) sub('s$', '', unit) else unit),
    display_title(x),
    if (!is.null(how)) paste0(toupper(substr(how, 1, 1)), substring(how, 2)) )

  # Rows 2 to k + 1 and columns 1 to k, as print.dist() shows the triangle,
  # with the upper triangle blank.
  k   <- if (size < 2) 0L else print_n_rows(size - 1L, n)
  lab <- attr(x, 'Labels', exact = TRUE)
  if (is.null(lab)) lab <- as.character(seq_len(size))
  df    <- data.frame()
  blank <- NULL
  if (k > 0) {
    idx   <- seq_len(k + 1L)
    m     <- as.matrix(x)[idx, idx, drop = FALSE][-1L, -(k + 1L), drop = FALSE]
    blank <- upper.tri(m)
    df    <- as.data.frame(unname(m), stringsAsFactors = FALSE)
    names(df)    <- lab[seq_len(k)]
    rownames(df) <- make.unique(lab[idx[-1L]])
  }

  cat(
    print_table_lines(
      df, head, n = Inf, width = width, digits = digits, types = FALSE,
      rows = unit, extra_rows = max(0L, size - 1L - k), blank = blank,
      attrs = print_attrs(x, c('Size', 'Labels', 'Diag', 'Upper', 'class', 'call')) ),
    sep = '\n' )
  invisible(x)
}

