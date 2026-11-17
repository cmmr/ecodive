# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# The shared result of ecodive's permutation tests on a distance matrix
# (mrpp(), anosim(), mantel(), permdisp()): a list of class
# c(<method>, "ecodive_test").
#
# Core fields, first and in this order:
#   statistic  Named numeric, e.g. c(delta = 0.412).
#   p_value    Permutation p-value, or NA.
#   z          Standardized effect size of the statistic against its
#              permuted values, as summary(vegan::permustats()):
#              (statistic - mean(perm)) / sd(perm). NA when undefined.
#   n_perm     Number of permutations.
# Then each method's own fields, and last:
#   permutation  How samples were permuted: 'free', 'within strata', or
#                'supplied matrix'.
#
# Attributes:
#   <name>_perm  The permuted statistics, named after `statistic`, in
#                lower case (delta_perm, r_perm, f_perm), as permanova()'s
#                f_perm.
#   title        The test's name, printed as the heading.
#
# The class names start with "ecodive_": vegan registers S3 methods for
# classes "mrpp", "anosim", and "mantel".
#
# Each method can add lines to the printout with a test_details() method,
# which returns a character vector (one element per line).



#------------------------------------------------------------------------------
# Builds a test result.
#
#   class        The method's class, e.g. 'ecodive_mrpp'.
#   title        The test's name, printed as the heading.
#   statistic    Named numeric scalar, e.g. c(delta = 0.412).
#   p_value, z, n_perm   The core fields.
#   ...          The method's own fields, named, in order.
#   perm         The permuted statistics (numeric vector, length n_perm).
#                Stored as attribute <tolower(names(statistic))>_perm.
#   permutation  'free', 'within strata', or 'supplied matrix'.
#------------------------------------------------------------------------------
new_ecodive_test <- function (
    class, title, statistic, p_value, z, n_perm, ..., perm,
    permutation = 'free' ) {

  stopifnot(
    is.character(class), length(class) >= 1,
    is.numeric(statistic), length(statistic) == 1, !is.null(names(statistic)) )

  res <- list(
    statistic = statistic,
    p_value   = as.double(p_value),
    z         = as.double(z),
    n_perm    = as.integer(n_perm),
    ...,
    permutation = permutation )

  attr(res, paste0(tolower(names(statistic)), '_perm')) <- as.double(perm)
  attr(res, 'title') <- title
  class(res) <- c(class, 'ecodive_test')

  res
}



#------------------------------------------------------------------------------
# Checks the `groups` argument of a test (permanova(), mrpp(), anosim(),
# permdisp(), dist_test()): one value per sample (a vector or factor,
# matched to the dist's `labels` by name when both are named, else in order),
# or, with `data`, the name of one of its columns (or a vector with one value
# per row of `data`). Every test gives the same messages.
#
# Returns list(data, col): pass `data` on to model_design(), which aligns it
# with the samples (a vector becomes a one-column data frame, `group`), and
# take the groups from column `col` of its result. test_group_values() does
# both.
#------------------------------------------------------------------------------
test_groups <- function (groups, data, n, labels, strata = NULL) {

  if (is.null(groups))
    stop('`groups` is required: one value per sample, or a column name of `data`.', call. = FALSE)

  if (!is.atomic(groups) || !is.null(dim(groups)))
    stop(
      '`groups` must be a vector or factor with one value per sample, ',
      'or the name of a column of `data`.', call. = FALSE )


  #--------------------------------------------------------
  # A column of `data`.
  #--------------------------------------------------------
  if (!is.null(data)) {

    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)

    if (is.character(groups) && length(groups) == 1 && !is.na(groups) &&
        groups %in% colnames(data))
      return (list(data = data, col = groups))

    # A vector parallel to the rows of `data`: a new column.
    if (length(groups) == nrow(data) && !(is.character(groups) && length(groups) == 1)) {
      data <- as.data.frame(data, stringsAsFactors = FALSE)
      col  <- make.unique(c(names(data), 'group'))[[ncol(data) + 1L]]
      data[[col]] <- groups
      return (list(data = data, col = col))
    }

    if (is.character(groups) && length(groups) == 1)
      stop(
        '`groups` must be the name of a column of `data`, but ',
        name_list(groups), ' is not one.', call. = FALSE )

    stop(
      'With `data`, `groups` must be the name of a column of `data`, ',
      'or have one value per row of `data` (', nrow(data), ').', call. = FALSE )
  }


  #--------------------------------------------------------
  # One value per sample.
  #--------------------------------------------------------
  if (is.character(groups) && length(groups) == 1 && n != 1)
    stop('`groups` names a column, but `data` is not given.', call. = FALSE)

  nms   <- names(groups)
  named <- !is.null(labels) && !is.null(nms) && !anyNA(nms) && !anyDuplicated(nms)

  if (named && !all(labels %in% nms)) {
    miss <- labels[!labels %in% nms]
    stop(
      sprintf('%d sample%s in `dm` %s no value in `groups`: %s',
        length(miss), if (length(miss) == 1) '' else 's',
        if (length(miss) == 1) 'has' else 'have', name_list(miss) ),
      call. = FALSE )
  }

  if (!named && length(groups) != n)
    stop(sprintf('`groups` has %d values, but `dm` has %d samples.', length(groups), n), call. = FALSE)

  # Values for other samples are ignored, even NAs.
  used <- if (named) groups[match(labels, nms)] else groups
  if (anyNA(used))
    stop('`groups` must not contain NAs.', call. = FALSE)

  if (!is.null(strata) && !(is.atomic(strata) && is.null(dim(strata)) && length(strata) == length(groups)))
    stop('Without `data`, `strata` must have one value per value of `groups`.', call. = FALSE)

  list(data = groups, col = 'group')
}


#------------------------------------------------------------------------------
# test_groups(), then each sample's group in the order of the dist's samples.
# Returns list(groups (the values, with unused factor levels dropped),
# strata (0-based codes, or NULL), col).
#------------------------------------------------------------------------------
test_group_values <- function (groups, data, n, labels, strata = NULL) {

  g   <- test_groups(groups, data, n, labels, strata)

  # `~ 1`: aligned rows and strata, without a model matrix (which fails for a
  # one-level factor).
  des <- model_design(g$data, ~ 1, labels, n, strata, allow_empty = TRUE)
  val <- des$data[[g$col]]

  if (anyNA(val))
    stop('Missing values (NA) in `data` column ', name_list(g$col), '.', call. = FALSE)

  # A factor's NA level (from addNA()) is a group like any other, as in
  # model.matrix(); anyNA() is FALSE for it. Name it "NA", so that it can
  # label results.
  if (is.factor(val) && anyNA(levels(val))) {
    lev <- levels(val)
    lev[is.na(lev)] <- make.unique(c(lev[!is.na(lev)], "NA"))[[sum(!is.na(lev)) + 1L]]
    levels(val) <- lev
  }

  list(groups = val, strata = des$strata, col = g$col)
}



#------------------------------------------------------------------------------
# z-score of a statistic against its permuted values, as
# summary(vegan::permustats()). NA with fewer than 2 permutations or when they
# are all equal.
#------------------------------------------------------------------------------
perm_z <- function (statistic, perm) {
  if (length(perm) < 2 || is.na(statistic)) return (NA_real_)
  z <- (statistic - mean(perm)) / stats::sd(perm)
  if (is.finite(z)) z else NA_real_
}



#------------------------------------------------------------------------------
# Lines a method adds to the printout, after the core results.
#------------------------------------------------------------------------------
test_details <- function (x, digits) UseMethod('test_details')

#' @export
test_details.default <- function (x, digits) character(0)



#' @export
print.ecodive_test <- function (x, digits = max(3L, getOption('digits') - 3L), ...) {

  fmt <- function (v) if (is.na(v)) 'NA' else formatC(v, digits = digits, format = 'g')

  title <- attr(x, 'title', exact = TRUE)
  if (is.null(title)) title <- class(x)[[1]]
  cat(title, '\n\n', sep = '')

  stat <- paste0(names(x$statistic), ' = ', fmt(x$statistic[[1]]))
  if (x$n_perm > 0)
    stat <- paste0(
      stat, ', p = ', if (is.na(x$p_value)) 'NA' else format.pval(x$p_value, digits = digits),
      ', z = ', fmt(x$z) )
  cat(stat, '\n', sep = '')

  perm <- x$permutation
  cat('Permutations: ', x$n_perm, if (!is.null(perm)) paste0(' (', perm, ')'), '\n', sep = '')

  extra <- test_details(x, digits)
  if (length(extra)) cat('\n', paste0(extra, '\n'), sep = '')

  invisible(x)
}



# One row of the core fields, with the same columns for every test, so that
# results from different tests can be rbind()-ed.
#' @export
as.data.frame.ecodive_test <- function (x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    test        = sub('^ecodive_', '', class(x)[[1]]),
    stat_name   = names(x$statistic),
    statistic   = unname(x$statistic),
    p_value     = x$p_value,
    z           = x$z,
    n_perm      = x$n_perm,
    row.names   = row.names,
    stringsAsFactors = FALSE )
}
