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
# permdisp(), dist_test()), and `strata`, by per_sample()'s one rule:
#   * with `data`: the name of one of its columns, or a vector with one value
#     per row of `data`;
#   * without `data`: one value per sample, matched to the dist's `labels`
#     by name when both are named, else in order.
# Every test gives the same messages.
#
# Returns list(data, col, strata):
#   data    A data frame with one row per sample, in the dist's order: the
#           columns of `data` (if given), and the groups in column `col`.
#           Its row names are the samples' (when they have unique names).
#   strata  The strata values in the same order, or NULL.
# Pass `data` and `strata` on to model_design(); test_group_values() gives
# the groups and strata codes directly.
#------------------------------------------------------------------------------
test_groups <- function (groups, data, n, labels, strata = NULL) {

  if (is.null(groups))
    stop('`groups` is required: one value per sample, or a column name of `data`.', call. = FALSE)

  if (!is.atomic(groups) || !is.null(dim(groups)))
    stop(
      '`groups` must be a vector or factor with one value per sample, ',
      'or the name of a column of `data`.', call. = FALSE )

  idx <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data))
      stop('`data` must be a data.frame.', call. = FALSE)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    dn   <- if (.row_names_info(data) > 0) rownames(data) else NULL
    idx  <- match_rows(labels, n, dn, nrow(data))
  }

  val <- per_sample(groups, 'groups', data, idx, labels, n)
  src <- if (!is.null(data) && is.character(groups) && length(groups) == 1 &&
             groups %in% names(data)) groups

  # Values for samples not in `dm` were left out, even NAs.
  if (anyNA(val))
    stop(
      if (is.null(src)) '`groups` must not contain NAs.'
      else paste0('Missing values (NA) in `data` column ', name_list(src), '.'),
      call. = FALSE )

  # A factor's NA level (from addNA()) is a group like any other, as in
  # model.matrix(); anyNA() is FALSE for it. Name it "NA", so that it can
  # label results.
  if (is.factor(val) && anyNA(levels(val))) {
    lev <- levels(val)
    lev[is.na(lev)] <- make.unique(c(lev[!is.na(lev)], 'NA'))[[sum(!is.na(lev)) + 1L]]
    levels(val) <- lev
  }

  strata <- per_sample(strata, 'strata', data, idx, labels, n)

  out <- if (is.null(data)) data.frame(row.names = seq_len(n)) else data[idx, , drop = FALSE]
  col <- if (!is.null(src)) src else make.unique(c(names(out), 'group'))[[ncol(out) + 1L]]
  out[[col]] <- val
  rownames(out) <- if (!is.null(labels) && !anyNA(labels) && !anyDuplicated(labels)) labels else NULL

  list(data = out, col = col, strata = strata)
}


#------------------------------------------------------------------------------
# test_groups(), then each sample's group and stratum, in the order of the
# dist's samples. Returns list(groups (the values, with unused factor levels
# dropped), strata (0-based codes, or NULL), col).
#------------------------------------------------------------------------------
test_group_values <- function (groups, data, n, labels, strata = NULL) {

  g   <- test_groups(groups, data, n, labels, strata)
  val <- g$data[[g$col]]
  if (is.factor(val)) val <- droplevels(val)

  list(groups = val, strata = strata_codes(g$strata), col = g$col)
}


# 0-based integer codes of strata values (already in the samples' order), or
# NULL. A factor's NA level is a stratum like any other.
strata_codes <- function (strata) {
  if (is.null(strata)) return (NULL)
  if (anyNA(strata)) stop('`strata` must not contain NAs.', call. = FALSE)
  as.integer(factor(strata, exclude = NULL)) - 1L
}



#------------------------------------------------------------------------------
# Permutation p-values: the proportion of the permutations, counting the
# observed one, with a statistic at least as extreme as the observed one,
# (count + 1) / (n_perm + 1). As vegan, values within `tol` of the observed
# one count as ties; give `tol` on the statistic's scale when it has units.
# Permuted values that are NA or NaN never count.
#
#   stat         The observed statistics (length m).
#   perm         The permuted ones: a vector (m = 1) or n_perm x m matrix.
#   alternative  'greater' (larger is more extreme), 'less', or 'two.sided'
#                (larger in absolute value).
#------------------------------------------------------------------------------
perm_p <- function (
    stat, perm, alternative = 'greater', tol = sqrt(.Machine$double.eps) ) {

  perm   <- as.matrix(perm)
  n_perm <- nrow(perm)
  if (n_perm == 0) return (rep(NA_real_, length(stat)))

  lim <- rep(switch(alternative, less = stat + tol, two.sided = abs(stat) - tol, stat - tol), each = n_perm)
  hit <- switch(alternative,
    greater   = perm >= lim,
    less      = perm <= lim,
    two.sided = abs(perm) >= lim )

  p <- (colSums(hit, na.rm = TRUE) + 1) / (n_perm + 1)
  p[is.na(stat)] <- NA_real_
  unname(p)
}



#------------------------------------------------------------------------------
# z-scores of statistics against their permuted values, as
# summary(vegan::permustats()): (stat - mean(perm)) / sd(perm). `perm` is a
# vector (one statistic) or an n_perm x m matrix. NA with fewer than 2
# permutations, or when the permuted values are all equal (an sd that is
# only rounding noise would give a huge, meaningless z).
#------------------------------------------------------------------------------
perm_z <- function (stat, perm) {

  perm <- as.matrix(perm)
  if (nrow(perm) < 2) return (rep(NA_real_, length(stat)))

  vapply(seq_along(stat), function (j) {
    p <- perm[, j]
    s <- stats::sd(p)
    if (is.na(stat[[j]]) || !is.finite(s) || s <= 1e-12 * max(abs(p))) return (NA_real_)
    z <- (stat[[j]] - mean(p)) / s
    if (is.finite(z)) z else NA_real_
  }, numeric(1))
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
