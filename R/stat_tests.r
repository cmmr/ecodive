# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



# Statistical tests of a distance matrix against one variable. `groups` is
# the kinds of variable each method's `groups` takes: "categorical" (factor,
# character, or logical) and/or "numeric".
DIST_TESTS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                            id,             groups,              alt_ids
      Analysis of Similarities,                        anosim,         categorical,         analysisofsimilarity
      Mantel Test,                                     mantel,         numeric,             manteltest
      Partial Mantel Test,                             mantel_partial, numeric,             partialmantel mantelpartial
      Multi-Response Permutation Procedure,            mrpp,           categorical,         multiresponsepermutation
      Permutational Multivariate Analysis of Variance, permanova,      categorical numeric, adonis adonis2 npmanova permutationalmanova
      Homogeneity of Multivariate Dispersions,         permdisp,       categorical,         betadisper permdisp2 dispersion multivariatedispersion
  " )

  return (df)
})



#' Distance Test Wrapper Function
#'
#' Runs any of ecodive's statistical tests of a distance matrix by name, as
#' [dist_ord()] does for ordination methods: does one variable, such as
#' treatment groups, explain the differences between samples?
#'
#' @param dm  The distances to test: a `dist` object such as from
#'        [beta_div()], or a symmetric matrix of distances.
#'
#' @param test  The name of a test: `"anosim"`, `"mantel"`,
#'        `"mantel_partial"`, `"mrpp"`, `"permanova"`, or `"permdisp"`.
#'        Flexible matching is supported (see below). Programmatic access via
#'        `list_dist_tests()`.
#'
#' @param groups  The variable to test: a vector with one value per sample,
#'        matched to the samples in `dm` by name when both are named, else in
#'        the same order as `dm`. Or, when `data` is given, the name of one of
#'        its columns (or a vector with one value per row of `data`). Each
#'        method takes categorical groups, numeric values, or both: see
#'        `list_dist_tests()`.
#'
#' @param data  Optional sample metadata: a data frame with one row per
#'        sample, holding `groups` (and `strata`) as columns. Rows are matched
#'        to the samples in `dm` by row name. Default: `NULL`.
#'
#' @param ...  Further arguments for the method's function, such as
#'        `permutations`, `strata`, `seed`, or `cpus`, which every method
#'        takes, or a method's own, such as `weights` for [mrpp()] or
#'        `method` (the correlation) for [mantel()], or `given` for
#'        [mantel_partial()]. Each method keeps its own defaults. For
#'        `mantel` and `mantel_partial`, `groups` is passed as the function's
#'        `y`.
#'
#' @return The method's result: a list of class `"ecodive_test"` whose first
#'   fields are `statistic`, `p_value`, `z`, and `n_perm`.
#'   `as.data.frame()` turns it into one row with the same columns for every
#'   method, so results can be `rbind()`-ed. `dist_test(dm, 'mrpp',
#'   groups, data, ...)` is the same as `mrpp(dm, groups, data, ...)`.
#'
#' @details
#'
#' **List of Tests**
#'
#' | Option / Function Name | Test Name                                         | Groups                 |
#' | :--------------------- | :------------------------------------------------ | :--------------------- |
#' | `anosim`               | Analysis of Similarities (ANOSIM)                 | categorical            |
#' | `mantel`               | Mantel Test                                       | numeric                |
#' | `mantel_partial`       | Partial Mantel Test                               | numeric                |
#' | `mrpp`                 | Multi-Response Permutation Procedure (MRPP)       | categorical            |
#' | `permanova`            | Permutational Multivariate Analysis of Variance   | categorical or numeric |
#' | `permdisp`             | Homogeneity of Multivariate Dispersions (PERMDISP) | categorical            |
#'
#' Categorical groups are factors, character, or logical vectors; a method
#' that takes only categorical groups gives an error for numbers, and the
#' reverse. `mantel` correlates `dm` with the differences between samples'
#' values, and `mantel_partial` does so after removing the effect of
#' `given` (such as geographic distances); to compare two distance
#' matrices, call [mantel()] or [mantel_partial()] directly. To
#' test several variables at once, call [permanova()] with a `formula`,
#' either directly or as `dist_test(dm, 'permanova', data = meta,
#' formula = ~ a + b)`.
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters. E.g.
#' `test = 'PERMANOVA'` or `'perma'` selects `permanova`. Common other names
#' also work: `'adonis'` or `'adonis2'` for `permanova`; `'betadisper'` or
#' `'PERMDISP2'` for `permdisp`; `'Analysis of Similarities'` for
#' `anosim`; and `'mantel.partial'` or `'partial Mantel'` for
#' `mantel_partial`.
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::anosim()`, `vegan::mantel()`, and `vegan::mrpp()` share names
#' with ecodive's functions, and whichever package was attached last masks
#' the other's.
#' Write `ecodive::anosim()` and so on, or use `dist_test()`, which always
#' runs ecodive's, to be sure.
#'
#' @seealso [anosim()], [mantel()], [mantel_partial()], [mrpp()],
#'   [permanova()], [permdisp()], [list_dist_tests()]
#'
#' @export
#' @examples
#'     # 30 samples in three treatment groups. The 'high' group has more
#'     # of the first five genera.
#'     set.seed(1)
#'     counts <- matrix(rpois(30 * 20, 5), 30, 20)
#'     rownames(counts) <- paste0('S', 1:30)
#'     meta <- data.frame(
#'       row.names = rownames(counts),
#'       treatment = rep(c('ctrl', 'low', 'high'), each = 10),
#'       ph        = round(rnorm(30, 7, 0.3), 1) )
#'     high <- meta$treatment == 'high'
#'     counts[high, 1:5] <- counts[high, 1:5] + 5
#'
#'     dm <- bray(counts)
#'
#'     dist_test(dm, 'permanova', 'treatment', meta)
#'
#'     # Every test of treatment, as one table.
#'     tests <- list_dist_tests(val = 'id', groups = 'categorical')
#'     do.call(rbind, lapply(tests, function (m)
#'       as.data.frame(dist_test(dm, m, 'treatment', meta, cpus = 1)) ))
#'
#'     # A numeric variable.
#'     dist_test(dm, 'adonis2', 'ph', meta)
#'     dist_test(dm, 'mantel', 'ph', meta, method = 'spearman')
#'
dist_test <- function (dm, test, groups, data = NULL, ...) {

  test <- match_dist_test(test)
  args <- list(...)

  # No `groups`: e.g. permanova()'s `formula`, given in `...`.
  if (missing(groups))
    return (do.call(test$func, c(list(quote(dm)), list(data = data), args)))

  test_group_type(groups, data, test)

  # quote(), so that a large `dm` isn't copied into the call.
  y <- if (test$id %in% c('mantel', 'mantel_partial')) 'y' else 'groups'
  do.call(test$func, c(list(quote(dm)), stats::setNames(list(quote(groups)), y), list(data = data), args))
}


# Checks that `groups` (a vector, or a column of `data`) is a kind of
# variable that `test` (from match_dist_test()) tests. Anything else is
# left to the test's function to check.
test_group_type <- function (groups, data, test) {

  vals <- NULL
  what <- '`groups`'

  if (is.null(data)) {
    if (is.atomic(groups) && !(is.character(groups) && length(groups) == 1))
      vals <- groups
  }
  else if (is.character(groups) && length(groups) == 1 && !is.na(groups) &&
           groups %in% colnames(data)) {
    vals <- if (is.data.frame(data)) data[[groups]] else data[, groups]
    what <- paste0('`data` column ', name_list(groups))
  }
  else if (is.atomic(groups) && length(groups) == NROW(data)) {
    vals <- groups
  }

  if (is.null(vals)) return (invisible())

  type <- if (is.factor(vals) || is.character(vals) || is.logical(vals)) 'categorical'
          else if (is.numeric(vals)) 'numeric'
          else class(vals)[[1]]

  ok <- strsplit(test$groups, ' ', fixed = TRUE)[[1]]
  if (type %in% ok) return (invisible())

  kinds <- c(
    categorical = 'categorical groups (a factor, character, or logical vector)',
    numeric     = 'numeric values' )

  others <- list_dist_tests(val = 'id', groups = intersect(type, names(kinds)))

  stop(
    '`test = "', test$id, '"` tests ', paste(kinds[ok], collapse = ' or '),
    ', but ', what, ' is ', if (type %in% names(kinds)) type else paste0('of class ', type), '.',
    if (length(others)) paste0(' For ', type, ' values, use: ', paste(others, collapse = ', '), '.'),
    call. = FALSE )
}



#' Find and Browse Distance Tests
#'
#' Programmatic access to the list of available statistical tests of a
#' distance matrix, and their associated functions.
#'
#' @param test   The name of a test to search for. Supports partial
#'        matching. All non-alpha characters are ignored.
#'
#' @param val   Sets the return value for this function call. See "Value"
#'        section below. Default: `"data.frame"`
#'
#' @param nm   What value to use for the names of the returned object.
#'        Default is `"id"` when `val` is `"list"` or `"func"`, otherwise the
#'        default is `NA` (no name).
#'
#' @param groups   Filter by the kind of variable the test takes:
#'        `"categorical"` or `"numeric"`. Default: `NULL` (no filtering).
#'
#' @return
#'
#' **`match_dist_test()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Test name, e.g. `"Multi-Response Permutation Procedure"`
#' * `id` : Test ID - also the name of the function, e.g. `"mrpp"`
#' * `groups` : The kinds of variable the test takes, e.g.
#'   `"categorical"`, or `"categorical numeric"`.
#' * `func` : The function for this test, e.g. `ecodive::mrpp`
#' * `params` : Formal args for `func`, e.g. `c("dm", "groups", "data", ...)`
#'
#'
#' **`list_dist_tests()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_dist_test()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of test IDs.
#' * `val = "name"` : A character vector of test names.
#' * `val = "groups"` : A character vector of the kinds of variable each test takes.
#'
#' If `nm` is set, then the names of the vector or list will be the test ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname dist_tests
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available tests.
#'     list_dist_tests()
#'
#'     # The tests of a numeric variable.
#'     list_dist_tests(val = 'id', groups = 'numeric')
#'
#'     # Find the test for 'adonis2'.
#'     m <- match_dist_test('adonis2')
#'     m$id
#'
list_dist_tests <- function (
    val    = c('data.frame', 'list', 'func', 'id', 'name', 'groups'),
    nm     = c(NA, 'id', 'name'),
    groups = NULL ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }


  df <- DIST_TESTS

  if (!is.null(groups)) {

    kinds <- c('categorical', 'numeric')
    if (length(bad <- setdiff(groups, kinds)) > 0)
      stop(
        "Invalid value for `groups`: ", paste(collapse = ", ", bad), "\n",
        "Options are: NULL, ", paste(collapse = ", ", kinds), "." )

    has <- vapply(strsplit(df$groups, ' ', fixed = TRUE), function (g) any(groups %in% g), logical(1))
    df  <- df[has,, drop = FALSE]
  }


  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) dist_test_info(df[df$id == id,]))
    if (!is.na(nm)) names(result) <- df[[nm]]
    return (result)
  }

  if (val == 'func') {
    result <- mget(df$id, ENV)
    if (!is.na(nm)) names(result) <- df[[nm]]
    return (result)
  }

  result <- df[[val]]
  if (!is.na(nm)) names(result) <- df[[nm]]
  return (result)
}


#' @rdname dist_tests
#' @export
#'
match_dist_test <- function (test, groups = NULL) {

  if (!is.character(test)) stop('`test` must be a character vector')
  if (length(test) != 1)   stop('`test` must be length 1')
  if (is.na(test))         stop('`test` must not be NA')

  df <- list_dist_tests(groups = groups)

  # Exact match
  if (test %in% df$id)
    return (dist_test_info(df[df$id == test,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(test))
  if (nchar(needle) == 0) stop('`test` must contain alphabetic characters')

  match <- match_name(needle, 'DIST_TESTS', df$id)

  if (length(match) > 1) stop(
    '`test = "', test, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`test = "', test, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  dist_test_info(df[df$id == match,])
}


# One row of DIST_TESTS as a list, with the method's function and its args.
dist_test_info <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(info$id, ENV)
  info$params  <- names(formals(info$func))
  info
}
