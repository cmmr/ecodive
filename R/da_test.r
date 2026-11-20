# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



# Differential-abundance methods. Each is a function of a count table that
# returns a per-feature data frame whose `method` attribute is its id.
#
#   design        'groups': f(counts, groups, data = NULL, ...), one
#                 categorical variable; 'formula': f(counts, formula = NULL,
#                 data = NULL, groups = NULL, ...), a model, with `groups`
#                 its term of interest (see da_test()).
#   effect        The result's column of effect sizes (NA: none).
#   base          The base of its logs: 'e', '2', or 'norm' (natural logs
#                 for a log-ratio `norm`, none otherwise); NA: no effect.
#   zero_free_ok  Whether the input may be a zero-free table, such as
#                 replace_zeros()'s (FALSE: the method models counts).
#
# `effect`, `base` and `zero_free_ok` are for da_sensitivity().
DA_TESTS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                    id,            design,  effect,      base, zero_free_ok, alt_ids
      Analysis of Composition of Microbiomes,  ancom,         groups,  ,            ,     TRUE,         analysisofcompositionsofmicrobiomes
      ANCOM-BC,                                ancombc,       formula, log_fc,      e,    TRUE,         ancombiascorrected analysisofcompositionsofmicrobiomeswithbiascorrection
      ANCOM-BC2,                               ancombc2,      formula, log_fc,      e,    TRUE,
      Dirichlet-Multinomial t-test,            dirmult_ttest, groups,  log2_fc,     2,    FALSE,        dirichletmultinomialwelchttest aldex aldex2
      Wilcoxon Rank-Sum Test,                  da_wilcoxon,   groups,  median_diff, norm, TRUE,         wilcoxon wilcox wilcoxontest mannwhitney mannwhitneyu ranksum
      Welch t-test,                            da_ttest,      groups,  diff,        norm, TRUE,         ttest welch studentttest
      Linear Model,                            da_lm,         formula, estimate,    norm, TRUE,         lm linearregression ols
      LinDA,                                   linda,         formula, log2_fc,     2,    TRUE,
  " )

  df$zero_free_ok <- as.logical(df$zero_free_ok)
  df$base         <- as.character(df$base)
  return (df)
})



#' Differential Abundance Testing
#'
#' Runs any of ecodive's differential-abundance (DA) methods by name, as
#' [count_ord()] does for ordinations, with one convention for the groups
#' or model to compare.
#'
#' @param counts  A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} features), or any of the objects listed under **Input
#'        Types**.
#'
#' @param method   The name of a DA method: `"ancom"`, `"ancombc"`,
#'        `"ancombc2"`, `"da_lm"`, `"da_ttest"`, `"da_wilcoxon"`,
#'        `"dirmult_ttest"`, or `"linda"`. Flexible matching is supported
#'        (see below). Programmatic access via `list_da_tests()`.
#'
#' @param groups  The grouping to compare: the name of a column of `data`,
#'        or a vector with one value per sample. For a model method, the
#'        term of interest: alone, it is the model (`formula = ~ groups`);
#'        with `formula`, it names one of its terms (see below).
#'
#' @param data  A data frame of sample metadata, with a row per sample,
#'        matched to the samples of `counts` by name when both have names,
#'        else by position.
#'
#' @param formula  For a model method (`design = "formula"` in
#'        `list_da_tests()`), the model: a one-sided formula over the columns
#'        of `data`, such as `~ status + age`. A group method takes a
#'        formula of one variable, `~ status`, as `groups = "status"`.
#'
#' @param ...   Further arguments for the method's function, such as
#'        `margin`, `alpha`, `p_adjust`, `cpus`, `contrast` ([da_ttest()],
#'        [da_wilcoxon()], [dirmult_ttest()]), `norm` and `pseudocount`. Each
#'        method keeps its own defaults; see [ancom()], [ancombc()],
#'        [ancombc2()], [da_lm()], [da_ttest()], [da_wilcoxon()],
#'        [dirmult_ttest()], and [linda()]. `seed` is ignored by the methods
#'        that have no random parts (all but [dirmult_ttest()]).
#'
#' @return The method's result, unchanged: `da_test(counts, 'ancombc2',
#'   ...)` is identical to `ancombc2(counts, ...)`, attributes included.
#'   Every result is a data frame with a row per feature (and covariate, for
#'   the model methods) and a `method` attribute, the method's id.
#'
#' @details
#'
#' **List of DA Methods**
#'
#' | Option / Function Name | Method Name                            | Design  |
#' | :--------------------- | :------------------------------------- | :------ |
#' | `ancom`                | Analysis of Composition of Microbiomes (ANCOM) | groups |
#' | `ancombc`              | ANCOM with bias correction (ANCOM-BC)  | formula |
#' | `ancombc2`             | ANCOM-BC2                              | formula |
#' | `da_lm`                | Linear model per feature               | formula |
#' | `da_ttest`             | Welch t-test per feature               | groups  |
#' | `da_wilcoxon`          | Wilcoxon rank-sum test per feature     | groups  |
#' | `dirmult_ttest`        | Dirichlet-multinomial t-test (ALDEx2-style) | groups |
#' | `linda`                | Linear models for differential abundance (LinDA) | formula |
#'
#' **Groups and models.** *Group methods* compare the samples of one
#' categorical variable: `groups`, a column name of `data` or one value per
#' sample. *Model methods* fit a linear model, `formula`, of any number of
#' variables; `groups` names its term of interest (the one their plots and
#' post-hoc tests use by default), and `groups` alone is the same as
#' `formula = ~ groups`. So `da_test(counts, m, 'status', meta)` compares
#' the samples by `status` for every method `m`, and covariates go in
#' `formula`: `da_test(counts, 'da_lm', 'status', meta, ~ status + age)`.
#' A one-variable `formula` works for the group methods too; more than one
#' variable is an error for them.
#'
#' The methods' assumptions differ, and so do their results: see each
#' method's page. In benchmarks on real data (Nearing et al. 2022; Wirbel
#' et al. 2024; Pelto et al. 2025), the plain per-feature tests
#' ([da_wilcoxon()], [da_ttest()], [da_lm()]) controlled false discoveries
#' well, and the methods often disagree: a feature found by several methods
#' is a firmer finding than one found by one.
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters
#' (digits count only to tell `'ANCOM-BC2'` from `'ANCOM-BC'`). E.g.
#' `method = 'Wilcoxon'`, `'Mann-Whitney'` or `'rank sum'` selects
#' `da_wilcoxon`; `'t-test'` or `'Welch'` selects `da_ttest`; `'lm'` or
#' `'linear model'` selects `da_lm`; and `'ALDEx2'` selects
#' `dirmult_ttest`, an ALDEx2-style method (not ALDEx2 itself, whose
#' results agree with it in direction and rough size only).
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#'
#' **Same names in other packages**
#'
#' None known.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). Microbiome
#' differential abundance methods produce different results across 38
#' datasets. *Nature Communications*, 13, 342.
#' \doi{10.1038/s41467-022-28034-z}
#'
#' Wirbel, J., Essex, M., Forslund, S. K., & Zeller, G. (2024). A realistic
#' benchmark for differential abundance testing and confounder adjustment in
#' human microbiome studies. *Genome Biology*, 25, 247.
#' \doi{10.1186/s13059-024-03390-9}
#'
#' Pelto, J., Auranen, K., Kujala, J. V., & Lahti, L. (2025). Elementary
#' methods provide more replicable results in microbial differential
#' abundance analysis. *Briefings in Bioinformatics*, 26, bbaf130.
#' \doi{10.1093/bib/bbaf130}
#'
#' @seealso [ancom()], [ancombc()], [da_lm()], [da_ttest()],
#'   [da_wilcoxon()], [dirmult_ttest()], [linda()], [list_da_tests()];
#'   [replace_zeros()] for the zeros; [da_sensitivity()] to rerun a method
#'   over pseudocounts, zero replacements and scale uncertainties.
#'
#' @export
#' @examples
#'     grp <- c('oral', 'oral', 'other', 'other')
#'     da_test(ex_counts, 'da_wilcoxon', grp, pseudocount = 1)
#'
#'     # The same comparison by every method, with a model method's covariate.
#'     meta <- data.frame(site = grp, depth = rowSums(ex_counts))
#'     da_test(ex_counts + 1, 'da_lm', 'site', meta, ~ site + depth)
#'
#'     # Which features each method finds.
#'     x <- rbind(
#'       c(12, 11, 10, 10, 10, 10, 10), c( 9, 11, 12, 10, 10, 10, 10),
#'       c( 1, 11, 10, 11, 10,  5,  9), c(22, 21,  9, 10, 10, 10, 10),
#'       c(20, 22, 10, 10, 13, 10, 10), c(23, 21, 14, 10, 10, 10, 10) )
#'     colnames(x) <- paste0('b', 1:7)
#'     meta <- data.frame(trt = rep(c('treatment', 'placebo'), each = 3))
#'     ids  <- c('ancom', 'da_ttest', 'da_wilcoxon', 'da_lm')
#'     sapply(ids, function (m) {
#'       res <- da_test(x, m, 'trt', meta)
#'       if ('covariate' %in% names(res)) res <- res[res$covariate != '(Intercept)', ]
#'       res$reject })
#'
da_test <- function (counts, method, groups = NULL, data = NULL, formula = NULL, ...) {

  method <- match_da_test(method)

  args <- list(...)
  if (!'seed' %in% method$params) args$seed <- NULL

  if (method$design == 'groups') {
    if (!is.null(formula)) {
      if (!is.null(groups))
        stop(
          '`', method$id, '` compares one grouping: give `groups` or a ',
          'formula of one variable, not both.', call. = FALSE )
      groups <- da_formula_groups(formula, method$id)
    }
    args <- c(list(quote(counts), groups = groups, data = data), args)
  }
  else {
    args <- c(list(quote(counts), formula = formula, data = data, groups = groups), args)
  }

  # quote(counts), so that a large table isn't copied into the call.
  do.call(method$func, args)
}



# The variable of a one-variable formula, such as ~ status, for a group
# method.
da_formula_groups <- function (formula, id) {

  if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
    txt <- if (grepl('~', formula, fixed = TRUE)) formula else paste('~', formula)
    formula <- tryCatch(stats::as.formula(txt), error = function (e) NULL)
  }
  if (!inherits(formula, 'formula'))
    stop('`formula` must be a formula, such as `~ status`.', call. = FALSE)

  tt  <- stats::terms(formula)
  lab <- attr(tt, 'term.labels')
  if (length(formula) == 3L) lab <- c(all.vars(formula[[2L]]), lab)

  if (length(lab) != 1L)
    stop(
      '`', id, '` compares one grouping; use a model method such as ',
      '`da_lm` or `ancombc2` for `formula` (', length(lab), ' terms).', call. = FALSE )
  if (!identical(all.vars(tt), lab))
    stop(
      '`formula` for `', id, '` must be a single variable of `data`, such ',
      'as `~ status`.', call. = FALSE )

  lab
}



#' Find and Browse Differential Abundance Methods
#'
#' Programmatic access to the list of differential-abundance methods that
#' [da_test()] runs, and their associated functions.
#'
#' @param method   The name of a DA method to search for. Supports partial
#'        matching. All non-alpha characters are ignored.
#'
#' @param val   Sets the return value for this function call. See "Value"
#'        section below. Default: `"data.frame"`
#'
#' @param nm   What value to use for the names of the returned object.
#'        Default is `"id"` when `val` is `"list"` or `"func"`, otherwise the
#'        default is `NA` (no name).
#'
#' @return
#'
#' **`match_da_test()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Method name, e.g. `"Welch t-test"`
#' * `id` : Method ID - also the name of the function, e.g. `"da_ttest"`
#' * `design` : `"groups"` (one grouping, given as `groups`) or `"formula"`
#'   (a model, given as `formula` and/or `groups`)
#' * `effect` : The result's column of effect sizes, e.g. `"log2_fc"`, or
#'   `NA` (ANCOM gives none)
#' * `base` : The base of the effect sizes' logs: `"e"`, `"2"`, or `"norm"`
#'   (natural logs for a log-ratio `norm`, no logs otherwise); `NA` without
#'   effect sizes
#' * `zero_free_ok` : Whether the method may be given a table with its zeros
#'   replaced ([replace_zeros()]): `FALSE` for the methods that model counts
#' * `func` : The function for this method, e.g. `ecodive::da_ttest`
#' * `params` : Formal args for `func`, e.g. `c("counts", "groups", ...)`
#'
#'
#' **`list_da_tests()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_da_test()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of method IDs.
#' * `val = "name"` : A character vector of method names.
#'
#' If `nm` is set, then the names of the vector or list will be the method ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname da_tests
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available methods.
#'     list_da_tests()
#'
#'     # Find the method for 'Mann-Whitney'.
#'     m <- match_da_test('Mann-Whitney')
#'     m$id
#'
list_da_tests <- function (
    val = c('data.frame', 'list', 'func', 'id', 'name'),
    nm  = c(NA, 'id', 'name') ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }

  df <- DA_TESTS

  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) da_test_info(df[df$id == id,]))
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


#' @rdname da_tests
#' @export
#'
match_da_test <- function (method) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')
  if (is.na(method))         stop('`method` must not be NA')

  df <- DA_TESTS

  # Exact match
  if (method %in% df$id)
    return (da_test_info(df[df$id == method,]))

  # An exact name with its digits, which the letters-only matching below
  # drops: 'ANCOM-BC2' is ancombc2, 'ANCOM-BC' is ancombc.
  digits <- gsub('[^a-z0-9]', '', tolower(method))
  keys   <- gsub('[^a-z0-9]', '', tolower(c(df$id, df$name)))
  hit    <- unique(c(df$id, df$id)[keys == digits])
  if (length(hit) == 1L)
    return (da_test_info(df[df$id == hit,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  # Methods of the other wrappers.
  if (needle %in% names(exact_names('DIST_TESTS'))) {
    id <- exact_names('DIST_TESTS')[[needle]]
    stop(
      '`method = "', method, '"`: ', id, '() is a distance test, not a ',
      'differential abundance test: use dist_test().', call. = FALSE )
  }
  for (reg in c('COUNT_ORDS', 'ORDINATIONS'))
    if (needle %in% names(exact_names(reg)))
      stop(
        '`method = "', method, '"` is an ordination, not a differential ',
        'abundance test: use ', if (reg == 'COUNT_ORDS') 'count_ord()' else 'dist_ord()',
        '.', call. = FALSE )

  match <- match_name(needle, 'DA_TESTS', df$id)

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  da_test_info(df[df$id == match,])
}


# One row of DA_TESTS as a list, with the method's function and its args.
da_test_info <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(info$id, ENV)
  info$params  <- names(formals(info$func))
  info
}


# A DA method's display name, from its id (the result's `method`
# attribute), for plots: "Welch t-test". Unknown ids are returned as given.
da_method_name <- function (id) {
  id <- paste0('', id)
  i  <- match(id, DA_TESTS$id)
  if (is.na(i)) id else DA_TESTS$name[[i]]
}
