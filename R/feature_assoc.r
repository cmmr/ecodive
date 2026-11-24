# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Associations between the features of a count table: methods whose result
# is a `dist` over features (vlr(), rho_p(), phi_s(), sparcc()). Facts about
# each method live here, keyed by the result's `method` attribute, never as
# flags on the result (plans/README.md, "Result attributes and plots"):
#
#   name      Display name: the result's title (see ?ecodive_display).
#   id        The `method` attribute of its results.
#   alt_ids   Other names, space-separated, for match_feature_assoc().
#   kind      'similarity' (larger is closer: rho_p, sparcc) or
#             'dissimilarity' (vlr, phi_s).
#   lower,    The range of its values. Numerical approximations may stray
#   upper     outside it (sparcc() warns when they do).
#   fisher_z  TRUE for Pearson-type correlations, for which Fisher's z
#             transform is meaningful.
#   sparse    TRUE when zeros mean "no edge was selected" (none yet).
#   strict    TRUE when its own cutoff keeps |r| > threshold rather than >=
#             (none yet).
#
# A row here means "over features": dist_clust_order() clusters a
# similarity as 1 - x, and plot.ecodive_dist() labels the axes "features"
# and draws a similarity on a diverging scale.
FEATURE_ASSOCS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                     id,     alt_ids, kind,          lower, upper, fisher_z, sparse, strict
      Variance of Log Ratios,   vlr,    logratiovariance,  dissimilarity,     0,   Inf,    FALSE,  FALSE,  FALSE
      Proportionality (rho_p),  rho_p,  rho rhop propr,    similarity,       -1,     1,    FALSE,  FALSE,  FALSE
      Proportionality (phi_s),  phi_s,  phi phis phs,      dissimilarity,     0,   Inf,    FALSE,  FALSE,  FALSE
      SparCC Correlation,       sparcc, fastspar,          similarity,       -1,     1,     TRUE,  FALSE,  FALSE
  " )

  df$alt_ids <- as.character(df$alt_ids)
  return (df)
})


# The FEATURE_ASSOCS row (a list) of association method `id`, or NULL.
feature_assoc_info <- function (id) {
  if (!is.character(id) || length(id) != 1L) return (NULL)
  i <- match(id, FEATURE_ASSOCS$id)
  if (is.na(i)) return (NULL)
  as.list(FEATURE_ASSOCS[i, , drop = FALSE])
}


# TRUE if `d` is a similarity over features, by its `method` attribute.
is_similarity <- function (d) {
  info <- feature_assoc_info(dist_metric_of(d))
  !is.null(info) && identical(info$kind, 'similarity')
}



#' Associations Between Features
#'
#' Runs any of ecodive's measures of association between the features of a
#' count table by name, as [count_ord()] does for ordinations. Every one
#' returns a `dist` over the features.
#'
#' @param counts  A numeric matrix of counts (samples \eqn{\times}
#'        features), or any of the objects listed under **Input Types**.
#'
#' @param method  The name of an association measure: `"vlr"`, `"rho_p"`,
#'        `"phi_s"` or `"sparcc"`. Flexible matching is supported (see
#'        below). Programmatic access via `list_feature_assocs()`.
#'
#' @param ...  Further arguments for the method's function, such as
#'        `margin`, `pseudocount`, `clust` and `cpus`, [rho_p()]'s and
#'        [phi_s()]'s `denominator`, or [sparcc()]'s `n_draws`,
#'        `permutations` and `seed`. Each method keeps its own defaults.
#'        `seed` is ignored by the methods that have no random parts (all
#'        but [sparcc()]). [sparcc()] takes no `pseudocount`: it handles
#'        zeros itself.
#'
#' @return The method's result: a `dist` over the features, of class
#'   `c("ecodive_dist", "dist")`, with attribute `method` set to the
#'   method's ID. `feature_assoc(counts, 'sparcc', ...)` is the same as
#'   `sparcc(counts, ...)`.
#'
#' @details
#'
#' **List of Association Measures**
#'
#' | Option / Function Name | Measure | Kind | Range |
#' | :--------------------- | :------ | :--- | :---- |
#' | `vlr`    | Variance of log ratios                    | dissimilarity | \eqn{[0, \infty)} |
#' | `rho_p`  | Proportionality \eqn{\rho_p}              | similarity    | \eqn{[-1, 1]} |
#' | `phi_s`  | Symmetric proportionality \eqn{\phi_s}    | dissimilarity | \eqn{[0, \infty)} |
#' | `sparcc` | SparCC correlation                        | similarity    | \eqn{[-1, 1]} |
#'
#' The similarities ([rho_p()] and [sparcc()]) can be made into networks
#' with [assoc_network()], or listed with [assoc_pairs()].
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters.
#' Other names also work: `'rho'` or `'propr'` for `rho_p`; `'phi'` or
#' `'phs'` for `phi_s`; `'FastSpar'` for `sparcc`; `'log ratio variance'`
#' for `vlr`. `'proportionality'` matches both `rho_p` and `phi_s`, and so
#' is an error.
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#' @inheritSection documentation Input Types
#'
#' @seealso [vlr()], [rho_p()], [phi_s()], [sparcc()],
#'   [list_feature_assocs()]; [assoc_network()]; [cross_assoc()] for the
#'   associations between the features of two tables.
#'
#' @export
#' @examples
#'     r <- feature_assoc(ex_counts, 'rho_p', pseudocount = 1)
#'     r
#'
#'     # The same arguments for every log-ratio measure.
#'     lapply(c('vlr', 'rho_p', 'phi_s'), function (m)
#'       range(feature_assoc(ex_counts, m, pseudocount = 1)) )
#'
feature_assoc <- function (counts, method, ...) {

  method <- match_feature_assoc(method)

  args <- list(...)
  if (!'seed' %in% method$params) args$seed <- NULL

  # quote(counts), so that a large table isn't copied into the call.
  do.call(method$func, c(list(quote(counts)), args))
}



#' Find and Browse Association Measures
#'
#' Programmatic access to the list of association measures that
#' [feature_assoc()] runs, and their associated functions.
#'
#' @param method   The name of an association measure to search for.
#'        Supports partial matching. All non-alpha characters are ignored.
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
#' **`match_feature_assoc()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Measure name, e.g. `"SparCC Correlation"`
#' * `id` : Measure ID - also the name of the function, e.g. `"sparcc"`
#' * `kind` : `"similarity"` (larger means more associated) or
#'   `"dissimilarity"`.
#' * `lower`, `upper` : The range of its values.
#' * `fisher_z` : Whether its values are Pearson-type correlations, which
#'   [assoc_network_diff()] can compare with Fisher's \eqn{z}.
#' * `sparse` : Whether its zeros mean "no edge" (a selected graph).
#' * `strict` : Whether its own cutoff keeps \eqn{|r| >} the cutoff,
#'   rather than \eqn{\ge}.
#' * `func` : The function for this measure, e.g. `ecodive::sparcc`
#' * `params` : Formal args for `func`, e.g. `c("counts", "margin", ...)`
#'
#'
#' **`list_feature_assocs()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_feature_assoc()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of measure IDs.
#' * `val = "name"` : A character vector of measure names.
#'
#' If `nm` is set, then the names of the vector or list will be the measure
#' ID (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names
#' will be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname feature_assocs
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available measures.
#'     list_feature_assocs()
#'
#'     # Find the measure for 'propr'.
#'     m <- match_feature_assoc('propr')
#'     m$id
#'
list_feature_assocs <- function (
    val = c('data.frame', 'list', 'func', 'id', 'name'),
    nm  = c(NA, 'id', 'name') ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }

  df <- FEATURE_ASSOCS

  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) feature_assoc_entry(df[df$id == id,]))
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


#' @rdname feature_assocs
#' @export
#'
match_feature_assoc <- function (method) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')
  if (is.na(method))         stop('`method` must not be NA')

  df <- FEATURE_ASSOCS

  # Exact match
  if (method %in% df$id)
    return (feature_assoc_entry(df[df$id == method,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  match <- match_name(needle, 'FEATURE_ASSOCS', df$id)

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  feature_assoc_entry(df[df$id == match,])
}


# One row of FEATURE_ASSOCS as a list, with the method's function and its args.
feature_assoc_entry <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(info$id, ENV)
  info$params  <- names(formals(info$func))
  info
}
