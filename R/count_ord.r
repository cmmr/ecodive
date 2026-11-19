# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



# Ordination methods that take a count table (`counts`) rather than
# distances. Methods that take distances are in ORDINATIONS (dist_ord.r).
COUNT_ORDS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                id,   alt_ids
      Correspondence Analysis,             ca,   correspondenceanalysis
      Canonical Correspondence Analysis,   cca,  canonicalcorrespondenceanalysis constrainedcorrespondenceanalysis
      Linear Discriminant Analysis,        lda,  lineardiscriminantanalysis
      Orthogonal Partial Least Squares Discriminant Analysis, oplsda, opls orthogonalplsda orthogonalpartialleastsquaresda orthogonalpartialleastsquaresdiscriminantanalysis
      Principal Component Analysis,        pca,  principalcomponents prcomp
      Partial Least Squares Discriminant Analysis, plsda, partialleastsquaresda partialleastsquaresdiscriminantanalysis plsdiscriminantanalysis
      Redundancy Analysis,                 rda,  redundancyanalysis tbrda transformationbasedrda
      Robust Aitchison PCA,                rpca, robustpca deicode gemelli
      Phylogenetic Robust Aitchison PCA,   rpca_phylo, phylogeneticrpca phylorpca phylogeneticdeicode
  " )

  return (df)
})



#' Count Table Ordination
#'
#' Runs any of ecodive's ordination methods that take a count table by name,
#' as [dist_ord()] does for the methods that take distances.
#'
#' @inheritParams documentation
#'
#' @param method   The name of an ordination method: `"ca"`, `"cca"`, `"lda"`,
#'        `"oplsda"`, `"pca"`, `"plsda"`, `"rda"`, `"rpca"`, or `"rpca_phylo"`.
#'        Flexible matching is supported (see below). Programmatic access via
#'        `list_count_ords()`.
#'
#' @param ...   Further arguments for the method's function, such as `n_axes`,
#'        `margin`, `cpus`, [pca()]'s `norm`, `pseudocount`, and `scale`, or
#'        [rda()]'s and [cca()]'s `groups`, `data`, and `formula`, or
#'        [lda()]'s, [plsda()]'s and [oplsda()]'s `groups` and `data`, or
#'        [rpca_phylo()]'s `tree`. Each
#'        method keeps its own defaults; see [ca()], [cca()], [lda()],
#'        [oplsda()], [pca()], [plsda()], [rda()], [rpca()], and [rpca_phylo()].
#'        `seed` is ignored by the methods that have no random parts (all
#'        but [cca()], whose adjusted \eqn{R^2} is found by permutation), so
#'        the same arguments can be given to every method.
#'
#' @return The method's result, as described in "Ordination result" below.
#'   `count_ord(counts, 'pca', ...)` is the same as `pca(counts, ...)`.
#'
#' @details
#'
#' **List of Ordination Methods**
#'
#' | Option / Function Name | Method Name                                |
#' | :--------------------- | :----------------------------------------- |
#' | `ca`                   | Correspondence Analysis (CA)               |
#' | `cca`                  | Canonical Correspondence Analysis (CCA), constrained |
#' | `lda`                  | Linear Discriminant Analysis (LDA), supervised |
#' | `oplsda`               | Orthogonal PLS Discriminant Analysis (OPLS-DA), supervised |
#' | `pca`                  | Principal Component Analysis (PCA)         |
#' | `plsda`                | Partial Least Squares Discriminant Analysis (PLS-DA), supervised |
#' | `rda`                  | Redundancy Analysis (RDA), constrained     |
#' | `rpca`                 | Robust Aitchison PCA                       |
#' | `rpca_phylo`           | Phylogenetic Robust Aitchison PCA (needs `tree`) |
#'
#' `rda` is constrained by sample metadata, given as further arguments:
#' `count_ord(counts, 'rda', 'treatment', meta)` is
#' `rda(counts, 'treatment', meta)`. Without them, it is a PCA. `cca` is
#' constrained the same way; without constraints, it is a correspondence
#' analysis. `lda`, `plsda` and `oplsda` need the groups to separate, given
#' the same way: `count_ord(counts, 'plsda', 'treatment', meta)` is
#' `plsda(counts, 'treatment', meta)`. `oplsda` takes `n_orthogonal` rather
#' than `n_axes`. `rpca_phylo` needs a phylogenetic tree:
#' `count_ord(counts, 'rpca_phylo', tree = tree)`.
#'
#' Methods that ordinate distances, such as [pcoa()] and [nmds()], are run
#' by [dist_ord()].
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters. E.g.
#' `method = 'PCA'` or `'principal comp'` selects `pca`. Common other names
#' also work: `'prcomp'` for `pca`; `'redundancy analysis'`, `'tbRDA'`, or
#' `'transformation-based RDA'` for `rda`; `'canonical correspondence
#' analysis'` or `'constrained correspondence analysis'` for `cca`;
#' `'correspondence analysis'` for `ca`; and `'robust PCA'`, `'DEICODE'`,
#' or `'gemelli'` for `rpca`; `'phylogenetic RPCA'`, `'phylo-RPCA'`, or
#' `'phylogenetic DEICODE'` for `rpca_phylo`; and `'PLS-DA'`, `'PLS'`, or `'partial least
#' squares discriminant analysis'` for `plsda`; and `'LDA'` or `'linear
#' discriminant analysis'` for `lda`; and `'OPLS-DA'`, `'OPLS'`, or
#' `'orthogonal PLS-DA'` for `oplsda`. `'CA'` selects `ca`, and
#' `'CCA'` selects `cca`.
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#'
#' **Same names in other packages**
#'
#' `vegan::pca()`, `vegan::rda()`, `vegan::cca()` and `vegan::ca()` (and
#' `ca::ca()`) share their names with ecodive's [pca()], [rda()], [cca()]
#' and [ca()], as `mixOmics::plsda()` and `caret::plsda()` do with
#' [plsda()] and `MASS::lda()` with [lda()], and whichever package was attached last masks the other's. Write `ecodive::pca()` and so on, or use
#' `count_ord()`, which
#' always runs ecodive's, to be sure.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [ca()], [cca()], [lda()], [oplsda()], [pca()], [plsda()], [rda()], [rpca()], [rpca_phylo()], [list_count_ords()]; [dist_ord()] for
#'   the methods that take distances.
#'
#' @export
#' @examples
#'     ord <- count_ord(ex_counts, 'pca', norm = 'clr', pseudocount = 1)
#'     head(ord)
#'     str(ord)
#'
#'     # The same arguments for every method (but lda, plsda and oplsda,
#'     # which need groups, and rpca_phylo, which needs a tree).
#'     ids <- setdiff(list_count_ords('id'), c('lda', 'oplsda', 'plsda', 'rpca_phylo'))
#'     lapply(ids, function (m) attr(count_ord(ex_counts, m, n_axes = 2), 'prop'))
#'
#'     ord <- count_ord(ex_counts, 'rpca_phylo', tree = ex_tree, n_axes = 2)
#'
#'     count_ord(ex_counts, 'plsda', c('a', 'a', 'b', 'b'), n_axes = 1)
#'
count_ord <- function (counts, method, ...) {

  method <- match_count_ord(method)

  args <- list(...)
  if (!'seed' %in% method$params) args$seed <- NULL

  # quote(counts), so that a large table isn't copied into the call.
  do.call(method$func, c(list(quote(counts)), args))
}



#' Find and Browse Count Table Ordination Methods
#'
#' Programmatic access to the list of ordination methods that [count_ord()]
#' runs, and their associated functions.
#'
#' @param method   The name of an ordination method to search for. Supports
#'        partial matching. All non-alpha characters are ignored.
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
#' **`match_count_ord()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Method name, e.g. `"Principal Component Analysis"`
#' * `id` : Method ID - also the name of the function, e.g. `"pca"`
#' * `func` : The function for this method, e.g. `ecodive::pca`
#' * `params` : Formal args for `func`, e.g. `c("counts", "margin", ...)`
#'
#'
#' **`list_count_ords()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_count_ord()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of method IDs.
#' * `val = "name"` : A character vector of method names.
#'
#' If `nm` is set, then the names of the vector or list will be the method ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname count_ords
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available methods.
#'     list_count_ords()
#'
#'     # Find the method for 'prcomp'.
#'     m <- match_count_ord('prcomp')
#'     m$id
#'
list_count_ords <- function (
    val = c('data.frame', 'list', 'func', 'id', 'name'),
    nm  = c(NA, 'id', 'name') ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }

  df <- COUNT_ORDS

  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) count_ord_info(df[df$id == id,]))
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


#' @rdname count_ords
#' @export
#'
match_count_ord <- function (method) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')
  if (is.na(method))         stop('`method` must not be NA')

  df <- COUNT_ORDS

  # Exact match
  if (method %in% df$id)
    return (count_ord_info(df[df$id == method,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  # Methods of the other wrappers.
  if (needle %in% names(exact_names('ORDINATIONS'))) {
    id <- exact_names('ORDINATIONS')[[needle]]
    stop(
      '`method = "', method, '"`: ', id, '() takes distances, not a count table. ',
      'Use dist_ord(dm, "', id, '"), or call ', id, '() directly.' )
  }
  if (needle %in% names(exact_names('DIST_TESTS')))
    stop('`method = "', method, '"` is a test, not an ordination: use dist_test().')

  match <- match_name(needle, 'COUNT_ORDS', df$id)

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  count_ord_info(df[df$id == match,])
}


# One row of COUNT_ORDS as a list, with the method's function and its args.
count_ord_info <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(info$id, ENV)
  info$params  <- names(formals(info$func))
  info
}
