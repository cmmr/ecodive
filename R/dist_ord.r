# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



# Ordination methods that take distances (a dist object or a symmetric
# matrix of beta diversity distances). Methods that take a count table are
# in COUNT_ORDS (count_ord.r).
ORDINATIONS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                id,   alt_ids
      Discriminant Canonical Analysis of Principal Coordinates, capda, capdiscriminant discriminantcap
      Constrained Analysis of Principal Coordinates, capscale, cap constrainedanalysisofprincipalcoordinates
      Distance-based Redundancy Analysis,  dbrda, distancebasedrda
      Isometric Feature Mapping,           isomap, isometricmapping
      Metric Multidimensional Scaling,     mds,  metricmds
      Non-metric Multidimensional Scaling, nmds, nonmetricmds metamds monomds isomds
      Principal Coordinates Analysis,      pcoa, pco cmdscale classicalmds
      t-distributed Stochastic Neighbor Embedding, tsne, rtsne bhtsne barneshuttsne
      Uniform Manifold Approximation and Projection, umap, uwot
  " )

  return (df)
})



#' Distance-based Ordination
#'
#' Runs any of ecodive's distance-based ordination methods by name, as
#' [beta_div()] does for beta diversity metrics and [dist_test()] does for
#' statistical tests. Compute the distances once with [beta_div()] (or a
#' metric's own function, such as [bray()]), then pass them to `dist_ord()`
#' and `dist_test()`.
#'
#' @param dm   The distances to ordinate: a `dist` object such as from
#'        [beta_div()], or a symmetric matrix of distances. Every method
#'        listed by `list_dist_ords()` takes distances.
#'
#' @param method   The name of an ordination method: `"capda"`, `"capscale"`,
#'        `"dbrda"`, `"isomap"`, `"mds"`, `"nmds"`, `"pcoa"`, `"tsne"`,
#'        or `"umap"`. Flexible matching is supported (see below).
#'        Programmatic access via `list_dist_ords()`.
#'
#' @param ...   Further arguments for the method's function, such as `n_axes`,
#'        `seed`, `cpus`, or `pcoa`'s `correction`. Each method keeps its own
#'        defaults; see [capscale()] and [dbrda()] (which also need `groups`,
#'        or `data` and `formula`, as in `dist_ord(dm, 'dbrda', 'treatment',
#'        meta)`), [capda()] (which needs `groups`, and `data` if `groups`
#'        names a column of it), [isomap()], [mds()], [nmds()], [pcoa()], [tsne()], and
#'        [umap()]. `seed` is ignored by the methods that have no random
#'        parts (`capscale`, `dbrda`, `isomap` and `pcoa`), so the same
#'        arguments can be given to every method.
#'
#' @return The method's result, as described in "Ordination result" below.
#'   `dist_ord(dm, 'pcoa', ...)` is the same as `pcoa(dm, ...)`.
#'
#' @details
#'
#' **List of Ordination Methods**
#'
#' | Option / Function Name | Method Name                                           |
#' | :--------------------- | :---------------------------------------------------- |
#' | `capda`                | Discriminant Canonical Analysis of Principal Coordinates (CAP), supervised |
#' | `capscale`             | Constrained Analysis of Principal Coordinates         |
#' | `dbrda`                | Distance-based Redundancy Analysis (db-RDA)           |
#' | `isomap`               | Isometric Feature Mapping (Isomap)                    |
#' | `mds`                  | Metric Multidimensional Scaling (MDS)                 |
#' | `nmds`                 | Non-metric Multidimensional Scaling (NMDS)            |
#' | `pcoa`                 | Principal Coordinates Analysis (PCoA)                 |
#' | `tsne`                 | t-distributed Stochastic Neighbor Embedding (t-SNE)   |
#' | `umap`                 | Uniform Manifold Approximation and Projection (UMAP)  |
#'
#' [pca()], [rda()], [cca()], [ca()] and [rpca()] ordinate a count table rather than
#' distances, so they are run by [count_ord()], or called directly. PCA is PCoA of Euclidean
#' distances: `dist_ord(euclidean(counts), 'pcoa')` gives the same scores as
#' `pca(counts)`, up to each axis's sign, but `pca()` also returns the
#' features' loadings and never forms the distance matrix. [permanova()] and
#' ecodive's other tests are run by [dist_test()].
#'
#' Default numbers of axes differ: `n_axes = 3` for the methods whose axes
#' have eigenvalues (`capscale`, `dbrda`, `isomap`, `pcoa`, and [pca()] and
#' [rpca()]; `capda` returns up to 3, and at most one fewer than the number
#' of groups), and `n_axes = 2` for those that fit a layout in a fixed number
#' of dimensions (`mds`, `nmds`, `tsne`, and `umap`).
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters. E.g.
#' `method = 'PCoA'` or `'principal coord'` selects `pcoa`. Common other names
#' also work: `'PCO'`, `'cmdscale'`, or `'classical MDS'` for `pcoa`;
#' `'MDS'` or `'metric MDS'` for `mds`; `'non-metric MDS'`, `'metaMDS'`,
#' `'monoMDS'`, or `'isoMDS'` for `nmds`; `'isometric mapping'` for `isomap`;
#' `'Rtsne'`, `'bhtsne'`, or `'Barnes-Hut t-SNE'` for `tsne`; and `'uwot'`
#' for `umap`. `'dbRDA'` or `'distance-based RDA'` select `dbrda`, and `'CAP'`
#' selects `capscale` (but `'CA'` names [ca()], which takes a count table).
#' `'CAPdiscrim'` (BiodiversityR's name) or `'discriminant CAP'` select
#' `capda`, Anderson & Willis's discriminant CAP.
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#'
#' **Same names in other packages**
#'
#' `ape::pcoa()`, `vegan::pca()`, `vegan::dbrda()`, `vegan::capscale()`,
#' `vegan::isomap()` (whose `k` is the number of neighbours, not axes), and
#' `uwot::umap()` share names with ecodive's functions, and whichever package was attached last
#' masks the other's. Write `ecodive::pcoa()` and so on, or use
#' `dist_ord()`, which always runs ecodive's, to be sure.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [capda()], [capscale()], [dbrda()], [isomap()], [mds()], [nmds()], [pcoa()],
#'   [tsne()], [umap()], [list_dist_ords()]; [count_ord()] for [pca()],
#'   [rda()], [cca()], [ca()] and [rpca()], which take a count table; [dist_test()] for
#'   ecodive's tests.
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     ord <- dist_ord(dm, 'pcoa')
#'     plot(ord, 'Species', iris)
#'     head(ord)
#'     str(ord)
#'
dist_ord <- function (dm, method, ...) {

  method <- match_dist_ord(method)

  args <- list(...)
  if (!'seed' %in% method$params) args$seed <- NULL

  # quote(dm), so that a large `dm` isn't copied into the call.
  do.call(method$func, c(list(quote(dm)), args))
}



#' Find and Browse Distance-based Ordination Methods
#'
#' Programmatic access to the list of ordination methods that [dist_ord()]
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
#' **`match_dist_ord()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Method name, e.g. `"Principal Coordinates Analysis"`
#' * `id` : Method ID - also the name of the function, e.g. `"pcoa"`
#' * `func` : The function for this method, e.g. `ecodive::pcoa`
#' * `params` : Formal args for `func`, e.g. `c("dm", "n_axes", "cpus")`
#'
#'
#' **`list_dist_ords()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_dist_ord()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of method IDs.
#' * `val = "name"` : A character vector of method names.
#'
#' If `nm` is set, then the names of the vector or list will be the method ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname dist_ords
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available methods.
#'     list_dist_ords()
#'
#'     # Find the method for 'MDS'.
#'     m <- match_dist_ord('MDS')
#'     m$id
#'
list_dist_ords <- function (
    val = c('data.frame', 'list', 'func', 'id', 'name'),
    nm  = c(NA, 'id', 'name') ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }

  df <- ORDINATIONS

  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }

  if (val == 'list') {
    result <- lapply(df$id, function (id) ord_info(df[df$id == id,]))
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


#' @rdname dist_ords
#' @export
#'
match_dist_ord <- function (method) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')
  if (is.na(method))         stop('`method` must not be NA')

  df <- ORDINATIONS

  # Exact match
  if (method %in% df$id)
    return (ord_info(df[df$id == method,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  # Methods of the other wrappers.
  if (needle %in% names(exact_names('COUNT_ORDS'))) {
    id <- exact_names('COUNT_ORDS')[[needle]]
    stop(
      '`method = "', method, '"`: ', id, '() takes a count table, not distances. ',
      'Use count_ord(counts, "', id, '"), or call ', id, '() directly.' )
  }
  if (needle %in% names(exact_names('DIST_TESTS')))
    stop('`method = "', method, '"` is a test, not an ordination: use dist_test().')

  match <- match_name(needle, 'ORDINATIONS', df$id)

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(df$id)) )

  ord_info(df[df$id == match,])
}


# One row of ORDINATIONS as a list, with the method's function and its args.
ord_info <- function (row) {
  info         <- as.list(row)
  info$alt_ids <- NULL
  info$func    <- get(info$id, ENV)
  info$params  <- names(formals(info$func))
  info
}
