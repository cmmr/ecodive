# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



# Ordination methods. `input` is what each method's function takes as `x`:
# "dist" for a dist object (or symmetric matrix) of beta diversity distances.
ORDINATIONS <- local({

  df <- read.table(
    header           = TRUE,
    sep              = ",",
    quote            = "",
    strip.white      = TRUE,
    stringsAsFactors = FALSE,
    na.strings       = "",
    text             = "
      name,                                id,   input, alt_ids
      Distance-based Redundancy Analysis,  dbrda, dist, distancebasedrda
      Non-metric Multidimensional Scaling, nmds, dist,  nonmetricmds metamds monomds
      Principal Coordinates Analysis,      pcoa, dist,  pco mds cmdscale classicalmds metricmds
      t-distributed Stochastic Neighbor Embedding, tsne, dist, rtsne bhtsne barneshuttsne
      Uniform Manifold Approximation and Projection, umap, dist, uwot
  " )

  return (df)
})

ORD_HAYSTACK <- make_haystack(ORDINATIONS)



#' Ordination Wrapper Function
#'
#' Runs any of ecodive's ordination methods by name, as [beta_div()] does for
#' beta diversity metrics.
#'
#' @param x   The data to ordinate: for `"pcoa"`, `"nmds"`, `"tsne"`, and `"umap"`, a
#'        `dist` object such as from [beta_div()], or a symmetric matrix of
#'        distances. See `list_ordinations()` for what each method takes.
#'
#' @param method   The name of an ordination method: `"nmds"`, `"pcoa"`,
#'        `"tsne"`, or `"umap"`. Flexible matching is supported (see below). Programmatic access via
#'        `list_ordinations()`.
#'
#' @param ...   Further arguments for the method's function, such as `k` or
#'        `cpus`. Each method keeps its own defaults; see [pcoa()],
#'        [nmds()], [tsne()], and [umap()].
#'
#' @return The method's result, as described in "Ordination result" below.
#'   `ordination(x, 'pcoa', ...)` is the same as `pcoa(x, ...)`.
#'
#' @details
#'
#' **List of Ordination Methods**
#'
#' | Option / Function Name | Method Name                                    |
#' | :--------------------- | :--------------------------------------------- |
#' | `dbrda`                | Distance-based Redundancy Analysis (db-RDA)    |
#' | `nmds`                 | Non-metric Multidimensional Scaling (NMDS)     |
#' | `pcoa`                 | Principal Coordinates Analysis (PCoA)          |
#' | `tsne`                 | t-distributed Stochastic Neighbor Embedding    |
#' | `umap`                 | Uniform Manifold Approximation and Projection  |
#'
#'
#' **Flexible name matching**
#'
#' Case insensitive and partial matching, ignoring anything but letters. E.g.
#' `method = 'PCoA'` or `'principal coord'` selects `pcoa`. Common other names
#' also work: `'MDS'`, `'PCO'`, `'cmdscale'`, or `'classical MDS'` for
#' `pcoa`, and `'non-metric MDS'`, `'metaMDS'`, or `'monoMDS'` for `nmds`.
#' `'dbRDA'` or `'distance-based RDA'` select `dbrda`, which also needs
#' `data` (see [dbrda()]).
#'
#' Finished code should always use the full primary option name to avoid
#' ambiguity with future additions to the list of methods.
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [nmds()], [tsne()], [umap()], [list_ordinations()]
#'
#' @export
#' @examples
#'     d <- bray(ex_counts)
#'
#'     ordination(d, 'pcoa')
#'
#'     ordination(d, 'NMDS', k = 2)
#'
ordination <- function (x, method, ...) {

  method <- match_ordination(method)

  method$func(x, ...)
}



#' Find and Browse Ordination Methods
#'
#' Programmatic access to the list of available ordination methods, and
#' their associated functions.
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
#' @param input   Filter by the kind of data the method takes: `"dist"`
#'        for distances. Default: `NULL` (no filtering).
#'
#' @return
#'
#' **`match_ordination()`**
#'
#' A `list` with the following elements.
#'
#' * `name` : Method name, e.g. `"Principal Coordinates Analysis"`
#' * `id` : Method ID - also the name of the function, e.g. `"pcoa"`
#' * `input` : The kind of data the method takes, e.g. `"dist"`.
#' * `func` : The function for this method, e.g. `ecodive::pcoa`
#' * `params` : Formal args for `func`, e.g. `c("dm", "k", "cpus")`
#'
#'
#' **`list_ordinations()`**
#'
#' The returned object's type and values are controlled with the `val` and
#' `nm` arguments.
#'
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_ordination()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of method IDs.
#' * `val = "name"` : A character vector of method names.
#' * `val = "input"` : A character vector of the kinds of data each method takes.
#'
#' If `nm` is set, then the names of the vector or list will be the method ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will
#' be applied to the `rownames()` property of the `data.frame`.
#'
#'
#' @rdname ordinations
#' @export
#'
#' @examples
#'
#'     # A data.frame of all available methods.
#'     list_ordinations()
#'
#'     # Find the method for 'MDS'.
#'     m <- match_ordination('MDS')
#'     m$id
#'
list_ordinations <- function (
    val   = c('data.frame', 'list', 'func', 'id', 'name', 'input'),
    nm    = c(NA, 'id', 'name'),
    input = NULL ) {

  val <- match.arg(val)

  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }


  df <- ORDINATIONS

  if (!is.null(input)) {

    if (length(bad <- setdiff(input, ORDINATIONS$input)) > 0)
      stop(
        "Invalid value for `input`: ", paste(collapse = ", ", bad), "\n",
        "Options are: NULL, ", paste(collapse = ", ", sort(unique(ORDINATIONS$input))), "." )

    df <- df[df$input %in% input,, drop = FALSE]
  }


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


#' @rdname ordinations
#' @export
#'
match_ordination <- function (method, input = NULL) {

  if (!is.character(method)) stop('`method` must be a character vector')
  if (length(method) != 1)   stop('`method` must be length 1')

  df <- list_ordinations(input = input)

  # Exact match
  if (method %in% df$id)
    return (ord_info(df[df$id == method,]))

  # Use partial matching
  needle <- gsub('[^a-z]', '', tolower(method))
  if (nchar(needle) == 0) stop('`method` must contain alphabetic characters')

  opts  <- ORD_HAYSTACK[unname(ORD_HAYSTACK) %in% df$id]
  match <- unique(unname(opts[startsWith(names(opts), needle)]))

  if (length(match) > 1) stop(
    '`method = "', method, '"` matches multiple options: ',
    paste(collapse = ', ', sort(match)) )

  if (length(match) == 0) stop(
    '`method = "', method, '"` does not match any option: ',
    paste(collapse = ', ', sort(unique(unname(opts)))) )

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
