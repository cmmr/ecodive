# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit




METRICS <- local({
  
  df <- read.table(
    header           = TRUE, 
    sep              = ",", 
    quote            = "", 
    strip.white      = TRUE, 
    stringsAsFactors = FALSE, 
    na.strings       = "", 
    text             = "
      name,                                          id,                        phylo, weighted, int_only, true_metric, div,   alt_ids
      Abundance-based Coverage Estimator (ACE),      ace,                       NO,    YES,      YES,      NA,          alpha, 
      Aitchison Distance,                            aitchison,                 NO,    YES,      NO,       YES,         beta,  
      Berger-Parker Index,                           berger,                    NO,    YES,      NO,       NA,          alpha, 
      Bhattacharyya Distance,                        bhattacharyya,             NO,    YES,      NO,       YES,         beta,  
      Bray-Curtis Dissimilarity,                     bray,                      NO,    YES,      NO,       NO,          beta,  
      Brillouin Index,                               brillouin,                 NO,    YES,      YES,      NA,          alpha, 
      Canberra Distance,                             canberra,                  NO,    YES,      NO,       YES,         beta,  
      Chao1,                                         chao1,                     NO,    YES,      YES,      NA,          alpha, 
      Chebyshev Distance,                            chebyshev,                 NO,    YES,      NO,       YES,         beta,  
      Chord Distance,                                chord,                     NO,    YES,      NO,       YES,         beta,  
      Clark's Divergence Distance,                   clark,                     NO,    YES,      NO,       YES,         beta,  
      Dice-Sorensen Dissimilarity,                   sorensen,                  NO,    NO,       NO,       NO,          beta,  
      Divergence,                                    divergence,                NO,    YES,      NO,       YES,         beta,  
      Euclidean Distance,                            euclidean,                 NO,    YES,      NO,       YES,         beta,  
      Faith's Phylogenetic Diversity,                faith,                     YES,   NO,       NO,       NA,          alpha, faithpd
      Fisher's Alpha,                                fisher,                    NO,    YES,      YES,      NA,          alpha, 
      Generalized UniFrac (GUniFrac),                generalized_unifrac,       YES,   YES,      NO,       YES,         beta,  gunifrac
      Gini-Simpson Index,                            simpson,                   NO,    YES,      NO,       NA,          alpha, 
      Gower Distance,                                gower,                     NO,    YES,      NO,       YES,         beta,  
      Hamming Distance,                              hamming,                   NO,    NO,       NO,       YES,         beta,  
      Hellinger Distance,                            hellinger,                 NO,    YES,      NO,       YES,         beta,  
      Horn-Morisita Dissimilarity,                   horn,                      NO,    YES,      NO,       NO,          beta,  
      Inverse Simpson Index,                         inv_simpson,               NO,    YES,      NO,       NA,          alpha, 
      Jaccard Distance,                              jaccard,                   NO,    NO,       NO,       YES,         beta,  
      Jensen-Shannon Distance,                       jensen,                    NO,    YES,      NO,       YES,         beta,  
      Jensen-Shannon Divergence (JSD),               jsd,                       NO,    YES,      NO,       YES,         beta,  
      Lorentzian Distance,                           lorentzian,                NO,    YES,      NO,       NO,          beta,  
      Manhattan Distance,                            manhattan,                 NO,    YES,      NO,       YES,         beta,  
      Margalef's Richness Index,                     margalef,                  NO,    YES,      YES,      NA,          alpha, 
      Matusita Distance,                             matusita,                  NO,    YES,      NO,       YES,         beta,  
      McIntosh Index,                                mcintosh,                  NO,    YES,      YES,      NA,          alpha, 
      Menhinick's Richness Index,                    menhinick,                 NO,    YES,      YES,      NA,          alpha, 
      Minkowski Distance,                            minkowski,                 NO,    YES,      NO,       YES,         beta,  
      Morisita Dissimilarity,                        morisita,                  NO,    YES,      YES,      NO,          beta,  
      Motyka Dissimilarity,                          motyka,                    NO,    YES,      NO,       NO,          beta,  
      Normalized Weighted UniFrac,                   normalized_unifrac,        YES,   YES,      NO,       YES,         beta,  nunifrac
      Observed Features,                             observed,                  NO,    NO,       NO,       NA,          alpha, otus asvs
      Otsuka-Ochiai Dissimilarity,                   ochiai,                    NO,    NO,       NO,       NO,          beta,  
      Probabilistic Symmetric Chi-Squared Distance,  psym_chisq,                NO,    YES,      NO,       NO,          beta,  
      Shannon Diversity Index,                       shannon,                   NO,    YES,      NO,       NA,          alpha, 
      Soergel Distance,                              soergel,                   NO,    YES,      NO,       YES,         beta,  
      Squared Chi-Squared Distance,                  squared_chisq,             NO,    YES,      NO,       NO,          beta,  
      Squared Chord Distance,                        squared_chord,             NO,    YES,      NO,       NO,          beta,  
      Squared Euclidean Distance,                    squared_euclidean,         NO,    YES,      NO,       NO,          beta,  
      Squares Richness Estimator,                    squares,                   NO,    YES,      YES,      NA,          alpha, 
      Topsoe Distance,                               topsoe,                    NO,    YES,      NO,       YES,         beta,  
      Unweighted UniFrac,                            unweighted_unifrac,        YES,   YES,      NO,       YES,         beta,  uunifrac
      Variance-Adjusted Weighted UniFrac,            variance_adjusted_unifrac, YES,   YES,      NO,       YES,         beta,  vunifrac
      Wave Hedges Distance,                          wave_hedges,               NO,    YES,      NO,       NO,          beta,  
      Weighted UniFrac,                              weighted_unifrac,          YES,   YES,      NO,       YES,         beta,  wunifrac
  ")
  
  # R 3.6.3 doesn't offer read.table(tryLogical) parameter
  for (i in c('int_only', 'phylo', 'weighted', 'true_metric'))
    df[[i]] <- unname(c('YES' = TRUE, 'NO' = FALSE, 'NA' = NA)[df[[i]]])

  return (df)
})

HAYSTACK <- local({
  result        <- character(0)
  names(result) <- character(0)
  alt_ids <- strsplit(METRICS$alt_ids, ' ', fixed = TRUE)
  for (i in seq_len(nrow(METRICS))) {
    ids <- c(METRICS$name[[i]], METRICS$id[[i]], na.omit(alt_ids[[i]]))
    ids <- unique(gsub('[^a-z]', '', tolower(ids)))
    ids <- ids[order(nchar(ids), decreasing = TRUE)]
    res <- rep_len(METRICS$id[[i]], length(ids))
    names(res) <- ids
    for (j in seq_along(res))
      if (!any(startsWith(names(result), names(res)[j])))
        result <- c(result, res[j])
  }
  return (result)
})

ENV <- environment()



#' Find and Browse Available Metrics
#' 
#' Programmatic access to the lists of available metrics, and their associated
#' functions.
#' 
#' @param metric   The name of an alpha/beta diversity metric to search for. 
#'        Supports partial matching. All non-alpha characters are ignored.
#'   
#' @param val   Sets the return value for this function call. See "Value" 
#'        section below. Default: `"data.frame"`
#'   
#' @param nm   What value to use for the names of the returned object.
#'        Default is `"id"` when `val` is `"list"` or `"func"`, otherwise the 
#'        default is `NA` (no name).
#'   
#' @param div,phylo,weighted,int_only,true_metric   Consider only metrics 
#'        matching specific criteria. For example, `div = "alpha"` will only 
#'        return alpha diversity metrics.
#'        Default: `NULL`
#' 
#' @return 
#' 
#' **`match_metric()`**
#' 
#' A `list` with the following elements.
#' 
#' * `name` : Metric name, e.g. `"Faith's Phylogenetic Diversity"`
#' * `id` : Metric ID - also the name of the function, e.g. `"faith"`
#' * `div` : Either `"alpha"` or `"beta"`.
#' * `phylo` : `TRUE` if metric requires a phylogenetic tree; `FALSE` otherwise.
#' * `weighted` : `TRUE` if metric takes relative abundance into account; `FALSE` if it only uses presence/absence.
#' * `int_only` : `TRUE` if metric requires integer counts; `FALSE` otherwise.
#' * `true_metric` : `TRUE` if metric is a true metric and satisfies the triangle inequality; `FALSE` if it is a non-metric dissimilarity; `NA` for alpha diversity metrics.
#' * `func` : The function for this metric, e.g. `ecodive::faith`
#' * `params` : Formal args for `func`, e.g. `c("counts", "tree", "cpus")`
#' 
#' 
#' **`list_metrics()`**
#' 
#' The returned object's type and values are controlled with the `val` and `nm` arguments.
#' 
#' * `val = "data.frame"` : The data.frame from which the below options are sourced.
#' * `val = "list"` : A list of objects as returned by `match_metric()` (above).
#' * `val = "func"` : A list of functions.
#' * `val = "id"` : A character vector of metric IDs.
#' * `val = "name"` : A character vector of metric names.
#' * `val = "div"` : A character vector `"alpha"` and/or `"beta"`.
#' * `val = "phylo"` : A logical vector indicating which metrics require a tree.
#' * `val = "weighted"` : A logical vector indicating which metrics take relative abundance into account (as opposed to just presence/absence).
#' * `val = "int_only"` : A logical vector indicating which metrics require integer counts.
#' * `val = "true_metric"` : A logical vector indicating which metrics are true metrics and satisfy the triangle inequality, which work better for ordinations such as PCoA.
#' 
#' If `nm` is set, then the names of the vector or list will be the metric ID
#' (`nm="id"`) or name (`nm="name"`). When `val="data.frame"`, the names will be
#' applied to the `rownames()` property of the `data.table`.
#' 
#' 
#' @rdname metrics
#' @export
#' 
#' @examples
#' 
#'     # A data.frame of all available metrics.
#'     head(list_metrics())
#'     
#'     # All alpha diversity function names.
#'     list_metrics('alpha', val = 'id')
#'     
#'     # Try to find a metric named 'otus'.
#'     m <- match_metric('otus')
#'     
#'     # The result is a list that includes the function.
#'     str(m)
#' 

list_metrics <- function (
    div = c(NA, 'alpha', 'beta'), 
    val = c('data.frame', 'list', 'func', 'id', 'name', 'div', 'phylo', 'weighted', 'int_only', 'true_metric'), 
    nm  = c(NA, 'id', 'name'), phylo = NULL, weighted = NULL, int_only = NULL, true_metric = NULL ) {
  
  div <- match.arg(div)
  val <- match.arg(val)
  
  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }
  
  
  #________________________________________________________
  # Subset the metrics according to `div`, `phylo`, etc.
  #________________________________________________________
  
  df <- METRICS
  
  for (k in c('div', 'phylo', 'weighted', 'int_only', 'true_metric')) {
    v <- get(k, inherits = FALSE)
    
    if (!is.null(v) && !is.na(v)) {
      
      if (length(bad <- setdiff(v, METRICS[[k]])) > 0)
        stop(
          "Invalid value for `", k, "`: ", paste(collapse = ", ", bad), "\n",
          "Options are: NULL, ", paste(collapse = ", ", sort(unique(METRICS[[k]]))), "." )
      
      df <- df[df[[k]] %in% v,, drop=FALSE]
    }
  }
  
  
  #________________________________________________________
  # Construct the result object.
  #________________________________________________________
  
  if (val == 'data.frame') {
    if (!is.na(nm)) rownames(df) <- df[[nm]]
    df$alt_ids <- NULL
    return (df)
  }
  
  if (val == 'list') {
    result <- apply(df, 1L, function (row) {
      row         <- as.list(row)
      row$alt_ids <- NULL
      row$func    <- get(row$id, ENV)
      row$params  <- names(formals(row$func))
      return(row)
    })
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


#' @rdname metrics
#' @export
#' 
match_metric <- function (
    metric, div = NULL, phylo = NULL, weighted = NULL, int_only = NULL, true_metric = NULL ) {
    
  if (!is.character(metric)) stop('`metric` must be a character vector')
  if (length(metric) != 1)   stop('`metric` must be length 1')
  
  
  df <- list_metrics(
    div         = div, 
    phylo       = phylo, 
    weighted    = weighted, 
    int_only    = int_only,
    true_metric = true_metric )
  
  
  # Exact match
  if (metric %in% df$id) {
    metric <- as.list(df[df$id == metric,])
  }
  
  # Use partial matching
  else {
    
    needle <- gsub('[^a-z]', '', tolower(metric))
    if (nchar(needle) == 0) stop('`metric` must contain alphabetic characters')
    
    opts  <- HAYSTACK[unname(HAYSTACK) %in% df$id]
    match <- unique(unname(opts[startsWith(names(opts), needle)]))
    
    if (length(match) > 1) stop(
      '`metric = "', metric, '"` matches multiple options: ', 
      paste(collapse = ', ', sort(match)) )
    
    if (length(match) == 0) stop(
      '`metric = "', metric, '"` does not match any option: \n', 
      paste(collapse = '\n', strwrap(paste(collapse = ', ', sort(unique(unname(opts)))))) )
    
    metric <- as.list(df[df$id == match,])
  }
  
  metric$alt_ids <- NULL
  metric$func    <- get(metric$id, ENV)
  metric$params  <- names(formals(metric$func))
  
  return (metric)
}


