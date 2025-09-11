# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit




METRICS <- local({
  
  df <- read.table(
    header           = TRUE, 
    sep              = ",", 
    quote            = "", 
    strip.white      = TRUE, 
    stringsAsFactors = FALSE, 
    na.strings       = "NA", 
    tryLogical       = FALSE, # R 3.6.3 doesn't offer this parameter
    text             = "
      name,                                          id,                        phylo, weighted, true_metric, div,   alt_ids
      Abundance-based Coverage Estimator (ACE),      ace,                       FALSE, TRUE,     NA,          alpha, 
      Aitchison Distance,                            aitchison,                 FALSE, TRUE,     TRUE,        beta,  
      Berger-Parker Index,                           berger,                    FALSE, TRUE,     NA,          alpha, 
      Bhattacharyya Distance,                        bhattacharyya,             FALSE, TRUE,     TRUE,        beta,  
      Bray-Curtis Dissimilarity,                     bray,                      FALSE, TRUE,     FALSE,       beta,  
      Brillouin Index,                               brillouin,                 FALSE, TRUE,     NA,          alpha, 
      Canberra Distance,                             canberra,                  FALSE, TRUE,     TRUE,        beta,  
      Chao1,                                         chao1,                     FALSE, TRUE,     NA,          alpha, 
      Chebyshev Distance,                            chebyshev,                 FALSE, TRUE,     TRUE,        beta,  
      Chord Distance,                                chord,                     FALSE, TRUE,     TRUE,        beta,  
      Clark's Divergence Distance,                   clark,                     FALSE, TRUE,     TRUE,        beta,  
      Dice-Sorensen Dissimilarity,                   sorensen,                  FALSE, FALSE,    FALSE,       beta,  
      Divergence,                                    divergence,                FALSE, TRUE,     TRUE,        beta,  
      Euclidean Distance,                            euclidean,                 FALSE, TRUE,     TRUE,        beta,  
      Faith's Phylogenetic Diversity,                faith,                     TRUE,  FALSE,    NA,          alpha, faithpd
      Fisher's Alpha,                                fisher,                    FALSE, TRUE,     NA,          alpha, 
      Generalized UniFrac (GUniFrac),                generalized_unifrac,       TRUE,  TRUE,     TRUE,        beta,  gunifrac
      Gini-Simpson Index,                            simpson,                   FALSE, TRUE,     NA,          alpha, 
      Gower Distance,                                gower,                     FALSE, TRUE,     TRUE,        beta,  
      Hamming Distance,                              hamming,                   FALSE, FALSE,    TRUE,        beta,  
      Hellinger Distance,                            hellinger,                 FALSE, TRUE,     TRUE,        beta,  
      Horn-Morisita Dissimilarity,                   horn,                      FALSE, TRUE,     FALSE,       beta,  
      Inverse Simpson Index,                         inv_simpson,               FALSE, TRUE,     NA,          alpha, 
      Jaccard Distance,                              jaccard,                   FALSE, FALSE,    TRUE,        beta,  
      Jensen-Shannon Distance,                       jensen,                    FALSE, TRUE,     TRUE,        beta,  
      Jensen-Shannon Divergence (JSD),               jsd,                       FALSE, TRUE,     TRUE,        beta,  
      Lorentzian Distance,                           lorentzian,                FALSE, TRUE,     FALSE,       beta,  
      Manhattan Distance,                            manhattan,                 FALSE, TRUE,     TRUE,        beta,  
      Margalef's Richness Index,                     margalef,                  FALSE, TRUE,     NA,          alpha, 
      Matusita Distance,                             matusita,                  FALSE, TRUE,     TRUE,        beta,  
      McIntosh Index,                                mcintosh,                  FALSE, TRUE,     NA,          alpha, 
      Menhinick's Richness Index,                    menhinick,                 FALSE, TRUE,     NA,          alpha, 
      Minkowski Distance,                            minkowski,                 FALSE, TRUE,     TRUE,        beta,  
      Morisita Dissimilarity,                        morisita,                  FALSE, TRUE,     FALSE,       beta,  
      Motyka Dissimilarity,                          motyka,                    FALSE, TRUE,     FALSE,       beta,  
      Normalized Weighted UniFrac,                   normalized_unifrac,        TRUE,  TRUE,     TRUE,        beta,  nunifrac
      Observed Features,                             observed,                  FALSE, FALSE,    NA,          alpha, otus asvs
      Otsuka-Ochiai Dissimilarity,                   ochiai,                    FALSE, FALSE,    FALSE,       beta,  
      Probabilistic Symmetric Chi-Squared Distance,  psym_chisq,                FALSE, TRUE,     FALSE,       beta,  
      Shannon Diversity Index,                       shannon,                   FALSE, TRUE,     NA,          alpha, 
      Soergel Distance,                              soergel,                   FALSE, TRUE,     TRUE,        beta,  
      Squared Chi-Squared Distance,                  squared_chisq,             FALSE, TRUE,     FALSE,       beta,  
      Squared Chord Distance,                        squared_chord,             FALSE, TRUE,     FALSE,       beta,  
      Squared Euclidean Distance,                    squared_euclidean,         FALSE, TRUE,     FALSE,       beta,  
      Squares Richness Estimator,                    squares,                   FALSE, TRUE,     NA,          alpha, 
      Topsoe Distance,                               topsoe,                    FALSE, TRUE,     TRUE,        beta,  
      Unweighted UniFrac,                            unweighted_unifrac,        TRUE,  TRUE,     TRUE,        beta,  uunifrac
      Variance-Adjusted Weighted UniFrac,            variance_adjusted_unifrac, TRUE,  TRUE,     TRUE,        beta,  vunifrac
      Wave Hedges Distance,                          wave_hedges,               FALSE, TRUE,     FALSE,       beta,  
      Weighted UniFrac,                              weighted_unifrac,          TRUE,  TRUE,     TRUE,        beta,  wunifrac
  ")
  
  for (i in c('phylo', 'weighted', 'true_metric'))
    df[[i]] <- unname(c('TRUE' = TRUE, 'FALSE' = FALSE, 'NA' = NA)[df[[i]]])

  return (df)
})

HAYSTACK <- local({
  result        <- character(0)
  names(result) <- character(0)
  alt_ids <- strsplit(METRICS$alt_ids, ' ', fixed = TRUE)
  for (i in seq_len(nrow(METRICS))) {
    ids <- c(METRICS$name[[i]], METRICS$id[[i]], alt_ids[[i]])
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
#' @param div,phylo,weighted,true_metric   Consider only metrics matching 
#'        specific criteria. For example, `div = "alpha"` will only return 
#'        alpha diversity metrics.
#'        Default: `div=NULL, phylo=NULL, weighted=NULL, true_metric=NULL`
#' 
#' @return 
#' 
#' **`match_metric()`**
#' 
#' A `list` with the following elements.
#' 
#' \describe{
#'   \item{`name` - }{ Metric name, e.g. `"Faith's Phylogenetic Diversity"` }
#'   \item{`id` - }{ Metric ID - also the name of the function, e.g. `"faith"` }
#'   \item{`div` - }{ Either `"alpha"` or `"beta"`. }
#'   \item{`phylo` - }{ `TRUE` if metric requires a phylogenetic tree; `FALSE` otherwise. }
#'   \item{`weighted` - }{ `TRUE` if metric takes relative abundance into account; `FALSE` if it only uses presence/absence. }
#'   \item{`true_metric` - }{ `TRUE` if metric satisfies the triangle inequality; `FALSE` if it is a non-metric dissimilarity; `NA` for alpha diversity metrics. }
#'   \item{`func` - }{ The function for this metric, e.g. `ecodive::faith` }
#'   \item{`params` - }{ Formal args for `func`, e.g. `c("counts", "tree", "cpus")` }
#' }
#' 
#' 
#' **`list_metrics()`**
#' 
#' The returned object's type and values are controlled with the `val` and `nm` arguments.
#' 
#' \describe{
#'   \item{`val = "data.frame"` - }{ The data.frame from which the below options are sourced. }
#'   \item{`val = "list"` - }{ A list of objects as returned by `match_metric()` (above). }
#'   \item{`val = "func"` - }{ A list of functions. }
#'   \item{`val = "id"` - }{ A character vector of metric IDs. }
#'   \item{`val = "name"` - }{ A character vector of metric names. }
#'   \item{`val = "div"` - }{ A character vector `"alpha"` and/or `"beta"`. }
#'   \item{`val = "phylo"` - }{ A logical vector indicating which metrics require a tree. }
#'   \item{`val = "weighted"` - }{ A logical vector indicating which metrics take relative abundance into account (as opposed to just presence/absence). }
#'   \item{`val = "true_metric"` - }{ A logical vector indicating which metrics satisfy the triangle inequality, which work better for ordinations such as PCoA. }
#' }
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
    val = c('data.frame', 'list', 'func', 'id', 'name', 'div', 'phylo', 'weighted', 'true_metric'), 
    nm  = c(NA, 'id', 'name'), phylo = NULL, weighted = NULL, true_metric = NULL ) {
  
  div <- match.arg(div)
  val <- match.arg(val)
  
  if (missing(nm)) { nm <- ifelse(val %in% c('list', 'func'), 'id', NA) }
  else             { nm <- match.arg(nm)                                }
  
  
  #________________________________________________________
  # Subset the metrics according to `div`, `phylo`, etc.
  #________________________________________________________
  
  df <- METRICS
  
  for (k in c('div', 'phylo', 'weighted', 'true_metric')) {
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
    metric, div = NULL, phylo = NULL, weighted = NULL, true_metric = NULL ) {
    
  if (!is.character(metric)) stop('`metric` must be a character vector')
  if (length(metric) != 1)   stop('`metric` must be length 1')
  
  
  df <- list_metrics(
    div         = div, 
    phylo       = phylo, 
    weighted    = weighted, 
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


