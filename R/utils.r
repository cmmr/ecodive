# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' API for browsing functions
#' 
#' Programmatic access to the lists of available metrics, and their associated
#' functions.
#' 
#' @param metric   The name, partial name, alternate name, or partial alternate
#'   name of an alpha/beta diversity metric.
#'   
#' @param evar   The name of the variable to display in an error message.
#'   
#' @param multiple   If `TRUE`, allow `length(metric)` to be longer than 1.
#' 
#' @return   `metrics` is a nested list of functions for alpha and beta
#'   diversity. `as_alpha_metric()` and `as_beta_metric()` return a validated
#'   metric name or throws an error.
#' 
#' @export
#' 
#' @examples
#'     names(metrics)
#'     
#'     names(metrics$alpha)
#'     
#'     names(formals(metrics$alpha$faith))
#'     
#'     metrics$alpha$faith(ex_counts, ex_tree)
#'     
#'     as_alpha_metric('otus')
#'     
#'     as_beta_metric('guni')
#' 

metrics <-list(
  
  alpha = list(
    'observed'    = observed, 
    'chao1'       = chao1, 
    'shannon'     = shannon, 
    'simpson'     = simpson, 
    'inv_simpson' = inv_simpson, 
    'faith'       = faith ),
  
  beta = list(
    'bray_curtis'               = bray_curtis, 
    'sorenson'                  = sorenson, 
    'canberra'                  = canberra, 
    'euclidean'                 = euclidean, 
    'gower'                     = gower, 
    'jaccard'                   = jaccard, 
    'kulczynski'                = kulczynski, 
    'manhattan'                 = manhattan, 
    'unweighted_unifrac'        = unweighted_unifrac, 
    'weighted_unifrac'          = weighted_unifrac, 
    'normalized_unifrac'        = normalized_unifrac, 
    'generalized_unifrac'       = generalized_unifrac, 
    'variance_adjusted_unifrac' = variance_adjusted_unifrac )
)


#' @rdname metrics
#' @export
as_alpha_metric <- function (metric, multiple = FALSE, evar = 'metric') {
  
  metric_match(
    metric   = metric,
    multiple = multiple,
    evar     = evar,
    primary  = names(metrics$alpha),
    mapped   = c(
      otus       = 'observed', 
      invsimpson = 'inv_simpson', 
      faith_pd   = 'faith' ))
}


#' @rdname metrics
#' @export
as_beta_metric <- function (metric, multiple = FALSE, evar = 'metric') {
  
  metric_match(
    metric   = metric,
    multiple = multiple,
    evar     = evar,
    primary  = names(metrics$beta),
    mapped   = c(
      uunifrac  = 'unweighted_unifrac', 
      u_unifrac = 'unweighted_unifrac', 
      wunifrac  = 'weighted_unifrac', 
      w_unifrac = 'weighted_unifrac', 
      nunifrac  = 'normalized_unifrac', 
      n_unifrac = 'normalized_unifrac', 
      gunifrac  = 'generalized_unifrac', 
      g_unifrac = 'generalized_unifrac', 
      vunifrac  = 'variance_adjusted_unifrac', 
      v_unifrac = 'variance_adjusted_unifrac' ))
}


metric_match <- function (metric, multiple, evar, primary, mapped) {
  
  if (!is.character(metric)) stop(evar, ' must be a character vector')
  if (any(is.na(metric)))    stop(evar, ' cannot be NA')
  
  if (isTRUE(multiple))
    if (length(metric) != 1) stop(evar, ' must be length 1')
  
  newick <- trimws(metric)
  if (any(nchar(metric) == 0)) stop(evar, ' cannot be ""')
  
  
  sapply(metric, USE.NAMES = FALSE, function (x) {
    
    x <- gsub('[^a-z]+', '_', tolower(x))
    
    i <- which(startsWith(primary, x))
    if (length(i) == 1) return (primary[[i]])
    
    if (length(i) > 1)
      stop(
        '`', evar, '` "', x, '" matches multiple options: ', 
        paste(collapse = ',', primary[i]) )
    
    
    i <- which(startsWith(names(mapped), x))
    if (length(i) == 1) return (mapped[[i]])
    
    stop(
      'invalid `', evar, ':` "', x, '"\n', 
      '`', evar, '` must be one of: ', paste(collapse = ',', primary))
  })
}



#' Number of CPU Cores
#' 
#' A thin wrapper around `parallely::availableCores()`. If the `parallely`
#' package is not installed, then it falls back to  
#' `parallel::detectCores(all.tests = TRUE, logical = TRUE)`. Returns `1` if
#' `pthread` support is unavailable or when the number of cpus cannot be
#' determined.
#' 
#' @return   A scalar integer, guaranteed to be at least `1`.
#' 
#' @importFrom parallel detectCores
#' @export
#' 
#' @examples
#'     n_cpus()
#' 
n_cpus <- function () {
  
  if (!n_cpus_cached) {
    n_cpus_cached <- 1L
    
    if (pthreads()) {
      
      if (nzchar(system.file(package = 'parallelly'))) {
        n <- do.call(`::`, list('parallelly', 'availableCores'))()
        n <- unname(n)
      }
      else {
        n <- parallel::detectCores(all.tests = TRUE, logical = TRUE) # nocov
      }
      
      if (isTRUE(n > 0))
        n_cpus_cached <- n
    }
  }
  
  return (n_cpus_cached)
}

n_cpus_cached <- 0


pthreads <- function () {
  .Call(C_pthreads)
}
