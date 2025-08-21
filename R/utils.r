# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' API for browsing functions
#' 
#' Programmatic access to the lists of available metrics, and their associated
#' functions.
#' 
#' @return   A nested list of functions for alpha and beta diversity.
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


as_alpha_metric <- function (metric) {
  
  metric_match(
    metric  = metric,
    primary = names(metrics$alpha),
    mapped  = c(
      otus       = 'observed', 
      invsimpson = 'inv_simpson', 
      faith_pd   = 'faith' ))
}


as_beta_metric <- function (metric) {
  
  metric_match(
    metric  = metric,
    primary = names(metrics$beta),
    mapped  = c(
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


metric_match <- function (metric, primary, mapped) {
  
  tryCatch({
      stopifnot(is.character(metric))
      stopifnot(length(metric) == 1)
      stopifnot(!is.na(metric))
      newick <- trimws(metric)
      stopifnot(nchar(metric) > 0)
    },
    
    error = function (e) 
      stop(e$message, '\n`metric` must be a character string.')
  )
  
  metric <- gsub('[^a-z]+', '_', tolower(metric))
  
  i <- which(startsWith(primary, metric))
  if (length(i) == 1) return (primary[[i]])
  
  if (length(i) > 1)
    stop(
      '`metric` "', metric, '" matches multiple options: ', 
      paste(collapse = ',', primary[i]) )
  
  
  i <- which(startsWith(names(mapped), metric))
  if (length(i) == 1) return (mapped[[i]])
  
  stop(
    'invalid `metric:` "', metric, '"\n', 
    '`metric` must be one of: ', paste(collapse = ',', primary))
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
