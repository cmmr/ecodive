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
    'bray_curtis'                 = bray_curtis, 
    'sorenson'                    = sorenson, 
    'canberra'                    = canberra, 
    'euclidean'                   = euclidean, 
    'gower'                       = gower, 
    'jaccard'                     = jaccard, 
    'kulczynski'                  = kulczynski, 
    'manhattan'                   = manhattan, 
    'unweighted_unifrac'          = unweighted_unifrac, 
    'weighted_unifrac'            = weighted_unifrac, 
    'weighted_normalized_unifrac' = weighted_normalized_unifrac, 
    'generalized_unifrac'         = generalized_unifrac, 
    'variance_adjusted_unifrac'   = variance_adjusted_unifrac )
)



#' Number of CPU Cores
#' 
#' A thin wrapper around 
#' `parallel::detectCores(all.tests = TRUE, logical = TRUE)` which falls back  
#' to `1` when the number of CPU cores cannot be detected, or when the system 
#' does not support `pthreads`. Consider using `parallely::availableCores()` 
#' in place of `n_cpus()` for more advanced interrogation of system resources.
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
      n <- detectCores(all.tests = TRUE, logical = TRUE)
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
