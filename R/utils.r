# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



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





