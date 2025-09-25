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



TRANSFORM_PCT   <- 1L
TRANSFORM_CLR   <- 2L
TRANSFORM_CHORD <- 3L

transform_pct <- function (counts, cpus = n_cpus()) {
  .Call(C_transform, counts, TRANSFORM_PCT, cpus, NULL)
}

transform_clr <- function (counts, pseudocount = NULL, cpus = n_cpus()) {
  if (is.null(pseudocount)) pseudocount <- min(counts[counts > 0])
  .Call(C_transform, counts, TRANSFORM_CLR, cpus, pseudocount)
}

transform_chord <- function (counts, cpus = n_cpus()) {
  .Call(C_transform, counts, TRANSFORM_CHORD, cpus, NULL)
}

