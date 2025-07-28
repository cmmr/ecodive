

#' Number of CPU cores.
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
