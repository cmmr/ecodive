# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


# Values worth computing once per session.
CACHE <- new.env(parent = emptyenv())



#' Number of CPU Cores
#' 
#' A thin wrapper around `parallely::availableCores()`. If the `parallely`
#' package is not installed, then it falls back to
#' `parallel::detectCores(all.tests = TRUE, logical = TRUE)`, limited to 2
#' under `R CMD check --as-cran` (when `_R_CHECK_LIMIT_CORES_` is set), as
#' `parallelly` does. Returns `1` if `pthread` support is unavailable or when
#' the number of cpus cannot be determined.
#'
#' Every function's `cpus` argument is capped at `n_cpus()`.
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
  
  # Not cached: parallelly::availableCores() follows options and environment
  # variables that can change during a session, such as `mc.cores`, which
  # the future package sets to 1 inside its workers.
  
  cpus <- 1L
  
  if (pthreads()) {
    
    if (nzchar(system.file(package = 'parallelly'))) {
      n <- do.call(`::`, list('parallelly', 'availableCores'))()
      n <- unname(n)
    }
    else { # nocov start
      # The machine's cores don't change, and counting them can run a
      # shell command (on Linux), so they are counted once.
      if (is.null(CACHE$cores))
        CACHE$cores <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
      n <- CACHE$cores
      # R CMD check --as-cran allows at most 2 cores, as parallelly knows.
      lim <- tolower(Sys.getenv('_R_CHECK_LIMIT_CORES_', ''))
      if (nzchar(lim) && lim != 'false') n <- min(n, 2L)
    } # nocov end

    if (isTRUE(n > 0))
      cpus <- as.integer(n)
  }
  
  return (cpus)
}


pthreads <- function () {
  .Call(C_pthreads)
}



#' Threads of the BLAS Library
#'
#' Get or set how many threads R's BLAS/LAPACK library uses, while R is
#' running.
#'
#' R computes matrix operations such as `%*%`, `crossprod()`, `solve()`,
#' `eigen()`, and `cmdscale()` with the BLAS/LAPACK library it is linked to.
#' Optimized libraries run each large operation on several threads.
#' Environment variables such as `OPENBLAS_NUM_THREADS` set how many before
#' the library loads; `blas_threads()` changes it afterward.
#'
#' The count can be read and set with OpenBLAS, Intel MKL, BLIS, and
#' FlexiBLAS (which passes it on to its backend). It can't with R's reference
#' BLAS (which is single-threaded), Apple's Accelerate, or any other library:
#' there `blas_threads()` returns `NA` and setting does nothing. It never
#' raises an error, warning, or message.
#'
#' Setting is limited to what is known to be safe. The count is only read,
#' not set:
#'
#' * With OpenBLAS's OpenMP builds before version 0.3.34, which take their
#'   count from OpenMP's, and OpenBLAS builds that don't report whether they
#'   use OpenMP, such as Ubuntu 20.04's.
#' * With FlexiBLAS, unless its backend is OpenBLAS (but not an OpenMP
#'   build), BLIS, or MKL, and only while that backend is in use.
#' * In a process forked from the one that loaded ecodive, e.g. by
#'   `parallel::mclapply()`.
#'
#' Otherwise, the count can be lowered and raised again, but not above the
#' highest count read from the library, so that it never has to start
#' threads it didn't already have: a larger `n` is reduced to that. When no
#' count was read, only `n = 1` is set. OpenMP's thread count is left as it
#' is.
#'
#' On Windows, the library is the one in `Rblas.dll` and `Rlapack.dll`, which
#' can be replaced with copies of an optimized library. If both are, each
#' copy has a count of its own: both are set, and the larger is returned.
#'
#' ecodive's own functions don't use the BLAS; their threads are set by
#' `cpus`.
#'
#' @param n   `NULL` to only read the count, or a whole number of at least 1
#'        to set it to. Any other value sets nothing.
#'
#' @return The thread count before the call: an integer of at least 1, or
#'         `NA` when it can't be read. Invisibly when `n` is given.
#'
#' @seealso [n_cpus()]
#'
#' @export
#'
#' @examples
#'     blas_threads()
#'
#'     # Use one thread for a matrix product, then restore the count.
#'     # Restoring NA (count unknown) sets nothing.
#'     old <- blas_threads(1)
#'     x   <- crossprod(matrix(runif(1e4), 100))
#'     blas_threads(old)
#'
blas_threads <- function (n = NULL) {

  # Rlapack.dll is only loaded when R first needs LAPACK. Load it now, so
  # that a copy of an optimized library there is found too.
  if (.Platform$OS.type == 'windows')
    tryCatch(La_version(), error = function (e) NULL, warning = function (w) NULL) # nocov

  if (is.null(n))
    return (.Call(C_blas_threads, NULL))

  valid <- is.numeric(n) && !is.object(n) && length(n) == 1L &&
    !is.na(n) && n >= 1 && n <= .Machine$integer.max && n == trunc(n)

  invisible(.Call(C_blas_threads, if (valid) as.integer(n) else NULL))
}





