# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Rarefy OTU counts.
#' 
#' Sub-sample OTU observations such that all samples have an equal number.
#' If called on data with non-integer abundances, values will be re-scaled to 
#' integers between 1 and `depth` such that they sum to `depth`.
#' 
#' @inherit documentation
#' 
#' @param depth   How many observations to keep per sample. When 
#'        `0 < depth < 1`, it is taken as the minimum percentage of the 
#'        dataset's observations to keep. Ignored when `n_samples` is 
#'        specified. Default: `0.1`
#'
#' @param n_samples   The number of samples to keep. When `0 < n_samples < 1`, 
#'        it is taken as the percentage of samples to keep. If negative, that 
#'        number of samples is dropped. If `0`, all samples are kept. 
#'        If `NULL`, then `depth` is used instead. Default: `NULL`
#'     
#' @param seed   An integer seed for randomizing which observations to keep or 
#'        drop. If you need to create different random rarefactions of the same 
#'        data, set the seed to a different number each time. Default: `0`
#'     
#' @param times   How many independent rarefactions to perform. If set, 
#'        `rarefy()` will return a list of matrices. The seeds for each matrix
#'        will be sequential, starting from `seed`. Default: `NULL`
#'     
#' @param drop   Drop rows and columns with zero observations after rarefying.
#'        Default: `TRUE`
#' 
#' 
#' @return A rarefied matrix. `Matrix` and `slam` objects will be returned with 
#'         the same type; otherwise a base R `matrix` will be returned.
#' 
#' @export
#' @examples
#'     # A 4-sample x 5-OTU matrix with samples in rows.
#'     counts <- matrix(c(0,0,0,0,0,8,9,10,5,5,5,5,2,0,0,0,6,5,7,0), 4, 5,
#'       dimnames = list(LETTERS[1:4], paste0('OTU', 1:5)))
#'     counts
#'     rowSums(counts)
#'     
#'     # Rarefy all samples to a depth of 13.
#'     # Note that sample 'A' has 0 counts and is dropped.
#'     r_mtx <- rarefy(counts, depth = 13, seed = 1)
#'     r_mtx
#'     rowSums(r_mtx)
#'     
#'     # Keep zero-sum rows and columns by setting `drop = FALSE`.
#'     rarefy(counts, depth = 13, drop = FALSE, seed = 1)
#'     
#'     # Rarefy to the depth of the 2nd most abundant sample (B, depth=22).
#'     rarefy(counts, n_samples = 2, seed = 1)
#'     
#'     # Perform 3 independent rarefactions.
#'     r_list <- rarefy(counts, depth = 13, times = 3, seed = 1)
#'     length(r_list)
#'     r_list[[1]]
#'     
#'     # The class of the input matrix is preserved.
#'     if (requireNamespace('Matrix', quietly = TRUE)) {
#'       counts_dgC <- Matrix::Matrix(counts, sparse = TRUE)
#'       class(counts_dgC)
#'       r_dgC <- rarefy(counts_dgC, depth = 13, seed = 1)
#'       class(r_dgC)
#'     }
#' 
rarefy <- function (
    counts, 
    depth     = 0.1, 
    n_samples = NULL, 
    seed      = 0, 
    times     = NULL, 
    drop      = TRUE, 
    margin    = 1L, 
    cpus      = n_cpus() ) {
  
  validate_args()
  assert_integer_counts()
  

  # Ensure counts are double
  if (is.matrix(counts)) {
    if (!is.double(counts)) counts[] <- as.double(counts)
  }
  else if (inherits(counts, "simple_triplet_matrix")) {
    if (!is.double(counts$v)) counts$v <- as.double(counts$v)
  }


  # Call C function
  if (is.null(times)) {
    result <- .Call(C_rarefy, counts, depth, n_samples, seed, margin, cpus)
  } else {
    seeds  <- ((seed + 2**31 - 1 + seq_len(times)) %% 2**32) - 2**31
    result <- lapply(seeds, function (seed) {
      .Call(C_rarefy, counts, depth, n_samples, seed, margin, cpus)
    })
  }
  
  
  # Drop zero-sum rows/cols
  if (drop) {
    
    # Generic function to drop zero-sum rows/cols from any matrix type.
    dropper <- function (m) {

      if (is.matrix(m)) {
        row_sums <- rowSums(m)
        col_sums <- colSums(m)
      } else if (inherits(m, "simple_triplet_matrix")) {
        row_sums <- slam::row_sums(m)
        col_sums <- slam::col_sums(m)
      } else {
        row_sums <- Matrix::rowSums(m)
        col_sums <- Matrix::colSums(m)
      }

      m[row_sums > 0, col_sums > 0, drop = FALSE]
    }
    
    if (is.null(times)) result <- dropper(result)
    else                result <- lapply(result, dropper)
  }
  
  return (result)
}
