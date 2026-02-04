# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Rarefy Observation Counts
#' 
#' Sub-sample observations from a feature table such that all samples have the 
#' same library size (depth). This is performed via random sampling without 
#' replacement.
#' 
#' @param counts  A numeric matrix or sparse matrix object (e.g., `dgCMatrix`).
#'        Counts must be integers.
#' 
#' @param depth   The number of observations to keep per sample. If `NULL` 
#'        (the default), a depth is auto-selected to maximize data retention.
#'     
#' @param seed    An integer seed for the random number generator. Providing 
#'        the same seed guarantees reproducible results. Default: `0`
#'     
#' @param times   The number of independent rarefactions to perform. If set, 
#'        returns a list of matrices. Seeds for subsequent iterations are 
#'        sequential (`seed`, `seed + 1`, ...). Default: `NULL`
#'     
#' @param drop    Logical. If `TRUE`, samples with fewer than `depth` 
#'        observations are discarded. If `FALSE`, they are kept with their 
#'        original counts. Default: `TRUE`
#'        
#' @param warn    Logical. If `TRUE`, emits a warning when samples are dropped 
#'        or returned unrarefied due to insufficient depth. 
#'        Default: `interactive()`
#' 
#' @return A rarefied matrix. The output class (`matrix`, `dgCMatrix`, etc.) 
#'         matches the input class.
#' 
#' @section Auto-Depth Selection:
#'   If `depth` is `NULL`, the function defaults to the highest depth that retains 
#'   at least 10% of the total observations in the dataset.
#' 
#' @section Dropping vs. Retaining Samples:
#'   If a sample has fewer observations than the specified `depth`:
#'   
#'   * `drop = TRUE` (Default): The sample is removed from the output matrix.
#'   * `drop = FALSE`: The sample is returned **unmodified** (with its original 
#'     counts). It is *not* rarefied or zeroed out.
#' 
#' @section Zero-Sum Features:
#'   Features (OTUs, ASVs, Genes) that lose all observations during rarefaction 
#'   are **always retained** as columns/rows of zeros. This ensures the output 
#'   matrix dimensions remain consistent with the input (barring dropped samples).
#' 
#' @export
#' @examples
#'     # A 4-sample x 5-OTU matrix with samples in rows.
#'     counts <- matrix(c(0,0,0,0,0,8,9,10,5,5,5,5,2,0,0,0,6,5,7,0), 4, 5,
#'       dimnames = list(LETTERS[1:4], paste0('OTU', 1:5)))
#'     counts
#'     rowSums(counts)
#'     
#'     # Rarefy all samples to a depth of 18.
#'     # Sample 'A' (13 counts) and 'D' (15 counts) will be dropped.
#'     r_mtx <- rarefy(counts, depth = 18)
#'     r_mtx
#'     rowSums(r_mtx)
#'     
#'     # Keep under-sampled samples by setting `drop = FALSE`.
#'     # Samples 'A' and 'D' are returned with their original counts.
#'     r_mtx <- rarefy(counts, depth = 18, drop = FALSE)
#'     r_mtx
#'     rowSums(r_mtx)
#'     
#'     # Perform 3 independent rarefactions.
#'     r_list <- rarefy(counts, times = 3)
#'     length(r_list)
#'     
#'     # Sparse matrices are supported and their class is preserved.
#'     if (requireNamespace('Matrix', quietly = TRUE)) {
#'       counts_dgC <- Matrix::Matrix(counts, sparse = TRUE)
#'       str(rarefy(counts_dgC))
#'     }
#' 
rarefy <- function (
    counts, 
    depth  = NULL, 
    seed   = 0, 
    times  = NULL, 
    drop   = TRUE, 
    margin = 1L, 
    cpus   = n_cpus(),
    warn   = interactive() ) {
  
  validate_args()
  assert_integer_counts()
  
  
  # Ensure counts are double
  if (is.matrix(counts)) {
    if (!is.double(counts)) counts[] <- as.double(counts)
  }
  else if (inherits(counts, "simple_triplet_matrix")) {
    if (!is.double(counts$v)) counts$v <- as.double(counts$v)
  }
  
  
  # Calculate library sizes first, as they are needed for both
  # auto-selection AND the warning check.
  if (is.null(depth) || isTRUE(warn)) {
    
    if (is.matrix(counts)) {
      if (margin == 1L) sums <- rowSums(counts)
      else              sums <- colSums(counts)
    }
    else if (inherits(counts, "simple_triplet_matrix")) {
      if (margin == 1L) sums <- slam::row_sums(counts)
      else              sums <- slam::col_sums(counts)
    }
    else {
      # Matrix package objects
      if (margin == 1L) sums <- Matrix::rowSums(counts)
      else              sums <- Matrix::colSums(counts)
    }
    
    
    if (is.null(depth)) {
      
      # Auto-select depth:
      # Select the lowest depth that still retains at least 10% 
      # of the total observations in the dataset.
      sums_sorted <- sort(sums)
      n_sams      <- length(sums_sorted)
      total_reads <- sum(sums_sorted)
      threshold   <- 0.1 * total_reads
      
      # Vectorized calculation:
      yields <- sums_sorted * seq.int(n_sams, 1L)
      
      # Find the first index (lowest depth) that meets the threshold
      idx <- which(yields >= threshold)[1L]
      
      if (!is.na(idx)) {
        depth <- as.integer(sums_sorted[idx])
      } else {
        # Fallback: Use the max depth (retains 1 sample, but 100% of its reads).
        # To reach the line, we need a distribution where NO depth meets the 10% 
        # retention threshold, such as a Zipfian distribution (1/rank) 
        # with N > 23,000 samples.
        depth <- as.integer(max(sums_sorted)) # nocov
      }
    }
    
    
    # Warning logic:
    # Check if any samples have insufficient depth and warn if requested.
    if (isTRUE(warn)) {
      
      n_insufficient <- sum(sums < depth)
      
      if (n_insufficient > 0) {
        
        action <- if (isTRUE(drop)) "dropped" else "returned unrarefied"
        
        warning(sprintf(
          "%d samples have fewer than %d observations and will be %s.",
          n_insufficient, depth, action
        ), call. = FALSE)
      }
    }
  }
  
  
  # Call C function
  if (is.null(times)) {
    result <- .Call(C_rarefy, counts, depth, seed, margin, cpus)
  } else {
    seeds  <- ((seed + 2**31 - 1 + seq_len(times)) %% 2**32) - 2**31
    result <- lapply(seeds, function (seed) {
      .Call(C_rarefy, counts, depth, seed, margin, cpus)
    })
  }
  
  
  # Drop samples with insufficient depth
  if (drop) {
    
    dropper <- function (m) {
      if (margin == 1L) {
        # Row samples
        if (is.matrix(m)) sums <- rowSums(m)
        else if (inherits(m, "simple_triplet_matrix")) sums <- slam::row_sums(m)
        else sums <- Matrix::rowSums(m)
        
        # We must keep samples where sum >= depth
        return(m[sums >= depth, , drop = FALSE])
        
      } else {
        # Col samples
        if (is.matrix(m)) sums <- colSums(m)
        else if (inherits(m, "simple_triplet_matrix")) sums <- slam::col_sums(m)
        else sums <- Matrix::colSums(m)
        
        return(m[, sums >= depth, drop = FALSE])
      }
    }
    
    if (is.null(times)) result <- dropper(result)
    else                result <- lapply(result, dropper)
  }
  
  return (result)
}
