# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit



#' Rarefy Observation Counts
#' 
#' Sub-sample observations from a feature table such that all samples have the 
#' same library size (depth). This is performed via random sampling without 
#' replacement.
#' 
#' @inherit documentation
#' 
#' @param counts  A numeric matrix or sparse matrix object (e.g., `dgCMatrix`).
#'        Counts must be integers.
#' 
#' @param depth   The number of observations to keep per sample: one number
#'        for every sample, or one per sample (matched to the samples by name
#'        when both are named, else in the samples' order), as
#'        [qmp_rarefy()] uses. If `NULL` (the default), a depth is
#'        auto-selected to maximize data retention.
#'     
#' @param seed    An integer seed for the random number generator. Providing 
#'        the same seed guarantees reproducible results. Default: `0`
#'     
#' @param times   The number of independent rarefactions to perform. If set, 
#'        returns a list of matrices. Seeds for subsequent iterations are 
#'        sequential (`seed`, `seed + 1`, ...). Default: `NULL`
#'     
#' @param drop    Logical. If `TRUE`, samples with fewer than `depth` 
#'        observations (their own depth, when there is one per sample) are
#'        discarded. If `FALSE`, they are kept with their 
#'        original counts. Default: `TRUE`
#'        
#' @param warn    Logical. If `TRUE`, emits a warning when samples are dropped 
#'        or returned unrarefied due to insufficient depth. 
#'        Default: `interactive()`
#' 
#' @param coverage  Rarefy to equal coverage instead of equal depth: a
#'        sample coverage greater than 0 and at most 1, such as `0.95`. Each
#'        sample is rarefied to the smallest depth whose expected coverage
#'        reaches it (see "Coverage-Based Rarefaction"). Give `depth` or
#'        `coverage`, not both. Default: `NULL` (rarefy to `depth`).
#' 
#' @seealso [qmp_rarefy()], which rarefies each sample to a depth set by
#'   its measured microbial load; [rarefaction_curve()], to choose a depth
#'   or a coverage.
#' 
#' @return A rarefied matrix. The output class (`matrix`, `dgCMatrix`, etc.) 
#'         matches the input class.
#' 
#' @section Auto-Depth Selection:
#'   If `depth` is `NULL`, the depth is the smallest library size at which
#'   the observations kept (that depth times the number of samples with at
#'   least that many) make up at least 10% of all the observations in the
#'   table. [beta_div_rarefied()] and [alpha_div_rarefied()] choose their
#'   default depth the same way.
#' 
#' @section Coverage-Based Rarefaction:
#'   Chao & Jost (2012) compare samples at equal **coverage** (completeness:
#'   the share of the community's individuals that belong to the features
#'   seen) rather than at equal depth. With `coverage`, each sample is
#'   rarefied to the smallest depth \eqn{d} whose expected coverage
#'   \eqn{C_d} (as in [rarefaction_curve()]'s `coverage` column) is at least
#'   `coverage`. At the sample's own total, \eqn{C_d} is the estimated
#'   coverage of the whole sample (Chao et al. 2014, Eq. 12). Samples whose
#'   whole-sample coverage is below `coverage` are dropped (or, with `drop =
#'   FALSE`, returned unrarefied). The depths come from the counts alone, so
#'   they are the same for every `seed`; each sample is then rarefied as
#'   with one `depth` per sample. `rowSums()` of the result gives the depths
#'   chosen.
#'   
#'   Any metric can then be compared at equal coverage, not only richness or
#'   Hill numbers.
#' 
#' @references
#' Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation:
#' standardizing samples by completeness rather than size. *Ecology*, 93,
#' 2533--2547. \doi{10.1890/11-1952.1}
#' 
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
#' R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
#' numbers: a framework for sampling and estimation in species diversity
#' studies. *Ecological Monographs*, 84, 45--67. \doi{10.1890/13-0133.1}
#' 
#' @section One Depth per Sample:
#'   With one `depth` per sample, each sample is rarefied to its own depth,
#'   and gives exactly what rarefying the table to that sample's depth would
#'   give it, with the same `seed`. `drop` and `warn` compare each sample with
#'   its own depth.
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
#'     # One depth per sample, matched by name.
#'     r_mtx <- rarefy(counts, depth = c(A = 10, B = 15, C = 20, D = 12))
#'     rowSums(r_mtx)
#'     
#'     # Rarefy each sample to an expected coverage of 80%.
#'     r_mtx <- rarefy(counts, coverage = 0.8)
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
    depth    = NULL, 
    seed     = 0L,
    times    = NULL, 
    drop     = TRUE, 
    margin   = 1L, 
    cpus     = n_cpus(),
    warn     = interactive(),
    coverage = NULL ) {
  
  validate_args()
  assert_integer_counts()
  
  if (!is.null(coverage) && !is.null(depth))
    stop('Give `depth` or `coverage`, not both.', call. = FALSE)
  
  
  # Ensure counts are double
  if (is.matrix(counts)) {
    if (!is.double(counts)) counts[] <- as.double(counts)
  }
  else if (inherits(counts, "simple_triplet_matrix")) {
    if (!is.double(counts$v)) counts$v <- as.double(counts$v)
  }
  
  
  # Coverage-based rarefaction (Chao & Jost 2012): each sample's depth is
  # the smallest whose expected coverage reaches `coverage`, or one more
  # than its total if even the whole sample's estimated coverage falls short
  # (so that it is dropped, or returned unrarefied, as for a depth).
  if (!is.null(coverage)) {
    if (any(stored_values(counts) < 0))
      stop('`counts` must not contain negative values.', call. = FALSE)
    depth <- .Call(C_coverage_depth, counts, margin, coverage, cpus)
    names(depth) <- count_samples(counts, margin)$labels
  }
  
  
  # One depth per sample: matched to the samples by name when both are
  # named, else in order. A single depth applies to every sample.
  else if (length(depth) > 1L) {
    info  <- count_samples(counts, margin)
    depth <- per_sample(
      depth, 'depth', NULL, NULL, info$labels, info$n, ref = '`counts`' )
    names(depth) <- info$labels
  }
  
  
  # Calculate library sizes first, as they are needed for both
  # auto-selection AND the warning check.
  if (is.null(depth) || isTRUE(warn)) {
    
    sums <- sample_sums(counts, margin)
    
    if (is.null(depth)) depth <- auto_depth(sums)
    
    
    # Warning logic:
    # Check if any samples have insufficient depth and warn if requested.
    if (isTRUE(warn)) {
      
      n_insufficient <- sum(sums < depth)
      
      if (n_insufficient > 0) {
        
        action <- if (isTRUE(drop)) "dropped" else "returned unrarefied"
        
        if (!is.null(coverage))
          warning(sprintf(
            "%d samples have an estimated coverage below %s and will be %s.",
            n_insufficient, format(coverage), action
          ), call. = FALSE)
        else if (length(depth) == 1L)
          warning(sprintf(
            "%d samples have fewer than %d observations and will be %s.",
            n_insufficient, depth, action
          ), call. = FALSE)
        else
          warning(sprintf(
            "%d samples have fewer observations than their `depth` and will be %s.",
            n_insufficient, action
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
        
        # We must keep samples where sum >= depth (each sample's own depth)
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



# Each sample's total count (library size): rows with margin = 1, columns
# with margin = 2, for any table class rarefy() takes.
sample_sums <- function (counts, margin) {
  if (is.matrix(counts)) {
    if (margin == 1L) rowSums(counts) else colSums(counts)
  }
  else if (inherits(counts, "simple_triplet_matrix")) {
    if (margin == 1L) slam::row_sums(counts) else slam::col_sums(counts)
  }
  else {
    # Matrix package objects
    if (margin == 1L) Matrix::rowSums(counts) else Matrix::colSums(counts)
  }
}


# The automatic rarefaction depth for library sizes `sums`: the smallest
# library size whose yield (that depth times the number of samples with at
# least that many observations) is at least 10% of all the observations.
auto_depth <- function (sums) {
  
  sums_sorted <- sort(sums)
  n_sams      <- length(sums_sorted)
  total_reads <- sum(sums_sorted)
  threshold   <- 0.1 * total_reads
  
  # Vectorized calculation:
  yields <- sums_sorted * seq.int(n_sams, 1L)
  
  # Find the first index (lowest depth) that meets the threshold
  idx <- which(yields >= threshold)[1L]
  
  if (!is.na(idx)) {
    as.integer(sums_sorted[idx])
  } else {
    # Fallback: Use the max depth (retains 1 sample, but 100% of its reads).
    # To reach the line, we need a distribution where NO depth meets the 10% 
    # retention threshold, such as a Zipfian distribution (1/rank) 
    # with N > 23,000 samples.
    as.integer(max(sums_sorted)) # nocov
  }
}
