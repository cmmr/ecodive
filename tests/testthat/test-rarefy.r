test_that("rarefaction", {
  
  # Basic return types
  expect_true(is.matrix(rarefy(counts)))
  expect_true(is.list(rarefy(counts, times = 2)))
  
  # Check reproducibility (seed) and default behavior
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 1, drop = FALSE, warn = FALSE)), 
    expected = c(0,0,0,0,0,4,3,10,5,4,5,3,2,0,0,0,6,5,5,0) )
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 2, drop = FALSE, warn = FALSE)), 
    expected = c(0,0,0,0,0,7,5,9,5,2,2,4,2,0,0,0,6,4,6,0) )
  
  # Test Auto-Depth Selection 
  # Default auto-selection (depth=13) retains all samples in this dataset 
  # (A=13, B=18, C=21, D=15), so NO warning is issued.
  expect_identical(
    object   = rarefy(counts, depth = NULL, warn = TRUE), 
    expected = rarefy(counts, depth = 13, warn = TRUE) 
  )
  
  # Test Auto-Depth Fallback (Branch coverage for validate_depth fallback)
  # Create a matrix where 10% yield is impossible unless we pick max depth
  # e.g., one huge sample and many tiny ones.
  # This SHOULD warn because the tiny samples (1, 1) are < 100.
  tiny_counts <- matrix(c(100, 1, 1), ncol=1)
  expect_warning(
    res <- rarefy(tiny_counts, depth = NULL, warn = TRUE), 
    regexp = "dropped"
  )
  expect_equal(sum(res), 100) # Should pick max depth (100) and drop others
  
  
  expect_silent(m <- rarefy(big_mtx, 15, warn = FALSE))
  expect_in(rowSums(m), range(rowSums(m)))
  
  # Test times parameter
  r_list <- rarefy(counts, times = 3, seed = 42)
  expect_length(r_list, 3)
  expect_true(is.matrix(r_list[[1]]))
  expect_false(identical(r_list[[1]], r_list[[2]]))
  
  # -----------------------------------------------------------------------
  # Test drop parameter & Warning Logic
  # -----------------------------------------------------------------------
  
  # Case 1: drop = FALSE (Retain Original Counts)
  # Rarefy to depth 15. Sample A (13) < 15.
  # They should be returned UNMODIFIED.
  expect_warning(
    r_drop_F <- rarefy(counts, depth = 15, drop = FALSE, seed = 1, warn = TRUE),
    regexp = "returned unrarefied"
  )
  
  expect_equal(dim(r_drop_F), c(4, 5))
  expect_equal(unname(rowSums(r_drop_F)['A']), 13) # Original count (13)
  expect_equal(unname(rowSums(r_drop_F)['D']), 15) # Rarefied to 15 (Exact match)
  expect_equal(unname(rowSums(r_drop_F)['B']), 15) # Rarefied to 15
  
  # Case 2: drop = TRUE (Remove Samples)
  # A (13) < 15. B, C, D >= 15. Only A drops.
  expect_warning(
    r_drop_T <- rarefy(counts, depth = 15, drop = TRUE, seed = 1, warn = TRUE),
    regexp = "will be dropped"
  )
  expect_equal(dim(r_drop_T), c(3, 5)) # B, C, D remain
  expect_false("A" %in% rownames(r_drop_T))
  expect_true("D" %in% rownames(r_drop_T))
  
  # Case 3: warn = FALSE (Silent operation)
  expect_silent(rarefy(counts, depth = 15, drop = TRUE, warn = FALSE))
  
  
  # Test drop with times > 1
  expect_warning(
    r_list_drop <- rarefy(counts, depth = 15, times = 2, drop = TRUE, seed = 1, warn = TRUE),
    regexp = "will be dropped"
  )
  expect_length(r_list_drop, 2)
  expect_equal(dim(r_list_drop[[1]]), c(3, 5))
  expect_false("A" %in% rownames(r_list_drop[[1]]))
  
  # Test error on non-integer counts
  expect_error(rarefy(counts * 1.5))
  
  
  # -----------------------------------------------------------------------
  # Matrix Types, Margin Logic & C Code Coverage
  # -----------------------------------------------------------------------
  
  # Ensure Matrix and slam packages are available for testing
  skip_if_not_installed('Matrix')
  skip_if_not_installed('slam')
  
  if (packageVersion("Matrix") >= "1.5") {
    COMPRESSED <- "sparseMatrix"
    TRIPLET    <- "TsparseMatrix"
    UNPACKED   <- "unpackedMatrix"
  } else {
    COMPRESSED <- "dgCMatrix"
    TRIPLET    <- "dgTMatrix"
    UNPACKED   <- "dgeMatrix"
  }
  
  # Test error on non-integer counts
  expect_error(rarefy(as(counts * 1.5, COMPRESSED)))
  expect_error(rarefy(slam::as.simple_triplet_matrix(counts * 1.5)))
  
  # Test with different matrix types
  counts_dgC  <- as(counts, COMPRESSED)
  counts_dgT  <- as(counts, TRIPLET)
  counts_dge  <- as(counts, UNPACKED)
  counts_slam <- slam::as.simple_triplet_matrix(counts)
  
  # rarefy with different matrix types
  expect_silent(rarefy(counts,      seed = 1))
  expect_silent(rarefy(counts_dgC,  seed = 1))
  expect_silent(rarefy(counts_dgT,  seed = 1))
  expect_silent(rarefy(counts_dge,  seed = 1))
  expect_silent(rarefy(counts_slam, seed = 1))
  
  # Test with integer class data
  counts_slam_int   <- counts_slam
  counts_slam_int$v <- as.integer(counts_slam_int$v)
  expect_equal(
    rarefy(counts_int, seed = 1),
    rarefy(counts,     seed = 1)
  )
  expect_equal(
    rarefy(counts_slam_int, seed = 1),
    rarefy(counts_slam,     seed = 1)
  )
  
  # Test margin = 2 (samples in columns)
  counts_t  <- t(counts)
  r_t_depth <- rarefy(counts_t, depth = 14, margin = 2L, seed = 1)
  expect_equal(colSums(r_t_depth), c(B=14, C=14, D=14)) # A (13) dropped implicit via helper.r values
  
  # Test Margin 2 Dropper Logic 
  # Sums: A=13, B=18, C=21, D=15. Depth 15.
  r_t_drop <- rarefy(counts_t, depth = 15, margin = 2L, drop = TRUE, warn = FALSE)
  expect_equal(ncol(r_t_drop), 3) 
  expect_false("A" %in% colnames(r_t_drop))
  expect_true("D" %in% colnames(r_t_drop))
  
  # Test different matrix types with margin = 2
  counts_t_dgC  <- as(counts_t, COMPRESSED)
  counts_t_dgT  <- as(counts_t, TRIPLET)
  counts_t_dge  <- as(counts_t, UNPACKED)
  counts_t_slam <- slam::as.simple_triplet_matrix(counts_t)
  
  expect_silent(rarefy(counts_t,      margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dgC,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dgT,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dge,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_slam, margin = 2L, seed = 1))
  
  
  # -----------------------------------------------------------------------
  # C Code Branch Coverage: Original Counts Retention
  # -----------------------------------------------------------------------
  
  # Test that samples with depth < target return ORIGINAL counts (not zeroed)
  r_depth_20 <- rarefy(counts, depth = 20, drop = FALSE, warn = FALSE)
  expect_equal(unname(rowSums(r_depth_20)), c(13, 18, 20, 15)) 
  
  # Test for triplet matrices (dgTMatrix, slam, dgCMatrix margin=1)
  r_dgT_20 <- rarefy(counts_dgT, depth = 20, drop = FALSE, warn = FALSE)
  expect_s4_class(r_dgT_20, TRIPLET)
  expect_equal(unname(Matrix::rowSums(r_dgT_20)), c(13, 18, 20, 15))
  
  r_slam_20 <- rarefy(counts_slam, depth = 20, drop = FALSE, warn = FALSE)
  expect_s3_class(r_slam_20, "simple_triplet_matrix")
  expect_equal(unname(slam::row_sums(r_slam_20)), c(13, 18, 20, 15))
  
  # Test for compressed sparse matrix (dgCMatrix margin=2)
  r_t_dgC_20 <- rarefy(counts_t_dgC, depth = 20, margin = 2L, drop = FALSE, warn = FALSE)
  expect_s4_class(r_t_dgC_20, COMPRESSED)
  expect_equal(unname(Matrix::colSums(r_t_dgC_20)), c(13, 18, 20, 15))
  
  
  # -----------------------------------------------------------------------
  # C Code Branch Coverage: Sparse Compaction & Multithreading
  # -----------------------------------------------------------------------
  
  # To trigger the "Allocate new vectors" block in compact_* functions,
  # we must ensure rarefaction turns some non-zero values into zeros.
  # Using depth=1 on a sample with multiple 1s guarantees this.
  
  s_counts <- matrix(c(1, 1, 1), nrow=1, dimnames=list("S1", c("O1", "O2", "O3")))
  
  # 1. dgCMatrix Compaction
  s_dgC <- as(s_counts, COMPRESSED)
  r_dgC <- rarefy(s_dgC, depth=1, seed=1, warn=FALSE)
  expect_s4_class(r_dgC, COMPRESSED)
  expect_equal(sum(r_dgC), 1)
  expect_equal(length(r_dgC@x), 1) # Explicit zeros should be gone
  
  # 2. dgTMatrix Compaction
  s_dgT <- as(s_counts, TRIPLET)
  r_dgT <- rarefy(s_dgT, depth=1, seed=1, warn=FALSE)
  expect_s4_class(r_dgT, TRIPLET)
  expect_equal(sum(r_dgT), 1)
  expect_equal(length(r_dgT@x), 1)
  
  # 3. simple_triplet_matrix Compaction
  s_slam <- slam::as.simple_triplet_matrix(s_counts)
  r_slam <- rarefy(s_slam, depth=1, seed=1, warn=FALSE)
  expect_s3_class(r_slam, "simple_triplet_matrix")
  expect_equal(sum(r_slam$v), 1)
  expect_equal(length(r_slam$v), 1)
  
  # -----------------------------------------------------------------------
  # Multithreading Coverage
  # -----------------------------------------------------------------------
  # Requires n_sams > 100 and cpus > 1 to trigger threaded blocks.
  # We use big_mtx (104 samples).
  # We use depth=1 to trigger zero-creation (compaction) inside threaded blocks.
  
  # Dense Threaded
  # NOTE: big_mtx contains 26 copies of 'counts'. Row A has 13 observations.
  # Using depth=15 would normally drop those rows if drop=TRUE (default).
  # We set drop=FALSE here to ensure we get back a full-sized matrix, 
  # verifying that the C threading code processed all rows correctly.
  expect_silent(r_big <- rarefy(big_mtx, depth=15, cpus=2, drop=FALSE, warn=FALSE))
  expect_equal(dim(r_big), dim(big_mtx))
  
  # Sparse Threaded + Compaction
  # dgCMatrix (margin 1, default) covers setup_triplet threading
  big_dgC <- as(big_mtx, COMPRESSED)
  r_big_dgC <- rarefy(big_dgC, depth=1, cpus=2, warn=FALSE)
  expect_equal(dim(r_big_dgC), dim(big_mtx))
  expect_true(length(r_big_dgC@x) < length(big_dgC@x)) # Zeros were removed
  
  # dgCMatrix (margin 2) covers rarefy_compressed threading
  # big_mtx is 104x5. Transpose to 5x104 to get >100 samples in columns.
  big_mtx_t <- t(big_mtx)
  big_dgC_t <- as(big_mtx_t, COMPRESSED)
  r_big_dgC_t <- rarefy(big_dgC_t, depth=1, margin=2, cpus=2, warn=FALSE)
  expect_equal(dim(r_big_dgC_t), dim(big_mtx_t))
  expect_true(length(r_big_dgC_t@x) < length(big_dgC_t@x))
  
  # dgTMatrix Threading
  big_dgT <- as(big_mtx, TRIPLET)
  r_big_dgT <- rarefy(big_dgT, depth=1, cpus=2, warn=FALSE)
  expect_equal(dim(r_big_dgT), dim(big_mtx))
  expect_true(length(r_big_dgT@x) < length(big_dgT@x))
  
  # slam Threading
  big_slam <- slam::as.simple_triplet_matrix(big_mtx)
  r_big_slam <- rarefy(big_slam, depth=1, cpus=2, warn=FALSE)
  expect_equal(dim(r_big_slam), dim(big_mtx))
  expect_true(length(r_big_slam$v) < length(big_slam$v))
  
})
