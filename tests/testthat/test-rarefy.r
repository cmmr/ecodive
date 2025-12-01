test_that("rarefaction", {
  
  # Basic return types
  
  expect_true(is.matrix(rarefy(counts)))
  expect_true(is.list(rarefy(counts, times = 2)))
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 1, drop = FALSE)), 
    expected = c(0,0,0,0,0,4,3,10,5,4,5,3,2,0,0,0,6,5,5,0) )
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 2, drop = FALSE)), 
    expected = c(0,0,0,0,0,7,5,9,5,2,2,4,2,0,0,0,6,4,6,0) )
  
  expect_identical(
    object   = rarefy(counts, n_samples = 0), 
    expected = rarefy(counts, n_samples = 4) )
  
  expect_identical(
    object   = rarefy(counts, n_samples = 0.5), 
    expected = rarefy(counts, n_samples = 2) )
  
  expect_identical(
    object   = rarefy(counts, n_samples = 3), 
    expected = rarefy(counts, n_samples = -1) )
  
  expect_identical(
    object   = rarefy(counts, depth = 0.1), 
    expected = rarefy(counts, depth = 13) )
  
  expect_silent(m <- rarefy(big_mtx, 15))
  expect_in(rowSums(m), range(rowSums(m)))
  
  # Test times parameter
  r_list <- rarefy(counts, times = 3, seed = 42)
  expect_length(r_list, 3)
  expect_true(is.matrix(r_list[[1]]))
  expect_false(identical(r_list[[1]], r_list[[2]]))
  
  # Test drop parameter
  
  # Rarefy to a high depth to create zero-sum rows/cols
  r_drop_F <- rarefy(counts, depth = 15, drop = FALSE, seed = 1)
  expect_equal(dim(r_drop_F), c(4, 5))
  expect_equal(unname(rowSums(r_drop_F)['A']), 0)
  
  r_drop_T <- rarefy(counts, depth = 15, drop = TRUE, seed = 1)
  expect_equal(dim(r_drop_T), c(3, 3))
  expect_false("A" %in% rownames(r_drop_T))
  expect_false("OTU1" %in% colnames(r_drop_T))
  
  # Test drop with times > 1
  r_list_drop <- rarefy(counts, depth = 15, times = 2, drop = TRUE, seed = 1)
  expect_length(r_list_drop, 2)
  expect_equal(dim(r_list_drop[[1]]), c(3, 3))
  expect_equal(dim(r_list_drop[[2]]), c(3, 3))
  expect_false("A" %in% rownames(r_list_drop[[1]]))
  
  # Test error on non-integer counts
  expect_error(rarefy(counts * 1.5))
  
  # Ensure Matrix and slam packages are available for testing
  skip_if_not_installed('Matrix')
  skip_if_not_installed('slam')
  
  # Test error on non-integer counts
  expect_error(rarefy(as(counts * 1.5, 'dgCMatrix')))
  expect_error(rarefy(slam::as.simple_triplet_matrix(counts * 1.5)))

  # Test with different matrix types
  counts_dgC  <- as(counts, 'dgCMatrix')
  counts_dgT  <- as(counts, 'TsparseMatrix')
  counts_dge  <- as(counts, 'unpackedMatrix')
  counts_slam <- slam::as.simple_triplet_matrix(counts)
  
  # rarefy with different matrix types
  expect_silent(rarefy(counts,      seed = 1))
  expect_silent(rarefy(counts_dgC,  seed = 1))
  expect_silent(rarefy(counts_dgT,  seed = 1))
  expect_silent(rarefy(counts_dge,  seed = 1))
  expect_silent(rarefy(counts_slam, seed = 1))
  
  # Test margin = 2 (samples in columns)
  counts_t  <- t(counts)
  r_t_depth <- rarefy(counts_t, depth = 14, margin = 2L, seed = 1)
  expect_equal(colSums(r_t_depth), c(B=14, C=14, D=14))
  
  # Test with n_samples
  r_t_nsamples <- rarefy(counts_t, n_samples = 2, margin = 2L, seed = 1)
  expect_equal(ncol(r_t_nsamples), 2)
  expect_equal(sum(colSums(r_t_nsamples) > 0), 2)
  
  # Test different matrix types with margin = 2
  counts_t_dgC  <- as(counts_t, 'dgCMatrix')
  counts_t_dgT  <- as(counts_t, 'TsparseMatrix')
  counts_t_dge  <- as(counts_t, 'unpackedMatrix')
  counts_t_slam <- slam::as.simple_triplet_matrix(counts_t)
  
  expect_silent(rarefy(counts_t,      margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dgC,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dgT,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_dge,  margin = 2L, seed = 1))
  expect_silent(rarefy(counts_t_slam, margin = 2L, seed = 1))
  
  # Test that samples with depth < target are zeroed out
  # This covers the `else if (depth < target)` branches in C
  r_depth_20 <- rarefy(counts, depth = 20, drop = FALSE)
  expect_equal(unname(rowSums(r_depth_20)), c(0, 0, 20, 0))
  
  # Test for triplet matrices (dgTMatrix, slam, dgCMatrix margin=1)
  # which use rarefy_triplet()
  r_dgT_20 <- rarefy(counts_dgT, depth = 20, drop = FALSE)
  expect_s4_class(r_dgT_20, "dgTMatrix")
  expect_equal(unname(Matrix::rowSums(r_dgT_20)), c(0, 0, 20, 0))
  
  r_slam_20 <- rarefy(counts_slam, depth = 20, drop = FALSE)
  expect_s3_class(r_slam_20, "simple_triplet_matrix")
  expect_equal(unname(slam::row_sums(r_slam_20)), c(0, 0, 20, 0))
  
  # Test for compressed sparse matrix (dgCMatrix margin=2)
  # which uses rarefy_compressed()
  r_t_dgC_20 <- rarefy(counts_t_dgC, depth = 20, margin = 2L, drop = FALSE)
  expect_s4_class(r_t_dgC_20, "dgCMatrix")
  expected_t_20 <- t(rarefy(counts, depth = 20, drop = FALSE))
  expect_equal(unname(Matrix::colSums(r_t_dgC_20)), c(0, 0, 20, 0))
  
})
