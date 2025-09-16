test_that("rarefaction", {
  
  expect_true(is.matrix(rarefy(counts)))
  expect_true(is.list(rarefy(counts, times = 2)))
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 1)), 
    expected = c(0,0,0,0,0,4,3,10,5,4,5,3,2,0,0,0,6,5,5,0) )
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 2)), 
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
  
  expect_error(rarefy(counts * 1.5))
})
