test_that("rarefaction", {
  
  expect_true(is.matrix(rarefy(counts)))
  expect_true(is.list(rarefy(counts, times = 2)))
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 1)), 
    expected = as.integer(c(0,0,5,2,6,0,4,4,0,5,0,3,5,0,5,0,10,3,0,0)) )
  
  expect_identical(
    object   = as.vector(rarefy(counts, seed = 2)), 
    expected = as.integer(c(0,0,5,2,6,0,7,2,0,4,0,5,2,0,6,0,9,4,0,0)) )
  
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
  expect_in(colSums(m), range(colSums(m)))
  
  expect_error(rarefy(counts * 1.5))
})
