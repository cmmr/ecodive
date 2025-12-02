test_that("ecomatrix.c parsing logic is covered", {

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

  # Define a test function that calls an internal C function
  # We will test an alpha diversity metric which uses the ecomatrix parsing
  test_parsing <- function(m, margin = 1L) {
    observed(m, margin = margin)
  }

  # Expected result: number of non-zero entries per sample
  expected_margin1 <- c(A = 3, B = 3, C = 3, D = 2)
  expected_margin2 <- c(OTU1 = 2, OTU2 = 4, OTU3 = 4, OTU4 = 3, OTU5 = 2)

  # === Test with base R matrix ===
  m_base <- counts
  expect_equal(test_parsing(m_base, margin = 1L), expected_margin1)
  expect_equal(test_parsing(t(m_base), margin = 2L), expected_margin1)

  # Test with logical matrix
  m_logical <- m_base > 0
  expect_equal(test_parsing(m_logical, margin = 1L), expected_margin1)
  expect_equal(test_parsing(t(m_logical), margin = 2L), expected_margin1)

  # Test with double matrix
  m_double <- m_base * 1.1
  expect_equal(test_parsing(m_double, margin = 1L), expected_margin1)
  expect_equal(test_parsing(t(m_double), margin = 2L), expected_margin1)


  # === Test with Matrix package types ===
  m_base_t <- t(m_base)
  m_dgC    <- as(m_base,   COMPRESSED)
  m_dgC_t  <- as(m_base_t, COMPRESSED)
  m_dgT    <- as(m_base,   TRIPLET)
  m_dgT_t  <- as(m_base_t, TRIPLET)
  m_dge    <- as(m_base,   UNPACKED)
  m_dge_t  <- as(m_base_t, UNPACKED)

  # dgCMatrix
  expect_equal(test_parsing(m_dgC, margin = 1L), expected_margin1, info = "dgCMatrix margin 1")
  expect_equal(test_parsing(m_dgC_t, margin = 2L), expected_margin1, info = "dgCMatrix margin 2")

  # dgTMatrix
  expect_equal(test_parsing(m_dgT, margin = 1L), expected_margin1, info = "dgTMatrix margin 1")
  expect_equal(test_parsing(m_dgT_t, margin = 2L), expected_margin1, info = "dgTMatrix margin 2")

  # dgeMatrix
  expect_equal(test_parsing(m_dge, margin = 1L), expected_margin1, info = "dgeMatrix margin 1")
  expect_equal(test_parsing(m_dge_t, margin = 2L), expected_margin1, info = "dgeMatrix margin 2")

  # Test with logical Matrix types
  m_logical_t <- t(m_logical)
  m_lgC       <- as(m_logical,   COMPRESSED)
  m_lgC_t     <- as(m_logical_t, COMPRESSED)
  m_lgT       <- as(m_logical,   TRIPLET)
  m_lge       <- as(m_logical,   UNPACKED)
  expect_equal(test_parsing(m_lgC, margin = 1L), expected_margin1, info = "lgCMatrix margin 1")
  expect_equal(test_parsing(m_lgC_t, margin = 2L), expected_margin1, info = "lgCMatrix margin 2")
  expect_equal(test_parsing(m_lgT, margin = 1L), expected_margin1, info = "lgTMatrix margin 1")
  expect_equal(test_parsing(m_lge, margin = 1L), expected_margin1, info = "lgeMatrix margin 1")


  # === Test with slam::simple_triplet_matrix ===
  m_slam   <- slam::as.simple_triplet_matrix(m_base)
  m_slam_t <- slam::as.simple_triplet_matrix(m_base_t)
  expect_equal(test_parsing(m_slam, margin = 1L), expected_margin1, info = "slam margin 1")
  expect_equal(test_parsing(m_slam_t, margin = 2L), expected_margin1, info = "slam margin 2")

  # Test with logical slam
  m_slam_logical <- slam::as.simple_triplet_matrix(m_logical)
  expect_equal(test_parsing(m_slam_logical, margin = 1L), expected_margin1, info = "slam logical margin 1")

  # Test with double slam
  m_slam_double <- slam::as.simple_triplet_matrix(m_double)
  expect_equal(test_parsing(m_slam_double, margin = 1L), expected_margin1, info = "slam double margin 1")


  # === Test sorting logic for sparse matrices ===
  # Create an unsorted dgTMatrix
  unsorted_m_dgT <- m_dgT
  set.seed(1)
  idx <- sample(length(unsorted_m_dgT@i))
  unsorted_m_dgT@i <- unsorted_m_dgT@i[idx]
  unsorted_m_dgT@j <- unsorted_m_dgT@j[idx]
  unsorted_m_dgT@x <- unsorted_m_dgT@x[idx]
  expect_equal(test_parsing(unsorted_m_dgT, margin = 1L), expected_margin1, info = "unsorted dgTMatrix")

  # Create an unsorted simple_triplet_matrix
  unsorted_m_slam <- m_slam
  unsorted_m_slam$i <- unsorted_m_slam$i[idx]
  unsorted_m_slam$j <- unsorted_m_slam$j[idx]
  unsorted_m_slam$v <- unsorted_m_slam$v[idx]
  expect_equal(test_parsing(unsorted_m_slam, margin = 1L), expected_margin1, info = "unsorted slam")

  # Test already sorted triplet
  sorted_m_dgT   <- as(m_dgC, TRIPLET) # dgC is sorted by column, then row
  sorted_m_dgT_t <- as(m_dgC_t, TRIPLET)
  expect_equal(test_parsing(sorted_m_dgT_t, margin = 2L), expected_margin1, info = "sorted dgTMatrix")


  # === Test error handling ===
  expect_error(test_parsing("not a matrix"))
  expect_error(test_parsing(m_base, margin = 3L))

  # Test non-numeric input
  char_matrix <- matrix(as.character(counts), nrow = 4)
  rownames(char_matrix) <- rownames(counts)
  expect_error(test_parsing(char_matrix))

  # === Test dimnames ===
  m_no_dimnames <- counts
  dimnames(m_no_dimnames) <- NULL
  res_no_dimnames <- test_parsing(m_no_dimnames)
  expect_null(names(res_no_dimnames))
  expect_equal(unname(res_no_dimnames), unname(expected_margin1))

  m_slam_no_dimnames <- slam::as.simple_triplet_matrix(m_no_dimnames)
  res_slam_no_dimnames <- test_parsing(m_slam_no_dimnames)
  expect_null(names(res_slam_no_dimnames))
  expect_equal(unname(res_slam_no_dimnames), unname(expected_margin1))

})
