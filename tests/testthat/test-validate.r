test_that("validation", {
  
  env <- new.env()
  
  
  
  
  # assert_integer_counts() ====
  
  env$counts <- counts_p
  expect_error(assert_integer_counts(env))



  # validate_pairs() ====

  env$counts <- matrix(nrow = 5, ncol = 0)
  expect_error(validate_pairs(env))

  env$counts <- counts

  env$pairs <- function (i,j) list()
  expect_error(validate_pairs(env))

  env$pairs <- TRUE; expect_error(validate_pairs(env))
  env$pairs <- 1.5;  expect_error(validate_pairs(env))
  env$pairs <- -1;   expect_error(validate_pairs(env))
  env$pairs <- 100;  expect_error(validate_pairs(env))

  env$pairs <- 2:5
  expect_silent(validate_pairs(env))

  env$pairs <- c(F,F,T,T,T,F,F,F,F,F)
  expect_silent(validate_pairs(env))




  # validate_alpha() ====

  env$alpha <- 1L
  expect_silent(validate_alpha(env))




  # validate_cpus() ====

  env$cpus <- 1
  expect_silent(validate_cpus(env))




  # validate_power() ====

  env$power <- 1L
  expect_silent(validate_power(env))



  
  # validate_norm() ====

  env$norm <- 'badoption'
  expect_error(validate_norm(env))

  env$norm <- 'chord'
  expect_silent(validate_norm(env))

  env$counts <- counts
  env$norm   <- 'percent'
  expect_silent(validate_norm(env))

  env$counts <- counts
  env$norm   <- 'binary'
  expect_silent(validate_norm(env))

  env$counts <- counts
  env$norm   <- 'clr'
  expect_silent(validate_norm(env))

  env$counts <- counts
  env$norm   <- 'none'
  expect_silent(validate_norm(env))



  # validate_tree() ====

  env$counts <- counts[1:3,,drop=FALSE]
  env$margin <- 1L
  tree2      <- tree
  tree2$edge.length <- as.integer(tree$edge.length * 100)
  tree2$edge <- matrix(
    data     = as.numeric(tree$edge),
    nrow     = nrow(tree$edge),
    ncol     = ncol(tree$edge),
    dimnames = dimnames(tree$edge) )
  env$tree <- tree2
  expect_silent(validate_tree(env))
  
  env$counts <- t(counts[1:3,,drop=FALSE])
  env$margin <- 2L
  expect_silent(validate_tree(env))

  
  env$counts <- counts[,1:4,drop=FALSE]
  env$margin <- 1L
  env$tree   <- tree
  expect_silent(validate_tree(env))
  
  env$counts <- t(counts[,1:4,drop=FALSE])
  env$margin <- 2L
  expect_silent(validate_tree(env))

  



  # validate_counts() ====

  env$tree   <- NULL
  env$counts <- counts
  env$margin <- 1L
  attr(env$counts, 'tree') <- tree
  expect_silent(validate_counts(env))

  skip_on_cran()
  skip_if_not_installed('rbiom')
  skip_if_not_installed('phyloseq')
  skip_if_not_installed('TreeSummarizedExperiment')

  hmp50               <- do.call(`::`, list('rbiom', 'hmp50'))
  convert_to_phyloseq <- do.call(`::`, list('rbiom', 'convert_to_phyloseq'))
  convert_to_TSE      <- do.call(`::`, list('rbiom', 'convert_to_TSE'))
  convert_to_SE       <- do.call(`::`, list('rbiom', 'convert_to_SE'))

  env$tree   <- NULL
  env$counts <- 1:10
  expect_silent(validate_counts(env))

  env$tree   <- NULL
  env$counts <- Matrix::t(hmp50$counts)
  expect_silent(validate_counts(env))

  env$tree   <- NULL
  env$counts <- hmp50
  expect_silent(validate_counts(env))

  env$tree   <- NULL
  env$counts <- convert_to_phyloseq(hmp50)
  expect_silent(validate_counts(env))

  env$tree   <- NULL
  env$counts <- convert_to_TSE(hmp50)
  expect_silent(validate_counts(env))

  env$tree   <- NULL
  env$counts <- convert_to_SE(hmp50)
  expect_silent(validate_counts(env))
  
  
  
  
  # validate_pseudocount() ====
  
  # Mock constants used in validate.r
  NORM_PERCENT <- 1L
  NORM_CLR     <- 2L
  
  # 1. Ignore when norm is NOT CLR
  # ----------------------------------------------------------------
  env$norm        <- NORM_PERCENT
  env$pseudocount <- 999
  env$counts      <- matrix(1:4, 2, 2)
  expect_silent(validate_pseudocount(env))
  expect_null(env$pseudocount)
  
  
  # 2. Explicit Pseudocount (Valid)
  # ----------------------------------------------------------------
  env$norm        <- NORM_CLR
  env$pseudocount <- 0.1
  expect_silent(validate_pseudocount(env))
  expect_equal(env$pseudocount, 0.1)
  
  
  # 3. Explicit Pseudocount (Invalid inputs)
  # ----------------------------------------------------------------
  env$norm <- NORM_CLR
  
  env$pseudocount <- -1
  expect_error(validate_pseudocount(env), "pseudocount` must be a single positive number")
  
  env$pseudocount <- 0
  expect_error(validate_pseudocount(env), "pseudocount` must be a single positive number")
  
  env$pseudocount <- c(1, 2)
  expect_error(validate_pseudocount(env), "pseudocount` must be a single positive number")
  
  env$pseudocount <- NA_real_
  expect_error(validate_pseudocount(env), "pseudocount` must be a single positive number")
  
  env$pseudocount <- "string"
  expect_error(validate_pseudocount(env), "pseudocount` must be a single positive number")
  
  
  # 4. Data Validation: Negative / Infinite / NA values in counts
  # ----------------------------------------------------------------
  env$norm        <- NORM_CLR
  env$pseudocount <- NULL
  
  env$counts <- matrix(c(1, 2, -1, 4), 2, 2)
  expect_error(validate_pseudocount(env), "contains negative values")
  
  env$counts <- matrix(c(1, 2, Inf, 4), 2, 2)
  expect_error(validate_pseudocount(env), "contains non-finite values")
  
  env$counts <- matrix(c(1, 2, NA, 4), 2, 2)
  expect_error(validate_pseudocount(env), "contains non-finite values")
  
  
  # 5. Auto-detection: Dense Matrix
  # ----------------------------------------------------------------
  env$norm <- NORM_CLR
  
  # Case: No zeros -> No change
  env$counts      <- matrix(1:4, 2, 2)
  env$pseudocount <- NULL
  expect_silent(validate_pseudocount(env))
  expect_null(env$pseudocount)
  
  # Case: Zeros present -> Warning + Default (min/2)
  env$counts      <- matrix(c(0, 10, 20, 30), 2, 2) # min non-zero is 10
  env$pseudocount <- NULL
  expect_warning(validate_pseudocount(env), "Zeros detected")
  expect_equal(env$pseudocount, 5) # 10 / 2
  
  # Case: All zeros -> Error
  env$counts      <- matrix(0, 2, 2)
  env$pseudocount <- NULL
  expect_error(validate_pseudocount(env), "`counts` contains only zeros")
  
  # Case: All same value -> Error
  env$counts      <- matrix(1, 2, 2)
  env$pseudocount <- NULL
  expect_error(validate_pseudocount(env), "`counts` contains a single unique value")
  
  
  # 6. Auto-detection: Sparse Matrix (slam::simple_triplet_matrix)
  # ----------------------------------------------------------------
  skip_if_not_installed('slam')
  
  # Case: Implicit zeros (v has length < nrow*ncol)
  # Matrix: 10  0
  #          0 20
  stm <- slam::simple_triplet_matrix(i=c(1,2), j=c(1,2), v=c(10,20), nrow=2, ncol=2)
  
  env$counts      <- stm
  env$pseudocount <- NULL
  expect_warning(validate_pseudocount(env), "Zeros detected")
  expect_equal(env$pseudocount, 5) # 10 / 2
  
  # Case: Explicit zeros in v
  stm_explicit <- slam::simple_triplet_matrix(i=c(1,1,2,2), j=c(1,2,1,2), v=c(10,0,0,20), nrow=2, ncol=2)
  
  env$counts      <- stm_explicit
  env$pseudocount <- NULL
  expect_warning(validate_pseudocount(env), "Zeros detected")
  expect_equal(env$pseudocount, 5)
  
  
  # 7. Auto-detection: Sparse Matrix (Matrix::dgCMatrix)
  # ----------------------------------------------------------------
  skip_if_not_installed('Matrix')
  
  # Case: Implicit zeros
  dgc <- Matrix::Matrix(c(10, 0, 0, 20), 2, 2, sparse=TRUE)
  
  env$counts      <- dgc
  env$pseudocount <- NULL
  expect_warning(validate_pseudocount(env), "Zeros detected")
  expect_equal(env$pseudocount, 5) # 10 / 2
  
  # Case: Fully dense sparse matrix (no zeros)
  dgc_full <- Matrix::Matrix(c(1:4), 2, 2, sparse=TRUE)
  
  env$counts      <- dgc_full
  env$pseudocount <- NULL
  expect_silent(validate_pseudocount(env))
  expect_null(env$pseudocount)
  
  
})
