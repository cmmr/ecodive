test_that("validation", {
  
  # validate_pairs() ==========================================================
  
  env <- new.env()
  
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
  
  env$pairs <- c(F, F, T, T, T, F)
  expect_silent(validate_pairs(env))
  
  
  
  
  # validate_alpha() ==========================================================
  
  env$alpha <- 1L
  validate_alpha(env)
  
  
  
  
  # validate_tree() ===========================================================
  
  env$counts <- counts[1:3,,drop=FALSE]
  tree2 <- tree
  tree2$edge.length <- as.integer(tree$edge.length * 100)
  tree2$edge <- matrix(
    data     = as.numeric(tree$edge), 
    nrow     = nrow(tree$edge), 
    ncol     = ncol(tree$edge),
    dimnames = dimnames(tree$edge) )
  env$tree <- tree2
  validate_tree(env)
  
  
  
  
  # validate_counts() =========================================================
  
  env$tree   <- NULL
  env$counts <- counts
  attr(env$counts, 'tree') <- tree
  expect_silent(validate_counts(env))
  
  skip_on_cran()
  skip_if_not_installed('rbiom')
  
  hmp50               <- do.call(`::`, list('rbiom', 'hmp50'))
  convert_to_phyloseq <- do.call(`::`, list('rbiom', 'convert_to_phyloseq'))
  convert_to_TSE      <- do.call(`::`, list('rbiom', 'convert_to_TSE'))
  convert_to_SE       <- do.call(`::`, list('rbiom', 'convert_to_SE'))
  
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
  
  
})
