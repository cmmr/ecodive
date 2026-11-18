#test_that("ecomatrix.c parsing logic is covered", {

  # Ensure Matrix and slam packages are available for testing
  exit_if_not(requireNamespace("Matrix", quietly=TRUE))
  exit_if_not(requireNamespace("slam",   quietly=TRUE))
  #skip_if_not_installed('Matrix')
  #skip_if_not_installed('slam')
  
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
  
  # Test with integer matrix
  expect_equal(test_parsing(counts_int, margin = 1L), expected_margin1)
  expect_equal(test_parsing(t(counts_int), margin = 2L), expected_margin1)


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
  m_slam_logical   <- m_slam
  m_slam_logical$v <- as.logical(m_slam_logical$v)
  expect_equal(test_parsing(m_slam_logical, margin = 1L), expected_margin1, info = "slam logical margin 1")

  # Test with double slam
  m_slam_double   <- m_slam
  m_slam_double$v <- as.double(m_slam_double$v)
  expect_equal(test_parsing(m_slam_double, margin = 1L), expected_margin1, info = "slam double margin 1")
  
  # Test with integer slam
  m_slam_int   <- m_slam
  m_slam_int$v <- as.integer(m_slam_int$v)
  expect_equal(test_parsing(m_slam_int, margin = 1L), expected_margin1, info = "slam int margin 1")


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


  # === Every input format gives identical results ===
  # Matrix and slam emit triplets column-major, which needs only the
  # scatter by sample; shuffled triplets also need the sort by otu.
  rand_tree <- function (tips) {
    nodes <- paste0(tips, ':', round(runif(length(tips), 0.1, 1), 3))
    while (length(nodes) > 1) {
      k     <- sample(length(nodes), 2)
      nodes <- c(nodes[-k], sprintf('(%s,%s):%.3f', nodes[k[1]], nodes[k[2]], runif(1, 0.1, 1)))
    }
    read_tree(paste0(sub(':[0-9.]+$', '', nodes), ';'))
  }
  shuffle <- function (x, idx) {
    if (inherits(x, 'dgTMatrix')) { x@i <- x@i[idx]; x@j <- x@j[idx]; x@x <- x@x[idx] }
    else                          { x$i <- x$i[idx]; x$j <- x$j[idx]; x$v <- x$v[idx] }
    x
  }
  alpha_ids <- list_metrics('alpha', val = 'id')
  beta_ids  <- list_metrics('beta',  val = 'id')
  all_divs  <- function (x, margin, tree) list(
    alpha = lapply(alpha_ids, function (id) alpha_div(x, id, tree = tree, margin = margin, cpus = 1L)),
    beta  = lapply(beta_ids,  function (id) beta_div( x, id, tree = tree, margin = margin, cpus = 1L, pseudocount = 0.5)) )

  set.seed(3)
  for (dims in list(c(12, 30, 0.3), c(25, 8, 0.7))) {
    n <- dims[[1]]; p <- dims[[2]]
    m <- matrix(rpois(n * p, 20) * rbinom(n * p, 1, dims[[3]]), n, p,
                dimnames = list(paste0('S', 1:n), paste0('OTU', 1:p)))
    m[,1] <- m[,1] + 1 # no empty samples
    storage.mode(m) <- 'integer'
    tr    <- rand_tree(colnames(m))
    ref   <- all_divs(m, 1L, tr)

    dgT   <- as(as(m,    COMPRESSED), TRIPLET)
    dgT_t <- as(as(t(m), COMPRESSED), TRIPLET)
    stm   <- slam::as.simple_triplet_matrix(m)
    stm_t <- slam::as.simple_triplet_matrix(t(m))
    stm_d <- stm; stm_d$v <- stm_d$v * 1.0 # double-valued slam
    idx   <- sample(length(stm$v))

    inputs <- list(
      matrix_2      = list(t(m),                 2L),
      dgC_1         = list(as(m,    COMPRESSED), 1L),
      dgC_2         = list(as(t(m), COMPRESSED), 2L),
      dgT_1         = list(dgT,                  1L),
      dgT_2         = list(dgT_t,                2L),
      dgT_1_shuf    = list(shuffle(dgT,   idx),  1L),
      dgT_2_shuf    = list(shuffle(dgT_t, idx),  2L),
      slam_1        = list(stm,                  1L),
      slam_2        = list(stm_t,                2L),
      slam_1_dbl    = list(stm_d,                1L),
      slam_1_shuf   = list(shuffle(stm,   idx),  1L),
      slam_2_shuf   = list(shuffle(stm_t, idx),  2L),
      slam_1_dblshf = list(shuffle(stm_d, idx),  1L) )

    for (nm in names(inputs))
      expect_identical(
        all_divs(inputs[[nm]][[1]], inputs[[nm]][[2]], tr), ref,
        info = paste(nm, n, 'x', p) )
  }


  # === Out-of-range sparse indices are errors ===
  bad <- m_slam;  bad$i[1] <- 99L
  expect_error(test_parsing(bad, margin = 1L), 'out of range', info = "slam i")
  bad <- m_slam;  bad$j[1] <- 0L
  expect_error(test_parsing(bad, margin = 1L), 'out of range', info = "slam j")
  bad <- m_dgT;   bad@i[1] <- -1L
  expect_error(test_parsing(bad, margin = 1L), 'out of range', info = "dgT i")
  bad <- m_dgT;   bad@j[1] <- 5L
  expect_error(test_parsing(bad, margin = 1L), 'out of range', info = "dgT j")
  bad <- m_dgC;   bad@i[1] <- 4L
  expect_error(test_parsing(bad, margin = 1L), 'out of range', info = "dgC i")
  bad <- m_dgC;   bad@p[3] <- bad@p[4] + 1L
  expect_error(test_parsing(bad, margin = 1L), 'column pointers', info = "dgC p")


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

  # Regression: integer NA in a sparse triplet matrix was read as
  # -2147483648, not NA. Both import orders: by sample and by feature.
  m_na <- matrix(rpois(200, 2) + 1L, 20, 10)
  m_na[3, 2] <- NA
  m_na_dbl <- m_na
  storage.mode(m_na_dbl) <- 'double'
  for (mg in 1:2) {
    x <- if (mg == 1) m_na else t(m_na)
    s <- slam::as.simple_triplet_matrix(x)
    expect_true(is.integer(s$v))
    expect_identical(bray(s, margin = mg), bray(m_na_dbl))
  }


  # Explicitly stored zeros are dropped, and repeated entries in a triplet
  # matrix are summed, as the Matrix package defines them: the same results
  # as the dense matrix.
  set.seed(7)
  dn  <- matrix(rpois(600, 1), 60, 10)
  dgc <- as(dn, COMPRESSED)
  dgc@x[1:3] <- 0
  dz  <- as.matrix(dgc)
  for (mg in 1:2) {
    x <- if (mg == 1) dgc else Matrix::t(dgc)
    expect_identical(bray(x, margin = mg), bray(dz))
    expect_identical(observed(x, margin = mg), observed(dz))
    expect_identical(pca(x, margin = mg, norm = 'binary', n_axes = 2, cpus = 1), pca(dz, norm = 'binary', n_axes = 2, cpus = 1))
  }

  dgt <- as(as(dn, COMPRESSED), TRIPLET)
  dup <- dgt
  dup@i <- c(dgt@i, dgt@i[1:2])
  dup@j <- c(dgt@j, dgt@j[1:2])
  dup@x <- c(dgt@x, 2, -dgt@x[2])  # the second repeat cancels its entry
  dd <- dn
  dd[cbind(dgt@i[1:2] + 1, dgt@j[1:2] + 1)] <- dd[cbind(dgt@i[1:2] + 1, dgt@j[1:2] + 1)] + c(2, -dgt@x[2])
  expect_equal(as.vector(bray(dup)), as.vector(bray(dd)), tolerance = 1e-14)
  expect_equal(as.vector(bray(Matrix::t(dup), margin = 2)), as.vector(bray(dd)), tolerance = 1e-14)


  # Invalid sparse objects are errors, not memory errors.
  bad <- as(dn, COMPRESSED)
  bad@p[3] <- 1000L
  expect_error(bray(bad, margin = 1), 'invalid column pointers')
  expect_error(bray(bad, margin = 2), 'invalid column pointers')
  bad <- as(dn, COMPRESSED)
  bad@i[1] <- 500L
  expect_error(bray(bad, margin = 2), 'out of range')
  s <- slam::as.simple_triplet_matrix(dn)
  s$i <- s$i[1:10]
  expect_error(bray(s), 'invalid')
  s <- slam::as.simple_triplet_matrix(dn)
  s$j <- as.double(s$j)
  expect_error(bray(s), 'invalid')

#})
