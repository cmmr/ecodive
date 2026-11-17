#test_that("normalize", {

  set.seed(1)
  x <- matrix(rpois(12 * 9, 2) * rbinom(12 * 9, 1, 0.6), 12, 9)
  x[, 1] <- x[, 1] + 1L  # no empty samples
  dimnames(x) <- list(paste0('S', 1:12), paste0('F', 1:9))



  # Values, against their definitions. ====

  pct <- x / rowSums(x)
  expect_equal(normalize(x, 'percent', cpus = 1), pct)
  expect_equal(normalize(x, 'chord',   cpus = 1), x / sqrt(rowSums(x^2)))
  expect_equal(normalize(x, 'binary',  cpus = 1), (x > 0) * 1)
  expect_equal(normalize(x, 'none',    cpus = 1), x * 1)

  clr <- log(x + 1)
  clr <- clr - rowMeans(clr)
  expect_equal(normalize(x, 'clr', pseudocount = 1, cpus = 1), clr)
  expect_equal(normalize(x, 'CLR', pseudocount = 1, cpus = 1), clr)

  rclr <- t(apply(x, 1, function (r) {
    nz <- r > 0
    r[nz] <- log(r[nz]) - mean(log(r[nz]))
    r }))
  expect_equal(normalize(x, 'rclr', cpus = 1), rclr)

  # Already proportions: left as they are, as in the norm = paths.
  expect_identical(normalize(pct, 'percent', cpus = 1), pct)

  # Integer and logical input.
  xi <- x
  storage.mode(xi) <- 'integer'
  expect_equal(normalize(xi, 'percent', cpus = 1), pct)
  expect_equal(normalize(x > 0, 'percent', cpus = 1), (x > 0) / rowSums(x > 0))



  # Exactly the values used by the functions' own `norm` argument. ====

  # pca()'s CLR path keeps the table sparse, with each sample's value for
  # zeros as an offset, so it sums in another order: equal to rounding.
  expect_equal(
    pca(x, norm = 'clr', pseudocount = 1, k = 2, cpus = 1),
    pca(normalize(x, 'clr', pseudocount = 1, cpus = 1), k = 2, cpus = 1),
    tolerance = 1e-12 )
  expect_identical(
    pca(x, norm = 'percent', k = 2, cpus = 1),
    pca(normalize(x, 'percent', cpus = 1), k = 2, cpus = 1) )
  expect_identical(
    bray(x, norm = 'percent', cpus = 1),
    bray(normalize(x, 'percent', cpus = 1), cpus = 1) )
  expect_identical(
    euclidean(x, norm = 'rclr', cpus = 1),
    euclidean(normalize(x, 'rclr', cpus = 1), cpus = 1) )
  expect_equal(
    aitchison(x, pseudocount = 1, cpus = 1),
    euclidean(normalize(x, 'clr', pseudocount = 1, cpus = 1), cpus = 1), tolerance = 1e-12 )

  ord <- pcoa(bray(x, cpus = 1), k = 2, cpus = 1)
  for (nm in c('percent', 'chord', 'rclr'))
    expect_equal(
      envfit(ord, counts = x, norm = nm, permutations = 19, cpus = 1),
      envfit(ord, counts = normalize(x, nm, cpus = 1), permutations = 19, cpus = 1),
      tolerance = 1e-12, info = nm )
  expect_equal(
    envfit(ord, counts = x, norm = 'clr', pseudocount = 1, permutations = 19, cpus = 1),
    envfit(ord, counts = normalize(x, 'clr', pseudocount = 1, cpus = 1), permutations = 19, cpus = 1),
    tolerance = 1e-10 )



  # Output form: orientation, dimnames, and class kept. ====

  # Samples in columns.
  expect_equal(normalize(t(x), 'percent', margin = 2, cpus = 1), t(pct))
  expect_equal(normalize(t(x), 'clr', margin = 2, pseudocount = 1, cpus = 1), t(clr))

  # No dimnames.
  expect_null(dimnames(normalize(unname(x), 'percent', cpus = 1)))

  if (requireNamespace('Matrix', quietly = TRUE)) {

    sp <- Matrix::Matrix(x, sparse = TRUE)
    r  <- normalize(sp, 'percent', cpus = 1)
    expect_inherits(r, 'dgCMatrix')
    expect_true(methods::validObject(r))
    expect_identical(dimnames(r), dimnames(x))
    expect_equal(as.matrix(r), pct)
    expect_identical(length(r@x), sum(x > 0))

    r2 <- normalize(Matrix::t(sp), 'rclr', margin = 2, cpus = 1)
    expect_inherits(r2, 'dgCMatrix')
    expect_true(methods::validObject(r2))
    expect_equal(as.matrix(r2), t(rclr))

    tr <- normalize(methods::as(sp, 'TsparseMatrix'), 'chord', cpus = 1)
    expect_inherits(tr, 'dgTMatrix')
    expect_true(methods::validObject(tr))
    expect_equal(as.matrix(tr), x / sqrt(rowSums(x^2)))

    de <- normalize(Matrix::Matrix(x, sparse = FALSE), 'binary', cpus = 1)
    expect_inherits(de, 'dgeMatrix')
    expect_equal(as.matrix(de), (x > 0) * 1)

    # CLR is dense.
    cl <- normalize(sp, 'clr', pseudocount = 1, cpus = 1)
    expect_true(is.matrix(cl) && is.double(cl))
    expect_equal(cl, clr)

    # The input isn't changed.
    expect_identical(sp, Matrix::Matrix(x, sparse = TRUE))
  }

  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(x)
    r  <- normalize(st, 'percent', cpus = 1)
    expect_inherits(r, 'simple_triplet_matrix')
    expect_identical(dimnames(r), dimnames(x))
    expect_equal(as.matrix(r), pct)
    expect_equal(normalize(st, 'clr', pseudocount = 1, cpus = 1), clr)
  }

  # The same for any number of cpus.
  big <- matrix(rpois(300 * 50, 1), 300, 50)
  big[, 1] <- big[, 1] + 1L
  for (nm in c('percent', 'clr', 'rclr'))
    expect_identical(
      suppressWarnings(normalize(big, nm, cpus = 1)),
      suppressWarnings(normalize(big, nm, cpus = 2)), info = nm )



  # Errors and warnings. ====

  expect_error(normalize(x), '`norm` is required')
  expect_error(normalize(x, 1), '`norm` must be')
  expect_error(normalize(x, 'log'))
  expect_error(normalize(x, 'percent', margin = 3), 'margin')
  expect_error(normalize(x, 'percent', cpus = 0), 'cpus')
  expect_warning(normalize(x, 'clr', cpus = 1), 'pseudocount')
  expect_silent(normalize(x, 'clr', pseudocount = 1, cpus = 1))

  # 'binary' keeps NA as NA, not a presence.
  xna <- matrix(c(NA, 2, 0, 3), 2)
  expect_identical(normalize(xna, 'binary', cpus = 1), matrix(c(NA, 1, 0, 1), 2))

  # An all-zero table has no stored values to look at.
  expect_identical(normalize(matrix(0, 2, 3), 'percent', cpus = 1), matrix(0, 2, 3))

#})
