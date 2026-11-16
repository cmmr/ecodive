#test_that("pca", {

  # Largest absolute difference between two matrices, relative to the
  # largest value of `y`, after matching each column's sign.
  axis_diff <- function (x, y) {
    x <- matrix(x, nrow(x))
    y <- matrix(y, nrow(y))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # Reproducible test data.
  set.seed(1)
  rand_counts <- function (n, p = 40) {
    m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
    m[, 1] <- m[, 1] + 1
    dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
    m
  }
  m150 <- rand_counts(150)
  m30  <- rand_counts(30, 90)



  # Same as prcomp() on the normalized table, up to sign. ====
  # Expected values: data-raw/expected/pca.r

  expected <- readRDS('expected/pca.rds')$cases

  for (cs in expected) {

    ord <- pca(
      counts      = cs$counts,
      norm        = cs$norm,
      pseudocount = cs$pseudocount,
      scale       = cs$scale,
      k           = cs$k,
      cpus        = 1 )

    ex   <- cs$expected
    k    <- cs$k
    feat <- rownames(ex$rotation)
    load <- attr(ord, 'loadings')

    expect_identical(dim(ord), c(nrow(cs$counts), k))
    expect_identical(dimnames(ord), list(rownames(cs$counts), paste0('PC', 1:k)))
    expect_identical(dimnames(load), list(colnames(cs$counts), paste0('PC', 1:k)))

    expect_true(axis_diff(ord, ex$x) < 1e-10)
    expect_true(axis_diff(load[feat, ], ex$rotation) < 1e-10)
    expect_equal(attr(ord, 'eig'),   ex$sdev[1:k]^2,                tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'),  ex$sdev[1:k]^2 / sum(ex$sdev^2), tolerance = 1e-10)
    expect_equal(attr(ord, 'trace'), sum(ex$sdev^2),                tolerance = 1e-10)

    # Loadings and scores are signed alike: scores = centered data %*% loadings.
    expect_true(all(sign(colSums(ord * ex$x)) == sign(colSums(load[feat, ] * ex$rotation))))

    # Constant features, left out with `scale = TRUE`, have no loading.
    expect_true(all(load[setdiff(colnames(cs$counts), feat), ] == 0))
  }



  # Same scores as pcoa() on the matching distances. ====

  for (m in list(m150, m30)) {

    ord <- pca(m)
    pco <- pcoa(euclidean(m))
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
    expect_equal(attr(ord, 'eig') * (nrow(m) - 1), attr(pco, 'eig'), tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'), attr(pco, 'prop'), tolerance = 1e-10)

    ord <- pca(m, norm = 'clr', pseudocount = 1, k = 4)
    pco <- pcoa(aitchison(m, pseudocount = 1), k = 4)
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
    expect_equal(attr(ord, 'prop'), attr(pco, 'prop'), tolerance = 1e-10)

    ord <- pca(m, norm = 'percent')
    pco <- pcoa(euclidean(m, norm = 'percent'))
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
  }

  # More axes than one pass over the data computes at a time.
  ord <- pca(m150, k = 9)
  pco <- pcoa(euclidean(m150), k = 9)
  expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)

  # The default pseudocount, with its warning, as aitchison() has.
  expect_warning(ord <- pca(m30, norm = 'clr'), 'pseudocount = 0.5')
  expect_equal(ord, pca(m30, norm = 'clr', pseudocount = 0.5))



  # Scores are the centered data times the loadings. ====

  ord  <- pca(m150, scale = TRUE, k = 4)
  x    <- scale(m150)
  expect_equal(unclass(ord), x %*% attr(ord, 'loadings'), check.attributes = FALSE)
  expect_equal(crossprod(attr(ord, 'loadings')), diag(4), check.attributes = FALSE)
  expect_equal(attr(ord, 'trace'), ncol(m150))



  # Every input format gives the same result. ====

  x <- pca(m150, norm = 'clr', pseudocount = 1)
  expect_identical(pca(t(m150), margin = 2, norm = 'clr', pseudocount = 1), x)

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(m150, sparse = TRUE)
    expect_identical(pca(sp, norm = 'clr', pseudocount = 1), x)
    expect_identical(pca(methods::as(sp, 'TsparseMatrix'), norm = 'clr', pseudocount = 1), x)
    expect_identical(pca(Matrix::t(sp), margin = 2, norm = 'clr', pseudocount = 1), x)
  }

  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(m150)
    expect_identical(pca(st, norm = 'clr', pseudocount = 1), x)
  }

  # Binary is presence/absence.
  expect_equal(pca(m30, norm = 'binary'), pca((m30 > 0) + 0))



  # Identical across runs and thread counts. ====

  for (m in list(m150, rand_counts(160, 400))) { # both solver sides
    x <- pca(m, norm = 'clr', pseudocount = 1, k = 5, cpus = 1)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, k = 5, cpus = 1), x)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, k = 5, cpus = 2), x)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, k = 5, cpus = 4), x)
  }

  # The shared ordination result format.
  expect_identical(
    names(attributes(x)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings') )
  expect_identical(attr(x, 'method'), 'pca')

  # Each axis's largest absolute score is positive.
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))



  # Edge cases. ====

  # Points on a plane in 3 dimensions have only two positive eigenvalues.
  pts <- cbind(sin(1:50), cos(1:50), sin(1:50) + cos(1:50))
  expect_warning(ord <- pca(pts), 'only 2 of the first 3 eigenvalues')
  expect_identical(ncol(ord), 2L)
  expect_identical(dim(attr(ord, 'loadings')), c(3L, 2L))

  # k is limited by the number of samples and of features.
  expect_identical(ncol(pca(m150[, 1:2])), 2L)
  expect_identical(ncol(pca(m150[1:3, ])), 2L)

  # Two samples: one axis, explaining everything.
  ord <- pca(m150[1:2, ])
  expect_identical(dim(ord), c(2L, 1L))
  expect_equal(attr(ord, 'prop'), 1)
  expect_equal(attr(ord, 'trace'), sum(apply(m150[1:2, ], 2, var)))

  # One sample: no axes.
  ord <- pca(m150[1, , drop = FALSE])
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(rownames(ord), 'S1')
  expect_identical(attr(ord, 'eig'), numeric(0))
  expect_identical(attr(ord, 'trace'), 0)
  expect_identical(dim(attr(ord, 'loadings')), c(40L, 0L))

  # No names.
  ord <- pca(unname(m30))
  expect_null(rownames(ord))
  expect_null(rownames(attr(ord, 'loadings')))
  expect_identical(colnames(ord), paste0('PC', 1:3))

  # Integer input.
  m <- m30
  storage.mode(m) <- 'integer'
  expect_identical(pca(m), pca(m30))

  # Bad input.
  expect_error(pca(m30, k = 0),        '`k` must be')
  expect_error(pca(m30, k = 1.5),      '`k` must be')
  expect_error(pca(m30, scale = NA),   '`scale` must be')
  expect_error(pca(m30, norm = 'xyz'))
  m[3, 3] <- NA
  expect_error(pca(m + 0), 'must not contain NA')

#})



#test_that("pcoa after moving its solver to lanczos.c", {

  # Same output as pcoa() had before. The last bits depend on the compiler
  # and its optimization level, so the expected values hold an optimized
  # (O2) and an unoptimized (O0, as from load_all()) build's results. On the
  # platform they were made on, the output must be bit-identical to one of
  # them; elsewhere, equal.
  # Expected values: data-raw/expected/pcoa_lanczos.r

  expected  <- readRDS('expected/pcoa_lanczos.rds')
  same_plat <- identical(expected$generated$O2$platform, R.version$platform)

  for (cs in expected$cases) {
    d   <- beta_div(cs$counts, cs$metric, pseudocount = 1)
    ord <- pcoa(d, k = cs$k, cpus = 2)
    expect_equal(ord, cs$expected$O2, tolerance = 1e-10)
    if (same_plat)
      expect_true(identical(ord, cs$expected$O2) || identical(ord, cs$expected$O0))
  }

#})
