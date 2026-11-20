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



  # Same as prcomp() on the normalized table, up to sign: stored results,
  # so at home only. ====

  if (at_home()) {

    # Expected values: data-raw/expected/pca.r
    expected <- readRDS('expected/pca.rds')$cases

    for (cs in expected) {

      ord <- pca(
        counts      = cs$counts,
        norm        = cs$norm,
        pseudocount = cs$pseudocount,
        scale       = cs$scale,
        n_axes           = cs$k,
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

    # Sparse input with samples in columns, constant features with
    # `scale = TRUE` on both solver sides, and k larger than the rank.
    for (cs in readRDS('expected/pca.rds')$more_cases) {

      counts <- cs$counts
      margin <- 1L
      if (cs$sparse) {
        if (!requireNamespace('Matrix', quietly = TRUE)) next
        counts <- Matrix::Matrix(t(counts), sparse = TRUE)
        margin <- 2L
      }

      run <- function () pca(
        counts      = counts,
        margin      = margin,
        norm        = cs$norm,
        pseudocount = cs$pseudocount,
        scale       = cs$scale,
        n_axes           = cs$k,
        cpus        = 2 )

      info <- paste(dim(cs$counts), cs$norm, cs$scale, cs$k, cs$sparse, collapse = ' ')
      if (is.null(cs$warn)) {
        ord <- run()
      } else {
        expect_warning(ord <- run(), cs$warn, info = info)
      }

      ex   <- cs$expected
      k    <- cs$k_ret
      feat <- rownames(ex$rotation)
      load <- attr(ord, 'loadings')

      expect_identical(dim(ord), c(nrow(cs$counts), k), info = info)
      expect_identical(dimnames(ord), list(rownames(cs$counts), paste0('PC', seq_len(k))), info = info)
      expect_identical(dimnames(load), list(colnames(cs$counts), paste0('PC', seq_len(k))), info = info)

      expect_true(axis_diff(ord, ex$x) < 1e-10, info = info)
      expect_true(axis_diff(load[feat, , drop = FALSE], ex$rotation) < 1e-10, info = info)
      expect_equal(attr(ord, 'eig'),   ex$sdev[1:k]^2,                  tolerance = 1e-10, info = info)
      expect_equal(attr(ord, 'prop'),  ex$sdev[1:k]^2 / sum(ex$sdev^2), tolerance = 1e-10, info = info)
      expect_equal(attr(ord, 'trace'), sum(ex$sdev^2),                  tolerance = 1e-10, info = info)
      expect_true(all(sign(colSums(ord * ex$x)) == sign(colSums(load[feat, , drop = FALSE] * ex$rotation))), info = info)
      expect_true(all(load[setdiff(colnames(cs$counts), feat), ] == 0), info = info)
    }
  }



  # Same scores as pcoa() on the matching distances. ====

  for (m in list(m150, m30)) {

    ord <- pca(m)
    pco <- pcoa(euclidean(m))
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
    expect_equal(attr(ord, 'eig') * (nrow(m) - 1), attr(pco, 'eig'), tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'), attr(pco, 'prop'), tolerance = 1e-10)

    ord <- pca(m, norm = 'clr', pseudocount = 1, n_axes = 4)
    pco <- pcoa(aitchison(m, pseudocount = 1), n_axes = 4)
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
    expect_equal(attr(ord, 'prop'), attr(pco, 'prop'), tolerance = 1e-10)

    ord <- pca(m, norm = 'percent')
    pco <- pcoa(euclidean(m, norm = 'percent'))
    expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)
  }

  # More axes than one pass over the data computes at a time.
  ord <- pca(m150, n_axes = 9)
  pco <- pcoa(euclidean(m150), n_axes = 9)
  expect_equal(unclass(ord), unclass(pco), tolerance = 1e-9, check.attributes = FALSE)

  # The default pseudocount, with its warning, as aitchison() has.
  expect_warning(ord <- pca(m30, norm = 'clr'), 'pseudocount = 0.5')
  expect_equal(ord, pca(m30, norm = 'clr', pseudocount = 0.5))



  # Scores are the centered data times the loadings. ====

  ord  <- pca(m150, scale = TRUE, n_axes = 4)
  x    <- scale(m150)
  expect_equal(unclass(ord), x %*% attr(ord, 'loadings'), check.attributes = FALSE)
  expect_equal(crossprod(attr(ord, 'loadings')), diag(4), check.attributes = FALSE)
  expect_equal(attr(ord, 'trace'), ncol(m150))



  # Every input format gives the same result. ====

  x <- pca(m150, norm = 'clr', pseudocount = 1)
  expect_identical(no_display(pca(t(m150), margin = 2, norm = 'clr', pseudocount = 1)), no_display(x))

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(m150, sparse = TRUE)
    expect_identical(pca(sp, norm = 'clr', pseudocount = 1), x)
    expect_identical(pca(methods::as(sp, 'TsparseMatrix'), norm = 'clr', pseudocount = 1), x)
    expect_identical(no_display(pca(Matrix::t(sp), margin = 2, norm = 'clr', pseudocount = 1)), no_display(x))
  }

  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(m150)
    expect_identical(pca(st, norm = 'clr', pseudocount = 1), x)
  }

  # Binary is presence/absence.
  expect_equal(no_display(pca(m30, norm = 'binary')), no_display(pca((m30 > 0) + 0)))



  # Identical across runs and thread counts. ====

  for (m in list(m150, rand_counts(160, 400))) { # both solver sides
    x <- pca(m, norm = 'clr', pseudocount = 1, n_axes = 5, cpus = 1)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, n_axes = 5, cpus = 1), x)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, n_axes = 5, cpus = 2), x)
    expect_identical(pca(m, norm = 'clr', pseudocount = 1, n_axes = 5, cpus = 4), x)
  }

  # The shared ordination result format.
  expect_identical(
    names(attributes(no_display(x))),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings', 'class') )
  expect_identical(class(x), c('ecodive_ord', 'matrix', 'array'))
  # print() lists the loadings by name instead of printing them (and adds a
  # line for each record of how it was made).
  expect_true(length(capture.output(print(x))) < nrow(x) + 7)
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

  # Shifting every feature changes nothing: PCA centers them. (Features
  # far from 0 relative to their spread used to lose most of their digits.)
  set.seed(5)
  xn <- matrix(rnorm(1200), 60, 20)
  p0 <- pca(xn, n_axes = 3, cpus = 1)
  for (off in c(1e4, 1e9)) {
    p1 <- pca(xn + off, n_axes = 3, cpus = 1)
    expect_equal(attr(p1, 'eig'), attr(p0, 'eig'), tolerance = 1e-6) # x + 1e9 is itself rounded
    expect_equal(unclass(p1)[, ], unclass(p0)[, ], tolerance = 1e-6)
  }
  xs <- xn
  xs[, 1] <- 5 + xs[, 1] * 1e-6
  expect_equal(attr(pca(xs, scale = TRUE, n_axes = 3), 'eig'), attr(pca(xn, scale = TRUE, n_axes = 3), 'eig'), tolerance = 1e-8)

  # Constant data: no axes, not rounding noise.
  expect_warning(ord <- pca(matrix(3, 10, 5)), 'only 0 of the first 3')
  expect_identical(dim(ord), c(10L, 0L))

  # NA is an error, whatever the normalization.
  m_na <- m30
  m_na[3, 3] <- NA
  expect_error(pca(m_na, norm = 'binary'), 'must not contain NA')
  expect_error(pca(m_na, scale = TRUE),    'must not contain NA')

  # Bad input.
  expect_error(pca(m30, n_axes = 0),        '`n_axes` must be')
  expect_error(pca(m30, n_axes = 1.5),      '`n_axes` must be')
  expect_error(pca(m30, scale = NA),   '`scale` must be')
  expect_error(pca(m30, norm = 'xyz'))
  m[3, 3] <- NA
  expect_error(pca(m + 0), 'must not contain NA')

#})



#test_that("pcoa after moving its solver to lanczos.c", {

  # Same output as pcoa() had before. (The last bits depend on the compiler
  # and its optimization level, so the comparison allows for rounding.)
  # Expected values: data-raw/expected/pcoa_lanczos.r. Stored results, so
  # at home only.

  if (at_home()) {
    expected <- readRDS('expected/pcoa_lanczos.rds')
    for (cs in expected$cases) {
      d   <- beta_div(cs$counts, cs$metric, pseudocount = 1)
      ord <- unclass(pcoa(d, n_axes = cs$k, cpus = 2))
      attr(ord, 'correction') <- attr(ord, 'constant') <- NULL  # added since
      expect_equal(no_display(ord), no_display(cs$expected), tolerance = 1e-12)
    }
  }

#})
