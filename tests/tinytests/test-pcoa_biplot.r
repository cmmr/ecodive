#test_that("pcoa_biplot", {

  # The same biplot from base R's matrix algebra, as scikit-bio computes it,
  # with each feature centered.
  direct <- function (ord, y) {
    y   <- scale(y[rownames(ord), , drop = FALSE], scale = FALSE)
    z   <- scale(unclass(ord))
    n   <- nrow(y)
    spc <- crossprod(y, z) / (n - 1)
    sqrt(n - 1) * spc %*% diag(attr(ord, 'eig')^-0.5, ncol(ord))
  }

  cpus <- if (at_home()) 4L else 2L

  set.seed(1)
  m <- matrix(rpois(40 * 15, 3) * rbinom(40 * 15, 1, 0.5), 40, 15)
  m[, 1] <- m[, 1] + 1
  dimnames(m) <- list(paste0('S', 1:40), paste0('OTU', 1:15))

  ord <- pcoa(bray(m), n_axes = 4)



  # The result. ====

  x <- pcoa_biplot(ord, m, cpus = 1)
  expect_identical(dim(x), c(15L, 4L))
  expect_identical(dimnames(x), list(colnames(m), colnames(ord)))
  expect_identical(attr(x, 'method'), 'pcoa_biplot')
  expect_false(inherits(x, 'ecodive_ord'))
  expect_identical(names(attributes(x)), c('dim', 'dimnames', 'method'))
  expect_equal(unclass(x), direct(ord, m), tolerance = 1e-12, check.attributes = FALSE)

  # The slope of each feature on each axis: cov(y, x) / var(x), so
  # Yc' x_j / eig_j.
  yc <- scale(m, scale = FALSE)
  expect_equal(
    unclass(x), crossprod(yc, unclass(ord)) %*% diag(1 / attr(ord, 'eig')),
    tolerance = 1e-12, check.attributes = FALSE )
  for (j in 1:4)
    expect_equal(
      unname(x[, j]), unname(coef(lm(m ~ ord[, j]))[2, ]),
      tolerance = 1e-10 )

  # The same direction as envfit()'s arrows, which are unit vectors.
  ef  <- envfit(ord, counts = m, permutations = 0, cpus = 1)
  dir <- x / sqrt(rowSums(x^2))
  ok  <- rowSums(x^2) > 0
  expect_equal(
    unname(dir[ok, ]), unname(as.matrix(ef[ok, colnames(ord)])),
    tolerance = 1e-10 )



  # Against base R, for each transformation. ====

  # Each transformation of a table, in base R.
  trans <- function (y, norm, pc = 1) switch(
    norm,
    none    = y,
    percent = y / rowSums(y),
    binary  = (y != 0) + 0,
    chord   = y / sqrt(rowSums(y^2)),
    clr     = t(apply(log(y + pc), 1, function (x) x - mean(x))),
    rclr    = t(apply(y, 1, function (x) {
      l <- log(x[x > 0])
      x[x > 0] <- l - mean(l)
      x })) )

  norms <- c('none', 'percent', 'binary', 'chord', 'clr', 'rclr')
  pc    <- function (norm, value = 1) if (norm == 'clr') value

  for (nm in norms)
    expect_equal(
      unclass(pcoa_biplot(ord, m, norm = nm, pseudocount = pc(nm))), direct(ord, trans(m, nm)),
      tolerance = 1e-12, check.attributes = FALSE, info = nm )
  expect_warning(pcoa_biplot(ord, m, norm = 'clr'), 'pseudocount')

  # A corrected ordination uses its corrected coordinates and eigenvalues.
  for (corr in c('lingoes', 'cailliez')) {
    oc <- pcoa(bray(m), n_axes = 4, correction = corr)
    expect_equal(unclass(pcoa_biplot(oc, m)), direct(oc, m), tolerance = 1e-12, check.attributes = FALSE)
  }

  # Variables with large means lose no precision: rounding error is
  # relative to their spread (a few units), not to their means (1e8).
  env <- cbind(a = 1e8 + m[, 2], b = -1e8 + m[, 3] / 7)
  expect_equal(
    unclass(pcoa_biplot(ord, env)), direct(ord, env),
    tolerance = 1e-14, check.attributes = FALSE )

  # Features with the same value in every sample, or all zeros, score
  # exactly 0.
  m2 <- cbind(m, const = 7, zero = 0)
  s2 <- pcoa_biplot(ord, m2)
  expect_identical(unname(s2['const', ]), rep(0, 4))
  expect_identical(unname(s2['zero', ]),  rep(0, 4))
  expect_identical(s2[1:15, ], x[1:15, ])



  # Every input type, margin, and thread count gives the same scores. ====

  expect_identical(pcoa_biplot(ord, m * 1L),            x)
  expect_identical(pcoa_biplot(ord, t(m), margin = 2),  x)
  expect_identical(pcoa_biplot(ord, as.data.frame(m)),  x)
  expect_identical(pcoa_biplot(ord, m[40:1, ]),         x)

  # Extra samples in `counts` are ignored, as in envfit().
  expect_identical(pcoa_biplot(ord, rbind(m, extra = 1:15)), x)

  # Threads.
  set.seed(4)
  big  <- matrix(rpois(40 * 300, 1), 40, 300, dimnames = list(rownames(m), paste0('F', 1:300)))
  seed <- .Random.seed
  xb   <- lapply(norms, function (nm) pcoa_biplot(ord, big, norm = nm, pseudocount = pc(nm, 0.5), cpus = 1))
  expect_identical(.Random.seed, seed)
  for (i in seq_along(norms)) {
    nm <- norms[[i]]
    expect_identical(pcoa_biplot(ord, big, norm = nm, pseudocount = pc(nm, 0.5), cpus = 2),    xb[[i]], info = nm)
    expect_identical(pcoa_biplot(ord, big, norm = nm, pseudocount = pc(nm, 0.5), cpus = cpus), xb[[i]], info = nm)
    expect_equal(unclass(xb[[i]]), direct(ord, trans(big, nm, 0.5)), tolerance = 1e-12, check.attributes = FALSE, info = nm)
  }

  # A slam simple_triplet_matrix (built by hand, without slam itself).
  stm <- function (y) {
    nz <- which(y != 0, arr.ind = TRUE)
    structure(
      list(i = unname(nz[, 1]), j = unname(nz[, 2]), v = y[nz], nrow = nrow(y), ncol = ncol(y),
           dimnames = dimnames(y)),
      class = 'simple_triplet_matrix' )
  }
  for (i in seq_along(norms)) {
    nm <- norms[[i]]
    expect_identical(pcoa_biplot(ord, stm(big), norm = nm, pseudocount = pc(nm, 0.5)), xb[[i]], info = nm)
  }

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(big, sparse = TRUE)
    for (i in seq_along(norms)) {
      nm <- norms[[i]]
      f  <- function (y, ...) pcoa_biplot(ord, y, norm = nm, pseudocount = pc(nm, 0.5), ...)
      expect_identical(f(sp),                                 xb[[i]], info = nm)
      expect_identical(f(Matrix::t(sp), margin = 2),          xb[[i]], info = nm)
      expect_identical(f(methods::as(sp, 'TsparseMatrix')),   xb[[i]], info = nm)
    }
  }



  # Edge cases. ====

  # Unnamed samples are matched by position.
  o  <- pcoa(bray(unname(m)), n_axes = 2)
  on <- structure(o, dimnames = list(rownames(m), colnames(o)))
  expect_identical(unname(pcoa_biplot(o, unname(m))), unname(pcoa_biplot(on, m)))
  expect_identical(dimnames(pcoa_biplot(o, unname(m))), list(NULL, c('PCoA1', 'PCoA2')))

  # No axes: one sample, or no positive eigenvalue.
  o1 <- pcoa(dist(m[1, , drop = FALSE]))
  s1 <- pcoa_biplot(o1, m[1, , drop = FALSE])
  expect_identical(dim(s1), c(15L, 0L))
  expect_identical(rownames(s1), colnames(m))
  s3 <- rownames(m)[1:3]
  o0 <- suppressWarnings(pcoa(as.dist(matrix(0, 3, 3, dimnames = list(s3, s3)))))
  expect_identical(dim(pcoa_biplot(o0, m[1:3, ])), c(15L, 0L))

  # Two samples.
  o2 <- pcoa(bray(m[1:2, ]))
  expect_equal(unclass(pcoa_biplot(o2, m[1:2, ])), direct(o2, m[1:2, ]), tolerance = 1e-12, check.attributes = FALSE)



  # Invalid input. ====

  expect_error(pcoa_biplot(unclass(ord)[, 1:2], m),     'result from pcoa')
  expect_error(pcoa_biplot(ord[, 1:2], m),              'result from pcoa')
  expect_error(pcoa_biplot(bray(m), m),                 'result from pcoa')
  expect_error(pcoa_biplot(structure(ord, method = 'pca'), m), 'result from pcoa')
  expect_error(pcoa_biplot(structure(ord, eig = NULL), m), 'result from pcoa')
  expect_error(pcoa_biplot(structure(ord, eig = c(1, NA, 1, 1)), m), 'result from pcoa')
  expect_error(pcoa_biplot(structure(ord, eig = c(1, 0, 1, 1)), m), 'result from pcoa')
  hand <- structure(matrix(1, 1, 1, dimnames = list('S1', 'PCoA1')), method = 'pcoa', eig = 1)
  expect_error(pcoa_biplot(hand, m[1, , drop = FALSE]), 'result from pcoa')

  expect_error(pcoa_biplot(ord, m[-1, ]),               'no sample in `counts`')
  expect_error(pcoa_biplot(ord, m[c(1, 1:40), ]),       'duplicate sample names')
  o_dup <- ord
  rownames(o_dup)[2] <- 'S1'
  expect_error(pcoa_biplot(o_dup, m),                   'duplicate sample names')
  expect_error(pcoa_biplot(o, unname(m)[-1, ]),         '39 samples')

  m_na <- m; m_na[2, 3] <- NA
  expect_error(pcoa_biplot(ord, m_na),                  'NA')
  expect_error(pcoa_biplot(ord, m_na * 1L),             'NA')
  m_na[2, 3] <- Inf
  expect_error(pcoa_biplot(ord, m_na),                  'infinite')
  expect_error(pcoa_biplot(ord, matrix(letters[1:80], 40, 2, dimnames = list(rownames(m), NULL))))
  expect_error(pcoa_biplot(ord, m, margin = 3),         'margin')
  expect_error(pcoa_biplot(ord, m, norm = 'xyz'))
  expect_error(pcoa_biplot(ord, m, cpus = 0),           'cpus')

  m_neg <- m; m_neg[2, 3] <- -1
  expect_error(pcoa_biplot(ord, m_neg, norm = 'clr', pseudocount = 1), 'negative')



  # Against scikit-bio: stored results, so at home only. ====

  if (at_home()) {

    # Expected results from data-raw/expected/pcoa_biplot.r.
    ex <- readRDS('expected/pcoa_biplot.rds')

    for (nm in names(ex$fits)) {

      fit <- ex$fits[[nm]]
      ord <- pcoa(fit$dist, n_axes = 5, cpus = 2)
      sgn <- sign(colSums(unclass(ord) * fit$coords))

      # scikit-bio's pcoa().
      expect_equal(attr(ord, 'eig'),  fit$eig,  tolerance = 1e-12, info = nm)
      expect_equal(attr(ord, 'prop'), fit$prop, tolerance = 1e-12, info = nm)
      expect_equal(
        sweep(unclass(ord), 2, sgn, '*'), fit$coords,
        tolerance = 1e-12, check.attributes = FALSE, info = nm )

      # Its pcoa_biplot(), for counts, proportions, and variables. Samples
      # are matched by name.
      for (lab in names(fit$biplot)) {

        ref    <- fit$biplot[[lab]]
        y      <- if (lab == 'env') ex$env else ex$counts
        y      <- y[c(10:1, 11:nrow(y)), ]
        scores <- pcoa_biplot(ord, y, norm = if (lab == 'percent') 'percent' else 'none', cpus = 2)

        expect_identical(dimnames(scores), list(colnames(y), colnames(ord)))
        expect_equal(
          sweep(unclass(scores), 2, sgn, '*'), ref,
          tolerance = 1e-10, check.attributes = FALSE, info = paste(nm, lab) )
      }
    }

    # An all-zero feature and one with the same count everywhere score
    # exactly 0.
    ord    <- pcoa(ex$fits$bray$dist, n_axes = 5)
    scores <- pcoa_biplot(ord, ex$counts)
    expect_identical(unname(scores['OTU11', ]), rep(0, 5))
    expect_identical(unname(scores['OTU12', ]), rep(0, 5))
  }

#})
