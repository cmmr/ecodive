#test_that("pcoa", {

  # Largest absolute difference between two coordinate matrices, relative to
  # the largest coordinate, after matching each axis's sign. `y` can have
  # extra axes.
  coord_diff <- function (x, y) {
    x <- matrix(x, nrow(x))
    y <- matrix(y, nrow(y))[, seq_len(ncol(x)), drop = FALSE]
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
  m30  <- rand_counts(30)



  # Same as cmdscale(), up to sign. ====

  for (d in list(bray(m150), euclidean(m150), jaccard(m150), bray(m30), bray(ecodive::ex_counts))) {

    k   <- min(3L, attr(d, 'Size') - 1L)
    ord <- pcoa(d, cpus = 1)
    cm  <- cmdscale(d, k = k, eig = TRUE)

    expect_identical(dim(ord), c(attr(d, 'Size'), k))
    expect_identical(dimnames(ord), list(attr(d, 'Labels'), paste0('PCoA', 1:k)))
    expect_true(coord_diff(ord, cm$points) < 1e-9)
    expect_equal(attr(ord, 'eig'),   cm$eig[1:k],     tolerance = 1e-10)
    expect_equal(attr(ord, 'trace'), sum(cm$eig),     tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'),  cm$eig[1:k] / sum(cm$eig), tolerance = 1e-10)
  }

  # More axes than one pass over `dm` computes at a time.
  d   <- bray(m150)
  ord <- pcoa(d, k = 6, cpus = 1)
  cm  <- cmdscale(d, k = 6, eig = TRUE)
  expect_true(coord_diff(ord, cm$points) < 1e-9)
  expect_equal(attr(ord, 'eig'), cm$eig[1:6], tolerance = 1e-10)

  # A symmetric matrix gives the same result as a dist.
  m <- as.matrix(d)
  expect_identical(pcoa(m), pcoa(d))
  
  # Labels come from row names, else column names.
  rownames(m) <- NULL
  expect_identical(rownames(pcoa(m)), colnames(m))



  # Same as ape and vegan, up to sign: stored results, so at home only. ====

  if (at_home()) {

    # Expected results from data-raw/expected/pcoa.r.
    ex <- readRDS('expected/pcoa.rds')
    ord <- pcoa(bray(ex$counts))

    # ape::pcoa(). Bray-Curtis has negative eigenvalues; prop matches
    # Relative_eig.
    expect_true(ex$ape$min_eig < 0)
    expect_true(coord_diff(ord, ex$ape$vectors) < 1e-9)
    expect_equal(attr(ord, 'eig'),   ex$ape$eig,          tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'),  ex$ape$relative_eig, tolerance = 1e-10)
    expect_equal(attr(ord, 'trace'), ex$ape$trace,        tolerance = 1e-10)

    # vegan::wcmdscale() and its "Proportion Explained".
    expect_true(coord_diff(ord, ex$vegan$points) < 1e-9)
    expect_equal(attr(ord, 'eig'),  ex$vegan$eig,  tolerance = 1e-10)
    expect_equal(attr(ord, 'prop'), ex$vegan$prop, tolerance = 1e-10)

    # More axes than positive eigenvalues, Euclidean distances, two and
    # three samples, duplicate samples, and 300 samples.
    for (cs in ex$cases) {

      d <- if (is.null(cs$pts)) bray(cs$counts) else dist(cs$pts)

      if (is.null(cs$warn)) {
        ord <- pcoa(d, k = cs$k, cpus = 2)
      } else {
        expect_warning(ord <- pcoa(d, k = cs$k, cpus = 2), cs$warn, info = cs$name)
      }

      k <- length(cs$ape$eig)
      expect_identical(dim(ord), c(attr(d, 'Size'), k), info = cs$name)
      expect_true(coord_diff(ord, cs$ape$vectors) < 1e-9, info = cs$name)
      expect_equal(attr(ord, 'eig'),   cs$ape$eig,          tolerance = 1e-10, info = cs$name)
      expect_equal(attr(ord, 'prop'),  cs$ape$relative_eig, tolerance = 1e-10, info = cs$name)
      expect_equal(attr(ord, 'trace'), cs$ape$trace,        tolerance = 1e-10, info = cs$name)
      expect_true(coord_diff(ord, cs$vegan$points) < 1e-9, info = cs$name)
      expect_equal(attr(ord, 'eig'),   cs$vegan$eig,        tolerance = 1e-10, info = cs$name)
      expect_equal(attr(ord, 'prop'),  cs$vegan$prop,       tolerance = 1e-10, info = cs$name)
    }
  }



  # Identical across runs and thread counts. ====

  d <- bray(m150)
  x <- pcoa(d, cpus = 1)
  expect_identical(pcoa(d, cpus = 1), x)
  expect_identical(pcoa(d, cpus = 2), x)
  expect_identical(pcoa(d, cpus = 4), x)

  # The shared ordination result format.
  expect_identical(
    names(attributes(x)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'correction', 'constant', 'class') )
  expect_identical(attr(x, 'correction'), 'none')
  expect_identical(attr(x, 'constant'), 0)
  expect_identical(attr(x, 'method'), 'pcoa')
  expect_identical(class(x), c('ecodive_ord', 'matrix', 'array'))
  expect_true(is.matrix(x) && is.numeric(x))

  # Subsetting and data frames give plain objects; print() names the
  # attributes instead of showing them.
  expect_identical(attributes(x[, 1:2]), list(dim = c(150L, 2L), dimnames = list(rownames(x), c('PCoA1', 'PCoA2'))))
  expect_identical(class(as.data.frame(x)), 'data.frame')
  expect_identical(dim(as.data.frame(x)), c(150L, 3L))
  out <- capture.output(print(x))
  expect_true(any(grepl('^pcoa: 150 samples, 3 axes; % of variance: ', out)))
  expect_true(any(grepl('^Attributes: method, eig, prop, trace, correction, constant$', out)))
  expect_identical(capture.output(print(x)), out)
  invisible(capture.output(v <- withVisible(print(x))))
  expect_identical(v$visible, FALSE)
  expect_identical(v$value, x)

  # Each axis's largest absolute coordinate is positive.
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))



  # Close and tied eigenvalues. ====

  # Points in 6 dimensions with chosen eigenvalues (squared spreads).
  config <- function (n, s) {
    set.seed(2)
    u <- qr.Q(qr(scale(matrix(rnorm(n * length(s)), n), scale = FALSE)))
    dist(u %*% diag(s))
  }

  # lambda3 and lambda4 within 0.02%.
  d   <- config(150, c(10, 5, 2, 2 * (1 - 1e-4), 1, 0.5))
  ord <- pcoa(d)
  cm  <- cmdscale(d, k = 3, eig = TRUE)
  expect_equal(attr(ord, 'eig'), c(100, 25, 4), tolerance = 1e-12)
  expect_true(coord_diff(ord, cm$points) < 1e-9)

  # lambda3 = lambda4 = lambda5: any basis of their space is right.
  d   <- config(150, c(10, 5, 2, 2, 2, 1))
  ord <- pcoa(d, k = 5)
  expect_equal(attr(ord, 'eig'), c(100, 25, 4, 4, 4), tolerance = 1e-12)
  expect_equal(crossprod(unclass(ord)), diag(c(100, 25, 4, 4, 4)), tolerance = 1e-10, check.attributes = FALSE)

  # Equidistant points: all n - 1 eigenvalues are 1/2.
  d   <- as.dist(matrix(1, 200, 200))
  ord <- pcoa(d)
  expect_equal(attr(ord, 'eig'), rep(0.5, 3), tolerance = 1e-12)
  expect_equal(attr(ord, 'trace'), 199 / 2)



  # Edge cases. ====

  # Points on a plane have only two positive eigenvalues.
  pts <- cbind(x = sin(1:120), y = cos(1:120 / 2) * 3)
  expect_warning(ord <- pcoa(dist(pts)), 'only 2 of the first 3 eigenvalues')
  expect_identical(ncol(ord), 2L)
  expect_true(coord_diff(ord, cmdscale(dist(pts), k = 2)) < 1e-9)

  # All distances zero.
  expect_warning(ord <- pcoa(as.dist(matrix(0, 5, 5))), 'only 0 of the first 3')
  expect_identical(dim(ord), c(5L, 0L))
  expect_identical(attr(ord, 'trace'), 0)

  # One, two, and three samples.
  ord <- pcoa(dist(matrix(1:3, 1, dimnames = list('A', NULL))))
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(rownames(ord), 'A')
  expect_identical(attr(ord, 'eig'), numeric(0))

  ord <- pcoa(dist(c(a = 1, b = 3)))
  expect_equal(unclass(ord), matrix(c(1, -1), 2, dimnames = list(c('a', 'b'), 'PCoA1')), check.attributes = FALSE)
  expect_identical(dimnames(ord), list(c('a', 'b'), 'PCoA1'))
  expect_equal(attr(ord, 'eig'), 2)
  expect_equal(attr(ord, 'prop'), 1)

  d   <- bray(counts[1:3, ])
  ord <- pcoa(d)
  expect_identical(ncol(ord), 2L)
  expect_true(coord_diff(ord, cmdscale(d, k = 2)) < 1e-9)

  # k larger than n - 1, and integer storage.
  d <- as.dist(matrix(c(0L, 3L, 4L, 3L, 0L, 5L, 4L, 5L, 0L), 3))
  expect_identical(pcoa(d, k = 10), pcoa(d, k = 2))
  expect_equal(attr(pcoa(d), 'eig'), cmdscale(d, k = 2, eig = TRUE)$eig[1:2])



  # Invalid input. ====

  d <- bray(counts)

  d_na <- d; d_na[2] <- NA
  expect_error(pcoa(d_na), 'NA')
  d_na[2] <- Inf
  expect_error(pcoa(d_na), 'infinite')

  expect_error(pcoa(d, k = 0))
  expect_error(pcoa(d, k = 1.5))
  expect_error(pcoa(d, k = NA))
  expect_error(pcoa(d, k = 'a'))
  expect_error(pcoa(d, cpus = 0))

  expect_error(pcoa(matrix(1:6, 2)),               'symmetric')
  expect_error(pcoa(matrix(c(0, 1, 2, 0), 2)),     'symmetric')
  expect_error(pcoa(matrix(letters[1:4], 2)),      'symmetric')
  expect_error(pcoa(structure(1:2, class = 'dist', Size = 3L)), 'valid dist')
  expect_error(pcoa(structure(1, class = 'dist', Size = -1L)),  'valid dist')
  # n * (n - 1) / 2 overflowed an integer from n = 46342.
  expect_error(pcoa(structure(1:2, class = "dist", Size = 46342L)), "valid dist")

  # Negative distances: an error, beyond rounding error, which is set to 0.
  d_neg <- d; d_neg[2] <- -0.1
  expect_error(pcoa(d_neg), 'must not contain negative distances')
  d_neg[2] <- -1e-15
  d_0 <- d; d_0[2] <- 0
  expect_identical(pcoa(d_neg, k = 2, cpus = 1), pcoa(d_0, k = 2, cpus = 1))

  # Distances too large to square.
  expect_error(pcoa(d * 1e160), 'too large to square')

  # A formula is not distances.
  expect_error(pcoa(d ~ x), 'not a formula')

#})
