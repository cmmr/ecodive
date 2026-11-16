#test_that("nmds", {

  # Test data and expected results, from data-raw/expected/nmds.r.
  ex <- readRDS('expected/nmds.rds')

  d_bray <- bray(ex$inputs$bray)
  d_jacc <- jaccard(ex$inputs$jacc)  # Many tied distances.
  d_pln  <- dist(ex$inputs$plane)    # Points near a plane: a clear solution.
  dists  <- list(bray = d_bray, jacc = d_jacc, plane = d_pln)

  # Largest absolute difference between two coordinate matrices, relative to
  # the largest coordinate, after matching each axis's sign.
  coord_diff <- function (x, y) {
    x <- matrix(x, nrow(x))
    y <- matrix(y, nrow(y))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # Root mean squared difference after a symmetric Procrustes rotation.
  pro_rmse <- function (x, y) procrustes_fit(unclass(x), unclass(y))$rmse

  # Run the C engine from given starting configurations.
  engine <- function (d, init, weakties = TRUE, stress = 1L, maxit = 200L, cpus = 1L) {
    n    <- attr(d, 'Size')
    k    <- dim(init)[2]
    prep <- .Call(C_nmds_prep, as.double(d), n)
    opts <- as.integer(c(maxit, weakties, stress, cpus, 0, 0, 1))
    res  <- .Call(C_nmds_run, prep, n, k, init, opts, c(1e-4, 1e-7, 0.999999))
    res$points <- array(res$points, dim(init))
    res
  }



  # Each start matches vegan::monoMDS() from the same configuration. ====

  # Built with the same compiler and optimization as vegan, the results are
  # identical. Elsewhere, last-bit differences in library functions such as
  # sqrt() can grow over hundreds of iterations, e.g. stopping an iteration
  # sooner or later, so long runs are only compared approximately.

  for (cs in ex$monomds) {

    d    <- dists[[cs$data]]
    init <- cs$init

    # Up to 25 iterations.
    res <- engine(d, init, cs$weakties, cs$stress, maxit = 25L)
    vg  <- cs$maxit25
    expect_equal(res$stress, vg$stress, tolerance = 1e-8)
    expect_equal(res$iters,  vg$iters)
    expect_equal(res$icause, vg$icause)
    for (i in seq_len(dim(init)[3]))
      expect_true(pro_rmse(res$points[, , i], vg$points[, , i]) < 1e-6)

    # Up to 200 iterations.
    res <- engine(d, init, cs$weakties, cs$stress, maxit = 200L)
    expect_equal(res$stress, cs$maxit200$stress, tolerance = 0.01)
  }

  # Few iterations: no chance for rounding differences to grow.
  vg  <- ex$monomds5
  res <- engine(d_bray, array(vg$init, c(dim(vg$init), 1)), maxit = 5L)
  expect_equal(res$stress, vg$stress, tolerance = 1e-12)
  expect_equal(res$iters,  vg$iters)
  expect_equal(res$icause, vg$icause)
  expect_equal(res$points[, , 1], vg$points, tolerance = 1e-12)



  # Same solution as vegan::metaMDS() and MASS::isoMDS(). ====

  ord <- nmds(d_pln)
  vg  <- ex$metamds_plane # set.seed(2); metaMDS(d_pln, trace = 0)
  expect_equal(attr(ord, 'stress'), vg$stress, tolerance = 1e-4)
  expect_true(coord_diff(ord, vg$points) < 0.01)

  iso <- ex$isomds_plane # isoMDS(d_pln, trace = FALSE)
  expect_equal(attr(ord, 'stress'), iso$stress / 100, tolerance = 0.01)
  expect_true(pro_rmse(ord, iso$points) < 0.01)

  # With a `maxdist` attribute: half-change scaling.
  d   <- d_pln / max(d_pln)
  attr(d, 'maxdist') <- 1
  ord <- nmds(d)
  vg  <- ex$metamds_hc # set.seed(2); metaMDS(d, trace = 0)
  expect_equal(attr(ord, 'stress'), vg$stress, tolerance = 1e-3)
  expect_true(coord_diff(ord, vg$points) < 0.05)



  # Return value. ====

  ord <- nmds(d_bray, cpus = 1)
  expect_identical(dim(ord), c(60L, 2L))
  expect_identical(dimnames(ord), list(attr(d_bray, 'Labels'), c('NMDS1', 'NMDS2')))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'stress', 'converged', 'tries', 'best_try') )
  expect_identical(attr(ord, 'method'), 'nmds')
  expect_identical(attr(ord, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(ord, 'prop'), c(NA_real_, NA_real_))
  expect_true(is.double(attr(ord, 'stress')) && attr(ord, 'stress') > 0)
  expect_true(is.logical(attr(ord, 'converged')))
  expect_true(attr(ord, 'tries') >= 20L)
  expect_true(attr(ord, 'best_try') %in% 0:attr(ord, 'tries'))

  # Centered, on principal axes, and scaled to the largest distance.
  x <- unclass(ord); attributes(x) <- list(dim = dim(x))
  expect_equal(colMeans(x), c(0, 0), tolerance = 1e-12)
  expect_equal(cor(x)[1, 2], 0, tolerance = 1e-10)
  expect_equal(max(dist(x)), max(d_bray), tolerance = 1e-12)
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))

  # k = 3, and a symmetric matrix input.
  ord3 <- nmds(d_bray, k = 3, trymax = 5)
  expect_identical(colnames(ord3), paste0('NMDS', 1:3))
  expect_identical(nmds(as.matrix(d_bray), trymax = 5), nmds(d_bray, trymax = 5))

  # Strong ties and stress formula 2.
  expect_true(attr(nmds(d_jacc, weakties = FALSE, trymax = 5), 'stress') > 0)
  d <- dist(c(0, 1, 2, 5, 11, 19)) # Two distances of 1: a tied pair.
  expect_true(attr(suppressWarnings(nmds(d, weakties = FALSE, trymax = 2)), 'stress') < 1e-3)
  expect_true(attr(nmds(d_bray, stress = 2, trymax = 5), 'stress') > 0)



  # Identical across runs and thread counts. ====

  x <- nmds(d_bray, cpus = 1)
  expect_identical(nmds(d_bray, cpus = 1), x)
  expect_identical(nmds(d_bray, cpus = 3), x)
  expect_identical(nmds(d_bray, cpus = 4), x)

  # Starts beyond `try`, run speculatively in batches of `cpus`.
  x <- nmds(d_bray, try = 2, trymax = 30, cpus = 1)
  expect_identical(nmds(d_bray, try = 2, trymax = 30, cpus = 4), x)

  # Stopping partway through a batch.
  x <- nmds(d_pln, try = 0, trymax = 10, cpus = 1)
  expect_true(attr(x, 'converged'))
  expect_true(attr(x, 'tries') < 4L)
  expect_identical(nmds(d_pln, try = 0, trymax = 10, cpus = 4), x)

  # R's random number generator is neither used nor changed.
  set.seed(5); r1 <- runif(1)
  set.seed(5); x  <- nmds(d_jacc, trymax = 3); r2 <- runif(1)
  expect_identical(r1, r2)

  # The seed changes the random starts.
  expect_false(identical(
    attr(nmds(d_bray, try = 5, trymax = 5, seed = 1), 'stress'),
    attr(nmds(d_bray, try = 5, trymax = 5, seed = 2), 'stress') ))

  # No random starts.
  ord <- nmds(d_bray, trymax = 0)
  expect_identical(attr(ord, 'tries'), 0L)
  expect_identical(attr(ord, 'best_try'), 0L)



  # Edge cases. ====

  # Perfect fit: points in a plane.
  pts <- cbind(x = sin(1:30), y = cos(1:30 / 2) * 3)
  expect_warning(ord <- nmds(dist(pts), trymax = 3), 'nearly')
  expect_true(attr(ord, 'stress') < 1e-3)

  # All distances tied.
  expect_warning(ord <- nmds(as.dist(matrix(1, 10, 10)), trymax = 2), 'nearly')
  expect_identical(dim(ord), c(10L, 2L))

  # Exactly equal configuration distances among tied dissimilarities: the
  # corners of a square.
  sq <- dist(cbind(c(0, 1, 0, 1), c(0, 0, 1, 1)))
  expect_warning(ord <- nmds(sq, trymax = 2), 'nearly')
  expect_equal(unname(as.matrix(dist(ord))), unname(as.matrix(sq)), tolerance = 1e-6)

  # One, two, three, and k + 1 samples.
  ord <- nmds(dist(matrix(1:3, 1, dimnames = list('A', NULL))))
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(rownames(ord), 'A')
  expect_identical(attr(ord, 'stress'), 0)
  expect_identical(names(attributes(ord))[3:5], c('method', 'eig', 'prop'))

  expect_warning(ord <- nmds(dist(c(a = 1, b = 3))), 'nearly')
  expect_identical(dimnames(ord), list(c('a', 'b'), 'NMDS1'))
  expect_equal(abs(unname(ord[, 1])), c(1, 1))

  expect_warning(ord <- nmds(bray(counts[1:3, ])), 'nearly')
  expect_identical(dim(ord), c(3L, 2L))

  expect_warning(ord <- nmds(bray(counts), k = 3), 'nearly')
  expect_identical(dim(ord), c(4L, 3L))

  # k larger than n - 1, and integer storage.
  d <- as.dist(matrix(c(0L, 3L, 4L, 3L, 0L, 5L, 4L, 5L, 0L), 3))
  expect_identical(
    suppressWarnings(nmds(d, k = 10)),
    suppressWarnings(nmds(d, k = 2)) )

  # Negative distances.
  d <- d_bray; d[1] <- -0.1
  expect_warning(nmds(d, trymax = 1), 'negative')

  # Half-change scaling's warnings.
  d <- d_bray
  attr(d, 'maxdist') <- 1e-10
  expect_warning(nmds(d, trymax = 1), 'dubious')
  attr(d, 'maxdist') <- 0.2
  expect_warning(nmds(d, trymax = 1), 'too few points')



  # Invalid input. ====

  d_na <- d_bray; d_na[2] <- NA
  expect_error(nmds(d_na), 'NA')
  d_na[2] <- Inf
  expect_error(nmds(d_na), 'infinite')
  expect_error(nmds(as.dist(matrix(0, 5, 5))), 'greater than 0')

  expect_error(nmds(d_bray, k = 0),         '`k`')
  expect_error(nmds(d_bray, k = 1.5),       '`k`')
  expect_error(nmds(d_bray, k = NA),        '`k`')
  expect_error(nmds(d_bray, try = -1),      '`try`')
  expect_error(nmds(d_bray, trymax = 'a'),  '`trymax`')
  expect_error(nmds(d_bray, maxit = 1:2),   '`maxit`')
  expect_error(nmds(d_bray, smin = -1),     '`smin`')
  expect_error(nmds(d_bray, sfgrmin = NA),  '`sfgrmin`')
  expect_error(nmds(d_bray, sratmax = 2),   '`sratmax`')
  expect_error(nmds(d_bray, weakties = NA), '`weakties`')
  expect_error(nmds(d_bray, stress = 3),    '`stress`')
  expect_error(nmds(d_bray, seed = 1.5),    'seed')
  expect_error(nmds(d_bray, cpus = 0),      'cpus')

  expect_error(nmds(matrix(1:6, 2)),           'symmetric')
  expect_error(nmds(matrix(c(0, 1, 2, 0), 2)), 'symmetric')
  expect_error(nmds(structure(1:2, class = 'dist', Size = 3L)), 'valid dist')

#})
