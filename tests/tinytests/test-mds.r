#test_that("mds", {

  home <- at_home()
  many <- if (home) 4L else 2L  # threads for direct .Call()s

  # Reproducible test data.
  set.seed(1)
  rand_counts <- function (n, p = 40) {
    m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
    m[, 1] <- m[, 1] + 1
    dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
    m
  }

  d_bray <- bray(rand_counts(60))
  d_jacc <- jaccard(rand_counts(50, 8))  # Many tied distances.
  d_big  <- bray(rand_counts(150))


  # The power of 2 that mds() reads dm at.
  dscale_of <- function (d) 2^-ceiling(log2(max(d)))

  # The C engine in metric mode, from given starting configurations.
  engine <- function (
      d, init = NULL, k = dim(init)[2], n_starts = 0L, max_iters = 300L, tol = 1e-6,
      cpus = 1L, seed = 0L, mem = 1024L, relax = 0L, check = 2^26,
      parmin = 2^18, dscale = 1 ) {
    .Call(
      C_nmds, as.double(d), attr(d, 'Size'), as.integer(k), init,
      as.integer(c(n_starts, seed, max_iters, cpus, 1L, 0L, mem, relax, 1L, check, parmin)),
      c(tol, dscale) )
  }

  # mds()'s PCoA start, as an n x k x 1 array in the units of dm * dscale.
  pcoa_init <- function (d, k) {
    n <- attr(d, 'Size')
    array(.Call(C_pcoa, as.double(d), n, as.integer(k), 1L)$points * dscale_of(d), c(n, k, 1))
  }

  # Stress as mds() reports it, and scikit-learn's normalized stress.
  stress_r <- function (d, x) {
    dx <- as.double(dist(x))
    sqrt(sum((dx - as.double(d))^2) / sum(as.double(d)^2))
  }
  skl_stress <- function (d, x) {
    dx <- as.double(dist(x))
    sqrt(sum((dx - as.double(d))^2) / sum(dx^2))
  }

  # The Guttman transform, in plain R.
  guttman_r <- function (x, d) {
    x <- matrix(unclass(x), nrow(x))
    n <- nrow(x); dx <- as.double(dist(x))
    R <- matrix(0, n, n); R[lower.tri(R)] <- ifelse(dx > 0, as.double(d) / dx, 0); R <- R + t(R)
    B <- -R; diag(B) <- rowSums(R)
    B %*% x / n
  }

  # Metric SMACOF in plain R, written from scikit-learn's description: the
  # per-iteration raw stress, and the final configuration.
  smacof_r <- function (d, init, max_iters = 300L, tol = 1e-6) {
    dv <- as.double(d)
    x  <- init; trace <- numeric(0); old <- NA
    for (it in seq_len(max_iters)) {
      x  <- guttman_r(x, d)
      dx <- as.double(dist(x))
      st <- sum((dx - dv)^2)
      trace <- c(trace, st)
      if (it > 1 && (old - st) / sum(dx^2) <= tol) break
      old <- st
    }
    list(x = x, trace = trace)
  }

  # Root mean squared difference after a Procrustes rotation, with both
  # configurations centered (not rescaled).
  pro_rmse <- function (x, y) {
    x  <- scale(unclass(x)[, , drop = FALSE], scale = FALSE)
    y  <- scale(unclass(y)[, , drop = FALSE], scale = FALSE)
    sv <- svd(crossprod(x, y))
    sqrt(sum((x - y %*% sv$v %*% t(sv$u))^2) / nrow(x))
  }

  rel_diff <- function (x, y) max(abs(x - y)) / max(abs(y))



  # Stored results from scikit-learn, so at home only. ====
  # Same iterations as scikit-learn's metric smacof(), from the same
  # starting configuration. From data-raw/expected/mds.r.

  if (home) {

    fx <- readRDS('expected/mds.rds')

    for (i in seq_len(nrow(fx$cases))) {

      cs  <- fx$cases[i, ]
      fv  <- fx$values[[i]]
      d   <- fx$dists[[cs$dist]]
      res <- engine(d, array(fv$init, c(dim(fv$init), 1)))

      expect_identical(res$iters, cs$n_iter, info = cs$case)
      expect_true(max(abs(res$trace[fv$trace$iter] / fv$trace$stress - 1)) < 1e-10, info = cs$case)
      expect_true(rel_diff(res$x_all[, , 1], fv$x) < 1e-8, info = cs$case)

      # The reported stress, normalized by the dissimilarities.
      expect_equal(res$stress, stress_r(d, res$x_all[, , 1]), tolerance = 1e-10)

      # Final quality: no worse than scikit-learn's own best, from 8 random
      # starts or from classical MDS, by the same (scale-safe) measure.
      if (cs$k >= 2) {
        ord <- mds(d, k = cs$k, cpus = 1)
        expect_true(attr(ord, 'stress') <= min(cs$skl_random, cs$skl_cmds) * (1 + 1e-3), info = cs$case)
      }
    }

    # The plain R version agrees too.
    d    <- fx$dists$bray30
    init <- fx$values[[1]]$init
    ref  <- smacof_r(d, init, max_iters = 60L)
    res  <- engine(d, array(init, c(dim(init), 1)), max_iters = 60L)
    expect_true(max(abs(res$trace / ref$trace - 1)) < 1e-10)
  }



  # The plain R version agrees, from a random start. ====

  init <- array(runif(50 * 2), c(50, 2, 1))
  ref  <- smacof_r(d_jacc, init[, , 1], max_iters = 40L)
  res  <- engine(d_jacc, init, max_iters = 40L)
  expect_identical(res$iters, length(ref$trace))
  expect_true(max(abs(res$trace[seq_len(res$iters)] / ref$trace - 1)) < 1e-10)
  expect_true(rel_diff(res$x_all[, , 1], ref$x) < 1e-8)



  # The returned configuration is a fixed point of the Guttman ====
  # transform, at the right scale; its stress is the one plain steps reach
  # from the same start. Tight tolerances take a while, so for k = 1, 2, 3
  # at home.

  for (k in if (home) 1:3 else 2L) {
    ord <- mds(d_bray, k = k, n_starts = 0, tol = 1e-12, max_iters = 1e5, cpus = 1)
    x   <- matrix(unclass(ord), 60)
    g   <- guttman_r(x, d_bray)
    expect_true(attr(ord, 'converged'))
    expect_true(sqrt(sum((x - g)^2) / sum(x^2)) < 1e-5, info = k)
    dx  <- as.double(dist(x))
    expect_equal(sum(dx * d_bray) / sum(dx^2), 1, tolerance = 1e-8, info = k)
    res <- engine(d_bray, pcoa_init(d_bray, k), tol = 1e-12, max_iters = 100000L,
                  dscale = dscale_of(d_bray))
    expect_equal(attr(ord, 'stress'), res$stress, tolerance = 1e-8, info = k)
  }

  # With default settings too.
  for (k in 1:3) {
    ord <- mds(d_bray, k = k, cpus = 1)
    x   <- matrix(unclass(ord), 60)
    dx  <- as.double(dist(x))
    expect_true(sqrt(sum((x - guttman_r(x, d_bray))^2) / sum(x^2)) < 1e-3, info = k)
    expect_equal(sum(dx * d_bray) / sum(dx^2), 1, tolerance = 1e-6, info = k)
  }



  # Euclidean distances are recovered. ====

  for (k in 1:3) {
    pts <- matrix(rnorm(40 * k), 40, k, dimnames = list(paste0('P', 1:40), NULL))
    d   <- dist(pts)
    ord <- mds(d, k = k, cpus = 1)
    expect_true(attr(ord, 'stress') < 1e-8)
    expect_true(pro_rmse(ord, pts) < 1e-6)
    expect_true(pro_rmse(ord, pcoa(d, k = k)) < 1e-6)
    expect_equal(as.double(dist(ord)), as.double(d), tolerance = 1e-8)
    # Every start reaches stress ~1e-16 here, so which one wins is a tie
    # broken by rounding, and it differs between platforms (musl, -O0).
    expect_true(attr(ord, 'best_start') %in% 0:4)
  }



  # Return value. ====

  ord <- mds(d_bray, cpus = 1)
  expect_identical(dim(ord), c(60L, 2L))
  expect_identical(dimnames(ord), list(attr(d_bray, 'Labels'), c('MDS1', 'MDS2')))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'stress', 'converged',
      'n_iters', 'n_starts', 'best_start', 'n_best', 'class') )
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(attr(ord, 'method'), 'mds')
  expect_identical(attr(ord, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(ord, 'prop'), c(NA_real_, NA_real_))
  expect_true(attr(ord, 'stress') > 0 && attr(ord, 'stress') < 1)
  expect_true(is.integer(attr(ord, 'n_iters')) && attr(ord, 'n_iters') >= 1L)
  expect_true(isTRUE(attr(ord, 'converged')))
  expect_identical(attr(ord, 'n_starts'), 4L)
  expect_equal(attr(ord, 'stress'), stress_r(d_bray, ord), tolerance = 1e-10)

  # scikit-learn's normalized stress, from the result.
  expect_true(skl_stress(d_bray, ord) > 0)

  # The best start is the one reported.
  init <- pcoa_init(d_bray, 2)
  res  <- engine(d_bray, init, n_starts = 4L, max_iters = 1000L, tol = 1e-7,
                 relax = 1L, dscale = dscale_of(d_bray))
  expect_identical(attr(ord, 'stress'), min(res$stress_all))
  expect_identical(attr(ord, 'best_start'), which.min(res$stress_all) - 1L)
  expect_identical(attr(ord, 'n_iters'), res$iters_all[[which.min(res$stress_all)]])
  expect_identical(attr(ord, 'n_best'), sum(res$stress_all <= min(res$stress_all) * 1.001))

  # Centered and on principal axes, in dm's units, signs fixed.
  x <- unclass(ord); attributes(x) <- list(dim = dim(x))
  expect_equal(colMeans(x), c(0, 0), tolerance = 1e-12)
  expect_equal(cor(x)[1, 2], 0, tolerance = 1e-10)
  expect_true(var(x[, 1]) >= var(x[, 2]))
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))
  expect_equal(pro_rmse(x, res$x_all[, , which.min(res$stress_all)] / dscale_of(d_bray)), 0, tolerance = 1e-10)

  # Same as through ordination(); k = 3, 4; matrix input.
  expect_identical(ordination(d_bray, 'mds', cpus = 1), ord)
  expect_identical(ordination(d_bray, 'MDS', cpus = 1), ord)
  ord3 <- mds(d_bray, k = 3, n_starts = 1)
  ord4 <- mds(d_bray, k = 4, n_starts = 1)
  expect_identical(colnames(ord4), paste0('MDS', 1:4))
  expect_true(attr(ord4, 'stress') < attr(ord3, 'stress'))
  expect_true(attr(ord3, 'stress') < attr(ord, 'stress'))
  expect_identical(mds(as.matrix(d_bray), n_starts = 1), mds(d_bray, n_starts = 1))

  # max_iters = 0, small max_iters, and tol = 0.
  ord0 <- mds(d_bray, max_iters = 0, n_starts = 0)
  expect_identical(attr(ord0, 'n_iters'), 0L)
  expect_equal(attr(ord0, 'stress'), stress_r(d_bray * dscale_of(d_bray), init[, , 1]), tolerance = 1e-10)
  ord <- mds(d_jacc, tol = 0, max_iters = 5000, n_starts = 0)
  expect_true(attr(ord, 'converged'))
  ord <- mds(d_jacc, max_iters = 2, n_starts = 0)
  expect_identical(attr(ord, 'n_iters'), 2L)
  expect_false(attr(ord, 'converged'))



  # Results scale with dm. ====

  for (p in c(-900, -1, 3, 400)) {
    x <- mds(d_bray * 2^p, n_starts = 1)
    y <- mds(d_bray, n_starts = 1)
    expect_identical(unclass(x)[, ], unclass(y)[, ] * 2^p)
    expect_identical(attr(x, 'stress'), attr(y, 'stress'))
  }
  for (s in c(1e-300, 1e-160, 1e100)) {
    x <- mds(d_bray * s, n_starts = 1)
    y <- mds(d_bray, n_starts = 1)
    expect_equal(unclass(x)[, ] / s, unclass(y)[, ], tolerance = 1e-6)
    expect_equal(attr(x, 'stress'), attr(y, 'stress'), tolerance = 1e-6)
  }



  # Identical across runs and thread counts (and, in test-nmds.r, ====
  # sessions).

  for (args in list(list(dm = d_bray), list(dm = d_jacc, k = 3))) {
    x <- do.call(mds, c(args, cpus = 1))
    expect_identical(do.call(mds, c(args, cpus = 1)), x)
    for (cpus in c(2, 3, 4, 100))
      expect_identical(do.call(mds, c(args, cpus = cpus)), x)
  }

  init <- pcoa_init(d_big, 2)
  r1   <- engine(d_big, init, n_starts = 3L, max_iters = 25L, relax = 1L)
  for (cpus in unique(c(2L, many))) for (mem in c(1L, 1024L)) for (check in c(2^26, 1000))
    expect_identical(
      engine(d_big, init, n_starts = 3L, max_iters = 25L, relax = 1L, cpus = cpus,
             mem = mem, check = check, parmin = 1000),
      r1 )

  set.seed(5); r1 <- runif(1)
  set.seed(5); x  <- mds(d_jacc, n_starts = 2); r2 <- runif(1)
  expect_identical(r1, r2)
  set.seed(5);  x1 <- mds(d_jacc, n_starts = 2)
  set.seed(99); x2 <- mds(d_jacc, n_starts = 2)
  expect_identical(x1, x2)

  r1 <- engine(d_bray, k = 2L, n_starts = 2L, seed = 1L)
  r2 <- engine(d_bray, k = 2L, n_starts = 2L, seed = 2L)
  expect_false(identical(r1$stress_all, r2$stress_all))



  # Edge cases. ====

  # All distances tied.
  ord <- mds(as.dist(matrix(1, 10, 10)))
  expect_identical(dim(ord), c(10L, 2L))
  expect_true(all(is.finite(ord)) && attr(ord, 'stress') < 1)

  # One nonzero distance among 20 samples: the stress is at most 1.
  d <- as.dist(matrix(0, 20, 20)); d[1] <- 1
  ord <- mds(d)
  expect_true(attr(ord, 'stress') <= 1)
  expect_equal(attr(ord, 'stress'), stress_r(d, ord), tolerance = 1e-10)

  # Duplicate samples: zero distances, pulled together.
  m   <- rand_counts(20)
  m   <- rbind(m, D1 = m[1, ], D2 = m[1, ])
  ord <- mds(bray(m))
  expect_true(all(is.finite(ord)))
  expect_true(max(abs(ord['D1', ] - ord[1, ])) < 0.05 * max(abs(ord)))

  # Zero, one, two, three samples; k larger than n - 1.
  ord <- mds(as.dist(matrix(numeric(0), 0, 0)))
  expect_identical(dim(ord), c(0L, 0L))
  ord <- mds(dist(matrix(1:3, 1, dimnames = list('A', NULL))))
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(attr(ord, 'method'), 'mds')
  expect_identical(attr(ord, 'stress'), 0)
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'stress', 'converged',
      'n_iters', 'n_starts', 'best_start', 'n_best', 'class') )
  ord <- mds(dist(c(a = 1, b = 4)))
  expect_identical(dimnames(ord), list(c('a', 'b'), 'MDS1'))
  expect_equal(abs(unname(ord[, 1])), c(1.5, 1.5))
  expect_equal(attr(ord, 'stress'), 0)
  d <- as.dist(matrix(c(0L, 3L, 4L, 3L, 0L, 5L, 4L, 5L, 0L), 3))
  ord <- mds(d, k = 1e12)
  expect_identical(dim(ord), c(3L, 2L))
  expect_equal(as.double(dist(ord)), c(3, 4, 5), tolerance = 1e-6)
  expect_identical(mds(d, cpus = 100), mds(d, cpus = 1))
  expect_identical(mds(d_bray, max_iters = 1e12, n_starts = 0), mds(d_bray, max_iters = 1e6, n_starts = 0))

  # A start whose points all collapse together: stress 1.
  d    <- as.dist(matrix(c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0), 4))
  init <- array(c(0, 1, 0, 1, 0, 1, 0, 1), c(4, 2, 1))
  res  <- engine(d, init)
  expect_identical(res$stress, 1)
  expect_true(all(is.finite(res$points)))



  # Invalid input. ====

  d_na <- d_bray; d_na[2] <- NA
  expect_error(mds(d_na), 'NA')
  d_na[2] <- Inf
  expect_error(mds(d_na), 'infinite')
  d_na[2] <- -0.1
  expect_error(mds(d_na), 'negative')
  expect_error(mds(as.dist(matrix(0, 5, 5))), 'greater than 0')
  expect_error(mds(d_bray * 1e151), 'too large')

  expect_error(mds(d_bray, k = 0),            '`k`')
  expect_error(mds(d_bray, k = 2.5),          '`k`')
  expect_error(mds(d_bray, k = c(2, 3)),      '`k`')
  expect_error(mds(d_bray, k = 'a'),          '`k`')
  expect_error(mds(d_bray, n_starts = -1),    '`n_starts`')
  expect_error(mds(d_bray, n_starts = 1001),  '`n_starts`')
  expect_error(mds(d_bray, n_starts = NA),    '`n_starts`')
  expect_error(mds(d_bray, max_iters = -1),    '`max_iters`')
  expect_error(mds(d_bray, max_iters = Inf),   '`max_iters`')
  expect_error(mds(d_bray, tol = -1),         '`tol`')
  expect_error(mds(d_bray, tol = NaN),        '`tol`')
  expect_error(mds(d_bray, seed = 'a'),       'seed')
  expect_error(mds(d_bray, cpus = 0),         'cpus')
  expect_error(mds(d_bray, ties = 'primary'), 'unused argument')

  expect_error(mds(matrix(1:6, 2)),           'symmetric')
  expect_error(mds(matrix(c(0, 1, 2, 0), 2)), 'symmetric')
  expect_error(mds('a'),                      'symmetric')
  expect_error(mds(NULL),                     'symmetric')
  expect_error(mds(structure(1:2, class = 'dist', Size = 3L)), 'valid dist')

#})
