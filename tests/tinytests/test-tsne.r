#test_that("tsne", {

  # Test data and expected results, from data-raw/expected/tsne.r.
  ex <- readRDS('expected/tsne.rds')
  d  <- dist(ex$pts)

  run <- function (a, n_iter = a$maxit, cpus = 1L)
    tsne(
      d, k = a$k, perplexity = a$perplexity, theta = a$theta, n_iter = n_iter,
      init = ex$inits[[a$init]], stop_lying_iter = a$stop_lying_iter,
      mom_switch_iter = a$mom_switch_iter, cpus = cpus )

  # Reproducible test data.
  set.seed(1)
  rand_counts <- function (n, p = 30) {
    m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.5), n, p)
    m[, 1] <- m[, 1] + 1
    dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
    m
  }
  m40 <- rand_counts(40)
  b40 <- bray(m40)



  # Same coordinates (up to each axis's sign) and cost as the reference,
  # from the same start: stored results, so at home only. ====

  # `y`'s axes flipped to match `x`'s signs.
  match_signs <- function (x, y) {
    x <- matrix(x, nrow(x))
    sweep(y, 2, sign(colSums(x * y)), '*')
  }

  if (at_home()) {

    # The first 20 iterations are compared everywhere. (Other math libraries
    # and fused multiply-adds change the result by up to 2e-7 there.)
    for (a in ex$cases) {
      x <- run(a, n_iter = 20L)
      expect_equal(matrix(x, nrow(x)), match_signs(x, a$Y20), tolerance = 1e-6)
      expect_equal(attr(x, 'kl_divergence'), a$cost20, tolerance = 1e-6)
    }

    # All iterations. t-SNE amplifies rounding differences, so after many
    # iterations it can settle in a different local minimum. Barnes-Hut
    # computes everything in the reference's order, which gives exactly its
    # result on the platform the expected results came from (Windows x86-64:
    # elsewhere, last-bit differences in exp() and log(), or fused
    # multiply-adds, lead to another minimum). Exact mode sums Q in
    # parallel, so only the quality of the result is compared: its cost,
    # and how many of each point's 10 nearest neighbours stay among its 10
    # nearest on the plot.
    nn_kept <- function (y) {
      dx <- as.matrix(d)
      dy <- as.matrix(dist(y))
      mean(sapply(seq_len(nrow(dx)), function (i)
        length(intersect(order(dx[i, ])[2:11], order(dy[i, ])[2:11])) / 10 ))
    }

    for (a in ex$cases) {
      x <- run(a)
      if (a$theta == 0) {
        expect_true(attr(x, 'kl_divergence') < 1.2 * a$cost)
        expect_true(nn_kept(x) > nn_kept(a$Y) - 0.05)
      }
      else if (.Platform$OS.type == 'windows' && identical(R.version$arch, 'x86_64')) {
        expect_equal(matrix(x, nrow(x)), match_signs(x, a$Y), tolerance = 1e-8)
        expect_equal(attr(x, 'kl_divergence'), a$cost, tolerance = 1e-8)
      }
    }

    # Barnes-Hut with early exaggeration on 600 points (a deep tree), in 3
    # dimensions (an octree), and at the largest perplexity, (n - 1) / 3.
    # The same comparisons.
    for (a in ex$more_cases) {

      pts  <- if (is.null(a$pts)) ex$pts else a$pts
      run2 <- function (n_iter) tsne(
        dist(pts), k = a$k, perplexity = a$perplexity, theta = a$theta,
        n_iter = n_iter, init = a$init, stop_lying_iter = a$stop_lying_iter,
        mom_switch_iter = a$mom_switch_iter, cpus = 2L )

      info <- paste(nrow(pts), a$k, a$perplexity)
      x    <- run2(20L)
      expect_identical(attr(x, 'perplexity'), a$perplexity, info = info)
      expect_equal(matrix(x, nrow(x)), match_signs(x, a$Y20), tolerance = 1e-6, info = info)
      expect_equal(attr(x, 'kl_divergence'), a$cost20, tolerance = 1e-6, info = info)

      if (.Platform$OS.type == 'windows' && identical(R.version$arch, 'x86_64')) {
        x <- run2(a$maxit)
        expect_equal(matrix(x, nrow(x)), match_signs(x, a$Y), tolerance = 1e-8, info = info)
        expect_equal(attr(x, 'kl_divergence'), a$cost, tolerance = 1e-8, info = info)
      }
    }
  }



  # Result format. ====

  x <- tsne(b40, perplexity = 5, n_iter = 50, cpus = 1)
  expect_true(is.matrix(x) && is.double(x))
  expect_identical(dim(x), c(40L, 2L))
  expect_identical(dimnames(x), list(paste0('S', 1:40), c('TSNE1', 'TSNE2')))
  expect_identical(
    names(attributes(x)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'kl_divergence', 'perplexity', 'init', 'class') )
  # Each axis's largest absolute coordinate is positive.
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))
  expect_identical(attr(x, 'method'), 'tsne')
  expect_identical(attr(x, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(x, 'prop'), c(NA_real_, NA_real_))
  expect_identical(attr(x, 'perplexity'), 5)
  expect_identical(attr(x, 'init'), 'pcoa')
  expect_true(attr(x, 'kl_divergence') > 0)
  expect_equal(unname(colMeans(x)), c(0, 0), tolerance = 1e-10)

  # A symmetric matrix works the same as a dist object.
  expect_identical(tsne(as.matrix(b40), perplexity = 5, n_iter = 50, cpus = 1), x)

  # Exact t-SNE allows more than 3 dimensions.
  x5 <- tsne(b40, k = 5, theta = 0, perplexity = 5, n_iter = 50, cpus = 1)
  expect_identical(colnames(x5), paste0('TSNE', 1:5))
  expect_identical(attr(x5, 'eig'), rep(NA_real_, 5))



  # Reproducibility. ====

  # Same seed, same result; another seed, another random start.
  xr <- tsne(b40, perplexity = 5, n_iter = 50, init = 'random', cpus = 1)
  expect_identical(attr(xr, 'init'), 'random')
  expect_identical(tsne(b40, perplexity = 5, n_iter = 50, init = 'random', cpus = 1), xr)
  expect_false(identical(tsne(b40, perplexity = 5, n_iter = 50, init = 'random', seed = 1, cpus = 1), xr))

  # The default PCoA start doesn't use the seed.
  expect_identical(tsne(b40, perplexity = 5, n_iter = 50, seed = 1, cpus = 1), x)

  # The seed has no effect when `init` is given.
  init <- pcoa(b40, k = 2)
  expect_identical(
    tsne(b40, perplexity = 5, n_iter = 50, init = init, seed = 1, cpus = 1),
    tsne(b40, perplexity = 5, n_iter = 50, init = init, seed = 2, cpus = 1) )

  # R's random number stream is untouched.
  set.seed(3); r1 <- runif(1)
  set.seed(3); invisible(tsne(b40, perplexity = 5, n_iter = 5, init = 'random', cpus = 1)); r2 <- runif(1)
  expect_identical(r1, r2)

  # Identical for any number of threads, in both modes.
  set.seed(4)
  pts  <- matrix(rnorm(1200 * 3), 1200, 3) + rep(c(0, 4), each = 600)
  d12  <- dist(pts)
  expect_identical(
    tsne(d12, perplexity = 20, n_iter = 30, cpus = 1),
    tsne(d12, perplexity = 20, n_iter = 30, cpus = 3) )
  d5 <- dist(pts[1:500, ])
  expect_identical(
    tsne(d5, theta = 0, perplexity = 20, n_iter = 10, cpus = 1),
    tsne(d5, theta = 0, perplexity = 20, n_iter = 10, cpus = 3) )




  # Starting coordinates. ====

  # "pcoa": PCoA axes scaled so that the first has standard deviation 1e-4,
  # with early exaggeration, unlike a matrix given as `init`.
  p  <- pcoa(b40, k = 2)
  p  <- unname(unclass(p)[, 1:2]) / sd(p[, 1]) * 1e-4
  xp <- tsne(b40, perplexity = 5, n_iter = 300, init = 'pcoa', cpus = 1)
  xm <- tsne(b40, perplexity = 5, n_iter = 300, init = p, cpus = 1)
  expect_identical(attr(xp, 'init'), 'pcoa')
  expect_identical(attr(xm, 'init'), 'matrix')
  expect_false(identical(c(xp), c(xm)))
  xm <- tsne(b40, perplexity = 5, n_iter = 300, init = p, stop_lying_iter = 250, mom_switch_iter = 250, cpus = 1)
  expect_identical(c(xp), c(xm))
  expect_identical(attr(xp, 'kl_divergence'), attr(xm, 'kl_divergence'))

  # Falls back to "random" without k positive PCoA axes.
  x <- tsne(b40, k = 5, theta = 0, perplexity = 5, n_iter = 20, init = 'pcoa', cpus = 1)
  expect_identical(attr(x, 'init'), 'pcoa')
  x <- tsne(bray(m40[1:4, ]), k = 4, theta = 0, n_iter = 20, init = 'pcoa', cpus = 1)
  expect_identical(attr(x, 'init'), 'random')
  expect_identical(c(x), c(tsne(bray(m40[1:4, ]), k = 4, theta = 0, n_iter = 20, init = 'random', cpus = 1)))

  # "auto" learning rate: max(n / exaggeration_factor, 200).
  expect_identical(
    tsne(b40, perplexity = 5, n_iter = 20, cpus = 1),
    tsne(b40, perplexity = 5, n_iter = 20, learning_rate = 200, cpus = 1) )
  expect_identical(
    tsne(d12, perplexity = 20, n_iter = 5, exaggeration_factor = 4, cpus = 1),
    tsne(d12, perplexity = 20, n_iter = 5, exaggeration_factor = 4, learning_rate = 300, cpus = 1) )
  expect_false(identical(
    tsne(d12, perplexity = 20, n_iter = 5, exaggeration_factor = 4, cpus = 1),
    tsne(d12, perplexity = 20, n_iter = 5, exaggeration_factor = 4, learning_rate = 200, cpus = 1) ))

  # A data frame works as a matrix.
  expect_identical(
    tsne(b40, perplexity = 5, n_iter = 20, init = as.data.frame(p), cpus = 1),
    tsne(b40, perplexity = 5, n_iter = 20, init = p, cpus = 1) )



  # Input similarities. ====

  # Each point's conditional distribution has the requested perplexity, and
  # P is symmetric and sums to 1.
  n   <- attr(b40, 'Size')
  aff <- .Call(C_tsne_affinity, as.double(b40), n, 5, TRUE, 1L)
  dm  <- unname(as.matrix(b40))
  P   <- matrix(0, n, n)
  P[lower.tri(P)] <- aff$P
  P   <- P + t(P)
  cond <- exp(-aff$beta * dm^2)
  diag(cond) <- 0
  cond <- cond / rowSums(cond)
  H <- -rowSums(ifelse(cond > 0, cond * log(cond), 0))
  expect_equal(H, rep(log(5), n), tolerance = 1e-4)
  expect_equal(sum(P), 1, tolerance = 1e-12)
  expect_equal(P, (cond + t(cond)) / (2 * n), tolerance = 1e-12)

  # Barnes-Hut: the same from each point's 15 nearest neighbors.
  aff  <- .Call(C_tsne_affinity, as.double(b40), n, 5, FALSE, 1L)
  rows <- rep(seq_len(n), diff(aff$row))
  Ps   <- matrix(0, n, n)
  Ps[cbind(rows, aff$col + 1L)] <- aff$val
  nn   <- t(apply(dm + diag(Inf, n), 1, function (r) order(r)[1:15]))
  cond <- matrix(0, n, n)
  for (i in seq_len(n)) {
    e <- exp(-aff$beta[[i]] * dm[i, nn[i, ]]^2)
    cond[i, nn[i, ]] <- e / sum(e)
  }
  expect_true(isSymmetric(Ps))
  expect_equal(sum(Ps), 1, tolerance = 1e-12)
  expect_equal(Ps, (cond + t(cond)) / (2 * n), tolerance = 1e-12)



  # Small inputs. ====

  # Default perplexity is reduced quietly; a given one with a warning.
  x <- tsne(bray(ecodive::ex_counts), n_iter = 50, cpus = 1)
  expect_identical(attr(x, 'perplexity'), 1)
  expect_warning(
    x <- tsne(bray(ecodive::ex_counts), perplexity = 30, n_iter = 50, cpus = 1),
    'reduced from 30 to 1' )
  expect_identical(attr(x, 'perplexity'), 1)

  one <- tsne(as.dist(matrix(0, 1, 1)), cpus = 1)
  expect_identical(dim(one), c(1L, 2L))
  expect_identical(attr(one, 'method'), 'tsne')
  expect_identical(attr(one, 'init'), 'none')
  expect_identical(attr(tsne(dist(1), init = matrix(5, 1, 2)), 'init'), 'matrix')
  expect_identical(dim(tsne(dist(matrix(numeric(0), 0, 2)), cpus = 1)), c(0L, 2L))

  # Distances of zero (duplicate samples) are fine.
  dup <- m40[c(1:20, 1:5), ]
  rownames(dup) <- NULL
  x   <- tsne(bray(dup), perplexity = 5, n_iter = 100, cpus = 1)
  expect_true(all(is.finite(x)))



  # Invalid input. ====

  expect_error(tsne(b40, k = 0),               '`k` must be a whole number')
  expect_error(tsne(b40, k = 4),               '`k` must be 1, 2, or 3')
  expect_error(tsne(b40, perplexity = 0),      '`perplexity` must be')
  expect_error(tsne(b40, perplexity = NA),     '`perplexity` must be')
  expect_error(tsne(b40, theta = 2),           '`theta` must be')
  expect_error(tsne(b40, n_iter = -1),          '`n_iter` must be')
  expect_error(tsne(b40, n_iter = 1.5),         '`n_iter` must be')
  expect_error(tsne(b40, learning_rate = 0),   '`learning_rate` must be "auto" or a number')
  expect_error(tsne(b40, learning_rate = 'a'), '`learning_rate` must be "auto" or a number')
  expect_error(tsne(b40, momentum = -1),       '`momentum` must be')
  expect_error(tsne(b40, final_momentum = 'a'), '`final_momentum` must be')
  expect_error(tsne(b40, mom_switch_iter = -1), '`mom_switch_iter` must be')
  expect_error(tsne(b40, exaggeration_factor = 0), '`exaggeration_factor` must be')
  expect_error(tsne(b40, stop_lying_iter = NA), '`stop_lying_iter` must be')
  expect_error(tsne(b40, init = matrix(0, 40, 3)), '`init` must be a numeric matrix with 40 rows and 2 columns')
  expect_error(tsne(b40, init = matrix(NA, 40, 2)), '`init` must be')
  expect_error(tsne(b40, init = 'spectral'), '`init` must be "random", "pcoa", or a numeric matrix')
  expect_error(tsne(b40, init = NULL), '`init` must be')
  expect_error(tsne(b40, seed = 1.5))
  expect_error(tsne(b40, cpus = 0))

  bad <- b40
  bad[[1]] <- NA
  expect_error(tsne(bad), 'must not contain NA')
  bad[[1]] <- -1
  expect_error(tsne(bad), 'negative')

  # Regressions: these crashed R. Below a perplexity of 1/3, Barnes-Hut had
  # no neighbors per point; with momentum >= 1 or huge steps, coordinates
  # became NaN and the Barnes-Hut tree subdivided without end.
  expect_error(tsne(b40, perplexity = 0.3),    '`perplexity` must be a number of at least 1')
  expect_error(tsne(b40, perplexity = 0.9, theta = 0), '`perplexity` must be a number of at least 1')
  expect_error(tsne(b40, momentum = 1),        '`momentum` must be a number from 0 to less than 1')
  expect_error(tsne(b40, final_momentum = 3),  '`final_momentum` must be a number from 0 to less than 1')
  expect_error(tsne(b40, learning_rate = Inf), '`learning_rate` must be "auto" or a number')
  expect_error(tsne(b40, perplexity = 5, n_iter = 50, learning_rate = 1e300), 't-SNE diverged')
  expect_error(tsne(b40, perplexity = 5, n_iter = 50, learning_rate = 1e300, theta = 0), 't-SNE diverged')
  expect_error(.Call(C_tsne_affinity, as.double(b40), 40L, 0.3, FALSE, 1L), 'perplexity of at least 1/3')
  # The smallest perplexity that `n` samples allow still runs.
  expect_identical(dim(tsne(dist(1:2), cpus = 1)), c(2L, 2L))

#})
