#test_that("nmds", {

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

  # Points near a plane, with noise: a clear NMDS solution.
  plane <- function (n) {
    pts <- cbind(runif(n), runif(n) * 0.6, rnorm(n, sd = 0.02))
    rownames(pts) <- paste0('P', seq_len(n))
    dist(pts)
  }

  d_bray <- bray(rand_counts(60))
  d_jacc <- jaccard(rand_counts(50, 8))  # Many tied distances.
  d_pln  <- plane(40)
  d_big  <- bray(rand_counts(150))       # Sorted in several buckets.


  # The power of 2 that nmds() reads dm at.
  dscale_of <- function (d) 2^-ceiling(log2(max(d)))

  # The C engine, from given starting configurations (n x k x m array).
  engine <- function (
      d, init = NULL, k = dim(init)[2], n_starts = 0L, max_iters = 300L, tol = 1e-6,
      ties = 'secondary', cpus = 1L, seed = 0L, mem = 1024L, relax = 0L,
      metric = 0L, check = 2^26, parmin = 2^18, dscale = 1 ) {
    .Call(
      C_nmds, as.double(d), attr(d, 'Size'), as.integer(k), init,
      as.integer(c(n_starts, seed, max_iters, cpus, 1L, ties == 'primary', mem,
                   relax, metric, check, parmin)),
      c(tol, dscale) )
  }

  # nmds()'s PCoA start, as an n x k x 1 array in the units of dm * dscale.
  pcoa_init <- function (d, k) {
    n <- attr(d, 'Size')
    array(.Call(C_pcoa, as.double(d), n, as.integer(k), 1L)$points * dscale_of(d), c(n, k, 1))
  }

  # Weighted pool-adjacent-violators, in plain R.
  pava_r <- function (y, w = rep(1, length(y))) {
    v <- y; ww <- w; len <- rep(1L, length(y)); top <- 0L
    for (i in seq_along(y)) {
      top <- top + 1L; v[top] <- y[i]; ww[top] <- w[i]; len[top] <- 1L
      while (top > 1L && v[top - 1L] >= v[top]) {
        sw <- ww[top - 1L] + ww[top]
        v[top - 1L]   <- (ww[top - 1L] * v[top - 1L] + ww[top] * v[top]) / sw
        ww[top - 1L]  <- sw
        len[top - 1L] <- len[top - 1L] + len[top]
        top <- top - 1L
      }
    }
    rep(v[seq_len(top)], len[seq_len(top)])
  }

  # Isotonic regression of distances `dx` on dissimilarities `dv` (both in
  # dist order), with zero dissimilarities fixed at 0.
  iso_r <- function (dv, dx, ties) {
    dhat <- numeric(length(dv))
    idx  <- which(dv != 0)
    if (ties == 'secondary') {
      o   <- idx[order(dv[idx])]
      g   <- cumsum(c(TRUE, diff(dv[o]) != 0))
      fit <- pava_r(as.double(tapply(dx[o], g, mean)), tabulate(g))
      dhat[o] <- fit[g]
    } else {
      o <- idx[order(dv[idx], dx[idx])]
      dhat[o] <- pava_r(dx[o])
    }
    dhat
  }

  # The Guttman transform of x for disparities dhat (dist order).
  guttman_r <- function (x, dhat) {
    n <- nrow(x)
    dx <- as.double(dist(x))
    R <- matrix(0, n, n); R[lower.tri(R)] <- ifelse(dx > 0, dhat / dx, 0); R <- R + t(R)
    B <- -R; diag(B) <- rowSums(R)
    B %*% x / n
  }

  # Non-metric SMACOF in plain R, written from scikit-learn's description:
  # the per-iteration raw stress, and the final configuration.
  smacof_r <- function (d, init, max_iters = 300L, tol = 1e-6, ties = 'secondary') {
    dv <- as.double(d); N <- length(dv)
    x  <- init; trace <- numeric(0); old <- NA
    for (it in seq_len(max_iters)) {
      dx   <- as.double(dist(x))
      dhat <- if (it == 1) dv else iso_r(dv, dx, ties)
      dhat <- dhat * sqrt(N / sum(dhat^2))
      x    <- guttman_r(x, dhat)
      dx   <- as.double(dist(x))
      st   <- sum((dx - dhat)^2)
      trace <- c(trace, st)
      if (it > 1 && (old - st) / sum(dx^2) <= tol) break
      old <- st
    }
    list(x = x, trace = trace)
  }

  # Kruskal's stress-1 in plain R.
  stress_r <- function (d, x, ties = 'secondary') {
    dx <- as.double(dist(x))
    sqrt(sum((dx - iso_r(as.double(d), dx, ties))^2) / sum(dx^2))
  }

  # How far x is from a fixed point of the Guttman transform (with
  # disparities fitted to x): the distance between the two
  # configurations, each centered and scaled to unit size.
  unit <- function (x) { x <- scale(x, scale = FALSE); x / sqrt(sum(x^2)) }
  stationarity <- function (d, x, ties) {
    x <- matrix(unclass(x), nrow(x))
    sqrt(sum((unit(x) - unit(guttman_r(x, iso_r(as.double(d), as.double(dist(x)), ties))))^2))
  }

  # Root mean squared difference after a Procrustes rotation, with
  # both configurations centered and scaled to unit size.
  pro_rmse <- function (x, y) {
    x <- unit(unclass(x)[, , drop = FALSE])
    y <- unit(unclass(y)[, , drop = FALSE])
    sqrt(max(0, 1 - sum(svd(crossprod(x, y))$d)^2) / nrow(x))
  }

  rel_diff <- function (x, y) max(abs(x - y)) / max(abs(y))



  # Stored results from scikit-learn, so at home only. ====
  # Same iterations as scikit-learn's non-metric smacof(), from the same
  # starting configuration. From data-raw/expected/nmds.r.

  if (home) {

    fx <- readRDS('expected/nmds.rds')

    for (i in seq_len(nrow(fx$cases))) {

      cs  <- fx$cases[i, ]
      fv  <- fx$values[[i]]
      d   <- fx$dists[[cs$dist]]
      res <- engine(d, array(fv$init, c(dim(fv$init), 1)))

      expect_identical(res$iters, cs$n_iter, info = cs$case)
      expect_true(max(abs(res$trace[fv$trace$iter] / fv$trace$stress - 1)) < 1e-10, info = cs$case)
      expect_true(rel_diff(res$x_all[, , 1], fv$x) < 1e-8, info = cs$case)

      # Final quality: no worse than scikit-learn's own best, from 8 random
      # starts or from classical MDS, with the same stress measure.
      if (cs$k >= 2) {
        ord <- nmds(d, n_axes = cs$k, ties = 'secondary', cpus = 1)
        expect_true(attr(ord, 'stress') <= min(cs$skl_random, cs$skl_cmds) * (1 + 1e-3), info = cs$case)
      }
    }

    # The plain R version agrees too.
    d    <- fx$dists$bray30
    init <- fx$values[[1]]$init
    ref  <- smacof_r(d, init, max_iters = 60L)
    res  <- engine(d, array(init, c(dim(init), 1)), max_iters = 60L)
    expect_true(max(abs(res$trace / ref$trace - 1)) < 1e-10)

    # Tied dissimilarities and zeros are recognized.
    d   <- fx$dists$ties50
    res <- engine(d, max_iters = 1L, k = 2L, n_starts = 1L)
    expect_identical(res$n_zero,   as.double(sum(d == 0)))
    expect_identical(res$n_groups, as.double(length(unique(d[d > 0]))))

    # The primary approach to ties matches plain R.
    init <- array(runif(50 * 2), c(50, 2, 1))
    res  <- engine(d, init, ties = 'primary', max_iters = 40L)
    ref  <- smacof_r(d, init[, , 1], ties = 'primary', max_iters = 40L)
    expect_identical(res$iters, length(ref$trace))
    expect_true(max(abs(res$trace[seq_len(res$iters)] / ref$trace - 1)) < 1e-9)
    expect_true(rel_diff(res$x_all[, , 1], ref$x) < 1e-8)
  }



  # Stored results from vegan::metaMDS(), which also uses the primary ====
  # approach to ties: the same solution quality, and the same stress-1
  # for the same configuration (when dm has no zeros), so at home only.
  # From data-raw/expected/nmds_vegan.r.

  if (home) {
    for (vg in readRDS('expected/nmds_vegan.rds')) {
      ord <- nmds(vg$dm, n_starts = 10, cpus = 2)
      expect_equal(vg$stress, stress_r(vg$dm, vg$points, 'primary'), tolerance = 1e-6)
      expect_true(attr(ord, 'stress') <= vg$stress * (1 + 1e-3) + 1e-4)
      expect_true(pro_rmse(ord, vg$points) < 0.01)
    }
  }



  # Tied dissimilarities and zeros are recognized. ====

  res <- engine(d_jacc, max_iters = 1L, k = 2L, n_starts = 1L)
  expect_identical(res$n_zero,   as.double(sum(d_jacc == 0)))
  expect_identical(res$n_groups, as.double(length(unique(d_jacc[d_jacc > 0]))))

  # Values within 1e-15 of the largest of a group's first value are tied,
  # as in scikit-learn, even when chains of such values cross the
  # boundaries of the parallel tie scan.
  set.seed(3)
  n     <- 400; N <- n * (n - 1) / 2
  chain <- 0.5 + (0:14999) * 3 * 2^-53
  v     <- c(runif(60000, 0.1, 0.4), chain, runif(N - 75001, 0.6, 0.99), 1)
  v     <- sample(v)
  d     <- structure(v, class = 'dist', Size = n, Diag = FALSE, Upper = FALSE)
  sv    <- sort(v)
  first <- -1; groups <- 0
  for (x in sv) if (x > first && x - first >= 1e-15) { groups <- groups + 1; first <- x }
  res   <- engine(d, k = 1L, n_starts = 1L, max_iters = 0L, cpus = 2L)
  expect_identical(res$n_groups, groups)
  expect_true(groups < N - 10000)



  # The primary approach to ties matches plain R. ====

  init <- array(runif(50 * 2), c(50, 2, 1))
  res  <- engine(d_jacc, init, ties = 'primary', max_iters = 40L)
  ref  <- smacof_r(d_jacc, init[, , 1], ties = 'primary', max_iters = 40L)
  expect_identical(res$iters, length(ref$trace))
  expect_true(max(abs(res$trace[seq_len(res$iters)] / ref$trace - 1)) < 1e-9)
  expect_true(rel_diff(res$x_all[, , 1], ref$x) < 1e-8)

  # Stress-1 as reported matches plain R, for both approaches; primary
  # stress is never above secondary for the same configuration.
  for (ties in c('primary', 'secondary')) {
    ord <- nmds(d_jacc, ties = ties, cpus = 1)
    expect_equal(attr(ord, 'stress'), stress_r(d_jacc, ord, ties), tolerance = 1e-10)
    expect_identical(attr(ord, 'ties'), ties)
  }
  expect_true(stress_r(d_jacc, ord, 'primary') <= stress_r(d_jacc, ord, 'secondary'))



  # The returned configuration is a fixed point of the Guttman ====
  # transform, and its stress is the one plain steps reach from the same
  # start. Tight tolerances take a while, so for k = 1, 2, 3 at home.

  for (k in if (home) 1:3 else 2L) for (ties in c('primary', 'secondary')) {
    ord <- nmds(d_bray, n_axes = k, ties = ties, n_starts = 0, tol = 1e-12, max_iters = 1e5, cpus = 1)
    expect_true(attr(ord, 'converged'))
    expect_true(stationarity(d_bray, ord, ties) < 1e-5, info = paste(k, ties))
    res <- engine(d_bray, pcoa_init(d_bray, k), ties = ties, tol = 1e-12,
                  max_iters = 100000L, dscale = dscale_of(d_bray))
    expect_equal(attr(ord, 'stress'), res$stress, tolerance = 1e-6, info = paste(k, ties))
  }

  # With default settings too.
  for (k in 1:3) {
    ord <- nmds(d_bray, n_axes = k, cpus = 1)
    expect_true(stationarity(d_bray, ord, 'primary') < 1e-3)
  }



  # Return value. ====

  ord <- nmds(d_bray, cpus = 1)
  expect_identical(dim(ord), c(60L, 2L))
  expect_identical(dimnames(ord), list(attr(d_bray, 'Labels'), c('NMDS1', 'NMDS2')))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'stress', 'converged',
      'n_iters', 'n_starts', 'best_start', 'n_best', 'ties', 'class') )
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(attr(ord, 'method'), 'nmds')
  expect_identical(attr(ord, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(ord, 'prop'), c(NA_real_, NA_real_))
  expect_true(attr(ord, 'stress') > 0 && attr(ord, 'stress') < 1)
  expect_true(is.integer(attr(ord, 'n_iters')) && attr(ord, 'n_iters') >= 1L)
  expect_true(isTRUE(attr(ord, 'converged')))
  expect_identical(attr(ord, 'n_starts'), 4L)
  expect_identical(attr(ord, 'ties'), 'primary')
  expect_equal(attr(ord, 'stress'), stress_r(d_bray, ord, 'primary'), tolerance = 1e-10)
  expect_identical(dist_ord(d_bray, 'nmds', cpus = 1), ord)

  # The best start is the one reported, and n_best counts the starts
  # within 0.1% of it.
  init <- pcoa_init(d_bray, 2)
  res  <- engine(d_bray, init, n_starts = 4L, ties = 'primary', max_iters = 1000L,
                 tol = 1e-7, relax = 1L, dscale = dscale_of(d_bray))
  expect_identical(attr(ord, 'stress'), min(res$stress_all))
  expect_identical(attr(ord, 'best_start'), which.min(res$stress_all) - 1L)
  expect_identical(attr(ord, 'n_iters'), res$iters_all[[which.min(res$stress_all)]])
  expect_identical(attr(ord, 'n_best'), sum(res$stress_all <= min(res$stress_all) * 1.001))

  # Centered, on principal axes, scaled to the RMS distance, signs fixed.
  x <- unclass(ord); attributes(x) <- list(dim = dim(x))
  expect_equal(colMeans(x), c(0, 0), tolerance = 1e-12)
  expect_equal(cor(x)[1, 2], 0, tolerance = 1e-10)
  expect_true(var(x[, 1]) >= var(x[, 2]))
  expect_equal(mean(dist(x)^2), mean(d_bray^2), tolerance = 1e-12)
  expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))
  expect_equal(pro_rmse(x, res$x_all[, , which.min(res$stress_all)]), 0, tolerance = 1e-7)

  # k = 3, 4 and 1, and a symmetric matrix input.
  ord3 <- nmds(d_bray, n_axes = 3, n_starts = 1)
  ord4 <- nmds(d_bray, n_axes = 4, n_starts = 1)
  expect_identical(colnames(ord4), paste0('NMDS', 1:4))
  expect_true(attr(ord4, 'stress') < attr(ord3, 'stress'))
  expect_true(attr(ord3, 'stress') < attr(ord, 'stress'))
  expect_identical(dim(nmds(d_pln, n_axes = 1, n_starts = 1)), c(40L, 1L))
  expect_identical(nmds(as.matrix(d_bray), n_starts = 1), nmds(d_bray, n_starts = 1))

  # max_iters = 0: the starting configurations themselves.
  ord0 <- nmds(d_bray, max_iters = 0, n_starts = 0)
  expect_identical(attr(ord0, 'n_iters'), 0L)
  expect_false(attr(ord0, 'converged'))
  expect_equal(attr(ord0, 'stress'), stress_r(d_bray, init[, , 1], 'primary'), tolerance = 1e-10)

  # Stopping at max_iters (its last step plain), or at tol = 0.
  ord <- nmds(d_bray, max_iters = 3, n_starts = 0)
  expect_identical(attr(ord, 'n_iters'), 3L)
  expect_false(attr(ord, 'converged'))
  ord <- nmds(d_pln, tol = 0, max_iters = 5000, n_starts = 0)
  expect_true(attr(ord, 'converged'))



  # Results don't depend on the scale of dm. ====

  for (p in c(-900, -1, 3, 400)) {
    x <- nmds(d_bray * 2^p, n_starts = 1)
    y <- nmds(d_bray, n_starts = 1)
    expect_identical(unclass(x)[, ], unclass(y)[, ] * 2^p)
    expect_identical(attr(x, 'stress'), attr(y, 'stress'))
  }
  for (s in c(1e-300, 1e-160, 1e100)) {
    x <- nmds(d_pln * s, n_starts = 1)
    y <- nmds(d_pln, n_starts = 1)
    expect_equal(unclass(x)[, ] / s, unclass(y)[, ], tolerance = 1e-6)
    expect_equal(attr(x, 'stress'), attr(y, 'stress'), tolerance = 1e-6)
  }



  # Identical across runs and thread counts (and, at home, sessions). ====

  for (args in list(
      list(dm = d_bray),                                  # several starts at once
      list(dm = d_jacc, ties = 'secondary', n_axes = 3) )) {
    x <- do.call(nmds, c(args, cpus = 1))
    expect_identical(do.call(nmds, c(args, cpus = 1)), x)
    for (cpus in c(2, 3, 4, 100))
      expect_identical(do.call(nmds, c(args, cpus = cpus)), x)
  }

  # Several starts at once, or one at a time when memory is short, with
  # threads within a start (passes in waves between interrupt checks,
  # parallel re-sorting of tie groups): the same results.
  init <- pcoa_init(d_big, 2)
  for (ties in c('primary', 'secondary')) {
    r1 <- engine(d_big, init, n_starts = 3L, max_iters = 25L, ties = ties, relax = 1L)
    for (cpus in unique(c(2L, many))) for (mem in c(1L, 1024L)) for (check in c(2^26, 1000))
      expect_identical(
        engine(d_big, init, n_starts = 3L, max_iters = 25L, ties = ties, relax = 1L,
               cpus = cpus, mem = mem, check = check, parmin = 1000),
        r1 )
  }

  # A new R session, on the library of this one, gives the same results
  # for nmds() and mds().
  if (home) {
    rscript <- file.path(R.home('bin'), 'Rscript')
    lib     <- dirname(getNamespaceInfo('ecodive', 'path'))
    out     <- tempfile(fileext = '.rds')
    src     <- tempfile(fileext = '.R')
    src_d   <- tempfile(fileext = '.rds')
    saveRDS(d_bray, src_d)
    writeLines(c(
      sprintf('.libPaths(c(%s, .libPaths()))', deparse(lib)),
      sprintf('d <- readRDS(%s)', deparse(src_d)),
      sprintf('saveRDS(list(getNamespaceInfo("ecodive", "path"), ecodive::nmds(d, cpus = 3), ecodive::mds(d, cpus = 3)), %s)', deparse(out)) ),
      src )
    status <- system2(rscript, c('--vanilla', shQuote(src)), stdout = FALSE, stderr = FALSE)
    expect_identical(status, 0L)
    child <- readRDS(out)
    expect_identical(normalizePath(child[[1]]), normalizePath(getNamespaceInfo('ecodive', 'path')))
    expect_identical(child[[2]], nmds(d_bray, cpus = 1))
    expect_identical(child[[3]], mds(d_bray, cpus = 1))
    unlink(c(out, src, src_d))
  }

  # R's random number generator is neither used nor changed.
  set.seed(5); r1 <- runif(1)
  set.seed(5); x  <- nmds(d_jacc, n_starts = 2); r2 <- runif(1)
  expect_identical(r1, r2)
  set.seed(5);  x1 <- nmds(d_jacc, n_starts = 2)
  set.seed(99); x2 <- nmds(d_jacc, n_starts = 2)
  expect_identical(x1, x2)
  rm('.Random.seed', envir = globalenv())
  x <- nmds(d_jacc, n_starts = 1)
  expect_false(exists('.Random.seed', envir = globalenv()))
  set.seed(1)

  # The seed changes the random starts.
  r1 <- engine(d_bray, k = 2L, n_starts = 2L, seed = 1L)
  r2 <- engine(d_bray, k = 2L, n_starts = 2L, seed = 2L)
  expect_false(identical(r1$stress_all, r2$stress_all))
  expect_identical(engine(d_bray, k = 2L, n_starts = 2L, seed = 1L), r1)



  # Edge cases. ====

  # Perfect fit: points in a plane.
  pts <- cbind(x = sin(1:30), y = cos(1:30 / 2) * 3)
  ord <- nmds(dist(pts), max_iters = 2000, tol = 1e-12)
  expect_true(attr(ord, 'stress') < 1e-3)

  # All distances tied.
  ord <- nmds(as.dist(matrix(1, 10, 10)))
  expect_identical(dim(ord), c(10L, 2L))
  expect_true(all(is.finite(ord)))
  expect_true(attr(ord, 'stress') >= 0 && attr(ord, 'stress') < 1)

  # Exactly equal configuration distances among tied dissimilarities: the
  # corners of a square.
  sq  <- dist(cbind(c(0, 1, 0, 1), c(0, 0, 1, 1)))
  ord <- nmds(sq, ties = 'secondary', max_iters = 1000, tol = 1e-12)
  expect_true(attr(ord, 'stress') < 1e-4)

  # Duplicate samples: zero distances, pulled together, for both
  # approaches to ties.
  m <- rand_counts(20)
  m <- rbind(m, D1 = m[1, ], D2 = m[1, ])
  for (ties in c('primary', 'secondary')) {
    ord <- nmds(bray(m), ties = ties)
    expect_true(all(is.finite(ord)))
    expect_true(max(abs(ord['D1', ] - ord[1, ])) < 0.05 * max(abs(ord)))
  }

  # Zero, one, two, three, and k + 1 samples.
  ord <- nmds(as.dist(matrix(numeric(0), 0, 0)))
  expect_identical(dim(ord), c(0L, 0L))
  expect_identical(attr(ord, 'stress'), 0)

  ord <- nmds(dist(matrix(1:3, 1, dimnames = list('A', NULL))))
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(rownames(ord), 'A')
  expect_identical(attr(ord, 'stress'), 0)
  expect_identical(attr(ord, 'n_best'), 1L)
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'stress', 'converged',
      'n_iters', 'n_starts', 'best_start', 'n_best', 'ties', 'class') )

  ord <- nmds(dist(c(a = 1, b = 3)))
  expect_identical(dimnames(ord), list(c('a', 'b'), 'NMDS1'))
  expect_equal(abs(unname(ord[, 1])), c(1, 1))
  expect_equal(attr(ord, 'stress'), 0)

  ord <- nmds(bray(counts[1:3, ]))
  expect_identical(dim(ord), c(3L, 2L))
  expect_true(attr(ord, 'stress') < 1e-3)

  ord <- nmds(bray(counts), n_axes = 3)
  expect_identical(dim(ord), c(4L, 3L))

  # k larger than n - 1 (even huge), and integer storage.
  d <- as.dist(matrix(c(0L, 3L, 4L, 3L, 0L, 5L, 4L, 5L, 0L), 3))
  expect_identical(nmds(d, n_axes = 10), nmds(d, n_axes = 2))
  expect_identical(nmds(d, n_axes = 1e12), nmds(d, n_axes = 2))

  # Huge max_iters still converges; cpus larger than n.
  expect_identical(nmds(d_pln, max_iters = 1e12), nmds(d_pln, max_iters = 1e6))
  expect_identical(nmds(d, cpus = 100), nmds(d, cpus = 1))
  expect_identical(nmds(d_bray, n_starts = 1, max_iters = 5, cpus = 1e5), nmds(d_bray, n_starts = 1, max_iters = 5, cpus = 1))

  # A PCoA axis without a positive eigenvalue gets random coordinates; so
  # does a start that isn't finite, or has no spread.
  d   <- as.dist(matrix(c(0, 1, 1, 3, 1, 0, 1, 1, 1, 1, 0, 1, 3, 1, 1, 0), 4))
  ord <- nmds(d, n_axes = 3, n_starts = 0)
  expect_true(all(apply(ord, 2, sd) > 0))
  res <- engine(d_bray, array(c(1:60, rep(5, 60)), c(60, 2, 1)), max_iters = 0L)
  expect_true(sd(res$x_all[, 2, 1]) > 0)
  res <- engine(d_bray, array(1, c(60, 2, 1)), max_iters = 0L)
  expect_true(all(apply(res$x_all[, , 1], 2, sd) > 0))
  res <- engine(d_bray, array(c(NaN, 1:119), c(60, 2, 1)), max_iters = 0L)
  expect_true(all(is.finite(res$x_all)))

  # A start whose points all collapse together.
  d    <- as.dist(matrix(c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0), 4))
  init <- array(c(0, 1, 0, 1, 0, 1, 0, 1), c(4, 2, 1))
  res  <- engine(d, init)
  expect_identical(res$stress, 1)
  expect_false(res$converged)
  expect_true(all(is.finite(res$points)))

  # Large tie groups that move a lot, for the primary approach.
  d   <- as.dist(matrix(rep(1:3, length.out = 900), 30, 30))
  ord <- nmds(d, n_starts = 2)
  expect_true(all(is.finite(ord)))



  # Invalid input. ====

  d_na <- d_bray; d_na[2] <- NA
  expect_error(nmds(d_na), 'NA')
  d_na[2] <- NaN
  expect_error(nmds(d_na), 'NA')
  d_na[2] <- Inf
  expect_error(nmds(d_na), 'infinite')
  d_na[2] <- -0.1
  expect_error(nmds(d_na), 'negative')
  expect_error(nmds(as.dist(matrix(0, 5, 5))), 'greater than 0')
  expect_error(nmds(d_bray * 1e151), 'too large')

  expect_error(nmds(d_bray, n_axes = 0),            '`n_axes`')
  expect_error(nmds(d_bray, n_axes = -1),           '`n_axes`')
  expect_error(nmds(d_bray, n_axes = 1.5),          '`n_axes`')
  expect_error(nmds(d_bray, n_axes = NA),           '`n_axes`')
  expect_error(nmds(d_bray, n_axes = Inf),          '`n_axes`')
  expect_error(nmds(d_bray, n_axes = 1:2),          '`n_axes`')
  expect_error(nmds(d_bray, n_axes = '2'),          '`n_axes`')
  expect_error(nmds(d_bray, n_starts = -1),    '`n_starts`')
  expect_error(nmds(d_bray, n_starts = 1001),  '`n_starts`')
  expect_error(nmds(d_bray, n_starts = 1e10),  '`n_starts`')
  expect_error(nmds(d_bray, n_starts = NULL),  '`n_starts`')
  expect_error(nmds(d_bray, max_iters = 1:2),   '`max_iters`')
  expect_error(nmds(d_bray, max_iters = -1),    '`max_iters`')
  expect_error(nmds(d_bray, max_iters = NaN),   '`max_iters`')
  expect_error(nmds(d_bray, tol = -1),         '`tol`')
  expect_error(nmds(d_bray, tol = NA),         '`tol`')
  expect_error(nmds(d_bray, tol = Inf),        '`tol`')
  expect_error(nmds(d_bray, tol = c(1, 2)),    '`tol`')
  expect_error(nmds(d_bray, tol = 'a'),        '`tol`')
  expect_error(nmds(d_bray, ties = 'weak'),    '`ties`')
  expect_error(nmds(d_bray, ties = NA),        '`ties`')
  expect_error(nmds(d_bray, ties = c('primary', 'secondary')), '`ties`')
  expect_error(nmds(d_bray, seed = 1.5),       'seed')
  expect_error(nmds(d_bray, seed = 2^40),      'seed')
  expect_error(nmds(d_bray, seed = NA),        'seed')
  expect_error(nmds(d_bray, cpus = 0),         'cpus')
  expect_error(nmds(d_bray, cpus = NA),        'cpus')
  expect_error(nmds(d_bray, cpus = 1:2),       'cpus')
  expect_error(nmds(d_bray, tries = 2),        'unused argument')
  expect_error(nmds(d_bray, maxit = 2),        'unused argument')

  expect_error(nmds(matrix(1:6, 2)),           'symmetric')
  expect_error(nmds(matrix(c(0, 1, 2, 0), 2)), 'symmetric')
  expect_error(nmds(matrix('a', 2, 2)),        'symmetric')
  expect_error(nmds(list(1, 2)),               'symmetric')
  expect_error(nmds(NULL),                     'symmetric')
  expect_error(nmds(structure(1:2, class = 'dist', Size = 3L)), 'valid dist')

#})
