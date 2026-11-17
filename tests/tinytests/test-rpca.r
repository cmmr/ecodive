#test_that("rpca", {

  # Threads for direct .Call()s and cpus arguments: on CRAN, at most 2
  # (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  # Largest absolute difference relative to the largest value of `y`, after
  # matching each column's sign.
  rel <- function (x, y) {
    x <- unname(as.matrix(x)); y <- unname(as.matrix(y))
    attributes(x) <- list(dim = dim(x))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # Distance between the column spaces of `a` and `b`: the largest
  # difference between their orthogonal projectors.
  subspace <- function (a, b) {
    qa <- qr.Q(qr(a)); qb <- qr.Q(qr(b))
    max(abs(tcrossprod(qa) - tcrossprod(qb)))
  }

  # rpca()'s C function, with its intermediate results.
  run <- function (counts, k = 3L, filters = c(0, 0, 0), step = 1e4, cpus = 1L)
    .Call(C_rpca, counts, 1L, as.integer(k), 5L, 1e-5,
          filters[1], filters[2], filters[3], TRUE, TRUE, step, as.integer(cpus))

  # Tables with planted groups.
  set.seed(1)
  group_counts <- function (sizes, p, depth = 2000, dropout = 0) {
    base <- rnorm(p, 0, 1.5)
    do.call(rbind, lapply(seq_along(sizes), function (g) {
      mu <- base + rnorm(p, 0, 1.5)
      t(replicate(sizes[g], {
        x <- rpois(p, depth * exp(rnorm(1, 0, 0.5)) * exp(mu + rnorm(p, 0, 0.5)) / sum(exp(mu)))
        x * (runif(p) >= dropout)
      }))
    }))
  }
  named <- function (m) {
    dimnames(m) <- list(paste0('S', seq_len(nrow(m))), paste0('F', seq_len(ncol(m))))
    m
  }



  # ====================================================================
  # gemelli 0.0.13 reference results (see data-raw/expected/rpca.r), stage
  # by stage.
  #
  # OptSpace is deterministic once its start is fixed, and its iterations
  # are equivariant under rotations of the initial factors, so rpca()
  # should give gemelli's results up to rounding error, and up to the
  # tolerance of the initial SVD (gemelli: ARPACK with tol = 0, i.e. to
  # machine precision; here: block Lanczos, to a residual of 1e-14 times
  # ||M||^2). The remaining differences come from computing the rclr as
  # log(x) - mean(log(x)) rather than log(x / sum(x)) - mean(...), from the
  # order of summation, and from the different SVD algorithms. On these
  # well-conditioned tables they are about 1e-13 (relative) at most; the
  # tolerance leaves a factor of 1000 for other compilers and platforms.
  # ====================================================================

  if (at_home()) {

    TOL      <- 1e-10
    expected <- readRDS('expected/rpca.rds')

    for (fx in expected) {

      name <- fx$name
      k    <- fx$k[1]
      res  <- run(fx$counts, k, fx$filters)
      n    <- length(res$samples)

      # Filtering keeps the same samples and features.
      expect_identical(res$samples,  as.integer(fx$e2e_samples), info = name)
      expect_identical(res$features, as.integer(fx$e2e_features), info = name)
      expect_identical(res$k, as.integer(k), info = name)

      # gemelli's sample loadings are u and its eigenvalues s; the C code's
      # scores are u * s and its eig s^2.
      u <- sweep(res$scores, 2, sqrt(res$eig), '/')
      expect_true(rel(u, fx$e2e_u) < TOL, info = name)
      expect_true(rel(res$loadings, fx$e2e_v) < TOL, info = name)
      expect_equal(sqrt(res$eig), c(fx$e2e_eig), tolerance = TOL, info = name)
      expect_equal(res$eig / sum(res$eig), c(fx$e2e_prop), tolerance = TOL, info = name)
      expect_equal(res$dist, c(fx$e2e_dist), tolerance = TOL, info = name)

      # Through rpca() itself, whose eig is s^2 / (n - 1).
      ord <- rpca(fx$counts, k = k, min_sample_count = fx$filters[1],
                  min_feature_count = fx$filters[2],
                  min_feature_frequency = fx$filters[3], distances = TRUE,
                  cpus = 1)
      expect_identical(rownames(ord), rownames(fx$counts)[fx$e2e_samples], info = name)
      expect_identical(rownames(attr(ord, 'loadings')), colnames(fx$counts)[fx$e2e_features], info = name)
      expect_true(rel(ord, fx$e2e_u %*% diag(c(fx$e2e_eig), k)) < TOL, info = name)
      expect_true(rel(attr(ord, 'loadings'), fx$e2e_v) < TOL, info = name)
      expect_equal(sqrt(attr(ord, 'eig') * (n - 1)), c(fx$e2e_eig), tolerance = TOL, info = name)
      expect_equal(attr(ord, 'prop'), c(fx$e2e_prop), tolerance = TOL, info = name)
      expect_equal(c(attr(ord, 'dm')), c(fx$e2e_dist), tolerance = TOL, info = name)

      # And with threads.
      expect_identical(rpca(fx$counts, k = k, min_sample_count = fx$filters[1],
                            min_feature_count = fx$filters[2],
                            min_feature_frequency = fx$filters[3],
                            distances = TRUE, cpus = max_cpus), ord, info = name)

      if (is.null(fx$rclr)) next

      # rclr values, in the order of the non-zero counts by sample.
      expect_equal(res$rclr, t(fx$rclr)[!is.na(t(fx$rclr))], tolerance = TOL, info = name)

      # Initial factors: the same subspaces (their bases are arbitrary), and
      # the same scaling.
      expect_true(subspace(res$u0, fx$u0) < TOL, info = name)
      expect_true(subspace(res$v0, fx$v0) < TOL, info = name)
      expect_equal(colSums(res$u0^2), colSums(fx$u0^2), tolerance = TOL, info = name)
      expect_equal(colSums(res$v0^2), colSums(fx$v0^2), tolerance = TOL, info = name)

      # The distortion after initialization and after each iteration.
      expect_equal(res$cost, c(fx$cost), tolerance = TOL, info = name)
      expect_identical(res$n_iters, 5L, info = name)

      # OptSpace's final factors, up to the rotation that the starting bases
      # leave free: the same completed matrix, the same singular values of
      # S, and the same distances between rows of U.
      expect_true(rel(res$u %*% res$s %*% t(res$v), fx$U %*% fx$S %*% t(fx$V)) < TOL, info = name)
      expect_equal(svd(res$s)$d, svd(fx$S)$d, tolerance = TOL, info = name)
      expect_true(subspace(res$u, fx$U) < TOL, info = name)
      expect_true(subspace(res$v, fx$V) < TOL, info = name)

      # The final PCA.
      expect_true(rel(u, fx$u) < TOL, info = name)
      expect_true(rel(res$loadings, fx$v) < TOL, info = name)
      expect_equal(sqrt(res$eig), c(fx$eig), tolerance = TOL, info = name)
      expect_equal(res$eig / sum(res$eig), c(fx$prop), tolerance = TOL, info = name)
    }


    # OptSpace's Grassmann penalty terms, which gemelli's default step size
    # leaves at zero for tables of this size, against gemelli's OptSpace
    # with step_size = 2.
    fx  <- expected$groups
    res <- run(fx$counts, 3, step = fx$pen_step[1])
    expect_equal(res$cost, c(fx$pen_cost), tolerance = TOL)
    expect_true(res$cost[6] != fx$cost[6])
    expect_true(rel(res$u %*% res$s %*% t(res$v), fx$pen_U %*% fx$pen_S %*% t(fx$pen_V)) < TOL)

    # Smaller steps make OptSpace diverge, in gemelli too (NaN).
    expect_error(run(fx$counts, 3, step = 0.5), 'diverged')


    # Planted groups: recovered, with the same separation as gemelli.
    ord <- rpca(fx$counts, cpus = 1)
    grp <- rep(1:3, c(7, 7, 6))

    # Between-group share of the variance on the first two axes.
    separation <- function (x) {
      x  <- as.matrix(x)[, 1:2]
      mu <- apply(x, 2, function (col) ave(col, grp))
      sum(sweep(mu, 2, colMeans(x))^2) / sum(sweep(x, 2, colMeans(x))^2)
    }
    gem <- fx$e2e_u %*% diag(c(fx$e2e_eig))
    expect_equal(separation(ord), separation(gem), tolerance = TOL)
    expect_true(separation(ord) > 0.9)

    # Every sample is nearest its own group's centroid.
    cen  <- apply(ord[, 1:2], 2, function (col) tapply(col, grp, mean))
    near <- apply(ord[, 1:2], 1, function (x) which.min(colSums((t(cen) - x)^2)))
    expect_identical(unname(near), grp)
  }



  # ====================================================================
  # Compositional scale invariance.
  #
  # Multiplying a sample's counts by c adds log(c) to each of its log
  # counts, which the centering removes exactly, but only up to rounding
  # error of order 1e-16 * |log(c)|. So the results agree to about
  # 1e-13 even for c = 1e-300 or 1e300.
  # ====================================================================

  m     <- named(group_counts(c(10, 10), 30, dropout = 0.3))
  ord   <- rpca(m, distances = TRUE, cpus = 1)
  scale <- 10^seq(-300, 300, length.out = nrow(m))

  for (scaled in list(m * scale, m / rowSums(m), m * 1e-300, m * 1e300)) {
    o2 <- rpca(scaled, distances = TRUE, cpus = 1)
    expect_true(rel(o2, ord) < 1e-9)
    expect_true(rel(attr(o2, 'loadings'), attr(ord, 'loadings')) < 1e-9)
    expect_equal(attr(o2, 'eig'), attr(ord, 'eig'), tolerance = 1e-9)
    expect_equal(c(attr(o2, 'dm')), c(attr(ord, 'dm')), tolerance = 1e-9)
  }

  # Subnormal counts too.
  expect_true(rel(rpca(m * 4.9e-324, cpus = 1), ord) < 1e-9)



  # ====================================================================
  # Relation to pca() and robust_aitchison(): with no zeros, and k at the
  # rank of the clr matrix, the completed matrix is the clr matrix itself,
  # so the scores and variances are pca()'s with norm = 'clr', and the
  # distances between scores are the (robust) Aitchison distances.
  # gemelli's distances (from OptSpace's U) are not.
  # ====================================================================

  x <- named(matrix(rpois(60, 20) + 1, 12))
  o <- rpca(x, k = 4, distances = TRUE, cpus = 1)
  p <- pca(x, norm = 'clr', k = 4, cpus = 1)
  expect_equal(c(dist(o)), c(robust_aitchison(x)), tolerance = 1e-12)
  expect_equal(c(dist(o)), c(aitchison(x)), tolerance = 1e-12)
  expect_true(max(abs(attr(o, 'dm') - robust_aitchison(x))) > 0.1)
  expect_true(rel(o, p) < 1e-10)
  expect_true(rel(attr(o, 'loadings'), attr(p, 'loadings')) < 1e-10)
  expect_equal(attr(o, 'eig'),   attr(p, 'eig'),   tolerance = 1e-10)
  expect_equal(attr(o, 'prop'),  attr(p, 'prop'),  tolerance = 1e-10)
  expect_equal(attr(o, 'trace'), attr(p, 'trace'), tolerance = 1e-10)



  # ====================================================================
  # Result format.
  # ====================================================================

  ord <- rpca(m, k = 2, distances = TRUE, cpus = 1)

  expect_true(inherits(ord, 'ecodive_ord') && is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings',
      'cost', 'n_iters', 'dm', 'class') )
  expect_identical(dimnames(ord), list(rownames(m), c('RPCA1', 'RPCA2')))
  expect_identical(attr(ord, 'method'), 'rpca')
  expect_identical(dimnames(attr(ord, 'loadings')), list(colnames(m), c('RPCA1', 'RPCA2')))
  expect_length(attr(ord, 'cost'), 6L)
  expect_identical(attr(ord, 'n_iters'), 5L)
  expect_true(all(diff(attr(ord, 'cost')) < 0))
  expect_true(all(attr(ord, 'prop') > 0) && sum(attr(ord, 'prop')) <= 1 + 1e-12)
  expect_equal(attr(ord, 'prop'), attr(ord, 'eig') / attr(ord, 'trace'), tolerance = 1e-14)

  # Scores times loadings is the double-centered completed matrix; eig is
  # each axis's variance.
  X <- unclass(ord) %*% t(attr(ord, 'loadings'))
  expect_true(max(abs(rowSums(X)), abs(colSums(X))) < 1e-10 * max(abs(X)))
  expect_equal(unname(apply(ord, 2, var)), attr(ord, 'eig'), tolerance = 1e-12)
  expect_equal(crossprod(attr(ord, 'loadings')), diag(2), tolerance = 1e-12, check.attributes = FALSE)

  # Each axis's largest score in absolute value is positive.
  expect_identical(ord_flip(unclass(ord)), unclass(ord))

  # The distances form a dist object that other ordinations take.
  d <- attr(ord, 'dm')
  expect_true(inherits(d, 'dist'))
  expect_identical(attr(d, 'Labels'), rownames(m))
  expect_identical(attr(d, 'Size'), nrow(m))
  expect_identical(dim(pcoa(d, k = 2, cpus = 1)), c(nrow(m), 2L))
  expect_null(attr(rpca(m, cpus = 1), 'dm'))

  # Unnamed tables get numbered samples and features.
  ord <- rpca(unname(m), cpus = 1)
  expect_identical(rownames(ord), as.character(seq_len(nrow(m))))
  expect_identical(rownames(attr(ord, 'loadings')), as.character(seq_len(ncol(m))))

  # max_iters and tol.
  o1 <- rpca(m, max_iters = 1, cpus = 1)
  expect_length(attr(o1, 'cost'), 2L)
  expect_identical(attr(o1, 'n_iters'), 1L)
  expect_length(attr(rpca(m, max_iters = 12, cpus = 1), 'cost'), 13L)
  o1 <- rpca(m, tol = 10, cpus = 1)
  expect_length(attr(o1, 'cost'), 2L)
  expect_identical(attr(o1, 'n_iters'), 1L)
  expect_identical(rpca(m, tol = 0, cpus = 1), rpca(m, cpus = 1))

  # Not one of ordination()'s methods: it takes counts.
  expect_error(ordination(dist(m), 'rpca'), 'call rpca\\(\\) directly')



  # ====================================================================
  # Input types and margins give identical results.
  # ====================================================================

  ord <- rpca(m, cpus = 1)
  mi  <- m; storage.mode(mi) <- 'integer'
  expect_identical(rpca(mi, cpus = 1), ord)
  expect_identical(rpca(t(m), margin = 2, cpus = 1), ord)
  expect_identical(rpca(as.data.frame(m), cpus = 1), ord)

  if (requireNamespace('Matrix', quietly = TRUE)) {

    sp <- methods::as(m, 'CsparseMatrix')
    expect_identical(rpca(sp, cpus = 1), ord)
    expect_identical(rpca(Matrix::t(sp), margin = 2, cpus = 1), ord)
    expect_identical(rpca(methods::as(sp, 'TsparseMatrix'), cpus = 1), ord)
    expect_identical(rpca(Matrix::Matrix(m, sparse = FALSE), cpus = 1), ord)
    expect_identical(rpca(methods::as(sp, 'RsparseMatrix'), cpus = 1), ord)

    # A dgTMatrix's repeated entries are added together.
    tm <- methods::as(sp, 'TsparseMatrix')
    i  <- which(tm@x > 1)[1:5]
    half <- tm@x[i] / 2
    tm2 <- tm
    tm2@x[i] <- half
    tm2@i <- c(tm2@i, tm@i[i]); tm2@j <- c(tm2@j, tm@j[i]); tm2@x <- c(tm2@x, half)
    expect_identical(rpca(tm2, cpus = 1), ord)

    # Explicit zeros are zeros.
    sp0 <- sp
    sp0@x[1:3] <- 0
    m0 <- m; m0[cbind(sp@i[1:3] + 1, rep(1, 3))] <- 0
    expect_identical(rpca(sp0, cpus = 1), rpca(m0, cpus = 1))

    # Repeated entries that overflow.
    # (Built slot by slot: sparseMatrix(repr = 'T') needs Matrix >= 1.3,
    # which R 3.6 doesn't have.)
    triplet <- function (i, j, x, dims) {
      tm <- methods::as(matrix(0, dims[1], dims[2]), 'CsparseMatrix')
      tm <- methods::as(tm, 'TsparseMatrix')
      tm@i <- as.integer(i - 1L); tm@j <- as.integer(j - 1L); tm@x <- as.double(x)
      tm
    }
    tm3 <- triplet(c(1, 1, 1:4), c(1, 1, 2:5), c(1e308, 1e308, 1:4), c(4, 5))
    expect_error(rpca(tm3), 'NA, NaN, or infinite')

    # A dgTMatrix in reverse order.
    rev_t <- triplet(rev(rep(1:300, each = 20)), rev(rep(1:20, 300)),
      rev(as.double(rep(1:7, length.out = 6000))), c(300, 20))
    expect_identical(rpca(rev_t, k = 2, cpus = 1), rpca(methods::as(rev_t, 'CsparseMatrix'), k = 2, cpus = 1))
  }

  # Objects that hold a table (here an rbiom-like object holding a data
  # frame, samples in columns).
  rb <- structure(list(counts = as.data.frame(t(m))), class = 'rbiom')
  expect_identical(rpca(rb, cpus = 1), ord)

  # A slam simple_triplet_matrix (built by hand, as slam may not be
  # installed): 1-based, unsorted, with double, integer or logical values.
  stm <- function (x, v = x[x != 0], perm = TRUE) {
    w <- which(x != 0, arr.ind = TRUE)
    o <- if (perm) rev(seq_len(nrow(w))) else seq_len(nrow(w))
    structure(
      list(i = w[o, 1], j = w[o, 2], v = v[o], nrow = nrow(x), ncol = ncol(x),
           dimnames = dimnames(x)),
      class = 'simple_triplet_matrix' )
  }
  expect_identical(rpca(stm(m), cpus = 1), ord)
  expect_identical(rpca(stm(t(m)), margin = 2, cpus = 1), ord)
  expect_identical(rpca(stm(mi, v = mi[mi != 0]), cpus = 1), ord)
  expect_error(rpca(stm(m > 0, v = rep(TRUE, sum(m > 0)))), 'no variation')
  bad_stm <- stm(m); bad_stm$i[1] <- 999L
  expect_error(rpca(bad_stm), 'out of range')
  bad_stm <- stm(m); bad_stm$i <- as.double(bad_stm$i)
  expect_error(rpca(bad_stm), 'invalid')
  bad_stm <- stm(m); bad_stm$v[1] <- -1
  expect_error(rpca(bad_stm), 'negative')
  bad_stm <- stm(m); bad_stm$v[1] <- Inf
  expect_error(rpca(bad_stm), 'NA, NaN, or infinite')



  # ====================================================================
  # Filters.
  # ====================================================================

  f <- m
  f[1:3, ] <- f[1:3, ] * (runif(length(f[1:3, ])) < 0.4)  # sparse samples
  f[, 1]   <- c(5, rep(0, nrow(f) - 1))                   # rare feature
  f[, 2]   <- 0                                           # absent feature

  ord <- rpca(f, cpus = 1)
  expect_identical(rownames(attr(ord, 'loadings')), colnames(f)[-2])

  # Samples by total count, features by total and by frequency.
  tot <- rowSums(f[, -2])
  ord <- rpca(f, min_sample_count = sort(tot)[3], cpus = 1)
  expect_identical(rownames(ord), rownames(f)[tot > sort(tot)[3]])
  ord <- rpca(f, min_feature_count = 5, cpus = 1)
  expect_identical(rownames(attr(ord, 'loadings')), colnames(f)[colSums(f) > 5])
  freq <- 100 * colSums(f > 0) / nrow(f)
  ord  <- rpca(f, min_feature_frequency = 50, cpus = 1)
  expect_identical(rownames(attr(ord, 'loadings')), colnames(f)[freq > 50])

  # Filtering everything away.
  expect_error(rpca(f, min_sample_count = 1e9), 'After filtering, 0 samples and 0 features')
  expect_error(rpca(f, min_feature_frequency = 100), 'After filtering')



  # ====================================================================
  # Edge cases.
  # ====================================================================

  # Samples with a single non-zero feature are dropped, with a warning.
  s1 <- m
  s1[2, ] <- 0; s1[2, 5] <- 7
  s1[4, ] <- 0; s1[4, 1] <- 1
  expect_warning(ord <- rpca(s1, cpus = 1), 'S2, S4')
  expect_identical(rownames(ord), rownames(m)[-c(2, 4)])
  expect_identical(suppressWarnings(rpca(s1, cpus = 1)), rpca(s1[-c(2, 4), ], cpus = 1))
  s1[1:10, ] <- 0
  for (i in 1:10) s1[i, i] <- i
  expect_warning(rpca(s1, cpus = 1), '\\.\\.\\. \\(5 more\\)')

  # All-zero samples and features are dropped silently.
  z <- rbind(m, Z = 0)
  expect_identical(rpca(z, cpus = 1), rpca(m, cpus = 1))

  # Samples whose non-zero counts are all equal are kept: their rclr is 0.
  eq <- m
  eq[1, ] <- (m[1, ] > 0) * 3
  ord <- rpca(eq, cpus = 1)
  expect_true(all(is.finite(ord)) && nrow(ord) == nrow(m))

  # Unless all of them are.
  expect_error(rpca((m > 0) + 0), 'no variation')

  # Or all samples have the same composition, so that nothing varies after
  # the robust clr (identical rows, or proportional ones).
  expect_error(rpca(matrix(rep(1:5, each = 6), 6)), 'no variation')
  expect_error(suppressWarnings(rpca(matrix(c(1, 2, 2, 4), 2))), 'no variation')
  same <- named(matrix(rep(c(3, 0, 5, 1, 9), each = 8), 8) * (1:8))
  expect_error(rpca(same), 'no variation')

  # Nearly so: the axis with the real (small) variation is kept, and the
  # others, at rounding-error level, are dropped with a warning.
  near <- matrix(rep(1:5, each = 6), 6)
  near[1, 1] <- 1.01
  expect_warning(o <- rpca(near, cpus = 1), 'only 1 of the first 3 eigenvalues')
  expect_identical(ncol(o), 1L)
  expect_true(attr(o, 'eig') > 1e-7)
  expect_identical(unname(which.max(abs(o[, 1]))), 1L)
  near[1, 1] <- 1 + 1e-3
  expect_warning(o <- rpca(near, cpus = 1), 'only 1 of the first 3 eigenvalues')
  expect_identical(unname(which.max(abs(o[, 1]))), 1L)

  # No zeros at all.
  expect_true(all(is.finite(rpca(m + 1, cpus = 1))))

  # k too large for the table is reduced, with a warning.
  small <- m[1:4, 1:6]
  expect_warning(ord <- rpca(small, k = 5, cpus = 1), '`k` reduced to 3')
  expect_identical(ncol(ord), 3L)
  expect_warning(ord <- rpca(m[, 1:3], k = 5, cpus = 1), '`k` reduced to 2')
  expect_identical(ncol(ord), 2L)

  # Smallest tables: 2 samples or 2 features, giving 1 axis.
  for (tbl in list(m[1:2, ], m[, 1:2], m[1:2, 1:2] + 1, m[1:3, ], m[, 1:4])) {
    ord <- suppressWarnings(rpca(tbl, cpus = 1))
    expect_true(all(is.finite(ord)) && ncol(ord) >= 1)
  }
  expect_identical(ncol(suppressWarnings(rpca(m[1:2, ], cpus = 1))), 1L)
  expect_error(rpca(m[1, , drop = FALSE]), 'After filtering, 1 sample and')
  expect_error(rpca(m[, 1, drop = FALSE]), 'After filtering')
  expect_error(rpca(m[0, ]), 'valid numeric matrix')
  expect_error(rpca(m[, 0]), 'valid numeric matrix')
  expect_error(rpca(m * 0), 'After filtering, 0 samples and 0 features')

  # When single-feature samples are why too few are left, the error says so.
  expect_error(rpca(diag(5) + 0), '5 samples with a single non-zero feature were removed')
  expect_error(rpca(rbind(diag(3), c(1, 1, 0)) + 0), '3 samples with a single')
  expect_error(rpca(rbind(c(1, 0, 0), c(1, 1, 0))), '1 sample with a single non-zero feature was removed')

  # k = 1.
  expect_identical(dim(rpca(m, k = 1, cpus = 1)), c(nrow(m), 1L))

  # Rank deficient: only two distinct samples, repeated. OptSpace's start
  # and fit then have directions without data; results are still finite,
  # and axes without variance are dropped with a warning.
  rd <- named(m[rep(1:2, 10), ])
  expect_warning(ord <- rpca(rd, k = 3, cpus = 1), 'eigenvalues are > 0')
  expect_true(all(is.finite(ord)))
  expect_true(ncol(ord) < 3)
  expect_true(max(abs(ord[seq(1, 20, 2), ] - ord[rep(1, 10), ])) < 1e-8 * max(abs(ord)))

  # Even rank 1 with k = 1 still works.
  r1 <- named(rbind(c(1, 2, 0, 4), c(2, 4, 0, 8), c(0, 3, 9, 1), c(0, 3, 9, 1)))
  expect_true(all(is.finite(suppressWarnings(rpca(r1, k = 1, cpus = 1)))))

  # Values that aren't counts.
  bad <- m; bad[1, 1] <- -1
  expect_error(rpca(bad), 'negative')
  for (v in list(NA, NaN, Inf, -Inf)) {
    bad <- m; bad[1, 1] <- v
    expect_error(rpca(bad), 'NA, NaN, or infinite', info = format(v))
  }
  bad <- mi; bad[1, 1] <- NA
  expect_error(rpca(bad), 'NA, NaN, or infinite')
  expect_error(rpca(matrix('a', 2, 2)), 'must be numeric')
  expect_error(rpca(data.frame(a = c('x', 'y'), b = 1:2)), 'must be numeric')



  # ====================================================================
  # Argument validation, with no warnings from the checks themselves.
  # ====================================================================

  for (k in list(0, 11, 1.5, NA, NaN, Inf, '3', c(2, 3), 1e300, -1e300, NULL, TRUE))
    expect_error(rpca(m, k = k), '`k` must be a whole number from 1 to 10')
  for (x in list(0, 1.5, NA, 2e6, 1e300, 'a'))
    expect_error(rpca(m, max_iters = x), '`max_iters` must be')
  for (x in list(-1, NA, Inf, 'a', c(1, 2), NULL)) {
    expect_error(rpca(m, tol = x), '`tol` must be')
    expect_error(rpca(m, min_sample_count = x), '`min_sample_count` must be')
    expect_error(rpca(m, min_feature_count = x), '`min_feature_count` must be')
    expect_error(rpca(m, min_feature_frequency = x), '`min_feature_frequency` must be')
  }
  expect_error(rpca(m, min_feature_frequency = 100.5), 'from 0 to 100')
  for (x in list(NA, 1, 'TRUE', c(TRUE, TRUE)))
    expect_error(rpca(m, distances = x), '`distances` must be TRUE or FALSE')
  expect_error(rpca(m, margin = 3), '`margin` must be 1 or 2')
  expect_error(rpca(m, cpus = 0), '`cpus` must be')
  expect_silent(rpca(m, k = 1e1, max_iters = 5, tol = 1e300,
                     min_sample_count = 1e300 * 0, cpus = 1))
  expect_identical(rpca(m, k = 3.0, cpus = 1), rpca(m, k = 3L, cpus = 1))



  # ====================================================================
  # Determinism: identical results for any number of threads, run to
  # run, and regardless of R's RNG.
  # ====================================================================

  # Big enough for threads, with the initial SVD's restarts in use.
  set.seed(2)
  big <- matrix(0, 400, 180)
  nz  <- sample(length(big), 45000)
  big[nz] <- rpois(length(nz), 3 + 4 * ((row(big)[nz] %% 2) == (col(big)[nz] %% 2))) + 1
  big <- named(big)

  ref <- rpca(big, distances = TRUE, cpus = 1)
  expect_identical(nrow(ref), 400L)
  for (cpus in 2:max_cpus)
    expect_identical(rpca(big, distances = TRUE, cpus = cpus), ref)

  # The initial SVD needed restarts here, and threads were used (up to the
  # number of blocks of work).
  res <- run(big, cpus = max_cpus)
  expect_true(res$lanczos_steps > 10)
  expect_true(if (pthreads()) res$threads > 1 else res$threads == 1)
  expect_identical(run(big, cpus = 1)$threads, 1L)
  expect_identical(run(m, cpus = max_cpus)$threads, 1L)

  # R's RNG is neither used nor changed.
  set.seed(99)
  seed <- .Random.seed
  o1 <- rpca(m, distances = TRUE, cpus = 2)
  expect_identical(.Random.seed, seed)
  set.seed(12345)
  expect_identical(rpca(m, distances = TRUE, cpus = 2), o1)

#})
