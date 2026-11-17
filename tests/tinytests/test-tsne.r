#test_that("tsne", {


  # Helpers. ====

  rows <- function (v, n) matrix(v, nrow = n, byrow = TRUE)
  rel  <- function (x, y) max(abs(x - y)) / max(abs(y))

  # Joint probabilities P: dense (K = 0) or K nearest neighbors (CSR).
  tsne_p <- function (d, n, perp, K = 0L, f32 = FALSE, scale = 1, cpus = 1L)
    .Call(C_tsne_p, as.double(d), as.integer(n), perp, as.integer(K), f32, scale, cpus)

  # C_tsne_run() from layout Y: `m` iterations (`exag_iter` of them with
  # exaggeration 12), uncentered. With m = 0, the KL divergence and
  # gradient at Y (exact for dense P, Barnes-Hut with `theta` for sparse).
  # `start`: 0 = Y as given, 1 = Y scaled as the PCoA start, 2 = random.
  tsne_run <- function (P, Y, m = 0L, theta = 0.5, exag_iter = 250L,
                        no_progress = 300L, lr = NULL, start = 0L, seed = 0L,
                        cpus = 1L, dims = dim(Y)) {
    if (is.null(lr)) lr <- max(dims[[1]] / 12 / 4, 50)
    iopt <- as.integer(c(dims, m, exag_iter, cpus, 50L, no_progress, 0L, start, seed))
    .Call(C_tsne_run, P, Y, iopt, c(theta, 12, lr, 1e-7, 0.5, 0.8))
  }

  # Exact gradient and KL divergence in plain R, from a full P matrix.
  ref_grad <- function (P, Y) {
    W <- 1 / (1 + as.matrix(dist(Y))^2)
    diag(W) <- 0
    Q <- W / sum(W)
    M <- (P - Q) * W
    G <- 4 * (rowSums(M) * Y - M %*% Y)
    pe <- P[row(P) != col(P)]
    q  <- pmax(Q[row(Q) != col(Q)], .Machine$double.eps)
    list(kl = sum(pe * log(pmax(pe, .Machine$double.eps) / q)), grad = unname(G))
  }

  # Full matrix from a packed lower triangle, or from CSR.
  full_dense <- function (p, n) { m <- matrix(0, n, n); m[lower.tri(m)] <- p; m + t(m) }
  full_csr   <- function (s, n) {
    m <- matrix(0, n, n)
    m[cbind(rep(seq_len(n), diff(s$ptr)), s$col + 1L)] <- s$val
    m
  }

  # Trustworthiness with m neighbors, as sklearn.manifold.trustworthiness().
  trust <- function (d, y, m = 5) {
    D <- as.matrix(d); diag(D) <- Inf
    E <- as.matrix(dist(y)); diag(E) <- Inf
    n <- nrow(D)
    rk <- t(apply(D, 1, function (v) rank(v, ties.method = 'first')))
    t  <- 0
    for (i in seq_len(n)) {
      r <- rk[i, order(E[i, ])[seq_len(m)]] - m
      t <- t + sum(r[r > 0])
    }
    1 - t * 2 / (n * m * (2 * n - 3 * m - 1))
  }

  # Leave-one-out 5-nearest-neighbor label accuracy on the plot.
  knn_acc <- function (y, labels, m = 5) {
    E <- as.matrix(dist(y)); diag(E) <- Inf
    mean(vapply(seq_along(labels), function (i) {
      tb <- table(labels[order(E[i, ])[seq_len(m)]])
      as.numeric(names(tb)[which.max(tb)]) == labels[[i]]
    }, logical(1)))
  }

  dist_of <- function (v) {
    n <- (1 + sqrt(1 + 8 * length(v))) / 2
    structure(v, Size = as.integer(n), Diag = FALSE, Upper = FALSE, class = 'dist')
  }

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  home <- at_home()
  many <- if (home) 4L else 2L


  set.seed(1)
  m60 <- matrix(rpois(60 * 25, 4), 60, 25, dimnames = list(paste0('S', 1:60), NULL))
  m60[1:20, 1:8] <- m60[1:20, 1:8] + 10
  d60 <- bray(m60)



  # Against scikit-learn 1.9.1 and Rtsne: stored results, so at home only. ====

  if (home) {

    # Frozen reference values, from data-raw/expected/tsne.r: a list of
    # numeric vectors per dataset.
    ref      <- readRDS('expected/tsne.rds')
    fx_blobs <- ref$blobs
    fx_bray  <- ref$bray


    # Affinities match scikit-learn's _joint_probabilities(_nn).

    # scikit-learn rounds squared distances to single precision (relative
    # error up to 2^-24, about 6e-8) before the binary search, which changes
    # p_ij by about that times beta * d^2, a few units for the neighbors that
    # matter. With the same rounding (f32 = TRUE), the results agree to the
    # fixtures' 10 digits; in double precision, to about 1e-7.

    for (fx in list(fx_blobs, fx_bray)) {

      d <- fx$dist
      n <- (1 + sqrt(1 + 8 * length(d))) / 2
      K <- min(n - 1, floor(3 * fx$perp_nn + 1))

      for (f32 in c(TRUE, FALSE)) {
        tol <- if (f32) 1e-9 else 1e-6
        P   <- tsne_p(d, n, fx$perp, 0L, f32)
        S   <- tsne_p(d, n, fx$perp_nn, K, f32)
        expect_true(rel(P, fx$P_dense) < tol)
        expect_identical(S$ptr, as.integer(fx$P_ptr))
        expect_identical(S$col, as.integer(fx$P_col))
        expect_true(rel(S$val, fx$P_val) < tol)
      }
    }


    # Objective and gradient match scikit-learn's _kl_divergence(_bh).

    # With the same affinities, the exact KL divergence and gradient agree to
    # the fixtures' 12 digits. scikit-learn's Barnes-Hut runs in single
    # precision, so agrees to about 1e-7 (relative; KL to ~1e-6 absolute).

    for (fx in list(fx_blobs, fx_bray)) {

      d <- fx$dist
      n <- (1 + sqrt(1 + 8 * length(d))) / 2
      K <- min(n - 1, floor(3 * fx$perp_nn + 1))
      P <- tsne_p(d, n, fx$perp, 0L, TRUE)
      S <- tsne_p(d, n, fx$perp_nn, K, TRUE)

      for (k in 1:3) {

        Y <- rows(fx[[paste0('Y', k)]], n)

        g <- tsne_run(P, Y)
        expect_equal(g$kl, fx[[paste0('exact', k, '_kl')]], tolerance = 1e-10)
        expect_true(rel(g$grad, rows(fx[[paste0('exact', k, '_grad')]], n)) < 1e-10)

        b <- tsne_run(S, Y, theta = 0.5)
        expect_true(abs(b$kl - fx[[paste0('bh', k, '_kl')]]) < 1e-5)
        expect_true(rel(b$grad, rows(fx[[paste0('bh', k, '_grad')]], n)) < 1e-6)
      }

      Y <- rows(fx$Y2, n)
      for (theta in c(0, 0.2, 0.8)) {
        b   <- tsne_run(S, Y, theta = theta)
        ref <- rows(fx[[sprintf('bh2_t%g_grad', theta)]], n)
        expect_true(rel(b$grad, ref) < 1e-6)
      }
    }


    # The optimizer follows scikit-learn's trajectory.

    # From the same start and affinities, positions agree to ~1e-12 after 10
    # iterations and ~1e-8 after 25. t-SNE's early iterations amplify
    # differences by about 2x per iteration (scikit-learn diverges from
    # itself at the same rate when its start is perturbed by 1e-15), so by
    # iteration ~45 rounding differences have grown to the size of the layout.

    for (fx in list(fx_blobs, fx_bray)) {
      d  <- fx$dist
      n  <- (1 + sqrt(1 + 8 * length(d))) / 2
      P  <- tsne_p(d, n, fx$perp, 0L, TRUE)
      Y0 <- rows(fx$Y0, n)
      expect_true(rel(tsne_run(P, Y0, 10L, no_progress = 250L)$points, rows(fx$traj10, n)) < 1e-10)
      expect_true(rel(tsne_run(P, Y0, 25L, no_progress = 250L)$points, rows(fx$traj25, n)) < 1e-6)

      # The second phase (no exaggeration, momentum 0.8) from a spread-out
      # layout, which isn't chaotic: 30 iterations agree to ~1e-12.
      r <- tsne_run(P, rows(fx$p2_start, n), 30L, exag_iter = 0L)
      expect_true(rel(r$points, rows(fx$p2_30, n)) < 1e-10)
      expect_identical(r$n_iter, 30L)
      expect_false(r$diverged)
    }


    # Final quality is as good as scikit-learn's and Rtsne's.

    # Same affinities (perplexity `qperp`), exact gradients. The fixtures
    # hold scikit-learn's KL divergence, trustworthiness (5 neighbors), and
    # 5-NN label accuracy for seeds 0-9 with a random start, and for the
    # PCoA start; and the same for Rtsne(theta = 0), seeds 0-9. Each run
    # ends in one of many local minima, and which one is chaotic (the last
    # bits of the arithmetic decide it), so medians over 10 random starts
    # are compared; and a single PCoA start is checked to be within the
    # range of scikit-learn's random starts, with a margin.

    check_quality <- function (fx, d) {

      q <- vapply(0:9, function (s) {
        o <- tsne(d, perplexity = fx$qperp, init = 'random', seed = s, cpus = 1)
        c(attr(o, 'kl_divergence'), trust(d, o), knn_acc(o, fx$labels))
      }, numeric(3))
      med <- apply(q, 1, stats::median)

      expect_true(med[[1]] <= 1.05 * min(stats::median(fx$q_random_kl), stats::median(fx$rtsne$kl)))
      expect_true(med[[2]] >= max(stats::median(fx$q_random_trust), stats::median(fx$rtsne$trust)) - 0.005)
      expect_true(min(q[3, ]) >= min(fx$q_random_knn, fx$rtsne$knn))

      o <- tsne(d, perplexity = fx$qperp, cpus = 1)
      expect_true(attr(o, 'kl_divergence') <= 1.1 * max(fx$q_random_kl, fx$q_pcoa_kl))
      expect_true(trust(d, o) >= min(fx$q_random_trust, fx$q_pcoa_trust) - 0.005)
      expect_true(knn_acc(o, fx$labels) >= min(fx$q_random_knn, fx$q_pcoa_knn))
    }

    check_quality(fx_blobs, dist_of(fx_blobs$dist))
    check_quality(fx_bray,  bray(rows(fx_bray$counts, fx_bray$counts_dim[[1]])))

    # Well-separated clusters stay separate.
    ord <- tsne(dist_of(fx_blobs$dist), perplexity = 10, cpus = 1)
    expect_identical(knn_acc(ord, fx_blobs$labels), 1)
    expect_identical(as.vector(table(fx_blobs$labels, kmeans(ord, unclass(ord)[c(1, 16, 31), ])$cluster) > 0),
                     as.vector(diag(3) > 0))

    # The scale of the distances doesn't matter (see below); far-off scales
    # still separate the clusters.
    d <- dist_of(fx_blobs$dist)
    for (s in c(1e-30, 1e100)) {
      expect_identical(knn_acc(tsne(d * s, perplexity = 10, cpus = 1), fx_blobs$labels), 1)
      expect_identical(knn_acc(tsne(d * s, perplexity = 10, init = 'random', cpus = 1), fx_blobs$labels), 1)
    }
  }



  # Affinities and gradients: identities within the package. ====

  set.seed(3)
  n  <- 150
  X  <- matrix(rnorm(n * 4), n) + matrix(rep(c(0, 6), each = n / 2 * 4), n, byrow = FALSE)
  dv <- as.vector(dist(X))
  P  <- tsne_p(dv, n, 10)
  S  <- tsne_p(dv, n, 10, 31L)

  # Joint probabilities sum to 1; every neighbor (K = n - 1) gives the
  # dense P (but for the dense P's floor at machine epsilon).
  expect_equal(2 * sum(P), 1, tolerance = 1e-12)
  expect_equal(sum(S$val), 1, tolerance = 1e-12)
  expect_equal(full_csr(tsne_p(dv, n, 10, n - 1L), n), full_dense(P, n), tolerance = 1e-12)

  # Each row's p_{j|i} has perplexity `perplexity`: with K = n - 1, P_cond
  # and its transpose are recovered from the symmetric P only in sum, so
  # check that rows are sorted by index and the neighbors are the nearest.
  D <- as.matrix(dist(X)); diag(D) <- Inf
  for (i in c(1L, 77L, n)) {
    cols <- S$col[(S$ptr[[i]] + 1L):S$ptr[[i + 1L]]] + 1L
    expect_false(is.unsorted(cols, strictly = TRUE))
    expect_true(all(order(D[i, ])[1:31] %in% cols))
  }

  # The exact gradient and KL divergence, as in plain R, for k = 1 to 5
  # (k > 3 runs the generic code), over several row blocks.
  for (k in 1:5) {
    Y <- matrix(rnorm(n * k), n)
    g <- tsne_run(P, Y)
    r <- ref_grad(full_dense(P, n), Y)
    expect_equal(g$kl, r$kl, tolerance = 1e-12)
    expect_true(rel(g$grad, r$grad) < 1e-12)
  }

  # Barnes-Hut converges to the exact gradient as theta -> 0: with
  # theta = 0 no cell is summarized, so it is exact (with the same sparse
  # P). The error grows with theta.
  Pf <- full_csr(S, n)
  Pd <- Pf[lower.tri(Pf)]
  for (k in 1:3) {
    Y   <- matrix(rnorm(n * k) * 5, n)
    ex  <- tsne_run(Pd, Y)
    r   <- ref_grad(Pf, Y)
    expect_true(rel(ex$grad, r$grad) < 1e-12)
    err <- vapply(c(0, 0.1, 0.3, 0.6, 1), function (theta)
      rel(tsne_run(S, Y, theta = theta)$grad, ex$grad), 0)
    expect_true(err[[1]] < 1e-12)
    expect_true(all(diff(err) > 0))
    expect_true(err[[4]] < 0.05)
    expect_equal(tsne_run(S, Y, theta = 0)$kl, ex$kl, tolerance = 1e-12)
  }

  # Barnes-Hut trees with identical points, and with points too close
  # together for the tree to separate them (they interact pairwise).
  S  <- tsne_p(as.vector(d60), 60L, 5, 16L)
  Pf <- full_csr(S, 60)
  Pd <- Pf[lower.tri(Pf)]
  for (k in 1:3) {
    Y <- matrix(rnorm(60 * k), 60)
    Y[2:6, ] <- Y[rep(1, 5), ]                       # Identical.
    Y[7, ]   <- Y[8, ] + 2^-60 * (1 + abs(Y[8, ]))   # Too close.
    Y[9, ]   <- 1e70                                 # Far out.
    b <- tsne_run(S, Y, theta = 0.5)
    e <- tsne_run(Pd, Y)
    expect_true(all(is.finite(b$grad)))
    expect_true(rel(tsne_run(S, Y, theta = 0)$grad, e$grad) < 1e-10)
    expect_true(rel(b$grad, e$grad) < 0.05)
  }

  # Random starts: normal, with mean 0 and sd 1e-4, from the seed; filled
  # column by column.
  P200 <- tsne_p(as.vector(dist(matrix(rnorm(200 * 3), 200))), 200L, 30)
  z <- tsne_run(P200, NULL, start = 2L, dims = c(200L, 50L))$points
  expect_true(abs(mean(z)) < 0.03e-4 && abs(sd(as.vector(z)) / 1e-4 - 1) < 0.03)
  z1 <- tsne_run(P200, NULL, start = 2L, seed = 7L, dims = c(200L, 1L))$points
  z2 <- tsne_run(P200, NULL, start = 2L, seed = 7L, dims = c(200L, 2L))$points
  expect_identical(z2[, 1], z1[, 1])
  expect_false(identical(z2, tsne_run(P200, NULL, start = 2L, seed = 8L, dims = c(200L, 2L))$points))



  # Return value. ====

  ord <- tsne(d60, cpus = 1)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(dim(ord), c(60L, 2L))
  expect_identical(dimnames(ord), list(paste0('S', 1:60), c('TSNE1', 'TSNE2')))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'kl_divergence',
      'perplexity', 'init', 'n_iters', 'learning_rate', 'gradient', 'class') )
  expect_identical(attr(ord, 'method'), 'tsne')
  expect_identical(attr(ord, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(ord, 'prop'), c(NA_real_, NA_real_))
  expect_true(attr(ord, 'kl_divergence') > 0 && attr(ord, 'kl_divergence') < 2)
  expect_identical(attr(ord, 'n_iters'), 1000L)
  expect_identical(attr(ord, 'perplexity'), 59 / 3)   # Reduced, silently.
  expect_identical(attr(ord, 'init'), 'pcoa')
  expect_identical(attr(ord, 'learning_rate'), 50)
  expect_identical(attr(ord, 'gradient'), 'exact')

  # Centered; each axis's largest absolute coordinate is positive.
  x <- unclass(ord); attributes(x) <- list(dim = dim(x))
  expect_equal(colMeans(x), c(0, 0), tolerance = 1e-12)
  expect_true(all(apply(x, 2, function (v) v[which(abs(v) >= max(abs(v)) * (1 - 1e-8))[[1]]] > 0)))

  # A symmetric matrix gives the same result; ordination() too.
  expect_identical(tsne(as.matrix(d60), cpus = 1), ord)
  expect_identical(ordination(d60, 'tsne', cpus = 1), ord)
  expect_identical(ordination(d60, 't-SNE', cpus = 1), ord)

  # The KL divergence is that of the result.
  sc60 <- 2^-floor(log2(max(d60)))
  P60  <- tsne_p(as.vector(d60), 60L, 59 / 3, scale = sc60)
  expect_equal(tsne_run(P60, x)$kl, attr(ord, 'kl_divergence'), tolerance = 1e-12)

  # Barnes-Hut, k = 1 and 3, and options.
  o <- tsne(d60, gradient = 'barnes_hut', cpus = 1)
  expect_identical(attr(o, 'gradient'), 'barnes_hut')
  expect_true(attr(o, 'kl_divergence') < 2)
  expect_identical(colnames(tsne(d60, k = 1, cpus = 1)), 'TSNE1')
  expect_identical(colnames(tsne(d60, k = 3, cpus = 1)), paste0('TSNE', 1:3))
  expect_identical(attr(tsne(d60, k = 3, gradient = 'barnes_hut', cpus = 1), 'gradient'), 'barnes_hut')

  o <- tsne(d60, perplexity = 5, learning_rate = 100, early_exaggeration = 4, theta = 0.3, max_iters = 300, cpus = 1)
  expect_identical(attr(o, 'perplexity'), 5)
  expect_identical(attr(o, 'learning_rate'), 100)
  expect_identical(attr(o, 'n_iters'), 300L)
  expect_identical(attr(tsne(d60, early_exaggeration = 1, max_iters = 2, cpus = 1), 'learning_rate'), 50)

  # Auto: exact up to 500 samples or for k > 3, else Barnes-Hut.
  expect_identical(attr(tsne(d60, k = 4, max_iters = 10, cpus = 1), 'gradient'), 'exact')
  set.seed(2)
  d501 <- dist(matrix(rnorm(501 * 3), 501))
  o <- tsne(d501, max_iters = 5, cpus = 1)
  expect_identical(attr(o, 'gradient'), 'barnes_hut')
  expect_identical(attr(o, 'learning_rate'), 50)
  expect_identical(attr(tsne(d501, max_iters = 5, early_exaggeration = 2, cpus = 1), 'learning_rate'), 501 / 2 / 4)

  # Fewer iterations than the early exaggeration phase.
  expect_identical(attr(tsne(d60, max_iters = 20, cpus = 1), 'n_iters'), 20L)

  # Starting layouts. A matrix start with the PCoA start's values gives
  # the PCoA result.
  y0 <- tsne_run(P60, pcoa_start(as.vector(d60), 60L, 2L, 1L), start = 1L)$points
  expect_equal(sqrt(mean((y0[, 1] - mean(y0[, 1]))^2)), 1e-4, tolerance = 1e-12)
  o <- tsne(d60, init = y0, cpus = 1)
  expect_identical(attr(o, 'init'), 'matrix')
  attr(o, 'init') <- 'pcoa'
  expect_identical(o, ord)
  expect_identical(tsne(d60, init = as.data.frame(y0), cpus = 1), tsne(d60, init = y0, cpus = 1))
  r1 <- tsne(d60, init = 'random', cpus = 1)
  expect_identical(attr(r1, 'init'), 'random')
  expect_false(identical(unclass(r1)[, ], unclass(ord)[, ]))
  expect_identical(tsne(d60, init = 'random', cpus = 1), r1)
  expect_false(identical(tsne(d60, init = 'random', seed = 1, cpus = 1), r1))
  # The seed doesn't matter for a PCoA start with every axis.
  expect_identical(tsne(d60, seed = 5, cpus = 1), ord)

  # PCoA with fewer positive axes than k: a random start instead.
  d <- dist(matrix(c(0, 1, 2, 4), 4))  # On a line: one axis.
  o <- tsne(d, k = 3, max_iters = 10, cpus = 1)
  expect_identical(attr(o, 'init'), 'random')
  expect_identical(o, tsne(d, k = 3, max_iters = 10, init = 'random', cpus = 1))
  expect_true(all(is.finite(o)) && all(apply(o, 2, sd) > 0))



  # Bit-identical for any number of threads, every path. ====

  # 200 samples are enough for every parallel step to use threads. The
  # attributes (KL divergence, iterations) are compared too.

  set.seed(4)
  n  <- 200
  db <- dist(matrix(rnorm(n * 5), n) + rep(c(0, 5, 10, 15), length.out = n))

  for (k in 1:3) {
    for (g in c('exact', 'barnes_hut')) {
      cps <- if (home && k == 2) c(2L, 3L, 4L, 1500L) else many
      x1  <- tsne(db, k = k, gradient = g, max_iters = 30, cpus = 1)
      for (cp in cps)
        expect_identical(tsne(db, k = k, gradient = g, max_iters = 30, cpus = cp), x1)
    }
  }

  dv <- as.vector(db)
  for (K in c(0L, 91L))
    expect_identical(tsne_p(dv, n, 30, K, cpus = many), tsne_p(dv, n, 30, K))

  # The whole schedule (both phases), with threads.
  for (g in c('exact', 'barnes_hut')) {
    x1 <- tsne(db, gradient = g, max_iters = 260, cpus = 1)
    expect_identical(tsne(db, gradient = g, max_iters = 260, cpus = many), x1)
    expect_identical(tsne(db, gradient = g, max_iters = 260, cpus = 1), x1)
  }

  # R's random number generator is neither used nor changed.
  set.seed(5); r1 <- runif(1)
  set.seed(5); x <- tsne(d60, init = 'random', max_iters = 30, cpus = 1); r2 <- runif(1)
  expect_identical(r1, r2)
  set.seed(6); x6 <- tsne(d60, init = 'random', max_iters = 30, cpus = 1)
  expect_identical(x6, x)

  # A fresh R session gives the same results: random and PCoA starts,
  # both gradients, k = 1 to 3, with threads in the child and none here.
  # The child loads this same installation of ecodive. (Nothing that feeds
  # the result uses R's long double arithmetic, whose precision varies,
  # e.g. under valgrind.)
  if (home) {
    fresh <- function (d, cpus) {
      runs <- list(ecodive::tsne(d, init = 'random', max_iters = 30, cpus = cpus))
      for (k in 1:3) for (g in c('exact', 'barnes_hut'))
        runs <- c(runs, list(ecodive::tsne(d, k = k, gradient = g, max_iters = 30, cpus = cpus)))
      runs
    }
    environment(fresh) <- globalenv()  # Don't save this test's objects.
    pkg  <- getNamespaceInfo('ecodive', 'path')
    f    <- gsub('\\', '/', tempfile(fileext = '.rds'), fixed = TRUE)
    saveRDS(list(d = db, fresh = fresh, lib = dirname(pkg)), f)
    code <- sprintf(paste(
      'x <- readRDS("%s"); .libPaths(c(x$lib, .libPaths()));',
      'loadNamespace("ecodive", lib.loc = x$lib);',
      'saveRDS(list(path = getNamespaceInfo("ecodive", "path"), runs = x$fresh(x$d, 3L)), "%s")' ),
      f, f )
    out <- system2(file.path(R.home('bin'), 'Rscript'), c('--vanilla', '-e', shQuote(code)),
                   stdout = FALSE, stderr = FALSE)
    expect_identical(out, 0L)  # (A failure here is reported, not skipped.)
    if (identical(out, 0L)) {
      child <- readRDS(f)
      expect_identical(normalizePath(child$path), normalizePath(pkg))
      expect_identical(child$runs, fresh(db, 1L))
    }
    unlink(f)
  }



  # Small and degenerate inputs. ====

  # No samples, or one: all at 0, with `k` columns.
  o <- tsne(dist(matrix(numeric(0), 0, 2)))
  expect_identical(dim(o), c(0L, 2L))
  expect_identical(attr(o, 'init'), 'none')
  o <- tsne(dist(matrix(1:2, 1, dimnames = list('A', NULL))), k = 3)
  expect_identical(dim(o), c(1L, 3L))
  expect_identical(unname(unclass(o)[1, ]), c(0, 0, 0))
  expect_identical(rownames(o), 'A')
  expect_identical(attr(o, 'eig'), rep(NA_real_, 3))
  expect_identical(attr(o, 'n_iters'), 0L)
  expect_identical(attr(o, 'kl_divergence'), 0)
  expect_identical(attr(o, 'perplexity'), 1)
  expect_identical(attr(o, 'init'), 'none')

  # 2, 3, and 4 samples: perplexity 1, `k` as given, finite results.
  for (n in 2:4) {
    d <- dist(matrix(c(0, 1, 3, 7)[seq_len(n)], n))
    for (g in c('exact', 'barnes_hut')) for (k in 1:3) {
      o <- tsne(d, k = k, gradient = g, cpus = 1)
      expect_true(all(is.finite(o)))
      expect_identical(attr(o, 'perplexity'), 1)
      expect_identical(dim(o), c(n, k))
    }
  }
  expect_identical(dim(tsne(d60, k = 70, max_iters = 2, cpus = 1)), c(60L, 70L))

  # Two samples: P and Q are both 1/2, so the gradient vanishes after the
  # early exaggeration phase and the run stops at the first check. A huge
  # max_iters is fine.
  o <- tsne(dist(1:2), max_iters = 1e12, cpus = 1)
  expect_identical(attr(o, 'n_iters'), 300L)
  expect_true(abs(attr(o, 'kl_divergence')) < 1e-12)

  # Duplicate samples (zero distances) stay together (up to rounding: their
  # affinities are summed in different orders); the rest is fine.
  dup <- rbind(m60, m60[1:5, ])
  rownames(dup) <- paste0('S', 1:65)
  for (g in c('exact', 'barnes_hut')) {
    o <- tsne(bray(dup), gradient = g, cpus = 1)
    expect_true(all(is.finite(o)))
    expect_true(max(abs(unclass(o)[61:65, ] - unclass(o)[1:5, ])) < 1e-3 * max(abs(o)))
    o <- tsne(bray(dup), gradient = g, init = 'random', cpus = 1)
    expect_true(all(is.finite(o)))
  }

  # All distances zero, or all equal: no structure, but finite. (With all
  # zero, the PCoA start has no positive eigenvalue, so is random.)
  for (v in c(0, 1)) {
    d <- dist_of(rep(v, 45))
    for (g in c('exact', 'barnes_hut')) {
      o <- tsne(d, gradient = g, cpus = 1)
      expect_true(all(is.finite(o)))
      if (v == 0) expect_identical(o, tsne(d, gradient = g, init = 'random', cpus = 1))
    }
  }

  # A runaway learning rate is an error, not a layout of huge or NaN values,
  # for both gradients.
  for (g in c('exact', 'barnes_hut')) {
    expect_error(tsne(d60, gradient = g, learning_rate = 1e300, cpus = 1), 'diverged')
    expect_error(tsne(d60, gradient = g, early_exaggeration = 1e305, cpus = 1), 'diverged')
    expect_error(tsne(d60, gradient = g, learning_rate = 1e300, max_iters = 20, cpus = 1), 'diverged')
  }

  # The scale of the distances doesn't matter: multiplying them by a power
  # of 2 gives exactly the same result, and by any other factor, the same
  # affinities to within the binary search's tolerance (P is invariant to
  # scale; its search starts at a precision of 1, so distances are scaled to
  # a largest one of 1 to 2 first). Far-off scales would otherwise give
  # uniform affinities.
  d   <- bray(m60[1:45, ])
  ref <- tsne(d, perplexity = 10, cpus = 1)
  for (e in c(-200, -40, 40, 200))
    expect_identical(tsne(d * 2^e, perplexity = 10, cpus = 1), ref)
  sc <- function (v) 2^-floor(log2(max(v)))
  P0 <- tsne_p(as.vector(d), 45L, 10, scale = sc(d))
  S0 <- tsne_p(as.vector(d), 45L, 4, 13L, scale = sc(d))
  for (s in 10^c(-300, -30, -8, 14, 76, 100)) {
    ds <- d * s
    P  <- tsne_p(as.vector(ds), 45L, 10, scale = sc(ds))
    S  <- tsne_p(as.vector(ds), 45L, 4, 13L, scale = sc(ds))
    expect_true(rel(P, P0) < 1e-4)
    expect_true(rel(S$val, S0$val) < 1e-4)
  }
  expect_true(all(is.finite(tsne(dist(matrix(rnorm(40), 20)) * 1e100, cpus = 1))))

  # Rows that scikit-learn's search can't calibrate in its 100 steps from a
  # precision of 1: three groups whose distances are about 1e-20 of the
  # largest (set by two outliers). The fallback search still separates them.
  set.seed(7)
  lab <- rep(1:3, each = 12)
  Xg  <- rbind(matrix(rnorm(36 * 2), 36) * 1e-20 + lab * 1e-19, c(1, 0), c(0, 1))
  o   <- tsne(dist(Xg), perplexity = 5, cpus = 1)
  expect_identical(knn_acc(unclass(o)[1:36, ], lab), 1)
  for (K in c(0L, 15L)) {
    P <- tsne_p(as.vector(dist(Xg)), 38L, 5, K)
    expect_true(all(is.finite(if (K) P$val else P)))
  }



  # Invalid input. ====

  expect_error(tsne(d60, k = 0),        '`k` must be a whole number of at least 1')
  expect_error(tsne(d60, k = 1.5),      '`k` must be a whole number')
  expect_error(tsne(d60, k = NA),       '`k` must be a whole number')
  expect_error(tsne(d60, k = 'a'),      '`k` must be a whole number')
  expect_error(tsne(d60, k = 1:2),      '`k` must be a whole number')
  expect_error(tsne(d60, k = 1e10),     '`k` is too large')
  expect_error(tsne(d60, perplexity = 0.5), '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, perplexity = -1),  '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, perplexity = NA),  '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, perplexity = Inf), '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, perplexity = c(5, 6)), '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, perplexity = '5'), '`perplexity` must be a number of at least 1')
  expect_error(tsne(d60, theta = -0.1), '`theta` must be a number from 0 to 1')
  expect_error(tsne(d60, theta = 1.1),  '`theta` must be a number from 0 to 1')
  expect_error(tsne(d60, theta = NaN),  '`theta` must be a number from 0 to 1')
  expect_error(tsne(d60, max_iters = 0),   '`max_iters` must be a whole number of at least 1')
  expect_error(tsne(d60, max_iters = 2.5), '`max_iters` must be a whole number of at least 1')
  expect_error(tsne(d60, max_iters = Inf), '`max_iters` must be a whole number of at least 1')
  expect_error(tsne(d60, early_exaggeration = 0.5), '`early_exaggeration` must be a number of at least 1')
  expect_error(tsne(d60, early_exaggeration = Inf), '`early_exaggeration` must be a number of at least 1')
  expect_error(tsne(d60, learning_rate = 0),      '`learning_rate` must be "auto" or a number greater than 0')
  expect_error(tsne(d60, learning_rate = -5),     '`learning_rate` must be "auto" or a number greater than 0')
  expect_error(tsne(d60, learning_rate = 'fast'), '`learning_rate` must be "auto" or a number greater than 0')
  expect_error(tsne(d60, learning_rate = NA),     '`learning_rate` must be "auto" or a number greater than 0')
  expect_error(tsne(d60, learning_rate = Inf),    '`learning_rate` must be "auto" or a number greater than 0')
  expect_error(tsne(d60, gradient = 'fast'), '`gradient` must be "auto", "exact", or "barnes_hut"')
  expect_error(tsne(d60, gradient = 1),      '`gradient` must be "auto", "exact", or "barnes_hut"')
  expect_error(tsne(d60, gradient = NA),     '`gradient` must be "auto", "exact", or "barnes_hut"')
  expect_error(tsne(d60, gradient = c('exact', 'barnes_hut')), '`gradient` must be')
  expect_error(tsne(d60, gradient = 'barnes_hut', k = 4), 'at most k = 3')
  expect_error(tsne(d60, init = 'pca'),       '`init` must be "pcoa", "random", or a numeric matrix')
  expect_error(tsne(d60, init = NA),          '`init` must be "pcoa", "random", or a numeric matrix')
  expect_error(tsne(d60, init = 1:120),       '`init` must be "pcoa", "random", or a numeric matrix')
  expect_error(tsne(d60, init = NULL),        '`init` must be "pcoa", "random", or a numeric matrix')
  expect_error(tsne(d60, init = c('pcoa', 'random')), '`init` must be')
  expect_error(tsne(d60, init = matrix(0, 60, 3)), '`init` must be a numeric matrix with 60 rows and 2 columns')
  expect_error(tsne(d60, init = matrix(0, 59, 2)), '`init` must be a numeric matrix with 60 rows and 2 columns')
  expect_error(tsne(d60, init = matrix(NA_real_, 60, 2)), '`init` must be a numeric matrix')
  expect_error(tsne(d60, init = matrix('a', 60, 2)), '`init` must be a numeric matrix')
  expect_error(tsne(d60, init = matrix(0, 60, 2)), 'must not have a column whose values are all the same')
  expect_error(tsne(d60, init = cbind(1:60, 5)),   'must not have a column whose values are all the same')
  expect_error(tsne(d60, seed = 1.5),  '`seed` must be an integer')
  expect_error(tsne(d60, seed = NA),   '`seed` must be an integer')
  expect_error(tsne(d60, cpus = 0),    '`cpus` must be an integer greater than 0')
  expect_error(tsne(d60, cpus = 1.5),  '`cpus` must be an integer greater than 0')
  expect_error(tsne(d60, cpus = NA),   '`cpus` must be an integer greater than 0')

  dna <- d60; dna[3] <- NA
  expect_error(tsne(dna), '`dm` must not contain NA, NaN, or infinite values')
  dna[3] <- Inf
  expect_error(tsne(dna), '`dm` must not contain NA, NaN, or infinite values')
  dna[3] <- NaN
  expect_error(tsne(dna), '`dm` must not contain NA, NaN, or infinite values')
  dna[3] <- -0.1
  expect_error(tsne(dna), '`dm` must not contain negative distances')
  dna[3] <- 1e200
  expect_error(tsne(dna), 'too large')

  expect_error(tsne(matrix(1:6, 2)),           'symmetric numeric matrix')
  expect_error(tsne(matrix(c(0, 1, 2, 0), 2)), 'symmetric numeric matrix')
  expect_error(tsne(matrix('a', 2, 2)),        'symmetric numeric matrix')
  expect_error(tsne(structure(1:2, Size = 3L, class = 'dist')), 'not a valid dist object')
  expect_error(tsne(list()),                   'symmetric numeric matrix')
  expect_error(tsne(NULL),                     'symmetric numeric matrix')

  # A perplexity too large for n is reduced, with a warning only if given.
  expect_warning(o <- tsne(d60, perplexity = 25, max_iters = 10, cpus = 1), 'reduced from 25 to 19.6667')
  expect_identical(attr(o, 'perplexity'), 59 / 3)
  expect_silent(tsne(d60, max_iters = 10, cpus = 1))

  # Integer storage is converted.
  di <- dist(matrix(c(1L, 5L, 2L, 8L, 9L), 5))
  storage.mode(di) <- 'integer'
  expect_true(all(is.finite(tsne(di, cpus = 1))))

#})
