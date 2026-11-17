#test_that("umap", {

  # Reference results, from data-raw/expected/umap.r.
  ref <- readRDS('expected/umap.rds')

  c1 <- ref$case1
  c2 <- ref$case2
  n1 <- attr(c1$d, 'Size')


  # Run the graph stage.
  graph <- function (d, k, local = 1, bw = 1, mix = 1, cpus = 1L)
    .Call(C_umap_graph, as.double(d), attr(d, 'Size'), as.integer(k), c(local, bw, mix), cpus)

  # Center, as uwot does, then match each axis's sign to `y`.
  center_like <- function (x, y) {
    x <- unname(scale(x, scale = FALSE))
    attributes(x) <- list(dim = dim(x))
    sweep(x, 2, sign(colSums(x * y)), '*')
  }

  # The per-vertex seeds that uwot's batch mode draws from R's RNG after
  # set.seed(seed): two per vertex per epoch.
  batch_seeds <- function (seed, n, n_epochs) {
    set.seed(seed)
    floor(runif(2 * n * n_epochs) * 4294967295)
  }

  # Layout quality, as in data-raw/expected/umap.r.

  # Fraction of each sample's m nearest neighbours in `d` that are among
  # its m nearest in the embedding `x`.
  nn_preservation <- function (d, x, m = 10) {
    dm <- as.matrix(d)
    de <- as.matrix(dist(x))
    mean(sapply(seq_len(nrow(dm)), function (i)
      length(intersect(order(dm[, i])[2:(m + 1)], order(de[, i])[2:(m + 1)])) / m ))
  }

  # Trustworthiness (Venna & Kaski 2001) of the embedding's m nearest
  # neighbours.
  trustworthiness <- function (d, x, m = 10) {
    dm <- as.matrix(d)
    de <- as.matrix(dist(x))
    n  <- nrow(dm)
    s  <- sum(sapply(seq_len(n), function (i) {
      rk <- integer(n)
      rk[order(dm[, i])] <- 0:(n - 1)
      nn <- order(de[, i])[2:(m + 1)]
      sum(pmax(rk[nn] - m, 0))
    }))
    1 - 2 * s / (n * m * (2 * n - 3 * m - 1))
  }

  # UMAP's fuzzy cross-entropy between the graph and the embedding's
  # similarities 1 / (1 + a d^(2b)), over all pairs.
  cross_entropy <- function (graph, x, ab) {
    n <- nrow(x)
    v <- matrix(0, n, n)
    v[cbind(graph$i + 1, rep(seq_len(n), diff(graph$p)))] <- graph$x
    q <- 1 / (1 + ab[[1]] * as.matrix(dist(x))^(2 * ab[[2]]))
    q <- pmin(pmax(q, 1e-12), 1 - 1e-12)
    ce <- -(v * log(q) + (1 - v) * log(1 - q))
    sum(ce[lower.tri(ce)])
  }

  # Adam's settings in uwot's batch mode.
  adam <- c(0.5, 0.9, 1e-7)

  # Large graphs' reference results are fingerprints of each part, as in
  # data-raw/expected/umap.r.
  fingerprint <- function (x) {
    x <- as.double(x)
    c(length(x), sum(x), sum(x * x), sum(x * cos(seq_along(x))))
  }

  graph_matches <- function (g, ref, info) {
    parts <- c(nn_idx = 'idx', nn_dist = 'dist', sigma = 'sigma', rho = 'rho', p = 'p', i = 'j', x = 'x')
    for (nm in names(parts))
      expect_equal(fingerprint(g[[parts[[nm]]]]), ref[[nm]], tolerance = 1e-12, info = paste(info, nm))
  }


  # Comparisons with uwot's stored results run at home only; the stored
  # distances are also used as inputs everywhere.
  home <- at_home()

  # Threads for the determinism checks below, which call the C code directly
  # (past validate_cpus()): at most 2 on CRAN.
  many <- if (home) 4L else 2L


  #--------------------------------------------------
  # Neighbours, smooth kNN, and fuzzy graph: exact.
  #--------------------------------------------------
  g1 <- graph(c1$d, c1$n_neighbors)

  if (home) {
    expect_identical(g1$idx,   c1$nn_idx)
    expect_identical(g1$dist,  c1$nn_dist)
    expect_identical(g1$sigma, c1$sigma)
    expect_identical(g1$rho,   c1$rho)
    expect_identical(g1$p,     c1$graph$p)
    expect_identical(g1$j,     c1$graph$i)
    expect_equal(g1$x, c1$graph$x, tolerance = 1e-14)  # exp() differs by an ulp across libms
  }
  expect_identical(g1$n_failures,   0L)
  expect_identical(g1$n_components, 1L)

  # The same for any number of threads.
  expect_identical(graph(c1$d, c1$n_neighbors, cpus = many), g1)

  # Ties, duplicate samples, partial local connectivity, mixed union.
  g2 <- graph(c2$d, c2$n_neighbors, local = c2$local_connectivity, mix = c2$set_op_mix_ratio)

  if (home) {
    expect_identical(g2$idx,   c2$nn_idx)
    expect_identical(g2$dist,  c2$nn_dist)
    expect_identical(g2$sigma, c2$sigma)
    expect_identical(g2$rho,   c2$rho)
    expect_identical(g2$p,     c2$graph$p)
    expect_identical(g2$j,     c2$graph$i)
    expect_equal(g2$x, c2$graph$x, tolerance = 1e-14)  # exp() differs by an ulp across libms
  }

  # Pure fuzzy intersection keeps only mutual neighbours.
  g3 <- graph(c2$d, c2$n_neighbors, mix = 0)
  expect_true(length(g3$x) < length(graph(c2$d, c2$n_neighbors)$x))
  expect_true(all(g3$x > 0))


  # Other numbers of neighbours, bandwidths and local connectivity; and
  # Bray-Curtis distances between 600 samples.
  if (home) {
    for (ref_g in ref$more$graphs) {
      a <- ref_g$args
      g <- graph(
        c1$d, ref_g$n_neighbors,
        local = if (is.null(a$local_connectivity)) 1 else a$local_connectivity,
        bw    = if (is.null(a$bandwidth)) 1 else a$bandwidth )
      graph_matches(g, ref_g, paste(ref_g$n_neighbors, unlist(a), collapse = ' '))
    }
    g600 <- graph(bray(ref$more$counts600), 15L, cpus = 2L)
    graph_matches(g600, ref$more$graph600, 'bray 600')
  }


  #--------------------------------------------------
  # Curve parameters.
  #--------------------------------------------------
  # ecodive's Levenberg-Marquardt fit against uwot's nls() fit of the same
  # curve, which stops at a relative tolerance of 1e-5.
  if (home) {
    expect_equal(umap_ab(1, 0.01),    c1$ab,      tolerance = 1e-6)
    expect_equal(umap_ab(1.5, 0.3),   c1$sgd2$ab, tolerance = 1e-6)
  }
  expect_error(umap_ab(1, 1e10))


  #--------------------------------------------------
  # Spectral start: the normalized Laplacian's
  # eigenvectors.
  #--------------------------------------------------
  sp <- .Call(C_umap_spectral, g1, 2L, 1L)
  expect_true(sp$converged)
  if (home) {
    expect_equal(sp$vectors, c1$spectral$vectors, tolerance = 1e-8)
    expect_equal(sp$values,  c1$spectral$values,  tolerance = 1e-8)
  }
  expect_identical(.Call(C_umap_spectral, g1, 2L, many), sp)


  #--------------------------------------------------
  # Layout: from the same start and seeds, uwot's batch
  # mode (batch = TRUE) gives the same result.
  # Bit-identical with the same compiler and math
  # library (Windows x86_64, where the fixture was
  # made). Elsewhere, powf() and FMA contraction can
  # differ in the last bit: a 5-epoch layout must still
  # agree to 1e-5, but over 60-100 epochs SGD amplifies
  # such differences (to a mean relative difference of
  # 0.6 on arm64), so long layouts are compared only on
  # the same platform; the quality checks below run
  # everywhere. Results across `cpus` are always
  # exactly identical.
  #--------------------------------------------------
  same_platform <- home && .Platform$OS.type == 'windows' && R.version$arch == 'x86_64'
  short_tol     <- if (same_platform) 0 else 1e-5

  s1    <- c1$sgd1
  seeds <- batch_seeds(s1$seed, n1, 5)
  xs    <- .Call(C_umap_sgd, g1, s1$init, 5L, c(c1$ab, 1, 1, 5, adam), seeds, 1L)
  expect_identical(.Call(C_umap_sgd, g1, s1$init, 5L, c(c1$ab, 1, 1, 5, adam), seeds, many), xs)
  if (home) expect_equal(center_like(xs, s1$short), s1$short, tolerance = short_tol)

  seeds <- batch_seeds(s1$seed, n1, s1$n_epochs)
  x1    <- .Call(C_umap_sgd, g1, s1$init, as.integer(s1$n_epochs), c(c1$ab, 1, 1, 5, adam), seeds, 1L)
  expect_identical(.Call(C_umap_sgd, g1, s1$init, as.integer(s1$n_epochs), c(c1$ab, 1, 1, 5, adam), seeds, many), x1)
  x1 <- center_like(x1, s1$embedding)
  if (same_platform) expect_equal(x1, s1$embedding, tolerance = 0)

  s2    <- c1$sgd2
  a2    <- s2$args
  g12   <- graph(c1$d, a2$n_neighbors)
  opt2  <- c(s2$ab, a2$repulsion_strength, a2$learning_rate, a2$negative_sample_rate, adam)
  seeds <- batch_seeds(s2$seed, n1, a2$n_epochs)
  x2    <- .Call(C_umap_sgd, g12, s2$init, as.integer(a2$n_epochs), opt2, seeds, 2L)
  x2    <- center_like(x2, s2$embedding)
  if (same_platform) expect_equal(x2, s2$embedding, tolerance = 0)

  # 3 dimensions, on case 2's Bray-Curtis graph.
  s3    <- ref$more$sgd3
  seeds <- batch_seeds(s3$seed, 40L, 5L)
  x3    <- .Call(C_umap_sgd, g2, s3$init, 5L, c(s3$ab, 1, 1, 5, adam), seeds, 1L)
  if (home) expect_equal(center_like(x3, s3$epochs5), s3$epochs5, tolerance = short_tol)
  seeds <- batch_seeds(s3$seed, 40L, 50L)
  x3    <- .Call(C_umap_sgd, g2, s3$init, 50L, c(s3$ab, 1, 1, 5, adam), seeds, 2L)
  if (same_platform) expect_equal(center_like(x3, s3$epochs50), s3$epochs50, tolerance = 0)

  # umap() with init = matrix runs the same layout, from its own seeds.
  u2 <- do.call(umap, c(list(dm = c1$d, k = 3, init = s2$init, seed = 9), a2))
  x9 <- .Call(C_umap_sgd, g12, s2$init, as.integer(a2$n_epochs), opt2, 9L, 1L)
  x9 <- center_like(x9, unclass(u2)[, ])
  expect_equal(unname(unclass(u2)[, ]), x9, tolerance = 0)
  expect_identical(attr(u2, 'init'), 'matrix')

  # A constant column would leave every sample at the same place.
  expect_error(umap(c1$d, init = matrix(0, n1, 2)), 'all the same')
  expect_error(umap(c1$d, init = cbind(seq_len(n1), 1)), 'all the same')


  #--------------------------------------------------
  # Quality: from uwot's spectral start, as good as
  # uwot's batch mode at preserving neighbours, and
  # with a similar cross-entropy.
  #--------------------------------------------------
  qb <- c1$quality['batch', ]
  if (home) {
    uq <- umap(c1$d, init = c1$init_q, min_dist = 0.01)
    expect_true(nn_preservation(c1$d, uq) > qb[['preservation']] - 0.05)
    expect_true(trustworthiness(c1$d, uq) > qb[['trustworthiness']] - 0.02)
    expect_true(cross_entropy(c1$graph, uq, c1$ab) < qb[['cross_entropy']] * 1.1)
  }


  #--------------------------------------------------
  # Defaults.
  #--------------------------------------------------
  u <- umap(c1$d)
  if (home) expect_true(nn_preservation(c1$d, u) > qb[['preservation']] - 0.1)

  expect_identical(dim(u), c(n1, 2L))
  expect_identical(colnames(u), c('UMAP1', 'UMAP2'))
  expect_identical(
    names(attributes(u)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'a', 'b', 'init', 'n_epochs', 'n_neighbors', 'class') )
  expect_identical(attr(u, 'method'), 'umap')
  expect_identical(attr(u, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(u, 'prop'), c(NA_real_, NA_real_))
  expect_identical(attr(u, 'init'), 'spectral')
  expect_identical(attr(u, 'n_epochs'), 500L)
  expect_equal(c(a = attr(u, 'a'), b = attr(u, 'b')), umap_ab(1, 0.1))
  expect_equal(unname(colMeans(u)), c(0, 0))

  # Deterministic, for any number of cpus, and R's RNG untouched.
  set.seed(1)
  rs <- .Random.seed
  expect_identical(umap(c1$d, cpus = 1), u)
  expect_identical(umap(c1$d, cpus = 4), u)
  expect_identical(.Random.seed, rs)
  expect_false(identical(umap(c1$d, seed = 1), u))

  # With enough samples for several threads to share the work.
  set.seed(2)
  d5 <- dist(matrix(rnorm(700 * 4), 700))
  u5 <- umap(d5, n_epochs = 30, cpus = 1)
  expect_identical(umap(d5, n_epochs = 30, cpus = 3), u5)
  expect_identical(umap(d5, n_epochs = 30, cpus = 8), u5)

  # Through ordination(), and from a matrix.
  expect_identical(ordination(c1$d, 'umap'), u)
  expect_identical(match_ordination('UMAP')$id, 'umap')
  m <- as.matrix(c2$d)
  dimnames(m) <- list(paste0('S', 1:40), paste0('S', 1:40))
  um <- umap(m, n_neighbors = 10, n_epochs = 20)
  expect_identical(rownames(um), paste0('S', 1:40))


  #--------------------------------------------------
  # Other starts.
  #--------------------------------------------------
  expect_identical(attr(umap(c1$d, init = 'pcoa',   n_epochs = 5), 'init'), 'pcoa')
  expect_identical(attr(umap(c1$d, init = 'random', n_epochs = 5), 'init'), 'random')

  # Three separate groups: unconnected graph, PCoA start.
  grp <- cbind(0:9, (0:9)^2 %% 7)
  far <- dist(rbind(grp, grp + 100, cbind(grp[, 1] + 1000, grp[, 2])))
  uf  <- umap(far, n_neighbors = 5, n_epochs = 10)
  expect_identical(attr(uf, 'init'), 'pcoa')
  expect_true(all(is.finite(uf)))

  # The PCoA start scaled to standard deviation 1: uwot's "spca" fallback.
  uf <- umap(dist(ref$more$far), n_neighbors = 5, n_epochs = 0)
  expect_identical(attr(uf, 'init'), 'pcoa')
  if (home) {
    y <- unname(scale(ref$more$far0, scale = FALSE)[, ])
    expect_true(max(abs(center_like(unclass(uf)[, ], y) - y)) / max(abs(y)) < 1e-10)
  }

  # No optimization: the (centered) start.
  u0 <- umap(c1$d, init = s1$init, n_epochs = 0)
  expect_equal(abs(unname(unclass(u0)[, ])), abs(unname(scale(s1$init, scale = FALSE)[, ])))


  #--------------------------------------------------
  # Small inputs.
  #--------------------------------------------------
  d3 <- dist(matrix(c(0, 1, 3), 3))
  expect_identical(attr(umap(d3, n_epochs = 10), 'init'), 'random')   # k > n - 2
  expect_identical(dim(umap(d3, k = 1, n_epochs = 10)), c(3L, 1L))
  expect_identical(dim(umap(dist(1:2), n_epochs = 10)), c(2L, 2L))
  expect_identical(dim(umap(dist(1))), c(1L, 2L))
  expect_identical(attr(umap(dist(1)), 'init'), 'none')
  expect_identical(attr(umap(dist(1), init = matrix(5, 1, 2)), 'init'), 'matrix')
  expect_identical(dim(umap(structure(numeric(0), Size = 0L, class = 'dist'))), c(0L, 2L))


  #--------------------------------------------------
  # Invalid arguments.
  #--------------------------------------------------
  expect_error(umap(c1$d, k = 0))
  expect_error(umap(c1$d, n_neighbors = 1))
  expect_error(umap(c1$d, min_dist = -1))
  expect_error(umap(c1$d, spread = 0))
  expect_error(umap(c1$d, n_epochs = -1))
  expect_error(umap(c1$d, learning_rate = 0))
  expect_error(umap(c1$d, set_op_mix_ratio = 2))
  expect_error(umap(c1$d, local_connectivity = 0.5))
  expect_error(umap(c1$d, init = 'pca'))
  expect_error(umap(c1$d, init = matrix(0, 3, 2)))
  expect_error(umap(c1$d, seed = 0.5))
  expect_error(umap(c1$d, cpus = 0))
  expect_error(umap(c(1, NA, 3)))
  expect_error(umap(dist(c(1, 2, NA))))
  bad <- c1$d
  bad[1] <- -1
  expect_error(umap(bad))
  expect_error(umap(c1$d, negative_sample_rate = Inf), '`negative_sample_rate` must be')
  expect_error(umap(c1$d, n_epochs = Inf), '`n_epochs` must be')

  # Regressions: a tiny negative_sample_rate made the number of negative
  # samples (size_t) -1, a loop without end; a diverged layout failed in
  # which.max() instead of with a clear error.
  un <- umap(c2$d, n_neighbors = 10, n_epochs = 5, negative_sample_rate = 1e-9)
  expect_true(all(is.finite(un)))
  expect_error(umap(c2$d, n_neighbors = 10, n_epochs = 5, learning_rate = 1e300), 'UMAP diverged')

#})
