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

  # Fraction of each sample's 10 nearest neighbours in `d` that are among
  # its 10 nearest in the embedding.
  nn_preservation <- function (d, x, m = 10) {
    dm <- as.matrix(d)
    de <- as.matrix(dist(x))
    mean(sapply(seq_len(nrow(dm)), function (i)
      length(intersect(order(dm[, i])[2:(m + 1)], order(de[, i])[2:(m + 1)])) / m ))
  }


  #--------------------------------------------------
  # Neighbours, smooth kNN, and fuzzy graph: exact.
  #--------------------------------------------------
  g1 <- graph(c1$d, c1$n_neighbors)

  expect_identical(g1$idx,   c1$nn_idx)
  expect_identical(g1$dist,  c1$nn_dist)
  expect_identical(g1$sigma, c1$sigma)
  expect_identical(g1$rho,   c1$rho)
  expect_identical(g1$p,     c1$graph$p)
  expect_identical(g1$j,     c1$graph$i)
  expect_identical(g1$x,     c1$graph$x)
  expect_identical(g1$n_failures,   0L)
  expect_identical(g1$n_components, 1L)

  # The same for any number of threads.
  expect_identical(graph(c1$d, c1$n_neighbors, cpus = 4L), g1)

  # Ties, duplicate samples, partial local connectivity, mixed union.
  g2 <- graph(c2$d, c2$n_neighbors, local = c2$local_connectivity, mix = c2$set_op_mix_ratio)

  expect_identical(g2$idx,   c2$nn_idx)
  expect_identical(g2$dist,  c2$nn_dist)
  expect_identical(g2$sigma, c2$sigma)
  expect_identical(g2$rho,   c2$rho)
  expect_identical(g2$p,     c2$graph$p)
  expect_identical(g2$j,     c2$graph$i)
  expect_identical(g2$x,     c2$graph$x)

  # Pure fuzzy intersection keeps only mutual neighbours.
  g3 <- graph(c2$d, c2$n_neighbors, mix = 0)
  expect_true(length(g3$x) < length(graph(c2$d, c2$n_neighbors)$x))
  expect_true(all(g3$x > 0))


  #--------------------------------------------------
  # Curve parameters.
  #--------------------------------------------------
  expect_equal(umap_ab(1, 0.01),    c1$ab,      tolerance = 1e-10)
  expect_equal(umap_ab(1.5, 0.3),   c1$sgd2$ab, tolerance = 1e-10)
  expect_error(umap_ab(1, 1e10))


  #--------------------------------------------------
  # Spectral start: the normalized Laplacian's
  # eigenvectors.
  #--------------------------------------------------
  sp <- .Call(C_umap_spectral, g1, 2L, 1L)
  expect_true(sp$converged)
  expect_equal(sp$vectors, c1$spectral$vectors, tolerance = 1e-8)
  expect_equal(sp$values,  c1$spectral$values,  tolerance = 1e-8)
  expect_identical(.Call(C_umap_spectral, g1, 2L, 4L), sp)


  #--------------------------------------------------
  # Layout: from the same start and per-epoch seeds,
  # uwot's single-threaded SGD gives the same result.
  # Bit-identical with the same compiler and math
  # library; other platforms' powf() can differ in the
  # last bit, and SGD amplifies that.
  #--------------------------------------------------
  same_platform <- .Platform$OS.type == 'windows' && R.version$arch == 'x86_64'
  layout_tol    <- if (same_platform) 0 else 0.5

  s1 <- c1$sgd1
  x1 <- .Call(C_umap_sgd, g1, s1$init, as.integer(s1$n_epochs), c(c1$ab, 1, 1, 5), s1$epoch_seeds)
  x1 <- center_like(x1, s1$embedding)
  expect_equal(x1, s1$embedding, tolerance = layout_tol)

  s2 <- c1$sgd2
  a2 <- s2$args
  g12 <- graph(c1$d, a2$n_neighbors)
  x2  <- .Call(
    C_umap_sgd, g12, s2$init, as.integer(a2$n_epochs),
    c(s2$ab, a2$repulsion_strength, a2$learning_rate, a2$negative_sample_rate),
    s2$epoch_seeds )
  x2 <- center_like(x2, s2$embedding)
  expect_equal(x2, s2$embedding, tolerance = layout_tol)

  # umap() with init = matrix runs the same layout, from its own seeds.
  u2 <- do.call(umap, c(list(dm = c1$d, k = 3, init = s2$init, seed = 9), a2))
  x9 <- .Call(
    C_umap_sgd, g12, s2$init, as.integer(a2$n_epochs),
    c(s2$ab, a2$repulsion_strength, a2$learning_rate, a2$negative_sample_rate), 9L )
  x9 <- center_like(x9, unclass(u2)[, ])
  expect_equal(unname(unclass(u2)[, ]), x9, tolerance = 0)
  expect_identical(attr(u2, 'init'), 'matrix')


  #--------------------------------------------------
  # Defaults: as good as uwot's at preserving
  # neighbours.
  #--------------------------------------------------
  u <- umap(c1$d)
  expect_true(nn_preservation(c1$d, u) > c1$preservation - 0.05)

  expect_identical(dim(u), c(n1, 2L))
  expect_identical(colnames(u), c('UMAP1', 'UMAP2'))
  expect_identical(
    names(attributes(u)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'a', 'b', 'init', 'n_epochs') )
  expect_identical(attr(u, 'method'), 'umap')
  expect_identical(attr(u, 'eig'),  c(NA_real_, NA_real_))
  expect_identical(attr(u, 'prop'), c(NA_real_, NA_real_))
  expect_identical(attr(u, 'init'), 'spectral')
  expect_identical(attr(u, 'n_epochs'), 500L)
  expect_equal(c(a = attr(u, 'a'), b = attr(u, 'b')), c1$ab)
  expect_equal(unname(colMeans(u)), c(0, 0))

  # Deterministic, for any number of cpus, and R's RNG untouched.
  set.seed(1)
  rs <- .Random.seed
  expect_identical(umap(c1$d, cpus = 1), u)
  expect_identical(umap(c1$d, cpus = 4), u)
  expect_identical(.Random.seed, rs)
  expect_false(identical(umap(c1$d, seed = 1), u))

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

#})
