#test_that("isomap", {

  # Threads for direct .Call()s: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  # isomap() itself caps `cpus` with n_cpus().
  max_cpus <- if (at_home()) 4L else 2L

  # Largest absolute difference between two coordinate matrices, relative to
  # the largest coordinate, after matching each axis's sign.
  coord_diff <- function (x, y) {
    x <- matrix(x, nrow(x))
    y <- matrix(y, nrow(y))
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # Largest elementwise relative difference; zeros must match exactly.
  rel_diff <- function (x, y) {
    if (any(x[y == 0] != 0)) return (Inf)
    nz <- y != 0
    if (!any(nz)) return (0)
    max(abs(x[nz] - y[nz]) / abs(y[nz]))
  }

  # The coordinates alone.
  vals <- function (x) { x <- unclass(x); attributes(x) <- list(dim = dim(x)); x }

  # join = c(bytes, ratio, floor); isomap() uses c(16 * 2^20, 10, 65536).
  geodesics <- function (d, n_neighbors = 0L, radius = NA_real_, cpus = 1L, join = c(16 * 2^20, 10, 65536))
    .Call(C_isomap, as.double(d), attr(d, 'Size'), as.integer(n_neighbors),
          as.double(radius), as.integer(cpus), as.double(join), TRUE)



  # scikit-learn reference results (see data-raw/expected/isomap.r). ====

  if (at_home()) {

    expected <- readRDS('expected/isomap.rds')

    for (fx in expected) {

      name  <- fx$name
      knn   <- !is.na(fx$n_neighbors)
      split <- fx$components > 1

      # Neighbor graph, including joining edges, and components.
      g <- if (knn) geodesics(fx$dm, n_neighbors = fx$n_neighbors)
           else     geodesics(fx$dm, radius = fx$radius)
      e <- unname(g$edges)
      expect_identical(e[order(e[, 1], e[, 2]), , drop = FALSE], fx$edges, info = name)
      expect_identical(g$components, fx$components, info = name)
      expect_identical(g$membership, fx$membership, info = name)

      # Edge weights are the distances.
      dmat <- as.matrix(fx$dm)
      expect_identical(g$weights, dmat[g$edges], info = name)

      # Geodesic distances.
      if (!is.null(fx$geo)) {
        expect_true(rel_diff(g$geo, fx$geo) < 1e-12, info = name)
        expect_identical(g$max, max(g$geo), info = name)
      }

      # Embedding and eigenvalues.
      run <- function ()
        if (knn) isomap(fx$dm, n_axes = fx$k, n_neighbors = fx$n_neighbors, cpus = 1)
        else     isomap(fx$dm, n_axes = fx$k, radius = fx$radius, cpus = 1)

      if (split) { expect_warning(ord <- run(), 'connected components', info = name) }
      else       { ord <- run()                                                  }

      expect_identical(dim(ord), dim(fx$embedding), info = name)
      expect_true(coord_diff(ord, fx$embedding) < 1e-9, info = name)
      expect_true(rel_diff(attr(ord, 'eig'), fx$eig) < 1e-10, info = name)
      expect_identical(attr(ord, 'components'), fx$components, info = name)

      # Same signs as scikit-learn, whose largest |value| on each axis is
      # positive, except on axes whose largest |values| are tied, where
      # isomap() makes the first of them positive, as pcoa() does.
      emb  <- fx$embedding
      tied <- apply(emb, 2, function (v) sum(abs(v) >= max(abs(v)) * (1 - 1e-8)) > 1)
      same <- sign(colSums(vals(ord) * emb)) == 1
      expect_true(all(same[!tied]), info = name)
      expect_identical(ord_flip(vals(ord)), vals(ord), info = name)
    }
  }



  # Shared ordination result format. ====

  set.seed(1)
  m <- matrix(rpois(120 * 30, 2), 120, 30, dimnames = list(paste0('S', 1:120), NULL))
  d <- bray(m)
  x <- isomap(d, n_neighbors = 10, cpus = 1)

  expect_identical(class(x), c('ecodive_ord', 'matrix', 'array'))
  expect_true(is.double(x))
  expect_identical(
    names(attributes(no_display(x))),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace',
      'n_neighbors', 'radius', 'components', 'membership', 'class') )
  expect_identical(dimnames(x), list(attr(d, 'Labels'), c('Isomap1', 'Isomap2', 'Isomap3')))
  expect_identical(attr(x, 'method'), 'isomap')
  expect_identical(attr(x, 'n_neighbors'), 10L)
  expect_identical(attr(x, 'radius'), NA_real_)
  expect_identical(attr(x, 'components'), 1L)
  expect_identical(attr(x, 'membership'), setNames(rep(1L, 120), attr(d, 'Labels')))
  expect_equal(attr(x, 'prop'), attr(x, 'eig') / attr(x, 'trace'))
  expect_identical(ord_flip(vals(x)), vals(x))

  # Through dist_ord().
  expect_identical(dist_ord(d, 'isomap', n_neighbors = 10, cpus = 1), x)
  expect_identical(dist_ord(d, 'isometric mapping', n_neighbors = 10, seed = 1, cpus = 1), x)

  # A symmetric matrix gives the same result as a dist.
  expect_identical(no_display(isomap(as.matrix(d), n_neighbors = 10, cpus = 1)), no_display(x))

  # Radius neighbors.
  r <- isomap(d, radius = 0.6, cpus = 1)
  expect_identical(class(r), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(attr(r, 'n_neighbors'), NA_integer_)
  expect_identical(attr(r, 'radius'), 0.6)

  # The embedding is pcoa() of the geodesic distances.
  g  <- geodesics(d, n_neighbors = 10)
  pc <- pcoa(structure(g$geo, Size = 120L, class = 'dist'), n_axes = 3, cpus = 1)
  expect_identical(vals(pc), vals(x))
  expect_identical(attr(pc, 'eig'), attr(x, 'eig'))
  expect_identical(attr(pc, 'trace'), attr(x, 'trace'))



  # With every pair linked, Euclidean distances are their own geodesics,
  # and isomap() is pcoa(). ====

  set.seed(2)
  pts <- matrix(rnorm(130 * 4), 130, 4)
  d   <- dist(pts)
  g   <- geodesics(d, n_neighbors = 129)
  expect_identical(g$geo, as.vector(d))
  expect_identical(nrow(g$edges), 8385L)

  iso <- isomap(d, n_axes = 3, n_neighbors = 129, cpus = 1)
  pc  <- pcoa(d, n_axes = 3, cpus = 1)
  expect_identical(vals(iso), vals(pc))
  expect_identical(attr(iso, 'eig'),   attr(pc, 'eig'))
  expect_identical(attr(iso, 'prop'),  attr(pc, 'prop'))
  expect_identical(attr(iso, 'trace'), attr(pc, 'trace'))

  # So do an infinite radius, and a radius of 0, which leaves every sample
  # on its own until they are all joined pairwise.
  rad <- isomap(d, n_axes = 3, radius = Inf, cpus = 1)
  expect_identical(vals(rad), vals(iso))
  expect_warning(rad <- isomap(d, n_axes = 3, radius = 0, cpus = 1), '130 connected components')
  expect_identical(vals(rad), vals(iso))
  expect_identical(attr(rad, 'membership'), 1:130)

  # Many components. Each thread that scans for the closest pairs between
  # them has a table of 16 bytes per pair of components; extra threads are
  # only used while all tables fit in the byte limit, but one always is.
  set.seed(4)
  d <- dist(matrix(rnorm(150 * 2), 150, 2) + 20 * rep(1:15, 10))
  g <- geodesics(d, n_neighbors = 3)
  expect_identical(g$components, 15L)
  for (cpus in 2:max_cpus) for (bytes in c(0, 2.5 * 15 * 14 / 2 * 16))
    expect_identical(geodesics(d, n_neighbors = 3, cpus = cpus, join = c(bytes, 10, 65536)), g)

  # Joining is refused when its edges would be more than `ratio` times the
  # neighbor graph's own and more than `floor`. Here 105 joining edges.
  n_own <- nrow(g$edges) - 105
  ratio <- function (r) c(16 * 2^20, r / n_own, 0)
  expect_identical(geodesics(d, n_neighbors = 3, join = ratio(105.5)), g)
  expect_error(geodesics(d, n_neighbors = 3, join = ratio(104.5)), '15 connected components')
  expect_identical(geodesics(d, n_neighbors = 3, join = c(16 * 2^20, 0, 105)), g)
  expect_error(geodesics(d, n_neighbors = 3, join = c(16 * 2^20, 0, 104)), 'add 105 edges')

  # isomap()'s limits: 10 times, and a floor of 65,536 edges, which up to 362
  # samples on their own stay within (130 did, above).
  d363 <- dist(seq_len(363))
  expect_error(isomap(d363, radius = 0, cpus = 1), '363 connected components.*Increase `radius`')
  pairs <- dist(rep(seq_len(400) * 10, each = 2) + rep(0:1, 400))
  expect_error(isomap(pairs, n_neighbors = 1, cpus = 1), '400 connected components.*Increase `n_neighbors`')
  expect_identical(attr(suppressWarnings(isomap(pairs, n_neighbors = 2, cpus = 1)), 'components'), 1L)

  # Bray-Curtis isn't metric: there, geodesics can be shorter than `dm`.
  d <- bray(rbind(c(1, 0), c(1, 1), c(0, 1)))
  expect_equal(as.vector(d), c(1/3, 1, 1/3))
  expect_equal(geodesics(d, n_neighbors = 2)$geo, c(1/3, 2/3, 1/3))



  # A long ecological gradient: Isomap keeps its order; PCoA bends it. ====

  set.seed(3)
  pos    <- sort(runif(80, 0, 20))
  opt    <- seq(-2, 22, by = 1)
  counts <- matrix(rpois(80 * 25, 50 * exp(-outer(pos, opt, '-')^2 / 4.5)), 80, 25)
  counts[, 1] <- counts[, 1] + (rowSums(counts) == 0)
  d      <- bray(counts)
  expect_true(mean(as.vector(d) > 0.99) > 0.25) # saturated distances

  rho_iso  <- abs(cor(isomap(d, n_axes = 1, n_neighbors = 6, cpus = 1)[, 1], pos, method = 'spearman'))
  rho_pcoa <- abs(cor(pcoa(d, n_axes = 1, cpus = 1)[, 1], pos, method = 'spearman'))
  expect_true(rho_iso > 0.99)
  expect_true(rho_iso > rho_pcoa + 0.05)



  # Identical across runs, thread counts, and R's RNG state. ====

  # Duplicate samples (5 or 6 copies each of 20 profiles), tied distances,
  # disconnected graphs, and radius graphs. At least 100 samples each, so
  # that threads are used.
  cases <- list(
    list(d = bray(big_mtx),    args = list(n_neighbors = 3)),
    list(d = bray(big_mtx),    args = list(n_neighbors = 40)),
    list(d = bray(m),          args = list(n_neighbors = 4, n_axes = 3)),
    list(d = bray(m),          args = list(radius = 0.55)),
    list(d = dist(round(pts)), args = list(n_neighbors = 2)) )

  for (case in cases) {

    run <- function (cpus) suppressWarnings(do.call(isomap, c(list(case$d, cpus = cpus), case$args)))
    nbr <- case$args[names(case$args) %in% c('n_neighbors', 'radius')]

    x <- run(1)
    for (cpus in 2:max_cpus)
      expect_identical(run(cpus), x)

    g <- do.call(geodesics, c(list(case$d), nbr, cpus = 1))
    expect_identical(do.call(geodesics, c(list(case$d), nbr, cpus = max_cpus)), g)
  }

  # More threads than samples (isomap() caps them at n_cpus()).
  case <- cases[[1]]
  expect_identical(
    suppressWarnings(isomap(case$d, n_neighbors = 3, cpus = 109)),
    suppressWarnings(isomap(case$d, n_neighbors = 3, cpus = 1)) )

  # 5 or 6 identical copies of each of 20 samples, 3 neighbors each: each
  # copy links to the first 3 others, so each profile is one component.
  x <- suppressWarnings(isomap(bray(big_mtx), n_neighbors = 3, cpus = 1))
  expect_identical(attr(x, 'components'), 20L)
  expect_identical(unname(attr(x, 'membership')), rep(1:20, length.out = 104))

  # R's random number generator is neither used nor disturbed.
  set.seed(10); seed <- .Random.seed
  x <- isomap(bray(m), n_neighbors = 10, cpus = 1)
  expect_identical(.Random.seed, seed)
  set.seed(11)
  expect_identical(isomap(bray(m), n_neighbors = 10, cpus = 1), x)

  # Identical in a fresh R session, loading ecodive from the same library
  # as this session.
  if (at_home()) {
    lib     <- dirname(getNamespaceInfo('ecodive', 'path'))
    # Forward slashes: Windows' backslashes would be escapes in the script.
    f_in    <- gsub('\\', '/', tempfile(fileext = '.rds'), fixed = TRUE)
    f_out   <- gsub('\\', '/', tempfile(fileext = '.rds'), fixed = TRUE)
    saveRDS(bray(m), f_in)
    script  <- sprintf(paste(
      ".libPaths(c('%s', .libPaths()));",
      "ns <- loadNamespace('ecodive');",
      "stopifnot(normalizePath(getNamespaceInfo(ns, 'path')) == normalizePath('%s'));",
      "saveRDS(ns$isomap(readRDS('%s'), n_neighbors = 10, cpus = 2), '%s')" ),
      lib, getNamespaceInfo('ecodive', 'path'), f_in, f_out )
    rscript <- file.path(R.home('bin'), 'Rscript')
    status  <- system2(rscript, c('--vanilla', '-e', shQuote(script)), stdout = FALSE, stderr = FALSE)
    expect_identical(status, 0L)
    expect_identical(readRDS(f_out), x)
    unlink(c(f_in, f_out))
  }



  # Edge cases. ====

  # No samples, one sample.
  ord <- isomap(as.dist(matrix(numeric(0), 0, 0)))
  expect_identical(dim(ord), c(0L, 0L))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(attr(ord, 'components'), 0L)

  ord <- isomap(dist(matrix(1:3, 1, dimnames = list('A', NULL))))
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(rownames(ord), 'A')
  expect_identical(attr(ord, 'eig'), numeric(0))
  expect_identical(attr(ord, 'membership'), c(A = 1L))
  expect_identical(attr(ord, 'n_neighbors'), 0L)
  expect_identical(attr(isomap(dist(1), radius = 1), 'radius'), 1)

  # Two samples: the default n_neighbors is reduced quietly.
  ord <- isomap(dist(c(a = 1, b = 3)))
  expect_equal(vals(ord), matrix(c(1, -1), 2))
  expect_identical(dimnames(ord), list(c('a', 'b'), 'Isomap1'))
  expect_identical(attr(ord, 'n_neighbors'), 1L)

  # Three and four samples.
  d <- dist(matrix(c(0, 0, 3, 0, 0, 4), 3, byrow = TRUE))
  expect_identical(vals(isomap(d, n_neighbors = 2)), vals(pcoa(d, n_axes = 2)))

  d <- dist(c(0, 1, 10, 11))
  expect_warning(ord <- isomap(d, n_axes = 1, n_neighbors = 1), '2 connected components')
  expect_identical(unname(attr(ord, 'membership')), c(1L, 1L, 2L, 2L))
  expect_equal(abs(vals(ord)[, 1]), c(5.5, 4.5, 4.5, 5.5))

  # n_neighbors larger than n - 1.
  d <- bray(counts[1:10, ])
  expect_warning(ord <- isomap(d, n_neighbors = 10),  'more than the 9 other samples')
  expect_identical(ord, isomap(d, n_neighbors = 9))
  expect_warning(ord <- isomap(d, n_neighbors = 1e10), '1e\\+10 is more than')
  expect_identical(attr(ord, 'n_neighbors'), 9L)

  # k larger than n - 1.
  expect_identical(
    no_display(suppressWarnings(isomap(d, n_axes = 1e10, n_neighbors = 9))),
    no_display(suppressWarnings(isomap(d, n_axes = 9,    n_neighbors = 9))) )

  # All distances zero: no axes.
  expect_warning(ord <- isomap(as.dist(matrix(0, 5, 5))), 'only 0 of the first 3')
  expect_identical(dim(ord), c(5L, 0L))
  expect_identical(attr(ord, 'trace'), 0)
  expect_identical(attr(ord, 'prop'), numeric(0))

  # All distances equal: ties go to the earlier samples.
  d <- as.dist(matrix(1, 6, 6))
  g <- geodesics(d, n_neighbors = 2)
  expect_identical(
    unname(g$edges),
    matrix(c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,  2L, 3L, 4L, 5L, 6L, 3L, 4L, 5L, 6L), ncol = 2) )
  expect_identical(g$geo, as.vector(as.dist(outer(1:6, 1:6, function (i, j) ifelse(i <= 2 | j <= 2, 1, 2)))))
  expect_identical(isomap(d, n_neighbors = 2, cpus = 2), isomap(d, n_neighbors = 2, cpus = 1))

  # Integer storage.
  d <- as.dist(matrix(c(0L, 3L, 4L, 3L, 0L, 5L, 4L, 5L, 0L), 3))
  expect_identical(isomap(d, n_neighbors = 2), isomap(d / 1, n_neighbors = 2))



  # Invalid input. ====

  d <- bray(counts[1:10, ])

  d_na <- d; d_na[2] <- NA
  expect_error(isomap(d_na), 'NA, NaN, or infinite')
  d_na[2] <- NaN
  expect_error(isomap(d_na), 'NA, NaN, or infinite')
  d_na[2] <- Inf
  expect_error(isomap(d_na), 'NA, NaN, or infinite')
  d_na[2] <- -Inf
  expect_error(isomap(d_na), 'NA, NaN, or infinite')
  expect_error(isomap(d_na, radius = 1), 'NA, NaN, or infinite')
  d_na[2] <- -0.1
  expect_error(isomap(d_na), 'negative')
  expect_error(isomap(d_na, radius = 1), 'negative')

  # Distances too large to square are rejected outright.
  expect_error(isomap(d * 1e160), 'too large')

  # Each distance is at most 1e150, but a geodesic sums several of them,
  # and the geodesics are held to the same limit. In a chain of samples,
  # each `step` from the next and 0.99e150 from all others, with one neighbor
  # each, the geodesics run along the chain: up to 1.8e152 for 200 samples
  # 0.9e150 apart, and 1.2e150 for 3 samples 0.6e150 apart.
  chain <- function (n, step) {
    m <- matrix(0.99e150, n, n)  # under valgrind, R parses 1e150 2 ulp high
    m[abs(row(m) - col(m)) == 1] <- step
    diag(m) <- 0
    as.dist(m)
  }
  expect_error(isomap(chain(200, 0.9e150), n_neighbors = 1, cpus = 1), 'exceed 1e150')
  expect_error(isomap(chain(3,   0.6e150), n_neighbors = 1, cpus = 1), 'exceed 1e150')
  expect_error(isomap(chain(3,   0.6e150), radius = 0.6e150, cpus = 1), 'exceed 1e150')
  g <- geodesics(chain(3, 0.5e150), n_neighbors = 1)
  expect_identical(g$geo, c(0.5e150, 1e150, 0.5e150))
  expect_identical(g$max, 1e150)

  for (bad in list(0, -1, 1.5, NA, NaN, Inf, 'a', c(1, 2), NULL, TRUE))
    expect_error(isomap(d, n_axes = bad), '`n_axes` must be')
  for (bad in list(0, -1, 1.5, NA, NaN, Inf, 'a', c(1, 2), NULL, TRUE))
    expect_error(isomap(d, n_neighbors = bad), '`n_neighbors` must be')
  for (bad in list(-1, NA, NaN, 'a', c(1, 2), TRUE))
    expect_error(isomap(d, radius = bad), '`radius` must be')
  for (bad in list(0, -1, 1.5, NA, 'a', c(1, 2)))
    expect_error(isomap(d, cpus = bad), '`cpus`')

  expect_error(isomap(d, n_neighbors = 3, radius = 0.5), 'not both')

  # n_neighbors = NULL means "not given", as in scikit-learn.
  expect_identical(isomap(d, n_neighbors = NULL, radius = 0.5), isomap(d, radius = 0.5))

  # Huge whole numbers are reduced without stray warnings (such as R's
  # "probable complete loss of accuracy in modulus").
  warnings_of <- function (expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function (cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart('muffleWarning') })
    w
  }
  w <- warnings_of(ord <- isomap(d, n_axes = 1e300, n_neighbors = 9))
  expect_identical(w, warnings_of(isomap(d, n_axes = 9, n_neighbors = 9)))
  expect_identical(no_display(ord), no_display(suppressWarnings(isomap(d, n_axes = 9, n_neighbors = 9))))
  w <- warnings_of(ord <- isomap(d, n_neighbors = 1e300))
  expect_identical(length(w), 1L)
  expect_true(grepl('more than the 9 other samples', w[1]))
  expect_identical(attr(ord, 'n_neighbors'), 9L)

  expect_error(isomap(matrix(1:6, 2)),           'symmetric')
  expect_error(isomap(matrix(c(0, 1, 2, 0), 2)), 'symmetric')
  expect_error(isomap(matrix(letters[1:4], 2)),  'symmetric')
  expect_error(isomap(list(1, 2)),               'symmetric')
  expect_error(isomap(structure(1:2, class = 'dist', Size = 3L)), 'valid dist')

#})
