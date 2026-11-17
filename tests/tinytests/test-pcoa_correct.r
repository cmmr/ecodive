#test_that("pcoa corrections", {

  # Largest absolute difference between two coordinate matrices, relative to
  # the largest coordinate, after matching each axis's sign.
  coord_diff <- function (x, y) {
    x <- matrix(x, nrow(x))
    y <- matrix(y, nrow(y))[, seq_len(ncol(x)), drop = FALSE]
    y <- sweep(y, 2, sign(colSums(x * y)), '*')
    max(abs(x - y)) / max(abs(y))
  }

  # Each correction computed in full from its formula, with base R's eigen().
  dense <- function (d, correction) {
    D  <- as.matrix(d)
    n  <- nrow(D)
    J  <- diag(n) - 1 / n
    B  <- -0.5 * J %*% (D^2) %*% J
    B2 <- -0.5 * J %*% D %*% J
    if (correction == 'lingoes') {
      const <- -min(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
      D2    <- D^2 + 2 * const
    } else {
      M     <- rbind(cbind(matrix(0, n, n), 2 * B), cbind(-diag(n), -4 * B2))
      ev    <- eigen(M, only.values = TRUE)$values
      const <- max(Re(ev[abs(Im(ev)) < 1e-8 * max(Mod(ev))]))
      D2    <- (D + const)^2
    }
    diag(D2) <- 0
    Bc <- -0.5 * J %*% D2 %*% J
    e  <- eigen(Bc, symmetric = TRUE)
    list(const = const, eig = e$values, vec = e$vectors, trace = sum(diag(Bc)))
  }

  # The smallest eigenvalue of the corrected matrix, built from pcoa()'s own
  # constant, relative to the largest.
  corrected_min <- function (d, o) {
    D <- as.matrix(d)
    D <- if (attr(o, 'correction') == 'lingoes') D^2 + 2 * attr(o, 'constant') else (D + attr(o, 'constant'))^2
    diag(D) <- 0
    J <- diag(nrow(D)) - 1 / nrow(D)
    e <- eigen(-0.5 * J %*% D %*% J, symmetric = TRUE, only.values = TRUE)$values
    min(e) / max(e)
  }

  # Without the correction's attributes.
  strip <- function (x) {
    attr(x, 'correction') <- NULL
    attr(x, 'constant')   <- NULL
    x
  }

  set.seed(1)
  rand_counts <- function (n, p = 40) {
    m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
    m[, 1] <- m[, 1] + 1
    dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
    m
  }

  cpus <- if (at_home()) 4L else 2L

  m60  <- rand_counts(60)
  d100 <- bray(rand_counts(100, 30))
  ds   <- list(
    bray    = bray(m60),                              # many negative eigenvalues
    jaccard = jaccard(rand_counts(30)),
    ties    = jaccard(rand_counts(80, 6) > 0),        # heavily tied distances
    unifrac = unweighted_unifrac(ecodive::ex_counts, tree = ecodive::ex_tree),
    tri     = as.dist(matrix(c(0, 1, 3, 1, 0, 1, 3, 1, 0), 3)) ) # breaks the triangle inequality



  # Against the formulas. ====

  for (nm in names(ds)) for (corr in c('lingoes', 'cailliez')) {

    d   <- ds[[nm]]
    n   <- attr(d, 'Size')
    k   <- min(4L, n - 1L)
    ord <- suppressWarnings(pcoa(d, k = k, correction = corr, cpus = 1))
    ref <- dense(d, corr)
    kk  <- ncol(ord)
    inf <- paste(nm, corr)

    # pcoa()'s constant makes the distances Euclidean: the corrected
    # matrix's smallest eigenvalue is 0.
    expect_true(attr(ord, 'constant') > 0, info = inf)
    expect_true(abs(corrected_min(d, ord)) < 1e-10, info = inf)

    # The reference comes from the platform's BLAS/LAPACK, so it is used only
    # if it passes the same test: no negative eigenvalue, and the smallest
    # non-trivial one 0. (On R-hub, OpenBLAS 0.3.20 has given a wrong
    # Cailliez constant on some runners' CPUs; pcoa() itself uses neither.)
    ref_ok <- ref$const > 0 &&
      min(ref$eig) > -1e-12 * ref$eig[[1]] &&
      abs(sort(ref$eig)[[2]]) < 1e-10 * ref$eig[[1]]
    if (ref_ok) {
      expect_equal(attr(ord, 'constant'), ref$const,               tolerance = 1e-10, info = inf)
      expect_equal(attr(ord, 'eig'),      ref$eig[1:kk],           tolerance = 1e-10, info = inf)
      expect_equal(attr(ord, 'trace'),    ref$trace,               tolerance = 1e-10, info = inf)
      expect_equal(attr(ord, 'prop'),     ref$eig[1:kk] / ref$trace, tolerance = 1e-10, info = inf)
      expect_true(coord_diff(ord, ref$vec %*% diag(sqrt(pmax(ref$eig, 0)))) < 1e-8, info = inf)
    }

    expect_identical(attr(ord, 'correction'), corr)
    expect_identical(
      names(attributes(ord)),
      c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'correction', 'constant', 'class') )
  }

  # The same constant as cmdscale(add = TRUE), base R's Cailliez correction,
  # where cmdscale()'s constant (from the platform's LAPACK) makes the
  # distances Euclidean, as above.
  for (d in list(ds$bray, ds$ties, d100, ds$jaccard)) {
    oC <- pcoa(d, correction = 'cailliez', cpus = cpus)
    cm <- cmdscale(d, k = 3, add = TRUE, eig = TRUE)
    cm_ok <- abs(corrected_min(d, structure(0, correction = 'cailliez', constant = cm$ac))) < 1e-10
    if (cm_ok) {
      expect_equal(attr(oC, 'constant'), cm$ac, tolerance = 1e-10)
      expect_true(coord_diff(oC, cm$points) < 1e-8)
    }
  }


  # Lingoes keeps the eigenvectors: same axes as no correction, each
  # eigenvalue plus c1, and trace plus (n - 1) c1.
  d  <- ds$bray
  o0 <- pcoa(d, k = 5, cpus = 1)
  o1 <- pcoa(d, k = 5, correction = 'lingoes', cpus = 1)
  c1 <- attr(o1, 'constant')
  expect_equal(attr(o1, 'eig'), attr(o0, 'eig') + c1, tolerance = 1e-12)
  expect_equal(attr(o1, 'trace'), attr(o0, 'trace') + 59 * c1, tolerance = 1e-12)
  expect_equal(
    sweep(unclass(o1), 2, sqrt(attr(o1, 'eig')), '/'),
    sweep(unclass(o0), 2, sqrt(attr(o0, 'eig')), '/'),
    tolerance = 1e-12, check.attributes = FALSE )

  # Axes whose eigenvalue was negative are kept, but for the most negative,
  # which becomes 0 and, with k = n - 1, is dropped as expected: no warning.
  # With fewer axes than the correction can make positive, there is one.
  d <- ds$jaccard
  expect_warning(o0 <- pcoa(d, k = 29), 'eigenvalues are > 0')
  expect_silent(o1 <- pcoa(d, k = 29, correction = 'lingoes'))
  expect_silent(o2 <- pcoa(d, k = 29, correction = 'cailliez'))
  expect_true(ncol(o0) < 28)
  expect_identical(c(ncol(o1), ncol(o2)), c(28L, 28L))
  expect_equal(attr(o1, 'eig'), dense(d, 'lingoes')$eig[1:28],  tolerance = 1e-10)
  expect_equal(attr(o2, 'eig'), dense(d, 'cailliez')$eig[1:28], tolerance = 1e-10)
  pts <- dist(cbind(x = sin(1:40), y = cos(1:40 / 2) * 3))
  expect_warning(pcoa(pts, k = 39, correction = 'lingoes'), 'only 2 of the first 39')

  # Through ordination().
  expect_identical(ordination(ds$bray, 'pcoa', correction = 'lingoes'), pcoa(ds$bray, correction = 'lingoes'))



  # The corrected matrix has no negative eigenvalue. ====

  # Samples on a circular gradient: an exactly symmetric design, on which
  # the Cailliez iteration can settle on a root below the largest, which
  # the check afterwards catches.
  circular <- function (n, s, w) {
    th <- 2 * pi * (0:(n - 1)) / n
    mu <- 2 * pi * (0:(s - 1)) / s
    m  <- outer(th, mu, function (a, b) round(100 * exp((cos(a - b) - 1) / w^2)))
    dimnames(m) <- list(paste0('S', 1:n), paste0('F', 1:s))
    jaccard(m)
  }

  d <- circular(48, 10, 0.2)
  o <- pcoa(d, correction = 'cailliez', cpus = cpus)
  expect_equal(attr(o, 'constant'), cmdscale(d, add = TRUE)$ac, tolerance = 1e-10)
  expect_true(corrected_min(d, o) > -1e-10)

  for (d in c(ds, list(d100)))
    for (corr in c('lingoes', 'cailliez')) {
      o <- pcoa(d, correction = corr, cpus = cpus)
      expect_true(corrected_min(d, o) > -1e-10)
    }



  # No negative eigenvalues: no correction, same results. ====

  for (d in list(euclidean(m60), dist(c(a = 1, b = 3)), as.dist(matrix(1, 50, 50)),
                 dist(cbind(sin(1:40), cos(1:40))))) {
    o0 <- suppressWarnings(pcoa(d, cpus = 1))
    expect_identical(attr(o0, 'correction'), 'none')
    expect_identical(attr(o0, 'constant'), 0)
    for (corr in c('lingoes', 'cailliez')) {
      o <- suppressWarnings(pcoa(d, correction = corr, cpus = 1))
      expect_identical(attr(o, 'constant'), 0)
      expect_identical(attr(o, 'correction'), corr)
      expect_identical(strip(o), strip(o0))
    }
  }

  # One sample, and all distances zero.
  for (corr in c('lingoes', 'cailliez')) {
    o <- pcoa(dist(matrix(1:3, 1, dimnames = list('A', NULL))), correction = corr)
    expect_identical(dim(o), c(1L, 0L))
    expect_identical(attr(o, 'constant'), 0)
    expect_identical(attr(o, 'correction'), corr)
    expect_warning(o <- pcoa(as.dist(matrix(0, 5, 5)), correction = corr), 'only 0')
    expect_identical(attr(o, 'trace'), 0)
    expect_identical(attr(o, 'constant'), 0)
  }



  # Identical across runs, threads, and R's random number state. ====

  d <- d100
  for (corr in c('none', 'lingoes', 'cailliez')) {
    set.seed(2)
    seed <- .Random.seed
    x    <- pcoa(d, correction = corr, cpus = 1)
    expect_identical(.Random.seed, seed)
    expect_identical(pcoa(d, correction = corr, cpus = 1), x)
    expect_identical(pcoa(d, correction = corr, cpus = 2), x)
    if (at_home()) expect_identical(pcoa(d, correction = corr, cpus = 4), x)
    set.seed(3)
    expect_identical(pcoa(d, correction = corr, cpus = cpus), x)
    expect_true(all(apply(x, 2, function (v) v[which.max(abs(v))] > 0)))
  }



  # Any scale of distances: a power of two apart, exactly the same,
  # scaled. ====

  d <- ds$bray
  for (corr in c('none', 'lingoes', 'cailliez')) {
    r <- pcoa(d, correction = corr, cpus = 1)
    for (e in c(-100, 100)) {
      o <- pcoa(d * 2^e, correction = corr, cpus = 1)
      expect_identical(as.vector(o) * 2^-e,     as.vector(r), info = e)
      expect_identical(attr(o, 'eig') * 2^-(2 * e),   attr(r, 'eig'), info = e)
      expect_identical(attr(o, 'trace') * 2^-(2 * e), attr(r, 'trace'), info = e)
      expect_identical(attr(o, 'prop'), attr(r, 'prop'), info = e)
      expect_identical(
        attr(o, 'constant') * 2^-(if (corr == 'cailliez') e else 2 * e),
        attr(r, 'constant'), info = e )
    }
  }



  # Invalid input. ====

  d <- ds$bray
  expect_error(pcoa(d, correction = 'xyz'),                    '`correction` must be')
  expect_error(pcoa(d, correction = NA),                       '`correction` must be')
  expect_error(pcoa(d, correction = NA_character_),            '`correction` must be')
  expect_error(pcoa(d, correction = 1),                        '`correction` must be')
  expect_error(pcoa(d, correction = NULL),                     '`correction` must be')
  expect_error(pcoa(d, correction = ''),                       '`correction` must be')
  expect_error(pcoa(d, correction = c('lingoes', 'cailliez')), '`correction` must be')
  expect_identical(pcoa(d, correction = 'ling'), pcoa(d, correction = 'lingoes'))

  d_neg <- d; d_neg[2] <- -0.1
  expect_error(pcoa(d_neg, correction = 'lingoes'),  'negative distances')
  expect_error(pcoa(d_neg, correction = 'cailliez'), 'negative distances')

  d_na <- d; d_na[2] <- NA
  expect_error(pcoa(d_na, correction = 'cailliez'), 'NA')
  d_na[2] <- Inf
  expect_error(pcoa(d_na, correction = 'lingoes'), 'infinite')



  # Against ape and vegan, and hard symmetric cases: stored results and
  # inputs, so at home only. ====

  if (at_home()) {

    # Expected results from data-raw/expected/pcoa_correct.r.
    ex <- readRDS('expected/pcoa_correct.rds')

    for (cs in ex$cases) {
      for (corr in c('lingoes', 'cailliez')) {

        o   <- pcoa(cs$dist, correction = corr, cpus = 2)
        ref <- cs[[corr]]
        inf <- paste(cs$name, corr)

        # ape::pcoa(): Corr_eig and Rel_corr_eig.
        expect_equal(attr(o, 'eig'),  ref$ape$eig,     tolerance = 1e-10, info = inf)
        expect_equal(attr(o, 'prop'), ref$ape$rel_eig, tolerance = 1e-10, info = inf)
        expect_true(coord_diff(o, ref$ape$vectors) < 1e-8, info = inf)

        # vegan::wcmdscale(add = ) and its "Proportion Explained".
        expect_equal(attr(o, 'constant'), ref$vegan$constant, tolerance = 1e-10, info = inf)
        expect_equal(attr(o, 'eig'),      ref$vegan$eig,      tolerance = 1e-10, info = inf)
        expect_equal(attr(o, 'prop'),     ref$vegan$prop,     tolerance = 1e-10, info = inf)
        expect_true(coord_diff(o, ref$vegan$points) < 1e-8, info = inf)
      }
    }

    # "Ring" distances, g points on a circle whose distances depend only on
    # their separation, through a profile chosen (data-raw/expected/
    # cailliez_rings.r) to put the mode that sets the Cailliez constant far
    # up B's spectrum; and larger circular gradients.
    ring <- function (h) {
      g <- 2L * length(h)
      K <- abs(outer(0:(g - 1), 0:(g - 1), '-'))
      K <- pmin(K, g - K)
      D <- matrix(0, g, g)
      D[K > 0] <- h[K[K > 0]]
      as.dist(D)
    }

    hard <- list(
      circular(120, 15, 0.2), circular(200, 20, 0.15),
      ring(ex$rings[['120']]), ring(ex$rings[['200']]) )

    for (i in seq_along(hard)) {
      d <- hard[[i]]
      o <- pcoa(d, correction = 'cailliez', cpus = cpus)
      expect_equal(attr(o, 'constant'), cmdscale(d, add = TRUE)$ac, tolerance = 1e-10, info = i)
      expect_true(corrected_min(d, o) > -1e-10, info = i)
    }

    # A Cailliez iteration long enough to restart.
    set.seed(5)
    d <- jaccard(rand_counts(150))
    o <- pcoa(d, correction = 'cailliez', cpus = cpus)
    expect_equal(attr(o, 'constant'), cmdscale(d, add = TRUE)$ac, tolerance = 1e-10)
  }

#})
