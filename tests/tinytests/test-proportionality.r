#test_that("proportionality", {

  # vlr(), rho_p() and phi_s(): scikit-bio 0.7.4's pairwise_vlr() and vlr()
  # (data-raw/expected/proportionality.py), base-R versions of the published
  # formulas (propr is GPL and could not be installed for fixtures),
  # identities, threading, input types and argument checks.

  ref <- readRDS(file.path('expected', 'proportionality.rds'))

  plain <- function (d) as.vector(unclass(d))
  maxabs <- function (a, b) max(abs(a - b))

  # The condensed vector (as `dist` stores it) of f(i, j) over features.
  condensed <- function (D, f) {
    unlist(lapply(seq_len(D - 1L), function (i)
      vapply((i + 1L):D, function (j) f(i, j), numeric(1)) ))
  }

  # The formulas, from the logs and a reference per sample.
  base_vlr <- function (L) condensed(ncol(L), function (i, j) var(L[, i] - L[, j]))
  base_rho <- function (z) condensed(ncol(z), function (i, j)
    1 - var(z[, i] - z[, j]) / (var(z[, i]) + var(z[, j])) )
  base_phi <- function (z) condensed(ncol(z), function (i, j)
    var(z[, i] - z[, j]) / var(z[, i] + z[, j]) )


  # -----------------------------------------------------------------------
  # scikit-bio
  # -----------------------------------------------------------------------

  for (nm in names(ref)) {
    r  <- ref[[nm]]
    x  <- r$counts
    pc <- if (r$pc > 0) r$pc else NULL

    if (r$pc > 0) expect_warning(vlr(x, cpus = 1L), 'Zeros detected', info = nm)
    v <- vlr(x, pseudocount = pc, cpus = 1L)
    expect_true(maxabs(plain(v), r$pairwise) < 1e-14, info = nm)
    expect_true(maxabs(plain(v), r$vlr)      < 1e-14, info = nm)
    if (nm != 'nearprop') {   # scikit-bio's covariance route cancels there
      expect_equal(plain(v), r$pairwise, tolerance = 1e-12, info = nm)
      expect_equal(plain(v), r$vlr,      tolerance = 1e-12, info = nm)
    }
  }

  # scikit-bio's own test values.
  expect_equal(plain(vlr(cbind(c(1, 1, 2), c(3, 5, 8)))), 0.0655828061998637, tolerance = 1e-14)
  expect_equal(sum(vlr(rbind(c(1, 1, 2), c(3, 5, 8), c(13, 21, 55)))), 0.2857382286903922, tolerance = 1e-14)


  # -----------------------------------------------------------------------
  # The formulas, every denominator
  # -----------------------------------------------------------------------

  for (nm in c('small', 'zeros', 'wide', 'prop', 'nearprop')) {
    r  <- ref[[nm]]
    x  <- r$counts
    L  <- log(x + r$pc)
    D  <- ncol(x)
    pc <- r$pc
    tol <- if (nm == 'nearprop') 1e-10 else 1e-12   # rounding in base R's var()

    expect_equal(plain(vlr(x, pseudocount = pc)), base_vlr(L), tolerance = 1e-12, info = nm)

    dens <- list(NULL, 1L, D, c(2L, 4L, 5L), seq_len(D))
    for (den in dens) {
      info <- paste(nm, paste(den, collapse = ' '))
      g    <- if (is.null(den)) rowMeans(L) else rowMeans(L[, den, drop = FALSE])
      z    <- L - g
      keep <- if (length(den) == 1L) -den else seq_len(D)

      rho <- rho_p(x, denominator = den, pseudocount = pc)
      phi <- phi_s(x, denominator = den, pseudocount = pc)
      expect_equal(plain(rho), base_rho(z[, keep]), tolerance = tol, info = info)
      expect_equal(plain(phi), base_phi(z[, keep]), tolerance = tol, info = info)

      # phi_s = (1 - rho_p) / (1 + rho_p); 1 - rho_p = vlr / (v_i + v_j)
      expect_equal(plain(phi), plain((1 - rho) / (1 + rho)), tolerance = 1e-10, info = info)
      vz <- apply(z[, keep], 2, var)
      vv <- condensed(length(vz), function (i, j) vz[i] + vz[j])
      vl <- plain(vlr(x[, keep], pseudocount = pc))
      expect_equal(1 - plain(rho), vl / vv, tolerance = 1e-10, info = info)

      # vlr is the same with any log-ratio transform.
      expect_equal(base_vlr(z[, keep]), base_vlr(L[, keep]), tolerance = 1e-10, info = info)
    }
  }

  # The clr is the default, and a reference set of every feature is the clr.
  x <- ref$small$counts
  expect_identical(plain(rho_p(x)), plain(rho_p(x, denominator = NULL)))
  expect_equal(plain(rho_p(x, denominator = 1:6)), plain(rho_p(x)), tolerance = 1e-14)


  # -----------------------------------------------------------------------
  # Exactly proportional features; scale invariance
  # -----------------------------------------------------------------------

  set.seed(1)
  b <- rpois(15, 30) + 1
  x <- cbind(a = b, b2 = 2 * b, b8 = 8 * b, c = rpois(15, 30) + 1, d = rpois(15, 30) + 1)
  m <- function (d) as.matrix(d)
  expect_equal(m(vlr(x))['a', 'b2'], 0, tolerance = 1e-14)
  expect_equal(m(vlr(x))['b2', 'b8'], 0, tolerance = 1e-14)
  expect_true(abs(m(rho_p(x))['a', 'b8'] - 1) < 1e-14)
  expect_true(abs(m(phi_s(x))['a', 'b8']) < 1e-14)
  expect_true(abs(m(rho_p(x, denominator = 'c'))['a', 'b2'] - 1) < 1e-14)

  y <- x * c(10, rep(1, 13), 0.001)     # samples 1 and 15 rescaled
  for (f in c('vlr', 'rho_p', 'phi_s'))
    expect_equal(plain(get(f)(y)), plain(get(f)(x)), tolerance = 1e-12, info = f)
  expect_equal(plain(rho_p(y, denominator = 4)), plain(rho_p(x, denominator = 4)), tolerance = 1e-12)

  # Reciprocal features: rho_p -1, phi_s infinite (or huge, from rounding).
  r2 <- cbind(p = c(1, 2, 4, 8), q = c(8, 4, 2, 1))
  expect_equal(plain(rho_p(r2, denominator = NULL)), -1, tolerance = 1e-14)
  expect_true(plain(phi_s(r2)) > 1e14)
  # (so are the two features of a two-feature reference: z_2 = -z_4)
  ph <- as.matrix(phi_s(ref$small$counts, denominator = c(2, 4)))
  expect_true(ph[2, 4] > 1e10)
  expect_equal(as.matrix(rho_p(ref$small$counts, denominator = c(2, 4)))[2, 4], -1, tolerance = 1e-12)

  # Two features constant relative to the reference: rho_p is NaN.
  k <- cbind(a = c(1, 2, 3), b = c(1, 2, 3), c = c(2, 4, 6), d = c(5, 1, 2))
  expect_true(is.nan(m(rho_p(k, denominator = c('a', 'b')))['a', 'b']))


  # -----------------------------------------------------------------------
  # Threads, input types, orientation
  # -----------------------------------------------------------------------

  set.seed(2)
  big <- matrix(rpois(150 * 300, 0.8), 150, 300,
    dimnames = list(paste0('s', 1:150), paste0('f', 1:300)) )
  big[, 1] <- big[, 1] + 1L
  sp <- Matrix::Matrix(big, sparse = TRUE)
  sl <- slam::as.simple_triplet_matrix(big)
  for (f in c('vlr', 'rho_p', 'phi_s')) {
    fn <- get(f)
    d1 <- fn(big, pseudocount = 0.5, cpus = 1L)
    for (cpus in c(2L, 3L, 4L, 8L))
      expect_identical(fn(big, pseudocount = 0.5, cpus = cpus), d1, info = paste(f, cpus))
    expect_identical(fn(sp, pseudocount = 0.5, cpus = 4L), d1, info = f)
    expect_identical(fn(sl, pseudocount = 0.5, cpus = 2L), d1, info = f)
    expect_identical(fn(t(big), margin = 2L, pseudocount = 0.5, cpus = 2L), d1, info = f)
    expect_identical(fn(Matrix::t(sp), margin = 2L, pseudocount = 0.5, cpus = 4L), d1, info = f)
  }
  d1 <- rho_p(big, denominator = 1, pseudocount = 0.5, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(rho_p(sp, denominator = 'f1', pseudocount = 0.5, cpus = cpus), d1)

  # Integer, logical and data-frame tables.
  x <- ref$small$counts
  expect_identical(vlr(x), vlr(x * 1.0))
  expect_identical(plain(rho_p(as.data.frame(x))), plain(rho_p(x)))
  expect_silent(vlr(matrix(TRUE, 3, 3)))
  expect_identical(plain(vlr(matrix(TRUE, 3, 3))), c(0, 0, 0))


  # -----------------------------------------------------------------------
  # Zeros
  # -----------------------------------------------------------------------

  x <- ref$zeros$counts
  expect_error(vlr(x, pseudocount = 0), 'has zeros')
  expect_error(rho_p(sp, pseudocount = 0), 'has zeros')
  expect_error(vlr(x, pseudocount = -1), 'non-negative')
  expect_error(vlr(x, pseudocount = c(1, 2)), 'non-negative')
  expect_warning(phi_s(x), 'Zeros detected')

  # replace_zeros() output goes straight in, with no pseudocount.
  z <- replace_zeros(x, 'mult')
  expect_silent(r <- rho_p(z))
  expect_identical(r, rho_p(z, pseudocount = 0))
  expect_equal(plain(r), base_rho(log(z) - rowMeans(log(z))), tolerance = 1e-12)
  expect_equal(plain(vlr(z)), base_vlr(log(z)), tolerance = 1e-12)

  # propr's zero handling (zeros alone replaced by 1).
  x1 <- x; x1[x1 == 0] <- 1
  expect_silent(r <- rho_p(x1))
  expect_equal(plain(r), base_rho(log(x1) - rowMeans(log(x1))), tolerance = 1e-12)


  # -----------------------------------------------------------------------
  # Results' structure
  # -----------------------------------------------------------------------

  x <- ref$small$counts
  colnames(x) <- paste0('F', 1:6)
  r <- rho_p(x, denominator = 'F3')
  expect_identical(attr(r, 'Labels'), paste0('F', c(1:2, 4:6)))
  expect_identical(attr(r, 'Size'), 5L)
  expect_identical(attr(r, 'method'), 'rho_p')
  expect_inherits(r, 'dist')
  expect_identical(rho_p(x, denominator = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)), rho_p(x, denominator = c('F1', 'F3')))
  expect_identical(phi_s(x, denominator = c('F4', 'F2')), phi_s(x, denominator = c(2, 4)))
  expect_identical(attr(vlr(x), 'method'), 'vlr')
  expect_identical(attr(phi_s(x), 'method'), 'phi_s')
  expect_identical(attr(vlr(x), 'Labels'), colnames(x))

  # Unnamed features: no labels, except for the positions kept by an alr.
  x <- unname(ref$small$counts)
  expect_null(attr(vlr(x), 'Labels'))
  expect_identical(attr(rho_p(x, denominator = 2), 'Labels'), as.character(c(1, 3:6)))
  expect_identical(rho_p(x, denominator = 2), rho_p(x, denominator = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)))
  expect_identical(rho_p(x, denominator = c(5, 1)), rho_p(x, denominator = c(1, 5)))

  # Usable as a dissimilarity.
  x <- ref$wide$counts
  expect_identical(length(stats::hclust(vlr(x))$order), 60L)
  expect_identical(nrow(pcoa(phi_s(x), n_axes = 2)), 60L)
  expect_identical(dim(as.matrix(rho_p(x))), c(60L, 60L))


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  x <- ref$small$counts
  colnames(x) <- paste0('F', 1:6)
  expect_error(vlr(x[1, , drop = FALSE]), 'at least 2 samples')
  expect_error(rho_p(x[, 1, drop = FALSE]), 'at least 2 features')
  expect_error(rho_p(x[, 1:2], denominator = 1), 'at least 2 features besides')
  expect_silent(rho_p(x[, 1:3], denominator = 1))
  expect_error(rho_p(x, denominator = 'F9'), 'not in `counts`')
  expect_error(phi_s(x, denominator = 7), 'outside 1 to 6')
  expect_error(phi_s(x, denominator = c(1, 1)), 'more than once')
  expect_error(phi_s(x, denominator = NA), 'NA')
  expect_error(phi_s(x, denominator = character(0)), 'no features')
  expect_error(rho_p(unname(x), denominator = 'F1'), 'no feature names')
  expect_error(vlr(x, margin = 3), 'margin')
  expect_error(vlr(x, cpus = 0), 'cpus')
  xn <- x; xn[1, 1] <- NA
  expect_error(vlr(xn), 'NA')
  expect_error(vlr(-x), 'negative')
  expect_error(vlr(x * 0), 'only zeros')

#})
