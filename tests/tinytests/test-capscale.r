#test_that("capscale", {

  # Inputs and reference results, written by data-raw/expected/capscale.r.
  ref <- readRDS(file.path('expected', 'capscale.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case_formula <- function (cs) {
    rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
    stats::as.formula(paste('~', paste(rhs, collapse = ' + ')))
  }

  fit_case <- function (cs, cpus = 1, fn = capscale) {
    suppressMessages(fn(cs$dm, data = cs$data, formula = case_formula(cs), n_axes = cs$k, cpus = cpus))
  }

  case <- function (name) ref$cases[[match(name, vapply(ref$cases, `[[`, '', 'name'))]]



  # Same as vegan::capscale(), up to the sign of each axis: stored results,
  # so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      exp <- cs$expected
      k   <- ncol(exp$sites)
      lbl <- cs$name

      if (is.null(cs$warn)) {
        ord <- fit_case(cs)
      } else {
        expect_warning(ord <- fit_case(cs), cs$warn, info = lbl)
      }

      expect_identical(dim(ord), c(attr(cs$dm, 'Size'), k), info = lbl)
      expect_identical(dimnames(ord), dimnames(exp$sites), info = lbl)

      s <- sign(colSums(unclass(ord) * exp$sites))
      expect_equal(sweep(unclass(ord), 2, s, '*'), exp$sites, check.attributes = FALSE, tolerance = 1e-10, info = lbl)

      expect_equal(attr(ord, 'eig'),   exp$eig,   tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'),  exp$prop,  tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'trace'), exp$trace, tolerance = 1e-10, info = lbl)
      expect_equal(
        attr(ord, 'inertia'),
        c(conditional = exp$conditional, constrained = exp$constrained,
          unconstrained = exp$unconstrained, imaginary = exp$imaginary),
        tolerance = 1e-10, info = lbl )
      expect_equal(attr(ord, 'constrained_eig'), exp$constrained_eig, tolerance = 1e-10, info = lbl)

      if (length(exp$r_squared) == 0) {
        expect_identical(attr(ord, 'r_squared'),     0, info = lbl)
        expect_identical(attr(ord, 'adj_r_squared'), 0, info = lbl)
      } else {
        expect_equal(attr(ord, 'r_squared'),     exp$r_squared,     tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'adj_r_squared'), exp$adj_r_squared, tolerance = 1e-10, info = lbl)
      }

      bp <- attr(ord, 'biplot')
      if (is.null(exp$biplot)) {
        expect_identical(nrow(bp), 0L, info = lbl)
      } else {
        expect_identical(dimnames(bp), dimnames(exp$biplot), info = lbl)
        expect_equal(sweep(bp, 2, s, '*'), exp$biplot, tolerance = 1e-10, info = lbl)
      }

      cn <- attr(ord, 'centroids')
      if (is.null(exp$centroids)) {
        expect_null(cn, info = lbl)
      } else {
        expect_identical(dimnames(cn), dimnames(exp$centroids), info = lbl)
        expect_equal(sweep(cn, 2, s, '*'), exp$centroids, tolerance = 1e-10, info = lbl)
      }
    }
  }



  # Result format. ====

  cs   <- case('bray, factor + continuous')
  d    <- cs$dm
  meta <- cs$data
  ord  <- fit_case(cs)

  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(
    names(attributes(no_display(ord))),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'inertia',
      'constrained_eig', 'r_squared', 'adj_r_squared', 'biplot', 'centroids', 'class') )
  expect_identical(attr(ord, 'method'), 'capscale')
  expect_identical(colnames(ord), c('CAP1', 'CAP2', 'CAP3'))
  expect_identical(names(attr(ord, 'inertia')), c('conditional', 'constrained', 'unconstrained', 'imaginary'))

  # Bray-Curtis distances have negative eigenvalues: imaginary inertia.
  ine <- attr(ord, 'inertia')
  expect_true(ine[['imaginary']] < 0)
  expect_equal(attr(ord, 'trace'), sum(ine))
  real <- attr(ord, 'trace') - ine[['imaginary']]
  expect_equal(attr(ord, 'r_squared'), ine[['constrained']] / real)
  expect_equal(attr(ord, 'adj_r_squared'), 1 - (1 - attr(ord, 'r_squared')) * 19 / 15)
  expect_equal(attr(ord, 'prop'), attr(ord, 'eig') / real)

  # The total inertia is dbrda()'s (and pcoa()'s, as no distance is 4 or more).
  db <- fit_case(cs, fn = dbrda)
  expect_equal(attr(ord, 'trace'), attr(db, 'trace'))
  expect_equal(attr(ord, 'trace'), attr(pcoa(d, cpus = 1), 'trace'))
  expect_true(max(abs(attr(ord, 'eig') - attr(db, 'eig'))) > 1e-3)



  # The same as the dense formula: RDA of G's positive part. ====

  dense_capscale <- function (d, x, z = NULL) {
    n <- attr(d, 'Size')
    A <- as.matrix(d) ^ 2
    if (max(d) >= 4) A <- A / (n - 1)
    J <- diag(n) - 1 / n
    G <- -0.5 * J %*% A %*% J
    e <- eigen(G, symmetric = TRUE)
    tol <- max(sqrt(.Machine$double.eps), sqrt(.Machine$double.eps) * e$values[[1]])
    pos <- e$values > tol
    Gp  <- e$vectors[, pos] %*% (e$values[pos] * t(e$vectors[, pos]))
    cen <- function (m) sweep(m, 2, colMeans(m))
    Qz  <- if (is.null(z)) matrix(0, n, 0) else qr.Q(qr(cen(z)))
    xr  <- cen(x)
    xr  <- xr - Qz %*% crossprod(Qz, xr)
    Qx  <- qr.Q(qr(xr))
    ec  <- eigen(crossprod(Qx, Gp %*% Qx), symmetric = TRUE)
    u   <- Qx %*% ec$vectors
    wa  <- (diag(n) - tcrossprod(Qz)) %*% Gp %*% u %*% diag(1 / ec$values, length(ec$values))
    list(
      sites = wa * sqrt(sqrt((n - 1) * sum(diag(G)))),
      eig   = ec$values,
      trace = sum(diag(G)),
      imag  = sum(e$values[e$values < -tol]),
      cond  = sum(diag(crossprod(Qz, Gp %*% Qz))) )
  }

  check_dense <- function (ord, ref, lbl) {
    k <- sum(startsWith(colnames(ord), 'CAP'))  # constrained axes
    o <- unclass(ord)[, seq_len(k), drop = FALSE]
    r <- ref$sites[, seq_len(k), drop = FALSE]
    s <- sign(colSums(o * r))
    expect_equal(sweep(o, 2, s, '*'), r, check.attributes = FALSE, tolerance = 1e-8, info = lbl)
    expect_equal(attr(ord, 'eig')[seq_len(k)], ref$eig[seq_len(k)], tolerance = 1e-8, info = lbl)
    expect_equal(attr(ord, 'trace'), ref$trace,           tolerance = 1e-8, info = lbl)
    expect_equal(attr(ord, 'inertia')[['imaginary']],   ref$imag, tolerance = 1e-8, info = lbl)
    expect_equal(attr(ord, 'inertia')[['conditional']], ref$cond, tolerance = 1e-8, info = lbl)
  }

  x <- stats::model.matrix(~ Management + A1, meta)[, -1]
  check_dense(ord, dense_capscale(d, x), 'bray')

  z  <- stats::model.matrix(~ A1, meta)[, -1, drop = FALSE]
  xp <- stats::model.matrix(~ Management, meta)[, -1]
  op <- capscale(d, data = meta, formula = ~ Management + Condition(A1), cpus = 1)
  check_dense(op, dense_capscale(d, xp, z), 'bray, partial')

  # Distances of 4 or more are rescaled, as in dbrda().
  dm <- d * 10
  check_dense(capscale(dm, 'Management', meta, cpus = 1), dense_capscale(dm, xp), 'rescaled')

  # A few negative eigenvalues among many positive ones (eigenvectors of the
  # few).
  cs <- case('few negative, partial')
  fz <- stats::model.matrix(~ z, cs$data)[, -1, drop = FALSE]
  fx <- stats::model.matrix(~ g + x, cs$data)[, -1]
  check_dense(fit_case(cs), dense_capscale(cs$dm, fx, fz), 'few negative')



  # Euclidean distances: the same as dbrda(), apart from the names. ====

  for (cs in ref$cases[vapply(ref$cases, `[[`, NA, 'euclidean')]) {
    oc <- fit_case(cs)
    od <- fit_case(cs, fn = dbrda)
    expect_identical(sub('^CAP', 'dbRDA', colnames(oc)), colnames(od), info = cs$name)
    expect_equal(unclass(oc)[, ], unclass(od)[, ], check.attributes = FALSE, tolerance = 1e-10, info = cs$name)
    for (a in c('eig', 'prop', 'trace', 'constrained_eig', 'r_squared', 'adj_r_squared'))
      expect_equal(attr(oc, a), attr(od, a), tolerance = 1e-10, info = paste(cs$name, a))
    expect_identical(attr(oc, 'inertia')[['imaginary']], 0, info = cs$name)
    expect_equal(attr(oc, 'inertia')[1:3], attr(od, 'inertia'), tolerance = 1e-10, info = cs$name)
  }



  # Ways to give the model, as dbrda() and permanova(). ====

  cs   <- case('bray, factor only')
  d    <- cs$dm
  meta <- cs$data
  one  <- capscale(d, data = meta, formula = ~ Management, cpus = 1)

  dmx <- d
  expect_identical(no_display(capscale(d, 'Management', meta, cpus = 1)), no_display(one))
  expect_identical(no_display(capscale(dmx ~ Management, meta, cpus = 1)), no_display(one))
  expect_identical(no_display(capscale(as.matrix(d), 'Management', meta, cpus = 1)), no_display(one))
  expect_identical(no_display(capscale(d, 'Management', meta[rev(seq_len(nrow(meta))), ], cpus = 1)), no_display(one))
  expect_equal(unclass(capscale(d, meta$Management, cpus = 1))[, ], unclass(one)[, ])

  # dist_ord() runs it, by its name or another.
  expect_identical(no_display(dist_ord(d, 'capscale', 'Management', meta, cpus = 1)), no_display(one))
  expect_identical(no_display(dist_ord(d, 'capscale', 'Management', meta, seed = 3, cpus = 1)), no_display(one))
  expect_identical(match_dist_ord('CAP')$id, 'capscale')
  expect_identical(match_dist_ord('constrained analysis of principal coordinates')$id, 'capscale')

  expect_error(capscale(d),                          '`groups` is required')
  expect_error(capscale(d, meta$Managment, meta),    '`groups` is NULL')
  expect_error(capscale(d, meta, ~ A1),               'capscale\\(dm, data = meta')
  expect_error(capscale(d, 'Management', meta, n_axes = 0), '`n_axes`')
  expect_error(capscale(d, 'Management', meta, cpus = 0))

  # Aliased columns are dropped with a message.
  dat <- cbind(meta, A1x2 = 2 * meta$A1)
  expect_message(capscale(d, data = dat, formula = ~ A1 + A1x2, cpus = 1), 'linear combinations')

  # Only conditions.
  oc <- capscale(d, data = meta, formula = ~ Condition(A1), cpus = 1)
  expect_true(all(startsWith(colnames(oc), 'MDS')))
  expect_equal(attr(oc, 'trace'), sum(attr(oc, 'inertia')))
  expect_identical(attr(oc, 'r_squared'), 0)



  # Deterministic, for any number of threads, on both paths. ====

  set.seed(1)
  n      <- 150
  cnt    <- matrix(rpois(n * 30, 2), n, dimnames = list(paste0('S', 1:n), NULL))
  env    <- data.frame(row.names = rownames(cnt), g = rep(c('a', 'b', 'c'), 50), x = rnorm(n))
  for (dd in list(bray(cnt), euclidean(cnt))) {
    o1 <- capscale(dd, data = env, formula = ~ g + x, n_axes = 4, cpus = 1)
    expect_identical(capscale(dd, data = env, formula = ~ g + x, n_axes = 4, cpus = 1),        o1)
    expect_identical(capscale(dd, data = env, formula = ~ g + x, n_axes = 4, cpus = 2),        o1)
    expect_identical(capscale(dd, data = env, formula = ~ g + x, n_axes = 4, cpus = max_cpus), o1)
  }



  # Edge cases. ====

  # Two samples.
  ord <- suppressWarnings(capscale(dist(c(a = 1, b = 3)), c(a = 1, b = 2), cpus = 1))
  expect_true(all(attr(ord, 'inertia')[1:3] >= 0))
  expect_identical(attr(ord, 'inertia')[['imaginary']], 0)

  # All distances 0.
  z0  <- as.dist(matrix(0, 5, 5))
  ord <- suppressWarnings(capscale(z0, 1:5, cpus = 1))
  expect_identical(ncol(ord), 0L)
  expect_identical(attr(ord, 'r_squared'), NA_real_)

#})
