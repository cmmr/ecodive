#test_that("dmm", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L


  # The model in base R, independently of src/dmm.c: the negative log
  # posterior (with the multinomial coefficients and the prior's constants)
  # and the Laplace approximation, from the Dirichlet parameters `alpha`
  # (types x features), weights `w` and posteriors `z` (samples x types).
  dmm_nle <- function (X, alpha, w, eta = 0.1, nu = 0.1) {
    J  <- rowSums(X)
    ll <- sapply(seq_len(nrow(alpha)), function (k) {
      a <- alpha[k, ]
      A <- sum(a)
      lgamma(A) - lgamma(A + J) + rowSums(lgamma(sweep(X, 2, a, '+'))) - sum(lgamma(a)) })
    ll <- sweep(matrix(ll, nrow(X)), 2, log(w), '+')
    m  <- apply(ll, 1, max)
    lse <- m + log(rowSums(exp(ll - m)))
    mc <- lgamma(J + 1) - rowSums(lgamma(X + 1))
    prior <- sum(eta * log(nu) - lgamma(eta) + eta * log(alpha) - nu * alpha)
    -sum(lse + mc) - prior
  }

  dmm_laplace <- function (X, alpha, w, z, nle, eta = 0.1, nu = 0.1) {
    J  <- rowSums(X)
    K  <- nrow(alpha)
    S  <- ncol(alpha)
    ld <- (K - 1) * log(nrow(X)) - sum(log(w))
    for (k in seq_len(K)) {
      a  <- alpha[k, ]
      A  <- sum(a)
      zk <- z[, k]
      Xa <- sweep(X, 2, a, '+')
      t  <- colSums(zk * (digamma(Xa) - matrix(digamma(a), nrow(X), S, byrow = TRUE))) -
            sum(zk * (digamma(A + J) - digamma(A)))
      u  <- colSums(zk * (trigamma(Xa) - matrix(trigamma(a), nrow(X), S, byrow = TRUE)))
      cc <- sum(zk * (trigamma(A + J) - trigamma(A)))
      H  <- cc * outer(a, a)
      diag(H) <- diag(H) - a * t - a^2 * u + nu * a
      ld <- ld + as.numeric(determinant(H, logarithm = TRUE)$modulus)
    }
    M <- K * S + K - 1
    nle + ld / 2 - M / 2 * log(2 * pi)
  }

  rel <- function (x, y) max(abs(x - y) / abs(y))



  # Results from DirichletMultinomial (data-raw/expected/dmm.r). ====

  if (at_home()) {

    ref <- readRDS(file.path('expected', 'dmm.rds'))

    for (case in names(ref)) {

      X <- ref[[case]]$counts

      for (f in ref[[case]]$fits) {

        K    <- f$k
        info <- paste(case, K)
        o    <- order(-f$pi)
        z0   <- f$z[, o, drop = FALSE]
        a0   <- t(f$alpha[, o, drop = FALSE])

        # From DirichletMultinomial's posteriors as the starting weights,
        # the same fit, to the precision of its optimizer.
        d   <- ecodive:::dmm_core(X, K, init = f$z, cpus = 1L)
        sel <- attr(d, 'selection')
        a   <- attr(d, 'alpha')
        expect_true(rel(sel$nle, f$gof[['NLE']]) < 1e-8, info = info)
        expect_true(sel$nle <= f$gof[['NLE']] + 1e-9 * abs(f$gof[['NLE']]), info = info)
        expect_true(rel(sel$aic, f$gof[['AIC']]) < 1e-8, info = info)
        expect_true(rel(sel$bic, f$gof[['BIC']]) < 1e-8, info = info)
        expect_true(max(abs(a - a0) / a0) < 1e-2, info = info)
        expect_true(median(abs(a - a0) / a0) < 1e-4, info = info)
        expect_equal(attr(d, 'types')$weight, f$pi[o], tolerance = 1e-4, info = info)
        expect_true(max(abs(attr(d, 'posterior') - z0)) < 1e-3, info = info)
        gap <- apply(z0, 1, function (p) if (length(p) > 1) -diff(sort(p, decreasing = TRUE)[1:2]) else 1)
        expect_identical(
          as.integer(d$type)[gap > 1e-3], max.col(z0, 'first')[gap > 1e-3], info = info )

        # For one type, the Laplace approximation is DirichletMultinomial's.
        if (K == 1) {
          expect_true(rel(sel$laplace, f$gof[['Laplace']]) < 1e-6, info = info)
          expect_true(max(abs(a - a0) / a0) < 1e-3, info = info)
        }

        # The criteria agree with the model computed in base R.
        nle <- dmm_nle(X, a, attr(d, 'types')$weight)
        expect_equal(sel$nle, nle, tolerance = 1e-10, info = info)
        expect_equal(
          sel$laplace,
          dmm_laplace(X, a, attr(d, 'types')$weight, attr(d, 'posterior'), sel$nle),
          tolerance = 1e-8, info = info )
        expect_equal(sel$aic, sel$nle + K * ncol(X) + K - 1, tolerance = 1e-12, info = info)
        expect_equal(
          sel$bic, sel$nle + (K * ncol(X) + K - 1) / 2 * log(nrow(X)),
          tolerance = 1e-12, info = info )

        # The fit is a stationary point: the M-step's gradient is ~0 for the
        # final posteriors, up to EM's tolerance.
        for (k in seq_len(K)) {
          ak <- a[k, ]
          A  <- sum(ak)
          zk <- attr(d, 'posterior')[, k]
          J  <- rowSums(X)
          g  <- -ak * (colSums(zk * (digamma(sweep(X, 2, ak, '+')) -
                                     matrix(digamma(ak), nrow(X), ncol(X), byrow = TRUE))) -
                       sum(zk * (digamma(A + J) - digamma(A)))) + 0.1 * ak - 0.1
          expect_true(max(abs(g)) < 1e-2, info = paste(info, k))
        }
      }
    }

    # From its own k-means start, dmm() finds DirichletMultinomial's fits of
    # the simulated types, and recovers them.
    X   <- ref$sim$counts
    d   <- dmm(X, 1:4, cpus = 1L)
    sel <- attr(d, 'selection')
    dmn <- vapply(ref$sim$fits, function (f) f$gof[['NLE']], 0)
    expect_true(all(abs(sel$nle[1:3] - dmn[1:3]) / dmn[1:3] < 1e-8))
    expect_identical(sel$n_types[sel$selected], 3L)
    tab <- table(d$type, ref$sim$groups) > 0
    expect_true(all(rowSums(tab) == 1) && all(colSums(tab) == 1))
  }



  # A known two-type mixture is recovered. ====

  set.seed(2)
  p1 <- c(30, 10, 5, 2, 1, 1, 1)
  p2 <- c(1, 2, 5, 10, 30, 5, 5)
  draw <- function (p, n) t(sapply(seq_len(n), function (i)
    stats::rmultinom(1, sample(300:900, 1), stats::rgamma(length(p), p))[, 1]))
  X2 <- rbind(draw(p1, 70), draw(p2, 40))
  dimnames(X2) <- list(sprintf('s%03d', seq_len(nrow(X2))), paste0('t', 1:7))
  g2 <- rep(1:2, c(70, 40))

  d2 <- dmm(X2, 1:4, cpus = 1L)
  expect_inherits(d2, c('ecodive_dmm', 'ecodive_types', 'data.frame'))
  expect_identical(names(d2), c('sample', 'type', 'probability'))
  expect_identical(rownames(d2), rownames(X2))
  expect_identical(d2$sample, rownames(X2))
  expect_identical(levels(d2$type), c('1', '2'))
  expect_identical(as.integer(d2$type), g2)
  expect_true(all(d2$probability > 0.99))
  sel <- attr(d2, 'selection')
  expect_identical(names(sel), c('n_types', 'laplace', 'aic', 'bic', 'nle', 'n_iters', 'converged', 'selected'))
  expect_identical(sel$n_types, 1:4)
  expect_identical(sel$selected, c(FALSE, TRUE, FALSE, FALSE))
  expect_true(all(sel$converged))
  expect_identical(attr(d2, 'criterion'), 'laplace')
  expect_identical(attr(d2, 'method'), 'dmm')

  ty <- attr(d2, 'types')
  expect_identical(names(ty), c('type', 'n_samples', 'top_feature', 'weight', 'precision'))
  expect_identical(ty$n_samples, c(70L, 40L))
  expect_identical(ty$top_feature, c('t1', 't5'))
  expect_equal(ty$weight, c(70, 40) / 110, tolerance = 1e-6)
  expect_equal(sum(ty$weight), 1, tolerance = 1e-12)
  al <- attr(d2, 'alpha')
  expect_identical(dimnames(al), list(c('1', '2'), paste0('t', 1:7)))
  expect_equal(ty$precision, unname(rowSums(al)), tolerance = 1e-12)
  pr <- attr(d2, 'profiles')
  expect_equal(pr, al / rowSums(al), tolerance = 1e-12)
  expect_equal(unname(rowSums(pr)), c(1, 1), tolerance = 1e-12)
  # The Dirichlet means are close to the generating proportions.
  expect_true(max(abs(pr[1, ] - p1 / sum(p1))) < 0.05)
  expect_true(max(abs(pr[2, ] - p2 / sum(p2))) < 0.05)
  post <- attr(d2, 'posterior')
  expect_identical(dimnames(post), list(rownames(X2), c('1', '2')))
  expect_equal(unname(rowSums(post)), rep(1, 110), tolerance = 1e-12)
  expect_equal(d2$probability, post[cbind(1:110, as.integer(d2$type))])

  # The selection's criteria, and the base-R model.
  for (K in 1:4) {
    dk <- dmm(X2, K, cpus = 1L)
    s  <- attr(dk, 'selection')
    sk <- sel[K, , drop = FALSE]
    sk$selected <- TRUE
    rownames(sk) <- NULL
    expect_equal(s, sk, info = K)
    expect_equal(s$nle, dmm_nle(X2, attr(dk, 'alpha'), attr(dk, 'types')$weight), tolerance = 1e-10)
    expect_equal(
      s$laplace,
      dmm_laplace(X2, attr(dk, 'alpha'), attr(dk, 'types')$weight, attr(dk, 'posterior'), s$nle),
      tolerance = 1e-8 )
  }

  # Choosing by AIC or BIC.
  da <- dmm(X2, 1:4, criterion = 'aic', cpus = 1L)
  expect_identical(attr(da, 'criterion'), 'aic')
  expect_identical(which(attr(da, 'selection')$selected), which.min(sel$aic))
  db <- dmm(X2, 1:4, criterion = 'b', cpus = 1L)
  expect_identical(attr(db, 'criterion'), 'bic')
  expect_identical(which(attr(db, 'selection')$selected), which.min(sel$bic))
  expect_identical(attr(db, 'selection')[, 1:7], sel[, 1:7])

  # A weaker prior changes the fit, and the prior's constants.
  dp <- dmm(X2, 2, eta = 1, nu = 0.01, cpus = 1L)
  expect_false(isTRUE(all.equal(attr(dp, 'alpha'), attr(dmm(X2, 2, cpus = 1L), 'alpha'))))
  expect_equal(
    attr(dp, 'selection')$nle,
    dmm_nle(X2, attr(dp, 'alpha'), attr(dp, 'types')$weight, eta = 1, nu = 0.01),
    tolerance = 1e-10 )
  expect_equal(
    attr(dp, 'selection')$laplace,
    dmm_laplace(X2, attr(dp, 'alpha'), attr(dp, 'types')$weight, attr(dp, 'posterior'),
                attr(dp, 'selection')$nle, eta = 1, nu = 0.01),
    tolerance = 1e-8 )



  # Threads, table types and invariances. ====

  X3 <- rbind(X2, draw(c(2, 2, 20, 20, 2, 1, 1), 40))
  rownames(X3) <- sprintf('s%03d', seq_len(nrow(X3)))
  d3 <- dmm(X3, 1:4, n_starts = 2L, cpus = 1L)
  for (cpus in 2:max_cpus)
    expect_identical(dmm(X3, 1:4, n_starts = 2L, cpus = cpus), d3, info = cpus)

  # Sparse tables, and samples in columns.
  expect_identical(dmm(Matrix::Matrix(X3, sparse = TRUE), 1:4, n_starts = 2L, cpus = 1L), d3)
  expect_identical(dmm(slam::as.simple_triplet_matrix(X3), 1:4, n_starts = 2L, cpus = 1L), d3)
  expect_identical(dmm(t(X3), 1:4, n_starts = 2L, margin = 2L, cpus = max_cpus), d3)
  X3d <- X3
  storage.mode(X3d) <- 'double'
  expect_identical(dmm(X3d, 1:4, n_starts = 2L, cpus = 1L), d3)

  # The same seed gives the same result; the start is drawn from it.
  expect_identical(dmm(X3, 3, seed = 5L, cpus = 1L), dmm(X3, 3, seed = 5L, cpus = 1L))

  # Feature order: the same fit, features permuted.
  pf <- c(7, 3, 1, 6, 2, 5, 4)
  df <- dmm(X3[, pf], 3, cpus = 1L)
  d33 <- dmm(X3, 3, cpus = 1L)
  expect_equal(attr(df, 'alpha'), attr(d33, 'alpha')[, pf], tolerance = 1e-6)
  expect_equal(attr(df, 'posterior'), attr(d33, 'posterior'), tolerance = 1e-6)
  expect_equal(attr(df, 'selection')$laplace, attr(d33, 'selection')$laplace, tolerance = 1e-9)

  # Sample order: from the same starting weights, the same fit, permuted.
  ps <- c(seq(150, 1, by = -2), seq(1, 150, by = 2))
  z0 <- attr(d33, 'posterior')
  sa <- ecodive:::dmm_core(X3, 3, init = z0, cpus = 1L)
  sb <- ecodive:::dmm_core(X3[ps, ], 3, init = z0[ps, ], cpus = 1L)
  expect_equal(attr(sb, 'posterior'), attr(sa, 'posterior')[ps, ], tolerance = 1e-6)
  expect_equal(attr(sb, 'alpha'), attr(sa, 'alpha'), tolerance = 1e-6)
  expect_identical(sb$type, sa$type[ps])
  # And from k-means starts, the same types (these are well separated).
  tab <- table(dmm(X3[ps, ], 3, cpus = 1L)$type, d33$type[ps]) > 0
  expect_true(all(rowSums(tab) == 1) && all(colSums(tab) == 1))



  # Edge cases and errors. ====

  # A feature absent from every sample gets a small Dirichlet parameter.
  X4 <- cbind(X2, t8 = 0L)
  d4 <- dmm(X4, 2, cpus = 1L)
  expect_true(all(attr(d4, 'alpha')[, 't8'] < 0.01))
  expect_identical(as.integer(d4$type), g2)

  # Two and three samples.
  d5 <- dmm(X2[c(1, 80), ], cpus = 1L)
  expect_identical(attr(d5, 'selection')$n_types, 1:2)
  expect_identical(nrow(d5), 2L)
  d6 <- dmm(X2[c(1, 2, 80), ], n_types = 1:3, cpus = 1L)
  expect_identical(attr(d6, 'selection')$n_types, 1:3)
  expect_true(all(is.finite(attr(d6, 'selection')$nle)))

  # The default 1:7 is cut to the number of samples.
  expect_identical(attr(dmm(X2[1:5, ], cpus = 1L), 'selection')$n_types, 1:5)

  # Unnamed tables number their samples and features.
  d7 <- dmm(unname(X2), 2, cpus = 1L)
  expect_identical(d7$sample, as.character(1:110))
  expect_identical(colnames(attr(d7, 'alpha')), as.character(1:7))

  expect_error(dmm(X2, 0:2), 'between 1 and 110')
  expect_error(dmm(X2, 111), 'between 1 and 110')
  expect_error(dmm(X2, 1.5), 'whole numbers')
  expect_error(dmm(X2, NA), 'whole numbers')
  expect_error(dmm(X2, 'a'), 'whole numbers')
  expect_error(dmm(X2, criterion = 'x'), 'should be one of')
  expect_error(dmm(X2, n_starts = 0), '`n_starts` must be')
  expect_error(dmm(X2, max_iters = 1.5), '`max_iters` must be')
  expect_error(dmm(X2, tol = -1), '`tol` must be')
  expect_error(dmm(X2, eta = 0), '`eta` must be')
  expect_error(dmm(X2, nu = c(1, 2)), '`nu` must be')
  expect_error(dmm(X2, seed = 1.5), '`seed` must be')
  expect_error(dmm(X2[1, , drop = FALSE]), 'at least 2 samples')
  expect_error(dmm(X2[, 1, drop = FALSE]), 'at least 2 features')
  expect_error(dmm(rbind(X2, s999 = 0L)), "positive total.*s999")
  expect_error(dmm(replace(X2, 1, -1L)), 'negative')
  expect_error(dmm(replace(X2, 1, NA)), 'NA')
  expect_error(dmm(X2[c(1, 1, 1, 1), ], 2), 'fewer distinct samples')
  expect_warning(dmm(X2 / 7, 1, cpus = 1L), 'not whole numbers')
  expect_warning(dmm(X3, 4, max_iters = 1L, cpus = 1L), 'did not converge')
  expect_identical(
    suppressWarnings(attr(dmm(X3, 4, max_iters = 1L, cpus = 1L), 'selection')$converged), FALSE )



  # print() and plot(). ====

  out <- capture.output(print(d2))
  expect_identical(out[[1]], '# Dirichlet-Multinomial Mixtures')
  expect_identical(out[[2]], '# 110 samples in 2 types, chosen by the Laplace approximation from 1 to 4')
  expect_true(any(out == '# Types:'))
  expect_true(any(out == '# Selection:'))
  expect_true(any(grepl('^# Samples: a data.frame: 110 . 3$', out)))
  expect_true(any(grepl('^# Attributes: types, profiles, selection, posterior, alpha, criterion, method$', out)))
  out <- capture.output(print(db))
  expect_true(any(grepl('chosen by the BIC from 1 to 4', out)))
  out <- capture.output(print(dmm(X2, 2, cpus = 1L)))
  expect_identical(out[[2]], '# 110 samples in 2 types')
  expect_false(any(out == '# Selection:'))

  grDevices::pdf(NULL)
  expect_identical(plot(d2), d2)
  expect_identical(plot(d2, 'selection'), d2)
  expect_identical(plot(db, 'selection', main = 'Mine'), db)
  expect_identical(plot(d2, n_features = 3), d2)
  expect_error(plot(d2, 'scores'), 'Only valencia')
  grDevices::dev.off()

#})
