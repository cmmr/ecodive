#test_that("vegdist metrics", {

  # Fixtures: data-raw/expected/vegdist_metrics.r (vegan 2.7-8, output
  # only; abdiv 0.2.0; SciPy's Mahalanobis).

  fx  <- readRDS(file.path('expected', 'vegdist_metrics.rds'))
  tol <- 1e-12

  # Equal within `tol` (relative to values above 1), with NaN in the same
  # places.
  same <- function (a, b, tol = 1e-12) {
    a <- as.vector(a); b <- as.vector(b)
    identical(is.na(a), is.na(b)) &&
      isTRUE(all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol * pmax(1, abs(b[!is.na(b)]))))
  }

  # vegan's NA for pairs with an empty sample, where ecodive's differs on
  # purpose (documented): chisq_distance() gives finite values.
  empty_pair <- function (x) {
    e <- rowSums(x) == 0
    as.vector(stats::as.dist(outer(e, e, `|`)))
  }



  # vegan's vegdist(), every table and method ====

  for (tab in names(fx$tables)) {
    x <- fx$tables[[tab]]
    for (nm in names(fx$vegan[[tab]])) {
      binary <- endsWith(nm, '_binary')
      f      <- sub('_binary$', '', nm)
      d      <- if (binary) do.call(f, list(x, norm = 'binary', cpus = 1L))
                else        do.call(f, list(x, cpus = 1L))
      ref    <- fx$vegan[[tab]][[nm]]
      tl     <- if (f == 'mountford') 1e-9 else tol
      if (f == 'chisq_distance') {
        ep <- empty_pair(x)
        expect_true(all(is.finite(d[ep])), info = paste(tab, nm))
        d <- d[!ep]; ref <- ref[!ep]
      }
      expect_true(same(d, ref, tl), info = paste(tab, nm))
    }
  }

  # soergel() is vegdist()'s quantitative "jaccard".
  expect_true(same(soergel(fx$tables$dune), fx$vegan$dune$soergel))

  # Mahalanobis: vegan and SciPy.
  for (tab in names(fx$maha)) {
    d <- mahalanobis_distance(fx$maha[[tab]], cpus = 1L)
    expect_true(same(d, fx$vegan_maha[[tab]], 1e-10), info = tab)
    expect_true(same(d, fx$scipy_maha[[tab]], 1e-10), info = tab)
  }
  expect_true(same(
    mahalanobis_distance(fx$maha$counts > 2, norm = 'binary'),
    fx$vegan_maha$counts_binary, 1e-10 ))



  # abdiv ====

  for (tab in names(fx$abdiv)) {
    x <- fx$tables[[tab]]
    for (f in names(fx$abdiv[[tab]]))
      expect_true(same(do.call(f, list(x)), fx$abdiv[[tab]][[f]]), info = paste(tab, f))
  }



  # Protocol grids ====

  pa_pair <- function (a, b, j, S = a + b - j) {
    x <- matrix(0, 2, S)
    x[1, seq_len(a)] <- 1
    if (b > 0) x[2, (a - j) + seq_len(b)] <- 1
    x
  }

  # Mountford: vegan stops its iterations at about 1e-10.
  g <- fx$mountford_grid
  m <- mapply(function (a, b, j) as.vector(mountford(pa_pair(a, b, j, a + b - j + 1))), g$a, g$b, g$j)
  expect_true(same(m, g$vegan, 1e-9))
  expect_true(all(is.nan(m[g$a == 0 | g$b == 0])))
  expect_true(all(m[g$a > 0 & g$b > 0 & g$j == 0] == log(2)))
  expect_true(all(m[g$j > 0 & (g$j == g$a | g$j == g$b)] == 0))

  # Mountford is log(2) - M, M the root, beyond where vegan converges.
  for (r in list(c(574, 571, 566), c(1616, 2493, 1572), c(3, 3000, 2))) {
    M <- log(2) - as.vector(mountford(pa_pair(r[1], r[2], r[3])))
    expect_true(abs(exp(-(r[2] - r[3]) * M) + exp(-(r[1] - r[3]) * M) - exp(-(r[1] + r[2] - r[3]) * M) - 1) < 1e-13)
  }

  # Raup-Crick: vegan, and phyper()'s upper tail at full precision.
  g <- fx$raup_grid
  r <- mapply(function (S, a, b, j) as.vector(raup_crick(pa_pair(a, b, j, S))), g$S, g$a, g$b, g$j)
  expect_true(same(r, g$vegan, 1e-12))
  p <- stats::phyper(g$j - 1, g$a, g$S - g$a, g$b, lower.tail = FALSE)
  expect_true(all(abs(r - p) <= 1e-12 * p + 1e-300))

  # An empty feature changes S, and so the result.
  x <- fx$tables$ex_counts
  expect_false(isTRUE(all.equal(raup_crick(x), raup_crick(cbind(x, 0)))))
  pa <- x > 0
  expect_true(same(
    as.matrix(raup_crick(x))[2, 1],
    stats::phyper(sum(pa[1, ] & pa[2, ]) - 1, sum(pa[1, ]), ncol(x) - sum(pa[1, ]), sum(pa[2, ]), lower.tail = FALSE) ))



  # Identities in plain R ====

  x <- fx$tables$dune
  pa <- (x > 0) * 1
  A <- rowSums(pa); J <- tcrossprod(pa)
  kb <- 1 - (J / A + t(J / A)) / 2
  expect_true(same(kulczynski(x, norm = 'binary'), stats::as.dist(kb)))
  expect_true(same(binomial_deviance(x, norm = 'binary'), log(2) * hamming(x)))
  expect_true(same(alt_gower(x, norm = 'binary'), jaccard(x)))

  y <- fx$maha$varechem
  w <- scale(y, scale = FALSE) %*% solve(chol(stats::cov(y)))
  expect_true(same(mahalanobis_distance(y), stats::dist(w), 1e-10))

  z <- fx$tables$sparse
  z <- z[rowSums(z) > 0, colSums(z) > 0]
  ch <- sqrt(sum(z)) * sweep(z / rowSums(z), 2, sqrt(colSums(z)), '/')
  expect_true(same(chisq_distance(z), stats::dist(ch)))
  expect_true(same(euclidean(z, norm = 'chisq'), chisq_distance(z)))
  expect_true(same(normalize(z, 'chisq'), ch))
  expect_identical(
    as.vector(chisq_distance(fx$tables$sparse)),
    as.vector(euclidean(fx$tables$sparse, norm = 'chisq')) )

  cao_r <- function (x, y) {
    k <- x > 0 | y > 0
    x <- pmax(x[k], 0.1); y <- pmax(y[k], 0.1); n <- x + y
    mean(log(n / 2) - (x * log(y) + y * log(x)) / n)
  }
  e <- fx$tables$ex_counts
  expect_true(same(as.matrix(cao(e))[3, 2], cao_r(e[2, ], e[3, ])))
  # Values below 0.1 are raised to it.
  expect_true(same(cao(rbind(c(0.05, 1, 0), c(0, 2, 3))), cao_r(c(0.05, 1, 0), c(0, 2, 3))))



  # Wrapper, aliases and registry ====

  e <- fx$tables$ex_counts
  for (a in list(
      c('binomial', 'binomial_deviance'), c('chao', 'chao_jaccard'),
      c('mahalanobis', 'mahalanobis_distance'), c('raup', 'raup_crick'),
      c('altGower', 'alt_gower'), c('chisq', 'chisq_distance'),
      c('kulczynski', 'kulczynski'), c('mountford', 'mountford'),
      c('cao', 'cao'), c('CYd', 'cao'), c('chi-square', 'chisq_distance') )) {
    expect_identical(match_metric(a[1], div = 'beta')$id, a[2], info = a[1])
  }
  expect_identical(match_metric('chao', div = 'alpha')$id, 'chao1')
  # chao_jaccard's 'chao' alias must not hide chao1 from prefix matching.
  expect_identical(match_metric('c', div = 'alpha')$id, 'chao1')
  expect_identical(match_metric('chaoj')$id, 'chao_jaccard')
  expect_identical(attr(beta_div(e, 'binomial'), 'method'), 'binomial_deviance')
  expect_identical(as.vector(beta_div(e, 'raup')), as.vector(raup_crick(e)))
  expect_identical(as.vector(beta_div(e, 'kulczynski', norm = 'binary')),
                   as.vector(kulczynski(e, norm = 'binary')))
  expect_identical(attr(chao_jaccard(e), 'display')$title, 'Chao-Jaccard Dissimilarity')
  expect_identical(attr(chisq_distance(e), 'method'), 'chisq_distance')
  expect_true(all(c('mountford', 'raup_crick') %in% list_metrics('beta', 'id', weighted = FALSE)))
  expect_true(all(c('mahalanobis_distance', 'chisq_distance') %in% list_metrics('beta', 'id', true_metric = TRUE)))

  # A brute-force search for triangle-inequality violations: the metrics
  # registered as not metric have one.
  set.seed(1)
  viol <- function (d) {
    m <- as.matrix(d); n <- nrow(m)
    for (i in 1:n) for (j in 1:n) for (k in 1:n)
      if (m[i, k] > m[i, j] + m[j, k] + 1e-12) return (TRUE)
    FALSE
  }
  r <- matrix(rpois(15 * 8, 2), 15, 8)
  for (f in c('kulczynski', 'binomial_deviance', 'cao', 'mountford'))
    expect_true(viol(do.call(f, list(r))), info = f)
  expect_true(viol(alt_gower(rbind(c(0, 0, 0, 1), c(1, 0, 1, 1), c(0, 5, 0, 0)))))
  expect_false(viol(chisq_distance(r)))



  # Threads, pairs, sparse input ====

  x  <- fx$tables$sparse
  xm <- x[, 1:10]
  xm <- xm[rowSums(xm) > 0, ]
  xm <- xm + matrix(rpois(length(xm), 3), nrow(xm))
  fns <- c('kulczynski', 'mountford', 'raup_crick', 'binomial_deviance',
           'chao_jaccard', 'cao', 'alt_gower', 'chisq_distance', 'mahalanobis_distance')
  for (f in fns) {
    xx <- if (f == 'mahalanobis_distance') xm else x
    d1 <- do.call(f, list(xx, cpus = 1L))
    for (cp in c(2L, 4L, 7L))
      expect_identical(do.call(f, list(xx, cpus = cp)), d1, info = paste(f, cp))
    expect_identical(do.call(f, list(t(xx), margin = 2L, cpus = 2L)), d1, info = f)
    pr <- c(1L, 5L, 9L)
    dp <- do.call(f, list(xx, pairs = pr, cpus = 2L))
    expect_identical(as.vector(dp)[pr], as.vector(d1)[pr], info = f)
    expect_true(all(is.na(as.vector(dp)[-pr])), info = f)
    if (requireNamespace('Matrix', quietly = TRUE)) {
      expect_identical(as.vector(do.call(f, list(Matrix::Matrix(xx, sparse = TRUE), cpus = 2L))),
                       as.vector(d1), info = f)
    }
    if (requireNamespace('slam', quietly = TRUE)) {
      expect_identical(as.vector(do.call(f, list(slam::as.simple_triplet_matrix(xx), cpus = 2L))),
                       as.vector(d1), info = f)
    }
  }



  # Errors ====

  neg <- e; neg[1, 1] <- -1
  for (f in setdiff(fns, 'mahalanobis_distance'))
    expect_error(do.call(f, list(neg)), 'negative', info = f)
  expect_error(chao_jaccard(e + 0.5), 'whole numbers')
  expect_silent(cao(e + 0.5))
  expect_error(mahalanobis_distance(e), 'more samples than features')
  expect_error(mahalanobis_distance(cbind(xm, xm[, 1] + xm[, 2])), 'collinear')
  expect_error(mahalanobis_distance(cbind(xm, 3)), 'collinear')
  expect_error(mahalanobis_distance(xm, norm = 'percent'), 'collinear')
  expect_error(normalize(neg, 'chisq'), 'chi-square')
  expect_error(random_forest(e, 1:nrow(e), norm = 'chisq'), 'chisq')

#})
