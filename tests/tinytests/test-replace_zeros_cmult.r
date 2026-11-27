#test_that("replace_zeros_czm, _gbm, _sq, _bl", {

  ref <- readRDS(file.path('expected', 'replace_zeros_cmult.rds'))

  plain  <- function (x) { attributes(x) <- list(dim = dim(x)); x }
  maxrel <- function (a, b) max(abs(plain(a) - b) / abs(b))

  fns <- list(czm = replace_zeros_czm, gbm = replace_zeros_gbm,
              sq = replace_zeros_sq, bl = replace_zeros_bl)

  run <- function (cs, x, ...) {
    args <- list(x, frac = cs$frac, adjust = cs$adjust, ...)
    if (cs$method == 'czm') args$threshold <- cs$threshold
    suppressWarnings(do.call(fns[[cs$method]], args))  # probe_real isn't whole
  }


  # -----------------------------------------------------------------------
  # zCompositions 1.6.2 cmultRepl(output = "prop", z.delete = FALSE)
  # -----------------------------------------------------------------------

  expect_true(length(ref$cases) >= 60)
  for (cs in ref$cases) {
    x <- ref$tables[[cs$table]]
    r <- run(cs, x)
    expect_true(maxrel(r, cs$result) < 1e-13, info = cs$id)
    expect_identical(attr(r, 'method'), cs$method, info = cs$id)
    expect_identical(attr(r, 'n_replaced'), as.double(sum(x == 0)), info = cs$id)

    # Features x samples gives the transpose.
    r2 <- run(cs, t(x), margin = 2L)
    expect_identical(plain(r2), t(plain(r)), info = cs$id)

    # Sparse tables give the same bits.
    if (requireNamespace('Matrix', quietly = TRUE)) {
      expect_identical(run(cs, Matrix::Matrix(x, sparse = TRUE)), r, info = cs$id)
      expect_identical(plain(run(cs, Matrix::Matrix(t(x), sparse = TRUE), margin = 2L)),
                       plain(r2), info = cs$id)
    }
    if (requireNamespace('slam', quietly = TRUE))
      expect_identical(run(cs, slam::as.simple_triplet_matrix(x)), r, info = cs$id)
  }

  # The adjustment was exercised, for every method.
  for (m in names(fns)) {
    n_adj <- vapply(ref$cases, function (cs)
      if (cs$method == m && cs$adjust) attr(run(cs, ref$tables[[cs$table]]), 'n_adjusted') else 0,
      numeric(1) )
    expect_true(sum(n_adj) > 0, info = m)
  }

  # ex_counts, with its names.
  for (m in names(fns)) {
    cs <- Filter(function (cs) cs$table == 'ex_counts' && cs$method == m &&
                   cs$adjust && cs$frac == 0.65, ref$cases)[[1]]
    r <- fns[[m]](ecodive::ex_counts)
    expect_identical(dimnames(r), dimnames(ecodive::ex_counts), info = m)
    expect_true(maxrel(r, cs$result) < 1e-13, info = m)
  }


  # -----------------------------------------------------------------------
  # The formulas, written out
  # -----------------------------------------------------------------------

  X <- rbind(c(10, 0, 5, 3, 2), c(4, 6, 0, 0, 10), c(1, 20, 7, 8, 4), c(30, 2, 3, 0, 15))
  n <- rowSums(X); x <- X / n; N <- 4; D <- 5
  Tl <- (matrix(colSums(X), N, D, byrow = TRUE) - X) / (sum(n) - n)
  by_hand <- function (imp) {
    r <- x
    for (i in 1:N) {
      z <- X[i, ] == 0
      r[i, z]  <- imp[i, z]
      r[i, !z] <- x[i, !z] * (1 - sum(imp[i, z]))
    }
    r
  }
  s_bl  <- rep(D, N)
  s_sq  <- sqrt(n)
  s_gbm <- 1 / exp(rowMeans(log(Tl)))
  expect_equal(plain(replace_zeros_bl(X,  adjust = FALSE)), by_hand(s_bl  * Tl / (n + s_bl)),  tolerance = 1e-14)
  expect_equal(plain(replace_zeros_sq(X,  adjust = FALSE)), by_hand(s_sq  * Tl / (n + s_sq)),  tolerance = 1e-14)
  expect_equal(plain(replace_zeros_gbm(X, adjust = FALSE)), by_hand(s_gbm * Tl / (n + s_gbm)), tolerance = 1e-14)
  expect_equal(plain(replace_zeros_czm(X, adjust = FALSE)), by_hand(matrix(0.325 / n, N, D)), tolerance = 1e-14)
  expect_equal(plain(replace_zeros_czm(X, frac = 0.2, threshold = 1, adjust = FALSE)),
               by_hand(matrix(0.2 / n, N, D)), tolerance = 1e-14)

  # adjust: above the feature's smallest observed proportion (strictly), an
  # imputed value becomes frac times it. Here 0.5 * 1 / 8 = 1/16, the
  # smallest proportion of feature 2: not adjusted.
  Xt <- rbind(c(4, 0, 4), c(8, 1, 7), c(2, 3, 3))
  r <- replace_zeros_czm(Xt, frac = 0.5, threshold = 1)
  expect_identical(r[1, 2], 1 / 16)
  expect_identical(attr(r, 'n_adjusted'), 0)
  expect_error(replace_zeros_czm(Xt, threshold = 1.1), 'at most 1')
  r <- replace_zeros_czm(rbind(c(1, 0, 1), c(8, 1, 7), c(2, 3, 3)), frac = 0.5, threshold = 1)
  expect_identical(attr(r, 'n_adjusted'), 1)
  expect_equal(r[1, 2], 0.5 * 1 / 16, tolerance = 1e-15)


  # -----------------------------------------------------------------------
  # Properties
  # -----------------------------------------------------------------------

  set.seed(1)
  m <- matrix(rpois(30 * 20, 20), 30, 20, dimnames = list(paste0('S', 1:30), paste0('F', 1:20)))
  m[runif(length(m)) < 0.4] <- 0
  stopifnot(all(rowSums(m) > 0), all(colSums(m > 0) >= 2))

  for (id in names(fns)) {
    r <- fns[[id]](m)
    expect_true(is.matrix(r) && is.double(r), info = id)
    expect_identical(dimnames(r), dimnames(m), info = id)
    expect_identical(names(attributes(r)), c('dim', 'dimnames', 'method', 'n_adjusted', 'n_replaced'), info = id)
    expect_identical(attr(r, 'method'), id, info = id)
    expect_identical(attr(r, 'n_replaced'), as.double(sum(m == 0)), info = id)
    expect_true(attr(r, 'n_adjusted') >= 0 && attr(r, 'n_adjusted') <= sum(m == 0), info = id)
    expect_identical(attr(fns[[id]](m, adjust = FALSE), 'n_adjusted'), 0, info = id)
    expect_true(all(r > 0), info = id)
    expect_true(max(abs(rowSums(r) - 1)) < 1e-15, info = id)

    # Non-zero parts keep their ratios.
    for (i in c(1, 17, 30)) {
      nz <- which(m[i, ] > 0)
      expect_equal(unname(r[i, nz] / r[i, nz[1]]), unname(m[i, nz] / m[i, nz[1]]), tolerance = 1e-14, info = id)
    }

    # Same inputs, same bits; and through the wrapper.
    expect_identical(fns[[id]](m), r, info = id)
    expect_identical(replace_zeros(m, id), r, info = id)
    expect_identical(replace_zeros(t(m), id, margin = 2L, frac = 0.5),
                     fns[[id]](t(m), margin = 2L, frac = 0.5), info = id)

    # Input types.
    mi <- m; storage.mode(mi) <- 'integer'
    expect_identical(fns[[id]](mi), r, info = id)
    expect_identical(fns[[id]](as.data.frame(m)), r, info = id)
    expect_null(dimnames(fns[[id]](unname(m))), info = id)
    expect_identical(plain(fns[[id]](unname(m))), plain(r), info = id)
    if (requireNamespace('Matrix', quietly = TRUE))
      expect_identical(fns[[id]](methods::as(Matrix::Matrix(m, sparse = TRUE), 'TsparseMatrix')), r, info = id)

    # No zeros: proportions.
    expect_equal(plain(fns[[id]](m + 1)), plain((m + 1) / rowSums(m + 1)), tolerance = 1e-15, info = id)
    expect_identical(attr(fns[[id]](m + 1), 'n_replaced'), 0, info = id)

    # Downstream: no pseudocount needed.
    expect_silent(aitchison(r))
    expect_silent(ilr(r))
  }

  # A sample's other samples set its prior: shallow samples get larger
  # values, and so does a feature that is abundant elsewhere.
  r <- replace_zeros_gbm(X, adjust = FALSE)
  expect_true(r[2, 3] > r[2, 4])  # feature 3 is more abundant elsewhere


  # -----------------------------------------------------------------------
  # Edge cases and errors
  # -----------------------------------------------------------------------

  # A feature non-zero in a single sample: GBM refuses, naming it; the
  # others handle it.
  X1 <- X; X1[, 2] <- c(0, 0, 5, 0)
  colnames(X1) <- paste0('f', 1:5)
  expect_error(replace_zeros_gbm(X1), "1 feature that is non-zero in only one sample \\('f2'\\).*use `sq` or `bl`")
  X2 <- X1; X2[, 4] <- c(0, 0, 0, 9)
  expect_error(replace_zeros_gbm(unname(X2)), "2 features that are non-zero in only one sample \\('2', '4'\\)")
  for (id in c('sq', 'bl', 'czm')) {
    r <- fns[[id]](X1)
    expect_true(all(r > 0) && max(abs(rowSums(r) - 1)) < 1e-15, info = id)
  }

  # Two samples: GBM refuses any zero (each zero feature is a single), the
  # others work.
  expect_error(replace_zeros_gbm(X[1:2, ]), 'non-zero in only one sample')
  expect_true(all(replace_zeros_bl(X[1:2, ]) > 0))

  # A feature that is zero everywhere: the Bayesian methods refuse it, czm
  # replaces it.
  Xz <- X1; Xz[, 2] <- 0
  for (id in c('gbm', 'sq', 'bl'))
    expect_error(fns[[id]](Xz), "1 feature that is zero in every sample \\('f2'\\)", info = id)
  r <- replace_zeros_czm(Xz)
  expect_equal(unname(r[, 2]), 0.325 / rowSums(Xz), tolerance = 1e-15)

  # One sample: only czm.
  for (id in c('gbm', 'sq', 'bl')) {
    expect_error(fns[[id]](X[1, , drop = FALSE]), 'at least 2 samples', info = id)
    expect_error(fns[[id]](c(1, 0, 3)), 'at least 2 samples', info = id)
  }
  v <- replace_zeros_czm(c(a = 1, b = 0, c = 3))
  expect_identical(dim(v), c(3L, 1L))
  expect_equal(as.vector(v), c((1 - 0.325 / 4) / 4, 0.325 / 4, (1 - 0.325 / 4) * 3 / 4), tolerance = 1e-15)

  # czm refuses zeros that would take up the whole sample.
  Xs <- rbind(a = c(1, 0, 0, 0, 0), b = c(3, 1, 1, 1, 1))
  expect_error(replace_zeros_czm(Xs, threshold = 1, frac = 0.9, adjust = FALSE),
               "negative proportions \\(or zeros\\) in 1 sample \\('a'\\).*smaller `frac` or `threshold`")
  expect_silent(replace_zeros_czm(Xs, threshold = 1, frac = 0.2, adjust = FALSE))

  # Counts are expected.
  for (id in names(fns))
    expect_warning(fns[[id]](m / 7), 'not whole numbers', info = id)

  me <- m; me[c(3, 9), ] <- 0
  for (id in names(fns))
    expect_error(fns[[id]](me), "2 samples with no non-zero values \\('S3', 'S9'\\)", info = id)
  mn <- m; mn[2, 2] <- -1
  expect_error(replace_zeros_gbm(mn), 'negative')
  mn[2, 2] <- NA
  expect_error(replace_zeros_sq(mn), 'NA')
  expect_error(replace_zeros_bl(m, margin = 3), 'margin')

  for (bad in list(0, 1, -1, 1.5, NA, NA_real_, 'a', c(0.1, 0.2), numeric(0)))
    for (id in names(fns))
      expect_error(fns[[id]](m, frac = bad), '`frac` must be a single number between 0 and 1', info = id)
  for (bad in list(0, -1, 1.01, NA, 'a', c(0.1, 0.2)))
    expect_error(replace_zeros_czm(m, threshold = bad), '`threshold` must be a single number above 0 and at most 1')
  expect_silent(replace_zeros_czm(m, threshold = 1))
  for (bad in list(NA, 'yes', 1, c(TRUE, FALSE)))
    for (id in names(fns))
      expect_error(fns[[id]](m, adjust = bad), '`adjust` must be TRUE or FALSE', info = id)


  # -----------------------------------------------------------------------
  # The wrapper and names
  # -----------------------------------------------------------------------

  mz <- function (x) match_zero_replacement(x)$id
  for (x in c('gbm', 'GBM', 'cmultRepl', 'geometric Bayesian-multiplicative', 'g'))
    expect_identical(mz(x), 'gbm', info = x)
  for (x in c('sq', 'SQ', 'square root', 'sqrt BM', 'Square-Root Bayesian-Multiplicative'))
    expect_identical(mz(x), 'sq', info = x)
  for (x in c('bl', 'BL', 'Bayes-Laplace', 'Laplace', 'b'))
    expect_identical(mz(x), 'bl', info = x)
  for (x in c('czm', 'CZM', 'count zero multiplicative', 'cz', 'count'))
    expect_identical(mz(x), 'czm', info = x)
  expect_error(mz('c'), 'matches multiple options: czm, gbm')  # cmultRepl
  for (x in character(0))
    expect_identical(mz(x), 'czm', info = x)

  expect_identical(match_zero_replacement('czm')$params,
                   c('counts', 'margin', 'frac', 'threshold', 'adjust'))
  for (id in c('gbm', 'sq', 'bl'))
    expect_identical(match_zero_replacement(id)$params, c('counts', 'margin', 'frac', 'adjust'), info = id)
  expect_identical(match_zero_replacement('gbm')$name, 'Geometric Bayesian-Multiplicative')

  expect_identical(replace_zeros(counts = m, method = 'cmultRepl'), replace_zeros_gbm(m))
  expect_identical(replace_zeros(m, 'czm', threshold = 0.3), replace_zeros_czm(m, threshold = 0.3))
  expect_error(replace_zeros(m, 'gbm', threshold = 0.3), 'unused argument')

#})
