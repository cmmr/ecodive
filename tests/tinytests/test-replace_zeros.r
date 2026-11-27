#test_that("replace_zeros", {

  ref <- readRDS(file.path('expected', 'replace_zeros.rds'))

  plain <- function (x) { attributes(x) <- list(dim = dim(x)); x }


  # -----------------------------------------------------------------------
  # scikit-bio 0.7.4 multi_replace(): bit-identical
  # -----------------------------------------------------------------------

  for (cs in ref$skbio) {
    r <- replace_zeros_mult(cs$x, delta = cs$delta)
    expect_identical(plain(r), cs$result, info = cs$id)
    expect_identical(attr(r, 'n_replaced'), as.double(sum(cs$x == 0)), info = cs$id)

    # Features x samples gives the transpose.
    r2 <- replace_zeros_mult(t(cs$x), margin = 2L, delta = cs$delta)
    expect_identical(plain(r2), t(cs$result), info = cs$id)
  }

  # The plan's worked example, rounded.
  r <- replace_zeros_mult(rbind(c(1, 2, 3, 0, 5), c(1, 0, 0, 4, 5), c(1, 2, 3, 4, 5)))
  expect_equal(round(plain(r), 6), rbind(
    c(.087273, .174545, .261818, .04, .436364),
    c(.092, .04, .04, .368, .46),
    c(.066667, .133333, .2, .266667, .333333) ))
  expect_equal(attr(r, 'delta'), 0.04, tolerance = 1e-15)
  expect_identical(attr(r, 'n_replaced'), 3)

  # zCompositions 1.6.2 multRepl() (output only), given the same delta.
  for (id in names(ref$zcomp))
    expect_equal(plain(replace_zeros_mult(ref$skbio[[id]]$x, delta = ref$skbio[[id]]$delta)),
                 ref$zcomp[[id]], tolerance = 1e-14, info = id)


  # -----------------------------------------------------------------------
  # Properties
  # -----------------------------------------------------------------------

  set.seed(1)
  m <- matrix(rpois(30 * 20, 2), 30, 20, dimnames = list(paste0('S', 1:30), paste0('F', 1:20)))
  m[m > 0 & runif(length(m)) < 0.3] <- 0
  m[rowSums(m) == 0, 1] <- 1
  stopifnot(any(m == 0))

  r <- replace_zeros_mult(m)
  expect_true(is.matrix(r) && is.double(r))
  expect_identical(dimnames(r), dimnames(m))
  expect_identical(names(attributes(r)), c('dim', 'dimnames', 'method', 'delta', 'n_replaced'))
  expect_identical(attr(r, 'method'), 'mult')
  expect_equal(attr(r, 'delta'), 1 / 400, tolerance = 1e-15)
  expect_identical(attr(r, 'n_replaced'), as.double(sum(m == 0)))
  expect_true(all(r > 0))
  expect_true(max(abs(rowSums(r) - 1)) < 1e-15)
  expect_true(all(r[m == 0] == attr(r, 'delta')))

  # Non-zero parts keep their ratios.
  for (i in 1:30) {
    nz <- which(m[i, ] > 0)
    expect_equal(unname(r[i, nz] / r[i, nz[1]]), unname(m[i, nz] / m[i, nz[1]]), tolerance = 1e-14)
  }

  # Counts and proportions give the same result.
  expect_equal(replace_zeros_mult(m / rowSums(m)), r, tolerance = 1e-15)


  # -----------------------------------------------------------------------
  # Input types
  # -----------------------------------------------------------------------

  mi <- m; storage.mode(mi) <- 'integer'
  expect_identical(replace_zeros_mult(mi), r)
  expect_identical(replace_zeros_mult(as.data.frame(m)), r)
  ml <- m > 0
  expect_identical(replace_zeros_mult(ml), replace_zeros_mult(ml * 1))

  # Unnamed input gives an unnamed result.
  expect_null(dimnames(replace_zeros_mult(unname(m))))
  expect_identical(plain(replace_zeros_mult(unname(m))), plain(r))

  # A vector is one sample: a features x 1 matrix.
  v <- replace_zeros_mult(c(a = 1, b = 0, c = 3))
  expect_identical(dim(v), c(3L, 1L))
  expect_equal(as.vector(v), c(8 / 9 * 1 / 4, 1 / 9, 8 / 9 * 3 / 4), tolerance = 1e-15)

  # Sparse tables give the same (dense) result.
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(m, sparse = TRUE)
    expect_identical(replace_zeros_mult(sp), r)
    r2 <- replace_zeros_mult(Matrix::t(sp), margin = 2L)
    expect_identical(plain(r2), t(plain(r)))
    expect_identical(dimnames(r2), rev(dimnames(r)))
    expect_identical(replace_zeros_mult(methods::as(sp, 'TsparseMatrix')), r)
    # Explicitly stored zeros are zeros.
    sp0 <- sp; sp0@x[1] <- 0
    m0 <- m; m0[cbind(sp@i[1] + 1L, 1L)] <- 0
    expect_identical(replace_zeros_mult(sp0), replace_zeros_mult(m0))
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(m)
    expect_identical(replace_zeros_mult(st), r)
  }


  # -----------------------------------------------------------------------
  # delta
  # -----------------------------------------------------------------------

  r1 <- replace_zeros_mult(m, delta = 1e-3)
  expect_identical(attr(r1, 'delta'), 1e-3)
  expect_true(all(r1[m == 0] == 1e-3))

  # One per sample: in order, or by name.
  d <- 0.5 / rowSums(m)
  rs <- replace_zeros_mult(m, delta = d)
  expect_identical(attr(rs, 'delta'), d)
  for (i in c(1, 7, 30))
    expect_identical(rs[i, ], replace_zeros_mult(m[i, , drop = FALSE], delta = d[[i]])[1, ])
  expect_identical(replace_zeros_mult(m, delta = rev(d)), rs)
  expect_identical(replace_zeros_mult(m, delta = unname(d)), rs)
  expect_identical(attr(replace_zeros_mult(unname(m), delta = unname(d)), 'delta'), unname(d))
  # One value for one sample is a scalar.
  expect_identical(attr(replace_zeros_mult(m[1, , drop = FALSE], delta = c(S1 = 0.01)), 'delta'), 0.01)

  # Too large for a sample's number of zeros: scikit-bio's error, naming the
  # samples. Exactly 1 (which would leave zeros) is an error too.
  z <- rbind(a = c(0, 0, 1, 1), b = c(1, 1, 1, 0))
  expect_error(replace_zeros_mult(z, delta = 0.6), "negative proportions \\(or zeros\\) in 1 sample \\('a'\\).*smaller `delta`")
  expect_error(replace_zeros_mult(z, delta = 0.5), "in 1 sample \\('a'\\)")
  expect_silent(replace_zeros_mult(z, delta = 0.49))
  expect_error(replace_zeros_mult(z, delta = c(0.6, 0.9)), "in 1 sample \\('a'\\)")
  expect_silent(replace_zeros_mult(z, delta = c(0.4, 0.9)))
  expect_error(replace_zeros_mult(unname(z), delta = 0.9), "in 1 sample \\('1'\\)")
  expect_error(replace_zeros_mult(rbind(z, z), delta = 0.9), "in 2 samples \\('a', 'a'\\)")

  for (bad in list(0, -0.1, 1, 2, NA_real_, Inf, NaN, 'a', TRUE, numeric(0), matrix(0.1, 1, 1)))
    expect_error(replace_zeros_mult(m, delta = bad), 'between 0 and 1', info = deparse(bad))
  expect_error(replace_zeros_mult(m, delta = c(0.1, 0.2)), '`delta` has 2 values, but `counts` has 30 samples')
  expect_error(replace_zeros_mult(m, delta = c(x = 0.1, S2 = 0.2)), '`delta`')


  # -----------------------------------------------------------------------
  # Edge cases and errors
  # -----------------------------------------------------------------------

  # A feature that is zero everywhere is replaced like any zero.
  mz <- m; mz[, 4] <- 0
  expect_true(all(replace_zeros_mult(mz)[, 4] == attr(r, 'delta')))

  # No zeros: proportions.
  expect_equal(plain(replace_zeros_mult(m + 1)), plain((m + 1) / rowSums(m + 1)), tolerance = 1e-15)
  expect_identical(attr(replace_zeros_mult(m + 1), 'n_replaced'), 0)

  # One feature.
  expect_identical(plain(replace_zeros_mult(matrix(c(2, 5), 2))), matrix(1, 2, 1))

  me <- m; me[c(3, 9), ] <- 0
  expect_error(replace_zeros_mult(me), "2 samples with no non-zero values \\('S3', 'S9'\\)")
  expect_error(replace_zeros_mult(unname(me)[1:3, ]), "1 sample with no non-zero values \\('3'\\)")
  mn <- m; mn[2, 2] <- -1
  expect_error(replace_zeros_mult(mn), 'negative')
  mn[2, 2] <- NA
  expect_error(replace_zeros_mult(mn), 'NA')
  mn[2, 2] <- Inf
  expect_error(replace_zeros_mult(mn), 'infinite')
  expect_error(replace_zeros_mult(matrix('a', 2, 2)), 'numeric')
  expect_error(replace_zeros_mult(m, margin = 3), 'margin')


  # -----------------------------------------------------------------------
  # Downstream: no pseudocount needed (item 14's log-ratio functions)
  # -----------------------------------------------------------------------

  expect_silent(a <- alr(r, 'F1'))
  expect_identical(attr(a, 'pseudocount'), 0)
  expect_equal(plain(a), plain(log(r[, -1] / r[, 1])), tolerance = 1e-14)
  expect_silent(i <- ilr(r))
  expect_identical(attr(i, 'pseudocount'), 0)
  expect_silent(cl <- normalize(r, 'clr'))
  expect_equal(plain(cl), plain(log(r) - rowMeans(log(r))), tolerance = 1e-14)
  expect_silent(ad <- aitchison(r))
  expect_equal(as.vector(ad), as.vector(dist(cl)), tolerance = 1e-12)
  expect_silent(pca(r, norm = 'clr', cpus = 1))

  # Orientation is kept from end to end.
  r2 <- replace_zeros_mult(t(m), margin = 2L)
  expect_identical(plain(ilr(r2, margin = 2L)), t(plain(i)))
  expect_identical(plain(alr(r2, 'F1', margin = 2L)), t(plain(a)))


  # -----------------------------------------------------------------------
  # The wrapper
  # -----------------------------------------------------------------------

  expect_identical(replace_zeros(m, 'mult'), replace_zeros_mult(m))
  expect_identical(replace_zeros(m, 'mult', delta = 1e-3), replace_zeros_mult(m, delta = 1e-3))
  expect_identical(replace_zeros(t(m), 'mult', margin = 2L), replace_zeros_mult(t(m), margin = 2L))
  expect_identical(replace_zeros(counts = m, method = 'mult'), replace_zeros_mult(m))
  expect_error(replace_zeros(m), '`method` is required')
  expect_error(replace_zeros(m, 'mult', seed = 1), 'unused argument')

  mz <- function (x) match_zero_replacement(x)$id
  for (x in c('mult', 'MULT', 'm', 'mu', 'multiplicative', 'Multiplicative Replacement',
              'multRepl', 'mult_repl', 'multi_replace', 'multiplicative_replacement',
              'mult replacement'))
    expect_identical(mz(x), 'mult', info = x)

  info <- match_zero_replacement('mult')
  expect_identical(info$func, replace_zeros_mult)
  expect_identical(info$params, c('counts', 'margin', 'delta'))
  expect_identical(info$name, 'Multiplicative Replacement')
  expect_null(info$alt_ids)

  df <- list_zero_replacements()
  expect_true(is.data.frame(df))
  expect_identical(names(df), c('name', 'id'))
  ids <- c('mult', 'czm', 'gbm', 'sq', 'bl')
  expect_identical(df$id, ids)
  expect_identical(rownames(list_zero_replacements(nm = 'id')), ids)
  expect_identical(list_zero_replacements('id'), ids)
  expect_identical(list_zero_replacements('name')[[1]], 'Multiplicative Replacement')
  expect_identical(list_zero_replacements('id', nm = 'name')[1:2],
                   c('Multiplicative Replacement' = 'mult', 'Count Zero Multiplicative' = 'czm'))
  expect_identical(list_zero_replacements('func'), list(
    mult = replace_zeros_mult, czm = replace_zeros_czm, gbm = replace_zeros_gbm,
    sq = replace_zeros_sq, bl = replace_zeros_bl ))
  expect_identical(list_zero_replacements('list')$mult, match_zero_replacement('mult'))
  expect_identical(names(formals(list_zero_replacements)), names(formals(list_count_ords)))

  expect_error(replace_zeros(m, 'xyz'),  'does not match any option: bl, czm, gbm, mult, sq')
  expect_error(replace_zeros(m, 123),    'character')
  expect_error(replace_zeros(m, c('mult', 'mult')), 'length 1')
  expect_error(replace_zeros(m, '-'),    'alphabetic')
  expect_error(replace_zeros(m, NA_character_), 'NA')
  expect_error(list_zero_replacements(val = 'xyz'))

  # Every method follows the shared shape: replace_zeros_<id>(counts, margin,
  # ...), a dense proportion matrix with `method` first and `n_replaced`
  # last.
  for (id in list_zero_replacements('id')) {
    f <- get(paste0('replace_zeros_', id))
    expect_identical(names(formals(f))[1:2], c('counts', 'margin'), info = id)
    z <- replace_zeros(m, id)
    expect_identical(attr(z, 'method'), id, info = id)
    a <- names(attributes(z))
    expect_identical(a[c(1:3, length(a))], c('dim', 'dimnames', 'method', 'n_replaced'), info = id)
    expect_true(all(z > 0) && max(abs(rowSums(z) - 1)) < 1e-14, info = id)
  }


  # -----------------------------------------------------------------------
  # Names: the planned registry (item 15) is now the registry (item 55)
  # -----------------------------------------------------------------------

  # The planned registry: every code and alias selects its own row, and each
  # bare code matches only its own row, exactly or by prefix. Item 55 built
  # it as planned.
  planned <- data.frame(
    stringsAsFactors = FALSE,
    name    = c('Multiplicative Replacement', 'Count Zero Multiplicative',
                'Geometric Bayesian-Multiplicative', 'Square-Root Bayesian-Multiplicative',
                'Bayes-Laplace Bayesian-Multiplicative'),
    id      = c('mult', 'czm', 'gbm', 'sq', 'bl'),
    alt_ids = c('multiplicative multreplacement multirepl multireplace',
                'countzeromultiplicative', 'geometricbayesianmultiplicative cmultrepl',
                'squareroot sqrtbm', 'bayeslaplace laplace') )
  hay <- ecodive:::make_haystack(planned)
  for (i in seq_len(nrow(planned)))
    for (a in c(planned$id[[i]], strsplit(planned$alt_ids[[i]], ' ')[[1]]))
      expect_identical(unique(unname(hay[startsWith(names(hay), a)])), planned$id[[i]], info = a)
  expect_identical(ZERO_REPLACEMENTS <- ecodive:::ZERO_REPLACEMENTS, planned)
  for (i in seq_len(nrow(planned)))
    for (a in c(planned$id[[i]], strsplit(planned$alt_ids[[i]], ' ')[[1]]))
      expect_identical(match_zero_replacement(a)$id, planned$id[[i]], info = a)

  # The current registry's ids and aliases begin no name of the other
  # registries' entries exactly (where the other wrappers look).
  ours <- names(ecodive:::exact_names('ZERO_REPLACEMENTS'))
  for (reg in c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'DA_TESTS', 'FEATURE_ASSOCS'))
    expect_identical(intersect(c(ours, planned$id), names(ecodive:::exact_names(reg))), character(0), info = reg)


  # -----------------------------------------------------------------------
  # The shared input checks, for the methods to come
  # -----------------------------------------------------------------------

  chk <- function (counts, margin = 1L, ...) ecodive:::zero_replacement_input(...)
  info <- chk(m)
  expect_identical(info[c('D', 'n')], list(D = 20L, n = 30L))
  expect_identical(info$features, colnames(m))
  expect_identical(info$samples, rownames(m))
  expect_identical(info$sample_zeros, unname(as.integer(rowSums(m == 0))))
  expect_identical(info$feature_nonzeros, unname(as.integer(colSums(m > 0))))
  expect_identical(chk(t(m), 2L), info)

  mz <- m; mz[, c(2, 5)] <- 0
  expect_silent(chk(mz))
  expect_error(chk(mz, empty_features_ok = FALSE),
               "2 features that are zero in every sample \\('F2', 'F5'\\)")
  expect_warning(chk(m / 3, whole = TRUE), 'not whole numbers')
  expect_silent(chk(m, whole = TRUE))
  sp <- matrix(byrow = TRUE, nrow = 5, dimnames = list(letters[1:5], paste0('f', 1:5)), data = c(
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 1,
    1, 0, 0, 1, 1,
    1, 1, 0, 1, 1,
    1, 1, 1, 1, 1 ))
  expect_warning(chk(sp, mostly_zero = TRUE),
                 "1 sample \\('a'\\) and 1 feature \\('f3'\\) are at least 80% zeros, so most of their values will be imputed")
  expect_warning(chk(sp[c('b', 'c', 'd'), ], mostly_zero = TRUE), "^1 feature \\('f3'\\) is at least 80% zeros")
  expect_silent(chk(sp))

#})
