#test_that("plsda", {

  # Inputs and reference results, written by data-raw/expected/plsda.r.
  ref <- readRDS(file.path('expected', 'plsda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  fit_case <- function (cs, cpus = 1, counts = cs$counts, ...)
    plsda(counts, cs$groups, norm = cs$norm, pseudocount = cs$pseudocount,
          scale = cs$scale, n_axes = cs$k, cpus = cpus, ...)

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }

  # `x` with each column's sign matched to `y`'s.
  match_signs <- function (x, y) {
    x <- unclass(x)
    attributes(x) <- list(dim = dim(x))
    sweep(x, 2, sign(colSums(x * y)), '*')
  }

  # VIP straight from the formula (Wold et al. 2001; Mehmood et al. 2012),
  # for the first a components, a = 1..k.
  vip_formula <- function (ord) {
    w   <- attr(ord, 'weights')
    ssy <- colSums(unclass(ord)^2) * colSums(attr(ord, 'y_loadings')^2)
    p   <- nrow(w)
    vapply(seq_len(ncol(w)), function (a) {
      ws <- w[, seq_len(a), drop = FALSE]
      sqrt(p * colSums(t(ws^2) * ssy[seq_len(a)]) / sum(ssy[seq_len(a)]))
    }, numeric(p))
  }



  # Same as scikit-learn's PLSRegression and mixOmics::plsda(), up to the
  # sign of each component: stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      lbl <- cs$name
      ord <- fit_case(cs)
      k   <- cs$k
      sk  <- cs$expected$sklearn

      expect_identical(dim(ord), c(nrow(cs$counts), as.integer(k)), info = lbl)
      expect_identical(rownames(ord), rownames(cs$counts), info = lbl)
      expect_identical(colnames(ord), paste0('PLSDA', seq_len(k)), info = lbl)

      expect_equal(match_signs(ord, sk$scores), sk$scores, tolerance = 1e-10, info = lbl)
      s <- sign(colSums(unclass(ord) * sk$scores))
      expect_equal(unname(sweep(attr(ord, 'weights'),    2, s, '*')), unname(sk$weights),    tolerance = 1e-10, info = lbl)
      expect_equal(unname(sweep(attr(ord, 'loadings'),   2, s, '*')), unname(sk$loadings),   tolerance = 1e-10, info = lbl)
      expect_equal(unname(sweep(attr(ord, 'y_loadings'), 2, s, '*')), unname(sk$y_loadings), tolerance = 1e-10, info = lbl)

      # scikit-learn's defaults: NIPALS converges at once with two groups.
      skd <- cs$expected$sklearn_default
      if (!is.null(skd))
        expect_equal(match_signs(ord, skd$scores), skd$scores, tolerance = 1e-10, info = lbl)

      mx <- cs$expected$mixomics
      if (!is.null(mx)) {
        expect_equal(match_signs(ord, mx$scores), mx$scores, tolerance = 1e-10, info = lbl)
        expect_equal(match_signs(attr(ord, 'weights'), mx$weights), mx$weights, tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'prop'), mx$prop, tolerance = 1e-10, info = lbl)
        if (cs$scale) {
          expect_equal(unname(attr(ord, 'vip')), mx$vip, tolerance = 1e-10, info = lbl)
        } else {
          # mixOmics weights each component by its squared correlations with
          # the dummy columns, which is the same only when they are scaled.
          w   <- attr(ord, 'weights')
          y   <- outer(as.integer(cs$groups), seq_along(levels(cs$groups)), '==') + 0
          ssy <- colSums(stats::cor(y, unclass(ord)[, ])^2)
          mv  <- vapply(seq_len(k), function (a)
            sqrt(nrow(w) * colSums(t(w[, seq_len(a), drop = FALSE]^2) * ssy[seq_len(a)]) / sum(ssy[seq_len(a)])),
            numeric(nrow(w)) )
          expect_equal(unname(mv), mx$vip, tolerance = 1e-10, info = lbl)
          expect_false(isTRUE(all.equal(unname(attr(ord, 'vip')), mx$vip)), info = lbl)
        }
      }

      # eig and trace: prop * trace, and the variance of the normalized
      # (scaled) table.
      x <- cs$x
      if (cs$scale) {
        x <- x[, apply(x, 2, stats::sd) > 0, drop = FALSE]
        x <- scale(x)
      }
      expect_equal(attr(ord, 'trace'), sum(apply(x, 2, stats::var)), tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'eig'), attr(ord, 'prop') * attr(ord, 'trace'), tolerance = 1e-12, info = lbl)

      expect_equal(attr(ord, 'vip'), vip_formula(ord), tolerance = 1e-12, check.attributes = FALSE, info = lbl)
    }


    # Sparse input: the same results.
    cs  <- case('big_site')
    ord <- fit_case(cs)
    if (requireNamespace('Matrix', quietly = TRUE)) {
      expect_equal(fit_case(cs, counts = Matrix::Matrix(cs$counts, sparse = TRUE)), ord, tolerance = 1e-12)
      expect_equal(fit_case(cs, counts = methods::as(Matrix::Matrix(cs$counts, sparse = TRUE), 'TsparseMatrix')), ord, tolerance = 1e-12)
    }
    if (requireNamespace('slam', quietly = TRUE))
      expect_equal(fit_case(cs, counts = slam::as.simple_triplet_matrix(cs$counts)), ord, tolerance = 1e-12)
    expect_equal(fit_case(cs, counts = t(cs$counts), margin = 2L), ord, tolerance = 1e-12)
  }



  # Identical for any number of threads. ====

  for (nm in c('big_site', 'dune_clr')) {
    cs  <- case(nm)
    one <- fit_case(cs, cpus = 1)
    for (cpus in unique(c(2L, max_cpus)))
      expect_identical(fit_case(cs, cpus = cpus), one, info = paste(nm, cpus))
  }



  # The result. ====

  cs  <- case('dune_mgmt')
  ord <- fit_case(cs)
  lev <- levels(cs$groups)

  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings', 'weights',
      'y_loadings', 'vip', 'groups', 'centroids', 'class') )
  expect_identical(attr(ord, 'method'), 'plsda')
  expect_identical(dimnames(attr(ord, 'loadings')), list(colnames(cs$counts), colnames(ord)))
  expect_identical(dimnames(attr(ord, 'weights')),  list(colnames(cs$counts), colnames(ord)))
  expect_identical(dimnames(attr(ord, 'vip')),      list(colnames(cs$counts), colnames(ord)))
  expect_identical(dimnames(attr(ord, 'y_loadings')), list(lev, colnames(ord)))
  expect_identical(dimnames(attr(ord, 'centroids')),  list(lev, colnames(ord)))
  expect_identical(unname(attr(ord, 'groups')), cs$groups)
  expect_identical(names(attr(ord, 'groups')), rownames(cs$counts))

  # Unit weights; orthogonal scores; each score flipped so that its largest
  # absolute value is positive.
  expect_equal(unname(colSums(attr(ord, 'weights')^2)), rep(1, 3), tolerance = 1e-12)
  tt <- crossprod(unclass(ord)[, ])
  expect_equal(tt[upper.tri(tt)], rep(0, 3), tolerance = 1e-10)
  expect_true(all(apply(unclass(ord)[, ], 2, function (v) v[which.max(abs(v))] > 0)))

  # Centroids: the mean score of each group.
  expect_equal(
    attr(ord, 'centroids'),
    as.matrix(stats::aggregate(unclass(ord)[, ], list(cs$groups), mean)[, -1]),
    check.attributes = FALSE, tolerance = 1e-12 )

  # Mean squared VIP is 1.
  expect_equal(unname(colMeans(attr(ord, 'vip')^2)), rep(1, 3), tolerance = 1e-12)

  print_out <- capture.output(print(ord))
  expect_true(any(grepl('plsda: 20 samples, 3 axes; % of variance', print_out, fixed = TRUE)))


  # Groups given as a column of `data`, a factor, characters, numbers, or
  # through count_ord(): all the same.
  meta <- data.frame(row.names = rownames(cs$counts), mgmt = cs$groups, num = as.integer(cs$groups))
  expect_identical(plsda(cs$counts, 'mgmt', meta, n_axes = 3, cpus = 1), ord)
  expect_identical(plsda(cs$counts, meta$mgmt, n_axes = 3, cpus = 1), ord)
  expect_identical(count_ord(cs$counts, 'plsda', 'mgmt', meta, n_axes = 3, cpus = 1), ord)
  by_chr <- plsda(cs$counts, as.character(cs$groups), n_axes = 3, cpus = 1)
  expect_equal(unclass(by_chr)[, ], unclass(ord)[, ])
  by_num <- plsda(cs$counts, 'num', meta, n_axes = 3, cpus = 1)
  expect_equal(unclass(by_num)[, ], unclass(ord)[, ])
  expect_identical(levels(attr(by_num, 'groups')), c('1', '2', '3', '4'))

  # Metadata rows are matched by name.
  shuffled <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
  expect_identical(plsda(cs$counts, 'mgmt', shuffled, n_axes = 3, cpus = 1), ord)

  # Unused factor levels are dropped.
  extra <- factor(cs$groups, levels = c(lev, 'unused'))
  expect_identical(plsda(cs$counts, extra, n_axes = 3, cpus = 1), ord)

  # Default: two components.
  expect_identical(colnames(plsda(cs$counts, cs$groups, cpus = 1)), c('PLSDA1', 'PLSDA2'))
  expect_equal(unclass(plsda(cs$counts, cs$groups, cpus = 1))[, ], unclass(ord)[, 1:2])


  # Constant features with scale = TRUE: left out (weights, loadings and VIP
  # of 0), and the rest as without them.
  cs2 <- case('dune_const')
  ord2 <- fit_case(cs2)
  expect_identical(unname(attr(ord2, 'weights')[c('Const', 'Zero'), ]),  matrix(0, 2, 3))
  expect_identical(unname(attr(ord2, 'loadings')[c('Const', 'Zero'), ]), matrix(0, 2, 3))
  expect_identical(unname(attr(ord2, 'vip')[c('Const', 'Zero'), ]),      matrix(0, 2, 3))
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ], tolerance = 1e-12)
  expect_equal(attr(ord2, 'trace'), attr(ord, 'trace'), tolerance = 1e-12)
  # VIP counts every feature, so the others' rise by sqrt(32 / 30).
  expect_equal(attr(ord2, 'vip')[colnames(cs$counts), ], attr(ord, 'vip') * sqrt(32 / 30), tolerance = 1e-12)



  # Edge cases. ====

  m <- matrix(c(1, 5, 2, 8, 3, 1, 4, 0, 2, 7, 6, 5), 4, 3, dimnames = list(LETTERS[1:4], c('x', 'y', 'z')))

  # Two samples, two groups: one component.
  o <- plsda(m[1:2, ], c('a', 'b'), n_axes = 3, cpus = 1)
  expect_identical(dim(o), c(2L, 1L))
  expect_equal(attr(o, 'prop'), 1)
  expect_equal(attr(o, 'vip')[, 1], sqrt(3 * attr(o, 'weights')[, 1]^2))

  # A singleton group.
  o <- plsda(m, c('a', 'a', 'b', 'c'), n_axes = 2, cpus = 1)
  expect_identical(dim(o), c(4L, 2L))
  expect_true(all(is.finite(unclass(o))))

  # More components than the table can fit: fewer, with a warning.
  dup <- cbind(x = m[, 'x'], x2 = m[, 'x'] * 2)
  expect_warning(o <- plsda(dup, c('a', 'a', 'b', 'b'), n_axes = 3, cpus = 1), 'only 1 of the first 2 components')
  expect_identical(ncol(o), 1L)
  expect_equal(attr(o, 'prop'), 1)

  # n_axes is capped at n - 1 and at the number of features, silently.
  expect_identical(ncol(plsda(m, c('a', 'b', 'a', 'b'), n_axes = 10, cpus = 1)), 3L)
  expect_identical(ncol(plsda(m[, 1:2], c('a', 'b', 'a', 'b'), n_axes = 10, cpus = 1)), 2L)

  # One feature.
  o <- plsda(m[, 1, drop = FALSE], c('a', 'b', 'a', 'b'), cpus = 1)
  expect_identical(dim(o), c(4L, 1L))
  expect_equal(unname(attr(o, 'weights')[1, 1]), 1)

  # Every feature constant, with scale = TRUE: no components.
  expect_warning(o <- plsda(matrix(3, 4, 2), c('a', 'b', 'a', 'b'), cpus = 1), 'only 0 of the first 2')
  expect_identical(dim(o), c(4L, 0L))
  expect_identical(attr(o, 'trace'), 0)

  # scale = FALSE: the unscaled table and dummy matrix.
  o <- plsda(m, c('a', 'b', 'a', 'b'), scale = FALSE, n_axes = 2, cpus = 1)
  expect_equal(attr(o, 'trace'), sum(apply(m, 2, stats::var)))
  w1 <- crossprod(scale(m, scale = FALSE), scale(outer(c(1, 2, 1, 2), 1:2, '==') + 0, scale = FALSE))
  w1 <- svd(w1)$u[, 1]
  expect_equal(abs(unname(attr(o, 'weights')[, 1])), abs(w1), tolerance = 1e-12)

  # Unnamed input.
  o <- plsda(unname(m), c('a', 'b', 'a', 'b'), cpus = 1)
  expect_null(rownames(o))
  expect_null(rownames(attr(o, 'loadings')))
  expect_null(names(attr(o, 'groups')))



  # Invalid input. ====

  expect_error(plsda(m),                                   '`groups` is required')
  expect_error(plsda(m, NULL),                             '`groups` is required')
  expect_error(plsda(m, rep('a', 4)),                      'at least two different values')
  expect_error(plsda(m, c('a', 'b', NA, 'b')),             'must not contain NAs')
  expect_error(plsda(m, c('a', 'b')),                      'groups')
  expect_error(plsda(m, 'nope', data.frame(g = 1:4)),      'nope')
  expect_error(plsda(m, c('a', 'b', 'a', 'b'), scale = NA), '`scale` must be TRUE or FALSE')
  expect_error(plsda(m, c('a', 'b', 'a', 'b'), n_axes = 0), '`n_axes` must be a whole number')
  expect_error(plsda(m, c('a', 'b', 'a', 'b'), norm = 'nope'))
  expect_error(plsda(m, c('a', 'b', 'a', 'b'), cpus = 0), 'cpus')
  bad <- m; bad[1, 1] <- NA
  expect_error(plsda(bad, c('a', 'b', 'a', 'b'), cpus = 1), 'must not contain NA')

#})
