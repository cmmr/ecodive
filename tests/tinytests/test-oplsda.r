#test_that("oplsda", {

  # Inputs and reference results, written by data-raw/expected/oplsda.r.
  ref <- readRDS(file.path('expected', 'oplsda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  fit_case <- function (cs, cpus = 1, counts = cs$counts, ...)
    oplsda(counts, cs$groups, norm = cs$norm, pseudocount = cs$pseudocount,
           scale = cs$scale, n_orthogonal = cs$k, cpus = cpus, ...)

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }

  # `x` (a vector or matrix) with each column's sign matched to `y`'s.
  match_signs <- function (x, y) {
    x <- as.matrix(unclass(x))
    attributes(x) <- list(dim = dim(x))
    y <- as.matrix(y)
    sweep(x, 2, sign(colSums(x * y)), '*')
  }



  # Same as pyopls and ropls::opls(), up to the sign of each component:
  # stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      lbl <- cs$name
      ord <- fit_case(cs)
      k   <- cs$k
      py  <- cs$expected$pyopls

      expect_identical(dim(ord), c(nrow(cs$counts), as.integer(k + 1)), info = lbl)
      expect_identical(rownames(ord), rownames(cs$counts), info = lbl)
      expect_identical(colnames(ord), c('p1', paste0('o', seq_len(k))), info = lbl)

      sc <- unclass(ord)[, , drop = FALSE]
      w  <- attr(ord, 'weights')
      ld <- attr(ord, 'loadings')

      # pyopls's OPLS: the orthogonal components.
      expect_equal(match_signs(sc[, -1], py$t_ortho), py$t_ortho, tolerance = 1e-10, info = lbl)
      s <- sign(colSums(sc[, -1, drop = FALSE] * py$t_ortho))
      expect_equal(unname(sweep(w[, -1, drop = FALSE],  2, s, '*')), py$w_ortho, tolerance = 1e-10, info = lbl)
      expect_equal(unname(sweep(ld[, -1, drop = FALSE], 2, s, '*')), py$p_ortho, tolerance = 1e-10, info = lbl)

      # PLSRegression(1) on pyopls's filtered table: the predictive one.
      expect_equal(match_signs(sc[, 1], py$t_pred), py$t_pred, tolerance = 1e-10, info = lbl)
      s <- sign(sum(sc[, 1] * py$t_pred))
      expect_equal(unname(w[, 1]  * s), drop(py$w_pred), tolerance = 1e-10, info = lbl)
      expect_equal(unname(ld[, 1] * s), drop(py$p_pred), tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'y_loading') * s, py$q_pred, tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'r2y'), py$r2y, tolerance = 1e-10, info = lbl)

      # ropls::opls(): scores, loadings, weights, VIP, R2X and R2Y.
      rp <- cs$expected$ropls
      if (!is.null(rp)) {
        expect_equal(match_signs(ord, rp$scores), rp$scores, tolerance = 1e-10, info = lbl)
        expect_equal(match_signs(ld, rp$loadings), rp$loadings, tolerance = 1e-10, info = lbl)
        expect_equal(match_signs(w, rp$weights), rp$weights, tolerance = 1e-10, info = lbl)
        expect_equal(unname(attr(ord, 'vip')), rp$vip, tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'prop'), rp$r2x_exact, tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'r2y'), rp$r2y_exact, tolerance = 1e-10, info = lbl)
        # Its rounded tables.
        expect_equal(signif(attr(ord, 'prop'), 3), rp$r2x[seq_len(k + 1)], info = lbl)
        expect_equal(signif(attr(ord, 'r2x'), 3), rp$r2x_cum, info = lbl)
        expect_equal(signif(attr(ord, 'r2y'), 3), rp$r2y_cum, info = lbl)
      }

      # eig, trace and r2x: prop * trace, the variance of the normalized
      # (scaled) table, and sum(prop).
      x <- cs$x
      if (cs$scale) {
        x <- x[, apply(x, 2, stats::sd) > 0, drop = FALSE]
        x <- scale(x)
      }
      expect_equal(attr(ord, 'trace'), sum(apply(x, 2, stats::var)), tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'eig'), attr(ord, 'prop') * attr(ord, 'trace'), tolerance = 1e-12, info = lbl)
      expect_equal(attr(ord, 'r2x'), sum(attr(ord, 'prop')), info = lbl)

      # The scores are orthogonal, and unit-length weights.
      cp <- crossprod(sc)
      expect_equal(cp[upper.tri(cp)], rep(0, k * (k + 1) / 2), tolerance = 1e-10 * max(cp), info = lbl)
      expect_equal(unname(colSums(w^2)), rep(1, k + 1), tolerance = 1e-12, info = lbl)

      # VIP from the formula.
      expect_equal(attr(ord, 'vip'), sqrt(nrow(ld)) * abs(ld[, 1]) / sqrt(sum(ld[, 1]^2)), tolerance = 1e-12, info = lbl)
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


    # With no orthogonal components: a one-component plsda(), up to scaling
    # of the Y loading (plsda() scales each dummy column; here one column).
    cs <- case('dune_two')
    o0 <- oplsda(cs$counts, cs$groups, n_orthogonal = 0, cpus = 1)
    pl <- plsda(cs$counts, cs$groups, n_axes = 1, cpus = 1)
    expect_identical(colnames(o0), 'p1')
    expect_equal(unclass(o0)[, 1], unclass(pl)[, 1], tolerance = 1e-12)
    expect_equal(attr(o0, 'weights')[, 1], attr(pl, 'weights')[, 1], tolerance = 1e-12)
    expect_equal(attr(o0, 'loadings')[, 1], attr(pl, 'loadings')[, 1], tolerance = 1e-12)
    expect_equal(attr(o0, 'eig'), attr(pl, 'eig'), tolerance = 1e-12)

    # The predictive component's fitted values are those of a
    # one-component PLS on the filtered table.
    ord <- fit_case(cs)
    yc  <- scale(as.integer(cs$groups) == 2)[, 1]
    fit <- unclass(ord)[, 1] * attr(ord, 'y_loading')
    expect_equal(attr(ord, 'r2y'), 1 - sum((yc - fit)^2) / sum(yc^2), tolerance = 1e-12)
  }



  # Identical for any number of threads. ====

  for (nm in c('big_site', 'dune_clr')) {
    cs  <- case(nm)
    ord <- fit_case(cs, cpus = 1)
    for (cpus in 2:max_cpus)
      expect_identical(fit_case(cs, cpus = cpus), ord, info = paste(nm, cpus))
  }



  # count_ord() runs it. ====

  cs   <- case('dune_two')
  meta <- data.frame(mgmt = cs$groups, row.names = rownames(cs$counts))
  ord  <- oplsda(cs$counts, 'mgmt', meta, n_orthogonal = 2, cpus = 1)
  expect_identical(count_ord(cs$counts, 'oplsda', 'mgmt', meta, n_orthogonal = 2, cpus = 1), ord)
  expect_identical(count_ord(cs$counts, 'OPLS-DA', 'mgmt', meta, n_orthogonal = 2, cpus = 1), ord)
  expect_identical(count_ord(cs$counts, 'opls', 'mgmt', meta, n_orthogonal = 2, cpus = 1), ord)



  # Attributes, groups and centroids. ====

  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'r2x', 'r2y',
      'loadings', 'weights', 'y_loading', 'vip', 'groups', 'centroids', 'class') )
  expect_identical(attr(ord, 'method'), 'oplsda')
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(names(attr(ord, 'vip')), colnames(cs$counts))
  expect_identical(dimnames(attr(ord, 'loadings')), list(colnames(cs$counts), c('p1', 'o1', 'o2')))
  expect_identical(dimnames(attr(ord, 'weights')),  list(colnames(cs$counts), c('p1', 'o1', 'o2')))
  expect_identical(attr(ord, 'groups'), stats::setNames(cs$groups, rownames(cs$counts)))
  lev <- levels(cs$groups)
  cen <- attr(ord, 'centroids')
  expect_identical(dimnames(cen), list(lev, c('p1', 'o1', 'o2')))
  for (g in lev)
    expect_equal(cen[g, ], colMeans(unclass(ord)[cs$groups == g, , drop = FALSE]))
  # The groups separate on p1 alone: the orthogonal centroids are 0.
  expect_equal(unname(cen[, -1]), matrix(0, 2, 2), tolerance = 1e-10)
  # VIP: the mean of the squares is 1.
  expect_equal(mean(attr(ord, 'vip')^2), 1)

  # Groups from a column name, a character vector, or numbers.
  meta$chr <- as.character(cs$groups)
  meta$num <- as.integer(cs$groups)
  by_chr <- oplsda(cs$counts, 'chr', meta, n_orthogonal = 2, cpus = 1)
  expect_equal(unclass(by_chr)[, ], unclass(ord)[, ])
  by_num <- oplsda(cs$counts, 'num', meta, n_orthogonal = 2, cpus = 1)
  expect_equal(unclass(by_num)[, ], unclass(ord)[, ])
  expect_identical(levels(attr(by_num, 'groups')), c('1', '2'))

  # Swapping the groups gives the same components (signs follow the
  # scores), and flips the Y loading's sign only if the scores' rule does.
  swapped <- factor(cs$groups, levels = rev(lev))
  sw <- oplsda(cs$counts, swapped, n_orthogonal = 2, cpus = 1)
  expect_equal(unclass(sw)[, ], unclass(ord)[, ], tolerance = 1e-12)
  expect_equal(attr(sw, 'y_loading'), -attr(ord, 'y_loading'), tolerance = 1e-12)
  expect_equal(attr(sw, 'r2y'), attr(ord, 'r2y'), tolerance = 1e-12)

  # Metadata rows are matched by name.
  shuffled <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
  expect_identical(oplsda(cs$counts, 'mgmt', shuffled, n_orthogonal = 2, cpus = 1), ord)

  # Unused factor levels are dropped.
  extra <- factor(cs$groups, levels = c(lev, 'unused'))
  expect_identical(oplsda(cs$counts, extra, n_orthogonal = 2, cpus = 1), ord)

  # Default: one orthogonal component.
  o1 <- oplsda(cs$counts, cs$groups, cpus = 1)
  expect_identical(colnames(o1), c('p1', 'o1'))
  expect_equal(unclass(o1)[, 'o1'], unclass(ord)[, 'o1'])


  # Constant features with scale = TRUE: left out (weights, loadings and VIP
  # of 0), and the rest as without them.
  cs2  <- case('dune_const')
  ord2 <- fit_case(cs2)
  expect_identical(unname(attr(ord2, 'weights')[c('Const', 'Zero'), ]),  matrix(0, 2, 3))
  expect_identical(unname(attr(ord2, 'loadings')[c('Const', 'Zero'), ]), matrix(0, 2, 3))
  expect_identical(unname(attr(ord2, 'vip')[c('Const', 'Zero')]),        c(0, 0))
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ], tolerance = 1e-12)
  expect_equal(attr(ord2, 'trace'), attr(ord, 'trace'), tolerance = 1e-12)
  expect_equal(attr(ord2, 'r2y'), attr(ord, 'r2y'), tolerance = 1e-12)
  # VIP counts every feature, so the others' rise by sqrt(32 / 30).
  expect_equal(attr(ord2, 'vip')[colnames(cs$counts)], attr(ord, 'vip') * sqrt(32 / 30), tolerance = 1e-12)



  # Edge cases. ====

  m <- matrix(c(1, 5, 2, 8, 3, 1, 4, 0, 2, 7, 6, 5), 4, 3, dimnames = list(LETTERS[1:4], c('x', 'y', 'z')))

  # Two samples: the predictive component only.
  o <- oplsda(m[1:2, ], c('a', 'b'), n_orthogonal = 3, cpus = 1)
  expect_identical(dim(o), c(2L, 1L))
  expect_equal(attr(o, 'prop'), 1)
  expect_equal(attr(o, 'r2y'), 1)

  # A singleton group.
  o <- oplsda(m, c('a', 'a', 'a', 'b'), cpus = 1)
  expect_identical(dim(o), c(4L, 2L))
  expect_true(all(is.finite(unclass(o))))

  # n_orthogonal is capped at n - 2 and at the number of features less one,
  # silently.
  expect_identical(ncol(oplsda(m, c('a', 'b', 'a', 'b'), n_orthogonal = 10, cpus = 1)), 3L)
  expect_identical(ncol(oplsda(m[, 1:2], c('a', 'b', 'a', 'b'), n_orthogonal = 10, cpus = 1)), 2L)
  expect_identical(ncol(oplsda(m[1:3, ], c('a', 'b', 'a'), n_orthogonal = 10, cpus = 1)), 2L)

  # No variation orthogonal to the groups: fewer, with a warning.
  dup <- cbind(x = m[, 'x'], x2 = m[, 'x'] * 2)
  expect_warning(o <- oplsda(dup, c('a', 'a', 'b', 'b'), n_orthogonal = 1, cpus = 1),
                 'only 0 of the first 1 orthogonal components')
  expect_identical(colnames(o), 'p1')
  expect_equal(attr(o, 'prop'), 1)

  # One feature.
  o <- oplsda(m[, 1, drop = FALSE], c('a', 'b', 'a', 'b'), cpus = 1)
  expect_identical(dim(o), c(4L, 1L))
  expect_equal(abs(unname(attr(o, 'weights')[1, 1])), 1)
  expect_equal(unname(attr(o, 'vip')), 1)

  # Every feature constant, or no covariance with the groups: an error.
  expect_error(oplsda(matrix(3, 4, 2), c('a', 'b', 'a', 'b'), cpus = 1), 'no covariance with `groups`')
  flat <- cbind(u = c(1, 1, 2, 2), v = c(5, 5, 9, 9))
  expect_error(oplsda(flat, c('a', 'b', 'b', 'a'), cpus = 1), 'no covariance with `groups`')

  # scale = FALSE: the unscaled table and indicator.
  o <- oplsda(m, c('a', 'b', 'a', 'b'), scale = FALSE, cpus = 1)
  expect_equal(attr(o, 'trace'), sum(apply(m, 2, stats::var)))
  w1 <- crossprod(scale(m, scale = FALSE), c(-1, 1, -1, 1))[, 1]
  expect_equal(abs(unname(attr(o, 'weights')[, 1])), unname(abs(w1)) / sqrt(sum(w1^2)), tolerance = 1e-12)
  yc <- c(-1, 1, -1, 1) / 2
  expect_equal(attr(o, 'y_loading')^2 * sum(unclass(o)[, 1]^2), attr(o, 'r2y') * sum(yc^2), tolerance = 1e-12)

  # Unnamed input.
  o <- oplsda(unname(m), c('a', 'b', 'a', 'b'), cpus = 1)
  expect_null(rownames(o))
  expect_null(rownames(attr(o, 'loadings')))
  expect_null(names(attr(o, 'vip')))
  expect_null(names(attr(o, 'groups')))



  # Invalid input. ====

  expect_error(oplsda(m),                                   '`groups` is required')
  expect_error(oplsda(m, NULL),                             '`groups` is required')
  expect_error(oplsda(m, rep('a', 4)),                      'must have two different values')
  expect_error(oplsda(m, c('a', 'b', 'c', 'a')),
               '`groups` has 3 different values, but oplsda() takes exactly two; use plsda() for more.', fixed = TRUE)
  expect_error(oplsda(m, c('a', 'b', NA, 'b')),             'must not contain NAs')
  expect_error(oplsda(m, c('a', 'b')),
               '`groups` has 2 values, but `counts` has 4 samples.', fixed = TRUE)
  expect_error(oplsda(m, 'g', data.frame(g = c('a', 'b', 'a'))),
               '`data` has 3 rows, but `counts` has 4 samples.', fixed = TRUE)
  expect_error(oplsda(m, 'nope', data.frame(g = 1:4)),      'nope')
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), scale = NA), '`scale` must be TRUE or FALSE')
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), n_orthogonal = -1), '`n_orthogonal` must be a whole number')
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), n_orthogonal = 1.5), '`n_orthogonal` must be a whole number')
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), n_orthogonal = NA), '`n_orthogonal` must be a whole number')
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), norm = 'nope'))
  expect_error(oplsda(m, c('a', 'b', 'a', 'b'), cpus = 0), 'cpus')
  bad <- m; bad[1, 1] <- NA
  expect_error(oplsda(bad, c('a', 'b', 'a', 'b'), cpus = 1), 'must not contain NA')

#})
