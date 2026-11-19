#test_that("lda", {

  # Inputs and reference results, written by data-raw/expected/lda.r.
  ref <- readRDS(file.path('expected', 'lda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  fit_case <- function (cs, cpus = 1, counts = cs$counts, ...)
    lda(counts, cs$groups, norm = cs$norm, pseudocount = cs$pseudocount,
        prior = cs$prior, cpus = cpus, ...)

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }

  # `x` (without attributes) with each column's sign matched to `y`'s.
  match_signs <- function (x, y) {
    x <- unclass(x)
    attributes(x) <- list(dim = dim(x))
    sweep(x, 2, sign(colSums(x * y)), '*')
  }



  # Same as MASS::lda() and scikit-learn's LinearDiscriminantAnalysis, up to
  # the sign of each axis (and, for scikit-learn, a constant factor):
  # stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      lbl <- cs$name
      n   <- nrow(cs$x)
      g   <- nlevels(cs$groups)
      ms  <- cs$expected$mass
      sk  <- cs$expected$sklearn
      r   <- length(ms$svd)
      ord <- fit_case(cs, n_axes = r)

      expect_identical(dim(ord), c(n, r), info = lbl)
      expect_identical(rownames(ord), rownames(cs$counts), info = lbl)
      expect_identical(colnames(ord), paste0('LD', seq_len(r)), info = lbl)

      # MASS: scores, scaling, eig = svd^2, proportion of trace, means, prior.
      expect_equal(match_signs(ord, ms$scores), ms$scores, tolerance = 1e-10, info = lbl)
      s <- sign(colSums(unclass(ord) * ms$scores))
      expect_equal(unname(sweep(attr(ord, 'scaling'), 2, s, '*')), ms$scaling, tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'eig'), ms$svd^2, tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'), ms$svd^2 / sum(ms$svd^2), tolerance = 1e-10, info = lbl)
      expect_equal(unname(attr(ord, 'means')), ms$means, tolerance = 1e-12, info = lbl)
      expect_equal(unname(attr(ord, 'prior')), ms$prior, tolerance = 1e-15, info = lbl)

      # scikit-learn: the same axes, larger by sqrt(n / (n - g)).
      f <- sqrt(n / (n - g))
      expect_equal(match_signs(ord, sk$scores) * f, sk$scores, tolerance = 1e-10, info = lbl)
      expect_equal(unname(sweep(attr(ord, 'scaling'), 2, s * f, '*')),
                   unname(sweep(sk$scalings, 2, sign(colSums(sk$scores * ms$scores)), '*')),
                   tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'), sk$evr[seq_len(r)], tolerance = 1e-10, info = lbl)
      expect_equal(unname(attr(ord, 'means')), sk$means, tolerance = 1e-12, info = lbl)

      # The scores are the table, minus prior %*% means, times scaling.
      xbar <- drop(attr(ord, 'prior') %*% attr(ord, 'means'))
      expect_equal(unname(xbar), sk$xbar, tolerance = 1e-12, info = lbl)
      expect_equal(matrix(unclass(ord), n), sweep(cs$x, 2, xbar) %*% attr(ord, 'scaling'),
                   tolerance = 1e-12, check.attributes = FALSE, info = lbl)

      # Leave-one-out classification (the core of capda()): MASS's CV = TRUE
      # for every number of leading columns.
      pr  <- unname(attr(ord, 'prior'))
      loo <- ecodive:::lda_loo(cs$x, cs$groups, pr, 1L)
      expect_identical(loo$class, ms$loo_class, info = lbl)
      expect_equal(loo$posterior, ms$loo_post, tolerance = 1e-10, info = lbl)
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

  big <- local({
    set.seed(3)
    n <- 700
    m <- matrix(rpois(n * 150, 3), n, 150, dimnames = list(paste0('S', 1:n), paste0('F', 1:150)))
    g <- factor(sample(letters[1:5], n, TRUE))
    m[g == 'b', 1:5] <- m[g == 'b', 1:5] + 2
    list(m = m, g = g)
  })
  for (nm in c('big_site', 'iris', 'big')) {
    if (nm == 'big') { x <- big$m; gr <- big$g } else { x <- case(nm)$counts; gr <- case(nm)$groups }
    one <- lda(x, gr, cpus = 1)
    lo1 <- ecodive:::lda_loo(x + 0, gr, as.vector(table(gr)) / length(gr), 1L)
    for (cpus in unique(c(2L, max_cpus))) {
      expect_identical(lda(x, gr, cpus = cpus), one, info = paste(nm, cpus))
      expect_identical(ecodive:::lda_loo(x + 0, gr, as.vector(table(gr)) / length(gr), cpus), lo1, info = paste(nm, cpus))
    }
  }



  # The result. ====

  cs  <- case('iris')
  ord <- fit_case(cs)
  lev <- levels(cs$groups)

  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'scaling', 'means', 'prior',
      'groups', 'class') )
  expect_identical(attr(ord, 'method'), 'lda')
  expect_identical(colnames(ord), c('LD1', 'LD2'))
  expect_identical(dimnames(attr(ord, 'scaling')), list(colnames(cs$counts), colnames(ord)))
  expect_identical(dimnames(attr(ord, 'means')),   list(lev, colnames(cs$counts)))
  expect_identical(names(attr(ord, 'prior')), lev)
  expect_equal(unname(attr(ord, 'prior')), rep(1/3, 3))
  expect_identical(unname(attr(ord, 'groups')), cs$groups)
  expect_identical(names(attr(ord, 'groups')), rownames(cs$counts))

  # Pooled within-group covariance of the scores is the identity; each axis
  # flipped so that its largest absolute score is positive.
  sc <- unclass(ord)[, ]
  wc <- crossprod(sc - apply(sc, 2, function (v) ave(v, cs$groups))) / (150 - 3)
  expect_equal(wc, diag(2), tolerance = 1e-12, check.attributes = FALSE)
  expect_true(all(apply(sc, 2, function (v) v[which.max(abs(v))] > 0)))

  # eig: the between- over within-group mean square of each axis's scores.
  bm <- apply(sc, 2, function (v) sum(table(cs$groups) * (tapply(v, cs$groups, mean) - mean(v))^2) / 2)
  expect_equal(unname(attr(ord, 'eig')), unname(bm), tolerance = 1e-10)

  print_out <- capture.output(print(ord))
  expect_true(any(grepl('lda: 150 samples, 2 axes', print_out, fixed = TRUE)))


  # Groups given as a column of `data`, a factor, characters, numbers, or
  # through count_ord(): all the same. A data frame is a table too.
  meta <- data.frame(row.names = rownames(cs$counts), sp = cs$groups, num = as.integer(cs$groups))
  expect_identical(lda(cs$counts, 'sp', meta, cpus = 1), ord)
  expect_identical(lda(cs$counts, meta$sp, cpus = 1), ord)
  expect_identical(lda(as.data.frame(cs$counts), meta$sp, cpus = 1), ord)
  expect_identical(count_ord(cs$counts, 'lda', 'sp', meta, cpus = 1), ord)
  expect_identical(count_ord(cs$counts, 'Linear Discriminant Analysis', 'sp', meta, cpus = 1), ord)
  by_chr <- lda(cs$counts, as.character(cs$groups), cpus = 1)
  expect_equal(unclass(by_chr)[, ], unclass(ord)[, ])
  by_num <- lda(cs$counts, 'num', meta, cpus = 1)
  expect_equal(unclass(by_num)[, ], unclass(ord)[, ])
  expect_identical(levels(attr(by_num, 'groups')), c('1', '2', '3'))

  # Metadata rows are matched by name; unused levels are dropped.
  shuffled <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
  expect_identical(lda(cs$counts, 'sp', shuffled, cpus = 1), ord)
  expect_identical(lda(cs$counts, factor(cs$groups, levels = c(lev, 'x')), cpus = 1), ord)

  # n_axes: fewer axes, the same values; more than there are, capped.
  o1 <- lda(cs$counts, cs$groups, n_axes = 1, cpus = 1)
  expect_identical(colnames(o1), 'LD1')
  expect_equal(unclass(o1)[, 1], unclass(ord)[, 1])
  expect_equal(attr(o1, 'prop'), attr(ord, 'prop')[1])
  expect_identical(ncol(lda(cs$counts, cs$groups, n_axes = 10, cpus = 1)), 2L)

  # Priors: named in any order, or unnamed in level order.
  pri <- case('iris_prior')$prior
  op  <- lda(cs$counts, cs$groups, prior = pri, cpus = 1)
  expect_identical(lda(cs$counts, cs$groups, prior = rev(pri), cpus = 1), op)
  expect_identical(lda(cs$counts, cs$groups, prior = unname(pri), cpus = 1), op)
  expect_identical(attr(op, 'prior'), pri)
  expect_false(isTRUE(all.equal(unclass(op)[, ], unclass(ord)[, ])))

  # Two groups: one axis, the default; prop is 1.
  o <- fit_case(case('iris_two'))
  expect_identical(dim(o), c(100L, 1L))
  expect_equal(attr(o, 'prop'), 1)

  # Four groups: three axes by default.
  expect_identical(ncol(fit_case(case('big_site'))), 3L)

  # Input that is already an ordination: lda() on PCoA axes (as capda()).
  pc <- pcoa(bray(case('big_site')$counts), n_axes = 4, cpus = 1)
  o  <- lda(pc, case('big_site')$groups, cpus = 1)
  expect_identical(dim(o), c(151L, 3L))
  expect_identical(rownames(attr(o, 'scaling')), paste0('PCoA', 1:4))



  # Edge cases. ====

  m <- matrix(c(1, 5, 2, 8, 3, 1, 4, 0, 2, 7, 6, 5, 3, 3, 9, 1, 2, 4), 6, 3,
              dimnames = list(LETTERS[1:6], c('x', 'y', 'z')))

  # A singleton group.
  o <- lda(m[, 1:2], c('a', 'a', 'b', 'b', 'b', 'c'), cpus = 1)
  expect_identical(dim(o), c(6L, 2L))
  expect_true(all(is.finite(unclass(o))))

  # One feature: one axis, the standardized feature.
  o <- lda(m[, 1, drop = FALSE], c('a', 'a', 'a', 'b', 'b', 'b'), cpus = 1)
  expect_identical(dim(o), c(6L, 1L))
  expect_equal(abs(unname(attr(o, 'scaling')[1, 1])),
               1 / sqrt(sum((m[, 1] - ave(m[, 1], rep(1:2, each = 3)))^2) / 4))

  # Groups with identical means: no discriminant axes, with a warning.
  e <- cbind(x = c(0, 2, 0, 2), y = c(0, 0, 2, 2))
  expect_warning(o <- lda(rbind(e, e + 0), c('a', 'a', 'a', 'a', 'b', 'b', 'b', 'b'), cpus = 1),
                 'no discriminant axes')
  expect_identical(dim(o), c(8L, 0L))

  # Singular: more features than samples minus groups; a constant feature;
  # collinear features (relative abundances); n <= g.
  expect_error(lda(cbind(m, w = 1:6), c('a', 'a', 'b', 'b', 'c', 'c'), cpus = 1), 'singular')
  expect_error(lda(cbind(m[, 1:2], k = 1), c('a', 'a', 'a', 'b', 'b', 'b'), cpus = 1), 'singular')
  expect_error(lda(cbind(m[, 1:2], k = rep(c(1, 2), each = 3)), c('a', 'a', 'a', 'b', 'b', 'b'), cpus = 1), 'singular')
  expect_error(lda(m, c('a', 'a', 'a', 'b', 'b', 'b'), norm = 'percent', cpus = 1), 'singular')
  expect_error(lda(m[1:2, ], c('a', 'b'), cpus = 1), 'more samples than groups')
  expect_error(lda(case('iris')$counts, case('iris')$groups, norm = 'clr', pseudocount = 1, cpus = 1),
               'pcoa\\(\\) or pca\\(\\)')

  # Unnamed input.
  o <- lda(unname(m[, 1:2]), c('a', 'a', 'a', 'b', 'b', 'b'), cpus = 1)
  expect_null(rownames(o))
  expect_null(rownames(attr(o, 'scaling')))
  expect_null(colnames(attr(o, 'means')))
  expect_null(names(attr(o, 'groups')))


  # lda_loo(): NA where a fit leaving a sample out is singular.
  x  <- cbind(a = c(1, 2, 3, 4, 5, 6, 7, 8), b = c(2, 1, 4, 3, 6, 5, 8, 7), c = 1)
  gr <- factor(rep(c('p', 'q'), each = 4))
  lo <- ecodive:::lda_loo(x, gr, c(0.5, 0.5), 1L)
  expect_true(all(!is.na(lo$class[, 1:2])))
  expect_true(all(is.na(lo$class[, 3])))
  expect_true(all(is.na(lo$posterior)))
  # A sample alone in its group: NA when left out (as MASS), the others not.
  gr <- factor(c('p', 'p', 'p', 'q', 'q', 'q', 'q', 'r'))
  lo <- ecodive:::lda_loo(x[, 1:2], gr, c(0.4, 0.4, 0.2), 1L)
  expect_identical(lo$class[8, ], c(NA_integer_, NA_integer_))
  expect_true(all(is.na(lo$posterior[8, ])))
  expect_true(all(!is.na(lo$class[1:7, ])))



  # Invalid input. ====

  g6 <- c('a', 'a', 'a', 'b', 'b', 'b')
  m2 <- m[, 1:2]
  expect_error(lda(m2),                                    '`groups` is required')
  expect_error(lda(m2, rep('a', 6)),                       'at least two different values')
  expect_error(lda(m2, c('a', 'b', NA, 'b', 'a', 'a')),    'must not contain NAs')
  expect_error(lda(m2, c('a', 'b')),
               '`groups` has 2 values, but `counts` has 6 samples.', fixed = TRUE)
  expect_error(lda(m2, g6, n_axes = 0),                    '`n_axes` must be a whole number')
  expect_error(lda(m2, g6, norm = 'nope'))
  expect_error(lda(m2, g6, cpus = 0),                      'cpus')
  expect_error(lda(m2, g6, prior = c(0.5, 0.6)),           'add up to 1')
  expect_error(lda(m2, g6, prior = c(1, 0)),               'positive')
  expect_error(lda(m2, g6, prior = 1),                     '2 positive numbers')
  expect_error(lda(m2, g6, prior = c(a = 0.5, c = 0.5)),   'names must be the groups')
  expect_error(lda(m2, g6, prior = 'a'),                   'positive numbers')
  bad <- m2; bad[1, 1] <- NA
  expect_error(lda(bad, g6, cpus = 1),                     'must not contain NA')

#})
