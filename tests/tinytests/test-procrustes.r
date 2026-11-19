#test_that("procrustes", {

  ex <- readRDS('expected/procrustes.rds')
  cf <- ex$configs
  X  <- cf$X
  n  <- nrow(X)

  target_of <- function (cs) {
    x <- cf[[cs$target]]
    if (!is.null(cs$tk)) x <- x[, seq_len(cs$tk), drop = FALSE]
    x
  }



  # Same results as vegan::procrustes(): vegan's procrustes(X, Y) is
  # ecodive's procrustes(Y, X). ====

  if (at_home()) {
    for (cs in ex$cases) {

      res <- procrustes(
        cf[[cs$ord]], target_of(cs), axes = cs$axes,
        symmetric = cs$symmetric, dilate = cs$dilate, cpus = 1 )

      info <- paste(cs$ord, cs$target, cs$symmetric, cs$dilate)
      tol  <- 1e-10

      # Absolute differences, relative to the coordinates' size: an exact
      # fit has ss and residuals of 0, up to rounding.
      size <- max(abs(cs$X))
      expect_equal(unname(unclass(res)[, ]), cs$Yrot, tolerance = tol, info = info)
      expect_equal(unname(attr(res, 'target')), cs$X, tolerance = tol, info = info)
      expect_true(abs(attr(res, 'ss') - cs$ss) <= tol * max(cs$ss, size^2), info = info)
      expect_true(max(abs(attr(res, 'residuals') - cs$residuals)) <= tol * size, info = info)
      expect_true(max(abs(attr(res, 'translation') - cs$translation)) <= tol * max(1, abs(cs$translation)), info = info)
      expect_equal(attr(res, 'dilation'), cs$scale, tolerance = tol, info = info)

      # vegan's rotation has one row per axis of ord, before padding; the
      # zero-padded axes' rows of ecodive's are arbitrary.
      rot <- unname(attr(res, 'rotation'))[seq_len(nrow(cs$rotation)), , drop = FALSE]
      if (cs$ord != 'X')  # an exact match's rotation is any for a zero singular value
        expect_true(max(abs(rot - cs$rotation)) <= tol, info = info)
    }
  }



  # Same statistic as scipy.spatial.procrustes(): the symmetric m^2. ====

  if (at_home() && !is.null(ex$scipy)) {
    pairs <- list(YX = list('Y', 1:3, X), RX = list('R', NULL, X[, 1:2]), NX = list('N', NULL, X[, 1:2]))
    for (nm in names(pairs)) {
      p   <- pairs[[nm]]
      res <- procrustes(cf[[p[[1]]]], p[[3]], axes = p[[2]], symmetric = TRUE)
      expect_true(abs(attr(res, 'ss') - ex$scipy[[nm]]) <= 1e-10 * max(ex$scipy[[nm]], 1e-6), info = nm)
    }
  }



  # Same statistics as vegan::protest(), given the same permutations:
  # every permuted r, and so the p-value. ====

  if (at_home()) {
    for (cs in ex$tests) {
      res  <- procrustes_test(
        cf[[cs$ord]], target_of(cs), axes = cs$axes,
        permutations = ex$perms, cpus = 1 )
      info <- paste(cs$ord, cs$target)
      expect_equal(res$statistic[['r']], cs$t0, tolerance = 1e-10, info = info)
      expect_equal(res$ss, cs$ss, tolerance = 1e-10, info = info)
      expect_identical(res$p_value, cs$signif, info = info)
      expect_equal(attr(res, 'r_perm'), cs$t, tolerance = 1e-10, info = info)
    }
  }



  # Result format: procrustes(). ====

  Y   <- cf$Y
  res <- procrustes(Y, X, cpus = 1)

  expect_inherits(res, 'ecodive_ord')
  expect_identical(dim(res), c(n, 2L))
  expect_identical(rownames(res), rownames(X))
  expect_identical(colnames(res), c('Dim1', 'Dim2'))
  expect_identical(attr(res, 'method'), 'procrustes')
  expect_identical(names(attributes(no_steps(res))), c(
    'dim', 'dimnames', 'method', 'eig', 'prop', 'ss', 'rotation', 'translation',
    'dilation', 'residuals', 'symmetric', 'target', 'class' ))
  expect_identical(attr(res, 'eig'), c(NA_real_, NA_real_))
  expect_identical(dim(attr(res, 'rotation')), c(2L, 2L))
  expect_identical(names(attr(res, 'residuals')), rownames(X))
  expect_false(attr(res, 'symmetric'))
  expect_equal(crossprod(attr(res, 'rotation')), diag(2), check.attributes = FALSE)
  expect_equal(attr(res, 'ss'), sum(attr(res, 'residuals')^2))
  expect_stdout(print(res), 'procrustes: 24 samples, 2 axes')

  # The fitted coordinates: dilation * ord %*% rotation + translation fits
  # target, with the same residuals.
  fitted <- attr(res, 'dilation') * Y %*% attr(res, 'rotation')
  fitted <- sweep(fitted, 2, attr(res, 'translation'), '+')
  expect_equal(sqrt(rowSums((X[, 1:2] - fitted)^2)), attr(res, 'residuals'), check.attributes = FALSE)

  # No dilation.
  expect_identical(attr(procrustes(Y, X, dilate = FALSE), 'dilation'), 1)

  # Symmetric: between 0 and 1, the same either way round, and zero
  # translation.
  s1 <- procrustes(Y, X, symmetric = TRUE)
  s2 <- procrustes(X, Y, symmetric = TRUE)
  expect_equal(attr(s1, 'ss'), attr(s2, 'ss'))
  expect_true(attr(s1, 'ss') > 0 && attr(s1, 'ss') < 1)
  expect_identical(unname(attr(s1, 'translation')), c(0, 0))

  # Column names come from target (or ord, for target's padded axes).
  o <- Y; colnames(o) <- c('NMDS1', 'NMDS2')
  t <- X; colnames(t) <- c('PCoA1', 'PCoA2', 'PCoA3')
  expect_identical(colnames(procrustes(o, t)), c('PCoA1', 'PCoA2'))
  expect_identical(colnames(procrustes(t, o, axes = 1:3)), c('NMDS1', 'NMDS2', 'PCoA3'))

  # A data frame works as a matrix.
  expect_identical(procrustes(as.data.frame(Y), X), res)



  # Exact matches: ss = 0, after any rotation, reflection, dilation, and
  # translation. ====

  expect_true(attr(procrustes(X, X), 'ss') < 1e-20)
  expect_true(attr(procrustes(cf$R, X[, 1:2]), 'ss') < 1e-20)
  expect_true(attr(procrustes(cf$R, X[, 1:2], dilate = FALSE), 'ss') > 0.1)
  expect_equal(attr(procrustes(cf$R, X[, 1:2]), 'dilation'), 1 / 2.5)
  expect_true(procrustes_test(X, X, permutations = 9)$statistic[['r']] > 1 - 1e-12)



  # Samples are matched by name; mismatched samples are an error. ====

  shuffled <- Y[rev(seq_len(n)), ]
  expect_equal(procrustes(Y, X[rev(seq_len(n)), ]), res)
  coords <- function (o) { m <- unclass(o); attributes(m) <- list(dim = dim(o)); m }
  expect_equal(coords(procrustes(shuffled, X)), coords(res)[rev(seq_len(n)), ])
  expect_equal(attr(procrustes(shuffled, X), 'ss'), attr(res, 'ss'))
  expect_equal(coords(procrustes(unname(Y), unname(X))), coords(res))
  expect_error(procrustes(Y[-1, ], X), 'no samples named "18"')
  expect_error(procrustes(Y, X[-1, ]), '`target` has no samples named "18"')
  expect_error(procrustes(unname(Y)[-1, ], X), '`ord` has 23 samples, but `target` has 24')
  y2 <- Y; rownames(y2)[2] <- rownames(y2)[1]
  expect_error(procrustes(y2, X), 'duplicate sample names')



  # dist input: the same as running pcoa() first. ====

  d1 <- dist(X)
  d2 <- dist(cf$N)
  expect_equal(no_steps(procrustes(d1, d1)), no_steps(procrustes(pcoa(d1), pcoa(d1))))
  expect_equal(no_steps(procrustes(d1, d2, n_axes = 2)), no_steps(procrustes(pcoa(d1, n_axes = 2), pcoa(d2, n_axes = 2))))
  expect_equal(no_steps(procrustes(d1, X)), no_steps(procrustes(pcoa(d1), X)))
  expect_equal(
    no_steps(procrustes_test(d1, d2, n_axes = 2, permutations = 19)),
    no_steps(procrustes_test(pcoa(d1, n_axes = 2), pcoa(d2, n_axes = 2), permutations = 19)) )



  # Result format: procrustes_test(). ====

  pt <- procrustes_test(Y, X, permutations = 99, cpus = 1)

  expect_identical(class(pt), c('ecodive_procrustes_test', 'ecodive_test'))
  expect_identical(names(pt), c(
    'statistic', 'p_value', 'z', 'n_perm', 'ss', 'n', 'n_axes', 'procrustes', 'permutation' ))
  expect_identical(names(pt$statistic), 'r')
  expect_identical(pt$n, n)
  expect_identical(pt$n_axes, 2L)
  expect_identical(pt$n_perm, 99L)
  expect_identical(pt$permutation, 'free')
  expect_identical(length(attr(pt, 'r_perm')), 99L)
  expect_equal(pt$statistic[['r']], sqrt(1 - pt$ss))
  expect_identical(no_steps(pt$procrustes), no_steps(procrustes(Y, X, symmetric = TRUE)))
  expect_equal(pt$z, (pt$statistic[['r']] - mean(attr(pt, 'r_perm'))) / sd(attr(pt, 'r_perm')))
  expect_identical(pt$p_value, 0.01)
  expect_stdout(print(pt), 'PROTEST')
  expect_stdout(print(pt), 'Samples: 24, axes: 2')
  row <- as.data.frame(pt)
  expect_identical(row$test, 'procrustes_test')
  expect_identical(row$stat_name, 'r')

  # The same whichever is rotated (symmetric statistic), given the same
  # permutations of the samples.
  expect_equal(
    procrustes_test(X, Y, permutations = ex$perms)$statistic,
    procrustes_test(Y, X, permutations = ex$perms)$statistic )

  # No permutations.
  pt0 <- procrustes_test(Y, X, permutations = 0)
  expect_identical(pt0$p_value, NA_real_)
  expect_identical(pt0$n_perm, 0L)
  expect_identical(attr(pt0, 'r_perm'), numeric(0))

  # Strata.
  strata <- rep(1:2, n / 2)
  ps <- procrustes_test(cf$N, X[, 1:2], permutations = 99, strata = strata)
  expect_identical(ps$permutation, 'within strata')
  expect_identical(procrustes_test(cf$N, X[, 1:2], permutations = ex$perms[1:5, ])$permutation, 'supplied matrix')
  expect_error(
    procrustes_test(Y, X, permutations = ex$perms, strata = strata),
    'Give `strata` or a `permutations` matrix, not both' )



  # Reproducible, and identical for any number of threads. ====

  a <- procrustes_test(cf$N, X[, 1:2], permutations = 199, seed = 3, cpus = 1)
  expect_identical(procrustes_test(cf$N, X[, 1:2], permutations = 199, seed = 3, cpus = 1), a)
  for (cpus in c(2, 4, 8))
    expect_identical(procrustes_test(cf$N, X[, 1:2], permutations = 199, seed = 3, cpus = cpus), a, info = cpus)
  expect_false(identical(
    attr(procrustes_test(cf$N, X[, 1:2], permutations = 199, seed = 4), 'r_perm'),
    attr(a, 'r_perm') ))
  expect_identical(
    procrustes_test(cf$N, X[, 1:2], permutations = 99, strata = strata, cpus = 1),
    procrustes_test(cf$N, X[, 1:2], permutations = 99, strata = strata, cpus = 4) )
  expect_identical(procrustes(d1, d2, n_axes = 2, cpus = 1), procrustes(d1, d2, n_axes = 2, cpus = 4))



  # Edge cases. ====

  # Two and three samples.
  two <- matrix(c(0, 1, 0, 0), 2)
  expect_true(attr(procrustes(two, two * 3), 'ss') < 1e-20)
  three <- matrix(c(0, 1, 0, 0, 0, 1), 3)
  p3 <- procrustes_test(three, three[c(2, 1, 3), ], permutations = 9)
  expect_true(is.finite(p3$statistic) && p3$p_value > 0 && p3$p_value <= 1)

  # One axis: a reflection at most.
  one <- matrix(c(1, 4, 2, 8), 4)
  expect_true(attr(procrustes(-one, one), 'ss') < 1e-20)
  expect_equal(unname(attr(procrustes(-one, one), 'rotation')), matrix(-1))

  # Duplicate points (ties) are fine.
  dup <- rbind(X[1:5, 1:2], X[1:5, 1:2])
  rownames(dup) <- NULL
  expect_true(is.finite(procrustes_test(dup, dup[, 2:1], permutations = 9)$statistic))

  # No variation.
  flat <- matrix(1, 4, 2)
  expect_error(procrustes(flat, matrix(1:8, 4)), '`ord` has no variation')
  expect_error(procrustes(matrix(1:8, 4), flat), '`target` has no variation')
  expect_error(procrustes(flat[1, , drop = FALSE], flat[1, , drop = FALSE]), 'at least 2 samples')



  # Invalid arguments. ====

  expect_error(procrustes('a', X), '`ord` must be a numeric matrix')
  expect_error(procrustes(Y, list(1)), '`target` must be a numeric matrix')
  expect_error(procrustes(data.frame(a = letters[1:3]), X), '`ord` must be a numeric matrix or data frame')
  y_na <- Y; y_na[1, 1] <- NA
  expect_error(procrustes(y_na, X), '`ord` must not contain NA')
  expect_error(procrustes(Y, X, axes = 4), '`axes` must be distinct whole numbers between 1 and 3')
  expect_error(procrustes(Y, X, axes = c(1, 1)), '`axes` must be distinct')
  expect_error(procrustes(Y, X, axes = 'Dim1'), '`axes` must be distinct')
  expect_error(procrustes(Y, X, symmetric = NA), '`symmetric` must be TRUE or FALSE')
  expect_error(procrustes(Y, X, dilate = 'yes'), '`dilate` must be TRUE or FALSE')
  expect_error(procrustes(Y, X, cpus = 0))
  expect_error(procrustes_test(Y, X, seed = 1.5))
  expect_error(procrustes_test(Y, X, permutations = -1))
  expect_error(procrustes_test(Y, X, strata = 1:3), 'strata')
  bad_d <- d1; bad_d[1] <- NA
  expect_error(procrustes(bad_d, d2), '`ord` must not contain NA')

#})
