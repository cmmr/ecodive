#test_that("mantel", {

  ex <- readRDS('expected/mantel.rds')
  x  <- ex$x
  n  <- attr(x, 'Size')



  # Same statistics as vegan::mantel() and mantel.partial(), given the same
  # permutations: stored results, so at home only. Every permuted statistic,
  # and so every p-value, matches. ====

  if (at_home()) {
    for (cs in ex$cases) {

      res <- mantel(
        x, ex[[cs$y]], method = cs$method, permutations = ex$perms,
        given = if (!is.null(cs$z)) ex[[cs$z]], cpus = 1 )

      info <- paste(cs$method, cs$y, cs$z)
      expect_equal(res$statistic[["r"]], cs$statistic, tolerance = 1e-12, info = info)
      expect_identical(res$p_value, cs$signif, info = info)
      expect_equal(attr(res, 'r_perm'), cs$perm, tolerance = 1e-12, info = info)
    }
  }



  # Result format. ====

  res <- mantel(x, ex$y, permutations = 99, cpus = 1)

  expect_identical(class(res), c('ecodive_mantel', 'ecodive_test'))
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm',
    'method', 'alternative', 'partial', 'n', 'permutation' ))
  expect_identical(names(res$statistic), 'r')
  expect_identical(res$method, 'pearson')
  expect_identical(res$alternative, 'greater')
  expect_false(res$partial)
  expect_identical(res$n, n)
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  expect_identical(attr(res, 'title'), 'Mantel test (Pearson)')
  expect_identical(length(attr(res, 'r_perm')), 99L)
  expect_true(res$p_value >= 1 / 100 && res$p_value <= 1)

  r_perm <- attr(res, 'r_perm')
  expect_equal(res$z, (res$statistic[["r"]] - mean(r_perm)) / sd(r_perm))

  res0 <- mantel(x, ex$y, permutations = 0)
  expect_identical(res0$p_value, NA_real_)
  expect_identical(res0$z, NA_real_)
  expect_identical(res0$n_perm, 0L)
  expect_identical(attr(res0, 'r_perm'), numeric(0))

  res <- mantel(x, ex$y, given = ex$z, method = 'kendall', permutations = 9,
                strata = rep(1:2, n / 2), alternative = 'less')
  expect_true(res$partial)
  expect_identical(res$permutation, 'within strata')
  expect_identical(res$alternative, 'less')
  expect_identical(attr(res, 'title'), 'Partial Mantel test (Kendall)')
  expect_identical(
    mantel(x, ex$y, permutations = ex$perms[1:3, ])$permutation,
    'supplied matrix' )

  details <- test_details.ecodive_mantel(res, 4)
  expect_true(is.character(details) && length(details) == 3)
  expect_true(any(grepl("Kendall's", details)))

  # The shared printout and data frame row.
  expect_stdout(print(res), 'Partial Mantel test \\(Kendall\\)')
  expect_stdout(print(res), "Kendall's rank correlation tau-b")
  row <- as.data.frame(res)
  expect_identical(row$test, 'mantel')
  expect_identical(row$stat_name, 'r')
  expect_identical(row$statistic, res$statistic[['r']])
  expect_identical(row$p_value, res$p_value)



  # The statistic is stats::cor() of the distances, for every method, with
  # and without ties (dist(age) is heavily tied). ====

  for (method in c('pearson', 'spearman', 'kendall')) {
    for (y in list(ex$y, ex$a)) {
      expect_equal(
        mantel(x, y, method = method, permutations = 0)$statistic[["r"]],
        stats::cor(as.vector(x), as.vector(y), method = method),
        tolerance = 1e-12, info = method )
    }

    # Partial: as mantel.partial() computes it.
    r  <- function (a, b) stats::cor(as.vector(a), as.vector(b), method = method)
    xy <- r(x, ex$y); xz <- r(x, ex$z); yz <- r(ex$y, ex$z)
    expect_equal(
      mantel(x, ex$y, given = ex$z, method = method, permutations = 0)$statistic[["r"]],
      (xy - xz * yz) / sqrt(1 - xz^2) / sqrt(1 - yz^2),
      tolerance = 1e-12, info = method )
  }

  # Kendall, both counting paths: a target with more than 4096 distinct
  # values (merge sort), and an x with ties within its tie groups.
  set.seed(2)
  m  <- 100
  xb <- dist(round(rnorm(m), 1))      # many ties
  yb <- dist(rnorm(m))                # no ties
  zb <- dist(sample(1:4, m, TRUE))    # few values
  for (pair in list(list(xb, yb), list(yb, xb), list(xb, zb), list(yb, zb))) {
    expect_equal(
      mantel(pair[[1]], pair[[2]], method = 'kendall', permutations = 0)$statistic[["r"]],
      stats::cor(as.vector(pair[[1]]), as.vector(pair[[2]]), method = 'kendall'),
      tolerance = 1e-12 )
  }



  # Each permutation relabels dm's samples, as in vegan. ====

  set.seed(3)
  perms <- t(replicate(5, sample(n)))
  dmx   <- as.matrix(x)

  for (method in c('pearson', 'spearman', 'kendall')) {
    res <- mantel(x, ex$a, method = method, permutations = perms, cpus = 1)
    ref <- apply(perms, 1, function (p)
      stats::cor(as.vector(stats::as.dist(dmx[p, p])), as.vector(ex$a), method = method))
    expect_equal(attr(res, 'r_perm'), ref, tolerance = 1e-12, info = method)
  }

  # The identity gives back the observed statistic.
  res <- mantel(x, ex$y, permutations = matrix(seq_len(n), 2, n, byrow = TRUE))
  expect_equal(attr(res, 'r_perm'), rep(res$statistic[["r"]], 2), tolerance = 1e-14)
  expect_identical(res$p_value, 1)



  # P-values for each alternative. ====

  eps <- sqrt(.Machine$double.eps)
  for (alt in c('greater', 'less', 'two.sided')) {
    res <- mantel(x, ex$y, given = ex$z, permutations = 49, alternative = alt, cpus = 1)
    rp  <- attr(res, 'r_perm')
    hit <- switch(alt,
      greater   = rp >= res$statistic[["r"]] - eps,
      less      = rp <= res$statistic[["r"]] + eps,
      two.sided = abs(rp) >= abs(res$statistic[["r"]]) - eps )
    expect_identical(res$p_value, (sum(hit) + 1) / 50, info = alt)
  }



  # Reproducible: the same seed gives the same permutations, for any number
  # of threads (more than 2 only at home). ====

  cpus <- if (at_home()) c(1, 2, 3, 4) else c(1, 2)
  strata <- rep(1:3, length.out = n)

  for (method in c('pearson', 'spearman', 'kendall')) {
    res <- lapply(cpus, function (cp)
      mantel(x, ex$a, method = method, permutations = 199, strata = strata, cpus = cp))
    for (r in res[-1]) expect_identical(r, res[[1]], info = method)
  }

  a <- mantel(x, ex$y, permutations = 29, seed = 1)
  b <- mantel(x, ex$y, permutations = 29, seed = 1)
  d <- mantel(x, ex$y, permutations = 29, seed = 2)
  expect_identical(a, b)
  expect_false(identical(attr(a, 'r_perm'), attr(d, 'r_perm')))



  # Strata: samples are only shuffled within them. With y = 1 between
  # strata and 0 within, every within-strata permutation keeps the sum over
  # pairs, so every permuted statistic is the observed one. ====

  grp <- rep(c('A', 'B', 'C'), each = n / 3)
  ys  <- stats::as.dist(outer(grp, grp, "!=") * 1)
  res <- mantel(x, ys, permutations = 50, strata = grp)
  expect_equal(attr(res, 'r_perm'), rep(res$statistic[["r"]], 50), tolerance = 1e-12)

  free <- mantel(x, ys, permutations = 50)
  expect_true(stats::sd(attr(free, 'r_perm')) > 0)

  # Named strata are matched to the samples.
  labs <- attr(x, 'Labels')
  res2 <- mantel(x, ys, permutations = 50, strata = rev(stats::setNames(grp, labs)))
  expect_identical(res2, res)



  # Inputs: matched by sample name; numeric vectors; matrices. ====

  y   <- ex$y
  ym  <- as.matrix(y)
  shf <- sample(n)

  ref <- mantel(x, y, permutations = 19)
  expect_identical(mantel(x, ym, permutations = 19), ref)
  expect_identical(mantel(x, stats::as.dist(ym[shf, shf]), permutations = 19), ref)
  # dm's order is the one permuted: y follows it.
  expect_identical(
    mantel(as.matrix(x)[shf, shf], y, permutations = 19),
    mantel(as.matrix(x)[shf, shf], ym[shf, shf], permutations = 19) )

  # A numeric vector: its absolute differences, as dist() computes them.
  age <- stats::setNames(seq_len(n) %% 7 + 0.5, labs)
  ref <- mantel(x, stats::dist(age), permutations = 19)
  expect_identical(mantel(x, age, permutations = 19), ref)
  expect_identical(mantel(x, rev(age), permutations = 19), ref)
  expect_identical(mantel(x, unname(age), permutations = 19), ref)

  # y and given perfectly correlated: no partial correlation.
  expect_warning(
    res <- mantel(x, age, given = rev(age), permutations = 9),
    "perfectly correlated" )
  expect_identical(res$statistic[["r"]], NA_real_)

  # dm and given perfectly correlated: the same (it gave r = 6.6e-10).
  expect_warning(
    res <- mantel(x, age, given = x, permutations = 9),
    "perfectly correlated" )
  expect_identical(res$statistic[["r"]], NA_real_)



  # `data`: y, given and strata as its column names, with rows matched to
  # the samples by name (extra rows ignored). ====

  meta <- data.frame(
    row.names = labs,
    age = seq_len(n) %% 7 + 0.5,
    ph  = round(seq_len(n) / 3, 1) %% 2,
    blk = rep(1:3, length.out = n),
    grp = rep(c('a', 'b'), length.out = n),
    stringsAsFactors = FALSE )
  v <- function (col) stats::setNames(meta[[col]], labs)

  ref <- mantel(x, v('age'), permutations = 19)
  expect_identical(mantel(x, 'age', meta, permutations = 19), ref)

  extra <- data.frame(row.names = 'X', age = 99, ph = 0, blk = 1, grp = 'a')
  shuffled <- rbind(meta[rev(seq_len(n)), ], extra)
  expect_identical(mantel(x, 'age', data = shuffled, permutations = 19), ref)

  expect_identical(
    mantel(x, 'age', meta, given = 'ph', permutations = 19),
    mantel(x, v('age'), given = v('ph'), permutations = 19) )
  expect_identical(
    mantel(x, 'age', meta, strata = 'blk', permutations = 19),
    mantel(x, v('age'), strata = v('blk'), permutations = 19) )
  expect_identical(
    mantel(x, ex$y, meta, given = 'ph', strata = meta$blk, permutations = 19),
    mantel(x, ex$y, given = v('ph'), strata = v('blk'), permutations = 19) )

  # Without row names, rows are in dm's order.
  expect_identical(mantel(x, 'age', data = `rownames<-`(meta, NULL), permutations = 19), ref)

  na_meta <- meta
  na_meta$age[[3]] <- NA
  expect_error(mantel(x, 'age'), '`y` names a column, but `data` is not given')
  expect_error(mantel(x, ex$y, given = 'ph'), '`given` names a column, but `data` is not given')
  expect_error(mantel(x, 'nope', meta), '`y` is not a column of `data`')
  expect_error(mantel(x, 'grp', meta), 'not numeric.*permanova')
  expect_error(mantel(x, 'age', na_meta), 'Missing values')
  expect_error(mantel(x, 'age', meta[-1, ]), 'no row in `data`')
  expect_error(mantel(x, 'age', list(age = 1)), '`data` must be a data.frame')
  expect_error(mantel(x), '`y` is required')



  # Errors and warnings. ====

  expect_error(mantel(x, 1:5), '`y` has 5 samples, but `dm` has 30')
  expect_error(mantel(x, c(age[-1], Z = 1)), '`y` has no samples named "S01"')
  expect_error(mantel(x, matrix(1:6, 2, 3)), '`y` must be a dist object')
  expect_error(mantel(x, c(NA, age[-1])), '`y` must not contain NA')
  expect_error(mantel(x, y, given = matrix('a', 2, 2)), '`given` must be')
  expect_error(mantel(stats::dist(1:2), stats::dist(1:2)), 'at least 3 samples')
  expect_error(mantel(x, y, method = 'foo'))
  expect_error(mantel(x, y, alternative = 'foo'))
  expect_error(mantel(x, y, permutations = -1), '`permutations` must be')
  expect_error(mantel(x, y, permutations = perms, strata = grp), 'not both')
  expect_error(mantel(x, y, strata = grp[-1]), 'one value per sample')
  expect_error(mantel(x, y, strata = c(NA, grp[-1])), 'must not contain NA')
  expect_error(mantel(x, y, seed = 1.5), '`seed`')
  expect_error(mantel(x, y, cpus = 0), '`cpus`')

  expect_warning(
    res <- mantel(x, stats::dist(rep(1, n)), permutations = 9),
    'same distance between every pair' )
  expect_identical(res$statistic[["r"]], NA_real_)
  expect_identical(res$p_value, NA_real_)
  expect_identical(attr(res, 'r_perm'), rep(NA_real_, 9))

#})
