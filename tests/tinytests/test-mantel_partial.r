#test_that("mantel_partial", {

  ex <- readRDS('expected/mantel.rds')
  x  <- ex$x
  n  <- attr(x, 'Size')



  # Same statistics as vegan::mantel.partial(), given the same permutations:
  # stored results, so at home only. Every permuted statistic, and so every
  # p-value, matches. (test-mantel.r has mantel()'s.) ====

  if (at_home()) {
    for (cs in ex$cases) {
      if (is.null(cs$z)) next

      res <- mantel_partial(
        x, ex[[cs$y]], given = ex[[cs$z]], method = cs$method,
        permutations = ex$perms, cpus = 1 )

      info <- paste(cs$method, cs$y, cs$z)
      expect_equal(res$statistic[["r"]], cs$statistic, tolerance = 1e-12, info = info)
      expect_identical(res$p_value, cs$signif, info = info)
      expect_equal(attr(res, 'r_perm'), cs$perm, tolerance = 1e-12, info = info)
    }
  }



  # Result format. ====

  res <- mantel_partial(x, ex$y, ex$z, permutations = 99, cpus = 1)

  expect_identical(class(res), c('ecodive_mantel', 'ecodive_test'))
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm',
    'method', 'alternative', 'partial', 'n', 'permutation' ))
  expect_identical(names(res$statistic), 'r')
  expect_identical(res$method, 'pearson')
  expect_identical(res$alternative, 'greater')
  expect_true(res$partial)
  expect_identical(res$n, n)
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  expect_identical(attr(res, 'display')$title, 'Partial Mantel test (Pearson)')
  expect_identical(length(attr(res, 'r_perm')), 99L)
  expect_true(res$p_value >= 1 / 100 && res$p_value <= 1)

  r_perm <- attr(res, 'r_perm')
  expect_equal(res$z, (res$statistic[["r"]] - mean(r_perm)) / sd(r_perm))

  res0 <- mantel_partial(x, ex$y, ex$z, permutations = 0)
  expect_identical(res0$p_value, NA_real_)
  expect_identical(res0$z, NA_real_)
  expect_identical(res0$n_perm, 0L)
  expect_identical(attr(res0, 'r_perm'), numeric(0))

  res <- mantel_partial(x, ex$y, given = ex$z, method = 'kendall', permutations = 9,
                        strata = rep(1:2, n / 2), alternative = 'less')
  expect_true(res$partial)
  expect_identical(res$permutation, 'within strata')
  expect_identical(res$alternative, 'less')
  expect_identical(attr(res, 'display')$title, 'Partial Mantel test (Kendall)')
  expect_identical(
    mantel_partial(x, ex$y, ex$z, permutations = ex$perms[1:3, ])$permutation,
    'supplied matrix' )

  details <- test_details.ecodive_mantel(res, 4)
  expect_true(is.character(details) && length(details) == 3)
  expect_true(any(grepl("Kendall's", details)))

  # The shared printout, and a data frame row named for the function (and
  # its dist_test() id).
  expect_stdout(print(res), 'Partial Mantel test \\(Kendall\\)')
  expect_stdout(print(res), "Kendall's rank correlation tau-b")
  row <- as.data.frame(res)
  expect_identical(row$test, 'mantel_partial')
  expect_identical(row$stat_name, 'r')
  expect_identical(row$statistic, res$statistic[['r']])
  expect_identical(row$p_value, res$p_value)
  expect_identical(names(row), c('test', 'stat_name', 'statistic', 'p_value', 'z', 'n_perm'))
  expect_identical(as.data.frame(mantel(x, ex$y, permutations = 9))$test, 'mantel')



  # The statistic is the partial correlation of stats::cor() of the
  # distances, as mantel.partial() computes it, for every method, with and
  # without ties (dist(age) is heavily tied). ====

  for (method in c('pearson', 'spearman', 'kendall')) {
    r <- function (a, b) stats::cor(as.vector(a), as.vector(b), method = method)
    for (yz in list(list(ex$y, ex$z), list(ex$a, ex$y), list(ex$y, ex$a))) {
      y  <- yz[[1]]
      z  <- yz[[2]]
      xy <- r(x, y); xz <- r(x, z); yz <- r(y, z)
      expect_equal(
        mantel_partial(x, y, z, method = method, permutations = 0)$statistic[["r"]],
        (xy - xz * yz) / sqrt(1 - xz^2) / sqrt(1 - yz^2),
        tolerance = 1e-12, info = method )
    }
  }

  # Each permutation relabels dm's samples, with y and given held fixed, as
  # in vegan.
  set.seed(3)
  perms <- t(replicate(5, sample(n)))
  dmx   <- as.matrix(x)
  for (method in c('pearson', 'spearman', 'kendall')) {
    r   <- function (a, b) stats::cor(as.vector(a), as.vector(b), method = method)
    yz  <- r(ex$a, ex$z)
    res <- mantel_partial(x, ex$a, ex$z, method = method, permutations = perms, cpus = 1)
    ref <- apply(perms, 1, function (p) {
      xp <- stats::as.dist(dmx[p, p])
      xy <- r(xp, ex$a); xz <- r(xp, ex$z)
      (xy - xz * yz) / sqrt(1 - xz^2) / sqrt(1 - yz^2)
    })
    expect_equal(attr(res, 'r_perm'), ref, tolerance = 1e-12, info = method)
  }



  # P-values for each alternative. ====

  eps <- sqrt(.Machine$double.eps)
  for (alt in c('greater', 'less', 'two.sided')) {
    res <- mantel_partial(x, ex$y, given = ex$z, permutations = 49, alternative = alt, cpus = 1)
    rp  <- attr(res, 'r_perm')
    hit <- switch(alt,
      greater   = rp >= res$statistic[["r"]] - eps,
      less      = rp <= res$statistic[["r"]] + eps,
      two.sided = abs(rp) >= abs(res$statistic[["r"]]) - eps )
    expect_identical(res$p_value, (sum(hit) + 1) / 50, info = alt)
  }



  # Reproducible: the same seed gives the same permutations, for any number
  # of threads (more than 2 only at home). ====

  cpus   <- if (at_home()) c(1, 2, 3, 4) else c(1, 2)
  strata <- rep(1:3, length.out = n)

  for (method in c('pearson', 'spearman', 'kendall')) {
    res <- lapply(cpus, function (cp)
      mantel_partial(x, ex$a, ex$z, method = method, permutations = 199,
                     strata = strata, cpus = cp))
    for (r in res[-1]) expect_identical(r, res[[1]], info = method)
  }

  a <- mantel_partial(x, ex$y, ex$z, permutations = 29, seed = 1)
  b <- mantel_partial(x, ex$y, ex$z, permutations = 29, seed = 1)
  d <- mantel_partial(x, ex$y, ex$z, permutations = 29, seed = 2)
  expect_identical(a, b)
  expect_false(identical(attr(a, 'r_perm'), attr(d, 'r_perm')))



  # Inputs: matched by sample name; numeric vectors; matrices. ====

  labs <- attr(x, 'Labels')
  y    <- ex$y
  zm   <- as.matrix(ex$z)
  shf  <- sample(n)

  ref <- mantel_partial(x, y, ex$z, permutations = 19)
  expect_identical(mantel_partial(x, as.matrix(y), zm, permutations = 19), ref)
  expect_identical(mantel_partial(x, y, stats::as.dist(zm[shf, shf]), permutations = 19), ref)

  # A numeric vector: its absolute differences, as dist() computes them.
  age <- stats::setNames(seq_len(n) %% 7 + 0.5, labs)
  ref <- mantel_partial(x, y, stats::dist(age), permutations = 19)
  expect_identical(mantel_partial(x, y, age, permutations = 19), ref)
  expect_identical(mantel_partial(x, y, rev(age), permutations = 19), ref)
  expect_identical(mantel_partial(x, y, unname(age), permutations = 19), ref)

  # y and given perfectly correlated: no partial correlation.
  expect_warning(
    res <- mantel_partial(x, age, given = rev(age), permutations = 9),
    "perfectly correlated" )
  expect_identical(res$statistic[["r"]], NA_real_)
  expect_identical(res$p_value, NA_real_)

  # dm and given perfectly correlated: the same (it gave r = 6.6e-10).
  expect_warning(
    res <- mantel_partial(x, age, given = x, permutations = 9),
    "perfectly correlated" )
  expect_identical(res$statistic[["r"]], NA_real_)

  # given constant: no correlation with anything.
  expect_warning(
    res <- mantel_partial(x, y, given = stats::dist(rep(1, n)), permutations = 9),
    'same distance between every pair' )
  expect_identical(res$statistic[["r"]], NA_real_)
  expect_identical(attr(res, 'r_perm'), rep(NA_real_, 9))



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

  ref <- mantel_partial(x, v('age'), v('ph'), permutations = 19)
  expect_identical(mantel_partial(x, 'age', 'ph', meta, permutations = 19), ref)
  expect_identical(mantel_partial(x, 'age', given = 'ph', data = meta, permutations = 19), ref)
  expect_identical(mantel_partial(x, 'age', meta$ph, meta, permutations = 19), ref)

  extra <- data.frame(row.names = 'X', age = 99, ph = 0, blk = 1, grp = 'a')
  shuffled <- rbind(meta[rev(seq_len(n)), ], extra)
  expect_identical(mantel_partial(x, 'age', 'ph', data = shuffled, permutations = 19), ref)

  expect_identical(
    mantel_partial(x, ex$y, 'ph', meta, strata = 'blk', permutations = 19),
    mantel_partial(x, ex$y, v('ph'), strata = v('blk'), permutations = 19) )
  expect_identical(
    mantel_partial(x, 'age', ex$z, meta, strata = meta$blk, permutations = 19),
    mantel_partial(x, v('age'), ex$z, strata = v('blk'), permutations = 19) )

  na_meta <- meta
  na_meta$ph[[3]] <- NA
  expect_error(mantel_partial(x, ex$y, 'ph'), '`given` names a column, but `data` is not given')
  expect_error(mantel_partial(x, ex$y, 'nope', meta), '`given` is not a column of `data`')
  expect_error(mantel_partial(x, ex$y, 'grp', meta), 'not numeric.*permanova')
  expect_error(mantel_partial(x, ex$y, 'ph', na_meta), 'Missing values')



  # dist_test(): `groups` is passed as `y`, and `given` through `...`. ====

  expect_identical(
    dist_test(x, 'mantel_partial', 'age', meta, given = 'ph', permutations = 19, cpus = 1),
    mantel_partial(x, 'age', given = 'ph', data = meta, permutations = 19, cpus = 1) )
  expect_identical(
    dist_test(x, 'mantel_partial', v('age'), given = ex$z, method = 'spearman',
              permutations = 19, cpus = 1),
    mantel_partial(x, v('age'), ex$z, method = 'spearman', permutations = 19, cpus = 1) )
  expect_error(dist_test(x, 'mantel_partial', 'grp', meta, given = 'ph'), "column 'grp' is categorical")
  expect_error(dist_test(x, 'mantel_partial', 'age', meta), '`given` is required')
  expect_identical(match_dist_test('mantel.partial')$id,  'mantel_partial')
  expect_identical(match_dist_test('mantel_partial')$id,  'mantel_partial')
  expect_identical(match_dist_test('partial Mantel')$id,  'mantel_partial')
  expect_identical(match_dist_test('Partial Mantel Test')$id, 'mantel_partial')
  expect_identical(match_dist_test('mantelpartial')$id,   'mantel_partial')
  expect_identical(match_dist_test('mantel')$id,          'mantel')
  expect_identical(match_dist_test('mantel_partial')$params, names(formals(mantel_partial)))
  expect_error(match_dist_test('mantel_partial', groups = 'categorical'), 'does not match any option')



  # Errors. ====

  expect_error(mantel_partial(x, y), '`given` is required.*use mantel\\(\\)')
  expect_error(mantel_partial(x, y, NULL), '`given` is required')
  expect_error(mantel_partial(x, given = ex$z), '`y` is required')
  expect_error(mantel_partial(x, y, 1:5), '`given` has 5 samples, but `dm` has 30')
  expect_error(mantel_partial(x, y, c(age[-1], Z = 1)), '`given` has no samples named "S01"')
  expect_error(mantel_partial(x, y, matrix('a', 2, 2)), '`given` must be')
  expect_error(mantel_partial(x, y, c(NA, age[-1])), '`given` must not contain NA')
  expect_error(mantel_partial(x, 1:5, ex$z), '`y` has 5 samples, but `dm` has 30')
  expect_error(mantel_partial(stats::dist(1:2), stats::dist(1:2), stats::dist(1:2)), 'at least 3 samples')
  expect_error(mantel_partial(x, y, ex$z, method = 'foo'))
  expect_error(mantel_partial(x, y, ex$z, alternative = 'foo'))
  expect_error(mantel_partial(x, y, ex$z, permutations = -1), '`permutations` must be')
  expect_error(mantel_partial(x, y, ex$z, permutations = perms, strata = strata), 'not both')
  expect_error(mantel_partial(x, y, ex$z, seed = 1.5), '`seed`')
  expect_error(mantel_partial(x, y, ex$z, cpus = 0), '`cpus`')

#})
