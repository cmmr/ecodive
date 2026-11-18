#test_that("mrpp", {

  ex   <- readRDS('expected/mrpp.rds')
  d    <- bray(ex$counts)
  meta <- ex$meta
  n    <- attr(d, 'Size')



  # Same results as vegan::mrpp() and summary(vegan::permustats()), given
  # the same permutations: stored results, so at home only. Every permuted
  # delta, and so the p-value and z, match. ====

  if (at_home()) {
    for (cs in ex$cases) {

      info <- paste(cs$groups, cs$weight_type)
      res  <- mrpp(
        d, cs$groups, meta, weights = cs$weight_type,
        permutations = ex$perms, cpus = 1 )

      expect_equal(res$statistic[['delta']],          cs$delta,   tolerance = 1e-12, info = info)
      expect_equal(res$expected_delta, cs$E_delta, tolerance = 1e-12, info = info)
      expect_equal(res$A,              cs$A,       tolerance = 1e-12, info = info)
      expect_equal(unname(res$class_delta), unname(cs$classdelta), tolerance = 1e-12, info = info)
      expect_identical(names(res$class_delta), names(cs$classdelta), info = info)
      expect_identical(res$class_n, cs$n, info = info)
      expect_equal(attr(res, 'delta_perm'), cs$boot_deltas, tolerance = 1e-12, info = info)
      expect_identical(res$p_value, cs$p, info = info)
      expect_equal(res$z, cs$z, tolerance = 1e-10, info = info)
    }

    # Free permutation within blocks: the same p-value up to sampling error.
    res <- mrpp(d, 'trt', meta, permutations = 9999, strata = 'block', cpus = 1)
    expect_true(abs(res$p_value - ex$strata_p) < 0.01)
  }



  # The same as a plain R implementation, given the same permutations. ====

  ref_delta <- function (g, wt) {
    m  <- as.matrix(d)
    gs <- unique(g)
    cd <- vapply(gs, function (x) {
      i <- which(g == x)
      if (length(i) < 2) NA_real_ else mean(m[i, i][lower.tri(diag(length(i)))]) }, 1)
    ng <- vapply(gs, function (x) sum(g == x), 1)
    w  <- switch(wt, ng, ng - 1, ng * (ng - 1) / 2)
    ok <- !is.na(cd)
    sum(cd[ok] * w[ok]) / sum(w[ok])
  }

  perms <- ex$perms[1:25, ]
  for (grp in c('trt', 'uneven')) {
    for (wt in 1:3) {
      g   <- meta[[grp]]
      res <- mrpp(d, g, weights = wt, permutations = perms, cpus = 1)
      dp  <- apply(perms, 1, function (p) ref_delta(g[p], wt))
      del <- ref_delta(g, wt)
      expect_equal(res$statistic[['delta']], del, tolerance = 1e-12)
      expect_equal(attr(res, 'delta_perm'), dp, tolerance = 1e-12)
      expect_equal(res$z, (del - mean(dp)) / sd(dp), tolerance = 1e-10)
      expect_identical(res$p_value, (1 + sum(dp <= del + sqrt(.Machine$double.eps))) / 26)
      expect_equal(res$expected_delta, mean(d), tolerance = 1e-14)
      expect_equal(res$A, 1 - del / mean(d), tolerance = 1e-12)
    }
  }

  # Group weights by name.
  expect_identical(mrpp(d, meta$uneven, weights = 'n-1', permutations = 9, cpus = 1)$weights, 'n-1')
  expect_identical(
    mrpp(d, meta$uneven, weights = 'n(n-1)', permutations = 9, cpus = 1),
    mrpp(d, meta$uneven, weights = 3, permutations = 9, cpus = 1) )
  expect_identical(
    mrpp(d, meta$uneven, weights = 'n(n-1)/2', permutations = 9, cpus = 1),
    mrpp(d, meta$uneven, weights = 3, permutations = 9, cpus = 1) )
  expect_error(mrpp(d, meta$trt, weights = 'n^2'), 'weights')
  expect_error(mrpp(d, meta$trt, weights = 4), 'weights')



  # Result format. ====

  res <- mrpp(d, 'trt', meta, permutations = 99, cpus = 1)

  expect_inherits(res, 'ecodive_mrpp')
  expect_inherits(res, 'ecodive_test')
  expect_identical(
    names(res),
    c('statistic', 'p_value', 'z', 'n_perm', 'expected_delta', 'A',
      'class_delta', 'class_n', 'weights', 'pairwise', 'p_adjust', 'permutation') )
  expect_identical(class(res), c('ecodive_mrpp', 'ecodive_test'))
  expect_identical(names(res$statistic), 'delta')
  expect_identical(attr(res, 'title'), 'Multi-response permutation procedure (MRPP)')

  # One row of the shared core fields.
  df <- as.data.frame(res)
  expect_identical(names(df), c('test', 'stat_name', 'statistic', 'p_value', 'z', 'n_perm'))
  expect_identical(df$test, 'mrpp')
  expect_identical(df$stat_name, 'delta')
  expect_identical(df$statistic, res$statistic[['delta']])
  expect_identical(df$n_perm, 99L)
  expect_identical(nrow(rbind(df, as.data.frame(mrpp(d, 'uneven', meta, permutations = 9, cpus = 1)))), 2L)
  expect_identical(names(res$class_delta), levels(meta$trt))
  expect_identical(res$class_n, c(ctrl = 10L, low = 10L, high = 10L))
  expect_identical(length(attr(res, 'delta_perm')), 99L)
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  expect_true(res$p_value >= 1 / 100 && res$p_value <= 1)
  expect_true(res$z < 0)  # A real treatment effect.
  out <- capture.output(print(res))
  expect_identical(out[[1]], 'Multi-response permutation procedure (MRPP)')
  expect_true(any(grepl('^delta = .*, p = .*, z = -', out)))
  expect_true(any(grepl('Permutations: 99 (free)', out, fixed = TRUE)))
  expect_true(any(grepl('agreement A', out)))
  expect_true(any(grepl('weights: n)', out, fixed = TRUE)))

  # A group of one sample has no within-group distances.
  res <- mrpp(d, 'uneven', meta, permutations = 9, cpus = 1)
  expect_true(is.na(res$class_delta[['d']]))
  expect_false(is.na(res$statistic[['delta']]))

  # No permutations: no p-value or z, but A.
  res <- mrpp(d, 'trt', meta, permutations = 0, cpus = 1)
  expect_true(is.na(res$p_value) && is.na(res$z))
  expect_false(is.na(res$A))
  expect_identical(attr(res, 'delta_perm'), numeric(0))
  expect_false(any(grepl('p = ', capture.output(print(res)), fixed = TRUE)))

  # One permutation: a p-value, but no z.
  res <- mrpp(d, 'trt', meta, permutations = 1, cpus = 1)
  expect_false(is.na(res$p_value))
  expect_true(is.na(res$z))



  # Deterministic: the same for any number of cpus, and from run to run. ====

  r1 <- mrpp(d, 'trt', meta, permutations = 199, seed = 7, cpus = 1)
  r2 <- mrpp(d, 'trt', meta, permutations = 199, seed = 7, cpus = 2)
  expect_identical(r1, r2)
  expect_identical(r1, mrpp(d, 'trt', meta, permutations = 199, seed = 7, cpus = 1))

  # The seed changes the permutations, not the statistics.
  r3 <- mrpp(d, 'trt', meta, permutations = 199, seed = 8, cpus = 1)
  expect_identical(r1[c('statistic', 'expected_delta', 'A', 'class_delta')],
                   r3[c('statistic', 'expected_delta', 'A', 'class_delta')])
  expect_false(identical(attr(r1, 'delta_perm'), attr(r3, 'delta_perm')))

  # More threads, at home only. Large enough for run_parallel() to start
  # threads.
  if (at_home()) {
    big <- bray(big_mtx)
    bg  <- rep(c('a', 'b', 'c', 'd'), 26)
    r1  <- mrpp(big, bg, permutations = 999, cpus = 1)
    r4  <- mrpp(big, bg, permutations = 999, cpus = 4)
    expect_identical(r1, r4)
    r1  <- mrpp(big, bg, permutations = 999, strata = rep(1:8, 13), cpus = 1)
    r4  <- mrpp(big, bg, permutations = 999, strata = rep(1:8, 13), cpus = 4)
    expect_identical(r1, r4)
  }



  # Strata. ====

  # Groups that are constant within strata can't change: every permuted
  # delta equals the observed one.
  res <- mrpp(d, 'block', meta, permutations = 49, strata = 'block', cpus = 1)
  expect_identical(res$p_value, 1)
  expect_identical(attr(res, 'delta_perm'), rep(res$statistic[['delta']], 49))
  expect_true(is.na(res$z))
  expect_identical(res$permutation, 'within strata')

  # Permutations stay within strata.
  res <- mrpp(d, meta$trt, permutations = 19, strata = meta$block, seed = 1, cpus = 1)
  expect_identical(res, mrpp(d, 'trt', meta, permutations = 19, strata = 'block', seed = 1, cpus = 1))

  expect_error(mrpp(d, 'trt', meta, permutations = ex$perms, strata = 'block'), 'not both')
  expect_error(mrpp(d, 'trt', meta, strata = 1:3), 'strata')
  expect_error(mrpp(d, meta$trt, strata = 1:3), 'strata')
  expect_error(mrpp(d, 'trt', meta, strata = c(NA, meta$block[-1])), 'NA')



  # Matching groups to samples. ====

  ref <- mrpp(d, 'trt', meta, permutations = 19, cpus = 1)

  # Rows of `data` by name ...
  shuffled <- meta[rev(seq_len(nrow(meta))), ]
  expect_identical(mrpp(d, 'trt', shuffled, permutations = 19, cpus = 1), ref)

  # ... or by position, without row names.
  unnamed <- meta
  rownames(unnamed) <- NULL
  expect_identical(mrpp(d, 'trt', unnamed, permutations = 19, cpus = 1), ref)

  # A vector, unnamed or named (in any order), factor or character.
  expect_identical(mrpp(d, meta$trt, permutations = 19, cpus = 1), ref)
  named <- setNames(meta$trt, rownames(meta))
  expect_identical(mrpp(d, rev(named), permutations = 19, cpus = 1), ref)
  chr <- mrpp(d, as.character(meta$trt), permutations = 19, cpus = 1)
  expect_identical(chr$statistic[['delta']], ref$statistic[['delta']])  # Levels sort differently.
  expect_identical(sort(chr$class_delta), sort(ref$class_delta))

  # A symmetric matrix works like a dist.
  expect_identical(mrpp(as.matrix(d), 'trt', meta, permutations = 19, cpus = 1), ref)

  # Rows for other samples are ignored, even with NAs, and so are their
  # factor levels.
  extra <- rbind(meta, meta[1, ])
  rownames(extra)[nrow(extra)] <- 'X'
  extra$trt <- factor(extra$trt, c(levels(meta$trt), 'only_x'))
  extra$trt[nrow(extra)] <- 'only_x'
  expect_identical(mrpp(d, 'trt', extra, permutations = 19, cpus = 1), ref)
  extra$trt[nrow(extra)] <- NA
  expect_identical(mrpp(d, 'trt', extra, permutations = 19, cpus = 1), ref)

  expect_error(mrpp(d, 'trt', meta[-1, ]), "1 sample in `dm` has no row in `data`: 'S01'")
  expect_error(mrpp(d, 'nothere', meta), "'nothere' is not one")
  expect_identical(mrpp(d, meta$trt, meta, permutations = 19, cpus = 1), mrpp(d, meta$trt, permutations = 19, cpus = 1))
  expect_error(mrpp(d, meta$trt[-1], meta), 'one value per row')
  expect_error(mrpp(d, 'trt'), '`data` is not given')
  expect_error(mrpp(d, meta$trt[-1]), 'has 29 values')
  expect_error(mrpp(d), '`groups` is required')
  na <- meta
  na$trt[3] <- NA
  expect_error(mrpp(d, 'trt', na), 'Missing values')
  expect_error(mrpp(d, na$trt), 'NAs')
  expect_error(mrpp(d, list(1, 2)), 'vector')



  # Other arguments and degenerate input. ====

  expect_error(mrpp(d, 'trt', meta, permutations = -1), 'permutations')
  expect_error(mrpp(d, 'trt', meta, permutations = ex$perms[, -1]), '30 columns')
  expect_error(mrpp(d, 'trt', meta, seed = 1.5), 'seed')
  expect_error(mrpp(d, 'trt', meta, cpus = 0), 'cpus')

  bad_d <- d
  bad_d[1] <- NA
  expect_error(mrpp(bad_d, 'trt', meta), 'NA')
  bad_d[1] <- -1
  expect_error(mrpp(bad_d, 'trt', meta), 'negative')
  expect_error(mrpp(as.dist(matrix(0, 1, 1)), 1), 'at least 2')

  # Every sample in its own group.
  expect_warning(res <- mrpp(d, rownames(meta), permutations = 9, cpus = 1), 'No group')
  expect_true(is.na(res$statistic[['delta']]) && is.na(res$A) && is.na(res$p_value) && is.na(res$z))
  expect_true(all(is.na(res$class_delta)))

  # One group.
  expect_warning(res <- mrpp(d, rep('a', n), permutations = 9, cpus = 1), 'one group')
  expect_equal(res$statistic[['delta']], res$expected_delta, tolerance = 1e-14)
  expect_identical(res$p_value, 1)

  # All distances 0.
  expect_warning(res <- mrpp(d * 0, 'trt', meta, permutations = 9, cpus = 1), 'All distances are 0')
  expect_true(is.na(res$A))
  expect_identical(res$statistic[['delta']], 0)

  # Two samples.
  d2  <- as.dist(matrix(c(0, 1, 1, 0), 2))
  expect_warning(res <- mrpp(d2, c('a', 'a'), permutations = 3, cpus = 1), 'one group')
  expect_identical(res$statistic[['delta']], 1)
  expect_warning(res <- mrpp(d2, c('a', 'b'), permutations = 3, cpus = 1), 'No group')
  expect_true(is.na(res$statistic[['delta']]))

#})
