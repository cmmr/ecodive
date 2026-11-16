#test_that("permdisp", {

  ex   <- readRDS('expected/permdisp.rds')
  meta <- ex$meta
  dms  <- list(bray = bray(ex$counts), euclidean = stats::dist(ex$counts))



  # The same results as betadisper() + permutest(), given the same
  # permutations: stored results, so at home only. ====

  if (at_home()) {
    for (cs in ex$cases) {

      info <- paste(cs$dist, cs$type, cs$bias, cs$groups)
      res  <- permdisp(
        dms[[cs$dist]], meta[[cs$groups]], type = cs$type, bias_adjust = cs$bias,
        pairwise = TRUE, permutations = ex$perms, cpus = 1 )

      tbl <- unname(as.matrix(res))

      # Medians: betadisper()'s optimizer stops early, so compare with it
      # loosely, and with the exact spatial medians to rounding.
      tol <- if (cs$type == 'median') 1e-3 else 1e-10

      expect_equal(unname(attr(res, 'distances')), cs$distances, tolerance = tol, info = info)
      expect_equal(attr(res, 'dispersions'), c(cs$dispersions), tolerance = tol, info = info)
      expect_equal(tbl[, 1:5], cs$table[, 1:5], tolerance = tol, info = info)
      expect_equal(unname(attr(res, 'f_perm')[, 1]), cs$f_perm, tolerance = tol, info = info)

      pw <- attr(res, 'pairwise')
      expect_equal(pw$t, cs$t_obs, tolerance = tol, info = info)
      expect_equal(pw[['Pr(>|t|)']], unname(cs$pairwise$observed), tolerance = tol, info = info)
      if (!is.null(cs$t_perm))
        expect_equal(unname(attr(res, 't_perm')), cs$t_perm, tolerance = tol, info = info)

      if (is.null(cs$exact)) {
        # Exact: identical p-values.
        expect_identical(tbl[, 6], cs$table[, 6], info = info)
        expect_identical(pw[['Pr(perm)']], unname(cs$pairwise$permuted), info = info)
      }
      else {
        expect_equal(unname(attr(res, 'distances')), cs$exact, tolerance = 1e-10, info = info)
      }
    }

    # Free permutation within strata: the reference's p-value, up to
    # sampling error.
    res <- permdisp(dms$bray, 'half', meta, type = 'centroid', permutations = 9999, strata = 'block', cpus = 2)
    expect_true(abs(res[['Pr(>F)']][[1]] - ex$strata_p) < 0.01)
  }



  # Identities that hold without the references. ====

  d  <- dms$bray
  dv <- as.matrix(d)
  g  <- meta$trt

  # Centroid distances from principal coordinates, real minus imaginary.
  G  <- -dv^2 / 2
  G  <- G - rowMeans(G)
  G  <- t(t(G) - colMeans(G))
  e  <- eigen(G, symmetric = TRUE)
  ok <- abs(e$values) > max(sqrt(.Machine$double.eps), sqrt(.Machine$double.eps) * e$values[[1]])
  X  <- e$vectors[, ok] %*% diag(sqrt(abs(e$values[ok])))
  sg <- sign(e$values[ok])
  z2 <- numeric(nrow(X))
  for (lv in levels(g)) {
    i <- g == lv
    C <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]))
    z2[i] <- C^2 %*% sg
  }
  res <- permdisp(d, g, type = 'centroid', permutations = 0, cpus = 1)
  expect_equal(unname(attr(res, 'distances')), sqrt(pmax(z2, 0)), tolerance = 1e-10)

  # The ANOVA is lm()'s.
  z <- attr(res, 'distances')
  a <- stats::anova(stats::lm(z ~ g))
  expect_equal(res$Df, a$Df)
  expect_equal(res[['Sum Sq']], a[['Sum Sq']], tolerance = 1e-12)
  expect_equal(res$F[[1]], a[['F value']][[1]], tolerance = 1e-12)
  expect_equal(unname(attr(res, 'dispersions')), as.vector(tapply(z, g, mean)), tolerance = 1e-14)

  # Euclidean distances: the spatial median minimizes the summed distance,
  # and the median of the distances themselves works.
  de  <- dms$euclidean
  res <- permdisp(de, g, permutations = 0, cpus = 1)
  zm  <- attr(res, 'distances')
  zc  <- attr(permdisp(de, g, type = 'centroid', permutations = 0, cpus = 1), 'distances')
  for (lv in levels(g))
    expect_true(sum(zm[g == lv]) <= sum(zc[g == lv]) * (1 + 1e-12))

  # The identity permutation: the residuals' group means are all 0, so F is
  # 0, while the pairwise t statistics (of the distances) are the observed
  # ones.
  id  <- matrix(seq_len(36), 5, 36, byrow = TRUE)
  res <- permdisp(d, g, permutations = id, pairwise = TRUE, cpus = 1)
  expect_equal(unname(attr(res, 'f_perm')[, 1]), rep(0, 5), tolerance = 1e-12)
  expect_equal(unname(attr(res, 't_perm')), matrix(attr(res, 'pairwise')$t, 5, 3, byrow = TRUE), tolerance = 1e-12)
  expect_identical(res[['Pr(>F)']][[1]], 1 / 6)
  expect_identical(attr(res, 'pairwise')[['Pr(perm)']], c(1, 1, 1))

  # bias_adjust scales each distance by sqrt(n / (n - 1)).
  r0 <- permdisp(d, g, type = 'centroid', permutations = 0, cpus = 1)
  r1 <- permdisp(d, g, type = 'centroid', bias_adjust = TRUE, permutations = 0, cpus = 1)
  ng <- as.vector(table(g))[as.integer(g)]
  expect_equal(attr(r1, 'distances'), attr(r0, 'distances') * sqrt(ng / (ng - 1)), tolerance = 1e-14)
  expect_true(grepl('small groups', attr(r1, 'heading')))



  # Result format. ====

  res <- permdisp(d, g, permutations = 99, cpus = 1)

  expect_inherits(res, 'anova')
  expect_inherits(res, 'data.frame')
  expect_identical(names(res), c('Df', 'Sum Sq', 'Mean Sq', 'F', 'N.Perm', 'Pr(>F)'))
  expect_identical(rownames(res), c('Groups', 'Residuals'))
  expect_identical(res$Df, c(2, 33))
  expect_identical(res$N.Perm, c(99, NA))
  expect_identical(names(attr(res, 'distances')), rownames(meta))
  expect_identical(names(attr(res, 'dispersions')), levels(g))
  expect_identical(dim(attr(res, 'f_perm')), c(99L, 1L))
  expect_identical(colnames(attr(res, 'f_perm')), 'Groups')
  expect_identical(attr(res, 'type'), 'median')
  expect_null(attr(res, 'pairwise'))
  expect_null(attr(res, 't_perm'))
  expect_true(res[['Pr(>F)']][[1]] >= 1 / 100)
  expect_true(any(grepl('spatial medians', attr(res, 'heading'))))
  expect_true(any(grepl('Number of permutations: 99', attr(res, 'heading'))))

  res <- permdisp(d, g, type = 'centroid', pairwise = TRUE, permutations = 19, cpus = 1)
  pw  <- attr(res, 'pairwise')
  expect_identical(rownames(pw), c('ctrl-low', 'ctrl-high', 'low-high'))
  expect_identical(names(pw), c('t', 'Df', 'Pr(>|t|)', 'Pr(perm)'))
  expect_identical(pw$Df, c(24, 22, 20))
  expect_identical(dim(attr(res, 't_perm')), c(19L, 3L))
  expect_identical(colnames(attr(res, 't_perm')), rownames(pw))
  expect_true(any(grepl('centroids', attr(res, 'heading'))))

  # No permutations: no p-values.
  res <- permdisp(d, g, permutations = 0, pairwise = TRUE, cpus = 1)
  expect_true(is.na(res[['Pr(>F)']][[1]]))
  expect_identical(dim(attr(res, 'f_perm')), c(0L, 1L))
  expect_true(all(is.na(attr(res, 'pairwise')[['Pr(perm)']])))

  # A group of one sample: no t statistics for its pairs.
  res <- permdisp(d, meta$solo, type = 'centroid', pairwise = TRUE, permutations = 9, cpus = 1)
  pw  <- attr(res, 'pairwise')
  solo <- grepl('d', rownames(pw))
  expect_true(all(is.na(pw$t[solo])) && all(!is.na(pw$t[!solo])))
  expect_true(all(is.na(attr(res, 't_perm')[, solo])))
  expect_identical(unname(attr(res, 'distances')[36]), 0)
  r1 <- permdisp(d, meta$solo, type = 'centroid', bias_adjust = TRUE, permutations = 0, cpus = 1)
  expect_identical(unname(attr(r1, 'distances')[36]), 0)



  # Deterministic: the same for any number of cpus, and from run to run. ====

  for (type in c('median', 'centroid')) {
    r1 <- permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1)
    r2 <- permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 2)
    expect_identical(r1, r2, info = type)
    expect_identical(r1, permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1))
  }

  # Large enough for run_parallel() to start threads.
  big <- bray(big_mtx + seq_len(104) %% 7)
  bg  <- rep(c('a', 'b', 'c', 'd'), 26)
  for (type in c('median', 'centroid')) {
    r1 <- permdisp(big, bg, type = type, pairwise = TRUE, permutations = 999, cpus = 1)
    r2 <- permdisp(big, bg, type = type, pairwise = TRUE, permutations = 999, cpus = 2)
    expect_identical(r1, r2, info = type)
    if (at_home()) {
      for (cpus in c(3L, 4L, 8L)) {
        expect_identical(r1, permdisp(
          big, bg, type = type, pairwise = TRUE, permutations = 999, cpus = cpus ),
          info = paste(type, cpus) )
      }
    }
  }

  # The seed changes the permutations, not the statistics.
  r2 <- permdisp(d, g, permutations = 99, seed = 2, cpus = 1)
  r3 <- permdisp(d, g, permutations = 99, seed = 3, cpus = 1)
  expect_identical(r2[, 1:5], r3[, 1:5])
  expect_false(identical(attr(r2, 'f_perm'), attr(r3, 'f_perm')))



  # Strata. ====

  # One stratum is free permutation.
  expect_identical(
    permdisp(d, g, strata = rep(1, 36), permutations = 49, cpus = 1)[, 1:6],
    permdisp(d, g, permutations = 49, cpus = 1)[, 1:6] )

  # Permuting only within the groups leaves the residuals' group means at 0:
  # every permuted F is 0.
  res <- permdisp(d, g, type = 'centroid', strata = g, permutations = 19, cpus = 1)
  expect_equal(as.vector(attr(res, 'f_perm')), rep(0, 19), tolerance = 1e-10)
  expect_true(grepl('within strata', attr(res, 'heading')))

  expect_error(permdisp(d, g, permutations = ex$perms, strata = meta$block), 'not both')
  expect_error(permdisp(d, g, strata = 1:3), 'strata')



  # Groups and metadata. ====

  ref <- permdisp(d, g, permutations = 19, cpus = 1)

  # A column of `data`, matched by row name.
  expect_identical(permdisp(d, 'trt', meta, permutations = 19, cpus = 1), ref)
  expect_identical(permdisp(d, 'trt', meta[36:1, ], permutations = 19, cpus = 1), ref)

  # Named groups are matched by name; strings and factors are the same.
  sh <- setNames(as.character(g), rownames(meta))[36:1]
  expect_identical(permdisp(d, sh, permutations = 19, cpus = 1)[, 1:6], ref[, 1:6])
  expect_identical(attr(permdisp(d, sh, permutations = 19, cpus = 1), 'distances'), attr(ref, 'distances'))

  # A symmetric matrix works like a dist.
  expect_identical(permdisp(as.matrix(d), g, permutations = 19, cpus = 1), ref)

  # Rows for other samples are ignored, with their levels.
  extra <- rbind(meta, meta[1, ])
  rownames(extra)[37] <- 'X'
  extra$trt <- factor(extra$trt, c(levels(g), 'only_x'))
  extra$trt[37] <- 'only_x'
  expect_identical(permdisp(d, 'trt', extra, permutations = 19, cpus = 1), ref)

  # Strata as a column.
  expect_identical(
    permdisp(d, 'trt', meta, strata = 'block', permutations = 19, cpus = 1),
    permdisp(d, g, strata = meta$block, permutations = 19, cpus = 1) )



  # Edge cases and errors. ====

  expect_warning(r <- permdisp(d, rep('a', 36), permutations = 9), 'only one group')
  expect_true(is.na(r$F[[1]]) && is.na(r[['Pr(>F)']][[1]]))

  d2 <- stats::dist(matrix(c(0, 1, 3), 3))
  expect_warning(r <- permdisp(d2, c('a', 'b', 'c'), permutations = 9), 'only one sample')
  expect_true(is.na(r$F[[1]]))
  expect_identical(unname(attr(r, 'distances')), c(0, 0, 0))

  # Two samples per group, in one dimension.
  r <- permdisp(stats::dist(c(0, 1, 5, 7)), c('a', 'a', 'b', 'b'), permutations = 0)
  expect_equal(unname(attr(r, 'distances')), c(0.5, 0.5, 1, 1), tolerance = 1e-14)

  # Non-metric distances can put a sample closer to its centre on the
  # imaginary axes than on the real ones: its distance is then 0.
  dn <- stats::as.dist(matrix(c(
    0, 1, 1, 2,
    1, 0, 3, 2,
    1, 3, 0, 2,
    2, 2, 2, 0 ), 4))
  for (type in c('centroid', 'median'))
    expect_warning(r <- permdisp(dn, c('a', 'a', 'a', 'b'), type = type, permutations = 0), 'negative')
  expect_identical(unname(attr(r, 'distances')[c(1, 4)]), c(0, 0))

  expect_error(permdisp(d), '`groups` is required')
  expect_error(permdisp(d, g[-1]), '35 rows')
  expect_error(permdisp(d, c(NA, as.character(g[-1]))), 'NA')
  expect_error(permdisp(d, 'nothere', meta), 'name of a column')
  expect_error(permdisp(d, g, meta), 'name of a column')
  expect_error(permdisp(d, matrix(1, 36, 2)), 'vector or factor')
  expect_error(permdisp(d, g, type = 'mean'), "'arg' should be one of")
  expect_error(permdisp(d, g, bias_adjust = NA), 'bias_adjust')
  expect_error(permdisp(d, g, pairwise = 'yes'), 'pairwise')
  expect_error(permdisp(d, g, permutations = -1), 'permutations')
  expect_error(permdisp(d, g, seed = 1.5), 'seed')
  expect_error(permdisp(d, g, cpus = 0), 'cpus')
  expect_error(permdisp(stats::dist(1), 1), 'at least 2')
  expect_error(permdisp(d, g, permutations = ex$perms[, -1]), 'columns')

#})
