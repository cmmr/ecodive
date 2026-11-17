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

      tbl <- unname(as.matrix(res$table))

      # Medians: betadisper()'s optimizer stops early, so compare with it
      # loosely, and with the exact spatial medians to rounding.
      tol <- if (cs$type == 'median') 1e-3 else 1e-10

      expect_equal(unname(res$distances), cs$distances, tolerance = tol, info = info)
      expect_equal(res$dispersions, c(cs$dispersions), tolerance = tol, info = info)
      expect_equal(tbl[, 1:5], cs$table[, 1:5], tolerance = tol, info = info)
      expect_equal(unname(attr(res, 'f_perm')), cs$f_perm, tolerance = tol, info = info)

      pw <- res$pairwise
      expect_equal(pw$t, cs$t_obs, tolerance = tol, info = info)
      expect_equal(pw$p_value, unname(cs$pairwise$observed), tolerance = tol, info = info)
      if (!is.null(cs$t_perm))
        expect_equal(unname(attr(res, 't_perm')), cs$t_perm, tolerance = tol, info = info)

      if (is.null(cs$exact)) {
        # Exact: identical p-values.
        expect_identical(tbl[, 6], cs$table[, 6], info = info)
        expect_identical(pw$p_value_perm, unname(cs$pairwise$permuted), info = info)
      }
      else {
        expect_equal(unname(res$distances), cs$exact, tolerance = 1e-10, info = info)
      }
    }

    # Free permutation within strata: the reference's p-value, up to
    # sampling error.
    res <- permdisp(dms$bray, 'half', meta, type = 'centroid', permutations = 9999, strata = 'block', cpus = 2)
    expect_true(abs(res$p_value - ex$strata_p) < 0.01)
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
  expect_equal(unname(res$distances), sqrt(pmax(z2, 0)), tolerance = 1e-10)

  # The ANOVA is lm()'s.
  z <- res$distances
  a <- stats::anova(stats::lm(z ~ g))
  expect_equal(res$table$Df, a$Df)
  expect_equal(res$table[['Sum Sq']], a[['Sum Sq']], tolerance = 1e-12)
  expect_equal(res$statistic[['F']], a[['F value']][[1]], tolerance = 1e-12)
  expect_equal(unname(res$dispersions), as.vector(tapply(z, g, mean)), tolerance = 1e-14)

  # Euclidean distances: the spatial median minimizes the summed distance,
  # where the centroid minimizes the summed squared distance.
  de  <- dms$euclidean
  zm  <- permdisp(de, g, type = 'median', permutations = 0, cpus = 1)$distances
  zc  <- permdisp(de, g, type = 'centroid', permutations = 0, cpus = 1)$distances
  for (lv in levels(g))
    expect_true(sum(zm[g == lv]) <= sum(zc[g == lv]) * (1 + 1e-12))

  # The identity permutation: the residuals' group means are all 0, so F is
  # 0, while the pairwise t statistics (of the distances) are the observed
  # ones.
  id  <- matrix(seq_len(36), 5, 36, byrow = TRUE)
  res <- permdisp(d, g, type = 'median', permutations = id, pairwise = TRUE, cpus = 1)
  expect_equal(unname(attr(res, 'f_perm')), rep(0, 5), tolerance = 1e-12)
  expect_equal(unname(attr(res, 't_perm')), matrix(res$pairwise$t, 5, 3, byrow = TRUE), tolerance = 1e-12)
  expect_identical(res$p_value, 1 / 6)
  expect_identical(res$pairwise$p_value_perm, c(1, 1, 1))

  # bias_adjust scales each distance by sqrt(n / (n - 1)).
  r0 <- permdisp(d, g, type = 'centroid', permutations = 0, cpus = 1)
  r1 <- permdisp(d, g, type = 'centroid', bias_adjust = TRUE, permutations = 0, cpus = 1)
  ng <- as.vector(table(g))[as.integer(g)]
  expect_equal(r1$distances, r0$distances * sqrt(ng / (ng - 1)), tolerance = 1e-14)
  expect_true(r1$bias_adjust)
  expect_true(any(grepl('small groups', test_details.ecodive_permdisp(r1, 4))))



  # Result format. ====

  # The default is the centroid.
  res <- permdisp(d, g, permutations = 99, cpus = 1)
  expect_identical(res, permdisp(d, g, type = 'centroid', permutations = 99, cpus = 1))

  expect_identical(class(res), c('ecodive_permdisp', 'ecodive_test'))
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm', 'table', 'distances',
    'dispersions', 'type', 'bias_adjust', 'pairwise', 'permutation' ))
  expect_identical(names(res$statistic), 'F')
  expect_identical(res$statistic[['F']], res$table$F[[1]])
  expect_identical(res$p_value, res$table[['Pr(>F)']][[1]])
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  f <- attr(res, 'f_perm')
  expect_true(is.double(f) && is.null(dim(f)) && length(f) == 99)
  expect_equal(res$z, (res$statistic[['F']] - mean(f)) / stats::sd(f), tolerance = 1e-14)
  expect_identical(attr(res, 'title'), 'Homogeneity of multivariate dispersions (PERMDISP)')

  tbl <- res$table
  expect_inherits(tbl, 'anova')
  expect_identical(names(tbl), c('Df', 'Sum Sq', 'Mean Sq', 'F', 'N.Perm', 'Pr(>F)'))
  expect_identical(rownames(tbl), c('Groups', 'Residuals'))
  expect_identical(tbl$Df, c(2, 33))
  expect_identical(tbl$N.Perm, c(99, NA))
  expect_identical(names(res$distances), rownames(meta))
  expect_identical(names(res$dispersions), levels(g))
  expect_identical(res$type, 'centroid')
  expect_false(res$bias_adjust)
  expect_null(res$pairwise)
  expect_null(attr(res, 't_perm'))
  expect_true(res$p_value >= 1 / 100)
  expect_true(any(grepl('centroids', test_details.ecodive_permdisp(res, 4))))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl('Homogeneity of multivariate dispersions', out, fixed = TRUE)))
  expect_true(any(grepl('Permutations: 99 (free)', out, fixed = TRUE)))
  expect_true(any(grepl('centroids', out)))

  res <- permdisp(d, g, type = 'median', pairwise = TRUE, permutations = 19, cpus = 1)
  expect_identical(res$type, 'median')
  pw  <- res$pairwise
  expect_identical(rownames(pw), c('ctrl-low', 'ctrl-high', 'low-high'))
  expect_identical(names(pw), c('t', 'Df', 'p_value', 'p_value_perm'))
  expect_identical(pw$Df, c(24, 22, 20))
  expect_identical(dim(attr(res, 't_perm')), c(19L, 3L))
  expect_identical(colnames(attr(res, 't_perm')), rownames(pw))
  out <- test_details.ecodive_permdisp(res, 4)
  expect_true(any(grepl('spatial medians', out)) && any(grepl('ctrl-high', out)))

  out <- utils::capture.output(print(res))
  expect_true(any(grepl('Pairwise t tests', out, fixed = TRUE)))

  # One row for comparing tests.
  df <- as.data.frame(res)
  expect_identical(df$test, 'permdisp')
  expect_identical(df$stat_name, 'F')
  expect_identical(df$statistic, res$statistic[['F']])
  expect_identical(df$p_value, res$p_value)

  # No permutations: no p-values.
  res <- permdisp(d, g, permutations = 0, pairwise = TRUE, cpus = 1)
  expect_true(is.na(res$p_value))
  expect_identical(attr(res, 'f_perm'), numeric(0))
  expect_true(is.na(res$z))
  expect_true(all(is.na(res$pairwise$p_value_perm)))

  # A group of one sample: no t statistics for its pairs.
  res <- permdisp(d, meta$solo, type = 'centroid', pairwise = TRUE, permutations = 9, cpus = 1)
  pw  <- res$pairwise
  solo <- grepl('d', rownames(pw))
  expect_true(all(is.na(pw$t[solo])) && all(!is.na(pw$t[!solo])))
  expect_true(all(is.na(attr(res, 't_perm')[, solo])))
  expect_identical(unname(res$distances[36]), 0)
  r1 <- permdisp(d, meta$solo, type = 'centroid', bias_adjust = TRUE, permutations = 0, cpus = 1)
  expect_identical(unname(r1$distances[36]), 0)



  # Deterministic: the same for any number of cpus, and from run to run. ====

  for (type in c('median', 'centroid')) {
    r1 <- permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1)
    r2 <- permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 2)
    expect_identical(r1, r2, info = type)
    expect_identical(r1, permdisp(d, g, type = type, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1))
  }

  # Large enough for threads: the medians' Weiszfeld steps and the
  # permutations (n = 400), and, at home, the centroids' sums (n = 3000).
  sim <- function (n) {
    m <- matrix(((seq_len(n * 12) * 7919) %% 23) * (seq_len(n * 12) %% 3 > 0), n, 12)
    m[, 1] <- m[, 1] + 1
    bray(m)
  }
  big <- sim(400)
  bg  <- rep(c('a', 'b', 'c', 'd'), 100)
  cpus_list <- if (at_home()) c(2L, 3L, 4L, 8L) else 2L
  for (type in c('median', 'centroid')) {
    r1 <- permdisp(big, bg, type = type, pairwise = TRUE, permutations = 2999, cpus = 1)
    for (cpus in cpus_list)
      expect_identical(r1, permdisp(
        big, bg, type = type, pairwise = TRUE, permutations = 2999, cpus = cpus ),
        info = paste(type, cpus) )
  }
  if (at_home()) {
    big <- sim(3000)
    bg  <- rep(c('a', 'b', 'c'), 1000)
    r1  <- permdisp(big, bg, type = 'centroid', permutations = 99, cpus = 1)
    for (cpus in c(2L, 4L, 8L))
      expect_identical(r1, permdisp(big, bg, type = 'centroid', permutations = 99, cpus = cpus))
  }

  # The seed changes the permutations, not the statistics.
  for (type in c('centroid', 'median')) {
    r2 <- permdisp(d, g, type = type, permutations = 99, seed = 2, cpus = 1)
    r3 <- permdisp(d, g, type = type, permutations = 99, seed = 3, cpus = 1)
    expect_identical(r2$table[, 1:5], r3$table[, 1:5], info = type)
    expect_false(identical(attr(r2, 'f_perm'), attr(r3, 'f_perm')), info = type)
  }



  # Strata. ====

  # One stratum is free permutation.
  for (type in c('centroid', 'median'))
    expect_identical(
      permdisp(d, g, type = type, strata = rep(1, 36), permutations = 49, cpus = 1)[1:4],
      permdisp(d, g, type = type, permutations = 49, cpus = 1)[1:4], info = type )

  # Permuting only within the groups leaves the residuals' group means at 0:
  # every permuted F is 0.
  res <- permdisp(d, g, type = 'centroid', strata = g, permutations = 19, cpus = 1)
  expect_equal(as.vector(attr(res, 'f_perm')), rep(0, 19), tolerance = 1e-10)
  expect_identical(res$permutation, 'within strata')
  expect_identical(permdisp(d, g, permutations = ex$perms)$permutation, 'supplied matrix')

  expect_error(permdisp(d, g, permutations = ex$perms, strata = meta$block), 'not both')
  expect_error(permdisp(d, g, strata = 1:3), 'strata')



  # Groups and metadata. ====

  ref <- permdisp(d, g, permutations = 19, cpus = 1)

  # A column of `data`, matched by row name.
  expect_identical(permdisp(d, 'trt', meta, permutations = 19, cpus = 1), ref)
  expect_identical(permdisp(d, 'trt', meta[36:1, ], permutations = 19, cpus = 1), ref)

  # Named groups are matched by name; strings and factors are the same.
  sh <- setNames(as.character(g), rownames(meta))[36:1]
  expect_identical(permdisp(d, sh, permutations = 19, cpus = 1)$table, ref$table)
  expect_identical(permdisp(d, sh, permutations = 19, cpus = 1)$distances, ref$distances)

  # A symmetric matrix works like a dist.
  expect_identical(permdisp(as.matrix(d), g, permutations = 19, cpus = 1), ref)
  expect_identical(
    permdisp(as.matrix(d), g, type = 'median', permutations = 19, cpus = 1),
    permdisp(d, 'trt', meta[36:1, ], type = 'median', permutations = 19, cpus = 1) )

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
  expect_true(is.na(r$statistic[['F']]) && is.na(r$p_value))

  d2 <- stats::dist(matrix(c(0, 1, 3), 3))
  expect_warning(r <- permdisp(d2, c('a', 'b', 'c'), permutations = 9), 'only one sample')
  expect_true(is.na(r$statistic[['F']]))
  expect_identical(unname(r$distances), c(0, 0, 0))

  # Two samples per group, in one dimension.
  for (type in c('centroid', 'median')) {
    # Each group's samples are equally far from its centre: F is infinite.
    expect_warning(
      r <- permdisp(stats::dist(c(0, 1, 5, 7)), c('a', 'a', 'b', 'b'), type = type, permutations = 0),
      'no residual variation', info = type )
    expect_equal(unname(r$distances), c(0.5, 0.5, 1, 1), tolerance = 1e-14, info = type)
    expect_identical(r$statistic[['F']], Inf, info = type)
  }

  # Non-metric distances can put a sample closer to its centre on the
  # imaginary axes than on the real ones: its distance is then 0.
  dn <- stats::as.dist(matrix(c(
    0, 1, 1, 2,
    1, 0, 3, 2,
    1, 3, 0, 2,
    2, 2, 2, 0 ), 4))
  for (type in c('centroid', 'median'))
    expect_warning(r <- permdisp(dn, c('a', 'a', 'a', 'b'), type = type, permutations = 0), 'negative')
  expect_identical(unname(r$distances[c(1, 4)]), c(0, 0))

  expect_error(permdisp(d), '`groups` is required')
  expect_error(permdisp(d, g[-1]), 'has 35 values')
  expect_error(permdisp(d, c(NA, as.character(g[-1]))), 'NA')
  expect_error(permdisp(d, 'nothere', meta), 'name of a column')
  expect_identical(permdisp(d, g, meta, permutations = 9, cpus = 1), permdisp(d, g, permutations = 9, cpus = 1))
  expect_error(permdisp(d, g[-1], meta), 'one value per row')
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
