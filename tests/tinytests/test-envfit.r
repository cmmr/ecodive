#test_that("envfit", {

  ex <- readRDS('expected/envfit.rds')



  # Same arrows, r2, centroids as vegan::envfit(), and identical p-values
  # given the same permutations: stored results, so at home only. ====

  if (at_home()) {
    for (nm in names(ex$cases)) {

      cs    <- ex$cases[[nm]]
      data  <- ex[[cs$data]]
      perms <- if (isTRUE(cs$args$na.rm)) ex$perms_na else ex$perms

      input <- if (cs$data == 'counts') list(counts = data) else list(data = data)
      res   <- do.call(envfit, c(
        list(ex[[cs$ord]]), input,
        list(axes = cs$axes, permutations = perms, cpus = 1), cs$args ))

      vec <- is.na(res$level)

      if (!is.null(cs$r)) {
        # vegan gives p = 1 for constant variables; envfit() gives NA.
        const <- is.na(res$p_value[vec])
        expect_identical(res$variable[vec], names(cs$r), info = nm)
        expect_equal(unname(as.matrix(res[vec, -c(1:2, ncol(res) - 2:0)][!const, ])),
                     unname(cs$arrows[!const, , drop = FALSE]), tolerance = 1e-10, info = nm)
        expect_equal(res$r2[vec], unname(cs$r), tolerance = 1e-10, info = nm)
        expect_identical(res$p_value[vec][!const], unname(cs$pvals[!const]), info = nm)
        expect_true(all(cs$pvals[const] == 1), info = nm)
      }

      if (!is.null(cs$f_r)) {
        lev <- res[!vec, ]
        expect_identical(rownames(lev), rownames(cs$centroids), info = nm)
        expect_equal(unname(as.matrix(lev[, -c(1:2, ncol(lev) - 2:0)])),
                     unname(cs$centroids), tolerance = 1e-10, info = nm)
        r2 <- tapply(lev$r2, factor(lev$variable, unique(lev$variable)), `[`, 1)
        p  <- tapply(lev$p_value,  factor(lev$variable, unique(lev$variable)), `[`, 1)
        expect_equal(as.vector(r2), unname(cs$f_r), tolerance = 1e-10, info = nm)
        expect_identical(as.vector(p), unname(cs$f_pvals), info = nm)
      }
    }
  }



  # Result format. ====

  ord <- ex$pco
  env <- ex$env

  res <- envfit(ord, env, permutations = 99, cpus = 1)

  expect_inherits(res, 'ecodive_envfit')
  expect_inherits(res, 'data.frame')
  expect_identical(names(res), c('variable', 'level', 'PCoA1', 'PCoA2', 'PCoA3', 'r2', 'p_value', 'z'))
  expect_identical(attr(res, 'n_perm'), 99L)
  expect_identical(
    rownames(res),
    c('N', 'P', 'K', 'Ca', 'pH', 'Humdepth', 'grazedno', 'grazedyes',
      'soilhigh', 'soillow', 'soilmid', 'zoneA', 'zoneB', 'zoneC', 'zoneD') )
  expect_identical(res$variable[1:8], c('N', 'P', 'K', 'Ca', 'pH', 'Humdepth', 'grazed', 'grazed'))
  expect_true(all(is.na(res$level[1:6])))
  expect_identical(res$level[7:8], c('no', 'yes'))

  # Unit arrows; r2 in [0, 1]; p in [1/100, 1].
  expect_equal(rowSums(as.matrix(res[1:6, 3:5])^2), rep(1, 6), tolerance = 1e-12, check.attributes = FALSE)
  expect_true(all(res$r2 >= 0 & res$r2 <= 1))
  expect_true(all(res$p_value >= 1 / 100 & res$p_value <= 1))

  rp <- attr(res, 'r2_perm')
  expect_identical(dim(rp), c(99L, 9L))
  expect_identical(colnames(rp), c('N', 'P', 'K', 'Ca', 'pH', 'Humdepth', 'grazed', 'soil', 'zone'))
  expect_equal(
    res$p_value[c(1:6, 7, 9, 12)],
    (colSums(rp >= rep(res$r2[c(1:6, 7, 9, 12)] - sqrt(.Machine$double.eps), each = 99)) + 1) / 100,
    check.attributes = FALSE )

  # z: r2 against its permutations.
  zc <- c(1:6, 7, 9, 12)
  expect_equal(
    res$z[zc], (res$r2[zc] - colMeans(rp)) / apply(rp, 2, stats::sd),
    tolerance = 1e-10, check.attributes = FALSE )

  # Centroids are the levels' mean coordinates.
  expect_equal(
    unname(unlist(res['grazedyes', 3:5])),
    unname(colMeans(ord[env$grazed == 'yes', ])), tolerance = 1e-12 )

  # No permutations: no p-values.
  res0 <- envfit(ord, env, permutations = 0, cpus = 1)
  expect_true(all(is.na(res0$p_value)))
  expect_true(all(is.na(res0$z)))
  expect_identical(attr(res0, 'n_perm'), 0L)
  expect_identical(dim(attr(res0, 'r2_perm')), c(0L, 9L))
  expect_identical(res0$r2, res$r2)

  # axes: indices or names; NULL = all.
  a12 <- envfit(ord, env, axes = 1:2, permutations = 9, cpus = 1)
  expect_identical(names(a12)[3:4], c('PCoA1', 'PCoA2'))
  expect_identical(a12, envfit(ord, env, axes = c('PCoA1', 'PCoA2'), permutations = 9, cpus = 1))
  a3 <- envfit(ord, env, axes = 3, permutations = 9, cpus = 1)
  expect_identical(names(a3)[3], 'PCoA3')
  expect_equal(abs(a3$PCoA3[1:6]), rep(1, 6))

  # A plain matrix without names; a data frame `ord`.
  m <- unname(ord)
  r <- envfit(m, unname(as.matrix(env[1:3])), permutations = 9, cpus = 1)
  expect_identical(names(r), c('variable', 'level', 'Dim1', 'Dim2', 'Dim3', 'r2', 'p_value', 'z'))
  expect_identical(r$variable, c('V1', 'V2', 'V3'))
  expect_equal(
    unname(as.matrix(envfit(as.data.frame(ord), env, permutations = 9, cpus = 1)[3:7])),
    unname(as.matrix(envfit(ord, env, permutations = 9, cpus = 1)[3:7])) )

  # A single vector or factor.
  g <- envfit(ord, setNames(env$grazed, rownames(env)), permutations = 9, cpus = 1)
  expect_identical(rownames(g), c('groupno', 'groupyes'))
  expect_equal(unique(g$r2), res$r2[7])



  # Matching samples by name. ====

  shuf <- env[c(24:1), ]
  expect_identical(
    envfit(ord, shuf, permutations = 49, cpus = 1)[, -1],
    envfit(ord, env,  permutations = 49, cpus = 1)[, -1] )

  # Extra rows in `data` are ignored.
  extra <- rbind(env, env[1, ])
  rownames(extra)[25] <- 'other'
  expect_identical(
    envfit(ord, extra, permutations = 49, cpus = 1),
    envfit(ord, env,   permutations = 49, cpus = 1) )

  # Row names that are only numbers (as here) count as names only when
  # they all match; otherwise rows are matched by position.
  expect_error(envfit(ord, env[-1, ], permutations = 9), '`data` has 23 rows, but `ord` has 24 samples')
  named <- ord
  rownames(named) <- paste0('S', rownames(ord))
  envn <- env
  rownames(envn) <- rownames(named)
  expect_error(envfit(named, envn[-1, ], permutations = 9), "1 sample in `ord` has no row in `data`: 'S18'")
  expect_error(envfit(named, envn[-1:-2, ], permutations = 9), '2 samples in `ord` have no row in `data`')



  # Count tables: dense, sparse (Matrix, slam), samples in columns. ====

  cnt <- ex$counts
  r1  <- envfit(ord, counts = cnt, permutations = 99, cpus = 1)

  expect_identical(r1$variable, colnames(cnt))
  expect_identical(rownames(r1), colnames(cnt))
  expect_identical(dim(attr(r1, 'r2_perm')), c(99L, ncol(cnt)))

  # NA in counts is an error, also where 'binary' would have made it a presence.
  cna <- cnt
  cna[2, 2] <- NA
  for (nm in c('none', 'binary', 'percent'))
    expect_error(envfit(ord, counts = cna, norm = nm, permutations = 9, cpus = 1), '`counts` must not contain NA', info = nm)

  # A constant taxon: r2 = 0, no arrow, p = NA.
  const <- which(apply(cnt, 2, function (x) length(unique(x)) == 1))
  expect_true(length(const) > 0)
  expect_true(all(r1$r2[const] == 0 & is.na(r1$p_value[const]) & is.na(r1$z[const])))
  expect_true(all(as.matrix(r1[const, 3:5]) == 0))

  # Same as a data frame of the same columns.
  expect_equal(
    unname(as.matrix(r1[3:7])),
    unname(as.matrix(envfit(ord, as.data.frame(cnt), permutations = 99, cpus = 1)[3:7])),
    tolerance = 1e-12 )

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(cnt, sparse = TRUE)
    expect_equal(envfit(ord, counts = sp, permutations = 99, cpus = 1), r1, tolerance = 1e-12)
    expect_equal(envfit(ord, counts = Matrix::t(sp), margin = 2, permutations = 99, cpus = 1), r1, tolerance = 1e-12)
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(cnt)
    expect_equal(envfit(ord, counts = st, permutations = 99, cpus = 1), r1, tolerance = 1e-12)
  }

  # Samples in other orders, and extra samples, are matched by name.
  cnt2 <- rbind(cnt[24:1, ], extra = seq_len(ncol(cnt)))
  expect_equal(envfit(ord, counts = cnt2, permutations = 99, cpus = 1), r1, tolerance = 1e-12)
  expect_equal(envfit(ord, counts = t(cnt2), margin = 2, permutations = 99, cpus = 1), r1, tolerance = 1e-12)

  # norm = 'percent' equals fitting the proportions.
  expect_equal(
    envfit(ord, counts = cnt, norm = 'percent', permutations = 99, cpus = 1),
    envfit(ord, counts = cnt / rowSums(cnt), permutations = 99, cpus = 1), tolerance = 1e-10 )

  # norm = 'clr': zeros are log(pseudocount) - mean log.
  clr <- log(cnt + 1)
  clr <- clr - rowMeans(clr)
  expect_equal(
    envfit(ord, counts = cnt, norm = 'clr', pseudocount = 1, permutations = 99, cpus = 1)[-2],
    envfit(ord, counts = clr, permutations = 99, cpus = 1)[-2], tolerance = 1e-10 )

  expect_error(envfit(ord, env, norm = 'percent'), 'only apply to `counts`')
  expect_error(envfit(ord, env, margin = 2), 'only apply to `counts`')
  expect_error(envfit(ord, env, pseudocount = 1), 'only apply to `counts`')
  cntn <- cnt
  rownames(cntn) <- rownames(named)
  expect_error(envfit(named, counts = cntn[-1, ]), 'has no sample in `counts`')



  # data and counts together: metadata rows, then feature rows, with the
  # same permutations as each alone.
  both <- envfit(ord, env, cnt, permutations = 99, cpus = 1)
  md   <- envfit(ord, env, permutations = 99, cpus = 1)
  expect_identical(nrow(both), nrow(md) + nrow(r1))
  expect_identical(rownames(both), c(rownames(md), rownames(r1)))
  expect_equal(unname(as.matrix(both[3:8])), unname(rbind(as.matrix(md[3:8]), as.matrix(r1[3:8]))))
  expect_identical(
    attr(both, 'r2_perm'), cbind(attr(md, 'r2_perm'), attr(r1, 'r2_perm')) )
  expect_identical(attr(both, 'n_perm'), 99L)

  # na.rm with both: the samples left out of `data` are left out of
  # `counts` too.
  bna <- envfit(ord, ex$env_na, cnt, na.rm = TRUE, permutations = 9, cpus = 1)
  ok  <- stats::complete.cases(ex$env_na)
  expect_equal(
    unname(as.matrix(bna[bna$variable %in% colnames(cnt), 3:6])),
    unname(as.matrix(envfit(ord[ok, ], counts = cnt[ok, ], permutations = 9, cpus = 1)[3:6])) )



  # Deterministic: the same for any number of cpus, and from run to run. ====

  r_a <- envfit(ord, env, permutations = 199, seed = 4, cpus = 1)
  expect_identical(r_a, envfit(ord, env, permutations = 199, seed = 4, cpus = 1))
  expect_identical(r_a, envfit(ord, env, permutations = 199, seed = 4, cpus = 2))
  expect_false(identical(
    attr(r_a, 'r2_perm'),
    attr(envfit(ord, env, permutations = 199, seed = 5, cpus = 1), 'r2_perm') ))

  if (at_home()) {
    set.seed(2)
    big  <- matrix(rpois(300 * 400, 0.4), 300, 400)
    bord <- matrix(rnorm(600), 300)
    b1 <- envfit(bord, counts = big, permutations = 999, cpus = 1)
    for (cp in c(3, 4, 7))
      expect_identical(b1, envfit(bord, counts = big, permutations = 999, cpus = cp), info = cp)
    bf <- data.frame(g = sample(letters[1:5], 300, TRUE), x = rnorm(300))
    expect_identical(
      envfit(bord, bf, permutations = 999, cpus = 1),
      envfit(bord, bf, permutations = 999, cpus = 4) )
  }



  # A supplied permutation matrix: row k permutes the samples, as
  # vegan's P[perm[k, ], ]. ====

  set.seed(7)
  pm <- t(replicate(19, sample(24)))
  r_m <- envfit(ord, env, permutations = pm, cpus = 1)
  x   <- scale(ord, scale = FALSE)
  for (k in 1:3) {
    v  <- env$N[pm[k, ]]
    r2 <- summary(stats::lm(v ~ x))$r.squared
    expect_equal(unname(attr(r_m, 'r2_perm')[k, 'N']), r2, tolerance = 1e-10)
  }



  # Strata: permutations stay within blocks. ====

  blk <- rep(1:4, each = 6)
  s1  <- envfit(ord, env, strata = blk, permutations = 99, cpus = 1)
  envb <- cbind(env, blk = blk)
  expect_identical(s1, envfit(ord, envb, strata = 'blk', permutations = 99, cpus = 1))

  # A factor constant within strata is never moved, so its r2 never changes.
  envc <- env
  envc$blkf <- factor(blk)
  sc <- envfit(ord, envc, strata = blk, permutations = 99, cpus = 1)
  expect_true(all(abs(attr(sc, 'r2_perm')[, 'blkf'] - sc$r2[sc$variable == 'blkf'][1]) < 1e-12))
  expect_identical(sc$p_value[sc$variable == 'blkf'][1], 1)

  expect_error(envfit(ord, env, strata = blk[-1]), 'one value per row')
  expect_error(envfit(ord, env, strata = 'nope'), 'no column')
  expect_error(envfit(ord, env, strata = blk, permutations = pm), 'not both')



  # Missing values. ====

  expect_error(envfit(ord, ex$env_na), 'na.rm')
  rn <- envfit(ord, ex$env_na, na.rm = TRUE, permutations = 9, cpus = 1)
  expect_identical(nrow(rn), nrow(res))
  ord_na <- ord
  ord_na[5, 1] <- NA
  expect_error(envfit(ord_na, counts = cnt), 'na.rm')
  expect_identical(
    envfit(ord_na, counts = cnt, na.rm = TRUE, permutations = 9, cpus = 1),
    envfit(ord[-5, ], counts = cnt, permutations = 9, cpus = 1) )



  # Edge cases. ====

  # Rank-deficient coordinates: the aliased axis gets 0.
  rd <- cbind(ord[, 1:2], sum = ord[, 1] + ord[, 2])
  r  <- envfit(rd, env[1:3], permutations = 9, cpus = 1)
  expect_true(all(r$sum == 0))
  expect_equal(r$r2, envfit(ord[, 1:2], env[1:3], permutations = 9, cpus = 1)$r2, tolerance = 1e-10)

  # A one-level factor: r2 = 0, p = NA.
  one <- envfit(ord, data.frame(f = rep('a', 24), row.names = rownames(ord)), permutations = 9, cpus = 1)
  expect_identical(one$r2, 0)
  expect_true(is.na(one$p_value))

  expect_error(envfit(ord), 'or both')
  expect_error(envfit(matrix(1, 5, 2), data.frame(x = 1:5)), 'same coordinates')
  expect_error(envfit('a', env), 'numeric matrix')
  expect_error(envfit(ord, env, axes = 4), 'axes')
  expect_error(envfit(ord, env, axes = 'PCoA9'), 'not in')
  expect_error(envfit(ord, env, permutations = -1), 'permutations')
  expect_error(envfit(ord, env, na.rm = NA), 'na.rm')
  expect_error(envfit(ord, data.frame(x = I(list(1)))[rep(1, 24), , drop = FALSE]))
  expect_error(envfit(ord[1, , drop = FALSE], env[1, ]), 'At least 2')

#})
