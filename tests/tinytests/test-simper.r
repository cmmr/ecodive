#test_that("simper", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  cols <- c('pair', 'taxon', 'average', 'sd', 'ratio', 'mean_a', 'mean_b',
            'cumsum', 'p_value', 'p_adjusted')

  set.seed(4)
  cnt <- matrix(rpois(18 * 15, 3) * rbinom(18 * 15, 1, 0.5), 18, 15,
    dimnames = list(paste0('S', 1:18), paste0('T', 1:15)) )
  cnt[, 1] <- cnt[, 1] + 1L
  cnt[, 'T7'] <- 0L
  grp <- rep(c('a', 'b', 'c'), c(7, 6, 5))
  cnt[grp == 'c', 2] <- cnt[grp == 'c', 2] + 6L


  # Stored results from vegan::simper() and summary() of it, on each pair's
  # samples with the induced permutations. p-values follow Clarke's "at
  # least as large" rule, where vegan counts only larger permuted averages:
  # they are checked against a direct computation of the rule on the same
  # permutations, and are never smaller than vegan's. ====

  # Clarke's p-values for one pair, directly: each sample pair's taxon
  # contributions, and for each permutation (permuted sample a takes the
  # label lab[p[a]]) the mean over the pairs whose labels differ.
  clarke_p <- function (x, lab, perms) {
    m  <- nrow(x)
    ij <- which(lower.tri(diag(m)), arr.ind = TRUE)
    tt <- rowSums(x)
    cc <- abs(x[ij[, 1], , drop = FALSE] - x[ij[, 2], , drop = FALSE]) /
          (tt[ij[, 1]] + tt[ij[, 2]])
    avg <- function (l) colMeans(cc[l[ij[, 1]] != l[ij[, 2]], , drop = FALSE])
    obs <- avg(lab)
    cnt <- numeric(ncol(x))
    for (b in seq_len(nrow(perms))) {
      pa  <- avg(lab[perms[b, ]])
      cnt <- cnt + (pa >= obs - 1e-7 * obs)
    }
    p <- (cnt + 1) / (nrow(perms) + 1)
    p[obs == 0] <- NA_real_
    stats::setNames(p, colnames(x))
  }

  if (at_home()) {

    ref <- readRDS(file.path('expected', 'simper.rds'))

    for (cs in ref$cases) {

      lbl <- cs$name
      res <- simper(cs$counts, cs$groups, norm = cs$norm,
                    permutations = cs$perms, cpus = 1)
      ex  <- cs$expected

      expect_inherits(res, c('ecodive_simper', 'data.frame'), info = lbl)
      expect_identical(names(res), cols, info = lbl)
      expect_identical(unique(res$pair), unique(ex$pair), info = lbl)
      expect_identical(attr(res, 'permutation'), 'supplied matrix', info = lbl)
      expect_identical(attr(res, 'n_perm'), nrow(cs$perms), info = lbl)

      for (pr in unique(ex$pair)) {
        r <- res[res$pair == pr, ]
        e <- ex[ex$pair == pr, ]
        info <- paste(lbl, pr)

        # The same order of taxa, by decreasing average (ties in the order
        # of the columns).
        expect_identical(r$taxon, e$taxon, info = info)

        for (v in c('average', 'sd', 'ratio', 'mean_a', 'mean_b', 'cumsum'))
          expect_equal(r[[v]], e[[v]], tolerance = 1e-10, info = paste(info, v))

        ab  <- strsplit(pr, '-', fixed = TRUE)[[1]]
        idx <- which(as.character(cs$groups) %in% ab)
        x   <- as.matrix(cs$counts)[idx, , drop = FALSE]
        if (identical(cs$norm, 'percent')) x <- x / rowSums(x)
        cp  <- clarke_p(x, as.character(cs$groups)[idx],
                        ecodive:::induced_perms(cs$perms, idx))
        expect_equal(r$p_value, unname(cp[r$taxon]), info = info)
        ok <- !is.na(e$p_value)
        expect_identical(is.na(r$p_value), is.na(e$p_value), info = info)
        expect_true(all(r$p_value[ok] >= e$p_value[ok]), info = info)
        expect_equal(attr(res, 'overall')[[pr]], e$overall[[1]], tolerance = 1e-10, info = info)
        expect_equal(r$p_adjusted, p.adjust(r$p_value, 'holm'), info = info)
      }

      # Results are identical for any number of threads.
      for (cpus in 2:max_cpus)
        expect_identical(
          simper(cs$counts, cs$groups, norm = cs$norm,
                 permutations = cs$perms, cpus = cpus),
          res, info = paste(lbl, cpus) )

      # The blocked permutations stay within the blocks; `strata` gives
      # random ones within them too.
      if (!is.null(cs$blocks)) {
        st <- simper(cs$counts, cs$groups, strata = cs$blocks,
                     permutations = 99, cpus = 1)
        expect_identical(attr(st, 'permutation'), 'within strata', info = lbl)
        expect_equal(st$average, res$average, info = lbl)
        expect_identical(simper(cs$counts, cs$groups, strata = cs$blocks,
                                permutations = 99, cpus = max_cpus), st, info = lbl)
      }
    }
  }



  # The averages sum to the mean between-group Bray-Curtis
  # dissimilarity. ====

  res <- simper(cnt, grp, permutations = 99, cpus = 1)
  bm  <- as.matrix(bray(cnt, cpus = 1))

  for (pr in c('a-b', 'a-c', 'b-c')) {
    ab <- strsplit(pr, '-')[[1]]
    r  <- res[res$pair == pr, ]
    mb <- mean(bm[grp == ab[[1]], grp == ab[[2]]])
    expect_equal(sum(r$average), mb, tolerance = 1e-12, info = pr)
    expect_equal(attr(res, 'overall')[[pr]], sum(r$average), info = pr)
    expect_equal(r$cumsum, cumsum(r$average) / sum(r$average), info = pr)
    expect_equal(r$cumsum[[nrow(r)]], 1, info = pr)
    expect_false(is.unsorted(rev(r$average)), info = pr)
    expect_identical(nrow(r), ncol(cnt), info = pr)
    expect_equal(r$mean_a, unname(colMeans(cnt[grp == ab[[1]], r$taxon])), info = pr)
    expect_equal(r$mean_b, unname(colMeans(cnt[grp == ab[[2]], r$taxon])), info = pr)

    # The breakdown straight from its definition.
    i  <- which(grp == ab[[1]])
    j  <- which(grp == ab[[2]])
    cc <- do.call(rbind, lapply(i, function (a) do.call(rbind, lapply(j, function (b)
      abs(cnt[a, ] - cnt[b, ]) / sum(cnt[a, ] + cnt[b, ]) ))))
    expect_equal(r$average, unname(colMeans(cc)[r$taxon]), tolerance = 1e-12, info = pr)
    expect_equal(r$sd, unname(apply(cc, 2, sd)[r$taxon]), tolerance = 1e-12, info = pr)
    expect_equal(r$ratio, r$average / r$sd, info = pr)
  }

  expect_identical(attr(res, 'n_perm'), 99L)
  expect_identical(attr(res, 'p_adjust'), 'holm')
  expect_identical(attr(res, 'permutation'), 'free')
  expect_identical(unique(res$pair), c('a-b', 'a-c', 'b-c'))
  expect_true(all(res$p_value > 0 & res$p_value <= 1, na.rm = TRUE))


  # A taxon absent everywhere: average 0, ratio NaN, and no p-value.
  t7 <- res[res$taxon == 'T7', ]
  expect_identical(t7$average, rep(0, 3))
  expect_identical(t7$sd, rep(0, 3))
  expect_true(all(is.nan(t7$ratio)))
  expect_true(all(is.na(t7$p_value)))
  expect_true(all(is.na(t7$p_adjusted)))

  # The taxon raised in group c ranks first against it, with a small p.
  expect_identical(res$taxon[res$pair == 'a-c'][[1]], 'T2')
  expect_true(res$p_value[res$pair == 'a-c'][[1]] <= 0.05)


  # p_adjust: within each pair, leaving out NAs.
  for (m in c('none', 'BH', 'bonferroni')) {
    r2 <- simper(cnt, grp, p_adjust = m, permutations = 99, cpus = 1)
    expect_identical(r2$p_value, res$p_value, info = m)
    expect_identical(attr(r2, 'p_adjust'), m, info = m)
    for (pr in unique(r2$pair))
      expect_equal(r2$p_adjusted[r2$pair == pr],
                   p.adjust(r2$p_value[r2$pair == pr], m), info = m)
  }
  expect_identical(attr(simper(cnt, grp, p_adjust = 'bonf', permutations = 0, cpus = 1), 'p_adjust'), 'bonferroni')


  # Seeds and threads: reproducible, and the same for any `cpus`. ====

  expect_identical(simper(cnt, grp, permutations = 199, seed = 3, cpus = 1),
                   simper(cnt, grp, permutations = 199, seed = 3, cpus = max_cpus))
  expect_false(identical(
    simper(cnt, grp, permutations = 199, seed = 3, cpus = 1)$p_value,
    simper(cnt, grp, permutations = 199, seed = 4, cpus = 1)$p_value ))

  st <- rep(1:3, 6)
  expect_identical(simper(cnt, grp, strata = st, permutations = 199, cpus = 1),
                   simper(cnt, grp, strata = st, permutations = 199, cpus = max_cpus))

  # A larger table, so that the observed sums split over several tasks.
  set.seed(5)
  big <- matrix(rpois(60 * 700, 1) * rbinom(60 * 700, 1, 0.2), 60, 700)
  big[, 1] <- big[, 1] + 1L
  bg  <- rep(c('u', 'v'), 30)
  b1  <- simper(big, bg, permutations = 49, cpus = 1)
  expect_identical(simper(big, bg, permutations = 49, cpus = max_cpus), b1)
  expect_identical(b1$taxon[1:3], as.character(order(-b1$average[order(as.integer(b1$taxon))]))[1:3])
  expect_equal(sum(b1$average), mean(as.matrix(bray(big, cpus = 1))[bg == 'u', bg == 'v']), tolerance = 1e-12)


  # Strata that hold one group each: no permutation changes the groups, so
  # every permuted average ties the observed one, and every p-value is 1.
  s2 <- simper(cnt[1:13, ], grp[1:13], strata = grp[1:13], permutations = 19, cpus = 1)
  expect_true(all(s2$p_value == 1 | is.na(s2$p_value)))

  # Two samples per group: a third of all permutations only swap samples
  # within the groups, and reproduce the observed averages, so no p-value
  # is below about 1/3 (vegan's rule would give the top taxon 1/(n+1)).
  s4 <- simper(cnt[c(1, 2, 8, 9), ],
               grp[c(1, 2, 8, 9)], permutations = 999, cpus = 1)
  expect_true(min(s4$p_value, na.rm = TRUE) > 0.25)

  # A supplied matrix is cut down to each pair, as in the pairwise tests.
  set.seed(6)
  pm <- t(replicate(30, sample(18)))
  sm <- simper(cnt, grp, permutations = pm, cpus = 1)
  for (pr in c('a-b', 'b-c')) {
    ab  <- strsplit(pr, '-')[[1]]
    idx <- which(grp %in% ab)
    one <- simper(cnt[idx, ], grp[idx], permutations = ecodive:::induced_perms(pm, idx), cpus = 1)
    expect_identical(sm$p_value[sm$pair == pr], one$p_value, info = pr)
  }



  # Inputs: sparse tables, margin, data, names. ====

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(cnt, sparse = TRUE)
    expect_equal(simper(sp, grp, permutations = 99, cpus = 1), res, tolerance = 1e-12)
    expect_equal(simper(Matrix::t(sp), grp, margin = 2, permutations = 99, cpus = 1), res, tolerance = 1e-12)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_equal(simper(slam::as.simple_triplet_matrix(cnt), grp, permutations = 99, cpus = 1), res, tolerance = 1e-12)

  expect_equal(simper(t(cnt), grp, margin = 2, permutations = 99, cpus = 1), res)
  storage.mode(cnt) <- 'double'
  expect_equal(simper(cnt, grp, permutations = 99, cpus = 1), res)

  meta <- data.frame(row.names = rev(rownames(cnt)), g = rev(grp), s = rev(st))
  expect_identical(simper(cnt, 'g', meta, permutations = 99, cpus = 1), res)
  expect_identical(simper(cnt, 'g', meta, strata = 's', permutations = 99, cpus = 1),
                   simper(cnt, grp, strata = st, permutations = 99, cpus = 1))
  expect_identical(simper(cnt, factor(grp, levels = c('b', 'a', 'c', 'd')), permutations = 0, cpus = 1)$pair[[1]], 'b-a')

  un <- cnt
  dimnames(un) <- NULL
  expect_identical(simper(un, grp, permutations = 0, cpus = 1)$taxon,
                   sub('T', '', simper(cnt, grp, permutations = 0, cpus = 1)$taxon))

  # norm: the same as normalizing first.
  pc <- cnt / rowSums(cnt)
  expect_equal(simper(cnt, grp, norm = 'percent', permutations = 99, cpus = 1),
               simper(pc, grp, permutations = 99, cpus = 1), tolerance = 1e-12)
  hc <- sqrt(pc)
  expect_equal(simper(cnt, grp, norm = 'hellinger', permutations = 99, cpus = 1),
               simper(hc, grp, permutations = 99, cpus = 1), tolerance = 1e-12)
  bc <- (cnt > 0) + 0
  expect_equal(simper(cnt, grp, norm = 'binary', permutations = 99, cpus = 1),
               simper(bc, grp, permutations = 99, cpus = 1), tolerance = 1e-12)



  # Edge cases. ====

  # No permutations: NA p-values.
  r0 <- simper(cnt, grp, permutations = 0, cpus = 1)
  expect_true(all(is.na(r0$p_value)))
  expect_true(all(is.na(r0$p_adjusted)))
  expect_identical(r0$average, res$average)
  expect_identical(attr(r0, 'n_perm'), 0L)

  # Two groups: one pair.
  r2 <- simper(cnt[1:13, ], grp[1:13], permutations = 9, cpus = 1)
  expect_identical(unique(r2$pair), 'a-b')
  expect_identical(nrow(r2), ncol(cnt))

  # A group of one sample.
  g1 <- c(rep('a', 17), 'z')
  r1 <- simper(cnt, g1, permutations = 99, cpus = 1)
  expect_identical(unique(r1$pair), 'a-z')
  expect_true(all(is.finite(r1$sd)))
  expect_true(all(r1$p_value >= 1 / 100 | is.na(r1$p_value)))

  # Two groups of one sample: no sd; their Bray-Curtis broken down.
  r11 <- simper(cnt[1:2, ], c('p', 'q'), permutations = 9, cpus = 1)
  expect_true(all(is.na(r11$sd)))
  expect_true(all(is.na(r11$ratio)))
  expect_equal(sum(r11$average), as.numeric(bray(cnt[1:2, ], cpus = 1)))

  # A taxon with the same abundance everywhere contributes nothing.
  cst <- cnt
  cst[, 'T3'] <- 4
  rc <- simper(cst, grp, permutations = 19, cpus = 1)
  expect_identical(rc$average[rc$taxon == 'T3'], rep(0, 3))
  expect_true(all(is.na(rc$p_value[rc$taxon == 'T3'])))

  # A factor's NA level is a group.
  expect_identical(unique(simper(cnt, addNA(factor(ifelse(grp == 'c', NA, grp))), permutations = 0, cpus = 1)$pair),
                   c('a-b', 'a-NA', 'b-NA'))



  # Errors. ====

  expect_error(simper(cnt), '`groups` is required')
  expect_error(simper(cnt, rep('a', 18)), 'at least two different values')
  expect_error(simper(cnt, grp[-1]), 'groups')
  expect_error(simper(cnt, c(grp[-1], NA)), 'NA')
  expect_error(simper(cnt, grp, norm = 'clr'), "Bray-Curtis can't use")
  expect_error(simper(cnt, grp, norm = 'rclr'), "Bray-Curtis can't use")
  expect_error(simper(cnt, grp, norm = 'nope'))
  expect_error(simper(cnt, grp, p_adjust = 'nope'), 'p.adjust')
  expect_error(simper(cnt, grp, p_adjust = 'h'), 'p.adjust')
  expect_error(simper(cnt, grp, permutations = -1), 'permutations')
  expect_error(simper(cnt, grp, permutations = matrix(1:3, 1)), 'permutations')
  expect_error(simper(cnt, grp, permutations = pm, strata = st), 'not both')
  expect_error(simper(cnt, grp, seed = 1.5), 'seed')
  expect_error(simper(cnt, grp, cpus = 0), 'cpus')

  neg <- cnt; neg[1, 1] <- -1
  expect_error(simper(neg, grp, cpus = 1), 'negative')
  nas <- cnt; nas[1, 1] <- NA
  expect_error(simper(nas, grp, cpus = 1), 'NAs')
  emp <- cnt; emp[c(2, 5), ] <- 0
  expect_error(simper(emp, grp, cpus = 1), 'Samples with no counts: S2, S5')

#})
