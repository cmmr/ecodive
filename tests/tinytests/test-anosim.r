#test_that("anosim", {

  ex   <- readRDS('expected/anosim.rds')
  d    <- bray(ex$counts)
  meta <- ex$meta
  n    <- nrow(meta)


  # ANOSIM written out with base R, as in vegan::anosim(): sample i of
  # permutation `take` gets the group of sample take[i].
  anosim_r <- function (dm, groups, perms = NULL) {
    x  <- as.vector(dm)
    rk <- rank(x)
    N  <- attr(dm, 'Size')
    ij <- which(lower.tri(diag(N)), arr.ind = TRUE)
    one <- function (g) {
      within <- g[ij[, 1]] == g[ij[, 2]]
      (mean(rk[!within]) - mean(rk[within])) / (length(x) / 2)
    }
    g <- as.integer(factor(groups))
    c(one(g), if (!is.null(perms)) apply(perms, 1L, function (p) one(g[p])))
  }



  # Same statistics as the reference, given the same permutations: stored
  # results, so at home only. Every p-value matches exactly; R and every
  # permuted R to the last bit or so (R's mean() rounds in long double). ====

  if (at_home()) {
    for (nm in names(ex$cases)) {
      cs  <- ex$cases[[nm]]
      res <- anosim(cs$dm, cs$groups, permutations = ex$perms, cpus = 1)
      expect_equal(res$statistic[['R']], cs$statistic, tolerance = 1e-14, info = nm)
      expect_identical(res$p_value, cs$signif, info = nm)
      expect_equal(attr(res, 'r_perm'), cs$perm, tolerance = 1e-14, info = nm)
      rk <- res$ranks
      expect_identical(rownames(rk), cs$ranks$class, info = nm)
      expect_equal(rk$pairs[!is.na(cs$ranks$pairs)], cs$ranks$pairs[!is.na(cs$ranks$pairs)], info = nm)
      expect_equal(rk$mean_rank, cs$ranks$mean_rank, tolerance = 1e-14, info = nm)
    }

    # The same p-value as the reference with strata, up to sampling error.
    res <- anosim(d, meta$site, permutations = 9999, strata = meta$block, cpus = 1)
    expect_true(abs(res$p_value - ex$strata_p) < 0.02)
  }



  # Same statistics as ANOSIM written out in base R, on stored inputs. ====

  perms <- ex$perms[1:49, ]
  for (dm in list(d, round(d, 1), round(d * 0 + 0.3, 1))) {
    for (g in list(meta$trt, meta$site, meta$uneq)) {
      ref <- anosim_r(dm, g, perms)
      res <- anosim(dm, g, permutations = perms, cpus = 1)
      expect_equal(res$statistic[['R']], ref[[1]], tolerance = 1e-14)
      expect_equal(attr(res, 'r_perm'), ref[-1], tolerance = 1e-14)
      expect_identical(
        res$p_value,
        (1 + sum(ref[-1] >= ref[[1]] - sqrt(.Machine$double.eps))) / 50 )
    }
  }

  # More samples than one tile of rows (512) in C, with ties; an odd number,
  # so that columns are left over from the pairs.
  set.seed(1)
  nb  <- 1101
  db  <- stats::dist(round(matrix(stats::runif(nb * 2), nb), 1))
  gb  <- sample(c('a', 'b', 'c'), nb, TRUE)
  pb  <- t(replicate(3, sample(nb)))
  ref <- anosim_r(db, gb, pb)
  res <- anosim(db, gb, permutations = pb, cpus = 2)
  expect_true(max(abs(c(res$statistic[['R']], attr(res, 'r_perm')) - ref)) < 1e-12)

  # Mean ranks of each class.
  rk  <- rank(as.vector(d))
  ij  <- which(lower.tri(diag(n)), arr.ind = TRUE)
  gi  <- as.character(meta$trt[ij[, 1]])
  gj  <- as.character(meta$trt[ij[, 2]])
  cls <- factor(ifelse(gi == gj, gi, 'Between'), c('Between', levels(meta$trt)))
  res <- anosim(d, meta$trt, permutations = 0)$ranks
  expect_identical(rownames(res), levels(cls))
  expect_equal(res$pairs, as.vector(table(cls)))
  expect_equal(res$mean_rank, as.vector(tapply(rk, cls, mean)), tolerance = 1e-14)

  # Zero and negative zero are tied.
  z0 <- d
  z0[1:5]   <- 0
  z0[6:10]  <- -0
  z0[11:20] <- 0.5
  expect_equal(anosim(z0, meta$trt, permutations = 0)$statistic[['R']], anosim_r(z0, meta$trt)[[1]], tolerance = 1e-14)

  # R is 1 when every within distance is below every between distance, and
  # -1 the other way round.
  g2  <- rep(1:2, each = 4)
  sep <- stats::dist(g2 + seq(0, 0.01, length.out = 8))
  expect_equal(anosim(sep, g2, permutations = 0)$statistic[['R']], 1)
  expect_equal(anosim(stats::dist(c(0, 5, 2.5)), c(1, 1, 2), permutations = 0)$statistic[['R']], -1)



  # Result format. ====

  res <- anosim(d, 'trt', meta, permutations = 99, cpus = 1)

  expect_inherits(res, 'ecodive_anosim')
  expect_inherits(res, 'ecodive_test')
  expect_identical(names(res), c('statistic', 'p_value', 'z', 'n_perm', 'ranks', 'permutation'))
  expect_identical(names(res$statistic), 'R')
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  expect_identical(length(attr(res, 'r_perm')), 99L)
  expect_identical(attr(res, 'title'), 'Analysis of similarities (ANOSIM)')
  expect_identical(names(res$ranks), c('pairs', 'mean_rank'))
  expect_true(res$p_value >= 1 / 100)
  expect_equal(res$z, (res$statistic[['R']] - mean(attr(res, 'r_perm'))) / stats::sd(attr(res, 'r_perm')))
  expect_identical(anosim(d, meta$trt, permutations = ex$perms)$permutation, 'supplied matrix')
  expect_identical(anosim(d, meta$trt, permutations = 9, strata = meta$block)$permutation, 'within strata')

  # The printout, with the mean ranks; one row of the shared columns.
  out <- capture.output(print(res))
  expect_true(any(grepl('Analysis of similarities (ANOSIM)', out, fixed = TRUE)))
  expect_true(any(grepl('R = ', out, fixed = TRUE)))
  expect_true(any(grepl('Mean rank', out, fixed = TRUE)))
  expect_true(any(grepl('^Between', out)))
  df <- as.data.frame(res)
  expect_identical(df$test, 'anosim')
  expect_identical(df$stat_name, 'R')
  expect_identical(df$statistic, res$statistic[['R']])
  expect_identical(df$p_value, res$p_value)

  # No permutations: no p-value or z.
  res <- anosim(d, meta$trt, permutations = 0)
  expect_true(is.na(res$p_value) && is.na(res$z))
  expect_identical(res$n_perm, 0L)
  expect_identical(attr(res, 'r_perm'), numeric(0))
  expect_true(is.na(anosim(d, meta$trt, permutations = 1)$z))

  # Every distance the same: R is 0, every permutation ties.
  res <- anosim(d * 0 + 1, meta$trt, permutations = 19)
  expect_identical(res$statistic, c(R = 0))
  expect_identical(res$p_value, 1)
  expect_true(is.na(res$z))



  # Deterministic: the same for any number of cpus, and from run to run. ====

  r1 <- anosim(d, meta$trt, permutations = 199, seed = 7, cpus = 1)
  r2 <- anosim(d, meta$trt, permutations = 199, seed = 7, cpus = 2)
  expect_identical(r1, r2)
  expect_identical(r1, anosim(d, meta$trt, permutations = 199, seed = 7, cpus = 1))

  # Large enough for run_parallel() to start threads.
  big <- bray(big_mtx)
  bg  <- rep(c('a', 'b', 'c', 'd'), 26)
  r1  <- anosim(big, bg, permutations = 999, cpus = 1)
  r2  <- anosim(big, bg, permutations = 999, cpus = 2)
  expect_identical(r1, r2)
  if (at_home()) {
    expect_identical(r1, anosim(big, bg, permutations = 999, cpus = 4))
    expect_identical(r1, anosim(big, bg, permutations = 999, cpus = 7))
  }

  # The seed changes the permutations, not the statistic.
  r2 <- anosim(d, meta$trt, permutations = 99, seed = 2, cpus = 1)
  r3 <- anosim(d, meta$trt, permutations = 99, seed = 3, cpus = 1)
  expect_identical(r2$statistic[['R']], r3$statistic[['R']])
  expect_false(identical(attr(r2, 'r_perm'), attr(r3, 'r_perm')))



  # Strata. ====

  # Groups that are constant within strata can't change: every permuted R
  # equals the observed one.
  res <- anosim(d, 'site', meta, permutations = 99, strata = 'site', cpus = 1)
  expect_identical(res$p_value, 1)
  expect_identical(attr(res, 'r_perm'), rep(res$statistic[['R']], 99))
  expect_identical(
    anosim(d, 'trt', meta, permutations = 19, strata = 'block', cpus = 1),
    anosim(d, meta$trt, permutations = 19, strata = meta$block, cpus = 1) )

  expect_error(anosim(d, meta$trt, permutations = ex$perms, strata = meta$block), 'not both')
  expect_error(anosim(d, meta$trt, strata = 1:3), 'strata')
  expect_error(anosim(d, 'trt', meta, strata = 1:3), 'strata')
  expect_error(anosim(d, meta$trt, strata = c(NA, meta$block[-1])), 'NA')



  # Groups and metadata. ====

  ref <- anosim(d, meta$trt, permutations = 19, cpus = 1)

  # A column of `data`, matched by name ...
  shuffled <- meta[rev(seq_len(n)), ]
  expect_identical(anosim(d, 'trt', shuffled, permutations = 19, cpus = 1), ref)

  # ... or by position, without row names.
  unnamed <- meta
  rownames(unnamed) <- NULL
  expect_identical(anosim(d, 'trt', unnamed, permutations = 19, cpus = 1), ref)

  # A vector parallel to the rows of `data`.
  expect_identical(anosim(d, shuffled$trt, shuffled, permutations = 19, cpus = 1), ref)

  # A named vector, matched by name; character, logical and numeric groups.
  expect_identical(anosim(d, setNames(meta$trt, rownames(meta))[n:1], permutations = 19, cpus = 1), ref)
  expect_identical(anosim(d, as.character(meta$trt), permutations = 19, cpus = 1)$statistic[['R']], ref$statistic[['R']])
  expect_identical(anosim(d, as.integer(meta$trt) * 10, permutations = 19, cpus = 1)$statistic[['R']], ref$statistic[['R']])
  expect_identical(
    anosim(d, meta$site == 'A', permutations = 19, cpus = 1)$statistic[['R']],
    anosim(d, meta$site, permutations = 19, cpus = 1)$statistic[['R']] )

  # A symmetric matrix works like a dist; a data matrix like a data frame.
  expect_identical(anosim(as.matrix(d), meta$trt, permutations = 19, cpus = 1), ref)
  mm <- as.matrix(meta[c('trt', 'site')])
  expect_identical(anosim(d, 'site', mm, permutations = 19, cpus = 1)$statistic[['R']], anosim(d, meta$site, permutations = 0)$statistic[['R']])

  # Rows for other samples are ignored, even with NAs, and so are groups
  # found only in them.
  extra <- rbind(meta, meta[1, ])
  rownames(extra)[nrow(extra)] <- 'X'
  extra$trt <- factor(extra$trt, c(levels(meta$trt), 'only_x'))
  extra$trt[nrow(extra)] <- NA
  res <- anosim(d, 'trt', extra, permutations = 19, cpus = 1)
  expect_identical(rownames(res$ranks), c('Between', levels(meta$trt)))
  extra$trt[nrow(extra)] <- 'only_x'
  expect_identical(anosim(d, 'trt', extra, permutations = 19, cpus = 1), ref)

  # Many groups don't make a big model matrix.
  pairs2 <- rep(seq_len(n / 2), each = 2)
  res <- anosim(d, pairs2, permutations = 9, cpus = 1)
  expect_equal(res$statistic[['R']], anosim_r(d, pairs2)[[1]], tolerance = 1e-14)

  expect_error(anosim(d), '`groups` is required')
  expect_error(anosim(d, 'trt'), '`data` is not given')
  expect_error(anosim(d, meta$trt[-1]), '29 values')
  expect_error(anosim(d, setNames(meta$trt, rownames(meta))[-1]), "1 sample in `dm` has no value in `groups`: 'S01'")
  expect_error(anosim(d, 'nothere', meta), 'name of a column')
  expect_error(anosim(d, meta$trt[-1], meta), 'one value per row')
  expect_error(anosim(d, meta[-1, ], 'trt'), 'vector or factor')
  expect_error(anosim(d, 'trt', meta[-1, ]), "1 sample in `dm` has no row in `data`: 'S01'")
  expect_error(anosim(d, 'trt', list(trt = meta$trt)), 'data.frame')
  expect_error(anosim(d, c(NA, meta$trt[-1])), 'NA')
  expect_error(anosim(d, rep('a', n)), 'at least 2 different')
  expect_error(anosim(d, seq_len(n)), 'replicates')



  # Other arguments. ====

  expect_error(anosim(d, meta$trt, permutations = -1), 'permutations')
  expect_error(anosim(d, meta$trt, permutations = 1.5), 'permutations')
  expect_error(anosim(d, meta$trt, permutations = ex$perms[, -1]), '30 columns')
  bad <- ex$perms
  bad[2, 1:2] <- 1L
  expect_error(anosim(d, meta$trt, permutations = bad), 'Row 2')
  expect_error(anosim(d, meta$trt, seed = 1.5), 'seed')
  expect_error(anosim(d, meta$trt, cpus = 0), 'cpus')

  bad_d <- d
  bad_d[1] <- NA
  expect_error(anosim(bad_d, meta$trt), 'NA')
  bad_d[1] <- -1
  expect_error(anosim(bad_d, meta$trt), 'negative')
  expect_error(anosim(stats::dist(1), 1), 'at least 2 samples')
  expect_error(anosim(d ~ trt, meta), 'formula')

  # Two samples in one group, one in another.
  tiny <- stats::dist(c(0, 1, 5))
  expect_equal(anosim(tiny, c(1, 1, 2), permutations = 0)$statistic[['R']], 1)

#})
