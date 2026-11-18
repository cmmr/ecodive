#test_that("pairwise", {

  # pairwise = TRUE in permanova(), anosim() and mrpp(): each pair of groups
  # tested on its own, with adjusted p-values. permdisp()'s pairwise t tests
  # are in test-permdisp.r.

  ex   <- readRDS('expected/pairwise.rds')
  meta <- ex$meta
  dms  <- list(bray = bray(ex$counts), euclidean = stats::dist(ex$counts))
  g    <- meta$grp
  labs <- c('a-b', 'a-c', 'a-d', 'b-c', 'b-d', 'c-d')

  # The same pairs, as each test computes them by hand on the pair's
  # samples. `idx` are the pair's samples.
  pair_dm <- function (d, idx) stats::as.dist(as.matrix(d)[idx, idx])



  # The same results as vegan on each pair's samples, given the
  # permutations the full matrix induces on them: stored results, so at
  # home only. ====

  if (at_home()) {
    for (cs in ex$cases) {

      info <- paste(cs$dist, cs$perms)
      d    <- dms[[cs$dist]]
      pm   <- ex$perms[[cs$perms]]
      tbl  <- cs$table

      res <- permanova(d, g, pairwise = TRUE, permutations = pm, cpus = 1)$pairwise
      expect_identical(rownames(res), labs, info = info)
      expect_identical(res$n, as.integer(tbl$n), info = info)
      expect_equal(res$F,  tbl$F,   tolerance = 1e-10, info = info)
      expect_equal(res$R2, tbl$R2,  tolerance = 1e-10, info = info)
      expect_equal(res$z,  tbl$F_z, tolerance = 1e-8,  info = info)
      expect_identical(res$p_value, tbl$F_p, info = info)

      res <- anosim(d, g, pairwise = TRUE, permutations = pm, cpus = 1)$pairwise
      expect_equal(res$R, tbl$R,   tolerance = 1e-10, info = info)
      expect_equal(res$z, tbl$R_z, tolerance = 1e-8,  info = info)
      expect_identical(res$p_value, tbl$R_p, info = info)

      res <- mrpp(d, g, pairwise = TRUE, permutations = pm, cpus = 1)$pairwise
      expect_equal(res$delta, tbl$delta1,   tolerance = 1e-10, info = info)
      expect_equal(res$A,     tbl$A1,       tolerance = 1e-10, info = info)
      expect_equal(res$z,     tbl$delta1_z, tolerance = 1e-8,  info = info)
      expect_identical(res$p_value, tbl$delta1_p, info = info)

      res <- mrpp(d, g, weights = 'n(n-1)', pairwise = TRUE, permutations = pm, cpus = 1)$pairwise
      expect_equal(res$delta, tbl$delta3, tolerance = 1e-10, info = info)
      expect_equal(res$A,     tbl$A3,     tolerance = 1e-10, info = info)
      expect_identical(res$p_value, tbl$delta3_p, info = info)

      # The induced permutations, as the fixture wrote them with a loop.
      if (!is.null(cs$induced)) for (k in seq_along(labs)) {
        idx <- which(g %in% strsplit(labs[[k]], '-')[[1]])
        expect_identical(ecodive:::induced_perms(pm, idx), unname(cs$induced[[k]]), info = labs[[k]])
      }
    }

    # RVAideMemoire's pairwise.perm.manova(), with its own random
    # permutations: the same p-values, up to sampling error, and the same
    # Holm adjustment.
    res <- permanova(dms$bray, g, pairwise = TRUE, permutations = 9999)$pairwise
    expect_true(all(abs(res$p_value - ex$rva$none) < 0.015))
    expect_true(all(abs(res$p_adjusted - ex$rva$holm) < 0.03))
    expect_equal(unname(stats::p.adjust(ex$rva$none, "holm")), unname(ex$rva$holm))
  }



  # Identities that hold without the references. ====

  d  <- dms$bray
  pm <- ex$perms$free

  # Each pair is the test run by hand on the pair's samples, with its own
  # seed and the permutations the full matrix induces.
  tests <- list(
    permanova = function (...) permanova(...),
    anosim    = function (...) anosim(...),
    mrpp      = function (...) mrpp(...) )

  for (nm in names(tests)) {
    f <- tests[[nm]]
    for (k in c(1, 6)) {
      idx  <- which(g %in% strsplit(labs[[k]], '-')[[1]])
      full <- f(d, g, pairwise = TRUE, permutations = 49, seed = 3, cpus = 1)
      one  <- f(pair_dm(d, idx), g[idx], permutations = 49, seed = ecodive:::pair_seed(3L, k), cpus = 1)
      expect_identical(full$pairwise$p_value[[k]], one$p_value, info = nm)
      expect_identical(full$pairwise$z[[k]], one$z, info = nm)
      expect_equal(full$pairwise[[1]][[k]], one$statistic[[1]], tolerance = 1e-14, info = nm)

      full <- f(d, g, pairwise = TRUE, permutations = pm, cpus = 1)
      one  <- f(pair_dm(d, idx), g[idx], permutations = ecodive:::induced_perms(pm, idx), cpus = 1)
      expect_identical(full$pairwise$p_value[[k]], one$p_value, info = nm)
    }

    # The overall result is unchanged.
    r0 <- f(d, g, permutations = 49, seed = 3, cpus = 1)
    r1 <- f(d, g, pairwise = TRUE, permutations = 49, seed = 3, cpus = 1)
    r1$pairwise <- r1$p_adjust <- r0$pairwise <- r0$p_adjust <- NULL
    expect_identical(r0, r1, info = nm)

    # Pairs have their own permutations: equal-sized pairs don't share them.
    expect_true(ecodive:::pair_seed(3L, 1L) != ecodive:::pair_seed(3L, 2L))
  }

  # Every p.adjust() method, matched as p.adjust() matches it.
  for (m in c(stats::p.adjust.methods, 'bonf', 'hoc', 'f')) {
    res  <- anosim(d, g, pairwise = TRUE, p_adjust = m, permutations = 49, cpus = 1)
    full <- match.arg(m, stats::p.adjust.methods)
    expect_identical(res$p_adjust, full, info = m)
    expect_identical(res$pairwise$p_adjusted, stats::p.adjust(res$pairwise$p_value, full), info = m)
  }
  res <- mrpp(d, g, pairwise = TRUE, p_adjust = 'none', permutations = 49, cpus = 1)
  expect_identical(res$pairwise$p_adjusted, res$pairwise$p_value)

  # The pairs of one factor in a formula, with strata cut down to each pair.
  r1 <- permanova(d, data = meta, formula = ~ grp, strata = 'block', pairwise = TRUE, permutations = 49, cpus = 1)
  r2 <- permanova(d, 'grp', meta, strata = 'block', pairwise = TRUE, permutations = 49, cpus = 1)
  expect_identical(r1$pairwise, r2$pairwise)
  idx <- which(g %in% c('a', 'b'))
  one <- permanova(
    pair_dm(d, idx), g[idx], strata = meta$block[idx], permutations = 49,
    seed = ecodive:::pair_seed(0L, 1L), cpus = 1 )
  expect_identical(r1$pairwise$p_value[[1]], one$p_value)

  # Within strata, the induced permutations stay within strata.
  pb  <- ex$perms$blocks
  b   <- meta$block
  for (k in seq_along(labs)) {
    idx <- which(g %in% strsplit(labs[[k]], '-')[[1]])
    ip  <- ecodive:::induced_perms(pb, idx)
    expect_true(all(apply(ip, 1, function (r) all(sort(r) == seq_along(idx)) && all(b[idx][r] == b[idx]))))
  }

  # A factor's level order sets the pairs' order; its unused levels are
  # dropped.
  gf  <- factor(g, c('d', 'c', 'b', 'a', 'z'))
  res <- mrpp(d, gf, pairwise = TRUE, permutations = 0, cpus = 1)
  expect_identical(rownames(res$pairwise), c('d-c', 'd-b', 'd-a', 'c-b', 'c-a', 'b-a'))
  expect_true(all(is.na(res$pairwise$p_value)))

  # The induced permutations of the identity are the identity; a cycle
  # through samples outside the pair is skipped.
  expect_identical(ecodive:::induced_perms(matrix(1:5, 1), c(2L, 4L)), matrix(1:2, 1))
  expect_identical(ecodive:::induced_perms(matrix(c(2L, 3L, 4L, 5L, 1L), 1), c(1L, 4L)), matrix(2:1, 1))
  expect_identical(ecodive:::induced_perms(matrix(integer(0), 0, 5), 1:2), matrix(integer(0), 0, 2))

  # dist_subset() matches as.matrix().
  dv <- ecodive:::dist_vector(d)
  sb <- ecodive:::dist_subset(dv, c(2L, 5L, 9L, 34L))
  expect_identical(c(sb), c(pair_dm(d, c(2, 5, 9, 34))))
  expect_identical(attr(sb, 'Labels'), rownames(meta)[c(2, 5, 9, 34)])



  # Result format. ====

  res <- permanova(d, g, pairwise = TRUE, permutations = 19, cpus = 1)
  pw  <- res$pairwise
  expect_inherits(pw, 'data.frame')
  expect_identical(names(pw), c('F', 'R2', 'z', 'n', 'n_perm', 'p_value', 'p_adjusted'))
  expect_identical(rownames(pw), labs)
  expect_identical(pw$n, c(20L, 19L, 17L, 17L, 15L, 14L))
  expect_identical(pw$n_perm, rep(19L, 6))
  expect_identical(res$p_adjust, 'holm')
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm', 'r2', 'table', 'by', 'pairwise',
    'p_adjust', 'permutation' ))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("Pairwise comparisons (p-values adjusted by Holm's method):", out, fixed = TRUE)))
  expect_true(any(grepl('^c-d ', out)))

  res <- anosim(d, g, pairwise = TRUE, permutations = 19, cpus = 1)
  expect_identical(names(res$pairwise), c('R', 'z', 'n', 'n_perm', 'p_value', 'p_adjusted'))
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm', 'ranks', 'pairwise', 'p_adjust', 'permutation' ))

  res <- mrpp(d, g, pairwise = TRUE, p_adjust = 'BH', permutations = 19, cpus = 1)
  expect_identical(names(res$pairwise), c('delta', 'A', 'z', 'n', 'n_perm', 'p_value', 'p_adjusted'))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl('adjusted by Benjamini-Hochberg', out, fixed = TRUE)))
  out <- utils::capture.output(print(mrpp(d, g, pairwise = TRUE, p_adjust = 'none', permutations = 9, cpus = 1)))
  expect_true(any(grepl('(p-values not adjusted)', out, fixed = TRUE)))

  # Without pairwise: NULL fields, and nothing printed.
  for (f in tests) {
    res <- f(d, g, permutations = 9, cpus = 1)
    expect_null(res$pairwise)
    expect_null(res$p_adjust)
    expect_true('pairwise' %in% names(res))
    expect_false(any(grepl('Pairwise', utils::capture.output(print(res)))))
  }

  # Through dist_test().
  r1 <- dist_test(d, 'adonis2', g, pairwise = TRUE, p_adjust = 'BY', permutations = 19, cpus = 1)
  r2 <- permanova(d, g, pairwise = TRUE, p_adjust = 'BY', permutations = 19, cpus = 1)
  expect_identical(r1, r2)
  expect_identical(
    dist_test(d, 'anosim', 'grp', meta, pairwise = TRUE, permutations = 19, cpus = 1),
    anosim(d, 'grp', meta, pairwise = TRUE, permutations = 19, cpus = 1) )



  # Pairs that can't be tested: NA, with a warning. ====

  gs <- g
  gs[1:2] <- c('s1', 's2')  # two groups of one sample each

  expect_warning(
    res <- anosim(d, gs, pairwise = TRUE, permutations = 19, cpus = 1),
    'Pair s1-s2 could not be tested: Every group has only one sample' )
  pw <- res$pairwise
  expect_true(all(is.na(unlist(pw['s1-s2', c('R', 'z', 'p_value', 'p_adjusted')]))))
  expect_identical(pw['s1-s2', 'n'], 2L)
  expect_identical(pw['s1-s2', 'n_perm'], 19L)
  # A group of one sample against one with replicates is tested.
  expect_true(all(!is.na(pw[rownames(pw) != 's1-s2', 'p_value'])))
  expect_identical(pw$p_adjusted, stats::p.adjust(pw$p_value, 'holm'))

  expect_warning(
    res <- permanova(d, gs, pairwise = TRUE, permutations = 19, cpus = 1),
    'Pair s1-s2: The model has no residual degrees of freedom' )
  expect_true(is.na(res$pairwise['s1-s2', 'p_value']))
  expect_true(all(!is.na(res$pairwise[rownames(res$pairwise) != 's1-s2', 'p_value'])))

  expect_warning(
    res <- mrpp(d, gs, pairwise = TRUE, permutations = 19, cpus = 1),
    'Pair s1-s2: No group has two or more samples' )
  expect_true(is.na(res$pairwise['s1-s2', 'delta']))

  # A matrix of permutations: the same pair, cut down.
  expect_warning(
    res <- anosim(d, gs, pairwise = TRUE, permutations = pm, cpus = 1),
    'Pair s1-s2 could not be tested' )
  expect_identical(res$pairwise['s1-s2', 'n_perm'], nrow(pm))

  # Two groups: one pair, adjusted to itself.
  res <- permanova(d, g == 'd', pairwise = TRUE, permutations = 19, cpus = 1)
  expect_identical(rownames(res$pairwise), 'FALSE-TRUE')
  expect_identical(res$pairwise$p_adjusted, res$pairwise$p_value)



  # Deterministic: the same for any number of cpus, and from run to run. ====

  for (nm in names(tests)) {
    f  <- tests[[nm]]
    r1 <- f(d, g, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1)
    for (cpus in if (at_home()) c(2, 3, 4, 8) else 2) {
      r2 <- f(d, g, pairwise = TRUE, permutations = 199, seed = 7, cpus = cpus)
      expect_identical(r1, r2, info = paste(nm, cpus))
    }
    expect_identical(r1, f(d, g, pairwise = TRUE, permutations = 199, seed = 7, cpus = 1), info = nm)
    expect_false(identical(r1$pairwise, f(d, g, pairwise = TRUE, permutations = 199, seed = 8, cpus = 1)$pairwise))
  }



  # Invalid arguments. ====

  for (f in tests) {
    expect_error(f(d, g, pairwise = NA),  '`pairwise` must be TRUE or FALSE')
    expect_error(f(d, g, pairwise = 'yes'), '`pairwise` must be TRUE or FALSE')
    expect_error(f(d, g, p_adjust = 'BH'), '`p_adjust` is only used with `pairwise = TRUE`')
    expect_error(f(d, g, pairwise = FALSE, p_adjust = 'holm'), '`p_adjust` is only used')
    expect_error(f(d, g, pairwise = TRUE, p_adjust = 'h'), "`p_adjust` must be one of 'holm'")
    expect_error(f(d, g, pairwise = TRUE, p_adjust = 'tukey'), '`p_adjust` must be one of')
    expect_error(f(d, g, pairwise = TRUE, p_adjust = NA), '`p_adjust` must be one of')
    expect_error(f(d, g, pairwise = TRUE, p_adjust = c('holm', 'BH')), '`p_adjust` must be one of')
  }
  expect_error(permdisp(d, g, p_adjust = 'BH'), '`p_adjust` is only used')
  expect_error(permdisp(d, g, pairwise = TRUE, p_adjust = 'x'), '`p_adjust` must be one of')

  # permanova(): one categorical variable.
  meta$ph <- seq_len(nrow(meta))
  expect_error(permanova(d, meta$ph, pairwise = TRUE), 'one categorical variable')
  expect_error(permanova(d, data = meta, formula = ~ grp + block, pairwise = TRUE), 'one categorical variable')
  expect_error(permanova(d, data = meta, formula = ~ grp + Condition(block), pairwise = TRUE), 'no `Condition\\(\\)`')

#})
