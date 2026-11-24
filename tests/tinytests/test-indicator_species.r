#test_that("indicator_species", {

  # indicator_species() against indicspecies::multipatt() and strassoc(),
  # and labdsv::indval() (outputs only: they are GPL). Fixture:
  # data-raw/expected/indicator_species.r. The permutations there are
  # ecodive's own (seed `ref$seed`), so p-values must be identical, with
  # random and with supplied permutations.

  max_cpus <- if (at_home()) 4L else 2L

  ref  <- readRDS(file.path('expected', 'indicator_species.rds'))
  x    <- ref$x
  grp  <- ref$grp
  lev  <- sort(unique(grp))
  perm_mat <- function (n, strata = NULL)
    t(.Call(C_perm_indices, n, ref$n_perm, strata, ref$seed))

  # NA where either is NA; the largest absolute difference elsewhere.
  same_na  <- function (a, b) identical(is.na(a), is.na(b))
  max_diff <- function (a, b) { ok <- !is.na(a) & !is.na(b); max(abs(a[ok] - b[ok]), 0) }


  # Every multipatt() case. ====

  case_args <- function (nm) {
    f    <- strsplit(nm, '/', fixed = TRUE)[[1]]
    args <- list(counts = x, groups = grp, index = f[[1]])
    switch(
      f[[2]],
      pa          = { args$counts <- ref$xpa },
      rel         = { args$counts <- ref$xrel },
      duleg       = { args$max_order <- 1L },
      strata      = { args$strata <- ref$strata },
      negative    = { args$counts <- ref$x_vary;   args$negative <- TRUE },
      negative_pa = { args$counts <- ref$xpa_vary; args$negative <- TRUE },
      max2        = { args$max_order <- 2L },
      restcomb    = { args$combinations <- c(1, 3, 5, 6) },
      restcomb_all  = { args$combinations <- c('forest', 'marsh', 'forest+marsh+meadow') },
      restcomb_max2 = { args$combinations <- c(1, 2, 4, 7); args$max_order <- 2L },
      four        = { args$counts <- ref$x4; args$groups <- ref$grp4 },
      four_min2   = { args$counts <- ref$x4; args$groups <- ref$grp4
                      args$min_order <- 2L; args$max_order <- 3L },
      four_restcomb = { args$counts <- ref$x4; args$groups <- ref$grp4
                        args$combinations <- c(2, 5, 7, 11, 14) } )
    args
  }

  for (nm in names(ref$cases)) {

    cs   <- ref$cases[[nm]]
    args <- c(case_args(nm), list(p_adjust = 'none', seed = ref$seed, permutations = ref$n_perm, cpus = 1L))
    d    <- do.call(indicator_species_details, args)
    r    <- d$result
    ok   <- !cs$near_tie

    # Every combination's value.
    expect_identical(colnames(d$strength), cs$comb, info = nm)
    expect_true(same_na(d$strength, cs$str), info = nm)
    expect_true(max_diff(d$strength, cs$str) < 1e-12, info = nm)

    # The best combination, its value, A and B, and the p-value.
    expect_true(same_na(r$stat, cs$stat), info = nm)
    expect_true(max_diff(r$stat, cs$stat) < 1e-12, info = nm)
    expect_identical(r$combination[ok], cs$best[ok], info = nm)
    # With restcomb, multipatt() tests a best combination of all the groups
    # (and gives 1); ecodive never does (NA), as multipatt() without it.
    all_k <- r$n_groups %in% length(unique(args$groups))
    expect_true(all(is.na(r$p_value[all_k])), info = nm)
    if (grepl('restcomb', nm)) cs$p_value[all_k & cs$p_value %in% 1] <- NA
    expect_identical(r$p_value[ok], cs$p_value[ok], info = nm)
    if (!is.null(cs$a)) {
      expect_true(max_diff(r$a[ok], cs$a[ok]) < 1e-12 && same_na(r$a[ok], cs$a[ok]), info = nm)
      expect_true(max_diff(r$b[ok], cs$b[ok]) < 1e-12 && same_na(r$b[ok], cs$b[ok]), info = nm)
    }
    else {
      expect_true(all(is.na(r$a)) && all(is.na(r$b)), info = nm)
    }

    # The same permutations, supplied as a matrix.
    args$permutations <- perm_mat(nrow(args$counts), if (!is.null(args$strata)) ref$strata - 1L)
    args$strata <- NULL
    r2 <- do.call(indicator_species, args)
    expect_identical(r2$p_value, r$p_value, info = nm)
    expect_identical(attr(r2, 'permutation'), 'supplied matrix', info = nm)
  }

  # Near ties (broken at random by multipatt()) are rare in these data.
  near <- unlist(lapply(ref$cases, `[[`, 'near_tie'))
  expect_true(sum(near) <= 5)


  # strassoc(): single groups, every index, and A and B. ====

  ix <- c(IndVal = 'indval', IndVal.g = 'indval_g', r = 'r', r.g = 'r_g')
  for (f in names(ix)) {
    d <- indicator_species_details(x, grp, index = ix[[f]], max_order = 1L, permutations = 0, cpus = 1L)
    # strassoc() gives 0 where the index is undefined (F09 is absent
    # everywhere, F10 the same everywhere), and multipatt() NaN.
    sa <- ref$strassoc[[f]]
    expect_identical(colnames(d$strength), lev, info = f)
    expect_true(all(sa[is.na(d$strength)] == 0), info = f)
    expect_identical(unname(which(is.na(d$strength[, 1]))), if (f %in% c('r', 'r.g')) 9:10 else 9L, info = f)
    expect_true(max_diff(d$strength, sa) < 1e-12, info = f)
    if (f %in% c('IndVal', 'IndVal.g')) {
      r  <- d$result
      at <- cbind(seq_len(nrow(r)), match(r$combination, lev))
      Aj <- ref$strassoc[[if (f == 'IndVal') 'A' else 'A.g']]
      expect_true(max_diff(r$a, Aj[at]) < 1e-12 && same_na(r$a, Aj[at]), info = f)
      expect_true(max_diff(r$b, ref$strassoc$B[at]) < 1e-12, info = f)
    }
  }


  # labdsv::indval(): Dufrene & Legendre's IndVal (A B) is index =
  # 'indval_g', single groups, squared. ====

  lv <- ref$labdsv
  d  <- indicator_species_details(x[, lv$features], grp, max_order = 1L, permutations = 0, cpus = 1L)
  expect_true(max(abs(d$strength^2 - lv$indval)) < 1e-12)
  expect_true(max(abs(d$result$stat^2 - lv$indcls)) < 1e-12)
  expect_identical(d$result$combination, lev[lv$maxcls])


  # Ties are exact, as in multipatt(): the probe has one permutation, whose
  # value is below the observed one by about d / 32 (from 1e-15 to 1e-2),
  # so p = 1/2 for every d. ====

  pg <- rep(1:2, c(4, 4))
  ds <- 10^-(15:2)
  px <- sapply(ds, function (d) c(5, 4, 3, 2, 2 - d, 0, 0, 0))
  pp <- matrix(c(1, 2, 3, 5, 4, 6, 7, 8), 1)
  for (f in names(ix)) {
    r <- indicator_species(px, pg, index = ix[[f]], max_order = 1L, permutations = pp, p_adjust = 'none', cpus = 1L)
    expect_identical(r$p_value, ref$probe_ties[, f], info = f)
    expect_true(all(r$p_value == 0.5), info = f)
  }


  # Determinism: any cpus, the same seed, dense = sparse, margin. ====

  set.seed(11)
  y  <- matrix(rpois(60 * 80, 1) * rbinom(60 * 80, 1, 0.3), 60, 80,
               dimnames = list(paste0('S', 1:60), paste0('T', 1:80)))
  yg <- sample(rep(c('p', 'q', 'r', 's'), c(10, 15, 15, 20)))
  y[yg == 'p', 1:5] <- y[yg == 'p', 1:5] + 2

  for (index in ix) {
    r1 <- indicator_species(y, yg, index = index, permutations = 299, seed = 3, cpus = 1L)
    for (cpus in 2:max_cpus)
      expect_identical(indicator_species(y, yg, index = index, permutations = 299, seed = 3, cpus = cpus), r1, info = index)
    expect_identical(indicator_species(y, yg, index = index, permutations = 299, seed = 3, cpus = 1L), r1)
    expect_identical(indicator_species(t(y), yg, margin = 2L, index = index, permutations = 299, seed = 3, cpus = 1L), r1)
    expect_identical(indicator_species(Matrix::Matrix(y, sparse = TRUE), yg, index = index, permutations = 299, seed = 3, cpus = 1L), r1)
    expect_identical(indicator_species(slam::as.simple_triplet_matrix(y), yg, index = index, permutations = 299, seed = 3, cpus = 1L), r1)
  }
  r3 <- indicator_species(y, yg, permutations = 299, seed = 4, cpus = 1L)
  r0 <- indicator_species(y, yg, permutations = 299, seed = 3, cpus = 1L)
  expect_false(identical(r3$p_value, r0$p_value))
  expect_identical(r3$stat, r0$stat)

  # Logical tables are presence-absence.
  expect_identical(
    indicator_species(y > 0, yg, permutations = 9, cpus = 1L),
    indicator_species((y > 0) * 1, yg, permutations = 9, cpus = 1L) )

  # The groups' order: a factor's levels.
  fg <- factor(yg, levels = c('s', 'r', 'q', 'p'))
  rf <- indicator_species(y, fg, permutations = 0, cpus = 1L)
  expect_identical(attr(rf, 'groups'), c('s', 'r', 'q', 'p'))
  expect_equal(rf$stat, indicator_species(y, yg, permutations = 0, cpus = 1L)$stat, tolerance = 1e-12)
  expect_identical(rf$combination[[1]], 'p')
  expect_identical(rf$combination_id[[1]], 4L)


  # Groups from `data`, strata from `data`, and a data frame's row order. ====

  meta <- data.frame(g = yg, b = rep(1:3, 20), row.names = rownames(y))
  r4   <- indicator_species(y, 'g', meta, strata = 'b', permutations = 99, cpus = 1L)
  expect_identical(r4, indicator_species(y, yg, strata = rep(1:3, 20), permutations = 99, cpus = 1L))
  expect_identical(attr(r4, 'permutation'), 'within strata')
  expect_identical(r4, indicator_species(y, 'g', meta[60:1, ], strata = 'b', permutations = 99, cpus = 1L))


  # Invariance: with single groups, 'indval_g' doesn't change when every
  # sample of one group is duplicated (its mean and fidelity don't). ====

  dup <- c(seq_len(60), which(yg == 'q'))
  a1  <- indicator_species(y, yg, max_order = 1L, permutations = 0, cpus = 1L)
  a2  <- indicator_species(y[dup, ], yg[dup], max_order = 1L, permutations = 0, cpus = 1L)
  expect_equal(a2$stat, a1$stat, tolerance = 1e-12)
  expect_identical(a2$combination, a1$combination)


  # Combinations: names in any order, ids, the order filter, ids of four
  # groups, and group names with '+'. ====

  cb <- c('q+p', 'r', 's+r+q')
  r5 <- indicator_species(y, yg, combinations = cb, permutations = 0, cpus = 1L)
  expect_identical(attr(r5, 'combinations'), c(r = 3L, 'p+q' = 5L, 'q+r+s' = 14L))
  expect_identical(r5, indicator_species(y, yg, combinations = c(14, 5, 3, 5), permutations = 0, cpus = 1L))
  expect_identical(
    attr(indicator_species(y, yg, combinations = cb, max_order = 2, permutations = 0, cpus = 1L), 'combinations'),
    c(r = 3L, 'p+q' = 5L) )
  all4 <- attr(indicator_species(y, yg, permutations = 0, cpus = 1L), 'combinations')
  expect_identical(unname(all4), 1:15)
  expect_identical(names(all4)[c(1, 5, 10, 11, 15)], c('p', 'p+q', 'r+s', 'p+q+r', 'p+q+r+s'))
  for (id in 1:15)
    expect_equal(ecodive:::indsp_rank(ecodive:::indsp_unrank(id, 4L), 4L), id)
  expect_identical(ecodive:::indsp_unrank(1023, 10L), 1:10)
  expect_equal(ecodive:::indsp_rank(c(2L, 9L, 10L), 10L),
               55 + which(apply(utils::combn(10, 3), 2, function (s) all(s == c(2, 9, 10)))))
  pg2 <- rep(c('a+b', 'a', 'b'), each = 4)
  r6  <- indicator_species(diag(12)[, 1:3] + 1, pg2, combinations = c('a+b+a', 'a+b'), permutations = 0, cpus = 1L)
  expect_identical(attr(r6, 'combinations'), c('a+b' = 2L, 'a+a+b' = 4L))
  r7  <- indicator_species(y, yg, index = 'r', permutations = 0, cpus = 1L)
  expect_identical(length(attr(r7, 'combinations')), 14L)


  # Edge cases. ====

  ez <- cbind(
    only_a = c(3, 2, 4, 0, 0, 0, 0, 0, 0),
    zero   = 0,
    const  = 2,
    every  = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
    single = c(0, 0, 0, 0, 0, 0, 0, 0, 5) )
  eg <- rep(c('a', 'b', 'c'), each = 3)
  e1 <- indicator_species(ez, eg, permutations = 99, p_adjust = 'none', cpus = 1L)
  expect_identical(e1$combination, c('a', NA, 'a+b+c', 'a+b+c', 'c'))
  expect_identical(is.na(e1$p_value), c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(e1$stat[1:3], c(1, NA, 1))
  expect_identical(e1$reject, !is.na(e1$p_value) & e1$p_value <= 0.05)
  expect_identical(e1$p_adjusted, e1$p_value)
  e2 <- indicator_species(ez, eg, index = 'r', negative = TRUE, permutations = 99, cpus = 1L)
  expect_identical(e2$combination[2:3], c(NA_character_, NA_character_))
  expect_true(all(is.na(e2$a)))
  expect_identical(e2$combination[[5]], 'c')

  # Holm over the tested features only.
  e3 <- indicator_species(ez, eg, permutations = 99, cpus = 1L)
  expect_equal(e3$p_adjusted, stats::p.adjust(e1$p_value, 'holm'))

  # A group of one sample, two groups, no permutations.
  g1 <- c(rep('a', 8), 'b')
  expect_identical(indicator_species(ez, g1, permutations = 9, cpus = 1L)$combination[[5]], 'b')
  e4 <- indicator_species(ez, rep(c('a', 'b'), c(4, 5)), permutations = 0, cpus = 1L)
  expect_true(all(is.na(e4$p_value)) && !any(e4$reject))
  expect_identical(attr(e4, 'n_perm'), 0L)

  # Many groups: only with fewer combinations.
  mg <- as.character(1:21)
  mx <- matrix(rpois(21 * 3, 2), 21, 3)
  expect_error(indicator_species(mx, mg, permutations = 0, cpus = 1L), 'more than 2\\^20')
  expect_identical(length(attr(indicator_species(mx, mg, max_order = 2, permutations = 0, cpus = 1L), 'combinations')), 231L)
  expect_error(indicator_species(matrix(1, 31, 2), as.character(1:31)), 'at most 30')


  # Errors. ====

  expect_error(indicator_species(ez, eg, combinations = 'a+z'), '"a\\+z" is not a group')
  expect_error(indicator_species(ez, eg, combinations = 8), 'from 1 to 7')
  expect_error(indicator_species(ez, eg, combinations = list(1)), '`combinations` must be')
  expect_error(indicator_species(ez, eg, combinations = 'a+b+c', index = 'r'), 'None of `combinations`')
  expect_error(indicator_species(ez, eg, negative = TRUE), 'only for')
  expect_error(indicator_species(ez, eg, negative = NA), '`negative`')
  expect_error(indicator_species(ez, eg, index = 'cos'), '`index` must be')
  expect_identical(attr(indicator_species(ez, eg, index = 'IndVal.g', permutations = 0), 'index'), 'indval_g')
  expect_identical(attr(indicator_species(ez, eg, index = 'R.G', permutations = 0), 'index'), 'r_g')
  expect_error(indicator_species(ez, eg, min_order = 3, max_order = 2), '`min_order` must not')
  expect_error(indicator_species(ez, eg, max_order = 4), 'at most the number of groups')
  expect_error(indicator_species(ez, eg, min_order = 0), '`min_order` must be')
  expect_error(indicator_species(ez, eg, index = 'r', min_order = 3), "can't be tested")
  expect_error(indicator_species(ez, eg, alpha = 1), '`alpha`')
  expect_error(indicator_species(ez, eg, p_adjust = 'nope'), '`p_adjust`')
  expect_error(indicator_species(ez, rep('a', 9)), 'at least two')
  expect_error(indicator_species(ez, eg, strata = rep(1:3, 3), permutations = matrix(1:9, 1)), 'not both')
  bad <- ez; bad[1, 1] <- -1
  expect_error(indicator_species(bad, eg), 'negative')
  bad[1, 1] <- NA
  expect_error(indicator_species(bad, eg))


  # da_test(): the registry and its aliases. ====

  expect_identical(
    da_test(y, 'indicator_species', yg, permutations = 49, seed = 2, cpus = 1L),
    indicator_species(y, yg, permutations = 49, seed = 2, cpus = 1L) )
  expect_identical(
    da_test(y, 'indval', 'g', meta, index = 'r_g', permutations = 49, cpus = 1L),
    indicator_species(y, 'g', meta, index = 'r_g', permutations = 49, cpus = 1L) )
  expect_identical(match_da_test('multipatt')$id, 'indicator_species')
  expect_identical(match_da_test('Indicator Species Analysis')$id, 'indicator_species')
  expect_error(da_sensitivity(y, 'indicator_species', yg, zeros = TRUE), 'zero')


  # Attributes, printing, the plot. ====

  expect_inherits(r1, c('ecodive_indicator_species', 'data.frame'))
  expect_identical(names(r1), c('feature', 'combination', 'combination_id', 'n_groups', 'stat',
                                'a', 'b', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(
    names(attributes(r1)),
    c('names', 'row.names', 'class', 'method', 'index', 'groups', 'combinations',
      'alpha', 'p_adjust', 'n_perm', 'permutation', 'means') )
  expect_identical(attr(r1, 'method'), 'indicator_species')
  expect_identical(attr(r1, 'p_adjust'), 'holm')
  expect_identical(attr(r1, 'permutation'), 'free')
  expect_identical(dim(attr(r1, 'means')), c(4L, 80L))
  expect_equal(attr(r1, 'means')['p', ], colMeans(y[yg == 'p', ]))
  expect_true(is.integer(r1$combination_id) && is.integer(r1$n_groups))

  out <- utils::capture.output(print(r1))
  expect_true(any(grepl('Indicator species (r.g)', out, fixed = TRUE)))
  out <- utils::capture.output(print(e1))
  expect_true(any(grepl('Indicator species (IndVal.g)', out, fixed = TRUE)))

  grDevices::pdf(NULL)
  rp <- indicator_species(y, yg, p_adjust = 'none', permutations = 99, cpus = 1L)
  expect_true(any(rp$reject))
  v <- withVisible(plot(rp))
  expect_identical(v$value, rp)
  expect_false(v$visible)
  expect_identical(plot(rp, 'heatmap', n_features = 3, main = 'x', sub = 'y'), rp)
  expect_silent(plot(e1))
  expect_silent(plot(e1, 'heatmap'))
  expect_silent(plot(e2))
  expect_error(plot(rp, n_features = -1), '`n_features`')
  expect_error(plot(rp, 'pie'))
  grDevices::dev.off()

#})
