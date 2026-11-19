#test_that("capda", {

  # Inputs and reference results, written by data-raw/expected/capda.r.
  ref <- readRDS(file.path('expected', 'capda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  # `x` (without attributes) with each column's sign matched to `y`'s.
  match_signs <- function (x, y) {
    x <- unclass(x)
    attributes(x) <- list(dim = dim(x))
    sweep(x, 2, sign(colSums(x * y)), '*')
  }

  attrs <- c(
    'dim', 'dimnames', 'method', 'eig', 'prop', 'n_pcoa_axes', 'success',
    'loo_table', 'loo_class', 'trace_stat', 'p_value', 'n_perm', 'trace_perm',
    'scaling', 'groups', 'class' )



  # Same as BiodiversityR::CAPdiscrim(): m, the classification success for
  # every m, the leave-one-out classes, delta^2, Pillai's trace and the
  # scores. Stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      lbl <- cs$name
      ex  <- cs$expected
      n   <- length(cs$groups)
      g   <- nlevels(cs$groups)

      # mmax = 10, as CAPdiscrim()'s default.
      ord <- capda(cs$dm, cs$groups, max_pcoa_axes = 10, permutations = 0, cpus = 1)
      m   <- attr(ord, 'n_pcoa_axes')
      suc <- attr(ord, 'success')
      k   <- min(g - 1L, m)

      expect_identical(m, as.integer(ex$m), info = lbl)
      expect_identical(suc$n_pcoa_axes, 1:10, info = lbl)
      expect_equal(suc$overall * 100, ex$by_m[1:10], tolerance = 1e-12, info = lbl)
      expect_equal(unlist(suc[m, levels(cs$groups)]) * 100, ex$percent_level, tolerance = 1e-12, info = lbl)
      expect_identical(as.character(attr(ord, 'loo_class')), ex$cv, info = lbl)
      expect_identical(
        unclass(attr(ord, 'loo_table')),
        unclass(table(group = cs$groups, loo_class = factor(ex$cv, levels(cs$groups)))), info = lbl )
      expect_identical(dim(ord), c(n, min(k, 3L)), info = lbl)
      expect_equal(attr(ord, 'eig'), ex$delta2[seq_len(ncol(ord))], tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'trace_stat'), ex$pillai, tolerance = 1e-10, info = lbl)
      all_ax <- capda(cs$dm, cs$groups, n_pcoa_axes = m, n_axes = 99, permutations = 0)
      expect_equal(sum(attr(all_ax, 'eig')), ex$pillai, tolerance = 1e-10, info = lbl)
      x <- ex$x[, seq_len(ncol(ord)), drop = FALSE]
      expect_equal(match_signs(ord, x), x, tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'), attr(ord, 'eig') / ex$pillai, tolerance = 1e-10, info = lbl)

      # The default search, up to 50 axes: every positive PCoA axis here.
      # Leave-one-out refits agree for every m; CAPdiscrim() agrees but where
      # MASS's posteriors underflow (mite, sample 67, m >= 20).
      ord <- capda(cs$dm, cs$groups, permutations = 0, cpus = 1)
      suc <- attr(ord, 'success')
      ok  <- if (lbl == 'mite_shrub') 1:19 else seq_along(ex$by_m)
      expect_identical(suc$n_pcoa_axes, seq_along(ex$by_m), info = lbl)
      expect_equal(suc$overall * 100, ex$by_m_exact, tolerance = 1e-12, info = lbl)
      expect_equal(suc$overall[ok] * 100, ex$by_m[ok], tolerance = 1e-12, info = lbl)
      expect_identical(attr(ord, 'n_pcoa_axes'), which.max(ex$by_m_exact), info = lbl)
    }


    # A singleton group: as MASS::lda(CV = TRUE), whose class for it is NA,
    # counted as misclassified.
    s   <- ref$single
    ord <- capda(s$dm, s$groups, permutations = 0, cpus = 1)
    suc <- attr(ord, 'success')
    hit <- !is.na(s$cv) & s$cv == as.character(s$groups)
    expect_equal(suc$overall, colMeans(hit), tolerance = 1e-15)
    expect_identical(suc$a, rep(0, 14))
    m <- attr(ord, 'n_pcoa_axes')
    expect_identical(m, which.max(colSums(hit)))
    expect_identical(as.character(attr(ord, 'loo_class')), s$cv[, m])
    expect_true(is.na(attr(ord, 'loo_class')[[1]]))
    expect_identical(colnames(attr(ord, 'loo_table')), c('a', 'b', 'c', NA))
    expect_identical(sum(attr(ord, 'loo_table')), 20L)


    # The permutation test, with the references' permutations: Pillai's trace
    # (MANOVA) of each permutation, and the p-value (also adonis2()'s).
    p   <- ref$perm
    ord <- capda(p$dm, p$groups, n_pcoa_axes = p$m, permutations = p$perms, cpus = 1)
    expect_equal(attr(ord, 'trace_stat'), p$trace_obs, tolerance = 1e-10)
    expect_equal(attr(ord, 'trace_perm'), p$trace_perm, tolerance = 1e-10)
    expect_identical(attr(ord, 'p_value'), p$p_value)
    expect_identical(attr(ord, 'p_value'), p$adonis_p)
    expect_identical(attr(ord, 'n_perm'), 199L)

    # Within strata: 9999 random permutations each, so the p-values agree
    # statistically.
    ord <- capda(p$dm, p$groups, n_pcoa_axes = p$m, strata = p$strata, permutations = 9999, seed = 3)
    se  <- sqrt(2 * p$strata_p * (1 - p$strata_p) / 9999)
    expect_true(abs(attr(ord, 'p_value') - p$strata_p) < 4 * se)
  }



  # The same scores as lda() on pcoa()'s axes. ====

  set.seed(1)
  n   <- 36
  cts <- matrix(rpois(n * 15, 3), n, 15, dimnames = list(sprintf('S%02d', 1:n), NULL))
  grp <- factor(rep(c('x', 'y', 'z'), each = 12), c('z', 'x', 'y'))
  cts[grp == 'x', 1:3] <- cts[grp == 'x', 1:3] + 3L
  cts[grp == 'y', 4:5] <- cts[grp == 'y', 4:5] + 2L
  dm  <- bray(cts)
  n_pos <- sum(suppressWarnings(attr(pcoa(dm, n_axes = n - 1L), 'eig')) > 0)  # positive PCoA axes
  meta <- data.frame(row.names = rev(rownames(cts)), trt = rev(grp), blk = rep(1:4, 9))

  for (M in c(1L, 2L, 5L, 12L)) {
    o  <- capda(dm, grp, n_pcoa_axes = M, permutations = 0)
    l  <- lda(pcoa(dm, n_axes = M), grp)
    expect_equal(unname(unclass(o)[, ]), unname(unclass(l)[, ]), tolerance = 1e-12, info = M)
    expect_equal(unname(attr(o, 'scaling')), unname(attr(l, 'scaling')), tolerance = 1e-12, info = M)
    f <- attr(l, 'eig') * 2
    expect_equal(attr(o, 'eig'), f / (f + n - 3), tolerance = 1e-12, info = M)
    expect_identical(rownames(attr(o, 'scaling')), paste0('PCoA', seq_len(M)))
    expect_identical(nrow(attr(o, 'success')), 1L)
    expect_identical(attr(o, 'success')$n_pcoa_axes, M)
  }



  # The result. ====

  ord <- capda(dm, 'trt', meta, permutations = 99)
  m   <- attr(ord, 'n_pcoa_axes')
  suc <- attr(ord, 'success')
  expect_identical(names(attributes(ord)), attrs)
  expect_identical(attr(ord, 'method'), 'capda')
  expect_identical(colnames(ord), c('CAP1', 'CAP2'))
  expect_identical(rownames(ord), rownames(cts))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_identical(names(suc), c('n_pcoa_axes', 'overall', 'z', 'x', 'y'))
  expect_identical(suc$n_pcoa_axes, seq_len(min(50L, n - 3L - 1L, n_pos)))
  expect_identical(m, which.max(suc$overall))
  expect_equal(suc$overall[[m]], mean(attr(ord, 'loo_class') == attr(ord, 'groups')))
  expect_identical(names(attr(ord, 'loo_class')), rownames(cts))
  expect_identical(levels(attr(ord, 'loo_class')), c('z', 'x', 'y'))
  expect_identical(attr(ord, 'groups'), setNames(grp, rownames(cts)))
  expect_identical(names(dimnames(attr(ord, 'loo_table'))), c('group', 'loo_class'))
  expect_true(all(attr(ord, 'eig') > 0 & attr(ord, 'eig') < 1))
  expect_equal(sum(attr(ord, 'prop')), 1)
  expect_equal(attr(ord, 'trace_stat'), sum(attr(ord, 'eig')), tolerance = 1e-10)
  expect_identical(attr(ord, 'n_perm'), 99L)
  expect_identical(length(attr(ord, 'trace_perm')), 99L)
  expect_identical(attr(ord, 'p_value'), 0.01)
  expect_identical(dim(attr(ord, 'scaling')), c(m, 2L))
  expect_identical(unname(apply(unclass(ord)[, ], 2, function (v) sign(v[which.max(abs(v))]))), c(1, 1))

  # The same from a symmetric matrix, and through dist_ord().
  expect_identical(capda(as.matrix(dm), 'trt', meta, permutations = 99), ord)
  expect_identical(dist_ord(dm, 'capda', 'trt', meta, permutations = 99), ord)
  expect_identical(dist_ord(dm, 'capda', grp, n_pcoa_axes = 3, permutations = 0),
                   capda(dm, grp, n_pcoa_axes = 3, permutations = 0))

  # max_pcoa_axes; n_axes.
  o <- capda(dm, grp, max_pcoa_axes = 4, permutations = 0)
  expect_identical(attr(o, 'success')$n_pcoa_axes, 1:4)
  expect_identical(attr(o, 'success'), attr(ord, 'success')[1:4, ], info = 'success is nested')
  expect_identical(nrow(attr(capda(dm, grp, max_pcoa_axes = 1e6, permutations = 0), 'success')), as.integer(min(n - 4L, n_pos)))
  expect_identical(ncol(capda(dm, grp, n_axes = 1, permutations = 0)), 1L)
  expect_identical(ncol(capda(dm, grp, n_axes = 9, permutations = 0)), 2L)

  # Without permutations.
  o <- capda(dm, grp, permutations = 0)
  expect_identical(attr(o, 'p_value'), NA_real_)
  expect_identical(attr(o, 'n_perm'), 0L)
  expect_identical(attr(o, 'trace_perm'), numeric(0))
  expect_identical(unclass(o)[, ], unclass(ord)[, ])

  # Two groups: one axis.
  g2 <- rep(c('a', 'b'), each = n / 2)
  o  <- capda(dm, g2, permutations = 9)
  expect_identical(colnames(o), 'CAP1')
  expect_identical(length(attr(o, 'eig')), 1L)
  expect_identical(attr(o, 'prop'), 1)

  # Groups with the same means: no axes, with a warning.
  x <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(2, 0), c(-2, 0), c(0, 2), c(0, -2))
  expect_warning(o <- capda(dist(x), rep(c('a', 'b'), each = 4), permutations = 9), 'same means')
  expect_identical(dim(o), c(8L, 0L))
  expect_true(attr(o, 'trace_stat') < 1e-12)



  # Permutations: strata, seeds, and any number of threads. ====

  # Strata that are the groups: every permutation keeps the groups.
  o <- capda(dm, grp, n_pcoa_axes = 3, strata = grp, permutations = 19)
  expect_identical(attr(o, 'p_value'), 1)
  expect_equal(attr(o, 'trace_perm'), rep(attr(o, 'trace_stat'), 19), tolerance = 1e-12)
  expect_identical(capda(dm, 'trt', meta, n_pcoa_axes = 3, strata = 'blk', permutations = 19),
                   capda(dm, grp, n_pcoa_axes = 3, strata = rev(meta$blk), permutations = 19))

  a <- capda(dm, grp, n_pcoa_axes = 3, permutations = 99, seed = 1)
  b <- capda(dm, grp, n_pcoa_axes = 3, permutations = 99, seed = 2)
  expect_identical(capda(dm, grp, n_pcoa_axes = 3, permutations = 99, seed = 1), a)
  expect_false(identical(attr(a, 'trace_perm'), attr(b, 'trace_perm')))
  expect_identical(unclass(a)[, ], unclass(b)[, ])

  # A supplied matrix: row 1 the identity gives the observed trace.
  pm <- rbind(seq_len(n), t(replicate(9, sample(n))))
  o  <- capda(dm, grp, n_pcoa_axes = 3, permutations = pm)
  expect_equal(attr(o, 'trace_perm')[[1]], attr(o, 'trace_stat'), tolerance = 1e-14)
  expect_identical(attr(o, 'n_perm'), 10L)

  big_n <- 300
  set.seed(2)
  bc  <- matrix(rpois(big_n * 30, 2), big_n, 30)
  bg  <- sample(letters[1:5], big_n, TRUE)
  bc[bg == 'a', 1:4] <- bc[bg == 'a', 1:4] + 1L
  bdm <- bray(bc)
  r1  <- capda(bdm, bg, permutations = 999, cpus = 1)
  for (cpus in 2:max_cpus)
    expect_identical(capda(bdm, bg, permutations = 999, cpus = cpus), r1, info = cpus)
  expect_identical(nrow(attr(r1, 'success')), 50L)



  # Invalid input. ====

  expect_error(capda(dm),                         '`groups` is required')
  expect_error(capda(dm, rep('a', n)),            'at least two different values')
  expect_error(capda(dm, c(NA, grp[-1])),         'must not contain NAs')
  expect_error(capda(dm, grp[-1]),                'groups')
  expect_error(capda(dm, grp, n_pcoa_axes = 0),   '`n_pcoa_axes` must be a whole number')
  expect_error(capda(dm, grp, n_pcoa_axes = 1.5), '`n_pcoa_axes` must be a whole number')
  expect_error(capda(dm, grp, n_pcoa_axes = 34),  'at most the number of samples minus the number of groups: 36 - 3 = 33')
  expect_error(capda(dm, grp, max_pcoa_axes = 'a'), '`max_pcoa_axes` must be a whole number')
  expect_error(capda(dm, grp, n_axes = 0),        '`n_axes` must be a whole number')
  expect_error(capda(dm, grp, permutations = -1), '`permutations` must be')
  expect_error(capda(dm, grp, seed = 1.5))
  expect_error(capda(dm, grp, cpus = 0))
  expect_error(capda(dm, grp, strata = meta$blk, permutations = pm), 'not both')
  expect_error(capda(bray(cts[1:4, ]), c('a', 'b', 'c', 'c')), 'at least two more samples than groups')
  expect_true(is.matrix(capda(bray(cts[1:4, ]), c('a', 'b', 'c', 'c'), n_pcoa_axes = 1, permutations = 0)))

  # Fewer positive PCoA axes than asked for.
  x <- cbind(1:10, 0)
  expect_error(capda(dist(x), rep(1:2, 5), n_pcoa_axes = 2), 'Only 1 PCoA axes have positive eigenvalues')

  # A singular fit: the groups don't vary along the first PCoA axis.
  x <- rbind(c(0, 0), c(1, 0), c(2, 0), c(0, 5), c(1, 5), c(2, 5))
  expect_error(capda(dist(x), rep(1:2, each = 3)), 'first 1 PCoA axes is singular')
  expect_error(capda(dist(x), rep(1:2, each = 3), n_pcoa_axes = 2), 'first 2 PCoA axes is singular.*smaller')

#})
