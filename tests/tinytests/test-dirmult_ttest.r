#test_that("dirmult_ttest", {

  # dirmult_ttest(): scikit-bio 0.7.4's dirmult_ttest() on the same Gamma
  # draws, each draw's Welch statistics, SciPy's t distribution, and
  # scikit-bio's run with its own draws (data-raw/expected/dirmult_ttest.py);
  # the Gamma, normal and uniform samplers by distribution; threading,
  # input types, contrasts and argument checks.

  ref <- readRDS(file.path('expected', 'dirmult_ttest.rds'))

  given <- function (r, contrast = r$contrast, ...)
    suppressWarnings(dirmult_ttest_details(
      ref$tables[[r$table]], r$groups, contrast = contrast,
      n_draws = r$n_draws, gamma_draws = r$gamma, cpus = 1L, ... ))

  maxrel <- function (a, b, floor = 1e-300) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.nan(a), is.nan(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }

  cols <- c('t', 'log2_fc', 'ci_lower', 'ci_upper', 'p_value', 'p_adjusted')

  rd <- function (n, dist, param = 0, seed = 0L, stream = 0)
    .Call(C_rng_dist, as.integer(n), as.integer(dist), as.double(param),
          as.integer(seed), as.double(stream))


  # -----------------------------------------------------------------------
  # scikit-bio, on the same draws
  # -----------------------------------------------------------------------

  for (cs in names(ref$cases)) {
    r <- ref$cases[[cs]]
    for (pa in names(r$results)) {
      res <- given(r, p_adjust = pa)$result
      e   <- r$results[[pa]]
      for (k in cols)
        expect_true(maxrel(res[[k]], e[[k]]) < 1e-12, info = paste(cs, pa, k))
      expect_identical(res$reject, e$reject, info = paste(cs, pa))
    }

    # Each draw's statistics. The bounds are differences, so are compared
    # on the scale of their parts.
    d <- given(r)$draws
    for (k in c('diff', 'se', 'dof', 't', 'p_tail'))
      expect_true(maxrel(d[[k]], r$draws[[k]]) < 1e-12, info = paste(cs, k))
    for (k in c('lower', 'upper'))
      expect_true(
        maxrel(d[[k]], r$draws[[k]], floor = max(abs(r$draws$diff), na.rm = TRUE)) < 1e-12,
        info = paste(cs, k) )
  }

  # The NaN pattern of a side with one sample, as scikit-bio's.
  r <- ref$cases$single
  expect_warning(
    res <- dirmult_ttest_details(
      ref$tables$single, r$groups, contrast = r$contrast, n_draws = r$n_draws,
      gamma_draws = r$gamma, cpus = 1L )$result,
    'single sample' )
  expect_true(all(is.finite(res$log2_fc)))
  for (k in setdiff(cols, 'log2_fc')) expect_true(all(is.nan(res[[k]])), info = k)
  expect_false(any(res$reject))

  # The default contrast is the first group against the rest.
  r <- ref$cases$three_dflt
  expect_identical(given(r)$result, given(r, contrast = c('ctl', NA))$result)
  expect_identical(attr(given(r)$result, 'contrast'), c(treatment = 'ctl', reference = NA))

  # conf_level and alpha: each draw's interval at another level, and the
  # envelope over the draws.
  r  <- ref$cases$wide
  x  <- given(r, conf_level = 0.9, alpha = 0.2)
  q  <- qt(0.95, r$draws$dof) * r$draws$se
  expect_true(maxrel(x$draws$lower, r$draws$diff - q, floor = 1) < 1e-12)
  expect_true(maxrel(x$draws$upper, r$draws$diff + q, floor = 1) < 1e-12)
  expect_true(maxrel(x$result$ci_lower, apply(r$draws$diff - q, 2, min) / log(2)) < 1e-12)
  expect_true(maxrel(x$result$ci_upper, apply(r$draws$diff + q, 2, max) / log(2)) < 1e-12)
  out <- (x$result$ci_lower > 0) | (x$result$ci_upper < 0)
  expect_identical(x$result$reject, x$result$p_adjusted <= 0.2 & out)
  expect_identical(x$result[, c('t', 'log2_fc', 'p_value')], given(r)$result[, c('t', 'log2_fc', 'p_value')])


  # -----------------------------------------------------------------------
  # SciPy's t distribution: R's pt() and qt()
  # -----------------------------------------------------------------------

  td <- ref$tdist
  a  <- td$what == 'stdtr'
  expect_true(maxrel(pt(td$x[a], td$df[a]), td$value[a]) < 1e-12)
  expect_true(maxrel(qt(td$x[!a], td$df[!a]), td$value[!a]) < 1e-12)


  # -----------------------------------------------------------------------
  # scikit-bio with its own draws: Monte Carlo error only
  # -----------------------------------------------------------------------

  s   <- ref$sampled
  res <- dirmult_ttest(ref$tables[[s$table]], s$groups, contrast = s$contrast, cpus = 2L)
  expect_true(max(abs(res$log2_fc - s$result$log2_fc)) < 1e-2)
  expect_true(max(abs(res$ci_lower - s$result$ci_lower)) < 2e-2)
  expect_true(max(abs(res$ci_upper - s$result$ci_upper)) < 2e-2)
  expect_true(max(abs(res$t - s$result$t) / abs(s$result$t)) < 0.1)
  expect_identical(res$reject, s$result$reject)


  # -----------------------------------------------------------------------
  # The samplers
  # -----------------------------------------------------------------------

  ks <- function (x, ...) suppressWarnings(stats::ks.test(x, ...)$p.value)

  u <- rd(2e5, 1)
  expect_true(all(u > 0 & u < 1))
  expect_true(abs(mean(u) - 0.5) < 4 * sqrt(1 / 12 / 2e5))
  expect_true(ks(u, 'punif') > 1e-3)

  z <- rd(2e5, 2, seed = 3L, stream = 5)
  expect_true(abs(mean(z)) < 4 * sqrt(1 / 2e5))
  expect_true(abs(var(z) - 1) < 4 * sqrt(2 / 2e5))
  expect_true(ks(z, 'pnorm') > 1e-3)

  for (a in c(0.5, 1, 7.5, 1000)) {
    g  <- rd(1e5, 3, a, seed = 7L, stream = a)
    lg <- rd(1e5, 4, a, seed = 8L, stream = a)
    expect_true(abs(mean(g) - a) < 4 * sqrt(a / 1e5), info = a)
    expect_true(abs(var(g) / a - 1) < 4 * sqrt((2 + 6 / a) / 1e5), info = a)
    expect_true(ks(g, 'pgamma', shape = a) > 1e-3, info = a)
    expect_true(ks(exp(lg), 'pgamma', shape = a) > 1e-3, info = a)
    expect_true(abs(mean(lg) - digamma(a)) < 4 * sqrt(trigamma(a) / 1e5), info = a)
  }

  # Tiny shapes stay finite on the log scale (G itself would underflow).
  lg <- rd(1e5, 4, 1e-3, seed = 2L)
  expect_true(all(is.finite(lg)))
  expect_true(abs(mean(lg) - digamma(1e-3)) < 4 * sqrt(trigamma(1e-3) / 1e5))

  # The same seed and stream give the same draws; others don't.
  expect_identical(rd(100, 3, 2.5, 1L, 9), rd(100, 3, 2.5, 1L, 9))
  expect_false(identical(rd(100, 3, 2.5, 1L, 9), rd(100, 3, 2.5, 1L, 10)))
  expect_false(identical(rd(100, 3, 2.5, 1L, 9), rd(100, 3, 2.5, 2L, 9)))


  # -----------------------------------------------------------------------
  # Threads, seeds and input types
  # -----------------------------------------------------------------------

  set.seed(1)
  big <- matrix(rpois(40 * 300, rep(rlnorm(300, 1, 1.5), each = 40)), 40, 300,
    dimnames = list(paste0('s', 1:40), paste0('f', 1:300)) )
  big[, 300] <- 0L                   # a feature absent everywhere
  bg  <- rep(c('u', 'v', 'w'), c(15, 13, 12))

  d1 <- dirmult_ttest_details(big, bg, contrast = c('v', NA), n_draws = 20L, cpus = 1L)
  for (cpus in c(2L, 3L, 4L, 8L))
    expect_identical(
      dirmult_ttest_details(big, bg, contrast = c('v', NA), n_draws = 20L, cpus = cpus),
      d1, info = cpus )

  sp <- Matrix::Matrix(big, sparse = TRUE)
  sl <- slam::as.simple_triplet_matrix(big)
  f  <- function (x, ...) dirmult_ttest_details(x, bg, contrast = c('v', NA), n_draws = 20L, ...)
  expect_identical(f(sp, cpus = 4L), d1)
  expect_identical(f(sl, cpus = 2L), d1)
  expect_identical(f(t(big), margin = 2L, cpus = 2L), d1)
  expect_identical(f(Matrix::t(sp), margin = 2L, cpus = 4L), d1)
  expect_identical(f(big * 1.0, cpus = 1L), d1)

  # Seeds.
  r0 <- dirmult_ttest(big, bg, n_draws = 10L, cpus = 2L)
  expect_identical(dirmult_ttest(big, bg, n_draws = 10L, cpus = 2L), r0)
  r1 <- dirmult_ttest(big, bg, n_draws = 10L, seed = 1L, cpus = 2L)
  expect_false(identical(r1$log2_fc, r0$log2_fc))
  expect_identical(attr(r1, 'seed'), 1L)
  expect_identical(dirmult_ttest(big, bg, n_draws = 10L, seed = 1, cpus = 2L), r1)

  # More draws extend the same sequence: draw d doesn't depend on n_draws.
  a <- dirmult_ttest_details(big, bg, n_draws = 5L, cpus = 1L)$draws
  b <- dirmult_ttest_details(big, bg, n_draws = 8L, cpus = 1L)$draws
  expect_identical(a$diff, b$diff[1:5, ])

  # Samples outside the contrast draw nothing that matters: the pair's
  # results don't change when the other group's samples are relabelled.
  bg2 <- replace(bg, bg == 'w', 'z')
  expect_identical(
    unclass(dirmult_ttest(big, bg,  contrast = c('u', 'v'), n_draws = 5L)),
    unclass(dirmult_ttest(big, bg2, contrast = c('u', 'v'), n_draws = 5L)) )


  # -----------------------------------------------------------------------
  # The result
  # -----------------------------------------------------------------------

  res <- dirmult_ttest(big, bg, contrast = c('w', 'u'), n_draws = 10L, cpus = 2L)
  expect_true(inherits(res, 'ecodive_dirmult_ttest') && is.data.frame(res))
  expect_identical(names(res), c('feature', cols, 'reject'))
  expect_identical(res$feature, colnames(big))
  expect_identical(attr(res, 'contrast'), c(treatment = 'w', reference = 'u'))
  expect_identical(attr(res, 'n'), c(treatment = 12L, reference = 15L))
  expect_identical(attr(res, 'n_draws'), 10L)
  expect_identical(attr(res, 'pseudocount'), 0.5)
  expect_identical(attr(res, 'alpha'), 0.05)
  expect_identical(attr(res, 'conf_level'), 0.95)
  expect_identical(attr(res, 'p_adjust'), 'holm')
  expect_identical(attr(res, 'seed'), 0L)
  expect_true(all(res$ci_lower <= res$log2_fc & res$log2_fc <= res$ci_upper))
  expect_identical(res$p_adjusted, p.adjust(res$p_value, 'holm'))

  # The reversed contrast mirrors the differences, from the same draws.
  rev <- dirmult_ttest(big, bg, contrast = c('u', 'w'), n_draws = 10L, cpus = 2L)
  expect_true(max(abs(rev$log2_fc + res$log2_fc)) < 1e-12)
  expect_true(max(abs(rev$ci_lower + res$ci_upper)) < 1e-12)
  expect_true(max(abs(rev$p_value - res$p_value)) < 1e-12)

  # n_draws = 1, and no feature names.
  one <- dirmult_ttest(unname(big), bg, n_draws = 1L, cpus = 1L)
  expect_identical(one$feature, as.character(1:300))
  expect_true(all(is.finite(one$t)))

  # p_adjust by prefix and 'fdr'.
  expect_identical(attr(dirmult_ttest(big, bg, n_draws = 2L, p_adjust = 'bonf'), 'p_adjust'), 'bonferroni')
  x <- dirmult_ttest(big, bg, n_draws = 2L, p_adjust = 'fdr')
  expect_identical(x$p_adjusted, p.adjust(x$p_value, 'BH'))
  x <- dirmult_ttest(big, bg, n_draws = 2L, p_adjust = 'none')
  expect_identical(x$p_adjusted, x$p_value)

  # Groups: a factor keeps its levels' order (the first level is the
  # default treatment); numbers; a column of `data`; names matching.
  fg <- factor(bg, levels = c('w', 'u', 'v', 'unused'))
  x  <- dirmult_ttest(big, fg, n_draws = 3L)
  expect_identical(attr(x, 'contrast'), c(treatment = 'w', reference = NA))
  expect_identical(attr(x, 'n'), c(treatment = 12L, reference = 28L))
  expect_identical(unclass(x), unclass(dirmult_ttest(big, bg, contrast = c('w', NA), n_draws = 3L)))
  expect_identical(
    unclass(dirmult_ttest(big, match(bg, c('w', 'u', 'v')) * 10, contrast = c(30, 10), n_draws = 3L)),
    { y <- unclass(dirmult_ttest(big, bg, contrast = c('v', 'w'), n_draws = 3L)); attr(y, 'contrast')[] <- c('30', '10'); y } )
  meta <- data.frame(row.names = rev(rownames(big)), grp = rev(bg))
  expect_identical(dirmult_ttest(big, 'grp', meta, n_draws = 3L), dirmult_ttest(big, bg, n_draws = 3L))
  expect_identical(
    dirmult_ttest(big, setNames(rev(bg), rev(rownames(big))), n_draws = 3L),
    dirmult_ttest(big, bg, n_draws = 3L) )


  # -----------------------------------------------------------------------
  # Argument checks
  # -----------------------------------------------------------------------

  expect_error(dirmult_ttest(big), '`groups` is required')
  expect_error(dirmult_ttest(big, rep('a', 40)), 'at least two different')
  expect_error(dirmult_ttest(big, replace(bg, 3, NA)), 'NA')
  expect_error(dirmult_ttest(big, bg[-1]), 'groups')
  expect_error(dirmult_ttest(big, bg, contrast = c('u', 'x')), '"x" is not one of the groups')
  expect_error(dirmult_ttest(big, fg, contrast = c('unused', NA)), '"unused" is not one of the groups')
  expect_error(dirmult_ttest(big, bg, contrast = c('u', 'u')), 'must be different')
  expect_error(dirmult_ttest(big, bg, contrast = 'u'), 'two groups')
  expect_error(dirmult_ttest(big, bg, contrast = c(NA, 'u')), 'two groups')
  expect_error(dirmult_ttest(big, bg, contrast = list('u', 'v')), 'two groups')
  expect_warning(dirmult_ttest(big, replace(bg, 1, 'solo'), contrast = c('solo', 'u'), n_draws = 2L), 'treatment side .* single sample')
  expect_warning(dirmult_ttest(big, replace(bg, 1, 'solo'), contrast = c('u', 'solo'), n_draws = 2L), 'reference side .* single sample')

  for (v in list(0, -1, NA, Inf, 'a', c(1, 2)))
    expect_error(dirmult_ttest(big, bg, pseudocount = v), 'pseudocount')
  for (v in list(0, -1, 1.5, NA, Inf, 'a', c(1, 2)))
    expect_error(dirmult_ttest(big, bg, n_draws = v), 'n_draws')
  for (arg in c('alpha', 'conf_level'))
    for (v in list(0, 1, -0.1, NA, 'a', c(0.1, 0.2)))
      expect_error(
        do.call(dirmult_ttest, setNames(list(big, bg, v), c('counts', 'groups', arg))),
        arg, info = arg )
  expect_error(dirmult_ttest(big, bg, p_adjust = 'nope'), 'p_adjust')
  expect_error(dirmult_ttest(big, bg, seed = 1.5), 'seed')
  expect_error(dirmult_ttest(big, bg, cpus = 0), 'cpus')

  bad <- big; bad[1, 1] <- -1L
  expect_error(dirmult_ttest(bad, bg), 'negative')
  bad <- big * 1.0; bad[1, 1] <- NA
  expect_error(dirmult_ttest(bad, bg), 'NA')
  expect_error(dirmult_ttest(big[, 1, drop = FALSE], bg), 'at least 2 features')
  expect_error(dirmult_ttest(matrix('a', 4, 3), c(1, 1, 2, 2)), 'numeric')

#})
