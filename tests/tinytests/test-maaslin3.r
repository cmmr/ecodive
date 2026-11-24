#test_that("maaslin3", {

  # maaslin3(): maaslin3 1.4.0's maaslin3() (data-raw/expected/maaslin3.r)
  # with its median comparisons' rnorm() draws regenerated here; glmer()'s
  # best fits for the mixed logistic models; threading, input types, the
  # wrappers, edge cases and argument checks.

  ref <- readRDS(file.path('expected', 'maaslin3.rds'))

  # Relative difference, with an absolute floor for values near 0.
  maxrel <- function (a, b, floor = 1e-8) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.na(a), is.na(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }

  # The median comparisons' draws, as maaslin3 makes them after
  # set.seed(seed): for each model with a median comparison (abundance
  # first), each covariate in turn, 10,000 draws of every feature without an
  # error, when the covariate's median is defined.
  draws <- function (d, seed, models) {
    set.seed(seed)
    out <- list()
    for (m in models) {
      tab <- d[[m]]
      out[[m]] <- list()
      for (k in colnames(tab$coef)) {
        ok <- is.na(tab$error[, k])
        p  <- tab$p[ok, k]
        if (is.na(stats::median(tab$coef[ok, k][!is.na(p) & p < 0.95]))) next
        out[[m]][[k]] <- matrix(stats::rnorm(sum(ok) * 10000), sum(ok))
      }
    }
    out
  }

  run <- function (cs, sims = NULL, cpus = 1L) {
    tb <- ref$tables[[cs$table]]
    suppressWarnings(do.call(maaslin3_details, c(
      list(tb$x, cs$formula, tb$meta, cpus = cpus, sims = sims), cs$eco )))
  }


  # -----------------------------------------------------------------------
  # maaslin3 1.4.0
  # -----------------------------------------------------------------------

  for (nm in names(ref$cases)) {

    cs    <- ref$cases[[nm]]
    mixed <- grepl('|', cs$formula, fixed = TRUE)
    small <- isTRUE(cs$eco$small_random_effects)
    mods  <- c(
      if (!isFALSE(cs$eco$median_comparison_abundance)) 'abundance',
      if (isTRUE(cs$eco$median_comparison_prevalence)) 'prevalence' )

    d0  <- run(cs)
    sims <- draws(d0, cs$seed, mods)
    e   <- run(cs, sims)$result

    key <- function (r) paste(r$feature, r$name)
    ra  <- cs$abundance[match(paste(e$feature, e$covariate), key(cs$abundance)), ]
    rp  <- cs$prevalence[match(paste(e$feature, e$covariate), key(cs$prevalence)), ]
    expect_false(anyNA(ra$feature) || anyNA(rp$feature), info = nm)
    expect_identical(nrow(e), nrow(cs$abundance), info = nm)

    # Abundance: lm() to rounding; lmerTest's REML to lme4's optimizer
    # tolerance. The median comparison's p-values come from multcomp's
    # numerical t integral (relative error up to about 1e-8).
    tl <- if (mixed) c(coef = 1e-3, se = 1e-3, p = 1e-3) else c(coef = 1e-10, se = 1e-10, p = 1e-7)
    expect_true(maxrel(e$coef_abundance, ra$coef, 1e-6) < tl[['coef']], info = nm)
    expect_true(maxrel(e$se_abundance, ra$stderr) < tl[['se']], info = nm)
    expect_true(maxrel(e$null_abundance, ra$null_hypothesis, 1e-6) < tl[['coef']], info = nm)
    expect_true(maxrel(e$p_abundance, ra$pval_individual, 1e-12) < tl[['p']], info = nm)
    expect_identical(e$flag_abundance, ra$error, info = nm)
    expect_identical(e$n_nonzero, as.integer(ra$N_not_zero), info = nm)

    if (!mixed || small) {
      # Prevalence: glm()'s iterations, to rounding (without the
      # augmentation, a separated feature's 25 iterations that never
      # converge amplify it).
      tp <- if (isFALSE(cs$eco$augment)) 1e-5 else 1e-8
      expect_true(maxrel(e$coef_prevalence, rp$coef, 1e-6) < tp, info = nm)
      expect_true(maxrel(e$se_prevalence, rp$stderr) < tp, info = nm)
      expect_true(maxrel(e$null_prevalence, rp$null_hypothesis, 1e-6) < 1e-8, info = nm)
      expect_true(maxrel(e$p_prevalence, rp$pval_individual, 1e-12) < 1e-7, info = nm)
      # q-values and joint p-values: with a random intercept in the
      # abundance model, to lme4's tolerance.
      tq <- if (mixed) 1e-3 else 1e-7
      expect_true(maxrel(e$q_abundance, ra$qval_individual, 1e-12) < tq, info = nm)
      expect_true(maxrel(e$q_prevalence, rp$qval_individual, 1e-12) < tq, info = nm)
      expect_true(maxrel(e$p_value, ra$pval_joint, 1e-12) < tq, info = nm)
      expect_true(maxrel(e$p_adjusted, ra$qval_joint, 1e-12) < tq, info = nm)
      expect_identical(e$flag_prevalence, rp$error, info = nm)
    }
    else {
      # Prevalence: glmer()'s best fit (three optimizers), to its tolerance.
      g <- cs$glmer[match(paste(e$feature, e$covariate), key(cs$glmer)), ]
      ok <- !is.na(g$coef)
      expect_true(max(abs(e$coef_prevalence[ok] - g$coef[ok])) < 1e-5, info = nm)
      expect_true(maxrel(e$se_prevalence[ok], g$stderr[ok]) < 1e-4, info = nm)
      expect_identical(is.na(e$flag_prevalence), is.na(rp$error), info = nm)
    }

    # Decisions: maaslin3's significant associations (q_joint <= alpha, or
    # a model's own q <= alpha without an error).
    al <- cs$eco$alpha
    sig_a <- (!is.na(ra$qval_joint) & ra$qval_joint <= al) |
      (!is.na(ra$qval_individual) & ra$qval_individual <= al & is.na(ra$error))
    sig_p <- (!is.na(rp$qval_joint) & rp$qval_joint <= al) |
      (!is.na(rp$qval_individual) & rp$qval_individual <= al &
         (is.na(rp$error) | grepl('^<4 average', rp$error)))
    expect_identical(e$reject, sig_a | sig_p, info = nm)
  }


  # -----------------------------------------------------------------------
  # The median comparison with PCG's draws: maaslin3's p-values to Monte
  # Carlo error; reproducible by seed.
  # -----------------------------------------------------------------------

  cs <- ref$cases$three
  tb <- ref$tables$three
  e  <- run(cs)$result
  ra <- cs$abundance[match(paste(e$feature, e$covariate), paste(cs$abundance$feature, cs$abundance$name)), ]
  # (Relative differences of -log p: the Monte Carlo error of the
  # covariances is amplified in the tails.)
  ok <- !is.na(ra$pval_individual)
  expect_true(maxrel(-log(e$p_abundance[ok]), -log(ra$pval_individual[ok]), 0.1) < 0.02)
  expect_identical(e$coef_abundance, run(cs)$result$coef_abundance)
  e1 <- maaslin3(tb$x, ~ trt + age + sex, tb$meta, seed = 1L, cpus = 1L)
  e2 <- maaslin3(tb$x, ~ trt + age + sex, tb$meta, seed = 2L, cpus = 1L)
  expect_false(identical(e1$p_abundance, e2$p_abundance))
  expect_identical(e1$coef_abundance, e2$coef_abundance)
  expect_identical(e1, maaslin3(tb$x, ~ trt + age + sex, tb$meta, seed = 1L, cpus = 1L))
  expect_true(maxrel(-log(e1$p_abundance), -log(e2$p_abundance), 0.1) < 0.02)

  # The kernel: covariances from given draws, as R computes them.
  set.seed(3)
  b <- rnorm(9); s <- runif(9, 0.1, 1); use <- c(rep(TRUE, 7), FALSE, FALSE)
  z <- matrix(rnorm(9 * 500), 9)
  mc <- .Call(ecodive:::C_maaslin3_median, b, s, use, 500L, 0L, 1L, z, 0.3, 1L)
  xs <- b + s * z
  med <- apply(xs[use, ], 2, stats::median)
  expect_equal(mc$medians, med, tolerance = 1e-14)
  expect_equal(mc$cov, apply(xs, 1, function (v) stats::cov(v, med)), tolerance = 1e-10)
  # PCG draws: standard normals, the same for any number of threads.
  mc1 <- .Call(ecodive:::C_maaslin3_median, b, s, use, 20000L, 5L, 1L, NULL, 0.3, 1L)
  mc4 <- .Call(ecodive:::C_maaslin3_median, b, s, use, 20000L, 5L, 1L, NULL, 0.3, 4L)
  expect_identical(mc1, mc4)
  set.seed(4)
  rmed <- apply(b[use] + s[use] * matrix(rnorm(7 * 20000), 7), 2, stats::median)
  expect_true(abs(mean(mc1$medians) - mean(rmed)) < 4 * stats::sd(rmed) / sqrt(10000))
  expect_true(abs(stats::sd(mc1$medians) / stats::sd(rmed) - 1) < 0.03)


  # -----------------------------------------------------------------------
  # Threads, input types, the wrappers
  # -----------------------------------------------------------------------

  tb <- ref$tables$mixed
  fm <- ~ grp + age
  r1 <- maaslin3(tb$x, fm, tb$meta, n_sims = 500L, cpus = 1L)
  for (cp in c(2L, 4L))
    expect_identical(maaslin3(tb$x, fm, tb$meta, n_sims = 500L, cpus = cp), r1, info = cp)
  rm1 <- maaslin3(tb$x, ~ grp + age + (1 | id), tb$meta, n_sims = 500L, cpus = 1L)
  expect_identical(maaslin3(tb$x, ~ grp + age + (1 | id), tb$meta, n_sims = 500L, cpus = 4L), rm1)
  expect_identical(maaslin3(t(tb$x), fm, tb$meta, margin = 2L, n_sims = 500L, cpus = 2L), r1)
  xd <- tb$x; storage.mode(xd) <- 'double'
  expect_identical(maaslin3(xd, fm, tb$meta, n_sims = 500L, cpus = 1L), r1)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(maaslin3(Matrix::Matrix(tb$x, sparse = TRUE), fm, tb$meta, n_sims = 500L, cpus = 3L), r1)
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(maaslin3(slam::as.simple_triplet_matrix(tb$x), fm, tb$meta, n_sims = 500L, cpus = 1L), r1)
  # Proportions: the same as the counts (percent normalization).
  expect_equal(maaslin3(tb$x / rowSums(tb$x), fm, tb$meta, n_sims = 500L, cpus = 1L)$coef_abundance,
               r1$coef_abundance, tolerance = 1e-12)

  expect_identical(da_test(tb$x, 'maaslin3', 'grp', tb$meta, fm, n_sims = 500L, cpus = 1L),
                   maaslin3(tb$x, fm, tb$meta, groups = 'grp', n_sims = 500L, cpus = 1L))
  expect_identical(da_test(tb$x, 'MaAsLin 3', 'grp', tb$meta, n_sims = 500L, cpus = 1L),
                   maaslin3(tb$x, groups = 'grp', data = tb$meta, n_sims = 500L, cpus = 1L))
  expect_identical(match_da_test('maaslin')$id, 'maaslin3')
  expect_identical(match_da_test('MaAsLin3')$id, 'maaslin3')
  expect_error(match_da_test('MaAsLin2'), 'MaAsLin 3')
  expect_identical(list_da_tests()$effect[list_da_tests()$id == 'maaslin3'], 'coef_abundance')

  # da_sensitivity(): maaslin3 has nothing to vary.
  expect_error(da_sensitivity(tb$x, 'maaslin3', 'grp', tb$meta, pseudocounts = TRUE), 'no `pseudocount`')
  expect_error(da_sensitivity(tb$x, 'maaslin3', 'grp', tb$meta, zeros = TRUE), 'presence')
  expect_error(da_sensitivity(tb$x, 'maaslin3', 'grp', tb$meta, scale_sds = TRUE), 'no `scale_sd`')


  # -----------------------------------------------------------------------
  # Edge cases
  # -----------------------------------------------------------------------

  tb <- ref$tables$two
  x  <- tb$x
  x[, 1] <- 0L                      # absent everywhere
  x[, 2] <- x[, 2] + 1L             # present everywhere
  x[, 3] <- c(5L, rep(0L, 39))      # in one sample
  x[, 4] <- c(rep(3L, 10), rep(0L, 30))   # only in group a
  r <- maaslin3(x, ~ grp, tb$meta, n_sims = 500L, cpus = 1L)
  expect_identical(nrow(r), ncol(x))
  expect_identical(r$flag_abundance[1:4], c(
    'All linear values are the same', NA, 'All linear values are the same',
    'contrasts can be applied only to factors with 2 or more levels' ))
  expect_identical(r$flag_prevalence[1:2], rep('All logistic values are the same', 2))
  expect_true(all(is.nan(c(r$coef_abundance[c(1, 3, 4)], r$coef_prevalence[1:2]))))
  expect_true(is.finite(r$coef_prevalence[3]) && is.finite(r$coef_abundance[2]))
  expect_identical(r$n_nonzero[1:4], c(0L, 40L, 1L, 10L))
  # The joint p-value of a feature with one model is that model's.
  expect_identical(r$p_value[2], r$p_abundance[2])
  expect_identical(r$p_value[3], r$p_prevalence[3])
  expect_true(is.nan(r$p_value[1]))

  # Separation: the augmentation keeps the prevalence estimate finite.
  xs <- tb$x; xs[1:20, 5] <- 0L; xs[21:40, 5] <- xs[21:40, 5] + 1L
  r <- maaslin3(xs, ~ grp, tb$meta, n_sims = 500L, cpus = 1L)
  expect_true(is.finite(r$coef_prevalence[5]) && abs(r$coef_prevalence[5]) < 10)
  r <- maaslin3(xs, ~ grp, tb$meta, augment = FALSE, warn_prevalence = FALSE, n_sims = 500L, cpus = 1L)
  expect_true(r$coef_prevalence[5] > 20)

  # A rank-deficient design: a warning, and the aliased column flagged.
  md <- tb$meta; md$z <- rnorm(40); md$z2 <- 2 * md$z
  expect_warning(r <- maaslin3(tb$x, ~ grp + z + z2, md, standardize = FALSE, n_sims = 500L, cpus = 1L), 'rank deficient')
  expect_true(all(r$flag_abundance[r$covariate == 'z2'] == 'Fitting error (NA p-value returned from fitting procedure)'))
  expect_true(all(is.finite(r$coef_abundance[r$covariate == 'grpb'])))

  # Few samples per random intercept group: the prevalence rows' flag.
  tb <- ref$tables$small
  r <- maaslin3(tb$x, ~ trt + age + (1 | id), tb$meta, n_sims = 500L, cpus = 1L)
  expect_true(all(grepl('^<4 average', r$flag_prevalence) | !is.na(r$flag_prevalence)))
  expect_identical(attr(r, 'random'), 'id')
  rs <- maaslin3(tb$x, ~ trt + age + (1 | id), tb$meta, small_random_effects = TRUE, n_sims = 500L, cpus = 1L)
  expect_false(any(grepl('^<4 average', rs$flag_prevalence)))
  expect_identical(rs$coef_abundance, r$coef_abundance)

  tb <- ref$tables$mixed
  x  <- tb$x
  pm <- tb$meta
  expect_error(maaslin3(x[c(1, 2, 30), ], ~ grp + age, pm[c(1, 2, 30), ]), 'more samples than design columns')
  mna <- pm; mna$age[3] <- NA
  expect_error(maaslin3(x, ~ age, mna), 'NA')
  expect_error(maaslin3(x, ~ age, transform(pm, age = 1)), 'standardized')
  expect_error(maaslin3(x, ~ 1, pm), 'besides the intercept')
  expect_error(maaslin3(-x, ~ grp, pm), 'negative')
  expect_error(maaslin3(x, ~ grp, pm, norm = 'clr'), 'norm')
  expect_error(maaslin3(x, ~ grp, pm, transform = 'plog'), 'transform')
  expect_error(maaslin3(x, ~ grp, pm, norm = 'rclr', warn_prevalence = FALSE), 'rclr values are logs')
  expect_error(maaslin3(x, ~ grp, pm, norm = 'rclr', transform = 'none'), 'warn_prevalence = FALSE')
  expect_error(maaslin3(x, ~ grp, pm, transform = 'none'), 'warn_prevalence = FALSE')
  expect_error(maaslin3(x, ~ grp, pm, small_random_effects = TRUE), 'random intercept')
  expect_error(maaslin3(x, ~ grp + (1 | age), pm), 'fewer values than there are samples')
  expect_error(maaslin3(x, ~ grp + (age | id), pm), 'one random intercept')
  expect_error(maaslin3(x, ~ grp, pm, augment = NA), '`augment`')
  expect_error(maaslin3(x, ~ grp, pm, median_comparison_abundance_threshold = -1), 'threshold')
  expect_error(maaslin3(x, ~ grp, pm, n_sims = 1), 'n_sims')
  expect_error(maaslin3(x, ~ grp, pm, alpha = 2), 'alpha')
  expect_error(maaslin3(x, ~ grp, pm, p_adjust = 'nope'))
  expect_error(maaslin3(x), 'formula')


  # -----------------------------------------------------------------------
  # Attributes, print and plot
  # -----------------------------------------------------------------------

  r <- maaslin3(x, ~ grp + age, pm, groups = 'grp', n_sims = 500L, cpus = 1L)
  expect_identical(class(r), c('ecodive_maaslin3', 'data.frame'))
  expect_identical(names(r), c(
    'feature', 'covariate', 'coef_abundance', 'se_abundance', 'dof_abundance',
    'null_abundance', 'p_abundance', 'q_abundance', 'coef_prevalence',
    'se_prevalence', 'null_prevalence', 'p_prevalence', 'q_prevalence',
    'n_nonzero', 'flag_abundance', 'flag_prevalence', 'p_value', 'p_adjusted',
    'reject', 'reject_abundance', 'reject_prevalence' ))
  expect_identical(unique(r$covariate), c('grpb', 'age'))
  expect_identical(attr(r, 'method'), 'maaslin3')
  expect_identical(attr(r, 'log_base'), 2)
  expect_identical(attr(r, 'alpha'), 0.05)
  expect_identical(attr(r, 'p_adjust'), 'holm')
  expect_identical(attr(r, 'groups'), 'grp')
  expect_identical(names(attr(r, 'mean')), colnames(x))
  expect_equal(unname(attr(r, 'prevalence')), unname(colMeans(x > 0)))
  expect_true(all(r$null_prevalence == 0))
  expect_true(all(r$reject == (r$p_adjusted <= 0.05 & !is.na(r$p_adjusted)) | r$reject_abundance | r$reject_prevalence))
  expect_identical(attr(maaslin3(x, ~ grp, pm, norm = 'rclr', transform = 'none',
                                 warn_prevalence = FALSE, n_sims = 500L), 'log_base'), exp(1))
  expect_null(attr(r, 'display'))
  expect_true(any(grepl('MaAsLin 3', utils::capture.output(print(r)))))
  # subtract_median: the coefficients minus the median, the same p-values.
  rsub <- maaslin3(x, ~ grp + age, pm, subtract_median = TRUE, n_sims = 500L, cpus = 1L)
  expect_equal(rsub$coef_abundance, r$coef_abundance - r$null_abundance, tolerance = 1e-14)
  expect_identical(rsub$p_abundance, r$p_abundance)
  # Without a median comparison: tested against 0.
  r0 <- maaslin3(x, ~ grp + age, pm, median_comparison_abundance = FALSE, cpus = 1L)
  expect_true(all(r0$null_abundance == 0))
  expect_equal(r0$p_abundance, 2 * stats::pt(-abs(r0$coef_abundance / r0$se_abundance), r0$dof_abundance))

  expect_identical(ecodive:::da_plot_title(r), 'Volcano plot of MaAsLin 3')
  grDevices::pdf(NULL)
  expect_identical(withVisible(plot(r))$visible, FALSE)
  expect_silent(plot(r))
  expect_silent(plot(r, 'age', model = 'prevalence', labels = 0, main = 'x'))
  for (ty in c('ma', 'forest', 'rank')) {
    expect_silent(plot(r, type = ty))
    expect_silent(plot(r, type = ty, model = 'prevalence'))
  }
  expect_silent(plot(rm1, type = 'forest', n_features = 5, conf_level = 0.9))
  expect_identical(ecodive:::da_plot_title(r, 'forest'), 'Forest plot of MaAsLin 3')
  expect_error(plot(r, 'nope'), '`covariate`')
  expect_error(plot(r, model = 'nope'))
  expect_error(plot(r, type = 'effect'))
  expect_error(plot(r, type = 'forest', conf_level = 1), 'conf_level')
  grDevices::dev.off()

#})
