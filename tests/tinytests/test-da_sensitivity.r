#test_that("da_sensitivity", {

  # da_sensitivity(): every run equals the direct call with that setting;
  # ANCOM-BC's pseudocount runs are the fits behind its `passed_ss`; the
  # runs, the table and the summary; threading; errors; plots.

  set.seed(3)
  n <- 16L; D <- 12L
  x <- matrix(rpois(n * D, 12), n, D, dimnames = list(paste0('S', 1:n), paste0('F', 1:D)))
  x[sample(length(x), 30)] <- 0L
  x[1:8, 1:3] <- x[1:8, 1:3] * 6L
  x[9:16, 4]  <- x[9:16, 4] * 6L
  x[, 5]      <- c(3L, 4L, 2L, 5L, 3L, 2L, 4L, 3L, 0L, 0L, 1L, 0L, 0L, 2L, 0L, 0L)
  meta <- data.frame(
    row.names = rownames(x),
    grp = rep(c('a', 'b'), each = 8),
    age = c(31, 45, 52, 38, 29, 61, 44, 50, 36, 58, 41, 47, 33, 55, 39, 60),
    stringsAsFactors = FALSE )
  xz <- replace_zeros(x, 'mult')

  # The rows of run `k` (non-intercept), as the method gives them.
  rows_of <- function (res, k) res[res$run == k, , drop = FALSE]
  fit_rows <- function (f) if (is.null(f$covariate)) f else f[f$covariate != '(Intercept)', , drop = FALSE]
  same_rows <- function (res, k, f, eff) {
    r <- rows_of(res, k)
    g <- fit_rows(f)
    identical(r$feature, as.character(g$feature)) &&
      identical(r$reject, g$reject) &&
      identical(r$effect, if (is.na(eff)) rep(NA_real_, nrow(g)) else g[[eff]]) &&
      identical(r$p_adjusted, if (is.null(g$p_adjusted)) rep(NA_real_, nrow(g)) else g$p_adjusted)
  }


  # -----------------------------------------------------------------------
  # Each run is the direct call with its setting.
  # -----------------------------------------------------------------------

  pcs <- c(0.1, 0.5, 1)

  # ANCOM
  res <- suppressWarnings(da_sensitivity(x, 'ancom', 'grp', meta, pseudocounts = TRUE, zeros = 'mult', cpus = 1))
  fits <- attr(res, 'fits')
  expect_identical(length(fits), 5L)
  expect_identical(fits[[1]], suppressWarnings(ancom(x, 'grp', meta, cpus = 1)))
  for (k in 1:3) {
    expect_identical(fits[[k + 1]], ancom(x, 'grp', meta, pseudocount = pcs[[k]], cpus = 1), info = k)
    expect_true(same_rows(res, k, fits[[k + 1]], NA), info = k)
  }
  expect_identical(fits[[5]], ancom(xz, 'grp', meta, cpus = 1))
  expect_true(same_rows(res, 4L, fits[[5]], NA))
  expect_true(all(is.na(res$effect)) && all(is.na(res$p_adjusted)) && all(is.na(res$covariate)))
  expect_identical(attr(res, 'effect_name'), NA_character_)
  expect_null(attr(res, 'covariate'))

  # ANCOM-BC2: a model method, with a covariate.
  res <- suppressWarnings(da_sensitivity(
    x, 'ancombc2', 'grp', meta, ~ grp + age, pseudocounts = TRUE, zeros = TRUE, cpus = 1 ))
  fits <- attr(res, 'fits')
  base <- suppressWarnings(ancombc2(x, ~ grp + age, meta, groups = 'grp', cpus = 1))
  expect_identical(fits[[1]], base)
  for (k in 1:3) {
    expect_identical(
      fits[[k + 1]],
      ancombc2(x, ~ grp + age, meta, groups = 'grp', pseudocount = pcs[[k]], sensitivity = FALSE, cpus = 1),
      info = k )
    expect_true(same_rows(res, k, fits[[k + 1]], 'log_fc'), info = k)
  }
  expect_identical(
    fits[[5]],
    suppressWarnings(ancombc2(xz, ~ grp + age, meta, groups = 'grp', sensitivity = FALSE, cpus = 1)) )
  expect_true(same_rows(res, 0L, base, 'log_fc'))
  expect_identical(unique(res$covariate), c('grpb', 'age'))
  expect_identical(attr(res, 'covariate'), 'grpb')
  expect_identical(attr(res, 'log_base'), exp(1))
  expect_identical(attr(res, 'alpha'), 0.05)

  # The Dirichlet-multinomial t-test: the prior and the scale model, with
  # the base run's seed in every run.
  res <- da_sensitivity(
    x, 'dirmult_ttest', 'grp', meta, n_draws = 16L, seed = 7L,
    pseudocounts = c(0.1, 1), scale_sds = c(0, 0.5), cpus = 1 )
  fits <- attr(res, 'fits')
  expect_identical(fits[[1]], dirmult_ttest(x, 'grp', meta, n_draws = 16L, seed = 7L, cpus = 1))
  expect_identical(fits[[2]], dirmult_ttest(x, 'grp', meta, n_draws = 16L, seed = 7L, pseudocount = 0.1, cpus = 1))
  expect_identical(fits[[3]], dirmult_ttest(x, 'grp', meta, n_draws = 16L, seed = 7L, pseudocount = 1, cpus = 1))
  expect_identical(fits[[4]], dirmult_ttest(x, 'grp', meta, n_draws = 16L, seed = 7L, scale_sd = 0, cpus = 1))
  expect_identical(fits[[5]], dirmult_ttest(x, 'grp', meta, n_draws = 16L, seed = 7L, scale_sd = 0.5, cpus = 1))
  for (k in 0:4) expect_true(same_rows(res, k, fits[[k + 1]], 'log2_fc'), info = k)
  expect_identical(attr(res, 'log_base'), 2)
  expect_identical(res$parameter, rep(c('base', 'pseudocount', 'pseudocount', 'scale_sd', 'scale_sd'), each = D))

  # da_lm, da_ttest, da_wilcoxon (base pseudocount given; the zeros runs
  # drop it).
  res <- da_sensitivity(x, 'da_lm', 'grp', meta, ~ grp + age, pseudocount = 0.5, pseudocounts = 1, zeros = 'mult', cpus = 1)
  fits <- attr(res, 'fits')
  expect_identical(fits[[1]], da_lm(x, ~ grp + age, meta, 'grp', pseudocount = 0.5, cpus = 1))
  expect_identical(fits[[2]], da_lm(x, ~ grp + age, meta, 'grp', pseudocount = 1, cpus = 1))
  expect_identical(fits[[3]], da_lm(xz, ~ grp + age, meta, 'grp', cpus = 1))
  for (k in 0:2) expect_true(same_rows(res, k, fits[[k + 1]], 'estimate'), info = k)

  for (m in c('da_ttest', 'da_wilcoxon')) {
    eff <- list_da_tests('list')[[m]]$effect
    f   <- get(m)
    res <- da_sensitivity(x, m, 'grp', meta, pseudocount = 1, norm = 'log', pseudocounts = c(0.5, 2), zeros = TRUE, cpus = 1)
    fits <- attr(res, 'fits')
    expect_identical(fits[[1]], f(x, 'grp', meta, pseudocount = 1, norm = 'log', cpus = 1), info = m)
    expect_identical(fits[[2]], f(x, 'grp', meta, pseudocount = 0.5, norm = 'log', cpus = 1), info = m)
    expect_identical(fits[[3]], f(x, 'grp', meta, pseudocount = 2, norm = 'log', cpus = 1), info = m)
    expect_identical(fits[[4]], f(xz, 'grp', meta, norm = 'log', cpus = 1), info = m)
    for (k in 0:3) expect_true(same_rows(res, k, fits[[k + 1]], eff), info = paste(m, k))
  }

  # linda (pseudocounts with zero_handling = 'pseudocount'; the zeros run
  # gives it proportions without zeros).
  res <- da_sensitivity(x, 'linda', 'grp', meta, ~ grp + age, zero_handling = 'pseudocount', pseudocount = 0.5,
                        pseudocounts = c(1, 2), zeros = 'mult', cpus = 1)
  fits <- attr(res, 'fits')
  expect_identical(fits[[1]], linda(x, ~ grp + age, meta, 'grp', zero_handling = 'pseudocount', pseudocount = 0.5, cpus = 1))
  expect_identical(fits[[2]], linda(x, ~ grp + age, meta, 'grp', zero_handling = 'pseudocount', pseudocount = 1, cpus = 1))
  expect_identical(fits[[3]], linda(x, ~ grp + age, meta, 'grp', zero_handling = 'pseudocount', pseudocount = 2, cpus = 1))
  expect_identical(fits[[4]], linda(xz, ~ grp + age, meta, 'grp', zero_handling = 'pseudocount', cpus = 1))
  for (k in 0:3) expect_true(same_rows(res, k, fits[[k + 1]], 'log2_fc'), info = k)
  expect_identical(attr(res, 'log_base'), 2)
  expect_identical(attr(res, 'covariate'), 'grpb')
  expect_error(da_sensitivity(x, 'linda', 'grp', meta, zero_handling = 'imputation', pseudocounts = TRUE), 'imputation')

  # A partial argument name is completed, so the run's setting replaces it.
  res <- da_sensitivity(x, 'da_ttest', 'grp', meta, pseudo = 1, pseudocounts = 2, cpus = 1)
  expect_identical(attr(res, 'fits')[[1]], da_ttest(x, 'grp', meta, pseudocount = 1, cpus = 1))
  expect_identical(attr(res, 'fits')[[2]], da_ttest(x, 'grp', meta, pseudocount = 2, cpus = 1))

  # Margin 2 (features as rows), with zeros.
  res <- da_sensitivity(t(x), 'da_ttest', 'grp', meta, margin = 2L, pseudocount = 1, zeros = 'mult', cpus = 1)
  expect_identical(attr(res, 'fits')[[2]], da_ttest(replace_zeros(t(x), 'mult', margin = 2L), 'grp', meta, margin = 2L, cpus = 1))
  expect_identical(res$effect, da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocount = 1, zeros = 'mult', cpus = 1)$effect)


  # -----------------------------------------------------------------------
  # ANCOM-BC's own sensitivity analysis.
  # -----------------------------------------------------------------------

  for (m in c('ancombc', 'ancombc2')) {
    f   <- get(m)
    own <- suppressWarnings(f(x, data = meta, groups = 'grp', sensitivity = pcs, cpus = 1))
    # Reused from the base run (sensitivity = pcs) ...
    r1 <- suppressWarnings(da_sensitivity(x, m, 'grp', meta, sensitivity = pcs, pseudocounts = TRUE, cpus = 1))
    # ... or fitted by da_sensitivity() (base without one).
    r2 <- suppressWarnings(da_sensitivity(x, m, 'grp', meta, sensitivity = FALSE, pseudocounts = TRUE, cpus = 1))
    expect_identical(attr(r1, 'fits')[-1], attr(own, 'sensitivity')$fits, info = m)
    expect_identical(attr(r2, 'fits')[-1], attr(own, 'sensitivity')$fits, info = m)
    keep <- own$covariate != '(Intercept)'
    expect_identical(attr(r1, 'summary')$passed, own$passed_ss[keep], info = m)
    expect_identical(attr(r2, 'summary')$passed, own$passed_ss[keep], info = m)
    expect_identical(attr(r1, 'summary')$robust, own$reject_robust[keep], info = m)
    # The table's `reject` is the method's own, not reject_robust.
    expect_identical(r1$reject[r1$run == 0L], own$reject[keep], info = m)
  }
  expect_true(any(!attr(r1, 'summary')$passed))  # some decisions change
  # ancombc2()'s default (sensitivity = TRUE) and pseudocounts = TRUE.
  own <- suppressWarnings(ancombc2(x, data = meta, groups = 'grp', cpus = 1))
  res <- suppressWarnings(da_sensitivity(x, 'ancombc2', 'grp', meta, pseudocounts = TRUE, cpus = 1))
  expect_identical(attr(res, 'fits')[-1], attr(own, 'sensitivity')$fits)
  expect_identical(attr(res, 'summary')$passed, own$passed_ss[own$covariate != '(Intercept)'])
  # Runs other than the base have no sensitivity analysis of their own.
  res <- suppressWarnings(da_sensitivity(x, 'ancombc2', 'grp', meta, zeros = 'mult', cpus = 1))
  expect_null(attr(attr(res, 'fits')[[2]], 'sensitivity'))
  expect_false(is.null(attr(attr(res, 'fits')[[1]], 'sensitivity')))


  # -----------------------------------------------------------------------
  # Runs, table and summary.
  # -----------------------------------------------------------------------

  res  <- da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = TRUE, zeros = TRUE, pseudocount = 1, cpus = 1)
  runs <- attr(res, 'runs')
  nz   <- length(list_zero_replacements('id'))
  expect_identical(nrow(runs), 4L + nz)
  expect_identical(runs$run, 0:(3L + nz))
  expect_identical(runs$parameter, c('base', rep('pseudocount', 3), rep('zeros', nz)))
  expect_identical(runs$setting[1:4], c('base', 'pseudocount = 0.1', 'pseudocount = 0.5', 'pseudocount = 1'))
  expect_identical(runs$zeros[5:(4 + nz)], list_zero_replacements('id'))
  expect_identical(runs$setting[[5]], 'zeros = mult')
  expect_identical(names(res), c('feature', 'covariate', 'run', 'parameter', 'setting', 'effect', 'p_adjusted', 'reject'))
  expect_identical(nrow(res), D * (4L + nz))
  expect_identical(class(res), c('ecodive_da_sensitivity', 'data.frame'))
  expect_identical(attr(res, 'method'), 'da_ttest')
  expect_identical(attr(res, 'effect_name'), 'diff')

  s <- attr(res, 'summary')
  expect_identical(names(s), c('feature', 'covariate', 'n_runs', 'n_reject', 'prop_reject', 'passed', 'robust', 'sign_stable'))
  expect_identical(s$feature, colnames(x))
  nrej <- as.vector(tapply(res$reject, factor(res$feature, colnames(x)), sum))
  expect_identical(s$n_reject, nrej)
  expect_identical(s$n_runs, rep(4L + nz, D))
  expect_identical(s$prop_reject, nrej / (4L + nz))
  expect_identical(s$passed, nrej == 0L | nrej == 4L + nz)
  expect_identical(s$robust, s$passed & res$reject[res$run == 0L])
  sgn <- tapply(res$effect, factor(res$feature, colnames(x)), function (e) length(unique(sign(e[is.finite(e)]))) == 1L)
  expect_identical(s$sign_stable, as.vector(sgn))
  expect_true(any(s$robust))

  # Crossed: the zero treatments (alternatives) x the scale_sds.
  res  <- da_sensitivity(
    x, 'dirmult_ttest', 'grp', meta, n_draws = 8L, pseudocounts = c(0.2, 0.7),
    scale_sds = c(0, 0.3, 1), crossed = TRUE, cpus = 1 )
  runs <- attr(res, 'runs')
  expect_identical(nrow(runs), 1L + 2L * 3L)
  expect_identical(runs$parameter, c('base', rep('crossed', 6)))
  expect_identical(runs$pseudocount, c(NA, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7))
  expect_identical(runs$scale_sd, c(NA, 0, 0.3, 1, 0, 0.3, 1))
  expect_identical(runs$setting[[3]], 'pseudocount = 0.2, scale_sd = 0.3')
  expect_identical(
    attr(res, 'fits')[[7]],
    dirmult_ttest(x, 'grp', meta, n_draws = 8L, pseudocount = 0.7, scale_sd = 1, cpus = 1) )
  res  <- da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocount = 1, pseudocounts = c(0.5, 2), zeros = 'mult', crossed = TRUE, cpus = 1)
  expect_identical(attr(res, 'runs')$parameter, c('base', 'pseudocount', 'pseudocount', 'zeros'))
  res  <- da_sensitivity(x, 'dirmult_ttest', 'grp', meta, n_draws = 8L, scale_sds = TRUE, crossed = TRUE, cpus = 1)
  expect_identical(attr(res, 'runs')$scale_sd, c(NA, 0, 0.25, 0.5, 1))
  expect_identical(attr(res, 'runs')$parameter, c('base', rep('scale_sd', 4)))


  # -----------------------------------------------------------------------
  # Identical for any cpus.
  # -----------------------------------------------------------------------

  for (m in c('ancombc2', 'dirmult_ttest', 'da_lm', 'ancom')) {
    grids <- switch(m,
      dirmult_ttest = list(pseudocounts = TRUE, scale_sds = c(0, 1)),
      list(pseudocounts = TRUE, zeros = TRUE) )
    run <- function (cpus) suppressWarnings(do.call(
      da_sensitivity,
      c(list(x, m, 'grp', meta, cpus = cpus), if (m == 'dirmult_ttest') list(n_draws = 16L), grids) ))
    r1 <- run(1L)
    expect_identical(run(2L), r1, info = m)
    expect_identical(run(4L), r1, info = m)
  }


  # -----------------------------------------------------------------------
  # Warnings: those of the base run as usual, the others once each.
  # -----------------------------------------------------------------------

  w <- character(0)
  res <- withCallingHandlers(
    da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = c(1, 2), zeros = 'mult', cpus = 1),
    warning = function (e) { w <<- c(w, conditionMessage(e)); invokeRestart('muffleWarning') })
  expect_identical(length(w), 1L)
  expect_true(grepl('Zeros detected', w[[1]]))
  xp <- x / rowSums(x)
  expect_warning(
    da_sensitivity(xp, 'da_ttest', 'grp', meta, pseudocounts = TRUE, cpus = 1),
    'not whole numbers' )


  # -----------------------------------------------------------------------
  # Errors.
  # -----------------------------------------------------------------------

  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta), 'at least one of')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = FALSE, zeros = NULL), 'at least one of')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, scale_sds = TRUE), 'has no `scale_sd` argument.*dirmult_ttest')
  expect_error(da_sensitivity(x, 'dirmult_ttest', 'grp', meta, zeros = TRUE), 'models counts')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = -1), 'positive numbers')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = 0), 'positive numbers')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = c(1, NA)), 'positive numbers')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = c(1, 1)), 'repeat')
  expect_error(da_sensitivity(x, 'dirmult_ttest', 'grp', meta, scale_sds = -0.5), 'non-negative')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, zeros = 'nonesuch'), 'does not match')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, zeros = 1), 'names of zero replacement')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, zeros = c('mult', 'multRepl')), 'twice')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocounts = 1, crossed = NA), '`crossed`')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, norm = 'rclr', pseudocounts = 1), 'adds no pseudocount')
  expect_error(da_sensitivity(x, 'da_ttest', 'grp', meta, norm = 'percent', pseudocounts = 1), 'adds no pseudocount')
  expect_error(da_sensitivity(x, 'bray', 'grp', meta, pseudocounts = 1))
  expect_error(
    da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocount = 1, zeros = 'mult', contrast = c('a', 'zz')),
    'contrast' )


  # -----------------------------------------------------------------------
  # Printing and plots.
  # -----------------------------------------------------------------------

  res <- da_sensitivity(x, 'da_ttest', 'grp', meta, pseudocount = 1, pseudocounts = c(0.5, 2), zeros = TRUE, cpus = 1)
  out <- capture.output(res)
  expect_true(any(grepl('Sensitivity of Welch t-test', out, fixed = TRUE)))
  expect_true(length(out) < 30L)

  expect_identical(da_sens_plot_title(res), 'Sensitivity of Welch t-test')
  grDevices::pdf(NULL)
  expect_identical(withVisible(plot(res)), list(value = res, visible = FALSE))
  expect_silent(plot(res, 'traces'))
  expect_silent(plot(res, 'traces', parameter = 'zeros'))
  expect_silent(plot(res, 'tiles', parameter = 'zeros', n_features = 2))
  expect_error(plot(res, 'traces', parameter = 'scale_sd'), '`parameter` must be one of')
  expect_error(plot(res, parameter = 'nonesuch'), '`parameter` must be one of')
  expect_error(plot(res, n_features = 0), 'at least 1')

  r2 <- suppressWarnings(da_sensitivity(x, 'ancombc2', 'grp', meta, ~ grp + age, pseudocounts = TRUE, cpus = 1))
  expect_silent(plot(r2))
  expect_silent(plot(r2, 'traces', covariate = 'age'))
  expect_error(plot(r2, covariate = '(Intercept)'), '`covariate` must be one of')

  r3 <- da_sensitivity(x, 'dirmult_ttest', 'grp', meta, n_draws = 8L, pseudocounts = c(0.2, 0.7), scale_sds = c(0, 1), crossed = TRUE, cpus = 1)
  expect_silent(plot(r3))
  expect_silent(plot(r3, 'traces'))
  expect_silent(plot(r3, 'traces', parameter = 'scale_sd'))

  r4 <- suppressWarnings(da_sensitivity(x, 'ancom', 'grp', meta, pseudocounts = TRUE, cpus = 1))
  expect_silent(plot(r4))
  expect_error(plot(r4, 'traces'), 'gives none')
  expect_identical(da_sens_plot_title(r4), 'Sensitivity of Analysis of Composition of Microbiomes')
  grDevices::dev.off()

#})
