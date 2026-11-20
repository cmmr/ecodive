#test_that("linda", {

  # linda(): MicrobiomeStat 1.4's linda() (data-raw/expected/linda.r) on
  # counts, proportions, factors, numeric and logical covariates, an
  # interaction, both zero handlings and the adaptive choice; modeest's
  # shorth() and meanshift(), LinDA's mode; threading, input types, the
  # wrappers, invariances, edge cases and argument checks.

  ref <- readRDS(file.path('expected', 'linda.rds'))

  # Relative difference, with an absolute floor for values near 0.
  maxrel <- function (a, b, floor = 1e-8) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.nan(a), is.nan(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }


  # -----------------------------------------------------------------------
  # MicrobiomeStat
  # -----------------------------------------------------------------------

  for (cs in names(ref$cases)) {

    r  <- ref$cases[[cs]]
    tb <- ref$tables[[r$table]]
    zh <- switch(r$zero, 'pseudo-count' = 'pseudocount', proportion = 'adaptive', r$zero)

    e <- linda(
      tb$counts, r$formula, tb$meta, zero_handling = zh, pseudocount = 0.5,
      outlier_pct = if (r$winsor) 0.03 else 0, p_adjust = 'BH', cpus = 1L )

    # Mixed models: lme4's optimizer stops short of the exact REML optimum
    # (see the lme4 section below), so compare at its tolerance.
    mixed <- grepl('|', r$formula, fixed = TRUE)
    tl <- if (mixed) c(fc = 1e-3, se = 1e-4, p = 1e-3, dof = 2e-3, bias = 1e-5, cov = 1e-4) else
                     c(fc = 1e-9, se = 1e-10, p = 1e-9, dof = 0, bias = 1e-12, cov = 1e-10)
    fl <- if (mixed) 1e-3 else 1e-8

    o <- r$result
    key <- paste(e$feature, e$covariate)
    i <- match(paste(o$feature, o$covariate), key)
    expect_false(anyNA(i), info = cs)
    expect_identical(nrow(e), nrow(o), info = cs)
    e <- e[i, ]

    expect_identical(attr(e, 'zero_handling'), if (r$zero == 'proportion') 'half_min' else r$chosen, info = cs)
    expect_true(maxrel(e$log2_fc, o$log2_fc, fl) < tl[['fc']], info = cs)
    expect_true(maxrel(e$se, o$se) < tl[['se']], info = cs)
    expect_true(maxrel(e$t, o$t, fl) < tl[['fc']], info = cs)
    expect_true(maxrel(e$p_value, o$p_value, fl) < tl[['p']], info = cs)
    expect_true(maxrel(e$p_adjusted, o$p_adjusted, fl) < tl[['p']], info = cs)
    expect_identical(e$reject, o$reject, info = cs)
    if (mixed) expect_true(maxrel(e$dof, o$dof) < tl[['dof']], info = cs)
    else       expect_identical(e$dof, o$dof, info = cs)

    # Steps 1 to 4: the biased coefficients; step 5: the mode.
    bias <- attr(e, 'bias')
    expect_true(maxrel(e$log2_fc + bias[e$covariate], o$biased, fl) < tl[['fc']], info = cs)
    expect_true(max(abs(bias - r$bias)) < tl[['bias']], info = cs)

    # Covariances, and each feature's residual standard deviation.
    feats <- unique(o$feature)
    if (length(r$cov12)) {
      v <- unique(o$covariate)
      expect_true(maxrel(attr(e, 'vcov')[feats, v[[1]], v[[2]]], r$cov12) < tl[['cov']], info = cs)
    }
    if (length(r$sig))
      expect_true(maxrel(attr(e, 'resid_sd')[feats], r$sig) < 1e-10, info = cs)
    if (mixed) expect_identical(attr(e, 'random'), 'subj', info = cs)
  }


  # -----------------------------------------------------------------------
  # Random intercepts: lme4's lmer() and lmerTest's degrees of freedom
  # -----------------------------------------------------------------------

  for (k in seq_along(ref$lmm)) {
    m  <- ref$lmm[[k]]
    Y  <- matrix(m$y)
    # At lme4's estimate of theta: the same REML fit.
    at <- .Call(C_feature_lmm, m$X, Y, m$groups, m$theta)
    expect_true(maxrel(at$coef, m$coef) < 1e-10, info = k)
    expect_true(maxrel(at$se, m$se) < 1e-10, info = k)
    expect_true(maxrel(at$df, m$df) < 1e-5, info = k)     # lmerTest: numerical derivatives
    expect_true(abs(at$deviance - m$dev) < 1e-9 * abs(m$dev), info = k)
    # The exact optimum: never worse than lme4's, and close to it.
    op <- .Call(C_feature_lmm, m$X, Y, m$groups, NULL)
    expect_true(op$deviance <= m$dev + 1e-10 * abs(m$dev), info = k)
    expect_true(abs(op$theta - m$theta) < 1e-3 * (1 + m$theta), info = k)
    expect_true(maxrel(op$coef, m$coef, 1e-3) < 1e-5, info = k)
    expect_true(maxrel(op$df, m$df) < 1e-3, info = k)
  }


  # -----------------------------------------------------------------------
  # The mode: modeest's shorth() and meanshift()
  # -----------------------------------------------------------------------

  for (k in seq_along(ref$modes)) {
    m <- ref$modes[[k]]
    v <- m$values
    s <- ecodive:::linda_shorth(v)
    expect_true(abs(s - m$shorth) <= 1e-13 * max(1, abs(m$shorth)), info = k)
    md <- .Call(C_linda_mode, matrix(v), s, stats::bw.nrd0(v), sqrt(.Machine$double.eps), 1000L)$mode
    expect_true(abs(md - m$mode) <= 1e-12 * max(abs(m$mode), abs(v)), info = k)
  }


  # -----------------------------------------------------------------------
  # Threads, input types, the wrappers
  # -----------------------------------------------------------------------

  tb   <- ref$tables$sparse
  x    <- tb$counts
  meta <- tb$meta
  fo   <- ~ g + x + c2
  base <- linda(x, fo, meta, groups = 'g', pseudocount = 0.5, cpus = 1L)
  for (cpus in c(2L, 3L, 4L, 8L))
    expect_identical(linda(x, fo, meta, groups = 'g', pseudocount = 0.5, cpus = cpus), base, info = cpus)

  # A larger table, so that the threads really share the work.
  set.seed(1)
  big  <- matrix(rnbinom(60 * 700, mu = 8, size = 0.5), 60, 700)
  bmet <- data.frame(g = rep(c('a', 'b'), 30), z = rnorm(60))
  rb   <- linda(big, ~ g + z, bmet, pseudocount = 0.5, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(linda(big, ~ g + z, bmet, pseudocount = 0.5, cpus = cpus), rb, info = cpus)
  rbi <- linda(big, ~ g + z, bmet, zero_handling = 'imputation', cpus = 1L)
  expect_identical(linda(big, ~ g + z, bmet, zero_handling = 'imputation', cpus = 4L), rbi)

  fn <- function (counts, ...) linda(counts, fo, meta, groups = 'g', pseudocount = 0.5, cpus = 2L, ...)
  tx <- t(x)
  expect_identical(fn(tx, margin = 2L), base)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    expect_identical(fn(Matrix::Matrix(x, sparse = TRUE)), base)
    expect_identical(fn(Matrix::Matrix(tx, sparse = TRUE), margin = 2L), base)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(fn(slam::as.simple_triplet_matrix(x)), base)
  storage.mode(x) <- 'integer'
  expect_identical(fn(x), base)

  expect_identical(
    da_test(x, 'linda', 'g', meta, fo, pseudocount = 0.5, cpus = 1L), base )
  expect_identical(da_test(x, 'LinDA', 'g', meta, pseudocount = 0.5, cpus = 1L),
                   linda(x, groups = 'g', data = meta, pseudocount = 0.5, cpus = 1L))
  expect_identical(linda(x, groups = meta$g, pseudocount = 0.5, cpus = 1L)$log2_fc,
                   linda(x, ~ g, meta, pseudocount = 0.5, cpus = 1L)$log2_fc)
  rs <- linda(x, 'g + x + c2', meta, groups = 'g', pseudocount = 0.5, cpus = 1L)
  expect_identical(format(attr(rs, 'formula')), format(fo))
  attr(rs, 'formula') <- fo
  expect_identical(rs, base)
  expect_identical(match_da_test('linda')$id, 'linda')
  expect_identical(match_da_test('LinDA')$effect, 'log2_fc')
  expect_identical(match_da_test('linear')$id, 'da_lm')


  # Mixed models: threads, input types, wrappers, formula forms.
  tb   <- ref$tables$subj
  xs   <- tb$counts
  fm   <- ~ g + x + (1 | subj)
  rm1  <- linda(xs, fm, tb$meta, groups = 'g', pseudocount = 0.5, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(linda(xs, fm, tb$meta, groups = 'g', pseudocount = 0.5, cpus = cpus), rm1, info = cpus)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(linda(Matrix::Matrix(t(xs), sparse = TRUE), fm, tb$meta, groups = 'g',
                           pseudocount = 0.5, margin = 2L, cpus = 2L), rm1)
  expect_identical(da_test(xs, 'linda', 'g', tb$meta, fm, pseudocount = 0.5, cpus = 1L), rm1)
  rs <- linda(xs, 'g + x + (1|subj)', tb$meta, groups = 'g', pseudocount = 0.5, cpus = 1L)
  expect_identical(format(attr(rs, 'formula')), '~g + x + (1 | subj)')
  expect_identical(rs$log2_fc, rm1$log2_fc)
  expect_identical(attr(rm1, 'random'), 'subj')
  expect_identical(names(attr(rm1, 'random_sd')), colnames(xs))
  expect_true(all(attr(rm1, 'random_sd') >= 0))
  expect_identical(dim(attr(rm1, 'vcov')), c(ncol(xs), 3L, 3L))
  expect_true(all(rm1$dof > 0 & rm1$dof <= nrow(xs) - 3 + 1e-8))
  # The term of interest is a fixed effect.
  expect_error(linda(xs, fm, tb$meta, groups = 'subj'), 'term in `formula`')
  # Unsupported random effects, and bad groupings.
  expect_error(linda(xs, ~ g + (x | subj), tb$meta), 'random intercept')
  expect_error(linda(xs, ~ g + (1 | subj) + (1 | c2), tb$meta), 'random intercept')
  expect_error(linda(xs, ~ g + (1 | subj:g), tb$meta), 'random intercept')
  expect_error(linda(xs, ~ (1 | subj), tb$meta), 'fixed effect')
  expect_error(linda(xs, ~ g + (1 | nope), tb$meta), 'nope')
  mm <- tb$meta; mm$subj[2] <- NA
  expect_error(linda(xs, fm, mm), 'NA')
  mm <- tb$meta; mm$subj <- 'one'
  expect_error(linda(xs, fm, mm), 'at least two')
  mm <- tb$meta; mm$subj <- seq_len(nrow(mm))
  expect_error(linda(xs, fm, mm), 'fewer values')
  mm <- tb$meta; mm$x2 <- 2 * mm$x
  expect_error(suppressWarnings(linda(xs, ~ g + x + x2 + (1 | subj), mm, standardize = FALSE)), 'full rank')


  # -----------------------------------------------------------------------
  # Invariances
  # -----------------------------------------------------------------------

  # Proportions without zeros, without winsorization: multiplying a sample
  # by a constant changes nothing (the clr cancels it).
  set.seed(2)
  p  <- matrix(rgamma(20 * 15, 2), 20, 15) / 7
  pm <- data.frame(g = rep(c('a', 'b'), 10), z = rnorm(20))
  a  <- linda(p, ~ g + z, pm, outlier_pct = 0, cpus = 1L)
  p2 <- p; p2[3, ] <- p2[3, ] * 13.7; p2[8, ] <- p2[8, ] / 3.3
  b  <- linda(p2, ~ g + z, pm, outlier_pct = 0, cpus = 1L)
  expect_identical(attr(a, 'zero_handling'), 'half_min')
  expect_true(maxrel(a$log2_fc, b$log2_fc) < 1e-12)
  expect_true(maxrel(a$se, b$se) < 1e-12)

  # standardize = FALSE: coefficients per unit, the same up to the scale.
  a <- linda(p, ~ z, pm, outlier_pct = 0, cpus = 1L)
  b <- linda(p, ~ z, pm, outlier_pct = 0, standardize = FALSE, cpus = 1L)
  expect_true(maxrel(a$log2_fc, b$log2_fc * sd(pm$z)) < 1e-10)
  expect_true(maxrel(a$p_value, b$p_value) < 1e-10)
  expect_true(abs(mean(attr(a, 'design')[, 'z'])) < 1e-12 && abs(sd(attr(a, 'design')[, 'z']) - 1) < 1e-12)
  expect_identical(unname(attr(b, 'design')[, 'z']), pm$z)


  # -----------------------------------------------------------------------
  # Zeros
  # -----------------------------------------------------------------------

  tb <- ref$tables$two
  x  <- tb$counts
  # The shared rule: half the smallest non-zero count, with a warning.
  expect_warning(r0 <- linda(x, ~ g, tb$meta, zero_handling = 'pseudocount', cpus = 1L), 'pseudocount = 0.5')
  expect_identical(attr(r0, 'pseudocount'), 0.5)
  expect_identical(r0, linda(x, ~ g, tb$meta, zero_handling = 'pseudocount', pseudocount = 0.5, cpus = 1L))
  expect_error(linda(x, ~ g, tb$meta, zero_handling = 'pseudocount', pseudocount = 0), 'pseudocount = 0')
  # No zeros: nothing added.
  r1 <- linda(x + 1L, ~ g, tb$meta, zero_handling = 'pseudocount', cpus = 1L)
  expect_identical(attr(r1, 'pseudocount'), 0)
  expect_identical(r1, linda(x + 1L, ~ g, tb$meta, zero_handling = 'pseudocount', pseudocount = 0, cpus = 1L))
  # Imputation adds nothing.
  ri <- linda(x, ~ g, tb$meta, zero_handling = 'imputation', cpus = 1L)
  expect_identical(attr(ri, 'pseudocount'), 0)
  expect_identical(attr(ri, 'zero_handling'), 'imputation')
  # A feature of zeros (counts): fitted.
  xz <- x; xz[, 5] <- 0L
  expect_true(all(is.finite(linda(xz, ~ g, tb$meta, zero_handling = 'imputation', cpus = 1L)$log2_fc)))
  expect_true(all(is.finite(linda(xz, ~ g, tb$meta, pseudocount = 0.5, cpus = 1L)$log2_fc)))
  # ... but not with proportions.
  pz <- p; pz[, 4] <- 0; pz[1, 1] <- 0
  expect_error(linda(pz, ~ g, pm), 'no non-zero value')
  # Every sample needs counts.
  x0 <- x; x0[2, ] <- 0L
  expect_error(linda(x0, ~ g, tb$meta), 'positive total')


  # -----------------------------------------------------------------------
  # Edge cases
  # -----------------------------------------------------------------------

  # Two features (the reference's `small` case covers the numbers).
  r <- linda(ref$tables$small$counts, ~ g, ref$tables$small$meta, pseudocount = 0.5, cpus = 1L)
  expect_identical(nrow(r), 2L)

  # All coefficients equal: each is its own mode.
  xe <- outer(c(10, 20, 15, 40, 30, 25), c(1, 2, 3, 4, 5))
  r  <- linda(xe, groups = rep(c('a', 'b'), 3), outlier_pct = 0, cpus = 1L)
  expect_true(all(abs(r$log2_fc) < 1e-12))

  # A rank-deficient design: a warning, and NaN for the aliased column.
  md <- data.frame(g = rep(c('a', 'b'), 10), z = rnorm(20))
  md$z2 <- 2 * md$z
  expect_warning(r <- linda(p, ~ g + z + z2, md, outlier_pct = 0, standardize = FALSE, cpus = 1L), 'rank deficient')
  expect_true(all(is.nan(r$log2_fc[r$covariate == 'z2'])))
  expect_true(all(is.finite(r$log2_fc[r$covariate == 'g' | r$covariate == 'gb'])))
  expect_identical(unique(r$dof), 17)

  expect_error(linda(p[1:3, ], ~ g + z, pm[1:3, ]), 'more samples than design columns')
  expect_error(linda(p, ~ g, data.frame(g = rep('a', 20))), 'single value')
  expect_error(linda(p, ~ z, data.frame(z = rep(1, 20))), 'standardized')
  mna <- pm; mna$z[3] <- NA
  expect_error(linda(p, ~ z, mna), 'NA')
  expect_error(linda(p, ~ log(z), pm), 'standardize = FALSE')
  expect_error(linda(p, ~ 1, pm), 'besides the intercept')
  expect_error(linda(p[, 1, drop = FALSE], ~ g, pm), 'two features')
  expect_error(linda(-p, ~ g, pm), 'negative')
  expect_error(linda(p, ~ g, pm, zero_handling = 'none'))
  expect_error(linda(p, ~ g, pm, zero_handling = NA), 'zero_handling')
  expect_error(linda(p, ~ g, pm, corr_cut = 1), 'corr_cut')
  expect_error(linda(p, ~ g, pm, outlier_pct = -0.1), 'outlier_pct')
  expect_error(linda(p, ~ g, pm, standardize = NA), 'standardize')
  expect_error(linda(p, ~ g, pm, pseudocount = -1), 'pseudocount')
  expect_error(linda(p, ~ g, pm, alpha = 2), 'alpha')
  expect_error(linda(p, ~ g, pm, p_adjust = 'nope'))
  expect_error(linda(p), 'formula')


  # -----------------------------------------------------------------------
  # Attributes, print and plot
  # -----------------------------------------------------------------------

  tb <- ref$tables$conf
  r  <- linda(tb$counts, ~ g + c1 + c2, tb$meta, groups = 'c2', pseudocount = 0.5, cpus = 1L)
  expect_identical(class(r), c('ecodive_linda', 'data.frame'))
  expect_identical(names(r), c('feature', 'covariate', 'log2_fc', 'se', 't', 'dof',
                               'p_value', 'p_adjusted', 'reject'))
  expect_identical(unique(r$covariate), c('gtrt', 'c1', 'c2north', 'c2south'))
  expect_identical(attr(r, 'method'), 'linda')
  expect_identical(attr(r, 'log_base'), 2)
  expect_identical(attr(r, 'alpha'), 0.05)
  expect_identical(attr(r, 'p_adjust'), 'holm')
  expect_identical(attr(r, 'groups'), 'c2')
  expect_identical(attr(r, 'standardize'), TRUE)
  expect_identical(names(attr(r, 'bias')), c('gtrt', 'c1', 'c2north', 'c2south'))
  expect_identical(dim(attr(r, 'vcov')), c(ncol(tb$counts), 5L, 5L))
  expect_identical(names(attr(r, 'mean')), colnames(tb$counts))
  expect_true(all(abs(attr(r, 'mean')) < 50))
  expect_null(attr(r, 'display'))
  expect_true(any(grepl('LinDA', utils::capture.output(print(r)))))

  expect_identical(ecodive:::da_plot_title(r), 'Volcano plot of LinDA')
  grDevices::pdf(NULL)
  expect_identical(withVisible(plot(r))$visible, FALSE)
  expect_silent(plot(r))
  expect_silent(plot(r, 'c1', labels = 0, main = 'x'))
  expect_error(plot(r, 'nope'), '`covariate`')
  rg <- linda(tb$counts, ~ g + c1, tb$meta, pseudocount = 0.5, cpus = 1L)
  for (ty in c('effect', 'ma', 'forest')) expect_silent(plot(rg, type = ty))
  expect_silent(plot(rg, 'c1', type = 'forest', n_features = 5, conf_level = 0.9))
  expect_identical(ecodive:::da_plot_title(rg, 'forest'), 'Forest plot of LinDA')
  expect_error(plot(rg, 'c1', type = 'effect'), 'two-level')
  expect_error(plot(rg, type = 'forest', conf_level = 1), 'conf_level')
  expect_error(plot(rg, type = 'nope'))
  grDevices::dev.off()

#})
