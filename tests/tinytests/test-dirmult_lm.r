#test_that("dirmult_lm", {

  # dirmult_lm(): ALDEx3 1.4.0's aldex() on the same Dirichlet and scale
  # draws (linear models with classical, HC0 and HC3 standard errors, three
  # p.adjust methods; random intercepts with lme4 / lmerTest); ALDEx3 and
  # ALDEx2 with their own draws, by distribution
  # (data-raw/expected/dirmult_lm.r); the scale models; agreement with
  # dirmult_ttest(); threading, input types, da_test(), argument checks,
  # print and plot.

  ref  <- readRDS(file.path('expected', 'dirmult_lm.rds'))
  x    <- ref$counts
  meta <- ref$meta
  G    <- ref$gamma_draws
  S    <- ref$n_draws
  n    <- nrow(x)
  D    <- ncol(x)

  # Largest difference relative to each value, or to `floor` (a fraction
  # of the largest value) for values near 0.
  maxrel <- function (a, b, floor = 1e-3) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.nan(a), is.nan(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) return (0)
    f  <- floor * max(abs(b[ok]))
    max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), f))
  }

  # The result's column as a P x D matrix (ALDEx3's layout).
  pd <- function (res, k) matrix(res[[k]], nrow = length(unique(res$covariate)))

  run <- function (cs, ...) {
    args <- list(
      x, cs$formula, meta, n_draws = S, se_type = cs$se_type,
      p_adjust = cs$p_adjust, gamma_draws = G, cpus = 1L, ... )
    if (cs$scale == 'given')     args$scale_draws <- cs$ref$logScale
    else if (cs$scale == 'load') { args$load <- 'load'; args$load_sd <- 0.4 }
    else                         { args$scale_model <- cs$scale; args$scale_sd <- 0.5 }
    do.call(dirmult_lm_details, args)
  }


  # -----------------------------------------------------------------------
  # ALDEx3, on the same draws
  # -----------------------------------------------------------------------

  for (nm in names(ref$cases)) {
    cs <- ref$cases[[nm]]
    if (grepl('lme4', nm)) next
    e <- run(cs)
    r <- e$result
    expect_true(maxrel(pd(r, 'estimate'),   cs$ref$estimate)  < 1e-10, info = nm)
    expect_true(maxrel(pd(r, 'se'),         cs$ref$std.error) < 1e-10, info = nm)
    expect_true(maxrel(pd(r, 'p_value'),    cs$ref$p.val)     < 1e-10, info = nm)
    expect_true(maxrel(pd(r, 'p_adjusted'), cs$ref$p.val.adj) < 1e-10, info = nm)
    expect_identical(r$reject, as.vector(cs$ref$p.val.adj) <= 0.05, info = nm)
    expect_true(all(r$dof == n - ncol(attr(r, 'design'))), info = nm)

    if (!is.null(cs$ref$draw_estimate)) {
      expect_true(maxrel(e$draws$estimate, cs$ref$draw_estimate)  < 1e-10, info = nm)
      expect_true(maxrel(e$draws$se,       cs$ref$draw_std.error) < 1e-10, info = nm)
    }
  }

  # A random intercept: lme4 / lmerTest stop short of the exact REML
  # optimum, so agreement is to their optimizer's tolerance.
  for (nm in grep('lme4', names(ref$cases), value = TRUE)) {
    cs <- ref$cases[[nm]]
    e  <- run(cs)
    r  <- e$result
    expect_identical(attr(r, 'random'), 'id', info = nm)
    expect_true(maxrel(e$draws$estimate, cs$ref$draw_estimate)  < 1e-4, info = nm)
    expect_true(maxrel(e$draws$se,       cs$ref$draw_std.error) < 1e-4, info = nm)
    expect_true(maxrel(e$draws$df[, , 1:4], cs$ref$draw_df)     < 1e-2, info = nm)
    expect_true(maxrel(pd(r, 'estimate'),   cs$ref$estimate)  < 1e-4, info = nm)
    expect_true(maxrel(pd(r, 'se'),         cs$ref$std.error) < 1e-4, info = nm)
    expect_true(maxrel(pd(r, 'p_value'),    cs$ref$p.val)     < 1e-3, info = nm)
    expect_true(maxrel(pd(r, 'p_adjusted'), cs$ref$p.val.adj) < 1e-3, info = nm)
    expect_identical(r$reject, as.vector(cs$ref$p.val.adj) <= 0.05, info = nm)
  }


  # lm() with many design columns (more than the kernel keeps on the stack),
  # classical and HC3 standard errors, on given draws without a scale.
  set.seed(2)
  nb <- 40L; Db <- 70L
  xb <- matrix(rpois(nb * Db, 30), nb, Db)
  mb <- data.frame(f = factor(rep(letters[1:18], length.out = nb)), z = rnorm(nb))
  Gb <- array(rgamma(nb * Db * 2, 3), c(nb, Db, 2))
  Xb <- model.matrix(~ f + z, mb)
  Yb <- log2(Gb[, , 2]); Yb <- Yb - rowMeans(Yb)
  fb <- stats::lm.fit(Xb, Yb)
  rb <- fb$residuals
  s2 <- colSums(rb^2) / (nb - ncol(Xb))
  Ab <- solve(crossprod(Xb), t(Xb))
  hb <- rowSums(Xb * t(Ab))
  for (se in c('classical', 'hc0', 'hc3')) {
    e <- dirmult_lm_details(xb, ~ f + z, mb, n_draws = 2L, gamma_draws = Gb, scale_sd = 0, se_type = se, cpus = 2L)
    want <- switch(se,
      classical = sqrt(outer(diag(solve(crossprod(Xb))), s2)),
      hc0       = sqrt(Ab^2 %*% rb^2),
      hc3       = sqrt(Ab^2 %*% (rb^2 / (1 - hb)^2)) )
    expect_true(maxrel(e$draws$estimate[, , 2], fb$coefficients) < 1e-10, info = se)
    expect_true(maxrel(e$draws$se[, , 2], want) < 1e-10, info = se)
  }


  # -----------------------------------------------------------------------
  # The scale models' draws
  # -----------------------------------------------------------------------

  lG  <- log(G)
  lse <- apply(lG, c(1, 3), function (v) { m <- max(v); m + log(sum(exp(v - m))) })
  lc  <- sweep(lG, c(1, 3), lse) / log(2)      # n x D x S log2 compositions
  X2  <- model.matrix(~ grp + age, meta)

  e <- dirmult_lm_details(x, ~ grp + age, meta, n_draws = S, gamma_draws = G, cpus = 1L)
  expect_identical(dim(e$draws$shift), c(ncol(X2), S))
  expect_true(max(abs(e$draws$scale - (-apply(lc, c(1, 3), mean) + X2 %*% e$draws$shift))) < 1e-12)

  e <- dirmult_lm_details(x, ~ grp + age, meta, n_draws = S, gamma_draws = G, scale_model = 'tss', cpus = 1L)
  expect_true(max(abs(e$draws$scale - X2 %*% e$draws$shift)) < 1e-12)

  e <- dirmult_lm_details(x, ~ grp, meta, n_draws = S, gamma_draws = G, load = 'load', load_sd = 0.4, cpus = 1L)
  expect_identical(dim(e$draws$shift), c(n, S))
  expect_true(max(abs(e$draws$scale - (log2(meta$load) + e$draws$shift))) < 1e-12)
  expect_identical(attr(e$result, 'scale_model'), 'load')
  expect_identical(attr(e$result, 'load_sd'), rep(0.4, n))
  expect_identical(attr(e$result, 'scale_sd'), NA_real_)

  # The shifts: N(0, scale_sd^2), the same for any table (one stream per
  # draw, after the samples'), and none with scale_sd = 0.
  sh <- dirmult_lm_details(x, ~ grp + age, meta, n_draws = 4000L, scale_sd = 0.7, cpus = 1L)$draws$shift
  expect_true(suppressWarnings(stats::ks.test(as.vector(sh), 'pnorm', 0, 0.7)$p.value) > 1e-3)
  expect_true(abs(stats::cor(sh[2, ], sh[3, ])) < 0.1)
  sh2 <- dirmult_lm_details(x + 5L, ~ grp + age, meta, n_draws = 10L, scale_sd = 0.7, cpus = 1L)$draws$shift
  expect_identical(sh2, sh[, 1:10])
  expect_true(all(dirmult_lm_details(x, ~ grp, meta, n_draws = 4L, scale_sd = 0, cpus = 1L)$draws$shift == 0))

  # Draw d doesn't depend on n_draws.
  a <- dirmult_lm_details(x, ~ grp, meta, n_draws = 6L, cpus = 1L)$draws
  b <- dirmult_lm_details(x, ~ grp, meta, n_draws = 3L, cpus = 1L)$draws
  expect_identical(a$estimate[, , 1:3], b$estimate)


  # -----------------------------------------------------------------------
  # ALDEx3 and ALDEx2 with their own draws, by distribution
  # -----------------------------------------------------------------------

  own <- ref$aldex3_own
  e   <- dirmult_lm_details(ref$deep, ~ grp + age, meta, n_draws = 2000L, p_adjust = 'BH')
  r   <- e$result
  sdv <- apply(e$draws$estimate, c(1, 2), stats::sd)
  # Mean estimates: Monte Carlo error of the scale model's shifts (sd 0.5 /
  # sqrt(2000) each, per design column).
  expect_true(max(abs(pd(r, 'estimate') - own$estimate)) < 0.05)
  expect_true(max(abs(log(sdv / own$estimate_sd))) < 0.1)
  expect_true(max(abs(pd(r, 'se') / own$std.error - 1)) < 0.02)
  expect_true(stats::cor(-log10(r$p_value), -log10(as.vector(own$p.val))) > 0.95)

  if (!is.null(ref$aldex2_glm)) {
    g  <- ref$aldex2_glm
    r0 <- dirmult_lm(ref$deep, ~ grp + age, meta, scale_sd = 0, n_draws = 500L, se_type = 'classical')
    for (cv in c('grpb', 'grpc', 'age')) {
      i   <- r0$covariate == cv
      est <- r0$estimate[i]
      expect_true(stats::cor(est, g[, paste0(cv, ':Est')]) > 0.9999, info = cv)
      expect_true(max(abs(est - g[, paste0(cv, ':Est')])) < 0.01, info = cv)
      expect_true(max(abs(r0$se[i] / g[, paste0(cv, ':SE')] - 1)) < 0.01, info = cv)
      expect_true(stats::cor(-log10(r0$p_value[i]), -log10(g[, paste0(cv, ':pval')])) > 0.999, info = cv)
    }
  }


  # -----------------------------------------------------------------------
  # dirmult_ttest(): the same draws, the same differences
  # -----------------------------------------------------------------------

  for (seed in 0:1) {
    lm <- dirmult_lm(x, ~ grp2, meta, scale_sd = 0, n_draws = 24L, seed = seed)
    tt <- dirmult_ttest(x, 'grp2', meta, contrast = c('trt', 'ctl'), n_draws = 24L, seed = seed)
    expect_true(max(abs(lm$estimate[lm$covariate == 'grp2trt'] - tt$log2_fc)) < 1e-12)
    expect_true(max(abs(attr(lm, 'resid_sd') - tt$pooled_sd)) < 1e-12)
    expect_true(max(abs(attr(lm, 'mean') - tt$mean)) < 1e-12)
  }

  # Classical standard errors and no scale model: each draw's test is the
  # pooled-variance t-test of the clr values (log2).
  e  <- dirmult_lm_details(x, ~ grp2, meta, scale_sd = 0, n_draws = 3L, se_type = 'classical')
  tt <- dirmult_ttest_details(x, 'grp2', meta, contrast = c('trt', 'ctl'), n_draws = 3L)$draws
  n1 <- sum(meta$grp2 == 'trt'); n2 <- n - n1
  se_pooled <- tt$pooled_sd * sqrt(1 / n1 + 1 / n2) / log(2)
  expect_true(maxrel(e$draws$se[2, , ], t(se_pooled)) < 1e-12)


  # -----------------------------------------------------------------------
  # Threads, input types, da_test(), seeds
  # -----------------------------------------------------------------------

  r1 <- dirmult_lm(x, ~ grp + age, meta, n_draws = 40L, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(dirmult_lm(x, ~ grp + age, meta, n_draws = 40L, cpus = cpus), r1, info = cpus)
  expect_identical(dirmult_lm(t(x), ~ grp + age, meta, margin = 2L, n_draws = 40L, cpus = 1L), r1)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(dirmult_lm(Matrix::Matrix(x, sparse = TRUE), ~ grp + age, meta, n_draws = 40L, cpus = 1L), r1)
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(dirmult_lm(slam::as.simple_triplet_matrix(x), ~ grp + age, meta, n_draws = 40L, cpus = 1L), r1)
  xd <- x; storage.mode(xd) <- 'double'
  expect_identical(dirmult_lm(xd, ~ grp + age, meta, n_draws = 40L, cpus = 1L), r1)
  expect_identical(da_test(x, 'dirmult_lm', formula = ~ grp + age, data = meta, n_draws = 40L, cpus = 1L), r1)
  expect_identical(da_test(x, 'ALDEx3', 'grp', meta, n_draws = 8L), dirmult_lm(x, groups = 'grp', data = meta, n_draws = 8L))
  expect_false(identical(dirmult_lm(x, ~ grp + age, meta, n_draws = 40L, seed = 1L, cpus = 1L), r1))

  m1 <- dirmult_lm(x, ~ grp2 + (1 | id), meta, n_draws = 8L, cpus = 1L)
  expect_identical(dirmult_lm(x, ~ grp2 + (1 | id), meta, n_draws = 8L, cpus = 4L), m1)
  expect_identical(attr(m1, 'se_type'), 'classical')
  expect_identical(attr(m1, 'formula'), ~ grp2 + (1 | id), info = 'formula keeps (1 | id)')

  # groups alone; values as groups.
  g1 <- dirmult_lm(x, groups = 'grp2', data = meta, n_draws = 8L)
  g0 <- dirmult_lm(x, ~ grp2, meta, n_draws = 8L)
  expect_identical(attr(g1, 'groups'), 'grp2')
  expect_identical(lapply(g1, identity), lapply(g0, identity), info = 'groups')
  expect_identical(attr(g1, 'resid_sd'), attr(g0, 'resid_sd'))
  expect_equal(attr(g1, 'formula'), ~ grp2, ignore_attr = TRUE)
  g2 <- dirmult_lm(x, groups = meta$grp2, n_draws = 8L)
  expect_identical(g2$estimate, g1$estimate)


  # -----------------------------------------------------------------------
  # Results' shape and attributes
  # -----------------------------------------------------------------------

  expect_identical(names(r1), c('feature', 'covariate', 'estimate', 'se', 'dof', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(nrow(r1), D * 4L)
  expect_identical(unique(r1$covariate), c('(Intercept)', 'grpb', 'grpc', 'age'))
  expect_identical(r1$feature[1:5], rep(colnames(x)[1:2], c(4, 1)))
  expect_identical(class(r1), c('ecodive_dirmult_lm', 'data.frame'))
  expect_identical(attr(r1, 'method'), 'dirmult_lm')
  expect_identical(attr(r1, 'scale_model'), 'clr')
  expect_identical(attr(r1, 'scale_sd'), 0.5)
  expect_identical(attr(r1, 'se_type'), 'hc3')
  expect_identical(attr(r1, 'n_draws'), 40L)
  expect_identical(attr(r1, 'p_adjust'), 'holm')
  expect_identical(attr(r1, 'log_base'), 2)
  expect_identical(names(attr(r1, 'mean')), colnames(x))
  expect_identical(names(attr(r1, 'resid_sd')), colnames(x))
  expect_identical(list_da_tests('list')$dirmult_lm$effect, 'estimate')
  expect_identical(match_da_test('ALDEx3')$id, 'dirmult_lm')
  expect_identical(match_da_test('aldex 3')$id, 'dirmult_lm')
  expect_identical(match_da_test('ALDEx2')$id, 'dirmult_ttest')
  expect_identical(match_da_test('aldex')$id, 'dirmult_ttest')
  expect_identical(match_da_test('Dirichlet-Multinomial Linear Model')$id, 'dirmult_lm')

  # p_adjust = 'none': the adjusted p-values are the p-values.
  r0 <- dirmult_lm(x, ~ grp, meta, n_draws = 8L, p_adjust = 'none')
  expect_identical(r0$p_adjusted, r0$p_value)

  # An intercept-free model and an intercept-only one.
  expect_identical(unique(dirmult_lm(x, ~ 0 + grp, meta, n_draws = 4L)$covariate), c('grpa', 'grpb', 'grpc'))
  expect_identical(unique(dirmult_lm(x, ~ 1, meta, n_draws = 4L)$covariate), '(Intercept)')

  # Leverage 1 (a group of one sample) with HC3: NaN, with a warning.
  m2 <- meta; m2$grp[1] <- 'z'
  expect_warning(r <- dirmult_lm(x, ~ grp, m2, n_draws = 4L), 'leverage 1')
  expect_true(all(is.nan(r$se)) && !any(r$reject))
  expect_silent(r <- dirmult_lm(x, ~ grp, m2, n_draws = 4L, se_type = 'classical'))
  expect_true(all(is.finite(r$se)))


  # -----------------------------------------------------------------------
  # Errors
  # -----------------------------------------------------------------------

  expect_error(dirmult_lm(x, ~ grp, meta, pseudocount = 0), '`pseudocount`')
  expect_error(dirmult_lm(x, ~ grp, meta, n_draws = 0), '`n_draws`')
  expect_error(dirmult_lm(x, ~ grp, meta, alpha = 1), '`alpha`')
  expect_error(dirmult_lm(x, ~ grp, meta, se_type = 'hc9'), '`se_type`')
  expect_error(dirmult_lm(x, ~ grp, meta, scale_model = 'alr'), '`scale_model`')
  expect_error(dirmult_lm(x, ~ grp, meta, scale_sd = -1), '`scale_sd`')
  expect_error(dirmult_lm(x, ~ grp, meta, p_adjust = 'nope'))
  expect_error(dirmult_lm(x[, 1, drop = FALSE], ~ grp, meta), 'at least 2 features')
  expect_error(dirmult_lm(-x, ~ grp, meta), 'negative')
  expect_error(dirmult_lm(x[1:3, ], ~ grp + age, meta[1:3, ]), 'more samples than design columns')
  m3 <- meta; m3$age2 <- 2 * m3$age
  expect_error(dirmult_lm(x, ~ age + age2, m3), 'rank deficient')
  expect_error(dirmult_lm(x, ~ grp, meta, load = 'load'), '`load_sd` is required')
  expect_error(dirmult_lm(x, ~ grp, meta, load = 'load', load_sd = 0.3, scale_sd = 1), 'only without `load`')
  expect_error(dirmult_lm(x, ~ grp, meta, load = 'load', load_sd = 0.3, scale_model = 'tss'), 'only without `load`')
  expect_error(dirmult_lm(x, ~ grp, meta, load_sd = 0.3), 'only for `load`')
  expect_error(dirmult_lm(x, ~ grp, meta, load = c(-1, meta$load[-1]), load_sd = 0.3), 'load')
  expect_error(dirmult_lm(x, ~ grp, meta, load = c(NA, meta$load[-1]), load_sd = 0.3), 'load')
  expect_error(dirmult_lm(x, ~ grp, meta, load = 'load', load_sd = -1), '`load_sd`')
  expect_error(dirmult_lm(x, ~ grp, meta, load = 'load', load_sd = c(1, 2)), '`load_sd`')
  expect_error(dirmult_lm(x, ~ grp + (1 | id), meta, se_type = 'hc3'), 'random intercept')
  expect_error(dirmult_lm(x, ~ grp + (age | id), meta), 'random intercept')
  m4 <- meta; m4$one <- 'a'
  expect_error(dirmult_lm(x, ~ grp + (1 | one), m4), 'at least two values')
  expect_error(dirmult_lm(x, ~ grp), '`data`')

  # load_sd per sample, or a column of `data`.
  m5 <- meta; m5$lsd <- seq(0.1, 0.5, length.out = n)
  a <- dirmult_lm(x, ~ grp, m5, load = 'load', load_sd = 'lsd', n_draws = 4L)
  b <- dirmult_lm(x, ~ grp, m5, load = m5$load, load_sd = m5$lsd, n_draws = 4L)
  expect_identical(a, b)
  expect_identical(attr(a, 'load_sd'), m5$lsd)


  # -----------------------------------------------------------------------
  # Print and plot
  # -----------------------------------------------------------------------

  expect_identical(ecodive:::da_plot_title(r1), 'Volcano plot of Dirichlet-Multinomial Linear Model')
  expect_identical(ecodive:::dirmult_lm_scale_label(r1), 'scale model: clr, scale_sd = 0.5')
  expect_identical(ecodive:::dirmult_lm_scale_label(a), 'scale model: measured load')
  expect_true(any(grepl('Dirichlet-Multinomial Linear Model', capture.output(print(r1)))))
  pdf(NULL)
  expect_identical(withVisible(plot(r1))$visible, FALSE)
  expect_silent(plot(r1))
  expect_silent(plot(r1, 'age', labels = 0, main = 'x'))
  for (ty in c('effect', 'ma', 'forest', 'rank')) expect_silent(plot(r1, type = ty))
  expect_silent(plot(r1, type = 'forest', n_features = 5, conf_level = 0.9))
  expect_error(plot(r1, 'age', type = 'effect'), 'two-level')
  expect_error(plot(r1, 'nope'), '`covariate`')
  expect_error(plot(r1, type = 'forest', conf_level = 1), 'conf_level')
  expect_silent(plot(m1, type = 'forest'))
  dev.off()

#})
