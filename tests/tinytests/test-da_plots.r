#test_that("da_plots", {

  # The per-feature quantities behind the effect, MA and forest plots of DA
  # results: dirmult_ttest()'s pooled_sd, mean and effect (against an R
  # reference from each draw's values, and ALDEx3 1.4.0's aldex.effect() on
  # the same draws, data-raw/expected/da_plots.r), ANCOM-BC's resid_sd and
  # mean attributes (against lm() and normalize()), da_wilcoxon()'s mean;
  # then the plots themselves.

  maxrel <- function (a, b, floor = 1e-300) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.nan(a), is.nan(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }

  set.seed(1)
  x <- matrix(
    rpois(18 * 40, rep(rlnorm(40, 2, 1.2), each = 18)), 18, 40,
    dimnames = list(paste0('s', 1:18), paste0('f', 1:40)) )
  x[1:6, 1:5] <- x[1:6, 1:5] * 6L
  grp <- rep(c('a', 'b', 'c'), each = 6)


  # -----------------------------------------------------------------------
  # dirmult_ttest(): an R reference from each draw's values
  # -----------------------------------------------------------------------

  for (ct in list(c('a', 'b'), c('b', NA), c('c', 'a'))) {
    for (sc in c(0, 0.5)) {
      d <- dirmult_ttest_details(x, grp, contrast = ct, scale_sd = sc, n_draws = 25L, cpus = 1L)
      r <- d$result
      info <- paste(c(ct, sc), collapse = ' ')

      # Each draw's pooled sd and mean, from the clr of the draws.
      expect_true(maxrel(r$pooled_sd, colMeans(d$draws$pooled_sd) / log(2)) < 1e-13, info = info)
      expect_true(maxrel(r$mean,      colMeans(d$draws$mean) / log(2))      < 1e-13, info = info)
      expect_true(maxrel(r$effect,    colMeans(d$draws$diff / d$draws$pooled_sd)) < 1e-13, info = info)
    }
  }

  # Each draw's pooled sd and mean are those of the clr values of given
  # Gamma draws, computed in R.
  G  <- array(stats::rgamma(18 * 40 * 3, shape = 2), c(18, 40, 3))
  d  <- dirmult_ttest_details(x, grp, contrast = c('a', 'c'), n_draws = 3L, gamma_draws = G, cpus = 1L)
  s1 <- grp == 'a'
  s2 <- grp == 'c'
  for (k in 1:3) {
    z  <- log(G[, , k])
    z  <- z - rowMeans(z)
    q  <- colSums(scale(z[s1, ], scale = FALSE)^2) + colSums(scale(z[s2, ], scale = FALSE)^2)
    expect_true(maxrel(d$draws$pooled_sd[k, ], sqrt(q / (sum(s1) + sum(s2) - 2))) < 1e-13, info = k)
    expect_true(maxrel(d$draws$mean[k, ], colMeans(z[s1 | s2, ])) < 1e-13, info = k)
  }

  # The scale model moves effect, but not pooled_sd or mean.
  r0 <- dirmult_ttest(x, grp, scale_sd = 0,   n_draws = 20L, cpus = 1L)
  r5 <- dirmult_ttest(x, grp, scale_sd = 0.5, n_draws = 20L, cpus = 1L)
  expect_identical(r0$pooled_sd, r5$pooled_sd)
  expect_identical(r0$mean, r5$mean)
  expect_false(identical(r0$effect, r5$effect))

  # Identical for any number of threads, and for every input type.
  r1 <- dirmult_ttest(x, grp, contrast = c('b', NA), scale_sd = 0.5, n_draws = 20L, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(dirmult_ttest(x, grp, contrast = c('b', NA), scale_sd = 0.5, n_draws = 20L, cpus = cpus), r1, info = cpus)
  sp <- Matrix::Matrix(x, sparse = TRUE)
  expect_identical(dirmult_ttest(sp, grp, contrast = c('b', NA), scale_sd = 0.5, n_draws = 20L, cpus = 2L), r1)
  expect_identical(dirmult_ttest(t(x), grp, contrast = c('b', NA), margin = 2L, scale_sd = 0.5, n_draws = 20L, cpus = 2L), r1)

  # One sample a side: no pooled sd.
  r <- suppressWarnings(dirmult_ttest(x[1:2, ], c('a', 'b'), n_draws = 4L, cpus = 1L))
  expect_true(all(is.nan(r$pooled_sd)) && all(is.nan(r$effect)))
  expect_true(all(is.finite(r$mean)))
  # One sample on one side: the other side's variance alone.
  r <- suppressWarnings(dirmult_ttest(x[1:4, ], c('a', 'b', 'b', 'b'), n_draws = 4L, cpus = 1L))
  expect_true(all(is.finite(r$pooled_sd)) && all(is.finite(r$effect)))


  # -----------------------------------------------------------------------
  # ALDEx3 1.4.0's aldex.effect() and ALDEx2 1.34.0's aldex.effect()
  # -----------------------------------------------------------------------

  ref <- readRDS(file.path('expected', 'da_plots.rds'))

  # ALDEx3, on its own Dirichlet draws (given as Gamma draws: proportions),
  # with and without a scale model whose shifts are ecodive's.
  for (cs in names(ref$aldex3)) {
    e <- ref$aldex3[[cs]]
    d <- dirmult_ttest_details(
      ref$counts, ref$groups, contrast = c('b', 'a'), scale_sd = e$scale_sd,
      n_draws = dim(e$gamma)[[3]], gamma_draws = e$gamma, seed = e$seed, cpus = 1L )
    r <- d$result
    expect_true(maxrel(r$log2_fc,   e$estimate)  < 1e-10, info = cs)
    expect_true(maxrel(r$pooled_sd, e$pooled_sd) < 1e-10, info = cs)
    expect_true(maxrel(r$effect,    e$cohens_d)  < 1e-10, info = cs)
  }

  # ALDEx3 with its own draws, on a deep table: Monte Carlo error only.
  r <- dirmult_ttest(ref$deep, ref$groups, contrast = c('b', 'a'), n_draws = 256L, cpus = 2L)
  e <- ref$aldex3_deep
  expect_true(maxrel(r$pooled_sd, e$pooled_sd) < 0.01)
  expect_true(max(abs(r$effect - e$cohens_d) / pmax(abs(e$cohens_d), 1)) < 0.02)

  # ALDEx2's effect statistics, from random pairs and medians: the same
  # ordering, and close in size (outputs only).
  e <- ref$aldex2
  expect_true(stats::cor(r$mean, e$rab.all, method = 'spearman') > 0.99)
  expect_true(stats::cor(r$log2_fc, e$diff.btw) > 0.99)
  expect_true(stats::cor(r$pooled_sd, e$diff.win, method = 'spearman') > 0.9)
  expect_true(stats::cor(r$effect, e$effect) > 0.99)


  # -----------------------------------------------------------------------
  # ANCOM-BC's resid_sd and mean
  # -----------------------------------------------------------------------

  meta <- data.frame(grp = grp, age = c(30, 41, 52, 33, 45, 50, 38, 29, 61, 44, 35, 48, 57, 40, 31, 46, 53, 36))
  X    <- stats::model.matrix(~ grp + age, meta)

  # ANCOM-BC: the initial fit's residuals, about each sample's mean
  # residual (theta). theta is orthogonal to the design, so lm() of the
  # logs less theta has the same residuals.
  f  <- ancombc(x, ~ grp + age, meta, pseudocount = 1, cpus = 1L)
  L  <- log(x + 1)
  th <- rowMeans(L - X %*% qr.coef(qr(X), L))
  sg <- apply(L - th, 2, function (y) stats::sigma(stats::lm(y ~ X - 1)))
  expect_true(maxrel(attr(f, 'resid_sd'), sg) < 1e-10)
  expect_identical(names(attr(f, 'resid_sd')), colnames(x))
  expect_true(maxrel(attr(f, 'mean'), colMeans(normalize(x, 'clr', pseudocount = 1))) < 1e-12)
  expect_identical(names(attr(f, 'mean')), colnames(x))

  # ANCOM-BC2: the final fit, of the centred logs less each sample's
  # sampling fraction.
  f2 <- ancombc2(x, ~ grp + age, meta, pseudocount = 1, sensitivity = FALSE, cpus = 1L)
  L2 <- scale(L, scale = FALSE) - attr(f2, 'sample_fractions')
  sg <- apply(L2, 2, function (y) stats::sigma(stats::lm(y ~ X - 1)))
  expect_true(maxrel(attr(f2, 'resid_sd'), sg) < 1e-10)
  expect_true(maxrel(attr(f2, 'mean'), attr(f, 'mean')) < 1e-12)

  # For one two-level factor alone, the pooled within-group sd.
  ab <- grp != 'c'
  f1 <- ancombc(x[ab, ], groups = grp[ab], pseudocount = 1, cpus = 1L)
  L1 <- log(x[ab, ] + 1)
  r1 <- L1 - rowMeans(L1 - stats::model.matrix(~ grp[ab]) %*% qr.coef(qr(stats::model.matrix(~ grp[ab])), L1))
  z  <- apply(r1, 2, function (y) {
    a <- y[grp[ab] == 'a']; b <- y[grp[ab] == 'b']
    sqrt((sum((a - mean(a))^2) + sum((b - mean(b))^2)) / (length(y) - 2))
  })
  expect_true(maxrel(attr(f1, 'resid_sd'), z) < 1e-10)

  # Zeros as missing values (ANCOM-BC2, pseudocount 0): the rclr's mean,
  # each feature over its observed samples; resid_sd over them too.
  xz <- x
  xz[sample(length(xz), 60)] <- 0L
  xz[xz[, 1] == 0, 1] <- 1L
  fz <- ancombc2(xz, ~ grp, meta, pseudocount = 0, sensitivity = FALSE, cpus = 1L)
  lz <- log(xz)
  lz[xz == 0] <- NA
  rc <- lz - rowMeans(lz, na.rm = TRUE)
  expect_true(maxrel(attr(fz, 'mean'), colMeans(rc, na.rm = TRUE)) < 1e-12)
  rs <- attr(fz, 'resid_sd')
  expect_true(all(is.finite(rs) | is.nan(rs)))
  expect_true(any(is.finite(rs)))

  # The same for any number of threads; the sensitivity fits carry them.
  big <- cbind(x, x[, 40:1] + 1L, x + 2L)
  colnames(big) <- paste0('f', seq_len(ncol(big)))
  g1 <- ancombc2(big, ~ grp + age, meta, sensitivity = FALSE, cpus = 1L)
  for (cpus in c(2L, 4L)) {
    g <- ancombc2(big, ~ grp + age, meta, sensitivity = FALSE, cpus = cpus)
    expect_identical(attr(g, 'resid_sd'), attr(g1, 'resid_sd'), info = cpus)
    expect_identical(attr(g, 'mean'), attr(g1, 'mean'), info = cpus)
  }
  fs <- ancombc2(x, ~ grp, meta, cpus = 1L)
  for (r in attr(fs, 'sensitivity')$fits)
    expect_true(!is.null(attr(r, 'resid_sd')) && !is.null(attr(r, 'mean')))


  # -----------------------------------------------------------------------
  # da_wilcoxon()'s mean: the same as da_ttest()'s
  # -----------------------------------------------------------------------

  w <- da_wilcoxon(x, grp, contrast = c('a', 'c'), pseudocount = 1, cpus = 1L)
  tt <- da_ttest(x, grp, contrast = c('a', 'c'), pseudocount = 1, cpus = 1L)
  expect_identical(names(w), c('feature', 'w', 'auc', 'median_diff', 'mean', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(w$mean, tt$mean)
  cl <- normalize(x, 'clr', pseudocount = 1)
  expect_true(maxrel(w$mean, colMeans(cl[grp != 'b', ])) < 1e-13)
  wr <- da_wilcoxon(x, grp, norm = 'rclr', cpus = 1L)
  tr <- da_ttest(x, grp, norm = 'rclr', cpus = 1L)
  expect_identical(wr$mean, tr$mean)
  for (cpus in c(2L, 4L))
    expect_identical(da_wilcoxon(x, grp, contrast = c('a', 'c'), pseudocount = 1, cpus = cpus), w, info = cpus)


  # -----------------------------------------------------------------------
  # The plots
  # -----------------------------------------------------------------------

  grDevices::pdf(NULL)

  invisibly <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }

  dt <- dirmult_ttest(x, grp, contrast = c('a', 'b'), n_draws = 16L, cpus = 1L)
  rl <- da_lm(x, ~ grp + age, meta, groups = 'grp', pseudocount = 1, cpus = 1L)
  res <- list(dirmult_ttest = dt, ancombc = f, ancombc2 = fs, da_ttest = tt, da_lm = rl)
  for (nm in names(res)) {
    r <- res[[nm]]
    for (ty in c('volcano', 'effect', 'ma', 'forest')) {
      expect_true(invisibly(r, type = ty), info = paste(nm, ty))
      expect_true(invisibly(r, type = ty, labels = 0, main = 'x', sub = 'y'), info = paste(nm, ty))
      expect_identical(
        da_plot_title(r, ty),
        sprintf(
          '%s of %s',
          c(volcano = 'Volcano plot', effect = 'Effect plot', ma = 'MA plot', forest = 'Forest plot')[[ty]],
          da_method_name(nm) ),
        info = paste(nm, ty) )
    }
    expect_true(invisibly(r, type = 'forest', n_features = 3), info = nm)
    expect_true(invisibly(r, type = 'forest', n_features = Inf), info = nm)
    expect_error(plot(r, type = 'forest', n_features = 0), '`n_features` must be at least 1', info = nm)
    expect_error(plot(r, type = 'nope'), '`type` must be', info = nm)
  }
  expect_identical(da_plot_title(fs, 'effect'), 'Effect plot of ANCOM-BC2')
  expect_identical(da_plot_title(dt), 'Volcano plot of Dirichlet-Multinomial t-test')

  # Wilcoxon: volcano and MA plots only.
  for (ty in c('volcano', 'ma')) expect_true(invisibly(w, type = ty), info = ty)
  expect_error(plot(w, type = 'effect'), '`type` must be "volcano" or "ma"')

  # Model results: the covariate stays second; the effect plot needs a
  # 0/1 design column.
  expect_true(invisibly(f, 'age', 'ma'))
  expect_true(invisibly(f, 'grpb', type = 'forest', conf_level = 0.9))
  expect_true(invisibly(rl, 'age', 'forest'))
  expect_error(plot(f, 'age', 'effect'), 'needs a two-level comparison.*`age`')
  expect_error(plot(rl, 'age', type = 'effect'), 'needs a two-level comparison')
  expect_error(plot(f, '(Intercept)', type = 'effect'), 'needs a two-level comparison')
  expect_error(plot(f, type = 'forest', conf_level = 1), '`conf_level` must be')
  expect_true(invisibly(fz, type = 'effect'))     # NaN dispersions left out

  # Axis labels follow log_base (and the norm, for the plain tests).
  expect_identical(lfc_label(dt), 'Log2 fold change')
  expect_identical(lfc_label(f), 'Log fold change')
  expect_identical(da_effect_label(tt, 'Difference in means'), 'Log fold change (difference in means, clr)')
  expect_identical(da_norm_units(tt), 'clr, natural log')
  tp <- da_ttest(x, grp, norm = 'percent', cpus = 1L)
  expect_identical(da_effect_label(tp, 'Difference in means'), 'Difference in means (percent)')
  expect_identical(da_norm_units(tp), 'percent')
  expect_true(invisibly(tp, type = 'effect'))
  expect_identical(da_plot_contrast(dt), 'a (n = 6) vs b (n = 6)')
  expect_identical(da_plot_contrast(da_ttest(x, grp, pseudocount = 1, cpus = 1L)), 'a (n = 6) vs all others (n = 12)')
  expect_identical(plot_pct(0.95), '95%')
  expect_identical(
    da_plot_sub('grpb', c(TRUE, FALSE, NA), robust = c(TRUE, FALSE, FALSE), extra = 'top 2 shown'),
    'grpb: 1 of 3 features differentially abundant, 1 robust to the pseudocount (filled); top 2 shown' )

  # A result made before these columns existed.
  old <- dt
  old$pooled_sd <- NULL
  expect_error(plot(old, 'effect'), 'no `pooled_sd` column')
  expect_true(invisibly(old, 'ma'))
  old <- f
  attr(old, 'mean') <- NULL
  expect_error(plot(old, type = 'ma'), 'no `mean` attribute')

  # The forest plot restores the margins it widened.
  mai <- graphics::par('mai')
  plot(dt, 'forest')
  expect_identical(graphics::par('mai'), mai)

  grDevices::dev.off()

#})
