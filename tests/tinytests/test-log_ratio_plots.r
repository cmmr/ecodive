#test_that("log_ratio_plots", {

  # log_ratio_features() against Qurro 0.9.0's autoselection, and the log
  # ratios of its features against qarcoal() and the sample plot's
  # computeBalance() (data-raw/expected/log_ratio_plots.r, .py, .js); the
  # rankings that each result defines; the ecodive_log_ratio class; and the
  # rank and log-ratio plots, drawn on a null device.

  ref <- readRDS(file.path('expected', 'log_ratio_plots.rds'))

  plain <- function (v) { attributes(v) <- NULL; v }


  # -----------------------------------------------------------------------
  # Qurro's autoselection and log ratios
  # -----------------------------------------------------------------------

  max_q <- 0
  max_b <- 0
  for (cs in ref$cases) {

    # Qurro's rank plot keeps every ranked feature, in the ranking's order.
    expect_identical(cs$qurro_features, names(cs$ranking), info = cs$id)

    # The same features, in the same order: Qurro returns the top n in
    # ascending order, ecodive the largest first.
    sel <- log_ratio_features(cs$ranking, cs$n)
    expect_identical(sel$numerator, rev(cs$top), info = cs$id)
    expect_identical(sel$denominator, cs$bot,    info = cs$id)

    # Their log ratio: NA exactly where Qurro drops a sample.
    lr <- suppressWarnings(log_ratio(ref$tables[[cs$table]], sel$numerator, sel$denominator))
    expect_identical(names(lr), names(cs$qarcoal), info = cs$id)
    expect_identical(is.na(plain(lr)), is.na(plain(cs$qarcoal)), info = cs$id)
    expect_identical(is.na(plain(lr)), is.na(plain(cs$balance)), info = cs$id)
    max_q <- max(max_q, abs(plain(lr) - plain(cs$qarcoal)), na.rm = TRUE)
    max_b <- max(max_b, abs(plain(lr) - plain(cs$balance)), na.rm = TRUE)
  }
  # qarcoal()'s log(N / D) differs only where numpy's log() and C's differ
  # by an ulp; computeBalance()'s log(N) - log(D) rounds differently.
  expect_true(max_q < 1e-16)
  expect_true(max_b < 1e-14)


  # The rankings that ecodive's results define are those Qurro was given.
  x    <- ref$tables$da
  meta <- ref$meta
  rk   <- function (id) ref$cases[[id]]$ranking
  same <- function (fit, id, ...) {
    s <- rank_source(fit, ...)
    e <- s$effect[is.finite(s$effect)]
    f <- s$feature[is.finite(s$effect)]
    expect_identical(f, names(rk(id)), info = id)
    expect_true(max(abs(e - rk(id)) / pmax(abs(rk(id)), 1e-12)) < 1e-10, info = id)
    n <- ref$cases[[id]]$n
    expect_identical(
      log_ratio_features(fit, n, ...),
      list(numerator = rev(ref$cases[[id]]$top), denominator = ref$cases[[id]]$bot),
      info = id )
  }

  fit_dm <- dirmult_ttest(x, meta$grp, contrast = c('a', 'b'), n_draws = 64L, seed = 1L, cpus = 1L)
  same(fit_dm, 'dirmult_ttest_n5')
  fit_bc <- ancombc2(x, ~ grp + age, meta, groups = 'grp')
  same(fit_bc, 'ancombc2_grpb_n5')
  same(fit_bc, 'ancombc2_age_n5', covariate = 'age')
  same(da_ttest(x, meta$grp, pseudocount = 0.5), 'da_ttest_n5')
  same(da_wilcoxon(x, meta$grp, pseudocount = 0.5), 'da_wilcoxon_n5')
  fit_lm <- da_lm(x, ~ grp + age, meta, groups = 'grp', pseudocount = 0.5)
  same(fit_lm, 'da_lm_n5')
  fit_ld <- linda(x, ~ grp + age, meta, groups = 'grp')
  same(fit_ld, 'linda_n5')
  ord <- pca(x, norm = 'clr', pseudocount = 0.5)
  same(ord, 'pca_ax1_n5')
  same(ord, 'pca_ax2_n10', axes = 2L)
  same(rpca(x), 'rpca_n5')

  # Each method's own effect column.
  expect_identical(rank_source(fit_dm)$effect, fit_dm$log2_fc)
  expect_identical(rank_source(fit_bc, 'age')$effect, fit_bc$log_fc[fit_bc$covariate == 'age'])
  expect_identical(rank_source(fit_bc)$what, 'grpb')
  expect_identical(rank_source(fit_lm)$what, 'grpb')
  expect_identical(rank_source(ord, axes = 2L)$effect, unname(attr(ord, 'loadings')[, 2]))
  expect_identical(rank_source(ord, axes = 2L)$axis, 'PC2')

  # Threads and sparse input change nothing.
  for (k in c(2L, 4L))
    expect_identical(
      log_ratio_features(dirmult_ttest(x, meta$grp, contrast = c('a', 'b'), n_draws = 64L, seed = 1L, cpus = k), 5),
      log_ratio_features(fit_dm, 5), info = k )
  if (requireNamespace('Matrix', quietly = TRUE)) {
    xs  <- Matrix::Matrix(x, sparse = TRUE)
    sel <- log_ratio_features(fit_dm, 5)
    expect_identical(do.call(log_ratio, c(list(xs), sel)), do.call(log_ratio, c(list(x), sel)))
    expect_identical(
      log_ratio_features(pca(xs, norm = 'clr', pseudocount = 0.5), 5),
      log_ratio_features(ord, 5) )
  }


  # -----------------------------------------------------------------------
  # log_ratio_features(): exact cases and errors
  # -----------------------------------------------------------------------

  # Ties keep the features' order; NA, NaN and Inf are left out.
  e <- c(a = 1, b = 0, c = 1, d = -2, e = NA, f = 0, g = -2, h = Inf)
  expect_identical(log_ratio_features(e, 1), list(numerator = 'c', denominator = 'd'))
  expect_identical(log_ratio_features(e, 2), list(numerator = c('c', 'a'), denominator = c('d', 'g')))
  expect_identical(log_ratio_features(e, 3), list(numerator = c('c', 'a', 'f'), denominator = c('d', 'g', 'b')))
  expect_identical(log_ratio_features(c(3, 1, 2, 4), 1), list(numerator = 4L, denominator = 2L))
  expect_identical(log_ratio_features(e, 2L), log_ratio_features(e, 2.0))

  expect_error(log_ratio_features(e, 4), '`n` must be at most 3: 6 features are ranked')
  expect_error(log_ratio_features(e, 0),   '`n` must be a whole number of at least 1')
  expect_error(log_ratio_features(e, 2.5), '`n` must be a whole number of at least 1')
  expect_error(log_ratio_features(e, NA),  '`n` must be a whole number of at least 1')
  expect_error(log_ratio_features(e, '1'), '`n` must be a whole number of at least 1')
  expect_error(log_ratio_features(e, 1, covariate = 'a'), '`covariate` is for differential')
  expect_error(log_ratio_features(e, 1, axes = 2), '`axes` is for ordinations')
  expect_error(log_ratio_features(fit_dm, 1, covariate = 'a'), 'dirmult_ttest\\(\\) results have none')
  expect_error(log_ratio_features(fit_bc, 1, covariate = 'x'), '`covariate` must be one of')
  expect_error(log_ratio_features(fit_dm, 1, axes = 1), '`axes` is for ordinations')
  expect_error(log_ratio_features(ord, 1, covariate = 'a'), 'give `axes` for an ordination')
  expect_error(log_ratio_features(ord, 1, axes = 9), '`axes` must be an axis number from 1 to')
  expect_error(log_ratio_features(pcoa(bray(x)), 1), 'no feature loadings')
  expect_error(log_ratio_features(ancom(x, meta$grp, pseudocount = 0.5), 1), 'ancom\\(\\) results have no signed effect')
  expect_error(log_ratio_features(lefse(x, meta$grp, cpus = 1L), 1), 'lefse\\(\\) results have no signed effect')
  expect_error(log_ratio_features(iris, 1), '`x` must be the result of')
  expect_error(log_ratio_features('a', 1), '`x` must be the result of')

  # ANCOM-BC's NaN estimates are left out.
  fit_nan <- fit_bc
  fit_nan$log_fc[fit_nan$covariate == 'grpb'][1:3] <- NaN
  s <- log_ratio_features(fit_nan, 18)
  expect_false(any(c('f1', 'f2', 'f3') %in% unlist(s)))
  expect_error(log_ratio_features(fit_nan, 19), '`n` must be at most 18: 37 features are ranked')

  # rpca_phylo(): the tips only, with a message.
  op <- rpca_phylo(ecodive::ex_counts, ecodive::ex_tree, n_axes = 2)
  expect_message(s <- log_ratio_features(op, 2), 'Ranking the 6 tips only')
  expect_true(all(unlist(s) %in% colnames(ecodive::ex_counts)))
  ld <- attr(op, 'loadings')[, 1]
  ld <- ld[names(ld) %in% attr(op, 'tree')$tip.label]
  expect_identical(s, log_ratio_features(ld, 2))
  expect_error(suppressMessages(log_ratio_features(op, 4)), '`n` must be at most 3')
  expect_identical(do.call(log_ratio, c(list(ecodive::ex_counts + 1), s)), log_ratio(ecodive::ex_counts + 1, s$numerator, s$denominator))


  # -----------------------------------------------------------------------
  # The ecodive_log_ratio class
  # -----------------------------------------------------------------------

  lr <- log_ratio(iris[, 1:4], 'Petal.Width', c('Sepal.Length', 'Sepal.Width'))
  bl <- balance(iris[, 1:4], 3:4, 1:2)
  expect_identical(class(lr), 'ecodive_log_ratio')
  expect_identical(class(bl), 'ecodive_log_ratio')
  expect_true(is.numeric(lr))
  expect_identical(
    names(attributes(lr)), c('method', 'numerator', 'denominator', 'pseudocount', 'class'))

  # Subsetting, sorting and arithmetic give plain (named) numbers.
  v <- setNames(plain(lr), names(lr))
  expect_identical(lr[1:2], v[1:2])
  expect_identical(sort(lr), sort(v))
  expect_identical(rev(lr), rev(v))
  expect_identical(-lr, -v)
  expect_identical(2 * lr, 2 * v)
  expect_identical(lr - lr, v - v)
  expect_identical(lr - mean(lr), v - mean(v))
  expect_identical(abs(lr), abs(v))
  expect_identical(round(lr, 2), round(v, 2))
  expect_identical(lr > -2, v > -2)
  expect_identical(lr + bl, v + setNames(plain(bl), names(bl)))
  expect_identical(as.numeric(lr), plain(lr))
  d <- data.frame(id = seq_along(lr))
  d$lr <- lr
  expect_true(is.numeric(d$lr))
  expect_identical(stats::coef(stats::lm(lr ~ Species, iris)), stats::coef(stats::lm(v ~ Species, iris)))

  # print(): a header, then the values, never the attributes.
  out <- utils::capture.output(print(lr[1:3]))
  expect_false(any(grepl('attr\\(', utils::capture.output(print(lr)))))
  out <- utils::capture.output(res <- withVisible(print(lr)))
  expect_identical(out[1:3], c(
    'Log ratio of 1 over 2 features; pseudocount 0',
    "Numerator:   'Petal.Width'",
    "Denominator: 'Sepal.Length', 'Sepal.Width'" ))
  expect_false(res$visible)
  expect_identical(res$value, lr)
  expect_identical(
    utils::capture.output(print(bl))[1],
    'Balance of 2 over 2 features; pseudocount 0' )
  x0 <- rbind(a = c(1, 0, 2), b = c(0, 0, 3))
  colnames(x0) <- c('f1', 'f2', 'f3')
  z  <- suppressWarnings(log_ratio(x0, 1:2, 3))
  expect_identical(utils::capture.output(print(z))[4], '1 sample is NA (a zero sum)')
  expect_identical(utils::capture.output(print(log_ratio(x0 + 1, 1, 2:3)))[1],
    'Log ratio of 1 over 2 features; pseudocount 0')


  # -----------------------------------------------------------------------
  # Rank plots
  # -----------------------------------------------------------------------

  grDevices::pdf(NULL)
  set.seed(1)
  seed <- .Random.seed

  invisibly <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }

  # The layout: ascending, ties in feature order, non-finite left out.
  lay <- rank_layout(e, names(e), reject = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), labels = 10)
  expect_identical(lay$x, 1:6)
  expect_identical(lay$feature, c('d', 'g', 'b', 'f', 'a', 'c'))
  expect_identical(lay$effect, c(-2, -2, 0, 0, 1, 1))
  expect_identical(lay$class, c('Rejected', 'None', 'None', 'None', 'Rejected', 'None'))
  expect_identical(lay$col, c('grey30', 'grey85', 'grey85', 'grey85', 'grey30', 'grey85'))
  expect_identical(lay$labelled, c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(attr(lay, 'omitted'), 2L)
  expect_identical(rank_omitted(lay), '2 features with no estimate not shown')

  # A selection replaces the decisions; labels go to its largest |effect|.
  lay <- rank_layout(e, names(e), reject = rep(TRUE, 8), num = c(1, 3), den = 4, labels = 2)
  expect_identical(lay$class, c('Denominator', 'None', 'None', 'None', 'Numerator', 'Numerator'))
  expect_identical(lay$col, c('royalblue3', 'grey85', 'grey85', 'grey85', 'firebrick', 'firebrick'))
  expect_identical(lay$labelled, c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))

  # No decisions (ordinations): all mid-grey, the largest |loading| named.
  lay <- rank_layout(c(p = 0.1, q = -0.5, r = 0.3), c('p', 'q', 'r'), labels = 1, name_any = TRUE)
  expect_identical(lay$col, rep('grey50', 3))
  expect_identical(lay$labelled, c(TRUE, FALSE, FALSE))
  expect_identical(rank_layout(1:3, 1:3, labels = 5)$labelled, rep(FALSE, 3))

  # Every DA result with a signed effect, and ordinations with loadings.
  sel <- log_ratio_features(fit_dm, 3)
  expect_true(invisibly(fit_dm, 'rank'))
  expect_true(invisibly(fit_dm, type = 'rank', numerator = sel$numerator, denominator = sel$denominator))
  expect_true(invisibly(fit_dm, type = 'rank', labels = 0, main = 'x', ylim = c(-5, 5)))
  expect_true(invisibly(fit_bc, type = 'rank'))
  expect_true(invisibly(fit_bc, 'age', 'rank', numerator = 'f1', denominator = 2:3))
  expect_true(invisibly(fit_nan, type = 'rank'))
  expect_true(invisibly(ancombc(x, ~ grp, meta, pseudocount = 0.5), type = 'rank'))
  expect_true(invisibly(da_ttest(x, meta$grp, pseudocount = 0.5), 'rank'))
  expect_true(invisibly(da_wilcoxon(x, meta$grp, pseudocount = 0.5), 'rank'))
  expect_true(invisibly(fit_lm, type = 'rank'))
  expect_true(invisibly(fit_lm, 'age', type = 'rank', numerator = sel$numerator))
  expect_true(invisibly(fit_ld, type = 'rank', denominator = c(TRUE, rep(FALSE, 39))))
  expect_true(invisibly(ord, type = 'rank'))
  expect_true(invisibly(ord, type = 'rank', axes = 2, numerator = 'f1', denominator = 'f2'))
  expect_true(invisibly(ord, type = 'rank', loadings = 0, legend = NULL, numerator = 'f1'))
  expect_true(invisibly(rpca(x), type = 'rank'))
  expect_true(invisibly(op, type = 'rank', numerator = 'n7'))
  expect_true(invisibly(rda(iris[, 1:4], data = iris, formula = ~ Species), type = 'rank'))

  # Other types of an ordination are still plot.default()'s.
  expect_true(invisibly(ord, type = 'n'))
  expect_true(invisibly(ord, type = 'p', axes = c(2, 1)))
  expect_true(invisibly(ord))

  # Titles and labels.
  expect_identical(da_plot_title(fit_dm, 'rank'), 'Rank plot of Dirichlet-Multinomial t-test')
  expect_identical(da_plot_title(fit_bc, 'rank'), 'Rank plot of ANCOM-BC2')

  # Errors.
  expect_error(plot(pcoa(bray(x)), type = 'rank'), 'no feature loadings')
  expect_error(plot(fit_dm, 'rank', numerator = 'f1', denominator = c('f2', 'f1')),
    "must not share features, but they share 1 feature: 'f1'")
  expect_error(plot(fit_dm, 'rank', numerator = 'zz'), "`numerator` names a feature not in the result: 'zz'")
  expect_error(plot(ord, type = 'rank', denominator = 99), '`denominator` has an index outside 1 to 40 \\(the result has 40 features\\)')
  expect_error(plot(fit_dm, type = 'x'), '`type` must be')
  expect_error(plot(fit_bc, type = 'x'), '`type` must be')
  expect_error(plot(fit_dm, 'volcano', numerator = 'f1'), 'are for `type = "rank"`')
  expect_error(plot(ord, numerator = 'f1'), 'are for `type = "rank"`')
  expect_error(plot(ord, type = 'rank', axes = 9), '`axes` must be an axis number from 1 to')


  # -----------------------------------------------------------------------
  # Log-ratio plots
  # -----------------------------------------------------------------------

  expect_true(invisibly(lr))
  expect_true(invisibly(lr, 'Species', iris))
  expect_true(invisibly(lr, 'Sepal.Length', iris))
  expect_true(invisibly(lr, iris$Species, col = c('red', 'green', 'blue'), pch = 1:3, legend = NULL))
  expect_true(invisibly(lr, iris$Species == 'setosa', main = 'x', ylab = 'y'))
  expect_true(invisibly(bl, iris$Species))
  sp <- factor(iris$Species, levels = c('virginica', 'setosa', 'versicolor', 'other'))
  expect_true(invisibly(lr, sp))

  # Undefined values are left out and counted.
  xz <- rbind(s1 = c(1, 2, 3, 4), s2 = c(0, 0, 3, 4), s3 = c(5, 1, 0, 0), s4 = c(2, 2, 2, 2))
  colnames(xz) <- c('a', 'b', 'c', 'd')
  z  <- suppressWarnings(log_ratio(xz, c('a', 'b'), c('c', 'd')))
  expect_true(invisibly(z, c('x', 'x', 'y', 'y')))
  expect_true(invisibly(z, c(1, 2, 3, 4)))
  expect_identical(
    log_ratio_plot_sub(z),
    "Numerator: 'a', 'b'; denominator: 'c', 'd'; 2 samples with an undefined log ratio not shown" )
  expect_identical(log_ratio_ylab(z), 'log((a + b) / (c + d))')
  expect_identical(log_ratio_ylab(lr), 'log(Petal.Width / (Sepal.Length + Sepal.Width))')
  expect_identical(log_ratio_ylab(bl), 'Balance of Petal.Length + Petal.Width vs Sepal.Length + Sepal.Width')
  zb <- suppressWarnings(balance(xz, 1, 2:4, pseudocount = 0))
  expect_identical(
    log_ratio_plot_sub(zb),
    "Numerator: 'a'; denominator: 'b', 'c', 'd'; 2 samples with an undefined balance not shown" )
  w  <- matrix(1:80, 2, 40, dimnames = list(NULL, paste0('f', 1:40)))
  expect_identical(log_ratio_ylab(log_ratio(w, 1:4, 5:40)), 'log(4 features / 36 features)')
  expect_identical(
    log_ratio_plot_sub(log_ratio(w, 1:4, 5:40)),
    "Numerator: 'f1', 'f2', 'f3', 'f4'; denominator: 'f5', 'f6', 'f7', 'f8', 'f9', ..." )
  expect_error(plot(suppressWarnings(log_ratio(xz[2:3, ], 1:2, 3:4))), 'no defined values')

  # The spread within a box: by value, from the golden ratio, no RNG.
  sp <- log_ratio_spread(c(1, 1, 1, 2), c(3, 1, 2, 5))
  expect_identical(sp$y, c(3, 1, 2, 5))
  expect_equal(sp$x, c(
    1 + 0.6 * ((3 * 0.618034) %% 1 - 0.5), 1 + 0.6 * ((1 * 0.618034) %% 1 - 0.5),
    1 + 0.6 * ((2 * 0.618034) %% 1 - 0.5), 2 + 0.6 * ((1 * 0.618034) %% 1 - 0.5) ))
  expect_true(all(abs(log_ratio_spread(rep(1, 1000), seq_len(1000))$x - 1) <= 0.3))

  # Nothing touched R's random number generator.
  expect_identical(.Random.seed, seed)

  grDevices::dev.off()

#})
