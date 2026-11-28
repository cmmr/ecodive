#test_that("fastancom", {

  # fastancom(): fastANCOM 0.0.4's fastANCOM() (data-raw/expected/fastancom.r)
  # on two groups, a numeric variable, a covariate, zeros, pseudocounts,
  # proportions, few features, reference-set ties and sizes, and 600 random
  # tables for the decision rules; an independent R computation of the
  # log-ratio tests and effects, for more than two groups; ancom()'s W;
  # threading, input types, the wrapper, invariances, edge cases and argument
  # checks.

  ref <- readRDS(file.path('expected', 'fastancom.rds'))

  # Relative difference, with an absolute floor for values near 0.
  maxrel <- function (a, b, floor = 1e-8) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.na(a), is.na(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    if (!identical(a[!ok & !is.na(a)], b[!ok & !is.na(b)])) return (Inf)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }


  # -----------------------------------------------------------------------
  # fastANCOM
  # -----------------------------------------------------------------------

  for (cs in names(ref$cases)) {

    r <- ref$cases[[cs]]
    a <- r$args
    v <- names(r$meta)[[1]]
    f <- if (cs == 'covariate') ~ grp + z else stats::as.formula(paste('~', v))

    e <- fastancom(
      r$counts, f, r$meta, groups = v, pseudocount = a$pseudo, alpha = a$sig,
      detect_rate = a$detect, ref_rate = a$ref, p_adjust = 'BH', cpus = 1L )

    expect_identical(e$feature, colnames(r$counts), info = cs)
    expect_identical(e$w, r$W, info = cs)
    expect_identical(e$reject, r$REJECT, info = cs)
    expect_true(maxrel(attr(e, 'coef')[, 2], r$beta) < 1e-10, info = cs)
    expect_true(maxrel(e$log2_fc, r$log2FC) < 1e-10, info = cs)
    expect_true(maxrel(e$se, r$SD) < 1e-10, info = cs)
    expect_true(maxrel(e$statistic, r$log2FC / r$SD) < 1e-10, info = cs)

    dof <- nrow(r$counts) - ncol(attr(e, 'design'))
    expect_identical(attr(e, 'dof'), dof, info = cs)

    if (cs != 'covariate') {
      expect_true(maxrel(e$p_value, r$pval) < 1e-10, info = cs)
      expect_true(maxrel(e$p_adjusted, r$qval) < 1e-10, info = cs)
    }
    else {
      # fastANCOM tests on n - 2 degrees of freedom whatever the covariates;
      # fastancom() on n - p: the same statistics, on dof.
      expect_identical(dof, 27L)
      p <- 2 * stats::pt(-abs(r$log2FC / r$SD), dof)
      expect_true(maxrel(e$p_value, p) < 1e-10)
      expect_true(maxrel(e$p_adjusted, stats::p.adjust(p, 'BH')) < 1e-10)
      expect_true(maxrel(r$pval, 2 * stats::pt(-abs(r$log2FC / r$SD), nrow(r$counts) - 2)) < 1e-10)
      W <- rowSums(2 * stats::pt(-abs(r$stat), dof) < 0.05, na.rm = TRUE)
      expect_identical(e$w, as.integer(W))
    }
  }

  # The decision rules, on fastANCOM's W of 600 random tables.
  for (i in seq_along(ref$decisions)) {
    d   <- ref$decisions[[i]]
    up  <- ecodive:::fastancom_detect(d$W, d$K)
    det <- if (is.na(up)) logical(d$K) else d$W >= up
    expect_identical(det, d$detected, info = i)
    dec <- ecodive:::fastancom_reject(d$W, d$K, d$detect_rate)
    expect_identical(dec$reject, d$reject, info = i)
    expect_identical(dec$reject, d$W >= dec$cutoff, info = i)
  }

  # Below 10 features, or below max(W) / K = 0.26, nothing is detected.
  expect_true(is.na(ecodive:::fastancom_detect(c(8L, 8L, 0:6), 9L)))
  expect_true(is.na(ecodive:::fastancom_detect(c(2L, rep(0L, 9)), 10L)))
  expect_false(is.na(ecodive:::fastancom_detect(c(9L, 9L, rep(0L, 8)), 10L)))


  # -----------------------------------------------------------------------
  # An independent computation: each log ratio's own linear model.
  # -----------------------------------------------------------------------

  set.seed(3)
  n <- 18; K <- 25
  x <- matrix(rpois(n * K, 20), n, K, dimnames = list(NULL, sprintf('f%02d', 1:K)))
  x[sample(length(x), 30)] <- 0
  meta <- data.frame(
    status = rep(c('lo', 'mid', 'hi'), each = 6),
    age    = round(rnorm(n, 40, 10)) )
  x[meta$status == 'mid', 1:3] <- x[meta$status == 'mid', 1:3] * 4
  x[meta$status == 'hi', 4:6]  <- x[meta$status == 'hi', 4:6] * 4

  e <- fastancom(x, ~ status + age, meta, groups = 'status', pseudocount = 0.5,
                 ref_rate = 0.12, cpus = 1L)

  cmps <- c('status', 'statuslo', 'statusmid', 'statusmid - statuslo')
  expect_identical(unique(e$covariate), cmps)
  expect_identical(e$covariate, rep(cmps, K))
  expect_identical(e$feature, rep(colnames(x), each = 4))

  L   <- log2(x + 0.5)
  mf  <- data.frame(meta, status = factor(meta$status, c('hi', 'lo', 'mid')))
  X   <- stats::model.matrix(~ status + age, mf)
  dof <- n - ncol(X)
  fit <- function (y) stats::lm.fit(X, y)
  con <- list(
    statuslo = c(0, 1, 0, 0), statusmid = c(0, 0, 1, 0),
    'statusmid - statuslo' = c(0, -1, 1, 0) )

  for (cm in cmps) {
    r <- e[e$covariate == cm, ]
    W <- integer(K)
    for (k in 1:(K - 1)) for (j in (k + 1):K) {
      f1 <- fit(L[, k] - L[, j])
      s2 <- sum(f1$residuals^2) / dof
      if (cm == 'status') {
        f0 <- stats::lm.fit(X[, -(2:3)], L[, k] - L[, j])
        Fs <- ((sum(f0$residuals^2) - sum(f1$residuals^2)) / 2) / s2
        p  <- stats::pf(Fs, 2, dof, lower.tail = FALSE)
      }
      else {
        V  <- s2 * chol2inv(qr.R(qr(X)))
        lam <- con[[cm]]
        t  <- sum(lam * f1$coefficients) / sqrt(drop(lam %*% V %*% lam))
        p  <- 2 * stats::pt(-abs(t), dof)
      }
      if (p < 0.05) { W[k] <- W[k] + 1L; W[j] <- W[j] + 1L }
    }
    expect_identical(r$w, W, info = cm)

    # The reference set: the floor(0.12 * 25) = 3 smallest W, ties in order.
    R <- order(W)[1:3]
    expect_identical(attr(e, 'reference')[[cm]], colnames(x)[R], info = cm)
    lr <- L - rowMeans(L[, R])
    if (cm == 'status') {
      Fk <- vapply(1:K, function (k) {
        f1 <- fit(lr[, k]); f0 <- stats::lm.fit(X[, -(2:3)], lr[, k])
        ((sum(f0$residuals^2) - sum(f1$residuals^2)) / 2) / (sum(f1$residuals^2) / dof)
      }, numeric(1))
      expect_true(all(is.nan(r$log2_fc)) && all(is.nan(r$se)), info = cm)
      expect_true(maxrel(r$statistic, Fk) < 1e-10, info = cm)
      expect_true(maxrel(r$p_value, stats::pf(Fk, 2, dof, lower.tail = FALSE)) < 1e-10, info = cm)
    }
    else {
      lam <- con[[cm]]
      est <- se <- numeric(K)
      for (k in 1:K) {
        f1 <- fit(lr[, k])
        V  <- sum(f1$residuals^2) / dof * chol2inv(qr.R(qr(X)))
        est[k] <- sum(lam * f1$coefficients)
        se[k]  <- sqrt(drop(lam %*% V %*% lam))
      }
      expect_true(maxrel(r$log2_fc, est) < 1e-10, info = cm)
      expect_true(maxrel(r$se, se) < 1e-10, info = cm)
      expect_true(maxrel(r$p_value, 2 * stats::pt(-abs(est / se), dof)) < 1e-10, info = cm)
    }
    expect_identical(r$p_adjusted, stats::p.adjust(r$p_value, 'holm'), info = cm)
    expect_identical(r$reject, ecodive:::fastancom_reject(W, K, 0.7)$reject, info = cm)
  }

  # W is ancom()'s with sig_test = 'anova' and no adjustment.
  an <- ancom(x, meta$status, pseudocount = 0.5, sig_test = 'anova', p_adjust = 'none')
  e3 <- fastancom(x, groups = meta$status, pseudocount = 0.5, cpus = 1L)
  expect_identical(e3$w[e3$covariate == 'group'], an$W)
  two <- rep(c('a', 'b'), 9)
  expect_identical(
    fastancom(x, groups = two, pseudocount = 0.5)$w,
    ancom(x, two, pseudocount = 0.5, sig_test = 'anova', p_adjust = 'none')$W )


  # -----------------------------------------------------------------------
  # Threads, input types, the wrapper, invariances
  # -----------------------------------------------------------------------

  r  <- ref$cases$manual
  fn <- function (counts, ...) fastancom(counts, groups = 'group', data = r$meta, pseudocount = 0.5, ...)
  base <- fn(r$counts, cpus = 1L)
  for (cp in 2:4) expect_identical(fn(r$counts, cpus = cp), base, info = cp)
  expect_identical(fn(t(r$counts), margin = 2L), base)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    expect_identical(fn(Matrix::Matrix(r$counts, sparse = TRUE)), base)
    expect_identical(fn(Matrix::Matrix(t(r$counts), sparse = TRUE), margin = 2L), base)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(fn(slam::as.simple_triplet_matrix(r$counts)), base)

  e3b <- fastancom(x, ~ status + age, meta, groups = 'status', pseudocount = 0.5, ref_rate = 0.12, cpus = 4L)
  expect_identical(e3b, e)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(
      fastancom(Matrix::Matrix(x, sparse = TRUE), ~ status + age, meta, groups = 'status',
                pseudocount = 0.5, ref_rate = 0.12, cpus = 3L), e )

  expect_identical(da_test(r$counts, 'fastancom', 'group', r$meta, pseudocount = 0.5), base)
  expect_identical(da_test(r$counts, 'fastANCOM', 'group', r$meta, pseudocount = 0.5), base)

  # A zero-free table without a pseudocount: scaling a sample cancels out.
  z  <- ref$cases$nopseudo
  e0 <- fastancom(z$counts, groups = 'grp', data = z$meta, pseudocount = 0)
  expect_identical(fastancom(z$counts, groups = 'grp', data = z$meta), e0)
  zs <- z$counts
  zs[3, ] <- zs[3, ] * 7.5
  zs[9, ] <- zs[9, ] / 3
  e1 <- fastancom(zs, groups = 'grp', data = z$meta, pseudocount = 0)
  expect_identical(e1$w, e0$w)
  expect_identical(e1$reject, e0$reject)
  expect_true(maxrel(e1$log2_fc, e0$log2_fc) < 1e-12)
  expect_true(maxrel(e1$se, e0$se) < 1e-12)

  # The shared pseudocount rule.
  zz <- ref$cases$zeros
  half <- min(zz$counts[zz$counts > 0]) / 2
  expect_warning(ez <- fastancom(zz$counts, groups = 'grp', data = zz$meta), 'Zeros detected')
  expect_identical(attr(ez, 'pseudocount'), half)
  expect_identical(ez, fastancom(zz$counts, groups = 'grp', data = zz$meta, pseudocount = half))


  # -----------------------------------------------------------------------
  # Edge cases
  # -----------------------------------------------------------------------

  g2 <- rep(c('a', 'b'), each = 5)
  set.seed(8)
  y  <- matrix(rpois(10 * 12, 30) + 1, 10, 12)

  # Two features.
  e2 <- fastancom(y[, 1:2], groups = g2)
  expect_identical(nrow(e2), 2L)
  expect_true(all(e2$w %in% 0:1) && !any(e2$reject))

  # One reference feature: its own log2 fold change is 0, with no p-value.
  e1r <- fastancom(y, groups = g2, ref_rate = 0.01)
  i   <- match(attr(e1r, 'reference')[[1]], e1r$feature)
  expect_identical(length(i), 1L)
  expect_identical(e1r$log2_fc[[i]], 0)
  expect_identical(e1r$se[[i]], 0)
  expect_true(is.nan(e1r$p_value[[i]]) && is.nan(e1r$p_adjusted[[i]]))

  # All W tied (every feature the same): no tests, no rejections.
  same <- matrix(rep(y[, 1], 12), 10, 12)
  es   <- fastancom(same, groups = g2)
  expect_identical(es$w, integer(12))
  expect_false(any(es$reject))
  expect_true(all(es$log2_fc == 0))

  # A constant feature: exactly no residual, so its log ratios are those of
  # the other feature.
  yc <- y
  yc[, 1] <- 20
  ec <- fastancom(yc, groups = g2)
  expect_identical(unname(attr(ec, 'coef')[1, ]), c(log2(20), 0))
  expect_true(all(is.finite(ec$se[-1])))
  yc[, 2] <- 40
  ec2 <- fastancom(yc, groups = g2)
  expect_true(ec2$w[[1]] <= ec$w[[1]])
  expect_identical(attr(ec2, 'coef')[2, 2], 0)

  # A numeric variable of interest, and character / factor / logical groups.
  en <- fastancom(y, groups = seq(1, 10))
  expect_identical(unique(en$covariate), 'group')
  ef <- fastancom(y, groups = factor(g2, c('b', 'a')))
  expect_identical(unique(ef$covariate), 'groupa')
  el <- fastancom(y, groups = g2 == 'b')
  expect_identical(unique(el$covariate), 'groupTRUE')
  expect_identical(el$w, fastancom(y, groups = g2)$w)

  # A term of several columns that isn't a factor: the global test only.
  ep <- fastancom(y, ~ poly(t, 2), data.frame(t = 1:10))
  expect_identical(unique(ep$covariate), 'poly(t, 2)')
  expect_true(all(is.nan(ep$log2_fc)))

  # Errors.
  m2 <- data.frame(g = g2, a = 1:10, b = 2 * (1:10))
  expect_error(fastancom(y, ~ g + a, m2), '`groups` is required')
  expect_error(fastancom(y, ~ g - 1, m2), 'intercept')
  expect_error(fastancom(y, ~ g + a + b, m2, groups = 'g'), 'rank deficient')
  expect_error(fastancom(y[c(1, 2, 6), ], ~ g + a, m2[c(1, 2, 6), ], groups = 'g'), 'more samples')
  expect_error(fastancom(y, ~ g, data.frame(g = c(NA, g2[-1]))), 'Missing values')
  expect_error(fastancom(y, ~ g + a, m2, groups = 'b'), 'term in `formula`')
  y0 <- y; y0[1, 1] <- 0
  expect_error(fastancom(y0, groups = g2, pseudocount = 0), 'pseudocount = 0')
  expect_error(fastancom(y, groups = g2, pseudocount = -1), 'pseudocount')
  expect_error(fastancom(y, groups = g2, pseudocount = NA), 'pseudocount')
  expect_error(fastancom(y, groups = g2, alpha = 1), '`alpha`')
  expect_error(fastancom(y, groups = g2, detect_rate = 0), '`detect_rate`')
  expect_error(fastancom(y, groups = g2, detect_rate = c(0.5, 0.6)), '`detect_rate`')
  expect_error(fastancom(y, groups = g2, ref_rate = 0.6), '`ref_rate`')
  expect_error(fastancom(y, groups = g2, ref_rate = 0), '`ref_rate`')
  expect_error(fastancom(y, groups = g2, p_adjust = 'nope'))
  expect_error(fastancom(y[, 1, drop = FALSE], groups = g2), 'at least 2 features')
  expect_error(fastancom(y), 'groups')
  expect_error(fastancom(y, groups = rep('a', 10)), 'single value')


  # -----------------------------------------------------------------------
  # Attributes, print and plot
  # -----------------------------------------------------------------------

  expect_identical(class(e), c('ecodive_fastancom', 'data.frame'))
  expect_identical(names(e), c('feature', 'covariate', 'log2_fc', 'se', 'statistic',
                               'w', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(attr(e, 'method'), 'fastancom')
  expect_identical(attr(e, 'groups'), 'status')
  expect_identical(attr(e, 'log_base'), 2)
  expect_identical(attr(e, 'alpha'), 0.05)
  expect_identical(attr(e, 'detect_rate'), 0.7)
  expect_identical(attr(e, 'ref_rate'), 0.12)
  expect_identical(attr(e, 'p_adjust'), 'holm')
  expect_identical(attr(e, 'pseudocount'), 0.5)
  expect_identical(attr(e, 'dof'), 14L)
  expect_identical(names(attr(e, 'cutoff')), cmps)
  expect_identical(names(attr(e, 'reference')), cmps)
  expect_identical(dim(attr(e, 'coef')), c(25L, 4L))
  expect_identical(names(attr(e, 'mean')), colnames(x))
  expect_true(maxrel(attr(e, 'mean'), colMeans(L)) < 1e-12)
  expect_null(attr(e, 'display'))
  expect_true(any(grepl('fastANCOM', utils::capture.output(print(e)))))
  expect_identical(match_da_test('fastancom')$name, 'fastANCOM')

  expect_identical(ecodive:::da_plot_title(base), 'Volcano plot of fastANCOM')
  grDevices::pdf(NULL)
  expect_identical(withVisible(plot(base))$visible, FALSE)
  expect_silent(plot(base))
  for (ty in c('volcano', 'ma', 'forest', 'rank')) expect_silent(plot(base, type = ty))
  expect_silent(plot(base, type = 'forest', n_features = 5, conf_level = 0.9))
  expect_silent(plot(e))
  expect_silent(plot(e, 'statusmid - statuslo', type = 'volcano', labels = 0, main = 'x'))
  expect_error(plot(e, type = 'volcano'), 'global test')
  expect_error(plot(e, 'nope'), '`covariate`')
  expect_error(plot(base, type = 'forest', conf_level = 1), 'conf_level')
  expect_error(plot(base, type = 'nope'))
  grDevices::dev.off()

#})
