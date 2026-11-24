#test_that("da_test", {

  # da_test(), list_da_tests(), match_da_test(); da_wilcoxon(), da_ttest()
  # and da_lm() against base R's wilcox.test(), t.test() and lm(), feature
  # by feature, on the same normalized values; the OLS kernel against lm()
  # and the sandwich formulas; threading, input types, edge cases and
  # argument checks.

  maxrel <- function (a, b, floor = 1e-300) {
    a <- as.vector(a); b <- as.vector(b)
    if (!identical(is.na(a), is.na(b))) return (Inf)
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), floor)))
  }

  set.seed(1)
  n <- 24L; D <- 30L
  x <- matrix(rpois(n * D, 6), n, D, dimnames = list(paste0('S', 1:n), paste0('F', 1:D)))
  x[sample(length(x), 120)] <- 0L
  x[, 1:4] <- x[, 1:4] + rep(c(0L, 8L), each = 12)
  meta <- data.frame(
    row.names = rownames(x),
    grp  = rep(c('a', 'b'), each = 12),
    g3   = rep(c('lo', 'mid', 'hi'), 8),
    age  = round(rnorm(n, 40, 10)),
    sex  = rep(c('f', 'm', 'm', 'f'), 6),
    stringsAsFactors = FALSE )
  pc <- 0.5

  # The values each test sees: a normalized table, NA where missing (rclr).
  nx <- function (x, norm, pc = 0.5) {
    y <- switch(norm, none = x + 0, normalize(x, norm, pseudocount = if (norm %in% c('clr', 'log')) pc))
    if (norm == 'rclr') y[x == 0] <- NA
    y
  }


  # -----------------------------------------------------------------------
  # da_wilcoxon() = wilcox.test()
  # -----------------------------------------------------------------------

  # With ties and both sides under 50, wilcox.test() gives the exact
  # p-value given the ties from R 4.6.0 on, and the normal approximation
  # before. On older R, the reference for those is this one: the
  # distribution of the doubled rank sum of every n1-subset, counted by
  # adding one rank at a time; twice the smaller tail.
  exact_ties_p <- function (u, v) {
    a  <- as.integer(round(2 * rank(c(u, v))))
    n1 <- length(u)
    S  <- sum(a)
    cnt <- matrix(0, n1 + 1L, S + 1L)
    cnt[1, 1] <- 1
    for (ai in a) for (k in n1:1)
      cnt[k + 1, (ai + 1):(S + 1)] <- cnt[k + 1, (ai + 1):(S + 1)] + cnt[k, 1:(S + 1 - ai)]
    d <- cnt[n1 + 1, ]; s <- 0:S; obs <- sum(a[seq_len(n1)])
    min(1, 2 * min(sum(d[s <= obs]), sum(d[s >= obs])) / sum(d))
  }
  new_wilcox <- getRversion() >= '4.6.0'

  ref_w <- function (y, a, b) t(vapply(seq_len(ncol(y)), function (f) {
    u <- y[a, f]; v <- y[b, f]; u <- u[!is.na(u)]; v <- v[!is.na(v)]
    if (!length(u) || !length(v)) return (c(NaN, NaN, NaN))
    w <- suppressWarnings(stats::wilcox.test(u, v))
    p <- w$p.value
    if (!new_wilcox && length(u) < 50 && length(v) < 50 && anyDuplicated(c(u, v)))
      p <- exact_ties_p(u, v)
    c(w$statistic, p, stats::median(u) - stats::median(v))
  }, numeric(3)))

  for (norm in c('clr', 'percent', 'rclr', 'log', 'none')) {
    r <- da_wilcoxon(x, 'grp', meta, norm = norm, pseudocount = pc, cpus = 1L)
    e <- ref_w(nx(x, norm), meta$grp == 'a', meta$grp == 'b')
    expect_identical(r$w, e[, 1], info = norm)
    expect_true(maxrel(r$p_value, e[, 2]) < 1e-12, info = norm)
    expect_true(maxrel(r$median_diff, e[, 3]) < 1e-12, info = norm)
    if (norm != 'rclr') expect_true(maxrel(r$auc, e[, 1] / (12 * 12)) < 1e-15, info = norm)
    expect_true(maxrel(r$p_adjusted, stats::p.adjust(e[, 2], 'holm')) < 1e-12, info = norm)
  }

  # The exact branch (5 vs 6, no ties), the normal one (60 vs 60), and ties.
  set.seed(2)
  for (sz in list(c(5, 6), c(60, 60), c(3, 49), c(49, 50))) {
    xx <- matrix(rexp(sum(sz) * 6) * 100, sum(sz), 6)
    g  <- rep(c('a', 'b'), sz)
    r  <- da_wilcoxon(xx, g, norm = 'none', cpus = 1L)
    e  <- ref_w(xx, g == 'a', g == 'b')
    expect_identical(r$w, e[, 1], info = paste(sz, collapse = ' '))
    expect_true(maxrel(r$p_value, e[, 2]) < 1e-12, info = paste(sz, collapse = ' '))
  }
  for (sz in list(c(10, 12), c(49, 49), c(49, 50), c(3, 40), c(1, 6))) {
    xx <- matrix(sample(0:3, sum(sz) * 5, TRUE), sum(sz), 5)
    g  <- rep(c('a', 'b'), sz)
    r  <- da_wilcoxon(xx, g, norm = 'none', cpus = 1L)
    e  <- ref_w(xx, g == 'a', g == 'b')
    expect_identical(r$w, e[, 1], info = paste(sz, collapse = ' '))
    expect_true(maxrel(r$p_value, e[, 2]) < 1e-12, info = paste(sz, collapse = ' '))
  }

  # The tie distribution itself, against every subset.
  set.seed(5)
  for (it in 1:20) {
    sz <- sample(2:7, 2, TRUE)
    xx <- matrix(sample(0:2, sum(sz), TRUE) + 0, ncol = 1)
    xx[1] <- 3                                     # not all tied
    g  <- rep(c('a', 'b'), sz)
    rk <- rank(xx)
    Ws <- utils::combn(sum(sz), sz[[1]], function (i) sum(rk[i]))
    ob <- sum(rk[g == 'a'])
    e  <- min(1, 2 * min(mean(Ws <= ob + 1e-9), mean(Ws >= ob - 1e-9)))
    expect_true(abs(da_wilcoxon(xx, g, norm = 'none')$p_value - e) < 1e-12, info = it)
  }


  # -----------------------------------------------------------------------
  # da_ttest() = t.test()
  # -----------------------------------------------------------------------

  ref_t <- function (y, a, b, ...) t(vapply(seq_len(ncol(y)), function (f) {
    u <- y[a, f]; v <- y[b, f]; u <- u[!is.na(u)]; v <- v[!is.na(v)]
    tt <- stats::t.test(u, v, ...)
    c(tt$statistic, tt$parameter, tt$p.value, tt$conf.int, tt$estimate[[1]] - tt$estimate[[2]],
      mean(c(u, v)), sqrt(((length(u) - 1) * stats::var(u) + (length(v) - 1) * stats::var(v)) /
                            (length(u) + length(v) - 2)))
  }, numeric(8)))

  cols <- c('t', 'dof', 'p_value', 'ci_lower', 'ci_upper', 'diff', 'mean', 'pooled_sd')
  for (norm in c('clr', 'percent', 'rclr', 'log')) {
    for (ve in c(FALSE, TRUE)) {
      r <- da_ttest(x, 'grp', meta, norm = norm, pseudocount = pc, var_equal = ve,
                    conf_level = 0.9, cpus = 1L)
      e <- ref_t(nx(x, norm), meta$grp == 'a', meta$grp == 'b', var.equal = ve, conf.level = 0.9)
      for (k in seq_along(cols))
        expect_true(maxrel(r[[cols[[k]]]], e[, k], floor = 1e-8) < 1e-12, info = paste(norm, ve, cols[[k]]))
    }
  }

  # Contrasts: reversed, one group against the rest, two of three groups.
  y <- nx(x, 'clr')
  r <- da_ttest(x, 'g3', meta, contrast = c('mid', 'hi'), pseudocount = pc, cpus = 1L)
  e <- ref_t(y, meta$g3 == 'mid', meta$g3 == 'hi')
  expect_true(maxrel(r$t, e[, 1]) < 1e-12 && maxrel(r$p_value, e[, 3]) < 1e-12)
  expect_identical(attr(r, 'n'), c(treatment = 8L, reference = 8L))
  expect_identical(attr(r, 'contrast'), c(treatment = 'mid', reference = 'hi'))
  r <- da_ttest(x, 'g3', meta, contrast = c('mid', NA), pseudocount = pc, cpus = 1L)
  e <- ref_t(y, meta$g3 == 'mid', meta$g3 != 'mid')
  expect_true(maxrel(r$t, e[, 1]) < 1e-12)
  r <- da_ttest(x, 'g3', meta, pseudocount = pc, cpus = 1L)   # 'hi' (sorted first) vs the rest
  e <- ref_t(y, meta$g3 == 'hi', meta$g3 != 'hi')
  expect_true(maxrel(r$t, e[, 1]) < 1e-12)
  r1 <- da_wilcoxon(x, 'grp', meta, contrast = c('b', 'a'), pseudocount = pc, cpus = 1L)
  r2 <- da_wilcoxon(x, 'grp', meta, contrast = c('a', 'b'), pseudocount = pc, cpus = 1L)
  expect_identical(r1$w, 144 - r2$w)
  expect_identical(r1$median_diff, -r2$median_diff)
  expect_true(maxrel(r1$p_value, r2$p_value) < 1e-12)


  # -----------------------------------------------------------------------
  # da_lm() = summary(lm()), confint(lm())
  # -----------------------------------------------------------------------

  ref_lm <- function (y, fml, data) do.call(rbind, lapply(seq_len(ncol(y)), function (f) {
    data$y <- y[, f]
    m  <- stats::lm(stats::update(fml, y ~ .), data)
    cf <- stats::coef(m)
    s  <- matrix(NaN, length(cf), 4, dimnames = list(names(cf), NULL))
    sc <- summary(m)$coefficients
    s[rownames(sc), ] <- sc
    ci <- matrix(NaN, length(cf), 2)
    ci[!is.na(cf), ] <- stats::confint(m)[!is.na(cf), ]
    cbind(cf, s[, 2:4, drop = FALSE], ci, sigma = summary(m)$sigma, df = m$df.residual)
  }))

  m2 <- meta
  m2$grp <- factor(m2$grp); m2$g3 <- factor(m2$g3, levels = c('hi', 'lo', 'mid')); m2$sex <- factor(m2$sex)
  for (fml in list(~ grp, ~ grp + age + g3, ~ g3 * age, ~ 0 + g3 + sex, ~ age)) {
    for (norm in c('clr', 'percent', 'log')) {
      r <- da_lm(x, fml, meta, norm = norm, pseudocount = pc, conf_level = 0.9, cpus = 1L)
      e <- ref_lm(nx(x, norm), fml, m2)
      e[, 5:6] <- do.call(rbind, lapply(seq_len(D), function (f) {   # conf_level 0.9
        m2$y <- nx(x, norm)[, f]; stats::confint(stats::lm(stats::update(fml, y ~ .), m2), level = 0.9) }))
      info <- paste(format(fml), norm)
      expect_identical(r$covariate, rep(colnames(attr(r, 'design')), D), info = info)
      expect_identical(unname(r$covariate), rep(rownames(e)[seq_len(ncol(attr(r, 'design')))], D), info = info)
      k <- c('estimate', 'se', 't', 'p_value', 'ci_lower', 'ci_upper')
      for (j in seq_along(k))
        expect_true(maxrel(r[[k[[j]]]], e[, j], floor = 1e-8) < 1e-10, info = paste(info, k[[j]]))
      P <- ncol(attr(r, 'design'))
      expect_true(maxrel(attr(r, 'resid_sd'), e[seq(1, nrow(e), P), 'sigma']) < 1e-10, info = info)
      expect_identical(unname(attr(r, 'dof')), unname(e[seq(1, nrow(e), P), 'df']) + 0, info = info)
      expect_true(maxrel(attr(r, 'mean'), colMeans(nx(x, norm))) < 1e-13, info = info)
    }
  }

  # rclr: each feature fitted to the samples it was observed in, with its
  # own degrees of freedom and, where a level is missing, aliased
  # coefficients (lm()'s NA, here NaN).
  xr <- x
  xr[meta$g3 == 'mid', 5] <- 0L      # F5: no 'mid' samples
  xr[, 6] <- 0L; xr[1:3, 6] <- 1:3   # F6: three samples
  xr[, 7] <- 0L                      # F7: none
  r <- da_lm(xr, ~ g3 + age, meta, norm = 'rclr', cpus = 1L)
  y <- nx(xr, 'rclr')
  X <- attr(r, 'design')
  e <- do.call(rbind, lapply(seq_len(D), function (f) {
    ok <- !is.na(y[, f])
    if (!any(ok)) return (matrix(NaN, ncol(X), 4))
    m  <- stats::lm(y[ok, f] ~ X[ok, , drop = FALSE] - 1)
    cf <- stats::coef(m)
    s  <- matrix(NaN, length(cf), 3)
    sc <- summary(m)$coefficients
    s[!is.na(cf), ] <- sc[, 2:4]
    cbind(cf, s)
  }))
  e[is.na(e)] <- NaN
  expect_true(maxrel(r$estimate, e[, 1], floor = 1e-8) < 1e-10)
  expect_true(maxrel(r$se, e[, 2], floor = 1e-8) < 1e-10)
  expect_true(maxrel(r$p_value, e[, 4], floor = 1e-8) < 1e-10)
  expect_true(is.nan(r$estimate[r$feature == 'F5' & r$covariate == 'g3mid']))
  expect_identical(unname(attr(r, 'dof')[c('F6', 'F7')]), c(NaN, NaN))   # 3 samples, 4 columns
  expect_false(any(r$reject[r$feature %in% c('F6', 'F7')]))

  # A rank-deficient design warns, and its aliased column is NaN as lm()'s NA.
  m3 <- meta; m3$age2 <- 2 * m3$age
  expect_warning(r <- da_lm(x, ~ age + age2, m3, pseudocount = pc, cpus = 1L), 'rank deficient')
  expect_true(all(is.nan(r$estimate[r$covariate == 'age2'])))
  e <- ref_lm(nx(x, 'clr'), ~ age + age2, m3)
  expect_true(maxrel(r$estimate[r$covariate != 'age2'], e[rownames(e) != 'age2', 1]) < 1e-10)

  # One categorical variable: da_lm()'s coefficient is the t-test's
  # difference (pooled variance), with the same p-value.
  r1 <- da_lm(x, groups = 'grp', data = meta, pseudocount = pc, cpus = 1L)
  r2 <- da_ttest(x, 'grp', meta, contrast = c('b', 'a'), var_equal = TRUE, pseudocount = pc, cpus = 1L)
  expect_true(maxrel(r1$estimate[r1$covariate == 'grpb'], r2$diff) < 1e-10)
  expect_true(maxrel(r1$p_value[r1$covariate == 'grpb'], r2$p_value) < 1e-10)
  expect_identical(attr(r1, 'groups'), 'grp')
  expect_identical(format(attr(r1, 'formula')), '~grp')

  # `groups` as values, alone: the variable `group`.
  r3 <- da_lm(x, groups = meta$grp, pseudocount = pc, cpus = 1L)
  expect_identical(unique(r3$covariate), c('(Intercept)', 'groupb'))
  expect_identical(r3$estimate, r1$estimate)


  # -----------------------------------------------------------------------
  # The OLS kernel: classical, HC0 and HC3 standard errors, missing values
  # and aliased columns.
  # -----------------------------------------------------------------------

  set.seed(3)
  X <- cbind(1, rnorm(30), rep(0:1, 15), rnorm(30))
  Y <- matrix(rnorm(30 * 5), 30, 5)
  Y[c(2, 9), 2] <- NaN
  Y[X[, 3] == 1, 3] <- NaN            # column 3 aliased for this response
  for (st in 0:2) {
    k <- .Call(C_feature_lm, X, Y, 1e-7, st)
    for (j in 1:5) {
      ok <- !is.na(Y[, j])
      Xo <- X[ok, , drop = FALSE]
      m  <- stats::lm(Y[ok, j] ~ Xo - 1)
      cf <- unname(stats::coef(m)); cf[is.na(cf)] <- NaN
      expect_true(maxrel(k$coef[j, ], cf) < 1e-11, info = paste(st, j))
      expect_identical(k$rank[[j]], m$rank, info = paste(st, j))
      expect_true(maxrel(k$rss[[j]], sum(stats::resid(m)^2)) < 1e-11, info = paste(st, j))
      keep <- !is.na(stats::coef(m))
      if (st == 0) {
        se <- rep(NaN, 4); se[keep] <- summary(m)$coefficients[, 2]
      }
      else {
        Xk <- Xo[, keep, drop = FALSE]
        e  <- stats::resid(m)
        h  <- stats::hatvalues(m)
        w  <- if (st == 1) e^2 else e^2 / (1 - h)^2
        B  <- solve(crossprod(Xk))
        V  <- B %*% crossprod(Xk * w, Xk) %*% B
        se <- rep(NaN, 4); se[keep] <- sqrt(diag(V))
        expect_true(maxrel(k$resid[ok, j], unname(e), floor = 1e-8) < 1e-10, info = paste(st, j))
      }
      expect_true(maxrel(k$se[j, ], se) < 1e-10, info = paste(st, j))
    }
  }


  # -----------------------------------------------------------------------
  # Invariance, threads, input types
  # -----------------------------------------------------------------------

  # Multiplying a sample by a constant changes nothing (clr, percent, log;
  # without zeros, so without a pseudocount).
  xi <- x + 1L
  xs <- xi; xs[3, ] <- xs[3, ] * 7
  for (norm in c('clr', 'percent', 'log')) {
    a <- da_ttest(xi, 'grp', meta, norm = norm, cpus = 1L)
    b <- da_ttest(xs, 'grp', meta, norm = norm, cpus = 1L)
    expect_true(maxrel(a$t, b$t, floor = 1e-10) < 1e-10, info = norm)
    a <- da_wilcoxon(xi, 'grp', meta, norm = norm, cpus = 1L)
    b <- da_wilcoxon(xs, 'grp', meta, norm = norm, cpus = 1L)
    expect_identical(a$w, b$w, info = norm)
  }
  xp <- x / rowSums(x)
  a <- da_lm(x,  ~ grp + age, meta, norm = 'percent', cpus = 1L)
  b <- da_lm(xp, ~ grp + age, meta, norm = 'percent', cpus = 1L)
  expect_true(maxrel(a$estimate, b$estimate, floor = 1e-10) < 1e-10)
  a <- da_wilcoxon(x * 3, 'grp', meta, norm = 'rclr', cpus = 1L)
  b <- da_wilcoxon(x,     'grp', meta, norm = 'rclr', cpus = 1L)
  expect_true(maxrel(a$median_diff, b$median_diff, floor = 1e-10) < 1e-10)

  # Identical for any number of threads, with enough work to use them.
  set.seed(4)
  big  <- matrix(rpois(100 * 4000, 2), 100, 4000)
  bgrp <- rep(c('u', 'v'), 50)
  bm   <- data.frame(g = bgrp, z = rnorm(100))
  bf   <- ~ g + z   # one formula: the result keeps it, and its environment
  for (f in list(
      function (cpus) da_wilcoxon(big, bgrp, pseudocount = 1, cpus = cpus),
      function (cpus) da_ttest(big, bgrp, norm = 'rclr', cpus = cpus),
      function (cpus) da_lm(big, bf, bm, norm = 'rclr', cpus = cpus),
      function (cpus) da_lm(big, bf, bm, pseudocount = 1, cpus = cpus) )) {
    r1 <- f(1L)
    expect_identical(f(2L), r1)
    expect_identical(f(4L), r1)
  }

  # Dense = Matrix = slam = margin 2.
  tx <- t(x)
  fa <- ~ grp + age
  for (fn in list(
      function (m, ...) da_wilcoxon(m, meta$grp, pseudocount = pc, cpus = 1L, ...),
      function (m, ...) da_ttest(m, meta$grp, norm = 'rclr', cpus = 1L, ...),
      function (m, ...) da_lm(m, fa, meta, norm = 'log', pseudocount = pc, cpus = 1L, ...) )) {
    r <- fn(x)
    expect_identical(fn(tx, margin = 2L), r)
    if (requireNamespace('Matrix', quietly = TRUE)) {
      expect_identical(fn(Matrix::Matrix(x, sparse = TRUE)), r)
      expect_identical(fn(Matrix::Matrix(tx, sparse = TRUE), margin = 2L), r)
    }
    if (requireNamespace('slam', quietly = TRUE))
      expect_identical(fn(slam::as.simple_triplet_matrix(x)), r)
  }


  # -----------------------------------------------------------------------
  # Edge cases
  # -----------------------------------------------------------------------

  # A side with one sample: Welch's t is NaN, with a warning; Student's is
  # defined; Wilcoxon's works.
  g1 <- c('a', rep('b', 23))
  expect_warning(r <- da_ttest(x, g1, pseudocount = pc, cpus = 1L), 'single sample')
  expect_true(all(is.nan(r$t)) && all(is.nan(r$p_value)) && !any(r$reject))
  r <- da_ttest(x, g1, pseudocount = pc, var_equal = TRUE, cpus = 1L)
  e <- ref_t(nx(x, 'clr'), g1 == 'a', g1 == 'b', var.equal = TRUE)
  expect_true(maxrel(r$t, e[, 1]) < 1e-12)
  r <- da_wilcoxon(x, g1, pseudocount = pc, cpus = 1L)
  expect_true(maxrel(r$p_value, ref_w(nx(x, 'clr'), g1 == 'a', g1 == 'b')[, 2]) < 1e-12)

  # A constant feature: NaN, never rejected (t.test() stops instead).
  xc <- x; xc[, 2] <- 5L
  r <- da_ttest(xc, 'grp', meta, norm = 'none', cpus = 1L)
  expect_true(is.nan(r$t[[2]]) && is.nan(r$p_value[[2]]) && !r$reject[[2]])
  expect_identical(r$p_adjusted[-2], stats::p.adjust(r$p_value[-2], 'holm'))
  r <- da_wilcoxon(xc, 'grp', meta, norm = 'none', cpus = 1L)
  expect_true(r$p_value[[2]] == 1 && !r$reject[[2]])   # every rank sum the same
  r <- da_lm(xc, ~ grp, meta, norm = 'none', cpus = 1L)
  expect_true(all(abs(r$estimate[r$feature == 'F2'] - c(5, 0)) < 1e-12))
  expect_true(all(r$se[r$feature == 'F2'] < 1e-12))

  # An all-zero feature with rclr: no values, NaN.
  xz <- x; xz[, 3] <- 0L
  r <- da_wilcoxon(xz, 'grp', meta, norm = 'rclr', cpus = 1L)
  expect_true(is.nan(r$w[[3]]) && is.nan(r$p_value[[3]]) && !r$reject[[3]])
  r <- da_ttest(xz, 'grp', meta, norm = 'rclr', cpus = 1L)
  expect_true(is.nan(r$t[[3]]) && is.nan(r$mean[[3]]))

  # Negative values only with norm = 'none'.
  r <- da_ttest(scale(x), 'grp', meta, norm = 'none', cpus = 1L)
  expect_true(all(is.finite(r$t)))
  expect_error(da_ttest(scale(x), 'grp', meta), 'negative')

  # No feature names: their indices.
  expect_identical(da_wilcoxon(unname(x), meta$grp, pseudocount = 1)$feature, as.character(1:D))

  # Argument checks.
  gna <- meta$grp; gna[2] <- NA
  expect_error(da_ttest(x, gna, pseudocount = 1), 'NA')
  expect_error(da_ttest(x, pseudocount = 1), '`groups` is required')
  expect_error(da_ttest(x, 'grp', meta, norm = 'nope'), 'arg')
  expect_error(da_ttest(x, 'grp', meta, norm = 1), '`norm`')
  expect_error(da_ttest(x, 'grp', meta, p_adjust = 'nope', pseudocount = 1), 'p_adjust')
  for (v in list(0, 1, NA, 'a', c(0.1, 0.2))) {
    expect_error(da_wilcoxon(x, 'grp', meta, alpha = v, pseudocount = 1), '`alpha`')
    expect_error(da_ttest(x, 'grp', meta, conf_level = v, pseudocount = 1), '`conf_level`')
    expect_error(da_lm(x, ~ grp, meta, conf_level = v, pseudocount = 1), '`conf_level`')
  }
  expect_error(da_ttest(x, 'grp', meta, var_equal = NA, pseudocount = 1), '`var_equal`')
  expect_error(da_ttest(x, 'grp', meta, contrast = c('a', 'z'), pseudocount = 1), 'not one of')
  expect_error(da_ttest(x, rep('a', n), pseudocount = 1), 'two different')
  expect_error(da_lm(x, pseudocount = 1), '`formula` or `groups` is required')
  expect_error(da_lm(x, ~ grp, pseudocount = 1), '`data` is required')
  expect_error(da_lm(x, ~ grp + nope, meta, pseudocount = 1), 'nope')
  expect_error(da_lm(x, ~ grp, meta, groups = 'age', pseudocount = 1), 'term in `formula`')
  expect_error(da_lm(x, ~ grp, meta, groups = meta$grp, pseudocount = 1), 'term of interest')
  expect_error(da_lm(x[1:3, ], ~ age + g3, meta, pseudocount = 1), 'more samples than design columns')
  xna <- x + 0; xna[1, 1] <- NA
  expect_error(da_ttest(xna, 'grp', meta), 'NA')
  expect_warning(da_ttest(x, 'grp', meta), 'pseudocount')


  # -----------------------------------------------------------------------
  # Attributes, print and plot
  # -----------------------------------------------------------------------

  r <- da_ttest(x, 'grp', meta, pseudocount = pc, cpus = 1L)
  expect_identical(class(r), c('ecodive_da_ttest', 'data.frame'))
  expect_identical(names(r), c('feature', 'diff', 'ci_lower', 'ci_upper', 't', 'dof',
                               'pooled_sd', 'mean', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(attr(r, 'method'), 'da_ttest')
  expect_identical(attr(r, 'log_base'), exp(1))
  expect_identical(attr(r, 'alpha'), 0.05)
  expect_identical(attr(r, 'conf_level'), 0.95)
  expect_identical(attr(r, 'norm'), 'clr')
  expect_identical(attr(r, 'pseudocount'), pc)
  expect_true(is.na(attr(da_ttest(x, 'grp', meta, norm = 'percent'), 'log_base')))
  expect_null(attr(r, 'display'))
  expect_true(any(grepl('Welch t-test', utils::capture.output(print(r)))))

  r <- da_wilcoxon(x, 'grp', meta, pseudocount = pc, cpus = 1L)
  expect_identical(names(r), c('feature', 'w', 'auc', 'median_diff', 'mean', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(class(r), c('ecodive_da_wilcoxon', 'data.frame'))
  expect_true(any(grepl('Wilcoxon Rank-Sum Test', utils::capture.output(print(r)))))

  rl <- da_lm(x, ~ grp + g3, meta, groups = 'g3', pseudocount = pc, cpus = 1L)
  expect_identical(names(rl), c('feature', 'covariate', 'estimate', 'se', 't', 'ci_lower',
                                'ci_upper', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(attr(rl, 'method'), 'da_lm')
  expect_identical(names(attr(rl, 'resid_sd')), colnames(x))
  expect_true(any(grepl('Linear Model', utils::capture.output(print(rl)))))

  expect_identical(ecodive:::da_plot_title(r),  'Volcano plot of Wilcoxon Rank-Sum Test')
  expect_identical(ecodive:::da_plot_title(rl), 'Volcano plot of Linear Model')
  grDevices::pdf(NULL)
  expect_identical(withVisible(plot(r))$visible, FALSE)
  expect_identical(plot(da_ttest(x, 'grp', meta, pseudocount = pc)), da_ttest(x, 'grp', meta, pseudocount = pc))
  expect_silent(plot(rl))
  expect_silent(plot(rl, 'grpb', labels = 0, main = 'x'))
  expect_error(plot(rl, 'nope'), '`covariate`')
  grDevices::dev.off()


  # -----------------------------------------------------------------------
  # da_test(): the method's own result
  # -----------------------------------------------------------------------

  xi <- x + 1L
  expect_identical(da_test(xi, 'ancom', 'grp', meta, cpus = 1L), ancom(xi, 'grp', meta, cpus = 1L))
  expect_identical(da_test(xi, 'ancombc', 'grp', meta, ~ grp + age, cpus = 1L),
                   ancombc(xi, ~ grp + age, meta, groups = 'grp', cpus = 1L))
  expect_identical(da_test(x, 'ancombc2', formula = ~ grp + age, data = meta, cpus = 1L),
                   ancombc2(x, ~ grp + age, meta, cpus = 1L))
  expect_identical(da_test(xi, 'dirmult_ttest', 'grp', meta, n_draws = 8L, seed = 3L, cpus = 1L),
                   dirmult_ttest(xi, 'grp', meta, n_draws = 8L, seed = 3L, cpus = 1L))
  expect_identical(da_test(xi, 'dirmult_lm', 'grp', meta, ~ grp + age, n_draws = 8L, seed = 3L, cpus = 1L),
                   dirmult_lm(xi, ~ grp + age, meta, groups = 'grp', n_draws = 8L, seed = 3L, cpus = 1L))
  expect_identical(da_test(x, 'linda', 'grp', meta, ~ grp + age, pseudocount = 1, cpus = 1L),
                   linda(x, ~ grp + age, meta, groups = 'grp', pseudocount = 1, cpus = 1L))
  expect_identical(da_test(x, 'lefse', 'grp', meta, n_boots = 4L, seed = 2L, cpus = 1L),
                   lefse(x, 'grp', meta, n_boots = 4L, seed = 2L, cpus = 1L))
  expect_identical(da_test(x, 'da_wilcoxon', 'grp', meta, pseudocount = 1, cpus = 1L),
                   da_wilcoxon(x, 'grp', meta, pseudocount = 1, cpus = 1L))
  expect_identical(da_test(x, 'da_ttest', meta$grp, norm = 'percent', seed = 9L, cpus = 1L),
                   da_ttest(x, meta$grp, norm = 'percent', cpus = 1L))
  expect_identical(da_test(x, 'lm', 'grp', meta, ~ grp + age, pseudocount = 1, cpus = 1L),
                   da_lm(x, ~ grp + age, meta, groups = 'grp', pseudocount = 1, cpus = 1L))
  expect_identical(da_test(counts = x, method = 'da_lm', groups = 'grp', data = meta, pseudocount = 1),
                   da_lm(x, data = meta, groups = 'grp', pseudocount = 1))

  # A one-variable formula is `groups` for a group method; more is an error.
  expect_identical(da_test(x, 'da_ttest', data = meta, formula = ~ g3, pseudocount = 1, cpus = 1L),
                   da_ttest(x, 'g3', meta, pseudocount = 1, cpus = 1L))
  expect_identical(da_test(xi, 'ancom', data = meta, formula = 'grp', cpus = 1L),
                   ancom(xi, 'grp', meta, cpus = 1L))
  expect_error(da_test(xi, 'ancom', data = meta, formula = ~ grp + age), 'compares one grouping')
  expect_error(da_test(xi, 'ancom', 'grp', meta, ~ grp), 'not both')
  expect_error(da_test(xi, 'ancom', data = meta, formula = ~ log(age)), 'single variable')
  expect_error(da_test(xi, 'da_ttest', meta$grp, n_draws = 3), 'unused argument')

  # ancombc()'s `groups` alone is the model ~ groups.
  a <- ancombc(xi, data = meta, groups = 'grp', cpus = 1L)
  b <- ancombc(xi, ~ grp, meta, cpus = 1L)
  expect_identical(lapply(a, identity), lapply(b, identity))
  expect_identical(attr(a, 'coef'), attr(b, 'coef'))
  expect_identical(format(attr(a, 'formula')), '~grp')
  expect_identical(attr(a, 'groups'), 'grp')
  c2 <- ancombc(xi, groups = meta$grp, cpus = 1L)
  expect_identical(unique(c2$covariate), c('(Intercept)', 'groupb'))
  expect_identical(c2$log_fc, a$log_fc)


  # -----------------------------------------------------------------------
  # The registry: list_da_tests(), match_da_test(), aliases
  # -----------------------------------------------------------------------

  df <- list_da_tests()
  expect_identical(names(df), c('name', 'id', 'design', 'effect', 'base', 'zero_free_ok', 'signed'))
  expect_identical(df$id[!df$signed], c('ancom', 'lefse', 'indicator_species'))
  expect_identical(sort(df$id), sort(c('ancom', 'ancombc', 'ancombc2', 'dirmult_ttest', 'dirmult_lm',
                                       'da_wilcoxon', 'da_ttest', 'da_lm', 'linda', 'lefse',
                                       'indicator_species')))
  expect_true(all(df$design %in% c('groups', 'formula')) && is.logical(df$zero_free_ok))
  expect_identical(rownames(list_da_tests(nm = 'id')), df$id)
  expect_identical(names(list_da_tests('func')), df$id)
  expect_identical(list_da_tests('func')$da_lm, da_lm)
  expect_identical(list_da_tests('list')$da_ttest$params, names(formals(da_ttest)))
  expect_identical(list_da_tests('name', 'id')[['da_ttest']], 'Welch t-test')
  for (id in df$id) expect_identical(match_da_test(id)$id, id, info = id)

  # Each name and alias selects its method.
  for (i in seq_len(nrow(DA_TESTS))) {
    nms <- c(DA_TESTS$name[[i]], DA_TESTS$id[[i]], stats::na.omit(strsplit(DA_TESTS$alt_ids[[i]], ' ')[[1]]))
    for (nm in nms) expect_identical(match_da_test(nm)$id, DA_TESTS$id[[i]], info = nm)
  }
  expect_identical(match_da_test('ANCOM-BC2')$id, 'ancombc2')
  expect_identical(match_da_test('ANCOM-BC')$id,  'ancombc')
  expect_identical(match_da_test('ANCOM')$id,     'ancom')
  expect_identical(match_da_test('ALDEx2')$id,    'dirmult_ttest')
  expect_identical(match_da_test('ALDEx3')$id,    'dirmult_lm')
  expect_identical(match_da_test('aldex')$id,     'dirmult_ttest')
  expect_error(match_da_test('dirmult'), 'multiple')
  expect_identical(match_da_test('Mann-Whitney U')$id, 'da_wilcoxon')
  expect_identical(match_da_test('t-test')$id,    'da_ttest')
  expect_identical(match_da_test('linear')$id,    'da_lm')
  expect_error(match_da_test('anc'), 'multiple')
  expect_error(match_da_test('zzz'), 'does not match')
  expect_error(match_da_test('---'), 'alphabetic')
  expect_error(match_da_test(1), 'character')
  expect_error(match_da_test(NA_character_), 'NA')

  # No name or alias of one wrapper is a name or alias of another, and each
  # wrapper sends the others' methods to their wrapper.
  regs <- c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'ZERO_REPLACEMENTS')
  for (reg in regs)
    expect_identical(intersect(names(exact_names('DA_TESTS')), names(exact_names(reg))), character(0), info = reg)
  expect_error(match_da_test('permanova'), 'dist_test')
  expect_error(match_da_test('pca'), 'count_ord')
  expect_error(match_da_test('pcoa'), 'dist_ord')
  expect_error(match_count_ord('wilcoxon'), 'da_test')
  expect_error(match_dist_ord('ancombc'), 'da_test')
  expect_error(match_dist_test('aldex2'), 'da_test')
  expect_error(da_test(x, 'permanova', 'grp', meta), 'dist_test')

  expect_identical(ecodive:::da_method_name('da_lm'), 'Linear Model')
  expect_identical(ecodive:::da_method_name('zzz'), 'zzz')

#})
