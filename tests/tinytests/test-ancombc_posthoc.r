#test_that("ancombc_posthoc", {

  # ancombc_global(), ancombc_pairwise(), ancombc_dunnett(), ancombc_trend()
  # and structural_zeros(): scikit-bio 0.7.4's global_test(),
  # pairwise_test(), dunnett_test(), trend_test() and struc_zero()
  # (data-raw/expected/ancombc.py), including the bootstraps on scikit-bio's
  # own null draws; ANCOMBC 2.8.1's deterministic post-hoc results where its
  # EM converges; the samplers by distribution; the QP solver against brute
  # force; threading, seeds and argument checks.

  ref <- readRDS(file.path('expected', 'ancombc.rds'))

  fit_of <- function (cs) {
    s  <- ref$skbio[[cs]]
    tb <- ref$tables[[s$table]]
    fn <- if (s$method == 'ancombc') ancombc else ancombc2
    args <- list(
      tb$counts, s$formula, tb$data, groups = s$groups,
      pseudocount = s$pseudocount, max_iters = s$max_iters, tol = s$tol,
      alpha = s$alpha, p_adjust = s$p_adjust, sensitivity = FALSE, cpus = 1L )
    if (s$method == 'ancombc2') args$var_quantile <- s$var_quantile
    do.call(fn, args)
  }

  maxabs <- function (a, b) { ok <- is.finite(a) & is.finite(b); max(c(0, abs(a[ok] - b[ok]))) }
  maxrel <- function (a, b) {
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300)))
  }
  same_na <- function (a, b) identical(is.finite(a), is.finite(b))


  # -----------------------------------------------------------------------
  # scikit-bio
  # -----------------------------------------------------------------------

  cases <- names(ref$skbio)[!vapply(ref$skbio, function (s) is.null(s$groups), logical(1))]
  expect_true(length(cases) >= 12L)

  for (cs in cases) {
    s   <- ref$skbio[[cs]]
    fit <- fit_of(cs)

    # Global test.
    g <- ancombc_global(fit, cpus = 1L)
    expect_inherits(g, c('ecodive_ancombc_global', 'data.frame'), info = cs)
    expect_true(same_na(g$w, s$global$w), info = cs)
    expect_true(maxrel(g$w, s$global$w) < 1e-7, info = cs)
    expect_true(maxrel(g$p_value, s$global$p) < 1e-6, info = cs)
    expect_true(maxrel(g$p_adjusted, s$global$q) < 1e-6, info = cs)
    expect_identical(g$reject, s$global$reject, info = cs)

    # Pairwise: rows in scikit-bio's order.
    r <- ancombc_pairwise(fit, cpus = 1L)
    e <- s$pairwise
    expect_identical(nrow(r), nrow(e), info = cs)
    expect_true(maxabs(r$log_fc, e$log_fc) < 1e-7, info = cs)
    expect_true(same_na(r$se, e$se) && maxrel(r$se, e$se) < 1e-8, info = cs)
    expect_equal(r$w, r$log_fc / r$se, info = cs)
    expect_true(maxrel(r$p_value, e$p) < 1e-6, info = cs)
    expect_true(maxrel(r$p_adjusted, e$q) < 1e-6, info = cs)
    expect_identical(r$reject, e$reject, info = cs)

    # Dunnett's: the level-against-reference rows of the pairwise table.
    d <- ancombc_dunnett(fit, n_bootstraps = 20L, cpus = 1L)
    k <- length(unique(d$comparison))
    m <- length(unique(r$comparison))
    first <- rep(seq_len(k), nrow(d) / k) + rep(seq(0, by = m, length.out = nrow(d) / k), each = k)
    expect_identical(d$comparison, r$comparison[first], info = cs)
    expect_identical(d$log_fc, r$log_fc[first], info = cs)
    expect_identical(d$se, r$se[first], info = cs)

    # Dunnett's and the trend test on scikit-bio's own null draws.
    if (!is.null(s$dunnett_draws)) {
      r <- ancombc_dunnett_details(fit, n_bootstraps = 20L, cpus = 1L, null_draws = s$dunnett_draws)
      expect_true(maxrel(r$p_value, s$dunnett$p) < 1e-6, info = cs)
      expect_true(maxrel(r$p_adjusted, s$dunnett$q) < 1e-6, info = cs)
      expect_identical(r$reject, s$dunnett$reject, info = cs)
    }
    tr <- if (!is.null(s$trend_draws))
            ancombc_trend_details(fit, n_bootstraps = 20L, cpus = 1L, null_draws = s$trend_draws)
          else ancombc_trend(fit, n_bootstraps = 20L, cpus = 1L)
    expect_true(same_na(tr$w, s$trend$w), info = cs)
    # SLSQP stops short of the exact solution (by up to 2e-4 in w; the
    # brute-force checks below show ecodive's is exact).
    expect_true(maxabs(tr$w, s$trend$w) < 1e-3, info = cs)
    if (!is.null(s$trend_draws)) {
      expect_identical(tr$p_value, s$trend$p, info = cs)
      expect_true(maxrel(tr$p_adjusted, s$trend$q) < 1e-8, info = cs)
      expect_identical(tr$reject, s$trend$reject, info = cs)
    }
  }

  # scikit-bio's documented values (ancombc2() docstring).
  fit <- fit_of('doc3')
  g   <- ancombc_global(fit)
  expect_equal(round(g$w, 3), c(1.786, 60.287, 0.04, 4.78, 1.385, 25.906, 0.865, 1.589))
  expect_identical(which(g$reject), c(2L, 6L))
  pw <- ancombc_pairwise(fit)
  expect_identical(pw$comparison[1:3], c('moderate - mild', 'severe - mild', 'severe - moderate'))
  expect_identical(paste(pw$feature, pw$comparison)[pw$reject],
    c('F02 moderate - mild', 'F02 severe - mild', 'F06 severe - mild', 'F06 severe - moderate'))
  expect_equal(round(pw$p_adjusted[pw$reject], 3), c(0.002, 0.001, 0.006, 0.023))
  tr <- ancombc_trend(fit)
  expect_equal(round(tr$w, 3), c(0.513, 1.880, 0.031, 0.833, 0.485, 1.729, 0.232, 0.539))
  expect_identical(tr$pattern[c(2, 6)], c('increasing', 'decreasing'))

  # The trend test's p-values with ecodive's draws, in distribution: 4,000
  # draws, against scikit-bio's 20,000.
  big <- ref$skbio$doc3$trend_big
  tr  <- ancombc_trend(fit, n_bootstraps = 4000L, seed = 5L, cpus = 2L)
  se  <- sqrt(pmax(big * (1 - big), 1e-4) * (1 / 4000 + 1 / 20000))
  expect_true(all(abs(tr$p_value - big) < 5 * se))


  # -----------------------------------------------------------------------
  # The samplers: Dunnett's null maxima are |t| (or |N(0, 1)|) values
  # -----------------------------------------------------------------------

  B <- 50000L
  for (dof in list(NULL, 3, 12.5, NaN)) {
    for (w in c(0.5, 1.5, 2.5)) {
      ex <- .Call(C_ancombc_dunnett, matrix(w), if (!is.null(dof)) dof, B, 3L, NULL, 1L)
      nu <- if (is.null(dof)) Inf else if (is.nan(dof)) 999 else dof
      p  <- 2 * stats::pt(-w, nu)
      expect_true(abs(ex / B - p) < 5 * sqrt(p * (1 - p) / B), info = paste(dof, w))
    }
  }


  # -----------------------------------------------------------------------
  # The QP solver: the trend test's constrained estimates are the exact
  # solutions, as brute force over the active sets finds them.
  # -----------------------------------------------------------------------

  brute <- function (b, V, C) {
    P <- solve((V + t(V)) / 2)
    G <- (2 * P + t(2 * P)) / 2 + diag(1e-10, length(b))
    d <- 2 * P %*% b
    best <- NULL; fbest <- Inf
    m <- nrow(C)
    for (s in 0:(2^m - 1)) {
      act <- which(bitwAnd(s, 2^(0:(m - 1))) > 0)
      if (length(act)) {
        A <- C[act, , drop = FALSE]
        K <- rbind(cbind(G, -t(A)), cbind(A, matrix(0, length(act), length(act))))
        sol <- tryCatch(solve(K, c(d, rep(0, length(act)))), error = function (e) NULL)
        if (is.null(sol)) next
        x <- sol[seq_along(b)]
        if (any(sol[-seq_along(b)] < -1e-9)) next
      }
      else x <- solve(G, d)
      if (any(C %*% x < -1e-9)) next
      f <- 0.5 * sum(x * (G %*% x)) - sum(d * x)
      if (f < fbest) { fbest <- f; best <- x }
    }
    best
  }

  for (cs in c('sim_bc2', 'sim_bc2_noint')) {
    fit <- fit_of(cs)
    tr  <- ancombc_trend(fit, n_bootstraps = 1L, cpus = 1L)
    gi  <- if (cs == 'sim_bc2') 2:3 else 1:3
    cf  <- attr(fit, 'coef')[, gi]
    vc  <- attr(fit, 'vcov')[, gi, gi]
    pat <- attr(tr, 'patterns')
    for (f in seq_len(nrow(cf))) {
      x <- brute(cf[f, ], vc[f, , ], pat[[tr$pattern[f]]])
      expect_true(max(abs(attr(tr, 'coef')[f, ] - x)) < 1e-9 * max(1, abs(x)), info = paste(cs, f))
    }
  }

  # More constraints than coefficients (a custom pattern): an umbrella
  # with its peak at the second level, and an extra bound.
  umb <- rbind(c(1, 0, 0), c(-1, 1, 0), c(0, -1, 1), c(0, 0, -1), c(1, 1, 1))
  x4   <- ref$tables$sparse$counts
  m4   <- ref$tables$sparse$data
  m4$trt4 <- rep(c('a', 'b', 'c', 'd'), length.out = nrow(x4))
  fit4 <- ancombc2(x4, ~ trt4 + age, m4, groups = 'trt4', sensitivity = FALSE, cpus = 1L)
  tr4  <- ancombc_trend(fit4, patterns = list(umbrella = umb), n_bootstraps = 10L, cpus = 1L)
  cf4  <- attr(fit4, 'coef')[, 2:4]
  vc4  <- attr(fit4, 'vcov')[, 2:4, 2:4]
  ok   <- is.finite(tr4$w)
  for (f in which(ok)) {
    x <- brute(cf4[f, ], vc4[f, , ], umb)
    expect_true(max(abs(attr(tr4, 'coef')[f, ] - x)) < 1e-8 * max(1, abs(x)), info = f)
  }
  expect_identical(attr(tr4, 'nodes'), c(umbrella = 3L))


  # -----------------------------------------------------------------------
  # Threads, seeds, groups, labels
  # -----------------------------------------------------------------------

  fit <- fit_of('sim_bc2')
  for (f in list(ancombc_global, ancombc_pairwise, ancombc_dunnett, ancombc_trend)) {
    r1 <- f(fit, cpus = 1L)
    for (cpus in c(2L, 4L)) expect_identical(f(fit, cpus = cpus), r1)
  }
  d1 <- ancombc_dunnett(fit, n_bootstraps = 300L, seed = 1L, cpus = 1L)
  expect_identical(ancombc_dunnett(fit, n_bootstraps = 300L, seed = 1L, cpus = 3L), d1)
  t1 <- ancombc_trend(fit, n_bootstraps = 300L, seed = 1L, cpus = 1L)
  t2 <- ancombc_trend(fit, n_bootstraps = 300L, seed = 2L, cpus = 1L)
  expect_identical(ancombc_trend(fit, n_bootstraps = 300L, seed = 1L, cpus = 4L), t1)
  expect_false(identical(t1$p_value, t2$p_value))
  expect_identical(t1$w, t2$w)
  expect_identical(attr(t1, 'seed'), 1L)
  expect_identical(attr(t1, 'n_bootstraps'), 300L)

  # A feature's bootstrap doesn't depend on the other features: feature f
  # uses stream f, so the first 5 features alone give the same counts.
  gi   <- 2:3
  pt   <- unname(attr(t1, 'patterns'))
  call <- function (rows) .Call(
    C_ancombc_trend, unname(attr(fit, 'coef')[rows, gi, drop = FALSE]),
    unname(attr(fit, 'var')[rows, gi, drop = FALSE]),
    unname(attr(fit, 'vcov')[rows, gi, gi, drop = FALSE]), NULL, pt,
    c(1L, 1L), 300L, 1L, NULL, 1L )
  D    <- nrow(attr(fit, 'coef'))
  full <- call(seq_len(D))
  part <- call(1:5)
  expect_identical(full$counts[seq_len(D)], as.integer(round(t1$p_value * 300)))
  expect_identical(part$counts[1:5], full$counts[1:5])
  expect_identical(part$W, full$W[1:5])

  # Labels and the grouping.
  expect_identical(unique(ancombc_dunnett(fit)$comparison), c('oral - gut', 'skin - gut'))
  expect_identical(unique(ancombc_pairwise(fit)$comparison), c('oral - gut', 'skin - gut', 'skin - oral'))
  fit0 <- fit_of('sim_bc2_noint')
  expect_identical(unique(ancombc_dunnett(fit0)$comparison), c('gut', 'oral', 'skin'))
  expect_identical(unique(ancombc_pairwise(fit0)$comparison), c('gut', 'oral', 'skin', 'oral - gut', 'skin - gut', 'skin - oral'))
  tb  <- ref$tables$sim
  fnog <- ancombc2(tb$counts, ~ site + ph + smoker, tb$data, sensitivity = FALSE, cpus = 1L)
  expect_error(ancombc_global(fnog), '`groups` is required')
  expect_identical(ancombc_global(fnog, groups = 'site'), ancombc_global(fit))
  expect_identical(ancombc_pairwise(fnog, groups = 'site'), ancombc_pairwise(fit))

  # alpha and p_adjust: the fit's by default, or given.
  g  <- ancombc_global(fit, alpha = 0.2, p_adjust = 'BH')
  expect_identical(attr(g, 'p_adjust'), 'BH')
  expect_identical(g$p_adjusted, stats::p.adjust(g$p_value, 'BH'))
  expect_identical(g$reject, g$p_adjusted <= 0.2)
  expect_identical(attr(ancombc_global(fit), 'alpha'), attr(fit, 'alpha'))
  pw <- ancombc_pairwise(fit, p_adjust = 'none')
  expect_true(all(pw$p_adjusted >= pw$p_value - 1e-15))

  # Non-estimable features.
  fs <- fit_of('sparse_bc2')
  bad <- c('F09', 'F10')
  for (r in list(ancombc_global(fs), ancombc_trend(fs))) {
    k <- r$feature %in% bad
    expect_true(all(is.nan(r$w[k])) && all(r$p_value[k] == 1) && !any(r$reject[k]))
  }
  for (r in list(ancombc_pairwise(fs), ancombc_dunnett(fs))) {
    k <- r$feature %in% bad
    expect_true(all(is.nan(r$log_fc[k])) && all(r$p_value[k] == 1) && !any(r$reject[k]))
  }

  # ANCOM-BC (z tests, chi-squared global test).
  fb <- fit_of('sim_bc')
  expect_identical(attr(ancombc_dunnett(fb), 'method'), 'ancombc')
  expect_true(all(is.finite(ancombc_trend(fb)$w)))


  # -----------------------------------------------------------------------
  # ANCOMBC 2.8.1: deterministic post-hoc results, where its EM converges
  # -----------------------------------------------------------------------

  n_post <- 0L
  for (cs in names(ref$ancombc)) {
    a <- ref$ancombc[[cs]]
    if (is.null(a$global)) next
    fit <- fit_of(cs)
    g   <- ancombc_global(fit)
    expect_true(maxrel(g$w, a$global$W) < 1e-4, info = cs)
    expect_true(maxrel(g$p_value, a$global$p_val) < 1e-3, info = cs)
    expect_identical(g$reject, a$global$diff_abn, info = cs)
    if (is.null(a$pairwise)) next
    n_post <- n_post + 1L

    pw <- ancombc_pairwise(fit)
    m  <- length(unique(pw$comparison))
    cols <- grep('^lfc_', names(a$pairwise), value = TRUE)
    expect_identical(length(cols), m, info = cs)
    lf <- as.vector(t(as.matrix(a$pairwise[cols])))
    expect_true(maxabs(pw$log_fc, lf) < 1e-4, info = cs)
    qa <- as.vector(t(as.matrix(a$pairwise[sub('^lfc_', 'q_', cols)])))
    expect_true(maxrel(pw$p_adjusted, qa) < 1e-3, info = cs)
    expect_identical(pw$reject, as.vector(t(as.matrix(a$pairwise[sub('^lfc_', 'diff_', cols)]))), info = cs)

    tr <- ancombc_trend(fit)
    expect_true(maxabs(tr$w, a$trend$W) < 1e-4, info = cs)
  }
  expect_true(n_post >= 1L)


  # -----------------------------------------------------------------------
  # structural_zeros()
  # -----------------------------------------------------------------------

  for (z in ref$struc) {
    tb  <- ref$tables[[z$table]]
    res <- structural_zeros(tb$counts, z$groups, tb$data, neg_lb = z$neg_lb)
    info <- paste(z$table, z$neg_lb)
    # scikit-bio sorts the groups even for a pandas Categorical; ecodive keeps
    # a factor's level order (the `levels` table's `dose`).
    expect_identical(sort(colnames(res)), sort(colnames(z$zero)), info = info)
    if (is.factor(tb$data[[z$groups]]))
      expect_identical(colnames(res), levels(tb$data[[z$groups]]), info = info)
    expect_identical(unname(res[, colnames(z$zero)]), unname(z$zero), info = info)
    expect_identical(rownames(res), colnames(tb$counts), info = info)
  }

  tb <- ref$tables$struc
  x  <- tb$counts
  grp <- tb$data$grouping
  r0 <- structural_zeros(x, grp)
  expect_true(is.logical(r0) && is.matrix(r0))
  expect_identical(colnames(r0), c('placebo', 'treatment'))
  expect_identical(unname(which(r0)), c(5L, 9L))   # f4 in placebo, f2 in treatment
  expect_identical(structural_zeros(Matrix::Matrix(x, sparse = TRUE), grp), r0)
  expect_identical(structural_zeros(slam::as.simple_triplet_matrix(x), grp), r0)
  expect_identical(structural_zeros(t(x), grp, margin = 2L), r0)
  expect_identical(structural_zeros(x, 'grouping', tb$data), r0)
  fg <- factor(grp, levels = c('treatment', 'placebo', 'unused'))
  expect_identical(colnames(structural_zeros(x, fg)), c('treatment', 'placebo'))
  expect_identical(structural_zeros(x, fg)[, c('placebo', 'treatment')], r0)
  expect_identical(structural_zeros(x * 0.5, grp), r0)


  # -----------------------------------------------------------------------
  # The pseudocount sensitivity analysis
  # -----------------------------------------------------------------------

  tb  <- ref$tables$doc3
  f0  <- ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', sensitivity = FALSE, cpus = 1L)
  fs  <- ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', cpus = 1L)
  fits <- attr(fs, 'sensitivity')$fits
  same <- function (main, rerun) Reduce(`&`, lapply(fits, function (f) rerun(f)$reject == main))

  for (test in list(ancombc_global, ancombc_pairwise, ancombc_dunnett)) {
    r0 <- test(f0, cpus = 1L)
    rs <- test(fs, cpus = 1L)
    expect_identical(names(rs), c(names(r0), 'passed_ss', 'reject_robust'))
    expect_identical(rs$reject, r0$reject)
    expect_identical(rs$passed_ss, same(r0$reject, function (f) test(f, cpus = 1L)))
    expect_identical(rs$reject_robust, rs$reject & rs$passed_ss)
    expect_identical(attr(rs, 'sensitivity'), c(0.1, 0.5, 1))
  }
  expect_false(any(ancombc_pairwise(fs)$reject_robust))   # as the primary results

  # The trend test uses the global test's decisions, as ANCOMBC does.
  tr <- ancombc_trend(fs, cpus = 1L)
  g  <- ancombc_global(fs, cpus = 1L)
  expect_identical(tr$passed_ss, g$passed_ss)
  expect_identical(tr$reject_robust, tr$reject & g$passed_ss)

  # A grouping given at post-hoc time works on the stored fits too.
  fn <- ancombc2(tb$counts, ~ status + age, tb$data, cpus = 1L)
  expect_identical(ancombc_global(fn, groups = 'status'), ancombc_global(fs))

  grDevices::pdf(NULL)
  for (r in list(g, tr, ancombc_pairwise(fs), ancombc_dunnett(fs))) expect_silent(plot(r))
  grDevices::dev.off()


  # -----------------------------------------------------------------------
  # Argument checks
  # -----------------------------------------------------------------------

  fit <- fit_of('doc3')
  expect_error(ancombc_global(list()), 'result of ancombc')
  expect_error(ancombc_global(data.frame(a = 1)), 'result of ancombc')
  expect_error(ancombc_global(fit, groups = 'age'), 'categorical')
  expect_error(ancombc_global(fit, groups = 'nope'), 'term in the fit')
  expect_error(ancombc_global(fit, groups = 1), 'name')
  expect_error(ancombc_global(fit, alpha = 1), 'between 0 and 1')
  expect_error(ancombc_global(fit, p_adjust = 'nope'), 'p_adjust')
  expect_error(ancombc_dunnett(fit, n_bootstraps = 0), '`n_bootstraps`')
  expect_error(ancombc_dunnett(fit, seed = 1.5), '`seed`')
  expect_error(ancombc_trend(fit, patterns = list(diag(3))), 'one column per')
  expect_error(ancombc_trend(fit, patterns = list(a = diag(2)), nodes = 3), '`nodes`')
  expect_error(ancombc_trend(fit, patterns = 'up'), 'patterns')
  m2 <- ref$tables$doc3$data
  m2$two <- rep(c('u', 'v', 'u'), 5)
  f2 <- ancombc2(ref$tables$doc3$counts, ~ two + status, m2, sensitivity = FALSE, cpus = 1L)
  expect_error(ancombc_global(f2, groups = 'two'), 'three')
  expect_error(structural_zeros(x), 'groups')
  expect_error(structural_zeros(x, grp, neg_lb = NA), '`neg_lb`')
  xna <- x; xna[1, 1] <- NA
  expect_error(structural_zeros(xna, grp), 'NA')

#})
