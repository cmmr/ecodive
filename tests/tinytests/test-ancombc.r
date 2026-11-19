#test_that("ancombc", {

  # ancombc() and ancombc2(): scikit-bio 0.7.4's ancombc() and ancombc2()
  # (data-raw/expected/ancombc.py), including their biases, sampling
  # fractions and fitted quantities; ANCOMBC 2.8.1's on the cases where its
  # EM ends at the same estimate (data-raw/expected/ancombc.r); threading,
  # input types,
  # invariances and argument checks.

  ref <- readRDS(file.path('expected', 'ancombc.rds'))

  run <- function (cs, ...) {
    s  <- ref$skbio[[cs]]
    tb <- ref$tables[[s$table]]
    fn <- if (s$method == 'ancombc') ancombc else ancombc2
    args <- list(
      tb$counts, s$formula, tb$data, groups = s$groups,
      pseudocount = s$pseudocount, max_iters = s$max_iters, tol = s$tol,
      alpha = s$alpha, p_adjust = s$p_adjust, sensitivity = FALSE, cpus = 1L )
    if (s$method == 'ancombc2') args$var_quantile <- s$var_quantile
    do.call(fn, utils::modifyList(args, list(...)))
  }

  # The result's column as a features x covariates matrix.
  as_mat <- function (res, col) {
    cov <- unique(res$covariate)
    matrix(res[[col]], ncol = length(cov), byrow = TRUE, dimnames = list(NULL, cov))
  }

  # Largest |a - b|, and largest |a - b| / |b|, over finite pairs; and
  # whether the non-finite values agree.
  maxabs <- function (a, b) { ok <- is.finite(a) & is.finite(b); max(c(0, abs(a[ok] - b[ok]))) }
  maxrel <- function (a, b) {
    ok <- is.finite(a) & is.finite(b)
    max(c(0, abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300)))
  }
  same_na <- function (a, b) identical(is.finite(a), is.finite(b))
  plain   <- function (r) { r <- unclass(r)[seq_along(r)]; as.data.frame(r, stringsAsFactors = FALSE) }


  # -----------------------------------------------------------------------
  # scikit-bio: the table, the biases, the sampling fractions and the fit
  # -----------------------------------------------------------------------

  for (cs in names(ref$skbio)) {
    s   <- ref$skbio[[cs]]
    res <- run(cs)
    cov <- colnames(s$log_fc)        # patsy's order, R's names

    expect_inherits(res, c('ecodive_ancombc', 'data.frame'), info = cs)
    expect_identical(sort(unique(res$covariate)), sort(cov), info = cs)

    for (col in c('log_fc', 'se', 'p_value', 'p_adjusted', 'reject')) {
      m <- as_mat(res, col)[, cov, drop = FALSE]
      e <- s[[col]]
      info <- paste(cs, col)
      if (col == 'reject') { expect_identical(unname(m), unname(e), info = info); next }
      expect_true(same_na(m, e), info = info)
      if (col == 'log_fc')     expect_true(maxabs(m, e) < 1e-7, info = info)
      if (col == 'se')         expect_true(maxrel(m, e) < 1e-8, info = info)
      if (col %in% c('p_value', 'p_adjusted'))
                               expect_true(maxrel(m, e) < 1e-6, info = info)
    }

    b <- attr(res, 'bias')[cov, , drop = FALSE]
    expect_true(maxabs(b[, 'delta_em'],  s$bias[, 'delta_em'])  < 1e-7, info = cs)
    expect_true(maxabs(b[, 'delta_wls'], s$bias[, 'delta_wls']) < 1e-7, info = cs)
    expect_true(maxrel(b[, 'var_delta'], s$bias[, 'var_delta']) < 1e-6, info = cs)
    # The EM's iterations (one case stops one iteration apart, at the
    # tolerance's edge).
    expect_true(all(abs(attr(res, 'n_iters')[cov] - s$n_iters) <= 1L), info = cs)

    if (!is.null(s$sample_fractions))
      expect_true(maxabs(attr(res, 'sample_fractions'), s$sample_fractions) < 1e-7, info = cs)

    expect_equal(res$w, res$log_fc / res$se, info = cs)

    # The fitted quantities (scikit-bio's _beta_hat, _var_hat, _vcov_hat),
    # stored for a few cases.
    if (!is.null(s$coef)) {
      expect_true(maxabs(attr(res, 'coef')[, cov], s$coef) < 1e-7, info = cs)
      expect_true(maxrel(attr(res, 'var')[, cov], s$var) < 1e-8, info = cs)
      design <- attr(res, 'design')
      gi <- which(attr(design, 'assign') == match(s$groups, attr(design, 'term.labels')))
      expect_true(maxrel(attr(res, 'vcov')[, gi, gi, drop = FALSE], s$vcov) < 1e-8, info = cs)
    }
    if (!is.null(s$groups)) {
      est <- attr(res, 'estimable')
      expect_identical(if (is.null(est)) rep(TRUE, nrow(attr(res, 'coef'))) else unname(apply(est, 1L, all)), s$estimable, info = cs)
      if (s$method == 'ancombc2') expect_identical(unname(attr(res, 'dof')), s$dof, info = cs)
      else                        expect_null(attr(res, 'dof'), info = cs)
    }
  }

  # scikit-bio's documented values (ancombc() and ancombc2() docstrings).
  res <- run('doc1')
  expect_equal(round(res$log_fc, 3), c(-0.045, 0.126, -0.873, 1.241, -0.202, 0.561, 1.005, -1.723, 0.141, -0.69, -0.873, 1.48, 0.157, 0.049))
  expect_identical(res$covariate[1:2], c('(Intercept)', 'statussick'))
  res <- run('doc3')
  expect_identical(which(res$reject), c(6L, 7L, 23L))
  expect_equal(round(res$se[res$feature == 'F06'], 3), c(0.522, 0.283, 0.343, 0.016))


  # -----------------------------------------------------------------------
  # ANCOMBC 2.8.1 (whose estimation code is 2.14.0's), on the cases where
  # its bias EM ends at scikit-bio's estimate (data-raw/expected/ancombc.r
  # keeps those): it minimises with nloptr's Nelder-Mead, not SciPy's, so
  # the results agree to about 1e-5.
  # -----------------------------------------------------------------------

  expect_true(length(ref$ancombc) >= 5L)
  for (cs in names(ref$ancombc)) {
    a   <- ref$ancombc[[cs]]
    res <- run(cs)
    cov <- unique(res$covariate)
    am  <- function (m) { m <- m[, cov, drop = FALSE]; storage.mode(m) <- "double"; unname(m) }
    expect_true(maxabs(as_mat(res, 'log_fc'), am(a$log_fc)) < 1e-5, info = cs)
    expect_true(maxrel(as_mat(res, 'se'), am(a$se)) < 1e-5, info = cs)
    expect_true(maxrel(as_mat(res, 'p_value'), am(a$p_value)) < 1e-3, info = cs)
    expect_identical(unname(as_mat(res, 'reject')), am(a$reject) == 1, info = cs)
    if (!is.null(attr(res, 'sample_fractions')))
      expect_true(maxabs(unname(attr(res, 'sample_fractions')), unname(a$sample_fractions)) < 1e-4, info = cs)
  }


  # -----------------------------------------------------------------------
  # Threads, input types and orientation: identical results
  # -----------------------------------------------------------------------

  tb   <- ref$tables$sim
  x    <- tb$counts
  meta <- tb$data
  sp   <- Matrix::Matrix(x, sparse = TRUE)
  sl   <- slam::as.simple_triplet_matrix(x)

  for (fn in c('ancombc', 'ancombc2')) {
    f  <- get(fn)
    pc <- if (fn == 'ancombc') 0.5 else 0
    r1 <- f(x, ~ site + ph + smoker, meta, groups = 'site', pseudocount = pc, cpus = 1L)
    for (cpus in c(2L, 3L, 8L))
      expect_identical(f(x, ~ site + ph + smoker, meta, groups = 'site', pseudocount = pc, cpus = cpus), r1, info = paste(fn, cpus))
    expect_identical(f(sp, ~ site + ph + smoker, meta, groups = 'site', pseudocount = pc, cpus = 4L), r1, info = fn)
    expect_identical(f(sl, ~ site + ph + smoker, meta, groups = 'site', pseudocount = pc, cpus = 2L), r1, info = fn)
    expect_identical(f(t(x), ~ site + ph + smoker, meta, groups = 'site', margin = 2L, pseudocount = pc, cpus = 2L), r1, info = fn)
    expect_identical(f(Matrix::t(sp), ~ site + ph + smoker, meta, groups = 'site', margin = 2L, pseudocount = pc, cpus = 4L), r1, info = fn)
    no_f <- function (r) {   # the formula's environment differs
      attr(r, 'formula') <- NULL
      sens <- attr(r, 'sensitivity')
      if (!is.null(sens)) {
        sens$fits <- lapply(sens$fits, function (g) { attr(g, 'formula') <- NULL; g })
        attr(r, 'sensitivity') <- sens
      }
      r
    }
    expect_identical(no_f(f(x, 'site + ph + smoker', meta, groups = 'site', pseudocount = pc, cpus = 1L)), no_f(r1), info = fn)
    expect_identical(no_f(f(x, '~ site + ph + smoker', meta, groups = 'site', pseudocount = pc, cpus = 1L)), no_f(r1), info = fn)
  }

  # Many features, zeros as missing values: several feature blocks per
  # thread, and the iterative fit.
  set.seed(7)
  big <- matrix(rpois(30 * 300, rep(rexp(300, 0.05), each = 30)), 30, 300)
  big[sample(length(big), 800)] <- 0L
  bmeta <- data.frame(g = rep(c('a', 'b', 'c'), 10), z = rnorm(30))
  r1 <- ancombc2(big, ~ g + z, bmeta, groups = 'g', sensitivity = FALSE, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(ancombc2(big, ~ g + z, bmeta, groups = 'g', sensitivity = FALSE, cpus = cpus), r1, info = cpus)
  expect_identical(ancombc2(Matrix::Matrix(big, sparse = TRUE), ~ g + z, bmeta, groups = 'g', sensitivity = FALSE, cpus = 3L), r1)


  # -----------------------------------------------------------------------
  # Each sample's scale is absorbed by its sampling fraction
  # -----------------------------------------------------------------------

  x  <- ref$tables$prop$counts
  m  <- ref$tables$prop$data
  sc <- x * exp(seq(-2, 3, length.out = nrow(x)))
  for (fn in c('ancombc', 'ancombc2')) {
    a <- get(fn)(x,  ~ grp, m, sensitivity = FALSE, cpus = 1L)
    b <- get(fn)(sc, ~ grp, m, sensitivity = FALSE, cpus = 1L)
    k <- a$covariate != '(Intercept)'
    # Rescaling shifts every slope by one constant, which the bias EM
    # absorbs to its tolerance (1e-5).
    expect_equal(b$log_fc[k], a$log_fc[k], tolerance = 1e-4, info = fn)
    expect_equal(b$se, a$se, tolerance = 1e-6, info = fn)
  }
  # With zeros as missing values, the sample effects are fitted only to
  # ANCOMBC's tolerance (0.01), so only a common scale is absorbed exactly
  # (scikit-bio's results move as much for per-sample scales).
  x  <- ref$tables$sparse$counts
  m  <- ref$tables$sparse$data
  a  <- ancombc2(x, ~ trt + age, m, sensitivity = FALSE, cpus = 1L)
  b  <- ancombc2(x * 7, ~ trt + age, m, sensitivity = FALSE, cpus = 1L)
  k  <- a$covariate != '(Intercept)'
  expect_equal(b$log_fc[k], a$log_fc[k], tolerance = 1e-8)
  expect_equal(b$p_value, a$p_value, tolerance = 1e-8)


  # -----------------------------------------------------------------------
  # The result's layout, covariates and reference levels
  # -----------------------------------------------------------------------

  x    <- ref$tables$doc2$counts + 1L
  meta <- ref$tables$doc2$data
  res  <- ancombc(x, ~ status + age, meta, groups = 'status', cpus = 1L)

  expect_identical(names(res), c('feature', 'covariate', 'log_fc', 'se', 'w', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(res$feature, rep(colnames(x), each = 4L))
  expect_identical(res$covariate, rep(c('(Intercept)', 'statusmoderate', 'statussevere', 'age'), 8L))
  expect_equal(res$w, res$log_fc / res$se)
  expect_identical(attr(res, 'method'), 'ancombc')
  expect_identical(attr(res, 'groups'), 'status')
  expect_identical(attr(res, 'p_adjust'), 'holm')
  expect_identical(attr(res, 'pseudocount'), 0)
  expect_null(attr(res, 'var_quantile'))
  expect_null(attr(res, 'sample_fractions'))
  expect_identical(dim(attr(res, 'coef')), c(8L, 4L))
  expect_identical(dim(attr(res, 'vcov')), c(8L, 4L, 4L))
  expect_identical(dimnames(attr(res, 'coef')), list(colnames(x), c('(Intercept)', 'statusmoderate', 'statussevere', 'age')))
  for (p in 1:4) expect_identical(attr(res, 'vcov')[, p, p], attr(res, 'var')[, p])
  expect_true(isSymmetric(attr(res, 'vcov')[3, , ]))
  expect_identical(unname(attr(res, 'coef')[, 2]), unname(as_mat(res, 'log_fc')[, 2]))
  expect_identical(unname(sqrt(attr(res, 'var')[, 2])), unname(as_mat(res, 'se')[, 2]))
  expect_identical(dim(attr(res, 'design')), c(15L, 4L))
  expect_identical(attr(attr(res, 'design'), 'assign'), c(0L, 1L, 1L, 2L))

  # Covariate order follows the formula.
  r2 <- ancombc(x, ~ age + status, meta, cpus = 1L)
  expect_identical(unique(r2$covariate), c('(Intercept)', 'age', 'statusmoderate', 'statussevere'))
  expect_equal(r2$log_fc[r2$covariate == 'age'], res$log_fc[res$covariate == 'age'], tolerance = 1e-9)

  # A factor keeps its levels' order; characters are sorted (C locale).
  m2 <- meta
  m2$status <- factor(m2$status, levels = c('severe', 'mild', 'moderate'))
  r3 <- ancombc(x, ~ status + age, m2, cpus = 1L)
  expect_identical(unique(r3$covariate), c('(Intercept)', 'statusmild', 'statusmoderate', 'age'))
  m2$status <- factor(m2$status, levels = c('mild', 'moderate', 'severe', 'unused'), ordered = TRUE)
  expect_identical(plain(ancombc(x, ~ status + age, m2, cpus = 1L)), plain(res), info = 'ordered, unused level')
  m2$status <- c('b', 'B', 'a')[as.integer(factor(meta$status))]
  expect_identical(unique(ancombc(x, ~ status, m2, cpus = 1L)$covariate), c('(Intercept)', 'statusa', 'statusb'))
  m2$flag <- rep(c(TRUE, FALSE, TRUE), 5)
  expect_identical(unique(ancombc(x, ~ flag, m2, cpus = 1L)$covariate), c('(Intercept)', 'flagTRUE'))

  # data matched by row name; a vector as `data`; unnamed counts.
  shuffled <- meta[15:1, ]
  expect_identical(plain(ancombc(x, ~ status + age, shuffled, groups = 'status', cpus = 1L)), plain(res))
  r4 <- ancombc(x, ~ group, setNames(meta$status, rownames(x)), cpus = 1L)
  expect_identical(r4$log_fc[r4$covariate != 'age'], ancombc(x, ~ status, meta, cpus = 1L)$log_fc)
  expect_identical(ancombc(unname(x), ~ status, meta, cpus = 1L)$feature, rep(as.character(1:8), each = 3L))

  # p_adjust and alpha.
  r5 <- ancombc(x, ~ status + age, meta, p_adjust = 'BH', alpha = 0.2, cpus = 1L)
  for (cv in unique(r5$covariate)) {
    k <- r5$covariate == cv
    expect_identical(r5$p_adjusted[k], stats::p.adjust(r5$p_value[k], 'BH'))
  }
  expect_identical(r5$reject, r5$p_adjusted <= 0.2)
  expect_identical(attr(ancombc(x, ~ status, meta, p_adjust = 'fdr', cpus = 1L), 'p_adjust'), 'fdr')
  expect_identical(ancombc(x, ~ status, meta, p_adjust = 'none', cpus = 1L)$p_adjusted,
                   ancombc(x, ~ status, meta, p_adjust = 'none', cpus = 1L)$p_value)


  # -----------------------------------------------------------------------
  # Zeros
  # -----------------------------------------------------------------------

  x0 <- ref$tables$doc2$counts
  expect_warning(r6 <- ancombc(x0, ~ status, meta, cpus = 1L), 'pseudocount = 0.5')
  expect_identical(attr(r6, 'pseudocount'), 0.5)
  expect_identical(plain(r6), plain(ancombc(x0 + 0.5, ~ status, meta, cpus = 1L)))
  expect_error(ancombc(x0, ~ status, meta, pseudocount = 0), 'zeros')
  expect_error(ancombc(x0, ~ status, meta, pseudocount = -1), 'non-negative')

  # ancombc2(): zeros are missing values by default.
  r7 <- ancombc2(x0, ~ status, meta, cpus = 1L)
  expect_identical(attr(r7, 'pseudocount'), 0)
  expect_identical(attr(r7, 'method'), 'ancombc2')
  expect_identical(attr(r7, 'var_quantile'), 0.05)
  expect_identical(length(attr(r7, 'sample_fractions')), 15L)
  expect_identical(names(attr(r7, 'sample_fractions')), rownames(x0))
  expect_identical(length(attr(r7, 'dof')), 8L)
  expect_error(ancombc2(x0, ~ status, meta, pseudocount = NULL), 'non-negative')
  x1 <- x0; x1[1, ] <- 0L
  expect_error(ancombc2(x1, ~ status, meta), 'only zeros')
  expect_true(all(is.finite(ancombc2(x1, ~ status, meta, pseudocount = 1, cpus = 1L)$log_fc)))

  # Non-estimable coefficients: NaN, p-value 1, not rejected.
  r8  <- run('sparse_bc2')
  est <- attr(r8, 'estimable')
  bad <- rownames(est)[!apply(est, 1L, all)]
  expect_identical(length(bad), 2L)
  k <- r8$feature %in% bad
  expect_true(all(is.nan(r8$log_fc[k]) & is.nan(r8$se[k]) & is.nan(r8$w[k])))
  expect_true(all(r8$p_value[k] == 1) && !any(r8$reject[k]))
  expect_true(all(is.finite(r8$log_fc[!k])))
  expect_true(all(is.finite(attr(r8, 'coef')[bad, ])))


  # -----------------------------------------------------------------------
  # The pseudocount sensitivity analysis
  # -----------------------------------------------------------------------

  tb  <- ref$tables$doc3
  f0  <- ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', sensitivity = FALSE, cpus = 1L)
  fs  <- ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', cpus = 1L)

  # On by default for ancombc2() (as ANCOMBC), off for ancombc(); the main
  # columns are the same either way, and two are added.
  expect_identical(formals(ancombc2)$sensitivity, TRUE)
  expect_identical(formals(ancombc)$sensitivity, FALSE)
  expect_identical(fs, ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', sensitivity = TRUE, cpus = 1L))
  expect_null(attr(ancombc(tb$counts + 1L, ~ status, tb$data, cpus = 1L), 'sensitivity'))
  expect_identical(names(f0), c('feature', 'covariate', 'log_fc', 'se', 'w', 'p_value', 'p_adjusted', 'reject'))
  expect_null(attr(f0, 'sensitivity'))
  expect_identical(names(fs), c(names(f0), 'passed_ss', 'reject_robust'))
  expect_identical(plain(fs)[names(f0)], plain(f0))
  expect_inherits(fs, c('ecodive_ancombc', 'data.frame'))

  # scikit-bio's docstring: none of the three significant results is robust.
  expect_identical(sum(fs$reject), 3L)
  expect_identical(sum(fs$reject_robust), 0L)

  # The rule: the decision is the same in the fits at 0.1, 0.5 and 1.
  man <- lapply(c(0.1, 0.5, 1), function (pc)
    ancombc2(tb$counts, ~ status + age, tb$data, groups = 'status', pseudocount = pc, sensitivity = FALSE, cpus = 1L) )
  passed <- Reduce(`&`, lapply(man, function (m) m$reject == f0$reject))
  expect_identical(fs$passed_ss, passed)
  expect_identical(fs$reject_robust, f0$reject & passed)
  sens <- attr(fs, 'sensitivity')
  expect_identical(sens$pseudocounts, c(0.1, 0.5, 1))
  no_f <- function (r) { attr(r, 'formula') <- NULL; r }   # its environment differs
  for (i in 1:3) expect_identical(no_f(sens$fits[[i]]), no_f(man[[i]]))

  # Custom pseudocounts, ancombc() too, threads.
  x1 <- ref$tables$sim$counts
  m1 <- ref$tables$sim$data
  f1 <- ancombc(x1, ~ site + ph, m1, pseudocount = 0.5, sensitivity = c(1, 2), cpus = 1L)
  m  <- lapply(c(1, 2), function (pc) ancombc(x1, ~ site + ph, m1, pseudocount = pc, cpus = 1L))
  expect_identical(f1$passed_ss, m[[1]]$reject == f1$reject & m[[2]]$reject == f1$reject)
  expect_identical(attr(f1, 'sensitivity')$pseudocounts, c(1, 2))
  expect_identical(ancombc(x1, ~ site + ph, m1, pseudocount = 0.5, sensitivity = c(1, 2), cpus = 4L), f1)
  expect_identical(attr(ancombc(x1, ~ site + ph, m1, pseudocount = 0.5, sensitivity = FALSE, cpus = 1L), 'sensitivity'), NULL)

  for (v in list('yes', c(0, 1), -1, NA, numeric(0), Inf))
    expect_error(ancombc2(tb$counts, ~ status, tb$data, sensitivity = v), '`sensitivity`')

  # The default pseudocounts assume counts: a warning for other values,
  # none when pseudocounts are given or the analysis is off.
  px <- ref$tables$prop$counts
  pm <- ref$tables$prop$data
  expect_warning(ancombc2(px, ~ grp, pm, cpus = 1L), 'not whole numbers')
  expect_silent(ancombc2(px, ~ grp, pm, sensitivity = c(1e-4, 1e-3), cpus = 1L))
  expect_silent(ancombc2(px, ~ grp, pm, sensitivity = FALSE, cpus = 1L))

  # Plots flag the robust decisions.
  grDevices::pdf(NULL)
  expect_silent(plot(fs))
  expect_silent(plot(f1, 'siteskin'))
  grDevices::dev.off()


  # -----------------------------------------------------------------------
  # Argument checks
  # -----------------------------------------------------------------------

  expect_error(ancombc(x), '`formula` is required')
  expect_error(ancombc(x, ~ status), '`data` is required')
  expect_error(ancombc(x, ~ status + nope, meta), "Not found in `data`: 'nope'")
  expect_error(ancombc(x, 3, meta), 'must be a formula')
  expect_error(ancombc(x, '~', meta), 'formula')
  m3 <- meta; m3$age[2] <- NA
  expect_error(ancombc(x, ~ age, m3), 'NA')
  m3$age[2] <- Inf
  expect_error(ancombc(x, ~ age, m3), 'Infinite')
  m3$one <- 'a'
  expect_error(ancombc(x, ~ one, m3), 'single value')
  expect_error(ancombc(x, ~ status, meta[-1, ]), 'no row')
  expect_error(ancombc(x[c(1, 2, 6, 11), ], ~ status + age, meta), 'more samples than design columns')
  m3 <- meta; m3$age2 <- m3$age * 2
  expect_warning(ancombc(x, ~ age + age2, m3, cpus = 1L), 'rank deficient')
  expect_error(ancombc(x, ~ status + age, meta, groups = 'age'), 'categorical')
  expect_error(ancombc(x, ~ status + age, meta, groups = 'nope'), 'column of `data`')
  m3 <- meta; m3$two <- rep(c('u', 'v', 'u'), 5)
  expect_error(ancombc(x, ~ status + age, m3, groups = 'two'), 'term in `formula`')
  expect_error(ancombc(x, ~ two + age, m3, groups = 'two'), 'three')
  expect_error(ancombc(x, ~ status, meta, groups = 1), 'name')
  for (arg in c('alpha')) for (v in list(0, 1, NA, 'a', c(0.1, 0.2)))
    expect_error(do.call(ancombc, setNames(list(x, ~ status, meta, v), c('counts', 'formula', 'data', arg))), 'between 0 and 1')
  expect_error(ancombc(x, ~ status, meta, tol = 0), '`tol`')
  expect_error(ancombc(x, ~ status, meta, max_iters = 0), '`max_iters`')
  expect_error(ancombc(x, ~ status, meta, max_iters = 1.5), '`max_iters`')
  expect_error(ancombc2(x, ~ status, meta, var_quantile = 1.5), '`var_quantile`')
  expect_error(ancombc2(x, ~ status, meta, var_quantile = NA), '`var_quantile`')
  expect_error(ancombc(x, ~ status, meta, p_adjust = 'nope'), 'p_adjust')
  expect_error(ancombc(x[, 1, drop = FALSE], ~ status, meta), 'at least 2 features')
  expect_error(ancombc(-x, ~ status, meta), 'negative')
  xna <- x; xna[1, 1] <- NA
  expect_error(ancombc(xna, ~ status, meta), 'NA')

#})
