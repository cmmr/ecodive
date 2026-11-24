#test_that("cross_assoc", {

  # cross_assoc(): correlations between the features of two tables, against
  # SciPy (BSD-3) and the output-only references stored by
  # data-raw/expected/cross_assoc.r (cor.test(), Hmisc, psych, ppcor,
  # microeco), and against base-R computations.

  f   <- readRDS(file.path('expected', 'cross_assoc.rds'))
  cnt <- f$counts
  env <- f$env
  rel <- cnt / rowSums(cnt)
  q   <- function (expr) suppressWarnings(expr)


  #--------------------------------------------------
  # SciPy and cor.test(exact = FALSE): values at 1e-12, p-values at 1e-10.
  #--------------------------------------------------
  for (m in c('spearman', 'pearson', 'kendall')) {
    res <- q(cross_assoc(cnt, env, method = m, cpus = 2L))
    r   <- as.matrix(res)
    pv  <- as.matrix(res, 'p_value')
    expect_identical(is.na(r), is.na(f$scipy[[m]]), info = m)
    expect_equal(r, f$scipy[[m]], tolerance = 1e-12, info = m)
    expect_equal(pv, f$scipy[[paste0(m, '_p')]], tolerance = 1e-10, info = m)
    expect_equal(r, f$cortest[[m]], tolerance = 1e-12, info = m)
    expect_equal(pv, f$cortest[[paste0(m, '_p')]], tolerance = 1e-10, info = m)
    n <- as.matrix(res, 'n')
    storage.mode(n) <- 'double'
    expect_equal(n, f$scipy$n, info = m)
  }

  # The table itself.
  res <- q(cross_assoc(cnt, env))
  expect_identical(class(res), c('ecodive_cross_assoc', 'data.frame'))
  expect_identical(names(res),
    c('x', 'y', 'association', 'statistic', 'p_value', 'p_adjusted', 'reject', 'n'))
  expect_identical(nrow(res), 25L * 6L)
  expect_identical(res$x, rep(colnames(cnt), each = 6))
  expect_identical(res$y, rep(names(env), times = 25))
  expect_identical(attr(res, 'method'), 'spearman')
  expect_identical(attr(res, 'alpha'), 0.05)
  expect_identical(attr(res, 'p_adjust'), 'BH')
  expect_identical(attr(res, 'n_samples'), 40L)
  expect_identical(attr(res, 'x_features'), colnames(cnt))
  expect_identical(attr(res, 'y_features'), names(env))
  expect_null(attr(res, 'n_perm'))
  expect_null(attr(res, 'given'))
  expect_equal(res$p_adjusted, p.adjust(res$p_value, 'BH'))
  expect_identical(res$reject, res$p_adjusted <= 0.05)
  df <- res$n - 2
  expect_equal(res$statistic, res$association * sqrt(df / (1 - res$association^2)))
  expect_true(res$reject[res$x == 'Taxon01' & res$y == 'ph'])
  expect_identical(sum(res$reject, na.rm = TRUE), 1L)
  expect_true(all(is.na(res$association[res$y == 'const' | res$x == 'Taxon25'])))

  # Constant columns warn, as cor() does.
  expect_warning(cross_assoc(cnt, env), 'feature of `x` and 1 column of `y`')

  # No `n` column when every pair has the same samples.
  r2 <- cross_assoc(cnt[, -25], env[, c('ph', 'od')])
  expect_false('n' %in% names(r2))


  #--------------------------------------------------
  # Hmisc (p-values only where Hmisc reports them precisely), psych,
  # microeco.
  #--------------------------------------------------
  h   <- f$hmisc
  sub <- function (m) m[rownames(h$r), colnames(h$r)]
  expect_equal(sub(as.matrix(res)), h$r, tolerance = 1e-12)
  hn <- sub(as.matrix(res, 'n')); storage.mode(hn) <- 'double'; storage.mode(h$n) <- 'double'
  expect_equal(hn, h$n)
  big <- h$P > 1e-8
  expect_equal(sub(as.matrix(res, 'p_value'))[big], h$P[big], tolerance = 1e-10)
  expect_true(all(abs(sub(as.matrix(res, 'p_value')) - h$P) < 1e-14))

  ps <- f$psych
  r2 <- cross_assoc(cnt[, -25], env[, names(env) != 'const'])
  expect_equal(as.matrix(r2), ps$r, tolerance = 1e-12)
  expect_equal(as.matrix(r2, 'p_value'), ps$p, tolerance = 1e-10)
  expect_equal(as.matrix(r2, 'p_adjusted'), ps$p.adj, tolerance = 1e-10)

  me <- f$microeco
  r3 <- cross_assoc(cnt[, -25], env[, c('ph', 'temp', 'od')])
  k  <- match(paste(me$Taxa, me$Env), paste(r3$x, r3$y))
  expect_false(anyNA(k))
  expect_equal(r3$association[k], me$Correlation, tolerance = 1e-12)
  # Its cor.test() is exact (AS 89) for pairs without ties.
  tied <- vapply(seq_len(nrow(me)), function (i)
    anyDuplicated(rel[, me$Taxa[[i]]]) > 0 || anyDuplicated(env[[me$Env[[i]]]]) > 0, logical(1))
  expect_true(sum(!tied) >= 1 && sum(tied) > 60)
  expect_equal(r3$p_value[k][tied], me$Pvalue[tied], tolerance = 1e-10)
  expect_equal(p.adjust(me$Pvalue, 'fdr'), me$AdjPvalue)


  #--------------------------------------------------
  # Partial correlations: ppcor, and base R.
  #--------------------------------------------------
  vars <- c('ph', 'temp', 'od')
  for (m in c('spearman', 'pearson')) {
    rp <- q(cross_assoc(cnt, env[vars], given = f$z, method = m))
    pp <- f$ppcor[[m]]
    expect_equal(as.matrix(rp)[1:24, ], pp$estimate, tolerance = 1e-12, info = m)
    expect_equal(as.matrix(rp, 'p_value')[1:24, ], pp$p.value, tolerance = 1e-10, info = m)
    expect_equal(as.matrix(rp, 'statistic')[1:24, ], pp$statistic, tolerance = 1e-10, info = m)
    expect_identical(attr(rp, 'given'), c('age', 'sex'))
  }

  # A formula over `data`, and names of its columns, are the same.
  dat <- cbind(env, f$z)
  a <- q(cross_assoc(cnt, vars, data = dat, given = ~ age + sex))
  b <- q(cross_assoc(cnt, vars, data = dat, given = c('age', 'sex')))
  d <- q(cross_assoc(cnt, env[vars], given = f$z))
  expect_identical(a, b)
  expect_identical(a, d)

  # With NAs in `y` and a factor covariate, against base R.
  z <- data.frame(age = f$z$age, site = factor(rep(c('a', 'b', 'c'), length.out = 40)))
  rp <- q(cross_assoc(cnt, env[c('salt', 'depth')], data = z, given = ~ age + site, method = 'pearson'))
  ref <- sapply(c('salt', 'depth'), function (v) {
    ok <- !is.na(env[[v]])
    ry <- resid(lm(env[[v]][ok] ~ z$age[ok] + z$site[ok]))
    sapply(1:24, function (i) cor(resid(lm(rel[ok, i] ~ z$age[ok] + z$site[ok])), ry))
  })
  expect_equal(unname(as.matrix(rp)[1:24, ]), unname(ref), tolerance = 1e-10)
  expect_identical(unique(rp$n[rp$y == 'salt']), 37L)
  dfp <- rp$n - 2 - 3
  expect_equal(rp$p_value, 2 * pt(-abs(rp$statistic), dfp))

  # NA in `given`: those samples are left out, with a warning.
  z2 <- f$z
  z2$age[1:2] <- NA
  expect_warning(rp <- cross_assoc(cnt[, 1:5], env['ph'], given = z2), '2 samples have NA')
  expect_identical(attr(rp, 'n_samples'), 38L)
  expect_equal(rp$association,
    q(cross_assoc(cnt[-(1:2), 1:5], env[-(1:2), 'ph', drop = FALSE], given = f$z[-(1:2), ]))$association)


  #--------------------------------------------------
  # Sparse and dense input; margins; threads; seeds.
  #--------------------------------------------------
  base <- q(cross_assoc(cnt, env, cpus = 1L))
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(cnt, sparse = TRUE)
    expect_equal(q(cross_assoc(sp, env)), base, tolerance = 1e-14)
    expect_equal(q(cross_assoc(Matrix::t(sp), env, margin = 2L)), base, tolerance = 1e-14)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_equal(q(cross_assoc(slam::as.simple_triplet_matrix(cnt), env)), base, tolerance = 1e-14)
  expect_identical(q(cross_assoc(t(cnt), env, margin = 2L)), base)
  expect_identical(q(cross_assoc(as.data.frame(cnt), env)), base)

  for (m in c('spearman', 'pearson', 'kendall')) {
    one <- q(cross_assoc(cnt, env, method = m, permutations = 49L, cpus = 1L))
    for (cp in c(2L, 3L, 4L))
      expect_identical(q(cross_assoc(cnt, env, method = m, permutations = 49L, cpus = cp)), one, info = m)
  }
  one <- q(cross_assoc(cnt, env[vars], given = f$z, permutations = 99L, cpus = 1L))
  expect_identical(q(cross_assoc(cnt, env[vars], given = f$z, permutations = 99L, cpus = 4L)), one)
  expect_identical(q(cross_assoc(cnt, env[vars], given = f$z, permutations = 99L, cpus = 1L)), one)
  expect_false(identical(
    q(cross_assoc(cnt, env[vars], permutations = 99L, seed = 1L))$p_value,
    q(cross_assoc(cnt, env[vars], permutations = 99L, seed = 2L))$p_value ))

  # Dense input with zeros equals cor(method = 'spearman') of the counts.
  expect_equal(
    as.matrix(q(cross_assoc(cnt, env['od'], norm = 'none'))),
    q(cor(cnt, env['od'], method = 'spearman')), tolerance = 1e-12 )


  #--------------------------------------------------
  # Permutations: a supplied matrix against a base-R loop.
  #--------------------------------------------------
  set.seed(3)
  P <- t(replicate(99, sample(40)))
  y <- env[vars]
  for (m in c('spearman', 'pearson', 'kendall')) {
    rp  <- cross_assoc(rel[, 1:8], y, method = m, norm = 'none', permutations = P)
    obs <- cor(rel[, 1:8], y, method = m)
    hit <- Reduce(`+`, lapply(1:99, function (b)
      abs(cor(rel[, 1:8], y[P[b, ], ], method = m)) >= abs(obs) - sqrt(.Machine$double.eps) ))
    expect_identical(as.matrix(rp, 'p_value'), (hit + 1) / 100, info = m)
    expect_identical(attr(rp, 'n_perm'), 99L)
  }
  # With `given`: the residuals of y are permuted (ter Braak).
  rp  <- cross_assoc(rel[, 1:8], y, norm = 'none', given = f$z, permutations = P)
  Z   <- cbind(rank(f$z$age), f$z$sex)
  rx  <- apply(rel[, 1:8], 2, function (v) resid(lm(rank(v) ~ Z)))
  ry  <- apply(y, 2, function (v) resid(lm(rank(v) ~ Z)))
  obs <- cor(rx, ry)
  hit <- Reduce(`+`, lapply(1:99, function (b)
    abs(cor(rx, ry[P[b, ], ])) >= abs(obs) - sqrt(.Machine$double.eps) ))
  expect_identical(as.matrix(rp, 'p_value'), (hit + 1) / 100)

  # Random permutations: p in [1/(B+1), 1]; strata keep their samples.
  rp <- q(cross_assoc(cnt, env, permutations = 199L))
  ok <- !is.na(rp$association)
  expect_true(all(rp$p_value[ok] >= 1 / 200 & rp$p_value[ok] <= 1))
  expect_true(all(is.na(rp$p_value[!ok])))
  expect_identical(rp$p_value[rp$x == 'Taxon01' & rp$y == 'ph'], 1 / 200)
  st <- rep(1:2, each = 20)
  rs <- q(cross_assoc(cnt, env, permutations = 19L, strata = st))
  expect_true(all(rs$p_value >= 1 / 20, na.rm = TRUE))
  expect_error(cross_assoc(cnt, env, permutations = P), 'without NAs')

  # With NAs: each set of samples (y's columns with the same NAs) gets its
  # own permutations of those samples, from `seed` for the first set and
  # pair_seed(seed, k) for set k + 1.
  y  <- env[c('ph', 'salt', 'od', 'depth')]
  rp <- cross_assoc(rel[, 1:6], y, norm = 'none', permutations = 49L, seed = 5L)
  sets <- list(c('ph', 'od'), 'salt', 'depth')
  seeds <- c(5L, pair_seed(5L, 1L), pair_seed(5L, 2L))
  for (k in 1:3) {
    v    <- sets[[k]]
    rows <- which(!is.na(y[[v[[1]]]]))
    pk   <- .Call(C_perm_indices, length(rows), 49L, NULL, seeds[[k]])
    yy   <- as.matrix(y[rows, v, drop = FALSE])
    obs  <- cor(rel[rows, 1:6], yy, method = 'spearman')
    hit  <- Reduce(`+`, lapply(1:49, function (b)
      abs(cor(rel[rows, 1:6], yy[pk[, b], , drop = FALSE], method = 'spearman')) >=
        abs(obs) - sqrt(.Machine$double.eps) ))
    expect_identical(as.matrix(rp, 'p_value')[, v, drop = FALSE], (hit + 1) / 50, info = k)
  }


  #--------------------------------------------------
  # Invariances and edge cases.
  #--------------------------------------------------
  # Monotone transforms of y's columns change neither Spearman nor Kendall.
  y2 <- data.frame(ph = exp(env$ph), od = env$od^3, row.names = rownames(env))
  for (m in c('spearman', 'kendall'))
    expect_equal(
      cross_assoc(cnt[, 1:10], y2, method = m)$association,
      cross_assoc(cnt[, 1:10], env[c('ph', 'od')], method = m)$association,
      tolerance = 1e-14, info = m )

  # Swapping x and y transposes the result.
  a <- cross_assoc(rel[, 1:5], env[c('ph', 'od')], norm = 'none')
  b <- cross_assoc(env[c('ph', 'od')], rel[, 1:5], norm = 'none')
  expect_equal(as.matrix(b), t(as.matrix(a)), tolerance = 1e-14)

  # n = 3.
  r3 <- cross_assoc(cnt[1:3, 1:4], env[1:3, c('ph', 'od')], method = 'pearson')
  expect_true(all(r3$n == 3L) || is.null(r3$n))
  expect_equal(as.matrix(r3),
    cor(cnt[1:3, 1:4] / rowSums(cnt[1:3, 1:4]), env[1:3, c('ph', 'od')]), tolerance = 1e-12)
  expect_equal(r3$p_value, 2 * pt(-abs(r3$statistic), 1))

  # Too few samples: NA.
  y3 <- env['ph']; y3$ph[3:40] <- NA
  r4 <- q(cross_assoc(cnt[, 1:3], cbind(y3, od = env$od)))
  expect_true(all(is.na(r4$association[r4$y == 'ph'])))
  expect_false(anyNA(r4$association[r4$y == 'od']))
  expect_identical(cross_assoc(cnt[, 1:3], env['od'], min_samples = 41L)$association, rep(NA_real_, 3))

  # An all-NA column, and |r| = 1 (p = 0).
  r5 <- cross_assoc(rel[, 1:2], data.frame(a = NA_real_, b = rel[, 1]), norm = 'none')
  expect_true(all(is.na(r5$association[r5$y == 'a'])))
  expect_identical(r5$association[[2]], 1)
  expect_identical(r5$p_value[[2]], 0)
  expect_identical(r5$statistic[[2]], Inf)

  # Logical columns are 0/1; a plain vector; columns of `data`.
  yl <- data.frame(hi = env$ph > 0)
  expect_equal(cross_assoc(cnt[, 1:3], yl)$association,
               cross_assoc(cnt[, 1:3], data.frame(hi = as.numeric(env$ph > 0)))$association)
  expect_identical(cross_assoc(cnt[, 1:3], setNames(env$ph, rownames(env)))$y, rep('y', 3))
  expect_identical(cross_assoc(cnt[, 1:3], 'ph', data = env)$association,
                   cross_assoc(cnt[, 1:3], env['ph'])$association)

  # Unnamed features and variables: their positions, as text.
  un <- cross_assoc(unname(cnt[, 1:2]), unname(as.matrix(env[c('ph', 'od')])))
  expect_identical(un$x, c('1', '1', '2', '2'))
  expect_identical(un$y, c('1', '2', '1', '2'))
  expect_identical(attr(un, 'y_features'), c('1', '2'))

  # Matching by name.
  shuffled <- env[40:1, c('ph', 'od')]
  expect_identical(cross_assoc(cnt[, 1:3], shuffled), cross_assoc(cnt[, 1:3], env[c('ph', 'od')]))
  bad <- env[c('ph', 'od')]; rownames(bad)[1] <- 'nope'
  expect_error(cross_assoc(cnt[, 1:3], bad), "1 sample in `x` has no row in `y`: 'S01'")

  # Norms: clr equals cor() of normalize()'s clr.
  expect_equal(
    as.matrix(cross_assoc(cnt[, 1:24], env['od'], norm = 'clr', pseudocount = 1, method = 'pearson')),
    cor(normalize(cnt[, 1:24], 'clr', pseudocount = 1), env['od']), tolerance = 1e-12 )
  y_cnt <- cnt[, 1:3] + 1L
  expect_equal(
    as.matrix(cross_assoc(rel[, 4:6], y_cnt, norm = 'none', y_norm = 'clr', method = 'pearson')),
    cor(rel[, 4:6], normalize(y_cnt, 'clr')), tolerance = 1e-12 )

  # as.matrix() after subsetting rows (attributes dropped).
  s <- res[res$y %in% c('ph', 'od'), ]
  expect_equal(as.matrix(s)[, c('ph', 'od')], as.matrix(res)[, c('ph', 'od')])


  #--------------------------------------------------
  # Errors.
  #--------------------------------------------------
  expect_error(cross_assoc(cnt, data.frame(site = letters[1:40])),
               "`y` column 'site' is not numeric; use da_test\\(\\) for groups")
  expect_error(cross_assoc(cnt, data.frame(site = factor(1:40))), 'not numeric')
  expect_error(cross_assoc(cnt, env, method = 'kendall', given = f$z), 'not available for Kendall')
  expect_error(cross_assoc(cnt, env, method = 'rho_p'), "'spearman', 'pearson', or 'kendall'")
  expect_error(cross_assoc(cnt, env, min_samples = 2), 'at least 3')
  expect_error(cross_assoc(cnt, env, p_adjust = 'nope'), 'p_adjust')
  expect_error(cross_assoc(cnt, env, alpha = 2), 'alpha')
  expect_error(cross_assoc(cnt, 'ph'), '`data` is not given')
  expect_error(cross_assoc(cnt, 'nope', data = env), "'nope', not column of `data`")
  expect_error(cross_assoc(cnt, env[1:30, ]), '10 samples in `x` have no row in `y`')
  expect_error(cross_assoc(unname(cnt), unname(as.matrix(env[1:30, ]))), '`y` has 30 rows, but `x` has 40 samples')
  expect_error(cross_assoc(cnt, env, given = ~ age), '`data` is not given')
  expect_error(cross_assoc(cnt, env, y_norm = 'clr'), 'must not contain NAs')
  x_na <- cnt; x_na[1, 1] <- NA
  expect_error(cross_assoc(x_na, env), '`x` must not contain NA')
  expect_error(cross_assoc(cnt, env, data = 1), '`data` must be a data frame')
  expect_error(plot(res, stars = 2), '`stars`')


  #--------------------------------------------------
  # Printing and plotting.
  #--------------------------------------------------
  out <- capture.output(print(res))
  expect_true(any(grepl("Spearman Correlation: 25 x 6 features, 1 pair with BH-adjusted p <= 0.05", out, fixed = TRUE)))
  rp <- cross_assoc(cnt[, 1:5], env[vars], given = f$z, method = 'pearson')
  expect_true(any(grepl('Partial Pearson Correlation', capture.output(print(rp)), fixed = TRUE)))

  expect_identical(cross_assoc_name(res), 'Spearman Correlation')
  expect_identical(cross_assoc_name(rp), 'Partial Pearson Correlation')
  expect_identical(cross_assoc_name(q(cross_assoc(cnt, env, method = 'kendall'))), "Kendall's tau-b")
  grDevices::pdf(NULL)
  expect_identical(plot(res), res)
  expect_identical(plot(res, n_x = 5, n_y = 3, clust = 'none', stars = NULL, key = FALSE), res)
  expect_identical(plot(rp, main = 'Mine', xlab = 'Variables'), rp)
  expect_identical(plot(res[res$y == 'ph', ]), res[res$y == 'ph', ])
  grDevices::dev.off()

#})
