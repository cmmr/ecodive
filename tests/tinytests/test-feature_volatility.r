# feature_volatility(), against q2-longitudinal 2026.10.0's and
# q2-sample-classifier 2026.10.0's own functions (fixtures from
# data-raw/expected/feature_volatility.py, with scikit-learn 1.9.1) and
# direct R implementations.

E    <- readRDS(file.path('expected', 'feature_volatility.rds'))
meta <- E$meta
cnt  <- E$counts
fv_details <- ecodive:::feature_volatility_details

max_cpus <- if (at_home()) 4L else 2L

rel_diff <- function (a, b) max(abs(a - b) / pmax(abs(b), 1e-300))



#--------------------------------------------------
# numpy's RandomState(rs).permutation(n), which train_test_split() draws.
#--------------------------------------------------
for (p in E$perm)
  expect_identical(
    as.integer(.Call(ecodive:::C_rf_rng_test, 5L, as.double(p$rs), as.double(p$n), p$n)),
    p$perm, info = paste(p$rs, p$n) )



#--------------------------------------------------
# q2's pipeline (without RFE), given its random_state: the same split,
# forest importances and held-out predictions, exactly; the same top
# features and statistics.
#--------------------------------------------------
for (cs in E$cases) {

  info <- paste('random_state', cs$random_state)
  fv <- fv_details(cnt, 'subject', 'state', 'group', data = meta, n_features = Inf,
                   norm = 'none', numpy_random_state = cs$random_state, cpus = 1L)
  smp <- rownames(cnt)

  expect_identical(smp[fv$split$train], cs$train, info = info)
  expect_identical(smp[fv$split$test], cs$test, info = info)

  # Importances, in q2's order.
  ft <- fv$features
  expect_identical(ft$feature, names(cs$importance), info = info)
  expect_equal(ft$importance, unname(cs$importance), tolerance = 1e-12, info = info)

  # The held-out samples' predictions and accuracy.
  pr <- fv$predictions
  expect_identical(sort(pr$sample), sort(cs$test), info = info)
  expect_identical(pr$sample, smp[sort(match(cs$test, smp))], info = info)
  expect_equal(pr$predicted[match(cs$test, pr$sample)], cs$predicted, tolerance = 1e-12, info = info)
  expect_equal(pr$observed, meta[pr$sample, 'state'], info = info)
  pf <- fv$performance
  expect_equal(pf[['mse']], cs$accuracy[['mse']], tolerance = 1e-12, info = info)
  expect_equal(pf[['r']], cs$accuracy[['r']], tolerance = 1e-12, info = info)
  expect_equal(pf[['r']]^2, cs$accuracy[['r2']], tolerance = 1e-12, info = info)
  expect_equal(pf[['slope']], cs$accuracy[['slope']], tolerance = 1e-12, info = info)
  expect_equal(pf[['intercept']], cs$accuracy[['intercept']], tolerance = 1e-12, info = info)

  # _importance_filtering() and _summarize_feature_stats(): all features, and the top 10.
  for (k in c('all', '10')) {
    ref <- cs[[paste0('stats_', k)]]
    res <- if (k == 'all') ft else
      fv_details(cnt, 'subject', 'state', 'group', data = meta, n_features = 10,
                 norm = 'none', numpy_random_state = cs$random_state, cpus = 1L)$features
    expect_identical(names(res), names(ref), info = info)
    expect_identical(res$feature, ref$feature, info = paste(info, k))
    for (col in names(ref)[-1])
      expect_true(rel_diff(res[[col]], ref[[col]]) < 1e-12, info = paste(info, k, col))
  }
}

# q2's actual pipeline also runs RFECV: the features it keeps are among the
# 10 most important of the same forest on every feature.
fv3 <- fv_details(cnt, 'subject', 'state', data = meta, n_features = 10, norm = 'none',
                  numpy_random_state = 3L, cpus = 1L)
expect_true(all(names(E$rfe) %in% fv3$features$feature))



#--------------------------------------------------
# The default (PCG) result.
#--------------------------------------------------
fv <- feature_volatility(cnt, 'subject', 'state', 'group', data = meta, n_features = 8, cpus = 1L)
expect_inherits(fv, 'ecodive_feature_volatility')
expect_identical(names(fv), c('features', 'means', 'abundances', 'samples', 'predictions', 'performance', 'forest'))
expect_identical(names(fv$features), c(
  'feature', 'importance', 'cumulative_avg_decrease', 'cumulative_avg_increase', 'net_avg_change',
  'global_variance', 'global_mean', 'global_median', 'global_sd', 'global_cv' ))
expect_identical(nrow(fv$features), 8L)
expect_identical(attr(fv, 'method'), 'feature_volatility')
expect_identical(attr(fv, 'conf_level'), 0.95)
expect_identical(attr(fv, 'n_resamples'), 1000L)
expect_identical(attr(fv, 'test_size'), 0.2)
expect_inherits(fv$forest, 'ecodive_random_forest')
expect_identical(attr(fv$forest, 'n_trees'), 100L)
expect_identical(attr(fv$forest, 'n_folds'), 0L)
expect_null(fv$forest$holdout)

# The forest's importances, ranked; the top 8 kept.
imp <- fv$forest$importance
expect_identical(fv$features$feature, imp$feature[1:8])
expect_identical(fv$features$importance, imp$impurity[1:8])
expect_false(is.unsorted(-imp$impurity))
expect_identical(nrow(imp), 25L)

# The split: ceiling(0.2 * 78) = 16 held out.
s <- fv$samples
expect_identical(names(s), c('sample', 'subject', 'group', 'state', 'set'))
expect_identical(s$sample, rownames(cnt))
expect_identical(s$subject, meta$subject)
expect_identical(s$state, meta$state)
expect_identical(sum(s$set == 'test'), 16L)
expect_identical(sort(fv$forest$predictions$sample), sort(s$sample[s$set == 'train']))
expect_identical(fv$predictions$sample, s$sample[s$set == 'test'])

# Held-out predictions: the kept forest's.
te <- s$set == 'test'
expect_identical(fv$predictions$predicted, unname(predict(fv$forest, cnt[te, ], cpus = 1L)))
obs <- fv$predictions$observed; pred <- fv$predictions$predicted
expect_equal(fv$performance[['mse']], mean((obs - pred)^2), tolerance = 1e-14)
expect_equal(fv$performance[['r_squared']], 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2), tolerance = 1e-14)
expect_equal(fv$performance[['r']], cor(obs, pred), tolerance = 1e-12)
lf <- coef(lm(pred ~ obs))
expect_equal(fv$performance[['slope']], unname(lf[[2]]), tolerance = 1e-12)
expect_equal(fv$performance[['intercept']], unname(lf[[1]]), tolerance = 1e-12)
expect_identical(fv$performance[['oob_r_squared']], fv$forest$performance[['oob_r_squared']])

# Relative abundances of the kept features.
rel <- cnt / rowSums(cnt)
expect_equal(fv$abundances, rel[, fv$features$feature], tolerance = 1e-14)

# The statistics, as R computes them.
for (f in fv$features$feature) {
  v  <- rel[, f]
  sm <- tapply(v, meta$state, mean)
  d  <- diff(sm[order(as.numeric(names(sm)))])
  r  <- fv$features[fv$features$feature == f, ]
  expect_equal(r$cumulative_avg_decrease, sum(d[d < 0]), tolerance = 1e-12, info = f)
  expect_equal(r$cumulative_avg_increase, sum(d[d > 0]), tolerance = 1e-12, info = f)
  expect_equal(r$net_avg_change, unname(sm[['28']] - sm[['0']]), tolerance = 1e-12, info = f)
  expect_equal(r$global_variance, var(v), tolerance = 1e-12, info = f)
  expect_equal(r$global_mean, mean(v), tolerance = 1e-14, info = f)
  expect_identical(r$global_median, median(v), info = f)
  expect_equal(r$global_sd, sd(v), tolerance = 1e-12, info = f)
  expect_equal(r$global_cv, sd(v) / mean(v), tolerance = 1e-12, info = f)
}

# The chart's means: volatility()'s numbers per feature (its own intervals
# use other streams).
m <- fv$means
expect_identical(names(m), c('feature', 'group', 'state', 'mean', 'lower', 'upper', 'n'))
expect_identical(nrow(m), 8L * 10L)
for (f in fv$features$feature[1:3]) {
  vo <- volatility(unname(rel[, f]), meta$subject, meta$state, meta$group)
  mf <- m[m$feature == f, -1L]
  rownames(mf) <- NULL
  expect_identical(mf[, c('group', 'state', 'mean', 'n')], vo$means[, c('group', 'state', 'mean', 'n')], info = f)
}
expect_true(all(m$lower <= m$mean & m$mean <= m$upper))

# Vega's bootstrapCI() in R, on the same uniforms: cell k of feature j
# (both 1-based) draws from stream 2^52 + (j - 1) * cells + k - 1.
boot_ci <- function (v, u, nres, alpha) {
  n  <- length(v)
  mu <- colSums(matrix(v[floor(u * n) + 1], n, nres)) / n
  stats::quantile(sort(mu), c(alpha / 2, 1 - alpha / 2), type = 7, names = FALSE)
}
for (row in c(1L, 7L, 34L, 80L)) {
  j <- match(m$feature[[row]], fv$features$feature)
  k <- row - (j - 1L) * 10L
  i <- which(meta$group == m$group[[row]] & meta$state == m$state[[row]])
  u <- .Call(ecodive:::C_rng_dist, as.integer(length(i) * 1000L), 1L, 0, 0L, 2^52 + (j - 1) * 10 + k - 1)
  expect_equal(c(m$lower[[row]], m$upper[[row]]), boot_ci(rel[i, m$feature[[row]]], u, 1000L, 0.05),
               tolerance = 1e-14, info = row)
}

# Display.
expect_identical(attr(fv, 'display')$title, 'Feature Volatility')
expect_identical(attr(fv, 'display')$xlab, 'state')
expect_identical(attr(fv, 'display')$legend, 'group')



#--------------------------------------------------
# Determinism: any cpus, input types and seeds.
#--------------------------------------------------
for (cpus in 2:max_cpus)
  expect_identical(feature_volatility(cnt, 'subject', 'state', 'group', data = meta, n_features = 8, cpus = cpus),
                   fv, info = cpus)
fp <- feature_volatility(cnt, 'subject', 'state', data = meta, n_features = 5,
                         importance = 'permutation', n_trees = 30, cpus = 1L)
for (cpus in 2:max_cpus)
  expect_identical(feature_volatility(cnt, 'subject', 'state', data = meta, n_features = 5,
                                      importance = 'permutation', n_trees = 30, cpus = cpus),
                   fp, info = cpus)
expect_identical(fp$features$importance, fp$forest$importance$permutation[1:5])
expect_false(is.unsorted(-fp$forest$importance$permutation))

# The same samples and abundances from a transposed table; dense tables
# give identical results.
ft <- feature_volatility(t(cnt), 'subject', 'state', 'group', data = meta, margin = 2L,
                         n_features = 8, cpus = 1L)
expect_identical(ft, fv)
st <- structure(
  list(i = row(cnt)[cnt > 0], j = col(cnt)[cnt > 0], v = as.double(cnt[cnt > 0]),
       nrow = nrow(cnt), ncol = ncol(cnt), dimnames = dimnames(cnt)),
  class = 'simple_triplet_matrix' )
fs <- feature_volatility(st, 'subject', 'state', 'group', data = meta, n_features = 25, cpus = 2L)
fd <- feature_volatility(cnt, 'subject', 'state', 'group', data = meta, n_features = 25, cpus = 2L)
expect_identical(fs$samples, fd$samples)
o <- match(fd$features$feature, fs$features$feature)
expect_equal(fs$abundances[, o], fd$abundances, tolerance = 1e-15)
expect_equal(fs$features[o, -2], fd$features[, -2], tolerance = 1e-15, check.attributes = FALSE)
expect_true(cor(fs$features$importance[o], fd$features$importance) > 0.99)

# Another seed: another split.
f7 <- feature_volatility(cnt, 'subject', 'state', 'group', data = meta, n_features = 8, seed = 7L, cpus = 1L)
expect_false(identical(f7$samples$set, fv$samples$set))
expect_identical(sum(f7$samples$set == 'test'), 16L)



#--------------------------------------------------
# Options.
#--------------------------------------------------
# No hold-out: every sample trains; accuracy only out of bag.
f0 <- feature_volatility(cnt, 'subject', 'state', data = meta, test_size = 0, n_features = Inf,
                         n_trees = 20, cpus = 1L)
expect_null(f0$predictions)
expect_true(all(f0$samples$set == 'train'))
expect_true(all(is.na(f0$performance[1:5])) && is.finite(f0$performance[['oob_r_squared']]))
expect_identical(nrow(f0$features), 25L)
expect_identical(f0$forest$predictions$sample, rownames(cnt))
expect_identical(names(f0$samples), c('sample', 'subject', 'group', 'state', 'set'))
expect_identical(levels(f0$samples$group), 'all')
expect_identical(levels(f0$means$group), 'all')
expect_null(attr(f0, 'display')$legend)

# n_features larger than the table; no subjects; states and groups as vectors.
f1 <- feature_volatility(cnt, states = meta$state, groups = meta$group, n_features = 1000,
                         n_trees = 20, cpus = 1L)
expect_identical(nrow(f1$features), 25L)
expect_identical(names(f1$samples), c('sample', 'group', 'state', 'set'))
expect_null(attr(f1, 'display')$xlab)

# Unnamed samples and features.
cu <- unname(cnt)
f2 <- feature_volatility(cu, states = meta$state, n_features = 3, n_trees = 20, cpus = 1L)
expect_identical(f2$samples$sample, as.character(seq_len(nrow(cnt))))
expect_true(all(grepl('^F[0-9]+$', f2$features$feature)))
expect_equal(unname(f2$abundances), unname(rel[, as.integer(sub('F', '', f2$features$feature))]),
             tolerance = 1e-14)

# Forest settings through `...`; norm.
f4 <- feature_volatility(cnt, 'subject', 'state', data = meta, n_features = 5, n_trees = 20,
                         max_features = 'sqrt', min_samples_leaf = 2L, n_folds = 3L,
                         keep_forest = FALSE, norm = 'clr', pseudocount = 0.5, cpus = 1L)
expect_identical(attr(f4$forest, 'max_features'), 5L)
expect_identical(attr(f4$forest, 'min_samples_leaf'), 2L)
expect_identical(attr(f4$forest, 'n_folds'), 3L)
expect_null(f4$forest$forest)
expect_true(is.finite(f4$forest$performance[['r_squared']]))
expect_equal(f4$abundances, rel[, f4$features$feature], tolerance = 1e-14)

# Interval settings.
f5 <- feature_volatility(cnt, 'subject', 'state', data = meta, n_features = 2, n_trees = 20,
                         conf_level = 0.8, n_resamples = 50L, cpus = 1L)
expect_identical(attr(f5, 'conf_level'), 0.8)
expect_identical(attr(f5, 'n_resamples'), 50L)



#--------------------------------------------------
# Errors.
#--------------------------------------------------
fvq <- function (...) feature_volatility(cnt, 'subject', 'state', data = meta, n_trees = 5, cpus = 1L, ...)
expect_error(fvq(n_features = 0), '`n_features` must be a whole number')
expect_error(fvq(n_features = 2.5), '`n_features` must be a whole number')
expect_error(fvq(n_features = NA), '`n_features` must be a whole number')
expect_error(fvq(test_size = 1), '`test_size` must be a number')
expect_error(fvq(test_size = -0.1), '`test_size` must be a number')
expect_error(fvq(test_size = 0.99), 'leaves fewer than two')
expect_error(fvq(importance = 'both'), "`importance` must be 'impurity' or 'permutation'")
expect_error(fvq(conf_level = 1), '`conf_level`')
expect_error(fvq(n_resamples = 0), '`n_resamples`')
expect_error(fvq(class_weight = 'balanced'), "Unknown argument 'class_weight'")
expect_error(fvq(folds = 1), 'Unknown argument')
expect_error(fvq(max_features = 0), '`max_features`')
expect_error(feature_volatility(cnt, 'subject', data = meta), '`states` is required')
m2 <- meta; m2$state[[3]] <- NA
expect_error(feature_volatility(cnt, 'subject', 'state', data = m2), 'Missing values')
m2 <- meta; m2$state <- as.character(m2$state)
expect_error(feature_volatility(cnt, 'subject', 'state', data = m2), 'must be numeric')
expect_error(feature_volatility(cnt, states = rep(1, nrow(cnt))), 'at least two distinct')
expect_error(feature_volatility(cnt, 'state', 'state', data = meta), 'different columns')
expect_error(feature_volatility(cnt[1, , drop = FALSE], states = 1), 'at least two samples')
expect_error(feature_volatility(cnt, states = meta$state, seed = 'a'))
expect_error(feature_volatility(cnt, states = meta$state[-1]))



#--------------------------------------------------
# print() and plot().
#--------------------------------------------------
out <- capture.output(print(fv))
expect_identical(out[[1]], 'Feature Volatility')
expect_true(grepl('trained on 62 of 78 samples; on the 16 held out', out[[2]]))
expect_true(any(grepl('The 8 most important of 25 features', out)))
expect_true(any(grepl('F01', out)))
expect_true(grepl('trained on 78 of 78 samples; OOB', capture.output(print(f0))[[2]]))

had_seed <- exists('.Random.seed', envir = globalenv())
seed0    <- if (had_seed) get('.Random.seed', envir = globalenv())
grDevices::pdf(NULL)
par0 <- graphics::par(no.readonly = TRUE)
drawn <- function (x, ...) {
  v <- withVisible(plot(x, ...))
  identical(v$value, x) && !v$visible
}
expect_true(drawn(fv))
expect_true(drawn(fv, feature = 'F03', n_bars = 3, spaghetti = FALSE, limits = FALSE,
                  legend = NULL, col = c('red', 'blue'), main = 'x', ylim = c(0, 1)))
expect_true(drawn(fv, feature = 8))
expect_true(drawn(f0))
expect_true(drawn(f1))
par1 <- graphics::par(no.readonly = TRUE)
expect_identical(par1[setdiff(names(par1), c('usr', 'xaxp', 'yaxp'))], par0[setdiff(names(par0), c('usr', 'xaxp', 'yaxp'))])
expect_error(plot(fv, feature = 'nope'), '`feature` must be the name')
expect_error(plot(fv, feature = 9), '`feature` must be the name')
expect_error(plot(fv, n_bars = 0), '`n_bars`')
expect_error(plot(fv, legend = 'middle'), '`legend`')
grDevices::dev.off()
expect_identical(exists('.Random.seed', envir = globalenv()), had_seed)
if (had_seed) expect_identical(get('.Random.seed', envir = globalenv()), seed0)
