# neutral_model(): Sloan et al.'s (2006) neutral community model.
#
# References (data-raw/expected/neutral_model.py and .r, folded into
# expected/neutral_model.rds):
#   * neufit (neufit.py, BSD-2), run as is on five tables: its m (lmfit, so
#     only to its stopping tolerance), R^2, standard error, predictions and
#     Wilson bands;
#   * a 40-digit mpmath reference for the least-squares m and its profile
#     interval (unbounded m);
#   * Burns et al.'s code, through iCAMP::snm() and tyRa's fit_sncm() (run,
#     never read): m, its confint() interval, R^2, RMSE, AIC and BIC of the
#     neutral, binomial and Poisson models, predictions and bands;
#   * statsmodels' and Hmisc's Wilson intervals at the edges.

E <- readRDS(file.path('expected', 'neutral_model.rds'))

rel <- function (a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
q   <- function (expr) suppressWarnings(expr)
A   <- function (x, nm) attr(x, nm, exact = TRUE)

if (requireNamespace('Matrix', quietly = TRUE)) {
  if (utils::packageVersion("Matrix") >= "1.4.2") {
    COMPRESSED <- "CsparseMatrix"
  } else {
    COMPRESSED <- "dgCMatrix"
  }
}


#--------------------------------------------------
# The 40-digit reference: m and its profile interval.
#--------------------------------------------------
for (i in seq_len(nrow(E$hp))) {
  h  <- E$hp[i, ]
  x  <- E$tables[[if (h$table == 'dm_pool') 'dm' else h$table]]
  r  <- q(if (h$table == 'dm_pool') neutral_model(x, pool = E$tables$pool) else neutral_model(x))
  expect_identical(nrow(r), h$K, info = h$table)
  expect_true(rel(A(r, 'm'), h$m) < 1e-10, info = h$table)
  expect_true(rel(A(r, 'm_lower'), h$m_lower) < 1e-10, info = h$table)
  expect_true(rel(A(r, 'm_upper'), h$m_upper) < 1e-10, info = h$table)
  expect_true(rel(stats::qf(0.95, 1, h$K - 1), h$cutoff) < 1e-12, info = h$table)
}


#--------------------------------------------------
# neufit.
#--------------------------------------------------
for (nm in names(E$neufit)) {

  nf  <- E$neufit[[nm]]
  tab <- if (nm == 'uneven') 'uneven_rarefied' else nm
  r   <- q(neutral_model(E$tables[[tab]]))
  i   <- match(nf$table$feature, r$feature)
  expect_false(anyNA(i), info = nm)

  expect_equal(A(r, 'n_individuals'), nf$n_reads, info = nm)
  expect_identical(A(r, 'n_samples'), nf$n_samples, info = nm)
  expect_true(rel(r$mean_abundance[i], nf$table$mean_abundance) < 1e-14, info = nm)
  expect_identical(r$occurrence[i], nf$table$occurrence, info = nm)

  # Given neufit's m: its predictions and Wilson bands (statsmodels).
  N <- nf$n_reads; n <- nf$n_samples
  pr <- ecodive:::neutral_curve(nf$table$mean_abundance, N, nf$m)
  expect_true(rel(pr, nf$table$predicted) < 1e-12, info = nm)
  w <- ecodive:::neutral_wilson(nf$table$predicted, n, 0.95)
  expect_true(max(abs(w$lower - nf$table$lower)) < 1e-14, info = nm)
  expect_true(max(abs(w$upper - nf$table$upper)) < 1e-14, info = nm)

  if (nm == 'sat') {
    # neufit stops at its bound, m = 1; the least-squares m is above 1.
    expect_true(abs(nf$m - 1) < 1e-6)
    expect_true(A(r, 'm') > 1.2)
    next
  }

  # Its own fit: lmfit stops at a relative tolerance of about 1e-6 in m.
  expect_true(rel(A(r, 'm'), nf$m) < 1e-5, info = nm)
  expect_true(rel(A(r, 'r_squared'), nf$r_squared) < 1e-8, info = nm)
  expect_true(rel(A(r, 'm_se'), nf$m_stderr) < 1e-3, info = nm)
  expect_true(max(abs(r$predicted[i] - nf$table$predicted)) < 1e-5, info = nm)
  expect_true(max(abs(r$lower[i] - nf$table$lower)) < 1e-5, info = nm)
  expect_true(max(abs(r$upper[i] - nf$table$upper)) < 1e-5, info = nm)
}


#--------------------------------------------------
# Burns et al.'s code: iCAMP::snm() and tyRa's fit_sncm().
#--------------------------------------------------
for (cs in names(E$burns)) {

  b    <- E$burns[[cs]]
  x    <- E$tables[[if (cs == 'dm_pool') 'dm' else cs]]
  pool <- if (cs == 'dm_pool') E$tables$pool
  r    <- q(neutral_model(x, pool = pool))

  ty <- b$tyra$stats
  td <- b$tyra$detail
  i  <- match(td$feature, r$feature)
  expect_false(anyNA(i), info = cs)
  expect_identical(nrow(td), nrow(r), info = cs)

  # tyRa: N is the mean depth; p and occurrence; m (nlsLM's tolerance);
  # m.ci, the lower half-width of confint()'s interpolated profile.
  expect_equal(A(r, 'n_individuals'), ty$N, info = cs)
  expect_identical(A(r, 'n_samples'), as.integer(ty$Samples), info = cs)
  expect_true(rel(r$mean_abundance[i], td$p) < 1e-14, info = cs)
  expect_true(rel(r$occurrence[i], td$freq) < 1e-14, info = cs)
  expect_true(rel(A(r, 'm'), ty$m) < 1e-5, info = cs)
  expect_true(rel(A(r, 'm') - A(r, 'm_lower'), ty$m.ci) < 1e-3, info = cs)
  expect_true(rel(A(r, 'r_squared'), ty$Rsqr) < 1e-6, info = cs)
  expect_true(rel(A(r, 'rmse'), ty$RMSE) < 1e-5, info = cs)
  expect_true(rel(A(r, 'r_squared_poisson'), ty$Rsqr.pois) < 1e-12, info = cs)
  expect_true(rel(A(r, 'rmse_poisson'), ty$RMSE.pois) < 1e-12, info = cs)
  expect_true(rel(A(r, 'aic_poisson'), ty$AIC.pois) < 1e-6, info = cs)
  expect_true(rel(A(r, 'bic_poisson'), ty$BIC.pois) < 1e-6, info = cs)

  # Given tyRa's m: predictions, bands (Hmisc's Wilson) and partitions, away
  # from the band's edges (where the last bit decides).
  N  <- ty$N; n <- as.integer(ty$Samples)
  pr <- ecodive:::neutral_curve(td$p, N, ty$m)
  expect_true(rel(pr, td$freq.pred) < 1e-12, info = cs)
  w  <- ecodive:::neutral_wilson(td$freq.pred, n, 0.95)
  expect_true(max(abs(w$lower - td$pred.lwr)) < 1e-14, info = cs)
  expect_true(max(abs(w$upper - td$pred.upr)) < 1e-14, info = cs)
  part <- ecodive:::neutral_partition(td$freq, w)
  tp   <- c('As predicted' = 'neutral', 'Above prediction' = 'above',
            'Below prediction' = 'below')[td$fit_class]
  away <- pmin(abs(td$freq - w$lower), abs(td$freq - w$upper)) > 1e-9
  expect_identical(unname(part[away]), unname(tp[away]), info = cs)

  # iCAMP keeps features whose p in the pool is 0, and rounds N: compare the
  # rest.
  if (cs %in% c('dm_pool', 'uneven')) next

  ic <- b$icamp$stats
  id <- b$icamp$detail
  i  <- match(id$feature, r$feature)
  expect_false(anyNA(i), info = cs)
  expect_true(rel(A(r, 'm'), ic$m) < 2e-5, info = cs)
  expect_true(rel(A(r, 'm_lower'), ic$m.ci.2.5) < 1e-3, info = cs)
  expect_true(rel(A(r, 'm_upper'), ic$m.ci.97.5) < 1e-3, info = cs)
  expect_true(rel(A(r, 'r_squared'), ic$Rsqr) < 1e-6, info = cs)
  expect_true(rel(A(r, 'rmse'), ic$RMSE) < 1e-5, info = cs)
  expect_true(rel(A(r, 'aic'), ic$AIC) < 1e-5, info = cs)
  expect_true(rel(A(r, 'bic'), ic$BIC) < 1e-5, info = cs)
  for (k in c('r_squared', 'rmse'))
    for (alt in c('binomial', 'poisson')) {
      key <- paste0(c(r_squared = 'Rsqr', rmse = 'RMSE')[[k]], '.', substr(alt, 1, 4))
      expect_true(rel(A(r, paste0(k, '_', alt)), ic[[key]]) < 1e-12, info = paste(cs, key))
    }
  for (k in c('aic', 'bic'))
    for (alt in c('binomial', 'poisson')) {
      key <- paste0(toupper(k), '.', substr(alt, 1, 4))
      expect_true(rel(A(r, paste0(k, '_', alt)), ic[[key]]) < 1e-5, info = paste(cs, key))
    }
  expect_true(rel(r$predicted_binomial[i], id$bino.pred) < 1e-12, info = cs)
  expect_true(rel(r$predicted_poisson[i], id$pois.pred) < 1e-12, info = cs)
  expect_true(rel(ecodive:::neutral_curve(id$p, ic$N, ic$m), id$freq.pred) < 1e-12, info = cs)
}

# iCAMP with unequal depths rounds N (to 956 here); given that N, its p
# differs only by the rounding, and its m is the least-squares m.
r  <- q(neutral_model(E$tables$uneven))
ic <- E$burns$uneven$icamp$stats
expect_equal(ic$N, round(A(r, 'n_individuals')))


#--------------------------------------------------
# The closed-form Gaussian log-likelihood equals a direct maximization.
#--------------------------------------------------
r  <- neutral_model(E$tables$sim)
K  <- nrow(r)
ll <- function (res, mu) {
  nll <- function (th) -sum(stats::dnorm(res, if (mu) th[[2]] else 0, exp(th[[1]]), log = TRUE))
  o   <- stats::optim(c(log(stats::sd(res)), mean(res))[seq_len(1 + mu)], nll,
                      method = 'BFGS', control = list(reltol = 1e-15, maxit = 1000))
  -o$value
}
aic <- function (ll) -2 * ll + 4
expect_true(rel(A(r, 'aic'), aic(ll(r$occurrence - r$predicted, FALSE))) < 1e-8)
expect_true(rel(A(r, 'aic_binomial'), aic(ll(r$occurrence - r$predicted_binomial, TRUE))) < 1e-8)
expect_true(rel(A(r, 'aic_poisson'), aic(ll(r$occurrence - r$predicted_poisson, TRUE))) < 1e-8)
expect_equal(A(r, 'bic') - A(r, 'aic'), 2 * log(K) - 4)
expect_equal(A(r, 'nm'), A(r, 'm') * 1000)
expect_equal(A(r, 'detection'), 1e-3)

# The binomial model is pbinom(1/N, N, p, lower.tail = FALSE).
expect_true(rel(r$predicted_binomial, stats::pbinom(1e-3, 1000, r$mean_abundance, lower.tail = FALSE)) < 1e-12)

# The m minimizes the residual sum of squares.
ss <- function (m) sum((r$occurrence - ecodive:::neutral_curve(r$mean_abundance, 1000, m)) ^ 2)
expect_true(ss(A(r, 'm')) <= ss(A(r, 'm') * (1 + 1e-6)))
expect_true(ss(A(r, 'm')) <= ss(A(r, 'm') * (1 - 1e-6)))

# The delta-method standard error.
m <- A(r, 'm'); h <- 1e-5 * m
J <- (ecodive:::neutral_curve(r$mean_abundance, 1000, m + h) -
      ecodive:::neutral_curve(r$mean_abundance, 1000, m - h)) / (2 * h)
expect_true(rel(A(r, 'm_se'), sqrt(ss(m) / (K - 1) / sum(J ^ 2))) < 1e-6)


#--------------------------------------------------
# The Wilson interval at the edges: statsmodels and Hmisc.
#--------------------------------------------------
for (ref in list(E$wilson_sm, E$wilson_hmisc)) {
  for (k in seq_len(nrow(ref))) {
    w <- ecodive:::neutral_wilson(ref$count[k] / ref$n[k], ref$n[k], ref$conf_level[k])
    expect_true(abs(w$lower - ref$lower[k]) < 1e-14, info = k)
    expect_true(abs(w$upper - ref$upper[k]) < 1e-14, info = k)
  }
}
w <- ecodive:::neutral_wilson(c(0, 1), 10, 0.95)
expect_identical(w$lower[[1]], 0)
expect_identical(w$upper[[2]], 1)


#--------------------------------------------------
# Draws from the model: m is recovered, and the interval covers it in most
# replicates (it assumes equal variances, so less often than 95%).
#--------------------------------------------------
set.seed(20261124)
hits <- 0L; ms <- numeric(0)
for (rep in 1:20) {
  p   <- stats::rgamma(60, 0.3); p <- p / sum(p)
  x   <- matrix(stats::rbeta(60 * 30, 1000 * 0.05 * p, 1000 * 0.05 * (1 - p)), 60)
  occ <- rowMeans(x > 1e-3)
  k   <- occ > 0
  o   <- order(p[k])
  f   <- q(ecodive:::neutral_fit(p[k][o], occ[k][o], 1000, 30L, 0.95))
  ms  <- c(ms, f$m)
  hits <- hits + (f$m_lower <= 0.05 && f$m_upper >= 0.05)
}
expect_true(hits >= 14L)
expect_true(abs(mean(ms) / 0.05 - 1) < 0.05)


#--------------------------------------------------
# Input types, margin, threads.
#--------------------------------------------------
x    <- E$tables$dm
base <- neutral_model(x, cpus = 1L)
expect_identical(neutral_model(t(x), margin = 2L, cpus = 1L), base)
expect_identical(neutral_model(x * 1.0, cpus = 1L), base)
if (requireNamespace('Matrix', quietly = TRUE)) {
  expect_identical(neutral_model(as(as(x * 1.0, 'CsparseMatrix'), COMPRESSED), cpus = 2L), base)
  expect_identical(neutral_model(Matrix::t(as(as(x * 1.0, 'CsparseMatrix'), COMPRESSED)), margin = 2L, cpus = 3L), base)
}
if (requireNamespace('slam', quietly = TRUE))
  expect_identical(neutral_model(slam::as.simple_triplet_matrix(x), cpus = 2L), base)

# Threads: a table large enough to start them, by group, with a pool.
set.seed(1)
big <- matrix(stats::rpois(300 * 1000, 2), 300, 1000)
big[, 1] <- big[, 1] + 1L
big[, 1] <- big[, 1] + (5000L - rowSums(big))   # equal depths
grp  <- rep(c('a', 'b', 'c'), 100)
pl   <- big[1:20, ] + 1L
b1   <- q(neutral_model(big, pool = pl, groups = grp, cpus = 1L))   # m > 1: Poisson counts
for (cpus in c(2L, 4L, 7L))
  expect_identical(q(neutral_model(big, pool = pl, groups = grp, cpus = cpus)), b1, info = cpus)
st1 <- .Call(ecodive:::C_neutral_stats, big, 1L, as.integer(factor(grp)) - 1L, 3L, 1L)
for (cpus in c(2L, 4L, 7L))
  expect_identical(.Call(ecodive:::C_neutral_stats, big, 1L, as.integer(factor(grp)) - 1L, 3L, cpus), st1)
expect_identical(unname(st1$total), unname(rowsum(big * 1.0, grp)))
expect_identical(unname(st1$occ), unname(rowsum((big > 0) * 1L, grp)))


#--------------------------------------------------
# The result.
#--------------------------------------------------
r <- neutral_model(E$tables$sim)
expect_true(inherits(r, c('ecodive_neutral_model', 'data.frame')))
expect_identical(
  names(r),
  c('feature', 'mean_abundance', 'occurrence', 'predicted', 'lower', 'upper',
    'partition', 'predicted_binomial', 'predicted_poisson') )
expect_identical(levels(r$partition), c('above', 'neutral', 'below'))
expect_false(is.unsorted(r$mean_abundance))
expect_identical(A(r, 'method'), 'neutral_model')
expect_identical(A(r, 'conf_level'), 0.95)
expect_identical(A(r, 'dropped'), character(0))
expect_true(all(r$lower <= r$predicted & r$predicted <= r$upper))
expect_identical(
  as.character(r$partition),
  ifelse(r$occurrence > r$upper, 'above', ifelse(r$occurrence < r$lower, 'below', 'neutral')) )

# Features without counts are dropped; unnamed features are numbered.
x <- E$tables$sim
x[, 1] <- x[, 1] + x[, 3] + x[, 7]
x[, c(3, 7)] <- 0L
r0 <- neutral_model(x)
expect_identical(A(r0, 'dropped'), colnames(x)[c(3, 7)])
expect_identical(nrow(r0), 98L)
r0 <- neutral_model(unname(x))
expect_identical(A(r0, 'dropped'), c('3', '7'))
expect_true(all(r0$feature %in% as.character(seq_len(100))))

# conf_level widens the interval and the band.
r9 <- neutral_model(E$tables$sim, conf_level = 0.99)
expect_true(A(r9, 'm_lower') < A(r, 'm_lower') && A(r9, 'm_upper') > A(r, 'm_upper'))
expect_true(all(r9$upper >= r$upper & r9$lower <= r$lower))
expect_identical(A(r9, 'm'), A(r, 'm'))


#--------------------------------------------------
# pool.
#--------------------------------------------------
x  <- E$tables$dm
pl <- E$tables$pool
rp <- neutral_model(x, pool = pl)
expect_identical(rp$mean_abundance, unname((colSums(pl) / sum(pl))[rp$feature]))
# Matched by name, in any order; the pool's other features count in its total.
expect_identical(neutral_model(x, pool = pl[, rev(seq_len(ncol(pl)))]), rp)
expect_identical(neutral_model(t(x), pool = t(pl), margin = 2L), rp)
expect_error(neutral_model(x, pool = pl[, !colnames(pl) %in% c('F1', 'F2')]), "2 features in `counts` are not in `pool`: 'F1', 'F2'")
# A feature absent from counts may be missing from the pool.
x2 <- x; x2[, 'F2'] <- x2[, 'F2'] + x2[, 'F1']; x2[, 'F1'] <- 0L
expect_silent(r2 <- neutral_model(x2, pool = pl[, colnames(pl) != 'F1']))
expect_true('F1' %in% A(r2, 'dropped'))
# By position without names.
pn <- pl[, colnames(x)]
ru <- neutral_model(unname(x), pool = unname(pn))
expect_identical(ru$mean_abundance, unname(colSums(pn) / sum(pn))[as.integer(ru$feature)])
expect_error(neutral_model(unname(x), pool = unname(pl)), '`pool` has 150 features, but `counts` has 110')
expect_error(neutral_model(x, pool = pl / 2), '`pool` must be counts')
expect_error(neutral_model(x, pool = pl * 0L), '`pool` is empty')
expect_error(neutral_model(x, pool = 'a'), '`pool`')


#--------------------------------------------------
# groups.
#--------------------------------------------------
x   <- E$tables$sim
grp <- rep(c('b', 'a'), 25)
rg  <- neutral_model(x, groups = grp)
expect_identical(levels(rg$group), c('a', 'b'))
expect_identical(names(rg)[[1]], 'group')
for (g in c('a', 'b')) {
  one <- neutral_model(x[grp == g, ])
  sub <- rg[rg$group == g, ]
  for (col in names(one))
    expect_identical(sub[[col]], one[[col]], info = paste(g, col))
  for (nm in c('m', 'm_lower', 'm_upper', 'm_se', 'nm', 'n_individuals', 'detection', 'n_samples',
               'r_squared', 'rmse', 'aic', 'bic', 'aic_binomial', 'aic_poisson'))
    expect_identical(unname(A(rg, nm)[[g]]), A(one, nm), info = paste(g, nm))
  expect_identical(A(rg, 'dropped')[[g]], A(one, 'dropped'))
}
expect_identical(names(A(rg, 'm')), c('a', 'b'))

# Factor levels keep their order; from `data`; NA groups are left out.
meta <- data.frame(hab = factor(grp, levels = c('b', 'a')), row.names = rownames(x))
rf   <- neutral_model(x, 'hab', meta)
expect_identical(levels(rf$group), c('b', 'a'))
expect_identical(unname(A(rf, 'm')), unname(A(rg, 'm')[c('b', 'a')]))
expect_identical(neutral_model(x, meta$hab), rf)
expect_identical(neutral_model(x[50:1, ], 'hab', meta), rf)
gna <- grp; gna[1:2] <- NA
rn  <- neutral_model(x, groups = gna)
expect_identical(unname(A(rn, 'n_samples')), c(24L, 24L))
expect_identical(rn$feature[rn$group == 'b'], neutral_model(x[grp == 'b', ][-1, ])$feature)

expect_error(neutral_model(x, groups = c('a', rep('b', 49))), "at least 2 samples in group 'a', but there is 1")
expect_error(neutral_model(x, groups = rep(NA, 50)), 'no groups')
expect_error(neutral_model(x, groups = 'hab'), 'names a column, but `data` is not given')
expect_error(neutral_model(x, groups = 'nope', data = meta), 'not one')
expect_error(neutral_model(x, groups = 1:3), '`groups`')


#--------------------------------------------------
# Warnings and errors.
#--------------------------------------------------
x <- E$tables$sim

expect_error(neutral_model(x[1, , drop = FALSE]), 'at least 2 samples, but there is 1')
expect_error(neutral_model(x / 1000), 'must be counts \\(whole numbers\\), not proportions')
expect_error(neutral_model(-x), 'negative')
y <- x * 1.0; y[1, 1] <- NA
expect_error(neutral_model(y), 'NA')
expect_error(neutral_model(x * 0L), 'Every sample in `counts` is empty')
expect_error(neutral_model(x, conf_level = 1), '`conf_level`')
expect_error(neutral_model(x, conf_level = c(0.9, 0.95)), '`conf_level`')
expect_error(neutral_model(x, data = 1:3), '`data` must be a data.frame')
expect_error(neutral_model(x, cpus = 0), '`cpus`')
expect_error(neutral_model(x, margin = 3), '`margin`')

# One feature with counts: no fit.
y <- matrix(0L, 5, 4); y[, 1] <- 100L
expect_error(neutral_model(y), 'at least 2 features with counts, but there is 1')

# Unequal depths: N is the mean, with a warning.
expect_warning(ru <- neutral_model(E$tables$uneven), 'assumes equal depths')
expect_equal(A(ru, 'n_individuals'), mean(rowSums(E$tables$uneven)))
expect_warning(neutral_model(E$tables$uneven, groups = rep(1:2, 15)), "groups '1', '2' have different totals")

# Empty samples are left out.
y <- x; y[2, ] <- 0L
expect_warning(ry <- neutral_model(y), "Left out 1 sample with no counts: 'sample_2'")
expect_identical(A(ry, 'n_samples'), 49L)
expect_identical(ry, neutral_model(x[-2, ]))

# m above 1 (near-saturated data).
expect_warning(rs <- neutral_model(E$tables$sat), 'above 1')
expect_true(A(rs, 'm') > 1)

# Every feature in every sample: R^2 undefined, m not identifiable.
y <- matrix(c(500L, 300L, 200L), 10, 3, byrow = TRUE)
w <- character(0)
ra <- withCallingHandlers(neutral_model(y), warning = function (e) {
  w <<- c(w, conditionMessage(e)); invokeRestart('muffleWarning') })
expect_true(any(grepl('R-squared is undefined', w)))
expect_true(any(grepl('not identifiable', w)))
expect_true(is.na(A(ra, 'r_squared')) && is.na(A(ra, 'm_lower')) && is.na(A(ra, 'aic')))
expect_true(all(ra$partition == 'neutral'))

# Two features: an interval that may be unbounded.
y <- matrix(0L, 5, 3); y[, 1] <- c(100L, 100L, 99L, 99L, 100L); y[, 2] <- c(0L, 0L, 1L, 1L, 0L)
r2 <- neutral_model(y)
expect_identical(nrow(r2), 2L)
expect_true(A(r2, 'm_lower') <= A(r2, 'm') && A(r2, 'm') <= A(r2, 'm_upper'))


#--------------------------------------------------
# Printing and plots.
#--------------------------------------------------
r   <- neutral_model(E$tables$sim)
out <- utils::capture.output(v <- withVisible(print(r)))
expect_false(v$visible)
expect_true(any(grepl('Neutral Community Model', out, fixed = TRUE)))
expect_true(any(grepl('m = 0.06176 (95% CI 0.05749 to 0.06642), Nm = 61.76, R2 = 0.955; 10 above, 88 neutral, 2 below', out, fixed = TRUE)))
out <- utils::capture.output(print(rg))
expect_true(any(grepl('^# a: m = ', out)) && any(grepl('^# b: m = ', out)))
out <- utils::capture.output(print(r[1:3, ]))
expect_true(any(grepl('Neutral Community Model', out, fixed = TRUE)))

drawn <- function (x, ...) {
  v <- withVisible(plot(x, ...))
  identical(v$value, x) && !v$visible
}
grDevices::pdf(NULL)
expect_true(drawn(r))
expect_true(drawn(r, 'binomial', legend = 'topleft'))
expect_true(drawn(r, legend = NULL, main = 'x', pch = 1))
expect_true(drawn(rg))
expect_true(drawn(rg, group = 'b'))
expect_true(drawn(rg, group = 2:1, type = 'binomial'))
expect_error(plot(r, 'nope'), '`type`')
expect_error(plot(r, group = 'a'), 'only for a result made with `groups`')
expect_error(plot(rg, group = 'z'), "`group` must name groups of `x`: 'a', 'b'")
if (!is.null(attr(r[1:2, ], 'm'))) {   # where `[` keeps the attributes
  expect_true(drawn(r[r$partition == 'above', ]))
  expect_true(drawn(rg[rg$group == 'b', ]))
}
y <- r; attr(y, 'm') <- NULL
expect_error(plot(y), 'lost the attributes')
expect_error(plot(r[0, ]), 'no rows')
expect_error(plot(r, legend = 'middle'), '`legend`')
grDevices::dev.off()

# The title: "Neutral Community Model" (with groups, and the group), from
# the recorded display list.
titles <- function (x, ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control('enable')
  plot(x, ...)
  txt <- unlist(lapply(grDevices::recordPlot()[[1]], function (op) {
    a <- unlist(op[[2]], recursive = TRUE)
    a[vapply(a, is.character, logical(1))]
  }))
  gsub('\n', ' ', unique(unlist(txt)))
}
expect_true('Neutral Community Model' %in% titles(r))
expect_true(all(c('Neutral Community Model: a', 'Neutral Community Model: b') %in% titles(rg)))
expect_true('m = 0.0618, Nm = 61.8, R2 = 0.955' %in% titles(r))
expect_true('x' %in% titles(r, main = 'x'))
expect_identical(ecodive:::neutral_model_label(r)[[1]], 'Neutral Community Model')
