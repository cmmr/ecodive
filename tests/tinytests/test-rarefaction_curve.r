#test_that("rarefaction curves", {

  # rarefaction_curve(): the closed forms against vegan::rarefy(), iNEXT's
  # sample coverage, a 50-digit mpmath computation of the published formulas
  # and a direct lchoose() computation; the resampled curves rebuilt from
  # ecodive's own rarefy() + alpha_div(); the depth grid against NumPy; the
  # median against pandas (q2-diversity's summary). And rarefy(coverage =).
  # Fixtures: data-raw/expected/rarefaction_curve.r.

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))
  exit_if_not(requireNamespace("slam",   quietly = TRUE))

  fx <- readRDS('expected/rarefaction_curve.rds')
  ex_counts <- ecodive::ex_counts
  ex_tree   <- ecodive::ex_tree

  # Relative difference, against max(|y|, 1).
  rdiff <- function (x, y) max(abs(x - y) / pmax(abs(y), 1))

  # The rows of `rc` for (sample, depth) keys of `ref`.
  rows_of <- function (rc, ref, metric = 'observed') {
    r <- rc[rc$metric == metric, ]
    m <- match(paste(ref$sample, ref$depth), paste(r$sample, r$depth))
    stopifnot(!anyNA(m))
    r[m, ]
  }

  cpus <- 2L



  # Closed forms. ====

  for (nm in names(fx$tables)) {
    x  <- fx$tables[[nm]]
    ex <- fx$exact[[nm]]
    rc <- rarefaction_curve(x, depths = sort(unique(ex$depth)), cpus = cpus)
    r  <- rows_of(rc, ex)

    # The published formulas at 50 digits (mpmath).
    expect_true(rdiff(r$value,    ex$S)        < 1e-12, info = nm)
    expect_true(rdiff(r$sd,       ex$sd)       < 1e-12, info = nm)
    expect_true(rdiff(r$coverage, ex$coverage) < 1e-12, info = nm)

    # vegan::rarefy(se = TRUE). vegan takes differences of lchoose() of
    # large numbers, so its own error grows with N: up to 1.4e-12 in the
    # richness and 1e-8 (relative) in the standard error of the deep
    # samples, where Heck et al.'s pair terms cancel, against the mpmath
    # reference. ecodive's are within 1e-12 of it (above).
    v  <- fx$vegan[[nm]]
    r  <- rows_of(rc, v)
    deep <- nm == 'deep'
    expect_true(rdiff(r$value, v$S) < if (deep) 1e-11 else 1e-12, info = nm)
    big <- v$se > 1e-3
    expect_true(max(abs(r$sd - v$se)[big] / v$se[big]) < if (deep) 1e-7 else 1e-10, info = nm)
    expect_true(max(abs(r$sd - v$se)) < if (deep) 1e-5 else 1e-7, info = nm)

    # iNEXT's expected richness and sample coverage (SC). iNEXT uses
    # lgamma() differences too.
    i <- fx$inext[[nm]]
    r <- rows_of(rc, i)
    expect_true(rdiff(r$value,    i$qD) < if (deep) 1e-10 else 1e-12, info = nm)
    expect_true(rdiff(r$coverage, i$SC) < if (deep) 1e-9  else 1e-11, info = nm)
  }

  # A direct computation with lchoose() on small samples.
  direct <- function (x, d) {
    x <- x[x > 0]; N <- sum(x)
    a  <- function (k) ifelse(N - k >= d, exp(lchoose(N - k, d) - lchoose(N, d)), 0)
    S  <- length(x) - sum(a(x))
    V  <- sum(a(x) * (1 - a(x)))
    for (i in seq_along(x)) for (j in seq_along(x)) if (i < j)
      V <- V + 2 * (a(x[i] + x[j]) - a(x[i]) * a(x[j]))
    f1 <- sum(x == 1); f2 <- sum(x == 2)
    C  <- if (d < N) 1 - sum(x / N * ifelse(N - x >= d, exp(lchoose(N - x, d) - lchoose(N - 1, d)), 0))
          else if (f1 == 0) 1
          else if (f2 > 0) 1 - f1 / N * (N - 1) * f1 / ((N - 1) * f1 + 2 * f2)
          else 1 - f1 / N * (N - 1) * (f1 - 1) / ((N - 1) * (f1 - 1) + 2)
    c(S = S, sd = sqrt(max(V, 0)), C = C)
  }
  set.seed(3)
  small <- matrix(rpois(5 * 12, 3), 5, 12, dimnames = list(sprintf('s%d', 1:5), NULL))
  small <- rbind(small, fx$tables$edge[, rep(1:6, 2)])
  rc <- rarefaction_curve(small, depths = 1:80, cpus = cpus)
  dd <- t(mapply(function (s, d) direct(small[s, ], d), rc$sample, rc$depth))
  expect_true(rdiff(rc$value,    dd[, 'S'])  < 1e-12)
  expect_true(max(abs(rc$sd - dd[, "sd"]))   < 1e-6)   # lchoose()'s own cancellation (1.8e-7 where sd is 0)
  expect_true(rdiff(rc$sd[dd[, 'sd'] > 0.01], dd[dd[, 'sd'] > 0.01, 'sd']) < 1e-10)
  expect_true(rdiff(rc$coverage, dd[, 'C'])  < 1e-12)

  # At each sample's total: the observed richness, sd 0, and Eq. 12.
  rc  <- rarefaction_curve(ex_counts, cpus = cpus)
  tot <- attr(rc, 'totals')
  at  <- rc[rc$depth == tot[rc$sample], ]
  expect_identical(at$sample, names(tot))
  expect_identical(at$value, unname(as.numeric(rowSums(ex_counts > 0))))
  expect_identical(at$sd, rep(0, 4))
  expect_true(rdiff(at$coverage, sapply(names(tot), function (s) direct(ex_counts[s, ], tot[[s]])[['C']])) < 1e-14)

  # At depth 1: one feature, sd exactly 0, coverage Simpson's sum.
  d1 <- rc[rc$depth == 1, ]
  expect_identical(d1$value, rep(1, 4))
  expect_identical(d1$sd, rep(0, 4))
  expect_true(rdiff(d1$coverage, apply(ex_counts, 1, function (x) sum(x * (x - 1)) / (sum(x) * (sum(x) - 1)))) < 1e-14)



  # Result layout. ====

  rc <- rarefaction_curve(ex_counts, c('shannon', 'observed'), n_depths = 5, times = 2, cpus = cpus)
  expect_true(inherits(rc, 'ecodive_rarefaction_curve'))
  expect_true(inherits(rc, 'data.frame'))
  expect_identical(names(rc), c('sample', 'metric', 'depth', 'value', 'sd', 'coverage'))
  expect_identical(class(rc$depth), 'integer')
  expect_identical(attr(rc, 'totals'), rowSums(ex_counts))
  expect_identical(attr(rc, 'analytic'), 'observed')
  expect_identical(attr(rc, 'times'), 2L)
  expect_identical(attr(rc, 'average'), 'mean')
  expect_null(attr(rc, 'display'))
  # Ordered by sample (input order), metric (as given), depth.
  expect_identical(unique(rc$sample), rownames(ex_counts))
  expect_identical(unique(rc$metric), c('shannon', 'observed'))
  expect_identical(order(match(rc$sample, rownames(ex_counts)), match(rc$metric, c('shannon', 'observed')), rc$depth), seq_len(nrow(rc)))
  # The coverage doesn't depend on the metric.
  expect_identical(rc$coverage[rc$metric == 'shannon'], rc$coverage[rc$metric == 'observed'])
  # Duplicate metrics (by ID) are dropped; aliases match.
  rc2 <- rarefaction_curve(ex_counts, c('Shannon', 'obs', 'shannon'), n_depths = 5, times = 2, cpus = cpus)
  expect_identical(unclass(rc2)[names(rc2)], unclass(rc)[names(rc)])
  # Printing.
  expect_true(any(grepl('Rarefaction curves', capture.output(print(rc)))))
  # Unnamed samples are labelled by position.
  u  <- unname(ex_counts)
  ru <- rarefaction_curve(u, cpus = cpus)
  expect_identical(unique(ru$sample), c('1', '2', '3', '4'))
  expect_identical(names(attr(ru, 'totals')), c('1', '2', '3', '4'))
  expect_identical(ru$value, rarefaction_curve(ex_counts, cpus = cpus)$value)



  # The depth grid. ====

  for (k in seq_along(fx$grid)) {
    g <- ecodive:::rc_grid(NULL, fx$grid_cases$steps[k], NULL, c(a = fx$grid_cases$max_depth[k]))
    expect_identical(as.integer(g), unique(fx$grid[[k]]), info = k)
    g <- ecodive:::rc_grid(NULL, fx$grid_cases$steps[k], fx$grid_cases$max_depth[k], c(a = 5))
    expect_identical(as.integer(g), unique(fx$grid[[k]]), info = k)
  }

  # NumPy's grid for (1011, 10), plus each sample's total; no rows above it.
  rc <- rarefaction_curve(ex_counts, n_depths = 10, cpus = cpus)
  g  <- unique(fx$grid[[2]])
  for (s in rownames(ex_counts)) {
    N <- sum(ex_counts[s, ])
    expect_identical(rc$depth[rc$sample == s], sort(unique(as.integer(c(g[g <= N], N)))))
  }

  # depths: sorted, duplicates dropped, totals added up to the largest.
  rc <- rarefaction_curve(ex_counts, depths = c(700, 10, 10, 400), cpus = cpus)
  expect_identical(rc$depth[rc$sample == 'Saliva'], c(10L, 345L))
  expect_identical(rc$depth[rc$sample == 'Gums'],   c(10L, 400L, 700L))
  expect_identical(rc$depth[rc$sample == 'Stool'],  c(10L, 400L, 615L))
  expect_identical(rc$depth[rc$sample == 'Nose'],   c(10L, 400L, 700L))

  # max_depth caps the grid, and a total above it isn't added.
  rc <- rarefaction_curve(ex_counts, max_depth = 500, n_depths = 3, cpus = cpus)
  expect_identical(rc$depth[rc$sample == 'Saliva'], c(1L, 250L, 345L))
  expect_identical(rc$depth[rc$sample == 'Gums'],   c(1L, 250L, 500L))
  rc <- rarefaction_curve(ex_counts, depths = c(100, 900), max_depth = 600, cpus = cpus)
  expect_identical(rc$depth[rc$sample == 'Nose'],  100L)
  expect_identical(rc$depth[rc$sample == 'Stool'], 100L)
  expect_identical(rc$depth[rc$sample == 'Saliva'], c(100L, 345L))

  # No rows at all.
  r0 <- rarefaction_curve(matrix(0L, 3, 4), 'shannon', cpus = cpus)
  expect_identical(dim(r0), c(0L, 6L))
  expect_identical(names(r0), c('sample', 'metric', 'depth', 'value', 'sd', 'coverage'))

  # An empty sample has no rows, but has its total.
  z  <- rbind(ex_counts, Empty = 0L)
  rz <- rarefaction_curve(z, c('observed', 'shannon'), n_depths = 4, times = 2, cpus = cpus)
  expect_false('Empty' %in% rz$sample)
  expect_identical(attr(rz, 'totals')[['Empty']], 0)
  expect_identical(rz$value, rarefaction_curve(ex_counts, c('observed', 'shannon'), n_depths = 4, max_depth = 1011, times = 2, cpus = cpus)$value)



  # Resampling. ====

  # Every rarefaction t is rarefy(counts, d, seed = seed + t - 1); at a
  # sample's own total, rarefy() returns the sample itself.
  mets  <- c('observed', 'shannon', 'faith', 'hill')
  stack <- function (x, depths, seed, times, metric, ...) {
    tot <- rowSums(as.matrix(x))
    lapply(seq_len(times), function (t) {
      do.call(rbind, lapply(depths, function (d) {
        keep <- tot >= d
        if (!any(keep)) return (NULL)
        v <- alpha_div(rarefy(x, d, seed = seed + t - 1, cpus = cpus), metric, tree = ex_tree, ...)
        data.frame(sample = names(v), depth = as.integer(d), value = unname(v), stringsAsFactors = FALSE)
      }))
    })
  }

  dps <- c(1L, 30L, 345L, 400L, 886L, 1000L)
  rc1 <- rarefaction_curve(ex_counts, mets, depths = dps, times = 1, analytic = FALSE,
                           tree = ex_tree, q = 1, seed = 7, cpus = cpus)
  for (m in mets) {
    st <- stack(ex_counts, sort(unique(c(dps, rowSums(ex_counts)))), 7, 1, m, q = 1)[[1]]
    r  <- rows_of(rc1, st[paste(st$sample, st$depth) %in% paste(rc1$sample, rc1$depth), ], m)
    expect_identical(nrow(r), sum(rc1$metric == m), info = m)
    st <- st[match(paste(r$sample, r$depth), paste(st$sample, st$depth)), ]
    expect_identical(r$value, st$value, info = m)
    expect_identical(r$sd, rep(NA_real_, nrow(r)), info = m)
  }

  # The mean: the sequential sum over rarefactions, divided by times.
  k   <- 4L
  rcm <- rarefaction_curve(ex_counts, c('shannon', 'observed'), depths = dps, times = k,
                           analytic = FALSE, seed = 3, cpus = cpus)
  rcd <- rarefaction_curve(ex_counts, c('shannon', 'observed'), depths = dps, times = k,
                           analytic = FALSE, seed = 3, average = 'median', cpus = cpus)
  for (m in c('shannon', 'observed')) {
    sts  <- stack(ex_counts, sort(unique(c(dps, rowSums(ex_counts)))), 3, k, m)
    r    <- rcm[rcm$metric == m, ]
    keys <- paste(r$sample, r$depth)
    M    <- sapply(sts, function (st) st$value[match(keys, paste(st$sample, st$depth))])
    acc  <- M[, 1]
    for (t in 2:k) acc <- acc + M[, t]
    expect_identical(r$value, acc / k, info = m)
    expect_true(max(abs(r$sd - apply(M, 1, stats::sd))) < 1e-12, info = m)
    rd <- rcd[rcd$metric == m, ]
    expect_identical(rd$value, apply(M, 1, stats::median), info = m)
    expect_true(max(abs(rd$sd - apply(M, 1, stats::sd))) < 1e-12, info = m)
  }
  expect_identical(attr(rcd, 'average'), 'median')
  expect_identical(attr(rcm, 'analytic'), character(0))

  # q2-diversity's per-sample summary: pandas' describe() 50% percentile.
  rq <- rarefaction_curve(fx$q2_table, 'shannon', depths = fx$q2_depths, max_depth = max(fx$q2_depths),
                          times = 7, average = 'med', seed = 0, cpus = cpus)
  r  <- rq[match(names(fx$q2_median), paste(rq$sample, rq$depth)), ]
  expect_identical(r$value, unname(fx$q2_median))

  # Nested draws: with the same seed, a smaller depth keeps a subset.
  x  <- matrix(rpois(30 * 40, rep(rexp(40, 1 / 8), each = 30)), 30, 40)
  x[, 1] <- x[, 1] + 50L
  tm <- function (m) as.matrix(m)
  for (cls in c('matrix', 'dgCMatrix', 'dgTMatrix', 'dgeMatrix', 'slam')) {
    for (margin in 1:2) {
      xx <- if (margin == 1) x else t(x)
      y  <- switch(
        cls,
        matrix    = xx,
        dgCMatrix = methods::as(Matrix::Matrix(xx, sparse = TRUE), 'CsparseMatrix'),
        dgTMatrix = methods::as(Matrix::Matrix(xx, sparse = TRUE), 'TsparseMatrix'),
        dgeMatrix = Matrix::Matrix(xx + 0, sparse = FALSE),
        slam      = slam::as.simple_triplet_matrix(xx) )
      for (s in c(0L, 9L)) {
        a <- tm(rarefy(y, 40,  seed = s, margin = margin, drop = FALSE))
        b <- tm(rarefy(y, 120, seed = s, margin = margin, drop = FALSE))
        expect_true(all(a <= b), info = paste(cls, margin, s))
        # And the same draws as the base matrix.
        expect_identical(unname(a), unname(tm(rarefy(xx, 40, seed = s, margin = margin, drop = FALSE))), info = cls)
      }
    }
  }

  # So each rarefaction's richness never decreases with depth.
  for (s in 0:4) {
    r <- rarefaction_curve(x, depths = seq(1, 300, by = 7), times = 1, analytic = FALSE, seed = s, cpus = cpus)
    expect_true(all(tapply(r$value, r$sample, function (v) all(diff(v) >= 0))), info = s)
  }

  # The resampled richness agrees with the closed form within 4 Monte Carlo
  # standard errors (exactly at depth 1 and at each total).
  ra <- rarefaction_curve(ex_counts, depths = c(1, 5, 20, 100, 300), times = 2000, analytic = FALSE, cpus = cpus)
  rc <- rarefaction_curve(ex_counts, depths = c(1, 5, 20, 100, 300), cpus = cpus)
  se <- ra$sd / sqrt(2000)
  expect_true(all(abs(ra$value - rc$value) <= 4 * se))
  expect_identical(ra$value[ra$sd == 0], rc$value[ra$sd == 0])
  expect_true(all(abs(ra$sd - rc$sd) <= 0.1 * rc$sd + 1e-12))



  # Margins, classes and threads. ====

  set.seed(38)
  big <- matrix(rpois(120 * 30, rep(rexp(30, 1 / 5), each = 120)), 120, 30,
                dimnames = list(sprintf('S%03d', 1:120), sprintf('F%02d', 1:30)))
  big[1:3, ] <- 0L; big[1, 1] <- 1L; big[2, 1:2] <- 1L   # tiny samples
  ref <- rarefaction_curve(big, c('observed', 'simpson'), n_depths = 7, times = 3,
                           seed = 2, cpus = 1L)
  for (cp in c(2L, 4L))
    expect_identical(rarefaction_curve(big, c('observed', 'simpson'), n_depths = 7,
                                       times = 3, seed = 2, cpus = cp), ref, info = cp)
  rfa <- rarefaction_curve(big, 'observed', n_depths = 7, times = 3, analytic = FALSE, seed = 2, cpus = 1L)
  expect_identical(rarefaction_curve(big, 'observed', n_depths = 7, times = 3, analytic = FALSE, seed = 2, cpus = 4L), rfa)

  same <- function (y, ...) {
    r <- rarefaction_curve(y, c('observed', 'simpson'), n_depths = 7, times = 3, seed = 2, cpus = cpus, ...)
    expect_identical(unclass(r)[names(r)], unclass(ref)[names(ref)])
    expect_identical(attr(r, 'totals'), attr(ref, 'totals'))
  }
  same(t(big), margin = 2L)
  same(Matrix::Matrix(big, sparse = TRUE))
  same(Matrix::t(Matrix::Matrix(big, sparse = TRUE)), margin = 2L)
  same(slam::as.simple_triplet_matrix(big))
  same(big + 0)

  # Arguments reach the metrics.
  r1 <- rarefaction_curve(ex_counts, 'hill', q = 0, n_depths = 4, times = 2, cpus = cpus)
  r2 <- rarefaction_curve(ex_counts, 'observed', n_depths = 4, times = 2, analytic = FALSE, cpus = cpus)
  expect_identical(r1$value, r2$value)
  r1 <- rarefaction_curve(ex_counts, 'chao1', bias_corrected = FALSE, n_depths = 4, times = 1, cpus = cpus)
  r2 <- rarefaction_curve(ex_counts, 'chao1', n_depths = 4, times = 1, cpus = cpus)
  expect_false(identical(r1$value, r2$value))



  # Errors. ====

  expect_error(rarefaction_curve(ex_counts, 'nope'))
  expect_error(rarefaction_curve(ex_counts, character(0)), '`metric` must be one or more')
  expect_error(rarefaction_curve(ex_counts, NA_character_), '`metric` must be one or more')
  expect_error(rarefaction_curve(ex_counts + 0.5))
  expect_error(rarefaction_curve(-ex_counts), 'negative')
  expect_error(rarefaction_curve(ex_counts, depths = c(0, 5)), '`depths` must be NULL or whole numbers')
  expect_error(rarefaction_curve(ex_counts, depths = c(NA, 5)), '`depths` must be NULL or whole numbers')
  expect_error(rarefaction_curve(ex_counts, depths = 1.5), '`depths` must be NULL or whole numbers')
  expect_error(rarefaction_curve(ex_counts, depths = 'a'), '`depths` must be NULL or whole numbers')
  expect_error(rarefaction_curve(ex_counts, n_depths = 1), '`n_depths` must be a whole number of at least 2')
  expect_error(rarefaction_curve(ex_counts, n_depths = 2.5), '`n_depths` must be a whole number of at least 2')
  expect_error(rarefaction_curve(ex_counts, max_depth = 0), '`max_depth` must be NULL or a whole number')
  expect_error(rarefaction_curve(ex_counts, max_depth = c(5, 6)), '`max_depth` must be NULL or a whole number')
  expect_error(rarefaction_curve(ex_counts, 'shannon', analytic = TRUE), 'only "observed" has; not "shannon"')
  expect_error(rarefaction_curve(ex_counts, analytic = NA), '`analytic` must be NULL, TRUE or FALSE')
  expect_error(rarefaction_curve(ex_counts, average = 'mode'), '`average = "mode"` must be one of')
  expect_error(rarefaction_curve(ex_counts, average = 'medoid'), 'only for beta diversity')
  expect_error(rarefaction_curve(ex_counts, times = 0), '`times` must be a whole number of at least 1')
  expect_error(rarefaction_curve(ex_counts, 'faith', n_depths = 3, times = 1))
  expect_identical(rarefaction_curve(ex_counts, analytic = TRUE), rarefaction_curve(ex_counts))



  # rarefy(coverage =). ====

  # Each sample's depth is the smallest whose coverage reaches the target.
  rc  <- rarefaction_curve(ex_counts, depths = seq_len(1011), cpus = cpus)
  for (cv in c(0.5, 0.9, 0.995, 0.9968)) {
    want <- sapply(rownames(ex_counts), function (s) {
      r <- rc[rc$sample == s, ]
      i <- which(r$coverage >= cv)
      if (length(i)) r$depth[min(i)] else NA_integer_
    })
    want <- want[!is.na(want)]
    got  <- rarefy(ex_counts, coverage = cv, seed = 4, cpus = cpus)
    expect_identical(rownames(got), names(want), info = cv)
    expect_identical(rowSums(got), stats::setNames(as.numeric(want), names(want)), info = cv)
  }

  # The same draws as rarefying to those depths.
  cv   <- 0.995
  got  <- rarefy(ex_counts, coverage = cv, seed = 4, cpus = cpus)
  dfull <- rowSums(ex_counts) + 1
  dfull[rownames(got)] <- rowSums(got)
  expect_identical(got, rarefy(ex_counts, depth = dfull, seed = 4, cpus = cpus))

  # Samples that can't reach it: dropped, or returned unrarefied.
  cv <- 0.9968
  kept <- rownames(rarefy(ex_counts, coverage = cv, cpus = cpus))
  lost <- setdiff(rownames(ex_counts), kept)
  expect_true(length(lost) > 0)
  nd <- rarefy(ex_counts, coverage = cv, drop = FALSE, cpus = cpus)
  expect_identical(nd[lost, ], ex_counts[lost, ] + 0)
  expect_warning(rarefy(ex_counts, coverage = cv, warn = TRUE, cpus = cpus), 'estimated coverage below 0.9968 and will be dropped')
  expect_warning(rarefy(ex_counts, coverage = cv, warn = TRUE, drop = FALSE, cpus = cpus), 'returned unrarefied')

  # Classes, margins, threads and times.
  ref <- rarefy(big, coverage = 0.8, seed = 1, cpus = 1L)
  expect_true(nrow(ref) < 120)
  for (cp in c(2L, 4L)) expect_identical(rarefy(big, coverage = 0.8, seed = 1, cpus = cp), ref)
  expect_identical(t(rarefy(t(big), coverage = 0.8, seed = 1, margin = 2L, cpus = cpus)), ref)
  expect_identical(as.matrix(rarefy(Matrix::Matrix(big, sparse = TRUE), coverage = 0.8, seed = 1, cpus = cpus)), ref)
  expect_identical(as.matrix(rarefy(slam::as.simple_triplet_matrix(big), coverage = 0.8, seed = 1, cpus = cpus)), ref)
  expect_identical(rarefy(big, coverage = 0.8, seed = 1, times = 2, cpus = cpus)[[1]], ref)

  # Errors.
  expect_error(rarefy(ex_counts, depth = 100, coverage = 0.9), 'not both')
  expect_error(rarefy(ex_counts, coverage = 0),        '`coverage` must be NULL or a single number')
  expect_error(rarefy(ex_counts, coverage = 1.1),      '`coverage` must be NULL or a single number')
  expect_error(rarefy(ex_counts, coverage = NA_real_), '`coverage` must be NULL or a single number')
  expect_error(rarefy(ex_counts, coverage = c(.5, .6)), '`coverage` must be NULL or a single number')
  expect_error(rarefy(ex_counts, coverage = '0.5'),    '`coverage` must be NULL or a single number')



  # plot(). ====

  grDevices::pdf(NULL)
  invisibly <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  rc <- rarefaction_curve(ex_counts, c('observed', 'shannon'), n_depths = 6, times = 3, cpus = cpus)
  grp <- c(Saliva = 'oral', Gums = 'oral', Nose = 'skin', Stool = 'gut')
  meta <- data.frame(site = grp, row.names = names(grp))
  expect_true(invisibly(rc))
  expect_true(invisibly(rc, 'shannon'))
  expect_true(invisibly(rc, 'Shannon Diversity Index', type = 'cov'))
  expect_true(invisibly(rc, type = 'completeness', depth = 500))
  expect_true(invisibly(rc, depth = 500, col = 'blue', lwd = 2))
  expect_true(invisibly(rc, groups = grp))
  expect_true(invisibly(rc, groups = unname(grp), legend = NULL))
  expect_true(invisibly(rc, 'shannon', groups = 'site', data = meta, type = 'coverage'))
  expect_true(invisibly(rc, groups = 'site', data = meta, col = c(oral = 'red', skin = 'blue', gut = 'green')))
  expect_true(invisibly(rc, groups = c(1, 1, 2, 3)))
  expect_true(invisibly(rc, main = 'Mine', xlab = 'Reads', log = 'x'))
  expect_true(invisibly(ru, groups = c('a', 'a', 'b', 'b')))
  expect_error(plot(rc, 'simpson'), 'is not in `x`')
  expect_error(plot(rc, type = 'co'), '`type` must be one of')
  expect_error(plot(rc, type = 'coverage', depth = 10), 'use type = "depth" or "completeness"')
  expect_error(plot(rc, depth = -1), '`depth` must be a single number')
  expect_error(plot(rc, groups = c('a', 'b')))
  expect_error(plot(rc[0, ]), 'nothing to draw')

  # Titles and labels.
  lab <- ecodive:::rc_plot_labels(rc, 'observed', 'depth')
  expect_identical(lab$main, 'Rarefaction Curves of Observed Features')
  expect_identical(lab$sub,  'Expected values')
  expect_identical(lab$ylab, 'Observed Features')
  lab <- ecodive:::rc_plot_labels(rc, 'shannon', 'coverage')
  expect_identical(lab$main, 'Rarefaction Curves of Shannon Diversity Index')
  expect_identical(lab$sub,  'Mean of 3 rarefactions')
  expect_identical(lab$xlab, 'Sample Coverage')
  lab <- ecodive:::rc_plot_labels(rcd, 'shannon', 'completeness')
  expect_identical(lab$main, 'Sample Completeness Curves')
  expect_identical(lab$ylab, 'Sample Coverage')
  expect_identical(ecodive:::rc_plot_labels(rcd, 'shannon', 'depth')$sub, 'Median of 4 rarefactions')
  expect_identical(ecodive:::rc_plot_labels(rc1, 'shannon', 'depth')$sub, 'Mean of 1 rarefaction')
  grDevices::dev.off()

#})
