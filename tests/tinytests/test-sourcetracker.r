#test_that("sourcetracker", {

  # sourcetracker(): SourceTracker2 2.0.1-dev's collapse_source_data(),
  # ConditionalProbability and draw schedule on fixed inputs (exact), and
  # its gibbs() over 40 seeds on simulated mixtures (by distribution; see
  # data-raw/expected/sourcetracker.py); rarefaction as rarefy(); threads,
  # input types, edge cases and argument checks.

  ref <- readRDS(file.path('expected', 'sourcetracker.rds'))

  maxrel <- function (a, b) {
    a <- as.vector(a); b <- as.vector(b)
    max(c(0, ifelse(b == 0, abs(a), abs(a - b) / abs(b))))
  }


  # -----------------------------------------------------------------------
  # Deterministic parts, against SourceTracker2
  # -----------------------------------------------------------------------

  # collapse_source_data(..., 'mean'), truncated, environments in C-locale
  # order ('C' < 'a' < 'b').
  cl <- ref$collapse
  d  <- sourcetracker_details(
    rbind(cl$table, 1), c(cl$env, NA), depth = NULL, n_burnin = 0L, n_restarts = 1L, cpus = 1L )
  expect_identical(rownames(d$collapsed), rownames(cl$expected))
  expect_identical(unname(d$collapsed), unname(cl$expected))
  expect_identical(d$sources, d$collapsed)   # no rarefaction
  expect_identical(colnames(d$result), c(rownames(cl$expected), 'Unknown'))

  # ConditionalProbability.calculate_cp_slice().
  for (i in seq_along(ref$cond_prob$cases)) {
    r   <- ref$cond_prob$cases[[i]]
    got <- .Call(
      C_st_cond_prob, ref$cond_prob$source, r$alpha1, r$alpha2, r$beta, r$n,
      r$xi, r$m_xiV, r$m_V, r$n_vnoti )
    expect_true(maxrel(got, r$expected) < 1e-10, info = paste('cond_prob case', i))
  }

  # gibbs_sampler()'s schedule: the passes followed by a draw.
  for (s in ref$schedule) {
    expect_identical(which(st_draw_passes(s$burnin, s$draws, s$delay) == 1L), s$passes)
    # ... and the sampler takes exactly those draws.
    d <- sourcetracker_details(
      rbind(a = c(5, 1), b = c(1, 5), s = c(3, 3)), c('a', 'b', NA), depth = NULL,
      n_restarts = 2L, n_draws = s$draws, n_burnin = s$burnin, thin = s$delay, cpus = 1L )
    expect_identical(dim(d$draws), c(3L, 2L * s$draws, 1L))
    expect_identical(attr(d$result, 'n_draws'), 2L * s$draws)
    expect_true(all(colSums(d$draws[, , 1]) == 6L))
  }


  # -----------------------------------------------------------------------
  # The whole sampler, by distribution: SourceTracker2's gibbs() over 40
  # seeds against ecodive over 20, on the same (already rarefied) tables.
  # -----------------------------------------------------------------------

  sim    <- ref$sim
  envs   <- rownames(sim$sources)
  counts <- rbind(sim$sources, sim$sinks)
  # In SourceTracker2's order (gibbs() takes the sources as given).
  labs   <- factor(c(envs, rep(NA, nrow(sim$sinks))), levels = envs)

  for (run in names(sim$runs)) {
    r     <- sim$runs[[run]]
    seeds <- 1:20
    fits  <- lapply(seeds, function (seed) sourcetracker(
      counts, labs, depth = NULL, n_restarts = r$n_restarts,
      n_draws = r$n_draws, n_burnin = r$n_burnin, thin = r$thin,
      loo = r$loo, features = TRUE, seed = seed ))
    P  <- simplify2array(lapply(fits, unclass))          # sinks x V x seeds
    m  <- apply(P, 1:2, mean)
    s  <- apply(P, 1:2, stats::sd)
    se <- sqrt(r$sd^2 / r$n_seeds + s^2 / length(seeds))
    expect_identical(dimnames(m), dimnames(r$mean), info = run)
    z  <- abs(m - r$mean) / pmax(se, 1e-3)
    expect_true(max(z) < 4.5, info = paste(run, 'means: max z', max(z)))
    expect_true(max(abs(m - r$mean)) < 0.01, info = paste(run, 'means'))

    # Spread between seeds, where there is any.
    big <- r$sd > 0.01
    ratio <- s[big] / r$sd[big]
    expect_true(all(ratio > 0.4 & ratio < 2.5), info = paste(run, 'between-seed sd'))

    # The sd attribute: SourceTracker2's mps times its number of draws.
    D   <- r$n_restarts * r$n_draws
    sds <- apply(simplify2array(lapply(fits, attr, 'sd')), 1:2, mean)
    big <- r$mps * D > 0.01
    ratio <- sds[big] / (r$mps[big] * D)
    expect_true(all(ratio > 0.75 & ratio < 1.33), info = paste(run, 'sd attribute'))

    # Feature tables: SourceTracker2's sums over draws are D times ecodive's
    # means. Each feature's reads are all assigned, every draw.
    ft  <- attr(fits[[1]], 'features')
    sk  <- if (r$loo) sim$sources else sim$sinks
    for (j in names(ft)) {
      expect_equal(colSums(ft[[j]]), sk[j, ], info = paste(run, j))
      expect_equal(rowSums(ft[[j]]) / sum(sk[j, ]), fits[[1]][j, ], info = paste(run, j))
      expect_identical(dimnames(ft[[j]]), dimnames(r$fas[[j]]), info = paste(run, j))
      expect_equal(colSums(r$fas[[j]]) / D, sk[j, ], info = paste(run, j))
    }
    Fm <- Reduce(`+`, lapply(fits, function (f) do.call(rbind, attr(f, 'features')))) / length(fits)
    Fr <- do.call(rbind, r$fas) / D
    expect_true(max(abs(Fm - Fr)) < 0.1 * max(sk), info = paste(run, 'features'))
    expect_true(stats::cor(as.vector(Fm), as.vector(Fr)) > 0.99, info = paste(run, 'features'))
  }

  # The true mixtures: SourceTracker is biased where sources overlap (it
  # underestimates the unknown, as Shenhav et al. 2019 found), so both are
  # up to about 0.1 off; ecodive is no farther off than SourceTracker2, and
  # pure sinks are recovered.
  r  <- sim$runs$default
  st <- sourcetracker(counts, labs, depth = NULL, n_restarts = 20L, seed = 1L)
  expect_true(max(abs(st - sim$truth)) < max(abs(r$mean - sim$truth)) + 0.02)
  expect_true(st['sink4', 'soil'] > 0.9 && st['sink5', 'Unknown'] > 0.85)

  # SourceTracker 1 (GPL, outputs only), when its run is in the fixture:
  # loosely.
  if (!is.null(sim$v1))
    expect_true(max(abs(r$mean - sim$v1$mean)) < 0.05)


  # -----------------------------------------------------------------------
  # Rarefaction: each sink as rarefy() rarefies it; each environment
  # rarefied from its pooled counts.
  # -----------------------------------------------------------------------

  x   <- rbind(sim$sources[c(1, 1, 2, 3), ], sim$sinks)
  lab <- factor(c(envs[c(1, 1, 2, 3)], rep(NA, 6)), levels = envs)
  d   <- sourcetracker_details(x, lab, depth = 1500L, sink_depth = 300L, n_burnin = 2L, seed = 7L)
  expect_identical(unname(d$sinks), unname(rarefy(x, 300L, seed = 7L)[5:10, ]))
  expect_true(all(rowSums(d$sources) == 1500))
  expect_true(all(d$sources <= d$collapsed))
  expect_identical(unname(d$collapsed), unname(sim$sources))
  d2 <- sourcetracker_details(x, lab, depth = 1500L, sink_depth = 300L, n_burnin = 2L, seed = 8L)
  expect_false(identical(d$sources, d2$sources))
  expect_false(identical(d$sinks, d2$sinks))
  # A sink as deep as sink_depth isn't rarefied; depth = NULL rarefies none.
  d3 <- sourcetracker_details(x, lab, depth = NULL, sink_depth = NULL, n_burnin = 2L)
  expect_identical(unname(d3$sinks), unname(x[5:10, ]))
  expect_identical(d3$sources, d3$collapsed)
  # Under loo, the sinks are the rarefied environments.
  d4 <- sourcetracker_details(x, lab, depth = 1500L, n_burnin = 2L, loo = TRUE, seed = 7L)
  expect_identical(d4$sinks, d4$sources)
  expect_identical(d4$sources, d$sources)


  # -----------------------------------------------------------------------
  # Determinism: any cpus, any input type, a sink's result without the
  # sinks after it.
  # -----------------------------------------------------------------------

  big <- do.call(rbind, rep(list(sim$sinks), 4))
  rownames(big) <- paste0('k', seq_len(nrow(big)))
  xb  <- rbind(sim$sources, big)
  lb  <- factor(c(envs, rep(NA, nrow(big))), levels = envs)
  run <- function (x = xb, cpus = 1L, features = TRUE, ...)
    sourcetracker(x, lb[seq_len(nrow(x))], depth = 1000L, sink_depth = 200L,
                  n_restarts = 3L, n_burnin = 10L, features = features, cpus = cpus, ...)
  a <- run()
  for (cpus in c(2L, 4L, 7L))
    expect_identical(run(cpus = cpus), a, info = paste('cpus', cpus))
  expect_identical(run(), a)
  expect_false(identical(run(seed = 1L), a))
  al <- run(loo = TRUE)
  for (cpus in c(2L, 4L))
    expect_identical(run(cpus = cpus, loo = TRUE), al, info = paste('loo cpus', cpus))

  # The first sinks alone: the same rows.
  part <- run(xb[1:10, ])
  expect_equal(unclass(part)[, ], unclass(a)[1:7, ], check.attributes = FALSE)
  expect_identical(attr(part, 'features'), attr(a, 'features')[1:7])
  expect_identical(attr(part, 'sd'), attr(a, 'sd')[1:7, ])

  # Sparse, transposed and integer tables.
  if (requireNamespace('Matrix', quietly = TRUE)) {
    expect_identical(run(Matrix::Matrix(xb, sparse = TRUE)), a)
    expect_identical(run(methods::as(Matrix::Matrix(xb, sparse = TRUE), 'TsparseMatrix')), a)
  }
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(run(slam::as.simple_triplet_matrix(xb)), a)
  expect_identical(sourcetracker(
    t(xb), lb, margin = 2L, depth = 1000L, sink_depth = 200L, n_restarts = 3L,
    n_burnin = 10L, features = TRUE, cpus = 1L ), a)
  xi <- xb; storage.mode(xi) <- 'integer'
  expect_identical(run(xi), a)


  # -----------------------------------------------------------------------
  # The result.
  # -----------------------------------------------------------------------

  expect_true(inherits(a, 'ecodive_source') && is.matrix(a) && is.numeric(a))
  expect_identical(dimnames(a), list(rownames(big), c(envs, 'Unknown')))
  expect_equal(unname(rowSums(a)), rep(1, nrow(big)))
  expect_identical(dim(attr(a, 'sd')), dim(a))
  expect_identical(dimnames(attr(a, 'sd')), dimnames(a))
  expect_identical(attr(a, 'n_draws'), 3L)
  expect_identical(attr(a, 'method'), 'sourcetracker')
  expect_true(all(attr(a, 'sd') >= 0))
  expect_identical(names(attr(a, 'features')), rownames(big))
  expect_identical(dimnames(attr(a, 'features')[[1]]), list(c(envs, 'Unknown'), colnames(xb)))
  expect_equal(colSums(attr(a, 'features')[[1]]), rarefy(xb, 200L)[4, ], check.attributes = FALSE)
  expect_null(attr(run(features = FALSE), 'features'))
  expect_identical(class(unclass(a)[, 1]), 'numeric')

  # Leave one out: environments x (environments, Unknown), 0 on the diagonal.
  expect_identical(dimnames(al), list(envs, c(envs, 'Unknown')))
  expect_identical(unname(diag(unclass(al)[, 1:3])), c(0, 0, 0))
  expect_identical(unname(diag(attr(al, 'sd')[, 1:3])), c(0, 0, 0))
  expect_equal(unname(rowSums(al)), rep(1, 3))
  expect_identical(unname(attr(al, 'features')[['soil']]['soil', ]), rep(0, ncol(xb)))

  # print() and plot().
  out <- utils::capture.output(print(a))
  expect_true(any(grepl('Source contributions (SourceTracker)', out, fixed = TRUE)))
  expect_true(any(grepl('Mean of 3 draws', out, fixed = TRUE)))
  expect_true(any(grepl('k24', out, fixed = TRUE)) || any(grepl('more sinks', out, fixed = TRUE)))
  grDevices::pdf(NULL)
  v <- withVisible(plot(a))
  expect_true(identical(v$value, a) && !v$visible)
  expect_true(identical(withVisible(plot(al, horiz = TRUE, legend = 'topright'))$visible, FALSE))
  expect_true(identical(withVisible(plot(a, legend = NULL, main = 'x', col = 1:4))$visible, FALSE))
  expect_error(plot(a, legend = 'nowhere'), '`legend`')
  grDevices::dev.off()
  # The title comes from the method's display name.
  expect_identical(source_tracker_name('sourcetracker'), 'SourceTracker')
  expect_null(source_tracker_name('nope'))


  # -----------------------------------------------------------------------
  # Edge cases.
  # -----------------------------------------------------------------------

  small <- rbind(
    a1 = c(50, 30, 20, 0, 0), a2 = c(40, 40, 20, 0, 0),
    b1 = c(0, 10, 20, 70, 0), b2 = c(0, 20, 10, 70, 0),
    s1 = c(25, 20, 20, 35, 0), s2 = c(0, 0, 0, 0, 100) )
  sl <- c('a', 'a', 'b', 'b', NA, NA)
  fit <- function (n_restarts = 4L, n_burnin = 20L, ...)
    sourcetracker(small, sl, depth = NULL, n_burnin = n_burnin, n_restarts = n_restarts, ...)

  # A sink with all its reads in a feature no source has: mostly unknown.
  f <- fit()
  expect_true(f['s2', 'Unknown'] > 0.95)
  expect_true(f['s1', 'Unknown'] < 0.5)

  # One environment.
  f1 <- sourcetracker(small[c(1, 2, 5, 6), ], c('a', 'a', NA, NA), depth = NULL, n_burnin = 5L)
  expect_identical(colnames(f1), c('a', 'Unknown'))
  expect_equal(unname(rowSums(f1)), c(1, 1))
  f1l <- sourcetracker(small[1:2, ], c('a', 'a'), depth = NULL, n_burnin = 5L, loo = TRUE)
  expect_identical(unname(unclass(f1l)[1, ]), c(0, 1))

  # A sample named in data; a factor's levels set the order.
  meta <- data.frame(env = sl, row.names = rownames(small))
  expect_identical(sourcetracker(small, 'env', meta, depth = NULL, n_burnin = 20L, n_restarts = 4L), f)
  ff <- fit(); expect_identical(colnames(ff), c('a', 'b', 'Unknown'))
  fr <- sourcetracker(small, factor(sl, levels = c('b', 'a')), depth = NULL, n_burnin = 20L, n_restarts = 4L)
  expect_identical(colnames(fr), c('b', 'a', 'Unknown'))
  expect_equal(unclass(fr)[, c('a', 'b', 'Unknown')], unclass(f)[, ], tolerance = 0.2, check.attributes = FALSE)
  # Names: matched by name when both are named.
  expect_identical(sourcetracker(small, setNames(rev(sl), rev(rownames(small))), depth = NULL, n_burnin = 20L, n_restarts = 4L), f)
  # Unnamed tables: rows named by position.
  fu <- sourcetracker(unname(small), sl, depth = NULL, n_burnin = 20L, n_restarts = 4L)
  expect_identical(rownames(fu), c('5', '6'))
  # Logical input.
  expect_silent(sourcetracker(small > 10, sl, depth = NULL, n_burnin = 2L))
  # loo ignores the sinks; with no sinks, loo still works.
  expect_identical(fit(loo = TRUE), sourcetracker(small[1:4, ], sl[1:4], depth = NULL, n_burnin = 20L, n_restarts = 4L, loo = TRUE))
  # alpha1 = 0 works while every environment has reads.
  expect_true(inherits(fit(alpha1 = 0), 'ecodive_source'))
  # Rarefied sources: depth equal to an environment's depth is fine.
  expect_true(inherits(sourcetracker(small, sl, depth = 100L, sink_depth = 100L, n_burnin = 2L), 'ecodive_source'))


  # -----------------------------------------------------------------------
  # Errors.
  # -----------------------------------------------------------------------

  expect_error(sourcetracker(small, sl, depth = 101L, sink_depth = NULL), '2 source environments have fewer than `depth` \\(101\\)')
  expect_error(sourcetracker(small, sl, depth = NULL, sink_depth = 101L), "2 sinks have fewer than `sink_depth` \\(101\\) reads: 's1', 's2'")
  expect_error(sourcetracker(rbind(small, s3 = 0), c(sl, NA), depth = NULL), "1 sink has no reads: 's3'")
  expect_error(sourcetracker(small, rep('a', 6), depth = NULL), 'no sinks')
  expect_error(sourcetracker(small, rep(NA, 6), depth = NULL), 'no source samples')
  expect_error(sourcetracker(small, depth = NULL), '`sources` is required')
  expect_error(sourcetracker(small, sl[1:5], depth = NULL), '`sources`')
  expect_error(sourcetracker(small, 'env', depth = NULL), 'names a column')
  expect_error(sourcetracker(small, c('Unknown', sl[-1]), depth = NULL), 'Unknown')
  expect_error(sourcetracker(small + 0.5, sl, depth = NULL), 'whole numbers')
  expect_error(sourcetracker(-small, sl, depth = NULL), 'negative')
  small_na <- small; small_na[1, 1] <- NA
  expect_error(sourcetracker(small_na, sl, depth = NULL), 'NA')
  zero_env <- rbind(small, z1 = 0)
  expect_error(sourcetracker(zero_env, c(sl, 'z'), depth = NULL, alpha1 = 0), "'z' has no reads")
  expect_error(sourcetracker(zero_env, c(sl, 'z'), depth = NULL, loo = TRUE), "'z' has no reads")
  expect_true(inherits(sourcetracker(zero_env, c(sl, 'z'), depth = NULL, n_burnin = 2L), 'ecodive_source'))
  for (bad in list(0, -1, 1.5, NA, 'a', c(1, 2)))
    expect_error(fit(n_restarts = bad), '`n_restarts`')
  expect_error(fit(n_restarts = 2^24), '`n_restarts`')
  expect_error(fit(n_draws = 0), '`n_draws`')
  expect_error(fit(n_burnin = -1), '`n_burnin`')
  expect_error(fit(thin = 0), '`thin`')
  expect_error(fit(alpha1 = -1), '`alpha1`')
  expect_error(fit(alpha2 = 0), '`alpha2`')
  expect_error(fit(beta = 0), '`beta`')
  expect_error(fit(beta = Inf), '`beta`')
  expect_error(fit(loo = NA), '`loo`')
  expect_error(fit(features = 'yes'), '`features`')
  expect_error(sourcetracker(small, sl, depth = 0), '`depth`')
  expect_error(sourcetracker(small, sl, depth = NULL, sink_depth = c(1, 2)), '`sink_depth`')
  expect_error(fit(seed = 1.5), '`seed`')
  expect_error(fit(cpus = 0), '`cpus`')
  expect_error(sourcetracker(small, sl, data = 1:6, depth = NULL), '`data`')

#})
