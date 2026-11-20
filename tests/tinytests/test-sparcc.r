#test_that("sparcc", {

  # sparcc(): the deterministic steps (given log fractions) against a base-R
  # implementation of Friedman & Alm (2012) (data-raw/expected/sparcc.r),
  # dense solves of the variance system, invariances, threading, the
  # Dirichlet draws against FastSpar 1.0.0 and SpiecEasi 1.1.1 (GPL; their
  # outputs only, data-raw/expected/sparcc/), p-values, clustering and
  # argument checks.

  ref <- readRDS(file.path('expected', 'sparcc.rds'))

  plain  <- function (d) as.vector(unclass(d))
  maxabs <- function (a, b) max(abs(a - b), na.rm = TRUE)

  # Packed pairs (i < j), in a dist's order.
  pairs_of <- function (D) {
    cbind(
      i = unlist(lapply(seq_len(D - 1L), function (i) rep(i, D - i))),
      j = unlist(lapply(seq_len(D - 1L), function (i) (i + 1L):D)) )
  }
  variation <- function (L) {
    pr <- pairs_of(ncol(L))
    vapply(seq_len(nrow(pr)), function (k) var(L[, pr[k, 1]] - L[, pr[k, 2]]), numeric(1))
  }

  # Basis variances by a dense solve, with the packed pairs `ex` excluded and
  # only the `active` features.
  dense_w <- function (t, D, ex, active) {
    pr   <- pairs_of(D)
    keep <- matrix(TRUE, D, D)
    keep[pr[ex, , drop = FALSE]] <- FALSE
    keep[pr[ex, 2:1, drop = FALSE]] <- FALSE
    diag(keep) <- FALSE
    T <- matrix(0, D, D)
    T[pr] <- t; T[pr[, 2:1]] <- t
    idx <- which(active)
    K <- keep[idx, idx]
    w <- rep(NA_real_, D)
    w[idx] <- solve(diag(rowSums(K), length(idx)) + K * 1, rowSums(T[idx, idx] * K))
    w
  }
  rho_from <- function (w, t, D) {
    pr <- pairs_of(D)
    wi <- w[pr[, 1]]; wj <- w[pr[, 2]]
    r  <- suppressWarnings((wi + wj - t) / (2 * sqrt(wi) * sqrt(wj)))
    r[!is.na(wi) & !is.na(wj) & !(wi > 0 & wj > 0)] <- NaN
    r
  }


  # -----------------------------------------------------------------------
  # The deterministic steps, against the base-R reference.
  # -----------------------------------------------------------------------

  for (nm in setdiff(names(ref), 'refs')) {
    cc <- ref[[nm]]
    ex <- cc$expected
    for (cpus in c(1L, 3L)) {
      got <- suppressWarnings(sparcc_details(
        cc$counts, log_fractions = cc$draws, max_exclusions = cc$max_exclusions,
        exclusion_threshold = cc$threshold, cpus = cpus ))
      info <- paste(nm, cpus)
      expect_identical(is.na(got$draw_rho), is.na(ex$draw_rho), info = info)
      expect_identical(is.nan(got$draw_rho), is.nan(ex$draw_rho), info = info)
      expect_true(maxabs(got$draw_rho, ex$draw_rho) < 1e-12, info = info)
      expect_true(maxabs(got$draw_variances, ex$draw_variances) < 1e-12 * max(abs(ex$draw_variances), na.rm = TRUE), info = info)
      expect_true(maxabs(plain(got$rho), ex$rho) < 1e-12, info = info)
      expect_identical(is.na(plain(got$rho)), is.na(ex$rho), info = info)
      expect_true(maxabs(unname(attr(got$rho, 'variances')), ex$variances) < 1e-12 * max(abs(ex$variances)), info = info)
      got_ex <- lapply(seq_along(cc$draws), function (d) {
        e <- got$exclusions[, d]; as.integer(e[!is.na(e)]) })
      expect_identical(got_ex, lapply(ex$exclusions, as.integer), info = info)
      expect_identical(got$n_exclusions, lengths(ex$exclusions), info = info)
      expect_identical(got$stop, ex$stop, info = info)
    }
  }

  # The fixtures cover ties (the first pair in packed order wins: features
  # 1, 2 and 3 are identical) and 20 exclusions among 30 features.
  expect_identical(as.integer(ref$ties$expected$exclusions[[1]][1:3]), c(1L, 2L, 8L))
  expect_true(all(lengths(ref$many$expected$exclusions) == 20L))
  expect_true(any(is.nan(ref$small$expected$draw_rho)))


  # -----------------------------------------------------------------------
  # Exclusion threshold 0: rounding decides which pairs whose correlation
  # is 0 are excluded, so features are dropped and systems become singular.
  # Each draw against dense solves with its own exclusions.
  # -----------------------------------------------------------------------

  set.seed(3)
  stops <- character(0)
  for (k in 1:30) {
    D <- 4L + k %% 4L
    n <- 9L
    L <- lapply(1:3, function (d) matrix(rnorm(n * D), n, D))
    got <- suppressWarnings(sparcc_details(
      matrix(1L, n, D), log_fractions = L, exclusion_threshold = 0, cpus = 1L ))
    stops <- c(stops, got$stop)
    for (d in seq_along(L)) {
      info <- paste(k, d)
      t  <- variation(L[[d]])
      e  <- got$exclusions[, d]; e <- e[!is.na(e)]
      pr <- pairs_of(D)
      # A feature is dropped when all its pairs are excluded.
      kept   <- setdiff(seq_len(nrow(pr)), e)
      active <- seq_len(D) %in% pr[kept, ]
      w_got  <- got$draw_variances[d, ]
      expect_identical(is.na(w_got), !active, info = info)
      expect_true(sum(active) >= 3L, info = info)
      # The variances that the last rho came from: for a stop at 3 features,
      # those before the last exclusion.
      if (got$stop[[d]] == 'three') {
        e0 <- e[-length(e)]
        a0 <- seq_len(D) %in% pr[setdiff(seq_len(nrow(pr)), e0), ]
        w  <- dense_w(t, D, e0, a0)
        w[!active] <- NA
      }
      else {
        w <- dense_w(t, D, e, active)
      }
      expect_true(maxabs(w_got, w) < 1e-10 * max(abs(w), na.rm = TRUE), info = info)
      r <- rho_from(w_got, t, D)
      expect_identical(is.na(got$draw_rho[d, ]), is.na(r), info = info)
      expect_true(maxabs(got$draw_rho[d, ], r) < 1e-12, info = info)
    }
  }
  expect_true(all(c('three', 'singular', 'threshold') %in% stops))


  # -----------------------------------------------------------------------
  # Closed form (no exclusions) equals solve(M, t_i).
  # -----------------------------------------------------------------------

  set.seed(4)
  L <- list(matrix(rnorm(15 * 9), 15, 9))
  got <- sparcc_details(matrix(1L, 15, 9), log_fractions = L, max_exclusions = 0L)
  t <- variation(L[[1]])
  expect_true(maxabs(got$draw_variances[1, ], dense_w(t, 9L, integer(0), rep(TRUE, 9))) < 1e-12)
  expect_identical(got$stop, 'max_exclusions')


  # -----------------------------------------------------------------------
  # Invariances of the deterministic steps.
  # -----------------------------------------------------------------------

  cc  <- ref$correlated
  D   <- ncol(cc$counts)
  run <- function (draws, counts = cc$counts) suppressWarnings(sparcc_details(
    counts, log_fractions = draws, max_exclusions = 20L, cpus = 2L ))
  base <- run(cc$draws)

  # Each sample's scale cancels: shifted log fractions give the same.
  shifted <- lapply(cc$draws, function (L) L + rnorm(nrow(L), sd = 3))
  expect_true(maxabs(plain(run(shifted)$rho), plain(base$rho)) < 1e-12)

  # The counts are not used with given draws (only their shape).
  expect_identical(plain(run(cc$draws, cc$counts * 10L)$rho), plain(base$rho))

  # Feature order: a permuted input gives the permuted output.
  perm <- c(7, 3, 12, 1, 5, 9, 2, 11, 4, 8, 10, 6)
  got  <- run(lapply(cc$draws, function (L) L[, perm]), cc$counts[, perm])
  m0   <- as.matrix(base$rho); m1 <- as.matrix(got$rho)
  expect_true(max(abs(m1 - m0[perm, perm])) < 1e-12)


  # -----------------------------------------------------------------------
  # The Dirichlet draws: threads, input types, seeds, blocks.
  # -----------------------------------------------------------------------

  set.seed(5)
  n <- 40; D <- 25
  basis <- matrix(rnorm(n * D), n, D)
  basis[, 2] <- basis[, 1] + rnorm(n, sd = 0.3)
  basis[, 4] <- -basis[, 3] + rnorm(n, sd = 0.3)
  x <- matrix(rpois(n * D, exp(basis + 3)), n, D,
    dimnames = list(paste0('s', 1:n), paste0('f', 1:D)) )
  x[1:5, 6] <- 0L

  r1 <- sparcc(x, n_draws = 20L, cpus = 1L)
  for (cpus in c(2L, 3L, 4L, 8L))
    expect_identical(sparcc(x, n_draws = 20L, cpus = cpus), r1, info = cpus)
  expect_identical(sparcc(Matrix::Matrix(x, sparse = TRUE), n_draws = 20L, cpus = 4L), r1)
  expect_identical(sparcc(slam::as.simple_triplet_matrix(x), n_draws = 20L, cpus = 2L), r1)
  expect_identical(sparcc(t(x), margin = 2L, n_draws = 20L, cpus = 2L), r1)
  expect_identical(sparcc(x * 1.0, n_draws = 20L), r1)
  expect_identical(sparcc(as.data.frame(x), n_draws = 20L), r1)
  expect_identical(sparcc(x, n_draws = 20L, seed = 0L), r1)
  expect_false(identical(plain(sparcc(x, n_draws = 20L, seed = 1L)), plain(r1)))

  # Draw d doesn't depend on n_draws: the first 20 of 30 draws are r1's.
  d20 <- sparcc_details(x, n_draws = 20L)
  d30 <- sparcc_details(x, n_draws = 30L)
  expect_identical(d30$draw_rho[1:20, ], d20$draw_rho)
  expect_identical(plain(d20$rho), plain(r1))

  # The medians are over the draws.
  expect_identical(plain(r1), unname(apply(d20$draw_rho, 2, stats::median, na.rm = TRUE)))
  expect_identical(unname(attr(r1, 'variances')), unname(apply(d20$draw_variances, 2, stats::median, na.rm = TRUE)))

  # Pairs stored in blocks (draws made again for each) give the same.
  for (bp in c(1, 7, 100, 299))
    expect_identical(
      sparcc_details(x, n_draws = 20L, block_pairs = bp, details = FALSE, cpus = 3L), r1, info = bp )
  rp <- sparcc(x, n_draws = 6L, permutations = 9L, cpus = 2L)
  expect_identical(
    sparcc_details(x, n_draws = 6L, permutations = 9L, block_pairs = 13, details = FALSE, cpus = 4L), rp )

  # The correlated pairs are found; the others are near 0.
  m <- as.matrix(r1)
  expect_true(m['f1', 'f2'] > 0.8)
  expect_true(m['f3', 'f4'] < -0.8)
  others <- m[5:D, 5:D][upper.tri(m[5:D, 5:D])]
  expect_true(mean(abs(others) < 0.3) > 0.9)


  # -----------------------------------------------------------------------
  # FastSpar and SpiecEasi, with the same settings: agreement to the noise
  # of the draws, about 0.002 for counts in the thousands (FastSpar prints 4
  # decimals). FastSpar makes at most D - 3 exclusions, which matters only
  # with fewer than max_exclusions + 3 features.
  # -----------------------------------------------------------------------

  refs <- ref$refs
  tb   <- refs$tables
  fs   <- refs$fastspar
  near <- function (d, m, tol) {
    a <- as.matrix(d); diag(a) <- 1
    max(abs(a - m[colnames(a), colnames(a)])) < tol
  }
  vnear <- function (d, v, tol) max(abs(attr(d, 'variances') - v) / v) < tol

  r <- sparcc(tb$tight, n_draws = 100L, max_exclusions = 10L, cpus = 4L)
  expect_true(near(r, fs$tight_i100_x10$cor, 0.02))
  expect_true(near(r, refs$spieceasi$tight_i100_x10, 0.02))
  expect_true(vnear(r, fs$tight_i100_x10$variances, 0.02))
  r <- sparcc(tb$tight, n_draws = 100L, max_exclusions = 17L, cpus = 4L)
  expect_true(near(r, fs$tight_i100_x17$cor, 0.02))
  expect_true(near(r, fs$tight_i100_x20$cor, 0.02))    # stops at D - 3 = 17
  expect_identical(fs$tight_i100_x20$cor, fs$tight_i100_x17$cor)
  expect_false(near(sparcc(tb$tight, n_draws = 100L, max_exclusions = 20L, cpus = 4L), fs$tight_i100_x20$cor, 0.02))

  r <- sparcc(tb$d8, n_draws = 200L, max_exclusions = 5L, cpus = 4L)
  expect_true(near(r, fs$d8_i200_x5$cor, 0.02))
  expect_true(vnear(r, fs$d8_i200_x5$variances, 0.02))
  r <- sparcc(tb$d8, n_draws = 200L, max_exclusions = 10L, cpus = 4L)
  expect_true(near(r, refs$spieceasi$d8_i200_x10, 0.02))    # no D - 3 limit
  expect_false(near(r, fs$d8_i200_x10$cor, 0.03))

  # Counts around 20: the medians of 1,000 draws differ by about 0.003 from
  # run to run (per-draw spread 0.06), and by as much between two seeds.
  r <- sparcc(tb$moderate, n_draws = 1000L, max_exclusions = 20L, cpus = 4L)
  a <- as.matrix(r); diag(a) <- 1
  d <- abs(a - fs$moderate_i1000_x20$cor)
  expect_true(mean(d[upper.tri(d)]) < 0.005 && max(d) < 0.05)
  expect_true(vnear(r, fs$moderate_i1000_x20$variances, 0.05))

  # p-values, from 99 null datasets each: the same pairs stand out.
  r  <- sparcc(tb$moderate, n_draws = 20L, max_exclusions = 10L, permutations = 99L, cpus = 4L)
  fp <- fs$moderate_i20_x10$p_value
  fp <- fp[lower.tri(fp)]
  p  <- attr(r, 'p_value')
  expect_true(stats::cor(p, fp, method = 'spearman') > 0.8)
  expect_identical(p[[1]], 0.01)
  expect_identical(fp[[1]], 0.01)
  expect_true(mean((p <= 0.05) == (fp <= 0.05)) > 0.9)


  # -----------------------------------------------------------------------
  # A tight posterior (large counts): the median over Dirichlet draws is
  # close to SparCC of the observed log fractions.
  # -----------------------------------------------------------------------

  set.seed(6)
  n <- 30; D <- 15
  basis <- matrix(rnorm(n * D, sd = 0.7), n, D)
  basis[, 2] <- basis[, 1] + rnorm(n, sd = 0.4)
  big <- matrix(rpois(n * D, exp(basis + 9)), n, D)
  r   <- sparcc(big, n_draws = 40L, cpus = 2L)
  obs <- sparcc_details(big, log_fractions = list(log(big)), cpus = 2L)
  expect_true(maxabs(plain(r), plain(obs$rho)) < 0.02)


  # -----------------------------------------------------------------------
  # Shuffled data: no strong correlations (the paper's null check).
  # -----------------------------------------------------------------------

  set.seed(7)
  n <- 200; D <- 30
  basis <- matrix(rnorm(n * D), n, D)
  basis[, 2] <- basis[, 1] + rnorm(n, sd = 0.3)
  y  <- matrix(rpois(n * D, exp(basis + 4)), n, D)
  ys <- apply(y, 2, sample)
  rs <- sparcc(ys, n_draws = 10L, cpus = 2L)
  expect_true(mean(abs(rs) > 0.3) < 0.005)
  expect_true(abs(as.matrix(sparcc(y, n_draws = 10L, cpus = 2L))[1, 2]) > 0.8)
  y <- y[1:50, ]


  # -----------------------------------------------------------------------
  # p-values.
  # -----------------------------------------------------------------------

  B  <- 19L
  dp <- sparcc_details(y[, 1:10], n_draws = 8L, permutations = B, cpus = 2L)
  p  <- attr(dp$rho, 'p_value')
  expect_identical(length(p), length(dp$rho))
  expect_true(all(p > 0 & p <= 1))
  expect_true(all(abs(p * (B + 1) - round(p * (B + 1))) < 1e-9))
  expect_identical(p, perm_p(plain(dp$rho), t(dp$null_rho), alternative = 'two.sided'))
  expect_identical(dim(dp$null_rho), c(45L, B))
  expect_true(p[[1]] == 1 / (B + 1))   # f1 ~ f2
  for (cpus in c(1L, 4L))
    expect_identical(
      sparcc(y[, 1:10], n_draws = 8L, permutations = B, cpus = cpus), dp$rho, info = cpus )
  # Each null dataset resamples each feature's counts: its own tables and
  # draws.
  expect_false(isTRUE(all.equal(dp$null_rho[, 1], dp$null_rho[, 2])))

  # Aligned with the triangle after clustering.
  rc <- sparcc(y[, 1:10], n_draws = 8L, permutations = B, clust = 'average', cpus = 2L)
  o  <- attr(rc, 'order')
  expect_true(!identical(o, 1:10))
  pm <- function (d, v) { m <- as.matrix(d); m[] <- NA; m[lower.tri(m)] <- v; m[upper.tri(m)] <- t(m)[upper.tri(m)]; m }
  expect_identical(pm(rc, attr(rc, 'p_value')), pm(dp$rho, p)[o, o])
  expect_identical(as.matrix(rc), as.matrix(dp$rho)[o, o])
  expect_identical(unname(attr(rc, 'variances')), unname(attr(dp$rho, 'variances'))[o])
  expect_identical(names(attr(rc, 'variances')), attr(rc, 'Labels'))

  # Clustered as the dissimilarity 1 - r.
  d1 <- 1 - dp$rho
  attributes(d1) <- list(Size = 10L, Labels = attr(dp$rho, 'Labels'), class = 'dist')
  expect_identical(o, stats::hclust(d1, 'average')$order)


  # -----------------------------------------------------------------------
  # The result: class, attributes, title, registry, printing and plots.
  # -----------------------------------------------------------------------

  expect_identical(class(r1), c('ecodive_dist', 'dist'))
  expect_identical(attr(r1, 'method'), 'sparcc')
  expect_identical(attr(r1, 'n_samples'), 40L)
  expect_identical(attr(r1, 'Labels'), colnames(x))
  expect_identical(names(attr(r1, 'variances')), colnames(x))
  expect_true(all(attr(r1, 'variances') > 0))
  expect_null(attr(r1, 'p_value'))
  expect_identical(attr(r1, 'display')$title, 'SparCC Correlation')
  expect_true(is_similarity(r1))
  expect_identical(dist_plot_name(r1)$unit, 'features')
  expect_false(is_similarity(vlr(x + 1L)))
  expect_true(is_similarity(rho_p(x + 1L)))
  expect_identical(attr(rho_p(x + 1L), 'display')$title, 'Proportionality (rho_p)')
  expect_identical(FEATURE_ASSOCS$id, c('vlr', 'rho_p', 'phi_s', 'sparcc'))
  expect_true(any(grepl('SparCC Correlation', utils::capture.output(print(r1)))))
  expect_true(any(grepl('25 features', utils::capture.output(print(r1)))))

  # Unnamed features: no labels.
  u <- sparcc(unname(x), n_draws = 5L)
  expect_null(attr(u, 'Labels'))
  expect_null(names(attr(u, 'variances')))

  # plot(): the heatmap, titled from the result.
  grDevices::pdf(NULL)
  expect_identical(plot(r1), r1)
  expect_identical(plot(rc), rc)
  grDevices::dev.off()

  # No registry has an alias starting with 'spar'.
  for (reg in c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'ZERO_REPLACEMENTS'))
    if (exists(reg)) expect_identical(grep('^spar', names(exact_names(reg)), value = TRUE), character(0), info = reg)

  # dist_reorder() reorders the attributes parallel to the triangle.
  d <- stats::dist(matrix(1:12, 4))
  attr(d, 'p_value') <- seq_along(d) / 10
  attr(d, 'other')   <- 1:6
  o2 <- c(3L, 1L, 4L, 2L)
  dr <- dist_reorder(d, o2)
  expect_identical(pm(dr, attr(dr, 'p_value')), pm(d, attr(d, 'p_value'))[o2, o2])
  expect_identical(attr(dr, 'other'), 1:6)


  # -----------------------------------------------------------------------
  # Warnings and errors.
  # -----------------------------------------------------------------------

  # Few samples: some variances not positive.
  expect_warning(sparcc(ecodive::ex_counts, n_draws = 20L), 'not positive')
  expect_warning(sparcc(x + 0.5, n_draws = 5L), 'whole numbers')
  expect_silent(sparcc(x, n_draws = 5L))

  expect_error(sparcc(x[, 1:3]), 'at least 4 features')
  expect_error(sparcc(x[1:2, ]), 'at least 3 samples')
  neg <- x; neg[1, 1] <- -1L
  expect_error(sparcc(neg), 'negative')
  na <- x * 1.0; na[1, 1] <- NA
  expect_error(sparcc(na), 'NA')
  expect_error(sparcc(matrix(letters[1:20], 4)), 'numeric')
  expect_error(sparcc(x, n_draws = 0), '`n_draws` must be a single positive integer')
  expect_error(sparcc(x, n_draws = 2.5), '`n_draws`')
  expect_error(sparcc(x, max_exclusions = -1), '`max_exclusions` must be a single non-negative integer')
  expect_error(sparcc(x, permutations = NA), '`permutations`')
  expect_error(sparcc(x, permutations = matrix(1:4, 2)), '`permutations`')
  expect_error(sparcc(x, exclusion_threshold = -0.1), '`exclusion_threshold`')
  expect_error(sparcc(x, exclusion_threshold = c(0.1, 0.2)), '`exclusion_threshold`')
  expect_error(sparcc(x, seed = 1.5), '`seed`')
  expect_error(sparcc(x, cpus = 0), '`cpus`')
  expect_error(sparcc(x, clust = 'foo'), '`clust`')
  expect_error(sparcc(x, margin = 3), 'margin')

  # max_exclusions = 0 and Inf thresholds run.
  expect_identical(sparcc_details(x, n_draws = 3L, max_exclusions = 0L)$n_exclusions, rep(0L, 3))
  expect_identical(sparcc_details(x, n_draws = 3L, exclusion_threshold = Inf)$n_exclusions, rep(0L, 3))
  expect_identical(
    suppressWarnings(sparcc(x[, 1:5], n_draws = 3L, max_exclusions = .Machine$integer.max)),
    suppressWarnings(sparcc(x[, 1:5], n_draws = 3L, max_exclusions = 10L)) )

#})
