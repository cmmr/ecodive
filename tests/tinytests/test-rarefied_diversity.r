#test_that("rarefied diversity", {

  # beta_div_rarefied() and alpha_div_rarefied(): every average is rebuilt
  # from ecodive's own rarefy() + beta_div() / alpha_div(), and the mean is
  # compared with vegan::avgdist() (data-raw/expected/rarefied_diversity.r).

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))
  exit_if_not(requireNamespace("slam",   quietly = TRUE))

  rel_diff <- function (x, y) max(abs(x - y)) / max(abs(y))

  # Attributes added by the rarefied functions, and the title.
  bare <- function (z) {
    for (a in c('depth', 'dropped', 'sd', 'medoid', 'display', 'method'))
      attr(z, a) <- NULL
    z
  }

  set.seed(11)
  x <- matrix(
    rpois(13 * 25, rep(rexp(25, 1 / 6), each = 13)), 13, 25,
    dimnames = list(sprintf('S%02d', 1:13), sprintf('F%02d', 1:25)) )
  x[13, ] <- 0L
  x[13, 1:3] <- 5L                 # 15 observations: dropped at depth 60
  x[1:12, 1] <- x[1:12, 1] + 40L   # every other sample has at least 60
  dd <- 60L
  ok <- rowSums(x) >= dd
  stopifnot(sum(!ok) == 1L)
  xk <- x[ok, ]

  k <- 7L
  stack_of <- function (metric, ...) {
    rs <- rarefy(x, dd, seed = 5L, times = k, warn = FALSE)
    lapply(rs, beta_div, metric, ...)
  }
  as_cols <- function (ds) vapply(ds, as.vector, numeric(length(ds[[1]])))


  # ----------------------------------------------------------------------
  # times = 1: exactly beta_div(rarefy()) and alpha_div(rarefy()).
  # ----------------------------------------------------------------------

  expect_warning(
    z <- beta_div_rarefied(x, 'bray', depth = dd, times = 1, seed = 3),
    '1 of the 13 samples has fewer than 60 observations and was dropped' )
  ref <- beta_div(rarefy(x, dd, seed = 3, warn = FALSE), 'bray')
  expect_identical(bare(z), bare(ref))
  expect_identical(attr(z, 'Labels'), rownames(xk))
  expect_identical(attr(z, 'method'), 'bray')
  expect_identical(attr(z, 'depth'), dd)
  expect_identical(attr(z, 'dropped'), 'S13')
  expect_identical(attr(z, 'display')$title, 'Bray-Curtis Dissimilarity (rarefied 1\u00d7)')

  a   <- suppressWarnings(alpha_div_rarefied(x, 'shannon', depth = dd, times = 1, seed = 3))
  ref <- alpha_div(rarefy(x, dd, seed = 3, warn = FALSE), 'shannon')
  expect_identical(as.vector(a), as.vector(ref))
  expect_identical(names(a), names(ref))
  expect_identical(attr(a, 'depth'), dd)
  expect_identical(attr(a, 'dropped'), 'S13')
  expect_null(attr(a, 'display'))


  # ----------------------------------------------------------------------
  # The mean: identical to Reduce('+', ...) / k; sd to apply(, 1, sd).
  # ----------------------------------------------------------------------

  for (m in c('bray', 'jaccard', 'aitchison')) {
    extra <- if (m == 'aitchison') list(pseudocount = 1) else list()
    ds <- do.call(stack_of, c(list(m), extra))
    z  <- suppressWarnings(do.call(
      beta_div_rarefied,
      c(list(x, m, depth = dd, times = k, seed = 5L, sd = TRUE), extra) ))
    expect_identical(as.vector(z), as.vector(Reduce('+', ds) / k), info = m)
    expect_identical(attr(z, 'Labels'), attr(ds[[1]], 'Labels'), info = m)
    expect_identical(attr(z, 'method'), m, info = m)
    M <- as_cols(ds)
    s <- attr(z, 'sd')
    expect_inherits(s, 'ecodive_dist')
    expect_true(rel_diff(as.vector(s), apply(M, 1, sd)) < 1e-12, info = m)
    expect_identical(attr(s, 'Labels'), attr(z, 'Labels'))
    expect_null(attr(s, 'method'))
    expect_identical(
      attr(s, 'display')$title,
      paste('SD of', attr(z, 'display')$title) )
  }

  # Weighted UniFrac, with a tree (helper.r's 4 x 5 table: depth 13).
  ds <- lapply(rarefy(counts, 13, seed = 2L, times = k), weighted_unifrac, tree = tree)
  z  <- beta_div_rarefied(counts, 'w_unifrac', depth = 13, times = k, seed = 2L, tree = tree, sd = TRUE)
  expect_identical(as.vector(z), as.vector(Reduce('+', ds) / k))
  expect_true(rel_diff(as.vector(attr(z, 'sd')), apply(as_cols(ds), 1, sd)) < 1e-12)
  expect_identical(attr(z, 'method'), 'weighted_unifrac')

  # Alpha: the sequential sum over k, and apply(, 1, sd).
  as_ <- lapply(rarefy(x, dd, seed = 5L, times = k, warn = FALSE), alpha_div, 'shannon')
  A   <- vapply(as_, as.vector, numeric(sum(ok)))
  a   <- suppressWarnings(alpha_div_rarefied(x, 'shannon', depth = dd, times = k, seed = 5L, sd = TRUE))
  expect_identical(as.vector(a), as.vector(Reduce('+', lapply(as_, as.vector)) / k))
  expect_identical(names(a), rownames(xk))
  expect_identical(as.vector(attr(a, 'sd')), apply(A, 1, sd))
  expect_identical(names(attr(a, 'sd')), rownames(xk))

  # times = 1: sd is NA.
  z <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = 1, sd = TRUE))
  expect_true(all(is.na(attr(z, 'sd'))))


  # ----------------------------------------------------------------------
  # Median and medoid.
  # ----------------------------------------------------------------------

  # An R version of q2-boots' _medoid(): Euclidean distances between the
  # iterations' condensed matrices, the first argmin of their sums.
  q2_medoid <- function (M) {
    n <- ncol(M)
    sums <- numeric(n)
    for (i in seq_len(n)) for (j in seq_len(i - 1L)) {
      dij <- sqrt(sum((M[, i] - M[, j])^2))
      sums[i] <- sums[i] + dij
      sums[j] <- sums[j] + dij
    }
    which.min(sums)
  }

  fit <- function (average, budget = 256 * 2^20, metric = 'bray', ...) {
    args <- list(
      counts = x, metric = metric, depth = dd, times = k, average = average,
      sd = TRUE, margin = 1L, norm = 'none', pseudocount = NULL, power = 1.5,
      alpha = 0.5, tree = NULL, pairs = NULL, clust = NULL, seed = 5L,
      cpus = 1L, budget = budget )
    extra <- list(...)
    args[names(extra)] <- extra
    suppressWarnings(do.call(rarefied_beta_fit, args))
  }

  for (m in c('bray', 'jaccard')) {
    M <- as_cols(stack_of(m))

    zm <- fit('median', metric = m)
    expect_true(rel_diff(as.vector(zm), apply(M, 1, median)) <= 1e-15, info = m)
    expect_true(rel_diff(as.vector(attr(zm, 'sd')), apply(M, 1, sd)) < 1e-12, info = m)
    expect_null(attr(zm, 'medoid'))

    zd <- fit('medoid', metric = m)
    expect_identical(attr(zd, 'medoid'), q2_medoid(M), info = m)
    expect_identical(as.vector(zd), M[, attr(zd, 'medoid')], info = m)
    expect_true(rel_diff(as.vector(attr(zd, 'sd')), apply(M, 1, sd)) < 1e-12, info = m)

    # Blocks of 3 pairs, and of 1: identical to one block.
    for (b in c(8 * k * 3, 1)) {
      expect_identical(fit('median', budget = b, metric = m), zm, info = paste(m, b))
      expect_identical(fit('medoid', budget = b, metric = m), zd, info = paste(m, b))
    }
  }

  # The exported function gives the same.
  expect_identical(
    suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, average = 'MEDI', sd = TRUE)),
    fit('median') )
  expect_identical(
    suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, average = 'medo', sd = TRUE)),
    fit('medoid') )

  # Alpha median.
  a <- suppressWarnings(alpha_div_rarefied(x, 'shannon', depth = dd, times = k, seed = 5L, average = 'median'))
  expect_identical(as.vector(a), apply(A, 1, median))

  # Medians of an even count and with NA, as median().
  R <- rbind(c(3, 1, 2, 4), c(1, 1, 1, 1), c(5, NA, 1, 2), c(2.5, -1, 1e300, 7))
  expect_identical(.Call(C_row_medians, R, 1L), apply(R, 1, median))
  R <- matrix(c(4, 2, 9, 1, 7), 1)
  expect_identical(.Call(C_row_medians, R, 1L), 4)
  R <- matrix(runif(250 * 9), 250)
  for (cpus in c(1L, 2L, 4L)) {
    expect_identical(.Call(C_row_medians, R, cpus), apply(R, 1, median))
    expect_identical(.Call(C_row_medians, R[, 1:8], cpus), apply(R[, 1:8], 1, median))
  }


  # ----------------------------------------------------------------------
  # Statistics: observed richness against its exact expectation under
  # rarefaction (Hurlbert 1971), within 4 Monte Carlo SEs.
  # ----------------------------------------------------------------------

  n_it <- 2000L
  a    <- suppressWarnings(alpha_div_rarefied(x, 'observed', depth = dd, times = n_it, sd = TRUE))
  ex   <- apply(xk, 1, function (v) {
    N <- sum(v)
    sum(1 - exp(lchoose(N - v[v > 0], dd) - lchoose(N, dd)))
  })
  se <- attr(a, 'sd') / sqrt(n_it)
  expect_true(all(abs(a - ex) <= 4 * se))


  # ----------------------------------------------------------------------
  # vegan::avgdist(): the mean over 5000 rarefactions agrees within 4
  # combined SEs in every cell.
  # ----------------------------------------------------------------------

  fx <- readRDS(file.path('expected', 'rarefied_diversity.rds'))
  z  <- beta_div_rarefied(fx$counts, 'bray', depth = fx$depth, times = fx$iterations, sd = TRUE)
  se <- sqrt((as.vector(attr(z, 'sd'))^2 + as.vector(fx$sd)^2) / fx$iterations)
  expect_identical(attr(z, 'Labels'), attr(fx$mean, 'Labels'))
  expect_true(all(abs(as.vector(z) - as.vector(fx$mean)) <= 4 * se + 1e-12))


  # ----------------------------------------------------------------------
  # Titles carry through ordinations and tests.
  # ----------------------------------------------------------------------

  z <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = 3))
  expect_identical(attr(z, 'display')$title, 'Bray-Curtis Dissimilarity (rarefied 3\u00d7)')
  expect_identical(
    attr(pcoa(z), 'display')$title,
    'Principal Coordinates Analysis of Bray-Curtis Dissimilarity (rarefied 3\u00d7)' )
  grp <- rep(c('a', 'b'), length.out = sum(ok))
  expect_identical(
    attr(permanova(z, grp, permutations = 9), 'display')$title,
    'Permutational multivariate analysis of variance (PERMANOVA) of Bray-Curtis Dissimilarity (rarefied 3\u00d7)' )
  pdf(NULL)
  expect_silent(plot(z))
  dev.off()
  expect_true(any(grepl('rarefied 3', capture.output(print(z)), fixed = TRUE)))


  # ----------------------------------------------------------------------
  # pairs, clust, norm, margin = 2, table classes, embedded trees.
  # ----------------------------------------------------------------------

  full <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L))
  prs  <- c(2L, 5L, 9L, 40L)
  for (avg in c('mean', 'median', 'medoid')) {
    zp <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, pairs = prs, average = avg))
    expect_identical(which(!is.na(zp)), prs, info = avg)
    if (avg != 'medoid')
      expect_identical(
        as.vector(zp)[prs],
        as.vector(suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, average = avg)))[prs],
        info = avg )
  }
  expect_identical(as.vector(full)[prs], as.vector(suppressWarnings(
    beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, pairs = prs)))[prs])
  expect_error(
    suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = 2, pairs = 67)),
    'between 1 and 66' )

  # clust: reorders the result and its sd alike.
  zc <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, sd = TRUE, clust = 'average'))
  zs <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, sd = TRUE))
  o  <- attr(zc, 'order')
  expect_identical(o, stats::hclust(stats::as.dist(as.matrix(zs)), 'average')$order)
  expect_identical(as.matrix(zc), as.matrix(zs)[o, o])
  expect_identical(as.matrix(attr(zc, 'sd')), as.matrix(attr(zs, 'sd'))[o, o])
  expect_identical(attr(attr(zc, 'sd'), 'order'), o)
  expect_identical(attr(zc, 'display')$clust, 'average')

  # norm, applied after each rarefaction.
  ds <- lapply(rarefy(x, dd, seed = 5L, times = k, warn = FALSE), beta_div, 'manhattan', norm = 'percent')
  zn <- suppressWarnings(beta_div_rarefied(x, 'manhattan', depth = dd, times = k, seed = 5L, norm = 'percent'))
  expect_identical(as.vector(zn), as.vector(Reduce('+', ds) / k))

  # margin = 2, and every table class: the same result.
  classes <- list(
    dgC  = function (m) as(m, 'CsparseMatrix'),
    dgT  = function (m) as(m, 'TsparseMatrix'),
    slam = slam::as.simple_triplet_matrix )
  for (avg in c('mean', 'median', 'medoid')) {
    zb <- suppressWarnings(beta_div_rarefied(x, 'bray', depth = dd, times = k, seed = 5L, average = avg, sd = TRUE))
    ab <- suppressWarnings(alpha_div_rarefied(x, 'shannon', depth = dd, times = k, seed = 5L, average = if (avg == 'medoid') 'mean' else avg))
    expect_identical(
      suppressWarnings(beta_div_rarefied(t(x), 'bray', depth = dd, times = k, seed = 5L, average = avg, sd = TRUE, margin = 2)),
      zb, info = avg )
    for (cl in names(classes)) {
      for (margin in 1:2) {
        m <- classes[[cl]](if (margin == 1L) x else t(x))
        expect_identical(
          suppressWarnings(beta_div_rarefied(m, 'bray', depth = dd, times = k, seed = 5L, average = avg, sd = TRUE, margin = margin)),
          zb, info = paste(cl, margin, avg) )
        expect_identical(
          suppressWarnings(alpha_div_rarefied(m, 'shannon', depth = dd, times = k, seed = 5L, average = if (avg == 'medoid') 'mean' else avg, margin = margin)),
          ab, info = paste(cl, margin, avg) )
      }
    }
  }

  # An rbiom-like object with its tree embedded (samples in columns).
  rb <- structure(list(counts = as(t(counts), 'CsparseMatrix'), tree = tree), class = 'rbiom')
  expect_identical(
    beta_div_rarefied(rb, 'unweighted_unifrac', depth = 13, times = k, seed = 1L),
    beta_div_rarefied(counts, 'unweighted_unifrac', depth = 13, times = k, seed = 1L, tree = tree) )
  expect_identical(
    alpha_div_rarefied(rb, 'faith', depth = 13, times = k, seed = 1L),
    alpha_div_rarefied(counts, 'faith', depth = 13, times = k, seed = 1L, tree = tree) )

  # Unnamed samples: dropped are positions.
  xu <- unname(x)
  z  <- suppressWarnings(beta_div_rarefied(xu, 'bray', depth = dd, times = 2))
  expect_identical(attr(z, 'dropped'), 13L)
  # None dropped: an empty vector.
  z  <- beta_div_rarefied(xk, 'bray', depth = dd, times = 2)
  expect_identical(attr(z, 'dropped'), character(0))

  # depth = NULL: rarefy()'s automatic depth.
  z <- beta_div_rarefied(counts, 'bray', times = 2)
  expect_identical(attr(z, 'depth'), 13L)
  expect_identical(bare(z), bare(beta_div_rarefied(counts, 'bray', depth = 13, times = 2)))

  # A warning from the metric is given once.
  w <- character(0)
  withCallingHandlers(
    beta_div_rarefied(xk, 'aitchison', depth = dd, times = 5),
    warning = function (cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart('muffleWarning') })
  expect_identical(length(w), 1L)
  expect_true(grepl('pseudocount', w))


  # ----------------------------------------------------------------------
  # The same for any cpus.
  # ----------------------------------------------------------------------

  big <- matrix(
    rpois(150 * 40, rep(rexp(40, 1 / 5), each = 150)) + 1L, 150, 40,
    dimnames = list(sprintf('S%03d', 1:150), sprintf('F%02d', 1:40)) )
  for (avg in c('mean', 'median', 'medoid')) {
    z1 <- beta_div_rarefied(big, 'bray', depth = 100, times = 12, average = avg, sd = TRUE, cpus = 1L)
    a1 <- alpha_div_rarefied(big, 'shannon', depth = 100, times = 12, average = if (avg == 'medoid') 'median' else avg, sd = TRUE, cpus = 1L)
    for (cpus in c(2L, 4L, 8L)) {
      expect_identical(beta_div_rarefied(big, 'bray', depth = 100, times = 12, average = avg, sd = TRUE, cpus = cpus), z1, info = paste(avg, cpus))
      expect_identical(alpha_div_rarefied(big, 'shannon', depth = 100, times = 12, average = if (avg == 'medoid') 'median' else avg, sd = TRUE, cpus = cpus), a1, info = paste(avg, cpus))
    }
  }
  M <- matrix(runif(200 * 30), 200)
  G <- matrix(0, 30, 30)
  for (cpus in c(2L, 4L))
    expect_identical(.Call(C_medoid_accum, M, G, cpus), .Call(C_medoid_accum, M, G, 1L))
  expect_equal(.Call(C_medoid_accum, M, G, 1L), unname(as.matrix(dist(t(M)))^2), tolerance = 1e-12)


  # ----------------------------------------------------------------------
  # Errors.
  # ----------------------------------------------------------------------

  expect_error(beta_div_rarefied(x + 0.5, 'bray'), 'whole numbers')
  expect_error(alpha_div_rarefied(x + 0.5, 'shannon'), 'whole numbers')
  expect_error(beta_div_rarefied(x, 'bray', depth = 0), 'single whole number of at least 1')
  expect_error(beta_div_rarefied(x, 'bray', depth = -5), 'single whole number of at least 1')
  expect_error(beta_div_rarefied(x, 'bray', depth = NA), 'single whole number of at least 1')
  expect_error(beta_div_rarefied(x, 'bray', depth = 2.5), 'single whole number of at least 1')
  expect_error(beta_div_rarefied(x, 'bray', depth = c(50, 60)), 'use rarefy')
  expect_error(beta_div_rarefied(x, 'bray', times = 0), '`times` must be a whole number')
  expect_error(beta_div_rarefied(x, 'bray', times = NA), '`times` must be a whole number')
  expect_error(beta_div_rarefied(x, 'bray', times = 1.5), '`times` must be a whole number')
  expect_error(beta_div_rarefied(x, 'bray', average = 'mode'), 'must be one of')
  expect_error(beta_div_rarefied(x, 'bray', average = 'me'), 'matches several')
  expect_error(alpha_div_rarefied(x, 'shannon', average = 'medoid'), 'only for beta diversity')
  expect_error(beta_div_rarefied(x, 'bray', sd = NA), '`sd` must be TRUE or FALSE')
  expect_error(beta_div_rarefied(x, 'bray', depth = 1000), 'Only 0 of the 13 samples have at least 1000')
  expect_error(
    suppressWarnings(beta_div_rarefied(x, 'bray', depth = max(rowSums(x)))),
    'at least two are needed' )
  expect_silent(alpha_div_rarefied(x[1, , drop = FALSE], 'observed', depth = dd, times = 2))
  expect_error(beta_div_rarefied(x, 'nonesuch'), 'does not match any option')
  expect_error(alpha_div_rarefied(x, 'bray'), 'does not match any option')
  expect_error(beta_div_rarefied(-x, 'bray'), 'negative')
  expect_error(beta_div_rarefied(x, 'bray', clust = 'nope'), '`clust` must be')
  expect_error(beta_div_rarefied(x, 'bray', seed = 1.5), '`seed` must be')

  # n_cpus() is held only while the functions run, errors or not.
  expect_null(CACHE$n_cpus_held)
  expect_identical(n_cpus(), parallel_n <- n_cpus())
  CACHE$n_cpus_held <- 1L
  expect_identical(n_cpus(), 1L)
  h <- hold_n_cpus()
  expect_identical(n_cpus(), 1L)
  release_n_cpus(h)
  expect_identical(CACHE$n_cpus_held, 1L)
  CACHE$n_cpus_held <- NULL
  expect_identical(n_cpus(), parallel_n)

#})
