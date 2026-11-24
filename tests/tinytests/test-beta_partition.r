#test_that("partition of beta diversity", {

  # Fixtures: data-raw/expected/beta_partition.r (betapart 1.6.1, output
  # only; abdiv 0.2.0).

  fx   <- readRDS(file.path('expected', 'beta_partition.rds'))
  tol  <- 1e-12
  fams <- c('sorensen', 'jaccard', 'bray', 'ruzicka')

  # Equal within `tol`, with NaN in the same places.
  same <- function (a, b) {
    a <- as.vector(a); b <- as.vector(b)
    identical(is.na(a), is.na(b)) &&
      isTRUE(all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tol))
  }



  # Pairwise: betapart's beta.pair() and beta.pair.abund() ====

  for (tab in names(fx$tables)) {
    x <- fx$tables[[tab]]
    for (fam in fams) {
      bp  <- beta_partition(x, fam, cpus = 1L)
      ref <- fx$betapart_pair[[tab]][[fam]]
      for (k in 1:3)
        expect_true(same(bp[[k]], ref[, k]), info = paste(tab, fam, k))
    }
  }



  # Pairwise: abdiv's component functions, run as the metrics ====

  for (tab in names(fx$tables)) {
    x <- fx$tables[[tab]]
    for (m in colnames(fx$abdiv[[tab]])) {
      d <- beta_div(x, m, cpus = 1L)
      expect_identical(attr(d, 'method'), m)
      expect_true(same(d, fx$abdiv[[tab]][, m]), info = paste(tab, m))
    }
  }



  # Multiple-site: betapart's beta.multi() and beta.multi.abund() ====

  for (tab in names(fx$tables)) {
    x <- fx$tables[[tab]]
    for (fam in fams) {
      v <- beta_partition_multi(x, fam, cpus = 1L)$value
      expect_true(same(v, fx$betapart_multi[[tab]][[fam]]), info = paste(tab, fam))
    }
  }



  # Draws: each is the value of some subset of n_sites samples ====

  sub_x <- fx$tables$sparse[1:9, ]
  for (fam in fams) {
    bm  <- beta_partition_multi(sub_x, fam, n_sites = 4, n_draws = 150, seed = 3, cpus = 1L)
    ref <- fx$betapart_subsets[[fam]]
    hit <- apply(bm$draws, 1L, function (r)
      any(apply(ref, 1L, function (s) all(abs(r - s) <= tol))) )
    expect_true(all(hit), info = fam)
    expect_identical(dim(bm$draws), c(150L, 3L))
    expect_equal(bm$mean, colMeans(bm$draws))
    expect_equal(bm$sd,   apply(bm$draws, 2L, sd))
    # Many different subsets are drawn.
    expect_true(nrow(unique(round(bm$draws, 10))) > 40, info = fam)
  }

  # And their means agree with beta.sample()'s, within 5 standard errors.
  x <- fx$tables$sparse
  for (fam in c('sorensen', 'bray')) {
    bm <- beta_partition_multi(x, fam, n_sites = 10, n_draws = 2000, seed = 7)
    se <- bm$sd * sqrt(2 / 2000)
    expect_true(all(abs(bm$mean - fx$sample_means[[fam]]) < 5 * se), info = fam)
  }

  # All the samples in every draw: the value of all the samples.
  bm <- beta_partition_multi(x, 'jaccard', n_sites = nrow(x), n_draws = 3)
  for (r in 1:3) expect_identical(unname(bm$draws[r, ]), unname(bm$value))

  # Two samples: the pairwise values.
  for (fam in fams) {
    v  <- beta_partition_multi(x[c(3, 9), ], fam)$value
    bp <- beta_partition(x[c(3, 9), ], fam)
    expect_equal(unname(v), unname(vapply(bp, as.vector, 0)), tolerance = tol, info = fam)
  }



  # The total is the base metric; nestedness is total - turnover ====

  base <- list(sorensen = sorensen, jaccard = jaccard, bray = bray, ruzicka = soergel)
  for (fam in fams) {
    bp <- beta_partition(x, fam)
    expect_identical(bp$total, base[[fam]](x), info = fam)
    expect_identical(as.vector(bp[[2]]), as.vector(bp$total) - as.vector(bp[[1]]), info = fam)
    expect_identical(bp[[1]], beta_div(x, attr(bp[[1]], 'method')), info = fam)
    expect_identical(bp[[2]], beta_div(x, attr(bp[[2]], 'method')), info = fam)
  }
  expect_identical(names(beta_partition(x, 'jaccard')), c('turnover', 'nestedness', 'total'))
  expect_identical(names(beta_partition(x, 'bray')),    c('balanced', 'gradient', 'total'))
  expect_identical(attr(beta_partition(x, 'Soergel'), 'index'), 'ruzicka')
  expect_identical(
    names(beta_partition_multi(x, 'ruzicka')$value), c('balanced', 'gradient', 'total') )



  # Edge cases ====

  e  <- fx$tables$edge
  at <- function (d, i, j) as.matrix(d)[i, j]
  st <- sorensen_turnover(e)
  sn <- sorensen_nestedness(e)
  expect_identical(at(st, 3, 4), 0)                     # nested: no turnover
  expect_true(at(sn, 3, 4) > 0)
  expect_identical(at(st, 7, 8), 1)                     # disjoint: all turnover
  expect_identical(at(sn, 7, 8), 0)
  expect_true(abs(at(sn, 9, 10)) < 1e-15)               # equal richness
  expect_identical(at(st, 5, 6), 0)                     # identical
  expect_identical(at(sn, 5, 6), 0)
  expect_true(is.nan(at(st, 1, 2)))                     # one empty sample
  expect_true(is.nan(at(sn, 1, 2)))
  expect_identical(at(sorensen(e), 1, 2), 1)
  expect_identical(at(bray_balanced(e), 11, 12), 1)     # all balanced
  expect_identical(at(bray_gradient(e), 11, 12), 0)
  expect_identical(at(bray_balanced(fx$tables$sparse), 9, 10), 0)   # all gradient
  expect_true(is.nan(at(soergel_gradient(e), 2, 1)))

  # Both samples empty.
  z <- matrix(0, 3, 4); z[3, 1] <- 1
  expect_true(is.nan(at(jaccard_turnover(z), 1, 2)))
  expect_true(is.nan(at(bray_gradient(z), 1, 2)))

  # The components are not metrics: a sample nested in two others.
  tri <- rbind(c(1, 0), c(1, 1), c(0, 1))
  expect_identical(as.vector(sorensen_turnover(tri)), c(0, 1, 0))



  # Arguments ====

  expect_error(beta_partition(x, 'horn'), '`index` must be one of')
  expect_error(beta_partition(x, c('bray', 'jaccard')), '`index` must be one of')
  expect_error(beta_partition(x, 'jaccard', norm = 'percent'), 'only used by the abundance')
  expect_error(beta_partition(x, 'bray', norm = 'clr'), 'negative values')
  expect_error(bray_gradient(x, norm = 'rclr'), 'negative values')
  expect_error(soergel_balanced(x, norm = 'cl'), 'negative values')
  expect_error(bray_balanced(x, norm = 'nonsense'))
  expect_error(sorensen_turnover(-x), 'negative values')
  expect_error(beta_partition(-x, 'bray'), 'negative values')
  expect_error(beta_partition_multi(x[1, , drop = FALSE]), 'at least two samples')
  expect_error(beta_partition(x[1, , drop = FALSE]), 'at least two samples')
  expect_error(beta_partition_multi(x, n_sites = 1), '`n_sites` must be NULL')
  expect_error(beta_partition_multi(x, n_sites = nrow(x) + 1), '`n_sites` must be NULL')
  expect_error(beta_partition_multi(x, n_sites = 2.5), '`n_sites` must be NULL')
  expect_error(beta_partition_multi(x, n_sites = 3, n_draws = 0), '`n_draws` must be')
  expect_error(beta_partition_multi(x, n_sites = 3, seed = 'a'), '`seed`')
  expect_silent(beta_partition_multi(x, n_sites = 2, n_draws = 5))

  # norm = 'percent' is the partition of relative abundances.
  xp <- x[-2, ] / rowSums(x[-2, ])
  expect_equal(beta_partition(x[-2, ], 'bray', norm = 'percent'), beta_partition(xp, 'bray'))
  expect_equal(bray_balanced(x[-2, ], norm = 'percent'), bray_balanced(xp))
  expect_equal(
    beta_partition_multi(x[-2, ], 'ruzicka', norm = 'percent'),
    beta_partition_multi(xp, 'ruzicka') )

  # Binary abundances are the incidence partition.
  expect_equal(
    as.vector(bray_balanced(x, norm = 'binary')), as.vector(sorensen_turnover(x)) )

  # Samples as columns.
  expect_identical(beta_partition(t(x), 'jaccard', margin = 2L), beta_partition(x, 'jaccard'))
  expect_identical(
    beta_partition_multi(t(x), 'bray', margin = 2L, n_sites = 5, n_draws = 20),
    beta_partition_multi(x, 'bray', n_sites = 5, n_draws = 20) )



  # pairs ====

  # (Pairs 40 to 77 have the empty sample 2, and NaN components.)
  bp <- beta_partition(x, 'bray', pairs = c(3, 85, 300))
  full <- beta_partition(x, 'bray')
  for (k in 1:3) {
    expect_identical(which(!is.na(bp[[k]])), c(3L, 85L, 300L))
    expect_identical(bp[[k]][c(3, 85, 300)], full[[k]][c(3, 85, 300)])
  }
  expect_identical(
    as.vector(soergel_gradient(x, pairs = 7:9))[7:9], as.vector(soergel_gradient(x))[7:9] )
  expect_true(all(is.na(beta_partition(x, pairs = integer(0))$turnover)))



  # clust: all three in the order of the total's clustering ====

  bc <- beta_partition(x[-2, ], 'jaccard', clust = 'average')
  expect_identical(bc$total, jaccard(x[-2, ], clust = 'average'))
  ord <- attr(bc$total, 'order')
  for (k in 1:2) {
    expect_identical(attr(bc[[k]], 'order'), ord)
    expect_identical(attr(bc[[k]], 'Labels'), attr(bc$total, 'Labels'))
    expect_identical(get_display(bc[[k]], 'clust'), 'average')
  }
  unclustered <- beta_partition(x[-2, ], 'jaccard')
  expect_equal(
    as.matrix(bc$turnover)[order(ord), order(ord)],
    as.matrix(unclustered$turnover), check.attributes = FALSE )
  expect_identical(get_display(beta_partition(x, clust = FALSE)$turnover, 'clust'), 'none')
  expect_error(beta_partition(x, clust = 'average', pairs = 1:3), 'cannot be used with `pairs`')
  expect_error(beta_partition(x, clust = 'nope'), '`clust` must be')



  # Threads, sparse input ====

  for (fam in fams) {
    r1 <- beta_partition(x, fam, cpus = 1L)
    for (cpus in c(2L, 4L, 7L))
      expect_identical(beta_partition(x, fam, cpus = cpus), r1, info = paste(fam, cpus))
    m1 <- beta_partition_multi(x, fam, n_sites = 7, n_draws = 250, seed = 11, cpus = 1L)
    for (cpus in c(2L, 4L, 7L))
      expect_identical(
        beta_partition_multi(x, fam, n_sites = 7, n_draws = 250, seed = 11, cpus = cpus), m1,
        info = paste(fam, cpus) )
  }
  for (m in c('sorensen_turnover', 'jaccard_nestedness', 'bray_gradient', 'soergel_balanced'))
    expect_identical(beta_div(x, m, cpus = 4L), beta_div(x, m, cpus = 1L), info = m)

  # Enough samples for the rows of the all-samples value to run on threads.
  set.seed(5)
  big <- matrix(rpois(160 * 25, 0.7), 160, 25)
  for (fam in fams) {
    v1 <- beta_partition_multi(big, fam, cpus = 1L)
    for (cpus in c(2L, 4L, 7L))
      expect_identical(beta_partition_multi(big, fam, cpus = cpus), v1, info = paste(fam, cpus))
    expect_identical(
      unname(beta_partition_multi(big, fam, n_sites = 160, n_draws = 1)$draws[1, ]),
      unname(v1$value), info = fam )
  }

  # Another seed, other draws.
  expect_false(identical(
    beta_partition_multi(x, n_sites = 5, n_draws = 20, seed = 1)$draws,
    beta_partition_multi(x, n_sites = 5, n_draws = 20, seed = 2)$draws ))

  if (requireNamespace('Matrix', quietly = TRUE)) {
    sx <- Matrix::Matrix(x, sparse = TRUE)
    for (fam in fams) {
      expect_identical(beta_partition(sx, fam), beta_partition(x, fam), info = fam)
      expect_identical(
        beta_partition_multi(sx, fam, n_sites = 6, n_draws = 30),
        beta_partition_multi(x, fam, n_sites = 6, n_draws = 30), info = fam )
    }
    expect_identical(bray_gradient(sx), bray_gradient(x))
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    sx <- slam::as.simple_triplet_matrix(x)
    expect_identical(beta_partition(sx, 'ruzicka'), beta_partition(x, 'ruzicka'))
    expect_identical(beta_partition_multi(sx, 'sorensen'), beta_partition_multi(x, 'sorensen'))
  }



  # Names, aliases and the registry ====

  expect_identical(attr(beta_div(x, 'betasim'), 'method'), 'sorensen_turnover')
  expect_identical(attr(beta_div(x, 'Simpson dissimilarity'), 'method'), 'sorensen_turnover')
  expect_identical(attr(beta_div(x, 'betajne'), 'method'), 'jaccard_nestedness')
  expect_identical(attr(beta_div(x, 'bcbal'), 'method'), 'bray_balanced')
  expect_identical(attr(beta_div(x, 'Ruzicka gradient'), 'method'), 'soergel_gradient')
  expect_identical(attr(beta_div(x, 'Bray-Curtis'), 'method'), 'bray')
  expect_identical(attr(beta_div(x, 'BC'), 'method'), 'bray')
  expect_identical(attr(beta_div(x, 'sorensen'), 'method'), 'sorensen')
  expect_identical(attr(beta_div(x, 'jaccard'), 'method'), 'jaccard')
  expect_identical(attr(beta_div(x, 'soergel'), 'method'), 'soergel')
  expect_identical(match_metric('simpson')$id, 'simpson')            # alpha
  expect_error(beta_div(x, 'ruzicka'), 'matches multiple options')  # no alias of soergel

  # Every alias of the new metrics names one entry, in every registry.
  new_ids <- c(
    'sorensen_turnover', 'sorensen_nestedness', 'jaccard_turnover', 'jaccard_nestedness',
    'bray_balanced', 'bray_gradient', 'soergel_balanced', 'soergel_gradient' )
  new_aliases <- c(
    'simpsondissimilarity', 'betasim', 'betasne', 'betajtu', 'betajne', 'bcbal',
    'braybal', 'bcgra', 'braygra', 'ruzickabalanced', 'ruzbal', 'ruzickagradient',
    'ruzgra', 'braycurtis', 'bc' )
  for (reg in c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'ZERO_REPLACEMENTS', 'DA_TESTS')) {
    ex <- exact_names(reg)
    for (a in c(new_aliases, gsub('_', '', new_ids))) {
      hits <- unique(unname(ex[names(ex) == a]))
      expect_true(length(hits) <= (reg == 'METRICS'), info = paste(reg, a))
    }
  }
  expect_true(all(new_ids %in% list_metrics('beta', 'id', true_metric = FALSE)))
  expect_true(all(new_ids %in% list_metrics('beta', 'id', phylo = FALSE)))
  expect_identical(
    sort(new_ids[list_metrics('beta', 'weighted', nm = 'id')[new_ids]]),
    c('bray_balanced', 'bray_gradient', 'soergel_balanced', 'soergel_gradient') )



  # Titles, print and plot ====

  expect_identical(display_title(sorensen_turnover(x)), 'Simpson Dissimilarity (Sorensen Turnover)')
  bp <- beta_partition(x, 'bray')
  expect_identical(display_title(bp$balanced), 'Bray-Curtis Balanced Variation')
  expect_identical(display_title(bp$total),    'Bray-Curtis Dissimilarity')

  out <- capture.output(bp)
  expect_true(grepl('Partition of Bray-Curtis Dissimilarity: 40 samples, 780 pairs', out[[1]]))
  expect_true(any(grepl('\\$gradient +Abundance gradient +mean', out)))
  out <- capture.output(beta_partition_multi(x, n_sites = 5, n_draws = 10))
  expect_true(grepl('Multiple-site partition of Sorensen Dissimilarity: 40 samples', out[[1]]))
  expect_true(any(grepl('10 draws of 5: mean \\(sd\\)', out)))
  out <- capture.output(beta_partition_multi(x, 'jaccard'))
  expect_true(any(grepl('^  Turnover', out)))

  # The plots' titles, from `index`.
  plot_text <- function (expr) {
    grDevices::pdf(NULL)
    grDevices::dev.control('enable')
    on.exit(grDevices::dev.off())
    res <- withVisible(expr)
    rec <- grDevices::recordPlot()
    list(res = res, txt = unlist(lapply(rec[[1]], function (op) {
      a <- op[[2]]
      unlist(lapply(a, function (v) if (is.character(v)) v))
    })))
  }

  p <- plot_text(plot(bp))
  expect_false(p$res$visible)
  expect_identical(p$res$value, bp)
  expect_true('Partition of Bray-Curtis Dissimilarity' %in% p$txt)

  p <- plot_text(plot(beta_partition(x), groups = rep(c('a', 'b'), 20), main = 'mine'))
  expect_true('mine' %in% p$txt)
  expect_true('Between groups' %in% p$txt)
  p <- plot_text(plot(beta_partition(x[-2, ], clust = 'average'), groups = rep(c('a', 'b'), length.out = 39)))
  expect_true('Partition of Sorensen Dissimilarity' %in% p$txt)
  expect_error(plot(beta_partition(x), groups = 1:3), '`groups` must have one value per sample')

  bm <- beta_partition_multi(x, 'ruzicka', n_sites = 5, n_draws = 30)
  p  <- plot_text(plot(bm))
  expect_false(p$res$visible)
  expect_true('Multiple-site partition of Ruzicka Dissimilarity' %in% p$txt)
  p  <- plot_text(plot(beta_partition_multi(x, 'jaccard'), legend = NULL))
  expect_true('Multiple-site partition of Jaccard Dissimilarity' %in% p$txt)
  # Draws with no spread.
  p  <- plot_text(plot(beta_partition_multi(x, n_sites = nrow(x), n_draws = 3)))
  expect_true('Multiple-site partition of Sorensen Dissimilarity' %in% p$txt)

#})
