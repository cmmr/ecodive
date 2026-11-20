#test_that("clust", {

  # The `clust` argument of the functions that return a distance matrix:
  # it reorders the samples by a clustering, and records itself and the
  # order in attributes.

  cpus <- 1L
  x    <- iris[, 1:4]
  rownames(x) <- paste0('s', 1:150)
  d0   <- bray(x, cpus = cpus)
  m0   <- as.matrix(d0)

  # Every metric function, and the wrapper, takes `clust`, before `cpus`.
  for (f in c(list_metrics('beta', 'func'), list(beta_div = beta_div, vlr = vlr, rho_p = rho_p, phi_s = phi_s))) {
    a <- names(formals(f))
    expect_true('clust' %in% a)
    expect_identical(a[[length(a)]], 'cpus')
    expect_null(formals(f)$clust)
  }

  # NULL, the default, changes nothing.
  expect_null(attr(d0, 'display')$clust)
  expect_null(attr(d0, 'order'))
  expect_identical(bray(x, clust = NULL, cpus = cpus), d0)

  # FALSE and 'none' keep the order, and record 'none'.
  for (cl in list(FALSE, 'none')) {
    d <- bray(x, clust = cl, cpus = cpus)
    expect_identical(attr(d, 'display')$clust, 'none')
    expect_null(attr(d, 'order'))
    expect_identical(as.vector(d), as.vector(d0))
    expect_identical(labels(d), labels(d0))
  }

  # Each linkage: the samples in hclust()'s leaf order, with their labels.
  for (cl in c('average', 'complete', 'single', 'ward.D', 'ward.D2', 'mcquitty', 'median', 'centroid')) {
    d <- bray(x, clust = cl, cpus = cpus)
    o <- attr(d, 'order')
    expect_identical(o, stats::hclust(stats::as.dist(m0), method = cl)$order, info = cl)
    expect_identical(attr(d, 'display')$clust, cl, info = cl)
    expect_identical(class(d), c('ecodive_dist', 'dist'), info = cl)
    expect_identical(labels(d), labels(d0)[o], info = cl)
    expect_equal(unname(as.matrix(d)), unname(m0[o, o]), info = cl)
    expect_identical(attr(d, 'Size'), 150L, info = cl)
  }

  # 'pcoa': along the first axis of pcoa().
  d <- bray(x, clust = 'pcoa', cpus = cpus)
  expect_identical(attr(d, 'order'), order(pcoa(d0, n_axes = 1, cpus = cpus)[, 1]))
  expect_equal(unname(as.matrix(d)), unname(m0[attr(d, 'order'), attr(d, 'order')]))

  # The same for any number of threads.
  expect_identical(bray(x, clust = 'pcoa', cpus = 1L), bray(x, clust = 'pcoa', cpus = 2L))

  # The order can be undone, and lines metadata up.
  d  <- bray(x, clust = 'average', cpus = cpus)
  o  <- attr(d, 'order')
  dr <- dist_reorder(d, order(o))
  expect_identical(as.vector(dr), as.vector(d0))
  expect_identical(labels(dr), labels(d0))
  expect_identical(rownames(x[o, ]), labels(d))
  # Tests match named samples to `data` by name, whatever their order.
  meta <- data.frame(Species = iris$Species, row.names = rownames(x))
  expect_equal(
    permanova(d, 'Species', meta, permutations = 0, cpus = cpus)$statistic,
    permanova(d0, 'Species', meta, permutations = 0, cpus = cpus)$statistic )

  # Analyses put a clustered matrix's samples back in their input order, so
  # that metadata lines up by position too, and results match.
  cases <- list(
    list(dd = bray(iris[, 1:4], clust = 'average', cpus = cpus), ref = bray(iris[, 1:4], cpus = cpus)),
    list(dd = d, ref = d0) )
  for (cs in cases) {
    dd  <- cs$dd
    ref <- cs$ref
    expect_equal(unclass(no_display(pcoa(dd, cpus = cpus)))[, ], unclass(no_display(pcoa(ref, cpus = cpus)))[, ])
    expect_identical(rownames(pcoa(dd, cpus = cpus)), labels(ref))
    expect_identical(
      permanova(dd, iris$Species, permutations = 0, cpus = cpus)$statistic,
      permanova(ref, iris$Species, permutations = 0, cpus = cpus)$statistic )
  }
  expect_identical(attr(pcoa(d, cpus = cpus), 'display')$title, 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity')

  # The heatmap finds the groups of a clustered matrix's samples by their
  # input positions.
  dc <- bray(iris[, 1:4], clust = 'average', cpus = cpus)
  g  <- dist_plot_groups(dc, 150L, labels(dc), 'Species', iris)
  expect_identical(g$f, iris$Species[attr(dc, 'order')])
  g  <- dist_plot_groups(d, 150L, labels(d), iris$Species, NULL)
  expect_identical(g$f, iris$Species[attr(d, 'order')])

  # Samples without names are labelled with their input positions.
  d <- bray(iris[, 1:4], clust = 'average', cpus = cpus)
  expect_identical(labels(d), as.character(attr(d, 'order')))
  expect_null(labels(bray(iris[, 1:4], cpus = cpus)))

  # The wrapper and the other functions pass it on.
  d <- beta_div(x, 'jaccard', clust = 'complete', cpus = cpus)
  expect_identical(attr(d, 'display')$clust, 'complete')
  expect_identical(attr(d, 'order'), stats::hclust(stats::as.dist(as.matrix(jaccard(x, cpus = cpus))), 'complete')$order)
  ex <- ecodive::ex_counts
  for (f in list(unweighted_unifrac, weighted_unifrac, normalized_unifrac, generalized_unifrac, variance_adjusted_unifrac)) {
    d <- f(ex, tree = ecodive::ex_tree, clust = 'average', cpus = cpus)
    expect_identical(attr(d, 'display')$clust, 'average')
    expect_identical(sort(attr(d, 'order')), 1:4)
  }
  for (f in list(hellinger, jensen, matusita, psym_chisq, squared_euclidean, topsoe)) {
    d <- f(x, clust = 'single', cpus = cpus)
    expect_identical(attr(d, 'display')$clust, 'single')
    expect_identical(class(d), c('ecodive_dist', 'dist'))
  }

  # Proportionality: over the features; rho_p as the dissimilarity 1 - rho.
  for (f in list(vlr, phi_s)) {
    d0f <- f(ex + 1, cpus = cpus)
    d   <- f(ex + 1, clust = 'average', cpus = cpus)
    expect_identical(attr(d, 'order'), stats::hclust(stats::as.dist(as.matrix(d0f)), 'average')$order)
    expect_identical(labels(d), labels(d0f)[attr(d, 'order')])
    expect_identical(attr(d, 'method'), attr(d0f, 'method'))
  }
  r0 <- rho_p(ex + 1, cpus = cpus)
  r  <- rho_p(ex + 1, clust = 'average', cpus = cpus)
  expect_identical(attr(r, 'order'), stats::hclust(stats::as.dist(1 - as.matrix(r0)), 'average')$order)
  expect_equal(unname(as.matrix(r)), unname(as.matrix(r0)[attr(r, 'order'), attr(r, 'order')]))

  # Two samples: nothing to cluster.
  d <- bray(x[1:2, ], clust = 'average', cpus = cpus)
  expect_identical(attr(d, 'order'), 1:2)

  # Errors.
  expect_error(bray(x, clust = 'ward'), '`clust` must be NULL, FALSE, or one of')
  expect_error(bray(x, clust = TRUE), '`clust` must be NULL, FALSE, or one of')
  expect_error(bray(x, clust = c('average', 'single')), '`clust` must be NULL, FALSE, or one of')
  expect_error(vlr(ex + 1, clust = 1), '`clust` must be NULL, FALSE, or one of')
  expect_error(bray(x, pairs = 1:10, clust = 'average'), 'cannot be used with `pairs`')
  expect_silent(bray(x, pairs = 1:10, clust = FALSE))



  # Each matrix's metric: its ID as `method`, and its display name as its
  # title (?ecodive_display). ====

  for (f in list_metrics('beta', 'func')) {
    id <- names(which(vapply(list_metrics('beta', 'func'), identical, TRUE, f)))
    d  <- if (id %in% list_metrics('beta', 'id', phylo = TRUE))
      f(ex + 1, tree = ecodive::ex_tree, cpus = cpus) else f(ex + 1, cpus = cpus)
    expect_identical(attr(d, 'method'), id, info = id)
    expect_identical(attr(d, 'display')$title, match_metric(id)$name, info = id)
  }
  d <- bray(x, norm = 'percent', cpus = cpus)
  expect_identical(names(attributes(d)), c('class', 'Size', 'Diag', 'Upper', 'Labels', 'method', 'display'))
  expect_identical(attr(beta_div(x, 'Bray-Curtis', cpus = cpus), 'method'), 'bray')
  # Clustering keeps them.
  d <- bray(x, clust = 'average', cpus = cpus)
  expect_identical(attr(d, 'display')$title, 'Bray-Curtis Dissimilarity')
  for (f in c('vlr', 'rho_p', 'phi_s')) {
    d <- get(f)(ex + 1, cpus = cpus)
    expect_identical(attr(d, 'method'), f)
    expect_true(is.character(attr(d, 'display')$title) && nzchar(attr(d, 'display')$title))
  }
  d <- attr(rpca(cnt <- round(x * 10), distances = TRUE, cpus = cpus), 'dm')
  expect_identical(attr(d, 'method'), 'rpca')
  expect_identical(attr(d, 'display')$title, 'Robust Aitchison Distance (RPCA)')



  # plot(): clusters by average linkage when `clust` isn't set; keeps the
  # order otherwise. ====

  expect_identical(dist_clust_order(d0, 'average'), stats::hclust(stats::as.dist(m0), 'average')$order)
  expect_identical(dist_plot_order_name(NULL), NULL)
  expect_identical(dist_plot_order_name('none'), NULL)
  expect_identical(dist_plot_order_name('average'), 'clustered (average linkage)')
  expect_identical(dist_plot_order_name('pcoa'), 'ordered by PCoA axis 1')

  grDevices::pdf(NULL)
  for (d in list(d0, bray(x, clust = FALSE), bray(x, clust = 'ward.D2'), bray(x, pairs = 1:10))) {
    v <- withVisible(plot(d, 'Species', meta))
    expect_false(v$visible)
    expect_identical(v$value, d)
  }
  grDevices::dev.off()

#})
