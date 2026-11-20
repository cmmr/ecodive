#test_that("display", {

  # The `display` list of results (?ecodive_display): a title for plot() and
  # print(), and a distance matrix's `clust`. Results record nothing else
  # about how they were made.

  cpus <- 1L
  x    <- iris[, 1:4]
  rownames(x) <- paste0('s', 1:150)
  meta <- data.frame(Species = iris$Species, row.names = rownames(x))

  # The helpers.
  y <- set_display(structure(1:3, a = 1), title = 'T', clust = 'none')
  expect_identical(names(attributes(y)), c('a', 'display'))
  expect_identical(attr(y, 'display'), list(title = 'T', clust = 'none'))
  y <- set_display(structure(y, b = 2), title = NULL)
  expect_identical(names(attributes(y)), c('a', 'b', 'display'))
  expect_identical(attr(y, 'display'), list(clust = 'none'))
  expect_null(attr(set_display(y, clust = NULL), 'display'))
  expect_identical(get_display(y, 'clust'), 'none')
  expect_null(get_display(1:3, 'clust'))
  expect_null(display_title(y))
  expect_null(display_title(set_display(1, title = '')))
  expect_identical(display_title(set_display(1, title = 'T')), 'T')
  expect_identical(plot_title(1, 'F'), 'F')
  expect_identical(title_of('PCoA', 'Thing'), 'PCoA of Thing')
  expect_identical(title_of('PCoA', NULL), 'PCoA')
  expect_identical(title_of('PCoA', NA_character_), 'PCoA')

  # Distances: `method` (the metric's ID) and the title.
  d <- bray(x, norm = 'percent', cpus = cpus)
  expect_identical(names(attributes(d)), c('class', 'Size', 'Diag', 'Upper', 'Labels', 'method', 'display'))
  expect_identical(attr(d, 'method'), 'bray')
  expect_identical(attr(d, 'display'), list(title = 'Bray-Curtis Dissimilarity'))
  expect_identical(stats::hclust(d)$dist.method, 'bray')

  # Clustering: `order` at the top level, `clust` in the display list.
  dc <- bray(x, clust = 'average', cpus = cpus)
  expect_identical(names(attributes(dc)), c('class', 'Size', 'Diag', 'Upper', 'Labels', 'method', 'order', 'display'))
  expect_identical(attr(dc, 'display'), list(title = 'Bray-Curtis Dissimilarity', clust = 'average'))

  # Arithmetic makes new values: no metric or title, but `clust` and `order`
  # still hold.
  for (s in list(sqrt(dc), 1 - dc, dc * 2, -dc, round(dc, 2))) {
    expect_null(attr(s, 'method'))
    expect_identical(attr(s, 'display'), list(clust = 'average'))
    expect_identical(attr(s, 'order'), attr(dc, 'order'))
    expect_inherits(s, 'ecodive_dist')
    expect_identical(dist_plot_name(s)$title, 'Distances')
  }
  expect_null(attr(sqrt(d), 'display'))
  expect_null(attr(sqrt(d), 'method'))
  # 1 - rho_p() is no longer taken for a similarity.
  r <- rho_p(ecodive::ex_counts + 1, cpus = cpus)
  expect_identical(dist_metric_of(r), 'rho_p')
  expect_identical(dist_metric_of(1 - r), '')

  # Ordinations: the method's name, then the distances' title.
  o <- pcoa(d, n_axes = 2, cpus = cpus)
  expect_identical(tail(names(attributes(o)), 1), 'display')
  expect_identical(attr(o, 'display'), list(title = 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity'))
  expect_null(attr(o, 'dist'))
  expect_null(attr(o, 'ord'))
  expect_true(any(capture.output(print(o)) == 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity'))
  expect_identical(attr(pcoa(dc, cpus = cpus), 'display')$title, 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity')
  for (m in list_dist_ords('id')) {
    args <- list(d, m, cpus = cpus)
    if ('groups' %in% names(formals(get(m)))) args <- c(args, list(groups = 'Species', data = meta))
    if (m == 'isomap') args$n_neighbors <- 50
    r <- do.call(dist_ord, args)
    expect_identical(
      attr(r, 'display')$title,
      paste(list_dist_ords()$name[list_dist_ords()$id == m], 'of Bray-Curtis Dissimilarity'), info = m )
  }

  # Without a title to name: plain matrices, stats::dist(), counts.
  expect_identical(attr(pcoa(as.matrix(d), cpus = cpus), 'display')$title, 'Principal Coordinates Analysis')
  expect_identical(attr(pcoa(stats::dist(x), cpus = cpus), 'display')$title, 'Principal Coordinates Analysis')
  expect_identical(attr(pca(x, cpus = cpus), 'display')$title, 'Principal Component Analysis')
  expect_identical(attr(procrustes(o, pca(x, n_axes = 2, cpus = cpus), cpus = cpus), 'display')$title, 'Procrustes Rotation')

  # Tests: the test's name, then the distances' title, also print()'s
  # heading; mantel() names both matrices.
  p <- permanova(d, 'Species', meta, permutations = 9, cpus = cpus)
  expect_identical(attr(p, 'display')$title,
    'Permutational multivariate analysis of variance (PERMANOVA) of Bray-Curtis Dissimilarity')
  expect_identical(capture.output(print(p))[[1]], attr(p, 'display')$title)
  expect_identical(test_title(p), attr(p, 'display')$title)
  for (t in list(
      anosim(d, 'Species', meta, permutations = 9, cpus = cpus),
      mrpp(d, 'Species', meta, permutations = 9, cpus = cpus),
      suppressWarnings(permdisp(d, 'Species', meta, permutations = 9, cpus = cpus)) ))
    expect_true(endsWith(attr(t, 'display')$title, ' of Bray-Curtis Dissimilarity'), info = class(t)[[1]])
  e <- euclidean(x[, 1:2], cpus = cpus)
  expect_identical(attr(mantel(d, e, permutations = 9, cpus = cpus), 'display')$title,
    'Mantel test (Pearson) of Bray-Curtis Dissimilarity and Euclidean Distance')
  expect_identical(attr(mantel(d, x[, 1], permutations = 9, cpus = cpus), 'display')$title,
    'Mantel test (Pearson) of Bray-Curtis Dissimilarity')
  expect_identical(attr(mantel(as.matrix(d), e, permutations = 9, cpus = cpus), 'display')$title, 'Mantel test (Pearson)')
  expect_identical(attr(permanova(as.matrix(d), 'Species', meta, permutations = 0, cpus = cpus), 'display')$title,
    'Permutational multivariate analysis of variance (PERMANOVA)')

  # envfit(): the ordination's title.
  expect_identical(attr(envfit(o, meta, permutations = 9, cpus = cpus), 'display')$title,
    'Fitted variables on Principal Coordinates Analysis of Bray-Curtis Dissimilarity')
  expect_identical(attr(envfit(unclass(o)[, 1:2], meta, permutations = 9, cpus = cpus), 'display')$title, 'Fitted variables')

  # A title set by the user is used.
  attr(o, 'display')$title <- 'Iris flowers'
  expect_identical(plot_title(o, 'x'), 'Iris flowers')

#})
