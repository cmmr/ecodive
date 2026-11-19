#test_that("steps", {

  # The records of how results were made (?ecodive_steps): one attribute
  # per step, `dist` and `ord`, carried forward through later steps.

  cpus <- 1L
  x    <- iris[, 1:4]
  rownames(x) <- paste0('s', 1:150)

  # A record: its name, then its settings; one line in print() and str().
  st <- ecodive_step('Thing', a = 1, b = NULL, c = 'z', settings = list(d = 1:10))
  expect_inherits(st, 'ecodive_step')
  expect_identical(names(st), c('name', 'a', 'c', 'd'))
  expect_identical(format(st), 'Thing (a = 1, c = z, d = integer [1:10])')
  expect_identical(capture.output(print(st)), 'Thing (a = 1, c = z, d = integer [1:10])')
  expect_identical(capture.output(str(st)), ' Thing (a = 1, c = z, d = integer [1:10])')
  expect_identical(format(ecodive_step('Bare')), 'Bare')
  expect_identical(step_name(st), 'Thing')
  expect_null(step_name(list(name = 'x')))
  expect_identical(step_title('PCoA', st), 'PCoA of Thing')
  expect_identical(step_title('PCoA', NULL), 'PCoA')

  # Distances: the `dist` record, and nothing else besides dist's own.
  d <- bray(x, norm = 'percent', cpus = cpus)
  expect_identical(names(attributes(d)), c('class', 'Size', 'Diag', 'Upper', 'Labels', 'dist'))
  rec <- attr(d, 'dist')
  expect_identical(rec$name, 'Bray-Curtis Dissimilarity')
  expect_identical(rec$metric, 'bray')
  expect_identical(rec$norm, 'percent')
  expect_identical(rec$pseudocount, 0)
  expect_null(rec$clust)

  # Ordinations: their `ord` record, after the matrix's `dist` record.
  o <- pcoa(d, n_axes = 2, cpus = cpus)
  expect_identical(attr(o, 'dist'), rec)
  oo <- attr(o, 'ord')
  expect_inherits(oo, 'ecodive_step')
  expect_identical(oo$name, 'Principal Coordinates Analysis')
  expect_identical(oo$method, 'pcoa')
  expect_identical(oo$n_axes, 2L)
  expect_identical(oo$correction, 'none')
  expect_null(oo$cpus)
  expect_identical(tail(names(attributes(o)), 2), c('dist', 'ord'))
  expect_true(any(grepl('"dist")= Bray-Curtis Dissimilarity (metric = bray', capture.output(str(o)), fixed = TRUE)))
  expect_true(any(grepl('"ord")= Principal Coordinates Analysis (method = pcoa', capture.output(str(o)), fixed = TRUE)))

  # The same by name, through dist_ord(), and with every method.
  expect_identical(attr(dist_ord(d, 'pcoa', n_axes = 2, cpus = cpus), 'ord'), oo)
  for (m in list_dist_ords('id')) {
    args <- list(d, m, cpus = cpus)
    if ('groups' %in% names(formals(get(m)))) args <- c(args, list(groups = 'Species', data = data.frame(Species = iris$Species, row.names = rownames(x))))
    if (m == 'isomap') args$n_neighbors <- 50
    r <- do.call(dist_ord, args)
    expect_identical(attr(r, 'dist'), rec, info = m)
    expect_identical(attr(r, 'ord')$method, m, info = m)
    expect_identical(attr(r, 'ord')$name, list_dist_ords()$name[list_dist_ords()$id == m], info = m)
  }
  expect_equal(attr(nmds(d, seed = 3, cpus = cpus), 'ord')$seed, 3)
  dbr <- dbrda(d, data = iris, formula = ~ Species + Sepal.Width, cpus = cpus)
  expect_identical(attr(dbr, 'ord')$formula, '~Species + Sepal.Width')

  # Ordinations of counts: an `ord` record only, with their settings.
  p <- pca(x, norm = 'clr', pseudocount = 1, scale = TRUE, cpus = cpus)
  expect_null(attr(p, 'dist'))
  expect_identical(attr(p, 'ord')$norm, 'clr')
  expect_identical(attr(p, 'ord')$scale, TRUE)
  expect_identical(attr(lda(x, iris$Species), 'ord')$name, 'Linear Discriminant Analysis')

  # Procrustes: what was rotated onto what.
  pr <- procrustes(o, pca(x, n_axes = 2, cpus = cpus), cpus = cpus)
  expect_identical(attr(pr, 'ord')$method, 'procrustes')
  expect_identical(attr(pr, 'ord')$ord, 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity')
  expect_identical(attr(pr, 'ord')$target, 'Principal Component Analysis')
  pt <- procrustes_test(o, pca(x, n_axes = 2, cpus = cpus), permutations = 9, cpus = cpus)
  expect_identical(attr(pt$procrustes, 'ord')$ord, 'Principal Coordinates Analysis of Bray-Curtis Dissimilarity')

  # Tests keep their distances' record; envfit() its ordination's.
  meta <- data.frame(Species = iris$Species, row.names = rownames(x))
  for (t in list(
      permanova(d, 'Species', meta, permutations = 9, cpus = cpus),
      anosim(d, 'Species', meta, permutations = 9, cpus = cpus),
      mrpp(d, 'Species', meta, permutations = 9, cpus = cpus),
      suppressWarnings(permdisp(d, 'Species', meta, permutations = 9, cpus = cpus)),
      mantel(d, euclidean(x[, 1:2], cpus = cpus), permutations = 9, cpus = cpus),
      mantel_partial(d, euclidean(x[, 1:2], cpus = cpus), given = euclidean(x[, 3:4], cpus = cpus), permutations = 9, cpus = cpus) ))
    expect_identical(attr(t, 'dist'), rec, info = class(t)[[1]])
  ef <- envfit(o, meta, permutations = 9, cpus = cpus)
  expect_identical(attr(ef, 'dist'), rec)
  expect_identical(attr(ef, 'ord'), oo)

  # Clustering adds to the record.
  dc <- bray(x, clust = 'average', cpus = cpus)
  expect_identical(attr(dc, 'dist')$clust, 'average')
  expect_identical(length(attr(dc, 'dist')$order), 150L)
  expect_identical(attr(pcoa(dc, cpus = cpus), 'dist')$clust, 'average')

  # Plain matrices and stats::dist() objects have no records.
  expect_null(attr(pcoa(stats::dist(x), cpus = cpus), 'dist'))
  expect_identical(attr(pcoa(stats::dist(x), cpus = cpus), 'ord')$method, 'pcoa')

  # Plot titles name both steps.
  expect_identical(test_title(permanova(d, 'Species', meta, permutations = 0, cpus = cpus)),
    'Permutational multivariate analysis of variance (PERMANOVA) of Bray-Curtis Dissimilarity')
  expect_identical(step_title(ord_method_name(attr(o, 'method')), attr(o, 'dist')),
    'Principal Coordinates Analysis of Bray-Curtis Dissimilarity')

#})
