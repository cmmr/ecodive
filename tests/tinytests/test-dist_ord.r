#test_that("ordination", {

  set.seed(1)
  m30 <- matrix(rpois(30 * 20, 3), 30, 20, dimnames = list(paste0('S', 1:30), NULL))
  d   <- bray(m30)


  # Same as calling the method's function. ====

  expect_identical(dist_ord(d, 'pcoa'),                 pcoa(d))
  expect_identical(dist_ord(d, 'pcoa', k = 2, cpus = 1), pcoa(d, k = 2, cpus = 1))
  expect_identical(dist_ord(d, 'nmds', n_starts = 1),    nmds(d, n_starts = 1))
  expect_identical(dist_ord(d, 'mds',  n_starts = 1),    mds(d, n_starts = 1))
  expect_identical(dist_ord(d, 'nmds', k = 3, seed = 2, n_starts = 2), nmds(d, k = 3, seed = 2, n_starts = 2))
  expect_identical(dist_ord(d, 'tsne', perplexity = 5, max_iters = 100), tsne(d, perplexity = 5, max_iters = 100))

  # Arguments the method doesn't take are an error, except `seed`, which is
  # dropped for methods without random parts.
  expect_error(dist_ord(d, 'pcoa', n_starts = 3), 'unused argument')
  expect_identical(dist_ord(d, 'pcoa', seed = 1, cpus = 1), pcoa(d, cpus = 1))
  expect_identical(dist_ord(dm = d, method = 'pcoa'), pcoa(d))


  # Every method returns the common attributes first, in the same order. ====

  # Sample metadata, for methods that need it.
  meta <- data.frame(row.names = rownames(m30), g = rep(c('a', 'b', 'c'), 10))
  args <- list(
    capscale = list(groups = 'g', data = meta),
    dbrda    = list(groups = 'g', data = meta) )

  for (id in list_dist_ords(val = 'id')) {
    # The same arguments, `seed` included, work for every method.
    ord <- do.call(dist_ord, c(list(d, id, k = 2, seed = 1, cpus = 1), args[[id]]))
    expect_true(is.matrix(ord) && is.double(ord))
    expect_identical(rownames(ord), attr(d, 'Labels'))
    expect_identical(names(attributes(ord))[1:5], c('dim', 'dimnames', 'method', 'eig', 'prop'))
    expect_identical(attr(ord, 'method'), id)
    expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
    expect_identical(ncol(ord), 2L)
    expect_identical(length(attr(ord, 'eig')),  ncol(ord))
    expect_identical(length(attr(ord, 'prop')), ncol(ord))
    expect_true(is.double(attr(ord, 'eig')) && is.double(attr(ord, 'prop')))
    expect_identical(names(attributes(unclass(ord)[, ])), c('dim', 'dimnames'))
  }


  # Flexible matching. ====

  m <- function (x) match_dist_ord(x)$id
  expect_identical(m('pcoa'),                'pcoa')
  expect_identical(m('PCoA'),                'pcoa')
  expect_identical(m('Principal coord'),     'pcoa')
  expect_identical(m('MDS'),                 'mds')
  expect_identical(m('metric MDS'),          'mds')
  expect_identical(m('pco'),                 'pcoa')
  expect_identical(m('classical MDS'),       'pcoa')
  expect_identical(m('cmdscale'),            'pcoa')
  expect_identical(m('NMDS'),                'nmds')
  expect_identical(m('non-metric MDS'),      'nmds')
  expect_identical(m('metaMDS'),             'nmds')
  expect_identical(m('monoMDS'),             'nmds')
  expect_identical(m('isoMDS'),              'nmds')
  expect_identical(m('dbRDA'),               'dbrda')
  expect_identical(m('distance-based RDA'),  'dbrda')
  expect_identical(m('capscale'),            'capscale')
  expect_identical(m('CAP'),                 'capscale')
  expect_identical(m('t-SNE'),               'tsne')
  expect_identical(m('Rtsne'),               'tsne')
  expect_identical(m('Barnes-Hut t-SNE'),    'tsne')
  expect_identical(m('Isomap'),              'isomap')
  expect_identical(m('isometric mapping'),   'isomap')
  expect_identical(match_dist_ord('nmds', input = 'dist')$id, 'nmds')

  info <- match_dist_ord('pcoa')
  expect_identical(info$func, pcoa)
  expect_identical(info$params, names(formals(pcoa)))
  expect_identical(info$input, 'dist')
  expect_null(info$alt_ids)


  # Listing. ====

  df <- list_dist_ords()
  expect_true(is.data.frame(df))
  expect_identical(names(df), c('name', 'id', 'input'))
  expect_identical(df$id, c('capscale', 'dbrda', 'isomap', 'mds', 'nmds', 'pcoa', 'tsne', 'umap'))
  expect_identical(rownames(list_dist_ords(nm = 'id')), c('capscale', 'dbrda', 'isomap', 'mds', 'nmds', 'pcoa', 'tsne', 'umap'))

  expect_identical(list_dist_ords(val = 'id', nm = 'name'),
    c('Constrained Analysis of Principal Coordinates' = 'capscale',
      'Distance-based Redundancy Analysis' = 'dbrda', 'Isometric Feature Mapping' = 'isomap',
      'Metric Multidimensional Scaling' = 'mds',
      'Non-metric Multidimensional Scaling' = 'nmds', 'Principal Coordinates Analysis' = 'pcoa',
      't-distributed Stochastic Neighbor Embedding' = 'tsne',
      'Uniform Manifold Approximation and Projection' = 'umap'))
  expect_identical(list_dist_ords(val = 'func'), list(capscale = capscale, dbrda = dbrda, isomap = isomap, mds = mds, nmds = nmds, pcoa = pcoa, tsne = tsne, umap = umap))
  expect_identical(list_dist_ords(val = 'list')$pcoa, match_dist_ord('pcoa'))
  expect_identical(list_dist_ords(val = 'input', input = 'dist'), rep('dist', 8))


  # Invalid input. ====

  expect_error(dist_ord(d, 'xyz'),   'does not match any option: capscale, dbrda, isomap, mds, nmds, pcoa, tsne, umap')
  expect_error(dist_ord(d, 'm'),     'matches multiple options: mds, nmds')
  expect_error(dist_ord(d, 123),     'character')
  expect_error(dist_ord(d, c('pcoa', 'nmds')), 'length 1')
  expect_error(dist_ord(d, '-'),     'alphabetic')
  expect_error(dist_ord(d, NA_character_), 'NA')
  expect_error(dist_ord(d, 'PCA'),       'call pca\\(\\) directly')
  expect_error(dist_ord(d, 'adonis2'),   'is a test, not an ordination: use dist_test')
  expect_error(dist_ord(d, 'betadisper'), 'use dist_test')
  expect_error(list_dist_ords(input = 'counts'), 'Options are: NULL, dist')
  expect_error(list_dist_ords(val = 'xyz'))

#})
