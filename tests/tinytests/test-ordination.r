#test_that("ordination", {

  set.seed(1)
  m30 <- matrix(rpois(30 * 20, 3), 30, 20, dimnames = list(paste0('S', 1:30), NULL))
  d   <- bray(m30)


  # Same as calling the method's function. ====

  expect_identical(ordination(d, 'pcoa'),                 pcoa(d))
  expect_identical(ordination(d, 'pcoa', k = 2, cpus = 1), pcoa(d, k = 2, cpus = 1))
  expect_identical(ordination(d, 'nmds', trymax = 3),      nmds(d, trymax = 3))
  expect_identical(ordination(d, 'nmds', k = 3, seed = 2, trymax = 2), nmds(d, k = 3, seed = 2, trymax = 2))

  # Arguments the method doesn't take.
  expect_error(ordination(d, 'pcoa', trymax = 3), 'unused argument')


  # Every method returns the common attributes first, in the same order. ====

  # Sample metadata, for methods that need it.
  meta <- data.frame(row.names = rownames(m30), g = rep(c('a', 'b', 'c'), 10))
  args <- list(dbrda = list(data = meta))

  for (id in list_ordinations(val = 'id')) {
    ord <- do.call(ordination, c(list(d, id, cpus = 1), args[[id]]))
    expect_true(is.matrix(ord) && is.double(ord))
    expect_identical(rownames(ord), attr(d, 'Labels'))
    expect_identical(names(attributes(ord))[1:5], c('dim', 'dimnames', 'method', 'eig', 'prop'))
    expect_identical(attr(ord, 'method'), id)
    expect_identical(length(attr(ord, 'eig')),  ncol(ord))
    expect_identical(length(attr(ord, 'prop')), ncol(ord))
  }


  # Flexible matching. ====

  m <- function (x) match_ordination(x)$id
  expect_identical(m('pcoa'),                'pcoa')
  expect_identical(m('PCoA'),                'pcoa')
  expect_identical(m('Principal coord'),     'pcoa')
  expect_identical(m('MDS'),                 'pcoa')
  expect_identical(m('pco'),                 'pcoa')
  expect_identical(m('classical MDS'),       'pcoa')
  expect_identical(m('cmdscale'),            'pcoa')
  expect_identical(m('NMDS'),                'nmds')
  expect_identical(m('non-metric MDS'),      'nmds')
  expect_identical(m('metaMDS'),             'nmds')
  expect_identical(m('monoMDS'),             'nmds')
  expect_identical(m('dbRDA'),               'dbrda')
  expect_identical(m('distance-based RDA'),  'dbrda')
  expect_identical(match_ordination('nmds', input = 'dist')$id, 'nmds')

  info <- match_ordination('pcoa')
  expect_identical(info$func, pcoa)
  expect_identical(info$params, names(formals(pcoa)))
  expect_identical(info$input, 'dist')
  expect_null(info$alt_ids)


  # Listing. ====

  df <- list_ordinations()
  expect_true(is.data.frame(df))
  expect_identical(names(df), c('name', 'id', 'input'))
  expect_identical(df$id, c('dbrda', 'nmds', 'pcoa', 'umap'))
  expect_identical(rownames(list_ordinations(nm = 'id')), c('dbrda', 'nmds', 'pcoa', 'umap'))

  expect_identical(list_ordinations(val = 'id', nm = 'name'),
    c('Distance-based Redundancy Analysis' = 'dbrda',
      'Non-metric Multidimensional Scaling' = 'nmds', 'Principal Coordinates Analysis' = 'pcoa',
      'Uniform Manifold Approximation and Projection' = 'umap'))
  expect_identical(list_ordinations(val = 'func'), list(dbrda = dbrda, nmds = nmds, pcoa = pcoa, umap = umap))
  expect_identical(list_ordinations(val = 'list')$pcoa, match_ordination('pcoa'))
  expect_identical(list_ordinations(val = 'input', input = 'dist'), c('dist', 'dist', 'dist', 'dist'))


  # Invalid input. ====

  expect_error(ordination(d, 'xyz'),   'does not match any option: dbrda, nmds, pcoa, umap')
  expect_error(ordination(d, 'm'),     'matches multiple options: nmds, pcoa')
  expect_error(ordination(d, 123),     'character')
  expect_error(ordination(d, c('pcoa', 'nmds')), 'length 1')
  expect_error(ordination(d, '-'),     'alphabetic')
  expect_error(list_ordinations(input = 'counts'), 'Options are: NULL, dist')
  expect_error(list_ordinations(val = 'xyz'))

#})
