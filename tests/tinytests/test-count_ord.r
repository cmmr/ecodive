#test_that("count_ord", {

  set.seed(1)
  m <- matrix(rpois(30 * 20, 3), 30, 20, dimnames = list(paste0('S', 1:30), paste0('F', 1:20)))
  m[, 1] <- m[, 1] + 1L


  # Same as calling the method's function. ====

  expect_identical(count_ord(m, 'pca', cpus = 1),                pca(m, cpus = 1))
  expect_identical(count_ord(m, 'pca', k = 2, norm = 'clr', pseudocount = 1, cpus = 1),
                   pca(m, k = 2, norm = 'clr', pseudocount = 1, cpus = 1))
  expect_identical(count_ord(m, 'rpca', k = 2, cpus = 1),        rpca(m, k = 2, cpus = 1))
  expect_identical(count_ord(counts = m, method = 'pca', cpus = 1), pca(m, cpus = 1))

  # Arguments the method doesn't take are an error, except `seed`, which is
  # dropped for methods without random parts.
  expect_error(count_ord(m, 'pca', n_starts = 3), 'unused argument')
  expect_identical(count_ord(m, 'pca', seed = 1, cpus = 1), pca(m, cpus = 1))

  # Every method returns the common attributes first, in the same order.
  for (id in list_count_ords(val = 'id')) {
    ord <- count_ord(m, id, k = 2, seed = 1, cpus = 1)
    expect_true(is.matrix(ord) && is.double(ord), info = id)
    expect_identical(rownames(ord), rownames(m), info = id)
    expect_identical(names(attributes(ord))[1:5], c('dim', 'dimnames', 'method', 'eig', 'prop'), info = id)
    expect_identical(attr(ord, 'method'), id, info = id)
    expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'), info = id)
    expect_identical(ncol(ord), 2L, info = id)
  }


  # Flexible matching. ====

  mo <- function (x) match_count_ord(x)$id
  expect_identical(mo('pca'),                 'pca')
  expect_identical(mo('PCA'),                 'pca')
  expect_identical(mo('principal comp'),      'pca')
  expect_identical(mo('prcomp'),              'pca')
  expect_identical(mo('RPCA'),                'rpca')
  expect_identical(mo('robust PCA'),          'rpca')
  expect_identical(mo('Robust Aitchison PCA'), 'rpca')
  expect_identical(mo('DEICODE'),             'rpca')
  expect_identical(mo('gemelli'),             'rpca')

  info <- match_count_ord('pca')
  expect_identical(info$func, pca)
  expect_identical(info$params, names(formals(pca)))
  expect_identical(info$name, 'Principal Component Analysis')
  expect_null(info$alt_ids)


  # Listing. ====

  df <- list_count_ords()
  expect_true(is.data.frame(df))
  expect_identical(names(df), c('name', 'id'))
  expect_identical(df$id, c('pca', 'rpca'))
  expect_identical(rownames(list_count_ords(nm = 'id')), c('pca', 'rpca'))
  expect_identical(list_count_ords(val = 'id', nm = 'name'),
    c('Principal Component Analysis' = 'pca', 'Robust Aitchison PCA' = 'rpca'))
  expect_identical(list_count_ords(val = 'func'), list(pca = pca, rpca = rpca))
  expect_identical(list_count_ords(val = 'list')$pca, match_count_ord('pca'))
  expect_identical(list_count_ords('name'), c('Principal Component Analysis', 'Robust Aitchison PCA'))


  # Invalid input, and methods of the other wrappers. ====

  expect_error(count_ord(m, 'xyz'),   'does not match any option: pca, rpca')
  expect_error(count_ord(m, 123),     'character')
  expect_error(count_ord(m, c('pca', 'rpca')), 'length 1')
  expect_error(count_ord(m, '-'),     'alphabetic')
  expect_error(count_ord(m, NA_character_), 'NA')
  expect_error(count_ord(m, 'PCoA'),     'takes distances, not a count table. Use dist_ord\\(dm, "pcoa"\\)')
  expect_error(count_ord(m, 'cmdscale'), 'dist_ord\\(dm, "pcoa"\\)')
  expect_error(count_ord(m, 'NMDS'),     'dist_ord\\(dm, "nmds"\\)')
  expect_error(count_ord(m, 'adonis2'),  'is a test, not an ordination: use dist_test')
  expect_error(list_count_ords(val = 'xyz'))

#})
