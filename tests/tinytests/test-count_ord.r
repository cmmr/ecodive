#test_that("count_ord", {

  set.seed(1)
  m <- matrix(rpois(30 * 20, 3), 30, 20, dimnames = list(paste0('S', 1:30), paste0('F', 1:20)))
  m[, 1] <- m[, 1] + 1L
  tr <- read_tree(paste0('((', paste0('F', 1:10, ':1', collapse = ','), '):0.5,(',
                         paste0('F', 11:20, ':0.5', collapse = ','), '):1);'))


  # Same as calling the method's function. ====

  expect_identical(count_ord(m, 'pca', cpus = 1),                pca(m, cpus = 1))
  expect_identical(count_ord(m, 'pca', n_axes = 2, norm = 'clr', pseudocount = 1, cpus = 1),
                   pca(m, n_axes = 2, norm = 'clr', pseudocount = 1, cpus = 1))
  expect_identical(count_ord(m, 'rpca', n_axes = 2, cpus = 1),        rpca(m, n_axes = 2, cpus = 1))
  expect_identical(count_ord(m, 'rpca_phylo', tree = tr, n_axes = 2, cpus = 1), rpca_phylo(m, tr, n_axes = 2, cpus = 1))
  expect_identical(count_ord(counts = m, method = 'pca', cpus = 1), pca(m, cpus = 1))
  expect_identical(count_ord(m, 'dpcoa', tree = tr, n_axes = 2, cpus = 1), dpcoa(m, tr, n_axes = 2, cpus = 1))

  # Arguments the method doesn't take are an error, except `seed`, which is
  # dropped for methods without random parts.
  expect_error(count_ord(m, 'pca', n_starts = 3), 'unused argument')
  expect_identical(count_ord(m, 'pca', seed = 1, cpus = 1), pca(m, cpus = 1))

  # Every method returns the common attributes first, in the same order.
  # lda() and plsda() need groups; the others run without them. oplsda()
  # needs two groups, and takes `n_orthogonal` rather than `n_axes`.
  # rpca_phylo() and dpcoa() need a tree.
  grp <- rep(c('a', 'b', 'c'), 10)
  for (id in list_count_ords(val = 'id')) {
    ord <- if (id %in% c('lda', 'plsda')) count_ord(m, id, grp, n_axes = 2, seed = 1, cpus = 1)
           else if (id == 'oplsda') count_ord(m, id, rep(c('a', 'b'), 15), n_orthogonal = 1, cpus = 1)
           else if (id %in% c('rpca_phylo', 'dpcoa')) count_ord(m, id, tree = tr, n_axes = 2, seed = 1, cpus = 1)
           else count_ord(m, id, n_axes = 2, seed = 1, cpus = 1)
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
  expect_identical(mo('rpca_phylo'),          'rpca_phylo')
  expect_identical(mo('rpca phylo'),          'rpca_phylo')
  expect_identical(mo('phylogenetic RPCA'),   'rpca_phylo')
  expect_identical(mo('phylo-RPCA'),          'rpca_phylo')
  expect_identical(mo('Phylogenetic DEICODE'), 'rpca_phylo')
  expect_identical(mo('Phylogenetic Robust Aitchison PCA'), 'rpca_phylo')
  expect_identical(mo('phylo'),               'rpca_phylo')
  expect_identical(mo('dpcoa'),               'dpcoa')
  expect_identical(mo('DPCoA'),               'dpcoa')
  expect_identical(mo('double principal coordinate analysis'), 'dpcoa')
  expect_identical(mo('Double Principal Coordinates Analysis'), 'dpcoa')
  expect_identical(mo('double'),              'dpcoa')
  expect_identical(mo('RDA'),                 'rda')
  expect_identical(mo('redundancy analysis'), 'rda')
  expect_identical(mo('tb-RDA'),              'rda')
  expect_identical(mo('transformation-based RDA'), 'rda')
  expect_identical(mo('CCA'),                 'cca')
  expect_identical(mo('canonical correspondence analysis'), 'cca')
  expect_identical(mo('constrained correspondence analysis'), 'cca')

  # 'ca' is a prefix of cca's 'canonical...' alias, but is ca's exact ID.
  expect_identical(mo('ca'),                  'ca')
  expect_identical(mo('CA'),                  'ca')
  expect_identical(mo('C.A.'),                'ca')
  expect_identical(mo('correspondence analysis'), 'ca')
  expect_identical(mo('Correspondence'),      'ca')
  expect_identical(mo('canonical'),           'cca')
  expect_identical(mo('cc'),                  'cca')
  expect_error(match_count_ord('c'), 'matches multiple options: ca, cca')
  expect_identical(mo('PLS-DA'),              'plsda')
  expect_identical(mo('pls'),                 'plsda')
  expect_identical(mo('partial least squares DA'), 'plsda')
  expect_identical(mo('Partial Least Squares Discriminant Analysis'), 'plsda')
  expect_identical(mo('PLS discriminant analysis'), 'plsda')
  expect_error(match_count_ord('p'), 'matches multiple options: pca, plsda')
  expect_identical(mo('OPLS-DA'),             'oplsda')
  expect_identical(mo('opls'),                'oplsda')
  expect_identical(mo('orthogonal PLS-DA'),   'oplsda')
  expect_identical(mo('orthogonal partial least squares DA'), 'oplsda')
  expect_identical(mo('Orthogonal Partial Least Squares Discriminant Analysis'), 'oplsda')
  expect_identical(mo('LDA'),                 'lda')
  expect_identical(mo('linear discriminant'), 'lda')
  expect_identical(mo('l'),                   'lda')
  expect_identical(count_ord(m, 'plsda', grp, cpus = 1), plsda(m, grp, cpus = 1))
  expect_identical(count_ord(m, 'pls-da', groups = grp, n_axes = 3, cpus = 1), plsda(m, grp, n_axes = 3, cpus = 1))

  info <- match_count_ord('pca')
  expect_identical(info$func, pca)
  expect_identical(info$params, names(formals(pca)))
  expect_identical(info$name, 'Principal Component Analysis')
  expect_null(info$alt_ids)


  # Listing. ====

  df <- list_count_ords()
  expect_true(is.data.frame(df))
  expect_identical(names(df), c('name', 'id'))
  expect_identical(df$id, c('ca', 'cca', 'dpcoa', 'lda', 'oplsda', 'pca', 'plsda', 'rda', 'rpca', 'rpca_phylo'))
  expect_identical(rownames(list_count_ords(nm = 'id')), c('ca', 'cca', 'dpcoa', 'lda', 'oplsda', 'pca', 'plsda', 'rda', 'rpca', 'rpca_phylo'))
  expect_identical(list_count_ords(val = 'id', nm = 'name'),
    c('Correspondence Analysis' = 'ca', 'Canonical Correspondence Analysis' = 'cca', 'Double Principal Coordinate Analysis' = 'dpcoa', 'Linear Discriminant Analysis' = 'lda', 'Orthogonal Partial Least Squares Discriminant Analysis' = 'oplsda', 'Principal Component Analysis' = 'pca', 'Partial Least Squares Discriminant Analysis' = 'plsda', 'Redundancy Analysis' = 'rda', 'Robust Aitchison PCA' = 'rpca', 'Phylogenetic Robust Aitchison PCA' = 'rpca_phylo'))
  expect_identical(list_count_ords(val = 'func'), list(ca = ca, cca = cca, dpcoa = dpcoa, lda = lda, oplsda = oplsda, pca = pca, plsda = plsda, rda = rda, rpca = rpca, rpca_phylo = rpca_phylo))
  expect_identical(list_count_ords(val = 'list')$pca, match_count_ord('pca'))
  expect_identical(list_count_ords('name'), c('Correspondence Analysis', 'Canonical Correspondence Analysis', 'Double Principal Coordinate Analysis', 'Linear Discriminant Analysis', 'Orthogonal Partial Least Squares Discriminant Analysis', 'Principal Component Analysis', 'Partial Least Squares Discriminant Analysis', 'Redundancy Analysis', 'Robust Aitchison PCA', 'Phylogenetic Robust Aitchison PCA'))


  # Invalid input, and methods of the other wrappers. ====

  expect_error(count_ord(m, 'xyz'),   'does not match any option: ca, cca, dpcoa, lda, oplsda, pca, plsda, rda, rpca')
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
