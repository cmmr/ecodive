#test_that("dist_test", {

  # Inputs only: no reference outputs, so every test here runs on CRAN.
  ex   <- readRDS('expected/permanova.rds')
  d    <- bray(ex$counts)
  meta <- ex$meta
  n    <- attr(d, 'Size')



  # The same as calling each method directly. ====

  for (m in c('anosim', 'mrpp', 'permanova', 'permdisp')) {

    f <- get(m, asNamespace('ecodive'))

    # A column of `data`, with strata.
    expect_identical(
      dist_test(d, m, 'trt', meta, permutations = 19, strata = 'block', cpus = 1),
      f(d, 'trt', meta, permutations = 19, strata = 'block', cpus = 1),
      info = m )

    # A vector, unnamed or named.
    expect_identical(
      dist_test(d, m, meta$trt, permutations = 19, seed = 3, cpus = 1),
      f(d, meta$trt, permutations = 19, seed = 3, cpus = 1),
      info = m )
    expect_identical(
      dist_test(d, m, setNames(meta$trt, rownames(meta)), permutations = 9, cpus = 1),
      f(d, setNames(meta$trt, rownames(meta)), permutations = 9, cpus = 1),
      info = m )

    # Errors come from the method.
    expect_error(dist_test(d, m, 'nothere', meta), "'nothere' is not one", info = m)
  }

  # A method's own arguments.
  expect_identical(
    dist_test(d, 'mrpp', 'trt', meta, weights = 'n-1', permutations = 9, cpus = 1),
    mrpp(d, 'trt', meta, weights = 'n-1', permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'permanova', 'trt', meta, by = 'onedf', permutations = 9, cpus = 1),
    permanova(d, 'trt', meta, by = 'onedf', permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'permdisp', 'trt', meta, type = 'median', permutations = 9, cpus = 1),
    permdisp(d, 'trt', meta, type = 'median', permutations = 9, cpus = 1) )

  # permanova()'s formula, without `groups`.
  expect_identical(
    dist_test(d, 'permanova', data = meta, formula = ~ trt + ph, permutations = 9, cpus = 1),
    permanova(d, data = meta, formula = ~ trt + ph, permutations = 9, cpus = 1) )
  expect_error(dist_test(d, 'mrpp', data = meta), '`groups` is required')

  # A symmetric matrix, as for the methods.
  expect_identical(
    no_display(dist_test(as.matrix(d), 'anosim', 'trt', meta, permutations = 9, cpus = 1)),
    no_display(anosim(d, 'trt', meta, permutations = 9, cpus = 1)) )

  # Results bind into one table.
  tbl <- do.call(rbind, lapply(list_dist_tests(val = 'id', groups = 'categorical'), function (m)
    as.data.frame(dist_test(d, m, 'trt', meta, permutations = 9, cpus = 1)) ))
  expect_identical(tbl$test, c('anosim', 'mrpp', 'permanova', 'permdisp'))
  expect_identical(names(tbl), c('test', 'stat_name', 'statistic', 'p_value', 'z', 'n_perm'))



  # The kinds of variable each method tests. ====

  # Numeric values: permanova only.
  expect_identical(
    dist_test(d, 'permanova', 'ph', meta, permutations = 9, cpus = 1),
    permanova(d, 'ph', meta, permutations = 9, cpus = 1) )
  for (m in c('anosim', 'mrpp', 'permdisp')) {
    expect_error(dist_test(d, m, 'ph', meta), "column 'ph' is numeric", info = m)
    expect_error(dist_test(d, m, meta$ph), '`groups` is numeric', info = m)
    expect_error(dist_test(d, m, meta$ph), 'use: mantel, mantel_partial, permanova', info = m)
  }

  # mantel(): numeric values only, passed as its `y`.
  expect_identical(
    dist_test(d, 'mantel', 'ph', meta, permutations = 19, strata = 'block', cpus = 1),
    mantel(d, 'ph', meta, permutations = 19, strata = 'block', cpus = 1) )
  expect_identical(
    dist_test(d, 'mantel', meta$ph, permutations = 19, seed = 3, cpus = 1),
    mantel(d, meta$ph, permutations = 19, seed = 3, cpus = 1) )
  expect_identical(
    dist_test(d, 'mantel', setNames(meta$ph, rownames(meta)), permutations = 9, cpus = 1),
    mantel(d, setNames(meta$ph, rownames(meta)), permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'Mantel test', 'depth', meta, alternative = 'two.sided', permutations = 9, cpus = 1),
    mantel(d, 'depth', meta, alternative = 'two.sided', permutations = 9, cpus = 1) )
  # mantel()'s own `method` (the correlation) passes through.
  sp <- dist_test(d, 'mantel', 'ph', meta, method = 'spearman', permutations = 9, cpus = 1)
  expect_identical(sp, mantel(d, 'ph', meta, method = 'spearman', permutations = 9, cpus = 1))
  expect_identical(sp$method, 'spearman')
  expect_error(dist_test(d, 'mantel', 'trt', meta), "column 'trt' is categorical")
  expect_error(dist_test(d, 'mantel', meta$trt), '`groups` is categorical')
  expect_error(
    dist_test(d, 'mantel', meta$trt),
    'use: anosim, mrpp, permanova, permdisp' )
  expect_error(dist_test(d, 'mantel', 'nothere', meta))

  # Character, logical and factor groups are categorical.
  expect_identical(
    dist_test(d, 'mrpp', as.character(meta$trt), permutations = 9, cpus = 1),
    mrpp(d, as.character(meta$trt), permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'anosim', meta$site == 'A', permutations = 9, cpus = 1),
    anosim(d, meta$site == 'A', permutations = 9, cpus = 1) )

  # Neither.
  dates <- as.Date('2020-01-01') + seq_len(n)
  expect_error(dist_test(d, 'permanova', dates), 'of class Date')

  # A vector parallel to the rows of `data`.
  expect_identical(
    dist_test(d, 'mrpp', meta$trt, meta, permutations = 9, cpus = 1),
    mrpp(d, meta$trt, meta, permutations = 9, cpus = 1) )
  expect_error(dist_test(d, 'mrpp', meta$ph, meta), '`groups` is numeric')



  # The registry. ====

  df <- list_dist_tests()
  expect_identical(names(df), c('name', 'id', 'groups'))
  expect_identical(df$id, c('anosim', 'mantel', 'mantel_partial', 'mrpp', 'permanova', 'permdisp'))
  expect_identical(list_dist_tests(val = 'id', groups = 'numeric'), c('mantel', 'mantel_partial', 'permanova'))
  expect_identical(
    list_dist_tests(val = 'id', groups = 'categorical'),
    c('anosim', 'mrpp', 'permanova', 'permdisp') )
  expect_identical(names(list_dist_tests(val = 'func')), df$id)
  expect_identical(list_dist_tests(val = 'func')$mrpp, mrpp)
  expect_identical(names(list_dist_tests(val = 'name', nm = 'id')), df$id)
  expect_identical(rownames(list_dist_tests(nm = 'id')), df$id)
  expect_error(list_dist_tests(groups = 'ordinal'), 'Invalid value for `groups`')

  info <- list_dist_tests(val = 'list')$permanova
  expect_identical(info$id, 'permanova')
  expect_identical(info$func, permanova)
  expect_identical(info$params, names(formals(permanova)))
  expect_identical(info$groups, 'categorical numeric')
  expect_identical(match_dist_test('permanova'), info)

  # Forgiving matching, and other packages' names.
  expect_identical(match_dist_test('PERMANOVA')$id,        'permanova')
  expect_identical(match_dist_test('adonis')$id,           'permanova')
  expect_identical(match_dist_test('adonis2')$id,          'permanova')
  expect_identical(match_dist_test('perma')$id,            'permanova')
  expect_identical(match_dist_test('betadisper')$id,       'permdisp')
  expect_identical(match_dist_test('PERMDISP2')$id,        'permdisp')
  expect_identical(match_dist_test('Analysis of Similarities')$id, 'anosim')
  expect_identical(match_dist_test('multi-response')$id,   'mrpp')
  expect_identical(match_dist_test('MRPP')$id,             'mrpp')
  expect_identical(match_dist_test('Mantel')$id,           'mantel')
  expect_identical(match_dist_test('mantel test')$id,      'mantel')
  expect_identical(match_dist_test('mantel')$params, names(formals(mantel)))
  expect_error(match_dist_test('mantel', groups = 'categorical'), 'does not match any option')
  expect_error(match_dist_test('perm'),   'matches multiple options')
  expect_error(match_dist_test('t-test'), 'use da_test')
  expect_error(match_dist_test('zzz'), 'does not match any option')
  expect_identical(match_dist_test('adonis', groups = 'numeric')$id, 'permanova')
  expect_error(match_dist_test('mrpp', groups = 'numeric'), 'does not match any option')
  expect_error(match_dist_test(1),        'character')
  expect_error(match_dist_test(c('a', 'b')), 'length 1')
  expect_error(match_dist_test(NA_character_), 'NA')
  expect_error(match_dist_test('!!'),     'alphabetic')

  expect_identical(
    dist_test(d, 'betadisper', 'trt', meta, permutations = 9, cpus = 1),
    permdisp(d, 'trt', meta, permutations = 9, cpus = 1) )



  # A factor's NA level (from addNA()) is a group, or a stratum, like any
  # other, as in model.matrix(). It crashed R. ====

  na_g  <- addNA(factor(rep(c('a', NA, 'b'), length.out = 12)))
  zz_g  <- factor(rep(c('a', 'zz', 'b'), length.out = 12), levels = c('a', 'b', 'zz'))
  na_s  <- addNA(factor(rep(c(1, 2, NA), each = 4)))
  zz_s  <- factor(rep(c(1, 2, 99), each = 4))
  d12   <- dist(c(1, 4, 2, 8, 5, 7, 3, 9, 6, 2, 8, 1))
  ab    <- rep(c('a', 'b'), 6)
  for (m in c('permanova', 'mrpp', 'anosim', 'permdisp')) {
    expect_identical(
      dist_test(d12, m, na_g, permutations = 19, cpus = 1)[1:4],
      dist_test(d12, m, zz_g, permutations = 19, cpus = 1)[1:4], info = m )
    expect_identical(
      dist_test(d12, m, ab, strata = na_s, permutations = 19, cpus = 1)[1:4],
      dist_test(d12, m, ab, strata = zz_s, permutations = 19, cpus = 1)[1:4], info = m )
  }
  x12 <- c(5, 1, 4, 2, 8, 3, 9, 6, 7, 2, 1, 4)
  expect_identical(
    mantel(d12, x12, strata = na_s, permutations = 19, cpus = 1)[1:4],
    mantel(d12, x12, strata = zz_s, permutations = 19, cpus = 1)[1:4] )
  expect_identical(
    permanova(d12, data = data.frame(g = na_g), formula = ~ g, permutations = 19, cpus = 1)[1:4],
    permanova(d12, data = data.frame(g = zz_g), formula = ~ g, permutations = 19, cpus = 1)[1:4] )
  ord12 <- pcoa(d12, n_axes = 1, cpus = 1)
  e_na  <- envfit(ord12, data.frame(g = na_g, row.names = NULL), strata = na_s, permutations = 19, cpus = 1)
  e_zz  <- envfit(ord12, data.frame(g = zz_g, row.names = NULL), strata = zz_s, permutations = 19, cpus = 1)
  expect_identical(e_na$r2, e_zz$r2)
  expect_identical(e_na$p_value, e_zz$p_value)
  expect_identical(nrow(e_na), 3L)  # a centroid for each level, NA's too

#})
