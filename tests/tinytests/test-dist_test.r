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
      dist_test(d, 'trt', meta, method = m, permutations = 19, strata = 'block', cpus = 1),
      f(d, 'trt', meta, permutations = 19, strata = 'block', cpus = 1),
      info = m )

    # A vector, unnamed or named.
    expect_identical(
      dist_test(d, meta$trt, method = m, permutations = 19, seed = 3, cpus = 1),
      f(d, meta$trt, permutations = 19, seed = 3, cpus = 1),
      info = m )
    expect_identical(
      dist_test(d, setNames(meta$trt, rownames(meta)), method = m, permutations = 9, cpus = 1),
      f(d, setNames(meta$trt, rownames(meta)), permutations = 9, cpus = 1),
      info = m )

    # Errors come from the method.
    expect_error(dist_test(d, 'nothere', meta, method = m), "'nothere' is not one", info = m)
  }

  # A method's own arguments.
  expect_identical(
    dist_test(d, 'trt', meta, method = 'mrpp', weights = 'n-1', permutations = 9, cpus = 1),
    mrpp(d, 'trt', meta, weights = 'n-1', permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'trt', meta, method = 'permanova', by = 'onedf', permutations = 9, cpus = 1),
    permanova(d, 'trt', meta, by = 'onedf', permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'trt', meta, method = 'permdisp', type = 'median', permutations = 9, cpus = 1),
    permdisp(d, 'trt', meta, type = 'median', permutations = 9, cpus = 1) )

  # permanova()'s formula, without `groups`.
  expect_identical(
    dist_test(d, data = meta, method = 'permanova', formula = ~ trt + ph, permutations = 9, cpus = 1),
    permanova(d, data = meta, formula = ~ trt + ph, permutations = 9, cpus = 1) )
  expect_error(dist_test(d, data = meta, method = 'mrpp'), '`groups` is required')

  # A symmetric matrix, as for the methods.
  expect_identical(
    dist_test(as.matrix(d), 'trt', meta, method = 'anosim', permutations = 9, cpus = 1),
    anosim(d, 'trt', meta, permutations = 9, cpus = 1) )

  # Results bind into one table.
  tbl <- do.call(rbind, lapply(list_dist_tests(val = 'id', groups = 'categorical'), function (m)
    as.data.frame(dist_test(d, 'trt', meta, m, permutations = 9, cpus = 1)) ))
  expect_identical(tbl$test, c('anosim', 'mrpp', 'permanova', 'permdisp'))
  expect_identical(names(tbl), c('test', 'stat_name', 'statistic', 'p_value', 'z', 'n_perm'))



  # The kinds of variable each method tests. ====

  # Numeric values: permanova only.
  expect_identical(
    dist_test(d, 'ph', meta, method = 'permanova', permutations = 9, cpus = 1),
    permanova(d, 'ph', meta, permutations = 9, cpus = 1) )
  for (m in c('anosim', 'mrpp', 'permdisp')) {
    expect_error(dist_test(d, 'ph', meta, method = m), "column 'ph' is numeric", info = m)
    expect_error(dist_test(d, meta$ph, method = m), '`groups` is numeric', info = m)
    expect_error(dist_test(d, meta$ph, method = m), 'use: mantel, permanova', info = m)
  }

  # mantel(): numeric values only, passed as its `y`.
  expect_identical(
    dist_test(d, 'ph', meta, method = 'mantel', permutations = 19, strata = 'block', cpus = 1),
    mantel(d, 'ph', meta, permutations = 19, strata = 'block', cpus = 1) )
  expect_identical(
    dist_test(d, meta$ph, method = 'mantel', permutations = 19, seed = 3, cpus = 1),
    mantel(d, meta$ph, permutations = 19, seed = 3, cpus = 1) )
  expect_identical(
    dist_test(d, setNames(meta$ph, rownames(meta)), method = 'mantel', permutations = 9, cpus = 1),
    mantel(d, setNames(meta$ph, rownames(meta)), permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, 'depth', meta, method = 'Mantel test', alternative = 'two.sided',
              permutations = 9, cpus = 1),
    mantel(d, 'depth', meta, alternative = 'two.sided', permutations = 9, cpus = 1) )
  expect_error(dist_test(d, 'trt', meta, method = 'mantel'), "column 'trt' is categorical")
  expect_error(dist_test(d, meta$trt, method = 'mantel'), '`groups` is categorical')
  expect_error(
    dist_test(d, meta$trt, method = 'mantel'),
    'use: anosim, mrpp, permanova, permdisp' )
  expect_error(dist_test(d, 'nothere', meta, method = 'mantel'))

  # Character, logical and factor groups are categorical.
  expect_identical(
    dist_test(d, as.character(meta$trt), method = 'mrpp', permutations = 9, cpus = 1),
    mrpp(d, as.character(meta$trt), permutations = 9, cpus = 1) )
  expect_identical(
    dist_test(d, meta$site == 'A', method = 'anosim', permutations = 9, cpus = 1),
    anosim(d, meta$site == 'A', permutations = 9, cpus = 1) )

  # Neither.
  dates <- as.Date('2020-01-01') + seq_len(n)
  expect_error(dist_test(d, dates, method = 'permanova'), 'of class Date')

  # A vector parallel to the rows of `data`.
  expect_identical(
    dist_test(d, meta$trt, meta, method = 'mrpp', permutations = 9, cpus = 1),
    mrpp(d, meta$trt, meta, permutations = 9, cpus = 1) )
  expect_error(dist_test(d, meta$ph, meta, method = 'mrpp'), '`groups` is numeric')



  # The registry. ====

  df <- list_dist_tests()
  expect_identical(names(df), c('name', 'id', 'groups'))
  expect_identical(df$id, c('anosim', 'mantel', 'mrpp', 'permanova', 'permdisp'))
  expect_identical(list_dist_tests(val = 'id', groups = 'numeric'), c('mantel', 'permanova'))
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
  expect_error(match_dist_test('t-test'), 'does not match any option')
  expect_identical(match_dist_test('adonis', groups = 'numeric')$id, 'permanova')
  expect_error(match_dist_test('mrpp', groups = 'numeric'), 'does not match any option')
  expect_error(match_dist_test(1),        'character')
  expect_error(match_dist_test(c('a', 'b')), 'length 1')
  expect_error(match_dist_test(NA_character_), 'NA')
  expect_error(match_dist_test('!!'),     'alphabetic')

  expect_identical(
    dist_test(d, 'trt', meta, method = 'betadisper', permutations = 9, cpus = 1),
    permdisp(d, 'trt', meta, permutations = 9, cpus = 1) )

#})
