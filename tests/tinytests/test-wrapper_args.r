#test_that("wrapper arguments don't shadow the methods' arguments", {

  # dist_ord(), count_ord(), dist_test(), replace_zeros(), da_test() and
  # feature_assoc() pass `...` on to
  # the function of
  # the method they run. R matches named arguments to the wrapper's own
  # arguments before `...`, exactly or by prefix, so a method argument that
  # equals or begins a wrapper argument would never reach the method: e.g.
  # `dist_ord(dm, 'capda', m = 3)` would set dist_ord()'s `method`. Every
  # method must therefore avoid its wrapper's argument names and their
  # prefixes, except those the wrapper passes through on purpose.
  #
  # Returns the method arguments that break the rule, as "method: arg".
  shadowed <- function (wrapper, funcs, passed = character(0)) {
    wargs <- setdiff(names(formals(wrapper)), '...')
    bad   <- character(0)
    for (id in names(funcs)) {
      args <- names(formals(funcs[[id]]))
      args <- setdiff(args[-1], c('...', passed))  # the first is the data, given by position
      hit  <- args[vapply(args, function (a) any(startsWith(wargs, a)), logical(1))]
      if (length(hit)) bad <- c(bad, paste0(id, ': ', hit))
    }
    bad
  }

  expect_identical(shadowed(dist_ord,  list_dist_ords('func')),  character(0))
  expect_identical(shadowed(count_ord, list_count_ords('func')), character(0))
  expect_identical(shadowed(dist_test, list_dist_tests('func'), passed = c('groups', 'data')), character(0))
  expect_identical(shadowed(replace_zeros, list_zero_replacements('func')), character(0))
  expect_identical(shadowed(da_test, list_da_tests('func'), passed = c('groups', 'data', 'formula')), character(0))
  expect_identical(shadowed(feature_assoc, list_feature_assocs('func')), character(0))

  # da_sensitivity() passes `...` to da_test() and takes its grids after
  # `...`, where they match only by exact name: no method may take an
  # argument named exactly `pseudocounts`, `zeros`, `scale_sds` or
  # `crossed`, or da_sensitivity() would capture it.
  grids <- setdiff(names(formals(da_sensitivity)), names(formals(da_test)))
  expect_identical(grids, c('pseudocounts', 'zeros', 'scale_sds', 'crossed'))
  expect_identical(names(formals(da_sensitivity))[1:3], c('counts', 'method', '...'))
  taken <- unlist(lapply(list_da_tests('func'), function (f) intersect(names(formals(f)), grids)))
  expect_identical(length(taken), 0L)

  # The wrappers' own arguments, which the rule protects. Changing them
  # changes what methods may take.
  expect_identical(names(formals(dist_ord)),  c('dm', 'method', '...'))
  expect_identical(names(formals(count_ord)), c('counts', 'method', '...'))
  expect_identical(names(formals(dist_test)), c('dm', 'test', 'groups', 'data', '...'))
  expect_identical(names(formals(replace_zeros)), c('counts', 'method', '...'))
  expect_identical(names(formals(da_test)), c('counts', 'method', 'groups', 'data', 'formula', '...'))
  expect_identical(names(formals(feature_assoc)), c('counts', 'method', '...'))

  # The check itself catches exact names and prefixes, and allows the
  # pass-through ones.
  f <- list(
    good   = function (dm, k, seed) NULL,
    prefix = function (dm, m, k) NULL,
    exact  = function (dm, method) NULL,
    data   = function (dm, groups, data, te) NULL )
  expect_identical(shadowed(dist_ord, f), c('prefix: m', 'exact: method'))
  expect_identical(shadowed(replace_zeros, list(z = function (counts, margin, co, delta) NULL)), 'z: co')
  expect_identical(shadowed(dist_test, f['data'], passed = c('groups', 'data')), 'data: te')
  expect_identical(
    shadowed(da_test, list(z = function (counts, groups, data, formula, form, me, g) NULL),
             passed = c('groups', 'data', 'formula')),
    c('z: form', 'z: me', 'z: g') )

  # Arguments that do pass through, by name, reach the method.
  dm <- bray(iris[1:30, 1:4])
  expect_identical(dist_ord(dm, 'pcoa', n_axes = 2, cpus = 1), pcoa(dm, n_axes = 2, cpus = 1))
  expect_identical(
    dist_test(dm, 'mantel', iris$Petal.Width[1:30], method = 'spearman', permutations = 9, cpus = 1),
    mantel(dm, iris$Petal.Width[1:30], method = 'spearman', permutations = 9, cpus = 1) )
  x <- matrix(c(0, 1, 2, 3, 0, 5), 2)
  expect_identical(replace_zeros(x, 'mult', delta = 0.01, margin = 2L),
                   replace_zeros_mult(x, delta = 0.01, margin = 2L))
  x <- matrix(1:24 + 0, 6, dimnames = list(NULL, paste0('f', 1:4)))
  expect_identical(da_test(x, 'da_ttest', rep(1:2, 3), conf_level = 0.9, p_adjust = 'BH', cpus = 1),
                   da_ttest(x, rep(1:2, 3), conf_level = 0.9, p_adjust = 'BH', cpus = 1))

#})
