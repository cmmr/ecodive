#test_that("wrapper arguments don't shadow the methods' arguments", {

  # dist_ord(), count_ord() and dist_test() pass `...` on to the function of
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

  # The wrappers' own arguments, which the rule protects. Changing them
  # changes what methods may take.
  expect_identical(names(formals(dist_ord)),  c('dm', 'method', '...'))
  expect_identical(names(formals(count_ord)), c('counts', 'method', '...'))
  expect_identical(names(formals(dist_test)), c('dm', 'test', 'groups', 'data', '...'))

  # The check itself catches exact names and prefixes, and allows the
  # pass-through ones.
  f <- list(
    good   = function (dm, k, seed) NULL,
    prefix = function (dm, m, k) NULL,
    exact  = function (dm, method) NULL,
    data   = function (dm, groups, data, te) NULL )
  expect_identical(shadowed(dist_ord, f), c('prefix: m', 'exact: method'))
  expect_identical(shadowed(dist_test, f['data'], passed = c('groups', 'data')), 'data: te')

  # Arguments that do pass through, by name, reach the method.
  dm <- bray(iris[1:30, 1:4])
  expect_identical(dist_ord(dm, 'pcoa', n_axes = 2, cpus = 1), pcoa(dm, n_axes = 2, cpus = 1))
  expect_identical(
    dist_test(dm, 'mantel', iris$Petal.Width[1:30], method = 'spearman', permutations = 9, cpus = 1),
    mantel(dm, iris$Petal.Width[1:30], method = 'spearman', permutations = 9, cpus = 1) )

#})
