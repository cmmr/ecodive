

test_that("list_methods", {
  
  df <- expect_silent(list_methods(nm = 'id'))
  expect_true(inherits(df, 'data.frame'))
  expect_identical(rownames(df), df$id)
  
  
  lst <- expect_silent(list_methods(val = 'list'))
  expect_true(inherits(lst, 'list'))
  expect_true(inherits(lst$faith, 'list'))
  expect_identical(lst$faith$func, faith)
  
  
  funcs <- expect_silent(list_methods(val = 'func'))
  expect_true(inherits(funcs, 'list'))
  expect_identical(funcs$faith, faith)
  
  
  vec <- expect_silent(list_methods(val = 'name', nm = 'id'))
  expect_true(inherits(vec, 'character'))
  expect_identical(names(vec), df$id)
  expect_identical(unname(vec), df$name)
  
  expect_error(list_methods(phylo = 4))
  
})


test_that("match_method", {
  
  expect_error(match_method('badoption'))   # not a valid method name
  expect_error(match_method('s'))           # matches shannon and simpson
  expect_error(match_method(3))             # stopifnot(is.character(method))
  expect_error(match_method(c('a', 'b')))   # stopifnot(length(method) == 1)
  expect_error(match_method(NA_character_)) # stopifnot(!is.na(method))
  expect_error(match_method('  '))          # stopifnot(nchar(method) > 0)
  
  m <- expect_silent(match_method('bray'))
  expect_identical(m$name,         "Bray-Curtis Dissimilarity")
  expect_identical(m$id,           "bray")
  expect_identical(m$div,          "beta")
  expect_identical(m$phylo,        FALSE)
  expect_identical(m$weighted,     TRUE)
  expect_identical(m$int_only,     FALSE)
  expect_identical(m$true_metric,  FALSE)
  expect_identical(m$func,         bray)
  expect_identical(m$params, names(formals(bray)))
  
})
