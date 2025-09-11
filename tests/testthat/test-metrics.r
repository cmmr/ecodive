

test_that("list_metrics", {
  
  df <- expect_silent(list_metrics(nm = 'id'))
  expect_true(inherits(df, 'data.frame'))
  expect_identical(rownames(df), df$id)
  
  
  lst <- expect_silent(list_metrics(val = 'list'))
  expect_true(inherits(lst, 'list'))
  expect_true(inherits(lst$faith, 'list'))
  expect_identical(lst$faith$func, faith)
  
  
  funcs <- expect_silent(list_metrics(val = 'func'))
  expect_true(inherits(funcs, 'list'))
  expect_identical(funcs$faith, faith)
  
  
  vec <- expect_silent(list_metrics(val = 'name', nm = 'id'))
  expect_true(inherits(vec, 'character'))
  expect_identical(names(vec), df$id)
  expect_identical(unname(vec), df$name)
  
  expect_error(list_metrics(phylo = 4))
  
})


test_that("match_metric", {
  
  expect_error(match_metric('badoption'))   # not a valid metric name
  expect_error(match_metric('s'))           # matches shannon and simpson
  expect_error(match_metric(3))             # stopifnot(is.character(metric))
  expect_error(match_metric(c('a', 'b')))   # stopifnot(length(metric) == 1)
  expect_error(match_metric(NA_character_)) # stopifnot(!is.na(metric))
  expect_error(match_metric('  '))          # stopifnot(nchar(metric) > 0)
  
  m <- expect_silent(match_metric('bray'))
  expect_identical(m$name,         "Bray-Curtis Dissimilarity")
  expect_identical(m$id,           "bray")
  expect_identical(m$div,          "beta")
  expect_identical(m$phylo,        FALSE)
  expect_identical(m$weighted,     TRUE)
  expect_identical(m$true_metric,  FALSE)
  expect_identical(m$func,         bray)
  expect_identical(m$params, names(formals(bray)))
  
})
