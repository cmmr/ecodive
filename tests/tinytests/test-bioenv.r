#test_that("bioenv", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L


  # A direct computation, with base R: every subset's distances (as
  # vegan's documentation describes them) and their correlation with dm.
  bioenv_ref <- function (dm, env, method = 'spearman', metric = 'euclidean',
                          max_vars = ncol(env), given = NULL) {
    x <- as.matrix(env)
    z <- if (metric %in% c('euclidean', 'mahalanobis')) scale(x)
         else apply(x, 2, function (v) (v - min(v)) / (max(v) - min(v)))
    dfun <- function (s) {
      zs <- z[, s, drop = FALSE]
      switch(metric,
        euclidean   = dist(zs),
        manhattan   = dist(zs, 'manhattan'),
        gower       = dist(zs, 'manhattan') / length(s),
        mahalanobis = dist(zs %*% solve(chol(stats::cov(zs)))) )
    }
    rc <- function (a, b) stats::cor(as.vector(a), as.vector(b), method = method)
    pc <- function (d) {
      if (is.null(given)) return (rc(dm, d))
      rxy <- rc(dm, d); rxz <- rc(dm, given); ryz <- rc(d, given)
      (rxy - rxz * ryz) / sqrt(1 - rxz^2) / sqrt(1 - ryz^2)
    }
    rows <- lapply(seq_len(max_vars), function (k) {
      cmb <- utils::combn(ncol(x), k, simplify = FALSE)
      r   <- vapply(cmb, function (s) pc(dfun(s)), 0)
      i   <- which.max(r)
      data.frame(size = k, variables = paste(colnames(x)[cmb[[i]]], collapse = ', '),
                 correlation = r[[i]], stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }

  plain <- function (b) { attributes(b) <- attributes(b)[c('names', 'row.names')]; class(b) <- 'data.frame'; b }


  # Stored results from vegan::bioenv() and scikit-bio's bioenv(), on
  # varespec and varechem. ====

  if (at_home()) {

    ref <- readRDS(file.path('expected', 'bioenv.rds'))

    for (cs in ref$cases) {
      info <- paste(cs$method, cs$metric, cs$max_vars, cs$partial, !is.null(cs$env))
      b <- bioenv(
        ref$dm, if (is.null(cs$env)) ref$env else cs$env,
        method = cs$method, metric = cs$metric, max_vars = cs$max_vars,
        given = if (cs$partial) ref$given, cpus = max_cpus )
      ex <- cs$expected

      expect_identical(b$size, ex$size, info = info)
      expect_identical(b$variables, ex$variables, info = info)

      # vegan's Mahalanobis whitening rounds tied distances apart (for one
      # variable, vegan has 254 distinct distances where the variable has
      # 252), which moves rank correlations by about 5e-5. bioenv() keeps
      # the ties: for one variable, its ranks are Euclidean's exactly.
      if (cs$metric == 'mahalanobis' && cs$method != 'pearson') {
        expect_equal(b$correlation[-1], ex$correlation[-1], tolerance = 1e-10, info = info)
        expect_equal(b$correlation[1], ex$correlation[1], tolerance = 5e-4, info = info)
        eu <- bioenv(ref$dm, ref$env, method = cs$method, max_vars = 1,
                     given = if (cs$partial) ref$given, cpus = 1)
        expect_identical(b$variables[1], eu$variables[1], info = info)
        expect_equal(b$correlation[1], eu$correlation[1], tolerance = 1e-14, info = info)
      }
      else {
        expect_equal(b$correlation, ex$correlation, tolerance = 1e-10, info = info)
      }
    }

    # scikit-bio, on its defaults.
    b <- bioenv(ref$dm, ref$env, cpus = max_cpus)
    expect_identical(b$size, ref$skbio$size)
    expect_identical(b$variables, ref$skbio$variables)
    expect_equal(b$correlation, ref$skbio$correlation, tolerance = 1e-10)
    expect_identical(attr(b, 'best'), c('N', 'P', 'Al', 'Mn', 'Baresoil'))
    expect_identical(attr(b, 'n_subsets'), 16383)

    # The ecodive argument forms give the same results: a matrix, the given
    # distances as a column of `data` or values per row, and rows in
    # another order.
    env2 <- cbind(ref$env, depth = seq_len(24))
    gd   <- dist(env2$depth)
    exp  <- bioenv(ref$dm, ref$env, given = gd, max_vars = 3, cpus = 1)
    expect_identical(bioenv(ref$dm, env2, given = 'depth', max_vars = 3, cpus = 1), exp)
    expect_identical(bioenv(ref$dm, ref$env, given = env2$depth, max_vars = 3, cpus = 1), exp)
    expect_identical(bioenv(as.matrix(ref$dm), ref$env, given = unname(as.matrix(gd)), max_vars = 3, cpus = 1), exp)
    o <- rev(seq_len(24))
    expect_identical(bioenv(ref$dm, env2[o, ], given = 'depth', max_vars = 3, cpus = 1), exp)
  }



  # A direct computation, for every option, including Kendall with
  # `given` (which vegan refuses). ====

  set.seed(11)
  n   <- 12
  env <- data.frame(
    row.names = paste0('S', 1:n),
    a = rnorm(n), b = runif(n), c = rpois(n, 4) + 0.5, d = rnorm(n, 10, 3) )
  cnt <- matrix(rpois(n * 8, 4), n, 8, dimnames = list(rownames(env), NULL))
  cnt[, 1] <- cnt[, 1] + round(3 * env$a - min(3 * env$a))
  dm  <- bray(cnt, cpus = 1)
  gv  <- dist(rnorm(n))

  if (at_home()) {
    for (method in c('spearman', 'pearson', 'kendall'))
      for (metric in c('euclidean', 'mahalanobis', 'manhattan', 'gower'))
        for (g in list(NULL, gv)) {
          info <- paste(method, metric, !is.null(g))
          b <- bioenv(dm, env, method = method, metric = metric, given = g, cpus = 2)
          r <- bioenv_ref(dm, env, method, metric, given = g)
          expect_identical(b$variables, r$variables, info = info)
          expect_equal(b$correlation, r$correlation, tolerance = 1e-10, info = info)
          expect_identical(attr(b, 'method'), method)
          expect_identical(attr(b, 'metric'), metric)
          expect_identical(attr(b, 'partial'), !is.null(g))
          expect_identical(attr(b, 'best'),
            strsplit(r$variables[which.max(r$correlation)], ', ')[[1]], info = info)
        }

    # 'gower' is 'manhattan' divided by the number of variables.
    expect_equal(plain(bioenv(dm, env, method = 'pearson', metric = 'gower', cpus = 1)),
                 plain(bioenv(dm, env, method = 'pearson', metric = 'manhattan', cpus = 1)),
                 tolerance = 1e-12)
  }



  # Results are identical for any number of threads. ====

  set.seed(5)
  n2   <- 40
  env2 <- as.data.frame(matrix(rnorm(n2 * 9), n2, 9, dimnames = list(NULL, letters[1:9])))
  dm2  <- dist(matrix(rnorm(n2 * 3), n2, 3) + as.matrix(env2[, 1:3]))
  for (method in c('spearman', 'pearson', 'kendall'))
    for (metric in c('euclidean', 'mahalanobis', 'manhattan')) {
      b1 <- bioenv(dm2, env2, method = method, metric = metric, given = dist(env2$i), cpus = 1)
      for (cpus in 2:max_cpus)
        expect_identical(
          bioenv(dm2, env2, method = method, metric = metric, given = dist(env2$i), cpus = cpus),
          b1, info = paste(method, metric, cpus) )
    }



  # Ties: the first subset found, in the order of the columns. ====

  # (The vegan fixtures above have a tied copy too.) `a` is the best
  # single variable here, and its copy ties with it.
  tied <- cbind(acopy = env$a, env[c('a', 'b', 'c')])
  b  <- bioenv(dm, tied, cpus = 1)
  b2 <- bioenv(dm, tied[c('a', 'b', 'c', 'acopy')], cpus = 1)
  expect_identical(bioenv(dm, env, max_vars = 1, cpus = 1)$variables, 'a')
  expect_identical(b$variables[1], 'acopy')
  expect_identical(b2$variables[1], 'a')
  expect_identical(b$correlation, b2$correlation)



  # Formulas and columns. ====

  b <- bioenv(dm, env, ~ log(c) + a, cpus = 1)
  r <- bioenv_ref(dm, data.frame(`log(c)` = log(env$c), a = env$a, check.names = FALSE))
  expect_identical(b$variables, r$variables)
  expect_equal(b$correlation, r$correlation, tolerance = 1e-12)
  expect_identical(bioenv(dm, env, '~ log(c) + a', cpus = 1), b)
  expect_identical(bioenv(dm, env, x ~ log(c) + a, cpus = 1), b)
  expect_identical(bioenv(dm, env, ~ . - d, cpus = 1), bioenv(dm, env[c('a', 'b', 'c')], cpus = 1))

  # Non-numeric columns are skipped by default, and the given column too.
  mixed <- cbind(env, f = rep(c('x', 'y'), n / 2), l = rep(c(TRUE, FALSE), n / 2),
                 g = as.double(1:n))
  expect_identical(bioenv(dm, mixed, given = 'g', cpus = 1),
                   bioenv(dm, env, given = dist(1:n), cpus = 1))
  expect_identical(bioenv(dm, mixed[c('a', 'b', 'c', 'd', 'f')], cpus = 1), bioenv(dm, env, cpus = 1))

  # max_vars.
  b <- bioenv(dm, env, max_vars = 2, cpus = 1)
  expect_identical(nrow(b), 2L)
  expect_identical(plain(b), plain(bioenv(dm, env, cpus = 1))[1:2, ])
  expect_identical(attr(b, 'n_subsets'), 10)
  expect_identical(bioenv(dm, env, max_vars = 10, cpus = 1), bioenv(dm, env, cpus = 1))

  # Unnamed samples and rows, matched by position.
  expect_identical(plain(bioenv(unname(as.matrix(dm)), unname(as.matrix(env)), cpus = 1))$correlation,
                   bioenv(dm, env, cpus = 1)$correlation)
  # Extra rows in `data` are ignored.
  extra <- rbind(env, data.frame(row.names = 'Z', a = NA, b = 1, c = 1, d = 1))
  expect_identical(bioenv(dm, extra[c(13, 12:1), ], cpus = 1), bioenv(dm, env, cpus = 1))

  # Three samples, the fewest.
  b <- bioenv(dist(1:3), data.frame(x = c(1, 2, 4), y = c(3, 1, 2)), cpus = 1)
  expect_identical(b$size, 1:2)
  expect_true(all(is.finite(b$correlation)))

  # Linearly dependent variables have no Mahalanobis distances.
  dep <- cbind(env, e = env$a + env$b)
  b <- bioenv(dm, dep[c('a', 'b', 'e')], metric = 'mahalanobis', cpus = 1)
  expect_identical(b$variables[3], NA_character_)
  expect_true(is.na(b$correlation[3]))
  expect_false(anyNA(b$correlation[1:2]))
  expect_true(all(is.na(bioenv(dm, dep[c('a', 'b', 'e')], metric = 'mahalanobis', max_vars = 3, cpus = 1)$correlation[3])))

  # Printed as a data frame.
  expect_true(any(grepl('correlation', capture.output(print(bioenv(dm, env, cpus = 1))))))



  # Errors. ====

  expect_error(bioenv(dm), '`data` is required')
  expect_error(bioenv(dm, NULL), '`data` is required')
  expect_error(bioenv(dm, list(a = 1)), '`data` must be a data.frame')
  expect_error(bioenv(dist(1:2), data.frame(a = 1:2)), 'at least 3 samples')
  expect_error(bioenv(dm, env[1:5, ]), 'no row in `data`')
  expect_error(bioenv(dm, env, method = 'nope'), 'should be one of')
  expect_error(bioenv(dm, env, metric = 'nope'), 'should be one of')
  expect_error(bioenv(dm, env, max_vars = 0), '`max_vars` must be a whole number')
  expect_error(bioenv(dm, env, max_vars = 1.5), '`max_vars` must be a whole number')
  expect_error(bioenv(dm, env, max_vars = NA), '`max_vars` must be a whole number')
  expect_error(bioenv(dm, env, cpus = 0), 'cpus')
  expect_error(bioenv(dm, mixed['f']), 'no numeric columns')
  expect_error(bioenv(dm, mixed, ~ a + f), "Not numeric: 'f'")
  expect_error(bioenv(dm, mixed, ~ a + l), "Not numeric: 'l'")
  expect_error(bioenv(dm, env, ~ a + nope), "Not found in `data`: 'nope'")
  expect_error(bioenv(dm, env, ~ a * b), 'cannot have interactions')
  expect_error(bioenv(dm, env, ~ a + Condition(b)), 'Condition')
  expect_error(bioenv(dm, env, ~ 1), 'no variables')
  expect_error(bioenv(dm, env, 'a + b'), '`formula` must be a formula')
  expect_error(bioenv(dm, env, 42), '`formula` must be a formula')
  na <- env; na$b[3] <- NA
  expect_error(bioenv(dm, na), "Missing values \\(NA\\) in 'b'")
  inf <- env; inf$b[3] <- Inf
  expect_error(bioenv(dm, inf), "Infinite values in 'b'")
  cst <- env; cst$b <- 2
  expect_error(bioenv(dm, cst), "cannot be standardized: 'b'")
  expect_error(bioenv(dm, env, given = 'nope'), '`given` is not a column of `data`')
  expect_error(bioenv(dm, mixed, given = 'f'), "column 'f' is not numeric")
  expect_error(bioenv(dm, na, ~ a, given = 'b'), 'Missing values')
  expect_error(bioenv(dm, env, given = letters[1:n]), '`given` must be distances')
  expect_error(bioenv(dm, env, given = 1:5), 'one value per row of `data`')
  expect_error(bioenv(dm, env, given = dist(1:5)), '`given` has 5 samples')
  expect_error(bioenv(dm, env, given = dist(rep(1, n))), '`given` has the same distance')
  expect_error(bioenv(dist(rep(1, n)), env), '`dm` has the same distance')
  expect_error(bioenv(dm, env, given = dm), 'perfectly correlated')
