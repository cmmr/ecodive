#test_that("varpart", {

  # Inputs and reference results, written by data-raw/expected/varpart.r.
  ref <- readRDS(file.path('expected', 'varpart.rds'))
  dat <- ref$data

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case_x <- function (cs) {
    y <- dat[[cs$x$table]]
    switch(
      cs$x$kind,
      hellinger = y,
      chisquare = y,
      bray      = bray(y),
      euclidean = euclidean(y, norm = 'hellinger') )
  }

  case_sets <- function (cs) {
    lapply(cs$sets, function (s) {
      if (!is.null(s$formula)) return (stats::as.formula(s$formula))
      t <- dat[[s$table]]
      if (!is.null(s$cols)) t <- t[, s$cols, drop = FALSE]
      as.data.frame(t)
    })
  }

  fit_case <- function (cs, cpus = 1L, ...) {
    args <- c(
      list(case_x(cs)), case_sets(cs),
      list(data = dat[[cs$env]], cpus = cpus),
      switch(
        cs$x$kind,
        hellinger = list(norm = 'hellinger'),
        chisquare = list(chisquare = TRUE, permutations = cs$permutations),
        list() ),
      list(...) )
    suppressWarnings(do.call(varpart, args))
  }

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }

  # vegan's "[a+d] = X1 | X3" -> c("[ad]", "X1 | X3").
  split_name <- function (s) {
    parts <- strsplit(s, ' = ', fixed = TRUE)
    list(
      fractions = gsub('+', '', vapply(parts, `[[`, '', 1L), fixed = TRUE),
      label     = trimws(vapply(parts, function (p) if (length(p) > 1) p[[2]] else '', '')) )
  }



  # Same as vegan::varpart(): stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      exp <- cs$expected
      lbl <- cs$name
      k   <- length(cs$sets)

      vp <- fit_case(cs)
      expect_inherits(vp, c('ecodive_varpart', 'data.frame'), info = lbl)

      # The Venn fractions, in vegan's order.
      expect_equal(vp$adj_r_squared, unname(exp$indfract), tolerance = 1e-10, info = lbl)
      expect_identical(rownames(vp), paste0('[', letters[seq_along(exp$indfract)], ']'), info = lbl)
      expect_identical(vp$testable, c(rep(TRUE, k), rep(FALSE, nrow(vp) - k)), info = lbl)
      expect_equal(sum(vp$adj_r_squared), 1, tolerance = 1e-12, info = lbl)
      expect_equal(sum(vp$r_squared), 1, tolerance = 1e-12, info = lbl)

      # The sets and their unions.
      un <- attr(vp, 'unions')
      fr <- exp$fract
      nm <- split_name(rownames(fr))
      expect_equal(un$adj_r_squared, fr$Adj.R.square, tolerance = 1e-10, info = lbl)
      expect_equal(un$r_squared, fr$R.square, tolerance = 1e-10, info = lbl)
      expect_identical(un$df, as.integer(fr$Df), info = lbl)
      expect_identical(gsub('+', '', un$fractions, fixed = TRUE), nm$fractions, info = lbl)
      # vegan calls the union of three or four sets "All".
      expect_identical(
        c(rownames(un)[-nrow(un)], if (k == 2) rownames(un)[[nrow(un)]] else 'All'),
        nm$label, info = lbl )

      # Each set given some of the others.
      cn <- attr(vp, 'conditional')
      if (k == 2) {
        expect_null(cn, info = lbl)
      }
      else {
        ex <- split_name(names(exp$contr))
        i  <- match(ex$label, rownames(cn))
        expect_false(anyNA(i), info = lbl)
        expect_identical(nrow(cn), length(exp$contr), info = lbl)
        expect_equal(cn$adj_r_squared[i], unname(exp$contr), tolerance = 1e-10, info = lbl)
        expect_identical(gsub('+', '', cn$fractions[i], fixed = TRUE), ex$fractions, info = lbl)
      }

      # Total inertia.
      expect_equal(attr(vp, 'total'), exp$SS.Y, tolerance = 1e-10, info = lbl)

      # vegan warns of collinear sets, and so does ecodive.
      if (any(grepl('collinearity', exp$warning)))
        expect_warning(do.call(varpart, c(
          list(case_x(cs)), case_sets(cs), list(data = dat[[cs$env]]),
          if (cs$x$kind == 'hellinger') list(norm = 'hellinger') )), 'collinear', info = lbl)
    }
  }



  # Identical for any cpus. ====

  for (nm in c('mite3_hel', 'mite4_bray', 'mite2_chi')) {
    cs <- case(nm)
    a  <- fit_case(cs, cpus = 1L)
    for (cpus in 2:max_cpus)
      expect_identical(fit_case(cs, cpus = cpus), a, info = paste(nm, cpus))
  }

  # Random permutations: the same for any cpus, and for the same seed.
  cs <- case('vare3_chi')
  cs$permutations <- 19L
  a <- fit_case(cs, cpus = 1L, seed = 3L)
  expect_identical(fit_case(cs, cpus = max_cpus, seed = 3L), a)
  expect_false(identical(fit_case(cs, cpus = 1L, seed = 4L)$adj_r_squared, a$adj_r_squared))
  expect_identical(attr(a, 'n_perm'), 19L)



  # Agreement with ecodive's own constrained ordinations. ====

  vare <- dat$varespec
  chem <- dat$varechem

  vp <- varpart(vare, ~ N + P, ~ K + Ca, ~ Al + Fe, data = chem, norm = 'hellinger')

  # A unique fraction is the partial RDA of its set given the others.
  r <- rda(vare, data = chem, formula = ~ N + P + Condition(K + Ca + Al + Fe), norm = 'hellinger')
  expect_equal(vp['[a]', 'adj_r_squared'], attr(r, 'adj_r_squared'), tolerance = 1e-12)
  expect_equal(vp['[a]', 'df'], 2L)

  # A union is the RDA of all its sets.
  r <- rda(vare, data = chem, formula = ~ N + P + K + Ca, norm = 'hellinger')
  un <- attr(vp, 'unions')
  expect_equal(un['X1+X2', 'adj_r_squared'], attr(r, 'adj_r_squared'), tolerance = 1e-12)
  expect_equal(un['X1+X2', 'r_squared'],     attr(r, 'r_squared'),     tolerance = 1e-12)
  expect_equal(attr(vp, 'total'), attr(r, 'trace') * (nrow(vare) - 1), tolerance = 1e-12)

  # A conditional fraction is a partial RDA.
  r <- rda(vare, data = chem, formula = ~ K + Ca + Condition(Al + Fe), norm = 'hellinger')
  expect_equal(attr(vp, 'conditional')['X2 | X3', 'adj_r_squared'], attr(r, 'adj_r_squared'), tolerance = 1e-12)

  # db-RDA of Bray-Curtis distances.
  dm <- bray(vare)
  vp <- varpart(dm, ~ N + P, ~ Al + Fe, data = chem)
  d  <- dbrda(dm, data = chem, formula = ~ Al + Fe + Condition(N + P))
  expect_equal(vp['[b]', 'adj_r_squared'], attr(d, 'adj_r_squared'), tolerance = 1e-12)
  expect_identical(attr(vp, 'method'), 'dbrda')

  # db-RDA of Euclidean distances of the normalized table is RDA.
  a <- varpart(euclidean(vare, norm = 'hellinger'), ~ N + P, ~ Al + Fe, data = chem)
  b <- varpart(vare, ~ N + P, ~ Al + Fe, data = chem, norm = 'hellinger')
  expect_equal(a$adj_r_squared, b$adj_r_squared, tolerance = 1e-12)
  expect_equal(attr(a, 'unions'), attr(b, 'unions'), tolerance = 1e-12)

  # CCA: each union's R-squared is cca()'s.
  vp <- varpart(vare, ~ N + P, ~ Al + Fe, data = chem, chisquare = TRUE, permutations = 9L)
  cc <- cca(vare, data = chem, formula = ~ N + P + Al + Fe, permutations = 9L)
  expect_equal(attr(vp, 'unions')['X1+X2', 'r_squared'], attr(cc, 'r_squared'), tolerance = 1e-12)
  expect_equal(attr(vp, 'total'), attr(cc, 'trace'), tolerance = 1e-12)
  expect_identical(attr(vp, 'method'), 'cca')
  expect_null(attr(vp, 'norm'))

  # Four sets: every union's R-squared is cca()'s, and the fractions add up.
  vp <- varpart(vare, ~ N + P, ~ K + Ca, ~ Al + Fe, ~ Humdepth + pH, data = chem,
                chisquare = TRUE, permutations = 19L)
  un <- attr(vp, 'unions')
  for (u in c('X1', 'X2+X4', 'X1+X3+X4')) {
    f  <- stats::as.formula(paste('~', paste(c('N + P', 'K + Ca', 'Al + Fe', 'Humdepth + pH')[
      as.integer(sub('X', '', strsplit(u, '+', fixed = TRUE)[[1]]))], collapse = ' + ')))
    cc <- cca(vare, data = chem, formula = f, permutations = 0L)
    expect_equal(un[u, 'r_squared'], attr(cc, 'r_squared'), tolerance = 1e-12, info = u)
  }
  expect_equal(sum(vp$adj_r_squared), 1, tolerance = 1e-12)
  expect_equal(attr(vp, 'total'), attr(cc, 'trace'), tolerance = 1e-12)



  # Inputs. ====

  base <- varpart(vare, ~ N + P, ~ pH, data = chem, norm = 'hellinger')

  # Every kind of set gives the same result, and they mix.
  m  <- as.matrix(chem[, c('N', 'P')])
  same <- function (vp) expect_equal(vp$adj_r_squared, base$adj_r_squared, tolerance = 1e-12)
  same(varpart(vare, chem[, c('N', 'P')], chem['pH'], norm = 'hellinger'))
  same(varpart(vare, m, ~ pH, data = chem, norm = 'hellinger'))
  same(varpart(vare, ~ N + P, chem$pH, data = chem, norm = 'hellinger'))
  same(varpart(vare, ~ N + P, stats::setNames(chem$pH, rownames(chem)), data = chem, norm = 'hellinger'))
  same(varpart(vare, unname(m), unname(chem$pH), norm = 'hellinger'))

  # Rows matched by name: shuffled, and with extra rows.
  shuf <- chem[rev(seq_len(nrow(chem))), ]
  same(varpart(vare, ~ N + P, ~ pH, data = shuf, norm = 'hellinger'))
  same(varpart(vare, m[rev(seq_len(nrow(m))), ], shuf['pH'], norm = 'hellinger'))
  extra <- rbind(chem, chem[1:2, ])
  rownames(extra) <- c(rownames(chem), 'x1', 'x2')
  same(varpart(vare, ~ N + P, extra['pH'], data = extra, norm = 'hellinger'))

  # Sparse tables, and samples in columns.
  same(varpart(Matrix::Matrix(vare, sparse = TRUE), ~ N + P, ~ pH, data = chem, norm = 'hellinger'))
  same(varpart(slam::as.simple_triplet_matrix(vare), ~ N + P, ~ pH, data = chem, norm = 'hellinger'))
  same(varpart(t(vare), ~ N + P, ~ pH, data = chem, norm = 'hellinger', margin = 2L))

  # Character and logical columns are coded as factors.
  meta <- data.frame(
    row.names = rownames(chem),
    g = rep(c('a', 'b', 'c'), 8),
    h = rep(c(TRUE, FALSE), 12) )
  a <- varpart(vare, meta['g'], meta['h'], norm = 'hellinger')
  b <- varpart(vare, ~ g, ~ h, data = meta, norm = 'hellinger')
  expect_equal(a$adj_r_squared, b$adj_r_squared, tolerance = 1e-12)
  expect_identical(attr(a, 'unions')$df, c(2L, 1L, 3L))

  # Named sets.
  vp <- varpart(vare, chem = ~ N + P, acidity = ~ pH, data = chem, norm = 'hellinger')
  expect_identical(vp$sets, c('chem', 'acidity', 'chem & acidity', 'residuals'))
  expect_identical(rownames(attr(vp, 'unions')), c('chem', 'acidity', 'chem+acidity'))
  expect_identical(attr(vp, 'set_labels'), c(chem = '~N + P', acidity = '~pH'))
  expect_identical(attr(vp, 'norm'), 'hellinger')
  expect_identical(attr(vp, 'n'), 24L)
  expect_identical(vp$df, c(2L, 1L, NA, NA))

  vp <- varpart(vare, a = ~ N, ~ P, b = ~ K, data = chem)
  expect_identical(vp$sets[1:3], c('a', 'X2', 'b'))
  expect_identical(attr(vp, 'norm'), 'none')

  # Collinear sets: shared columns move into the shared fraction, with a
  # warning; a set nested in another has no unique fraction.
  expect_warning(
    vp <- varpart(vare, ~ N + P + K, ~ N, data = chem, norm = 'hellinger'),
    'collinear' )
  expect_equal(vp['[b]', 'adj_r_squared'], 0, tolerance = 1e-12)
  expect_identical(vp['[b]', 'df'], 0L)

  # No residual degrees of freedom: NA.
  small <- vare[1:5, ]
  vp <- varpart(small, ~ N + P, ~ K + Ca, data = chem[1:5, ], norm = 'hellinger')
  expect_true(is.na(vp['[a]', 'adj_r_squared']))
  expect_false(anyNA(vp['[a]', 'r_squared']))



  # Printing. ====

  for (sets in list(list(~ N, ~ P), list(~ N, ~ P, ~ K), list(~ N, ~ P, ~ K, ~ Ca))) {
    vp  <- do.call(varpart, c(list(vare), sets, list(data = chem, norm = 'hellinger')))
    out <- capture.output(print(vp))
    k   <- length(sets)
    expect_true(any(grepl('Variation partitioning by RDA (hellinger): 24 samples', out, fixed = TRUE)), info = k)
    expect_true(any(grepl(sprintf('[%s]', letters[nrow(vp)]), out, fixed = TRUE)), info = k)
    expect_identical(any(grepl('Residuals', out, fixed = TRUE)), k < 4, info = k)
    expect_identical(any(grepl('+-', out, fixed = TRUE)), k < 4, info = k)
    capture.output(pv <- print(vp))
    expect_identical(pv, vp)
  }
  out <- capture.output(print(varpart(bray(vare), ~ N, ~ P, data = chem)))
  expect_true(any(grepl('by db-RDA', out)))
  out <- capture.output(print(varpart(vare, ~ N, ~ P, data = chem, chisquare = TRUE, permutations = 9)))
  expect_true(any(grepl('by CCA (adjusted R-squared from 9 permutations)', out, fixed = TRUE)))



  # Errors. ====

  expect_error(varpart(vare, ~ N, data = chem), 'two to four sets')
  expect_error(varpart(vare, ~ N, ~ P, ~ K, ~ Ca, ~ Mg, data = chem), 'two to four sets')
  expect_error(varpart(vare, a = ~ N, a = ~ P, data = chem), 'different names')
  expect_error(varpart(bray(vare), ~ N, ~ P, data = chem, norm = 'hellinger'), '`norm` applies')
  expect_error(varpart(bray(vare), ~ N, ~ P, data = chem, chisquare = TRUE), 'needs a count table')
  expect_error(varpart(vare, ~ N, ~ P, data = chem, chisquare = TRUE, norm = 'hellinger'), "must be 'none'")
  expect_error(varpart(vare, ~ N, ~ P, data = chem, chisquare = TRUE, permutations = 0), 'at least 1')
  expect_error(varpart(vare, ~ N, ~ P, data = chem, chisquare = NA), 'TRUE or FALSE')
  expect_error(varpart(vare ~ N, ~ P, data = chem), 'not a formula')
  expect_error(varpart(vare, ~ N + Condition(P), ~ K, data = chem), 'Condition')
  expect_error(varpart(vare, y ~ N, ~ K, data = chem), 'one-sided')
  expect_error(varpart(vare, ~ N, ~ K), '`data` is required')
  expect_error(varpart(vare, ~ N, ~ K, data = list(1)), 'data.frame')
  expect_error(varpart(vare, ~ N, ~ nothere, data = chem), 'Not found')
  expect_error(varpart(vare, ~ N, list(1), data = chem), 'must be a formula')
  expect_error(varpart(vare, ~ N, matrix('a', 24, 1), data = chem), 'must be numeric')
  expect_error(varpart(vare, ~ N, matrix(1, 23, 1), data = chem), '23 rows')
  expect_error(varpart(vare, ~ N, matrix(c(NA, 1:23), 24, 1), data = chem), 'missing values')
  expect_error(varpart(vare, ~ N, matrix(c(Inf, 1:23), 24, 1), data = chem), 'infinite')
  expect_error(varpart(vare, ~ N, data.frame(a = c(NA, 1:23)), data = chem), 'missing values')
  expect_error(varpart(vare, ~ N, data.frame(a = c(Inf, 1:23)), data = chem), 'infinite')
  expect_error(varpart(vare, ~ N, data.frame(row.names = 1:24), data = chem), 'no columns')
  expect_error(varpart(vare, ~ N, rep(1, 24), data = chem), 'no variation')
  expect_error(varpart(vare, ~ N, list2 = sum, data = chem), 'must be a formula')
  expect_error(varpart(vare[1:2, ], ~ N, ~ P, data = chem[1:2, ]), 'at least 3')
  expect_error(varpart(matrix(1, 5, 3), 1:5, 5:1), 'no variation')
  expect_error(varpart(vare, ~ N, ~ P, data = chem, norm = 'nope'))

  zero <- vare
  zero[3, ] <- 0
  expect_error(varpart(zero, ~ N, ~ P, data = chem, chisquare = TRUE, permutations = 9), 'no observations')

#})
