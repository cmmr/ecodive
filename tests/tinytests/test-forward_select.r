#test_that("forward_select", {

  # Inputs and reference results, written by data-raw/expected/forward_select.r.
  ref <- readRDS(file.path('expected', 'forward_select.rds'))
  dat <- ref$data

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case_x <- function (cs) {
    y <- dat[[cs$x$table]]
    switch(cs$x$kind, hellinger = y, bray = bray(y))
  }

  # Without the formula of the selected variables, whose environment is
  # the caller's.
  strip <- function (x) { attr(x, 'selected') <- NULL; x }

  fit_case <- function (cs, cpus = 1L, ...) {
    args <- c(
      list(case_x(cs), stats::as.formula(cs$formula), dat[[cs$env]], cpus = cpus),
      if (cs$x$kind == 'hellinger') list(norm = 'hellinger'),
      list(...) )
    do.call(forward_select, args)
  }



  # Same as vegan::ordiR2step(): stored results, so at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      lbl <- cs$name
      exp <- cs$expected
      fs  <- fit_case(cs, permutations = dat[[cs$permutations]], global_test = FALSE)
      k   <- length(exp$variable)

      # The same variables, in the same order, and the same stopping point.
      expect_identical(fs$variable, exp$variable, info = lbl)
      expect_identical(fs$df, as.integer(exp$df), info = lbl)
      expect_equal(fs$adj_r_squared_cum, exp$adj, tolerance = 1e-10, info = lbl)
      expect_equal(fs$F, exp$F, tolerance = 1e-10, info = lbl)
      expect_identical(fs$p_value, exp$p, info = lbl)
      expect_equal(attr(fs, 'global')$adj_r_squared, exp$global, tolerance = 1e-10, info = lbl)

      # The increments add up.
      expect_equal(cumsum(fs$adj_r_squared), fs$adj_r_squared_cum, tolerance = 1e-12, info = lbl)
      expect_equal(cumsum(fs$r_squared), fs$r_squared_cum, tolerance = 1e-12, info = lbl)
      expect_identical(
        attr(fs, 'selected'),
        if (k) stats::reformulate(exp$variable, env = environment(attr(fs, 'selected')))
        else stats::as.formula(~ 1, env = environment(attr(fs, 'selected'))), info = lbl )
    }


    # Why the reference cases stopped.
    cs <- ref$cases[[which(vapply(ref$cases, `[[`, '', 'name') == 'vare_hel')]]
    fs <- fit_case(cs, permutations = dat$perm_vare, global_test = FALSE)
    expect_identical(attr(fs, 'stop'), 'adjusted R-squared would exceed the global model\'s')
    expect_identical(attr(fs, 'next_best')$variable, 'K')
    expect_true(attr(fs, 'next_best')$adj_r_squared_cum > attr(fs, 'global')$adj_r_squared)
    expect_true(is.na(attr(fs, 'next_best')$p_value))

    cs <- ref$cases[[which(vapply(ref$cases, `[[`, '', 'name') == 'mite_inter')]]
    fs <- fit_case(cs, permutations = dat$perm_mite, global_test = FALSE)
    expect_identical(attr(fs, 'stop'), 'not significant')
    expect_identical(attr(fs, 'next_best')$variable, 'SubsDens:Topo')  # after its main effects
    expect_true(attr(fs, 'next_best')$p_value > 0.05)
    # That last candidate is the global model itself: the same terms give
    # the same numbers whatever order they were chosen in, so it can't
    # "exceed" the global adjusted R-squared by rounding (as it did on CI).
    expect_identical(attr(fs, 'next_best')$adj_r_squared_cum, attr(fs, 'global')$adj_r_squared)


    # adespatial::forward.sel(): the same order, R-squared, and F. Its
    # p-values are from its own permutations, so select here by adjusted
    # R-squared alone (permutations = 0) and compare its steps.
    fw  <- ref$forward_sel
    num <- dat$mite.vars[, vapply(dat$mite.vars, is.numeric, TRUE)]
    fs  <- forward_select(dat$mite, data = num, norm = 'hellinger', permutations = 0L)
    k   <- nrow(fw)
    expect_true(nrow(fs) >= k)
    expect_identical(fs$variable[1:k], as.character(fw$variable))
    expect_equal(fs$r_squared[1:k], fw$r2, tolerance = 1e-10)
    expect_equal(fs$r_squared_cum[1:k], fw$r2_cum, tolerance = 1e-10)
    expect_equal(fs$adj_r_squared_cum[1:k], fw$adj_cum, tolerance = 1e-10)
    expect_equal(fs$F[1:k], fw$F, tolerance = 1e-10)
    expect_true(all(is.na(fs$p_value)))
    expect_identical(attr(fs, 'n_perm'), 0L)


    # Same as ecodive's own methods: RDA of the Hellinger table equals
    # db-RDA of the Euclidean distances among its samples; each step's F
    # and p-value equal permanova()'s of the larger model by terms, with
    # the full permutation (no conditions); the global test equals
    # permanova()'s overall test.
    cs  <- ref$cases[[1]]
    a   <- fit_case(cs, permutations = dat$perm_mite)
    de  <- euclidean(dat$mite, norm = 'hellinger')
    b   <- forward_select(de, ~ ., dat$mite.vars, permutations = dat$perm_mite)
    expect_identical(a$variable, b$variable)
    expect_equal(a$adj_r_squared_cum, b$adj_r_squared_cum, tolerance = 1e-10)
    expect_equal(a$F, b$F, tolerance = 1e-10)
    expect_identical(a$p_value, b$p_value)
    expect_identical(attr(b, 'method'), 'dbrda')
    expect_identical(attr(a, 'method'), 'rda')
    expect_identical(attr(a, 'norm'), 'hellinger')

    pm <- permanova(de, data = dat$mite.vars, formula = ~ ., permutations = dat$perm_mite)
    expect_equal(attr(a, 'global')$F, pm$table$F[[1]], tolerance = 1e-10)
    expect_identical(attr(a, 'global')$p_value, pm$p_value)
    expect_equal(attr(a, 'global')$r_squared, pm$r2, tolerance = 1e-10)

    pm <- permanova(
      de, data = dat$mite.vars, by = 'terms', permutations = dat$perm_mite,
      formula = stats::reformulate(a$variable[1:3]) )
    expect_equal(a$F[[3]], pm$table$F[[3]], tolerance = 1e-10)  # the same residual

    # The first step is a test without conditions, as permanova()'s.
    pm <- permanova(de, data = dat$mite.vars, formula = ~ WatrCont, permutations = dat$perm_mite)
    expect_identical(a$p_value[[1]], pm$p_value)
  }



  # Results don't depend on the number of threads. ====

  cs <- ref$cases[[which(vapply(ref$cases, `[[`, '', 'name') == 'mite_cond')]]
  a  <- fit_case(cs, permutations = 49L, seed = 2L)
  for (cpus in 2:max_cpus)
    expect_identical(strip(fit_case(cs, cpus = cpus, permutations = 49L, seed = 2L)), strip(a), info = cpus)

  cs <- ref$cases[[which(vapply(ref$cases, `[[`, '', 'name') == 'mite_bray')]]
  a  <- fit_case(cs, permutations = 49L)
  for (cpus in 2:max_cpus)
    expect_identical(strip(fit_case(cs, cpus = cpus, permutations = 49L)), strip(a), info = cpus)
  expect_identical(attr(a, 'permutation'), 'free')
  expect_identical(attr(a, 'n_perm'), 49L)
  expect_true(all(a$p_value <= 0.05))
  expect_true(all(a$p_value >= 1 / 50))



  # A small, simulated data set. ====

  set.seed(1)
  n    <- 30
  meta <- data.frame(
    row.names = paste0('S', 1:n),
    ph   = runif(n, 5, 8),
    temp = runif(n, 10, 20),
    z1   = rnorm(n),
    z2   = rnorm(n),
    site = factor(rep(c('a', 'b', 'c'), 10)),
    blk  = rep(c('u', 'v'), each = 15) )
  counts <- matrix(rpois(n * 12, 3), n, 12, dimnames = list(rownames(meta), paste0('F', 1:12)))
  counts[, 1:4] <- counts[, 1:4] + round(10 * (meta$ph - 5))
  counts[, 5:8] <- counts[, 5:8] + round(2 * (meta$temp - 10))

  fs <- forward_select(counts, ~ ph + temp + z1 + z2 + site, meta, norm = 'hellinger', permutations = 99L)
  expect_inherits(fs, c('ecodive_forward_select', 'data.frame'))
  expect_identical(names(fs), c('variable', 'df', 'r_squared', 'r_squared_cum',
                                'adj_r_squared', 'adj_r_squared_cum', 'F', 'p_value'))
  expect_true(all(c('ph', 'temp') %in% fs$variable))
  expect_true(attr(fs, 'global')$p_value <= 0.05)
  expect_identical(attr(fs, 'n'), 30L)
  out <- capture.output(print(fs))
  expect_true(any(grepl('Forward selection by RDA (hellinger): 30 samples', out, fixed = TRUE)))
  expect_true(any(grepl('Selected: ~', out, fixed = TRUE)))

  # Default formula: every column of `data` except `strata`; samples are
  # matched by name.
  fs1 <- forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 99L)
  fs2 <- forward_select(counts, data = meta[rev(rownames(meta)), ], strata = 'blk',
                        norm = 'hellinger', permutations = 99L)
  expect_identical(fs1, forward_select(counts, ~ ph + temp + z1 + z2 + site, meta,
                                       norm = 'hellinger', permutations = 99L))
  expect_identical(attr(fs2, 'permutation'), 'within strata')
  expect_identical(fs2$variable, fs1$variable)
  expect_equal(fs2$adj_r_squared_cum, fs1$adj_r_squared_cum, tolerance = 1e-12)
  expect_false(identical(attr(fs2, 'next_best')$p_value, attr(fs1, 'next_best')$p_value))
  expect_identical(
    forward_select(counts, data = meta[, 1:5], strata = meta$blk, norm = 'hellinger', permutations = 99L),
    fs2 )

  # max_vars.
  fm <- forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 99L, max_vars = 1)
  expect_identical(nrow(fm), 1L)
  expect_identical(attr(fm, 'stop'), 'max_vars reached')
  expect_identical(fm$variable, fs1$variable[[1]])
  expect_identical(attr(fm, 'next_best')$variable, fs1$variable[[2]])

  # Seeds.
  expect_identical(forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 99L), fs1)
  expect_false(identical(
    attr(forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 99L, seed = 7L), 'next_best')$p_value,
    attr(fs1, 'next_best')$p_value ))

  # alpha = 1: every step that raises adjusted R-squared is kept.
  fa <- forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 19L, alpha = 1)
  fz <- forward_select(counts, data = meta[, 1:5], norm = 'hellinger', permutations = 0L)
  expect_identical(fa$variable, fz$variable)
  expect_equal(fa$F, fz$F, tolerance = 1e-12)

  # Nothing to find: the global test stops selection, with a warning.
  noise <- meta[, c('z1', 'z2')]
  set.seed(9)
  rnd <- matrix(rpois(n * 12, 5), n, 12)
  expect_warning(
    fn <- forward_select(rnd, data = noise, norm = 'hellinger', permutations = 99L),
    'global model is not significant' )
  expect_identical(nrow(fn), 0L)
  expect_identical(attr(fn, 'stop'), 'global test not significant')
  expect_identical(deparse(attr(fn, 'selected')), '~1')
  expect_identical(environment(attr(fn, 'selected')), environment())
  expect_true(attr(fn, 'global')$p_value > 0.05)
  expect_true(any(grepl('No variable selected', capture.output(print(fn)))))
  fn2 <- forward_select(rnd, data = noise, norm = 'hellinger', permutations = 99L, global_test = FALSE)
  expect_true(is.na(attr(fn2, 'global')$p_value))
  expect_true(attr(fn2, 'stop') %in% c('adjusted R-squared does not increase', 'not significant',
                                       'adjusted R-squared would exceed the global model\'s'))

  # A constant candidate is never chosen.
  fc <- forward_select(counts, ~ ph + k, cbind(meta, k = 1), norm = 'hellinger', permutations = 19L)
  expect_identical(fc$variable, 'ph')
  expect_identical(attr(fc, 'stop'), 'adjusted R-squared does not increase')

  # A candidate aliased with Condition() terms.
  expect_error(
    forward_select(counts, ~ ph + Condition(ph), meta, norm = 'hellinger'),
    'No candidate varies' )

  # Saturated global model: no adjusted R-squared, no global test.
  small <- meta[1:6, ]
  fsat <- suppressWarnings(
    forward_select(counts[1:6, ], ~ ph + temp + z1 + z2 + site, small, permutations = 9L) )
  expect_warning(
    forward_select(counts[1:6, ], ~ ph + temp + z1 + z2 + site, small, permutations = 9L, global_test = FALSE),
    'adjusted R-squared is undefined' )
  expect_true(is.na(attr(fsat, 'global')$p_value))
  expect_true(is.na(attr(fsat, 'global')$adj_r_squared))



  # Invalid arguments. ====

  expect_error(forward_select(counts, ~ ph, meta, alpha = 0), '`alpha`')
  expect_error(forward_select(counts, ~ ph, meta, alpha = 2), '`alpha`')
  expect_error(forward_select(counts, ~ ph, meta, max_vars = 0), '`max_vars`')
  expect_error(forward_select(counts, ~ ph, meta, max_vars = 1.5), '`max_vars`')
  expect_error(forward_select(counts, ~ ph, meta, global_test = NA), '`global_test`')
  expect_error(forward_select(counts, ~ ph), '`data` is required')
  expect_error(forward_select(counts, meta, ~ ph), 'Name the formula|Give `formula` before|`data` must be')
  expect_error(forward_select(counts, ~ ph, data = 1:3), '`data` must be a data.frame')
  expect_error(forward_select(counts, ~ nope, meta), 'Not found in `data`')
  expect_error(forward_select(counts[1:2, ], ~ ph, meta[1:2, ]), 'at least 3 samples')
  expect_error(forward_select(counts, ~ ph, meta, permutations = -1), '`permutations`')
  expect_error(forward_select(counts, ~ ph, meta, strata = 'blk', permutations = rbind(1:30, 30:1)),
               '`strata` or a `permutations` matrix')
  expect_error(forward_select(bray(counts), ~ ph, meta, norm = 'hellinger'), '`norm` applies')
  expect_error(forward_select(~ ph, ~ ph, meta), 'not a formula')
  expect_error(forward_select(matrix(1, n, 4, dimnames = list(rownames(meta), NULL)), ~ ph, meta),
               'no variation')
  expect_error(forward_select(counts, ~ Condition(ph), meta), 'no terms to test')
  expect_error(forward_select(counts, ~ ph, meta, cpus = 0), '`cpus`')

#})
