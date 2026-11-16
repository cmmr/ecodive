#test_that("dbrda", {

  # Inputs and reference results, written by data-raw/expected/dbrda.r.
  ref <- readRDS(file.path('expected', 'dbrda.rds'))

  case_formula <- function (cs) {
    rhs <- cs$constraints
    if (!is.null(cs$conditions)) rhs <- paste0(rhs, ' + Condition(', cs$conditions, ')')
    stats::as.formula(paste('~', rhs))
  }

  fit_case <- function (cs, cpus = 1) {
    dbrda(cs$dm, cs$data, case_formula(cs), k = cs$k, cpus = cpus)
  }



  # Same as the reference, up to the sign of each axis. ====

  for (cs in ref$cases) {

    ord <- fit_case(cs)
    exp <- cs$expected
    k   <- ncol(exp$sites)
    lbl <- cs$name

    expect_identical(dim(ord), c(attr(cs$dm, 'Size'), k), info = lbl)
    expect_identical(dimnames(ord), dimnames(exp$sites), info = lbl)

    # Flip axes to match the reference.
    s <- sign(colSums(unclass(ord) * exp$sites))
    expect_equal(sweep(unclass(ord), 2, s, '*'), exp$sites, check.attributes = FALSE, tolerance = 1e-10, info = lbl)

    expect_equal(attr(ord, 'eig'),             exp$eig,             tolerance = 1e-10, info = lbl)
    expect_equal(attr(ord, 'prop'),            exp$prop,            tolerance = 1e-10, info = lbl)
    expect_equal(attr(ord, 'trace'),           exp$trace,           tolerance = 1e-10, info = lbl)
    expect_equal(
      attr(ord, 'inertia'),
      c(conditional = exp$conditional, constrained = exp$constrained, unconstrained = exp$unconstrained),
      tolerance = 1e-10, info = lbl )
    expect_equal(attr(ord, 'constrained_eig'), exp$constrained_eig, tolerance = 1e-10, info = lbl)

    bp <- attr(ord, 'biplot')
    expect_identical(dimnames(bp), dimnames(exp$biplot), info = lbl)
    expect_equal(sweep(bp, 2, s, '*'), exp$biplot, tolerance = 1e-10, info = lbl)

    cn <- attr(ord, 'centroids')
    if (is.null(exp$centroids)) {
      expect_null(cn, info = lbl)
    } else {
      expect_identical(dimnames(cn), dimnames(exp$centroids), info = lbl)
      expect_equal(sweep(cn, 2, s, '*'), exp$centroids, tolerance = 1e-10, info = lbl)
    }
  }



  # Result format. ====

  cs  <- ref$cases[[1]]
  ord <- fit_case(cs)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'inertia',
      'constrained_eig', 'biplot', 'centroids') )
  expect_identical(attr(ord, 'method'), 'dbrda')

  # Inertia components add up.
  expect_equal(attr(ord, 'trace'), sum(attr(ord, 'inertia')))



  # Scale relative to pcoa(): the same, unless a distance is 4 or more. ====

  # Bray-Curtis: every distance is below 4.
  ds  <- ref$cases[[1]]$dm
  pcs <- pcoa(ds, k = 2, cpus = 1)
  ords <- dbrda_core(dist_vector(ds), NULL, k = 2, cpus = 1)
  expect_true(max(ds) < 4)
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace'), tolerance = 1e-12)
  expect_equal(attr(ords, 'eig'),   attr(pcs, 'eig'),   tolerance = 1e-10)

  # Euclidean distances between counts: dbrda's are pcoa's / (ns - 1).
  ds  <- ref$cases[[6]]$dm
  ns  <- attr(ds, 'Size')
  pcs <- pcoa(ds, k = 2, cpus = 1)
  ords <- dbrda_core(dist_vector(ds), NULL, k = 2, cpus = 1)
  expect_true(max(ds) >= 4)
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace') / (ns - 1), tolerance = 1e-12)
  expect_equal(attr(ords, 'eig'),   attr(pcs, 'eig')   / (ns - 1), tolerance = 1e-10)
  expect_equal(attr(ords, 'prop'),  attr(pcs, 'prop'),            tolerance = 1e-10)
  ords <- fit_case(ref$cases[[6]])
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace') / (ns - 1), tolerance = 1e-12)

  # A symmetric matrix gives the same result as a dist.
  expect_identical(dbrda(as.matrix(cs$dm), cs$data, case_formula(cs), cpus = 1), ord)

  # Rows of `data` are matched to the samples by name.
  shuffled <- cs$data[rev(seq_len(nrow(cs$data))), ]
  expect_identical(dbrda(cs$dm, shuffled, case_formula(cs), cpus = 1), ord)

  # The default formula uses every column of `data`; ordination() runs it.
  vars <- cs$data[, c('Management', 'A1')]
  expect_identical(dbrda(cs$dm, vars, cpus = 1), ord)
  expect_identical(ordination(cs$dm, 'dbrda', data = vars, cpus = 1), ord)
  expect_identical(match_ordination('distance-based RDA')$id, 'dbrda')

  # Not capscale, whose results differ for non-Euclidean distances.
  expect_error(match_ordination('capscale'), 'does not match')
  expect_error(match_ordination('CAP'),      'does not match')

  # `k` is validated.
  expect_error(dbrda(cs$dm, vars, k = 0),   '`k`')
  expect_error(dbrda(cs$dm, vars, k = 1.5), '`k`')



  # Deterministic, for any number of threads. ====

  for (cs in ref$cases[8:9]) {
    ord <- fit_case(cs, cpus = 1)
    expect_identical(fit_case(cs, cpus = 1), ord)
    expect_identical(fit_case(cs, cpus = 3), ord)
    expect_identical(fit_case(cs, cpus = 8), ord)
  }



  # Aliased (redundant) constraints are dropped with a message. ====

  cs  <- ref$cases[[1]]
  dat <- cbind(cs$data, A1x2 = 2 * cs$data$A1)
  expect_message(ord2 <- dbrda(cs$dm, dat, ~ Management + A1 + A1x2, cpus = 1), 'aliased')
  ord <- fit_case(cs)
  expect_identical(rownames(attr(ord2, 'biplot')), colnames(cs$x))
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ])

  x2 <- cbind(cs$x, const = 1)
  expect_message(ord2 <- dbrda_core(dist_vector(cs$dm), x2, k = 3, cpus = 1), 'aliased')
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ])

  # Constraints entirely within the conditions leave only MDS axes.
  expect_message(ord3 <- dbrda_core(dist_vector(cs$dm), cs$x[, 'A1', drop = FALSE], z = cs$x, k = 2, cpus = 1), 'aliased')
  expect_identical(colnames(ord3), c('MDS1', 'MDS2'))
  expect_equal(attr(ord3, 'inertia')[['constrained']], 0)



  # Edge cases. ====

  # No constraints: the same axes as PCoA, named MDS.
  d   <- cs$dm
  ord <- dbrda_core(dist_vector(d), NULL, k = 2, cpus = 1)
  pc  <- pcoa(d, k = 2, cpus = 1)
  expect_identical(colnames(ord), c('MDS1', 'MDS2'))
  expect_equal(attr(ord, 'eig'), attr(pc, 'eig'), tolerance = 1e-10)
  expect_equal(attr(ord, 'inertia')[['unconstrained']], attr(ord, 'trace'))

  # More axes than exist.
  ord <- dbrda_core(dist_vector(as.dist(as.matrix(d)[1:4, 1:4])), cs$x[1:4, 'A1', drop = FALSE], k = 10, cpus = 1)
  expect_true(ncol(ord) <= 3)

  # One sample.
  ord <- dbrda_core(dist_vector(as.dist(matrix(0, 1, 1))), matrix(1, 1, 1), k = 3, cpus = 1)
  expect_identical(dim(ord), c(1L, 0L))

  # Non-finite values.
  bad <- d
  bad[[3]] <- NA
  expect_error(dbrda_core(dist_vector(bad), cs$x, k = 3, cpus = 1), 'NA')
  xb <- cs$x
  xb[2, 1] <- NA
  expect_error(dbrda_core(dist_vector(d), xb, k = 3, cpus = 1), 'NA')

#})
