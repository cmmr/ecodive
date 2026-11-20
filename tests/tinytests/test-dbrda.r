#test_that("dbrda", {

  # Inputs and reference results, written by data-raw/expected/dbrda.r.
  ref <- readRDS(file.path('expected', 'dbrda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case_formula <- function (cs) {
    rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
    stats::as.formula(paste('~', paste(rhs, collapse = ' + ')))
  }

  fit_case <- function (cs, cpus = 1) {
    if (is.null(cs$dm)) cs$dm <- dist(cs$counts)
    suppressMessages(dbrda(cs$dm, data = cs$data, formula = case_formula(cs), n_axes = cs$k, cpus = cpus))
  }



  # Same as the reference, up to the sign of each axis: stored results, so
  # at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      exp <- cs$expected
      k   <- ncol(exp$sites)
      lbl <- cs$name

      # Euclidean distances between counts, or the distances themselves.
      if (is.null(cs$dm)) cs$dm <- dist(cs$counts)

      if (is.null(cs$warn)) {
        ord <- fit_case(cs)
      } else {
        expect_warning(ord <- fit_case(cs), cs$warn, info = lbl)
      }

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

      # vegan gives no R-squared without constraints; dbrda() gives 0.
      if (length(exp$r_squared) == 0) {
        expect_identical(attr(ord, 'r_squared'),     0, info = lbl)
        expect_identical(attr(ord, 'adj_r_squared'), 0, info = lbl)
      } else {
        expect_equal(attr(ord, 'r_squared'),     exp$r_squared,     tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'adj_r_squared'), exp$adj_r_squared, tolerance = 1e-10, info = lbl)
      }

      bp <- attr(ord, 'biplot')
      if (is.null(exp$biplot)) { # Conditions only: no constraints.
        expect_identical(nrow(bp), 0L, info = lbl)
      } else {
        expect_identical(dimnames(bp), dimnames(exp$biplot), info = lbl)
        expect_equal(sweep(bp, 2, s, '*'), exp$biplot, tolerance = 1e-10, info = lbl)
      }

      cn <- attr(ord, 'centroids')
      if (is.null(exp$centroids)) {
        expect_null(cn, info = lbl)
      } else {
        expect_identical(dimnames(cn), dimnames(exp$centroids), info = lbl)
        expect_equal(sweep(cn, 2, s, '*'), exp$centroids, tolerance = 1e-10, info = lbl)
      }
    }
  }



  # Result format. ====

  cs  <- ref$cases[[1]]
  ord <- fit_case(cs)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(no_display(ord))),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'inertia',
      'constrained_eig', 'r_squared', 'adj_r_squared', 'biplot', 'centroids', 'class') )
  expect_identical(attr(ord, 'method'), 'dbrda')

  # Inertia components add up.
  expect_equal(attr(ord, 'trace'), sum(attr(ord, 'inertia')))

  # R-squared is the constrained share of the inertia, as permanova()'s.
  expect_equal(attr(ord, 'r_squared'), attr(ord, 'inertia')[['constrained']] / attr(ord, 'trace'))
  pm <- permanova(cs$dm, data = cs$data, formula = ~ Management + A1, permutations = 0, cpus = 1)
  expect_equal(attr(ord, 'r_squared'), pm$r2, tolerance = 1e-12)

  # Adjusted for 4 constraint columns, 20 samples.
  expect_equal(attr(ord, 'adj_r_squared'), 1 - (1 - attr(ord, 'r_squared')) * 19 / 15)



  # Ways to give the model: as permanova(). ====

  d    <- cs$dm
  meta <- cs$data
  one  <- dbrda(d, data = meta, formula = ~ Management, cpus = 1)

  # `groups`: a column of `data`, a vector parallel to its rows, or a vector
  # with one value per sample.
  expect_identical(no_display(dbrda(d, 'Management', meta, cpus = 1)), no_display(one))
  expect_equal(unclass(dbrda(d, meta$Management, meta, cpus = 1))[, ], unclass(one)[, ])
  expect_equal(unclass(dbrda(d, setNames(meta$Management, rownames(meta)), cpus = 1))[, ], unclass(one)[, ])
  expect_equal(unclass(dbrda(d, meta$Management, cpus = 1))[, ], unclass(one)[, ])

  # vegan's style: the distances on the formula's left-hand side, with the
  # metadata second or named.
  dmx <- d
  expect_identical(no_display(dbrda(dmx ~ Management, meta, cpus = 1)), no_display(one))
  expect_identical(no_display(dbrda(dmx ~ Management, data = meta, cpus = 1)), no_display(one))
  expect_identical(no_display(dbrda(dmx ~ Management + A1, meta, cpus = 1)), no_display(ord))
  expect_error(dbrda(dmx ~ A1, meta, ~ A1), 'Name the formula argument')
  expect_error(dbrda(dmx ~ A1, data = meta, formula = ~ A1), 'not both')

  # dist_ord() passes the same arguments on.
  expect_identical(no_display(dist_ord(d, 'dbrda', 'Management', meta, cpus = 1)), no_display(one))
  expect_identical(dist_ord(d, 'dbrda', data = meta, formula = ~ Management, cpus = 1), one)
  expect_identical(match_dist_ord('distance-based RDA')$id, 'dbrda')

  # dist_ord() drops `seed` for dbrda, which has no random parts.
  expect_identical(no_display(dist_ord(d, 'dbrda', 'Management', meta, seed = 1, cpus = 1)), no_display(one))

  # A symmetric matrix gives the same result as a dist.
  expect_identical(no_display(dbrda(as.matrix(d), 'Management', meta, cpus = 1)), no_display(one))

  # Rows of `data` are matched to the samples by name.
  expect_identical(no_display(dbrda(d, 'Management', meta[rev(seq_len(nrow(meta))), ], cpus = 1)), no_display(one))

  # A model is required: no default formula, and an explicit (misspelled)
  # NULL `groups` is an error.
  expect_error(dbrda(d),                               '`groups` is required')
  expect_error(dbrda(d, data = meta),                  '`groups` is required')
  expect_error(dbrda(d, meta$Managment, meta),         '`groups` is NULL')
  expect_error(dbrda(d, NULL, meta),                   '`groups` is NULL')
  expect_error(dbrda(d, formula = ~ Management),       '`data` is required')
  expect_error(dbrda(d, meta, ~ A1),                   'Name the formula argument')
  expect_error(dbrda(d, 'Management', meta, formula = ~ A1), 'not both')
  expect_error(dbrda(d, meta),                         '`groups` must be a vector')
  expect_error(dbrda(d, 'Mangement', meta),            'not one')
  expect_error(dbrda(d, 'Management'),                 'not given')

  # Only conditions: the unconstrained ordination after removing them.
  oc <- dbrda(d, data = meta, formula = ~ Condition(A1), cpus = 1)
  expect_true(all(startsWith(colnames(oc), 'MDS')))
  expect_true(attr(oc, 'inertia')[['conditional']] > 0)
  expect_equal(attr(oc, 'inertia')[['constrained']], 0)
  expect_equal(attr(oc, 'trace'), sum(attr(oc, 'inertia')))
  expect_identical(attr(oc, 'r_squared'),     0)
  expect_identical(attr(oc, 'adj_r_squared'), 0)

  # `n_axes` is validated.
  expect_error(dbrda(d, 'Management', meta, n_axes = 0),   '`n_axes`')
  expect_error(dbrda(d, 'Management', meta, n_axes = 1.5), '`n_axes`')



  # Scale relative to pcoa(): the same, unless a distance is 4 or more. ====

  # Bray-Curtis: every distance is below 4.
  ds  <- ref$cases[[1]]$dm
  pcs <- pcoa(ds, n_axes = 2, cpus = 1)
  ords <- dbrda_core(dist_vector(ds), NULL, k = 2, cpus = 1)
  expect_true(max(ds) < 4)
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace'), tolerance = 1e-12)
  expect_equal(attr(ords, 'eig'),   attr(pcs, 'eig'),   tolerance = 1e-10)

  # Euclidean distances between counts: dbrda's are pcoa's / (ns - 1).
  ds  <- ref$cases[[6]]$dm
  ns  <- attr(ds, 'Size')
  pcs <- pcoa(ds, n_axes = 2, cpus = 1)
  ords <- dbrda_core(dist_vector(ds), NULL, k = 2, cpus = 1)
  expect_true(max(ds) >= 4)
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace') / (ns - 1), tolerance = 1e-12)
  expect_equal(attr(ords, 'eig'),   attr(pcs, 'eig')   / (ns - 1), tolerance = 1e-10)
  expect_equal(attr(ords, 'prop'),  attr(pcs, 'prop'),            tolerance = 1e-10)
  ords <- fit_case(ref$cases[[6]])
  expect_equal(attr(ords, 'trace'), attr(pcs, 'trace') / (ns - 1), tolerance = 1e-12)



  # Unconstrained axes are orthogonal to the constraints, the conditions and
  # the ones vector (the Lanczos basis is kept orthogonal to them). ====

  for (cs in ref$cases[8:9]) {
    ord <- fit_case(cs)
    mds <- unclass(ord)[, startsWith(colnames(ord), 'MDS'), drop = FALSE]
    des <- cbind(1, cs$x, cs$z)
    des <- des / rep(sqrt(colSums(des ^ 2)), each = nrow(des))
    mds <- mds / rep(sqrt(colSums(mds ^ 2)), each = nrow(mds))
    expect_true(ncol(mds) >= 1, info = cs$name)
    expect_true(max(abs(crossprod(des, mds))) < 1e-10, info = cs$name)
  }



  # Deterministic, for any number of threads. ====

  for (cs in ref$cases[8:9]) {
    ord <- fit_case(cs, cpus = 1)
    expect_identical(fit_case(cs, cpus = 1),        ord)
    expect_identical(fit_case(cs, cpus = 2),        ord)
    expect_identical(fit_case(cs, cpus = max_cpus), ord)
  }



  # Aliased (redundant) constraints are dropped with a message. ====

  cs  <- ref$cases[[1]]
  dat <- cbind(cs$data, A1x2 = 2 * cs$data$A1)
  expect_message(ord2 <- dbrda(cs$dm, data = dat, formula = ~ Management + A1 + A1x2, cpus = 1), 'linear combinations')
  ord <- fit_case(cs)
  expect_identical(rownames(attr(ord2, 'biplot')), colnames(cs$x))
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ])
  # ... and R-squared is adjusted for the columns kept.
  expect_equal(attr(ord2, 'adj_r_squared'), attr(ord, 'adj_r_squared'))

  x2 <- cbind(cs$x, const = 1)
  expect_message(ord2 <- dbrda_core(dist_vector(cs$dm), x2, k = 3, cpus = 1), 'linear combinations')
  expect_equal(unclass(ord2)[, ], unclass(ord)[, ])

  # Constraints entirely within the conditions leave only MDS axes.
  expect_message(ord3 <- dbrda_core(dist_vector(cs$dm), cs$x[, 'A1', drop = FALSE], z = cs$x, k = 2, cpus = 1), 'linear combinations')
  expect_identical(colnames(ord3), c('MDS1', 'MDS2'))
  expect_equal(attr(ord3, 'inertia')[['constrained']], 0)



  # R-squared without residual degrees of freedom. ====

  sat <- dbrda(cs$dm, factor(seq_len(20)), cpus = 1)
  expect_equal(attr(sat, 'r_squared'), 1)
  expect_identical(attr(sat, 'adj_r_squared'), NA_real_)



  # Edge cases. ====

  # No constraints: the same axes as PCoA, named MDS.
  d   <- cs$dm
  ord <- dbrda_core(dist_vector(d), NULL, k = 2, cpus = 1)
  pc  <- pcoa(d, n_axes = 2, cpus = 1)
  expect_identical(colnames(ord), c('MDS1', 'MDS2'))
  expect_equal(attr(ord, 'eig'), attr(pc, 'eig'), tolerance = 1e-10)
  expect_equal(attr(ord, 'inertia')[['unconstrained']], attr(ord, 'trace'))

  # More axes than exist.
  ord <- suppressWarnings(dbrda_core(dist_vector(as.dist(as.matrix(d)[1:4, 1:4])), cs$x[1:4, 'A1', drop = FALSE], k = 10, cpus = 1))
  expect_true(ncol(ord) <= 3)

  # Axes without positive eigenvalues are dropped with a warning.
  pts <- cbind(x = sin(1:40), y = cos(1:40 / 2) * 3)
  expect_warning(ord <- dbrda_core(dist_vector(dist(pts)), NULL, k = 3, cpus = 1), 'only 2 of the first 3')
  expect_identical(colnames(ord), c('MDS1', 'MDS2'))

  # Two samples: no tiny negative inertia from rounding.
  ord <- suppressWarnings(dbrda(dist(c(a = 1, b = 3)), c(a = 1, b = 2), cpus = 1))
  expect_true(all(attr(ord, 'inertia') >= 0))

  # One sample.
  ord <- dbrda_core(dist_vector(as.dist(matrix(0, 1, 1))), matrix(1, 1, 1), k = 3, cpus = 1)
  expect_identical(dim(ord), c(1L, 0L))
  expect_identical(attr(ord, 'r_squared'), NA_real_)

  # Non-finite values.
  bad <- d
  bad[[3]] <- NA
  expect_error(dbrda_core(dist_vector(bad), cs$x, k = 3, cpus = 1), 'NA')
  xb <- cs$x
  xb[2, 1] <- NA
  expect_error(dbrda_core(dist_vector(d), xb, k = 3, cpus = 1), 'NA')

#})
