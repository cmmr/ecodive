#test_that("rda", {

  # Inputs and reference results, written by data-raw/expected/rda.r.
  ref <- readRDS(file.path('expected', 'rda.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case_formula <- function (cs) {
    rhs <- c(cs$constraints, if (!is.null(cs$conditions)) paste0('Condition(', cs$conditions, ')'))
    if (length(rhs) == 0) return (NULL)
    stats::as.formula(paste('~', paste(rhs, collapse = ' + ')))
  }

  fit_case <- function (cs, cpus = 1, counts = cs$counts, ...) {
    f <- case_formula(cs)
    suppressMessages(if (is.null(f)) {
      rda(counts, norm = cs$norm, pseudocount = cs$pseudocount, scale = cs$scale,
          n_axes = cs$k, cpus = cpus, ...)
    } else {
      rda(counts, data = cs$data, formula = f, norm = cs$norm,
          pseudocount = cs$pseudocount, scale = cs$scale, n_axes = cs$k, cpus = cpus, ...)
    })
  }

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }



  # Same as vegan::rda(), up to the sign of each axis: stored results, so at
  # home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      exp <- cs$expected
      k   <- ncol(exp$sites)
      lbl <- cs$name

      if (is.null(cs$warn)) {
        ord <- fit_case(cs)
      } else {
        expect_warning(ord <- fit_case(cs), cs$warn, info = lbl)
      }

      expect_identical(dim(ord), c(nrow(cs$counts), k), info = lbl)
      expect_identical(dimnames(ord), dimnames(exp$sites), info = lbl)

      # Flip axes to match the reference.
      s <- sign(colSums(unclass(ord) * exp$sites))
      expect_equal(sweep(unclass(ord), 2, s, '*'), exp$sites, check.attributes = FALSE, tolerance = 1e-10, info = lbl)

      ld <- attr(ord, 'loadings')
      expect_identical(dimnames(ld), dimnames(exp$loadings), info = lbl)
      expect_equal(sweep(ld, 2, s, '*'), exp$loadings, tolerance = 1e-10, info = lbl)

      expect_equal(attr(ord, 'eig'),             exp$eig,             tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'),            exp$prop,            tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'trace'),           exp$trace,           tolerance = 1e-10, info = lbl)
      expect_equal(
        attr(ord, 'inertia'),
        c(conditional = exp$conditional, constrained = exp$constrained, unconstrained = exp$unconstrained),
        tolerance = 1e-10, info = lbl )
      expect_equal(attr(ord, 'constrained_eig'), exp$constrained_eig, tolerance = 1e-10, info = lbl)

      # vegan gives no R-squared without constraints; rda() gives 0.
      if (length(exp$r_squared) == 0) {
        expect_identical(attr(ord, 'r_squared'),     0, info = lbl)
        expect_identical(attr(ord, 'adj_r_squared'), 0, info = lbl)
      } else {
        expect_equal(attr(ord, 'r_squared'),     exp$r_squared,     tolerance = 1e-10, info = lbl)
        expect_equal(attr(ord, 'adj_r_squared'), exp$adj_r_squared, tolerance = 1e-10, info = lbl)
      }

      bp <- attr(ord, 'biplot')
      if (is.null(exp$biplot)) { # No constraints.
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


    # The permutation test of an RDA is PERMANOVA of Euclidean distances
    # between the normalized samples: vegan's anova.cca() pseudo-F and
    # p-value, for the same permutations.
    pr <- ref$perm
    de <- euclidean(pr$counts, norm = 'hellinger', cpus = 1)
    for (m in pr[c('full', 'partial')]) {
      pm <- permanova(de, data = pr$data, formula = stats::as.formula(m$formula),
                      permutations = pr$permutations, cpus = 1)
      expect_equal(unname(pm$statistic), m$f, tolerance = 1e-10, info = m$formula)
      expect_equal(pm$p_value, m$p, info = m$formula)
      ord <- rda(pr$counts, data = pr$data, formula = stats::as.formula(m$formula), norm = 'hellinger', cpus = 1)
      if (m$formula == '~ Management + A1')
        expect_equal(attr(ord, 'r_squared'), pm$r2, tolerance = 1e-10)
    }


    # scikit-bio's rda() in scaling 2, every axis. Its "eigvals" are
    # singular values s of the centred table: rda()'s eig is s^2 / (n - 1).
    sk  <- ref$skbio
    m   <- ncol(sk$sites)
    ord <- rda(sk$counts, data = sk$x, formula = ~ ManagementHF + ManagementNM + ManagementSF + A1,
               norm = 'hellinger', n_axes = m, cpus = 1)
    expect_identical(ncol(ord), m)
    s <- sign(colSums(unclass(ord) * sk$sites))
    expect_equal(sweep(unclass(ord), 2, s, '*'),         sk$sites,   check.attributes = FALSE, tolerance = 1e-10)
    expect_equal(sweep(attr(ord, 'loadings'), 2, s, '*'), sk$species, check.attributes = FALSE, tolerance = 1e-10)
    expect_equal(sweep(attr(ord, 'biplot'), 2, s, '*'),   sk$biplot,  check.attributes = FALSE, tolerance = 1e-10)
    expect_equal(attr(ord, 'eig'), unname(sk$eigvals ^ 2 / (nrow(sk$counts) - 1)), tolerance = 1e-10)
  }



  # Result format. ====

  cs  <- case('dune, hellinger, factor + numeric')
  ord <- fit_case(cs)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(no_steps(ord))),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'inertia',
      'constrained_eig', 'r_squared', 'adj_r_squared', 'biplot', 'centroids',
      'loadings', 'class') )
  expect_identical(attr(ord, 'method'), 'rda')
  expect_identical(colnames(ord), c('RDA1', 'RDA2', 'RDA3'))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))

  # Inertia components add up, and R-squared is the constrained share.
  expect_equal(attr(ord, 'trace'), sum(attr(ord, 'inertia')))
  expect_equal(attr(ord, 'r_squared'), attr(ord, 'inertia')[['constrained']] / attr(ord, 'trace'))
  expect_equal(attr(ord, 'adj_r_squared'), 1 - (1 - attr(ord, 'r_squared')) * 19 / 15)

  # Constrained axes, then PCs.
  cs6 <- case('dune, k > constrained axes')
  ord6 <- fit_case(cs6)
  expect_identical(colnames(ord6), c('RDA1', 'RDA2', 'PC1', 'PC2', 'PC3', 'PC4'))
  expect_identical(colnames(attr(ord6, 'loadings')), colnames(ord6))
  expect_true(all(attr(ord6, 'biplot')[, 3:6] == 0))



  # Ways to give the model, as dbrda(). ====

  counts <- cs$counts
  meta   <- cs$data
  one    <- rda(counts, data = meta, formula = ~ Management, norm = 'hellinger', cpus = 1)

  expect_identical(no_steps(rda(counts, 'Management', meta, norm = 'hellinger', cpus = 1)), no_steps(one))
  expect_equal(unclass(rda(counts, meta$Management, norm = 'hellinger', cpus = 1))[, ], unclass(one)[, ])

  # vegan's style: the table on the formula's left-hand side.
  expect_identical(rda(counts ~ Management, meta, norm = 'hellinger', cpus = 1), one)
  expect_identical(rda(counts ~ Management, data = meta, norm = 'hellinger', cpus = 1), one)
  expect_identical(rda(counts ~ Management + A1, meta, norm = 'hellinger', cpus = 1), ord)
  expect_error(rda(counts ~ A1, meta, ~ A1), 'Name the formula argument')
  expect_error(rda(counts ~ A1, data = meta, formula = ~ A1), 'not both')

  # count_ord() passes the same arguments on, and drops `seed`.
  expect_identical(no_steps(count_ord(counts, 'rda', 'Management', meta, norm = 'hellinger', cpus = 1)), no_steps(one))
  expect_identical(count_ord(counts, 'tbrda', data = meta, formula = ~ Management, norm = 'hellinger', seed = 1, cpus = 1), one)
  expect_error(dist_ord(bray(counts), 'rda'), 'count_ord')

  # Rows of `data` are matched to the samples by name.
  expect_identical(no_steps(rda(counts, 'Management', meta[rev(seq_len(nrow(meta))), ], norm = 'hellinger', cpus = 1)), no_steps(one))

  # Errors in the model.
  expect_error(rda(counts, meta$Managment, meta),         '`groups` is NULL')
  expect_error(rda(counts, NULL, meta),                   '`groups` is NULL')
  expect_error(rda(counts, formula = ~ Management),       '`data` is required')
  expect_error(rda(counts, meta, ~ A1),                   'Name the formula argument')
  expect_error(rda(counts, 'Management', meta, formula = ~ A1), 'not both')
  expect_error(rda(counts, 'Mangement', meta),            'not one')
  bad <- meta
  bad$A1[3] <- NA
  expect_error(rda(counts, 'A1', bad), 'NA')

  # Other arguments are validated.
  expect_error(rda(counts, 'A1', meta, n_axes = 0),   '`n_axes`')
  expect_error(rda(counts, 'A1', meta, scale = NA),   '`scale`')
  expect_error(rda(counts, 'A1', meta, norm = 'xyz'), 'should be one of')
  expect_error(rda(counts, 'A1', meta, cpus = 0),     '`cpus`')
  bc <- counts + 0
  bc[2, 2] <- NA
  expect_error(rda(bc, 'A1', meta), 'NA')
  expect_error(rda(-counts, 'A1', meta, norm = 'hellinger'), 'negative')



  # Unconstrained: a PCA. Same eig, prop and trace as pca(), and the same
  # directions, in vegan's scaling 2. ====

  for (sc in c(FALSE, TRUE)) {
    for (nm in c('none', 'hellinger', 'clr')) {
      lbl <- paste(nm, sc)
      pc  <- pca(counts, norm = nm, pseudocount = if (nm == 'clr') 1, scale = sc, n_axes = 4, cpus = 1)
      ru  <- rda(counts, norm = nm, pseudocount = if (nm == 'clr') 1, scale = sc, n_axes = 4, cpus = 1)
      n   <- nrow(counts)
      expect_identical(colnames(ru), colnames(pc), info = lbl)
      expect_equal(attr(ru, 'eig'),   attr(pc, 'eig'),   tolerance = 1e-10, info = lbl)
      expect_equal(attr(ru, 'prop'),  attr(pc, 'prop'),  tolerance = 1e-10, info = lbl)
      expect_equal(attr(ru, 'trace'), attr(pc, 'trace'), tolerance = 1e-12, info = lbl)
      const <- sqrt(sqrt((n - 1) * attr(pc, 'trace'))) / sqrt((n - 1) * attr(pc, 'eig'))
      s     <- sign(colSums(unclass(ru) * unclass(pc)))
      expect_equal(sweep(unclass(ru), 2, s, '*'), sweep(unclass(pc), 2, const, '*'),
                   check.attributes = FALSE, tolerance = 1e-8, info = lbl)
      # Species scores: pca()'s unit loadings times sqrt(eig / trace) * const.
      lc <- sqrt(attr(pc, 'eig') / attr(pc, 'trace')) * sqrt(sqrt((n - 1) * attr(pc, 'trace')))
      expect_equal(sweep(attr(ru, 'loadings'), 2, s, '*'), sweep(attr(pc, 'loadings'), 2, lc, '*'),
                   check.attributes = FALSE, tolerance = 1e-8, info = lbl)
      expect_equal(attr(ru, 'inertia')[['unconstrained']], attr(ru, 'trace'), info = lbl)
    }
  }



  # tb-RDA is db-RDA of the transformed table's Euclidean distances. Every
  # Hellinger distance is below 4, so dbrda()'s eigenvalues and inertia are
  # rda()'s times n - 1, and its scores rda()'s times (n - 1)^(1/4); prop,
  # R-squared and biplot scores are the same. ====

  hd <- euclidean(counts, norm = 'hellinger', cpus = 1)
  expect_true(max(hd) < 4)
  for (f in list(~ Management + A1, ~ Management + Condition(A1))) {
    o1 <- rda(counts, data = meta, formula = f, norm = 'hellinger', n_axes = 5, cpus = 1)
    o2 <- dbrda(hd, data = meta, formula = f, n_axes = 5, cpus = 1)
    lbl <- deparse(f)
    expect_equal(attr(o2, 'constrained_eig'), attr(o1, 'constrained_eig') * 19, tolerance = 1e-10, info = lbl)
    expect_equal(attr(o2, 'eig'),     attr(o1, 'eig') * 19,     tolerance = 1e-10, info = lbl)
    expect_equal(attr(o2, 'trace'),   attr(o1, 'trace') * 19,   tolerance = 1e-12, info = lbl)
    expect_equal(attr(o2, 'inertia'), attr(o1, 'inertia') * 19, tolerance = 1e-10, info = lbl)
    expect_equal(attr(o2, 'prop'),    attr(o1, 'prop'),         tolerance = 1e-10, info = lbl)
    expect_equal(attr(o2, 'r_squared'),     attr(o1, 'r_squared'),     tolerance = 1e-10, info = lbl)
    expect_equal(attr(o2, 'adj_r_squared'), attr(o1, 'adj_r_squared'), tolerance = 1e-10, info = lbl)
    expect_equal(unclass(o2)[, ], unclass(o1)[, ] * 19^0.25, check.attributes = FALSE, tolerance = 1e-8, info = lbl)
    expect_equal(attr(o2, 'biplot'), attr(o1, 'biplot'), check.attributes = FALSE, tolerance = 1e-8, info = lbl)
  }

  # Raw counts' Euclidean distances reach 4, and dbrda() then divides by
  # n - 1 as vegan does: the same eigenvalues as rda().
  ed <- euclidean(counts, cpus = 1)
  expect_true(max(ed) >= 4)
  o1 <- rda(counts, 'Management', meta, cpus = 1)
  o2 <- dbrda(ed, 'Management', meta, cpus = 1)
  expect_equal(attr(o2, 'eig'),   attr(o1, 'eig'),   tolerance = 1e-10)
  expect_equal(attr(o2, 'trace'), attr(o1, 'trace'), tolerance = 1e-12)
  expect_equal(unclass(o2)[, ], unclass(o1)[, ], check.attributes = FALSE, tolerance = 1e-8)



  # Sparse and dense input, and either margin, give the same result. ====

  for (nm in c('hellinger', 'clr', 'none')) {
    pc  <- if (nm == 'clr') 0.5
    d1  <- rda(counts, data = meta, formula = ~ Management + Condition(A1), norm = nm,
               pseudocount = pc, n_axes = 5, cpus = 1)
    tc  <- t(counts)
    expect_identical(no_steps(rda(tc, data = meta, formula = ~ Management + Condition(A1), margin = 2,
                         norm = nm, pseudocount = pc, n_axes = 5, cpus = 1)), no_steps(d1), info = nm)
    if (requireNamespace('Matrix', quietly = TRUE)) {
      sp <- Matrix::Matrix(tc, sparse = TRUE)
      expect_true(inherits(sp, 'dgCMatrix'))
      expect_equal(no_steps(rda(sp, data = meta, formula = ~ Management + Condition(A1), margin = 2,
                       norm = nm, pseudocount = pc, n_axes = 5, cpus = 1)), no_steps(d1), tolerance = 1e-12, info = nm)
    }
    if (requireNamespace('slam', quietly = TRUE)) {
      st <- slam::as.simple_triplet_matrix(counts)
      expect_equal(rda(st, data = meta, formula = ~ Management + Condition(A1),
                       norm = nm, pseudocount = pc, n_axes = 5, cpus = 1), d1, tolerance = 1e-12, info = nm)
    }
  }



  # Deterministic, for any number of threads. ====

  for (nm in c('n = 110, hellinger, partial', 'n = 110, none, partial', 'dune, scale')) {
    cs  <- case(nm)
    ord <- fit_case(cs, cpus = 1)
    expect_identical(fit_case(cs, cpus = 1),        ord, info = nm)
    expect_identical(fit_case(cs, cpus = 2),        ord, info = nm)
    expect_identical(fit_case(cs, cpus = max_cpus), ord, info = nm)
    if (at_home()) expect_identical(fit_case(cs, cpus = 3), ord, info = nm)
  }

  # A larger, sparse table, where the products run in several blocks.
  set.seed(7)
  nb   <- 300
  bigc <- matrix(rpois(nb * 150, 2) * rbinom(nb * 150, 1, 0.2), nb, 150,
                 dimnames = list(paste0('S', 1:nb), paste0('F', 1:150)))
  bigm <- data.frame(row.names = rownames(bigc), g = rep(letters[1:3], each = 100), x = rnorm(nb))
  ob   <- rda(bigc, data = bigm, formula = ~ g + Condition(x), norm = 'hellinger', n_axes = 6, cpus = 1)
  expect_identical(colnames(ob), c('RDA1', 'RDA2', paste0('PC', 1:4)))
  expect_identical(rda(bigc, data = bigm, formula = ~ g + Condition(x), norm = 'hellinger', n_axes = 6, cpus = max_cpus), ob)

  # Unconstrained axes are orthogonal to the constraints, the conditions and
  # the ones vector.
  des <- cbind(1, model.matrix(~ g + x, bigm)[, -1])
  des <- des / rep(sqrt(colSums(des ^ 2)), each = nrow(des))
  pcs <- unclass(ob)[, 3:6]
  pcs <- pcs / rep(sqrt(colSums(pcs ^ 2)), each = nrow(pcs))
  expect_true(max(abs(crossprod(des, pcs))) < 1e-10)



  # Aliased (redundant) constraints are dropped with a message. ====

  cs  <- case('dune, hellinger, factor + numeric')
  dat <- cbind(cs$data, A1x2 = 2 * cs$data$A1)
  expect_message(o2 <- rda(cs$counts, data = dat, formula = ~ Management + A1 + A1x2, norm = 'hellinger', cpus = 1), 'linear combinations')
  o1 <- fit_case(cs)
  expect_equal(unclass(o2)[, ], unclass(o1)[, ])
  expect_identical(rownames(attr(o2, 'biplot')), rownames(attr(o1, 'biplot')))



  # Edge cases. ====

  # A constant feature and an all-zero feature: loadings of 0 with scale.
  cc <- cbind(counts, Const = 3, Zero = 0)
  os <- rda(cc, 'Management', meta, scale = TRUE, cpus = 1)
  expect_true(all(attr(os, 'loadings')[c('Const', 'Zero'), ] == 0))
  expect_equal(attr(os, 'trace'), ncol(counts))
  expect_equal(unclass(os)[, ], unclass(rda(counts, 'Management', meta, scale = TRUE, cpus = 1))[, ], tolerance = 1e-12)

  # A singleton level.
  m1 <- meta
  m1$Management <- as.character(m1$Management)
  m1$Management[1] <- 'ZZ'
  o1 <- rda(counts, 'Management', m1, norm = 'hellinger', cpus = 1)
  expect_identical(ncol(o1), 3L)
  expect_true('ManagementZZ' %in% rownames(attr(o1, 'centroids')))
  expect_true(all(is.finite(unclass(o1))))

  # Three samples.
  o3 <- suppressWarnings(rda(counts[1:3, ], c(1, 2, 4), norm = 'hellinger', n_axes = 3, cpus = 1))
  expect_identical(colnames(o3), c('RDA1', 'PC1'))
  expect_equal(attr(o3, 'trace'), sum(attr(o3, 'eig')))

  # Two samples: a saturated model.
  o2 <- suppressWarnings(rda(counts[1:2, ], c(1, 2), cpus = 1))
  expect_identical(colnames(o2), 'RDA1')
  expect_equal(attr(o2, 'r_squared'), 1)
  expect_true(all(attr(o2, 'inertia') >= 0))

  # One sample: no axes.
  o1 <- rda(counts[1, , drop = FALSE], cpus = 1)
  expect_identical(dim(o1), c(1L, 0L))
  expect_identical(dim(attr(o1, 'loadings')), c(ncol(counts), 0L))

  # One feature: one constrained axis (the fitted values) and one
  # unconstrained (the residuals).
  expect_warning(of <- rda(counts[, 1, drop = FALSE], 'A1', meta, cpus = 1), 'only 2 of the first 3')
  expect_identical(colnames(of), c('RDA1', 'PC1'))
  expect_equal(sum(attr(of, 'eig')), attr(of, 'trace'))

  # No variation at all: no axes, with a warning.
  expect_warning(oz <- rda(matrix(1, 5, 3), c(1, 2, 3, 4, 6), cpus = 1), 'only 0 of the first 3')
  expect_identical(ncol(oz), 0L)
  expect_identical(attr(oz, 'trace'), 0)

  # Constraints entirely within the conditions leave only PCs.
  expect_message(oc <- rda(counts, data = meta, formula = ~ A1 + Condition(A1), n_axes = 2, cpus = 1), 'linear combinations')
  expect_identical(colnames(oc), c('PC1', 'PC2'))
  expect_equal(attr(oc, 'inertia')[['constrained']], 0)

  # Unnamed samples and features.
  ou <- rda(unname(counts), meta$Management, cpus = 1)
  expect_null(rownames(ou))
  expect_null(rownames(attr(ou, 'loadings')))

#})
