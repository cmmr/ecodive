#test_that("cca", {

  # Inputs and reference results, written by data-raw/expected/cca.r.
  ref <- readRDS(file.path('expected', 'cca.rds'))

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
      cca(counts, n_axes = cs$k, permutations = cs$permutations, cpus = cpus, ...)
    } else {
      cca(counts, data = cs$data, formula = f, n_axes = cs$k,
          permutations = cs$permutations, cpus = cpus, ...)
    })
  }

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }



  # Same as vegan::cca(), up to the sign of each axis: stored results, so at
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

      # vegan gives no R-squared without constraints; cca() gives 0. For
      # partial models, cca()'s adjusted R-squared is NA.
      if (length(exp$r_squared) == 0) {
        expect_identical(attr(ord, 'r_squared'),     0, info = lbl)
        expect_identical(attr(ord, 'adj_r_squared'), 0, info = lbl)
      } else {
        expect_equal(attr(ord, 'r_squared'), exp$r_squared, tolerance = 1e-10, info = lbl)
        if (length(exp$adj_r_squared) == 0) { # No residual inertia.
          expect_identical(attr(ord, 'adj_r_squared'), NA_real_, info = lbl)
        } else if (is.null(cs$conditions)) {
          expect_equal(attr(ord, 'adj_r_squared'), exp$adj_r_squared, tolerance = 1e-10, info = lbl)
        } else {
          expect_identical(attr(ord, 'adj_r_squared'), NA_real_, info = lbl)
        }
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


    # scikit-bio's cca() in scaling 2, every axis: sites ("wa"), species and
    # eigenvalues. Its biplot scores are unweighted correlations of the
    # weighted constraints with the axes, where vegan's (and cca()'s) are
    # weighted correlations, so only their signs' pattern is compared.
    sk  <- ref$skbio
    m   <- ncol(sk$sites)
    ord <- cca(sk$counts, data = sk$x, formula = ~ ManagementHF + ManagementNM + ManagementSF + A1,
               n_axes = m, permutations = 0, cpus = 1)
    expect_identical(ncol(ord), m)
    s <- sign(colSums(unclass(ord) * sk$sites))
    expect_equal(sweep(unclass(ord), 2, s, '*'),           sk$sites,   check.attributes = FALSE, tolerance = 1e-10)
    expect_equal(sweep(attr(ord, 'loadings'), 2, s, '*'), sk$species, check.attributes = FALSE, tolerance = 1e-10)
    expect_equal(attr(ord, 'eig'), unname(sk$eigvals), tolerance = 1e-10)
  }



  # Result format. ====

  cs  <- case('dune, factor + numeric')
  ord <- fit_case(cs)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'inertia',
      'constrained_eig', 'r_squared', 'adj_r_squared', 'biplot', 'centroids',
      'loadings', 'class') )
  expect_identical(attr(ord, 'method'), 'cca')
  expect_identical(colnames(ord), c('CCA1', 'CCA2', 'CCA3'))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))

  # Inertia components add up, and R-squared is the constrained share.
  expect_equal(attr(ord, 'trace'), sum(attr(ord, 'inertia')))
  expect_equal(attr(ord, 'r_squared'), attr(ord, 'inertia')[['constrained']] / attr(ord, 'trace'))

  # The total inertia is Pearson's chi-square over the grand total.
  x  <- cs$counts
  e  <- outer(rowSums(x), colSums(x)) / sum(x)
  expect_equal(attr(ord, 'trace'), sum((x - e) ^ 2 / e) / sum(x), tolerance = 1e-12)

  # Constrained axes, then CA axes.
  cs6 <- case('dune, k > constrained axes')
  ord6 <- fit_case(cs6)
  expect_identical(colnames(ord6), c('CCA1', 'CCA2', 'CA1', 'CA2', 'CA3', 'CA4'))
  expect_identical(colnames(attr(ord6, 'loadings')), colnames(ord6))
  expect_true(all(attr(ord6, 'biplot')[, 3:6] == 0))

  # "wa" scores: in scaling 2, each site score is the weighted average of
  # the species scores divided by the axis's eigenvalue (ter Braak's
  # transition formula), for constrained and CA axes.
  # (CA axes follow it for the residual table instead.)
  for (o in list(ord, ord6)) {
    j  <- grep('^CCA', colnames(o))
    wa <- ((x / rowSums(x)) %*% attr(o, 'loadings')[, j]) / rep(attr(o, 'eig')[j], each = nrow(x))
    expect_equal(unclass(o)[, j], wa, check.attributes = FALSE, tolerance = 1e-10)
  }



  # The adjusted R-squared, by permutation. ====

  meta <- cs$data
  a1 <- attr(cca(x, 'A1', meta, cpus = 1), 'adj_r_squared')
  expect_true(is.finite(a1) && a1 < attr(cca(x, 'A1', meta, cpus = 1), 'r_squared'))
  expect_identical(attr(cca(x, 'A1', meta, seed = 0L, cpus = 1), 'adj_r_squared'), a1)
  expect_false(identical(attr(cca(x, 'A1', meta, seed = 1L, cpus = 1), 'adj_r_squared'), a1))
  expect_identical(attr(cca(x, 'A1', meta, permutations = 0, cpus = 1), 'adj_r_squared'), NA_real_)
  expect_identical(attr(cca(x, data = meta, formula = ~ A1 + Condition(Use), cpus = 1), 'adj_r_squared'), NA_real_)
  expect_identical(attr(cca(x, cpus = 1), 'adj_r_squared'), 0)

  # From a supplied matrix: 1 - (1 - R2) / (1 - mean permuted R2), where a
  # permutation moves the samples' rows of the table against the
  # constraints.
  pm  <- case('dune, factor')$permutations[1:10, ]
  r2p <- vapply(seq_len(nrow(pm)), function (i)
    attr(cca(x[pm[i, ], ], meta$Management, permutations = 0, cpus = 1), 'r_squared'), 1)
  o   <- cca(x, 'Management', meta, permutations = pm, cpus = 1)
  expect_equal(attr(o, 'adj_r_squared'), 1 - (1 - attr(o, 'r_squared')) / (1 - mean(r2p)), tolerance = 1e-12)

  expect_error(cca(x, 'A1', meta, permutations = -1),             '`permutations`')
  expect_error(cca(x, 'A1', meta, permutations = matrix(1L, 2, 3)), 'columns')
  expect_error(cca(x, 'A1', meta, seed = 'a'),                    '`seed`')



  # Ways to give the model, as rda(). ====

  counts <- x
  one    <- cca(counts, data = meta, formula = ~ Management, cpus = 1)

  expect_identical(cca(counts, 'Management', meta, cpus = 1), one)
  expect_equal(unclass(cca(counts, meta$Management, cpus = 1))[, ], unclass(one)[, ])

  # vegan's style: the table on the formula's left-hand side.
  expect_identical(cca(counts ~ Management, meta, cpus = 1), one)
  expect_identical(cca(counts ~ Management, data = meta, cpus = 1), one)
  expect_error(cca(counts ~ A1, meta, ~ A1), 'Name the formula argument')
  expect_error(cca(counts ~ A1, data = meta, formula = ~ A1), 'not both')

  # count_ord() passes the same arguments on, including `seed`.
  expect_identical(count_ord(counts, 'cca', 'Management', meta, cpus = 1), one)
  expect_identical(count_ord(counts, 'canonical correspondence analysis', data = meta, formula = ~ Management, cpus = 1), one)
  expect_identical(
    count_ord(counts, 'cca', 'Management', meta, seed = 3L, cpus = 1),
    cca(counts, 'Management', meta, seed = 3L, cpus = 1) )
  expect_error(dist_ord(bray(counts), 'cca'), 'count_ord')

  # Rows of `data` are matched to the samples by name.
  expect_identical(cca(counts, 'Management', meta[rev(seq_len(nrow(meta))), ], cpus = 1), one)

  # Errors in the model.
  expect_error(cca(counts, meta$Managment, meta),         '`groups` is NULL')
  expect_error(cca(counts, formula = ~ Management),       '`data` is required')
  expect_error(cca(counts, meta, ~ A1),                   'Name the formula argument')
  expect_error(cca(counts, 'Management', meta, formula = ~ A1), 'not both')
  bad <- meta
  bad$A1[3] <- NA
  expect_error(cca(counts, 'A1', bad), 'NA')

  # Other arguments are validated, and so is the table.
  expect_error(cca(counts, 'A1', meta, n_axes = 0), '`n_axes`')
  expect_error(cca(counts, 'A1', meta, cpus = 0),   '`cpus`')
  bc <- counts + 0
  bc[2, 2] <- NA
  expect_error(cca(bc, 'A1', meta), 'NA')
  bc[2, 2] <- -1
  expect_error(cca(bc, 'A1', meta), 'negative')

  # Samples with no observations are an error, naming them.
  bc <- counts + 0
  bc[c(3, 5), ] <- 0
  expect_error(cca(bc, 'A1', meta), '2 samples with no observations \\(all zero\\), which CCA cannot place: 3, 5')
  rownames(bc) <- paste0('s', 1:20)
  expect_error(cca(bc), 's3, s5')
  expect_error(cca(bc[3, , drop = FALSE]), 'no observations')



  # Unconstrained: a correspondence analysis. Its eigenvalues are the
  # squared singular values of the chi-square residuals. ====

  qb <- (x / sum(x) - outer(rowSums(x), colSums(x)) / sum(x) ^ 2) /
        sqrt(outer(rowSums(x), colSums(x)) / sum(x) ^ 2)
  sv <- svd(qb)$d ^ 2
  ou <- cca(x, n_axes = 6, cpus = 1)
  expect_identical(colnames(ou), paste0('CA', 1:6))
  expect_equal(attr(ou, 'eig'), sv[1:6], tolerance = 1e-10)
  expect_equal(attr(ou, 'trace'), sum(sv), tolerance = 1e-12)
  expect_equal(attr(ou, 'inertia')[['unconstrained']], attr(ou, 'trace'))
  expect_identical(attr(ou, 'r_squared'), 0)
  expect_null(attr(ou, 'centroids'))

  # The site scores have weighted mean 0 and weighted variance 1.
  r <- rowSums(x) / sum(x)
  expect_equal(colSums(unclass(ou) * r), rep(0, 6), check.attributes = FALSE, tolerance = 1e-12)
  expect_equal(colSums(unclass(ou) ^ 2 * r), rep(1, 6), check.attributes = FALSE, tolerance = 1e-10)

  # Conditions only: the CA of what they leave.
  oc <- cca(x, data = meta, formula = ~ Condition(A1), n_axes = 3, cpus = 1)
  expect_identical(colnames(oc), paste0('CA', 1:3))
  expect_true(attr(oc, 'inertia')[['conditional']] > 0)



  # Sparse and dense input, and either margin, give the same result. ====

  d1 <- cca(counts, data = meta, formula = ~ Management + Condition(A1), n_axes = 5, cpus = 1)
  tc <- t(counts)
  expect_identical(cca(tc, data = meta, formula = ~ Management + Condition(A1), margin = 2,
                       n_axes = 5, cpus = 1), d1)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(tc, sparse = TRUE)
    expect_true(inherits(sp, 'dgCMatrix'))
    expect_equal(cca(sp, data = meta, formula = ~ Management + Condition(A1), margin = 2,
                     n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(counts)
    expect_equal(cca(st, data = meta, formula = ~ Management + Condition(A1),
                     n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
  }

  # Non-integer abundances work too, and scaling the table changes nothing.
  expect_equal(cca(counts / 7, data = meta, formula = ~ Management + Condition(A1), n_axes = 5, cpus = 1),
               d1, tolerance = 1e-12)



  # Deterministic, for any number of threads. ====

  for (nm in c('n = 110, partial', 'n = 110, 2 groups', 'varespec, 12 numeric')) {
    cs  <- case(nm)
    ord <- fit_case(cs, cpus = 1)
    expect_identical(fit_case(cs, cpus = 1),        ord, info = nm)
    expect_identical(fit_case(cs, cpus = 2),        ord, info = nm)
    expect_identical(fit_case(cs, cpus = max_cpus), ord, info = nm)
    if (at_home()) expect_identical(fit_case(cs, cpus = 3), ord, info = nm)
  }

  # A larger, sparse table, where the products run in several blocks, and
  # random permutations for the adjusted R-squared.
  set.seed(7)
  nb   <- 300
  bigc <- matrix(rpois(nb * 150, 2) * rbinom(nb * 150, 1, 0.2), nb, 150,
                 dimnames = list(paste0('S', 1:nb), paste0('F', 1:150)))
  bigc[, 1] <- bigc[, 1] + 1
  bigm <- data.frame(row.names = rownames(bigc), g = rep(letters[1:3], each = 100), x = rnorm(nb))
  ob   <- cca(bigc, data = bigm, formula = ~ g + Condition(x), n_axes = 6, cpus = 1)
  expect_identical(colnames(ob), c('CCA1', 'CCA2', paste0('CA', 1:4)))
  expect_identical(cca(bigc, data = bigm, formula = ~ g + Condition(x), n_axes = 6, cpus = max_cpus), ob)
  oa <- cca(bigc, data = bigm, formula = ~ g + x, n_axes = 3, permutations = 200, cpus = 1)
  expect_identical(cca(bigc, data = bigm, formula = ~ g + x, n_axes = 3, permutations = 200, cpus = max_cpus), oa)

  # Unconstrained axes are weighted-orthogonal to the constraints, the
  # conditions and the constant.
  rb  <- rowSums(bigc) / sum(bigc)
  des <- cbind(1, model.matrix(~ g + x, bigm)[, -1])
  ca_ <- unclass(ob)[, 3:6]
  wc  <- crossprod(des * rb, ca_) / sqrt(crossprod(des ^ 2, rb) %*% t(crossprod(ca_ ^ 2, rb)))
  expect_true(max(abs(wc)) < 1e-10)



  # Aliased (redundant) constraints are dropped with a message. ====

  cs  <- case('dune, factor + numeric')
  dat <- cbind(cs$data, A1x2 = 2 * cs$data$A1)
  expect_message(o2 <- cca(cs$counts, data = dat, formula = ~ Management + A1 + A1x2, cpus = 1), 'linear combinations')
  o1 <- cca(cs$counts, data = cs$data, formula = ~ Management + A1, cpus = 1)
  expect_equal(unclass(o2)[, ], unclass(o1)[, ])
  expect_equal(attr(o2, 'adj_r_squared'), attr(o1, 'adj_r_squared'))
  expect_identical(rownames(attr(o2, 'biplot')), rownames(attr(o1, 'biplot')))



  # Edge cases. ====

  # An all-zero feature has no weight: loadings of 0, nothing else changes.
  cz <- cbind(counts, Zero = 0)
  oz <- cca(cz, 'Management', meta, cpus = 1)
  expect_true(all(attr(oz, 'loadings')['Zero', ] == 0))
  expect_equal(unclass(oz)[, ], unclass(cca(counts, 'Management', meta, cpus = 1))[, ], tolerance = 1e-12)

  # A singleton level.
  m1 <- meta
  m1$Management <- as.character(m1$Management)
  m1$Management[1] <- 'ZZ'
  o1 <- cca(counts, 'Management', m1, cpus = 1)
  expect_identical(ncol(o1), 3L)
  expect_true('ManagementZZ' %in% rownames(attr(o1, 'centroids')))
  expect_true(all(is.finite(unclass(o1))))

  # Three samples.
  o3 <- suppressWarnings(cca(counts[1:3, ], c(1, 2, 4), n_axes = 3, cpus = 1))
  expect_identical(colnames(o3), c('CCA1', 'CA1'))
  expect_equal(attr(o3, 'trace'), sum(attr(o3, 'eig')))

  # Two samples: a saturated model.
  o2 <- suppressWarnings(cca(counts[1:2, ], c(1, 2), cpus = 1))
  expect_identical(colnames(o2), 'CCA1')
  expect_equal(attr(o2, 'r_squared'), 1)

  # One sample: no axes.
  o1 <- cca(counts[1, , drop = FALSE], cpus = 1)
  expect_identical(dim(o1), c(1L, 0L))
  expect_identical(dim(attr(o1, 'loadings')), c(ncol(counts), 0L))

  # A single non-zero cell per sample, all in one feature: no inertia.
  expect_warning(os <- cca(matrix(c(1, 2, 3, 0, 0, 0), 3), c(1, 2, 4), cpus = 1), 'only 0 of the first')
  expect_identical(ncol(os), 0L)
  expect_identical(attr(os, 'trace'), 0)

  # No variation at all: no axes, with a warning.
  expect_warning(oz <- cca(matrix(1, 5, 3), c(1, 2, 3, 4, 6), cpus = 1), 'only 0 of the first 3')
  expect_identical(ncol(oz), 0L)

  # Constraints entirely within the conditions leave only CA axes.
  expect_message(oc <- cca(counts, data = meta, formula = ~ A1 + Condition(A1), n_axes = 2, cpus = 1), 'linear combinations')
  expect_identical(colnames(oc), c('CA1', 'CA2'))
  expect_equal(attr(oc, 'inertia')[['constrained']], 0)

  # Unnamed samples and features.
  ou <- cca(unname(counts), meta$Management, cpus = 1)
  expect_null(rownames(ou))
  expect_null(rownames(attr(ou, 'loadings')))

#})
