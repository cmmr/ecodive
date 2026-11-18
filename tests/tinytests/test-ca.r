#test_that("ca", {

  # Inputs and reference results, written by data-raw/expected/ca.r.
  ref <- readRDS(file.path('expected', 'ca.rds'))

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  case <- function (name) {
    for (cs in ref$cases) if (cs$name == name) return (cs)
    stop('no case ', name)
  }

  # Flips the columns of `x` to match `y`'s signs.
  flip <- function (x, y) sweep(x, 2, sign(colSums(x * y)), '*')



  # Same as vegan::ca() (vegan::cca() without constraints), ade4::dudi.coa()
  # and scikit-bio's ca(), up to the sign of each axis: stored results, so
  # at home only. ====

  if (at_home()) {

    for (cs in ref$cases) {

      exp <- cs$expected
      k   <- ncol(exp$sites)
      lbl <- cs$name

      if (is.null(cs$warn)) {
        ord <- ca(cs$counts, n_axes = cs$k, cpus = 1)
      } else {
        expect_warning(ord <- ca(cs$counts, n_axes = cs$k, cpus = 1), cs$warn, info = lbl)
      }

      expect_identical(dim(ord), c(nrow(cs$counts), k), info = lbl)
      expect_identical(dimnames(ord), dimnames(exp$sites), info = lbl)

      s <- sign(colSums(unclass(ord) * exp$sites))
      expect_equal(sweep(unclass(ord), 2, s, '*'), exp$sites, check.attributes = FALSE, tolerance = 1e-10, info = lbl)

      ld <- attr(ord, 'loadings')
      expect_identical(dimnames(ld), dimnames(exp$loadings), info = lbl)
      expect_equal(sweep(ld, 2, s, '*'), exp$loadings, tolerance = 1e-10, info = lbl)

      expect_equal(attr(ord, 'eig'),   exp$eig,   tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'prop'),  exp$prop,  tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'trace'), exp$trace, tolerance = 1e-10, info = lbl)

      # ade4: every eigenvalue (their sum is the trace), $l1 is ca()'s
      # scores, and $co its loadings (of the features with observations).
      a <- cs$ade4
      expect_equal(attr(ord, 'eig'),   a$eig[seq_len(k)], tolerance = 1e-10, info = lbl)
      expect_equal(attr(ord, 'trace'), sum(a$eig),        tolerance = 1e-10, info = lbl)
      expect_equal(flip(unclass(ord)[, ], a$l1), a$l1, check.attributes = FALSE, tolerance = 1e-10, info = lbl)
      seen <- colSums(cs$counts) > 0
      expect_equal(flip(ld[seen, , drop = FALSE], a$co), a$co, check.attributes = FALSE, tolerance = 1e-10, info = lbl)
    }


    # scikit-bio's ca() in scaling 2, every axis.
    for (nm in names(ref$skbio)) {
      sk  <- ref$skbio[[nm]]
      m   <- ncol(sk$sites)
      ord <- ca(sk$counts, n_axes = m, cpus = 1)
      expect_identical(ncol(ord), m, info = nm)
      s <- sign(colSums(unclass(ord) * sk$sites))
      expect_equal(sweep(unclass(ord), 2, s, '*'),           sk$sites,   check.attributes = FALSE, tolerance = 1e-10, info = nm)
      expect_equal(sweep(attr(ord, 'loadings'), 2, s, '*'), sk$species, check.attributes = FALSE, tolerance = 1e-10, info = nm)
      expect_equal(attr(ord, 'eig'), unname(sk$eigvals), tolerance = 1e-10, info = nm)
    }
  }



  # Result format. ====

  x   <- case('dune')$counts
  ord <- ca(x, cpus = 1)
  expect_true(is.matrix(ord) && is.double(ord))
  expect_identical(
    names(attributes(ord)),
    c('dim', 'dimnames', 'method', 'eig', 'prop', 'trace', 'loadings', 'class') )
  expect_identical(attr(ord, 'method'), 'ca')
  expect_identical(colnames(ord), c('CA1', 'CA2', 'CA3'))
  expect_identical(rownames(ord), rownames(x))
  expect_identical(dimnames(attr(ord, 'loadings')), list(colnames(x), colnames(ord)))
  expect_identical(class(ord), c('ecodive_ord', 'matrix', 'array'))
  expect_equal(attr(ord, 'prop'), attr(ord, 'eig') / attr(ord, 'trace'))

  # The total inertia is Pearson's chi-square over the grand total, and the
  # eigenvalues are the squared singular values of the chi-square residuals.
  e  <- outer(rowSums(x), colSums(x)) / sum(x)
  expect_equal(attr(ord, 'trace'), sum((x - e) ^ 2 / e) / sum(x), tolerance = 1e-12)
  qb <- (x - e) / sqrt(e) / sqrt(sum(x))
  sv <- svd(qb)
  oa <- ca(x, n_axes = 19, cpus = 1)
  expect_equal(attr(oa, 'eig'), sv$d[1:19] ^ 2, tolerance = 1e-10)
  expect_equal(sum(attr(oa, 'eig')), attr(oa, 'trace'), tolerance = 1e-12)
  expect_equal(sum(attr(oa, 'prop')), 1, tolerance = 1e-12)

  # Scaling 2: site scores D_r^-1/2 U, loadings D_c^-1/2 V Sigma.
  r <- rowSums(x) / sum(x)
  cw <- colSums(x) / sum(x)
  expect_equal(flip(unclass(ord)[, ], sv$u[, 1:3] / sqrt(r)), sv$u[, 1:3] / sqrt(r),
               check.attributes = FALSE, tolerance = 1e-10)
  vs <- sweep(sv$v[, 1:3] / sqrt(cw), 2, sv$d[1:3], '*')
  expect_equal(flip(attr(ord, 'loadings'), vs), vs, check.attributes = FALSE, tolerance = 1e-10)

  # The site scores have weighted mean 0 and weighted variance 1, and are
  # the weighted averages of the species scores over sqrt(eig) squared
  # (ter Braak's transition formula).
  expect_equal(colSums(unclass(ord) * r), rep(0, 3), check.attributes = FALSE, tolerance = 1e-12)
  expect_equal(colSums(unclass(ord) ^ 2 * r), rep(1, 3), check.attributes = FALSE, tolerance = 1e-10)
  wa <- ((x / rowSums(x)) %*% attr(ord, 'loadings')) / rep(attr(ord, 'eig'), each = nrow(x))
  expect_equal(unclass(ord)[, ], wa, check.attributes = FALSE, tolerance = 1e-10)

  # Each axis's largest score is positive.
  for (j in 1:3) expect_true(unclass(ord)[which.max(abs(ord[, j])), j] > 0)



  # The same as cca() without constraints. ====

  for (nm in c('dune', 'varespec', 'n = 110')) {
    cs <- case(nm)
    o1 <- ca(cs$counts, n_axes = cs$k, cpus = 1)
    o2 <- cca(cs$counts, n_axes = cs$k, cpus = 1)
    expect_identical(unclass(o1)[, ], unclass(o2)[, ], info = nm)
    expect_identical(attr(o1, 'eig'),      attr(o2, 'eig'),      info = nm)
    expect_identical(attr(o1, 'prop'),     attr(o2, 'prop'),     info = nm)
    expect_identical(attr(o1, 'trace'),    attr(o2, 'trace'),    info = nm)
    expect_identical(attr(o1, 'loadings'), attr(o2, 'loadings'), info = nm)
  }



  # Ways to call it. ====

  counts <- x
  expect_identical(count_ord(counts, 'ca', cpus = 1), ord)
  expect_identical(count_ord(counts, 'CA', cpus = 1), ord)
  expect_identical(count_ord(counts, 'correspondence analysis', cpus = 1), ord)
  expect_identical(count_ord(counts, 'ca', seed = 3L, cpus = 1), ord) # seed is dropped
  expect_identical(count_ord(counts, 'ca', n_axes = 2, cpus = 1), ca(counts, n_axes = 2, cpus = 1))
  expect_error(dist_ord(bray(counts), 'ca'), 'count_ord')

  # Arguments and the table are validated.
  expect_error(ca(counts, n_axes = 0),   '`n_axes`')
  expect_error(ca(counts, n_axes = 1.5), '`n_axes`')
  expect_error(ca(counts, cpus = 0),     '`cpus`')
  expect_error(ca(counts, margin = 3),   '`margin`')
  expect_error(ca(counts, groups = 'Management'), 'unused argument')
  bc <- counts + 0
  bc[2, 2] <- NA
  expect_error(ca(bc), 'NA')
  bc[2, 2] <- -1
  expect_error(ca(bc), 'negative')
  bc[2, 2] <- Inf
  expect_error(ca(bc), 'infinite')

  # Samples with no observations are an error, naming them.
  bc <- counts + 0
  bc[c(3, 5), ] <- 0
  expect_error(ca(bc), '2 samples with no observations \\(all zero\\), which CA cannot place: 3, 5')
  rownames(bc) <- paste0('s', 1:20)
  expect_error(ca(bc), 's3, s5')
  expect_error(ca(unname(bc)), 'CA cannot place: 3, 5')
  expect_error(ca(bc[3, , drop = FALSE]), '1 sample with no observations')



  # Sparse and dense input, and either margin, give the same result. ====

  d1 <- ca(counts, n_axes = 5, cpus = 1)
  tc <- t(counts)
  expect_identical(ca(tc, margin = 2, n_axes = 5, cpus = 1), d1)
  if (requireNamespace('Matrix', quietly = TRUE)) {
    sp <- Matrix::Matrix(tc, sparse = TRUE)
    expect_true(inherits(sp, 'dgCMatrix'))
    expect_equal(ca(sp, margin = 2, n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
    expect_equal(ca(Matrix::t(sp), n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
  }
  if (requireNamespace('slam', quietly = TRUE)) {
    st <- slam::as.simple_triplet_matrix(counts)
    expect_equal(ca(st, n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
  }

  # Non-integer abundances work too, and scaling the table changes nothing.
  expect_equal(ca(counts / 7, n_axes = 5, cpus = 1), d1, tolerance = 1e-12)
  expect_equal(ca(as.data.frame(counts), n_axes = 5, cpus = 1), d1)



  # Deterministic, for any number of threads. ====

  for (nm in c('n = 110', 'varespec')) {
    cs  <- case(nm)
    o   <- ca(cs$counts, n_axes = cs$k, cpus = 1)
    expect_identical(ca(cs$counts, n_axes = cs$k, cpus = 1),        o, info = nm)
    expect_identical(ca(cs$counts, n_axes = cs$k, cpus = 2),        o, info = nm)
    expect_identical(ca(cs$counts, n_axes = cs$k, cpus = max_cpus), o, info = nm)
    if (at_home()) expect_identical(ca(cs$counts, n_axes = cs$k, cpus = 3), o, info = nm)
  }

  # A larger, sparse table, where the products run in several blocks.
  set.seed(7)
  nb   <- 300
  bigc <- matrix(rpois(nb * 150, 2) * rbinom(nb * 150, 1, 0.2), nb, 150,
                 dimnames = list(paste0('S', 1:nb), paste0('F', 1:150)))
  bigc[, 1] <- bigc[, 1] + 1
  ob <- ca(bigc, n_axes = 6, cpus = 1)
  expect_identical(ca(bigc, n_axes = 6, cpus = max_cpus), ob)
  if (at_home()) {
    eb <- outer(rowSums(bigc), colSums(bigc)) / sum(bigc)
    sb <- svd((bigc - eb) / sqrt(eb) / sqrt(sum(bigc)), nu = 0, nv = 0)$d
    expect_equal(attr(ob, 'eig'), sb[1:6] ^ 2, tolerance = 1e-10)
  }



  # Edge cases. ====

  # An all-zero feature has no weight: loadings of 0, nothing else changes.
  cz <- cbind(counts, Zero = 0)
  oz <- ca(cz, cpus = 1)
  expect_true(all(attr(oz, 'loadings')['Zero', ] == 0))
  expect_equal(unclass(oz)[, ], unclass(ord)[, ], tolerance = 1e-12)
  expect_equal(attr(oz, 'loadings')[colnames(counts), ], attr(ord, 'loadings'), tolerance = 1e-12)

  # More axes than n - 1: the n - 1 there are, without a warning; fewer
  # non-zero eigenvalues than that: a warning.
  expect_silent(o8 <- ca(counts[1:8, ], n_axes = 10, cpus = 1))
  expect_identical(colnames(o8), paste0('CA', 1:7))
  expect_warning(o4 <- ca(case('rank 3')$counts, n_axes = 8, cpus = 1), 'only 3 of the first 8')
  expect_identical(colnames(o4), paste0('CA', 1:3))

  # Three samples: two axes at most.
  o3 <- ca(counts[c(1, 7, 15), ], n_axes = 2, cpus = 1)
  expect_identical(colnames(o3), c('CA1', 'CA2'))
  expect_equal(attr(o3, 'trace'), sum(attr(o3, 'eig')), tolerance = 1e-12)
  expect_identical(colnames(ca(counts[c(1, 7, 15), ], cpus = 1)), c('CA1', 'CA2'))

  # Two samples: one axis.
  o2 <- ca(counts[1:2, ], cpus = 1)
  expect_identical(colnames(o2), 'CA1')
  expect_equal(attr(o2, 'eig'), attr(o2, 'trace'), tolerance = 1e-12)
  expect_equal(attr(o2, 'prop'), 1, tolerance = 1e-12)

  # One sample: no axes.
  o1 <- ca(counts[1, , drop = FALSE], cpus = 1)
  expect_identical(dim(o1), c(1L, 0L))
  expect_identical(dim(attr(o1, 'loadings')), c(ncol(counts), 0L))
  expect_identical(attr(o1, 'trace'), 0)

  # A single non-zero cell per sample, all in one feature: no inertia.
  expect_warning(os <- ca(matrix(c(1, 2, 3, 0, 0, 0), 3), cpus = 1), 'only 0 of the first 2')
  expect_identical(ncol(os), 0L)
  expect_identical(attr(os, 'trace'), 0)
  expect_identical(attr(os, 'prop'), numeric(0))

  # A single non-zero cell in the whole table: one sample, no axes; with
  # more samples, the others are empty.
  os <- ca(matrix(c(0, 5, 0), 1), cpus = 1)
  expect_identical(dim(os), c(1L, 0L))
  expect_error(ca(matrix(c(5, 0, 0, 0), 2), cpus = 1), '1 sample with no observations')

  # No variation at all: no axes, with a warning.
  expect_warning(oz <- ca(matrix(1, 5, 3), cpus = 1), 'only 0 of the first 3')
  expect_identical(ncol(oz), 0L)

  # Ties: a block-diagonal table has a repeated eigenvalue of 1 (the
  # table splits into separate groups); the eigenvalues are still exact.
  bd <- diag(3)[rep(1:3, each = 2), ] + 0
  ob <- ca(bd, n_axes = 2, cpus = 1)
  expect_equal(attr(ob, 'eig'), c(1, 1), tolerance = 1e-10)
  expect_equal(attr(ob, 'trace'), 2, tolerance = 1e-12)

  # Duplicate samples have identical scores.
  dd <- rbind(counts, counts[1, , drop = FALSE])
  od <- ca(dd, cpus = 1)
  expect_equal(unclass(od)[1, ], unclass(od)[21, ], tolerance = 1e-10)

  # Unnamed samples and features.
  ou <- ca(unname(counts), cpus = 1)
  expect_null(rownames(ou))
  expect_null(rownames(attr(ou, 'loadings')))
  expect_identical(colnames(ou), c('CA1', 'CA2', 'CA3'))

#})
