#test_that("assoc_pairs", {

  # assoc_pairs(): the pairs of a sparcc() or rho_p() result (or a
  # correlation matrix) at or above a cutoff, with p-values and their
  # adjustment; order, labels, clustering, and argument checks.

  set.seed(1)
  basis <- matrix(rnorm(60 * 30), 60, 30)
  basis[, 2] <- basis[, 1] + rnorm(60, sd = 0.3)
  basis[, 4] <- -basis[, 3] + rnorm(60, sd = 0.3)
  counts <- matrix(rpois(1800, exp(basis + 4)), 60, 30,
    dimnames = list(NULL, paste0('f', 1:30)) )
  top2 <- function (a) paste(a$from[1:2], a$to[1:2])

  r <- sparcc(counts, n_draws = 10L, cpus = 2L)
  m <- as.matrix(r)
  diag(m) <- 1

  # The pairs, strongest first, matching the matrix.
  a <- assoc_pairs(r)
  expect_identical(class(a), c('ecodive_assoc_pairs', 'data.frame'))
  expect_identical(names(a), c('from', 'to', 'association'))
  expect_identical(nrow(a), sum(abs(m[upper.tri(m)]) >= 0.3))
  expect_true(all(abs(a$association) >= 0.3))
  expect_false(is.unsorted(-abs(a$association)))
  expect_identical(a$association, m[cbind(a$from, a$to)])
  expect_true(all(match(a$from, colnames(counts)) < match(a$to, colnames(counts))))
  expect_identical(top2(a), c('f3 f4', 'f1 f2'))
  expect_identical(attr(a, 'method'), 'sparcc')
  expect_identical(attr(a, 'threshold'), 0.3)
  expect_null(attr(a, 'alpha'))
  expect_null(attr(a, 'p_adjust'))
  expect_identical(rownames(a), as.character(seq_len(nrow(a))))

  # threshold = 0 keeps every pair; a large one keeps none.
  expect_identical(nrow(assoc_pairs(r, threshold = 0)), 435L)
  none <- assoc_pairs(r, threshold = 2)
  expect_identical(nrow(none), 0L)
  expect_identical(names(none), c('from', 'to', 'association'))

  # The threshold is inclusive (|r| >= threshold).
  t2 <- abs(a$association[[2]])
  expect_identical(nrow(assoc_pairs(r, threshold = t2)), 2L)

  # Clustering doesn't change the pairs, nor which one is `from`.
  rc <- sparcc(counts, n_draws = 10L, clust = 'average', cpus = 2L)
  expect_false(identical(attr(rc, 'Labels'), colnames(counts)))
  expect_identical(assoc_pairs(rc), a)

  # Unnamed features: their input positions.
  u <- assoc_pairs(sparcc(unname(counts), n_draws = 10L, cpus = 2L))
  expect_identical(top2(u), c('3 4', '1 2'))
  u <- assoc_pairs(sparcc(unname(counts), n_draws = 10L, clust = 'average', cpus = 2L))
  expect_identical(top2(u), c('3 4', '1 2'))

  # Missing correlations are never kept.
  rn <- r
  rn[1] <- NA
  attributes(rn) <- attributes(r)
  expect_false(any(rn[1] %in% assoc_pairs(rn, threshold = 0)$association))
  expect_identical(nrow(assoc_pairs(rn, threshold = 0)), 434L)


  # -----------------------------------------------------------------------
  # p-values: adjusted over all pairs, then filtered by alpha.
  # -----------------------------------------------------------------------

  rp <- sparcc(counts, n_draws = 5L, permutations = 99L, cpus = 2L)
  p  <- attr(rp, 'p_value')
  a  <- assoc_pairs(rp, threshold = 0)
  expect_identical(names(a), c('from', 'to', 'association', 'p_value', 'p_adjusted'))
  mp <- as.matrix(rp); mp[] <- NA; mp[lower.tri(mp)] <- p
  mp[upper.tri(mp)] <- t(mp)[upper.tri(mp)]
  expect_identical(a$p_value, mp[cbind(a$from, a$to)])
  ma <- as.matrix(rp); ma[] <- NA; ma[lower.tri(ma)] <- p.adjust(p, 'BH')
  ma[upper.tri(ma)] <- t(ma)[upper.tri(ma)]
  expect_identical(a$p_adjusted, ma[cbind(a$from, a$to)])
  expect_identical(attr(a, 'p_adjust'), 'BH')

  # Adjusted over all 435 pairs, not just those kept.
  b <- assoc_pairs(rp, threshold = 0.5, p_adjust = 'bonferroni')
  expect_true(nrow(b) >= 2L)
  expect_identical(b$p_adjusted, pmin(1, b$p_value * 435))

  # alpha keeps the pairs at or below it, with the threshold too.
  for (al in c(0.05, 0.5, 1)) {
    k <- assoc_pairs(rp, alpha = al, p_adjust = 'none')
    expect_true(all(k$p_adjusted <= al & abs(k$association) >= 0.3), info = al)
    expect_identical(nrow(k), sum(p <= al & abs(rp) >= 0.3), info = al)
    expect_identical(attr(k, 'alpha'), al, info = al)
  }
  expect_identical(nrow(assoc_pairs(rp, threshold = 0, alpha = 1)), 435L)

  # Clustering keeps the p-values with their pairs.
  rpc <- sparcc(counts, n_draws = 5L, permutations = 99L, clust = 'average', cpus = 2L)
  expect_identical(assoc_pairs(rpc, threshold = 0), assoc_pairs(rp, threshold = 0))


  # -----------------------------------------------------------------------
  # Other inputs.
  # -----------------------------------------------------------------------

  # rho_p().
  x   <- counts + 1L
  rho <- rho_p(x)
  a   <- assoc_pairs(rho, threshold = 0.2)
  mr  <- as.matrix(rho)
  expect_identical(a$association, mr[cbind(a$from, a$to)])
  expect_identical(attr(a, 'method'), 'rho_p')

  # A correlation matrix.
  cm <- stats::cor(log(x))
  a  <- assoc_pairs(cm)
  expect_identical(a$association, cm[cbind(a$from, a$to)])
  expect_identical(nrow(a), sum(abs(cm[upper.tri(cm)]) >= 0.3))
  expect_null(attr(a, 'method'))
  expect_identical(assoc_pairs(unname(cm))[, 3], a[, 3])

  # Printing: a table with the method and cutoff in its header.
  out <- utils::capture.output(print(assoc_pairs(r)))
  expect_true(any(grepl('SparCC Correlation: pairs with |r| >= 0.3', out, fixed = TRUE)))
  out <- utils::capture.output(print(assoc_pairs(rp, alpha = 0.5)))
  expect_true(any(grepl('BH-adjusted p <= 0.5', out, fixed = TRUE)))
  out <- utils::capture.output(print(assoc_pairs(cm)))
  expect_true(any(grepl('Correlations: pairs with |r| >= 0.3', out, fixed = TRUE)))


  # -----------------------------------------------------------------------
  # Errors.
  # -----------------------------------------------------------------------

  expect_error(assoc_pairs(vlr(x)), 'vlr\\(\\), a dissimilarity')
  expect_error(assoc_pairs(phi_s(x)), 'phi_s\\(\\), a dissimilarity')
  expect_error(assoc_pairs(bray(x)), 'must be correlations')
  expect_error(assoc_pairs(stats::dist(x)), 'must be correlations')
  expect_error(assoc_pairs(x[1:30, ]), 'must be correlations')       # not symmetric
  expect_error(assoc_pairs(cm * 2), 'must be correlations')          # diagonal not 1
  expect_error(assoc_pairs(1:5), 'must be correlations')
  expect_error(assoc_pairs(r, alpha = 0.05), '`alpha` needs p-values')
  expect_error(assoc_pairs(r, threshold = -1), '`threshold`')
  expect_error(assoc_pairs(r, threshold = NA), '`threshold`')
  expect_error(assoc_pairs(r, threshold = c(0.1, 0.2)), '`threshold`')
  expect_error(assoc_pairs(rp, alpha = 0), '`alpha`')
  expect_error(assoc_pairs(rp, alpha = 1.5), '`alpha`')
  expect_error(assoc_pairs(rp, p_adjust = 'foo'), '`p_adjust`')

#})
