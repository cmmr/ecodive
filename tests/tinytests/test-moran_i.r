#test_that("moran_i", {

  ex <- readRDS('expected/moran_i.rds')

  styles <- c(r = 'row', b = 'binary', v = 'variance', o = 'none')



  # Same statistics as esda's Moran (esda 2.10.0), for every case and
  # transformation, and the same permuted statistics given the same
  # permutations: stored results, so at home only. ====

  if (at_home()) {
    for (cs in names(ex)) {
      k <- ex[[cs]]
      for (i in seq_len(nrow(k$esda))) {

        e     <- k$esda[i, ]
        style <- styles[[e$style]]
        info  <- paste(cs, style)

        # esda's one-tailed analytic p-values are in the direction of z.
        alt <- if (is.finite(e$z_norm) && e$z_norm < 0) 'less' else 'greater'

        res <- suppressWarnings(moran_i(
          k$x, k$w, style = style, alternative = alt,
          permutations = k$perms, cpus = 1 ))

        expect_equal(res$statistic[['I']], e$I, tolerance = 1e-12, info = info)
        expect_equal(res$expected, e$EI, tolerance = 1e-14, info = info)
        expect_equal(attr(res, 'i_perm'), k$sims[, e$style], tolerance = 1e-12, info = info)

        # With binary weights between every pair of the idw case's samples,
        # I is the same for every arrangement (esda's variances are 0 to
        # rounding error, its p-values 0): ecodive gives NA, and counts
        # every permuted I as a tie.
        if (cs == 'idw' && style == 'binary') {
          expect_identical(res$sd_norm, NA_real_, info = info)
          expect_identical(res$sd_rand, NA_real_, info = info)
          expect_identical(res$p_value_norm, NA_real_, info = info)
          expect_identical(res$p_value, 1, info = info)
          next
        }

        expect_equal(res$sd_norm, sqrt(e$VI_norm), tolerance = 1e-12, info = info)
        expect_equal(res$sd_rand, sqrt(e$VI_rand), tolerance = 1e-12, info = info)
        expect_equal(res$p_value_norm, e$p_norm, tolerance = 1e-10, info = info)
        expect_equal(res$p_value_rand, e$p_rand, tolerance = 1e-10, info = info)

        # esda's p_sim is from the smaller tail, as ecodive's p-value from
        # that side, unless a permuted I ties with the observed one (esda
        # counts ties in the upper tail only, ecodive within rounding).
        ip <- attr(res, 'i_perm')
        i0 <- res$statistic[['I']]
        if (!any(abs(ip - i0) < 1e-10)) {
          side <- if (sum(ip >= i0) <= length(ip) / 2) 'greater' else 'less'
          p    <- suppressWarnings(moran_i(
            k$x, k$w, style = style, alternative = side, permutations = k$perms, cpus = 1 ))$p_value
          expect_identical(p, e$p_sim, info = info)
        }

        # Two-tailed analytic p-values.
        res2 <- suppressWarnings(moran_i(
          k$x, k$w, style = style, alternative = 'two.sided', permutations = 0 ))
        expect_equal(res2$p_value_norm, e$p_norm2, tolerance = 1e-10, info = info)
        expect_equal(res2$p_value_rand, e$p_rand2, tolerance = 1e-10, info = info)

        # libpysal's sums of the transformed weights.
        sp <- ecodive:::moran_space(k$w, NULL)
        wt <- suppressWarnings(ecodive:::moran_weights(sp, style))
        expect_equal(c(wt$s0, wt$s1, wt$s2), c(e$s0, e$s1, e$s2), tolerance = 1e-12, info = info)
      }
    }
  }



  # ape's Moran.I() (ape 5.7.1): row-standardised weights, the standard
  # deviation under randomisation. spdep's moran.test() (spdep 1.3.1), for
  # each style, under randomisation and normality. Outputs only. ====

  if (at_home()) {
    sp_style <- c(W = 'row', B = 'binary', S = 'variance', U = 'none')
    for (cs in names(ex)) {
      k <- ex[[cs]]

      res <- suppressWarnings(moran_i(k$x, k$w, alternative = 'two.sided', permutations = 0))
      expect_equal(res$statistic[['I']], k$ape$observed, tolerance = 1e-12, info = cs)
      expect_equal(res$expected, k$ape$expected, tolerance = 1e-14, info = cs)
      expect_equal(res$sd_rand, k$ape$sd, tolerance = 1e-12, info = cs)
      expect_equal(res$p_value_rand, k$ape$p_two_sided, tolerance = 1e-10, info = cs)
      res <- suppressWarnings(moran_i(k$x, k$w, permutations = 0))
      expect_equal(res$p_value_rand, k$ape$p_greater, tolerance = 1e-10, info = cs)

      for (i in seq_len(nrow(k$spdep))) {
        s    <- k$spdep[i, ]
        info <- paste(cs, s$style, s$randomisation)
        res  <- suppressWarnings(moran_i(k$x, k$w, style = sp_style[[s$style]], permutations = 0))
        expect_equal(res$statistic[['I']], s$I, tolerance = 1e-12, info = info)
        expect_equal(res$expected, s$expectation, tolerance = 1e-14, info = info)
        if (cs == 'idw' && s$style == 'B') next  # variance 0 (see above)
        sd <- if (s$randomisation) res$sd_rand else res$sd_norm
        p  <- if (s$randomisation) res$p_value_rand else res$p_value_norm
        expect_equal(sd^2, s$variance, tolerance = 1e-12, info = info)
        expect_equal(p, s$p_greater, tolerance = 1e-10, info = info)
      }
    }
  }



  # Distances, coordinates and weights give the same results. ====

  k   <- ex$idw
  xy  <- k$xy
  rownames(xy) <- paste0('S', seq_len(nrow(xy)))
  x   <- setNames(k$x, rownames(xy))
  d   <- dist(xy)
  w   <- 1 / as.matrix(d)    # Inf on the diagonal, which is ignored.

  r_w <- moran_i(x, w, permutations = 99, cpus = 1)
  r_d <- moran_i(x, d, permutations = 99, cpus = 1)
  r_c <- moran_i(x, xy, permutations = 99, cpus = 1)
  r_f <- moran_i(x, as.data.frame(xy), permutations = 99, cpus = 1)
  meta <- data.frame(row.names = rownames(xy), east = xy[, 1], north = xy[, 2], val = k$x)
  r_m <- moran_i('val', c('east', 'north'), meta, permutations = 99, cpus = 1)

  for (r in list(r_d, r_c, r_f, r_m)) {
    expect_equal(r$statistic, r_w$statistic, tolerance = 1e-13)
    expect_equal(r$sd_norm, r_w$sd_norm, tolerance = 1e-13)
    expect_equal(r$sd_rand, r_w$sd_rand, tolerance = 1e-13)
    expect_equal(attr(r, 'i_perm'), attr(r_w, 'i_perm'), tolerance = 1e-13)
    expect_identical(r$p_value, r_w$p_value)
  }
  expect_identical(r_w$space_type, 'weights')
  expect_identical(r_d$space_type, 'distances')
  expect_identical(r_c$space_type, 'coordinates')
  expect_identical(r_m$space_type, 'coordinates')

  # As esda with the matrix from the fixture.
  if (at_home())
    expect_equal(r_w$statistic[['I']], k$esda$I[k$esda$style == 'r'], tolerance = 1e-12)



  # Samples are matched by name. ====

  o  <- c(5:30, 1:4)
  r1 <- moran_i(x[o], d, permutations = 0)
  expect_equal(r1$statistic, r_w$statistic, tolerance = 1e-13)
  r2 <- moran_i('val', d, meta[rev(seq_len(30)), ], permutations = 0)
  expect_equal(r2$statistic, r_w$statistic, tolerance = 1e-13)
  r3 <- moran_i(k$x, w[o, o], permutations = 0)   # unnamed x: in w's order
  expect_false(isTRUE(all.equal(r3$statistic, r_w$statistic)))
  r4 <- moran_i(x, w[o, o], permutations = 0)
  expect_equal(r4$statistic, r_w$statistic, tolerance = 1e-13)
  # data's rows for other samples are ignored.
  meta2 <- rbind(meta, data.frame(row.names = 'extra', east = 0, north = 0, val = NA))
  expect_equal(moran_i('val', d, meta2, permutations = 0)$statistic, r_w$statistic, tolerance = 1e-13)



  # Result format. ====

  res <- moran_i(x, d, permutations = 99, cpus = 1)

  expect_identical(class(res), c('ecodive_moran_i', 'ecodive_test'))
  expect_identical(names(res), c(
    'statistic', 'p_value', 'z', 'n_perm',
    'expected', 'sd_norm', 'p_value_norm', 'sd_rand', 'p_value_rand',
    'style', 'alternative', 'n', 'space_type', 'permutation' ))
  expect_identical(names(res$statistic), 'I')
  expect_identical(res$style, 'row')
  expect_identical(res$alternative, 'greater')
  expect_identical(res$n, 30L)
  expect_identical(res$n_perm, 99L)
  expect_identical(res$permutation, 'free')
  expect_identical(res$expected, -1 / 29)
  expect_identical(attr(res, 'display')$title, "Moran's I")
  expect_identical(length(attr(res, 'i_perm')), 99L)
  expect_identical(res$p_value, ecodive:::perm_p(res$statistic[[1]], attr(res, 'i_perm')))
  expect_identical(res$z, ecodive:::perm_z(res$statistic[[1]], attr(res, 'i_perm')))

  out <- capture.output(print(res))
  expect_identical(out[[1]], "Moran's I")
  expect_true(any(grepl('Expected I: ', out, fixed = TRUE)))
  expect_true(any(grepl('Under randomisation: sd = ', out, fixed = TRUE)))
  expect_true(any(grepl('Weights: row-standardised inverse distances', out, fixed = TRUE)))
  expect_true(any(grepl('Samples: 30', out, fixed = TRUE)))

  row <- as.data.frame(res)
  expect_identical(row$test, 'moran_i')
  expect_identical(row$stat_name, 'I')
  expect_identical(row$p_value, res$p_value)

  pdf(NULL)
  expect_identical(plot(res), res)
  dev.off()

  res0 <- moran_i(x, d, permutations = 0)
  expect_identical(res0$p_value, NA_real_)
  expect_identical(res0$z, NA_real_)
  expect_identical(attr(res0, 'i_perm'), numeric(0))
  expect_identical(res0$statistic, res$statistic)

  # The analytic p-values for each alternative.
  zr <- (res$statistic[[1]] - res$expected) / res$sd_rand
  expect_equal(moran_i(x, d, alternative = 'less', permutations = 0)$p_value_rand, pnorm(zr))
  expect_equal(res$p_value_rand, pnorm(zr, lower.tail = FALSE))
  expect_equal(
    moran_i(x, d, alternative = 'two.sided', permutations = 0)$p_value_rand,
    2 * pnorm(-abs(zr)) )



  # Random permutations: the p-value counts permuted values at least as
  # extreme, for each alternative. ====

  for (alt in c('greater', 'less', 'two.sided')) {
    res  <- moran_i(x, d, alternative = alt, permutations = 199, seed = 3, cpus = 1)
    i    <- res$statistic[[1]]
    ip   <- attr(res, 'i_perm')
    tol  <- sqrt(.Machine$double.eps)
    hit  <- switch(alt, greater = ip >= i - tol, less = ip <= i + tol, two.sided = abs(ip) >= abs(i) - tol)
    expect_identical(res$p_value, (sum(hit) + 1) / 200, info = alt)
  }

  # The same seed gives the same permutations; another seed, others.
  a <- moran_i(x, d, permutations = 99, seed = 1, cpus = 1)
  b <- moran_i(x, d, permutations = 99, seed = 1, cpus = 1)
  c <- moran_i(x, d, permutations = 99, seed = 2, cpus = 1)
  expect_identical(a, b)
  expect_false(identical(attr(a, 'i_perm'), attr(c, 'i_perm')))

  # A supplied permutation matrix: sample i takes the value of sample p[i].
  set.seed(1)
  pm  <- t(replicate(20, sample(30)))
  res <- moran_i(x, d, permutations = pm, cpus = 1)
  expect_identical(res$permutation, 'supplied matrix')
  for (j in c(1, 7, 20))
    expect_equal(attr(res, 'i_perm')[[j]], moran_i(unname(x)[pm[j, ]], d, permutations = 0)$statistic[[1]], tolerance = 1e-13)

  # Negative autocorrelation.
  if (at_home()) {
    k   <- ex$checker
    res <- moran_i(k$x, k$w, alternative = 'less', permutations = k$perms)
    expect_true(res$statistic[[1]] < res$expected)
    expect_identical(res$p_value, 1 / 200)
  }



  # Identical for any number of threads. ====

  cpus <- if (at_home()) c(1, 2, 3, 4, 8) else c(1, 2)
  set.seed(2)
  big  <- matrix(runif(400 * 2), ncol = 2)
  bx   <- big[, 1] + rnorm(400)
  ref  <- moran_i(bx, big, permutations = 999, seed = 5, cpus = 1)
  for (cp in cpus[-1])
    expect_identical(moran_i(bx, big, permutations = 999, seed = 5, cpus = cp), ref, info = cp)



  # Edge cases. ====

  # A constant x: NA, with a warning.
  expect_warning(
    res <- moran_i(rep(2, 30), d, permutations = 9),
    'same value for every sample' )
  expect_identical(res$statistic[['I']], NA_real_)
  expect_identical(res$p_value, NA_real_)
  expect_identical(res$p_value_norm, NA_real_)
  expect_identical(res$sd_rand, NA_real_)
  expect_identical(attr(res, 'i_perm'), rep(NA_real_, 9))

  # Islands: kept, with a warning; they count in n.
  if (at_home()) {
    k <- ex$island
    expect_warning(res <- moran_i(k$x, k$w, permutations = 0), '1 sample has no neighbours')
    expect_identical(res$n, 12L)
  }
  w0 <- matrix(0, 5, 5)
  w0[cbind(c(1, 2, 2, 3, 3, 4), c(2, 1, 3, 2, 4, 3))] <- 1   # sample 5 has none
  expect_warning(moran_i(1:5, w0, permutations = 0), '1 sample has no neighbours')

  # Three samples: no variance under randomisation.
  res <- moran_i(c(1, 3, 2), dist(c(0, 1, 3)), permutations = 0)
  expect_true(is.finite(res$statistic[[1]]))
  expect_true(is.finite(res$sd_norm))
  expect_identical(res$sd_rand, NA_real_)
  expect_identical(res$p_value_rand, NA_real_)

  # Asymmetric weights: the same as their symmetric part with style 'none'.
  set.seed(3)
  wa <- matrix(runif(36), 6); diag(wa) <- 0
  xa <- rnorm(6)
  expect_equal(
    moran_i(xa, wa, style = 'none', permutations = 0)$statistic,
    moran_i(xa, (wa + t(wa)) / 2, style = 'none', permutations = 0)$statistic )
  # The direct formula.
  z  <- xa - mean(xa)
  expect_equal(
    moran_i(xa, wa, style = 'none', permutations = 0)$statistic[[1]],
    6 / sum(wa) * sum(wa * outer(z, z)) / sum(z^2) )
  wr <- wa / rowSums(wa)
  expect_equal(
    moran_i(xa, wa, permutations = 0)$statistic[[1]],
    6 / sum(wr) * sum(wr * outer(z, z)) / sum(z^2) )
  # Logical weights are binary.
  expect_equal(
    moran_i(xa, wa > 0.5, style = 'none', permutations = 0)$statistic,
    moran_i(xa, (wa > 0.5) * 1, style = 'binary', permutations = 0)$statistic )


  # Errors.
  expect_error(moran_i(x), '`space` is required')
  expect_error(moran_i(space = d), '`x` is required')
  expect_error(moran_i(x[1:2], dist(1:2)), 'at least 3 samples')
  expect_error(moran_i(x, d, style = 'col'), "'arg' should be one of")
  expect_error(moran_i(x, d, alternative = 'up'), "'arg' should be one of")
  expect_error(moran_i(c(x[-1], S1 = NA), d), 'must not contain NA')
  expect_error(moran_i(letters[1:30], d), '`x` must be numeric')
  expect_error(moran_i('val', d), 'names a column, but `data` is not given')
  expect_error(moran_i('nope', d, meta), 'must be the name of a column of `data`')
  expect_error(moran_i('val', c('east', 'up'), meta), "not in `data`: 'up'")
  expect_error(moran_i('val', 'east'), '`space` names columns, but `data` is not given')
  meta$chr <- letters[1:30]
  expect_error(moran_i('chr', d, meta), "`data` column 'chr' is not numeric")
  expect_error(moran_i(1:3, dist(c(0, 0, 1))), 'at the same place')
  expect_error(moran_i(1:3, data.frame(a = c(0, 1, 1), b = c(0, 2, 2))), "samples '2' and '3' at the same place")
  expect_error(moran_i(1:4, matrix(-1, 4, 4)), 'finite and not negative')
  expect_error(moran_i(1:4, matrix(0, 4, 4)), 'Every spatial weight')
  expect_error(moran_i(1:4, matrix('a', 4, 4)), 'must be numeric')
  expect_error(moran_i(1:3, data.frame(a = c('x', 'y', 'z'))), 'must be numeric columns')
  expect_error(moran_i(1:3, list(1, 2, 3)), '`space` must be a matrix of spatial weights')
  expect_error(moran_i(1:30, d, permutations = -1), '`permutations` must be')
  expect_error(moran_i(unname(x[1:29]), d), '`x` has 29 values, but `space` has 30 samples')
  expect_error(moran_i(setNames(1:30, paste0('T', 1:30)), d), 'have no value')
  wn <- w; colnames(wn) <- rev(colnames(w))
  expect_error(moran_i(x, wn), 'different row and column names')

#})
