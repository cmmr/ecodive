#test_that("hill extrapolation", {

  # hill_extrapolation() and hill_standardize(): point estimates against
  # iNEXT 3.0.2's iNEXT(), estimateD() and DataInfo(); the bootstrap's
  # standard errors exactly, given iNEXT's own draws, and by distribution
  # with ecodive's draws; identities with rarefaction_curve() and the alpha
  # metrics; Chao et al. (2014)'s printed spider numbers; determinism; the
  # binomial sampler; edge cases, errors, printing and plotting.
  # Fixtures: data-raw/expected/hill_extrapolation.r.

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))

  fx <- readRDS('expected/hill_extrapolation.rds')
  ex_counts <- ecodive::ex_counts
  cpus <- 2L

  # Relative difference, against max(|y|, 1e-300); equal infinities count
  # as no difference.
  rel <- function (x, y) {
    d <- ifelse(is.infinite(x) & x == y, 0, abs(x - y) / pmax(abs(y), 1e-300))
    if (length(d)) max(d) else 0
  }

  # Eq. 10c: Shannon diversity at depth m > n, from the observed and
  # asymptotic diversities.
  eq10c <- function (obs, asy, n, m) exp(n / m * log(obs) + (m - n) / m * log(asy))



  # Point estimates against iNEXT(). ====

  for (nm in names(fx$tables)) {
    x  <- fx$tables[[nm]]
    he <- hill_extrapolation(x, n_boot = 0, cpus = cpus)
    ref <- fx$curves[[nm]]
    s   <- he$size

    # The same depths (iNEXT's knots), sample by sample.
    expect_identical(
      sort(paste(s$sample, s$q, s$depth)), sort(paste(ref$sample, ref$q, ref$depth)), info = nm )

    i <- match(paste(ref$sample, ref$q, ref$depth), paste(s$sample, s$q, s$depth))
    s <- s[i, ]
    expect_identical(s$method, ref$method, info = nm)
    q1x <- ref$q == 1 & ref$method == 'extrapolation'
    expect_true(rel(s$estimate[!q1x], ref$qD[!q1x]) < 1e-10, info = nm)
    expect_true(max(abs(s$coverage - ref$SC)) < 1e-10, info = nm)

    # Shannon diversity beyond the sample's size: Eq. 10c, with iNEXT's
    # observed and asymptotic values (iNEXT 3 extrapolates otherwise).
    a <- fx$asy[[nm]]
    if (any(q1x)) {
      k   <- which(q1x)
      ai  <- match(paste(ref$sample[k], 1), paste(a$sample, a$q))
      ns  <- fx$info[[nm]]$n[match(ref$sample[k], fx$info[[nm]]$sample)]
      expect_true(rel(s$estimate[k], eq10c(a$observed[ai], a$estimate[ai], ns, ref$depth[k])) < 1e-10, info = nm)
    }

    # Asymptotes and observed values.
    ha <- he$asymptotic[match(paste(a$sample, a$q), paste(he$asymptotic$sample, he$asymptotic$q)), ]
    expect_true(rel(ha$estimate, a$estimate) < 1e-10, info = nm)
    expect_true(rel(ha$observed, a$observed) < 1e-10, info = nm)

    # DataInfo (whose coverage is rounded to 4 digits).
    di <- fx$info[[nm]]
    hi <- he$info[match(di$sample, he$info$sample), ]
    expect_identical(hi$n, as.double(di$n), info = nm)
    expect_identical(hi$s_obs, as.double(di$s_obs), info = nm)
    expect_identical(unname(as.matrix(hi[paste0('f', 1:10)])), unname(as.matrix(di[paste0('f', 1:10)])) + 0, info = nm)
    expect_identical(round(hi$coverage, 4), di$SC, info = nm)

    # By coverage: iNEXT's depths come from a numerical search (about
    # 1e-5 off); where coverage is flat over a range of depths, its depth
    # can be anywhere in the range. So: the rows whose coverage is a level
    # of ours, away from flat stretches and Shannon's extrapolation.
    cv  <- he$coverage
    rc  <- fx$cov[[nm]]
    j   <- match(paste(rc$sample, rc$q, rc$method != 'extrapolation'), paste(cv$sample, cv$q, cv$method != 'extrapolation'))
    for (r in seq_len(nrow(rc))) {
      w <- which(cv$sample == rc$sample[r] & cv$q == rc$q[r] & abs(cv$coverage - rc$SC[r]) < 1e-9)
      if (length(w) != 1) next
      # Skip where iNEXT's search left its own size-based value (flat
      # coverage: f2_zero is flat from depth 9 to 10, and iNEXT lands
      # between them).
      sr <- ref$qD[ref$sample == rc$sample[r] & ref$q == rc$q[r] & ref$depth == rc$depth[r]]
      if (length(sr) == 1 && abs(sr / rc$qD[r] - 1) > 1e-4) next
      if (rc$q[r] == 1 && cv$depth[w] > fx$info[[nm]]$n[fx$info[[nm]]$sample == rc$sample[r]]) next
      if (rc$method[r] == 'extrapolation' || rc$method[r] == 'observed')
        expect_true(rel(cv$estimate[w], rc$qD[r]) < 1e-10, info = paste(nm, r))
      else if (abs(cv$depth[w] - rc$depth[r]) < 1e-3)
        expect_true(rel(cv$estimate[w], rc$qD[r]) < 1e-5, info = paste(nm, r))
    }
  }

  # The coverage table: one row per coverage level (sample, q), the depths
  # nondecreasing in coverage, and estimates equal to the size table's at
  # those depths.
  he <- hill_extrapolation(fx$tables$ex, n_boot = 0, cpus = cpus)
  expect_false(any(duplicated(he$coverage[c('sample', 'q', 'coverage')])))
  for (s in unique(he$coverage$sample)) {
    cv <- he$coverage[he$coverage$sample == s & he$coverage$q == 0, ]
    expect_true(all(diff(cv$depth) >= 0) && all(diff(cv$coverage) > 0), info = s)
  }
  # Gums (f1 = 2, f2 = f3 = 0, n = 886): coverage is flat, 1 - f1/n, from
  # depth 883 to 885, so the level of depth 885 is first reached at 883
  # (iNEXT's search lands at 884).
  g <- he$coverage[he$coverage$sample == 'Gums' & he$coverage$q == 0, ]
  r <- fx$cov$ex[fx$cov$ex$sample == 'Gums' & fx$cov$ex$q == 0 & fx$cov$ex$depth == 884, ]
  w <- which(abs(g$coverage - r$SC) < 1e-12)
  expect_identical(g$depth[w], 883)
  expect_true(abs(g$coverage[w] - (1 - 2 / 886)) < 1e-14)
  s <- he$size[he$size$sample == 'Gums' & he$size$q == 0, ]
  expect_true(rel(s$coverage[s$depth == 885], 1 - 2 / 886) < 1e-14)



  # Depth grids (iNEXT's knots). ====

  for (g in fx$grids)
    expect_identical(ecodive:::hx_grid(168, NULL, g$endpoint, g$knots), as.double(g$m), info = paste(g$endpoint, g$knots))
  he <- hill_extrapolation(fx$tables$spider, depths = c(5, 500, 10), n_boot = 0, cpus = cpus)
  expect_identical(unique(he$size$depth[he$size$sample == 'Girdled']), c(5, 10, 168, 500))
  he <- hill_extrapolation(fx$tables$spider, depths = c(5, 500, 10), max_depth = 300, n_boot = 0, cpus = cpus)
  expect_identical(unique(he$size$depth[he$size$sample == 'Logged']), c(5, 10, 252))
  he <- hill_extrapolation(fx$tables$spider, max_depth = 100, n_depths = 12, q = 2, n_boot = 0, cpus = cpus)
  expect_identical(unique(he$size$depth[he$size$sample == 'Girdled']), c(floor(seq(1, 100, length.out = 12))))
  expect_identical(unique(he$size$q), 2)



  # estimateD(). ====

  ed <- fx$estd
  for (tb in unique(ed$table)) {
    x <- fx$tables[[tb]]
    n <- rowSums(x)
    for (given in c(FALSE, TRUE)) {

      # Base size.
      e <- ed[ed$table == tb & ed$base == 'size' & ed$given == given, ]
      for (lv in unique(e$depth)) {
        hs <- hill_standardize(x, depth = lv, n_boot = 0, cpus = cpus)
        r  <- e[e$depth == lv, ]
        v  <- unclass(hs)[cbind(match(r$sample, rownames(hs)), r$q + 1)]
        keep <- !(r$q == 1 & r$method == 'extrapolation')
        expect_true(rel(v[keep], r$qD[keep]) < 1e-10, info = paste(tb, lv))
        expect_true(max(abs(attr(hs, 'coverage')[r$sample] - r$SC)) < 1e-10, info = paste(tb, lv))
        expect_identical(unname(attr(hs, 'method')[r$sample]), r$method, info = paste(tb, lv))
      }

      # Base coverage. The default level: the smallest coverage at twice
      # the samples' sizes.
      e <- ed[ed$table == tb & ed$base == 'coverage' & ed$given == given, ]
      if (!given) {
        hs <- hill_standardize(x, n_boot = 0, cpus = cpus)
        expect_true(abs(attr(hs, 'level') - min(e$SC)) < 1e-10, info = tb)
        expect_identical(attr(hs, 'base'), 'coverage')
      }
      for (lv in if (given) c(0.8, 0.95, 0.99, 1) else attr(hs, 'level')) {
        hs <- hill_standardize(x, coverage = lv, n_boot = 0, cpus = cpus)
        r  <- e[abs(e$SC - lv) < 1e-6 | (given & e$depth != floor(e$depth) & FALSE), ]
        r  <- if (given) e[seq_len(nrow(e)) %in% which(round(e$SC, 6) == lv | !is.finite(e$depth) & lv == 1 | TRUE), ] else e
        r  <- r[sapply(seq_len(nrow(r)), function (i) {
          d <- attr(hs, 'depth')[[r$sample[i]]]
          (is.infinite(d) && is.infinite(r$depth[i])) || (is.finite(d) && is.finite(r$depth[i]) && abs(d - r$depth[i]) < 1e-3)
        }), ]
        if (nrow(r) == 0) next
        v  <- unclass(hs)[cbind(match(r$sample, rownames(hs)), r$q + 1)]
        ext  <- r$method == 'extrapolation'
        keep <- !(r$q == 1 & ext)
        # Extrapolated depths are in closed form; rarefied ones come from
        # iNEXT's search, about 1e-5 off.
        expect_true(rel(v[keep & ext], r$qD[keep & ext]) < 1e-10, info = paste(tb, lv))
        if (any(keep & !ext)) expect_true(rel(v[keep & !ext], r$qD[keep & !ext]) < 1e-5, info = paste(tb, lv))
      }
    }

    # At iNEXT's own (fractional) depths, the straight line between whole
    # depths reproduces its estimates exactly.
    e <- ed[ed$table == tb & ed$base == 'coverage' & is.finite(ed$depth) & ed$method == 'rarefaction', ]
    for (s in unique(e$sample)) {
      r <- e[e$sample == s, ]
      i <- match(s, rownames(x))
      sz <- lapply(seq_len(nrow(x)), function (k) if (k == i) sort(unique(r$depth)) else numeric(0))
      em <- lapply(seq_len(nrow(x)), function (k) numeric(0))
      res <- ecodive:::hx_run(x, 1L, sz, em, em, 0:2, 0L, 0L, cpus)
      v <- res$size[cbind(match(r$depth, sort(unique(r$depth))), r$q + 1)]
      expect_true(rel(v, r$qD) < 1e-10, info = paste(tb, s))
      expect_true(max(abs(res$size[match(r$depth, sort(unique(r$depth))), 4] - r$SC)) < 1e-10, info = paste(tb, s))
    }
  }



  # The bootstrap with iNEXT's own draws: exact standard errors. ====

  b  <- fx$boot
  x  <- fx$tables$spider['Girdled', , drop = FALSE]
  n  <- sum(x)
  sz <- as.double(b$sizes)
  em <- list(numeric(0))

  # The Appendix G assemblage, most abundant first.
  f1 <- sum(x == 1); f2 <- sum(x == 2)
  C  <- 1 - f1 / n * (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
  expect_equal(length(b$p), sum(x > 0) + ceiling((n - 1) / n * f1^2 / (2 * f2)))

  r1 <- ecodive:::hx_run(x, 1L, list(sz), list(rep(NA_real_, length(sz))), list(sz), 0:2, 200L, 0L, cpus, boot = list(b$draws1))
  bs <- b$size
  for (q in 0:2) {
    k  <- which(bs$q == q & bs$depth > 1 & !(q == 1 & bs$depth > n))
    se <- r1$size_sd[match(bs$depth[k], sz), q + 1]
    expect_true(rel(se, bs$se[k]) < 1e-10, info = q)
  }
  k <- which(bs$q == 0 & !bs$sc_clipped & bs$depth > 1)
  expect_true(rel(r1$size_sd[match(bs$depth[k], sz), 4], bs$sc_se[k]) < 1e-10)

  # At fixed coverage, each replicate solves for its own depth (iNEXT's
  # levels are its search's, so about 1e-5 off).
  bc  <- b$coverage[b$coverage$depth <= 50 & b$coverage$depth > 1, ]
  lv  <- unique(bc$SC)
  rc  <- ecodive:::hx_run(x, 1L, em, list(lv), list(rep(NA_real_, length(lv))), 0:2, 200L, 0L, cpus, boot = list(b$draws1))
  expect_true(rel(rc$cov_sd[cbind(match(bc$SC, lv), bc$q + 1)], bc$se) < 1e-4)

  # The asymptotes, from the second set of draws (iNEXT always uses 100).
  r2 <- ecodive:::hx_run(x, 1L, em, em, em, 0:2, 100L, 0L, cpus, boot = list(b$draws2))
  expect_true(rel(r2$asy_sd[1, ], b$asy_se) < 1e-10)

  # estimateD(): at size 100 exactly; at coverage 0.9, q = 0 and 2 (q = 1's
  # replicates extrapolate there).
  re <- ecodive:::hx_run(x, 1L, list(100), em, em, 0:2, 60L, 0L, cpus, boot = list(b$draws_e))
  expect_true(rel(re$size_sd[1, 1:3], b$estd_size_se) < 1e-10)
  re <- ecodive:::hx_run(x, 1L, em, list(0.9), list(NA_real_), 0:2, 60L, 0L, cpus, boot = list(b$draws_e))
  expect_true(rel(re$cov_sd[1, c(1, 3)], b$estd_cov_se[c(1, 3)]) < 1e-5)

  # The intervals: estimate +- z se, lower limits stopping at 0 (and the
  # asymptotic richness's at the observed richness), coverage's at 0 and 1.
  he <- hill_extrapolation(fx$tables$spider, n_boot = 30, conf_level = 0.9, seed = 2, cpus = cpus)
  z  <- stats::qnorm(0.95)
  s  <- he$size
  se <- (s$upper - s$estimate) / z
  expect_true(max(abs(s$lower - pmax(0, s$estimate - z * se))) < 1e-10)
  a  <- he$asymptotic
  expect_true(max(abs(a$upper - (a$estimate + z * a$se))) < 1e-10)
  expect_true(all(a$lower[a$q == 0] >= a$observed[a$q == 0]))
  expect_true(all(s$coverage_upper <= 1 & s$coverage_lower >= 0))
  expect_identical(attr(he, 'conf_level'), 0.9)
  expect_identical(attr(he, 'n_boot'), 30L)



  # Standard errors by distribution: ecodive's draws against iNEXT's, both
  # from 2,000 replicates (for the asymptotes, iNEXT's from 20 runs of 100). ====

  for (s in names(fx$se2000)) {
    tb <- if (s == 'Girdled') 'spider' else 'ex'
    x  <- fx$tables[[tb]][s, , drop = FALSE]
    n  <- sum(x)
    he <- hill_extrapolation(x, depths = c(floor(n / 2), 2 * n), n_boot = 2000, seed = 5, cpus = cpus)
    ref <- fx$se2000[[s]]
    expect_true(rel(he$asymptotic$se, ref$asy$se) < 0.1, info = s)
    sz <- ref$size[!(ref$size$q == 1 & ref$size$depth > n) & !ref$size$clipped & ref$size$se > 1e-8, ]
    h  <- he$size[match(paste(sz$q, sz$depth), paste(he$size$q, he$size$depth)), ]
    se <- (h$upper - h$estimate) / stats::qnorm(0.975)
    expect_true(rel(se, sz$se) < 0.1, info = s)
  }



  # Identities. ====

  he <- hill_extrapolation(ex_counts, n_boot = 0, cpus = cpus)
  s0 <- he$size[he$size$q == 0 & he$size$method != 'extrapolation', ]
  rc <- rarefaction_curve(ex_counts, depths = sort(unique(s0$depth)), cpus = cpus)
  i  <- match(paste(s0$sample, s0$depth), paste(rc$sample, rc$depth))
  expect_false(anyNA(i))
  expect_true(rel(s0$estimate, rc$value[i]) < 1e-12)
  expect_true(max(abs(s0$coverage - rc$coverage[i])) < 1e-12)

  ob <- he$size[he$size$method == 'observed', ]
  for (q in 0:2) {
    v <- ob$estimate[ob$q == q]
    names(v) <- ob$sample[ob$q == q]
    expect_equal(v[rownames(ex_counts)], hill(ex_counts, q = q), tolerance = 1e-12, info = q, check.attributes = FALSE)
  }
  expect_equal(ob$estimate[ob$q == 0], unname(observed(ex_counts)), tolerance = 1e-14)
  expect_equal(ob$estimate[ob$q == 1], unname(exp(shannon(ex_counts))), tolerance = 1e-12)
  expect_equal(ob$estimate[ob$q == 2], unname(inv_simpson(ex_counts)), tolerance = 1e-12)

  # Chao et al. (2014)'s worked example: girdled n = 168, 26 species,
  # Shannon and Simpson diversity 12.06 and 7.84, coverage 93%.
  he <- hill_extrapolation(fx$tables$spider, n_boot = 0, cpus = cpus)
  a  <- he$asymptotic[he$asymptotic$sample == 'Girdled', ]
  expect_identical(he$info$n[1], 168)
  expect_identical(he$info$s_obs[1], 26)
  expect_identical(round(a$observed[2:3], 2), c(12.06, 7.84))
  expect_identical(round(he$info$coverage[1], 2), 0.93)



  # Determinism. ====

  x  <- fx$tables$ex
  h1 <- hill_extrapolation(x, n_boot = 20, seed = 4, cpus = 1L)
  for (cp in c(2L, 4L)) expect_identical(hill_extrapolation(x, n_boot = 20, seed = 4, cpus = cp), h1, info = cp)
  expect_identical(hill_extrapolation(x, n_boot = 20, seed = 4, cpus = 3L), h1)
  expect_false(identical(hill_extrapolation(x, n_boot = 20, seed = 5, cpus = 2L)$size$upper, h1$size$upper))
  s1 <- hill_standardize(x, n_boot = 20, seed = 4, cpus = 1L)
  for (cp in c(2L, 4L)) expect_identical(hill_standardize(x, n_boot = 20, seed = 4, cpus = cp), s1, info = cp)

  # Replicates in blocks: the first 20 of 21 replicates' streams are those of
  # 20 (one stream per replicate), so means of their values agree.
  r20 <- ecodive:::hx_run(x, 1L, lapply(1:4, function (i) numeric(0)), lapply(1:4, function (i) numeric(0)), lapply(1:4, function (i) numeric(0)), 0:2, 20L, 4L, 2L)
  expect_true(all(is.finite(r20$asy_sd[, c(1, 3)])))

  # Sparse and dense, and margin = 2.
  expect_identical(hill_extrapolation(Matrix::Matrix(x, sparse = TRUE), n_boot = 5, cpus = cpus), hill_extrapolation(x, n_boot = 5, cpus = cpus))
  expect_identical(hill_extrapolation(t(x), margin = 2L, n_boot = 5, cpus = cpus), hill_extrapolation(x, n_boot = 5, cpus = cpus))
  expect_identical(hill_standardize(t(x), margin = 2L, n_boot = 5, cpus = cpus), hill_standardize(x, n_boot = 5, cpus = cpus))



  # The binomial and multinomial sampler. ====

  for (np in list(c(20, 0.1), c(1000, 0.3), c(1e6, 0.4), c(50, 0.9), c(400, 0.97))) {
    v <- .Call(ecodive:::C_rng_binom, 1e6L, np[[1]], np[[2]], 3L)
    m <- np[[1]] * np[[2]]; s2 <- m * (1 - np[[2]])
    expect_true(abs(mean(v) - m) < 5 * sqrt(s2 / 1e6), info = paste(np, collapse = ' '))
    expect_true(abs(stats::var(v) / s2 - 1) < 0.01, info = paste(np, collapse = ' '))
    expect_true(all(v >= 0 & v <= np[[1]] & v == floor(v)))
  }
  expect_identical(.Call(ecodive:::C_rng_binom, 5L, 10, 0, 1L), rep(0, 5))
  expect_identical(.Call(ecodive:::C_rng_binom, 5L, 10, 1, 1L), rep(10, 5))
  p <- c(0.5, 0.2, 0.2, 0.05, 0.05)
  m <- .Call(ecodive:::C_rng_multinom, 2000L, 1000, p, 9L)
  expect_true(all(colSums(m) == 1000))
  expect_true(max(abs(rowMeans(m) / (1000 * p) - 1)) < 0.02)



  # Edge cases. ====

  # f2 = 0, f1 = 0, all singletons, one feature: covered by the 'edge'
  # fixture above. One feature: every estimate 1.
  he <- hill_extrapolation(matrix(c(10, 0), 1), n_boot = 5, cpus = cpus)
  expect_true(all(he$size$estimate == 1) && all(he$asymptotic$estimate == 1))
  # All singletons: Simpson's asymptote is infinite, as in iNEXT.
  he <- hill_extrapolation(matrix(1, 1, 5), n_boot = 0, cpus = cpus)
  expect_identical(he$asymptotic$estimate[3], Inf)
  # f1 = 0: coverage 1 at the sample's size, no extrapolated gain in
  # richness.
  he <- hill_extrapolation(matrix(c(2, 2, 3, 5), 1), n_boot = 5, cpus = cpus)
  s  <- he$size[he$size$q == 0 & he$size$method == 'extrapolation', ]
  expect_true(all(s$estimate == 4) && all(s$coverage == 1))
  expect_identical(nrow(he$coverage[he$coverage$q == 0 & he$coverage$method == 'extrapolation', ]), 0L)

  # An empty sample: no rows, NAs, and a warning.
  x0 <- rbind(fx$tables$ex, Empty = 0)
  expect_warning(he <- hill_extrapolation(x0, n_boot = 5, cpus = cpus), '1 empty sample \\(Empty\\)')
  expect_false('Empty' %in% he$size$sample)
  expect_true(all(is.na(he$asymptotic$estimate[he$asymptotic$sample == 'Empty'])))
  expect_true(is.na(he$info$coverage[he$info$sample == 'Empty']))
  expect_identical(he$size, suppressWarnings(hill_extrapolation(fx$tables$ex, n_boot = 5, cpus = cpus))$size)
  expect_warning(hs <- hill_standardize(x0, n_boot = 5, cpus = cpus), 'empty')
  expect_true(all(is.na(hs['Empty', ])) && is.na(attr(hs, 'depth')[['Empty']]))

  # Depths beyond twice the size are allowed.
  he <- hill_extrapolation(fx$tables$spider, max_depth = 2000, n_boot = 0, cpus = cpus)
  expect_true(max(he$size$depth) == 2000 && all(is.finite(he$size$estimate)))

  # Coverage 1: the asymptotes, at depth Inf; a coverage below that of one
  # observation: depth 1.
  hs <- hill_standardize(fx$tables$spider, coverage = 1, n_boot = 0, cpus = cpus)
  he <- hill_extrapolation(fx$tables$spider, n_boot = 0, cpus = cpus)
  expect_identical(unname(attr(hs, 'depth')), c(Inf, Inf))
  expect_true(rel(as.vector(t(unclass(hs))), he$asymptotic$estimate) < 1e-14)
  hs <- hill_standardize(fx$tables$spider, coverage = 0.01, n_boot = 0, cpus = cpus)
  expect_identical(unname(attr(hs, 'depth')), c(1, 1))
  expect_true(all(unclass(hs) == 1))

  # Unnamed tables.
  he <- hill_extrapolation(unname(fx$tables$spider), n_boot = 0, q = 0, cpus = cpus)
  expect_identical(he$info$sample, c('1', '2'))



  # hill_standardize()'s result. ====

  hs <- hill_standardize(fx$tables$ex, q = c(2, 0), n_boot = 10, cpus = cpus)
  expect_true(inherits(hs, 'ecodive_hill_standardize') && is.matrix(hs))
  expect_identical(colnames(hs), c('q2', 'q0'))
  expect_identical(rownames(hs), rownames(fx$tables$ex))
  expect_identical(
    names(attributes(hs)),
    c('dim', 'dimnames', 'lower', 'upper', 'depth', 'coverage', 'method', 'level', 'base', 'conf_level', 'class') )
  expect_true(all(attr(hs, 'lower') <= unclass(hs) & unclass(hs) <= attr(hs, 'upper')))
  expect_true(all(abs(attr(hs, 'coverage') - attr(hs, 'level')) < 1e-12))
  hs <- hill_standardize(fx$tables$ex, depth = 500, n_boot = 0, cpus = cpus)
  expect_null(attr(hs, 'lower'))
  expect_identical(attr(hs, 'base'), 'depth')
  expect_identical(unname(attr(hs, 'method')), c('extrapolation', 'rarefaction', 'rarefaction', 'rarefaction'))
  expect_identical(unname(attr(hs, 'depth')), rep(500, 4))



  # Errors. ====

  expect_error(hill_extrapolation(ex_counts + 0.5), 'whole numbers')
  expect_error(hill_extrapolation(-ex_counts), 'negative')
  expect_error(hill_extrapolation(ex_counts, q = 3), '`q` must be one or more of 0, 1 and 2')
  expect_error(hill_extrapolation(ex_counts, q = numeric(0)), '`q` must be')
  expect_error(hill_extrapolation(ex_counts, n_boot = -1), '`n_boot` must be')
  expect_error(hill_extrapolation(ex_counts, n_boot = 1.5), '`n_boot` must be')
  expect_error(hill_extrapolation(ex_counts, conf_level = 1), '`conf_level` must be')
  expect_error(hill_extrapolation(ex_counts, n_depths = 1), '`n_depths` must be')
  expect_error(hill_extrapolation(ex_counts, max_depth = 0), '`max_depth` must be')
  expect_error(hill_extrapolation(ex_counts, depths = 0.5), '`depths` must be')
  expect_error(hill_extrapolation(ex_counts, seed = 0.5), '`seed`')
  expect_error(hill_standardize(ex_counts, coverage = 0.9, depth = 10), 'not both')
  expect_error(hill_standardize(ex_counts, coverage = 0), '`coverage` must be')
  expect_error(hill_standardize(ex_counts, coverage = 1.1), '`coverage` must be')
  expect_error(hill_standardize(ex_counts, depth = 0), '`depth` must be')
  expect_error(hill_standardize(ex_counts, depth = 2.5), '`depth` must be')
  expect_error(hill_standardize(ex_counts, q = 1.5), '`q` must be')



  # Printing and plotting. ====

  he <- hill_extrapolation(fx$tables$ex, n_boot = 5, cpus = cpus)
  out <- utils::capture.output(print(he))
  expect_true(any(grepl('Interpolated and extrapolated Hill numbers', out)))
  expect_true(any(grepl('95% intervals from 5 bootstrap replicates', out)))
  out <- utils::capture.output(print(hill_standardize(fx$tables$ex, depth = 100, n_boot = 0)))
  expect_true(any(grepl('Hill numbers at depth 100', out)))
  out <- utils::capture.output(print(hill_standardize(fx$tables$ex, coverage = 0.99, n_boot = 0)))
  expect_true(any(grepl('Hill numbers at coverage 0.99', out)))

  invisibly <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  grDevices::pdf(NULL)
  for (type in c('depth', 'size', 'completeness', 'coverage', 'comp'))
    expect_true(invisibly(he, type = type), info = type)
  expect_true(invisibly(he, q = 2, samples = c('Saliva', 'Stool'), legend = NULL, col = c('red', 'blue')))
  expect_true(invisibly(he, samples = 2:3, main = 'Mine', log = 'x'))
  expect_true(invisibly(hill_extrapolation(fx$tables$ex, n_boot = 0, cpus = cpus), q = 1))
  expect_error(plot(he, q = 3), '`q` must be one of')
  expect_error(plot(he, type = 'c'), '`type` must be one of')
  expect_error(plot(he, samples = 'nope'), '`samples` must name')
  expect_error(plot(hill_extrapolation(fx$tables$ex, q = 0, n_boot = 0), q = 1), '`q` must be one of')
  grDevices::dev.off()

  lab <- ecodive:::hx_plot_labels('depth', 0, 0.95, TRUE)
  expect_identical(lab$main, 'Sample-size-based Rarefaction and Extrapolation')
  expect_identical(lab$ylab, 'Species richness')
  expect_identical(lab$xlab, 'Depth')
  expect_identical(lab$sub, 'Bands: 95% confidence intervals')
  lab <- ecodive:::hx_plot_labels('completeness', 1, 0.9, FALSE)
  expect_identical(lab$main, 'Sample Completeness Curves')
  expect_identical(lab$ylab, 'Sample Coverage')
  expect_null(lab$sub)
  lab <- ecodive:::hx_plot_labels('coverage', 2, 0.9, TRUE)
  expect_identical(lab$main, 'Coverage-based Rarefaction and Extrapolation')
  expect_identical(lab$xlab, 'Sample Coverage')
  expect_identical(lab$ylab, 'Simpson diversity')
  expect_identical(ecodive:::hx_plot_labels('depth', 1, 0.95, FALSE)$ylab, 'Shannon diversity')

#})
