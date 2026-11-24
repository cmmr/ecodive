#test_that("lefse", {

  # lefse(): LEfSe 1.1.2's tests and scores (with its own jitter and
  # subsamples injected), R's kruskal.test() and coin's wilcox_test(),
  # lefser 1.12.1; threads, input types, the wrapper, edge cases, errors,
  # printing and the plot. Fixture: data-raw/expected/lefse.r and .py.

  fx  <- readRDS(file.path('expected', 'lefse.rds'))
  rel <- function (a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

  run <- function (nm, details = TRUE, ...) {
    cs  <- fx$cases[[nm]]
    st  <- cs$settings
    inp <- fx$inputs[[st$input]]
    args <- list(
      inp$counts, inp$class, subgroups = inp$subclass, p_adjust = 'none',
      wilcoxon = st$wilc == 1, one_against_one = st$multiclass == 1,
      same_subgroups = st$same_subcl == 1, min_subgroup = st$min_c,
      n_boots = st$n_boots, boot_fraction = st$f, ... )
    if (details) do.call(lefse_details, c(args, list(jitter = cs$jitter, boots = cs$boots)))
    else do.call(lefse, args)
  }


  # -----------------------------------------------------------------------
  # LEfSe 1.1.2, with its jitter and subsamples.
  # -----------------------------------------------------------------------

  for (nm in names(fx$cases)) {
    cs  <- fx$cases[[nm]]
    inp <- fx$inputs[[cs$settings$input]]
    d   <- run(nm, cpus = 2L)
    r   <- d$result

    # The table scaled to 1e6, exactly as numerical_values().
    norm <- .Call(C_lefse_columns, inp$counts, 1L, 1e6, seq_len(ncol(inp$counts)) - 1L, 1L)
    expect_identical(unname(norm), unname(cs$norm), info = nm)

    # Class means, Kruskal-Wallis (LEfSe's, and R's on LEfSe's table).
    expect_true(rel(t(attr(r, 'means')), t(cs$means)) < 1e-12, info = nm)
    ok <- unname(!is.na(cs$kw_p))
    expect_identical(is.nan(r$p_value), !ok, info = nm)
    expect_true(rel(r$p_value[ok], cs$kw_p[ok]) < 1e-12, info = nm)
    expect_true(rel(r$p_value[ok], cs$kruskal[ok]) < 1e-12, info = nm)
    expect_identical(r$p_value, r$p_adjusted, info = nm)

    # Wilcoxon decisions, the features scored, and their scores.
    expect_identical(unname(r$wilcoxon), unname(cs$wilcoxon), info = nm)
    expect_identical(r$feature[d$features], names(cs$lda), info = nm)
    expect_true(rel(r$lda_score[d$features], cs$lda) < if (nm == 'wide') 1e-6 else 1e-10, info = nm)
    expect_identical(r$reject, r$feature %in% names(cs$lda)[cs$lda > 2], info = nm)
    expect_identical(is.nan(r$lda_score), !r$feature %in% names(cs$lda), info = nm)

    # enriched and log_max_mean, as LEfSe's output.
    m <- cs$means
    expect_identical(r$enriched, rownames(m)[max.col(t(m), 'first')], info = nm)
    expect_equal(r$log_max_mean, unname(log10(pmax(apply(m, 2, max), 1))), tolerance = 1e-12, info = nm)

    # With ecodive's own draws: the same decisions, close scores.
    r2 <- run(nm, details = FALSE, cpus = 2L)
    expect_identical(r2[, 1:7], r[, 1:7], info = nm)
    expect_identical(r2$reject, r$reject, info = nm)
    expect_true(max(abs(r2$lda_score - r$lda_score) / r$lda_score, na.rm = TRUE) < if (cs$settings$n_boots == 1) 0.2 else 0.1, info = nm)
  }

  # The 'low' case's rare features score below 2.
  r <- run('low')$result
  expect_true(any(r$wilcoxon %in% TRUE & !r$reject))

  # The Wilcoxon p-values are coin's (min_subgroup = 5: every pair tested).
  d  <- run('sub_boots', cpus = 1L)
  fi <- which(d$result$p_value < 0.05)
  expect_true(rel(d$wilcoxon, fx$coin[fi, c(1, 3, 2, 4)]) < 1e-10)

  # The rank-deficient case scores more features than a subsample has
  # samples.
  expect_true(length(run('wide')$features) > nrow(fx$cases$wide$boots))


  # -----------------------------------------------------------------------
  # lefser 1.12.1 (R's random numbers, MASS's default tol, signed scores
  # from mean(b) - mean(a)): for the features richer in the second class,
  # its scores and LEfSe's agree closely; for those richer in the first,
  # lefser's are negative and smaller.
  # -----------------------------------------------------------------------

  inp <- fx$inputs$two
  r   <- lefse(inp$counts, inp$class, p_adjust = 'none', cpus = 2L)
  lr  <- fx$lefser$scores
  s   <- r$lda_score[match(names(lr), r$feature)]
  e   <- r$enriched[match(names(lr), r$feature)]
  b   <- e == 'b' & names(lr) != 'F10'   # F10: jittered, differently in lefser
  expect_true(rel(s[b], lr[b]) < 0.03)
  expect_true(all(lr[e == 'a'] < 0) && all(-lr[e == 'a'] <= s[e == 'a'] + 0.05))
  expect_identical(sort(names(lr)), sort(r$feature[r$reject]))


  # -----------------------------------------------------------------------
  # Threads, input types, margin, the wrapper.
  # -----------------------------------------------------------------------

  x   <- fx$inputs$sub$counts
  grp <- fx$inputs$sub$class
  sg  <- fx$inputs$sub$subclass
  r1  <- lefse(x, grp, subgroups = sg, cpus = 1L)
  for (cp in c(2L, 3L, 4L, 8L))
    expect_identical(lefse(x, grp, subgroups = sg, cpus = cp), r1, info = cp)
  expect_identical(lefse(Matrix::Matrix(x, sparse = TRUE), grp, subgroups = sg, cpus = 2L), r1)
  expect_identical(lefse(slam::as.simple_triplet_matrix(x), grp, subgroups = sg, cpus = 2L), r1)
  expect_identical(lefse(t(x), grp, subgroups = sg, margin = 2L, cpus = 2L), r1)
  meta <- data.frame(g = grp, s = sg)
  expect_identical(lefse(x, 'g', meta, subgroups = 's', cpus = 2L), r1)
  expect_identical(da_test(x, 'lefse', 'g', meta, subgroups = 's', cpus = 2L), r1)
  expect_identical(da_test(x, 'LEfSe', formula = ~ g, data = meta, subgroups = 's', cpus = 2L), r1)
  expect_identical(match_da_test('lefser')$id, 'lefse')
  expect_identical(match_da_test('LEfSe')$id, 'lefse')
  info <- list_da_tests()[list_da_tests()$id == 'lefse', ]
  expect_identical(info$effect, 'lda_score')
  expect_true(is.na(info$base) && info$zero_free_ok && info$design == 'groups')

  # Rank-deficient bootstraps, threaded.
  w <- fx$inputs$wide
  rw <- lefse(w$counts, w$class, cpus = 1L)
  expect_identical(lefse(w$counts, w$class, cpus = 4L), rw)

  # The seed changes the scores only.
  r9 <- lefse(x, grp, subgroups = sg, seed = 9L, cpus = 2L)
  expect_identical(r9[, 1:7], r1[, 1:7])
  expect_false(identical(r9$lda_score, r1$lda_score))

  # A sample scaled by a constant changes nothing but rounding.
  x2 <- x
  x2[3, ] <- x2[3, ] * 7
  r2 <- lefse(x2, grp, subgroups = sg, cpus = 2L)
  expect_equal(r2$p_value, r1$p_value, tolerance = 1e-10)
  expect_equal(r2$lda_score, r1$lda_score, tolerance = 1e-8)
  expect_identical(r2$reject, r1$reject)

  # total = NULL takes the values as given: proportions, scaled by 1e6, are
  # the default.
  rp <- lefse(x / rowSums(x) * 1e6, grp, subgroups = sg, total = NULL, cpus = 2L)
  expect_equal(rp$lda_score, r1$lda_score, tolerance = 1e-8)
  expect_identical(attr(rp, 'total'), NA_real_)
  expect_identical(attr(r1, 'total'), 1e6)


  # -----------------------------------------------------------------------
  # The Wilcoxon step's rules.
  # -----------------------------------------------------------------------

  set.seed(5)
  y <- matrix(rpois(16 * 8, 30), 16, 8, dimnames = list(NULL, paste0('F', 1:8)))
  y[9:16, 1:3] <- y[9:16, 1:3] * 4
  g <- rep(c('a', 'b'), each = 8)

  # Under min_subgroup samples a class, only the medians are compared, and
  # the first comparison always passes (LEfSe's rule).
  r <- lefse(y, g, p_adjust = 'none', cpus = 1L)
  expect_true(all(r$wilcoxon[r$p_value < 0.05]))
  expect_true(all(is.na(r$wilcoxon[!(r$p_value < 0.05)])))
  # min_subgroup = 0: the rank-sum test decides.
  r0 <- lefse(y, g, p_adjust = 'none', min_subgroup = 0L, cpus = 1L)
  d0 <- lefse_details(y, g, p_adjust = 'none', min_subgroup = 0L, cpus = 1L)
  pw <- vapply(which(r0$p_value < 0.05), function (j)
    stats::wilcox.test(y[g == 'a', j] / rowSums(y)[g == 'a'] * 1e6, y[g == 'b', j] / rowSums(y)[g == 'b'] * 1e6,
                       exact = FALSE, correct = FALSE)$p.value, 0)
  expect_equal(as.vector(d0$wilcoxon), unname(pw), tolerance = 1e-12)
  expect_identical(r0$wilcoxon[r0$p_value < 0.05], pw < 0.1)
  # wilcoxon = FALSE skips it.
  rn <- lefse(y, g, p_adjust = 'none', wilcoxon = FALSE, cpus = 1L)
  expect_true(all(is.na(rn$wilcoxon)))
  expect_identical(rn$reject, rn$p_value < 0.05 & rn$lda_score > 2 & !is.na(rn$lda_score))

  # p_adjust: on the Kruskal-Wallis p-values, before the alpha cut.
  rh <- lefse(y, g, cpus = 1L)
  expect_identical(rh$p_adjusted, stats::p.adjust(rh$p_value, 'holm'))
  rb <- lefse(y, g, p_adjust = 'fdr', cpus = 1L)
  expect_identical(attr(rb, 'p_adjust'), 'fdr')
  expect_identical(rb$p_adjusted, stats::p.adjust(rb$p_value, 'BH'))

  # Three classes: one against all passes more than one against one.
  t1 <- run('three_all', details = FALSE, cpus = 1L)
  t2 <- run('three_one', details = FALSE, cpus = 1L)
  expect_true(all(t1$wilcoxon[t2$wilcoxon %in% TRUE]) && sum(t1$wilcoxon, na.rm = TRUE) > sum(t2$wilcoxon, na.rm = TRUE))


  # -----------------------------------------------------------------------
  # Edge cases.
  # -----------------------------------------------------------------------

  # A feature absent everywhere, one constant: NaN, never reported.
  z <- y
  z[, 7] <- 0
  z[, 8] <- 5
  r <- lefse(z, g, cpus = 1L)
  expect_true(is.nan(r$p_value[7]) && is.nan(r$kw_stat[7]) && !r$reject[7] && is.na(r$wilcoxon[7]))
  expect_identical(r$enriched[7], 'a')
  expect_identical(r$log_max_mean[7], 0)

  # Nothing passes: no scores, no error.
  set.seed(6)
  nul <- matrix(rpois(20 * 6, 50), 20, 6)
  r <- lefse(nul, rep(c('a', 'b'), 10), alpha = 1e-9, cpus = 1L)
  expect_true(!any(r$reject) && all(is.nan(r$lda_score)) && all(is.na(r$wilcoxon)))

  # One surviving feature; one bootstrap.
  one <- y[, c(1, 4:6)]
  r <- lefse(one, g, total = NULL, p_adjust = 'none', alpha = 0.01, n_boots = 1L, cpus = 1L)
  expect_identical(sum(!is.nan(r$lda_score)), 1L)
  expect_true(r$reject[1])

  # A class with one sample.
  g1 <- c('a', rep('b', 15))
  r <- lefse(y, g1, p_adjust = 'none', cpus = 2L)
  expect_true(inherits(r, 'ecodive_lefse') && !anyNA(r$reject))

  # A negative threshold reports every feature that passed both tests.
  r <- lefse(y, g, lda_threshold = -1, cpus = 1L)
  expect_identical(r$reject, r$wilcoxon %in% TRUE)

  # Factor levels set the classes' order.
  r <- lefse(y, factor(g, c('b', 'a')), cpus = 1L)
  expect_identical(attr(r, 'groups'), c('b', 'a'))
  expect_identical(rownames(attr(r, 'means')), c('b', 'a'))

  # Errors.
  expect_error(lefse(y, rep('a', 16)), 'at least two groups')
  expect_error(lefse(y, c(NA, g[-1])), 'NAs')
  expect_error(lefse(y, g, subgroups = c(NA, g[-1])), '`subgroups` must not contain NAs')
  expect_error(lefse(y, g, subgroups = 1:3), '`subgroups`')
  expect_error(lefse(y, g, boot_fraction = 0), '`boot_fraction`')
  expect_error(lefse(y, g, boot_fraction = 1.5), '`boot_fraction`')
  expect_error(lefse(y, g, boot_fraction = 0.1), 'more than the 2 groups')
  expect_error(lefse(y, g, n_boots = 0), '`n_boots`')
  expect_error(lefse(y, g, total = 0), '`total`')
  expect_error(lefse(y, g, alpha = 1), '`alpha`')
  expect_error(lefse(y, g, wilcoxon = NA), '`wilcoxon`')
  expect_error(lefse(y, g, p_adjust = 'nope'), '`p_adjust`')
  expect_error(lefse(-y, g), 'negative')
  expect_error(lefse(y, g, seed = 1.5), '`seed`')
  expect_error(lefse_details(y, g, jitter = matrix(1, 2, 2)), '`jitter`')


  # -----------------------------------------------------------------------
  # Attributes, printing, the plot.
  # -----------------------------------------------------------------------

  r <- r1
  expect_identical(names(r), c('feature', 'enriched', 'log_max_mean', 'kw_stat', 'p_value', 'p_adjusted', 'wilcoxon', 'lda_score', 'reject'))
  expect_identical(attr(r, 'method'), 'lefse')
  expect_identical(attr(r, 'groups'), c('a', 'b'))
  expect_identical(attr(r, 'subgroups'), c('s1', 's2'))
  expect_identical(attr(r, 'log_base'), NA_real_)
  expect_identical(attr(r, 'p_adjust'), 'holm')
  expect_identical(dim(attr(r, 'means')), c(2L, ncol(x)))
  expect_identical(class(r), c('ecodive_lefse', 'data.frame'))

  out <- utils::capture.output(print(r))
  expect_true(any(grepl('LEfSe: a vs b', out, fixed = TRUE)))

  grDevices::pdf(NULL)
  v <- withVisible(plot(r))
  expect_true(identical(v$value, r) && !v$visible)
  expect_identical(plot(r, n_features = 2, main = 'x', sub = 'y'), r)
  expect_silent(plot(run('three_all', details = FALSE, cpus = 1L)))
  expect_silent(plot(lefse(nul, rep(c('a', 'b'), 10), alpha = 1e-9, cpus = 1L)))
  expect_error(plot(r, n_features = -1), '`n_features`')
  grDevices::dev.off()
  expect_identical(paste('LDA scores of', da_method_name('lefse')), 'LDA scores of LEfSe')

#})
