#test_that("valencia", {

  # Threads: on CRAN, at most 2 (_R_CHECK_LIMIT_CORES_).
  max_cpus <- if (at_home()) 4L else 2L

  cents <- ecodive::valencia_centroids
  subs  <- rownames(cents)

  # Valencia.py's form of the Yue-Clayton similarity, in base R, over the
  # union of the taxa: sum(pq) / (sum((p - q)^2) + sum(pq)).
  yue_clayton <- function (counts, cents, totals = rowSums(counts)) {
    taxa <- union(colnames(counts), colnames(cents))
    P <- Q <- NULL
    P <- matrix(0, nrow(counts), length(taxa), dimnames = list(rownames(counts), taxa))
    Q <- matrix(0, nrow(cents),  length(taxa), dimnames = list(rownames(cents),  taxa))
    P[, colnames(counts)] <- counts / totals
    Q[, colnames(cents)]  <- cents
    t(apply(P, 1L, function (p) apply(Q, 1L, function (q) {
      pq <- sum(p * q)
      pq / (sum((p - q)^2) + pq) })))
  }

  # Read counts drawn around the centroids, with taxa they lack, and the
  # centroids' taxa in another order.
  set.seed(1)
  draw <- function (k, depth) {
    x <- stats::rmultinom(1, depth, cents[k, ] + 1e-4)[, 1]
    c(x, Extra_1 = stats::rpois(1, 20), Extra_2 = 0)
  }
  counts <- t(sapply(rep(seq_along(subs), each = 9), draw, depth = 5000))
  rownames(counts) <- sprintf('S%03d', seq_len(nrow(counts)))
  counts <- counts[, c(ncol(counts):1)]
  storage.mode(counts) <- 'integer'



  # Results from Valencia.py (data-raw/expected/valencia.r). ====

  if (at_home()) {

    ref <- readRDS(file.path('expected', 'valencia.rds'))

    for (case in names(ref)) {

      r   <- ref[[case]]
      cen <- if (is.null(r$centroids)) cents else r$centroids
      sub_ref <- ifelse(r$subCST == '', NA, r$subCST)
      cst_ref <- ifelse(r$CST == '', NA, r$CST)

      # With Valencia.py's read counts; and with the row sums where those are
      # the read counts.
      for (nm in c('totals', 'rowsum')) {

        if (nm == 'rowsum' && !all(r$totals == rowSums(r$counts))) next
        info <- paste(case, nm)
        tt   <- if (nm == 'totals') r$totals

        v  <- suppressWarnings(valencia(r$counts, centroids = cen, totals = tt))
        sc <- attr(v, 'scores')

        expect_identical(dimnames(sc), dimnames(r$scores), info = info)
        expect_identical(is.nan(sc), is.nan(r$scores), info = info)
        expect_true(max(abs(sc - r$scores), na.rm = TRUE) < 1e-12, info = info)
        expect_identical(as.character(v$type), sub_ref, info = info)
        expect_identical(as.character(v$cst),  cst_ref, info = info)
        expect_identical(is.nan(v$score), is.nan(r$score), info = info)
        expect_true(max(abs(v$score - r$score), na.rm = TRUE) < 1e-12, info = info)
      }
    }

    # Exact ties go to the first sub-CST, as pandas' idxmax() breaks them.
    r <- ref$ties
    v <- suppressWarnings(valencia(r$counts, centroids = r$centroids))
    sc <- attr(v, 'scores')
    expect_true(all(sc[, 'I-A'] == sc[, 'I-B']))
    expect_true(all(sc[, 'III-A'] == sc[, 'III-B']))
    expect_false(any(v$type %in% c('I-B', 'III-B')))
    expect_true(sum(v$type == 'I-A') >= 10 && sum(v$type == 'III-A') >= 10)

    # The centroids as a data frame read from VALENCIA's file (`sub_CST`).
    df <- data.frame(sub_CST = rownames(r$centroids), r$centroids, check.names = FALSE)
    expect_identical(suppressWarnings(valencia(r$counts, centroids = df)), v)

    # The warnings of the synthetic case: an all-zero sample and a sample
    # that shares no taxa with the centroids.
    r <- ref$synthetic
    expect_warning(valencia(r$counts, totals = r$totals), "no reads.*'zero'")
    expect_warning(valencia(r$counts, totals = r$totals), "no taxa with the centroids.*first, I-A.*'unmatched'")

    # Sparse tables give the same results as dense ones, and threads change
    # nothing.
    v <- suppressWarnings(valencia(r$counts, totals = r$totals, cpus = 1L))
    for (cpus in 2:max_cpus)
      expect_identical(suppressWarnings(valencia(r$counts, totals = r$totals, cpus = cpus)), v, info = cpus)
    expect_identical(suppressWarnings(valencia(t(r$counts), totals = r$totals, margin = 2L)), v)

    if (requireNamespace('Matrix', quietly = TRUE)) {
      sp <- Matrix::Matrix(r$counts, sparse = TRUE)
      expect_inherits(sp, 'dgCMatrix')
      expect_identical(suppressWarnings(valencia(sp, totals = r$totals, cpus = 2L)), v)
      expect_identical(suppressWarnings(valencia(Matrix::t(sp), totals = r$totals, margin = 2L)), v)
      expect_identical(suppressWarnings(valencia(methods::as(sp, 'TsparseMatrix'), totals = r$totals)), v)
    }
    if (requireNamespace('slam', quietly = TRUE)) {
      st <- slam::as.simple_triplet_matrix(r$counts)
      expect_identical(suppressWarnings(valencia(st, totals = r$totals)), v)
    }
  }



  # The result. ====

  v <- valencia(counts)

  expect_inherits(v, c('ecodive_valencia', 'ecodive_types', 'data.frame'))
  expect_identical(names(v), c('sample', 'type', 'cst', 'score', 'matched'))
  expect_identical(rownames(v), rownames(counts))
  expect_identical(v$sample, rownames(counts))
  expect_identical(levels(v$type), subs)
  expect_identical(levels(v$cst), c('I', 'II', 'III', 'IV-A', 'IV-B', 'IV-C', 'V'))
  expect_identical(as.character(v$type), rep(subs, each = 9))
  cst_of <- c('I', 'I', 'II', 'III', 'III', 'IV-A', 'IV-B', rep('IV-C', 5), 'V')
  expect_identical(as.character(v$cst), rep(cst_of, each = 9))
  expect_true(all(v$score > 0 & v$score <= 1))

  sc <- attr(v, 'scores')
  expect_identical(dimnames(sc), list(rownames(counts), subs))
  expect_equal(sc, yue_clayton(counts, cents), tolerance = 1e-13)
  expect_identical(v$score, sc[cbind(seq_len(nrow(sc)), as.integer(v$type))])
  expect_identical(as.integer(v$type), unname(apply(sc, 1L, which.max)))

  # matched: the share of reads in the centroids' taxa.
  inc <- colnames(counts) %in% colnames(cents)
  expect_equal(v$matched, unname(rowSums(counts[, inc]) / rowSums(counts)), tolerance = 1e-14)

  types <- attr(v, 'types')
  expect_identical(names(types), c(
    'type', 'n_samples', 'top_feature', 'cst', 'score_mean',
    'reference_score_mean', 'reference_score_sd' ))
  expect_identical(types$type, factor(subs, levels = subs))
  expect_identical(types$n_samples, rep(9L, 13))
  expect_identical(types$top_feature[1:7], c(
    'Lactobacillus_crispatus', 'Lactobacillus_crispatus', 'Lactobacillus_gasseri',
    'Lactobacillus_iners', 'Lactobacillus_iners', 'BVAB1', 'Gardnerella_vaginalis' ))
  expect_identical(as.character(types$cst), cst_of)
  expect_equal(types$score_mean, as.vector(tapply(v$score, v$type, mean)), tolerance = 1e-15)
  expect_identical(types$reference_score_mean[c(1, 13)], c(0.995678, 0.734973))
  expect_identical(types$reference_score_sd[c(1, 13)], c(0.007615, 0.140276))
  expect_identical(attr(v, 'profiles'), cents)
  expect_null(attr(v, 'selection'))
  expect_identical(attr(v, 'method'), 'valencia')
  expect_null(attr(v, 'display'))
  expect_identical(plot_title(v, community_types_row(v)$name), 'VALENCIA Community State Types')

  # The types as groups in a test (unused levels are dropped).
  pm <- permanova(bray(counts[1:36, ]), v$type[1:36], permutations = 9)
  expect_true(is.finite(pm$statistic))

  # Threads.
  big <- counts[rep(seq_len(nrow(counts)), 2), ]
  rownames(big) <- sprintf('B%03d', seq_len(nrow(big)))
  vb <- valencia(big, cpus = 1L)
  for (cpus in 2:max_cpus)
    expect_identical(valencia(big, cpus = cpus), vb, info = cpus)

  # Sparse input.
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(valencia(Matrix::Matrix(counts, sparse = TRUE)), v)



  # Hand-made cases. ====

  # A sample identical to a centroid scores 1, against that sub-CST.
  vc <- valencia(cents, totals = rep(1, 13))
  expect_identical(as.character(vc$type), subs)
  expect_identical(vc$score, rep(1, 13))
  vc <- valencia(cents)
  expect_identical(as.character(vc$type), subs)
  expect_equal(vc$score, rep(1, 13), tolerance = 1e-6)

  # Proportions give the counts' results.
  expect_equal(attr(valencia(counts / rowSums(counts)), 'scores'), sc, tolerance = 1e-14)

  # Feature order and a feature absent everywhere change nothing; samples
  # permute with the input.
  perm <- sample(ncol(counts))
  vp <- valencia(cbind(counts[, perm], Nothing = 0L))
  expect_equal(attr(vp, 'scores'), sc, tolerance = 1e-14)
  expect_identical(vp$type, v$type)
  o  <- sample(nrow(counts))
  vo <- valencia(counts[o, ])
  expect_identical(vo$type, v$type[o])
  expect_identical(attr(vo, 'scores'), sc[o, ])

  # One sample, as a named vector or a one-row table.
  v1 <- valencia(counts[5, ])
  expect_identical(v1$type, v$type[5])
  expect_identical(unname(attr(v1, 'scores')), unname(sc[5, , drop = FALSE]))
  expect_identical(valencia(counts[5, , drop = FALSE])$score, v$score[5])

  # totals: by name or position; above the sum lowers the proportions.
  tt <- rowSums(counts) + 1000
  vt <- valencia(counts, totals = tt)
  expect_equal(attr(vt, 'scores'), yue_clayton(counts, cents, tt), tolerance = 1e-13)
  expect_identical(valencia(counts, totals = tt[sample(length(tt))]), vt)
  expect_identical(valencia(counts, totals = unname(tt)), vt)
  expect_equal(vt$matched, unname(v$matched * rowSums(counts) / tt), tolerance = 1e-14)
  expect_identical(valencia(counts, totals = rowSums(counts)), v)

  # Ties: two identical centroids; the first wins.
  tc <- cents[c('III-B', 'III-A', 'I-A'), ]
  rownames(tc) <- c('B', 'A', 'C')
  tc['A', ] <- tc['B', ]
  vt <- valencia(counts[28:45, ], centroids = tc)
  expect_true(all(attr(vt, 'scores')[, 'A'] == attr(vt, 'scores')[, 'B']))
  expect_false(any(vt$type == 'A'))

  # Centroids of other names: their own CSTs, no reference scores.
  expect_identical(levels(vt$cst), c('B', 'A', 'C'))
  expect_true(all(is.na(attr(vt, 'types')$reference_score_mean)))
  expect_identical(valencia_colours(c('I-A', 'X', 'V', 'Y')),
                   c('#ff6868', plot_palette(2)[1], '#ffff71', plot_palette(2)[2]))



  # Edge cases and warnings. ====

  z <- rbind(counts[1:3, ], Z = 0L, U = 0L)
  z <- cbind(z, Mine = c(0L, 0L, 0L, 0L, 7L))
  expect_warning(vz <- valencia(z), "1 sample has no reads.*'Z'")
  expect_warning(valencia(z), "1 sample shares no taxa.*'U'")
  vz <- suppressWarnings(valencia(z))
  expect_true(is.na(vz$type[4]) && is.na(vz$cst[4]) && is.nan(vz$score[4]) && is.na(vz$matched[4]))
  expect_true(all(is.nan(attr(vz, 'scores')[4, ])))
  expect_identical(as.character(vz$type[5]), 'I-A')
  expect_identical(vz$score[5], 0)
  expect_identical(vz$matched[5], 0)
  expect_identical(attr(vz, 'types')$n_samples[1], 4L)
  out <- capture.output(print(vz))
  expect_identical(out[[2]], '# 5 samples in 13 types (1 with samples); 1 without a type')

  # A sample of zeros with a total scores 0 (no shared taxa).
  expect_warning(vz2 <- valencia(z[4, , drop = FALSE], totals = 10), 'no taxa')
  expect_identical(vz2$score, 0)

  # Missing key taxa.
  few <- counts[, setdiff(colnames(counts), c('Lactobacillus_iners', 'BVAB1'))]
  expect_warning(valencia(few), "no feature named 'Lactobacillus_iners', 'BVAB1', which dominate")
  expect_identical(valencia_key_taxa(cents), c(
    'Lactobacillus_crispatus', 'Lactobacillus_gasseri', 'Lactobacillus_iners', 'BVAB1',
    'Gardnerella_vaginalis', 'g_Streptococcus', 'g_Enterococcus', 'g_Bifidobacterium',
    'g_Staphylococcus', 'Lactobacillus_jensenii' ))

  # Errors.
  expect_error(valencia(unname(counts)), 'must name its features')
  dup <- counts
  colnames(dup)[2] <- colnames(dup)[1]
  expect_error(valencia(dup), "duplicated feature names: '")
  na <- counts
  na[1, 1] <- NA
  expect_error(valencia(na), 'NA')
  neg <- counts
  neg[1, 1] <- -1L
  expect_error(valencia(neg), 'negative')
  expect_error(valencia(counts, margin = 3), '`margin`')
  expect_error(valencia(counts, cpus = 0), '`cpus`')
  expect_error(valencia(counts, totals = 1:3), 'one value per sample')
  expect_error(valencia(counts, totals = rep(NA_real_, nrow(counts))), 'finite')
  expect_error(valencia(counts, totals = replace(tt, 1, 0)), "positive.*'S001'")
  expect_warning(valencia(counts, totals = replace(tt, 2, 10)), "smaller than the sum of counts for 1 sample: 'S002'")
  bad <- tt
  names(bad)[1] <- 'nope'
  expect_error(valencia(counts, totals = bad), 'S001')
  expect_error(valencia(counts, centroids = unname(cents)), 'unique row names')
  expect_error(valencia(counts, centroids = cents[, c(1, 1, 2)]), 'unique column names')
  expect_error(valencia(counts, centroids = 'x'), 'numeric matrix')
  expect_error(valencia(counts, centroids = -cents), 'non-negative')
  expect_error(valencia(counts, centroids = rbind(cents, Z = 0)), "centroids of zeros: 'Z'")



  # print() and plot(). ====

  out <- capture.output(print(v))
  expect_identical(out[[1]], '# VALENCIA Community State Types')
  expect_identical(out[[2]], sprintf('# %d samples in 13 types', nrow(counts)))
  expect_true(any(out == '# Types:'))
  expect_false(any(out == '# Selection:'))
  expect_true(any(grepl('^# Samples: a data.frame: 117 . 5$', out)))
  expect_true(any(grepl('^# Attributes: types, profiles, scores, method$', out)))
  out <- capture.output(print(valencia(counts[1:20, ])))
  expect_identical(out[[2]], '# 20 samples in 13 types (3 with samples)')

  grDevices::pdf(NULL)
  expect_identical(plot(v), v)
  expect_identical(plot(v, 'scores'), v)
  expect_identical(plot(v, 'scores', main = 'Mine', ylim = c(0.5, 1)), v)
  expect_identical(plot(vt, 'scores'), vt)
  expect_identical(plot(vz, 'scores'), vz)
  expect_error(plot(v, 'selection'), 'no criterion')
  grDevices::dev.off()

#})
