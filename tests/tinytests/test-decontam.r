#test_that("decontam", {

  # decontam() and decontam_not(): against decontam's isContaminant() and
  # isNotContaminant() (Artistic-2.0), stored by
  # data-raw/expected/decontam.r; the prevalence score against base R's
  # prop.test() and fisher.test(); and the rules for inputs, threads and
  # table classes.

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))
  exit_if_not(requireNamespace("slam",   quietly = TRUE))

  f <- readRDS(file.path('expected', 'decontam.rds'))
  q <- function (expr) suppressWarnings(expr)

  # The largest relative difference, where both are numbers; and whether the
  # NAs agree.
  rel_diff <- function (a, b) {
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) return (0)
    max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
  }
  same_na <- function (a, b) identical(is.na(a), is.na(b))
  strip   <- function (x) { attributes(x) <- attributes(x)[c('names', 'row.names', 'class')]; x }


  #--------------------------------------------------
  # isContaminant(): every method, batch combination and
  # threshold form; scores at 1e-12, calls and
  # prevalences identical.
  #--------------------------------------------------
  ref  <- f$contaminant
  keys <- unique(ref[, c('case', 'method', 'batch_combine', 'batched', 'threshold', 'normalize', 'use')])
  expect_true(nrow(keys) >= 100)

  worst <- 0
  for (i in seq_len(nrow(keys))) {
    k   <- keys[i, ]
    r   <- ref[ref$case == k$case & ref$method == k$method & ref$batch_combine == k$batch_combine &
                 ref$batched == k$batched & ref$threshold == k$threshold &
                 ref$normalize == k$normalize & ref$use == k$use, ]
    t   <- f$tables[[k$case]]
    X   <- t$X
    if (!k$normalize) { X <- X / rowSums(X) * 1000; X[!is.finite(X)] <- 0 }
    use <- strsplit(k$use, ' ')[[1]]
    thr <- as.numeric(strsplit(k$threshold, ' ')[[1]])
    res <- q(decontam(
      X, conc = if ('conc' %in% use) t$conc, neg = if ('neg' %in% use) t$neg,
      method = k$method, batch = if (k$batched) t$batch, batch_combine = k$batch_combine,
      threshold = thr, norm = if (k$normalize) 'percent' else 'none', cpus = 2L ))
    info <- paste(unlist(k), collapse = ' ')

    expect_identical(res$feature, r$feature, info = info)
    expect_identical(res$prev, r$prev, info = info)
    expect_identical(res$contaminant, r$contaminant, info = info)
    for (col in c('freq', 'p_freq', 'p_prev', 'p')) {
      expect_true(same_na(res[[col]], r[[col]]), info = paste(info, col))
      d <- rel_diff(res[[col]], r[[col]])
      worst <- max(worst, d)
      expect_true(d < 1e-12, info = paste(info, col))
    }
  }
  expect_true(worst < 1e-12)


  #--------------------------------------------------
  # isNotContaminant().
  #--------------------------------------------------
  ref <- f$not_contaminant
  for (nm in unique(ref$case)) for (thr in unique(ref$threshold)) {
    r   <- ref[ref$case == nm & ref$threshold == thr, ]
    t   <- f$tables[[nm]]
    res <- q(decontam_not(t$X, t$neg, threshold = as.numeric(thr), cpus = 2L))
    info <- paste(nm, thr)
    expect_identical(res$feature, r$feature, info = info)
    expect_identical(res$prev, r$prev, info = info)
    expect_identical(res$not_contaminant, r$not_contaminant, info = info)
    for (col in c('freq', 'p_prev', 'p')) {
      expect_true(same_na(res[[col]], r[[col]]), info = paste(info, col))
      expect_true(rel_diff(res[[col]], r[[col]]) < 1e-12, info = paste(info, col))
    }
    expect_identical(names(res), c(
      'feature', 'freq', 'prev', 'prev_neg', 'prev_pos', 'p_prev', 'p', 'not_contaminant' ))
    expect_identical(class(res), c('ecodive_decontam_not', 'data.frame'))
  }


  #--------------------------------------------------
  # The prevalence score on every 2 x 2 table with up to
  # 10 controls and 10 samples: each feature of a table
  # is one (a, c), present in `a` of the `n1` controls
  # and `c` of the `n2` samples. Against prop.test() and,
  # where it warns, fisher.test()'s mid-p (decontam's
  # rule); 0.5 when present everywhere, NA in under two.
  #--------------------------------------------------
  fisher_mid <- function (tab) {
    g <- stats::fisher.test(tab, alternative = 'greater')$p.value
    l <- stats::fisher.test(tab, alternative = 'less')$p.value
    g - (g + l - 1) / 2
  }
  n_warned <- 0L
  worst    <- 0
  for (n1 in c(1:4, 7, 10)) for (n2 in c(1:4, 6, 10)) {
    ac <- expand.grid(a = 0:n1, c = 0:n2)
    X  <- matrix(0L, n1 + n2, nrow(ac))
    for (j in seq_len(nrow(ac))) {
      if (ac$a[j] > 0) X[seq_len(ac$a[j]), j] <- 1L
      if (ac$c[j] > 0) X[n1 + seq_len(ac$c[j]), j] <- 1L
    }
    X   <- cbind(X, 1L)   # no empty samples
    neg <- rep(c(TRUE, FALSE), c(n1, n2))
    res <- q(decontam(X, neg = neg, cpus = 1L))
    for (j in seq_len(nrow(ac))) {
      a <- ac$a[j]; c <- ac$c[j]
      tab <- matrix(c(a, c, n1 - a, n2 - c), 2)
      exp <- if (a + c < 2) NA_real_
             else if (a == n1 && c == n2) 0.5
             else tryCatch(
               stats::prop.test(tab, alternative = 'greater')$p.value,
               warning = function (w) { n_warned <<- n_warned + 1L; fisher_mid(tab) } )
      got <- res$p_prev[[j]]
      expect_identical(is.na(got), is.na(exp), info = paste(n1, n2, a, c))
      if (!is.na(exp)) {
        d <- abs(got - exp) / max(abs(exp), 1e-300)
        worst <- max(worst, d)
        expect_true(d < 1e-12, info = paste(n1, n2, a, c))
      }
    }
  }
  expect_true(n_warned > 100L)   # both branches are covered
  expect_true(worst < 1e-12)


  #--------------------------------------------------
  # Threads, table classes, orientation.
  #--------------------------------------------------
  t <- f$tables$main
  X <- t$X
  if (packageVersion("Matrix") >= "1.5") {
    COMPRESSED <- "CsparseMatrix"; TRIPLET <- "TsparseMatrix"; UNPACKED <- "unpackedMatrix"
  } else {
    COMPRESSED <- "dgCMatrix"; TRIPLET <- "dgTMatrix"; UNPACKED <- "dgeMatrix"
  }
  Xd <- X; storage.mode(Xd) <- 'double'
  base <- decontam(X, t$conc, t$neg, batch = t$batch, cpus = 1L)
  for (cpus in c(2L, 4L, 7L))
    expect_identical(decontam(X, t$conc, t$neg, batch = t$batch, cpus = cpus), base, info = cpus)
  expect_identical(decontam(Xd, t$conc, t$neg, batch = t$batch), base)
  expect_identical(decontam(as(as(Xd, 'CsparseMatrix'), COMPRESSED), t$conc, t$neg, batch = t$batch, cpus = 3L), base)
  expect_identical(decontam(as(as(Xd, 'CsparseMatrix'), TRIPLET),    t$conc, t$neg, batch = t$batch), base)
  expect_identical(decontam(as(Matrix::Matrix(Xd, sparse = FALSE), UNPACKED), t$conc, t$neg, batch = t$batch), base)
  expect_identical(decontam(slam::as.simple_triplet_matrix(Xd), t$conc, t$neg, batch = t$batch, cpus = 2L), base)
  expect_identical(decontam(t(X), t$conc, t$neg, margin = 2L, batch = t$batch), base)
  expect_identical(decontam(slam::as.simple_triplet_matrix(t(Xd)), t$conc, t$neg, margin = 2, batch = t$batch), base)

  base_not <- decontam_not(X, t$neg, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(decontam_not(X, t$neg, cpus = cpus), base_not)
  expect_identical(decontam_not(as(as(Xd, 'CsparseMatrix'), COMPRESSED), t$neg), base_not)
  expect_identical(decontam_not(t(X), t$neg, margin = 2L), base_not)

  # A wide table, so that the threads have work.
  set.seed(3)
  W <- matrix(rpois(60 * 3000, 0.7), 60, 3000)
  cw <- exp(stats::runif(60, 0, 3)); nw <- rep(c(FALSE, TRUE), c(50, 10)); bw <- rep(1:3, 20)
  w1 <- decontam(W, cw, nw, batch = bw, cpus = 1L)
  for (cpus in c(2L, 4L))
    expect_identical(decontam(W, cw, nw, batch = bw, cpus = cpus), w1)
  expect_identical(decontam(as(as(W * 1.0, 'CsparseMatrix'), COMPRESSED), cw, nw, batch = bw, cpus = 4L), w1)


  #--------------------------------------------------
  # The result.
  #--------------------------------------------------
  expect_identical(class(base), c('ecodive_decontam', 'data.frame'))
  expect_identical(names(base), c(
    'feature', 'freq', 'prev', 'prev_neg', 'prev_pos', 'p_freq', 'p_prev', 'p', 'contaminant' ))
  expect_identical(attr(base, 'method'), 'decontam')
  expect_identical(attr(base, 'test'), 'combined')
  expect_identical(attr(base, 'threshold'), 0.1)
  expect_identical(attr(base, 'batch_combine'), 'minimum')
  expect_identical(attr(base, 'n_neg'), 8L)
  expect_identical(attr(base, 'n_pos'), 32L)
  pfb <- attr(base, 'p_freq_batch')
  ppb <- attr(base, 'p_prev_batch')
  expect_identical(dimnames(pfb), list(c('run1', 'run2'), colnames(X)))
  expect_identical(dim(ppb), c(2L, ncol(X)))
  m <- q(apply(pfb, 2, min, na.rm = TRUE)); m[is.infinite(m)] <- NA
  expect_equal(unname(m), base$p_freq)
  expect_identical(base$prev_neg, as.integer(colSums(X[t$neg, ] > 0)))
  expect_identical(base$prev_pos, as.integer(colSums(X[!t$neg, ] > 0)))
  expect_identical(base$prev, as.integer(colSums(X > 0)))
  expect_true(is.integer(base$prev_neg) && is.integer(base$prev))

  # Without batches, no per-batch attributes; without neg, NA counts.
  r <- decontam(X, t$conc)
  expect_identical(attr(r, 'test'), 'frequency')
  expect_null(attr(r, 'p_freq_batch'))
  expect_true(all(is.na(r$prev_neg)) && all(is.na(r$p_prev)))
  expect_identical(attr(r, 'n_neg'), 0L)
  expect_identical(attr(r, 'n_pos'), 40L)
  expect_identical(attr(decontam(X, neg = t$neg), 'test'), 'prevalence')
  r <- decontam(X, t$conc, t$neg, method = 'either', threshold = 0.2)
  expect_identical(attr(r, 'threshold'), c(0.2, 0.2))
  expect_true(all(is.na(r$p)))

  # The removal one-liner.
  expect_identical(dim(X[, !base$contaminant]), c(40L, sum(!base$contaminant)))

  # Features without names are numbered.
  expect_identical(decontam(unname(X), unname(t$conc))$feature, as.character(seq_len(ncol(X))))


  #--------------------------------------------------
  # Per-sample values: by name, from `data`.
  #--------------------------------------------------
  meta <- data.frame(dna = t$conc, blank = t$neg, run = t$batch, row.names = rownames(X))
  expect_identical(decontam(X, 'dna', 'blank', meta, batch = 'run'), base)
  expect_identical(decontam(X, 'dna', 'blank', meta[40:1, ], batch = 'run'), base)
  o <- 40:1
  expect_identical(decontam(X, t$conc[o], t$neg[o], batch = t$batch[o]), base)
  expect_identical(decontam_not(X, 'blank', meta), base_not)
  expect_error(decontam(X, 'nope', data = meta), 'column of `data`')
  expect_error(decontam(X, 'dna'), '`data` is not given')
  expect_error(decontam(X, t$conc[1:39]), 'sample')


  #--------------------------------------------------
  # Edge cases.
  #--------------------------------------------------
  e <- f$tables$edge
  expect_warning(r <- decontam(e$X, e$conc, e$neg), 'Removed 1 sample with zero total')
  expect_true(is.na(r$p_freq[10]) && is.na(r$p_prev[10]) && !r$contaminant[10])  # one sample
  expect_true(is.na(r$p_freq[11]) && r$prev[11] == 0L && !r$contaminant[11])     # none
  expect_identical(r$p_prev[12], 0.5)                                            # every sample
  expect_identical(attr(r, 'n_pos'), 11L)

  # All samples controls, or none: no prevalence scores.
  r <- decontam(X, neg = rep(TRUE, 40))
  expect_true(all(is.na(r$p_prev)) && !any(r$contaminant))
  r <- decontam(X, neg = rep(FALSE, 40))
  expect_true(all(is.na(r$p_prev)))

  # NA in neg: those samples are left out of both tests (but counted in
  # freq and prev).
  na_neg <- t$neg
  na_neg[c(3, 37)] <- NA
  r1 <- decontam(X, t$conc, na_neg)
  r2 <- decontam(X[-c(3, 37), ], t$conc[-c(3, 37)], t$neg[-c(3, 37)])
  for (col in c('p_freq', 'p_prev', 'p', 'contaminant', 'prev_neg', 'prev_pos'))
    expect_identical(r1[[col]], r2[[col]], info = col)
  expect_identical(r1$prev, base$prev)
  expect_identical(attr(r1, 'n_neg') + attr(r1, 'n_pos'), 38L)

  # Batches: at least 2 samples; a warning at 4 or fewer.
  b <- t$batch; b[1] <- 'solo'
  expect_error(decontam(X, t$conc, batch = b), 'at least 2 samples')
  b <- t$batch; b[1:4] <- 'small'
  expect_warning(decontam(X, t$conc, batch = b), 'very few')
  b <- t$batch; b[2] <- NA
  expect_error(decontam(X, t$conc, batch = b), 'NA')
  # Unused factor levels are dropped.
  expect_identical(
    decontam(X, t$conc, t$neg, batch = factor(t$batch, levels = c('run1', 'run2', 'run3'))), base )

  # 'product' with no scores in any batch gives 1, as decontam.
  r <- decontam(e$X[-3, ], e$conc[-3], method = 'frequency', batch_combine = 'product')
  expect_identical(r$p_freq[11], 1)
  r <- decontam(e$X[-3, ], e$conc[-3], method = 'frequency', batch_combine = 'fisher')
  expect_equal(r$p_freq[11], 0.5)

  # Invalid arguments.
  cz <- t$conc; cz[5] <- 0
  expect_error(decontam(X, cz), 'greater than 0')
  cz[5] <- NA
  expect_error(decontam(X, cz), 'greater than 0')
  expect_error(decontam(X, as.character(t$conc)), 'numeric')
  expect_error(decontam(X, neg = as.integer(t$neg)), 'logical')
  expect_error(decontam(X, t$conc, method = 'prevalence'), '`neg` is required')
  expect_error(decontam(X, neg = t$neg, method = 'frequency'), '`conc` is required')
  expect_error(decontam(X), '`neg` is required')
  expect_error(decontam(X, t$conc, threshold = 1.5), 'from 0 to 1')
  expect_error(decontam(X, t$conc, threshold = -0.1), 'from 0 to 1')
  expect_error(decontam(X, t$conc, threshold = NA), 'from 0 to 1')
  expect_error(decontam(X, t$conc, threshold = c(0.1, 0.5)), 'single number')
  expect_error(decontam(X, t$conc, t$neg, method = 'both', threshold = c(0.1, 0.2, 0.3)), 'one or two')
  expect_error(decontam(X, t$conc, method = 'nope'), '`method` must be one of')
  expect_error(decontam(X, t$conc, batch_combine = 'mean'), '`batch_combine`')
  expect_error(decontam(X, t$conc, norm = 'clr'), '"percent" or "none"')
  expect_error(decontam(-X, t$conc), 'negative')
  expect_error(decontam_not(X), '`neg` is required')
  expect_error(decontam_not(X, t$neg, threshold = 2), 'from 0 to 1')
  expect_identical(attr(decontam(X, t$conc, method = 'freq'), 'test'), 'frequency')
  expect_error(decontam(matrix(0, 3, 3), c(1, 2, 3)), 'empty')


  #--------------------------------------------------
  # Printing and plots.
  #--------------------------------------------------
  out <- utils::capture.output(print(base))
  expect_true(any(grepl('decontam (combined)', out, fixed = TRUE)))
  out <- utils::capture.output(print(base_not))
  expect_true(any(grepl('decontam_not (prevalence)', out, fixed = TRUE)))

  drawn <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  grDevices::pdf(NULL)
  expect_true(drawn(base))
  expect_true(drawn(base, 'prevalence'))
  expect_true(drawn(base, 'prevalence', legend = NULL, main = 'x'))
  expect_true(drawn(decontam(X, t$conc, t$neg, method = 'both', threshold = c(0.1, 0.5))))
  expect_true(drawn(base_not))
  expect_true(drawn(base_not, 'prevalence'))
  expect_true(drawn(decontam(X, t$conc), breaks = 10))
  expect_error(plot(decontam(X, t$conc), 'prevalence'), 'needs negative controls')
  expect_error(plot(base, 'nope'), '`type`')
  expect_error(plot(q(decontam(e$X[-3, 10:11], e$conc[-3]))), 'no scores')

  v <- withVisible(decontam_frequency_plot(X, c('F01', 'F20'), t$conc, t$neg))
  expect_false(v$visible)
  expect_identical(names(v$value), c('feature', 'sample', 'conc', 'neg', 'abundance'))
  expect_identical(nrow(v$value), 80L)
  expect_equal(v$value$abundance[1:40], unname(X[, 1] / rowSums(X)))
  d1 <- decontam_frequency_plot(X, 1, 'dna', 'blank', meta, log = FALSE, legend = NULL)
  d2 <- decontam_frequency_plot(t(X), 'F01', t$conc, margin = 2L, models = FALSE, norm = 'none')
  expect_identical(d1$abundance, v$value$abundance[1:40])
  expect_identical(d2$abundance, as.double(X[, 1]))
  expect_warning(decontam_frequency_plot(e$X, 1:3, e$conc, e$neg), 'Removed 1 sample')
  expect_error(decontam_frequency_plot(X, 'nope', t$conc), 'not in `counts`')
  expect_error(decontam_frequency_plot(X, 1), '`conc` is required')
  expect_error(decontam_frequency_plot(X, 1, t$conc, models = NA), '`models`')
  grDevices::dev.off()

#})
