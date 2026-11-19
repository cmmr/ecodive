#test_that("qmp_scale() and qmp_rarefy()", {

  exit_if_not(requireNamespace("Matrix", quietly = TRUE))
  exit_if_not(requireNamespace("slam",   quietly = TRUE))

  if (packageVersion("Matrix") >= "1.5") {
    COMPRESSED <- "CsparseMatrix"
    TRIPLET    <- "TsparseMatrix"
    UNPACKED   <- "unpackedMatrix"
  } else {
    COMPRESSED <- "dgCMatrix"
    TRIPLET    <- "dgTMatrix"
    UNPACKED   <- "dgeMatrix"
  }

  classes <- list(
    matrix = identity,
    dgC    = function (m) as(m, COMPRESSED),
    dgT    = function (m) as(m, TRIPLET),
    dge    = function (m) as(m, UNPACKED),
    slam   = slam::as.simple_triplet_matrix )

  # A 30 x 40 table with zeros and a few deep cells.
  x <- readRDS(file.path('expected', 'rarefy.rds'))$counts
  n <- rowSums(x)
  set.seed(11)
  L <- stats::setNames(exp(stats::rnorm(nrow(x), 23, 1)), rownames(x))

  # Base-R references.
  ref_scale <- function (x, L) x / rowSums(x) * L
  ref_target <- function (x, L) {
    t <- min(rowSums(x) / L) * L
    r <- round(t)
    t <- ifelse(abs(t - r) <= 1e-9 * pmax(1, r), r, floor(t))
    pmin(t, rowSums(x))
  }

  strip <- function (m) { attributes(m) <- attributes(m)[c('dim', 'dimnames')]; m }

  rel_diff <- function (a, b) max(abs(a - b) / pmax(abs(b), 1e-300))


  # ---------------------------------------------------------------------------
  # qmp_scale(): x / n * L, for every class and both margins.
  # ---------------------------------------------------------------------------

  expected <- ref_scale(x, L)
  for (cl in names(classes)) {
    for (margin in 1:2) {
      m <- if (margin == 1L) x else t(x)
      for (cpus in c(1L, 2L, 4L)) {
        q <- qmp_scale(classes[[cl]](m), L, margin = margin, cpus = cpus)
        expect_identical(class(q), class(classes[[cl]](m)), info = cl)
        qm <- as.matrix(q)
        if (margin == 2L) qm <- t(qm)
        expect_identical(dimnames(qm), dimnames(x), info = cl)
        expect_true(rel_diff(qm, expected) < 1e-15, info = paste(cl, margin, cpus))
        if (cpus == 1L) q1 <- q else expect_identical(q, q1, info = paste(cl, margin, cpus))
      }
    }
  }
  q <- qmp_scale(x, L)
  expect_equal(rowSums(q), L, tolerance = 1e-14)

  # Sparse stays sparse: zeros are not stored.
  q <- qmp_scale(classes$dgC(x), L)
  expect_identical(length(q@x), sum(x != 0))
  q <- qmp_scale(classes$slam(x), L)
  expect_identical(length(q$v), sum(x != 0))

  # Integer and logical input.
  xi <- x; storage.mode(xi) <- 'integer'
  expect_identical(qmp_scale(xi, L), qmp_scale(x, L))
  xl <- x > 0
  expect_equal(qmp_scale(xl, L), ref_scale(xl * 1, L), tolerance = 1e-15)

  # Proportions or any non-negative values: the same formula.
  xp <- x / rowSums(x) * 0.3
  expect_equal(qmp_scale(xp, L), expected, tolerance = 1e-14)

  # load by name, by position, and as a data column all agree.
  o <- rev(seq_along(L))
  expect_identical(qmp_scale(x, L[o]), qmp_scale(x, L))
  expect_identical(qmp_scale(x, unname(L)), qmp_scale(x, L))
  meta <- data.frame(cells = L[o], grp = 'a', row.names = names(L)[o])
  expect_identical(qmp_scale(x, 'cells', data = meta), qmp_scale(x, L))
  expect_identical(qmp_scale(x, L[o], data = meta), qmp_scale(x, L))
  meta2 <- data.frame(cells = unname(L))
  expect_identical(qmp_scale(x, 'cells', data = meta2), qmp_scale(x, L))
  # Extra rows in data are ignored.
  meta3 <- rbind(meta, data.frame(cells = 1, grp = 'b', row.names = 'other'))
  expect_identical(qmp_scale(x, 'cells', data = meta3), qmp_scale(x, L))
  # Unnamed samples: by position.
  expect_identical(qmp_scale(unname(x), L), unname(qmp_scale(x, L)))


  # ---------------------------------------------------------------------------
  # qmp_rarefy(): the deterministic parts exactly.
  # ---------------------------------------------------------------------------

  tgt <- ref_target(x, L)
  q   <- qmp_rarefy(x, L, seed = 5)
  expect_identical(attr(q, 'sampling_depth'), n / L)
  expect_identical(attr(q, 'target'), stats::setNames(as.integer(tgt), rownames(x)))
  expect_equal(rowSums(q), L, tolerance = 1e-14)
  expect_true(all(attr(q, 'target') <= n))

  # The sample with the smallest sampling depth keeps all its reads.
  i <- which.min(n / L)
  expect_identical(attr(q, 'target')[[i]], as.integer(n[[i]]))
  expect_true(rel_diff(q[i, ], expected[i, ]) < 1e-15)

  # The same as rarefy() to each sample's target, scaled by target and load.
  r <- rarefy(x, depth = tgt, seed = 5, drop = FALSE, warn = FALSE)
  expect_identical(unname(rowSums(r)), unname(tgt))
  expect_true(rel_diff(unclass(q)[, ], r / tgt * L) < 1e-15)

  # Bit-identical for cpus 1, 2 and 4, and for the same seed; differs for
  # another seed.
  for (cpus in c(2L, 4L))
    expect_identical(qmp_rarefy(x, L, seed = 5, cpus = cpus), q)
  expect_identical(qmp_rarefy(x, L, seed = 5), q)
  expect_false(identical(unclass(qmp_rarefy(x, L, seed = 6))[, ], unclass(q)[, ]))

  # Every class and margin (slam's draws are seeded by 1-based sample number,
  # as in rarefy(), so they differ from the others').
  for (cl in names(classes)) {
    for (margin in 1:2) {
      m  <- if (margin == 1L) x else t(x)
      for (cpus in c(1L, 2L, 4L)) {
        qc <- qmp_rarefy(classes[[cl]](m), L, margin = margin, seed = 5, cpus = cpus)
        if (cpus == 1L) q1 <- qc else expect_identical(qc, q1, info = paste(cl, margin, cpus))
      }
      expect_inherits(qc, class(classes[[cl]](m))[[1]], info = cl)
      expect_identical(attr(qc, 'target'), attr(q, 'target'), info = cl)
      expect_identical(attr(qc, 'sampling_depth'), attr(q, 'sampling_depth'), info = cl)
      qm <- strip(as.matrix(qc))
      if (margin == 2L) qm <- t(qm)
      expect_equal(rowSums(qm), L, tolerance = 1e-14, info = cl)
      if (cl != 'slam')
        expect_identical(qm, strip(q), info = paste(cl, margin))
      else {
        rs <- as.matrix(rarefy(
          classes$slam(m), depth = tgt, seed = 5, margin = margin, drop = FALSE, warn = FALSE ))
        if (margin == 2L) rs <- t(rs)
        expect_true(rel_diff(qm, rs / tgt * L) < 1e-15)
      }
    }
  }
  qs <- qmp_rarefy(classes$dgC(x), L, seed = 5)
  expect_true(length(qs@x) <= sum(x != 0))
  expect_false(any(qs@x == 0))

  # The draws are unbiased: the mean over many seeds approaches x / n * L.
  # Each r_ij is hypergeometric, so the z-scores of the means are ~N(0, 1).
  xs <- x[1:6, 1:8]
  Ls <- L[1:6]
  ts <- ref_target(xs, Ls)
  K  <- 400
  qs <- vapply(seq_len(K), function (s) unclass(qmp_rarefy(xs, Ls, seed = s))[, ], xs)
  mu <- ref_scale(xs, Ls)
  p  <- xs / rowSums(xs)
  ns <- rowSums(xs)
  v  <- ts * p * (1 - p) * (ns - ts) / pmax(ns - 1, 1) * (Ls / ts)^2
  z  <- (apply(qs, 1:2, mean) - mu) / sqrt(v / K)
  z  <- z[v > 0]
  expect_true(max(abs(z)) < 5)
  expect_true(abs(mean(z)) < 0.5)
  expect_true(stats::sd(z) > 0.6 && stats::sd(z) < 1.4)

  # Non-integer targets are floored; a hair below a whole number is that
  # number, so the smallest-depth sample keeps all its reads.
  y  <- matrix(c(3, 7, 10,  4, 4, 2,  5, 1, 1), 3, byrow = TRUE)
  Ly <- c(1, 0.8, 0.6)       # depths 20, 12.5, 11.67; min 11.67
  qy <- qmp_rarefy(y, Ly)
  expect_identical(attr(qy, 'target'), c(11L, 9L, 7L))
  expect_identical(attr(qy, 'sampling_depth'), rowSums(y) / Ly)
  # (n / L) * L can come out a hair below n; the target is still n.
  y2 <- matrix(c(49, 0, 0, 98), 2)
  L2 <- c(0.7, 0.7 * 2)
  expect_identical(attr(qmp_rarefy(y2, L2), 'target'), c(49L, 98L))
  found <- FALSE
  for (k in 1:2000) {
    nk <- 10 + k; Lk <- 0.1 * k + 0.07
    if (nk / Lk * Lk < nk) { found <- TRUE; break }
  }
  if (found) {
    yk <- matrix(c(nk, 0, 0, 3 * nk), 2)
    tk <- attr(qmp_rarefy(yk, c(Lk, 2 * Lk)), 'target')
    expect_identical(tk[[1]], as.integer(nk))
  }

  # The authors' example table (QMP.R), with cell counts per sample.
  a <- matrix(
    c(4,4,2,1, 8,5,2,0, 3,5,3,1, 10,8,3,0, 0,6,4,3), nrow = 5, byrow = TRUE,
    dimnames = list(paste('Sample', LETTERS[1:5]), paste0('taxa', 1:4)) )
  b <- c(10, 20, 34, 21, 12) * 1e5
  qa <- qmp_rarefy(a, b)
  expect_equal(unname(attr(qa, 'sampling_depth')), unname(rowSums(a)) / b)
  expect_identical(unname(attr(qa, 'target')), as.integer(ref_target(a, b)))
  expect_equal(unname(rowSums(qa)), b)


  # ---------------------------------------------------------------------------
  # The authors' QMP.R (data-raw/expected/qmp.r): its rarefaction is
  # phyloseq's, so the draws differ; everything else must agree.
  # ---------------------------------------------------------------------------

  fx <- readRDS(file.path('expected', 'qmp.rds'))

  n_rescued <- 0L
  for (nm in names(fx$cases)) {
    cs <- fx$cases[[nm]]
    q  <- qmp_rarefy(cs$x, cs$load)

    # Sampling depths: the script's n / L.
    expect_identical(unname(attr(q, 'sampling_depth')), unname(rowSums(cs$x)) / cs$load, info = nm)

    # Targets: the script asks phyloseq for s_min * L, which it truncates.
    # ecodive floors too, except where s_min * L is a hair below a whole
    # number, where the script loses a read.
    expect_equal(cs$asked, min(rowSums(cs$x) / cs$load) * cs$load, tolerance = 1e-15, info = nm)
    tg   <- unname(attr(q, 'target'))
    near <- abs(cs$asked - round(cs$asked)) <= 1e-9 * pmax(1, round(cs$asked))
    expect_identical(tg[!near], as.integer(cs$got[!near]), info = nm)
    expect_identical(tg[near], as.integer(round(cs$asked[near])), info = nm)
    expect_true(all(tg - cs$got %in% 0:1), info = nm)
    n_rescued <- n_rescued + sum(tg != cs$got)

    # Scaling: the script's output is its rarefied counts / their total *
    # load, which qmp_scale() of those counts gives.
    r <- cs$qmp / cs$load * cs$got
    expect_equal(r, round(r), tolerance = 1e-12, info = nm)
    expect_equal(unname(rowSums(round(r))), cs$got, info = nm)
    expect_true(rel_diff(strip(qmp_scale(round(r), cs$load)), cs$qmp) < 1e-15, info = nm)
    expect_equal(unname(rowSums(q)), cs$load, tolerance = 1e-14, info = nm)
  }
  # The smallest-depth samples of 'random' and 'hair' keep a read the script
  # loses, and 'hair''s second sample too.
  expect_identical(n_rescued, 3L)

  # The mean over many draws: ecodive's and the script's agree, as both
  # approach x / n * load.
  md <- fx$mean_draws
  ts <- ref_target(md$x, md$load)
  qs <- vapply(seq_len(md$K), function (s) strip(qmp_rarefy(md$x, md$load, seed = s)), md$x * 1)
  p  <- md$x / rowSums(md$x)
  ns <- rowSums(md$x)
  v  <- ts * p * (1 - p) * (ns - ts) / pmax(ns - 1, 1) * (md$load / ts)^2
  z  <- (apply(qs, 1:2, mean) - md$mean) / sqrt(2 * v / md$K)
  z  <- z[v > 0]
  expect_true(max(abs(z)) < 5)
  expect_true(stats::sd(z) > 0.6 && stats::sd(z) < 1.4)


  # ---------------------------------------------------------------------------
  # Errors.
  # ---------------------------------------------------------------------------

  expect_error(qmp_scale(x), '`load` is required')
  expect_error(qmp_rarefy(x), '`load` is required')
  expect_error(qmp_scale(x, NULL), '`load` is required')
  expect_error(qmp_scale(x, unname(L)[-1]), '`load` has 29 values, but `counts` has 30 samples')
  expect_error(qmp_scale(x, L[-1]), "1 sample in `counts` has no value in `load`: 'S01'")
  expect_error(qmp_scale(x, 5), '`load` has 1 value, but `counts` has 30 samples')
  bad <- L; bad[c(3, 7)] <- NA
  expect_error(qmp_scale(x, bad), "not for 2 samples: 'S03', 'S07'")
  bad <- L; bad[3] <- 0
  expect_error(qmp_scale(x, bad), "`load` must be finite and greater than 0, but is not for 1 sample: 'S03'")
  bad <- L; bad[4] <- -1
  expect_error(qmp_rarefy(x, bad), "not for 1 sample: 'S04'")
  bad <- L; bad[4] <- Inf
  expect_error(qmp_scale(x, bad), "not for 1 sample: 'S04'")
  expect_error(qmp_scale(x, as.character(L)), '`load` must be numeric')
  expect_error(qmp_scale(x, factor(L)), '`load` must be numeric')
  expect_error(qmp_scale(x, matrix(L)), '`load` must be a vector with one value per sample')
  expect_error(qmp_scale(x, stats::setNames(L, paste0('X', 1:30))), 'samples in `counts` have no value in `load`')
  expect_error(qmp_scale(x, 'cells'), '`load` names a column, but `data` is not given')
  expect_error(qmp_scale(x, 'nope', data = meta), "`load` must be the name of a column of `data`, but 'nope' is not one")
  expect_error(qmp_scale(x, 'grp', data = meta), "`data` column 'grp' must be numeric")
  meta_na <- meta; meta_na$cells[2] <- NA
  expect_error(qmp_scale(x, 'cells', data = meta_na), "`data` column 'cells' must be finite and greater than 0")
  expect_error(qmp_scale(x, 'cells', data = meta[-1, ]), "1 sample in `counts` has no row in `data`")
  expect_error(qmp_scale(x, 'cells', data = list(cells = L)), '`data` must be a data.frame')

  xe <- x; xe[c(2, 5), ] <- 0
  expect_error(qmp_scale(xe, L), "2 samples in `counts` are empty \\(no counts\\): 'S02', 'S05'")
  expect_error(qmp_rarefy(xe, L), "2 samples in `counts` are empty")
  expect_error(qmp_scale(classes$dgC(xe), L), "'S02', 'S05'")

  expect_error(qmp_rarefy(x + 0.5, L), '`counts` must be whole numbers')
  expect_error(qmp_rarefy(classes$dgC(x + 0.5), L), '`counts` must be whole numbers')
  expect_error(qmp_scale(-x, L), '`counts` must not contain negative values')
  xn <- x; xn[1] <- NA
  expect_error(qmp_scale(xn, L), '`counts` must not contain NA')
  expect_error(qmp_scale(x, L, margin = 3), '`margin` must be 1 or 2')
  expect_error(qmp_scale(x, L, cpus = 0), '`cpus` must be an integer')
  expect_error(qmp_rarefy(x, L, seed = 0.5), '`seed` must be an integer')

  # A load so small that the sample would keep no reads.
  Lz <- L; Lz[[9]] <- Lz[[9]] * 1e-6
  expect_error(qmp_rarefy(x, Lz), "1 sample would be rarefied to 0 reads \\('S09'\\)")

  # n = 1 and 2 samples.
  expect_equal(unclass(qmp_rarefy(x[1, , drop = FALSE], L[1]))[, ], expected[1, ], tolerance = 1e-15)
  expect_equal(qmp_scale(x[1:2, ], L[1:2]), expected[1:2, ], tolerance = 1e-15)

#})
