# Writes tests/tinytests/expected/longitudinal.rds: inputs and reference
# results for test-longitudinal.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/longitudinal.r      # writes q2's inputs
#     Q2L_DIR=/path/to/q2-longitudinal \
#       python data-raw/expected/longitudinal.py    # q2-longitudinal's outputs
#     Rscript data-raw/expected/longitudinal.r      # folds them in
#
# q2-longitudinal 2026.10.0 (BSD-3) is the reference; see the Python
# script for how its functions are run without a QIIME 2 install.
#
# Inputs (simulated, written to data-raw/expected/longitudinal/):
#
#   meta.csv     One row per sample (in a shuffled order): sample, subject,
#                state (days 0, 3, 7, 14, 28), group (A or B) and shannon
#                (hexadecimal floats; NA for one sample). Subjects S01-S12
#                (S01-S06 in group A) have a sample at every state, except
#                that S03 has none at day 7 (a skipped middle state), S09
#                none at day 0, and S05 two at day 14 (a replicate). S13 has
#                samples at days 0 and 7 only (two samples), S14 at day 3
#                only (one sample).
#   counts.csv   Samples x 10 features, integer counts, in meta.csv's order.
#                Feature F07 is 0 in every sample of S02 (constant within a
#                subject), and counts are small, so with ties.
#   rel.csv      counts / rowSums(counts), hexadecimal floats: q2's relative
#                frequencies (FeatureTable[RelativeFrequency]).
#   dm.csv       Bray-Curtis dissimilarities of rel.csv, a full matrix of
#                hexadecimal floats (sum |a - b| / sum (a + b) per pair).
#
# Fields of the .rds:
#
#   meta, counts, dm   The inputs (dm as a `dist`).
#   fd                 q2's _first_differences() of `shannon` (named list):
#                      drop (consecutive states, replicate_handling =
#                      'drop'), baseline (baseline = 0, 'drop'), baseline7
#                      (baseline = 7, 'drop'). Each a data frame: sample (the
#                      later sample, q2's index), difference, in q2's order.
#   fd_table           The same for each feature of rel.csv ('drop',
#                      consecutive): a list of data frames named by feature.
#   fdist              q2's _first_differences(distance_matrix = dm): drop,
#                      baseline (0). Data frames: sample, distance.
#   errors             q2's error messages (named character vector).
#   nmit               q2's _nmit() for corr_method kendall, pearson,
#                      spearman and dist_method fro, nuc: named list of
#                      matrices (`<corr>_<dist>`), with q2's labels (each
#                      subject's first sample) as dimnames, in q2's order.
#   nmit_ids           The subject of each of those labels.
#   versions           Of the Python packages used.

dir <- file.path('data-raw', 'expected', 'longitudinal')
dir.create(dir, showWarnings = FALSE)

hexf <- function (x) ifelse(is.na(x), 'NA', sprintf('%a', x))
unhex <- function (x) { x[x %in% c('NA', 'nan', '')] <- NA; as.numeric(x) }

if (!file.exists(file.path(dir, 'meta.csv'))) {

  set.seed(20261124)
  states <- c(0, 3, 7, 14, 28)
  rows   <- list()
  for (s in 1:12) {
    st <- states
    if (s == 3) st <- setdiff(st, 7)
    if (s == 9) st <- setdiff(st, 0)
    if (s == 5) st <- sort(c(st, 14))
    rows[[s]] <- data.frame(subject = sprintf('S%02d', s), state = st)
  }
  rows[[13]] <- data.frame(subject = 'S13', state = c(0, 7))
  rows[[14]] <- data.frame(subject = 'S14', state = 3)
  meta <- do.call(rbind, rows)
  meta$group <- ifelse(meta$subject %in% sprintf('S%02d', c(1:6, 13)), 'A', 'B')

  n <- nrow(meta); p <- 10
  base <- exp(rnorm(p, 2, 0.8))
  subj <- exp(matrix(rnorm(14 * p, 0, 0.5), 14, p))
  sid  <- as.integer(substr(meta$subject, 2, 3))
  mu   <- subj[sid, ] * rep(base, each = n) *
          exp(outer(meta$state / 28 * ifelse(meta$group == 'A', 1, -1), seq(-1, 1, length.out = p)))
  counts <- matrix(rnbinom(n * p, mu = mu, size = 2), n, p)
  counts[meta$subject == 'S02', 7] <- 0
  counts[rowSums(counts) == 0, 1] <- 1

  ord  <- sample(n)
  meta <- meta[ord, ]
  counts <- counts[ord, ]
  meta$sample <- sprintf('smp%02d', seq_len(n))
  dimnames(counts) <- list(meta$sample, sprintf('F%02d', 1:p))

  rel <- counts / rowSums(counts)
  meta$shannon <- -rowSums(ifelse(rel > 0, rel * log(rel), 0))
  meta$shannon[meta$subject == 'S08' & meta$state == 3] <- NA
  meta <- meta[, c('sample', 'subject', 'state', 'group', 'shannon')]
  rownames(meta) <- NULL

  bc <- matrix(0, n, n, dimnames = list(meta$sample, meta$sample))
  for (i in seq_len(n)) for (j in seq_len(n))
    bc[i, j] <- sum(abs(rel[i, ] - rel[j, ])) / sum(rel[i, ] + rel[j, ])

  m2 <- meta; m2$shannon <- hexf(m2$shannon)
  write.csv(m2, file.path(dir, 'meta.csv'), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(sample = meta$sample, counts, check.names = FALSE),
            file.path(dir, 'counts.csv'), row.names = FALSE, quote = FALSE)
  r2 <- rel; r2[] <- hexf(rel)
  write.csv(data.frame(sample = meta$sample, r2, check.names = FALSE),
            file.path(dir, 'rel.csv'), row.names = FALSE, quote = FALSE)
  b2 <- bc; b2[] <- hexf(bc)
  write.csv(data.frame(sample = meta$sample, b2, check.names = FALSE),
            file.path(dir, 'dm.csv'), row.names = FALSE, quote = FALSE)

  message('Inputs written; now run data-raw/expected/longitudinal.py.')
  quit(save = 'no')
}


#----------------------------------------------------------
# Fold q2-longitudinal's outputs in.
#----------------------------------------------------------
rd <- function (f) read.csv(file.path(dir, f), colClasses = 'character', check.names = FALSE)

meta <- rd('meta.csv')
meta$state   <- as.numeric(meta$state)
meta$shannon <- unhex(meta$shannon)

cc     <- rd('counts.csv')
counts <- as.matrix(cc[, -1]); storage.mode(counts) <- 'integer'
rownames(counts) <- cc$sample

dd <- rd('dm.csv')
dm <- matrix(unhex(as.matrix(dd[, -1])), nrow(dd), dimnames = list(dd$sample, dd$sample))
dm <- as.dist(dm)

pairs <- function (f, col) {
  x <- rd(f)
  data.frame(sample = x$sample, value = unhex(x[[col]]))
}

fd <- list(
  drop      = pairs('fd-drop.csv', 'difference'),
  baseline  = pairs('fd-baseline.csv', 'difference'),
  baseline7 = pairs('fd-baseline7.csv', 'difference') )

ft <- rd('fd-table.csv')
fd_table <- lapply(split(ft, factor(ft$feature, unique(ft$feature))), function (x)
  data.frame(sample = x$sample, value = unhex(x$difference)))

fdist <- list(
  drop     = pairs('fdist-drop.csv', 'distance'),
  baseline = pairs('fdist-baseline.csv', 'distance') )

er <- rd('errors.csv')
errors <- setNames(er$message, er$case)

nmit <- list()
for (cm in c('kendall', 'pearson', 'spearman')) for (dmeth in c('fro', 'nuc')) {
  x <- rd(sprintf('nmit-%s-%s.csv', cm, dmeth))
  nmit[[paste(cm, dmeth, sep = '_')]] <-
    matrix(unhex(as.matrix(x[, -1])), nrow(x), dimnames = list(x$id, x$id))
}
ids <- rd('nmit-ids.csv')
nmit_ids <- setNames(ids$subject, ids$id)

vv <- rd('versions.csv')
versions <- setNames(vv$version, vv$package)

saveRDS(
  list(meta = meta, counts = counts, dm = dm, fd = fd, fd_table = fd_table,
       fdist = fdist, errors = errors, nmit = nmit, nmit_ids = nmit_ids,
       versions = versions),
  file.path('tests', 'tinytests', 'expected', 'longitudinal.rds'),
  version = 2 )
