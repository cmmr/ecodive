# Writes tests/tinytests/expected/feature_volatility.rds: inputs and
# reference results for test-feature_volatility.r, so that the tests need no
# reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/feature_volatility.r     # writes q2's inputs
#     Q2L_DIR=/path/to/q2-longitudinal Q2SC_DIR=/path/to/q2-sample-classifier \
#       python -I data-raw/expected/feature_volatility.py   # q2's outputs
#     Rscript data-raw/expected/feature_volatility.r     # folds them in
#
# q2-longitudinal 2026.10.0 and q2-sample-classifier 2026.10.0 (both BSD-3)
# are the reference, with scikit-learn 1.9.1; see the Python script for how
# their functions are run without a QIIME 2 install.
#
# Inputs (simulated, written to data-raw/expected/feature_volatility/):
#
#   meta.csv     One row per sample (in a shuffled order): sample, subject
#                (S01-S16; S01-S08 in group A), state (days 0, 3, 7, 14,
#                28) and group. S04 has no sample at day 7 and S12 none at
#                day 28, so 78 samples.
#   counts.csv   Samples x 25 features (F01-F25), integer counts, in
#                meta.csv's order. F01-F06 rise or fall over time; every
#                feature is present in at least 20 samples.
#
# Fields of the .rds:
#
#   meta, counts   The inputs.
#   cases          One list per q2 run (named "rs<random_state>" without
#                  RFE, "rfe-rs<random_state>" with q2's RFECV), each:
#                  random_state; train and test (sample ids, in q2's order:
#                  the shuffle's); importance (named numeric, q2's
#                  feature_importances_ in its sorted order, from the forest
#                  of 100 trees on the training samples' counts: with RFE,
#                  on the features RFECV kept only); predicted
#                  (the test samples' predictions); accuracy (q2's
#                  _linear_regress(): mse, r, r2, slope, intercept); and
#                  stats_all / stats_10 (q2's _importance_filtering() with
#                  feature_count 'all' and 10, then the importances joined
#                  with _summarize_feature_stats() on relative frequencies:
#                  a data frame, feature first, q2's columns renamed as
#                  ecodive's). With RFE also rfe: ranking (RFECV's
#                  ranking_, named by feature) and scores (its cv_results_:
#                  n_features, mean, sd, fold_1 ... fold_5).
#   perm           numpy's RandomState(rs).permutation(n): a list of
#                  list(rs, n, perm) (0-based).
#   versions       Of the Python packages used.

dir <- file.path('data-raw', 'expected', 'feature_volatility')
dir.create(dir, showWarnings = FALSE)

hexf  <- function (x) ifelse(is.na(x), 'NA', sprintf('%a', x))
unhex <- function (x) { x[x %in% c('NA', 'nan', '')] <- NA; as.numeric(x) }

if (!file.exists(file.path(dir, 'meta.csv'))) {

  set.seed(20261126)
  states <- c(0, 3, 7, 14, 28)
  rows   <- list()
  for (s in 1:16) {
    st <- states
    if (s == 4)  st <- setdiff(st, 7)
    if (s == 12) st <- setdiff(st, 28)
    rows[[s]] <- data.frame(subject = sprintf('S%02d', s), state = st)
  }
  meta <- do.call(rbind, rows)
  meta$group <- ifelse(meta$subject %in% sprintf('S%02d', 1:8), 'A', 'B')

  n <- nrow(meta); p <- 25
  base  <- exp(rnorm(p, 2.5, 0.7))
  subj  <- exp(matrix(rnorm(16 * p, 0, 0.4), 16, p))
  sid   <- as.integer(substr(meta$subject, 2, 3))
  slope <- c(1.5, -1.5, 1.0, -1.0, 0.6, -0.6, rep(0, p - 6))
  mu    <- subj[sid, ] * rep(base, each = n) * exp(outer(meta$state / 28, slope))
  counts <- matrix(rnbinom(n * p, mu = mu, size = 3), n, p)
  stopifnot(all(colSums(counts > 0) >= 20), all(rowSums(counts) > 0))

  ord  <- sample(n)
  meta <- meta[ord, ]
  counts <- counts[ord, ]
  meta$sample <- sprintf('smp%02d', seq_len(n))
  dimnames(counts) <- list(meta$sample, sprintf('F%02d', 1:p))
  meta <- meta[, c('sample', 'subject', 'state', 'group')]
  rownames(meta) <- NULL

  write.csv(meta, file.path(dir, 'meta.csv'), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(sample = meta$sample, counts, check.names = FALSE),
            file.path(dir, 'counts.csv'), row.names = FALSE, quote = FALSE)

  message('Inputs written; now run data-raw/expected/feature_volatility.py.')
  quit(save = 'no')
}


#----------------------------------------------------------
# Fold q2's outputs in.
#----------------------------------------------------------
rd <- function (f) read.csv(file.path(dir, f), colClasses = 'character', check.names = FALSE)

meta <- rd('meta.csv')
meta$state <- as.numeric(meta$state)
rownames(meta) <- meta$sample

cc     <- rd('counts.csv')
counts <- as.matrix(cc[, -1]); storage.mode(counts) <- 'integer'
rownames(counts) <- cc$sample

# "key<TAB>values" files.
read_kv <- function (f) {
  ln <- readLines(file.path(dir, f))
  kv <- strsplit(ln, '\t', fixed = TRUE)
  structure(
    lapply(kv, function (x) if (length(x) < 2 || x[[2]] == '') character(0) else strsplit(x[[2]], ' ', fixed = TRUE)[[1]]),
    names = vapply(kv, `[[`, '', 1) )
}

STATS <- c(
  'Cumulative Avg Decrease'   = 'cumulative_avg_decrease',
  'Cumulative Avg Increase'   = 'cumulative_avg_increase',
  'Net Avg Change'            = 'net_avg_change',
  'Global Variance'           = 'global_variance',
  'Global Mean'               = 'global_mean',
  'Global Median'             = 'global_median',
  'Global Standard Deviation' = 'global_sd',
  'Global CV (%)'             = 'global_cv' )

read_stats <- function (kv, prefix) {
  df <- data.frame(feature = kv[[paste0(prefix, '_feature')]],
                   importance = unhex(kv[[paste0(prefix, '_importance')]]),
                   stringsAsFactors = FALSE)
  for (s in names(STATS)) df[[STATS[[s]]]] <- unhex(kv[[paste0(prefix, '_', s)]])
  df
}

cases <- list()
for (f in sort(list.files(dir, '^case-(rfe-)?rs[0-9]+[.]txt$'))) {
  kv <- read_kv(f)
  nm <- sub('^case-((rfe-)?rs[0-9]+)[.]txt$', '\\1', f)
  cases[[nm]] <- list(
    random_state = as.integer(kv$random_state),
    train        = kv$train,
    test         = kv$test,
    importance   = structure(unhex(kv$importance), names = kv$importance_features),
    predicted    = unhex(kv$predicted),
    accuracy     = structure(unhex(kv$accuracy), names = c('mse', 'r', 'r2', 'slope', 'intercept')),
    stats_all    = read_stats(kv, 'all'),
    stats_10     = read_stats(kv, 'top10') )
  if (!is.null(kv$rfe_ranking)) {
    sc <- data.frame(n_features = as.integer(kv$rfe_n_features),
                     mean = unhex(kv$rfe_mean), sd = unhex(kv$rfe_sd))
    for (k in 1:5) sc[[paste0('fold_', k)]] <- unhex(kv[[paste0('rfe_fold_', k)]])
    cases[[nm]]$rfe <- list(
      ranking = structure(as.integer(kv$rfe_ranking), names = kv$rfe_features),
      scores  = sc )
  }
}

kv   <- read_kv('perm.txt')
perm <- lapply(grep('^perm_', names(kv), value = TRUE), function (k) {
  a <- as.integer(strsplit(sub('^perm_', '', k), '_', fixed = TRUE)[[1]])
  list(rs = a[[1]], n = a[[2]], perm = as.integer(kv[[k]]))
})

kv <- read_kv('versions.txt')
versions <- vapply(kv, paste, '', collapse = ' ')

saveRDS(
  list(meta = meta, counts = counts, cases = cases, perm = perm,
       versions = versions),
  file.path('tests', 'tinytests', 'expected', 'feature_volatility.rds'),
  version = 2 )
message('Wrote tests/tinytests/expected/feature_volatility.rds')
