# Writes tests/tinytests/expected/random_forest.rds: inputs and reference
# results for test-random_forest.r, so that the tests need no reference
# package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/random_forest.r     # writes the inputs
#     Q2SC_DIR=/path/to/q2-sample-classifier \
#       python -I data-raw/expected/random_forest.py  # scikit-learn and q2
#     R_LIBS=<lib with randomForest and ranger> \
#       Rscript data-raw/expected/random_forest.r   # folds them in
#
# References:
#   * scikit-learn 1.9.1 (BSD-3) with numpy >= 2: trees, forests,
#     importances, OOB predictions, folds, roc_curve() / auc(), and
#     numpy's legacy RandomState draws (see the Python script).
#   * q2-sample-classifier 2026.10.0 (BSD-3): _fit_and_predict_cv() and its
#     ROC averages, run from a checkout without a QIIME 2 install.
#   * randomForest and ranger (GPL; outputs only, distributional checks):
#     OOB error over 20 seeds and permutation importances.
#
# Inputs (simulated with R's RNG, written to data-raw/expected/random_forest/;
# real numbers as hexadecimal floats):
#
#   bin.csv     60 samples x 40 features, label A or B (30 each, shuffled).
#               Counts with many zeros and ties; F07 is 3 in every sample
#               (a constant feature); F01-F03 are larger in B.
#   three.csv   45 samples x 30 features, label A, B or C (20, 15, 10).
#               Counts, with fewer zeros (every feature is non-zero in at
#               least 10 samples, so that q2's DictVectorizer keeps every
#               feature in every training fold).
#   reg.csv     50 samples x 30 features, y (an age-like number). Relative
#               abundances (proportions, hexadecimal), every feature
#               non-zero in at least 10 samples.
#
# Fields of the .rds:
#
#   bin, three, reg   The inputs: list(x (samples x features matrix, with
#                     dimnames), y (factor or numeric)).
#   cases             Named list, one per scikit-learn case: data ('bin',
#                     'three' or 'reg'), sparse (TRUE: CSC input), params
#                     (the estimator's arguments, ecodive's names),
#                     random_state, and the fitted forest: trees (offsets
#                     and the concatenated tree_ arrays: left, right,
#                     feature, threshold, value (k x nodes), impurity,
#                     n_node_samples, weighted_n_node_samples), inbag (n x
#                     trees bootstrap counts from estimators_samples_),
#                     importance (feature_importances_), pred
#                     (predict_proba() or predict() on the training table),
#                     oob (oob_decision_function_ or oob_prediction_),
#                     oob_score.
#   rng               numpy draws: random_sample (seed 0, 1000), randint
#                     (seed 1; n = 7, 60, 1000, 2^31 - 1; 1000 each),
#                     tree_seeds (randint(2^31 - 1) x 1000 from seed 7),
#                     our_rand_r (1000 from 12345), np_sum (arrays and
#                     their np.sum()).
#   folds             StratifiedKFold(5, shuffle, random_state = 3) on
#                     three's y and KFold(5, shuffle, random_state = 3) on
#                     reg: the test fold (1-based) of each sample.
#   roc               roc_curve() / auc() on stored scores (binary, with
#                     ties) and q2's per-class, micro and macro curves.
#   q2_cls, q2_reg    q2-sample-classifier's _fit_and_predict_cv() on three
#                     (RandomForestClassifier, 20 trees, random_state 5,
#                     cv 5) and reg (RandomForestRegressor): folds (test
#                     fold of each sample), predictions, probabilities,
#                     importances (mean over folds), scores (per fold),
#                     and q2's ROC (classification).
#   dist              Distributional references: scikit-learn's OOB
#                     accuracy and impurity importances over 20
#                     random_states on bin (500 trees); randomForest's and
#                     ranger's OOB error over 20 seeds and mean permutation
#                     importances on bin.
#   versions          Of the packages used.

dir <- file.path('data-raw', 'expected', 'random_forest')
dir.create(dir, showWarnings = FALSE)

hexf  <- function (x) ifelse(is.na(x), 'NA', sprintf('%a', x))
unhex <- function (x) {
  x[x %in% c('NA', 'nan', '')] <- NA
  x[x == 'inf'] <- 'Inf'
  x[x == '-inf'] <- '-Inf'
  as.numeric(x)
}

write_table <- function (x, y, file) {
  df <- data.frame(
    sample = rownames(x),
    y      = if (is.numeric(y)) hexf(y) else as.character(y),
    apply(x, 2, hexf), check.names = FALSE )
  utils::write.csv(df, file.path(dir, file), row.names = FALSE, quote = FALSE)
}

read_table <- function (file, numeric_y = FALSE) {
  df <- utils::read.csv(file.path(dir, file), colClasses = 'character', check.names = FALSE)
  x  <- sapply(df[-(1:2)], unhex)
  rownames(x) <- df$sample
  y  <- if (numeric_y) unhex(df$y) else factor(df$y)
  list(x = x, y = y)
}

feat_names <- function (p) sprintf('F%02d', seq_len(p))


if (!file.exists(file.path(dir, 'bin.csv'))) {

  set.seed(2026)

  # bin: two classes.
  n <- 60; p <- 40
  y <- sample(rep(c('A', 'B'), each = 30))
  mu <- exp(stats::rnorm(p, 1, 1))
  x  <- sapply(mu, function (m) stats::rpois(n, m))
  x[stats::runif(n * p) < 0.5] <- 0
  b  <- y == 'B'
  x[b, 1:3] <- x[b, 1:3] + stats::rpois(sum(b) * 3, 4)
  x[, 7] <- 3
  dimnames(x) <- list(sprintf('S%02d', seq_len(n)), feat_names(p))
  write_table(x, y, 'bin.csv')

  # three: three classes, every feature in >= 10 samples.
  n <- 45; p <- 30
  y <- sample(rep(c('A', 'B', 'C'), c(20, 15, 10)))
  mu <- exp(stats::rnorm(p, 1.5, 0.8))
  repeat {
    x <- sapply(mu, function (m) stats::rpois(n, m))
    x[stats::runif(n * p) < 0.3] <- 0
    if (all(colSums(x > 0) >= 10)) break
  }
  x[y == 'B', 4:5] <- x[y == 'B', 4:5] + 5
  x[y == 'C', 6]   <- x[y == 'C', 6] + 8
  dimnames(x) <- list(sprintf('S%02d', seq_len(n)), feat_names(p))
  write_table(x, y, 'three.csv')

  # reg: proportions, an age-like y.
  n <- 50; p <- 30
  age <- round(stats::runif(n, 0, 36), 1)
  repeat {
    x <- sapply(exp(stats::rnorm(p, 2, 1)), function (m) stats::rpois(n, m))
    x[stats::runif(n * p) < 0.3] <- 0
    if (all(colSums(x > 0) >= 10) && all(rowSums(x) > 0)) break
  }
  x[, 1] <- x[, 1] + stats::rpois(n, age)
  x[, 2] <- x[, 2] + stats::rpois(n, 36 - age)
  x <- x / rowSums(x)
  dimnames(x) <- list(sprintf('S%02d', seq_len(n)), feat_names(p))
  write_table(x, age, 'reg.csv')

  cat('Inputs written to', dir, '\n')
  quit(save = 'no')
}


#----------------------------------------------------------------------
# Fold the outputs in.
#----------------------------------------------------------------------

bin   <- read_table('bin.csv')
three <- read_table('three.csv')
reg   <- read_table('reg.csv', numeric_y = TRUE)

# Text outputs of the Python script: one file per case, "name<TAB>values"
# lines, values separated by spaces.
read_kv <- function (file) {
  lines <- readLines(file.path(dir, file))
  kv    <- strsplit(lines, '\t', fixed = TRUE)
  keys  <- vapply(kv, `[[`, '', 1)
  vals  <- lapply(kv, function (z) if (length(z) < 2 || z[[2]] == '') character(0) else strsplit(z[[2]], ' ', fixed = TRUE)[[1]])
  stats::setNames(vals, keys)
}

num <- function (kv, key) unhex(kv[[key]])
int <- function (kv, key) as.integer(kv[[key]])

case_names <- readLines(file.path(dir, 'cases.txt'))
cases <- list()
for (nm in case_names) {
  kv <- read_kv(paste0('case-', nm, '.txt'))
  data <- kv$data
  n  <- nrow(get(data)$x)
  k  <- int(kv, 'k')
  nt <- int(kv, 'n_trees')
  params <- list()
  for (key in grep('^param_', names(kv), value = TRUE)) {
    v <- kv[[key]]
    nmp <- sub('^param_', '', key)
    params[[nmp]] <- if (nmp == 'class_weight' && length(v) > 1)
      stats::setNames(unhex(v[c(FALSE, TRUE)]), v[c(TRUE, FALSE)])
    else if (grepl('^[-0-9.]+$', v[[1]]) && !grepl('[.]', v[[1]])) as.integer(v)
    else if (grepl('^[-0-9.]+$', v[[1]])) as.numeric(v)
    else v
  }
  cases[[nm]] <- list(
    data         = data,
    sparse       = kv$sparse == 'True',
    params       = params,
    random_state = int(kv, 'random_state'),
    trees = list(
      offsets                 = int(kv, 'offsets'),
      left                    = int(kv, 'left'),
      right                   = int(kv, 'right'),
      feature                 = int(kv, 'feature'),
      threshold               = num(kv, 'threshold'),
      value                   = matrix(num(kv, 'value'), nrow = k),
      impurity                = num(kv, 'impurity'),
      n_node_samples          = int(kv, 'n_node_samples'),
      weighted_n_node_samples = num(kv, 'weighted_n_node_samples') ),
    inbag      = matrix(int(kv, 'inbag'), nrow = n, ncol = nt),
    importance = num(kv, 'importance'),
    pred       = matrix(num(kv, 'pred'), nrow = k),
    oob        = matrix(num(kv, 'oob'), nrow = k),
    oob_score  = num(kv, 'oob_score') )
}

kv  <- read_kv('rng.txt')
rng <- list(
  random_sample = num(kv, 'random_sample'),
  randint       = list(
    `7`          = num(kv, 'randint_7'),
    `60`         = num(kv, 'randint_60'),
    `1000`       = num(kv, 'randint_1000'),
    `2147483647` = num(kv, 'randint_2147483647') ),
  tree_seeds    = num(kv, 'tree_seeds'),
  our_rand_r    = num(kv, 'our_rand_r'),
  np_sum        = lapply(grep('^np_sum_in_', names(kv), value = TRUE), function (key)
    list(x = num(kv, key), sum = num(kv, sub('_in_', '_out_', key)))) )

kv <- read_kv('folds.txt')
folds <- list(three = int(kv, 'three'), reg = int(kv, 'reg'))

kv  <- read_kv('roc.txt')
roc_keys <- unique(sub('_(fpr|tpr|thresholds|auc|truth|score)$', '', names(kv)))
roc <- lapply(stats::setNames(roc_keys, roc_keys), function (key) {
  out <- list()
  for (f in c('truth', 'score', 'fpr', 'tpr', 'thresholds', 'auc'))
    if (!is.null(kv[[paste0(key, '_', f)]])) out[[f]] <- num(kv, paste0(key, '_', f))
  out
})

read_q2 <- function (file) {
  kv <- read_kv(file)
  out <- list(
    folds       = int(kv, 'folds'),
    sample      = kv$sample,
    predictions = if (!is.null(kv$classes)) kv$predictions else num(kv, 'predictions'),
    importances = stats::setNames(num(kv, 'importances'), kv$importance_features),
    scores      = num(kv, 'scores') )
  if (!is.null(kv$classes)) {
    out$classes       <- kv$classes
    out$probabilities <- matrix(num(kv, 'probabilities'), ncol = length(kv$classes),
                                dimnames = list(kv$sample, kv$classes))
    for (cl in c(kv$classes, 'micro', 'macro'))
      out$roc[[cl]] <- list(
        fpr = num(kv, paste0('roc_', cl, '_fpr')),
        tpr = num(kv, paste0('roc_', cl, '_tpr')),
        auc = num(kv, paste0('roc_', cl, '_auc')) )
  }
  out
}
q2_cls <- read_q2('q2-cls.txt')
q2_reg <- read_q2('q2-reg.txt')

kv <- read_kv('dist.txt')
dist <- list(
  sklearn_oob_accuracy = num(kv, 'oob_accuracy'),
  sklearn_importance   = matrix(num(kv, 'importance'), nrow = ncol(bin$x)) )

versions <- read_kv('versions.txt')


# randomForest and ranger (GPL; run for their outputs only).
if (requireNamespace('randomForest', quietly = TRUE) &&
    requireNamespace('ranger', quietly = TRUE)) {

  x <- bin$x; y <- bin$y
  rf_err  <- numeric(20)
  rg_err  <- numeric(20)
  rf_perm <- matrix(0, ncol(x), 20)
  rg_perm <- matrix(0, ncol(x), 20)
  for (s in 1:20) {
    set.seed(s)
    fit <- randomForest::randomForest(x, y, ntree = 500, importance = TRUE)
    rf_err[[s]]   <- fit$err.rate[500, 'OOB']
    rf_perm[, s]  <- randomForest::importance(fit, type = 1, scale = FALSE)[, 1]
    fit <- ranger::ranger(x = x, y = y, num.trees = 500, importance = 'permutation', seed = s)
    rg_err[[s]]   <- fit$prediction.error
    rg_perm[, s]  <- fit$variable.importance
  }
  dist$randomForest_oob_error <- rf_err
  dist$ranger_oob_error       <- rg_err
  dist$randomForest_perm      <- rowMeans(rf_perm)
  dist$ranger_perm            <- rowMeans(rg_perm)
  versions$randomForest <- as.character(utils::packageVersion('randomForest'))
  versions$ranger       <- as.character(utils::packageVersion('ranger'))
} else {
  stop('randomForest and ranger are needed for the distributional references.')
}

saveRDS(
  list(
    bin = bin, three = three, reg = reg, cases = cases, rng = rng,
    folds = folds, roc = roc, q2_cls = q2_cls, q2_reg = q2_reg,
    dist = dist, versions = versions ),
  file.path('tests', 'tinytests', 'expected', 'random_forest.rds'),
  version = 2, compress = 'xz' )

cat('Wrote tests/tinytests/expected/random_forest.rds\n')
