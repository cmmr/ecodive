# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Reference values for lefse(), from LEfSe 1.1.2 (data-raw/expected/lefse.py),
# R's kruskal.test() and coin's wilcox_test() (the functions LEfSe calls),
# and Bioconductor lefser 1.12.1 (the newest that installs on R 4.3.3).
#
# Run three times:
#
#   Rscript data-raw/expected/lefse.r                    # writes the inputs
#   PYTHONHASHSEED=0 python data-raw/expected/lefse.py   # LEfSe's outputs
#   Rscript data-raw/expected/lefse.r                    # folds them into
#                                                        # tests/tinytests/expected/lefse.rds
#
# The second R run needs coin (1.4-3 or later) and lefser with its
# dependencies (SummarizedExperiment); ecodive itself is not used.

dir <- file.path('data-raw', 'expected', 'lefse')
dir.create(dir, showWarnings = FALSE, recursive = TRUE)


#------------------------------------------------------------------------------
# Inputs: samples sorted by class, then subclass, as LEfSe sorts them.
#------------------------------------------------------------------------------
make_input <- function (name, cls, sub, D, mu, eff) {
  n <- length(cls)
  x <- matrix(0L, n, D, dimnames = list(sprintf('S%02d', seq_len(n)), sprintf('F%02d', seq_len(D))))
  for (j in seq_len(D)) x[, j] <- rpois(n, mu[[j]] * eff[, j])
  utils::write.csv(x, file.path(dir, paste0(name, '_counts.csv')), row.names = FALSE, quote = FALSE)
  utils::write.csv(
    data.frame(class = cls, subclass = if (is.null(sub)) '' else sub),
    file.path(dir, paste0(name, '_meta.csv')), row.names = FALSE, quote = FALSE )
}

if (!file.exists(file.path(dir, 'cases.csv'))) {

  set.seed(20261124)

  # Two classes, 20 + 20, no subclasses. F01-F06 richer in b, F07-F09 in a;
  # F10 has few distinct values (jittered); F11 is all zeros; the rest are
  # null, some sparse.
  cls <- rep(c('a', 'b'), each = 20)
  D   <- 30
  mu  <- c(rep(40, 9), 0.8, 0, exp(seq(log(0.5), log(400), length.out = D - 11)))
  eff <- matrix(1, 40, D)
  eff[cls == 'b', 1:6] <- 3
  eff[cls == 'a', 7:9] <- 3
  eff[cls == 'b', 10]  <- 2.5
  make_input('two', cls, NULL, D, mu, eff)

  # Two classes with two subclasses each; a's s2 has 6 samples, under
  # min_c = 10 (compared by medians only). F01-F04 richer in b in every
  # subclass, F05-F06 only in b's s1, F07-F08 richer in a.
  cls <- rep(c('a', 'b'), c(18, 23))
  sub <- c(rep(c('s1', 's2'), c(12, 6)), rep(c('s1', 's2'), c(11, 12)))
  D   <- 25
  mu  <- c(rep(60, 8), exp(seq(log(1), log(300), length.out = D - 8)))
  eff <- matrix(1, 41, D)
  eff[cls == 'b', 1:4] <- 3
  eff[cls == 'b' & sub == 's1', 5:6] <- 4
  eff[cls == 'a', 7:8] <- 2.5
  make_input('sub', cls, sub, D, mu, eff)

  # The same subclass names in both classes (for -e 1): x and y, 12 each.
  cls <- rep(c('a', 'b'), each = 24)
  sub <- rep(rep(c('x', 'y'), each = 12), 2)
  D   <- 20
  mu  <- c(rep(50, 6), exp(seq(log(1), log(200), length.out = D - 6)))
  eff <- matrix(1, 48, D)
  eff[cls == 'b', 1:3] <- 3
  eff[cls == 'b' & sub == 'x', 4:5] <- 3     # only in the x subclasses
  eff[cls == 'a' & sub == 'y', 6]   <- 3
  make_input('same', cls, sub, D, mu, eff)

  # Three classes of 15. F01-F03 richer in a, F04-F06 in b, F07-F08 poorer
  # in c, F09 richer in b and c.
  cls <- rep(c('a', 'b', 'c'), each = 15)
  D   <- 25
  mu  <- c(rep(50, 9), exp(seq(log(1), log(300), length.out = D - 9)))
  eff <- matrix(1, 45, D)
  eff[cls == 'a', 1:3] <- 3
  eff[cls == 'b', 4:6] <- 3
  eff[cls == 'c', 7:8] <- 0.3
  eff[cls != 'a', 9]   <- 3
  make_input('three', cls, NULL, D, mu, eff)

  # Wide: 10 + 10 samples, 40 features of which most differ, so that the
  # features scored outnumber the 13 samples of a subsample.
  cls <- rep(c('a', 'b'), each = 10)
  D   <- 40
  mu  <- exp(seq(log(20), log(500), length.out = D))
  eff <- matrix(1, 20, D)
  eff[cls == 'b', 1:20]  <- 5
  eff[cls == 'a', 21:34] <- 5
  make_input('wide', cls, NULL, D, mu, eff)

  # Deep samples: F01 dominates; F02-F05 are rare (a few per million) and
  # richer in b, so that their scores fall below 2; F06-F09 are richer in b
  # at a few hundred per million.
  cls <- rep(c('a', 'b'), each = 20)
  D   <- 20
  mu  <- c(1e5, rep(0.6, 4), rep(30, 4), exp(seq(log(1), log(3000), length.out = D - 9)))
  eff <- matrix(1, 40, D)
  eff[cls == 'b', 2:5] <- 4
  eff[cls == 'b', 6:9] <- 3
  make_input('low', cls, NULL, D, mu, eff)

  cases <- data.frame(
    name       = c('two', 'sub', 'same', 'three_all', 'three_one', 'three_nowilc', 'wide', 'sub_boots', 'low'),
    input      = c('two', 'sub', 'same', 'three',     'three',     'three',        'wide', 'sub',       'low'),
    wilc       = c(1,     1,     1,      1,           1,           0,              1,      1,           1),
    multiclass = c(0,     0,     0,      0,           1,           0,              0,      0,           0),
    same_subcl = c(0,     0,     1,      0,           0,           0,              0,      0,           0),
    min_c      = c(10,    10,    10,     10,          10,          10,             10,     5,           10),
    n_boots    = c(30,    30,    30,     30,          30,          30,             30,     1,           30),
    f          = c(0.67,  0.67,  0.67,   0.67,        0.67,        0.67,           0.67,   0.5,         0.67) )
  utils::write.csv(cases, file.path(dir, 'cases.csv'), row.names = FALSE, quote = FALSE)

  message('Wrote the inputs; now run data-raw/expected/lefse.py.')
  quit(save = 'no')
}


#------------------------------------------------------------------------------
# Fold LEfSe's outputs, and R's, into the fixture.
#------------------------------------------------------------------------------
cases <- utils::read.csv(file.path(dir, 'cases.csv'), stringsAsFactors = FALSE)
rd    <- function (f, ...) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE, ...)

inputs <- list()
for (inp in unique(cases$input)) {
  meta <- rd(paste0(inp, '_meta.csv'), colClasses = 'character')
  inputs[[inp]] <- list(
    counts   = as.matrix(rd(paste0(inp, '_counts.csv'))),
    class    = meta$class,
    subclass = if (all(is.na(meta$subclass) | meta$subclass == '')) NULL else meta$subclass )
}

out <- list(inputs = inputs, cases = list())
for (i in seq_len(nrow(cases))) {
  cs   <- cases[i, ]
  nm   <- cs$name
  kw   <- rd(paste0(nm, '_kw.csv'))
  wil  <- rd(paste0(nm, '_wilcoxon.csv'))
  lda  <- rd(paste0(nm, '_lda.csv'))
  jit  <- if (nrow(lda)) as.matrix(rd(paste0(nm, '_jitter.csv'))) else NULL
  bts  <- if (nrow(lda)) as.matrix(utils::read.csv(file.path(dir, paste0(nm, '_boots.csv')), header = FALSE)) else NULL
  mns  <- rd(paste0(nm, '_means.csv'))
  norm <- as.matrix(rd(paste0(nm, '_norm.csv')))
  if (!is.null(jit)) dimnames(jit) <- NULL
  if (!is.null(bts)) dimnames(bts) <- NULL
  out$cases[[nm]] <- list(
    settings = as.list(cs),
    kw_p     = stats::setNames(kw$p, kw$feature),
    wilcoxon = stats::setNames(as.logical(wil$pass), wil$feature),
    lda      = stats::setNames(lda$score, lda$feature),
    jitter   = jit,
    boots    = if (!is.null(bts)) t(bts) + 1L,     # m x n_boots, 1-based
    means    = { m <- t(as.matrix(mns[, -1])); colnames(m) <- mns$feature; m },
    norm     = norm )
}

# R's own tests on the scaled tables: kruskal.test() for every feature, and
# coin's asymptotic wilcox_test() for every pair of subclasses of the 'sub'
# input (LEfSe tests a pair only with at least min_c samples a side; these
# are all of them).
suppressMessages(library(coin))
for (nm in names(out$cases)) {
  cs  <- out$cases[[nm]]
  inp <- out$inputs[[cs$settings$input]]
  g   <- factor(inp$class)
  out$cases[[nm]]$kruskal <- apply(cs$norm, 2, function (y) stats::kruskal.test(y ~ g)$p.value)
}
cs  <- out$cases$sub
inp <- out$inputs$sub
key <- paste(inp$class, inp$subclass)
grid <- expand.grid(a = c('a s1', 'a s2'), b = c('b s1', 'b s2'), stringsAsFactors = FALSE)
out$coin <- sapply(seq_len(nrow(grid)), function (r) apply(cs$norm, 2, function (y) {
  i <- key %in% c(grid$a[[r]], grid$b[[r]])
  d <- data.frame(x = y[i], y = factor(key[i] == grid$b[[r]]))
  if (length(unique(d$x)) < 2) return (NA_real_)
  coin::pvalue(coin::wilcox_test(x ~ y, data = d)) }))
out$coin_pairs <- grid

# lefser 1.12.1 on the 'two' input as counts per million, with its seed.
suppressMessages(library(lefser))
inp <- out$inputs$two
cpm <- t(inp$counts / rowSums(inp$counts) * 1e6)
se  <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(counts = cpm),
  colData = S4Vectors::DataFrame(GROUP = inp$class) )
set.seed(1)
lr <- suppressMessages(suppressWarnings(lefser::lefser(se, groupCol = 'GROUP', lda.threshold = 0)))
out$lefser <- list(
  version = as.character(utils::packageVersion('lefser')),
  scores  = stats::setNames(lr$scores, lr$Names) )

saveRDS(out, file.path('tests', 'tinytests', 'expected', 'lefse.rds'), version = 2)
message('Wrote tests/tinytests/expected/lefse.rds')
