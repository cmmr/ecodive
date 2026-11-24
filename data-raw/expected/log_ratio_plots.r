# Fixtures for test-log_ratio_plots.r: Qurro 0.9.0's autoselection and log
# ratios for rankings made by ecodive, so that the tests need no reference
# package. Three steps, from the package root, with ecodive installed:
#
#     Rscript data-raw/expected/log_ratio_plots.r inputs   # 1. the rankings
#     python  data-raw/expected/log_ratio_plots.py         # 2. Qurro's outputs
#     Rscript data-raw/expected/log_ratio_plots.r          # 3. the .rds
#
# Step 1 writes the count tables and, for each case, a ranking (one value
# per feature) and n to data-raw/expected/log_ratio_plots/: synthetic
# rankings with and without ties, and the effects and loadings of
# dirmult_ttest(), ancombc2(), da_ttest(), da_wilcoxon(), da_lm(), linda(),
# pca() and rpca() on one table. Non-finite effects are left out, as
# log_ratio_features() leaves them out, since Qurro can't take them.
#
# Step 2 (log_ratio_plots.py, with log_ratio_plots.js under Node.js) runs
# Qurro's own code: process_input() (which features and in which order the
# rank plot holds), the autoselection of its JavaScript filterFeatures()
# ("autoLiteralTop" / "autoLiteralBot"), qarcoal()'s log ratio, and the
# sample plot's computeBalance(). Every number is a hexadecimal float, read
# back exactly.
#
# Step 3 stores the tables, the cases and Qurro's results in
# tests/tinytests/expected/log_ratio_plots.rds.

dir  <- file.path('data-raw', 'expected', 'log_ratio_plots')
args <- commandArgs(trailingOnly = TRUE)

hexs <- function (v) ifelse(is.na(v), 'NA', sprintf('%a', v))
hex  <- function (s) suppressWarnings(as.numeric(s))  # 'NA' -> NA

write_table <- function (m, nm)
  write.csv(
    data.frame(sample = rownames(m), m, check.names = FALSE),
    file.path(dir, sprintf('table-%s.csv', nm)), row.names = FALSE )

read_table <- function (nm) {
  d <- read.csv(file.path(dir, sprintf('table-%s.csv', nm)), check.names = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE])
  rownames(m) <- d$sample
  storage.mode(m) <- 'integer'
  m
}


if (identical(args, 'inputs')) {

  library(ecodive)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Table 'syn': 20 samples x 30 features, counts with zeros (no feature
  # or sample all zero). Table 'da': 18 x 40, with a difference in group a.
  set.seed(20261124)
  syn <- matrix(
    rpois(20 * 30, rep(rlnorm(30, 1.5, 1.2), each = 20)), 20, 30,
    dimnames = list(paste0('S', 1:20), paste0('F', 1:30)) )
  syn[, colSums(syn) == 0] <- 1L

  set.seed(1)
  da <- matrix(
    rpois(18 * 40, rep(rlnorm(40, 2, 1.2), each = 18)), 18, 40,
    dimnames = list(paste0('s', 1:18), paste0('f', 1:40)) )
  da[1:6, 1:5] <- da[1:6, 1:5] * 6L
  da[, colSums(da) == 0] <- 1L
  meta <- data.frame(
    grp = rep(c('a', 'b', 'c'), each = 6), age = round(runif(18, 20, 60)),
    row.names = rownames(da) )

  write_table(syn, 'syn')
  write_table(da,  'da')

  cases <- list()
  add <- function (id, table, effect, n)
    for (k in n) cases[[length(cases) + 1L]] <<- list(
      id = sprintf('%s_n%d', id, k), table = table, n = k,
      feature = names(effect)[is.finite(effect)],
      effect = unname(effect[is.finite(effect)]) )

  # Synthetic rankings: no ties, ties at one decimal, and seven levels.
  set.seed(7)
  r <- rnorm(30)
  add('syn_norm',  'syn', setNames(r, colnames(syn)), c(1, 3, 5, 15))
  add('syn_round', 'syn', setNames(round(r, 1), colnames(syn)), c(2, 5, 9))
  add('syn_level', 'syn', setNames(as.double(sample(-3:3, 30, TRUE)), colnames(syn)), c(1, 4, 7, 12))
  add('syn_zero',  'syn', setNames(c(rep(0, 10), r[11:30]), colnames(syn)), c(3, 10))

  # Rankings that ecodive's results define.
  eff <- function (x, ...) { s <- ecodive:::rank_source(x, ...); setNames(s$effect, s$feature) }
  fit <- dirmult_ttest(da, meta$grp, contrast = c('a', 'b'), n_draws = 64L, seed = 1L, cpus = 1L)
  add('dirmult_ttest', 'da', eff(fit), c(1, 5, 20))
  fit <- ancombc2(da, ~ grp + age, meta, groups = 'grp')
  add('ancombc2_grpb', 'da', eff(fit), 5)
  add('ancombc2_age',  'da', eff(fit, 'age'), 5)
  add('da_ttest',      'da', eff(da_ttest(da, meta$grp, pseudocount = 0.5)), 5)
  add('da_wilcoxon',   'da', eff(da_wilcoxon(da, meta$grp, pseudocount = 0.5)), 5)
  add('da_lm',         'da', eff(da_lm(da, ~ grp + age, meta, groups = 'grp', pseudocount = 0.5)), 5)
  add('linda',         'da', eff(linda(da, ~ grp + age, meta, groups = 'grp')), 5)
  ord <- pca(da, norm = 'clr', pseudocount = 0.5)
  add('pca_ax1', 'da', eff(ord, axes = 1L), 5)
  add('pca_ax2', 'da', eff(ord, axes = 2L), c(3, 10))
  add('rpca',    'da', eff(rpca(da), axes = 1L), 5)

  # One line per case and feature: the ranking, in its order.
  rows <- do.call(rbind, lapply(cases, function (cs) data.frame(
    id = cs$id, table = cs$table, n = cs$n, feature = cs$feature,
    effect = hexs(cs$effect) )))
  write.csv(rows, file.path(dir, 'rankings.csv'), row.names = FALSE)
  saveRDS(meta, file.path(dir, 'meta.rds'), version = 2)

} else {

  rankings <- read.csv(file.path(dir, 'rankings.csv'), colClasses = 'character')
  qurro    <- read.csv(file.path(dir, 'qurro.csv'),    colClasses = 'character')
  values   <- read.csv(file.path(dir, 'values.csv'),   colClasses = 'character')

  ids   <- unique(rankings$id)
  cases <- lapply(ids, function (id) {
    r <- rankings[rankings$id == id, ]
    q <- qurro[qurro$id == id, ]
    v <- values[values$id == id, ]
    split_ids <- function (s) if (nzchar(s)) strsplit(s, ' ')[[1]] else character(0)
    list(
      id      = id,
      table   = r$table[[1]],
      n       = as.integer(r$n[[1]]),
      ranking = setNames(hex(r$effect), r$feature),
      # Qurro's rank plot features, in its dataset's order.
      qurro_features = split_ids(q$features),
      top     = split_ids(q$top),
      bot     = split_ids(q$bot),
      qarcoal = setNames(hex(v$qarcoal), v$sample),
      balance = setNames(hex(v$balance), v$sample) )
  })
  names(cases) <- ids

  out <- list(
    tables = list(syn = read_table('syn'), da = read_table('da')),
    meta   = readRDS(file.path(dir, 'meta.rds')),
    cases  = cases )
  saveRDS(out, file.path('tests', 'tinytests', 'expected', 'log_ratio_plots.rds'), version = 2)
}
