# Generates tests/tinytests/expected/pcoa_lanczos.rds: pcoa() results from
# before its block Lanczos solver moved to src/lanczos.c (shared with pca()).
# test-pca.r checks that pcoa() still returns exactly these.
#
# Install ecodive from commit 66d7a80 (the last one before the move) into a
# separate library, then run from the package root:
#
#   R CMD INSTALL --library=<lib> .        (at commit 66d7a80)
#   Rscript data-raw/expected/pcoa_lanczos.r <lib>

library(ecodive, lib.loc = commandArgs(TRUE)[1])

set.seed(1)

rand_counts <- function (n, p) {
  m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
  m[, 1] <- m[, 1] + 1L
  storage.mode(m) <- 'integer'
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
  m
}

m150 <- rand_counts(150, 40)
m300 <- rand_counts(300, 80)

# Distances come from ecodive's own (deterministic) beta diversity code.
cases <- list(
  list(counts = m150, metric = 'bray',      k = 3L),
  list(counts = m150, metric = 'bray',      k = 12L),
  list(counts = m150, metric = 'jaccard',   k = 6L),
  list(counts = m300, metric = 'euclidean', k = 3L),
  list(counts = m300, metric = 'aitchison', k = 5L) )

for (i in seq_along(cases)) {
  cs <- cases[[i]]
  d  <- beta_div(cs$counts, cs$metric, pseudocount = 1)
  cases[[i]]$expected <- pcoa(d, k = cs$k, cpus = 1)
}

saveRDS(
  file    = 'tests/tinytests/expected/pcoa_lanczos.rds',
  version = 2,
  object  = list(
    cases     = cases,
    generated = list(
      ecodive  = paste(packageVersion('ecodive'), '(commit 66d7a80)'),
      R        = R.version.string,
      platform = R.version$platform )))
