# Generates tests/tinytests/expected/pcoa_lanczos.rds: pcoa() results from
# before its block Lanczos solver moved to src/lanczos.c (shared with pca()).
# test-pca.r checks that pcoa() still returns exactly these.
#
# The last bit of a result depends on how the C code was compiled, so the
# results of two builds of commit 66d7a80 (the last one before the move) are
# kept: optimized (as installed) and unoptimized (as devtools::load_all()
# compiles). From the package root, with that commit's source in <src>:
#
#   R CMD INSTALL --library=<lib> <src>
#   Rscript data-raw/expected/pcoa_lanczos.r O2 <lib>
#   Rscript data-raw/expected/pcoa_lanczos.r O0 <src>
#
# Each run adds or replaces its build's results in the file.

args  <- commandArgs(TRUE)
build <- args[1]
path  <- args[2]

if (build == 'O2') library(ecodive, lib.loc = path)
if (build == 'O0') devtools::load_all(path, quiet = TRUE)

set.seed(1)

rand_counts <- function (n, p) {
  m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
  m[, 1] <- m[, 1] + 1L
  storage.mode(m) <- 'integer'
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
  m
}

m150 <- rand_counts(150, 40)
m200 <- rand_counts(200, 60)

# Distances come from ecodive's own (deterministic) beta diversity code.
cases <- list(
  list(counts = m150, metric = 'bray',      k = 3L),
  list(counts = m150, metric = 'bray',      k = 12L),
  list(counts = m150, metric = 'jaccard',   k = 6L),
  list(counts = m200, metric = 'euclidean', k = 3L),
  list(counts = m200, metric = 'aitchison', k = 5L) )

file <- 'tests/tinytests/expected/pcoa_lanczos.rds'
old  <- if (file.exists(file)) readRDS(file) else NULL

for (i in seq_along(cases)) {
  cs <- cases[[i]]
  d  <- beta_div(cs$counts, cs$metric, pseudocount = 1)
  cases[[i]]$expected <- old$cases[[i]]$expected
  cases[[i]]$expected[[build]] <- pcoa(d, k = cs$k, cpus = 1)
}

generated <- old$generated
generated[[build]] <- list(
  ecodive  = paste(packageVersion('ecodive'), '(commit 66d7a80)'),
  R        = R.version.string,
  platform = R.version$platform )

saveRDS(
  file    = file,
  version = 2,
  object  = list(cases = cases, generated = generated) )
