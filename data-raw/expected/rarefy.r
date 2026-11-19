# Writes tests/tinytests/expected/rarefy.rds: rarefy() results from ecodive
# 0.x as it was before rarefy(depth =) took one depth per sample (commit
# 4600cb8), so that test-rarefy.r can check that a single depth still gives
# bit-identical results after C_rarefy was generalised.
#
# Run from the package root, with that version of ecodive installed:
#
#     Rscript data-raw/expected/rarefy.r

library(ecodive)
library(Matrix)

gen <- function (n, d, seed) {
  set.seed(seed)
  m <- matrix(rpois(n * d, rexp(d, 1 / 8)), n, d)
  m[sample(length(m), length(m) %/% 3)] <- 0
  m[m > 0 & runif(length(m)) < 0.05] <- 250  # a few deep cells
  dimnames(m) <- list(sprintf('S%02d', seq_len(n)), sprintf('F%02d', seq_len(d)))
  m
}

x <- gen(30, 40, 1)

classes <- list(
  matrix = identity,
  dgC    = function (m) as(m, 'CsparseMatrix'),
  dgT    = function (m) as(m, 'TsparseMatrix'),
  dge    = function (m) as(as(m, 'unpackedMatrix'), 'generalMatrix'),
  slam   = slam::as.simple_triplet_matrix )

depths <- c(50L, 120L, as.integer(min(rowSums(x))))

res <- list()
for (cl in names(classes))
  for (margin in 1:2)
    for (depth in depths)
      for (seed in c(0L, 7L)) {
        m <- if (margin == 1L) x else t(x)
        r <- rarefy(classes[[cl]](m), depth = depth, seed = seed,
                    margin = margin, drop = FALSE, warn = FALSE, cpus = 1L)
        key <- sprintf('%s_m%d_d%d_s%d', cl, margin, depth, seed)
        res[[key]] <- as.matrix(r)
      }

saveRDS(
  list(counts = x, depths = depths, results = res),
  file.path('tests', 'tinytests', 'expected', 'rarefy.rds'),
  version = 2 )
