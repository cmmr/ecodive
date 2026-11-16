# Expected results for tests/tinytests/test-mrpp.r, from vegan::mrpp() and
# summary(vegan::permustats()). Run from the package root:
#
#     Rscript data-raw/expected/mrpp.r
#
# Writes tests/tinytests/expected/mrpp.rds.
#
# The tests give mrpp() the same permutation matrix as vegan, so every
# permuted delta, the p-value and the z-score must match, not just
# statistically.

library(ecodive)
library(vegan)

set.seed(1)
n <- 30

counts <- matrix(rpois(n * 25, 4) * rbinom(n * 25, 1, 0.5), n, 25)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:25))
storage.mode(counts) <- 'integer'

# Three groups with a treatment effect, and an unbalanced grouping with a
# group of one sample (which vegan leaves out of delta).
meta <- data.frame(
  row.names = rownames(counts),
  trt   = factor(rep(c('ctrl', 'low', 'high'), length.out = n), c('ctrl', 'low', 'high')),
  uneven = rep(c('a', 'b', 'c', 'd'), c(14, 9, 6, 1)),
  block = factor(rep(1:5, each = n / 5)) )
counts[meta$trt == 'high', 2:4] <- counts[meta$trt == 'high', 2:4] + 5L

d     <- bray(counts)
perms <- shuffleSet(n, 199)

cases <- list()
for (grp in c('trt', 'uneven')) {
  for (wt in 1:3) {
    m <- mrpp(d, meta[[grp]], permutations = perms, weight.type = wt)
    s <- summary(permustats(m))
    cases[[length(cases) + 1]] <- list(
      groups      = grp,
      weight_type = wt,
      delta       = m$delta,
      E_delta     = m$E.delta,
      A           = m$A,
      classdelta  = m$classdelta,
      n           = m$n,
      p           = m$Pvalue,
      boot_deltas = unname(m$boot.deltas),
      z           = unname(s$z) )
  }
}

# Free permutations within blocks: vegan's p-value from 9999 random
# permutations, for a statistical comparison.
set.seed(2)
strata_p <- mrpp(
  d, meta$trt, permutations = how(nperm = 9999, blocks = meta$block) )$Pvalue

expected <- list(

  versions = list(
    R       = R.version.string,
    vegan   = as.character(packageVersion('vegan')),
    permute = as.character(packageVersion('permute')) ),

  counts   = counts,
  meta     = meta,
  perms    = matrix(as.integer(perms), nrow(perms)),
  cases    = cases,
  strata_p = strata_p )

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/mrpp.rds', version = 2)
