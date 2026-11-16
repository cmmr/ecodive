# Expected results for tests/tinytests/test-anosim.r, from vegan::anosim().
# Run from the package root:
#
#     Rscript data-raw/expected/anosim.r
#
# Writes tests/tinytests/expected/anosim.rds.
#
# The tests give anosim() the same permutation matrix as vegan, so p-values
# (and every permuted R) must match exactly, not just statistically.

library(ecodive)
library(vegan)

set.seed(1)
n <- 30

counts <- matrix(rpois(n * 25, 4) * rbinom(n * 25, 1, 0.5), n, 25)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:25))
storage.mode(counts) <- 'integer'

# A treatment effect on the first few OTUs, so that some p-values are small.
meta <- data.frame(
  row.names = rownames(counts),
  trt   = factor(rep(c('ctrl', 'low', 'high'), length.out = n), c('ctrl', 'low', 'high')),
  site  = rep(c('A', 'B'), each = n / 2),
  uneq  = c('x', rep('y', 4), rep('z', 25)),  # a singleton group
  block = factor(rep(1:5, each = n / 5)) )
counts[meta$trt == 'high', 2:4] <- counts[meta$trt == 'high', 2:4] + 5L

d     <- bray(counts)
perms <- shuffleSet(n, 199)

# Mean rank of each class (Between, then each group), as anosim() reports.
class_means <- function (a) {
  m <- tapply(a$dis.rank, a$class.vec, mean)
  k <- tapply(a$dis.rank, a$class.vec, length)
  list(class = names(m), pairs = as.integer(k), mean_rank = as.numeric(m))
}

cases <- list()

add_case <- function (name, dm, groups) {
  a <- anosim(dm, groups, permutations = perms)
  cases[[name]] <<- list(
    dm        = dm,
    groups    = groups,
    statistic = a$statistic,
    signif    = a$signif,
    perm      = as.numeric(a$perm),
    ranks     = class_means(a) )
}

add_case('trt',        d,            meta$trt)
add_case('site_ties',  round(d, 1),  meta$site)   # heavy ties
add_case('uneq',       d,            meta$uneq)
add_case('all_tied',   d * 0 + 0.5,  meta$trt)    # every distance equal
add_case('jaccard',    jaccard(counts), meta$trt)

# Free permutations within blocks: anosim's p-value from 9999 random
# permutations, for a statistical comparison.
set.seed(2)
strata_p <- anosim(d, meta$site, permutations = 9999, strata = meta$block)$signif

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
saveRDS(expected, 'tests/tinytests/expected/anosim.rds', version = 2)
