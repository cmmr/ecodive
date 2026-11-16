# Expected results for tests/tinytests/test-permanova.r, from
# vegan::adonis2() and vegan::anova.cca() on vegan::dbrda(). Run from the
# package root:
#
#     Rscript data-raw/expected/permanova.r
#
# Writes tests/tinytests/expected/permanova.rds.
#
# The tests give permanova() the same permutation matrix as vegan, so p-values
# (and every permuted F) must match exactly, not just statistically.

library(ecodive)
library(vegan)

set.seed(1)
n <- 36

counts <- matrix(rpois(n * 30, 4) * rbinom(n * 30, 1, 0.5), n, 30)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:30))
storage.mode(counts) <- 'integer'

# A treatment effect on the first few OTUs, so that some p-values are small.
meta <- data.frame(
  row.names = rownames(counts),
  trt   = factor(rep(c('ctrl', 'low', 'high'), length.out = n), c('ctrl', 'low', 'high')),
  site  = factor(rep(c('A', 'B'), each = n / 2)),
  ph    = round(rnorm(n, 7, 0.5), 2),
  depth = round(runif(n, 1, 10), 1),
  block = factor(rep(1:6, each = n / 6)) )
counts[meta$trt == 'high', 2:4] <- counts[meta$trt == 'high', 2:4] + 6L

d     <- bray(counts)
perms <- shuffleSet(n, 199)  # adonis2() needs its "control" attribute

# A table as a plain list, without vegan's classes and attributes.
tbl <- function (x, f_perm = TRUE) {
  list(
    rows   = rownames(x),
    table  = unname(as.matrix(x)),
    f_perm = if (f_perm) unname(attr(x, 'F.perm')) )
}

cases <- list()

# adonis2: each `by` for several designs.
formulas <- c(
  '~ trt',
  '~ ph',
  '~ trt + ph',
  '~ trt * site',
  '~ site + trt * ph + depth' )

for (f in formulas) {
  for (by in c('overall', 'terms', 'margin', 'onedf')) {
    a <- adonis2(
      stats::as.formula(paste('d', f)), data = meta, permutations = perms,
      by = if (by == 'overall') NULL else by )
    cases[[length(cases) + 1]] <- c(
      list(formula = f, by = by, ref = 'adonis2'),
      tbl(a, f_perm = (f == '~ trt * site')) )
  }
}

# Partial models: anova(dbrda(...)) has no Total row and no R2 column.
for (f in c('~ trt + Condition(site)', '~ trt * ph + Condition(site + depth)')) {
  m <- dbrda(stats::as.formula(paste('d', f)), data = meta)
  for (by in c('overall', 'terms', 'margin')) {
    a <- anova(m, permutations = perms, by = if (by == 'overall') NULL else by)
    cases[[length(cases) + 1]] <- c(
      list(formula = f, by = by, ref = 'dbrda'),
      tbl(a, f_perm = FALSE) )
  }
}

# Free permutations within blocks: adonis2's p-value from 9999 random
# permutations, for a statistical comparison.
set.seed(2)
strata_p <- adonis2(d ~ ph, data = meta, permutations = 9999, strata = meta$block)[1, 'Pr(>F)']

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
saveRDS(expected, 'tests/tinytests/expected/permanova.rds', version = 2)
