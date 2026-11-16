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

# More cases. Each names its `set` of distances, metadata and permutations:
#   'meta2'   bray(counts), `meta2`, and `perms`.
#   'blocks'  bray(counts), `meta`, and `block_perms`: permutations within
#             blocks.
#   'n500'    bray(counts500), `meta500`, and `perms500`.
# By default, the reference is adonis2().

set.seed(3)

# Ordered, character and logical columns, and an unbalanced factor with an
# unused level.
meta2 <- data.frame(
  row.names = rownames(counts),
  moist = factor(sample(c('low', 'mid', 'high'), n, TRUE), c('low', 'mid', 'high'), ordered = TRUE),
  ph    = meta$ph,
  chr   = sample(c('a', 'b', 'c'), n, TRUE),
  lgl   = sample(c(TRUE, FALSE), n, TRUE),
  unbal = factor(rep(c('w', 'x', 'y', 'z'), c(20, 10, 4, 2)), c('w', 'x', 'unused', 'y', 'z')),
  stringsAsFactors = FALSE )

block_perms <- shuffleSet(n, 199, control = how(blocks = meta$block))

# 500 samples and one factor: permanova()'s fast path for a single factor.
counts500 <- matrix(rpois(500 * 10, 2), 500, 10)
counts500[, 1] <- counts500[, 1] + 1L
dimnames(counts500) <- list(sprintf('S%03d', 1:500), paste0('OTU', 1:10))
storage.mode(counts500) <- 'integer'
meta500 <- data.frame(
  row.names = rownames(counts500),
  g = factor(sample(c('a', 'b', 'c', 'd'), 500, TRUE, prob = c(0.5, 0.3, 0.15, 0.05))),
  x = round(runif(500), 3) )
counts500[meta500$g == 'a', 2] <- counts500[meta500$g == 'a', 2] + 1L
perms500 <- shuffleSet(500, 9)

sets <- list(
  meta2  = list(d = d,               data = meta2,   perms = perms),
  blocks = list(d = d,               data = meta,    perms = block_perms),
  n500   = list(d = bray(counts500), data = meta500, perms = perms500) )

more_cases <- list()
more_case  <- function (set, f, by, f_perm = FALSE) {
  s <- sets[[set]]
  a <- adonis2(
    stats::as.formula(paste('s$d', f)), data = s$data, permutations = s$perms,
    by = if (by == 'overall') NULL else by )
  more_cases[[length(more_cases) + 1]] <<- c(
    list(set = set, formula = f, by = by, ref = 'adonis2'),
    tbl(a, f_perm = f_perm) )
}

for (by in c('overall', 'terms', 'onedf')) more_case('meta2', '~ moist', by)
for (by in c('terms', 'onedf'))            more_case('meta2', '~ poly(ph, 2)', by)
more_case('meta2', '~ chr', 'overall', f_perm = TRUE)
more_case('meta2', '~ lgl', 'overall')
for (by in c('terms', 'margin'))           more_case('meta2', '~ chr + lgl', by)
for (by in c('overall', 'onedf'))          more_case('meta2', '~ unbal', by)
for (by in c('terms', 'margin'))           more_case('meta2', '~ moist * lgl + poly(ph, 2)', by)
more_case('blocks', '~ trt', 'overall', f_perm = TRUE)
for (by in c('terms', 'margin'))           more_case('blocks', '~ ph + trt', by)
more_case('n500', '~ g', 'overall', f_perm = TRUE)
more_case('n500', '~ g + x', 'terms')

expected <- list(

  versions = list(
    R       = R.version.string,
    vegan   = as.character(packageVersion('vegan')),
    permute = as.character(packageVersion('permute')) ),

  counts   = counts,
  meta     = meta,
  perms    = matrix(as.integer(perms), nrow(perms)),
  cases    = cases,
  strata_p = strata_p,

  # The cases above, and their inputs.
  meta2       = meta2,
  block_perms = matrix(as.integer(block_perms), nrow(block_perms)),
  counts500   = counts500,
  meta500     = meta500,
  perms500    = matrix(as.integer(perms500), nrow(perms500)),
  more_cases  = more_cases )

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/permanova.rds', version = 2)
