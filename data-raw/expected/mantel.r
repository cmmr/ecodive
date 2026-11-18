# Expected results for tests/tinytests/test-mantel.r and
# test-mantel_partial.r, from vegan::mantel() and vegan::mantel.partial().
# Run from the package root:
#
#     Rscript data-raw/expected/mantel.r
#
# Writes tests/tinytests/expected/mantel.rds.
#
# The tests give mantel() and mantel_partial() the same permutation matrix
# as vegan, so p-values (and every permuted statistic) must match, not just
# statistically.

library(vegan)

set.seed(1)
n <- 30

counts <- matrix(rpois(n * 25, 4) * rbinom(n * 25, 1, 0.5), n, 25)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:25))

# Environment: a pH gradient that shifts OTUs 2-5, so some tests are
# significant. Whole-number ages give many tied distances (as in rbiom's
# Mantel test of dist(age)).
env <- data.frame(
  row.names = rownames(counts),
  ph  = round(rnorm(n, 7, 0.5), 2),
  age = sample(20:30, n, replace = TRUE),
  lat = round(runif(n, 30, 35), 1),
  lon = round(runif(n, -100, -95), 1) )
shift <- round(12 * (env$ph - min(env$ph)) / diff(range(env$ph)))
counts[, 2:5] <- counts[, 2:5] + shift

x <- vegdist(counts, 'bray')                          # community
y <- dist(env[, c('ph', 'lat')])                      # environment
z <- dist(env[, c('lat', 'lon')])                     # geography
a <- dist(setNames(env$age, rownames(env)))           # ties

perms <- shuffleSet(n, 199)
attributes(perms) <- list(dim = dim(perms))

# One case: vegan's statistic, p-value and permuted statistics.
run <- function (ydis, zdis, method) {
  res <- if (is.null(zdis))
    mantel(x, ydis, method = method, permutations = perms)
  else
    mantel.partial(x, ydis, zdis, method = method, permutations = perms)
  list(statistic = res$statistic, signif = res$signif, perm = unname(res$perm))
}

cases <- list()
for (method in c('pearson', 'spearman', 'kendall')) {
  cases[[length(cases) + 1]] <- c(list(y = 'y', z = NULL,  method = method), run(y, NULL, method))
  cases[[length(cases) + 1]] <- c(list(y = 'a', z = NULL,  method = method), run(a, NULL, method))
  cases[[length(cases) + 1]] <- c(list(y = 'y', z = 'z',   method = method), run(y, z,    method))
  cases[[length(cases) + 1]] <- c(list(y = 'a', z = 'y',   method = method), run(a, y,    method))
}

saveRDS(
  object = list(
    x        = x,
    y        = y,
    z        = z,
    a        = a,
    perms    = perms,
    cases    = cases,
    versions = list(vegan = as.character(packageVersion('vegan'))) ),
  file     = 'tests/tinytests/expected/mantel.rds',
  version  = 2 )
