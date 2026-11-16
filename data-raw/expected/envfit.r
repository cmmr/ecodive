# Expected results for tests/tinytests/test-envfit.r, from vegan::envfit().
# Run from the package root:
#
#     Rscript data-raw/expected/envfit.r
#
# Writes tests/tinytests/expected/envfit.rds.
#
# The tests give envfit() the same permutation matrices as vegan, so p-values
# must match exactly, not just statistically.

library(vegan)

data(varespec, varechem, package = 'vegan')

set.seed(1)
n <- nrow(varespec)

# Two ordinations: PCoA axes (orthogonal), and a skewed copy of them whose
# axes are correlated, as NMDS or t-SNE axes can be.
pco  <- cmdscale(vegdist(varespec), k = 3)
colnames(pco) <- paste0('PCoA', 1:3)
skew <- pco %*% matrix(c(1, 0.6, -0.3, 0.2, 1, 0.5, 0, -0.4, 1), 3)
dimnames(skew) <- list(rownames(pco), paste0('Ax', 1:3))

env <- data.frame(
  varechem[, c('N', 'P', 'K', 'Ca', 'pH', 'Humdepth')],
  grazed = factor(rep(c('no', 'yes'), length.out = n)),
  soil   = cut(varechem$Baresoil, 3, labels = c('low', 'mid', 'high')),
  zone   = as.character(rep(c('A', 'B', 'C', 'D'), each = n / 4)),
  check.names = FALSE )
env$soil <- as.character(env$soil)  # character columns become factors

perms <- shuffleSet(n, 199)

# The same samples with NAs, for na.rm = TRUE.
env_na <- env
env_na$P[c(3, 10)]  <- NA
env_na$zone[17]     <- NA
keep_na  <- complete.cases(env_na)
perms_na <- shuffleSet(sum(keep_na), 199)

# Count tables, transformed as envfit()'s `norm` does.
counts  <- round(as.matrix(varespec) * 10)
storage.mode(counts) <- 'integer'
percent <- counts / rowSums(counts)
clr     <- log(counts + 1)
clr     <- clr - rowMeans(clr)

fit <- function (ord, P, choices, perms) {
  f <- envfit(ord, P, permutations = perms, choices = choices)
  list(
    arrows    = if (!is.null(f$vectors)) unclass(f$vectors$arrows),
    r         = if (!is.null(f$vectors)) f$vectors$r,
    pvals     = if (!is.null(f$vectors)) f$vectors$pvals,
    centroids = if (!is.null(f$factors)) unclass(f$factors$centroids),
    f_r       = if (!is.null(f$factors)) f$factors$r,
    f_pvals   = if (!is.null(f$factors)) f$factors$pvals )
}

cases <- list()

for (ord in c('pco', 'skew')) {
  for (ch in list(1:2, 1:3)) {
    lab <- paste(ord, paste(ch, collapse = ''))
    cases[[paste('env', lab)]] <- c(
      list(ord = ord, data = 'env', axes = ch, args = list()),
      fit(get(ord), env, ch, perms) )
    cases[[paste('counts', lab)]] <- c(
      list(ord = ord, data = 'counts', axes = ch, args = list()),
      fit(get(ord), counts, ch, perms) )
    cases[[paste('percent', lab)]] <- c(
      list(ord = ord, data = 'counts', axes = ch, args = list(norm = 'percent')),
      fit(get(ord), percent, ch, perms) )
    cases[[paste('clr', lab)]] <- c(
      list(ord = ord, data = 'counts', axes = ch, args = list(norm = 'clr', pseudocount = 1)),
      fit(get(ord), clr, ch, perms) )
  }
}

f <- envfit(pco, env_na, permutations = perms_na, na.rm = TRUE)
cases[['env_na']] <- list(
  ord = 'pco', data = 'env_na', axes = 1:2, args = list(na.rm = TRUE),
  arrows = unclass(f$vectors$arrows), r = f$vectors$r, pvals = f$vectors$pvals,
  centroids = unclass(f$factors$centroids), f_r = f$factors$r, f_pvals = f$factors$pvals )

out <- list(
  versions = list(
    R     = as.character(getRversion()),
    vegan = as.character(packageVersion('vegan')),
    permute = as.character(packageVersion('permute')) ),
  pco      = pco,
  skew     = skew,
  env      = env,
  env_na   = env_na,
  counts   = counts,
  perms    = unclass(perms)[seq_len(nrow(perms)), , drop = FALSE],
  perms_na = unclass(perms_na)[seq_len(nrow(perms_na)), , drop = FALSE],
  cases    = cases )

attributes(out$perms)    <- list(dim = dim(perms))
attributes(out$perms_na) <- list(dim = dim(perms_na))

saveRDS(out, 'tests/tinytests/expected/envfit.rds', version = 2)
