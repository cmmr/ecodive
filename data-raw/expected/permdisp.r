# Expected results for tests/tinytests/test-permdisp.r, from
# vegan::betadisper() and vegan::permutest(). Run from the package root:
#
#     Rscript data-raw/expected/permdisp.r
#
# Writes tests/tinytests/expected/permdisp.rds.
#
# The tests give permdisp() the same permutation matrix as vegan, so p-values
# must match exactly, not just statistically.
#
# betadisper(type = "median") stops optim()'s BFGS early (to ~1e-5), so for
# medians this also stores `exact`: the distances to spatial medians solved
# to full precision (Weiszfeld's algorithm, in vegan's own principal
# coordinates, separately on the real and imaginary axes), which permdisp()
# must match to rounding.

library(vegan)

set.seed(1)
n <- 36

counts <- matrix(rpois(n * 30, 4) * rbinom(n * 30, 1, 0.5), n, 30)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:30))
storage.mode(counts) <- 'integer'

# Unequal groups; 'high' is more dispersed.
meta <- data.frame(
  row.names = rownames(counts),
  trt   = factor(rep(c('ctrl', 'low', 'high'), c(14, 12, 10)), c('ctrl', 'low', 'high')),
  block = factor(rep(1:6, length.out = n)),
  solo  = c(rep(c('a', 'b', 'c'), c(15, 12, 8)), 'd'),
  half  = rep(c('x', 'y'), each = n / 2) )
hi <- meta$trt == 'high'
counts[hi, ] <- counts[hi, ] * matrix(rpois(sum(hi) * 30, 1), sum(hi), 30)
counts[rowSums(counts) == 0, 1] <- 1L

dists <- list(
  bray      = vegdist(counts, 'bray'),
  euclidean = vegdist(counts, 'euclidean') )

perms <- shuffleSet(n, 99)


# Spatial median to full precision, with Vardi & Zhang's step at data points.
weiszfeld <- function (X) {
  m <- colMeans(X)
  if (nrow(X) <= 2) return(m)
  if (ncol(X) == 1) return(median(X))
  for (it in 1:10000) {
    dd <- sqrt(rowSums(sweep(X, 2, m)^2))
    z  <- dd == 0
    w  <- 1 / dd[!z]
    T  <- colSums(X[!z, , drop = FALSE] * w) / sum(w)
    if (any(z)) {
      r <- sqrt(sum(((T - m) * sum(w))^2))
      g <- min(1, sum(z) / r)
      T <- (1 - g) * T + g * m
    }
    step <- sqrt(sum((T - m)^2))
    m    <- T
    if (step <= 1e-15 * mean(dd)) break
  }
  m
}

exact_median <- function (bd) {
  v   <- bd$vectors
  pos <- bd$eig > 0
  g   <- bd$group
  z2  <- numeric(nrow(v))
  for (lv in levels(g)) {
    i <- which(g == lv)
    for (s in c(TRUE, FALSE)) {
      ax <- if (s) pos else !pos
      if (!any(ax)) next
      X  <- v[i, ax, drop = FALSE]
      dq <- rowSums(sweep(X, 2, weiszfeld(X))^2)
      z2[i] <- z2[i] + if (s) dq else -dq
    }
  }
  sqrt(pmax(z2, 0))
}


cases <- list()

for (dn in names(dists)) {
  for (type in c('centroid', 'median')) {
    for (bias in c(FALSE, TRUE)) {
      for (grp in c('trt', 'solo')) {

        if (grp == 'solo' && (bias || dn == 'euclidean')) next

        bd <- betadisper(dists[[dn]], meta[[grp]], type = type, bias.adjust = bias)
        pt <- permutest(bd, permutations = perms, pairwise = TRUE)

        cases[[length(cases) + 1]] <- list(
          dist        = dn,
          type        = type,
          bias        = bias,
          groups      = grp,
          distances   = unname(bd$distances),
          dispersions = bd$group.distances,
          exact       = if (type == 'median' && !bias) exact_median(bd),
          table       = unname(as.matrix(pt$tab)),
          f_perm      = unname(pt$perm[, 1]),
          t_perm      = if (!bias) unname(pt$perm[, -1]),
          t_obs       = unname(pt$statistic[-1]),
          pairwise    = pt$pairwise )
      }
    }
  }
}


# Free permutations within blocks: permutest's p-value from 9999 random
# permutations, for a statistical comparison.
set.seed(2)
bd <- betadisper(dists$bray, meta$half, type = 'centroid')
strata_p <- permutest(bd, permutations = how(nperm = 9999, blocks = meta$block))$tab[1, 'Pr(>F)']

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
saveRDS(expected, 'tests/tinytests/expected/permdisp.rds', version = 2)
