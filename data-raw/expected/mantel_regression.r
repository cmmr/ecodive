# Expected results for tests/tinytests/test-mantel_regression.r, from
# ecodist::MRM() (ecodist 2.1.3) on the oribatid mites (vegan's `mite`,
# `mite.env` and `mite.xy`). Run from the package root:
#
#     Rscript data-raw/expected/mantel_regression.r
#
# Writes tests/tinytests/expected/mantel_regression.rds.
#
# ecodist is GPL: only its outputs are used (see plans/mantel_regression.md).
# MRM() has no argument for a permutation matrix, and draws its permutations
# in C from R's random number generator. From its outputs alone (the number
# of uniform draws it consumes, and which candidate shuffles reproduce its
# p-values), its permutations with `nperm` are nperm - 1 cumulative shuffles
# of the response's samples, each drawn as below; the observed samples count
# as the first of `nperm`. ecodist_perms() rebuilds them, one row per
# permutation of the response's samples, so that mantel_regression() given
# that matrix must reproduce MRM()'s p-values exactly, not just
# statistically.

library(ecodist)
library(vegan)  # the data

data(mite)
data(mite.env)
data(mite.xy)

ecodist_perms <- function (n, nperm, seed) {
  set.seed(seed)
  cum <- seq_len(n)
  P   <- matrix(0L, nperm - 1L, n)
  for (k in seq_len(nperm - 1L)) {
    p <- seq_len(n)
    for (i in n:2) {
      j    <- floor(runif(1) * (i - 1)) + 1
      tmp  <- p[i]; p[i] <- p[j]; p[j] <- tmp
    }
    cum     <- cum[p]
    P[k, ]  <- cum
  }
  P
}

n     <- nrow(mite)
nperm <- 200L
seed  <- 1L
perms <- ecodist_perms(n, nperm, seed)

x     <- vegdist(decostand(mite, 'hellinger'), 'euclidean')  # community
geo   <- dist(as.matrix(mite.xy))                            # m
dens  <- dist(mite.env$SubsDens)
water <- dist(mite.env$WatrCont)
# 0 for two samples with the same microtopography, 1 for different.
topo  <- as.dist(outer(mite.env$Topo, mite.env$Topo, '!=') * 1)
shrub <- as.dist(outer(mite.env$Shrub, mite.env$Shrub, '!=') * 1)


# One case: MRM()'s coefficients, p-values, R^2 and F.
run <- function (f, mrank) {
  set.seed(seed)
  res <- MRM(f, nperm = nperm, mrank = mrank)
  list(
    formula   = deparse(f),
    mrank     = mrank,
    estimate  = unname(res$coef[, 1]),
    p_value   = unname(res$coef[, 2]),
    r_squared = unname(res$r.squared[[1]]),
    r2_p      = unname(res$r.squared[[2]]),
    F         = unname(res$F.test[[1]]),
    F_p       = unname(res$F.test[[2]]) )
}

cases <- list(
  run(x ~ geo + dens + water, FALSE),
  run(x ~ geo + dens + water, TRUE),
  run(x ~ geo, FALSE),
  run(x ~ water + topo + shrub, FALSE),
  run(x ~ geo + water + topo, TRUE) )

saveRDS(
  object = list(
    x     = x,
    xy    = as.matrix(mite.xy),
    env   = mite.env[, c('SubsDens', 'WatrCont', 'Shrub', 'Topo')],
    perms = perms,
    cases = cases ),
  file     = 'tests/tinytests/expected/mantel_regression.rds',
  version  = 2,
  compress = 'xz' )
