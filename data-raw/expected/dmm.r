# Expected results for tests/tinytests/test-dmm.r, from Bioconductor's
# DirichletMultinomial 1.44.0 (Bioconductor 3.18, R 4.3.3, GSL 2.7.1). Run
# from the package root:
#
#     Rscript data-raw/expected/dmm.r
#
# Writes tests/tinytests/expected/dmm.rds.
#
# DirichletMultinomial is LGPL: it is only run here, to make fixtures (see
# plans/community_typing.md). Its source was not read.
#
# Cases:
#   sim:  120 samples x 30 features of multinomial counts (depths 500 to
#         3,000) from three Dirichlet types (50, 40 and 30 samples), with
#         precisions 40, 15 and 80.
#   dune: vegan's dune (20 sites x 30 species, abundance classes 0 to 9),
#         a small real table of whole numbers. (ecodive's ex_counts has only
#         6 features.)
#
# For each case and k = 1..4, dmn(counts, k, seed = 1): goodnessOfFit()
# (NLE, LogDet, Laplace, BIC, AIC), fitted() (features x k), mixturewt()$pi
# and mixture() (samples x k), in DirichletMultinomial's component order.
# The inputs are stored.

library(DirichletMultinomial)
data(dune, package = 'vegan')

set.seed(66)
S <- 30
n <- 120
rdir <- function (a) { g <- rgamma(length(a), a); g / sum(g) }
prof <- lapply(1:3, function (t) {
  p <- rexp(S)^2
  p[sample(S, 3)] <- p[sample(S, 3)] + 3
  p / sum(p) })
prec <- c(40, 15, 80)
grp  <- rep(1:3, c(50, 40, 30))
sim  <- t(sapply(grp, function (g) rmultinom(1, sample(500:3000, 1), rdir(prof[[g]] * prec[g]))))
dimnames(sim) <- list(sprintf('S%03d', 1:n), sprintf('F%02d', 1:S))
storage.mode(sim) <- 'integer'

dune <- as.matrix(dune)
storage.mode(dune) <- 'integer'

ref <- function (counts) {
  lapply(1:4, function (k) {
    f <- dmn(counts, k, seed = 1)
    list(
      k     = k,
      gof   = goodnessOfFit(f),
      alpha = unname(fitted(f)),
      pi    = unname(mixturewt(f)$pi),
      z     = unname(mixture(f)) )
  })
}

expected <- list(
  sim  = list(counts = sim,  groups = grp, fits = ref(sim)),
  dune = list(counts = dune, fits = ref(dune)) )

saveRDS(expected, 'tests/tinytests/expected/dmm.rds', version = 2)
