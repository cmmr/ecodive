# Writes tests/tinytests/expected/da_plots.rds: inputs and reference
# results for test-da_plots.r, so that the tests need no reference package.
#
# Run from the package root, with ecodive, ALDEx3 1.4.0 (CRAN, MIT) and
# ALDEx2 1.34.0 (Bioconductor 3.18, GPL >= 3: run for its outputs only)
# installed:
#
#     Rscript data-raw/expected/da_plots.r
#
# * ALDEx3's aldex.effect() (estimate, pooled.SD, cohens.d) for the binary
#   model ~ grp, with ALDEx3's own Dirichlet draws written out as the
#   Gamma draws that dirmult_ttest_details() takes (proportions), with no
#   scale model and with a scale model whose shifts are ecodive's: ALDEx3
#   is given the scale matrix -colMeans(logComp) + x * shift, which is
#   clr.sm()'s with the treatment's shift fixed. The two aldex() calls
#   start from the same R seed, so they draw the same compositions (checked).
# * ALDEx3's aldex.effect() with its own 2000 draws on a deep table (counts
#   of about 1e4, so the draws barely vary), and ALDEx2's aldex.effect()
#   (rab.all, diff.btw, diff.win, effect) with its own 1000 draws on it,
#   for comparisons by distribution.

library(ecodive)
library(ALDEx3)
library(ALDEx2)

stopifnot(
  packageVersion('ALDEx3') == '1.4.0',
  packageVersion('ALDEx2') == '1.34.0' )

out <- file.path('tests', 'tinytests', 'expected', 'da_plots.rds')


#--------------------------------------------------------
# Tables: samples x features, two groups of six.
#--------------------------------------------------------
set.seed(20261120)
n_f    <- 20L
groups <- rep(c('a', 'b'), each = 6)
mu     <- rlnorm(n_f, 4, 1)
fc     <- c(4, 3, 0.25, 0.5, rep(1, n_f - 4))
lam    <- outer(ifelse(groups == 'b', 1, 0), log(fc))
counts <- matrix(
  rpois(12 * n_f, exp(log(rep(mu, each = 12)) + lam + rnorm(12 * n_f, 0, 0.4))),
  12, n_f, dimnames = list(paste0('S', 1:12), paste0('F', 1:n_f)) )
storage.mode(counts) <- 'integer'

deep <- matrix(
  rpois(12 * n_f, 1e4 * exp(log(rep(mu / mean(mu), each = 12)) + lam + rnorm(12 * n_f, 0, 0.4))),
  12, n_f, dimnames = dimnames(counts) )
storage.mode(deep) <- 'integer'

meta <- data.frame(grp = factor(groups), row.names = rownames(counts))
xb   <- as.numeric(groups == 'b')   # model.matrix(~ grp)'s grpb column


#--------------------------------------------------------
# ALDEx3 on the same draws.
#--------------------------------------------------------
n_draws <- 32L

aldex3_case <- function (scale_sd, seed) {

  set.seed(1)
  m0 <- ALDEx3::aldex(
    t(counts), ~ grp, meta, nsample = n_draws, scale = ALDEx3::clr.sm,
    gamma = 0, n.cores = 1 )
  lc <- m0$logComp                                   # D x N x S, log2
  G  <- aperm(2^lc, c(2, 1, 3))                      # N x D x S
  dimnames(G) <- NULL

  # ecodive's shifts for these draws (natural logs), as log2.
  shift <- rep(0, n_draws)
  if (scale_sd > 0)
    shift <- ecodive:::dirmult_ttest_details(
      counts, groups, contrast = c('b', 'a'), scale_sd = scale_sd,
      n_draws = n_draws, gamma_draws = G, seed = seed, cpus = 1 )$shift / log(2)

  sc <- -colMeans(lc, dims = 1) + outer(xb, shift)  # N x S
  dimnames(sc) <- NULL

  set.seed(1)
  m <- ALDEx3::aldex(
    t(counts), ~ grp, meta, nsample = n_draws, scale = sc, n.cores = 1 )
  stopifnot(identical(unname(m$logComp), unname(lc)))

  e <- ALDEx3::aldex.effect(m, 'grpb')
  list(
    scale_sd  = scale_sd,
    seed      = seed,
    gamma     = G,
    estimate  = e$estimate,
    pooled_sd = e$pooled.SD,
    cohens_d  = e$cohens.d )
}

aldex3 <- list(
  none  = aldex3_case(0,   0L),
  scale = aldex3_case(0.5, 3L) )


#--------------------------------------------------------
# The deep table: ALDEx3's and ALDEx2's own draws.
#--------------------------------------------------------
set.seed(2)
m3 <- ALDEx3::aldex(
  t(deep), ~ grp, meta, nsample = 2000, scale = ALDEx3::clr.sm, gamma = 0,
  n.cores = 1 )
e3 <- ALDEx3::aldex.effect(m3, 'grpb')

set.seed(3)
cl <- ALDEx2::aldex.clr(
  as.data.frame(t(deep)), groups, mc.samples = 1000, denom = 'all',
  verbose = FALSE, useMC = FALSE )
e2 <- ALDEx2::aldex.effect(cl, verbose = FALSE, include.sample.summary = FALSE)

saveRDS(
  list(
    counts      = counts,
    deep        = deep,
    groups      = groups,
    aldex3      = aldex3,
    aldex3_deep = list(
      estimate  = e3$estimate,
      pooled_sd = e3$pooled.SD,
      cohens_d  = e3$cohens.d ),
    aldex2      = as.list(e2[, c('rab.all', 'diff.btw', 'diff.win', 'effect')]) ),
  out, version = 2 )

message('Wrote ', out)
