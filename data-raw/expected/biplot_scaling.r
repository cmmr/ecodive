# Writes tests/tinytests/expected/biplot_scaling.rds: reference biplot
# coordinates for the "Biplot scalings" section of test-plot.r, so that the
# tests need no reference package.
#
# Run from the package root:  Rscript data-raw/expected/biplot_scaling.r
#
# vegan, easyCODA and compositions are GPL: only their outputs are used,
# never their source. Generated with vegan 2.7.8 (permute 0.10.1),
# easyCODA 0.50 and compositions 2.0-9, under R 4.3.3.
#
# Holds:
#
#   clr      A zero-free 25 x 7 table (`counts`), and for its clr transform:
#            vegan::rda() of the clr values, `scores(scaling = 1)` and
#            `scores(scaling = 2)` of sites and species on all 6 axes, with
#            `tot_chi` (its total inertia); easyCODA's unweighted PCA()
#            of CLR(counts): `rowcoord`, `rowpcoord`, `colcoord`,
#            `colpcoord`, `colmass` and `sv`; and compositions'
#            princomp(acomp(counts)): `sdev`, `loadings` and `scores`.
#   rda      vegan::rda() of the Hellinger-transformed dune table on
#            `~ Management + A1 + Condition(Moisture)` (and unconstrained,
#            `rda_u`): scores of sites ("wa"), species, biplot ("bp") and
#            centroids ("cn") on all axes in scalings 1, 2 and 3.
#   cca      vegan::cca(dune ~ Management + A1, dune.env): the same.
#   ca       vegan::cca(dune): sites and species in scalings 1, 2 and 3.
#
# Each of rda, rda_u, cca and ca also holds `k`, the number of axes; `dune`
# and `dune_env` are their inputs (vegan's dune and dune.env).

lib <- Sys.getenv('ECODIVE_REF_LIB')   # where vegan 2.7.8 etc. are installed
if (nzchar(lib)) .libPaths(c(strsplit(lib, ':')[[1]], .libPaths()))

library(vegan)
library(easyCODA)
suppressPackageStartupMessages(library(compositions))

stopifnot(packageVersion('vegan') == '2.7.8')

data(dune, dune.env)
dune <- as.matrix(dune)


#--------------------------------------------------
# A zero-free composition.
#--------------------------------------------------
set.seed(20261124)
n <- 25; p <- 7
mu     <- c(5, 4, 4, 3, 3, 2, 1)
counts <- t(vapply(seq_len(n), function (i) {
  g <- if (i <= 12) c(0, 1, 0, -1, 0, 0.5, 0) else c(1, 0, -1, 0, 0.5, 0, 0)
  round(exp(mu + g + rnorm(p, 0, 0.4))) + 1
}, numeric(p)))
dimnames(counts) <- list(sprintf('S%02d', seq_len(n)), sprintf('F%d', seq_len(p)))
stopifnot(all(counts > 0))

clr_x <- log(counts) - rowMeans(log(counts))
v     <- vegan::rda(clr_x)
k     <- p - 1L
sc1   <- scores(v, choices = 1:k, display = c('sites', 'species'), scaling = 1)
sc2   <- scores(v, choices = 1:k, display = c('sites', 'species'), scaling = 2)

e  <- easyCODA::PCA(easyCODA::CLR(counts), nd = k, weight = FALSE, msg = FALSE)
pc <- compositions::princomp.acomp(compositions::acomp(counts))

clr <- list(
  counts = counts,
  vegan  = list(
    sites1 = unclass(sc1$sites), species1 = unclass(sc1$species),
    sites2 = unclass(sc2$sites), species2 = unclass(sc2$species),
    tot_chi = v$tot.chi ),
  easycoda = list(
    rowcoord  = unname(e$rowcoord[, 1:k]),  rowpcoord = unname(e$rowpcoord[, 1:k]),
    colcoord  = unname(e$colcoord[, 1:k]),  colpcoord = unname(e$colpcoord[, 1:k]),
    colmass   = e$colmass, sv = e$sv ),
  compositions = list(
    sdev     = unname(pc$sdev),
    loadings = unname(unclass(pc$loadings)),
    scores   = unname(unclass(pc$scores)) ))


#--------------------------------------------------
# vegan's scalings of rda(), cca() and ca().
#--------------------------------------------------
vegan_scores <- function (fit, display) {
  k   <- length(fit$CCA$eig) + length(fit$CA$eig)
  out <- list(k = k)
  for (s in 1:3) {
    sc <- scores(fit, choices = 1:k, display = display, scaling = s)
    if (!is.list(sc)) sc <- setNames(list(sc), display)
    out[[paste0('scaling', s)]] <- lapply(sc, function (m) unclass(m))
  }
  out
}

hel <- decostand(dune, 'hellinger')
dsp <- c('wa', 'species', 'bp', 'cn')

res <- list(
  clr      = clr,
  dune     = dune,
  dune_env = dune.env,
  rda   = vegan_scores(vegan::rda(hel ~ Management + A1 + Condition(Moisture), dune.env), dsp),
  rda_u = vegan_scores(vegan::rda(hel), c('sites', 'species')),
  cca   = vegan_scores(vegan::cca(dune ~ Management + A1, dune.env), dsp),
  ca    = vegan_scores(vegan::cca(dune), c('sites', 'species')) )

saveRDS(res, 'tests/tinytests/expected/biplot_scaling.rds', version = 2)
