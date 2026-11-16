# Expected results for tests/tinytests/test-pcoa.r, from ape::pcoa() and
# vegan::wcmdscale(). Run from the package root:
#
#     Rscript data-raw/expected/pcoa.r
#
# Writes tests/tinytests/expected/pcoa.rds.

library(ecodive)

set.seed(1)
counts <- matrix(rpois(150 * 40, 3) * rbinom(150 * 40, 1, 0.4), 150, 40)
counts[, 1] <- counts[, 1] + 1
dimnames(counts) <- list(paste0('S', 1:150), paste0('OTU', 1:40))
storage.mode(counts) <- 'integer'

d  <- bray(counts)
ap <- ape::pcoa(d)
vg <- vegan::wcmdscale(d, k = 3, eig = TRUE)

expected <- list(

  versions = list(
    R     = R.version.string,
    ape   = as.character(packageVersion('ape')),
    vegan = as.character(packageVersion('vegan')) ),

  # Input: Bray-Curtis distances between the rows of `counts`.
  counts = counts,

  # ape::pcoa(d): first three axes. Bray-Curtis has negative eigenvalues.
  ape = list(
    vectors      = matrix(as.vector(ap$vectors[, 1:3]), nrow(ap$vectors)),
    eig          = ap$values$Eigenvalues[1:3],
    relative_eig = ap$values$Relative_eig[1:3],
    min_eig      = min(ap$values$Eigenvalues),
    trace        = ap$trace ),

  # vegan::wcmdscale(d, k = 3, eig = TRUE), and its "Proportion Explained".
  vegan = list(
    points = matrix(as.vector(vg$points), nrow(vg$points)),
    eig    = vg$eig[1:3],
    prop   = unname(summary(vegan::eigenvals(vg))[2, 1:3]) )
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/pcoa.rds', version = 2)
