# Expected results for tests/tinytests/test-pcoa_correct.r: the Lingoes and
# Cailliez corrections from ape::pcoa() and vegan::wcmdscale(), whose
# definitions match pcoa()'s, and the "ring" distance profiles from
# data-raw/expected/cailliez_rings.r. Run from the package root:
#
#     Rscript data-raw/expected/pcoa_correct.r
#
# Writes tests/tinytests/expected/pcoa_correct.rds.

library(ecodive)

rand_counts <- function (n, p) {
  m <- matrix(rpois(n * p, 3) * rbinom(n * p, 1, 0.4), n, p)
  m[, 1] <- m[, 1] + 1L
  storage.mode(m) <- 'integer'
  dimnames(m) <- list(paste0('S', seq_len(n)), paste0('OTU', seq_len(p)))
  m
}

set.seed(1)
dists <- list(
  'bray, n = 60'       = bray(rand_counts(60, 40)),          # many negative eigenvalues
  'jaccard, ties'      = jaccard(rand_counts(80, 6) > 0),    # heavily tied distances
  'bray, n = 100'      = bray(rand_counts(100, 30)),
  'jaccard, n = 30'    = jaccard(rand_counts(30, 40)) )

# ape's and vegan's first three corrected axes.
refs <- function (d, corr) {

  ap <- ape::pcoa(d, correction = corr)
  vg <- vegan::wcmdscale(d, k = 3, add = corr, eig = TRUE)
  sm <- summary(vegan::eigenvals(vg))

  list(
    ape = list(
      eig      = ap$values$Corr_eig[1:3],
      rel_eig  = ap$values$Rel_corr_eig[1:3],
      vectors  = unname(ap$vectors.cor[, 1:3]) ),
    vegan = list(
      constant = vg$ac,
      eig      = vg$eig[1:3],
      prop     = unname(sm['Proportion Explained', 1:3]),
      points   = unname(vg$points[, 1:3]) ))
}

cases <- lapply(names(dists), function (nm) {
  d <- dists[[nm]]
  list(
    name     = nm,
    dist     = d,
    lingoes  = refs(d, 'lingoes'),
    cailliez = refs(d, 'cailliez') )
})


# Ring profiles, from their hexadecimal floats.
rings <- read.csv(
  file.path('data-raw', 'expected', 'pcoa_correct', 'cailliez_rings.csv'),
  stringsAsFactors = FALSE )
rings <- lapply(
  X   = split(rings$h, rings$g),
  FUN = function (h) as.numeric(h) )

stopifnot(identical(names(rings), c('120', '200')))


expected <- list(

  versions = list(
    R     = R.version.string,
    ape   = as.character(packageVersion('ape')),
    vegan = as.character(packageVersion('vegan')) ),

  # Distances, and ape's and vegan's results for each correction.
  cases = cases,

  # Profiles h[1..g/2] of the ring distances for g = 120 and 200 points.
  rings = rings
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/pcoa_correct.rds', version = 2)
