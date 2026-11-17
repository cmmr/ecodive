# Expected results for tests/tinytests/test-pcoa_biplot.r, from scikit-bio's
# pcoa() and pcoa_biplot(). There is no Python here: data-raw/expected/
# pcoa_biplot.py wrote the CSV files in data-raw/expected/pcoa_biplot/ (inputs
# and outputs, with 17 significant digits), which are committed, and this
# script turns them into an .rds file. Run from the package root:
#
#     Rscript data-raw/expected/pcoa_biplot.r
#
# Writes tests/tinytests/expected/pcoa_biplot.rds.

dir <- file.path('data-raw', 'expected', 'pcoa_biplot')

rd <- function (name)
  as.matrix(read.csv(file.path(dir, name), row.names = 1, check.names = FALSE))

fits <- list()

for (nm in c('bray', 'jaccard', 'euclidean')) {

  eig <- rd(sprintf('skbio_%s_eig.csv', nm))

  fit <- list(
    dist   = as.dist(rd(sprintf('skbio_%s_dist.csv', nm))),
    eig    = unname(eig[, 'eig']),
    prop   = unname(eig[, 'prop']),
    coords = unname(rd(sprintf('skbio_%s_coords.csv', nm))) )

  if (nm != 'jaccard')
    fit$biplot <- lapply(
      X   = c(counts = 'counts', percent = 'percent', env = 'env'),
      FUN = function (lab) rd(sprintf('skbio_%s_biplot_%s.csv', nm, lab)) )

  fits[[nm]] <- fit
}

expected <- list(

  versions = list(
    R        = R.version.string,
    skbio    = '0.7.4 (data-raw/expected/pcoa_biplot.py)' ),

  # Inputs: 24 samples x 12 features (OTU11 all zero, OTU12 constant), and
  # four continuous variables.
  counts = rd('skbio_counts.csv'),
  env    = rd('skbio_env.csv'),

  # Per distance: scikit-bio's distances, eigenvalues, proportions, and
  # coordinates (5 axes), and biplot scores for the counts, their
  # proportions ('percent'), and the variables.
  fits = fits
)

dir.create('tests/tinytests/expected', showWarnings = FALSE)
saveRDS(expected, 'tests/tinytests/expected/pcoa_biplot.rds', version = 2)
