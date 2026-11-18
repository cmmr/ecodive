# Expected results for tests/tinytests/test-procrustes.r, from
# vegan::procrustes() and vegan::protest() (GPL: outputs only), and from
# scipy.spatial.procrustes() (BSD-3). Run from the package root:
#
#     Rscript data-raw/expected/procrustes.r       # writes the inputs
#     python data-raw/expected/procrustes.py       # SciPy's disparities
#     Rscript data-raw/expected/procrustes.r       # folds them in
#
# Writes tests/tinytests/expected/procrustes.rds. Generated with vegan
# 2.7.8, permute 0.10.0 and SciPy 1.18.1.
#
# vegan's procrustes(X, Y) rotates Y onto X: ecodive's procrustes(Y, X).
# procrustes_test() is given the same permutation matrix as protest(), so
# every permuted statistic, and the p-value, must match.

library(vegan)

dir <- file.path('data-raw', 'expected', 'procrustes')
dir.create(dir, showWarnings = FALSE)

data(varespec)
d <- vegdist(varespec, 'bray')

# The target: a 3-axis PCoA. The ordination rotated onto it: a 2-axis
# NMDS of the same samples (so vegan pads it to 3 axes with zeros).
X <- cmdscale(d, k = 3)
set.seed(1)
Y <- metaMDS(d, k = 2, trace = 0)$points
attributes(Y) <- list(dim = dim(Y), dimnames = list(rownames(X), NULL))

# A reflected, rotated, dilated, and shifted copy of the first 2 axes of X
# (ss = 0 with dilation), and some unrelated noise (a p-value that isn't
# the smallest possible).
a   <- 0.7
rot <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
R   <- 2.5 * (X[, 1:2] %*% diag(c(1, -1))) %*% rot
R   <- sweep(R, 2, c(3, -1), '+')
N   <- matrix(round(rnorm(nrow(X) * 2), 3), nrow(X), 2, dimnames = list(rownames(X), NULL))

perms <- shuffleSet(nrow(X), 199)
attributes(perms) <- list(dim = dim(perms))

configs <- list(X = X, Y = Y, R = R, N = N)

# One case: ord rotated onto target, with axes (NULL: vegan's own default
# is the same when the widths agree).
fit <- function (ord, target, symmetric, scale) {
  v <- suppressWarnings(procrustes(configs[[target]], configs[[ord]], symmetric = symmetric, scale = scale))
  list(
    Yrot = unname(v$Yrot), X = matrix(v$X, nrow(v$X)), ss = v$ss, rotation = v$rotation,
    translation = drop(v$translation), scale = v$scale,
    residuals = unname(residuals(v)) )
}

cases <- list()
for (cs in list(
    list(ord = 'Y', target = 'X', axes = 1:3),           # widths 2 and 3
    list(ord = 'Y', target = 'X', axes = NULL, tk = 2),  # both 2
    list(ord = 'R', target = 'X', axes = NULL, tk = 2),  # reflected
    list(ord = 'X', target = 'X', axes = NULL),          # exact match
    list(ord = 'N', target = 'X', axes = NULL, tk = 2) ))
  for (symmetric in c(FALSE, TRUE))
    for (dilate in c(TRUE, FALSE)) {
      saved <- configs$X
      if (!is.null(cs$tk)) configs$X <- X[, seq_len(cs$tk), drop = FALSE]
      cases[[length(cases) + 1]] <- c(cs, list(symmetric = symmetric, dilate = dilate),
                                      fit(cs$ord, cs$target, symmetric, dilate))
      configs$X <- saved
    }

tests <- list()
for (cs in list(
    list(ord = 'Y', target = 'X', axes = 1:3),
    list(ord = 'N', target = 'X', axes = NULL, tk = 2),
    list(ord = 'X', target = 'Y', axes = 1:3) )) {
  tx <- configs[[cs$target]]
  if (!is.null(cs$tk)) tx <- tx[, seq_len(cs$tk), drop = FALSE]
  v <- suppressWarnings(protest(tx, configs[[cs$ord]], permutations = perms))
  tests[[length(tests) + 1]] <- c(cs, list(t0 = v$t0, ss = v$ss, signif = v$signif, t = unname(v$t)))
}


#--------------------------------------------------------
# Inputs for SciPy (same width: zero-padded as SciPy's
# docs say), and its disparities if already computed.
#--------------------------------------------------------
pad <- function (m, k) cbind(unname(m), matrix(0, nrow(m), k - ncol(m)))
pairs <- list(YX = list(pad(Y, 3), X), RX = list(R, X[, 1:2]), NX = list(N, X[, 1:2]))
for (nm in names(pairs)) {
  write.csv(format(pairs[[nm]][[1]], digits = 17), file.path(dir, paste0(nm, '_ord.csv')), row.names = FALSE)
  write.csv(format(pairs[[nm]][[2]], digits = 17), file.path(dir, paste0(nm, '_target.csv')), row.names = FALSE)
}

scipy_file <- file.path(dir, 'scipy_disparity.csv')
scipy <- if (file.exists(scipy_file)) {
  s <- read.csv(scipy_file, colClasses = c('character', 'numeric'))
  setNames(s$disparity, s$pair)
}

saveRDS(
  object = list(
    configs  = configs,
    perms    = perms,
    cases    = cases,
    tests    = tests,
    scipy    = scipy,
    versions = list(
      vegan   = as.character(packageVersion('vegan')),
      permute = as.character(packageVersion('permute')),
      scipy   = if (!is.null(scipy)) '1.18.1' ) ),
  file     = 'tests/tinytests/expected/procrustes.rds',
  version  = 2 )
