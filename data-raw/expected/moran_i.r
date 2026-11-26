# Writes tests/tinytests/expected/moran_i.rds: inputs and reference results
# for test-moran_i.r, so that the tests need no reference package.
#
# Run from the package root, after data-raw/expected/moran_i.py:
#
#     python data-raw/expected/moran_i.py
#     Rscript data-raw/expected/moran_i.r
#
# References:
#   * esda 2.10.0's Moran (BSD-3-Clause), from the hexadecimal floats that
#     moran_i.py writes to data-raw/expected/moran_i/: for each case and each
#     of esda's transformations "r", "b", "v" and "o", I, EI, VI_norm,
#     VI_rand, z_norm, z_rand, p_norm and p_rand (one-tailed, in the
#     direction of z), p_norm2 and p_rand2 (two-tailed), p_sim, EI_sim,
#     seI_sim and z_sim from the case's 199 permutations, and libpysal's s0,
#     s1 and s2; and the permuted statistics (`sims`, one column per
#     transformation), from the permutations in `perms` (1-based, one per
#     row).
#   * ape 5.7.1's Moran.I() (GPL-2 | GPL-3) and spdep 1.3.1's moran.test()
#     (GPL-2 | GPL-3), run as black boxes (outputs only; their code is not
#     read). Moran.I() with the weights as given (its outputs show that it
#     row-standardises them, and that its sd is the one under
#     randomisation); moran.test() with mat2listw() styles "W"
#     (row-standardised), "B" (of the binary matrix w != 0), "S"
#     (variance-stabilising) and "U" (scaled to sum 1, which gives the same I
#     and variance as the weights as given), under randomisation and
#     normality, with adjust.n = FALSE (so that an island counts in n, as in
#     esda) and zero.policy = TRUE.
#
# Each case holds x, w (the full weight matrix, zero diagonal), xy (the
# coordinates the weights came from, or NULL), perms, sims, esda (a data
# frame, one row per transformation), ape (list(observed, expected, sd,
# p_two_sided, p_greater)) and spdep (a data frame: style, randomisation,
# I, expectation, variance, p_greater).

library(ape)
library(spdep)

stopifnot(
  packageVersion('ape')   == '5.7.1',
  packageVersion('spdep') == '1.3.1' )

dir <- file.path('data-raw', 'expected', 'moran_i')

hex_matrix <- function (f) {
  m <- as.matrix(read.csv(file.path(dir, f), header = FALSE, colClasses = 'character'))
  matrix(as.numeric(m), nrow(m))
}

esda <- read.csv(file.path(dir, 'esda.csv'), colClasses = 'character')
for (col in setdiff(names(esda), c('case', 'style')))
  esda[[col]] <- as.numeric(esda[[col]])

cases <- list()

for (cs in unique(esda$case)) {

  x  <- hex_matrix(sprintf('x-%s.csv', cs))[, 1]
  w  <- hex_matrix(sprintf('w-%s.csv', cs))
  xy <- if (file.exists(file.path(dir, sprintf('xy-%s.csv', cs))))
          hex_matrix(sprintf('xy-%s.csv', cs))
  perms <- unname(as.matrix(read.csv(
    file.path(dir, sprintf('perms-%s.csv', cs)), header = FALSE )))
  storage.mode(perms) <- 'integer'

  sims <- read.csv(file.path(dir, sprintf('sims-%s.csv', cs)), colClasses = 'character')
  sims <- vapply(sims, as.numeric, numeric(nrow(sims)))

  a2 <- Moran.I(x, w)
  ag <- Moran.I(x, w, alternative = 'greater')

  sp <- do.call(rbind, lapply(c('W', 'B', 'S', 'U'), function (style) {
    # mat2listw()'s "B" keeps the matrix's own values, so it is given the
    # binary matrix.
    lw <- mat2listw(if (style == 'B') (w != 0) * 1 else w, style = style, zero.policy = TRUE)
    do.call(rbind, lapply(c(TRUE, FALSE), function (rnd) {
      mt <- moran.test(
        x, lw, randomisation = rnd, zero.policy = TRUE, adjust.n = FALSE,
        alternative = 'greater' )
      data.frame(
        style = style, randomisation = rnd,
        I           = unname(mt$estimate[[1]]),
        expectation = unname(mt$estimate[[2]]),
        variance    = unname(mt$estimate[[3]]),
        p_greater   = mt$p.value )
    }))
  }))

  cases[[cs]] <- list(
    x = x, w = w, xy = xy, perms = perms, sims = sims,
    esda  = esda[esda$case == cs, setdiff(names(esda), 'case')],
    ape   = list(
      observed = a2$observed, expected = a2$expected, sd = a2$sd,
      p_two_sided = a2$p.value, p_greater = ag$p.value ),
    spdep = sp )
}

saveRDS(cases, file.path('tests', 'tinytests', 'expected', 'moran_i.rds'), version = 2)
