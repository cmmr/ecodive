# Writes tests/tinytests/expected/coda_diagnostics.rds: inputs and
# reference results for test-coda_diagnostics.r, so that the tests need no
# reference package.
#
# Run from the package root:
#
#     Rscript data-raw/expected/coda_diagnostics.r
#
# References, run as black boxes (outputs only; their code is not read):
#   * compositions 2.0.9's variation() (GPL-2 | GPL-3): the variation matrix
#     tau, full and symmetric, of acomp(x) with robust = FALSE; and mvar(),
#     the total variance.
#   * zCompositions 1.6.2's zPatterns() (GPL-2 | GPL-3): each sample's
#     pattern id, with label = 0. Its ids are numbered its own way; the
#     tests compare the partition of the samples and the patterns' sizes.
#
# Cases (no zeros for compositions, which reads zeros as below a detection
# limit):
#   va_pos   a 15 x 7 table of positive values (Gamma draws), named.
#   va_pc    ex_counts + 0.5 (ecodive's example table, with a pseudocount
#            added by hand).
#   zp_*     count tables with zeros: a hand-built 6 x 4 table, 40 x 8
#            Poisson counts, and 100 x 15 counts with a range of means.

library(compositions)
library(zCompositions)

stopifnot(
  packageVersion('compositions')  == '2.0.9',
  packageVersion('zCompositions') == '1.6.2' )

ex_counts <- local({ e <- new.env(); load('data/ex_counts.rda', e); e$ex_counts })

set.seed(20261127)
va_pos <- matrix(
  rgamma(15 * 7, shape = rep(c(0.8, 2, 5, 1, 3, 10, 0.5), each = 15)), 15, 7,
  dimnames = list(sprintf('s%02d', 1:15), paste0('f', 1:7)) )
va_pc  <- ex_counts + 0.5

va_ref <- function (x) list(
  x         = x,
  variation = unclass(variation(acomp(x), robust = FALSE)),
  mvar      = mvar(acomp(x), robust = FALSE) )

zp_hand <- rbind(
  c(1, 0, 2, 3), c(0, 0, 1, 1), c(1, 0, 2, 3),
  c(5, 1, 0, 1), c(0, 0, 1, 1), c(1, 0, 4, 9) )
zp_pois <- matrix(rpois(40 * 8, 1.2), 40, 8)
zp_wide <- matrix(rpois(100 * 15, rep(rexp(15, 0.5), each = 100)), 100, 15)
zp_wide <- zp_wide[rowSums(zp_wide > 0) > 0, ]   # zPatterns() needs a non-zero per row

zp_ref <- function (x) list(
  x  = x,
  id = as.integer(as.character(zPatterns(x, label = 0, plot = FALSE, show.means = FALSE))) )

expected <- list(
  va = list(pos = va_ref(va_pos), pc = va_ref(va_pc)),
  zp = list(hand = zp_ref(zp_hand), pois = zp_ref(zp_pois), wide = zp_ref(zp_wide)) )

saveRDS(expected, 'tests/tinytests/expected/coda_diagnostics.rds', version = 2)
