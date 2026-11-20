# Writes the input tables and SpiecEasi's outputs for the reference
# comparisons of sparcc() into data-raw/expected/sparcc/. Then run
# data-raw/expected/sparcc_fastspar.sh (FastSpar's outputs) and
# data-raw/expected/sparcc.r (folds everything into the fixtures).
#
# Run from the package root, with SpiecEasi 1.1.1 installed (GPL-3; output
# only, installed from the jsDelivr mirror of github.com/zdk123/SpiecEasi):
#
#     Rscript data-raw/expected/sparcc_refs.r
#
# Tables are in FastSpar's BIOM tsv layout: features as rows, samples as
# columns, a "#OTU ID" header.

dir <- file.path('data-raw', 'expected', 'sparcc')
dir.create(dir, showWarnings = FALSE)

# Counts from log-normal "basis" abundances with three planted correlations
# (1 ~ 2 strongly positive, 3 ~ 4 strongly negative, 5 ~ 6 moderate).
make_table <- function (n, D, shift) {
  b <- matrix(rnorm(n * D), n, D)
  b[, 2] <- b[, 1] + rnorm(n, sd = 0.4)
  b[, 4] <- -b[, 3] + rnorm(n, sd = 0.4)
  b[, 6] <- 0.6 * b[, 5] + rnorm(n, sd = 0.8)
  matrix(
    rpois(n * D, exp(b + shift)), n, D,
    dimnames = list(sprintf('s%02d', 1:n), sprintf('otu%02d', 1:D)) )
}

write_biom <- function (x, file) {
  out <- cbind(`#OTU ID` = colnames(x), t(x))
  utils::write.table(out, file, sep = '\t', quote = FALSE, row.names = FALSE)
}

set.seed(20261121)
tables <- list(
  tight    = make_table(40, 20, 9),     # counts in the thousands
  moderate = make_table(60, 30, 2.5),   # counts around 10-30, some zeros
  d8       = make_table(40,  8, 8) )    # 8 features: past FastSpar's D - 3
for (nm in names(tables))
  write_biom(tables[[nm]], file.path(dir, paste0(nm, '.tsv')))


# SpiecEasi::sparcc(data, iter, inner_iter, th): samples as rows. Its draws
# use R's RNG, hence set.seed().
runs <- list(
  list(table = 'tight', iter = 100, inner_iter = 10),
  list(table = 'd8',    iter = 200, inner_iter = 10) )
for (r in runs) {
  set.seed(1)
  cor <- SpiecEasi::sparcc(tables[[r$table]], iter = r$iter, inner_iter = r$inner_iter, th = 0.1)$Cor
  dimnames(cor) <- list(colnames(tables[[r$table]]), colnames(tables[[r$table]]))
  out <- cbind(`#OTU ID` = rownames(cor), formatC(cor, digits = 17, format = 'g'))
  utils::write.table(
    out, file.path(dir, sprintf('spieceasi_%s_i%d_x%d.tsv', r$table, r$iter, r$inner_iter)),
    sep = '\t', quote = FALSE, row.names = FALSE )
}
