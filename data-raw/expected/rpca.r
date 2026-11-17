# Expected results for tests/tinytests/test-rpca.r, from gemelli 0.0.13's
# rpca(). Run from the package root:
#
#     Rscript data-raw/expected/rpca.r
#
# Converts the text files in data-raw/expected/rpca/, written by
# data-raw/expected/rpca.py, into tests/tinytests/expected/rpca.rds: a list
# with one element per table. There is no Python where the package is
# developed, so the committed text files are the reference; this script only
# reads them.
#
# Each text file is a series of blocks: a line "@ name nrow ncol", then one
# line per row of space-separated values ("NaN" for the rclr matrix's
# missing values). Each element of the list is a named list of those blocks
# as numeric matrices, with the counts given sample and feature names (S1,
# S2, ... and F1, F2, ..., as gemelli was given them), plus `name`:
#
#   counts, k, filters (min_sample_count, min_feature_count,
#   min_feature_frequency):  the input.
#
#   e2e_samples, e2e_features, e2e_u, e2e_v, e2e_eig, e2e_prop, e2e_dist:
#     gemelli.rpca.rpca() end to end: the kept samples and features
#     (1-based), the sample loadings u, the feature loadings v, the singular
#     values ("eigvals"), proportion_explained, and the distances (in the
#     order of a dist object).
#
#   For some tables, OptSpace's stages, on the unfiltered table: rclr (n x
#   p, NaN where missing), u0 and v0 (the initial factors after OptSpace's
#   scaling), cost (the fit error after initialization and after each
#   iteration), U, S, V (OptSpace's final factors), u, v, eig, prop (the
#   final PCA).
#
#   For "groups", also pen_step, pen_U, pen_S, pen_V, pen_cost: OptSpace
#   with step_size = pen_step, so that its Grassmann penalty is active.

read_fixture <- function (path) {

  lines <- readLines(path)
  lines <- lines[!startsWith(lines, '#')]
  out   <- list()

  for (h in which(startsWith(lines, '@ '))) {
    p  <- strsplit(lines[h], ' ', fixed = TRUE)[[1]]
    nr <- as.integer(p[3])
    nc <- as.integer(p[4])
    v  <- as.numeric(unlist(strsplit(lines[h + seq_len(nr)], ' ', fixed = TRUE)))
    stopifnot(length(v) == nr * nc)
    out[[p[2]]] <- matrix(v, nr, nc, byrow = TRUE)
  }

  dimnames(out$counts) <- list(
    paste0('S', seq_len(nrow(out$counts))),
    paste0('F', seq_len(ncol(out$counts))) )

  c(list(name = sub('\\.txt$', '', basename(path))), out)
}

names <- c('groups', 'sparse', 'filter', 'tall', 'wide', 'khigh', 'nozero')

expected <- lapply(names, function (nm)
  read_fixture(file.path('data-raw', 'expected', 'rpca', paste0(nm, '.txt'))) )
names(expected) <- names

saveRDS(expected, 'tests/tinytests/expected/rpca.rds', version = 2)
