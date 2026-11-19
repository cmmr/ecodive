# Expected results for tests/tinytests/test-rpca_phylo.r, from gemelli
# 0.0.13's phylogenetic_rpca(). Run from the package root:
#
#     Rscript data-raw/expected/rpca_phylo.r
#
# Converts the text files in data-raw/expected/rpca_phylo/, written by
# data-raw/expected/rpca_phylo.py, into
# tests/tinytests/expected/rpca_phylo.rds: a list with one element per case.
# There is no Python where the package is developed, so the committed text
# files are the reference; this script only reads them.
#
# Each text file is a series of blocks: a line "@ name nrow ncol", then one
# line per row of space-separated numbers (hexadecimal floats, integers, or
# "NaN" for the robust clr matrix's missing values); or a line "$ name n",
# then n lines of strings. Each element of the list is a named list of those
# blocks (numeric matrices and character vectors), plus `name`:
#
#   newick, sample_ids, feature_ids, counts (samples x features, with those
#   dimnames), k, filters (min_sample_count, min_feature_count,
#   min_feature_frequency), min_depth:  the input.
#
#   samples (the kept samples' ids), fids (the node ids, in gemelli's
#   column order), parent (each column's parent column, 1-based, 0 for the
#   root), lengths (each column's branch length), node_counts (kept samples
#   x nodes), flat (1 for samples whose values are all equal), rclr (the
#   robust clr matrix, NaN where missing):  gemelli's node table.
#
#   u, v, eig, prop, dist:  phylogenetic_rpca() end to end: the sample
#   loadings, the feature (node) loadings, the singular values
#   ("eigvals"), proportion_explained, and the distances (in the order of a
#   dist object).

read_fixture <- function (path) {

  lines <- readLines(path)
  lines <- lines[!startsWith(lines, '#')]
  out   <- list()

  for (h in which(startsWith(lines, '@ ') | startsWith(lines, '$ '))) {
    p <- strsplit(lines[h], ' ', fixed = TRUE)[[1]]
    if (p[1] == '$') {
      out[[p[2]]] <- lines[h + seq_len(as.integer(p[3]))]
      next
    }
    nr <- as.integer(p[3])
    nc <- as.integer(p[4])
    v  <- as.numeric(unlist(strsplit(lines[h + seq_len(nr)], ' ', fixed = TRUE)))
    stopifnot(length(v) == nr * nc, !anyNA(v[!is.nan(v)]))
    out[[p[2]]] <- matrix(v, nr, nc, byrow = TRUE)
  }

  dimnames(out$counts) <- list(out$sample_ids, out$feature_ids)

  c(list(name = sub('\\.txt$', '', basename(path))), out)
}

names <- c('demo', 'gemelli', 'gemelli_filter', 'random', 'random_depth',
           'single_root', 'star', 'flat')

expected <- lapply(names, function (nm)
  read_fixture(file.path('data-raw', 'expected', 'rpca_phylo', paste0(nm, '.txt'))) )
names(expected) <- names

saveRDS(expected, 'tests/tinytests/expected/rpca_phylo.rds', version = 2)
