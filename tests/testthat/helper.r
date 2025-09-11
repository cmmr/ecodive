
counts <- matrix(
  as.integer(c(0, 0, 5, 2, 6, 0, 8, 5, 0, 5, 0, 9, 5, 0, 7, 0, 10, 5, 0, 0)),
  nrow     = 5,
  dimnames = list(paste0('OTU', 1:5), LETTERS[1:4]) )

nz_counts <- counts + 1

counts_p <- apply(counts, 2L, function (x) x / sum(x))
nz_counts_p <- apply(nz_counts, 2L, function (x) x / sum(x))

tree_str <- "((OTU4:0.676,(('OTU2':0.548,OTU3:0.629):0.25,OTU1:0.751):0.276):0.056,OTU5:0.432);"
tree     <- read_tree(tree_str)

big_mtx <- matrix(
    as.integer(rep(c(0, 0, 5, 2, 6, 0, 8, 5, 0, 5, 0, 9, 5, 0, 7, 0, 10, 5, 0, 0), 26)),
    nrow     = 5,
    dimnames = list(paste0('OTU', 1:5), paste0('S', 1:104) ))
