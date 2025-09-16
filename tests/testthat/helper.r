
counts <- structure(
  .Data = c(0,0,0,0,0,8,9,10,5,5,5,5,2,0,0,0,6,5,7,0), 
  dim = 4:5, dimnames = list(LETTERS[1:4], paste0('OTU', 1:5)) )

big_mtx <- structure(
  .Data = rep(c(0,0,0,0,0,8,9,10,5,5,5,5,2,0,0,0,6,5,7,0), 26), 
  dim = c(104, 5), dimnames = list(paste0('S', 1:104), paste0('OTU', 1:5)) )

nz_counts <- counts + 1

counts_p    <- t(apply(counts,    1L, function (x) x / sum(x)))
nz_counts_p <- t(apply(nz_counts, 1L, function (x) x / sum(x)))

tree_str <- "((OTU4:0.676,(('OTU2':0.548,OTU3:0.629):0.25,OTU1:0.751):0.276):0.056,OTU5:0.432);"
tree     <- read_tree(tree_str)
