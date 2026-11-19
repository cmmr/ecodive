# Expected results for tests/tinytests/test-philr.r. Run from the package
# root, after data-raw/expected/philr.py:
#
#     Rscript data-raw/expected/philr.r
#
# Writes tests/tinytests/expected/philr.rds, so that the tests need no
# reference package.
#
# philr (Bioconductor, GPL-3) is output only: this script runs it, and its
# source was never read. Version 1.28.0 (Bioconductor 3.18, the newest that
# installs on R 4.3), with ape 5.x for the trees. Every weighting is run:
# part.weights 'uniform', 'gm.counts', 'anorm', 'anorm.x.gm.counts', 'enorm'
# and 'enorm.x.gm.counts', times ilr.weights 'uniform', 'blw', 'blw.sqrt'
# and 'mean.descendants'.
#
# scikit-bio 0.7.4's tree_basis() (BSD-3), from the text files that
# philr.py writes to data-raw/expected/philr/, is folded in as well.
#
# The result is a named list of cases, each a list of:
#
#   tree         the tree given to ecodive::philr(), as a plain phylo list
#                (no ape needed to use it);
#   counts       samples x features (features not in the tree's order);
#   pseudocount  the pseudocount (given to both);
#   ref_tree     the binary tree given to philr::philr(), when it differs
#                (resolved with ape::multi2di(random = FALSE), or pruned
#                with ape::keep.tip());
#   results      a list with one element per weighting run: part, ilr (the
#                names), coords (samples x coordinates, philr's columns),
#                p and ilr_weights (philr's return.all values);
#   dist         stats::dist() of the coordinates with the paper's weights
#                ('enorm.x.gm.counts', 'blw.sqrt'), where computed;
#   newick, tips, nodes, basis:  scikit-bio's tree_basis() (skbio cases).

suppressMessages({
  library(ape)
  library(philr)
})
stopifnot(packageVersion('philr') == '1.28.0')

PART <- c('uniform', 'gm.counts', 'anorm', 'anorm.x.gm.counts', 'enorm', 'enorm.x.gm.counts')
ILR  <- c('uniform', 'blw', 'blw.sqrt', 'mean.descendants')
ALL  <- expand.grid(part = PART, ilr = ILR, stringsAsFactors = FALSE)
PAPER <- data.frame(part = 'enorm.x.gm.counts', ilr = 'blw.sqrt')


plain_tree <- function (tr) {
  out <- list(
    edge        = unname(tr$edge),
    edge.length = tr$edge.length,
    tip.label   = tr$tip.label,
    Nnode       = tr$Nnode )
  if (!is.null(tr$node.label)) out$node.label <- tr$node.label
  structure(out, class = 'phylo')
}


run_case <- function (tree, counts, pc, runs = ALL, ref_tree = tree, dist = FALSE) {

  results <- lapply(seq_len(nrow(runs)), function (i) {
    r <- suppressWarnings(suppressMessages(philr::philr(
      counts, ref_tree, part.weights = runs$part[i], ilr.weights = runs$ilr[i],
      pseudocount = pc, return.all = TRUE )))
    list(
      part        = runs$part[i],
      ilr         = runs$ilr[i],
      coords      = unclass(r$x.ilrp),
      p           = r$p,
      ilr_weights = r$ilr.weights )
  })

  out <- list(tree = plain_tree(tree), counts = counts, pseudocount = pc, results = results)
  if (!identical(ref_tree, tree)) out$ref_tree <- plain_tree(ref_tree)

  if (dist) {
    r <- suppressWarnings(suppressMessages(philr::philr(
      counts, ref_tree, part.weights = 'enorm.x.gm.counts', ilr.weights = 'blw.sqrt',
      pseudocount = pc )))
    out$dist <- as.vector(stats::dist(r))
  }
  out
}


counts_for <- function (tips, n, zeros = 0, seed = 1) {
  set.seed(seed)
  D <- length(tips)
  m <- matrix(rpois(n * D, 30) * rbinom(n * D, 1, 1 - zeros), n, D)
  m[, 1] <- m[, 1] + 1
  dimnames(m) <- list(paste0('S', seq_len(n)), sample(tips))   # Not in tip order.
  m
}


expected <- list()


# Labelled internal nodes, zeros, zero-length tip and internal branches.
set.seed(20261119)
tr <- rtree(10)
tr <- makeNodeLabel(tr, method = 'number', prefix = 'node')
tr$edge.length[tr$edge[, 2] == 3]  <- 0
tr$edge.length[tr$edge[, 2] == 14] <- 0
x <- counts_for(tr$tip.label, 8, zeros = 0.15, seed = 1)
x <- x[, sample(ncol(x))]
expected$labelled <- run_case(tr, x, 0.5, dist = TRUE)


# Unlabelled, no zeros, no pseudocount; the smallest positive branch is
# internal, and two tip branches have length 0.
set.seed(20261120)
tr <- rtree(15)
tr$edge.length[which(tr$edge[, 2] > 15)[3]] <- 1e-3
tr$edge.length[which(tr$edge[, 2] <= 15)[c(2, 5)]] <- 0
x <- counts_for(tr$tip.label, 6, seed = 2)
expected$unlabelled <- run_case(tr, x, 0, dist = TRUE)


# Internal nodes numbered out of preorder: columns follow the numbers.
set.seed(20261121)
tr <- rtree(9)
perm <- c(1:9, 10L, 10L + sample(tr$Nnode - 1L))
tr$edge[] <- perm[tr$edge]
x <- counts_for(tr$tip.label, 5, seed = 3)
expected$renumbered <- run_case(tr, x, 1, runs = ALL[ALL$part %in% c('uniform', 'gm.counts'), ])


# Polytomies (an unrooted tree, with three children at its root, and two
# more) and tips that aren't features: philr is given the tree pruned with
# ape::keep.tip() and resolved with ape::multi2di(random = FALSE).
set.seed(20261122)
tr <- rtree(24)
tr <- di2multi(tr, tol = 0.12)
tr <- unroot(tr)
x  <- counts_for(sample(tr$tip.label, 20), 7, zeros = 0.1, seed = 4)
ref <- multi2di(keep.tip(tr, colnames(x)), random = FALSE)
stopifnot(is.binary(ref), is.rooted(ref), Ntip(ref) == 20)
expected$polytomy <- run_case(
  tr, x, 1, runs = rbind(ALL[ALL$part == 'uniform', ], PAPER, data.frame(part = 'anorm', ilr = 'mean.descendants')),
  ref_tree = ref, dist = TRUE )


# scikit-bio's tree_basis(), on trees read from Newick.
read_skbio <- function (name) {
  lines <- readLines(file.path('data-raw', 'expected', 'philr', paste0('skbio-', name, '.txt')))
  lines <- lines[!startsWith(lines, '#')]
  out   <- list()
  for (h in which(startsWith(lines, '@ ') | startsWith(lines, '$ '))) {
    p <- strsplit(lines[h], ' ', fixed = TRUE)[[1]]
    if (p[1] == '$') { out[[p[2]]] <- lines[h + seq_len(as.integer(p[3]))]; next }
    nr <- as.integer(p[3]); nc <- as.integer(p[4])
    v  <- as.numeric(unlist(strsplit(lines[h + seq_len(nr)], ' ', fixed = TRUE)))
    stopifnot(length(v) == nr * nc, !anyNA(v))
    out[[p[2]]] <- matrix(v, nr, nc, byrow = TRUE)
  }
  out
}

for (nm in c('small', 'random')) {
  sk <- read_skbio(nm)
  tr <- read.tree(text = sk$newick)
  x  <- counts_for(tr$tip.label, 6, zeros = 0.1, seed = 5)
  cs <- run_case(tr, x, 0.5, runs = rbind(ALL[ALL$part == 'uniform' & ALL$ilr == 'uniform', ], PAPER), dist = TRUE)
  expected[[paste0('skbio_', nm)]] <- c(cs, sk[c('newick', 'tips', 'nodes', 'basis')])
}


attr(expected, 'philr') <- as.character(packageVersion('philr'))
attr(expected, 'ape')   <- as.character(packageVersion('ape'))

saveRDS(expected, 'tests/tinytests/expected/philr.rds', version = 2)
