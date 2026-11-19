# Expected results for tests/tinytests/test-dpcoa.r. Run from the package
# root:
#
#     Rscript data-raw/expected/dpcoa.r
#
# Writes tests/tinytests/expected/dpcoa.rds, so that the tests need no
# reference package.
#
# ade4 (GPL) and phyloseq (AGPL) are output only: this script runs them, and
# their source was never read. ade4 1.7-24's dpcoa(df, dis, nf), with dis
# the square roots of the patristic distances (ape::cophenetic.phylo())
# for a tree; phyloseq 1.46.0's ordinate(physeq, "DPCoA") (Bioconductor
# 3.18, the newest that installs on R 4.3), which calls ade4::dpcoa().
#
# The result is a named list of cases, each a list of:
#
#   counts     samples x taxa, as given to ecodive::dpcoa();
#   tree       the tree given to ecodive::dpcoa(), as a plain phylo list
#              (no ape needed to use it), or
#   taxa_dm    the taxa dissimilarities given instead (a dist object);
#   used       the taxa given to ade4 (it refuses taxa absent from every
#              sample, which ecodive keeps);
#   li, dls    ade4's sample and taxa coordinates (plain matrices);
#   eig        ade4's eigenvalues (all of them);
#   rao        ade4's RaoDiv;
#   between    ade4's between-samples diversity (RaoDecodiv[1, 1]);
#   phyloseq   (case 'globalpatterns' only) the same from phyloseq's
#              ordinate(): li, dls, eig, rao.

suppressMessages({
  library(ape)
  library(ade4)
})
stopifnot(packageVersion('ade4') == '1.7.24')

NF <- 4L

plain_tree <- function (tr) {
  out <- list(
    edge        = unname(tr$edge),
    edge.length = tr$edge.length,
    tip.label   = tr$tip.label,
    Nnode       = tr$Nnode )
  structure(out, class = 'phylo')
}

plain <- function (x) { x <- as.matrix(x); x }

run_ade4 <- function (counts, dis) {
  used <- colnames(counts)[colSums(counts) > 0]
  dis  <- as.dist(as.matrix(dis)[used, used])
  r    <- ade4::dpcoa(as.data.frame(counts[, used, drop = FALSE]), dis, scannf = FALSE, nf = NF)
  list(
    used    = used,
    li      = plain(r$li),
    dls     = plain(r$dls),
    eig     = r$eig,
    rao     = r$RaoDiv,
    between = r$RaoDecodiv[1, 1] )
}

tree_case <- function (counts, tree) {
  pd <- sqrt(cophenetic.phylo(tree))
  c(list(counts = counts, tree = plain_tree(tree)), run_ade4(counts, as.dist(pd)))
}

dm_case <- function (counts, taxa_dm) {
  c(list(counts = counts, taxa_dm = taxa_dm), run_ade4(counts, taxa_dm))
}

cases <- list()


#----------------------------------------------------------
# ecodive's example data (4 samples, 6 taxa, and the
# 6-tip ex_tree).
#----------------------------------------------------------
load('data/ex_counts.rda')
load('data/ex_tree.rda')
cases$ex <- tree_case(ex_counts, ex_tree)


#----------------------------------------------------------
# A random, non-ultrametric tree with extra tips (pruned by
# ecodive; ade4 gets the distances among the taxa), the
# taxa in another order than the tips, and a taxon absent
# from every sample.
#----------------------------------------------------------
set.seed(20261119)
tr <- rtree(60)
taxa <- sample(tr$tip.label, 45)
X <- matrix(rpois(18 * 45, 3) * rbinom(18 * 45, 1, 0.4), 18, 45,
  dimnames = list(paste0('S', 1:18), taxa))
X[, 1] <- X[, 1] + 1L   # no empty samples
X[, 7] <- 0L            # absent from every sample
cases$rtree <- tree_case(X, tr)


#----------------------------------------------------------
# An ultrametric tree, made unrooted (3 children at the
# root) and with a polytomy, with proportions rather than
# counts and one sample holding a single taxon.
#----------------------------------------------------------
set.seed(7)
tr <- rcoal(25)
tr$edge.length[tr$edge[, 2] == 30] <- 0
tr <- di2multi(tr)        # a polytomy where that branch was
tr <- unroot(tr)
X <- matrix(runif(9 * 25) * rbinom(9 * 25, 1, 0.6), 9, 25,
  dimnames = list(paste0('S', 1:9), sample(tr$tip.label)))
X[, 2] <- X[, 2] + 0.5
X[9, ] <- 0; X[9, 3] <- 2.5
cases$unrooted <- tree_case(X, tr)
stopifnot(!is.rooted(tr), !is.binary(tr))


#----------------------------------------------------------
# ade4's own examples: humDNAm (taxa_dm, the square roots
# of its distances, as its help page does) and ecomor
# (dist.taxo() of a taxonomy).
#----------------------------------------------------------
data(humDNAm)
cases$humDNAm <- dm_case(t(as.matrix(humDNAm$samples)), sqrt(humDNAm$distances))

data(ecomor)
cases$ecomor <- dm_case(t(as.matrix(ecomor$habitat)), dist.taxo(ecomor$taxo))


#----------------------------------------------------------
# phyloseq's GlobalPatterns, cut to its 150 most abundant
# taxa: ade4, and phyloseq's ordinate(method = 'DPCoA').
#----------------------------------------------------------
suppressMessages(library(phyloseq))
stopifnot(packageVersion('phyloseq') == '1.46.0')
data(GlobalPatterns)
top <- names(sort(taxa_sums(GlobalPatterns), decreasing = TRUE))[1:150]
gp  <- prune_taxa(top, GlobalPatterns)
X   <- t(as(otu_table(gp), 'matrix'))
storage.mode(X) <- 'integer'
cases$globalpatterns <- tree_case(X, phy_tree(gp))

o <- ordinate(gp, method = 'DPCoA', scannf = FALSE, nf = NF)
cases$globalpatterns$phyloseq <- list(
  li  = plain(o$li),
  dls = plain(o$dls),
  eig = o$eig,
  rao = o$RaoDiv )


saveRDS(cases, 'tests/tinytests/expected/dpcoa.rds', version = 2)
cat('Wrote tests/tinytests/expected/dpcoa.rds:',
  format(file.size('tests/tinytests/expected/dpcoa.rds'), big.mark = ','), 'bytes\n')
