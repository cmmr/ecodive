# Expected results for tests/tinytests/test-phylo_structure.r. Run from the
# package root:
#
#     Rscript data-raw/expected/phylo_structure.r
#
# Writes tests/tinytests/expected/phylo_structure.rds, so that the tests
# need no reference package.
#
# picante 1.8.4 (GPL-2) is output only: this script runs it, and its source
# was never read. ape 5.8.1's cophenetic.phylo() gives picante its `dis`.
#
# picante's ses.mpd() and ses.mntd() draw their random communities with
# R's RNG, through its exported taxaShuffle() (taxa.labels) and
# randomizeMatrix() (richness, frequency, independentswap, trialswap). It
# was checked, from outputs only, that set.seed(s); ses.mpd(...) gives the
# same null means as set.seed(s) followed by `runs` calls of
# mpd(randomizeMatrix(samp, model), dis) (or mpd(samp, taxaShuffle(dis))),
# so the random communities of a seeded ses.mpd() call can be rebuilt here
# and handed to ecodive's internal hook (`null_tables`): ecodive's SES must
# then equal picante's. picante's sample.pool gives the same communities as
# its richness (checked too); phylogeny.pool can't be rebuilt, so it is
# compared statistically only.
#
# The result is a list:
#
#   tree, samp, dis   picante's `phylocom` example: its tree (32 tips, 25 of
#                     them in the table, as a plain phylo list), table
#                     (6 x 25) and cophenetic(tree) (32 x 32);
#   obs               picante's mpd() / mntd() of samp, unweighted and
#                     weighted (mpd, wmpd, mntd, wmntd);
#   big               a 300-tip rtree(), a sparse 40 x 250 table of its tips
#                     (an empty sample and a one-taxon sample), and picante's
#                     four values for it;
#   traits            a trait-distance taxa_dm for samp's taxa, and
#                     picante's four values with it;
#   hook              per model: the seed, the random tables of a 20-run
#                     ses call (each 6 x 32 for taxa.labels, else 6 x 25),
#                     and picante's ses.mpd() / ses.mntd() results,
#                     unweighted and weighted, with that seed;
#   stats             per model: picante's null means and sds over 9,999
#                     runs (statistical comparisons only).

suppressMessages({
  library(ape)
  library(picante)
})
stopifnot(packageVersion('picante') == '1.8.4')

plain_tree <- function (tr) {
  out <- list(
    edge        = unname(tr$edge),
    edge.length = tr$edge.length,
    tip.label   = tr$tip.label,
    Nnode       = tr$Nnode )
  structure(out, class = 'phylo')
}

four <- function (samp, dis) list(
  mpd   = unname(picante::mpd(samp, dis)),
  wmpd  = unname(picante::mpd(samp, dis, TRUE)),
  mntd  = unname(picante::mntd(samp, dis)),
  wmntd = unname(picante::mntd(samp, dis, TRUE)) )

plain_ses <- function (x) { x <- as.data.frame(x); rownames(x) <- NULL; x }


data(phylocom)
tree <- phylocom$phylo
samp <- phylocom$sample
dis  <- cophenetic(tree)
out  <- list(tree = plain_tree(tree), samp = samp, dis = dis, obs = four(samp, dis))


# A larger, sparse case.
set.seed(1)
tr <- rtree(300)
m  <- matrix(
  rpois(40 * 250, 0.3) * rbinom(40 * 250, 1, 0.3), 40, 250,
  dimnames = list(paste0('s', 1:40), sample(tr$tip.label, 250)) )
m[1, ] <- 0; m[1, 5] <- 2   # one taxon
m[2, ] <- 0                 # empty
storage.mode(m) <- 'double'
out$big <- list(tree = plain_tree(tr), counts = m, obs = four(m, cophenetic(tr)))


# Trait distances.
set.seed(2)
traits <- matrix(rnorm(25 * 3), 25, 3, dimnames = list(colnames(samp), NULL))
tdm    <- dist(traits)
out$traits <- list(taxa_dm = tdm, obs = four(samp, as.matrix(tdm)))


# Shared randomizations.
R <- 20L
rebuild <- function (model, seed) {
  set.seed(seed)
  lapply(seq_len(R), function (r) {
    if (model == 'taxa.labels') {
      d2  <- taxaShuffle(dis)
      tip <- rownames(dis)[match(colnames(samp), rownames(d2))]
      x   <- matrix(0, nrow(samp), nrow(dis), dimnames = list(rownames(samp), rownames(dis)))
      x[, tip] <- samp
      x
    }
    else {
      rm <- if (model == 'sample.pool') 'richness' else model
      randomizeMatrix(samp, rm, iterations = 100)
    }
  })
}

hook <- list()
for (model in c('taxa.labels', 'richness', 'frequency', 'sample.pool', 'independentswap', 'trialswap')) {
  seed   <- 100L + length(hook)
  tables <- rebuild(model, seed)
  ses    <- list()
  for (w in c(FALSE, TRUE)) {
    set.seed(seed)
    a <- ses.mpd(samp, dis, null.model = model, abundance.weighted = w, runs = R, iterations = 100)
    set.seed(seed)
    b <- ses.mntd(samp, dis, null.model = model, abundance.weighted = w, runs = R, iterations = 100)
    # The rebuilt tables give picante's own null means.
    v <- sapply(tables, function (x) picante::mpd(x, dis, w))
    stopifnot(isTRUE(all.equal(unname(rowMeans(v, na.rm = TRUE)), unname(a$mpd.rand.mean))))
    v <- sapply(tables, function (x) picante::mntd(x, dis, w))
    stopifnot(isTRUE(all.equal(unname(rowMeans(v, na.rm = TRUE)), unname(b$mntd.rand.mean))))
    ses[[if (w) 'weighted' else 'unweighted']] <- list(mpd = plain_ses(a), mntd = plain_ses(b))
  }
  hook[[model]] <- list(seed = seed, tables = tables, ses = ses)
}
out$hook <- hook


# Null summaries over 9,999 runs.
stats <- list()
for (model in c('taxa.labels', 'phylogeny.pool', 'richness', 'frequency', 'independentswap', 'trialswap')) {
  set.seed(7)
  a <- ses.mpd(samp, dis, null.model = model, runs = 9999)
  set.seed(8)
  b <- ses.mpd(samp, dis, null.model = model, abundance.weighted = TRUE, runs = 9999)
  set.seed(9)
  d <- ses.mntd(samp, dis, null.model = model, runs = 9999)
  stats[[model]] <- list(
    mpd   = plain_ses(a)[, c('mpd.rand.mean', 'mpd.rand.sd')],
    wmpd  = plain_ses(b)[, c('mpd.rand.mean', 'mpd.rand.sd')],
    mntd  = plain_ses(d)[, c('mntd.rand.mean', 'mntd.rand.sd')] )
}
out$stats <- stats

saveRDS(out, 'tests/tinytests/expected/phylo_structure.rds', version = 2)
