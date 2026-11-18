# Expected results for tests/tinytests/test-pairwise.r: pairwise = TRUE in
# permanova(), anosim() and mrpp(), from vegan's adonis2(), anosim() and
# mrpp() run on each pair of groups' samples (as pairwiseAdonis'
# pairwise.adonis2() and q2-diversity do). Run from the package root:
#
#     Rscript data-raw/expected/pairwise.r
#
# Writes tests/tinytests/expected/pairwise.rds.
#
# The tests give ecodive a permutation matrix of all the samples, and each
# pair is given the permutations it induces on the pair's samples (written
# out here independently of ecodive, as a plain loop), so p-values must
# match exactly, not just statistically. permdisp()'s pairwise t tests are
# checked against vegan::permutest() in data-raw/expected/permdisp.r.
#
# Also stores RVAideMemoire's pairwise.perm.manova() p-values (9999 random
# permutations, unadjusted and Holm-adjusted), to compare with up to
# sampling error. RVAideMemoire is GPL: only its outputs are used.

library(vegan)
library(RVAideMemoire)

set.seed(2)
n <- 34

counts <- matrix(rpois(n * 25, 4) * rbinom(n * 25, 1, 0.5), n, 25)
counts[, 1] <- counts[, 1] + 1L
dimnames(counts) <- list(sprintf('S%02d', 1:n), paste0('OTU', 1:25))
storage.mode(counts) <- 'integer'

# Four unequal groups, shuffled; 'd' differs on a few OTUs, 'c' a little.
meta <- data.frame(
  row.names = rownames(counts),
  grp   = sample(rep(c('a', 'b', 'c', 'd'), c(11, 9, 8, 6))),
  block = factor(rep(1:2, length.out = n)) )
counts[meta$grp == 'd', 2:4] <- counts[meta$grp == 'd', 2:4] + 6L
counts[meta$grp == 'c', 5]   <- counts[meta$grp == 'c', 5] + 3L

dists <- list(
  bray      = vegdist(counts, 'bray'),
  euclidean = vegdist(counts, 'euclidean') )

perms <- list(
  free   = shuffleSet(n, 199),
  blocks = shuffleSet(n, 199, control = how(blocks = meta$block)) )
attributes(perms$free)   <- list(dim = dim(perms$free))
attributes(perms$blocks) <- list(dim = dim(perms$blocks))


# The permutation of samples `idx` induced by permutation `p` of 1:n: each
# sample of idx goes to the first sample of idx after it on its cycle.
induce <- function (p, idx) {
  vapply(idx, function (i) {
    j <- p[[i]]
    while (!j %in% idx) j <- p[[j]]
    match(j, idx)
  }, integer(1))
}
induce_all <- function (pm, idx) t(apply(pm, 1, induce, idx = idx))

z_score <- function (stat, perm) (stat - mean(perm)) / sd(perm)


lev    <- sort(unique(meta$grp))
combos <- combn(lev, 2)
labels <- paste(combos[1, ], combos[2, ], sep = '-')

cases <- list()
for (d in names(dists)) for (pn in names(perms)) {

  dm  <- dists[[d]]
  dmm <- as.matrix(dm)
  out <- list(dist = d, perms = pn, induced = list())

  rows <- lapply(seq_len(ncol(combos)), function (k) {

    idx <- which(meta$grp %in% combos[, k])
    sub <- as.dist(dmm[idx, idx])
    g   <- meta$grp[idx]
    ip  <- induce_all(perms[[pn]], idx)

    # Induced permutations stay within the blocks.
    if (pn == 'blocks') {
      b <- meta$block[idx]
      stopifnot(all(apply(ip, 1, function (r) all(b[r] == b))))
    }

    a2 <- adonis2(sub ~ g, permutations = ip)
    fp <- attr(a2, 'F.perm')[, 1]
    an <- anosim(sub, g, permutations = ip)
    m1 <- mrpp(sub, g, permutations = ip, weight.type = 1)
    m3 <- mrpp(sub, g, permutations = ip, weight.type = 3)

    list(
      induced = ip,
      row = data.frame(
        n         = length(idx),
        F         = a2$F[[1]],
        R2        = a2$R2[[1]],
        F_z       = z_score(a2$F[[1]], fp),
        F_p       = a2$`Pr(>F)`[[1]],
        R         = an$statistic,
        R_z       = z_score(an$statistic, an$perm),
        R_p       = an$signif,
        delta1    = m1$delta,
        A1        = m1$A,
        delta1_z  = z_score(m1$delta, m1$boot.deltas),
        delta1_p  = m1$Pvalue,
        delta3    = m3$delta,
        A3        = m3$A,
        delta3_p  = m3$Pvalue ) )
  })

  out$table   <- do.call(rbind, lapply(rows, `[[`, 'row'))
  rownames(out$table) <- labels
  out$induced <- setNames(lapply(rows, `[[`, 'induced'), labels)
  # Keep the induced permutations of one case only: the tests check
  # ecodive's own against them.
  if (d != 'bray') out$induced <- NULL

  cases[[length(cases) + 1]] <- out
}

# RVAideMemoire: lower-triangle matrices of p-values, rows b, c, d and
# columns a, b, c, as pair "column-row".
rva <- lapply(c(none = 'none', holm = 'holm'), function (m) {
  set.seed(3)
  r <- pairwise.perm.manova(dists$bray, meta$grp, nperm = 9999, p.method = m)$p.value
  v <- r[lower.tri(r, diag = TRUE)]
  names(v) <- outer(rownames(r), colnames(r), function (a, b) paste(b, a, sep = '-'))[lower.tri(r, diag = TRUE)]
  v[labels]
})

expected <- list(
  counts = counts,
  meta   = meta,
  perms  = perms,
  cases  = cases,
  rva    = rva,
  rva_version = as.character(packageVersion('RVAideMemoire')),
  vegan  = as.character(packageVersion('vegan')) )

saveRDS(expected, 'tests/tinytests/expected/pairwise.rds', version = 2)
