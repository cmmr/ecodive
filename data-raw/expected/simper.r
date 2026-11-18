# Expected results for tests/tinytests/test-simper.r, from vegan::simper()
# and summary() of it. vegan is GPL: only its outputs are used. Run from the
# package root:
#
#     Rscript data-raw/expected/simper.r
#
# Writes tests/tinytests/expected/simper.rds.
#
# ecodive's simper() permutes only each pair's samples, as its pairwise
# tests do, while vegan::simper() permutes all the samples' groups. So
# vegan is run on each pair's samples on their own, given the permutations
# that a full permutation matrix induces on them (written out here
# independently of ecodive, as a plain loop, as in pairwise.r), and the
# tests give ecodive the full matrix: p-values must then match exactly.

library(vegan)

data(dune)
data(dune.env)

set.seed(8)


# A sparse table with three unequal groups, a taxon absent everywhere, one
# with the same abundance in every sample, and a few that differ by group.
n      <- 25
sparse <- matrix(rpois(n * 40, 3) * rbinom(n * 40, 1, 0.3), n, 40)
dimnames(sparse) <- list(sprintf('S%02d', 1:n), sprintf('OTU%02d', 1:40))
sgrp <- sample(rep(c('x', 'y', 'z'), c(11, 8, 6)))
sparse[, 'OTU05'] <- 0L
sparse[, 'OTU06'] <- 2L
sparse[sgrp == 'z', 1:3] <- sparse[sgrp == 'z', 1:3] + 5L
sparse[sgrp == 'y', 4]   <- sparse[sgrp == 'y', 4] + 4L
storage.mode(sparse) <- 'integer'
stopifnot(all(rowSums(sparse) > 0))


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

plain <- function (pm) { attributes(pm) <- list(dim = dim(pm)); pm }


#------------------------------------------------------------------
# One case: vegan::simper() on each pair of groups' samples.
#   counts  The table ecodive is given (with `norm`).
#   vcounts The table vegan is given (counts, normalized).
#------------------------------------------------------------------
simper_case <- function (name, counts, groups, perms, norm = 'none',
                         vcounts = counts, blocks = NULL) {

  groups <- as.character(groups)
  lev    <- sort(unique(groups))
  combos <- combn(lev, 2)

  pairs <- lapply(seq_len(ncol(combos)), function (k) {

    a   <- combos[1, k]
    b   <- combos[2, k]
    idx <- which(groups %in% c(a, b))
    ip  <- induce_all(perms, idx)

    if (!is.null(blocks)) {
      bl <- blocks[idx]
      stopifnot(all(apply(ip, 1, function (r) all(bl[r] == bl))))
    }

    v  <- simper(vcounts[idx, , drop = FALSE], groups[idx], permutations = ip)
    stopifnot(length(v) == 1)
    vn <- strsplit(names(v), '_', fixed = TRUE)[[1]]  # vegan's (a, b) order
    sm <- summary(v)[[1]]
    sp <- rownames(sm)

    # vegan's ava is the mean of the first group in its own order.
    swap <- !identical(vn, c(a, b))

    data.frame(
      pair    = paste(a, b, sep = '-'),
      taxon   = sp,
      average = sm$average,
      sd      = sm$sd,
      ratio   = sm$ratio,
      mean_a  = if (swap) sm$avb else sm$ava,
      mean_b  = if (swap) sm$ava else sm$avb,
      cumsum  = sm$cumsum,
      p_value = unname(v[[1]]$p[sp]),
      overall = v[[1]]$overall,
      stringsAsFactors = FALSE )
  })

  list(
    name     = name,
    counts   = counts,
    groups   = groups,
    norm     = norm,
    perms    = perms,
    blocks   = blocks,
    expected = do.call(rbind, pairs) )
}


dune_m  <- as.matrix(dune)
storage.mode(dune_m) <- 'double'
mgmt    <- as.character(dune.env$Management)
use     <- as.character(dune.env$Use)

cases <- list(

  simper_case(
    'dune_free', dune_m, mgmt,
    plain(shuffleSet(nrow(dune_m), 999)) ),

  simper_case(
    'dune_blocks', dune_m, mgmt,
    plain(shuffleSet(nrow(dune_m), 199, control = how(blocks = factor(use)))),
    blocks = use ),

  simper_case(
    'sparse_free', sparse, sgrp,
    plain(shuffleSet(n, 499)) ),

  simper_case(
    'sparse_percent', sparse, sgrp,
    plain(shuffleSet(n, 199)),
    norm = 'percent', vcounts = decostand(sparse, 'total') ) )


expected <- list(
  cases        = cases,
  vegan        = as.character(packageVersion('vegan')),
  permute      = as.character(packageVersion('permute')) )

saveRDS(expected, 'tests/tinytests/expected/simper.rds', version = 2)
expected$vegan
