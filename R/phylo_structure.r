# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' Mean Pairwise and Mean Nearest Taxon Distance
#'
#' How closely related the taxa in each sample are, along a phylogeny (or
#' any dissimilarity among the taxa): the mean phylogenetic distance between
#' pairs of the sample's taxa (MPD), or from each taxon to its nearest
#' relative in the sample (MNTD) (Webb 2000; Webb et al. 2002). MPD responds
#' to structure deep in the tree, MNTD to structure near its tips. The
#' `weighted_*` forms weight each taxon by its abundance. To test them
#' against a null model (NRI and NTI), see [mpd_ses()].
#'
#' @inheritParams documentation
#'
#' @param counts   A numeric matrix of non-negative abundances (samples
#'        \eqn{\times} taxa, see `margin`), or any of the other objects
#'        listed in "Input Types". Values need not be integers: `mpd()` and
#'        `mntd()` use only presence, and the weighted forms only relative
#'        abundances. The taxa must have names, to match them to the tips of
#'        `tree` or the labels of `taxa_dm`.
#'
#' @param tree   A `phylo` object with branch lengths (none negative or
#'        missing), whose tips include every taxon of `counts` (see
#'        [read_tree()]). Tips that aren't taxa of `counts` are ignored.
#'        Can be omitted if a tree is embedded with the `counts` object or
#'        as `attr(counts, 'tree')`. Distances are path lengths along the
#'        tree (patristic distances), as from `ape::cophenetic.phylo()`.
#'        The tree may be unrooted, or have polytomies or nodes with one
#'        child: none of these change path lengths.
#'
#' @param taxa_dm   Instead of `tree`: the dissimilarities among the taxa,
#'        as a `dist` object or symmetric matrix (picante's `dis`, which may
#'        be functional-trait distances, for instance). Labels are matched to
#'        the taxa of `counts`; without labels, rows are taken in the order
#'        of `counts`' taxa. An embedded tree is ignored when `taxa_dm` is
#'        given.
#'
#' @return A numeric vector with one value per sample, named by sample.
#'   A sample with fewer than two taxa is `NA`.
#'
#' @details
#'
#' For a sample with taxa \eqn{i = 1 \ldots S}, abundances \eqn{x_i}
#' (total \eqn{W}), and distances \eqn{d_{ij}}:
#'
#' * **MPD:** \eqn{\sum_{i<j} d_{ij} / (S(S-1)/2)}, the mean over
#'   pairs of distinct taxa.
#' * **Weighted MPD:** \eqn{\sum_{i,j} x_i x_j d_{ij} / W^2}, the
#'   mean distance between two individuals drawn at random, with
#'   replacement, so that pairs within a taxon count too (at distance 0).
#'   This is picante's and phylocom's definition.
#' * **MNTD:** \eqn{\sum_i m_i / S}, with \eqn{m_i = \min_{j \ne i}
#'   d_{ij}}.
#' * **Weighted MNTD:** \eqn{\sum_i x_i m_i / W}.
#'
#' The weighted forms are unchanged by scaling a sample, so counts and
#' proportions give the same values. "Weighted" here means weighting by
#' abundance; it is not the weighting of [weighted_unifrac()].
#'
#' These are also alpha diversity metrics: `alpha_div(counts, 'mpd', tree =
#' tree)`, and likewise `'weighted_mpd'`, `'mntd'` and `'weighted_mntd'`.
#'
#'
#' **Speed and memory**
#'
#' With a tree, no taxa \eqn{\times} taxa matrix is formed. Each sample's
#' taxa are sorted into the tree's preorder, and the tree's lowest common
#' ancestors of neighbours (found in constant time from a table built once)
#' give the part of the tree that spans them, with at most \eqn{2S - 1}
#' nodes; MPD and MNTD are one or two passes over it. Each sample takes
#' time proportional to \eqn{S \log S}, however large the tree. With
#' `taxa_dm`, each sample reads its \eqn{S^2} distances. Results are
#' identical for any number of `cpus`.
#'
#'
#' **Differences from picante**
#'
#' Values agree with picante 1.8.4's `mpd()` and `mntd()` (with `dis =
#' cophenetic(tree)` and `abundance.weighted` set to match) to about
#' \eqn{10^{-12}}. picante computes `cophenetic(tree)` first, a dense
#' matrix of every pair of tips (3.2 GB for 20,000 tips); `mpd()` reads the
#' tree directly. phylocom's `comstruct` gives MNTD 0, not NA, for a sample
#' with one taxon.
#'
#'
#' **Same names in other packages**
#'
#' `picante::mpd()` and `picante::mntd()` share these names (taking `samp`,
#' `dis` and `abundance.weighted`); whichever package was attached last
#' masks the other's. Write `ecodive::mpd()` to be sure.
#' `PhyloMeasures::mpd.query()` computes the unweighted MPD.
#'
#' @references
#' Webb, C. O. (2000). Exploring the phylogenetic structure of ecological
#' communities: an example for rain forest trees. *American Naturalist*,
#' 156(2), 145-155. \doi{10.1086/303378}
#'
#' Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
#' Phylogenies and community ecology. *Annual Review of Ecology and
#' Systematics*, 33, 475-505.
#' \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H.,
#' Ackerly, D. D., Blomberg, S. P., & Webb, C. O. (2010). Picante: R tools
#' for integrating phylogenies and ecology. *Bioinformatics*, 26(11),
#' 1463-1464. \doi{10.1093/bioinformatics/btq166}
#'
#' @inheritSection documentation Input Types
#'
#' @family Phylogenetic metrics
#' @seealso [mpd_ses()] for NRI and NTI; [faith()]; [alpha_div()];
#'   [read_tree()].
#'
#' @export
#' @examples
#'     mpd(ex_counts, ex_tree)
#'     weighted_mpd(ex_counts, ex_tree)
#'     mntd(ex_counts, ex_tree)
#'     weighted_mntd(ex_counts, ex_tree)
#'
#'     # The same through alpha_div().
#'     alpha_div(ex_counts, 'weighted_mntd', tree = ex_tree)
#'
#'     # Distances between (made-up) trait values instead of a tree.
#'     traits <- data.frame(
#'       row.names = colnames(ex_counts),
#'       gram_pos  = c(1, 0, 1, 0, 1, 1),
#'       size      = c(0.8, 1.5, 0.7, 1.2, 0.9, 1.0) )
#'     mpd(ex_counts, taxa_dm = dist(traits))
#'
mpd <- function (counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus()) {
  phylo_observed(counts, tree, taxa_dm, 0L, FALSE, margin, cpus)
}

#' @rdname mpd
#' @export
weighted_mpd <- function (counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus()) {
  phylo_observed(counts, tree, taxa_dm, 0L, TRUE, margin, cpus)
}

#' @rdname mpd
#' @export
mntd <- function (counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus()) {
  phylo_observed(counts, tree, taxa_dm, 1L, FALSE, margin, cpus)
}

#' @rdname mpd
#' @export
weighted_mntd <- function (counts, tree = NULL, taxa_dm = NULL, margin = 1L, cpus = n_cpus()) {
  phylo_observed(counts, tree, taxa_dm, 1L, TRUE, margin, cpus)
}



#' Net Relatedness and Nearest Taxon Indices (SES of MPD and MNTD)
#'
#' Tests whether the taxa in each sample are more closely related
#' (phylogenetic clustering) or more distantly related (overdispersion)
#' than in random communities: the standardized effect sizes of [mpd()]
#' and [mntd()] against a null model. Their negatives are Webb's Net
#' Relatedness Index (**NRI**, from MPD) and Nearest Taxon Index (**NTI**,
#' from MNTD), so positive values mean clustering. This is picante's
#' `ses.mpd()` and `ses.mntd()`, and phylocom's `comstruct`.
#'
#' @inheritParams mpd
#' @inheritParams documentation
#'
#' @param weighted   Weight each taxon by its abundance: [weighted_mpd()]
#'        and [weighted_mntd()] (picante's `abundance.weighted`, phylocom's
#'        `-a`). Default: `FALSE`
#'
#' @param null_model   How the random communities are made. One of
#'        `"taxa_labels"` (the default), `"richness"`, `"frequency"`,
#'        `"sample_pool"`, `"phylogeny_pool"`, `"independent_swap"`,
#'        `"trial_swap"`, or (for `mpd_ses()` only) `"exact"`; see "Null
#'        models". picante's spellings (`"taxa.labels"`, `"independentswap"`,
#'        ...) are accepted too.
#'
#' @param n_draws   The number of random communities per sample (picante's
#'        `runs`). Default: `999`
#'
#' @param n_swaps   For the two swap models: the number of swaps that make
#'        each random table from the observed one (picante's `iterations`).
#'        Ignored otherwise. Default: `1000`
#'
#' @param seed   An integer seed for the random communities. The same `seed`
#'        gives identical results, for any number of `cpus`. Default: `0`
#'
#' @return A data frame with one row per sample, of class
#'   `"ecodive_mpd_ses"` (`"ecodive_mntd_ses"`), with columns:
#'
#'   * `sample`: The sample's name.
#'   * `n_taxa`: Its number of taxa (non-zero values).
#'   * `obs`: Its MPD (MNTD).
#'   * `null_mean`, `null_sd`: The mean and standard deviation (divisor
#'     \eqn{n - 1}) of MPD (MNTD) over the random communities.
#'   * `z`: The standardized effect size, `(obs - null_mean) / null_sd`
#'     (picante's `mpd.obs.z`).
#'   * `nri` (`nti`): `-z`.
#'   * `rank`: The rank of `obs` among itself and the random values, ties
#'     averaged.
#'   * `p_value`: `rank / (n_draws + 1)`: the lower tail, so a small
#'     `p_value` means clustering, and one near 1 overdispersion (picante's
#'     `mpd.obs.p`).
#'   * `n_draws`: The random communities used: those where MPD (MNTD) is
#'     defined (with `"frequency"`, a random sample can have fewer than two
#'     taxa).
#'
#'   With `null_model = "exact"`, `null_mean` and `null_sd` are exact, and
#'   `rank`, `p_value` and `n_draws` are `NA`.
#'
#'   Attributes: `method` (`"mpd_ses"` or `"mntd_ses"`), `null_model`,
#'   `n_draws` (as given), `n_swaps` (the swap models only) and `weighted`.
#'
#' @details
#'
#' **Null models**
#'
#' | `null_model` | Each random community | Keeps |
#' |:-------------|:----------------------|:------|
#' | `"taxa_labels"` | the tips' labels shuffled: the table's taxa are moved to tips drawn at random from **every** tip of `tree` (or label of `taxa_dm`) | the table |
#' | `"richness"` | each sample's values moved to taxa drawn at random from all the table's taxa | each sample's richness and values |
#' | `"frequency"` | each taxon's values moved to samples drawn at random | each taxon's frequency and values |
#' | `"sample_pool"` | as `"richness"`, from the taxa present in at least one sample | richness |
#' | `"phylogeny_pool"` | as `"richness"`, from every tip of `tree` (or label of `taxa_dm`) | richness |
#' | `"independent_swap"` | `n_swaps` swaps of 2 \eqn{\times} 2 checkerboards (Gotelli 2000) | each sample's richness and each taxon's frequency |
#' | `"trial_swap"` | `n_swaps` attempted swaps (Miklós & Podani 2004) | as above |
#' | `"exact"` | (no draws) the exact mean and standard deviation of MPD over random sets of tips, as `"phylogeny_pool"` | richness |
#'
#' A sample's abundances go with its taxa: under the pool models and
#' `"richness"` they are given in random order to the drawn taxa, and in
#' the swaps each taxon keeps its values (as in picante and phylocom). They
#' matter only with `weighted = TRUE`.
#'
#' Which null model answers which question is discussed by Kembel (2009)
#' and Hardy (2008): `"taxa_labels"` and `"phylogeny_pool"` ask whether the
#' sample is a random draw from the tree, `"richness"` and `"sample_pool"`
#' whether it is a random draw from the taxa observed, and the swaps also
#' hold how common each taxon is. Because `"taxa_labels"` and
#' `"phylogeny_pool"` draw from **every tip**, a tree with many taxa that
#' could never occur in the study makes every sample look clustered: prune
#' the tree to the study's species pool first (picante gives the same
#' advice).
#'
#' **Exact moments** (`null_model = "exact"`, `mpd_ses()` only). For
#' \eqn{S} tips drawn at random from all \eqn{T} without replacement (the
#' `"phylogeny_pool"` model; per sample, also `"taxa_labels"`), with the
#' sample's abundances given to them at random, MPD is a sum over pairs of
#' positions. Its mean is the mean distance \eqn{\mu} between two distinct
#' tips, and its variance follows from three sums over the tips: of the
#' squared deviations \eqn{(d_{ab} - \mu)^2}, and of the squares of each
#' tip's total deviation \eqn{\sum_b (d_{ab} - \mu)}, weighted by how many
#' pairs of sample positions share two, one or no positions. With a tree,
#' all three take one pass over its branches, so no draws are needed (as
#' in Tsirogiannis, Sandel & Cheliotis 2012, for the unweighted MPD). Then
#' `z` is exact, but has no rank or p-value; `pnorm(z)` is a normal
#' approximation to the lower-tail p-value. No exact form is implemented
#' for MNTD.
#'
#'
#' **Speed and memory**
#'
#' Each random community's MPD (MNTD) is computed from the tree as in
#' [mpd()], in time proportional to \eqn{S \log S}, with no taxa
#' \eqn{\times} taxa matrix. Draws run in parallel, each from its own
#' random-number stream, and are summed in a fixed order, so results are
#' identical for any number of `cpus`. The swap models keep one samples
#' \eqn{\times} taxa map per thread (4 bytes per cell), and use fewer
#' threads when the maps would take more than 1 GB.
#'
#'
#' **Differences from picante**
#'
#' picante 1.8.4's `ses.mpd()` and `ses.mntd()` draw with R's random
#' numbers, ecodive with its own (PCG), so the null values differ draw by
#' draw; given the same random communities, the results agree to about
#' \eqn{10^{-12}}. Other differences:
#'
#' * picante's `"sample.pool"` makes the same communities as its
#'   `"richness"`, drawing from every column of the table, including taxa
#'   absent from every sample; ecodive's `"sample_pool"` draws only from
#'   taxa present in some sample. They agree when no taxon is absent from
#'   every sample.
#' * picante's shuffles for `"richness"` and `"frequency"` are not quite
#'   uniform (each position's chance of holding a value varies by a few
#'   percent); ecodive's are.
#' * With `"frequency"`, a random sample can have fewer than two taxa.
#'   picante leaves those draws out of the mean and standard deviation but
#'   still divides the rank by `runs + 1`; ecodive divides by the number of
#'   draws used plus one.
#' * When a sample has fewer than two taxa, picante gives `NaN` for
#'   `mpd.rand.mean` (all draws `NA`) and an arbitrary rank and p-value if
#'   it isn't `NA`; ecodive gives `NA`.
#' * Random values within \eqn{10^{-12}} (relative) of `obs` count as ties.
#' * Each swap-model draw is a chain of its own from the observed table, as
#'   in picante (phylocom instead continues one chain). When no swap is
#'   possible (every pair of samples nested), picante's independent swap
#'   never ends; ecodive warns and uses the observed table.
#'
#'
#' **Same names in other packages**
#'
#' picante calls these `ses.mpd()` and `ses.mntd()`; there are no
#' functions named `nri()` or `nti()`.
#'
#' @references
#' Webb, C. O., Ackerly, D. D., McPeek, M. A., & Donoghue, M. J. (2002).
#' Phylogenies and community ecology. *Annual Review of Ecology and
#' Systematics*, 33, 475-505.
#' \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). Phylocom: software
#' for the analysis of phylogenetic community structure and trait
#' evolution. *Bioinformatics*, 24(18), 2098-2100.
#' \doi{10.1093/bioinformatics/btn358}
#'
#' Kembel, S. W. (2009). Disentangling niche and neutral influences on
#' community assembly: assessing the performance of community phylogenetic
#' structure tests. *Ecology Letters*, 12(9), 949-960.
#' \doi{10.1111/j.1461-0248.2009.01354.x}
#'
#' Hardy, O. J. (2008). Testing the spatial phylogenetic structure of local
#' communities: statistical performances of different null models and test
#' statistics on a locally neutral community. *Journal of Ecology*, 96(5),
#' 914-926. \doi{10.1111/j.1365-2745.2008.01421.x}
#'
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence
#' patterns. *Ecology*, 81(9), 2606-2621.
#' \doi{10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2}
#'
#' Miklós, I., & Podani, J. (2004). Randomization of presence-absence
#' matrices: comments and new algorithms. *Ecology*, 85(1), 86-92.
#' \doi{10.1890/03-0101}
#'
#' Tsirogiannis, C., Sandel, B., & Cheliotis, D. (2012). Efficient
#' computation of popular phylogenetic tree measures. *Algorithms in
#' Bioinformatics* (WABI 2012), LNCS 7534, 30-43.
#' \doi{10.1007/978-3-642-33122-0_3}
#'
#' @inheritSection documentation Input Types
#'
#' @family Phylogenetic metrics
#' @seealso [mpd()], [mntd()]; [faith()].
#'
#' @export
#' @examples
#'     # Net Relatedness Index, against shuffled tip labels.
#'     res <- mpd_ses(ex_counts, ex_tree, n_draws = 99)
#'     res
#'     plot(res)
#'
#'     # Nearest Taxon Index, abundance-weighted, holding richness.
#'     mntd_ses(ex_counts, ex_tree, weighted = TRUE,
#'              null_model = 'richness', n_draws = 99)
#'
#'     # Exact null moments for MPD: no draws.
#'     mpd_ses(ex_counts, ex_tree, null_model = 'exact')
#'
mpd_ses <- function (
    counts, tree = NULL, taxa_dm = NULL, weighted = FALSE,
    null_model = 'taxa_labels', n_draws = 999L, n_swaps = 1000L,
    margin = 1L, seed = 0L, cpus = n_cpus() ) {

  phylo_ses(
    counts, tree, taxa_dm, 0L, weighted, null_model, n_draws, n_swaps,
    margin, seed, cpus )
}

#' @rdname mpd_ses
#' @export
mntd_ses <- function (
    counts, tree = NULL, taxa_dm = NULL, weighted = FALSE,
    null_model = 'taxa_labels', n_draws = 999L, n_swaps = 1000L,
    margin = 1L, seed = 0L, cpus = n_cpus() ) {

  phylo_ses(
    counts, tree, taxa_dm, 1L, weighted, null_model, n_draws, n_swaps,
    margin, seed, cpus )
}



#------------------------------------------------------------------------------
# Internals.
#------------------------------------------------------------------------------

# C codes of the null models (src/phylo_structure.c).
PHYLO_NULLS <- c(
  taxa_labels = 1L, richness = 2L, frequency = 3L, sample_pool = 4L,
  phylogeny_pool = 5L, independent_swap = 6L, trial_swap = 7L, exact = 8L )

# Models that draw from every tip of the tree (or label of taxa_dm).
PHYLO_ALL_TIPS <- c('taxa_labels', 'phylogeny_pool', 'exact')


phylo_null_model <- function (null_model, stat) {
  ok <- names(PHYLO_NULLS)
  if (stat == 1L) ok <- setdiff(ok, 'exact')
  key <- function (s) gsub('[._ ]', '', tolower(s))
  i <- if (is.character(null_model) && length(null_model) == 1 && !is.na(null_model))
    match(key(null_model), key(ok))
  if (is.null(i) || is.na(i)) {
    msg <- paste0("`null_model` must be one of ", name_list(ok, max = 10), ".")
    if (identical(key(null_model), 'exact'))
      msg <- paste(msg, 'No exact null moments are implemented for MNTD.')
    stop(msg, call. = FALSE)
  }
  ok[[i]]
}


phylo_observed <- function (counts, tree, taxa_dm, stat, weighted, margin, cpus) {
  res <- phylo_core(counts, tree, taxa_dm, stat, weighted, 'none', 0L, 0L, margin, 0L, cpus)
  stats::setNames(res$obs, res$samples)
}


phylo_ses <- function (
    counts, tree, taxa_dm, stat, weighted, null_model, n_draws, n_swaps,
    margin, seed, cpus, null_tables = NULL ) {

  if (!is.logical(weighted) || length(weighted) != 1 || is.na(weighted))
    stop('`weighted` must be TRUE or FALSE.', call. = FALSE)
  model <- phylo_null_model(null_model, stat)

  if (!is.numeric(n_draws) || length(n_draws) != 1 || !is.finite(n_draws) ||
      n_draws < 1 || n_draws %% 1 != 0 || n_draws > 1e8)
    stop('`n_draws` must be a whole number from 1 to 1e8.', call. = FALSE)
  if (!is.numeric(n_swaps) || length(n_swaps) != 1 || !is.finite(n_swaps) ||
      n_swaps < 0 || n_swaps %% 1 != 0 || n_swaps > .Machine$integer.max)
    stop('`n_swaps` must be a whole number of at least 0.', call. = FALSE)
  n_draws <- as.integer(n_draws)
  n_swaps <- as.integer(n_swaps)

  res <- phylo_core(
    counts, tree, taxa_dm, stat, weighted, model, n_draws, n_swaps,
    margin, seed, cpus, null_tables )

  if (isTRUE(res$nested) && model %in% c('independent_swap', 'trial_swap') &&
      is.null(null_tables))
    warning(
      'No swap is possible: every pair of samples is nested (one holds all ',
      'the taxa of the other), so every null table is the observed one.',
      call. = FALSE )

  exact <- model == 'exact'
  z     <- (res$obs - res$mean) / res$sd
  used  <- if (exact) rep(NA_integer_, length(res$obs)) else res$n_used
  p     <- if (exact) rep(NA_real_, length(res$obs)) else res$rank / (used + 1)

  df <- data.frame(
    sample    = if (is.null(res$samples)) as.character(seq_along(res$obs)) else res$samples,
    n_taxa    = res$n_taxa,
    obs       = res$obs,
    null_mean = res$mean,
    null_sd   = res$sd,
    z         = z,
    index     = -z,
    rank      = res$rank,
    p_value   = p,
    n_draws   = used,
    stringsAsFactors = FALSE )
  names(df)[[7]] <- if (stat == 0L) 'nri' else 'nti'

  method <- if (stat == 0L) 'mpd_ses' else 'mntd_ses'
  attr(df, 'method')     <- method
  attr(df, 'null_model') <- model
  attr(df, 'n_draws')    <- if (exact) NA_integer_ else n_draws
  if (model %in% c('independent_swap', 'trial_swap')) attr(df, 'n_swaps') <- n_swaps
  attr(df, 'weighted')   <- weighted
  class(df) <- c(paste0('ecodive_', method), 'data.frame')
  set_display(df, title = if (stat == 0L) 'Net Relatedness Index (NRI)' else 'Nearest Taxon Index (NTI)')
}


# The work shared by the six functions.
#
#   stat         0L MPD, 1L MNTD.
#   model        'none' (observed only) or a name of PHYLO_NULLS.
#   null_tables  Internal hook, for exact comparisons with other software:
#                a list of n_draws matrices (samples x taxa, rows in the
#                order of `counts`' samples, columns named by tips of `tree`
#                or labels of `taxa_dm`), used as the random communities
#                instead of drawing them. `model` is then only a label.
#
# Returns list(samples, n_taxa, obs, mean, sd, rank, n_used, nested).
phylo_core <- function (
    counts, tree, taxa_dm, stat, weighted, model, n_draws, n_swaps, margin,
    seed, cpus, null_tables = NULL ) {

  if (!is.null(tree) && !is.null(taxa_dm))
    stop('Give either `tree` or `taxa_dm`, not both.', call. = FALSE)
  use_dm <- !is.null(taxa_dm)

  validate_counts()
  validate_margin()
  validate_seed()
  validate_cpus()

  if (use_dm) tree <- NULL  # An embedded tree isn't used.

  switch(
    EXPR = get_matrix_package(counts),
    'base'   = { dims <- dim(counts);                 dn <- dimnames(counts) },
    'slam'   = { dims <- c(counts$nrow, counts$ncol); dn <- counts$dimnames  },
    'Matrix' = { dims <- counts@Dim;                  dn <- counts@Dimnames  } )
  if (is.null(dn)) dn <- list(NULL, NULL)
  if (margin == 2L) { dims <- rev(dims); dn <- rev(dn) }

  if (is.logical(counts) && is.matrix(counts)) storage.mode(counts) <- 'double'
  vals <- switch(get_matrix_package(counts), base = counts, slam = counts$v, Matrix = counts@x)
  if (!is.numeric(vals))
    stop('`counts` must be numeric.', call. = FALSE)
  remove('vals')

  n        <- as.integer(dims[[1]])
  D        <- as.integer(dims[[2]])
  samples  <- if (length(dn[[1]])) as.character(dn[[1]])
  features <- if (length(dn[[2]])) as.character(dn[[2]])

  supplied <- !is.null(null_tables)
  all_tips <- supplied || model %in% PHYLO_ALL_TIPS


  # The labels of the taxa that distances are needed for: the table's
  # taxa first (in order), then the other tips, for the models that use
  # them.
  epar <- echd <- elen <- dm <- NULL
  if (use_dm) {
    tdm    <- phylo_taxa_dm(taxa_dm, features, D, all_tips)
    dm     <- tdm$dm
    labels <- tdm$labels
  }
  else {
    if (is.null(tree))
      stop(
        '`tree` or `taxa_dm` is required: a phylo object whose tips are the ',
        'taxa of `counts` (or a phyloseq, rbiom or TreeSummarizedExperiment ',
        'object that holds one), or the dissimilarities among the taxa.',
        call. = FALSE )
    if (!inherits(tree, 'phylo') || is.null(tree$tip.label))
      stop('`tree` is not a valid phylo object.', call. = FALSE)
    if (is.null(features))
      stop(
        '`counts` must have taxa names, to match them to the tips of `tree`.',
        call. = FALSE )
    tips   <- as.character(tree$tip.label)
    labels <- if (all_tips) c(features, setdiff(tips, features)) else features
    if (is.null(tree$edge.length) || anyNA(tree$edge.length))
      stop('`tree` must have branch lengths, with none missing.', call. = FALSE)
    if (!is.numeric(tree$edge.length) || any(tree$edge.length < 0))
      stop('`tree` must not have negative branch lengths.', call. = FALSE)
    if (length(labels) >= 2L) {
      tr   <- philr_tree(tree, labels)
      epar <- tr$edge_par
      echd <- tr$edge_child
      elen <- tr$edge_len
    }
    else if (length(setdiff(features, tips)))
      philr_tree(tree, labels) # for its error message
  }
  T <- length(labels)


  # Supplied random communities: one table of n_draws x n rows, with the
  # T labels as columns.
  sup <- NULL
  if (supplied) {
    if (!is.list(null_tables) || length(null_tables) < 1L)
      stop('`null_tables` must be a list of tables.', call. = FALSE) # nocov
    n_draws <- length(null_tables)
    sup <- matrix(0, n * n_draws, T)
    for (r in seq_len(n_draws)) {
      x <- as.matrix(null_tables[[r]])
      j <- match(colnames(x), labels)
      if (nrow(x) != n || anyNA(j))
        stop('`null_tables[[', r, ']]` must have one row per sample and columns ',
             'named by tips.', call. = FALSE) # nocov
      sup[(r - 1L) * n + seq_len(n), j] <- x
    }
    model <- 'supplied'
  }

  code <- switch(model, none = 0L, supplied = 9L, PHYLO_NULLS[[model]])

  if (T < 2L) {
    # At most one taxon: every MPD and MNTD is NA.
    ntaxa <- if (n > 0) phylo_richness(counts, margin, n) else integer(0)
    na    <- rep(NA_real_, n)
    return (list(
      samples = samples, n_taxa = ntaxa, obs = na, mean = na, sd = na,
      rank = na, n_used = rep(0L, n), nested = TRUE ))
  }

  res <- .Call(
    C_phylo_structure, counts, margin, epar, echd, elen, dm, T,
    as.integer(stat), as.logical(weighted), code, as.integer(n_draws),
    as.integer(n_swaps), as.double(seed), as.integer(cpus), sup )

  c(list(samples = samples), res)
}


# Each sample's number of non-zero values.
phylo_richness <- function (counts, margin, n) {
  m <- as.matrix(counts)
  if (margin == 2L) m <- t(m)
  as.integer(rowSums(m > 0))
}


# taxa_dm as a packed lower triangle over its labels: the taxa of `counts`
# first, in order, then (with all_tips) its other labels.
phylo_taxa_dm <- function (taxa_dm, features, D, all_tips) {

  v <- tryCatch(
    dist_vector(taxa_dm),
    error = function (e)
      stop(gsub('`dm`', '`taxa_dm`', conditionMessage(e), fixed = TRUE), call. = FALSE) )

  if (anyNA(v$dm) || any(!is.finite(v$dm)) || any(v$dm < 0))
    stop('`taxa_dm` must have finite, non-negative values.', call. = FALSE)

  labels <- if (!is.null(v$labels)) as.character(v$labels)

  if (is.null(labels) || is.null(features)) {
    if (v$n != D)
      stop(
        '`taxa_dm` must have one row per taxon of `counts` (', D, '), or ',
        'labels to match them by.', call. = FALSE )
    return (list(dm = v$dm, labels = if (is.null(features)) as.character(seq_len(D)) else features))
  }

  if (anyDuplicated(labels))
    stop('`taxa_dm` has duplicated labels: ',
      dQuote(labels[anyDuplicated(labels)], FALSE), '.', call. = FALSE)
  if (anyDuplicated(features))
    stop('`counts` has duplicated taxa names: ',
      dQuote(features[anyDuplicated(features)], FALSE), '.', call. = FALSE)

  missing <- setdiff(features, labels)
  if (length(missing))
    stop(
      length(missing), ' tax', if (length(missing) > 1) 'a' else 'on', ' of `counts` ',
      if (length(missing) > 1) 'are' else 'is', ' not in `taxa_dm`: ',
      paste(dQuote(utils::head(missing, 5), FALSE), collapse = ', '),
      if (length(missing) > 5) ', ...', '.', call. = FALSE )

  keep <- if (all_tips) c(features, setdiff(labels, features)) else features
  i    <- match(keep, labels)
  if (identical(i, seq_len(v$n))) return (list(dm = v$dm, labels = keep))

  # Reorder (and subset): element (r, c), r > c, of the packed triangle.
  n   <- v$n
  K   <- length(keep)
  pos <- function (r, c) (c - 1) * n - (c - 1) * c / 2 + (r - c)
  rc  <- which(lower.tri(matrix(0, K, K)), arr.ind = TRUE)
  a   <- i[rc[, 1]]
  b   <- i[rc[, 2]]
  list(dm = v$dm[pos(pmax(a, b), pmin(a, b))], labels = keep)
}



#------------------------------------------------------------------------------
# plot()
#------------------------------------------------------------------------------

#' @rdname mpd_ses
#' @param x   A result of `mpd_ses()` or `mntd_ses()`.
#' @param groups   For `plot()`: optional groups of the samples (one value
#'        per row of `x`), which split and colour the dot chart.
#' @param ...   For `plot()`: graphical parameters for
#'        [graphics::dotchart()], such as `main`, `xlab` or `cex`.
#' @export
plot.ecodive_mpd_ses <- function (x, groups = NULL, ...) {
  phylo_ses_plot(x, groups, 'nri', 'Net Relatedness Index (NRI)', ...)
}

#' @rdname mpd_ses
#' @export
plot.ecodive_mntd_ses <- function (x, groups = NULL, ...) {
  phylo_ses_plot(x, groups, 'nti', 'Nearest Taxon Index (NTI)', ...)
}


phylo_ses_plot <- function (x, groups, col, fallback, ...) {

  args <- list(...)
  if (!col %in% names(x) || nrow(x) == 0)
    stop('`x` has no `', col, '` values to plot.', call. = FALSE)

  v   <- x[[col]]
  lab <- if ('sample' %in% names(x)) as.character(x$sample) else as.character(seq_along(v))
  p   <- if ('p_value' %in% names(x)) x$p_value else rep(NA_real_, length(v))
  sig <- (p < 0.05 | p > 0.95) %in% TRUE

  ord <- order(v, na.last = NA)
  if (!is.null(groups)) {
    if (length(groups) != nrow(x))
      stop('`groups` must have one value per row of `x`.', call. = FALSE)
    g   <- factor(groups)
    ord <- ord[order(g[ord], v[ord])]
  }
  if (length(ord) == 0)
    stop('`x` has no finite `', col, '` values to plot.', call. = FALSE)

  model <- attr(x, 'null_model', exact = TRUE)
  nd    <- attr(x, 'n_draws', exact = TRUE)
  sub   <- if (is.null(model)) NULL else if (identical(model, 'exact'))
    'Null model: exact moments' else sprintf('Null model: %s, %s draws', model, plot_num(nd, 6L))

  cols <- if (is.null(groups)) ifelse(sig[ord], PLOT_HIGHLIGHT, PLOT_MUTED)
          else plot_palette(nlevels(g))[as.integer(g[ord])]

  plot_do(
    graphics::dotchart, c(list(v[ord]), plot_dot_args(args)),
    list(
      labels = lab[ord], groups = if (is.null(groups)) NULL else droplevels(g[ord]),
      pch    = ifelse(sig[ord], 19, 1),
      color  = cols,
      xlim   = range(c(-1.96, 1.96, v[ord]), finite = TRUE),
      main   = plot_title(x, fallback),
      xlab   = paste0(toupper(col), ' (positive: clustered)') ))

  graphics::abline(v = c(-1.96, 0, 1.96), lty = c(2, 3, 2), col = PLOT_MUTED)
  graphics::title(sub = plot_arg(args, 'sub', sub))
  invisible(x)
}
