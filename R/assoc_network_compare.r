# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Comparing two association networks, as NetCoMi's netCompare() does,
# written from the papers cited and NetCoMi's documentation and outputs
# (no NetCoMi or igraph source was read). See plans/assoc_network.md.



#' Compare Two Association Networks
#'
#' Compares two groups' association networks made by [assoc_network()]:
#' the differences in each feature's centralities and in the networks'
#' global properties, how much their most central features and hubs
#' overlap (Jaccard index), and how similar their clusterings are (adjusted
#' Rand index), as NetCoMi's `netCompare()` does. Given the count table
#' and the groups, it also tests the differences by permuting the samples'
#' groups.
#'
#'
#' @param x,y  Two networks from [assoc_network()] over the same features
#'        (matched by name): `x` of the first group, `y` of the second.
#'
#' @param quantile  The "most central" features of a network are those
#'        above this quantile of a centrality (type 7, see
#'        [stats::quantile()]), over all the features. Default `0.75`, as
#'        NetCoMi's `jaccQuant`.
#'
#' @param permutations  How many permutations for the tests: of the
#'        features' cluster labels, for the adjusted Rand index; and, with
#'        `counts` and `groups`, of the samples' groups. Default `1000L`.
#'        `0L` for no tests.
#'
#' @param counts  The count table both networks were estimated from (all
#'        samples of both groups), or `NULL` (the default) for no
#'        group-permutation tests. Any of the objects listed under **Input
#'        Types**; with `margin = 2` (features as rows), give `margin` in
#'        `assoc_args`.
#'
#' @param groups  Each sample's group, two different values (or the name of
#'        a column of `data`). `x` was made from the samples of the first
#'        group (the first level of a factor, else the first in sorted
#'        order), `y` from those of the second.
#'
#' @param data  A data frame with the samples' metadata, when `groups` is a
#'        column name.
#'
#' @param assoc_args  A list of the arguments that made the associations of
#'        `x` and `y` (besides `counts`), for [feature_assoc()] with the
#'        networks' method: e.g. `list(n_draws = 20)` for [sparcc()]. The
#'        same arguments must be given again, as results don't record them.
#'
#' @param network_args  A list of the arguments that made `x` and `y` with
#'        [assoc_network()] (besides `x`), e.g. `list(threshold = 0.4)`.
#'
#' @param seed  An integer seed for the random permutations, and for the
#'        seeds of each permutation's association runs. The same `seed`
#'        gives the same results, for any number of `cpus`.
#'
#' @param cpus  How many threads to use, for the association runs and the
#'        networks' shortest paths. Results are the same for any number.
#'
#' @param p_adjust  How to adjust the centralities' permutation p-values
#'        for testing every feature: one of [stats::p.adjust.methods].
#'        Default `'BH'`.
#'
#' @return A list of class `ecodive_network_compare`:
#'
#'   * `nodes`: a data frame, one row per feature and centrality, with
#'     `feature`, `centrality` (`"degree"`, `"betweenness"`, `"closeness"`
#'     or `"eigenvector"`), `x`, `y`, `difference` (`x - y`), and, with the
#'     group-permutation test, `p_value` and `p_adjusted` (over the
#'     features, for each centrality).
#'   * `global`: a data frame, one row per global property of
#'     [assoc_network()], with `property`, `x`, `y`, `difference` and,
#'     with the group-permutation test, `p_value`.
#'   * `jaccard`: a data frame, one row per set (`"degree"`,
#'     `"betweenness"`, `"closeness"`, `"eigenvector"`: the most central
#'     features; `"hubs"`), with `jaccard` (the Jaccard index of the two
#'     networks' sets), and `p_less` and `p_greater`, the probabilities of
#'     a Jaccard index at most and at least as large by chance (see
#'     **Details**).
#'   * `ari`: the adjusted Rand index of the two clusterings, `value`, and
#'     its permutation `p_value`.
#'   * `clusters`: a data frame, one row per feature, with `feature`,
#'     `cluster_x`, `cluster_y`, `hub_x` and `hub_y`: each network's
#'     clusters and hubs.
#'   * `edges_x`, `edges_y`: the two networks' `edges`, which `plot()`
#'     draws ([plot.ecodive_network_compare()]).
#'
#'   Attributes: `method`, `quantile`, `permutations`, `groups` (the two
#'   group names, when given), and `display` (the title; see
#'   [ecodive_display]).
#'
#' @details
#'
#' **Most central features.** Those above the `quantile` quantile of a
#' centrality over all the features (0 outside the largest connected
#' component), as for hubs in [assoc_network()].
#'
#' **Jaccard index.** For two sets \eqn{A} and \eqn{B} of features, \eqn{J
#' = |A \cap B| / |A \cup B|} (0 when both are empty). `p_less` and
#' `p_greater` are \eqn{P(J \le j)} and \eqn{P(J \ge j)} under Real &
#' Vargas's (1996) model, as NetCoMi computes them: each of the \eqn{n =
#' |A \cup B|} features is in both sets with probability 1/3, so
#' \eqn{|A \cap B|} is Binomial(\eqn{n}, 1/3). (This depends only on the
#' union, not on how many features there are.)
#'
#' **Adjusted Rand index** (Hubert & Arabie 1985) of the two networks'
#' `cluster` columns as they are, so cluster 0 (features alone in their
#' cluster) is one cluster, as in NetCoMi. Its p-value tests the null
#' hypothesis ARI = 0: the cluster labels of `y` are permuted over the
#' features `permutations` times, and \eqn{p = (1 + \#\{|ARI_b| \ge
#' |ARI|\}) / (B + 1)}. This is two-sided, as NetCoMi's p-values are
#' (its documentation describes a one-sided test; its outputs, compared
#' over 20,000 permutations, are two-sided); NetCoMi adds no 1 to the
#' counts.
#'
#' **Group-permutation tests** (with `counts` and `groups`). Each
#' permutation shuffles the samples' groups, keeping the groups' sizes,
#' re-estimates both groups' associations with [feature_assoc()] (with
#' `assoc_args`) and rebuilds both networks with [assoc_network()] (with
#' `network_args`). A difference's p-value is the proportion of the
#' permutations, counting the observed one, with an absolute difference at
#' least as large. This runs the association method \eqn{2B} times: with
#' [sparcc()], lower its `n_draws` in `assoc_args`, or `permutations`.
#' Permutation \eqn{b}'s association runs get their own `seed`, drawn from
#' (`seed`, \eqn{b}), replacing any `seed` in `assoc_args`.
#'
#' **Names in NetCoMi.** `quantile` is `jaccQuant`, `permutations` is
#' `nPermRand` and `nPerm` (with `permTest = TRUE`), and `p_adjust` is
#' `adjust`. NetCoMi's results `jaccDeg`, `jaccBetw`, `jaccClose`,
#' `jaccEigen`, `jaccHub` and `randInd` are `jaccard` and `ari` here.
#'
#' **Differences from NetCoMi:** p-values are adjusted with
#' [stats::p.adjust()] only (NetCoMi's default is "adaptBH"); the ARI test's
#' p-value counts the observed value; the permutations come from ecodive's
#' own generator, so they match NetCoMi's in distribution only; no graphlet
#' distance; no comparison of the LCCs alone.
#'
#' @inheritSection documentation Input Types
#'
#' @references
#' Peschel, S., Müller, C. L., von Mutius, E., Boulesteix, A.-L., & Depner,
#' M. (2021). NetCoMi: network construction and comparison for microbiome
#' data in R. *Briefings in Bioinformatics*, 22, bbaa290.
#'
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
#' Classification*, 2, 193–218.
#'
#' Real, R., & Vargas, J. M. (1996). The probabilistic basis of Jaccard's
#' index of similarity. *Systematic Biology*, 45, 380–385.
#'
#' @seealso [assoc_network()], [assoc_network_diff()], [feature_assoc()].
#'
#' @export
#' @examples
#'     # Two groups of 60 samples of 15 taxa
#'     set.seed(1)
#'     make <- function (linked) {
#'       basis <- matrix(rnorm(60 * 15), 60, 15)
#'       basis[, linked] <- basis[, linked] + basis[, linked[1]]
#'       matrix(rpois(900, exp(basis + 4)), 60, 15)
#'     }
#'     counts <- rbind(make(1:5), make(6:10))
#'     colnames(counts) <- paste0('OTU', 1:15)
#'     groups <- rep(c('a', 'b'), each = 60)
#'
#'     nx <- assoc_network(sparcc(counts[groups == 'a', ], n_draws = 10))
#'     ny <- assoc_network(sparcc(counts[groups == 'b', ], n_draws = 10))
#'     cmp <- assoc_network_compare(nx, ny)
#'     cmp$jaccard
#'     cmp$ari
#'     plot(cmp)
#'
#'     # Group-permutation tests (slow: 2 x 19 SparCC runs)
#'     cmp <- assoc_network_compare(
#'       nx, ny, permutations = 19, counts = counts, groups = groups,
#'       assoc_args = list(n_draws = 10) )
#'     cmp$global
#'
assoc_network_compare <- function (
    x, y, quantile = 0.75, permutations = 1000L,
    counts = NULL, groups = NULL, data = NULL,
    assoc_args = list(), network_args = list(),
    p_adjust = 'BH', seed = 0L, cpus = n_cpus() ) {

  network_compare_core(
    x, y, quantile, permutations, counts, groups, data, assoc_args,
    network_args, p_adjust, seed, cpus )
}



# assoc_network_compare(), with two internal hooks: `ari_perms`, a matrix
# of permutations of the features (n x B, 1-based) for the ARI test, and
# `group_perms`, one of the samples (n_samples x B) for the group test,
# replacing the random ones (for exact comparisons with other software).
network_compare_core <- function (
    x, y, quantile = 0.75, permutations = 1000L,
    counts = NULL, groups = NULL, data = NULL,
    assoc_args = list(), network_args = list(),
    p_adjust = 'BH', seed = 0L, cpus = n_cpus(),
    ari_perms = NULL, group_perms = NULL ) {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  if (!inherits(x, 'ecodive_network') || !inherits(y, 'ecodive_network'))
    stop('`x` and `y` must be networks from assoc_network().', call. = FALSE)

  if (!is.numeric(quantile) || length(quantile) != 1L || is.na(quantile) || quantile < 0 || quantile > 1)
    stop('`quantile` must be a single number from 0 to 1.', call. = FALSE)

  if (!is.numeric(permutations) || length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 0 || permutations %% 1 != 0 || permutations > 1e8)
    stop('`permutations` must be a whole number from 0 to 1e8.', call. = FALSE)
  permutations <- as.integer(permutations)

  if (!is.list(assoc_args) || !is.list(network_args))
    stop('`assoc_args` and `network_args` must be lists.', call. = FALSE)

  p_adjust <- network_p_adjust(p_adjust)
  validate_seed()
  validate_cpus()

  mx <- attr(x, 'method', exact = TRUE)
  my <- attr(y, 'method', exact = TRUE)
  if (!identical(mx, my))
    stop('`x` and `y` must be made from associations of the same method.', call. = FALSE)

  fx <- x$nodes$feature
  fy <- y$nodes$feature
  if (length(fx) != length(fy) || !setequal(fx, fy) || anyDuplicated(fx))
    stop('`x` and `y` must be over the same features (matched by name).', call. = FALSE)
  iy <- match(fx, fy)
  ny <- y$nodes[iy, , drop = FALSE]
  nx <- x$nodes
  n  <- length(fx)

  cents <- c('degree', 'betweenness', 'closeness', 'eigenvector')


  #--------------------------------------------------------
  # Observed differences.
  #--------------------------------------------------------
  nodes <- data.frame(
    feature    = rep(fx, times = 4L),
    centrality = rep(cents, each = n),
    x          = unlist(nx[cents], use.names = FALSE),
    y          = unlist(ny[cents], use.names = FALSE),
    stringsAsFactors = FALSE )
  nodes$difference <- nodes$x - nodes$y

  props  <- names(x$global)
  global <- data.frame(
    property   = props,
    x          = unname(x$global),
    y          = unname(y$global[props]),
    stringsAsFactors = FALSE )
  global$difference <- global$x - global$y


  #--------------------------------------------------------
  # Jaccard indices of the most central sets, and of hubs.
  #--------------------------------------------------------
  sets <- c(cents, 'hubs')
  jac  <- lapply(sets, function (s) {
    if (s == 'hubs') {
      a <- nx$hub
      b <- ny$hub
    }
    else {
      a <- network_most_central(nx[[s]], quantile)
      b <- network_most_central(ny[[s]], quantile)
    }
    network_jaccard(a, b)
  })
  jaccard <- data.frame(
    set       = sets,
    jaccard   = vapply(jac, `[[`, numeric(1), 'jaccard'),
    p_less    = vapply(jac, `[[`, numeric(1), 'p_less'),
    p_greater = vapply(jac, `[[`, numeric(1), 'p_greater'),
    stringsAsFactors = FALSE )


  #--------------------------------------------------------
  # Adjusted Rand index, tested by permuting y's labels.
  #--------------------------------------------------------
  cx  <- nx$cluster
  cy  <- ny$cluster
  ari <- adjusted_rand(cx, cy)

  if (!is.null(ari_perms)) {
    pm <- ari_perms
  }
  else if (permutations > 0L) {
    pm <- .Call(C_perm_indices, n, permutations, NULL, seed)
  }
  else pm <- NULL

  ari_p <- if (is.null(pm) || is.na(ari)) NA_real_ else {
    perm <- apply(pm, 2L, function (p) adjusted_rand(cx, cy[p]))
    perm_p(ari, perm, 'two.sided')
  }


  #--------------------------------------------------------
  # Group-permutation tests.
  #--------------------------------------------------------
  grp_names <- NULL
  if (!is.null(counts)) {

    if (is.null(groups))
      stop('`counts` needs `groups`: the samples\' two groups.', call. = FALSE)
    if (is.null(mx))
      stop(
        'The group-permutation tests re-estimate the associations, so `x` ',
        'and `y` must be made from results of an association method, not ',
        'matrices.', call. = FALSE )

    gt <- network_group_test(
      x, y, mx, counts, groups, data, assoc_args, network_args,
      permutations, seed, group_perms, fx )

    nodes$p_value    <- gt$node_p
    nodes$p_adjusted <- unlist(lapply(split(gt$node_p, nodes$centrality)[cents],
                                      stats::p.adjust, method = p_adjust),
                               use.names = FALSE)
    global$p_value   <- gt$global_p
    grp_names        <- gt$levels
  }
  else if (!is.null(groups)) {
    stop('`groups` needs `counts`: the count table both networks were made from.', call. = FALSE)
  }

  res <- structure(
    list(
      nodes   = nodes,
      global  = global,
      jaccard = jaccard,
      ari     = c(value = ari, p_value = ari_p),
      clusters = data.frame(
        feature   = fx,
        cluster_x = nx$cluster,
        cluster_y = ny$cluster,
        hub_x     = nx$hub,
        hub_y     = ny$hub,
        stringsAsFactors = FALSE ),
      edges_x = x$edges,
      edges_y = y$edges ),
    class = 'ecodive_network_compare' )
  if (!is.null(mx)) attr(res, 'method') <- mx
  attr(res, 'quantile')     <- as.double(quantile)
  attr(res, 'permutations') <- if (is.null(pm)) 0L else ncol(pm)
  if (!is.null(grp_names)) attr(res, 'groups') <- grp_names
  set_display(res, title = title_of('Comparison', display_title(x)))
}



# The "most central" nodes: above the quantile of v, over all nodes.
network_most_central <- function (v, quantile) {
  above_quantile(v, quantile)
}



# The Jaccard index of two sets (logical vectors over the same features),
# and the probabilities of an index at most and at least as large by
# chance, as Real & Vargas (1996) derive them: each of the n features in
# either set is in both with probability 1/3, so the number shared is
# Binomial(n, 1/3), and J is that number over n. Two empty sets give J = 0
# and probabilities 1 (as NetCoMi).
network_jaccard <- function (a, b) {
  k <- sum(a & b)
  u <- sum(a | b)
  list(
    jaccard   = if (u == 0L) 0 else k / u,
    p_less    = stats::pbinom(k, u, 1/3),
    p_greater = stats::pbinom(k - 1, u, 1/3, lower.tail = FALSE) )
}



# The adjusted Rand index (Hubert & Arabie 1985) of two labelings.
adjusted_rand <- function (a, b) {
  n <- length(a)
  if (n < 2L) return (NA_real_)
  tab <- table(a, b)
  c2  <- function (v) sum(v * (v - 1) / 2)
  sij <- c2(tab)
  sa  <- c2(rowSums(tab))
  sb  <- c2(colSums(tab))
  exp <- sa * sb / (n * (n - 1) / 2)
  mx  <- (sa + sb) / 2
  if (mx == exp) return (if (sij == exp) 1 else NA_real_)
  (sij - exp) / (mx - exp)
}



# The group-permutation tests of assoc_network_compare(): returns
# list(node_p (4 x n, centrality-major, as `nodes`), global_p, levels).
network_group_test <- function (
    x, y, method, counts, groups, data, assoc_args, network_args,
    permutations, seed, group_perms, features ) {

  # The count table, and which margin holds the samples.
  margin <- assoc_args$margin
  if (is.null(margin)) margin <- 1L
  env <- environment()
  validate_counts(env)
  n   <- if (margin == 1L) nrow(counts) else ncol(counts)
  lab <- if (margin == 1L) rownames(counts) else colnames(counts)

  g <- test_group_values(groups, data, n, lab, ref = '`counts`')$groups
  g <- if (is.factor(g)) droplevels(g) else factor(g)
  if (nlevels(g) != 2L)
    stop('`groups` must have exactly two different values.', call. = FALSE)
  lev <- levels(g)
  gi  <- as.integer(g)

  info  <- match_feature_assoc(method)
  aargs <- assoc_args
  aargs$margin <- margin
  nargs <- network_args
  nargs$x <- NULL

  pm <- if (!is.null(group_perms)) group_perms
        else .Call(C_perm_indices, n, permutations, NULL, seed)
  B  <- ncol(pm)

  cents <- c('degree', 'betweenness', 'closeness', 'eigenvector')
  stats_of <- function (nx, ny) {
    iy <- match(features, ny$nodes$feature)
    c(unlist(nx$nodes[cents], use.names = FALSE) -
        unlist(ny$nodes[iy, cents], use.names = FALSE),
      unname(nx$global - ny$global[names(nx$global)]))
  }
  obs <- stats_of(x, y)

  sub <- function (idx) {
    if (margin == 1L) counts[idx, , drop = FALSE] else counts[, idx, drop = FALSE]
  }

  perm <- matrix(NA_real_, B, length(obs))
  for (b in seq_len(B)) {
    gb <- gi[pm[, b]]
    if ('seed' %in% info$params)
      aargs$seed <- .Call(C_pair_seed, as.integer(seed), as.integer(b))
    nets <- lapply(1:2, function (k) {
      a <- do.call(info$func, c(list(sub(which(gb == k))), aargs))
      do.call(assoc_network, c(list(a), nargs))
    })
    perm[b, ] <- stats_of(nets[[1]], nets[[2]])
  }

  p <- perm_p(obs, perm, 'two.sided')
  k <- 4L * length(features)
  list(node_p = p[seq_len(k)], global_p = p[-seq_len(k)], levels = lev)
}



#' @export
print.ecodive_network_compare <- function (x, ...) {

  cat(plot_title(x, 'Comparison of Association Networks'), '\n', sep = '')
  grp <- attr(x, 'groups', exact = TRUE)
  if (!is.null(grp)) cat('Groups: ', grp[[1]], ' (x) vs ', grp[[2]], ' (y)\n', sep = '')

  cat(sprintf(
    'Adjusted Rand index of the clusterings: %s (p = %s)\n',
    plot_num(x$ari[['value']]), plot_p(x$ari[['p_value']]) ))

  cat('Jaccard index of the most central features (above the ',
      format(attr(x, 'quantile', exact = TRUE)), ' quantile) and of hubs:\n', sep = '')
  j <- x$jaccard
  for (i in seq_len(nrow(j)))
    cat(sprintf(
      '  %-12s %s  (P(J <= j) = %s, P(J >= j) = %s)\n', j$set[[i]],
      plot_num(j$jaccard[[i]]), plot_p(j$p_less[[i]]), plot_p(j$p_greater[[i]]) ))

  cat('Global properties (x, y, difference', if (!is.null(x$global$p_value)) ', p-value', '):\n', sep = '')
  gl <- x$global
  for (i in seq_len(nrow(gl)))
    cat(sprintf(
      '  %-23s %10s %10s %10s%s\n', gl$property[[i]], plot_num(gl$x[[i]]),
      plot_num(gl$y[[i]]), plot_num(gl$difference[[i]]),
      if (is.null(gl$p_value)) '' else paste0('  ', plot_p(gl$p_value[[i]])) ))
  invisible(x)
}
