# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Association networks: the analysis layer that NetCoMi (Peschel et al.
# 2021) made standard, written from textbook formulas, the papers cited
# below and NetCoMi's documentation and outputs (NetCoMi and igraph are
# GPL: no source of either was read). See plans/assoc_network.md and
# src/assoc_network.c.



#' Association Network
#'
#' Makes a network of the features from their associations, such as
#' [sparcc()] correlations: an edge joins two features whose association is
#' strong (and significant) enough. Then measures each feature's centrality,
#' finds the hubs and clusters, and summarises the whole network, as
#' NetCoMi's `netConstruct()` and `netAnalyze()` do.
#'
#' @param x  Associations between features: a result of [sparcc()] or
#'        [rho_p()] (or of [feature_assoc()] with a similarity), or a
#'        symmetric matrix of correlations with 1 on the diagonal, such as
#'        from [stats::cor()].
#'
#' @param threshold  An edge needs an association of at least this in
#'        absolute value (\eqn{|r| \ge} `threshold`). `NULL` (the default)
#'        is the result's own cutoff when its method chooses one (for a
#'        method whose result is a selected graph, its non-zero entries;
#'        for one that returns a `threshold` attribute, that), and
#'        otherwise `0.3`, NetCoMi's default and the SparCC paper's cutoff.
#'        `0` makes every pair with a non-zero association an edge (with
#'        `alpha`, every significant pair).
#'
#' @param alpha  Also make edges only of pairs whose adjusted p-value is at
#'        most this. `NULL` (the default) for no p-value filter. Needs
#'        p-values: a [sparcc()] result made with `permutations`.
#'
#' @param p_adjust  How to adjust the p-values for testing every pair: one
#'        of [stats::p.adjust.methods]. Default `'BH'`. The adjustment is
#'        over all the pairs of `x`. Used only with `alpha`.
#'
#' @param sign  How an edge's association \eqn{r} becomes its length (a
#'        dissimilarity \eqn{d}) and weight (a similarity \eqn{s = 1 - d}):
#'        `'signed'` (the default): \eqn{d = \sqrt{(1 - r)/2}}, so that
#'        negative associations make long, weak edges; `'unsigned'`:
#'        \eqn{d = \sqrt{1 - r^2}}, so that strong negative associations
#'        make short, strong edges; `'positive'`: as `'signed'`, with the
#'        negative edges left out.
#'
#' @param clusters  How to find clusters: `'greedy'` (the default), fast
#'        greedy modularity optimisation (Clauset et al. 2004) on the edge
#'        weights \eqn{s}, or `'none'`.
#'
#' @param hubs  Which centralities make a hub: one or more of `'degree'`,
#'        `'betweenness'`, `'closeness'` and `'eigenvector'`. A hub is above
#'        the `hub_quantile` quantile of every one named. Default
#'        `'eigenvector'`.
#'
#' @param hub_quantile  The quantile a hub must be above. Default `0.95`.
#'
#' @param weighted_degree  If `TRUE`, a node's degree is the sum of its
#'        edges' weights \eqn{s} rather than their number. Default `FALSE`.
#'
#' @param cpus  How many threads to use for shortest paths (betweenness and
#'        closeness). Results are the same for any number.
#'
#' @return A list of class `ecodive_network`:
#'
#'   * `edges`: a data frame, one row per edge, with `from` and `to` (the
#'     two features, `from` first in the input table), `association`
#'     (\eqn{r}), `similarity` (\eqn{s}, the edge's weight), `distance`
#'     (\eqn{d}, its length), `sign` (`1` or `-1`), and `p_value` and
#'     `p_adjusted` when `x` has p-values. Ordered by `from`, then `to`.
#'   * `nodes`: a data frame, one row per feature (in the input table's
#'     order), with `feature`, `degree`, `betweenness`, `closeness`,
#'     `eigenvector` (computed on the largest connected component: 0 for
#'     nodes outside it), `cluster` (an integer, `0` for a feature alone in
#'     its cluster, as every feature without edges is), `hub` (logical) and
#'     `in_lcc` (logical: whether it is in the largest connected
#'     component).
#'   * `global`: a named numeric vector of the whole network's properties:
#'     `n_nodes`, `n_edges`, `density`, `positive_fraction`,
#'     `n_components`, `lcc_size`, `average_dissimilarity`,
#'     `average_path_length`, `clustering_coefficient` and `modularity`
#'     (see **Details**).
#'
#'   Attributes: `method` (the method that made `x`, such as `"sparcc"`;
#'   absent for a matrix), `n_samples` (when `x` has it), `threshold` (the
#'   cutoff used), `alpha` and `p_adjust` (when used), `sign`, `hubs` and
#'   `hub_quantile` (which define the `hub` column), and
#'   `display` (the title; see [ecodive_display]).
#'
#'   `as.data.frame()` gives the `nodes` table; `plot()` draws the network
#'   ([plot.ecodive_network()]).
#'
#' @details
#'
#' **Edges.** Features \eqn{i} and \eqn{j} are joined when \eqn{|r_{ij}|
#' \ge} `threshold` and, with `alpha`, their adjusted p-value is at most
#' `alpha`, as [assoc_pairs()] lists them. Associations slightly beyond
#' \eqn{\pm 1} (possible for [sparcc()]) are taken as \eqn{\pm 1}. As in
#' NetCoMi, whose network is a matrix of weights and whose graph of lengths
#' is a matrix of dissimilarities, a pair with weight \eqn{s = 0} (\eqn{r =
#' -1} with `'signed'`) is no edge, and an edge of length \eqn{d = 0}
#' (\eqn{r = 1}, or \eqn{|r| = 1} with `'unsigned'`) isn't used by shortest
#' paths.
#'
#' **Average dissimilarity** of a set of nodes: the mean of \eqn{d} over all
#' its pairs, counting 1 for the pairs that aren't edges (NetCoMi's
#' `avDiss`).
#'
#' **Centralities** (NetCoMi's definitions, with `centrLCC = TRUE`): computed
#' on the largest connected component (LCC), of \eqn{k} nodes (a tie for the
#' largest goes to the component with the first feature), and 0 for the
#' other nodes. \eqn{N} is the number of features.
#'
#' * Degree: the number of a node's neighbours over \eqn{N - 1}; with
#'   `weighted_degree`, the sum of its edges' weights \eqn{s}, not
#'   normalised.
#' * Betweenness: the number of shortest paths (on edge lengths \eqn{d})
#'   between pairs of other nodes that pass through the node, each pair's
#'   count divided by its number of shortest paths (Brandes 2001), over
#'   \eqn{(k - 1)(k - 2)/2}. Two paths are equally short when their lengths
#'   agree to a relative \eqn{10^{-12}}.
#' * Closeness: \eqn{\sum_u 1 / \delta(v, u)} over the other nodes of the
#'   LCC, over \eqn{k - 1}, where \eqn{\delta} is the shortest-path length
#'   divided by the LCC's average dissimilarity (NetCoMi's `sPathNorm =
#'   TRUE`).
#' * Eigenvector: the leading eigenvector of the LCC's matrix of edge
#'   weights \eqn{s}, scaled to a largest value of 1.
#'
#' **Hubs** are the nodes above the empirical `hub_quantile` quantile
#' ([stats::quantile()]'s default type 7) of each centrality in `hubs`, over
#' all \eqn{N} features (those outside the LCC count, with 0). Values equal
#' to the quantile but for rounding (a relative \eqn{10^{-12}}) are not
#' above it.
#'
#' **Clusters.** Fast greedy modularity optimisation (Clauset et al. 2004),
#' on the weights \eqn{s}: every feature starts in a cluster of its own;
#' then, one pair at a time, the two connected clusters whose merge most
#' increases the weighted modularity are merged (ties to the pair with the
#' first features), until no two clusters are connected; the partition with
#' the largest modularity is kept, except that a single cluster of every
#' feature is never kept (as igraph's `cluster_fast_greedy()`; the
#' partition one merge before it is). Clusters are numbered in the order of
#' their first features; a feature alone in its cluster, as every feature
#' without edges is, is in cluster `0`, as in NetCoMi.
#'
#' **Global properties** (of the whole network):
#'
#' * `n_nodes`, `n_edges`: the numbers of features and of edges.
#' * `density`: `n_edges` over the number of pairs.
#' * `positive_fraction`: the share of the associations that passed the
#'   filters that are positive. With `'positive'`, the negative ones count
#'   though they make no edges (as NetCoMi's positive edge percentage with
#'   `"signedPos"`). `NA` without edges.
#' * `n_components`: connected components, counting single features.
#' * `lcc_size`: nodes in the LCC.
#' * `average_dissimilarity`: the average dissimilarity of all features
#'   (see above).
#' * `average_path_length`: the mean shortest-path length over the pairs of
#'   connected features, divided by `average_dissimilarity` (NetCoMi's
#'   `sPathNorm = TRUE`); 1 without edges.
#' * `clustering_coefficient`: the mean of Barrat et al.'s (2004) weighted
#'   local clustering coefficient over the nodes with at least two
#'   neighbours; 0 when there are none.
#' * `modularity`: Newman's (2006) modularity of the clusters, unweighted,
#'   with cluster 0 as one group (as NetCoMi reports it); 0 with
#'   `clusters = 'none'` or without edges.
#'
#' **Size.** Shortest paths from every node take time
#' \eqn{O(N(E + N \log N))} for \eqn{E} edges, split over `cpus` threads;
#' memory is \eqn{O(N + E)} per thread. 3,000 features with 168,000 edges
#' take about 4 seconds on 4 threads.
#'
#' **Names in NetCoMi.** `threshold` is `thresh` (with `sparsMethod =
#' "threshold"`), `alpha` and `p_adjust` are `alpha` and `adjust` (with
#' `sparsMethod = "t-test"` or a permutation test), `sign` is `dissFunc`
#' (`"signed"`, `"unsigned"`, `"signedPos"`), `clusters = 'greedy'` is
#' `clustMethod = "cluster_fast_greedy"`, `hubs` and `hub_quantile` are
#' `hubPar` and `hubQuant`, and `weighted_degree` is `weightDeg`. NetCoMi
#' calls the global properties `nComp`, `avDiss`, `avPath`, `clustCoef`,
#' `modularity`, `density` and `pep` (a percentage).
#'
#' **Differences from NetCoMi** (1.2.0, whose outputs were compared on 160
#' networks to a relative \eqn{10^{-13}}):
#'
#' * p-values are adjusted with [stats::p.adjust()] only (NetCoMi's default
#'   is the local false discovery rate).
#' * Hubs use the empirical quantile; NetCoMi's log-normal fit
#'   (`lnormFit = TRUE`) isn't offered.
#' * NetCoMi switches to the weighted degree when every pair is an edge;
#'   `weighted_degree` is never changed here.
#' * A network without edges is an error in NetCoMi; here every centrality
#'   is 0.
#' * Two values that tie but for rounding are never one above the other
#'   (e.g. two eigenvector centralities of 1); NetCoMi's hubs can depend on
#'   that rounding.
#' * Only the global properties above, of the whole network: not natural
#'   connectivity, vertex or edge connectivity, graphlet measures, or the
#'   LCC's own properties.
#' * `plot()` draws a circular layout, with colour-blind-safer edge
#'   colours.
#'
#' @references
#' Peschel, S., Müller, C. L., von Mutius, E., Boulesteix, A.-L., & Depner,
#' M. (2021). NetCoMi: network construction and comparison for microbiome
#' data in R. *Briefings in Bioinformatics*, 22, bbaa290.
#'
#' Brandes, U. (2001). A faster algorithm for betweenness centrality.
#' *Journal of Mathematical Sociology*, 25, 163–177.
#'
#' Clauset, A., Newman, M. E. J., & Moore, C. (2004). Finding community
#' structure in very large networks. *Physical Review E*, 70, 066111.
#'
#' Barrat, A., Barthélemy, M., Pastor-Satorras, R., & Vespignani, A.
#' (2004). The architecture of complex weighted networks. *PNAS*, 101,
#' 3747–3752.
#'
#' Newman, M. E. J. (2006). Modularity and community structure in networks.
#' *PNAS*, 103, 8577–8582.
#'
#' @seealso [assoc_network_compare()], [assoc_network_diff()],
#'   [assoc_pairs()], [sparcc()], [rho_p()].
#'
#' @export
#' @examples
#'     # 60 samples of 20 taxa, with two groups of correlated taxa
#'     set.seed(1)
#'     basis <- matrix(rnorm(60 * 20), 60, 20)
#'     basis[, 2:5]   <- basis[, 2:5]   + basis[, 1]
#'     basis[, 12:14] <- basis[, 12:14] + basis[, 11]
#'     counts <- matrix(
#'       rpois(1200, exp(basis + 4)), 60, 20,
#'       dimnames = list(NULL, paste0('OTU', 1:20)) )
#'
#'     net <- assoc_network(sparcc(counts, n_draws = 20))
#'     net
#'     net$global
#'     plot(net)
#'
assoc_network <- function (
    x, threshold = NULL, alpha = NULL, p_adjust = 'BH', sign = 'signed',
    clusters = 'greedy', hubs = 'eigenvector', hub_quantile = 0.95,
    weighted_degree = FALSE, cpus = n_cpus() ) {

  #--------------------------------------------------------
  # Arguments.
  #--------------------------------------------------------
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0))
    stop('`threshold` must be NULL or a single non-negative number.', call. = FALSE)

  if (!is.null(alpha) && (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha > 1))
    stop('`alpha` must be NULL or a single number above 0 and at most 1.', call. = FALSE)

  p_adjust <- network_p_adjust(p_adjust)
  sign     <- network_choice(sign, c('signed', 'unsigned', 'positive'), 'sign')
  clusters <- network_choice(clusters, c('greedy', 'none'), 'clusters')

  cents <- c('degree', 'betweenness', 'closeness', 'eigenvector')
  if (!is.character(hubs) || length(hubs) < 1L || anyNA(hubs) || !all(hubs %in% cents))
    stop(
      '`hubs` must name one or more of ',
      paste0('"', cents, '"', collapse = ', '), '.', call. = FALSE )
  hubs <- unique(hubs)

  if (!is.numeric(hub_quantile) || length(hub_quantile) != 1L || is.na(hub_quantile) ||
      hub_quantile < 0 || hub_quantile > 1)
    stop('`hub_quantile` must be a single number from 0 to 1.', call. = FALSE)

  if (!is.logical(weighted_degree) || length(weighted_degree) != 1L || is.na(weighted_degree))
    stop('`weighted_degree` must be TRUE or FALSE.', call. = FALSE)

  validate_cpus()

  a <- network_input(x)

  if (!is.null(alpha) && is.null(a$p_value))
    stop(
      '`alpha` needs p-values, and `x` has none. ',
      'Run sparcc() with `permutations` to get them.', call. = FALSE )


  #--------------------------------------------------------
  # The cutoff: given, the method's own, or 0.3.
  #--------------------------------------------------------
  strict <- FALSE
  if (is.null(threshold)) {
    info <- feature_assoc_info(a$method)
    own  <- attr(x, 'threshold', exact = TRUE)
    if (!is.null(info) && isTRUE(info$sparse)) {
      threshold <- 0
      strict    <- TRUE
    }
    else if (!is.null(info) && is.numeric(own) && length(own) == 1L && is.finite(own)) {
      threshold <- own
      strict    <- isTRUE(info$strict)
    }
    else threshold <- 0.3
  }
  threshold <- as.double(threshold)


  #--------------------------------------------------------
  # Edges.
  #--------------------------------------------------------
  sel <- assoc_select(a, threshold, alpha, p_adjust, strict = strict)$pairs
  sel <- sel[sel$association != 0, , drop = FALSE]

  # The share of positive associations among those that pass, counted
  # before `positive` drops the negative ones (as NetCoMi's "signedPos").
  pos_frac <- if (nrow(sel)) mean(sel$association > 0) else NA_real_
  if (sign == 'positive') sel <- sel[sel$association > 0, , drop = FALSE]
  sel <- sel[order(sel$pos_from, sel$pos_to), , drop = FALSE]

  r <- pmin(pmax(sel$association, -1), 1)
  d <- if (sign == 'unsigned') sqrt(1 - r^2) else sqrt(0.5 * (1 - r))

  # An edge of weight 0 (r = -1 with 'signed') is no edge, as in NetCoMi,
  # whose adjacency matrix holds the weights.
  w0 <- d >= 1
  if (any(w0)) {
    sel <- sel[!w0, , drop = FALSE]
    r   <- r[!w0]
    d   <- d[!w0]
  }

  edges <- data.frame(
    from        = sel$from,
    to          = sel$to,
    association = sel$association,
    similarity  = 1 - d,
    distance    = d,
    sign        = ifelse(sel$association > 0, 1L, -1L),
    stringsAsFactors = FALSE )
  if (!is.null(a$p_value)) {
    edges$p_value    <- sel$p_value
    edges$p_adjusted <- sel$p_adjusted
  }
  rownames(edges) <- NULL


  #--------------------------------------------------------
  # The analysis.
  #--------------------------------------------------------
  labels <- character(a$n)
  labels[a$position] <- a$labels

  res <- network_analyze(
    n = a$n, from = sel$pos_from, to = sel$pos_to, r = r, d = d,
    clusters = clusters, hubs = hubs, hub_quantile = hub_quantile,
    weighted_degree = weighted_degree, cpus = cpus, positive_fraction = pos_frac )

  nodes <- data.frame(feature = labels, res$nodes, stringsAsFactors = FALSE)
  rownames(nodes) <- NULL

  net <- structure(
    list(edges = edges, nodes = nodes, global = res$global),
    class = 'ecodive_network' )
  if (!is.null(a$method)) attr(net, 'method') <- a$method
  ns <- attr(x, 'n_samples', exact = TRUE)
  if (is.numeric(ns) && length(ns) == 1L) attr(net, 'n_samples') <- ns
  attr(net, 'threshold') <- threshold
  if (!is.null(alpha)) {
    attr(net, 'alpha')    <- as.double(alpha)
    attr(net, 'p_adjust') <- p_adjust
  }
  attr(net, 'sign')         <- sign
  attr(net, 'hubs')         <- hubs
  attr(net, 'hub_quantile') <- as.double(hub_quantile)

  from <- if (inherits(x, 'dist')) display_title(x)
  set_display(net, title = title_of('Association Network', from))
}



# The associations of `x` for assoc_network(): as assoc_pairs_input(), with
# values beyond the range of a correlation refused.
network_input <- function (x) {

  if (inherits(x, 'dist')) {
    info <- feature_assoc_info(dist_metric_of(x))
    if (!is.null(info) && identical(info$kind, 'dissimilarity'))
      stop(
        '`x` is from ', info$id, '(), a dissimilarity, not an association ',
        'on a correlation scale. Use rho_p(), or convert it to a similarity ',
        'matrix with 1 on the diagonal, such as 1 / (1 + as.matrix(x)).', call. = FALSE )
  }

  a <- assoc_pairs_input(x)

  if (is.null(a$method)) {
    v <- a$values[!is.na(a$values)]
    if (any(abs(v) > 1 + 1e-8))
      stop('`x` must have values from -1 to 1.', call. = FALSE)
  }
  a
}



# The graph analysis of assoc_network(): n nodes, edges from[e] -- to[e]
# (1-based input positions) with associations r and lengths d. Returns
# list(nodes = data frame without `feature`, global = named vector).
network_analyze <- function (
    n, from, to, r, d, clusters, hubs, hub_quantile, weighted_degree, cpus,
    positive_fraction ) {

  m  <- length(from)
  s  <- 1 - d
  f0 <- as.integer(from - 1L)
  t0 <- as.integer(to - 1L)

  # Connected components (ids: each component's first node), and the
  # largest, LCC (ties: the one with the first node).
  comp     <- network_components(n, from, to)
  size     <- tabulate(comp, nbins = n)[comp]
  lcc_size <- max(size)
  in_lcc   <- comp == comp[which(size == lcc_size)[1L]]
  ein      <- in_lcc[from]
  k        <- lcc_size

  # Average dissimilarity: the mean over all pairs, with 1 for the pairs
  # that aren't edges (NetCoMi's avDiss, avDissIgnoreInf = FALSE); of the
  # whole network, and of the LCC.
  pairs  <- n * (n - 1) / 2
  av_all <- if (n > 1L) (sum(d) + (pairs - m)) / pairs else NA_real_
  kp     <- k * (k - 1) / 2
  av_lcc <- if (k > 1L) (sum(d[ein]) + (kp - sum(ein))) / kp else NA_real_

  deg <- bet <- clo <- eig <- numeric(n)
  apl <- 1   # NetCoMi's value for a network without edges

  if (m > 0L) {

    # Shortest paths on edge lengths d, from every node. An edge of length
    # 0 (r = 1) is left out, as in NetCoMi, whose igraph graph of
    # dissimilarities reads a 0 as "no edge".
    pe  <- d > 0
    p   <- .Call(C_network_paths, n, f0[pe], t0[pe], as.double(d[pe]), cpus)

    # Average path length: the mean shortest path between connected nodes,
    # in units of the average dissimilarity.
    if (sum(p[[4]]) > 0) apl <- (sum(p[[3]]) / sum(p[[4]])) / av_all

    # Centralities, on the LCC. Degree is normalised by the number of other
    # nodes in the whole network (as igraph's degree(normalized = TRUE) of
    # the whole network); weighted degree, the sum of the edge weights, is
    # not normalised (as igraph's strength()). Betweenness and closeness
    # are normalised by the LCC's size.
    lcc <- which(in_lcc)
    deg[lcc] <- if (weighted_degree)
      vapply(split(c(s, s), factor(c(from, to), levels = seq_len(n))), sum, numeric(1))[lcc]
    else
      tabulate(c(from, to), nbins = n)[lcc] / (n - 1)

    bet[lcc] <- if (k > 2L) p[[1]][lcc] * 2 / ((k - 1) * (k - 2)) else 0

    # Closeness from shortest paths normalised by the LCC's average
    # dissimilarity (NetCoMi's sPathNorm = TRUE).
    clo[lcc] <- av_lcc * p[[2]][lcc] / (k - 1)

    idx <- match(seq_len(n), lcc)
    eig[lcc] <- as.vector(.Call(
      C_network_eigen, k, idx[from[ein]] - 1L, idx[to[ein]] - 1L, as.double(s[ein]), 1L ))
  }

  # Hubs: above the quantile of every named centrality, over all nodes
  # (those outside the LCC count, with 0).
  cent <- list(degree = deg, betweenness = bet, closeness = clo, eigenvector = eig)
  hub  <- rep(m > 0L, n)
  for (h in hubs) {
    v   <- cent[[h]]
    hub <- hub & above_quantile(v, hub_quantile)
  }

  # Clusters: numbered in order of their first feature; 0 for features
  # alone in their cluster, including those without edges (as NetCoMi).
  has_edge <- tabulate(c(from, to), nbins = n) > 0L
  cluster  <- integer(n)
  if (clusters == 'greedy' && m > 0L) {
    mem  <- .Call(C_network_greedy, n, f0, t0, as.double(s))[[1]]
    many <- tabulate(mem + 1L, nbins = n)[mem + 1L] > 1L
    cluster[many] <- match(mem[many], unique(mem[many]))
  }

  # Modularity of the clusters, unweighted, with cluster 0 as one group
  # (as NetCoMi reports it).
  modularity <- if (m == 0L || clusters == 'none') 0 else {
    dg0 <- tabulate(c(from, to), nbins = n)
    cl  <- cluster
    tot <- vapply(split(dg0, cl), sum, numeric(1))
    sum(cl[from] == cl[to]) / m - sum((tot / (2 * m))^2)
  }

  # Barrat's weighted clustering coefficient, averaged over the nodes
  # with at least two neighbours (0 when there are none, as NetCoMi).
  cc <- if (m == 0L) 0 else {
    loc <- .Call(C_network_barrat, n, f0, t0, as.double(s))
    if (all(is.nan(loc))) 0 else mean(loc[!is.nan(loc)])
  }

  global <- c(
    n_nodes                = n,
    n_edges                = m,
    density                = if (n > 1L) m / pairs else NA_real_,
    positive_fraction      = positive_fraction,
    n_components           = sum(comp == seq_len(n)),
    lcc_size               = lcc_size,
    average_dissimilarity  = av_all,
    average_path_length    = apl,
    clustering_coefficient = cc,
    modularity             = modularity )

  list(
    nodes = data.frame(
      degree      = deg,
      betweenness = bet,
      closeness   = clo,
      eigenvector = eig,
      cluster     = cluster,
      hub         = hub,
      in_lcc      = in_lcc ),
    global = global )
}



# Which of `v` are above its quantile `q` (type 7) by more than rounding
# (a relative 1e-12), so that values equal but for rounding, such as two
# nodes' eigenvector centralities of 1, are never above one another.
above_quantile <- function (v, q) {
  lim <- stats::quantile(v, q, names = FALSE, type = 7)
  v - lim > 1e-12 * max(abs(v), 1e-300)
}



# Each node's component, as the smallest node in it (union-find).
network_components <- function (n, from, to) {
  parent <- seq_len(n)
  find <- function (v) {
    while (parent[v] != v) v <- parent[v]
    v
  }
  for (e in seq_along(from)) {
    a <- find(from[e])
    b <- find(to[e])
    if (a < b) parent[b] <- a
    else if (b < a) parent[a] <- b
    # Shorten the paths walked, so later finds are quick.
    parent[from[e]] <- parent[to[e]] <- min(a, b)
  }
  for (v in seq_len(n)) parent[v] <- parent[parent[v]]
  parent
}



network_p_adjust <- function (p_adjust) {
  if (!is.character(p_adjust) || length(p_adjust) != 1L || !p_adjust %in% stats::p.adjust.methods)
    stop(
      '`p_adjust` must be one of ',
      paste0('"', stats::p.adjust.methods, '"', collapse = ', '), '.', call. = FALSE )
  p_adjust
}

network_choice <- function (x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices)
    stop(
      '`', name, '` must be one of ',
      paste0('"', choices, '"', collapse = ', '), '.', call. = FALSE )
  x
}

#' @export
as.data.frame.ecodive_network <- function (x, ...) x$nodes



#' @export
print.ecodive_network <- function (x, ...) {

  g     <- x$global
  nodes <- x$nodes
  title <- plot_title(x, 'Association Network')

  cat(title, '\n', sep = '')
  cat(sprintf(
    '%d features, %d edges (%s positive) with |r| >= %s\n',
    as.integer(g[['n_nodes']]), as.integer(g[['n_edges']]),
    if (is.na(g[['positive_fraction']])) 'none'
    else paste0(round(100 * g[['positive_fraction']]), '%'),
    format(attr(x, 'threshold', exact = TRUE)) ))
  a <- attr(x, 'alpha', exact = TRUE)
  if (!is.null(a))
    cat(sprintf('  and %s-adjusted p <= %s\n', attr(x, 'p_adjust', exact = TRUE), format(a)))
  cat(sprintf(
    '%d components; the largest has %d features\n',
    as.integer(g[['n_components']]), as.integer(g[['lcc_size']]) ))
  k <- length(unique(nodes$cluster[nodes$cluster > 0L]))
  cat(sprintf('%d clusters, modularity %s\n', k, plot_num(g[['modularity']])))
  hubs <- nodes$feature[nodes$hub]
  cat('Hubs (', paste(attr(x, 'hubs', exact = TRUE), collapse = ', '), '): ',
      if (length(hubs)) paste(hubs, collapse = ', ') else 'none', '\n', sep = '')
  invisible(x)
}
