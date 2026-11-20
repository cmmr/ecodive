#test_that("assoc_network", {

  # assoc_network(), assoc_network_compare() and assoc_network_diff()
  # against NetCoMi 1.2.0 and igraph outputs (data-raw/expected/
  # assoc_network.r), against base-R implementations of the formulas, and
  # their edge cases, determinism and argument checks.

  ref <- readRDS(file.path('expected', 'assoc_network.rds'))

  same_partition <- function (a, b) {
    t <- table(a, b)
    all(rowSums(t > 0) == 1) && all(colSums(t > 0) == 1)
  }

  check_ref <- function (net, r, label, hubs = TRUE, degree = TRUE) {
    nd <- net$nodes
    if (degree) expect_equal(nd$degree, r$degree, tolerance = 1e-10, info = label)
    expect_equal(nd$betweenness, r$betweenness, tolerance = 1e-10, info = label)
    expect_equal(nd$closeness,   r$closeness,   tolerance = 1e-10, info = label)
    expect_equal(nd$eigenvector, r$eigenvector, tolerance = 1e-10, info = label)
    expect_true(same_partition(nd$cluster, r$cluster), info = label)
    expect_identical(nd$cluster == 0L, r$cluster == 0, info = label)
    if (hubs) expect_identical(sort(nd$feature[nd$hub]), sort(r$hubs), info = label)
    g <- net$global
    for (p in names(r$global))
      expect_equal(g[[p]], r$global[[p]], tolerance = 1e-10, info = paste(label, p))
  }


  #----------------------------------------------------------------
  # NetCoMi: 12 pairs of networks with varied settings.
  #----------------------------------------------------------------
  for (i in seq_along(ref$networks)) {
    cs <- ref$networks[[i]]
    nx <- do.call(assoc_network, c(list(cs$A1), cs$args))
    ny <- do.call(assoc_network, c(list(cs$A2), cs$args))
    check_ref(nx, cs$ref1, paste('network', i, 'x'))
    check_ref(ny, cs$ref2, paste('network', i, 'y'))

    cmp <- assoc_network_compare(nx, ny, permutations = 0L)
    j   <- cs$jaccard
    expect_equal(cmp$jaccard$jaccard,   unname(j[, 'jacc']),      tolerance = 1e-12, info = paste('jaccard', i))
    expect_equal(cmp$jaccard$p_less,    unname(j[, 'p.less']),    tolerance = 1e-12, info = paste('p_less', i))
    expect_equal(cmp$jaccard$p_greater, unname(j[, 'p.greater']), tolerance = 1e-12, info = paste('p_greater', i))
    expect_equal(cmp$ari[['value']], cs$ari, tolerance = 1e-12, info = paste('ari', i))
    expect_true(is.na(cmp$ari[['p_value']]))
  }

  # A pair with r = 1 (left out of shortest paths, as NetCoMi).
  net <- assoc_network(ref$edge_r1$A)
  check_ref(net, ref$edge_r1$ref, 'r = 1')
  expect_equal(net$nodes$eigenvector[1], net$nodes$eigenvector[2], tolerance = 1e-12)

  # A complete network: NetCoMi switches to weighted degree.
  net <- assoc_network(ref$complete$A, threshold = 0)
  check_ref(net, ref$complete$ref, 'complete', degree = FALSE)
  net <- assoc_network(ref$complete$A, threshold = 0, weighted_degree = TRUE)
  expect_equal(net$nodes$degree, ref$complete$ref$degree, tolerance = 1e-10)

  # A single edge.
  check_ref(assoc_network(ref$one_edge$A), ref$one_edge$ref, 'one edge')

  # NetCoMi's Jaccard probabilities, for sets of every size up to 6.
  g <- ref$jacc_grid
  for (r in seq_len(nrow(g))) {
    a <- seq_len(30) <= g[r, 'a']
    b <- seq_len(30) > 30 - g[r, 'b'] + g[r, 'k'] | seq_len(30) <= g[r, 'k']
    jr <- network_jaccard(a, b)
    expect_equal(c(jr$jaccard, jr$p_greater, jr$p_less), unname(g[r, 4:6]), tolerance = 1e-12)
  }

  # Fisher's z: NetCoMi computes 1 - pnorm(), which loses digits below
  # 1e-11; compared absolutely.
  d  <- ref$diff
  dd <- assoc_network_diff(d$A1, d$A2, n_samples = c(d$n1, d$n2))
  ij <- cbind(match(dd$feature1, rownames(d$p_value)), match(dd$feature2, rownames(d$p_value)))
  expect_equal(dd$p_value,    d$p_value[ij],    tolerance = 1e-14, scale = 1)
  expect_equal(dd$p_adjusted, d$p_adjusted[ij], tolerance = 1e-14, scale = 1)


  #----------------------------------------------------------------
  # igraph on hand-built graphs.
  #----------------------------------------------------------------
  for (nm in names(ref$igraph)) {
    h  <- ref$igraph[[nm]]
    n  <- max(h$edges)
    f0 <- as.integer(h$edges[, 1] - 1L)
    t0 <- as.integer(h$edges[, 2] - 1L)
    p  <- .Call(C_network_paths, n, f0, t0, 1 - h$s, 1L)
    expect_equal(p[[1]], h$betweenness, tolerance = 1e-10, info = nm)
    expect_equal(p[[2]], h$harmonic,    tolerance = 1e-10, info = nm)
    expect_equal(as.vector(.Call(C_network_eigen, n, f0, t0, h$s, 1L)), h$eigenvector, tolerance = 1e-10, info = nm)
    expect_equal(.Call(C_network_barrat, n, f0, t0, h$s), h$barrat, tolerance = 1e-12, info = nm)
    # The star's partitions tie but for rounding (Q = 0 and -7e-17): only
    # its modularity is compared.
    gr <- .Call(C_network_greedy, n, f0, t0, h$s)
    if (nm != 'star') expect_true(same_partition(gr[[1]], h$fg_membership), info = nm)
    expect_equal(max(gr[[2]]), h$fg_modularity, tolerance = 1e-10, info = nm)
  }


  #----------------------------------------------------------------
  # Base-R references: Floyd-Warshall paths, brute-force betweenness,
  # eigen(), Barrat's formula, Newman's modularity.
  #----------------------------------------------------------------
  base_network <- function (A, threshold) {
    n  <- nrow(A)
    E  <- abs(A) >= threshold & row(A) != col(A)
    D  <- ifelse(E, sqrt(0.5 * (1 - A)), Inf); diag(D) <- 0
    S  <- ifelse(E, 1 - sqrt(0.5 * (1 - A)), 0)
    for (k in 1:n) D <- pmin(D, outer(D[, k], D[k, ], '+'))
    comp <- apply(is.finite(D), 1, function (r) min(which(r)))
    lcc  <- comp == as.integer(names(which.max(table(factor(comp, unique(comp))))))
    k    <- sum(lcc)
    avl  <- (sum(ifelse(E, sqrt(0.5 * (1 - A)), 1)[lcc, lcc][upper.tri(diag(k))])) / (k * (k - 1) / 2)
    # Betweenness by counting shortest paths (no ties in random data).
    np <- matrix(0, n, n)
    for (s in 1:n) for (t in 1:n) np[s, t] <- if (s == t) 1 else NA
    bet <- numeric(n)
    L <- which(lcc)
    for (s in L) for (t in L) if (s < t) for (v in L) if (v != s && v != t)
      if (abs(D[s, v] + D[v, t] - D[s, t]) < 1e-12) bet[v] <- bet[v] + 1
    harm <- rowSums(ifelse(is.finite(D) & D > 0, 1 / D, 0))
    ev   <- abs(eigen(S[lcc, lcc], symmetric = TRUE)$vectors[, 1])
    eig  <- numeric(n); eig[lcc] <- ev / max(ev)
    list(
      degree      = ifelse(lcc, rowSums(E) / (n - 1), 0),
      betweenness = ifelse(lcc, bet * 2 / ((k - 1) * (k - 2)), 0),
      closeness   = ifelse(lcc, avl * harm / (k - 1), 0),
      eigenvector = eig,
      barrat      = vapply(1:n, function (i) {
        nb <- which(E[i, ]); if (length(nb) < 2) return (NaN)
        tot <- 0
        for (j in nb) for (h in nb) if (j != h && E[j, h]) tot <- tot + (S[i, j] + S[i, h]) / 2
        tot / (sum(S[i, ]) * (length(nb) - 1))
      }, numeric(1)),
      E = E )
  }

  set.seed(2)
  for (rep in 1:3) {
    b <- matrix(rnorm(40 * 14), 40, 14)
    b[, 2:4] <- b[, 2:4] + b[, 1]
    b[, 8:9] <- b[, 8:9] - b[, 7]
    A <- cor(b)
    net <- assoc_network(A, threshold = 0.25, clusters = 'greedy')
    br  <- base_network(A, 0.25)
    expect_equal(net$nodes$degree,      br$degree,      tolerance = 1e-12)
    expect_equal(net$nodes$betweenness, br$betweenness, tolerance = 1e-12)
    expect_equal(net$nodes$closeness,   br$closeness,   tolerance = 1e-12)
    expect_equal(net$nodes$eigenvector, br$eigenvector, tolerance = 1e-10)
    cc <- br$barrat
    expect_equal(net$global[['clustering_coefficient']], mean(cc[!is.nan(cc)]), tolerance = 1e-12)
    # Newman's modularity (unweighted) of the clusters found.
    E  <- br$E; m <- sum(E) / 2; kd <- rowSums(E)
    cl <- ifelse(net$nodes$cluster > 0, net$nodes$cluster, -seq_len(nrow(A)))
    Q  <- sum((E - outer(kd, kd) / (2 * m)) * outer(cl, cl, '==')) / (2 * m)
    expect_equal(net$global[['modularity']], Q, tolerance = 1e-12)
    # Greedy merges never leave a better neighbouring merge... at least the
    # clusters are connected and the modularity is positive.
    expect_true(Q > 0)
    # Edges.
    ex <- net$edges
    expect_identical(nrow(ex), as.integer(sum(br$E) / 2))
    expect_identical(ex$association, A[cbind(as.integer(ex$from), as.integer(ex$to))])
    expect_equal(ex$distance, sqrt(0.5 * (1 - ex$association)), tolerance = 1e-15)
    expect_equal(ex$similarity, 1 - ex$distance, tolerance = 1e-15)
    expect_identical(ex$sign, ifelse(ex$association > 0, 1L, -1L))
  }

  # Unsigned and positive lengths.
  A   <- cor(matrix(rnorm(300), 30, 10)) * 3; diag(A) <- 1; A <- pmin(pmax(A, -1), 1)
  un  <- assoc_network(A, sign = 'unsigned')
  expect_equal(un$edges$distance, sqrt(1 - un$edges$association^2), tolerance = 1e-15)
  ps  <- assoc_network(A, sign = 'positive')
  expect_true(all(ps$edges$association > 0))
  expect_identical(ps$global[['positive_fraction']], mean(un$edges$association > 0))


  #----------------------------------------------------------------
  # Hand-built graphs with known answers.
  #----------------------------------------------------------------
  graph_matrix <- function (n, edges, r = 0.8) {
    A <- diag(n)
    A[edges] <- A[edges[, 2:1]] <- r
    dimnames(A) <- list(paste0('n', 1:n), paste0('n', 1:n))
    A
  }
  # Star: the centre lies on every path.
  star <- assoc_network(graph_matrix(5, cbind(1, 2:5)))
  expect_equal(star$nodes$betweenness, c(1, 0, 0, 0, 0))
  expect_equal(star$nodes$degree, c(1, 0.25, 0.25, 0.25, 0.25))
  expect_identical(star$nodes$hub, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(star$global[['clustering_coefficient']], 0)
  # Path of 5: betweenness 0, 3/6, 4/6, 3/6, 0.
  path <- assoc_network(graph_matrix(5, cbind(1:4, 2:5)))
  expect_equal(path$nodes$betweenness, c(0, 3, 4, 3, 0) / 6)
  # Ring: all equal.
  ring <- assoc_network(graph_matrix(6, cbind(1:6, c(2:6, 1))))
  expect_equal(diff(range(ring$nodes$betweenness)), 0)
  expect_equal(ring$nodes$eigenvector, rep(1, 6))
  expect_false(any(ring$nodes$hub))
  # Two cliques joined by one edge: two clusters.
  cl <- assoc_network(graph_matrix(8, rbind(t(combn(1:4, 2)), t(combn(5:8, 2)), c(4, 5))))
  expect_identical(cl$nodes$cluster, rep(1:2, each = 4))
  expect_equal(cl$global[['modularity']], 2 * (6 / 13 - (13 / 26)^2))
  # Disconnected graph with isolated nodes: LCC only.
  dg <- assoc_network(graph_matrix(7, rbind(c(1, 2), c(2, 3), c(5, 6))))
  expect_identical(dg$nodes$in_lcc, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(dg$nodes$cluster[c(4, 7)], c(0L, 0L))
  expect_equal(dg$nodes$degree[5:6], c(0, 0))
  expect_identical(unname(dg$global[c('n_components', 'lcc_size', 'n_edges')]), c(4, 3, 3))
  # A tie for the largest component goes to the first.
  tie <- assoc_network(graph_matrix(4, rbind(c(3, 4), c(1, 2))))
  expect_identical(tie$nodes$in_lcc, c(TRUE, TRUE, FALSE, FALSE))


  #----------------------------------------------------------------
  # Edge cases.
  #----------------------------------------------------------------
  # No edges.
  none <- assoc_network(diag(5), threshold = 0.5)
  expect_identical(nrow(none$edges), 0L)
  expect_true(all(unlist(none$nodes[c('degree', 'betweenness', 'closeness', 'eigenvector', 'cluster')]) == 0))
  expect_false(any(none$nodes$hub))
  expect_equal(unname(none$global[c('modularity', 'clustering_coefficient', 'average_path_length',
                                    'average_dissimilarity', 'n_components', 'lcc_size')]), c(0, 0, 1, 1, 5, 1))
  expect_true(is.na(none$global[['positive_fraction']]))
  # Every pair an edge.
  full <- assoc_network(matrix(0.5, 4, 4) + diag(0.5, 4), threshold = 0.1)
  expect_identical(nrow(full$edges), 6L)
  expect_equal(full$global[['density']], 1)
  expect_equal(full$nodes$betweenness, rep(0, 4))
  # |r| = 1.
  one <- assoc_network(matrix(c(1, 1, 1, 1), 2, 2, dimnames = list(c('a', 'b'), c('a', 'b'))))
  expect_identical(one$edges$distance, 0)
  expect_identical(one$nodes$closeness, c(0, 0))
  # Two features.
  two <- assoc_network(matrix(c(1, -0.6, -0.6, 1), 2))
  expect_identical(two$edges$from, '1')
  expect_identical(two$edges$sign, -1L)
  expect_identical(two$nodes$eigenvector, c(1, 1))
  # clusters = 'none'.
  nc <- assoc_network(ref$networks[[2]]$A1, clusters = 'none')
  expect_true(all(nc$nodes$cluster == 0L))
  expect_identical(nc$global[['modularity']], 0)


  #----------------------------------------------------------------
  # Results from sparcc() and rho_p(): attributes, titles, clustering.
  #----------------------------------------------------------------
  set.seed(1)
  basis <- matrix(rnorm(50 * 12), 50, 12)
  basis[, 2:4] <- basis[, 2:4] + basis[, 1]
  counts <- matrix(rpois(600, exp(basis + 4)), 50, 12, dimnames = list(NULL, paste0('f', 1:12)))
  sp  <- sparcc(counts, n_draws = 10L, permutations = 19L, cpus = 2L)
  net <- assoc_network(sp)
  expect_identical(class(net), 'ecodive_network')
  expect_identical(names(net), c('edges', 'nodes', 'global'))
  expect_identical(attr(net, 'method'), 'sparcc')
  expect_identical(attr(net, 'n_samples'), 50L)
  expect_identical(attr(net, 'threshold'), 0.3)
  expect_identical(attr(net, 'display')$title, 'Association Network of SparCC Correlation')
  expect_identical(names(net$edges), c('from', 'to', 'association', 'similarity', 'distance', 'sign', 'p_value', 'p_adjusted'))
  expect_identical(net$nodes$feature, colnames(counts))
  expect_identical(as.data.frame(net), net$nodes)
  pa <- assoc_pairs(sp)
  expect_identical(nrow(net$edges), nrow(pa))
  expect_true(setequal(paste(net$edges$from, net$edges$to), paste(pa$from, pa$to)))
  expect_false(is.unsorted(match(net$edges$from, colnames(counts))))
  # A clustered result gives the same network.
  spc <- sparcc(counts, n_draws = 10L, permutations = 19L, cpus = 2L, clust = 'average')
  expect_identical(assoc_network(spc), net)
  # alpha.
  net_a <- assoc_network(sp, threshold = 0, alpha = 1)
  expect_identical(nrow(net_a$edges), 66L)
  expect_identical(attr(net_a, 'alpha'), 1)
  expect_identical(attr(net_a, 'p_adjust'), 'BH')
  # rho_p().
  rp <- rho_p(counts, pseudocount = 1)
  expect_identical(attr(assoc_network(rp), 'display')$title, 'Association Network of Proportionality (rho_p)')
  expect_identical(attr(assoc_network(as.matrix(rp) + diag(12)), 'display')$title, 'Association Network')
  # feature_assoc() gives the same.
  expect_identical(assoc_network(feature_assoc(counts, 'rho_p', pseudocount = 1)), assoc_network(rp))
  # A result's own sparsity: a method whose zeros mean "no edge".
  sp0 <- sp
  sp0[abs(sp0) < 0.2] <- 0
  local({
    old <- FEATURE_ASSOCS
    on.exit(assign('FEATURE_ASSOCS', old, envir = asNamespace('ecodive')))
    fa <- old; fa$sparse[fa$id == 'sparcc'] <- TRUE
    unlockBinding('FEATURE_ASSOCS', asNamespace('ecodive'))
    assign('FEATURE_ASSOCS', fa, envir = asNamespace('ecodive'))
    sn <- assoc_network(sp0)
    expect_identical(attr(sn, 'threshold'), 0)
    expect_identical(nrow(sn$edges), sum(sp0 != 0))
  })


  #----------------------------------------------------------------
  # Determinism: identical for any cpus.
  #----------------------------------------------------------------
  set.seed(3)
  big <- cor(matrix(rnorm(80 * 150), 80) %*% diag(150) + rnorm(80) %o% rnorm(150) * 0.4)
  dimnames(big) <- list(paste0('F', 1:150), paste0('F', 1:150))
  b1 <- assoc_network(big, threshold = 0.15, cpus = 1L)
  for (cpus in c(2L, 3L, 4L, 8L))
    expect_identical(assoc_network(big, threshold = 0.15, cpus = cpus), b1)
  expect_true(nrow(b1$edges) > 1000)


  #----------------------------------------------------------------
  # assoc_network_compare().
  #----------------------------------------------------------------
  nx <- assoc_network(ref$networks[[2]]$A1)
  ny <- assoc_network(ref$networks[[2]]$A2)
  cmp <- assoc_network_compare(nx, ny, permutations = 99L)
  expect_identical(class(cmp), 'ecodive_network_compare')
  expect_identical(names(cmp), c('nodes', 'global', 'jaccard', 'ari', 'clusters', 'edges_x', 'edges_y'))
  expect_identical(names(cmp$nodes), c('feature', 'centrality', 'x', 'y', 'difference'))
  expect_identical(cmp$nodes$x[cmp$nodes$centrality == 'closeness'], nx$nodes$closeness)
  expect_identical(cmp$global$difference, unname(nx$global - ny$global))
  expect_identical(cmp$clusters$cluster_y, ny$nodes$cluster)
  expect_identical(attr(cmp, 'display')$title, 'Comparison of Association Network')
  expect_identical(attr(cmp, 'permutations'), 99L)
  # The ARI test: |ARI| of permuted labels, as perm_p() on the internal
  # hook's permutations.
  pm  <- .Call(C_perm_indices, 15L, 99L, NULL, 0L)
  aperm <- apply(pm, 2, function (p) adjusted_rand(nx$nodes$cluster, ny$nodes$cluster[p]))
  expect_identical(cmp$ari[['p_value']], perm_p(cmp$ari[['value']], aperm, 'two.sided'))
  expect_identical(network_compare_core(nx, ny, permutations = 0L, ari_perms = pm)$ari, cmp$ari)
  # y's features in another order give the same.
  ny2 <- ny; ny2$nodes <- ny$nodes[15:1, ]
  expect_identical(assoc_network_compare(nx, ny2, permutations = 99L)$jaccard, cmp$jaccard)
  # ARI against base R's formula.
  ari_base <- function (a, b) {
    t <- table(a, b); c2 <- function (v) sum(choose(v, 2)); n <- length(a)
    e <- c2(rowSums(t)) * c2(colSums(t)) / choose(n, 2)
    (c2(t) - e) / ((c2(rowSums(t)) + c2(colSums(t))) / 2 - e)
  }
  expect_equal(cmp$ari[['value']], ari_base(nx$nodes$cluster, ny$nodes$cluster), tolerance = 1e-14)
  expect_identical(adjusted_rand(c(1, 1, 2, 2), c(5, 5, 7, 7)), 1)

  # Group-permutation tests, checked against the networks rebuilt by hand
  # for the same permutations.
  set.seed(4)
  mk <- function (linked) {
    b <- matrix(rnorm(30 * 8), 30, 8)
    b[, linked] <- b[, linked] + b[, linked[1]]
    matrix(rpois(240, exp(b + 3)), 30, 8)
  }
  cnt <- rbind(mk(1:3), mk(5:7))
  colnames(cnt) <- paste0('t', 1:8)
  grp <- rep(c('ctl', 'trt'), each = 30)
  ax  <- function (k) rho_p(cnt[k, ], pseudocount = 1)
  gx  <- assoc_network(ax(grp == 'ctl'), threshold = 0.4)
  gy  <- assoc_network(ax(grp == 'trt'), threshold = 0.4)
  gp  <- .Call(C_perm_indices, 60L, 19L, NULL, 5L)
  gc  <- network_compare_core(
    gx, gy, permutations = 19L, counts = cnt, groups = grp,
    assoc_args = list(pseudocount = 1), network_args = list(threshold = 0.4),
    seed = 5L, group_perms = gp )
  gr <- assoc_network_compare(
    gx, gy, permutations = 19L, counts = cnt, groups = grp,
    assoc_args = list(pseudocount = 1), network_args = list(threshold = 0.4), seed = 5L )
  expect_identical(gr, gc)
  stat <- function (a, b) c(
    unlist(a$nodes[c('degree', 'betweenness', 'closeness', 'eigenvector')], use.names = FALSE) -
      unlist(b$nodes[c('degree', 'betweenness', 'closeness', 'eigenvector')], use.names = FALSE),
    unname(a$global - b$global) )
  perm <- t(apply(gp, 2, function (p) {
    g <- grp[p]
    stat(assoc_network(ax(g == 'ctl'), threshold = 0.4), assoc_network(ax(g == 'trt'), threshold = 0.4))
  }))
  pv <- perm_p(stat(gx, gy), perm, 'two.sided')
  expect_identical(gc$nodes$p_value, pv[1:32])
  expect_identical(gc$global$p_value, pv[-(1:32)])
  expect_identical(
    gc$nodes$p_adjusted,
    unlist(lapply(split(pv[1:32], rep(1:4, each = 8)), p.adjust, 'BH'), use.names = FALSE) )
  expect_identical(attr(gc, 'groups'), c('ctl', 'trt'))
  # Groups by column name, a factor, and identical for any cpus.
  meta <- data.frame(g = factor(grp))
  for (cpus in c(1L, 3L))
    expect_identical(
      assoc_network_compare(
        gx, gy, permutations = 19L, counts = cnt, groups = 'g', data = meta,
        assoc_args = list(pseudocount = 1), network_args = list(threshold = 0.4),
        seed = 5L, cpus = cpus ),
      gr )
  # sparcc(): each permutation's own seed, identical for any cpus.
  set.seed(5)
  mk20 <- function (linked) {
    b <- matrix(rnorm(30 * 20), 30, 20)
    b[, linked] <- b[, linked] + b[, linked[1]]
    matrix(rpois(600, exp(b + 4)), 30, 20)
  }
  cnt20 <- rbind(mk20(1:4), mk20(5:8))
  colnames(cnt20) <- paste0('t', 1:20)
  rx <- sparcc(cnt20[grp == 'ctl', ], n_draws = 5L)
  ry <- sparcc(cnt20[grp == 'trt', ], n_draws = 5L)
  sx <- assoc_network(rx)
  sy <- assoc_network(ry)
  s1 <- assoc_network_compare(sx, sy, permutations = 3L, counts = cnt20, groups = grp,
                              assoc_args = list(n_draws = 5L), cpus = 1L)
  s4 <- assoc_network_compare(sx, sy, permutations = 3L, counts = cnt20, groups = grp,
                              assoc_args = list(n_draws = 5L), cpus = 4L)
  expect_identical(s1, s4)
  expect_true(all(s1$global$p_value >= 0.25))
  # Features as rows.
  expect_identical(
    assoc_network_compare(
      gx, gy, permutations = 19L, counts = t(cnt), groups = grp,
      assoc_args = list(pseudocount = 1, margin = 2L), network_args = list(threshold = 0.4), seed = 5L ),
    gr )


  #----------------------------------------------------------------
  # assoc_network_diff().
  #----------------------------------------------------------------
  dd <- assoc_network_diff(rx, ry)
  expect_identical(class(dd), c('ecodive_network_diff', 'data.frame'))
  expect_identical(names(dd), c('feature1', 'feature2', 'association_x', 'association_y', 'z', 'p_value', 'p_adjusted', 'reject'))
  expect_identical(nrow(dd), 190L)
  mx <- as.matrix(rx); my <- as.matrix(ry)
  ij <- cbind(dd$feature1, dd$feature2)
  z  <- (atanh(mx[ij]) - atanh(my[ij])) / sqrt(1 / 27 + 1 / 27)
  expect_equal(dd$z, z, tolerance = 1e-14)
  expect_equal(dd$p_value, 2 * pnorm(-abs(z)), tolerance = 1e-14)
  expect_identical(dd$p_adjusted, p.adjust(dd$p_value, 'BH'))
  expect_identical(dd$reject, dd$p_adjusted <= 0.05)
  expect_identical(attr(dd, 'display')$title, 'Differential Associations of SparCC Correlation')
  expect_identical(attr(dd, 'method'), 'sparcc')
  # y clustered (another order) gives the same.
  ryc <- sparcc(cnt20[grp == 'trt', ], n_draws = 5L, clust = 'average')
  expect_identical(assoc_network_diff(rx, ryc), dd)
  # n_samples overrides; |r| >= 1 gives NA with a warning.
  m1 <- as.matrix(rx); diag(m1) <- 1
  m2 <- as.matrix(ry); diag(m2) <- 1
  m1[1, 2] <- m1[2, 1] <- 1
  expect_warning(d2 <- assoc_network_diff(m1, m2, n_samples = c(30, 40)), 'correlation of 1')
  expect_true(is.na(d2$z[1]) && is.na(d2$p_value[1]))
  expect_equal(d2$z[2], (atanh(m1[1, 3]) - atanh(m2[1, 3])) / sqrt(1 / 27 + 1 / 37), tolerance = 1e-14)


  #----------------------------------------------------------------
  # Printing and plots.
  #----------------------------------------------------------------
  out <- capture.output(print(net))
  expect_identical(out[1], 'Association Network of SparCC Correlation')
  expect_true(any(grepl('^Hubs \\(eigenvector\\)', out)))
  expect_identical(capture.output(print(cmp))[1], 'Comparison of Association Network')
  expect_true(any(grepl('Differential Associations of SparCC Correlation', capture.output(print(dd)))))
  grDevices::pdf(NULL)
  expect_identical(plot(net), net)
  expect_identical(plot(net, curved = TRUE, isolated = FALSE, labels = FALSE, legend = NULL), net)
  expect_identical(plot(none), none)
  expect_identical(plot(gr), gr)
  expect_identical(plot(dd), dd)
  lay <- cbind(x = 1:12, y = 12:1); rownames(lay) <- colnames(counts)
  expect_identical(plot(net, layout = lay), net)
  grDevices::dev.off()
  pos <- network_plot_layout(NULL, nx$nodes$feature, nx$nodes$cluster, nx$nodes$degree,
                             c(nx$edges$from, nx$edges$to), TRUE)
  expect_equal(unname(pos[1, ]), c(0, 1))
  expect_identical(nrow(pos), 15L)
  linked <- nx$nodes$feature %in% c(nx$edges$from, nx$edges$to)
  expect_identical(rownames(pos)[seq_len(sum(linked))], rownames(pos)[rownames(pos) %in% nx$nodes$feature[linked]])


  #----------------------------------------------------------------
  # feature_assoc() and its registry.
  #----------------------------------------------------------------
  expect_identical(feature_assoc(counts, 'sparcc', n_draws = 5L, cpus = 2L),
                   sparcc(counts, n_draws = 5L, cpus = 2L))
  expect_identical(feature_assoc(counts, 'vlr', pseudocount = 1, seed = 3L), vlr(counts, pseudocount = 1))
  expect_identical(feature_assoc(t(counts), 'phi_s', margin = 2L, pseudocount = 1),
                   phi_s(counts, pseudocount = 1))
  expect_identical(list_feature_assocs('id'), c('vlr', 'rho_p', 'phi_s', 'sparcc'))
  expect_identical(list_feature_assocs('name')[[4]], 'SparCC Correlation')
  expect_identical(names(list_feature_assocs('func')), list_feature_assocs('id'))
  expect_identical(list_feature_assocs('list')$sparcc$fisher_z, TRUE)
  expect_identical(rownames(list_feature_assocs(nm = 'id')), list_feature_assocs('id'))
  expect_false('alt_ids' %in% names(list_feature_assocs()))
  for (nm in c('rho', 'propr', 'RhoP', 'rho_p'))      expect_identical(match_feature_assoc(nm)$id, 'rho_p', info = nm)
  for (nm in c('phi', 'phs', 'phi_s'))               expect_identical(match_feature_assoc(nm)$id, 'phi_s', info = nm)
  for (nm in c('SparCC', 'FastSpar', 'spar'))        expect_identical(match_feature_assoc(nm)$id, 'sparcc', info = nm)
  for (nm in c('vlr', 'log ratio variance', 'variance of log')) expect_identical(match_feature_assoc(nm)$id, 'vlr', info = nm)
  expect_identical(match_feature_assoc('sparcc')$params[1], 'counts')
  expect_error(match_feature_assoc('proportionality'), 'multiple')
  expect_error(match_feature_assoc('pearson'), 'does not match')
  expect_error(match_feature_assoc(1), 'character')
  expect_error(match_feature_assoc(c('a', 'b')), 'length 1')
  expect_error(match_feature_assoc(NA_character_), 'NA')
  expect_error(match_feature_assoc('--'), 'alphabetic')
  expect_null(feature_assoc_info(NULL))
  expect_null(feature_assoc_info('bray'))
  # The new aliases, and the assoc/net prefixes, are in no other registry.
  for (reg in c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'ZERO_REPLACEMENTS', 'DA_TESTS'))
    if (exists(reg))
      expect_identical(
        c(intersect(names(exact_names(reg)), c('rho', 'propr', 'rhop', 'phi', 'phs', 'phis', 'fastspar', 'logratiovariance')),
          grep('^(assoc|net)', names(exact_names(reg)), value = TRUE)),
        character(0), info = reg )


  #----------------------------------------------------------------
  # Argument checks.
  #----------------------------------------------------------------
  expect_error(assoc_network(vlr(counts + 1)), 'dissimilarity')
  expect_error(assoc_network(matrix(c(1, 2, 2, 1), 2)), 'from -1 to 1')
  expect_error(assoc_network(matrix(1, 1, 1)), 'correlations between features')
  expect_error(assoc_network(rp, alpha = 0.05), 'needs p-values')
  expect_error(assoc_network(rp, threshold = -1), '`threshold`')
  expect_error(assoc_network(rp, sign = 'both'), '`sign`')
  expect_error(assoc_network(rp, clusters = 'louvain'), '`clusters`')
  expect_error(assoc_network(rp, hubs = 'pagerank'), '`hubs`')
  expect_error(assoc_network(rp, hub_quantile = 2), '`hub_quantile`')
  expect_error(assoc_network(rp, weighted_degree = NA), '`weighted_degree`')
  expect_error(assoc_network(rp, p_adjust = 'xx'), '`p_adjust`')
  expect_error(assoc_network_compare(nx, rp), 'networks from assoc_network')
  expect_error(assoc_network_compare(nx, assoc_network(rp)), 'same method|same features')
  expect_error(assoc_network_compare(nx, ny, quantile = 2), '`quantile`')
  expect_error(assoc_network_compare(nx, ny, permutations = -1), '`permutations`')
  expect_error(assoc_network_compare(nx, ny, groups = grp), 'needs `counts`')
  expect_error(assoc_network_compare(gx, gy, counts = cnt), 'needs `groups`')
  expect_error(assoc_network_compare(nx, ny, counts = cnt, groups = grp), 'not matrices')
  expect_error(
    assoc_network_compare(gx, gy, counts = cnt, groups = rep(1:3, 20), permutations = 1L), 'exactly two')
  expect_error(assoc_network_diff(rp, rp), 'not Pearson-type')
  expect_error(assoc_network_diff(rx, m2), 'same method')
  expect_error(assoc_network_diff(m1, m2), 'n_samples')
  expect_error(assoc_network_diff(m1, m2[1:7, 1:7], n_samples = 30), 'same features')
  expect_error(assoc_network_diff(m1, m2, n_samples = c(1, 2, 3)), '`n_samples`')
  expect_error(assoc_network_diff(m1, m2, n_samples = 3), 'more than 3')
  expect_error(assoc_network_diff(rx, ry, alpha = 0), '`alpha`')

#})
