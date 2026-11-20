# Writes tests/tinytests/expected/assoc_network.rds: reference outputs for
# assoc_network(), assoc_network_compare() and assoc_network_diff() from
# NetCoMi and igraph (both GPL: run for their outputs only; no source of
# either was read).
#
# Run from the package root with NetCoMi 1.2.0 (the jsDelivr mirror of
# github.com/stefpeschel/NetCoMi@1.2.0) and igraph installed:
#
#     Rscript data-raw/expected/assoc_network.r
#
# NetCoMi 1.2.0 was installed on R 4.3.3 with stand-in packages (an empty
# function per imported name) for the imports that its correlation-input
# path never calls: mia, mixedCCA, phyloseq, SPRING, SpiecEasi,
# SummarizedExperiment, Biobase and qgraph. igraph 1.6.0, WGCNA 1.72-5,
# fdrtool, orca, pulsar and the rest were the real packages. (NetCoMi 1.3.0
# needs R >= 4.6; its NEWS lists no change to these computations.)
#
# Every network is built from a fixed correlation matrix
# (dataType = "correlation", sparsMethod = "threshold"), so the fixtures
# test only the network layer.

suppressMessages({
  library(NetCoMi)
  library(igraph)
})

mk <- function (n, D) {
  b <- matrix(rnorm(n * D), n, D)
  for (i in seq_len(D %/% 2)) {
    a <- sample(D, 2)
    b[, a[2]] <- b[, a[2]] + b[, a[1]] * runif(1, -1.5, 1.5)
  }
  m <- cor(b)
  dimnames(m) <- list(paste0('T', 1:D), paste0('T', 1:D))
  m
}

analyze <- function (A1, A2, thresh, diss, wdeg, hubs, hq) {
  net <- suppressMessages(netConstruct(
    data = A1, data2 = A2, dataType = 'correlation', sparsMethod = 'threshold',
    thresh = thresh, dissFunc = diss, verbose = 0 ))
  pr <- suppressMessages(netAnalyze(
    net, graphlet = FALSE, gcmHeat = FALSE, verbose = 0, weightDeg = wdeg,
    hubPar = hubs, hubQuant = hq, connectivity = FALSE ))
  list(net = net, props = pr)
}

ref_of <- function (pr, g) {
  ce <- pr$centralities
  gp <- pr$globalProps
  k  <- function (nm) gp[[paste0(nm, g)]]
  list(
    degree      = unname(ce[[paste0('degree',  g)]]),
    betweenness = unname(ce[[paste0('between', g)]]),
    closeness   = unname(ce[[paste0('close',   g)]]),
    eigenvector = unname(ce[[paste0('eigenv',  g)]]),
    cluster     = unname(pr$clustering[[paste0('clust', g)]]),
    hubs        = pr$hubs[[g]],
    global      = c(
      n_components = k('nComp'), average_dissimilarity = k('avDiss'),
      average_path_length = k('avPath'), clustering_coefficient = k('clustCoef'),
      modularity = k('modularity'), density = k('density'),
      positive_fraction = k('pep') / 100,
      lcc_size = pr$globalPropsLCC[[paste0('lccSize', g)]] ) )
}


#--------------------------------------------------------------------------
# 1. assoc_network(): 12 pairs of matrices, varied settings.
#--------------------------------------------------------------------------
set.seed(20261120)
sign_of <- c(signed = 'signed', unsigned = 'unsigned', signedPos = 'positive')
settings <- list(
  list(D =  8, thresh = 0.3, diss = 'signed',    wdeg = FALSE, hubs = 'eigenvector',                 hq = 0.95),
  list(D = 15, thresh = 0.3, diss = 'signed',    wdeg = FALSE, hubs = 'eigenvector',                 hq = 0.75),
  list(D = 15, thresh = 0.2, diss = 'unsigned',  wdeg = TRUE,  hubs = 'betweenness',                 hq = 0.90),
  list(D = 25, thresh = 0.2, diss = 'signedPos', wdeg = FALSE, hubs = c('degree', 'betweenness'),    hq = 0.90),
  list(D = 25, thresh = 0.4, diss = 'unsigned',  wdeg = FALSE, hubs = 'degree',                      hq = 0.90),
  list(D = 25, thresh = 0.3, diss = 'signed',    wdeg = TRUE,  hubs = c('closeness', 'eigenvector'), hq = 0.75),
  list(D = 30, thresh = 0.3, diss = 'signedPos', wdeg = TRUE,  hubs = 'betweenness',                 hq = 0.90),
  list(D = 30, thresh = 0.4, diss = 'signed',    wdeg = FALSE, hubs = 'closeness',                   hq = 0.95),
  list(D = 12, thresh = 0.25, diss = 'unsigned', wdeg = FALSE, hubs = c('degree', 'eigenvector'),    hq = 0.75),
  list(D = 20, thresh = 0.35, diss = 'signedPos', wdeg = FALSE, hubs = 'eigenvector',                hq = 0.95),
  list(D = 30, thresh = 0.2, diss = 'signed',    wdeg = FALSE, hubs = 'degree',                      hq = 0.95),
  list(D = 10, thresh = 0.3, diss = 'signed',    wdeg = FALSE, hubs = 'eigenvector',                 hq = 0.90) )

networks <- lapply(settings, function (s) {
  A1 <- mk(sample(25:50, 1), s$D)
  A2 <- mk(sample(25:50, 1), s$D)
  r  <- analyze(A1, A2, s$thresh, s$diss, s$wdeg, s$hubs, s$hq)
  list(
    A1 = A1, A2 = A2,
    args = list(threshold = s$thresh, sign = sign_of[[s$diss]],
                weighted_degree = s$wdeg, hubs = s$hubs, hub_quantile = s$hq),
    ref1 = ref_of(r$props, 1), ref2 = ref_of(r$props, 2),
    jaccard = local({
      cmp <- netCompare(r$props, permTest = FALSE, testRand = FALSE, gcd = FALSE, verbose = FALSE)
      rbind(
        degree = cmp$jaccDeg, betweenness = cmp$jaccBetw, closeness = cmp$jaccClose,
        eigenvector = cmp$jaccEigen, hubs = cmp$jaccHub )
    }),
    ari = suppressWarnings(netCompare(
      r$props, permTest = FALSE, testRand = FALSE, gcd = FALSE, verbose = FALSE)$randInd[['value']]) )
})


#--------------------------------------------------------------------------
# 2. Edge cases: a pair with r = 1, a complete network (thresh = 0, where
#    NetCoMi switches to weighted degree), a single edge.
#--------------------------------------------------------------------------
A <- networks[[2]]$A1
A['T2', ] <- A['T1', ]; A[, 'T2'] <- A[, 'T1']; A['T2', 'T2'] <- 1
r <- analyze(A, A, 0.3, 'signed', FALSE, 'eigenvector', 0.95)
edge_r1 <- list(A = A, ref = ref_of(r$props, 1))

r <- analyze(networks[[1]]$A1, networks[[1]]$A2, 0, 'signed', FALSE, 'eigenvector', 0.95)
complete <- list(A = networks[[1]]$A1, ref = ref_of(r$props, 1))

B <- diag(4); B[1, 2] <- B[2, 1] <- 0.5; dimnames(B) <- list(letters[1:4], letters[1:4])
r <- analyze(B, B, 0.3, 'signed', FALSE, 'eigenvector', 0.95)
one_edge <- list(A = B, ref = ref_of(r$props, 1))


#--------------------------------------------------------------------------
# 3. NetCoMi's Jaccard probabilities for sets of sizes a, b sharing k.
#--------------------------------------------------------------------------
calc_jacc <- get('.calcJaccard', asNamespace('NetCoMi'))
jacc_grid <- NULL
for (a in 0:6) for (b in 0:6) for (k in 0:min(a, b)) {
  A <- letters[seq_len(a)]
  B <- c(letters[seq_len(k)], LETTERS[seq_len(b - k)])
  jacc_grid <- rbind(jacc_grid, c(a = a, b = b, k = k, unlist(calc_jacc(A, B, sigTest = TRUE))))
}


#--------------------------------------------------------------------------
# 4. diffnet(): Fisher's z test, with sample sizes given.
#--------------------------------------------------------------------------
r  <- analyze(networks[[2]]$A1, networks[[2]]$A2, 0.3, 'signed', FALSE, 'eigenvector', 0.95)
dn <- diffnet(r$net, diffMethod = 'fisherTest', n1 = 40, n2 = 33, adjust = 'BH', verbose = FALSE)
diff <- list(A1 = networks[[2]]$A1, A2 = networks[[2]]$A2, n1 = 40, n2 = 33,
             p_value = dn$pvalsMat, p_adjusted = dn$pAdjustMat)


#--------------------------------------------------------------------------
# 5. igraph on hand-built graphs (a second check of the textbook
#    formulas): betweenness on lengths, harmonic closeness sums, eigenvector
#    centrality, Barrat's coefficient, and fast greedy clusters with their
#    weighted modularity.
#--------------------------------------------------------------------------
hand <- list(
  star  = rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5)),
  path  = rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5)),
  ring  = rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 6), c(6, 1)),
  cliques = rbind(t(combn(1:4, 2)), t(combn(5:8, 2)), c(4, 5)),
  random  = local({ set.seed(3); e <- t(combn(1:12, 2)); e[runif(nrow(e)) < 0.3, ] }) )
igraph_ref <- lapply(hand, function (el) {
  n <- max(el)
  set.seed(nrow(el))
  w <- round(runif(nrow(el), 0.1, 0.9), 3)     # similarities s
  d <- 1 - w                                   # lengths
  g <- graph_from_edgelist(el, directed = FALSE)
  sp <- distances(g, weights = d)
  fg <- cluster_fast_greedy(g, weights = w)
  list(
    edges = el, s = w,
    betweenness = betweenness(g, weights = d, normalized = FALSE),
    harmonic    = rowSums(ifelse(is.finite(sp) & sp > 0, 1 / sp, 0)),
    eigenvector = eigen_centrality(g, weights = w)$vector,
    barrat      = transitivity(g, type = 'barrat', weights = w),
    fg_membership = as.vector(membership(fg)),
    fg_modularity = max(fg$modularity) )
})


saveRDS(
  list(
    networks = networks, edge_r1 = edge_r1, complete = complete, one_edge = one_edge,
    jacc_grid = jacc_grid, diff = diff, igraph = igraph_ref,
    versions = c(NetCoMi = as.character(packageVersion('NetCoMi')),
                 igraph  = as.character(packageVersion('igraph'))) ),
  file.path('tests', 'tinytests', 'expected', 'assoc_network.rds'),
  version = 2 )
