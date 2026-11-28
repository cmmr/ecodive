# mpd(), mntd(), weighted_mpd(), weighted_mntd(), mpd_ses(), mntd_ses().
# Fixtures: data-raw/expected/phylo_structure.r (picante 1.8.4).

fx <- readRDS(file.path('expected', 'phylo_structure.rds'))

max_cpus <- if (at_home()) 4L else 2L

rel  <- function (x, y) {
  x <- unname(x); y <- unname(y)
  if (!identical(is.na(x), is.na(y))) return (Inf)
  ok <- !is.na(x)
  if (!any(ok)) return (0)
  max(abs(x[ok] - y[ok]) / pmax(abs(y[ok]), 1))
}

# Patristic distances between the tips of a phylo list, in plain R, by
# summing branch lengths along each path (no depth differences).
patristic <- function (tr) {
  ntip <- length(tr$tip.label)
  par  <- integer(ntip + tr$Nnode); len <- numeric(ntip + tr$Nnode)
  par[tr$edge[, 2]] <- tr$edge[, 1]; len[tr$edge[, 2]] <- tr$edge.length
  up <- lapply(seq_len(ntip), function (t) {
    nodes <- t; lens <- numeric(0)
    while (par[nodes[length(nodes)]] > 0) {
      v <- nodes[length(nodes)]
      nodes <- c(nodes, par[v]); lens <- c(lens, len[v])
    }
    list(nodes = nodes, lens = lens)
  })
  d <- matrix(0, ntip, ntip, dimnames = list(tr$tip.label, tr$tip.label))
  for (i in seq_len(ntip)) for (j in seq_len(ntip)) if (i < j) {
    a <- match(up[[i]]$nodes, up[[j]]$nodes)
    f <- which(!is.na(a))[1]
    d[i, j] <- d[j, i] <- sum(up[[i]]$lens[seq_len(f - 1)]) + sum(up[[j]]$lens[seq_len(a[f] - 1)])
  }
  d
}

# MPD and MNTD in plain R, from a distance matrix.
r_four <- function (x, d) {
  x <- as.matrix(x)
  one <- function (v, stat, w) {
    sp <- colnames(x)[v > 0]
    if (length(sp) < 2) return (NA_real_)
    dd <- d[sp, sp]; a <- v[v > 0]
    if (stat == 'mpd') {
      if (w) sum(outer(a, a) * dd) / sum(a)^2 else mean(dd[lower.tri(dd)])
    } else {
      diag(dd) <- Inf; m <- apply(dd, 1, min)
      if (w) sum(a * m) / sum(a) else mean(m)
    }
  }
  list(
    mpd   = apply(x, 1, one, 'mpd', FALSE),  wmpd  = apply(x, 1, one, 'mpd', TRUE),
    mntd  = apply(x, 1, one, 'mntd', FALSE), wmntd = apply(x, 1, one, 'mntd', TRUE) )
}

e_four <- function (x, tree = NULL, taxa_dm = NULL, ...) list(
  mpd   = mpd(x, tree, taxa_dm, ...),          wmpd  = weighted_mpd(x, tree, taxa_dm, ...),
  mntd  = mntd(x, tree, taxa_dm, ...),         wmntd = weighted_mntd(x, tree, taxa_dm, ...) )


#--------------------------------------------------------------------
# Observed values against picante.
#--------------------------------------------------------------------
samp <- fx$samp; tree <- fx$tree; dis <- fx$dis

got <- e_four(samp, tree)
for (k in names(got)) {
  expect_true(rel(got[[k]], fx$obs[[k]]) < 1e-12, info = k)
  expect_identical(names(got[[k]]), rownames(samp), info = k)
}
got <- e_four(samp, taxa_dm = dis)
for (k in names(got)) expect_true(rel(got[[k]], fx$obs[[k]]) < 1e-12, info = paste('dm', k))
got <- e_four(samp, taxa_dm = stats::as.dist(dis))
for (k in names(got)) expect_true(rel(got[[k]], fx$obs[[k]]) < 1e-12, info = paste('dist', k))

big <- fx$big
got <- e_four(big$counts, big$tree, cpus = 1L)
for (k in names(got)) expect_true(rel(got[[k]], big$obs[[k]]) < 1e-12, info = paste('big', k))
expect_true(all(is.na(sapply(got, `[`, 1:2))))   # one taxon; empty

got <- e_four(samp, taxa_dm = fx$traits$taxa_dm)
for (k in names(got)) expect_true(rel(got[[k]], fx$traits$obs[[k]]) < 1e-12, info = paste('traits', k))

# Against plain R on path sums.
pd  <- patristic(tree)
ref <- r_four(samp, pd)
got <- e_four(samp, tree)
for (k in names(got)) expect_true(rel(got[[k]], ref[[k]]) < 1e-12, info = paste('R', k))

# Deep tree with short tips: depths near 1e6, tip distances near 1e-3.
deep <- read_tree(paste0(
  '(((a:0.001,b:0.002):0.0005,(c:0.003,d:0.001):0.002):1000000,',
  '((e:0.002,f:0.004):1e-4,(g:0.001,(h:0.0007,i:0.0009):0.0002):3e-4):999999.5);' ))
x <- matrix(c(1, 2, 0, 3, 1, 1, 0, 2, 5,
              0, 1, 4, 1, 0, 2, 1, 1, 1), 2, 9, byrow = TRUE,
            dimnames = list(c('u', 'v'), letters[1:9]))
ref <- r_four(x, patristic(deep))
got <- e_four(x, deep)
for (k in names(got)) expect_true(rel(got[[k]], ref[[k]]) < 1e-10, info = paste('deep', k))
expect_true(max(abs(got$mntd - ref$mntd) / ref$mntd) < 1e-6)  # relative to the tiny values


#--------------------------------------------------------------------
# alpha_div() and METRICS.
#--------------------------------------------------------------------
expect_true(all(c('mpd', 'mntd', 'weighted_mpd', 'weighted_mntd') %in% list_metrics('alpha', 'id')))
expect_identical(alpha_div(samp, 'mpd', tree = tree), mpd(samp, tree))
expect_identical(alpha_div(samp, 'weighted_mntd', tree = tree), weighted_mntd(samp, tree))
expect_identical(alpha_div(samp, 'wmpd', taxa_dm = dis), weighted_mpd(samp, taxa_dm = dis))
expect_identical(match_metric('mnnd')$id, 'mntd')
expect_identical(match_metric('Mean Pairwise Distance')$id, 'mpd')
expect_true(all(list_metrics('alpha', 'phylo', nm = 'id')[c('mpd', 'mntd', 'weighted_mpd', 'weighted_mntd')]))
expect_identical(unname(list_metrics('alpha', 'weighted', nm = 'id')[c('mpd', 'mntd', 'weighted_mpd', 'weighted_mntd')]), c(FALSE, FALSE, TRUE, TRUE))


#--------------------------------------------------------------------
# Invariances and input types.
#--------------------------------------------------------------------
# Scaling a sample; proportions.
expect_true(rel(weighted_mpd(samp * 7.5, tree), weighted_mpd(samp, tree)) < 1e-14)
expect_true(rel(weighted_mntd(samp / rowSums(samp), tree), weighted_mntd(samp, tree)) < 1e-14)
expect_identical(mpd(samp > 0, tree), mpd(samp, tree))

# Polytomies, single-child nodes and extra tips don't change values.
t1 <- read_tree('((A:1,B:2,C:0.5):1,(D:1,(E:2):0.5):0.25,F:3);')
t2 <- read_tree('(((A:1,B:2):0,C:0.5):1,((D:1,E:2.5):0.25,F:3):0);')
t3 <- read_tree('((A:1,B:2,C:0.5,X:4):1,(D:1,(E:2):0.5,Y:0.1):0.25,F:3);')
y  <- matrix(c(1, 2, 0, 1, 5, 1,  0, 1, 1, 1, 0, 3,  2, 0, 0, 0, 1, 1), 3, 6, byrow = TRUE,
             dimnames = list(c('p', 'q', 'r'), LETTERS[1:6]))
a <- e_four(y, t1); b <- e_four(y, t2); cc <- e_four(y, t3)
for (k in names(a)) {
  expect_true(rel(a[[k]], b[[k]]) < 1e-14, info = k)
  expect_true(rel(a[[k]], cc[[k]]) < 1e-14, info = k)
}
ref <- r_four(y, patristic(t1))
for (k in names(a)) expect_true(rel(a[[k]], ref[[k]]) < 1e-14, info = paste('poly', k))

# Sparse, dense, transposed; any cpus.
for (k in names(fx$obs)) {
  f  <- switch(k, mpd = mpd, wmpd = weighted_mpd, mntd = mntd, wmntd = weighted_mntd)
  d0 <- f(big$counts, big$tree, cpus = 1L)
  expect_identical(f(t(big$counts), big$tree, margin = 2L, cpus = 1L), d0, info = k)
  for (cpus in 2:max_cpus) expect_identical(f(big$counts, big$tree, cpus = cpus), d0, info = k)
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(f(Matrix::Matrix(big$counts, sparse = TRUE), big$tree, cpus = 1L), d0, info = k)
  if (requireNamespace('slam', quietly = TRUE))
    expect_identical(f(slam::as.simple_triplet_matrix(big$counts), big$tree, cpus = 1L), d0, info = k)
}

# An embedded tree.
z <- samp; attr(z, 'tree') <- tree
expect_identical(mpd(z), mpd(samp, tree))
expect_identical(mpd(z, taxa_dm = dis), mpd(samp, taxa_dm = dis))


#--------------------------------------------------------------------
# SES given picante's random communities (internal hook).
#--------------------------------------------------------------------
for (model in names(fx$hook)) {
  h <- fx$hook[[model]]
  for (wn in c('unweighted', 'weighted')) for (st in 0:1) {
    pr  <- if (st == 0L) 'mpd' else 'mntd'
    p   <- h$ses[[wn]][[pr]]
    for (route in c('tree', 'dm')) {
      r <- ecodive:::phylo_ses(
        samp, if (route == 'tree') tree, if (route == 'dm') dis, st, wn == 'weighted',
        'taxa_labels', 20L, 0L, 1L, 0L, 2L, null_tables = h$tables )
      info <- paste(model, wn, pr, route)
      expect_true(rel(r$obs,       p[[paste0(pr, '.obs')]])       < 1e-12, info = info)
      expect_true(rel(r$null_mean, p[[paste0(pr, '.rand.mean')]]) < 1e-12, info = info)
      expect_true(rel(r$null_sd,   p[[paste0(pr, '.rand.sd')]])   < 1e-12, info = info)
      expect_true(rel(r$z,         p[[paste0(pr, '.obs.z')]])     < 1e-12, info = info)
      expect_true(rel(r$rank,      p[[paste0(pr, '.obs.rank')]])  < 1e-12, info = info)
      expect_true(rel(r$p_value,   p[[paste0(pr, '.obs.p')]])     < 1e-12, info = info)
      expect_identical(r$n_taxa, as.integer(p$ntaxa), info = info)
    }
  }
}


#--------------------------------------------------------------------
# Results: columns, attributes, title, determinism.
#--------------------------------------------------------------------
models <- c('taxa_labels', 'richness', 'frequency', 'sample_pool',
            'phylogeny_pool', 'independent_swap', 'trial_swap')

r <- mpd_ses(samp, tree, n_draws = 99L)
expect_true(inherits(r, 'ecodive_mpd_ses') && is.data.frame(r))
expect_identical(names(r), c('sample', 'n_taxa', 'obs', 'null_mean', 'null_sd', 'z', 'nri', 'rank', 'p_value', 'n_draws'))
expect_identical(r$sample, rownames(samp))
expect_identical(r$nri, -r$z)
expect_identical(r$p_value, r$rank / (r$n_draws + 1))
expect_identical(attr(r, 'method'), 'mpd_ses')
expect_identical(attr(r, 'null_model'), 'taxa_labels')
expect_identical(attr(r, 'n_draws'), 99L)
expect_identical(attr(r, 'weighted'), FALSE)
expect_null(attr(r, 'n_swaps'))
expect_identical(attr(r, 'display')$title, 'Net Relatedness Index (NRI)')
r2 <- mntd_ses(samp, tree, weighted = TRUE, null_model = 'independentswap', n_draws = 9L, n_swaps = 50L)
expect_identical(names(r2)[[7]], 'nti')
expect_identical(attr(r2, 'null_model'), 'independent_swap')
expect_identical(attr(r2, 'n_swaps'), 50L)
expect_identical(attr(r2, 'weighted'), TRUE)
expect_identical(attr(r2, 'display')$title, 'Nearest Taxon Index (NTI)')
expect_identical(attr(mpd_ses(samp, tree, null_model = 'taxa.labels', n_draws = 5L), 'null_model'), 'taxa_labels')
expect_identical(attr(mpd_ses(samp, tree, null_model = 'Phylogeny.Pool', n_draws = 5L), 'null_model'), 'phylogeny_pool')

# Identical for any cpus, with tree or taxa_dm, sparse or dense; the seed.
x <- big$counts[3:40, ]
for (model in models) for (st in c('mpd', 'mntd')) {
  f  <- if (st == 'mpd') mpd_ses else mntd_ses
  a  <- f(x, big$tree, weighted = TRUE, null_model = model, n_draws = 25L, n_swaps = 200L, seed = 3L, cpus = 1L)
  for (cpus in 2:max_cpus)
    expect_identical(f(x, big$tree, weighted = TRUE, null_model = model, n_draws = 25L, n_swaps = 200L, seed = 3L, cpus = cpus), a, info = paste(model, st, cpus))
  if (requireNamespace('Matrix', quietly = TRUE))
    expect_identical(f(Matrix::Matrix(x, sparse = TRUE), big$tree, weighted = TRUE, null_model = model, n_draws = 25L, n_swaps = 200L, seed = 3L, cpus = 2L), a, info = paste(model, st, 'sparse'))
  expect_identical(f(t(x), big$tree, weighted = TRUE, null_model = model, n_draws = 25L, n_swaps = 200L, seed = 3L, margin = 2L), a, info = paste(model, st, 'margin'))
  b <- f(x, big$tree, weighted = TRUE, null_model = model, n_draws = 25L, n_swaps = 200L, seed = 4L, cpus = 2L)
  expect_false(identical(a$null_mean, b$null_mean), info = paste(model, st, 'seed'))
}
# The taxa_dm route draws the same communities as the tree route.
for (model in models) {
  a <- mntd_ses(samp, tree,         null_model = model, n_draws = 30L, n_swaps = 100L, seed = 5L)
  b <- mntd_ses(samp, taxa_dm = dis, null_model = model, n_draws = 30L, n_swaps = 100L, seed = 5L)
  expect_true(rel(a$null_mean, b$null_mean) < 1e-12 && rel(a$z, b$z) < 1e-12, info = model)
  expect_identical(a$rank, b$rank, info = model)
}


#--------------------------------------------------------------------
# Exact moments.
#--------------------------------------------------------------------
# Against every ordered draw of S of the 7 tips of a small tree.
small <- read_tree('(((t1:1,t2:2.5):0.5,t3:1.25):1,((t4:0.5,t5:0.75):2,(t6:1,t7:3):0.25):0.5);')
sx <- matrix(c(3, 1, 0, 0, 0,  1, 2, 5, 0, 0,  1, 1, 1, 7, 0,  2, 1, 4, 1, 3), 4, 5, byrow = TRUE,
             dimnames = list(paste0('s', 1:4), paste0('t', 1:5)))
pd7 <- patristic(small)
perms <- function (v) if (length(v) == 1) list(v) else
  do.call(c, lapply(seq_along(v), function (i) lapply(perms(v[-i]), function (p) c(v[i], p))))
enum <- function (a, w) {
  S <- length(a); vals <- c()
  cb <- utils::combn(7, S)
  for (k in seq_len(ncol(cb))) for (p in perms(cb[, k])) {
    dd <- pd7[p, p]
    vals <- c(vals, if (w) sum(outer(a, a) * dd) / sum(a)^2 else mean(dd[lower.tri(dd)]))
  }
  c(mean(vals), sqrt(mean((vals - mean(vals))^2)))
}
for (w in c(FALSE, TRUE)) {
  e  <- t(sapply(1:4, function (i) enum(sx[i, sx[i, ] > 0], w)))
  r  <- mpd_ses(sx, small, weighted = w, null_model = 'exact')
  r2 <- mpd_ses(sx, taxa_dm = pd7, weighted = w, null_model = 'exact')
  expect_true(rel(r$null_mean, e[, 1]) < 1e-13 && rel(r$null_sd, e[, 2]) < 1e-13, info = w)
  expect_true(rel(r2$null_mean, e[, 1]) < 1e-13 && rel(r2$null_sd, e[, 2]) < 1e-13, info = w)
  expect_true(all(is.na(r$rank)) && all(is.na(r$p_value)) && all(is.na(r$n_draws)))
  expect_identical(r$z, (r$obs - r$null_mean) / r$null_sd)
  expect_identical(attr(r, 'n_draws'), NA_integer_)
}
expect_error(mntd_ses(samp, tree, null_model = 'exact'), 'No exact null moments are implemented for MNTD')


#--------------------------------------------------------------------
# Null distributions (statistical; at home only).
#--------------------------------------------------------------------
if (at_home()) {

  N  <- 9999L
  nm <- c(taxa.labels = 'taxa_labels', phylogeny.pool = 'phylogeny_pool',
          independentswap = 'independent_swap', trialswap = 'trial_swap')
  for (model in names(nm)) {
    s <- fx$stats[[model]]
    for (k in c('mpd', 'wmpd', 'mntd')) {
      f <- if (k == 'mntd') mntd_ses else mpd_ses
      r <- f(samp, tree, weighted = k == 'wmpd', null_model = nm[[model]], n_draws = N, seed = 11L)
      se <- sqrt((s[[k]][[2]]^2 + r$null_sd^2) / N)
      expect_true(max(abs(r$null_mean - s[[k]][[1]]) / se) < 4.5, info = paste(model, k))
      expect_true(max(abs(r$null_sd / s[[k]][[2]] - 1)) < 0.08, info = paste(model, k))
    }
  }

  # phylogeny_pool and taxa_labels against the exact moments.
  e <- mpd_ses(samp, tree, null_model = 'exact')
  for (model in c('phylogeny_pool', 'taxa_labels')) {
    r <- mpd_ses(samp, tree, null_model = model, n_draws = N, seed = 12L)
    expect_true(max(abs(r$null_mean - e$null_mean) / (e$null_sd / sqrt(N))) < 4.5, info = model)
    expect_true(max(abs(r$null_sd / e$null_sd - 1)) < 0.05, info = model)
  }
  ew <- mpd_ses(samp, tree, weighted = TRUE, null_model = 'exact')
  r  <- mpd_ses(samp, tree, weighted = TRUE, null_model = 'phylogeny_pool', n_draws = N, seed = 13L)
  expect_true(max(abs(r$null_mean - ew$null_mean) / (ew$null_sd / sqrt(N))) < 4.5)

  # richness and sample_pool (picante's shuffles for these aren't
  # uniform): the exact moments over the table's own taxa.
  dtab <- dis[colnames(samp), colnames(samp)]
  for (w in c(FALSE, TRUE)) {
    e <- mpd_ses(samp, taxa_dm = dtab, weighted = w, null_model = 'exact')
    for (model in c('richness', 'sample_pool')) {
      r <- mpd_ses(samp, tree, weighted = w, null_model = model, n_draws = N, seed = 14L)
      expect_true(max(abs(r$null_mean - e$null_mean) / (e$null_sd / sqrt(N))) < 4.5, info = paste(model, w))
      expect_true(max(abs(r$null_sd / e$null_sd - 1)) < 0.05, info = paste(model, w))
    }
  }

  # frequency: against a uniform shuffle of each taxon's column, in R.
  set.seed(1)
  M  <- 4000L
  v  <- replicate(M, {
    z <- apply(samp, 2, sample)
    dimnames(z) <- dimnames(samp)
    unlist(r_four(z, dtab)[c('mpd', 'mntd')])
  })
  rm <- rowMeans(v, na.rm = TRUE)
  rs <- apply(v, 1, stats::sd, na.rm = TRUE)
  a  <- mpd_ses(samp, tree, null_model = 'frequency', n_draws = N, seed = 15L)
  b  <- mntd_ses(samp, tree, null_model = 'frequency', n_draws = N, seed = 15L)
  se <- sqrt(rs^2 / M + c(a$null_sd, b$null_sd)^2 / N)
  expect_true(max(abs(c(a$null_mean, b$null_mean) - rm) / se) < 4.5)
}


#--------------------------------------------------------------------
# Edge cases.
#--------------------------------------------------------------------
e3 <- matrix(c(0, 0, 0, 0,  0, 2, 0, 0,  1, 1, 0, 0,  1, 1, 1, 1), 4, 4, byrow = TRUE,
             dimnames = list(c('S0', 'S1', 'S2', 'S4'), c('A', 'B', 'C', 'D')))
t4 <- read_tree('((A:1,B:2):1,(C:1,D:3):0.5);')
expect_identical(is.na(mpd(e3, t4)), c(S0 = TRUE, S1 = TRUE, S2 = FALSE, S4 = FALSE))
expect_identical(is.na(weighted_mntd(e3, t4)), c(S0 = TRUE, S1 = TRUE, S2 = FALSE, S4 = FALSE))
expect_equal(mpd(e3, t4)[['S2']], 3)
expect_equal(mntd(e3, t4)[["S4"]], mean(c(3, 3, 3.5, 4)))  # C is nearest A
for (model in c(models, 'exact')) {
  r <- suppressWarnings(mpd_ses(e3, t4, null_model = model, n_draws = 20L, n_swaps = 20L))
  expect_identical(r$n_taxa, c(0L, 1L, 2L, 4L), info = model)
  expect_true(all(is.na(r[1:2, c('obs', 'z', 'nri', 'rank', 'p_value')])), info = model)
}
r <- mpd_ses(e3, t4, n_draws = 20L)
expect_true(all(is.na(r$null_mean[1:2])) && all(r$n_draws[1:2] == 0L))
# Every draw of a sample with all the taxa is the observed one: sd 0.
r <- mpd_ses(e3, t4, null_model = 'richness', n_draws = 20L)
expect_identical(r$null_sd[[4]], 0)
expect_true(is.nan(r$z[[4]]))
expect_identical(r$rank[[4]], 11)   # ties averaged: 1 + 20 / 2
# n_draws = 1: no sd.
r <- mntd_ses(e3, t4, n_draws = 1L)
expect_true(all(is.na(r$null_sd)) && all(r$n_draws[3:4] == 1L))
# A zero-length tree.
z0 <- t4; z0$edge.length[] <- 0
expect_identical(unname(mpd(e3, z0)[3:4]), c(0, 0))
# Identical samples.
idn <- rbind(a = c(A = 1, B = 2, C = 0, D = 1), b = c(A = 1, B = 2, C = 0, D = 1))
r <- mpd_ses(idn, t4, null_model = 'taxa_labels', n_draws = 30L)
expect_identical(r$null_mean[[1]], r$null_mean[[2]])
# Nested tables: no swap is possible.
nest <- rbind(a = c(A = 1, B = 1, C = 1, D = 0), b = c(A = 1, B = 1, C = 0, D = 0))
expect_warning(r <- mpd_ses(nest, t4, null_model = 'independent_swap', n_draws = 5L), 'No swap is possible')
expect_identical(r$null_mean, r$obs)
expect_warning(mpd_ses(nest, t4, null_model = 'trial_swap', n_draws = 5L), 'nested')
# frequency: draws with fewer than two taxa are left out.
fq <- rbind(a = c(A = 1, B = 1, C = 0, D = 0), b = c(A = 0, B = 0, C = 1, D = 0), c = c(A = 0, B = 0, C = 0, D = 1))
r <- mpd_ses(fq, t4, null_model = 'frequency', n_draws = 200L)
expect_true(all(r$n_draws < 200L) && all(r$n_draws > 0L))
expect_identical(r$p_value, r$rank / (r$n_draws + 1))
# One taxon in all: every value NA.
r <- mpd_ses(e3[, 'A', drop = FALSE], read_tree('(A:1,B:2);'), null_model = 'richness', n_draws = 3L)
expect_true(all(is.na(r$obs)) && all(is.na(r$z)))
expect_identical(r$n_taxa, c(0L, 0L, 1L, 1L))

# Errors.
expect_error(mpd(samp), '`tree` or `taxa_dm` is required')
expect_error(mpd(samp, tree, taxa_dm = dis), 'not both')
bad <- samp; colnames(bad)[1] <- 'nope'
expect_error(mpd(bad, tree), 'not in `tree`')
expect_error(mpd(bad, taxa_dm = dis), 'not in `taxa_dm`')
nl <- tree; nl$edge.length <- NULL
expect_error(mpd(samp, nl), 'branch lengths')
ng <- tree; ng$edge.length[1] <- -1
expect_error(mpd(samp, ng), 'negative')
expect_error(mpd(samp, taxa_dm = -dis), "negative")
expect_error(mpd(unname(samp), tree), 'taxa names')
expect_error(mpd_ses(samp, tree, null_model = 'nope'), "`null_model` must be one of 'taxa_labels', 'richness', 'frequency', 'sample_pool', 'phylogeny_pool', 'independent_swap', 'trial_swap', 'exact'")
expect_error(mpd_ses(samp, tree, n_draws = 0), '`n_draws`')
expect_error(mpd_ses(samp, tree, n_swaps = -1), '`n_swaps`')
expect_error(mpd_ses(samp, tree, weighted = NA), '`weighted`')
expect_error(mpd_ses(samp, tree, seed = 1.5), '`seed`')
neg <- samp; neg[1, 1] <- -1
expect_error(mpd(neg, tree), 'negative')
# Unlabelled taxa_dm: rows in the order of the taxa.
expect_identical(mpd(samp, taxa_dm = unname(dis[colnames(samp), colnames(samp)])), mpd(samp, taxa_dm = dis))
expect_error(mpd(samp, taxa_dm = unname(dis)), 'one row per taxon')


#--------------------------------------------------------------------
# print() and plot().
#--------------------------------------------------------------------
r  <- mpd_ses(samp, tree, n_draws = 99L)
rt <- mntd_ses(samp, tree, n_draws = 99L)
out <- utils::capture.output(print(r))
expect_true(any(grepl('Net Relatedness Index (NRI)', out, fixed = TRUE)))
out <- utils::capture.output(print(rt[1:3, ]))
expect_true(any(grepl('Nearest Taxon Index (NTI)', out, fixed = TRUE)))

drawn <- function (x, ...) {
  v <- withVisible(plot(x, ...))
  identical(v$value, x) && !v$visible
}
grDevices::pdf(NULL)
expect_true(drawn(r))
expect_true(drawn(rt, groups = c('a', 'a', 'b', 'b', 'c', 'c'), main = 'x', cex = 0.8))
expect_true(drawn(mpd_ses(samp, tree, null_model = 'exact')))
expect_error(plot(r, groups = 1:2), 'one value per row')
expect_error(plot(r[0, ]), 'no `nri` values')
grDevices::dev.off()

titles <- function (x, ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control('enable')
  plot(x, ...)
  txt <- unlist(lapply(grDevices::recordPlot()[[1]], function (op) {
    a <- unlist(op[[2]], recursive = TRUE)
    a[vapply(a, is.character, logical(1))]
  }))
  gsub('\n', ' ', unique(unlist(txt)))
}
expect_true('Net Relatedness Index (NRI)' %in% titles(r))
expect_true('Nearest Taxon Index (NTI)' %in% titles(rt))
expect_true('Null model: taxa_labels, 99 draws' %in% titles(r))
y <- r; attr(y, 'display')$title <- 'Mine'
expect_true('Mine' %in% titles(y))
