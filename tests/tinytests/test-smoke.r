#test_that("smoke", {

  # Smoke tests: the only tests run on CRAN (see tests/tinytest.r). Every
  # exported function runs on small inputs and returns results of the right
  # shape, class and attributes, with finite values. Nothing here compares
  # numbers with stored results or with references computed by R (whose
  # BLAS/LAPACK and math library vary between machines); those tests run at
  # home and on GitHub CI, with NOT_CRAN=true.

  cpus <- 2L

  set.seed(1)
  counts <- matrix(rpois(20 * 12, 5) + 1L, 20, 12,
    dimnames = list(paste0('S', 1:20), paste0('F', 1:12)) )
  meta <- data.frame(
    row.names = rownames(counts),
    grp       = rep(c('a', 'b'), each = 10),
    x         = seq(0, 1, length.out = 20) )
  dm <- bray(counts, cpus = cpus)



  # Diversity metrics. ====

  for (m in list_metrics('alpha', 'id')) {
    a <- if (m == 'faith') alpha_div(ecodive::ex_counts, m, tree = ecodive::ex_tree) else alpha_div(counts, m, cpus = cpus)
    expect_true(is.numeric(a) && length(a) > 0, info = m)  # chao1 is NaN without singletons and doubletons
  }

  for (m in list_metrics('beta', 'id')) {
    b <- if (m %in% list_metrics('beta', 'id', phylo = TRUE))
      beta_div(ecodive::ex_counts, m, tree = ecodive::ex_tree) else beta_div(counts, m, cpus = cpus)
    expect_inherits(b, 'dist', info = m)
    expect_true(all(is.finite(b)), info = m)  # morisita can be < 0
  }

  expect_identical(match_metric('bray')$id, 'bray')

  r <- rarefy(counts, depth = 40, cpus = cpus)
  expect_equal(unname(rowSums(r)), rep(40, 20))

  for (nm in c('percent', 'chord', 'binary', 'clr', 'rclr')) {
    z <- normalize(counts, nm, cpus = cpus)
    expect_identical(dim(z), dim(counts), info = nm)
    expect_true(all(is.finite(as.matrix(z))), info = nm)
  }

  tr <- read_tree('((A:1,B:2):1,C:3);')
  expect_inherits(tr, 'phylo')
  expect_identical(sort(tr$tip.label), c('A', 'B', 'C'))



  # Ordinations. ====

  is_ord <- function (o, prefix, n = 20L) {
    is.matrix(o) && inherits(o, 'ecodive_ord') && nrow(o) == n &&
      ncol(o) >= 1 && startsWith(colnames(o)[[1]], prefix) &&
      identical(names(attributes(o))[3:5], c('method', 'eig', 'prop')) &&
      all(is.finite(unclass(o)))
  }

  expect_true(is_ord(pcoa(dm, cpus = cpus), 'PCoA'))
  expect_true(is_ord(pcoa(dm, correction = 'cailliez', cpus = cpus), 'PCoA'))
  expect_true(is_ord(mds(dm, n_starts = 1, max_iters = 50, cpus = cpus), 'MDS'))
  expect_true(is_ord(nmds(dm, n_starts = 1, max_iters = 50, cpus = cpus), 'NMDS'))
  expect_true(is_ord(tsne(dm, perplexity = 5, max_iters = 50, cpus = cpus), 'TSNE'))
  expect_true(is_ord(umap(dm, n_neighbors = 5, n_epochs = 20, cpus = cpus), 'UMAP'))
  expect_true(is_ord(isomap(dm, cpus = cpus), 'Isomap'))
  expect_true(is_ord(dbrda(dm, 'grp', meta, cpus = cpus), 'dbRDA'))
  expect_true(is_ord(capscale(dm, 'grp', meta, cpus = cpus), 'CAP'))
  expect_true(is_ord(pca(counts, cpus = cpus), 'PC'))
  expect_true(is_ord(rpca(counts, cpus = cpus), 'RPCA'))
  expect_true(is_ord(rda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus), 'RDA'))
  expect_true(is_ord(cca(counts, 'grp', meta, permutations = 9, cpus = cpus), 'CCA'))
  expect_true(is_ord(ca(counts, cpus = cpus), 'CA'))
  expect_true(is_ord(lda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus), 'LD'))
  expect_true(is_ord(plsda(counts, 'grp', meta, norm = 'clr', pseudocount = 1, cpus = cpus), 'PLSDA'))
  expect_true(is_ord(dist_ord(dm, 'pcoa', cpus = cpus), 'PCoA'))
  expect_true(all(c('pcoa', 'nmds', 'umap') %in% list_dist_ords('id')))
  expect_identical(match_dist_ord('cmdscale')$id, 'pcoa')
  expect_true(is_ord(count_ord(counts, 'pca', cpus = cpus), 'PC'))
  expect_identical(list_count_ords('id'), c('ca', 'cca', 'lda', 'pca', 'plsda', 'rda', 'rpca'))
  expect_identical(match_count_ord('prcomp')$id, 'pca')

  o  <- pcoa(dm, cpus = cpus)
  bp <- pcoa_biplot(o, counts, cpus = cpus)
  expect_identical(dim(bp), c(ncol(counts), ncol(o)))
  expect_true(all(is.finite(bp)))

  ef <- envfit(o, meta, permutations = 9, cpus = cpus)
  expect_inherits(ef, 'data.frame')
  expect_true(all(ef$p_value >= 0 & ef$p_value <= 1))

  pr <- procrustes(o, pcoa(dm, n_axes = 2, cpus = cpus), cpus = cpus)
  expect_true(is_ord(pr, 'PCoA'))
  expect_true(is.finite(attr(pr, 'ss')))

  expect_true(any(grepl('PCoA1', capture.output(print(o)))))

  sp <- simper(counts, 'grp', meta, permutations = 9, cpus = cpus)
  expect_inherits(sp, c('ecodive_simper', 'data.frame'))
  expect_identical(dim(sp), c(12L, 10L))
  expect_true(all(is.finite(sp$average)) && all(sp$p_value > 0 & sp$p_value <= 1))

  be <- bioenv(dm, cbind(meta, z = rep(1:4, 5)), cpus = cpus)
  expect_inherits(be, c('ecodive_bioenv', 'data.frame'))
  expect_identical(dim(be), c(2L, 3L))
  expect_true(all(is.finite(be$correlation)) && all(abs(be$correlation) <= 1))



  # Distance-based tests. ====

  tests <- list(
    permanova = permanova(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    mrpp      = mrpp(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    anosim    = anosim(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    permdisp  = permdisp(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    mantel    = mantel(dm, 'x', meta, permutations = 9, cpus = cpus),
    mantel_partial = mantel_partial(dm, 'x', rep(1:4, 5), meta, permutations = 9, cpus = cpus),
    procrustes_test = procrustes_test(o, pcoa(dm, n_axes = 2, cpus = cpus), permutations = 9, cpus = cpus),
    dist_test = dist_test(dm, 'permanova', 'grp', meta, permutations = 9, cpus = cpus) )

  for (nm in names(tests)) {
    t <- tests[[nm]]
    expect_inherits(t, 'ecodive_test', info = nm)
    expect_true(is.finite(t$statistic), info = nm)
    expect_true(t$p_value >= 0 && t$p_value <= 1, info = nm)
    expect_identical(nrow(as.data.frame(t)), 1L, info = nm)
    expect_true(length(capture.output(print(t))) > 0, info = nm)
    if (!is.null(t$pairwise)) {
      expect_inherits(t$pairwise, 'data.frame', info = nm)
      expect_identical(rownames(t$pairwise), 'a-b', info = nm)
      expect_true(all(t$pairwise$p_adjusted >= 0 & t$pairwise$p_adjusted <= 1), info = nm)
    }
  }
  for (nm in c('permanova', 'mrpp', 'anosim', 'permdisp'))
    expect_identical(tests[[nm]]$p_adjust, 'holm', info = nm)
  expect_true(all(c('permanova', 'mantel') %in% list_dist_tests()$id))
  expect_identical(match_dist_test('adonis2')$id, 'permanova')



  # Input checks. ====

  bad <- dm
  bad[1] <- -1
  expect_error(pcoa(bad), 'negative')
  bad[1] <- NA
  expect_error(permanova(bad, 'grp', meta), 'NA')
  expect_error(pcoa(dm, n_axes = 0), '`n_axes`')
  expect_error(dist_ord(dm, 'no_such_method'))
  expect_error(dist_test(dm, 'no_such_test', 'grp', meta))

  expect_true(n_cpus() >= 1)
  bt <- blas_threads()
  expect_true(is.integer(bt) && length(bt) == 1)

#})
