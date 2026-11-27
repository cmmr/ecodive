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
    expect_true(is.numeric(a) && length(a) > 0 && all(is.finite(a)), info = m)
  }

  # The alpha metrics ported from scikit-bio.
  expect_true(all(is.finite(dominance(counts, cpus = cpus))))
  expect_true(all(is.finite(gini(counts, cpus = cpus))))
  expect_true(all(is.finite(goods_coverage(counts, cpus = cpus))))
  expect_true(all(is.finite(heip(counts, cpus = cpus))))
  expect_true(all(is.finite(hill(counts, q = 1, cpus = cpus))))
  expect_true(all(is.finite(mcintosh_evenness(counts, cpus = cpus))))
  expect_true(all(is.finite(pielou(counts, cpus = cpus))))
  expect_true(all(is.finite(renyi(counts, q = Inf, cpus = cpus))))
  expect_true(all(is.finite(simpson_evenness(counts, cpus = cpus))))
  expect_true(all(is.finite(strong(counts, cpus = cpus))))
  expect_true(all(is.finite(tsallis(counts, q = 0.5, cpus = cpus))))

  for (m in list_metrics('beta', 'id')) {
    b <- if (m %in% list_metrics('beta', 'id', phylo = TRUE))
      beta_div(ecodive::ex_counts, m, tree = ecodive::ex_tree) else beta_div(counts, m, cpus = cpus)
    expect_inherits(b, 'dist', info = m)
    expect_true(all(is.finite(b)), info = m)  # morisita can be < 0
  }

  expect_identical(match_metric('bray')$id, 'bray')

  # Partition of beta diversity.
  bp <- beta_partition(counts, 'bray', cpus = cpus)
  expect_inherits(bp, 'ecodive_beta_partition')
  expect_true(all(vapply(bp, function (d) inherits(d, 'dist') && all(is.finite(d)), TRUE)))
  bm <- beta_partition_multi(counts, 'sorensen', n_sites = 5, n_draws = 20, cpus = cpus)
  expect_inherits(bm, 'ecodive_beta_partition_multi')
  expect_true(all(is.finite(bm$value)) && identical(dim(bm$draws), c(20L, 3L)))

  r <- rarefy(counts, depth = 40, cpus = cpus)
  expect_equal(unname(rowSums(r)), rep(40, 20))
  r <- rarefy(counts, depth = 31:50, cpus = cpus)
  expect_equal(unname(rowSums(r)), 31:50)

  # Diversity over repeated rarefactions.
  z <- beta_div_rarefied(counts, 'bray', depth = 40, times = 3, cpus = cpus)
  expect_true(inherits(z, 'ecodive_dist') && attr(z, 'Size') == 20L && all(is.finite(z)))
  z <- alpha_div_rarefied(counts, 'shannon', depth = 40, times = 3, cpus = cpus)
  expect_true(is.numeric(z) && length(z) == 20L && all(is.finite(z)))
  z <- alpha_div_rarefied(counts, 'hill', q = 1, depth = 40, times = 3, cpus = cpus)
  expect_true(is.numeric(z) && length(z) == 20L && all(is.finite(z)))

  # Rarefaction curves, and rarefying to equal coverage.
  z <- rarefaction_curve(counts, c('observed', 'shannon'), n_depths = 5, times = 2, cpus = cpus)
  expect_true(inherits(z, 'ecodive_rarefaction_curve') && all(is.finite(as.matrix(z[, 3:6]))))
  r <- rarefy(counts, coverage = 0.9, cpus = cpus)
  expect_true(is.matrix(r) && ncol(r) == 12L && all(rowSums(r) <= rowSums(counts)[rownames(r)]))

  # Absolute abundance.
  load <- seq(1, 2, length.out = 20) * 1e9
  z <- qmp_scale(counts, load, cpus = cpus)
  expect_identical(dim(z), dim(counts))
  expect_true(all(is.finite(z)) && all(z >= 0))
  z <- qmp_rarefy(counts, load, cpus = cpus)
  expect_identical(dim(z), dim(counts))
  expect_true(all(is.finite(z)) && all(z >= 0))
  expect_identical(length(attr(z, 'target')), 20L)
  expect_identical(length(attr(z, 'sampling_depth')), 20L)

  for (nm in c('percent', 'chord', 'binary', 'clr', 'rclr')) {
    z <- normalize(counts, nm, cpus = cpus)
    expect_identical(dim(z), dim(counts), info = nm)
    expect_true(all(is.finite(as.matrix(z))), info = nm)
  }

  # Log-ratio transforms.
  z <- alr(counts, denominator = 'F1', cpus = cpus)
  expect_identical(dim(z), c(20L, 11L))
  expect_true(all(is.finite(z)) && identical(attr(z, 'method'), 'alr'))
  sbp <- t(vapply(1:11, function (j) c(rep(0, j - 1), 1, rep(-1, 12 - j)), numeric(12)))
  for (b in list(NULL, sbp)) {
    z <- ilr(counts, basis = b, cpus = cpus)
    expect_identical(dim(z), c(20L, 11L))
    expect_true(all(is.finite(z)) && identical(dim(attr(z, 'basis')), c(11L, 12L)))
  }
  z <- philr(ecodive::ex_counts, ecodive::ex_tree, pseudocount = 1,
             part_weights = 'enorm.x.gm.counts', ilr_weights = 'blw.sqrt', cpus = cpus)
  expect_identical(dim(z), c(4L, 5L))
  expect_true(all(is.finite(z)) && identical(attr(z, 'method'), 'philr'))
  z <- philr_distance(ecodive::ex_counts, ecodive::ex_tree, pseudocount = 1, cpus = cpus)
  expect_true(inherits(z, 'dist') && attr(z, 'Size') == 4L && all(is.finite(z)))

  # Log ratios.
  z <- log_ratio(counts, c('F1', 'F2'), 3:5)
  expect_true(is.numeric(z) && length(z) == 20 && all(is.finite(z)))
  expect_identical(names(z), rownames(counts))
  z <- balance(counts, c('F1', 'F2'), 3:5)
  expect_true(is.numeric(z) && length(z) == 20 && all(is.finite(z)))
  expect_identical(attr(z, 'method'), 'balance')
  expect_true(inherits(z, 'ecodive_log_ratio') && !inherits(z[1:2], 'ecodive_log_ratio'))
  z <- log_ratio_features(pca(counts, cpus = cpus), n = 2)
  expect_true(is.character(z$numerator) && length(z$numerator) == 2 && length(z$denominator) == 2)

  # Proportionality.
  for (f in c('vlr', 'rho_p', 'phi_s')) {
    z <- get(f)(counts, cpus = cpus)
    expect_true(inherits(z, 'dist') && attr(z, 'Size') == 12L && length(z) == 66L, info = f)
    expect_true(all(is.finite(z)) && identical(attr(z, 'method'), f), info = f)
    expect_identical(attr(z, 'Labels'), colnames(counts), info = f)
  }
  z <- rho_p(counts, denominator = 'F1', cpus = cpus)
  expect_true(attr(z, 'Size') == 11L && all(is.finite(z)) && all(abs(z) <= 1))

  # Compositional diagnostics.
  z <- variation_array(counts, cpus = cpus)
  expect_true(inherits(z, 'ecodive_variation_array') && identical(dim(z), c(12L, 12L)) &&
              all(is.finite(z[upper.tri(z)])) && all(is.na(diag(z))) &&
              abs(sum(attr(z, 'clr_var')) - attr(z, 'total_var')) < 1e-10)
  z <- zero_pattern(counts - 1L)
  expect_true(inherits(z, 'ecodive_zero_pattern') && length(z) == 20L &&
              sum(attr(z, 'n_samples')) == 20L && ncol(attr(z, 'patterns')) == 12L)

  # Associations between two tables.
  z <- cross_assoc(counts[1:10, 1:5], data.frame(a = 1:10, b = (1:10)^2, c = c(5:1, 1:5)), cpus = cpus)
  expect_true(inherits(z, 'ecodive_cross_assoc') && nrow(z) == 15L && all(abs(z$association) <= 1) && all(z$p_value >= 0 & z$p_value <= 1))

  # Contaminant identification.
  dc_conc <- exp(seq(0, 3, length.out = 20))
  dc_neg  <- rep(c(FALSE, TRUE), c(16, 4))
  dc <- decontam(counts, dc_conc, dc_neg, cpus = cpus)
  expect_true(inherits(dc, 'ecodive_decontam') && nrow(dc) == 12L && is.logical(dc$contaminant))
  expect_true(all(dc$p >= 0 & dc$p <= 1) && identical(attr(dc, 'test'), 'combined'))
  z <- decontam_not(counts, dc_neg, cpus = cpus)
  expect_true(inherits(z, 'ecodive_decontam_not') && nrow(z) == 12L && is.logical(z$not_contaminant))

  # Sloan's neutral community model: a sparse table of equal depths (the
  # last feature tops each sample up), with no random numbers.
  nc <- outer(1:20, 1:30, function (i, j) (i * 7L + j * 13L) %% 11L + j %/% 3L)
  nc[outer(1:20, 1:30, function (i, j) (i * 7L + j * 13L) %% 11L) < (30:1) %/% 3L] <- 0L
  nc[, 30] <- nc[, 30] + (max(rowSums(nc)) - rowSums(nc)) + 5L
  dimnames(nc) <- list(rownames(counts), paste0('F', 1:30))
  nmod <- neutral_model(nc, 'grp', meta, cpus = cpus)
  expect_true(inherits(nmod, 'ecodive_neutral_model') && nrow(nmod) == 60L && is.factor(nmod$partition))
  expect_true(all(is.finite(attr(nmod, 'm'))) && all(attr(nmod, 'm') > 0) && identical(names(attr(nmod, 'm')), c('a', 'b')))

  # Longitudinal: 5 subjects at 4 states.
  lm_meta <- data.frame(
    row.names = rownames(counts), subj = rep(sprintf('m%d', 1:5), each = 4),
    day = rep(c(0, 3, 7, 14), 5), grp = rep(c('a', 'b'), each = 10) )
  shan <- shannon(counts)
  fdi <- first_differences(shan, 'subj', 'day', lm_meta)
  expect_true(inherits(fdi, 'ecodive_first_differences') && nrow(fdi) == 15L && all(is.finite(fdi$difference)))
  fdt <- first_distances(dm, 'subj', 'day', lm_meta, baseline = 0)
  expect_true(inherits(fdt, 'ecodive_first_distances') && nrow(fdt) == 15L && all(fdt$distance >= 0))
  vol <- volatility(shan, 'subj', 'day', 'grp', lm_meta, n_resamples = 50L)
  expect_true(inherits(vol, 'ecodive_volatility') && nrow(vol$means) == 8L && all(is.finite(vol$global)))
  fvo <- feature_volatility(counts, 'subj', 'day', 'grp', lm_meta, n_trees = 10L, n_resamples = 50L, cpus = cpus)
  expect_true(inherits(fvo, 'ecodive_feature_volatility') && nrow(fvo$predictions) == 4L &&
              all(is.finite(fvo$features$importance)) && nrow(fvo$means) == 8L * nrow(fvo$features))
  nm <- nmit(counts, 'subj', lm_meta, cpus = cpus)
  expect_true(inherits(nm, 'ecodive_dist') && attr(nm, 'Size') == 5L && all(is.finite(nm)) && all(nm >= 0))

  # Batch correction: MMUPHin and ComBat.
  bm <- batch_mmuphin(counts, 'grp', ~ x, meta, cpus = cpus)
  expect_true(inherits(bm, 'ecodive_batch') && identical(dim(bm), dim(counts)) && all(is.finite(bm)))
  expect_true(identical(dimnames(bm), dimnames(counts)) && all(bm == round(bm)) && nrow(attr(bm, 'batch_effects')) == 24L)
  bc <- batch_combat(counts, 'grp', ~ x, meta, norm = 'clr', cpus = cpus)
  expect_true(inherits(bc, 'ecodive_batch') && identical(dim(bc), dim(counts)) && all(is.finite(bc)))
  z <- batch_combat(counts, 'grp', data = meta, norm = 'log', parametric = FALSE, cpus = cpus)
  expect_true(inherits(z, 'ecodive_batch') && all(is.finite(z)) && is.logical(attr(z, 'adjusted')))

  # Associations: SparCC.
  z <- sparcc(counts[, 1:10], n_draws = 5L, permutations = 9L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_dist') && attr(z, 'Size') == 10L && length(z) == 45L)
  expect_true(all(is.finite(z)) && identical(attr(z, 'method'), 'sparcc'))
  expect_identical(attr(z, 'Labels'), colnames(counts)[1:10])
  expect_identical(names(attr(z, 'variances')), colnames(counts)[1:10])
  expect_true(length(attr(z, 'p_value')) == 45L && all(attr(z, 'p_value') > 0 & attr(z, 'p_value') <= 1))
  a <- assoc_pairs(z, threshold = 0, alpha = 1)
  expect_true(inherits(a, 'ecodive_assoc_pairs') && nrow(a) == 45L)
  expect_identical(names(a), c('from', 'to', 'association', 'p_value', 'p_adjusted'))
  expect_true(all(a$from %in% colnames(counts)) && all(is.finite(a$association)))

  # Associations by name, and networks.
  expect_identical(feature_assoc(counts, 'rho_p', cpus = cpus), rho_p(counts, cpus = cpus))
  expect_identical(list_feature_assocs('id'), c('vlr', 'rho_p', 'phi_s', 'sparcc'))
  expect_identical(match_feature_assoc('FastSpar')$id, 'sparcc')
  an <- assoc_network(z, threshold = 0.2, cpus = cpus)
  expect_true(inherits(an, 'ecodive_network') && nrow(an$nodes) == 10L)
  expect_identical(an$nodes$feature, colnames(counts)[1:10])
  expect_true(all(is.finite(an$global[c('n_nodes', 'n_edges', 'density', 'modularity')])))
  expect_true(all(is.finite(as.matrix(an$nodes[2:5]))) && is.logical(an$nodes$hub))
  z2 <- sparcc(counts[, 1:10], n_draws = 5L, seed = 1L, cpus = cpus)
  ac <- assoc_network_compare(an, assoc_network(z2, threshold = 0.2), permutations = 9L, cpus = cpus)
  expect_true(inherits(ac, 'ecodive_network_compare') && nrow(ac$jaccard) == 5L)
  expect_true(nrow(ac$nodes) == 40L && all(ac$jaccard$p_less >= 0 & ac$jaccard$p_less <= 1))
  ad <- assoc_network_diff(z, z2)
  expect_true(inherits(ad, 'ecodive_network_diff') && nrow(ad) == 45L)
  expect_true(is.logical(ad$reject) && all(ad$p_value[!is.na(ad$p_value)] <= 1))

  # Differential abundance.
  z <- ancom(counts, 'grp', meta, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancom') && nrow(z) == 12L)
  expect_identical(z$feature, colnames(counts))
  expect_true(is.integer(z$W) && all(z$W >= 0L & z$W <= 11L) && is.logical(z$reject))
  z <- ancom(counts, meta$grp, sig_test = 'kruskal', p_adjust = 'BH', cpus = cpus)
  expect_true(nrow(z) == 12L && !anyNA(z$reject))
  z <- dirmult_ttest(counts, 'grp', meta, n_draws = 8L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_dirmult_ttest') && nrow(z) == 12L)
  expect_identical(z$feature, colnames(counts))
  expect_true(all(is.finite(as.matrix(z[, 2:7]))) && is.logical(z$reject) && !anyNA(z$reject))
  expect_true(all(z$ci_lower <= z$log2_fc & z$log2_fc <= z$ci_upper))
  z <- dirmult_lm(counts, ~ grp + x, meta, n_draws = 8L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_dirmult_lm') && nrow(z) == 36L && all(is.finite(z$estimate)) && !anyNA(z$reject))
  z <- ancombc(counts, ~ grp + x, meta, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc') && nrow(z) == 36L)
  expect_identical(z$feature, rep(colnames(counts), each = 3L))
  expect_identical(unique(z$covariate), c('(Intercept)', 'grpb', 'x'))
  expect_true(all(is.finite(as.matrix(z[, 3:7]))) && is.logical(z$reject) && !anyNA(z$reject))

  # Source tracking.
  st <- sourcetracker(counts, ifelse(meta$x < 0.8, meta$grp, NA), depth = 20L, n_restarts = 1L, n_burnin = 5L, cpus = cpus)
  expect_true(inherits(st, 'ecodive_source') && identical(dim(st), c(4L, 3L)) && all(is.finite(st)))
  expect_true(all(abs(rowSums(st) - 1) < 1e-12) && identical(colnames(st), c('a', 'b', 'Unknown')))
  cz <- counts; cz[cz < 3] <- 0L
  z <- ancombc2(cz, 'grp', meta, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc') && nrow(z) == 24L && identical(attr(z, 'method'), 'ancombc2'))
  expect_true(all(z$p_value >= 0 & z$p_value <= 1) && !anyNA(z$reject))
  expect_true(length(attr(z, 'sample_fractions')) == 20L && all(is.finite(attr(z, 'sample_fractions'))))
  meta$g3 <- rep(c('a', 'b', 'c', 'd'), 5)
  fit <- ancombc2(counts, ~ g3 + x, meta, groups = 'g3', cpus = cpus)
  z <- ancombc2(counts, ~ g3 + x, meta, groups = 'g3', sensitivity = c(0.5, 1), cpus = cpus)
  expect_true(is.logical(z$passed_ss) && identical(z$reject_robust, z$reject & z$passed_ss))
  z <- ancombc_global(fit, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc_global') && nrow(z) == 12L && all(z$p_value >= 0 & z$p_value <= 1))
  z <- ancombc_pairwise(fit, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc_pairwise') && nrow(z) == 72L && !anyNA(z$reject))
  z <- ancombc_dunnett(fit, n_bootstraps = 20L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc_dunnett') && nrow(z) == 36L && !anyNA(z$reject))
  z <- ancombc_trend(fit, n_bootstraps = 20L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_ancombc_trend') && nrow(z) == 12L && all(is.finite(z$w)))
  z <- structural_zeros(cz, 'grp', meta)
  expect_true(is.logical(z) && identical(dim(z), c(12L, 2L)))
  z <- da_wilcoxon(cz, 'grp', meta, pseudocount = 1, cpus = cpus)
  expect_true(inherits(z, 'ecodive_da_wilcoxon') && nrow(z) == 12L && identical(attr(z, 'method'), 'da_wilcoxon'))
  expect_true(all(z$p_value >= 0 & z$p_value <= 1) && is.logical(z$reject) && !anyNA(z$reject))
  z <- da_ttest(cz, 'grp', meta, norm = 'percent', cpus = cpus)
  expect_true(inherits(z, 'ecodive_da_ttest') && nrow(z) == 12L && all(is.finite(as.matrix(z[, 2:10]))))
  expect_true(all(z$ci_lower <= z$diff & z$diff <= z$ci_upper))
  z <- da_lm(cz, ~ grp + x, meta, norm = 'rclr', cpus = cpus)
  expect_true(inherits(z, 'ecodive_da_lm') && nrow(z) == 36L && !anyNA(z$reject))
  expect_identical(unique(z$covariate), c('(Intercept)', 'grpb', 'x'))
  z <- lefse(cz, 'grp', meta, n_boots = 5L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_lefse') && nrow(z) == 12L && identical(attr(z, 'method'), 'lefse'))
  expect_true(all(z$p_value >= 0 & z$p_value <= 1) && is.logical(z$reject) && !anyNA(z$reject))
  expect_true(all(is.nan(z$lda_score) | z$lda_score >= 0))
  z <- da_test(counts, 'da_lm', 'grp', meta, cpus = cpus)
  expect_true(inherits(z, 'ecodive_da_lm') && nrow(z) == 24L)
  z <- suppressWarnings(linda(cz, ~ grp + x, meta, cpus = cpus))
  expect_true(inherits(z, 'ecodive_linda') && nrow(z) == 24L && !anyNA(z$reject))
  expect_true(all(is.finite(z$log2_fc)) && all(z$p_value >= 0 & z$p_value <= 1))
  z <- indicator_species(cz, 'grp', meta, permutations = 19, cpus = cpus)
  expect_true(inherits(z, 'ecodive_indicator_species') && nrow(z) == 12L && identical(attr(z, 'method'), 'indicator_species'))
  expect_true(all(is.na(z$p_value) | (z$p_value > 0 & z$p_value <= 1)) && is.logical(z$reject) && !anyNA(z$reject))
  z <- linda(counts, ~ grp + x + (1 | subj), data.frame(meta, subj = rep(sprintf('s%d', 1:5), 4)), cpus = cpus)
  expect_true(inherits(z, 'ecodive_linda') && nrow(z) == 24L && all(z$dof > 0) && !anyNA(z$reject))
  z <- maaslin3(cz, ~ grp + x, meta, n_sims = 200L, cpus = cpus)
  expect_true(inherits(z, 'ecodive_maaslin3') && nrow(z) == 24L && !anyNA(z$reject) && identical(attr(z, 'method'), 'maaslin3'))
  expect_true(all(z$p_value >= 0 & z$p_value <= 1, na.rm = TRUE) && any(is.finite(z$coef_abundance)))
  expect_true(is.data.frame(list_da_tests()) && 'da_ttest' %in% list_da_tests('id'))
  expect_identical(match_da_test('Mann-Whitney')$id, 'da_wilcoxon')
  z <- da_sensitivity(counts, 'da_ttest', 'grp', meta, pseudocounts = TRUE, cpus = cpus)
  expect_true(inherits(z, 'ecodive_da_sensitivity') && nrow(z) == 48L && identical(sort(unique(z$run)), 0:3))
  expect_true(is.data.frame(attr(z, 'summary')) && nrow(attr(z, 'summary')) == 12L && is.logical(attr(z, 'summary')$passed))

  # Zero replacement.
  cz <- counts; cz[cz < 4] <- 0L
  for (m in list_zero_replacements('id')) {
    z <- replace_zeros(cz, m)
    expect_identical(dim(z), dim(cz), info = m)
    expect_true(all(is.finite(z)) && all(z > 0) && identical(attr(z, 'method'), m), info = m)
  }
  z <- replace_zeros_mult(cz)
  expect_true(all(is.finite(z)) && all(z > 0) && identical(dimnames(z), dimnames(cz)))
  expect_true(is.finite(sum(alr(z, 'F1', cpus = cpus))))
  expect_true(all(is.finite(balance(z, 1:3, 4:6))))
  expect_true(all(is.finite(phi_s(z, cpus = cpus))))
  expect_true('mult' %in% list_zero_replacements('id'))
  expect_identical(match_zero_replacement('multRepl')$id, 'mult')
  z <- replace_zeros_czm(cz)
  expect_true(all(is.finite(z)) && all(z > 0) && identical(attr(z, 'method'), 'czm'))
  z <- replace_zeros_gbm(cz)
  expect_true(all(is.finite(z)) && all(z > 0) && max(abs(rowSums(z) - 1)) < 1e-12)
  z <- replace_zeros_sq(cz, adjust = FALSE)
  expect_true(all(is.finite(z)) && all(z > 0) && identical(attr(z, 'n_adjusted'), 0))
  z <- replace_zeros_bl(t(cz), margin = 2L)
  expect_true(all(is.finite(z)) && all(z > 0) && identical(dim(z), rev(dim(cz))))
  expect_identical(match_zero_replacement('cmultRepl')$id, 'gbm')

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
  expect_true(is_ord(capda(dm, 'grp', meta, permutations = 9, cpus = cpus), 'CAP'))
  expect_true(is_ord(pca(counts, cpus = cpus), 'PC'))
  expect_true(is_ord(rpca(counts, cpus = cpus), 'RPCA'))
  expect_true(is_ord(rpca_phylo(ecodive::ex_counts, ecodive::ex_tree, cpus = cpus), 'RPCA', n = 4L))
  expect_true(is_ord(dpcoa(ecodive::ex_counts, ecodive::ex_tree, cpus = cpus), 'DPCoA', n = 4L))
  expect_true(is_ord(rda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus), 'RDA'))
  expect_true(is_ord(cca(counts, 'grp', meta, permutations = 9, cpus = cpus), 'CCA'))
  expect_true(is_ord(ca(counts, cpus = cpus), 'CA'))
  expect_true(is_ord(lda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus), 'LD'))
  expect_true(is_ord(plsda(counts, 'grp', meta, norm = 'clr', pseudocount = 1, cpus = cpus), 'PLSDA'))
  expect_true(is_ord(oplsda(counts, 'grp', meta, norm = 'clr', pseudocount = 1, cpus = cpus), 'p1'))
  expect_true(is_ord(dist_ord(dm, 'pcoa', cpus = cpus), 'PCoA'))
  expect_true(all(c('pcoa', 'nmds', 'umap') %in% list_dist_ords('id')))
  expect_identical(match_dist_ord('cmdscale')$id, 'pcoa')
  expect_true(is_ord(count_ord(counts, 'pca', cpus = cpus), 'PC'))
  expect_identical(list_count_ords('id'), c('ca', 'cca', 'dpcoa', 'lda', 'oplsda', 'pca', 'plsda', 'rda', 'rpca', 'rpca_phylo'))
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

  mc <- mantel_correlogram(dm, cbind(meta$x, rep(1:4, 5)), permutations = 9, cpus = cpus)
  expect_inherits(mc, c('ecodive_mantel_correlogram', 'data.frame'))
  expect_identical(ncol(mc), 7L)
  expect_true(any(is.finite(mc$r)) && all(abs(mc$r) <= 1, na.rm = TRUE))
  expect_true(all(mc$p_adjusted >= 0 & mc$p_adjusted <= 1, na.rm = TRUE))

  mr <- mantel_regression(dm ~ x + grp, data = meta, z = rep(1:4, 5), permutations = 9, cpus = cpus)
  expect_inherits(mr, 'ecodive_mantel_regression')
  expect_identical(dim(mr$coefficients), c(4L, 5L))
  expect_true(all(is.finite(mr$coefficients$estimate)) && is.finite(mr$r_squared))
  expect_true(all(mr$coefficients$p_value > 0 & mr$coefficients$p_value <= 1))

  et <- enterotype(counts, n_types = 2:3, cpus = cpus)
  expect_inherits(et, c('ecodive_enterotype', 'ecodive_types', 'data.frame'))
  expect_identical(dim(et), c(20L, 3L))
  expect_true(is.factor(et$type) && all(is.finite(et$silhouette)))
  expect_true(all(is.finite(attr(et, 'selection')$ch)))

  vl <- valencia(round(ecodive::valencia_centroids * 1000), cpus = cpus)
  expect_inherits(vl, c('ecodive_valencia', 'ecodive_types', 'data.frame'))
  expect_identical(dim(vl), c(13L, 5L))
  expect_true(is.factor(vl$type) && all(vl$score > 0 & vl$score <= 1))

  dmx <- dmm(counts, n_types = 1:2, cpus = cpus)
  expect_inherits(dmx, c('ecodive_dmm', 'ecodive_types', 'data.frame'))
  expect_identical(dim(dmx), c(20L, 3L))
  expect_true(is.factor(dmx$type) && all(dmx$probability > 0 & dmx$probability <= 1))
  expect_true(all(is.finite(as.matrix(attr(dmx, 'selection')[, c('laplace', 'aic', 'bic', 'nle')]))))

  vp <- varpart(counts, ~ grp, ~ x, data = meta, norm = 'hellinger', cpus = cpus)
  expect_inherits(vp, c('ecodive_varpart', 'data.frame'))
  expect_identical(dim(vp), c(4L, 5L))
  expect_true(all(is.finite(vp$adj_r_squared)))
  expect_true(any(grepl('[a]', capture.output(print(vp)), fixed = TRUE)))

  fs <- forward_select(counts, ~ grp + x, meta, norm = 'hellinger', permutations = 9, alpha = 1, cpus = cpus)
  expect_inherits(fs, c('ecodive_forward_select', 'data.frame'))
  expect_identical(ncol(fs), 8L)
  expect_true(is.finite(attr(fs, 'global')$adj_r_squared))
  expect_true(any(grepl('Forward selection', capture.output(print(fs)), fixed = TRUE)))

  prc <- principal_response_curves(
    counts, rep(c('c', 'd1', 'd2', 'c'), 5), rep(1:5, each = 4), norm = 'hellinger',
    n_axes = 2, permutations = 9, cpus = cpus )
  expect_inherits(prc, 'ecodive_prc')
  expect_identical(dim(prc$coefficients), c(15L, 4L))
  expect_true(all(is.finite(prc$coefficients$PRC1)) && all(is.finite(prc$weights)))
  expect_true(all(prc$test$p_value > 0 & prc$test$p_value <= 1))
  expect_true(any(grepl('Principal Response Curves', capture.output(print(prc)), fixed = TRUE)))



  # Distance-based tests. ====

  tests <- list(
    permanova = permanova(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    mrpp      = mrpp(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    anosim    = anosim(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    permdisp  = permdisp(dm, 'grp', meta, pairwise = TRUE, permutations = 9, cpus = cpus),
    mantel    = mantel(dm, 'x', meta, permutations = 9, cpus = cpus),
    mantel_partial = mantel_partial(dm, 'x', rep(1:4, 5), meta, permutations = 9, cpus = cpus),
    moran_i   = moran_i('x', dm, meta, permutations = 9, cpus = cpus),
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



  # Plots, on a null device: each returns its argument invisibly. ====

  drawn <- function (x, ...) {
    v <- withVisible(plot(x, ...))
    identical(v$value, x) && !v$visible
  }
  grDevices::pdf(NULL)
  expect_true(drawn(dm))
  expect_true(drawn(dm, 'grp', meta))
  expect_true(drawn(o))
  expect_true(drawn(o, 'grp', meta))
  expect_true(drawn(o, 'x', meta, axes = c(1, 3)))
  expect_true(drawn(rda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus), 'grp', meta))
  expect_true(suppressWarnings(drawn(pca(counts, norm = 'clr', cpus = cpus), scaling = 'covariance', links = 1:3)))  # biplot scalings
  expect_true(drawn(rda(counts, norm = 'hellinger', cpus = cpus), scaling = 1))
  expect_true(drawn(prc))
  expect_true(drawn(lda(counts, 'grp', meta, norm = 'hellinger', cpus = cpus)))  # one axis
  expect_true(drawn(capda(dm, 'grp', meta, permutations = 0, cpus = cpus), 'grp', meta))  # one axis
  expect_true(drawn(pr))
  expect_true(drawn(ef))
  expect_true(drawn(ef, add = TRUE))
  expect_true(drawn(sp))
  expect_true(drawn(be))
  expect_true(drawn(an))
  expect_true(drawn(ac))
  expect_true(drawn(ad))
  expect_true(drawn(vp))
  expect_true(drawn(fs))
  expect_true(drawn(ancom(counts, 'grp', meta, cpus = cpus)))
  expect_true(drawn(dirmult_ttest(counts, 'grp', meta, n_draws = 8L, cpus = cpus)))
  expect_true(drawn(dirmult_lm(counts, ~ grp + x, meta, n_draws = 8L, cpus = cpus)))
  expect_true(drawn(sparcc(counts, n_draws = 5L, cpus = cpus)))
  expect_true(drawn(st))   # sourcetracker()
  expect_true(drawn(dc))   # decontam()
  expect_true(drawn(dc, 'prevalence'))
  expect_true(drawn(nmod))   # neutral_model()
  expect_true(drawn(fdi))    # first_differences()
  expect_true(drawn(fdt))    # first_distances()
  expect_true(drawn(vol))    # volatility()
  expect_true(drawn(fvo))    # feature_volatility()
  expect_true(drawn(bm))   # batch_mmuphin()
  expect_true(drawn(bc, type = 'means'))   # batch_combat()
  expect_true(is.data.frame(decontam_frequency_plot(counts, 'F1', dc_conc, dc_neg)))
  expect_true(drawn(fit))   # ancombc2()
  expect_true(drawn(ancombc_global(fit, cpus = cpus)))
  expect_true(drawn(ancombc_pairwise(fit, cpus = cpus)))
  expect_true(drawn(ancombc_dunnett(fit, n_bootstraps = 20L, cpus = cpus)))
  expect_true(drawn(ancombc_trend(fit, n_bootstraps = 20L, cpus = cpus)))
  expect_true(drawn(da_wilcoxon(counts, 'grp', meta, cpus = cpus)))
  expect_true(drawn(da_ttest(counts, 'grp', meta, cpus = cpus)))
  expect_true(drawn(da_lm(counts, ~ grp + x, meta, cpus = cpus)))
  expect_true(drawn(linda(counts, ~ grp + x, meta, cpus = cpus)))
  expect_true(drawn(maaslin3(cz, ~ grp + x, meta, n_sims = 200L, cpus = cpus)))
  for (nm in names(tests)) {
    expect_true(drawn(tests[[nm]]), info = nm)
    expect_true(drawn(tests[[nm]], type = 'null'), info = nm)
  }
  z <- withVisible(composition_barplot(counts, 'grp', meta, n_features = 5, heatmap = TRUE, cpus = cpus))
  expect_true(!z$visible && identical(dim(z$value$proportions), c(6L, 20L)) &&
              all(abs(colSums(z$value$proportions) - 1) < 1e-12) && all(is.finite(z$value$clr)))
  z <- withVisible(ternary_plot(counts, parts = list(1, 2, NA), groups = 'grp', data = meta,
                                center = TRUE, ellipses = 'confidence'))
  expect_true(!z$visible && nrow(z$value) == 20L && all(abs(rowSums(z$value[, 1:3]) - 1) < 1e-12) &&
              length(attr(z$value, 'ellipses')) == 2L)
  expect_true(drawn(variation_array(counts, cpus = cpus)))
  expect_true(drawn(zero_pattern(counts - 1L), n_patterns = 5, n_features = 8))
  grDevices::dev.off()



  # Prediction: random forests. ====

  rf <- random_forest(counts, 'grp', meta, n_trees = 25, n_folds = 2, cpus = cpus)
  expect_true(inherits(rf, 'ecodive_random_forest') && nrow(rf$predictions) == 20L &&
              nrow(rf$importance) == 12L && is.factor(predict(rf, counts, cpus = cpus)))
  rf <- random_forest(counts, 'x', meta, n_trees = 25, n_folds = 2, cpus = cpus)
  expect_true(identical(attr(rf, 'type'), 'regression') &&
              all(is.finite(predict(rf, counts, cpus = cpus))))



  # Spatial: Moran's eigenvector maps. ====

  mem <- dbmem(data.frame(row.names = rownames(counts), x = meta$x, y = rep(0:1, 10)), cpus = cpus)
  expect_true(inherits(mem, 'ecodive_dbmem') && nrow(mem) == 20L && ncol(mem) >= 1L &&
              all(is.finite(mem)) && all(attr(mem, 'moran_i') > 0))



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
  expect_true(n_cpus() <= max_cpus())
  bt <- blas_threads()
  expect_true(is.integer(bt) && length(bt) == 1)

#})
