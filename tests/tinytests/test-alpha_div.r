#test_that("alpha diversity", {
  
  
  
  # alpha_div wrapper ====
  
  expect_equal(alpha_div(counts, 'observed'), observed(counts))
  expect_equal(alpha_div(counts, 'otus'),     observed(counts))
  
  
  
  # Matrix with > 100 columns to trigger pthreading ====
  
  expect_silent(simpson(big_mtx))
  expect_silent(faith(big_mtx, tree))
  
  
  
  
  # ACE ====
  
  expect_equal( # vegan::estimateR(counts)['S.ACE',]
    current = ace(counts, cutoff = 5), 
    target  = c(A = 3, B = 3, C = 3, D = 2) )
  
  expect_error(ace(counts, cutoff = 'A'))  # stopifnot(is.numeric(cutoff))
  expect_error(ace(counts, cutoff = NULL)) # stopifnot(length(cutoff) == 1)
  expect_error(ace(counts, cutoff = NA))   # stopifnot(!is.na(cutoff))
  expect_error(ace(counts, cutoff = 0))    # stopifnot(cutoff > 0)
  expect_error(ace(counts, cutoff = 2.3))  # stopifnot(cutoff %% 1 == 0)
  
  
  
  # Berger ====
  
  expect_equal( # apply(counts, 1, tabula::index_berger)
    current = berger(counts), 
    target  = c(A = 0.461538461538462, B = 0.444444444444444, 
                 C = 0.428571428571429, D = 0.666666666666667 ))
  
  
  
  # Brillouin ====
  
  expect_equal( # apply(counts, 1, tabula::index_brillouin)
    current = brillouin(counts), 
    target  = c(A = 0.807097978290102, B = 0.900881045540206, 
                 C = 0.91741230069276,  D = 0.533824471198889 ))
  
  
  
  # Chao1 ====
  
  # Bias-corrected by default, as scikit-bio.
  expect_equal(
    current = chao1(counts), 
    target  = c(A = 3, B = 3, C = 3,  D = 2) )
  expect_equal(chao1(counts, bias_corrected = FALSE), chao1(counts))
  
  expect_equal(chao1(t(1:10)), 10)
  expect_equal(chao1(t(1:10), bias_corrected = FALSE), 10.5)
  expect_equal(chao1(t(rep(1:10, 2))), 20 + 2 / 6)
  expect_equal(chao1(t(rep(1:10, 2)), bias_corrected = FALSE), 21)
  
  # The classic form falls back on the bias-corrected one without
  # singletons or doubletons (formerly Inf and NaN).
  expect_identical(chao1(t(c(1, 1, 5)), bias_corrected = FALSE), 3 + 1)
  expect_identical(chao1(t(c(3, 4, 5)), bias_corrected = FALSE), 3)
  expect_error(chao1(counts, bias_corrected = NA), '`bias_corrected` must be TRUE or FALSE.')
  expect_error(chao1(counts, bias_corrected = 1),  '`bias_corrected` must be TRUE or FALSE.')
  expect_identical(
    alpha_div(counts_int, 'chao1', bias_corrected = FALSE), 
    chao1(counts_int, bias_corrected = FALSE) )
  
  
  
  # Faith's Phylogenetic Diversity ====
 
  expect_equal( # apply(counts, 1L, abdiv::faith_pd, tree)
    current = faith(counts, tree), 
    target  = c(A = 2.319, B = 2.191, C = 2.191, D = 1.759) )
  
  
  
  # Fisher ====
  
  expect_equal( # vegan::fisher.alpha(counts)
    current = fisher(counts, digits = 5), 
    target  = c(A = 1.22255, B = 1.02800, C = 0.95777, D = 0.61979) )
  
  expect_error(fisher(counts, digits = 'A'))  # stopifnot(is.numeric(digits))
  expect_error(fisher(counts, digits = NULL)) # stopifnot(length(digits) == 1)
  expect_error(fisher(counts, digits = NA))   # stopifnot(!is.na(digits))
  expect_error(fisher(counts, digits = -1))   # stopifnot(digits >= 0)
  expect_error(fisher(counts, digits = 100))  # stopifnot(digits <= 10)
  expect_error(fisher(counts, digits = 2.3))  # stopifnot(digits %% 1 == 0)
  
  
  
  # Inverse Simpson ====
 
  expect_equal( # vegan::diversity(counts, 'invsimpson')
    current = inv_simpson(counts), 
    target  = c(A = 2.6, B = 2.84210526315789, 
                 C = 2.84516129032258, D = 1.8 ))
  
  
  
  # Margalef ====
  
  expect_equal( # apply(counts, 1, tabula::index_margalef)
    current = margalef(counts), 
    target  = c(A = 0.77974249050256,  B = 0.691952512522387, 
                 C = 0.656917477506102, D = 0.369269373068855 ))
  
  
  
  # McIntosh ====
  
  expect_equal( # apply(counts, 1, tabula::index_mcintosh)
    current = mcintosh(counts), 
    target  = c(A = 0.525602129138804, B = 0.532291232744964, 
                 C = 0.520794263652204, D = 0.34327800805431 ))
  
  
  
  # Menhinick ====
  
  expect_equal( # apply(counts, 1, tabula::index_menhinick)
    current = menhinick(counts), 
    target  = c(A = 0.832050294337844, B = 0.707106781186548, 
                 C = 0.654653670707977, D = 0.516397779494322 ))
  
  
  
  # Observed Features ====
  
  expect_equal( # rowSums(counts > 0)
    current = observed(counts), 
    target  = c(A = 3, B = 3, C = 3, D = 2) )
  
  
  
  # Shannon ====
  
  expect_equal( # vegan::diversity(counts, 'shannon')
    current = shannon(counts), 
    target  = c(A = 1.01233083910317, B = 1.07204334357507,
                 C = 1.07101854240991, D = 0.636514168294813 ))
  
  
  
  # Simpson ====
  
  expect_equal( # vegan::diversity(counts, 'simpson')
    current = simpson(counts), 
    target  = c(A = 0.615384615384615, B = 0.648148148148148, 
                 C = 0.648526077097506, D = 0.444444444444444 ))
  
  
  
  # Squares Estimator ====
  
  expect_equal( # apply(counts, 1, tabula::index_squares)
    current = squares(counts), 
    target  = c(A = 3, B = 3, C = 3, D = 2))
  
#})



  # Metrics ported from scikit-bio (plans/alpha_metrics.md) ====
  #
  # Fixtures: data-raw/expected/alpha_metrics.r and .py (scikit-bio 0.7.4,
  # vegan 2.7.8 output only).

  fx <- readRDS(file.path('expected', 'alpha_metrics.rds'))

  # ecodive's value of one fixture row's metric, for a whole table.
  adiv_of <- function (x, metric, q, ...) {
    if (metric == 'chao1_classic') return (chao1(x, bias_corrected = FALSE, ...))
    if (is.na(q))                  return (alpha_div(x, metric, ...))
    alpha_div(x, metric, q = q, ...)
  }

  rel_err <- function (a, b) {
    d <- abs(a - b)
    ifelse(d == 0, 0, d / pmax(abs(b), 1e-300))
  }

  has_matrix <- requireNamespace('Matrix', quietly = TRUE)
  has_slam   <- requireNamespace('slam',   quietly = TRUE)

  sk  <- fx$skbio
  key <- unique(sk[c('table', 'metric', 'q')])

  for (i in seq_len(nrow(key))) {

    k    <- key[i, ]
    same_q <- (is.na(sk$q) & is.na(k$q)) | (!is.na(sk$q) & sk$q %in% k$q)
    rows   <- sk[sk$table == k$table & sk$metric == k$metric & same_q, ]
    x    <- fx$tables[[k$table]]
    info <- paste(k$table, k$metric, k$q)
    got  <- unname(adiv_of(x, k$metric, k$q))[rows$sample]
    nz   <- rowSums(x)[rows$sample] > 0

    # An empty sample is NA; scikit-bio gives nan (0 for chao1 and gini).
    expect_true(all(is.na(got[!nz])), info = info)
    expect_true(all(rel_err(got[nz], rows$value[nz]) < 1e-12), info = info)

    # The same from samples in columns, and from sparse tables.
    expect_identical(unname(adiv_of(t(x), k$metric, k$q, margin = 2L))[rows$sample], got, info = info)
    if (has_matrix)
      expect_equal(unname(adiv_of(Matrix::Matrix(x, sparse = TRUE), k$metric, k$q))[rows$sample], got, tolerance = 1e-14, info = info)
    if (has_slam)
      expect_equal(unname(adiv_of(slam::as.simple_triplet_matrix(x), k$metric, k$q))[rows$sample], got, tolerance = 1e-14, info = info)
  }

  # vegan (output only): Hill numbers, Renyi entropies, Pielou's J and
  # Simpson's dominance. vegan's Pielou is NaN for one feature; ecodive's
  # is 1, as scikit-bio's.
  vg <- fx$vegan
  for (i in seq_len(nrow(vg))) {
    got <- unname(adiv_of(fx$tables[[vg$table[i]]], vg$metric[i], vg$q[i]))[vg$sample[i]]
    if (is.nan(vg$value[i])) expect_identical(got, 1, info = vg$metric[i])
    else expect_true(rel_err(got, vg$value[i]) < 1e-12, info = paste(vg$table[i], vg$metric[i], vg$q[i]))
  }


  # Identities. ====

  x <- fx$tables$rand
  x <- x[rowSums(x) > 0, ]
  expect_identical(hill(x, q = 0), observed(x))
  expect_equal(hill(x, q = 1),   exp(shannon(x)),  tolerance = 1e-14)
  expect_equal(hill(x, q = 2),   inv_simpson(x),   tolerance = 1e-14)
  expect_equal(hill(x, q = Inf), 1 / berger(x),    tolerance = 1e-14)
  for (q in c(0, 0.5, 1, 2, 3, Inf)) {
    expect_equal(renyi(x, q = q), log(hill(x, q = q)), tolerance = 1e-13, info = q)
  }
  expect_equal(tsallis(x, q = 2),  simpson(x),               tolerance = 1e-14)
  expect_equal(tsallis(x, q = 0),  observed(x) - 1,          tolerance = 1e-14)
  expect_equal(tsallis(x, q = 1),  shannon(x),               tolerance = 1e-14)
  expect_equal(dominance(x),       1 - simpson(x),           tolerance = 1e-14)
  expect_equal(simpson_evenness(x), inv_simpson(x) / observed(x), tolerance = 1e-14)
  expect_equal(pielou(x),          shannon(x) / log(observed(x)), tolerance = 1e-14)
  expect_equal(heip(x),            (exp(shannon(x)) - 1) / (observed(x) - 1), tolerance = 1e-14)
  expect_equal(goods_coverage(x),  1 - rowSums(x == 1) / rowSums(x), tolerance = 1e-14)

  # The defaults: q = 2, as scikit-bio.
  expect_identical(hill(x),    hill(x, q = 2))
  expect_identical(renyi(x),   renyi(x, q = 2))
  expect_identical(tsallis(x), tsallis(x, q = 2))

  # The documented base R equivalents.
  r <- ecodive::ex_counts[1, ]
  p <- r[r > 0] / sum(r)
  expect_equal(unname(gini(ecodive::ex_counts)[1]), max(0, 1 - 2 * sum(cumsum(sort(r)) / sum(r)) / length(r)))
  s <- sort(r[r > 0], decreasing = TRUE)
  expect_equal(unname(strong(ecodive::ex_counts)[1]), max(cumsum(s) / sum(s) - seq_along(s) / length(s)))
  n <- sum(r > 0)
  expect_equal(unname(mcintosh_evenness(ecodive::ex_counts)[1]), sqrt(sum(r^2)) / sqrt((sum(r) - n + 1)^2 + n - 1))
  expect_equal(unname(hill(ecodive::ex_counts, q = 0.5)[1]), sum(p ** 0.5) ** (1 / (1 - 0.5)))


  # Proportions give the same values as counts, except for the metrics that
  # need integers. ====

  xp <- x / rowSums(x)
  for (m in c('pielou', 'heip', 'simpson_evenness', 'dominance', 'strong', 'gini', 'hill', 'renyi', 'tsallis'))
    expect_equal(alpha_div(xp, m), alpha_div(x, m), tolerance = 1e-13, info = m)
  for (m in c('goods_coverage', 'mcintosh_evenness'))
    expect_error(alpha_div(xp, m), 'integer', info = m)


  # Edge cases. ====

  e <- fx$tables$edge
  # One feature: pielou and heip 1, renyi and tsallis 0, hill 1.
  expect_identical(pielou(e)[[1]], 1)
  expect_identical(heip(e)[[1]], 1)
  for (q in c(0, 0.5, 1, 2, Inf)) {
    expect_identical(renyi(e, q = q)[[1]], 0, info = q)
    expect_identical(tsallis(e, q = q)[[1]], 0, info = q)
    expect_equal(hill(e, q = q)[[1]], 1, info = q)
  }
  # All singletons: Good's coverage 0.
  expect_identical(goods_coverage(e)[[2]], 0)
  # An empty sample is NA for every metric.
  for (m in list_metrics('alpha', 'id'))
    if (m != 'faith') expect_true(is.na(alpha_div(e, m)[[3]]), info = m)
  # tsallis(q = Inf) is 0.
  expect_true(all(tsallis(x, q = Inf) == 0))
  # The Gini index changes with absent features; the others don't.
  w <- fx$tables$edge_wide
  expect_true(all(gini(w)[-3] > gini(e)[-3]))
  for (m in c('pielou', 'heip', 'strong', 'goods_coverage', 'dominance', 'hill'))
    expect_identical(alpha_div(w, m), alpha_div(e, m), info = m)
  # Evenly spread over every feature: Gini 0, Strong 0, evenness 1.
  expect_identical(gini(e)[[8]], 0)
  expect_equal(strong(e)[[8]], 0)
  expect_equal(pielou(e)[[8]], 1)
  expect_equal(simpson_evenness(e)[[8]], 1)
  # One feature of K: Gini 1 - 2/K (by rectangles).
  expect_equal(gini(e)[[1]], 1 - 2/8)
  expect_equal(gini(w)[[1]], 1 - 2/20)
  # A single sample as a vector.
  expect_identical(pielou(c(1, 2, 3)), pielou(t(c(1, 2, 3))))
  expect_identical(gini(c(0, 2, 3)), gini(t(c(0, 2, 3))))


  # q. ====

  q_msg <- '`q` must be a single number of at least 0 (Inf is allowed).'
  for (f in list(hill, renyi, tsallis)) {
    expect_error(f(counts, q = -1),     q_msg, fixed = TRUE)
    expect_error(f(counts, q = NA),     q_msg, fixed = TRUE)
    expect_error(f(counts, q = c(1, 2)), q_msg, fixed = TRUE)
    expect_error(f(counts, q = 'a'),    q_msg, fixed = TRUE)
    expect_error(f(counts, q = NULL),   q_msg, fixed = TRUE)
    expect_silent(f(counts, q = Inf))
    expect_identical(f(counts, q = 1L), f(counts, q = 1))
  }
  expect_error(alpha_div(counts, 'hill', q = -1), q_msg, fixed = TRUE)
  expect_identical(alpha_div(counts, 'hill', q = 0.5), hill(counts, q = 0.5))
  expect_identical(alpha_div(counts, 'renyi', q = 3),  renyi(counts, q = 3))
  expect_identical(alpha_div(counts, 'tsallis', q = Inf), tsallis(counts, q = Inf))


  # alpha_div(): IDs and aliases. ====

  ids <- c(
    'dominance', 'gini', 'goods_coverage', 'heip', 'hill', 'mcintosh_evenness',
    'pielou', 'renyi', 'simpson_evenness', 'strong', 'tsallis' )
  expect_true(all(ids %in% list_metrics('alpha', 'id')))
  for (m in ids)
    expect_identical(alpha_div(counts_int, m), do.call(m, list(counts_int)), info = m)

  aliases <- c(
    pielou_e = 'pielou', pielou_j = 'pielou', heip_e = 'heip',
    simpson_e = 'simpson_evenness', mcintosh_e = 'mcintosh_evenness',
    goods = 'goods_coverage', good = 'goods_coverage', hill_number = 'hill',
    simpson_d = 'dominance', gini_index = 'gini', 'Gini Index' = 'gini',
    sobs = 'observed', observed_features = 'observed', enspie = 'inv_simpson',
    berger_parker_d = 'berger', mcintosh_d = 'mcintosh', brillouin_d = 'brillouin',
    fisher_alpha = 'fisher', 'Renyi Entropy' = 'renyi', renyi_entropy = 'renyi' )
  for (a in names(aliases))
    expect_identical(match_metric(a, div = 'alpha')$id, aliases[[a]], info = a)

  # 'gini' is the Gini index; the Gini-Simpson index is 'simpson'.
  expect_identical(match_metric('gini')$id, 'gini')
  expect_identical(match_metric('Gini-Simpson')$id, 'simpson')

  # One metric per call, always a named vector.
  expect_identical(alpha_div(counts_int, 'shannon'), shannon(counts_int))
  expect_identical(alpha_div(counts_int, 'hill', q = 1), hill(counts_int, q = 1))
  expect_error(alpha_div(counts_int, c('shannon', 'pielou')), '`metric` must be length 1', fixed = TRUE)
  expect_error(alpha_div(counts_int, character(0)), '`metric` must be length 1', fixed = TRUE)
  expect_error(alpha_div(counts_int, 1), '`metric` must be a character vector', fixed = TRUE)


  # Results are identical for any cpus (and big enough to run threads). ====

  set.seed(3)
  big <- matrix(rnbinom(300 * 60, mu = 10, size = 0.4), 300, 60)
  for (m in c(ids, 'chao1')) {
    ref <- alpha_div(big, m, cpus = 1)
    for (cpus in c(2L, 4L, 7L))
      expect_identical(alpha_div(big, m, cpus = cpus), ref, info = paste(m, cpus))
  }


  # alpha_div_rarefied() passes q and bias_corrected, one metric per call.
  # ====

  x  <- fx$tables$rand
  dd <- 30L
  a  <- suppressWarnings(alpha_div_rarefied(x, 'hill', q = 1, depth = dd, times = 3, seed = 2))
  r  <- rarefy(x, dd, seed = 2, times = 3, warn = FALSE)
  expect_equal(as.vector(a), Reduce(`+`, lapply(r, hill, q = 1)) / 3)
  a0 <- suppressWarnings(alpha_div_rarefied(x, 'chao1', bias_corrected = FALSE, depth = dd, times = 3, seed = 2))
  expect_equal(as.vector(a0), Reduce(`+`, lapply(r, chao1, bias_corrected = FALSE)) / 3)
  expect_error(alpha_div_rarefied(x, 'hill', q = -1, depth = dd), q_msg, fixed = TRUE)

  at <- suppressWarnings(alpha_div_rarefied(x, 'tsallis', q = 0.5, depth = dd, times = 3, seed = 2, sd = TRUE))
  expect_true(is.numeric(at) && is.null(dim(at)))
  expect_identical(names(attr(at, "sd")), names(at))
  expect_equal(as.vector(at), Reduce(`+`, lapply(r, tsallis, q = 0.5)) / 3)
  expect_error(
    suppressWarnings(alpha_div_rarefied(x, c('shannon', 'hill'), depth = dd)),
    '`metric` must be length 1', fixed = TRUE )
