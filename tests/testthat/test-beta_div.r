test_that("beta diversity", {
  
  
  
  # beta_div wrapper =========================================================
  
  expect_equal(beta_div(counts, 'manhattan'), manhattan(counts))
  expect_equal(metrics$beta$manhattan(counts), manhattan(counts))
  
  
  
  # Bray-Curtis / Sorenson ====================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'bray', FALSE))
    object   = as.vector(bray_curtis(counts, TRUE)), 
    expected = c(0.354838709677419, 0.352941176470588, 0.642857142857143, 
                 0.0769230769230769, 0.212121212121212, 0.222222222222222 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'bray', TRUE))
    object   = as.vector(bray_curtis(counts, FALSE)), 
    expected = c(0.333333333333333, 0.333333333333333, 0.6, 0, 0.2, 0.2) )
  
  expect_equal(sorenson(counts), bray_curtis(counts, FALSE))
  
  
  
  # Canberra ==================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'canberra', FALSE))
    object   = as.vector(canberra(counts, TRUE)), 
    expected = c(0.522727272727273, 0.519230769230769, 0.75, 
                 0.0751633986928105, 0.37037037037037, 0.350877192982456 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'canberra', TRUE))
    object   = as.vector(canberra(counts, FALSE)), 
    expected = c(0.5, 0.5, 0.75, 0, 0.333333333333333, 0.333333333333333) )
  
  
  
  # Euclidean =================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'euclidean', FALSE))
    object   = as.vector(euclidean(counts, weighted = 1)), 
    expected = c(8.30662386291807, 9.2736184954957, 11.8321595661992, 
                 2.23606797749979, 5.3851648071345, 7.07106781186548 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'euclidean', TRUE))
    object   = as.vector(euclidean(counts, weighted = 0)), 
    expected = c(1.4142135623731, 1.4142135623731, 1.73205080756888, 0, 1, 1) )
  
  
  
  # Gower =====================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'gower', FALSE))
    object   = as.vector(gower(counts, TRUE)), 
    expected = c(0.388571428571429, 0.408571428571429, 0.571428571428571, 
                 0.0771428571428571, 0.182857142857143, 0.22 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'gower', TRUE))
    object   = as.vector(gower(counts, FALSE, cpus = 1)), 
    expected = c(0.4, 0.4, 0.6, 0, 0.2, 0.2) )
  
  
  
  # Jaccard ===================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'jaccard', FALSE))
    object   = as.vector(jaccard(counts, TRUE)), 
    expected = c(0.523809523809524, 0.521739130434783, 0.782608695652174, 
                 0.142857142857143, 0.35, 0.363636363636364 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'jaccard', TRUE))
    object   = as.vector(jaccard(counts, FALSE)), 
    expected = c(0.5, 0.5, 0.75, 0, 0.333333333333333, 0.333333333333333) )
  
  
  
  # Kulcynski =================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'kulczynski', FALSE))
    object   = as.vector(kulczynski(counts, TRUE)), 
    expected = c(0.337606837606838, 0.315018315018315, 0.641025641025641, 
                 0.0714285714285714, 0.205555555555556, 0.2 ))
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'kulczynski', TRUE))
    object   = as.vector(kulczynski(counts, FALSE)), 
    expected = c(0.333333333333333, 0.333333333333333, 0.583333333333333, 
                 0, 0.166666666666667, 0.166666666666667 ))
  
  
  
  # Manhattan =================================================================
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'manhattan', FALSE))
    object   = as.vector(manhattan(counts, TRUE)), 
    expected = c(11, 12, 18, 3, 7, 8) )
  
  expect_equal( # as.vector(vegan::vegdist(t(counts), 'manhattan', TRUE))
    object   = as.vector(manhattan(counts, FALSE)), 
    expected = c(2, 2, 3, 0, 1, 1) )
  
  
  
  
  # # Hand columns to an abdiv function in pairs.
  # pairwise <- function (m, f, ...)
  #   apply(combn(seq_len(ncol(m)), 2), 2L, \(i) f(m[,i[1]], m[,i[2]], ...))
  
  
  # Unweighted UniFrac ========================================================
  
  expect_equal( # pairwise(counts, abdiv::unweighted_unifrac, tree)
    object   = as.vector(unweighted_unifrac(counts, tree)), 
    expected = c(0.426927101499826, 0.426927101499826, 0.577607254970352, 
                 0, 0.197170241898676, 0.197170241898676 ))
  
  
  
  # Weighted UniFrac ========================================================
  
  expect_equal( # pairwise(counts, abdiv::weighted_unifrac, tree)
    object   = as.vector(weighted_unifrac(counts, tree)), 
    expected = c(0.682012820512821, 0.641941391941392, 1.05051282051282, 
                 0.0899920634920635, 0.438388888888889, 0.528380952380952 ))
  
  
  
  # Normalized UniFrac ========================================================
  
  expect_equal( # pairwise(counts, abdiv::weighted_normalized_unifrac, tree)
    object   = as.vector(weighted_normalized_unifrac(counts, tree)), 
    expected = c(0.392778473738226, 0.378863215786652, 0.542965436810857, 
                 0.0479898763749635, 0.207216197053649, 0.254811004455059 ))
  
  
  
  # Generalized UniFrac ========================================================
  
  expect_equal( # pairwise(counts, abdiv::generalized_unifrac, tree)
    object   = as.vector(generalized_unifrac(counts, tree)), 
    expected = c(0.461633430665443, 0.452969548713512, 0.601617618207151, 
                 0.0514229205275909, 0.248010203981353, 0.294478054449811 ))
  
  
  
  # Variance Adjusted UniFrac =================================================
  
  expect_equal( # pairwise(counts, abdiv::variance_adjusted_unifrac, tree)
    object   = as.vector(variance_adjusted_unifrac(counts, tree)), 
    expected = c(0.432212413424998, 0.422907225603088, 0.565799257415752, 
                 0.0488081800681702, 0.213158661751972, 0.259184248897931 ))
  
  
  
  # Matrix with > 100 columns to trigger pthreading ===========================
  
  expect_silent(bray_curtis(big_mtx))
  expect_silent(unweighted_unifrac(big_mtx, tree))
  
})
