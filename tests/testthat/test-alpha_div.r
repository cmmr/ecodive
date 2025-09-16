test_that("alpha diversity", {
  
  
  
  # alpha_div wrapper ====
  
  expect_equal(alpha_div(counts, 'observed'), observed(counts))
  expect_equal(alpha_div(counts, 'otus'),     observed(counts))
  
  
  
  # Matrix with > 100 columns to trigger pthreading ====
  
  expect_silent(simpson(big_mtx))
  expect_silent(faith(big_mtx, tree))
  
  
  
  
  # ACE ====
  
  expect_equal( # vegan::estimateR(counts)['S.ACE',]
    object   = ace(counts, cutoff = 5), 
    expected = c(A = 3, B = 3, C = 3, D = 2) )
  
  expect_error(ace(counts, cutoff = 'A'))  # stopifnot(is.numeric(cutoff))
  expect_error(ace(counts, cutoff = NULL)) # stopifnot(length(cutoff) == 1)
  expect_error(ace(counts, cutoff = NA))   # stopifnot(!is.na(cutoff))
  expect_error(ace(counts, cutoff = 0))    # stopifnot(cutoff > 0)
  expect_error(ace(counts, cutoff = 2.3))  # stopifnot(cutoff %% 1 == 0)
  
  
  
  # Berger ====
  
  expect_equal( # apply(counts, 1, tabula::index_berger)
    object   = berger(counts), 
    expected = c(A = 0.461538461538462, B = 0.444444444444444, 
                 C = 0.428571428571429, D = 0.666666666666667 ))
  
  
  
  # Brillouin ====
  
  expect_equal( # apply(counts, 1, tabula::index_brillouin)
    object   = brillouin(counts), 
    expected = c(A = 0.807097978290102, B = 0.900881045540206, 
                 C = 0.91741230069276,  D = 0.533824471198889 ))
  
  
  
  # Chao1 ====
  
  expect_equal(
    object   = chao1(counts), 
    expected = c(A = 3, B = NaN, C = NaN,  D = NaN) )
  
  expect_equal(chao1(t(1:10)), 10.5)
  expect_equal(chao1(t(rep(1:10, 2))), 21)
  
  
  
  # Faith's Phylogenetic Diversity ====
 
  expect_equal( # apply(counts, 1L, abdiv::faith_pd, tree)
    object   = faith(counts, tree), 
    expected = c(A = 2.319, B = 2.191, C = 2.191, D = 1.759) )
  
  
  
  # Fisher ====
  
  expect_equal( # vegan::fisher.alpha(counts)
    object   = fisher(counts, digits = 5), 
    expected = c(A = 1.22255, B = 1.02800, C = 0.95777, D = 0.61979) )
  
  expect_error(fisher(counts, digits = 'A'))  # stopifnot(is.numeric(digits))
  expect_error(fisher(counts, digits = NULL)) # stopifnot(length(digits) == 1)
  expect_error(fisher(counts, digits = NA))   # stopifnot(!is.na(digits))
  expect_error(fisher(counts, digits = -1))   # stopifnot(digits >= 0)
  expect_error(fisher(counts, digits = 100))  # stopifnot(digits <= 10)
  expect_error(fisher(counts, digits = 2.3))  # stopifnot(digits %% 1 == 0)
  
  
  
  # Inverse Simpson ====
 
  expect_equal( # vegan::diversity(counts, 'invsimpson')
    object   = inv_simpson(counts), 
    expected = c(A = 2.6, B = 2.84210526315789, 
                 C = 2.84516129032258, D = 1.8 ))
  
  
  
  # Margalef ====
  
  expect_equal( # apply(counts, 1, tabula::index_margalef)
    object   = margalef(counts), 
    expected = c(A = 0.77974249050256,  B = 0.691952512522387, 
                 C = 0.656917477506102, D = 0.369269373068855 ))
  
  
  
  # McIntosh ====
  
  expect_equal( # apply(counts, 1, tabula::index_mcintosh)
    object   = mcintosh(counts), 
    expected = c(A = 0.525602129138804, B = 0.532291232744964, 
                 C = 0.520794263652204, D = 0.34327800805431 ))
  
  
  
  # Menhinick ====
  
  expect_equal( # apply(counts, 1, tabula::index_menhinick)
    object   = menhinick(counts), 
    expected = c(A = 0.832050294337844, B = 0.707106781186548, 
                 C = 0.654653670707977, D = 0.516397779494322 ))
  
  
  
  # Observed Features ====
  
  expect_equal( # rowSums(counts > 0)
    object   = observed(counts), 
    expected = c(A = 3, B = 3, C = 3, D = 2) )
  
  
  
  # Shannon ====
  
  expect_equal( # vegan::diversity(counts, 'shannon')
    object   = shannon(counts), 
    expected = c(A = 1.01233083910317, B = 1.07204334357507,
                 C = 1.07101854240991, D = 0.636514168294813 ))
  
  
  
  # Simpson ====
  
  expect_equal( # vegan::diversity(counts, 'simpson')
    object   = simpson(counts), 
    expected = c(A = 0.615384615384615, B = 0.648148148148148, 
                 C = 0.648526077097506, D = 0.444444444444444 ))
  
  
  
  # Squares Estimator ====
  
  expect_equal( # apply(counts, 1, tabula::index_squares)
    object   = squares(counts), 
    expected = c(A = 3, B = 3, C = 3, D = 2))
  
})
