test_that("beta diversity", {
  
  
  # # Hand columns to an abdiv function in pairs.
  # abdiv        <- function (x, metric, ...) apply(combn(seq_len(nrow(x)), 2), 2L, \(i) do.call(`::`, list('abdiv', metric))(x[i[1],], x[i[2],], ...))
  # parallelDist <- function (x, metric, ...) as.vector(parallelDist::parallelDist(x, metric, ...))
  # philentropy  <- function (x, metric, ...) suppressMessages(as.vector(as.dist(philentropy::distance(x, metric, ...))))
  # vegan        <- function (x, metric, ...) as.vector(vegan::vegdist(x, metric, ...))
  
  
  
  
  # beta_div wrapper ====
  
  expect_equal(beta_div(counts, 'ochiai'), ochiai(counts))
  expect_equal(beta_div(counts, 'Otsuka'), ochiai(counts))
  
  
  
  # Matrix with > 100 columns to trigger pthreading ====
  
  expect_silent(bray(big_mtx))
  expect_silent(unweighted_unifrac(big_mtx, tree))
  
  
  
  # Pairs != NULL ====
  
  expect_equal(
    object   = as.vector(bray(counts, pairs = 1:2)), 
    expected = c(0.354838709677419, 0.352941176470588, NA, NA, NA, NA) )
  
  expect_equal(
    object   = as.vector(unweighted_unifrac(counts, tree, pairs = 1:2)), 
    expected = c(0.426927101499826, 0.426927101499826, NA, NA, NA, NA) )
  
  
  
  # Pairs == integer(0) ====
  
  expect_equal(
    object   = as.vector(bray(counts, pairs = integer(0))), 
    expected = as.numeric(c(NA, NA, NA, NA, NA, NA)) )
  
  expect_equal(
    object   = as.vector(unweighted_unifrac(counts, tree, pairs = integer(0))), 
    expected = as.numeric(c(NA, NA, NA, NA, NA, NA)) )
  
  
  
  # Aitchison ====
  
  expect_equal( # vegan(counts + 1, 'aitchison')
    object   = as.vector(aitchison(counts, pseudocount = 1)), 
    expected = c(2.42489292341729,  2.48372587230057, 3.26493363948397, 
                 0.250928909965064, 1.6566104027154,  1.88303098843395 ))
  
  expect_equal( # clr euclidean == aitchison
    object   = as.vector(euclidean(counts, norm = 'clr')), 
    expected = as.vector(aitchison(counts)) )
  
  
  
  # Bhattacharyya ====
  
  expect_equal( # philentropy(counts_p, 'bhattacharyya')
    object   = as.vector(bhattacharyya(counts)), 
    expected = c(0.378456649039499,   0.364064988015199, 1.02706186684777, 
                 0.00210387145241093, 0.164142171328222, 0.203046152841888 ))
  
  
  
  # Bray-Curtis ====
  
  expect_equal( # vegan(counts, 'bray')
    object   = as.vector(bray(counts)), 
    expected = c(0.354838709677419,  0.352941176470588, 0.642857142857143, 
                 0.0769230769230769, 0.212121212121212, 0.222222222222222 ))
  
  expect_equal( # binary bray == sorensen
    object   = as.vector(bray(counts, norm = 'binary')), 
    expected = as.vector(sorensen(counts)) )
  
  
  
  # Canberra ====
  
  expect_equal( # philentropy(counts, 'canberra')
    object   = as.vector(canberra(counts)), 
    expected = c(2.09090909090909,  2.07692307692308, 3, 
                 0.225490196078431, 1.11111111111111, 1.05263157894737 ))
  
  
  
  # Chebyshev ====
  
  expect_equal( # philentropy(counts, 'chebyshev')
    object   = as.vector(chebyshev(counts)), 
    expected = c(8, 9, 10, 2, 5, 7) )
  
  
  
  # Chord ====
  
  expect_equal( # abdiv(counts, 'chord')
    object   = as.vector(chord(counts)), 
    expected = c(0.849787680450186, 0.815473426840873, 1.2022062234803, 
                 0.11819770141174,  0.490727607919993, 0.589602596121889 ))
  
  
  
  # Clark ====
  
  expect_equal( # philentropy(counts, 'clark')
    object   = as.vector(clark(counts)), 
    expected = c(1.41713247891999,  1.416304049194,   1.73205080756888, 
                 0.176742709581568, 1.00615390423749, 1.00138408370739 ))
  
  
  
  # Divergence ====
  
  expect_equal( # philentropy(counts_p, 'divergence')
    object   = as.vector(divergence(counts)), 
    expected = c(4.17558822616231,   4.16275577991495, 6.01020408163265, 
                 0.0290244021712553, 2.09652892561983, 2.15007351396765 ))
  
  
  
  # Euclidean ====
  
  expect_equal( # vegan(counts, 'euclidean')
    object   = as.vector(euclidean(counts)), 
    expected = c(8.30662386291807, 9.2736184954957, 11.8321595661992, 
                 2.23606797749979, 5.3851648071345,  7.07106781186548 ))
  
  
  
  # Generalized UniFrac ====
  
  expect_equal( # abdiv(counts, 'generalized_unifrac', tree)
    object   = as.vector(generalized_unifrac(counts, tree)), 
    expected = c(0.461633430665443,  0.452969548713512, 0.601617618207151, 
                 0.0514229205275909, 0.248010203981353, 0.294478054449811 ))
  
  
  
  # Gower ====
  
  expect_equal( # vegan(counts, 'gower'); cluster::daisy(counts, 'gower')
    object   = as.vector(gower(counts)), 
    expected = c(0.388571428571429,  0.408571428571429, 0.571428571428571, 
                 0.0771428571428571, 0.182857142857143, 0.22))
  
  
  
  # Hellinger ====
  
  expect_equal( # vegan(counts, 'hellinger')
    object   = as.vector(hellinger(counts)), 
    expected = c(0.793829122355658,  0.781222072639055, 1.13308654830978, 
                 0.0648330142150295, 0.550233834516,    0.606233340003081 ))
  
  
  
  # Hamming ====
  
  expect_equal( # abdiv(counts > 0, 'hamming')
    object   = as.vector(hamming(counts)), 
    expected = c(2, 2, 3, 0, 1, 1) )
  
  
  
  # Horn ====
  
  expect_equal( # vegan(counts, 'horn')
    object   = as.vector(horn(counts)), 
    expected = c(0.361702127659574,   0.333175355450237, 0.727272727272727, 
                 0.00698549167114448, 0.142857142857143, 0.195 ))
  
  
  
  # Jaccard ====
  
  expect_equal( # vegan(counts, 'jaccard', binary = TRUE)
    object   = as.vector(jaccard(counts)), 
    expected = c(0.5, 0.5, 0.75, 0, 0.333333333333333, 0.333333333333333) )
  
  
  
  # Jensen-Shannon distance ====
  
  expect_equal( # sqrt(philentropy(counts_p, 'jensen-shannon'))
    object   = as.vector(jensen(counts)), 
    expected = c(0.472459308330185,  0.464481366702569, 0.667264300672983, 
                 0.0458301797899792, 0.329728993801882, 0.364081336931884 ))
  
  
  
  # Jensen-Shannon Divergence (JSD) ====
  
  expect_equal( # philentropy(counts_p, 'jensen-shannon')
    object   = as.vector(jsd(counts)), 
    expected = c(0.223217798027837,   0.215742940013886, 0.445241646952604, 
                 0.00210040537958182, 0.108721209353602, 0.132555219902108 ))
  
  
  
  # Lorentzian ====
  
  expect_equal( # philentropy(counts, 'lorentzian')
    object   = as.vector(lorentzian(counts)), 
    expected = c(3.98898404656427, 4.0943445622221,  5.44241771052179, 
                 1.79175946922805, 2.89037175789616, 2.77258872223978))
  
  
  
  # Manhattan ====
  
  expect_equal( # vegan(counts, 'manhattan')
    object   = as.vector(manhattan(counts)), 
    expected = c(11, 12, 18, 3, 7, 8) )
  
  
  
  # Matusita ====
  
  expect_equal( # philentropy(counts_p, 'matusita')
    object   = as.vector(matusita(counts)), 
    expected = c(0.793829122355658,  0.781222072639055, 1.13308654830978, 
                 0.0648330142150295, 0.550233834516,    0.606233340003081 ))
  
  
  
  # Minkowski ====
  
  expect_equal( # philentropy(counts, 'minkowski', p = 1.5)
    object   = as.vector(minkowski(counts, power = 1.5)), 
    expected = c(8.87866660953074, 9.83182738119092, 13.417493822195, 
                 2.44726081477148, 5.81121051366182,  7.24976150371949 ))
  
  
  
  # Morisita ====
  
  expect_equal( # abdiv(counts, 'morisita')
    object   = as.vector(morisita(counts)), 
    expected = c( 0.273504273504273, 0.247613700168445,  0.700854700854701, 
                  -0.103733215287,    0.0713489409141583, 0.133709981167608 ))
  
  
  
  # Motyka ====
  
  expect_equal( # philentropy(counts, 'motyka')
    object   = as.vector(motyka(counts)), 
    expected = c(0.67741935483871,  0.676470588235294, 0.821428571428571, 
                 0.538461538461538, 0.606060606060606, 0.611111111111111))
  
  
  
  # Normalized UniFrac ====
  
  expect_equal( # abdiv(counts, 'weighted_normalized_unifrac', tree)
    object   = as.vector(normalized_unifrac(counts, tree)), 
    expected = c(0.392778473738226,  0.378863215786652, 0.542965436810857, 
                 0.0479898763749635, 0.207216197053649, 0.254811004455059 ))
  
  
  
  # Otsuka-Ochiai ====
  
  expect_equal( # parallelDist(counts, 'ochiai')
    object   = as.vector(ochiai(counts)), 
    expected = c(0.333333333333333, 0.333333333333333, 0.591751709536137, 
                 0,                 0.183503419072274, 0.183503419072274 ))
  
  
  
  # Probabilistic Symmetric Chi-Squared ====
  
  expect_equal( # philentropy(counts_p, 'prob_symm')
    object   = as.vector(psym_chisq(counts)), 
    expected = c(1.32239418236062,   1.27514231499051,  2.57142857142857, 
                 0.0167832167832168, 0.654545454545455, 0.801932367149758 ))
  
  
  
  # Soergel ====
  
  expect_equal( # philentropy(counts, 'soergel')
    object   = as.vector(soergel(counts)), 
    expected = c(0.523809523809524, 0.521739130434783, 0.782608695652174, 
                 0.142857142857143, 0.35,              0.363636363636364 ))
  
  expect_equal( # binary soergel == jaccard
    object   = as.vector(soergel(counts, norm = 'binary')), 
    expected = as.vector(jaccard(counts)) )
  
  
  
  # Dice-Sorensen ====
  
  expect_equal( # abdiv(counts, 'sorenson')
    object   = as.vector(sorensen(counts)), 
    expected = c(0.333333333333333, 0.333333333333333, 0.6, 0, 0.2, 0.2) )
  
  
  
  # Squared Chi-Squared ====
  
  expect_equal( # philentropy(counts_p, 'squared_chi')
    object   = as.vector(squared_chisq(counts)), 
    expected = c(0.66119709118031, 0.637571157495256, 1.28571428571429, 
                 0.00839160839160839, 0.327272727272727, 0.400966183574879 ))
  
  
  
  # Squared Chord ====
  
  expect_equal( # philentropy(counts_p, 'squared_chord')
    object   = as.vector(squared_chord(counts)), 
    expected = c(0.630164675499954,  0.610307926778461, 1.28388512596057, 
                 0.0042033197322062, 0.302757272646181, 0.367518862531291 ))
  
  
  
  # Squared Euclidean ====
  
  expect_equal( # philentropy(counts, 'squared_euclidean')
    object   = as.vector(squared_euclidean(counts)), 
    expected = c(69, 86, 140, 5, 29, 50) )
  
  
  
  # Topsoe ====
  
  expect_equal( # philentropy(nz_counts_p, 'topsoe')
    object   = as.vector(topsoe(nz_counts)), 
    expected = c(0.198353107654462,   0.197099053749219,  0.416832493477321, 
                 0.00318836636509382, 0.0936295595052999, 0.12396818587543 ))
  
  
  
  # Unweighted UniFrac ====
  
  expect_equal( # abdiv(counts, 'unweighted_unifrac', tree)
    object   = as.vector(unweighted_unifrac(counts, tree)), 
    expected = c(0.426927101499826, 0.426927101499826, 0.577607254970352, 
                 0,                 0.197170241898676, 0.197170241898676 ))
  
  
  
  # Variance Adjusted UniFrac ====
  
  expect_equal( # abdiv(counts, 'variance_adjusted_unifrac', tree)
    object   = as.vector(variance_adjusted_unifrac(counts, tree)), 
    expected = c(0.432212413424998,  0.422907225603088, 0.565799257415752, 
                 0.0488081800681702, 0.213158661751972, 0.259184248897931 ))
  
  
  
  # Wave Hedges ====
  
  expect_equal( # philentropy(counts, 'wavehedges')
    object   = as.vector(wave_hedges(counts)), 
    expected = c(2.16666666666667,  2.14285714285714, 3, 
                 0.396825396825397, 1.2,              1.1 ))
  
  
  
  # Weighted UniFrac ====
  
  expect_equal( # abdiv(counts, 'weighted_unifrac', tree)
    object   = as.vector(weighted_unifrac(counts, tree)), 
    expected = c(0.682012820512821,  0.641941391941392, 1.05051282051282, 
                 0.0899920634920635, 0.438388888888889, 0.528380952380952 ))
  
})
