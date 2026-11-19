# ecodive 2.4.0

* ecodive now covers a whole diversity analysis: besides alpha and beta
  diversity, it ordinates samples and tests what explains their distances
  (see below). New vignettes, `vignette("ordination")` and
  `vignette("stat_tests")`, walk through both, and the README and
  `vignette("ecodive")` give an overview of everything the package does. The
  package title is now "Parallel and Memory-Efficient Ecological Diversity
  Analysis".

* New `plot()` methods draw ecodive's results with base graphics (only
  'graphics' and 'grDevices', from R 3.6.0 on), taking their titles, axis
  labels, subtitles and legends from what each result holds, and computing
  nothing new:
  * Ordinations (`pcoa()`, `pca()`, `dbrda()`, `nmds()`, `procrustes()`, ...)
    draw two axes, labelled with their share of the variance, with the
    method's fit in the subtitle (stress, R-squared, ...). `plot(ord,
    "Species", iris)` colours the samples by a column of metadata, with a
    legend, by the same `groups` and `data` as the tests; supervised
    ordinations colour by their own groups. Loadings, constraint arrows,
    centroids and a Procrustes rotation's target are drawn on top.
  * `plot(envfit(...), add = TRUE)` adds the fitted variables to an
    ordination plot.
  * Permutation tests show the permutation distribution of their statistic,
    with the observed value marked; `mrpp()` and `permdisp()` show each
    group's within-group distance or dispersion, and `procrustes_test()` its
    Procrustes fit (each with `type = "null"` for the distribution).
  * `simper()` gives a dot chart of the taxa that contribute most to each
    pair's dissimilarity, `dirmult_ttest()` a volcano plot, `ancom()` a dot
    chart of W, `varpart()` a Venn diagram of the fractions (two to four
    sets), `forward_select()` the adjusted R-squared at each step, and
    `bioenv()` the best correlation for each number of variables.
  * Distance matrices give a heatmap, in the samples' order or sorted by
    `groups`, with a colour key: those of `beta_div()` and every metric
    function, of `vlr()`, `rho_p()` and `phi_s()` (over features), and
    `rpca()`'s `dm`. These now have class `c("ecodive_dist", "dist")`, so
    that `plot()` knows them; base R has no `plot()` method for `dist`
    objects, and every function for `dist` objects works on them as
    before. Code that checks `identical(class(d), "dist")` or
    `class(d) == "dist"` should use `inherits(d, "dist")` instead.
  * `ancombc()` and `ancombc2()` give a volcano plot of one covariate's log
    fold changes, `ancombc_pairwise()` and `ancombc_dunnett()` a heatmap of
    the log fold change of each feature in each comparison, and
    `ancombc_global()` and `ancombc_trend()` a dot chart of their
    statistics.

* Every function that returns a distance matrix, `beta_div()` and its 37
  metric functions, `vlr()`, `rho_p()` and `phi_s()`, takes `clust`, to
  cluster the samples before returning it: one of `stats::hclust()`'s
  linkage methods (`"average"`, `"complete"`, `"single"`, `"ward.D"`,
  `"ward.D2"`, `"mcquitty"`, `"median"`, `"centroid"`), or `"pcoa"` to order
  them along the first axis of `pcoa()`. Clustering **reorders the
  samples** of the result, with their labels (unnamed samples are labelled
  with their input positions); the matrix's `dist` record gives each one's
  position in the input as `order`, and the method as `clust`. ecodive's
  ordinations and tests put the samples back in their input order before
  analysing them, so metadata lines up whether matched by name or by
  position; for other packages, reorder metadata with
  `attr(d, "dist")$order`. `clust = NULL`, the default, leaves the
  matrix as it was, and `plot()` then clusters it by average linkage for
  the heatmap; `clust = FALSE` (recorded as `"none"`) keeps the input order
  in the plot too.

* Results record how they were made, one attribute per step, and the
  records travel with them through later steps (`?ecodive_steps`). A
  distance matrix has a `dist` record: the metric's display name (such as
  `"Bray-Curtis Dissimilarity"`, so that other packages can label it
  without ecodive's table of metrics), its ID, and the `norm`,
  `pseudocount`, `power`, `alpha` and `clust` it was made with. An
  ordination adds an `ord` record (its method, display name and settings,
  such as `n_axes` or `seed`) and keeps its distance matrix's `dist`; the
  permutation tests keep their distance matrix's `dist`, and `envfit()` its
  ordination's `dist` and `ord`. Each record is a small list that `str()`
  and `print()` show on one line, read as `attr(x, "dist")$metric`.
  `plot()` titles with them, such as "Principal Coordinates Analysis of
  Bray-Curtis Dissimilarity". `vlr()`, `rho_p()` and `phi_s()` move their
  `method`, `denominator` and `pseudocount` attributes into their `dist`
  record (`method` as `metric`).

* New `pcoa()` runs a principal coordinates analysis (classical
  multidimensional scaling) on a `dist` object, such as from `beta_div()`. It
  gives the same coordinates and eigenvalues as `cmdscale()`, `ape::pcoa()`,
  and `vegan::wcmdscale()` (up to each axis's sign), but computes only the
  first few axes (`n_axes = 3` by default), in parallel C, without forming an
  n-by-n matrix. Its results are identical for any number of `cpus`.
  `correction = "lingoes"` or `"cailliez"` adds the smallest constant to the
  squared distances, or to the distances, that leaves no negative eigenvalue,
  giving the same constants, eigenvalues and coordinates as `ape::pcoa()`,
  `vegan::wcmdscale(add = )` and `cmdscale(add = TRUE)`, again without an
  n-by-n (or 2n-by-2n) matrix: with 2,000 samples, `"cailliez"` takes 20
  seconds where `cmdscale(add = TRUE)` takes 2 minutes. Every result has
  `correction` and `constant` attributes.

* New `pcoa_biplot()` projects the features of a table (taxa counts, or
  environmental variables) onto the axes of a `pcoa()` result, for drawing
  arrows on a PCoA plot. It computes the same scores as scikit-bio's
  `pcoa_biplot()`, in parallel C, keeps sparse tables sparse, and takes every
  table type and `norm` that `pca()` does. `envfit()` gives the same
  directions with permutation tests.

* Every ordination function takes the number of axes to return as `n_axes`
  (3 by default for the methods whose axes have eigenvalues: `pcoa()`,
  `pca()`, `rpca()`, `dbrda()`, `capscale()` and `isomap()`; 2 for `mds()`,
  `nmds()`, `tsne()` and `umap()`). The name says what it counts; `k`, as
  `cmdscale()` and vegan call it, often means a number of neighbours in
  these methods.

* New `dist_ord()` runs any distance-based ordination method by name, as
  `beta_div()` does for beta diversity metrics and `dist_test()` does for
  statistical tests, with the method's own arguments: `dist_ord(dm, "nmds")`
  is `nmds(dm)`. `list_dist_ords()` and `match_dist_ord()` list and look up
  the methods. `count_ord()`, with `list_count_ords()` and
  `match_count_ord()`, does the same for the methods that take a count table
  rather than distances, `pca()`, `rda()`, `cca()`, `ca()`, `lda()`, `plsda()`, `oplsda()`, `rpca()` and `rpca_phylo()`: `count_ord(counts, "pca")` is
  `pca(counts)`. Each wrapper points to the other when given the other's
  method, and to `dist_test()` when given a test.
  Every ordination method returns the same kind of result: a sample-by-axis
  matrix (of class `"ecodive_ord"`, which prints its attributes' names
  rather than their contents) whose attributes start with `method`, `eig`,
  and `prop`. Every method rejects negative distances (beyond rounding
  error), and `seed` can be given to any of them through `dist_ord()`.

* New `nmds()` runs non-metric multidimensional scaling with the SMACOF
  algorithm of scikit-learn's non-metric `MDS()`, ported to parallel C. It
  starts from the PCoA solution and from `n_starts` random configurations,
  keeps the one with the lowest Kruskal stress-1, and centers it, rotates it
  to principal axes, and scales it to the distances. Tied distances use
  Kruskal's primary approach (as in `vegan::metaMDS()`), or scikit-learn's
  secondary one with `ties = "secondary"`. The distances are sorted once and
  shared by all starts; nothing n-by-n is formed. Each start runs at most
  `max_iters` iterations, or stops when the stress improves by less than
  `tol`. Random starts come from a `seed` argument, and results are
  identical for any number of `cpus`. The result's attributes are `stress`,
  `converged`, `n_iters`, `n_starts`, `best_start`, `n_best`, and `ties`.

* New `mds()` runs metric multidimensional scaling with the SMACOF algorithm
  of scikit-learn's metric `MDS()`, fitting the distances themselves rather
  than their ranks. It shares `nmds()`'s engine, starts, arguments (without
  `ties`) and result format, and has no limit on the number of samples.
  `dist_ord(dm, "mds")` and `"metric MDS"` select it; `"classical MDS"`,
  `"pco"` and `"cmdscale"` select `pcoa()`.

* New `isomap()` runs an Isomap (isometric feature mapping) ordination: it
  links each sample to its `n_neighbors` nearest neighbors (or to those
  within a `radius`), replaces each distance with the shortest path through
  that graph, and runs `pcoa()` on those geodesic distances, which keeps long
  ecological gradients straight where `pcoa()` bends them into an arch. A C
  port of scikit-learn's `Isomap(metric = "precomputed")`, it gives the same
  geodesic distances and embedding to rounding, breaks ties by sample order,
  and joins disconnected neighbor graphs at their closest samples, with a
  warning (`vegan::isomap()` stops instead; note that its `k` is the number
  of neighbors, while `isomap()`'s `n_axes` is the number of axes). The neighbor
  and shortest-path searches run in parallel without an n-by-n matrix, and
  results are identical for any number of `cpus`.

* New `tsne()` runs t-distributed stochastic neighbor embedding (t-SNE) on a
  `dist` object, as a port of scikit-learn's `TSNE` to parallel C.
  `gradient = "auto"` uses the exact gradient for up to 500 samples (or more
  than 3 dimensions) and the Barnes-Hut approximation (a quadtree in 2
  dimensions, an octree in 3) above that. It finds nearest neighbors exactly
  from the distances and never forms an n-by-n matrix. It starts from the
  PCoA solution by default, sets the learning rate with scikit-learn's rule
  (`learning_rate = "auto"`), stops early once the layout stops improving
  (at most `max_iters` iterations; the number run is the `n_iters`
  attribute), and lowers `perplexity` for small datasets. Its layouts match
  scikit-learn's and `Rtsne::Rtsne()`'s in quality. Random starts come from
  a `seed` argument, and results are identical for any number of `cpus`.

* New `umap()` runs a uniform manifold approximation and projection (UMAP)
  on a `dist` object, following `uwot::umap()` and its defaults: exact
  nearest neighbours, smooth-kNN calibration and fuzzy union (identical to
  uwot's), a spectral start from the normalized graph Laplacian, and the
  layout optimization of uwot's batch mode (`batch = TRUE`: per-epoch
  gradients, then Adam), which gives uwot's coordinates from the same start
  and seeds. Every stage runs in parallel C, the neighbours straight from
  the dist triangle without an n-by-n matrix. Random numbers come from a
  `seed` argument, and results are identical for any number of `cpus`.
  `min_dist` defaults to the canonical 0.1 (umap-learn's), not uwot's 0.01.

* New `permanova()` runs a permutational multivariate analysis of variance
  on a `dist` object and sample metadata. Like ecodive's other tests, it
  takes `groups`: one value per sample (categories or a numeric covariate),
  or a column name of `data`. For several variables, give `data` and a
  `formula` instead: factors, covariates, interactions, and `Condition()`
  terms. It gives the same table as `vegan::adonis2()` (Df, sums of squares,
  R2, F, and p-values), testing the model as a whole or each term
  sequentially, marginally, or by contrast, with optional `strata`. Given the
  same permutation matrix, its permuted F statistics and p-values are
  identical to `adonis2()`'s. It returns a list of class `"ecodive_test"`,
  like `mrpp()`: the model's F, p-value, z-score (as
  `summary(vegan::permustats())`) and number of permutations come first,
  whatever `by` is, then its R2 and `adonis2()`'s table, as `table`. The
  sums of squares come straight from the dist triangle, without the n-by-n
  matrices `adonis2()` builds for each permutation, and permutations run in
  parallel. On one thread it is 25 to 50 times faster than `adonis2()` with
  200 samples, and 100 to 300 times faster with 1,000.
  Permutations come from a `seed` argument, and results are identical for
  any number of `cpus`. Like `adonis2()`, it also takes the distances on the
  left-hand side of a formula: `permanova(dm ~ treatment + ph, data)`.

* New `mrpp()` runs a multi-response permutation procedure on a `dist`
  object and each sample's group, as `vegan::mrpp()`: the observed and
  expected weighted mean within-group distance (delta), the chance-corrected
  within-group agreement A, each group's mean distance, and a permutation
  p-value, with group weights `'n'`, `'n-1'` or `'n(n-1)'` (vegan's
  `weight.type` 1 to 3) and optional `strata`. It also gives the
  standardized effect size z of `summary(vegan::permustats())`. Given the
  same permutation matrix, its permuted deltas, p-values and z are those of
  `vegan::mrpp()`. Each permutation costs one comparison per pair of
  samples, read straight from the dist triangle; with 500 samples and 999
  permutations it is about 60 times faster than `vegan::mrpp()` on one
  thread. Results are identical for any number of `cpus`. Like ecodive's
  other permutation tests, it returns a list of class `"ecodive_test"`
  whose `statistic`, `p_value`, `z` and `n_perm` come first, and
  `as.data.frame()` turns it into one row with the same columns for every
  test.

* New `anosim()` runs an analysis of similarities (ANOSIM) on a `dist`
  object and each sample's group, given as a vector or as a column of sample
  metadata, with optional `strata`. Given the same permutation matrix, its R
  statistic, permuted R values and p-value match `vegan::anosim()`'s. It
  ranks the distances once, with a radix sort, and then each permutation
  needs only one pass over the ranks, without the n-by-n matrices
  `vegan::anosim()` builds; permutations run in parallel. With 300 samples
  and 999 permutations it takes under a tenth of a second, against about
  seven for `vegan::anosim()`. Results are identical for any number of `cpus`.

* New `permdisp()` tests the homogeneity of multivariate dispersions
  (PERMDISP), as `vegan::betadisper()` followed by `vegan::permutest()`: each
  sample's distance to its group's centroid (the default) or spatial
  median, each group's mean distance, and the ANOVA F with a permutation
  p-value, plus optional pairwise t tests and small-group bias adjustment.
  Distances to centroids come from the distances within groups alone, with
  no eigendecomposition: they match `betadisper()`'s to rounding error, and
  are about 100 times faster with 1,000 samples. `betadisper()` defaults to
  spatial medians instead: pass `type = "median"` to compare with it.
  `permdisp()` solves them exactly, where `betadisper()`'s optimizer stops
  at about 1e-5. Given the same permutation matrix, its p-values are
  identical to `permutest()`'s. Permutations come from a `seed` argument,
  can be restricted to `strata`, and results are identical for any number
  of `cpus`.

* New `mantel()` runs a Mantel test of the correlation between two distance
  matrices (Pearson, Spearman, or Kendall), and new `mantel_partial()` a
  partial Mantel test, which removes the effect of a third (`given`). They
  give the same statistics as `vegan::mantel()` and `vegan::mantel.partial()`,
  and, given the same permutation matrix, the same permuted statistics and
  p-values. The second and third matrices may also be numeric vectors, one
  value per sample (such as age), and samples are matched by name. They add
  lower-tail and two-sided tests (`alternative`). Each permutation is one pass
  over the dist triangle, with no n-by-n matrix (Kendall's tau-b by Knight's
  merge sort rather than `cor()`'s quadratic count), and permutations run in
  parallel. Permutations come from a `seed` argument, and results are
  identical for any number of `cpus`.

* New `dist_test(dm, test, groups, data, ...)` runs any of these tests by
  name, as `dist_ord()` does for ordination methods: `dist_test(dm, "mrpp",
  "treatment", meta)` is `mrpp(dm, "treatment", meta)`, and other arguments
  (such as `mantel()`'s `method`) pass through. It checks that the variable is
  of a kind the test takes: categorical groups for `anosim`, `mrpp` and
  `permdisp`, numeric values for `mantel` and `mantel_partial` (passed as
  their `y`), and either for `permanova`. `list_dist_tests()` and
  `match_dist_test()` list and look up the tests, with forgiving name matching
  (`"adonis2"` finds `permanova`, `"betadisper"` finds `permdisp`).

* `permanova()`, `anosim()`, `mrpp()` and `permdisp()` take `pairwise = TRUE`
  to also compare each pair of groups, as is usual after an overall test
  across several groups, and `p_adjust` to adjust the pairs' p-values for
  multiple testing: any of `stats::p.adjust.methods`, Holm's method by
  default. `permanova()`, `anosim()` and `mrpp()` rerun the test on each
  pair's samples, as q2-diversity's `beta_group_significance(pairwise=True)`
  and pairwiseAdonis's `pairwise.adonis2()` do, with the pair's statistics,
  z-score, sample size, and raw and adjusted p-values in a `pairwise` data
  frame. Given the same permutations, the pairs' statistics and p-values are
  identical to vegan's tests on each pair. Each pair's random permutations
  come from `seed` and the pair, so results are identical for any number of
  `cpus`; a supplied `permutations` matrix is cut down to each pair's
  samples. A pair that can't be tested, such as two groups of one sample,
  gives `NA`s with a warning. `permdisp()`'s pairwise t tests gain
  `p_adjusted` and `p_adjusted_perm` (each of its p-value columns adjusted
  across the pairs, which `vegan::permutest()` doesn't do), and `z`, `n` and
  `n_perm` columns, so its pairwise table now has more columns than before.

* One rule aligns every per-sample argument (`groups`, `strata`, and the
  Mantel tests' `y` and `given`) in `permanova()`, `mrpp()`, `anosim()`,
  `permdisp()`, `mantel()`, `mantel_partial()`, `envfit()` and `dist_test()`:
  with `data`, a column name or one value per row of `data`; without it, one
  value per sample, matched by name when named, else by position. The same
  checks and messages apply everywhere. As with `vegan::adonis2()`,
  `permanova()` needs `groups` or a `formula` (use `formula = ~ .` to test
  every column of `data`). `permutations` is limited to 1e8.

* New `dbrda()` runs a distance-based redundancy analysis (db-RDA): the axes
  of a `dist` object that sample metadata explain, given as a formula such as
  `~ treatment + Condition(site)`. It gives the same eigenvalues, inertia,
  site scores, biplot scores and centroids as `vegan::dbrda()` (up to each
  axis's sign; not `vegan::capscale()`, which differs for non-Euclidean
  distances), without forming an n-by-n matrix: the constrained axes come from
  a small eigenproblem, and unconstrained axes (`MDS1`, ...) are only computed
  when the constraints give fewer than `n_axes` axes. Its results are identical for
  any number of `cpus`. As in `vegan`, `dbrda(dm ~ treatment, data)` works
  too, and a formula of only `Condition()` terms gives the unconstrained
  ordination after removing their effect.

* `dbrda()` takes the same `groups`, `data` and `formula` as `permanova()`:
  one variable as `groups` (a vector or a column of `data`), or `data` and a
  `formula`, or `dm ~ terms`, so `dist_ord(dm, "dbrda", "treatment", meta)`
  works. Its results gain `r_squared` and `adj_r_squared`, as
  `vegan::RsquareAdj()`, including for partial models.

* New `capscale()` runs a constrained analysis of principal coordinates. It
  gives the same results as `vegan::capscale()` (and
  `phyloseq::ordinate(method = "CAP")`), with any imaginary inertia from
  negative eigenvalues reported in `inertia`. Euclidean distances take
  `dbrda()`'s fast path; otherwise a threaded eigendecomposition makes it
  about five times faster than `vegan::capscale()` on one thread. Also
  available as `dist_ord(dm, "capscale", ...)`.

* New `rda()` runs a redundancy analysis of a count table: the axes of the
  (normalized) table that sample metadata explain (`RDA1`, ...), then
  principal components of the rest (`PC1`, ...). `norm = 'hellinger'` or
  `'chord'` gives transformation-based RDA (tb-RDA), and `Condition()` terms
  partial RDA. It takes `groups`, `data` and `formula` as `dbrda()` does
  (and `counts ~ terms`, as `vegan::rda()`), and with no model it is a PCA.
  It gives the same site scores, species scores (`loadings`), eigenvalues,
  inertia, R-squared, biplot scores and centroids as `vegan::rda()` (up to
  each axis's sign), sharing `dbrda()`'s solver and `pca()`'s sparse table
  products, so the table is never made dense: with 3,000 samples and 2,000
  taxa it takes half a second where `vegan::rda()` takes a minute. Results
  are identical for any number of `cpus`. Also available as
  `count_ord(counts, "rda", ...)`. Its permutation test is `permanova()` of
  Euclidean distances of the normalized table.

* New `cca()` runs a canonical correspondence analysis (ter Braak 1986) of a
  table of abundances: the correspondence analysis axes, in the chi-square
  metric, that sample metadata explain (`CCA1`, ...), then correspondence
  analysis axes of the rest (`CA1`, ...). It takes `groups`, `data` and
  `formula` as `rda()` does, including `Condition()` for partial CCA and
  `counts ~ terms`, and with no model it is a correspondence analysis. It
  gives the same site scores, species scores (`loadings`), eigenvalues,
  inertia, R-squared, biplot scores and centroids as `vegan::cca()` (up to
  each axis's sign), and the same permutation-based adjusted R-squared as
  `vegan::RsquareAdj()` for the same permutations (`permutations` and
  `seed`), partial CCA included. It shares `rda()`'s solver, with the
  chi-square weights, and multiplies the table's non-zero values only, so a
  sparse table is never made dense: with 3,000 samples and 2,000 taxa it
  takes half a second where `vegan::cca()` takes a minute. Results are
  identical for any number of `cpus`. Also available as `count_ord(counts,
  "cca", ...)`.

* New `ca()` runs a correspondence analysis of a table of abundances: the
  unconstrained ordination in the chi-square metric, with species scores
  (`loadings`). It gives the same site scores, species scores, eigenvalues
  and inertia as `vegan::ca()` (`vegan::cca()` without constraints),
  `ade4::dudi.coa()` and scikit-bio's `ca()` (up to each axis's sign), and
  is the same as `cca()` with no model, whose sparse solver it shares.
  Results are identical for any number of `cpus`. Also available as
  `count_ord(counts, "ca")`; `"CA"` selects it, and `dist_ord(dm, "ca")`
  points to `count_ord()`.

* New `envfit()` fits sample metadata (`data`) and feature tables
  (`counts`, such as taxa abundances) onto any ordination, for biplots: an
  arrow and r-squared for each numeric variable or feature, the centroid of
  each factor level, and permutation p-values. It gives the same arrows,
  r-squared and centroids as `vegan::envfit()`, and the same p-values for the
  same permutation matrix. Each permutation costs a few multiply-adds per
  non-zero value of the table, with no refitting, and sparse tables are never
  made dense, so thousands of taxa with 999 permutations take well under a
  second (minutes in `vegan`). It fits every axis of the ordination by
  default (`vegan` fits the first two), and its results are identical for any
  number of `cpus`.

* New `procrustes()` rotates, reflects, and optionally dilates one
  ordination to fit another of the same samples, and reports the Procrustes
  sum of squares, m^2: how similar, say, Bray-Curtis and UniFrac PCoAs (or
  16S and metagenomic ones) of the same samples are. It gives the same
  rotation, translation, dilation, residuals and m^2 as `vegan::procrustes()`
  (and the same symmetric m^2 as SciPy's `scipy.spatial.procrustes()`), but
  takes its arguments as `procrustes(ord, target)`: `vegan::procrustes(target,
  ord)`. vegan's `scale` is called `dilate`. Inputs may be `dist` objects,
  ordinated with `pcoa()` first (`n_axes` axes), and samples are matched by
  name. The result is the rotated ordination, with the rest as attributes.
  New `procrustes_test()` runs the Procrustes randomization test (PROTEST),
  giving the same statistic as `vegan::protest()` and the same p-values for
  the same permutation matrix, in parallel C, with `strata`, and identical
  results for any number of `cpus`.

* New `plsda()` runs a partial least squares discriminant analysis
  (PLS-DA): the components of a count table that best separate known
  groups, from PLS regression on a dummy matrix of the groups (ported from
  scikit-learn's `PLSRegression`). It gives the same scores, weights,
  loadings and Y loadings as scikit-learn, and the same scores, explained
  variance and VIP as `mixOmics::plsda()` (up to each component's sign),
  and also returns each feature's VIP and each group's centroid. Each
  component's weights are computed directly rather than by NIPALS
  iteration, so they don't depend on a convergence tolerance. It takes
  `groups` and `data` as `permanova()` does, every `norm` (`'clr'` is the
  usual choice for microbiome counts), and sparse tables, which it never
  makes dense: 1,000 samples and 5,000 taxa take under a tenth of a second,
  where `mixOmics::plsda()` takes over a second. Results are identical for
  any number of `cpus`. `count_ord(counts, 'plsda', groups)` runs it too.
  Like every supervised method, it separates any groups it is given, so
  validate the separation (e.g. with `permanova()`).

* New `oplsda()` runs an orthogonal PLS discriminant analysis (OPLS-DA) of
  two groups: one predictive component (`p1`) that separates them and
  orthogonal components (`o1`, `o2`, ...) for the rest of the table's
  variation (ported from pyopls). It gives the same orthogonal scores,
  weights and loadings as pyopls, and the same scores, loadings, weights,
  R2X, R2Y and VIP as `ropls::opls()` (up to each component's sign), and
  also returns each group's centroid. It takes `groups`, `data`, every
  `norm` and every table type as `plsda()` does, and never makes a sparse
  table dense: 1,000 samples, 5,000 taxa and two orthogonal components
  take about a tenth of a second, where `ropls::opls()` takes 8 (including
  its cross-validation). Results are identical for any number of `cpus`.
  `count_ord(counts, 'oplsda', groups)` runs it too. For more than two
  groups, use `plsda()`.

* New `lda()` runs a linear discriminant analysis (Fisher's discriminant
  axes, ported from scikit-learn's `LinearDiscriminantAnalysis` with its
  `"svd"` solver). It gives the same scores, coefficients (`scaling`),
  group means and proportion of trace as `MASS::lda()` and its
  `predict()` (up to each axis's sign), with optional `prior`
  probabilities, and the same scores as scikit-learn up to a constant
  factor (scikit-learn divides the within-group covariance by n, not
  n - g). It takes `groups` and `data` as `plsda()` does, every `norm` and
  every table type, and stops with an error when the within-group
  covariance is singular, as it is for microbiome tables with more taxa
  than samples: run it on a few PCoA or PCA axes instead. 5,000 samples
  and 200 features take 0.2 seconds, where `MASS::lda()` takes 0.8.
  Results are identical for any number of `cpus`.
  `count_ord(counts, 'lda', groups)` runs it too.

* Ordinations of many axes are faster on several threads. The block
  Lanczos solver behind `pcoa()`, `pca()`, `dbrda()`, `capscale()`,
  `rda()`, `cca()`, `rpca()` and `umap()`'s spectral start now splits its
  own steps on the basis (products with small matrices, dot products and,
  for more than four axes, block Gram-Schmidt) among `cpus`, each as far
  as its work is worth it, as well as the products with the data. On 4
  threads, `pcoa()` of 50 axes takes 1.4 s instead of 2.5 for 1,000
  samples and 15 s instead of 27 for 5,000, and `dbrda()` of 40 axes for
  2,000 samples 3.3 s instead of 5.6. Results stay identical for any number
  of `cpus`; for up to four axes they are unchanged, and for more they
  differ from before only by rounding.

* New `capda()` runs Anderson & Willis's (2003) discriminant canonical
  analysis of principal coordinates (CAP, as in PRIMER-PERMANOVA+): `lda()`
  of the first m axes of `pcoa()`, for any distances, with m chosen by
  leave-one-out classification success (up to `max_pcoa_axes`, 50 by
  default; on ties, the smallest m), and a permutation test of the trace
  statistic (the sum of the squared canonical correlations), with `strata`
  or a supplied permutation matrix. It gives the same m, classification
  success, leave-one-out classes, squared canonical correlations and scores
  as `BiodiversityR::CAPdiscrim()`, whose leave-one-out posteriors underflow
  for samples far from every group, and which stops at a group of one
  sample (`capda()` counts such samples as misclassified). The success for
  every m considered and the leave-one-out table are returned. Results are
  identical for any number of `cpus`. `dist_ord(dm, 'capda', groups)` runs
  it too; `'CAP'` still selects `capscale()`, which is db-RDA. `lda()`'s
  error for a singular within-group covariance now suggests `capda()`.

* New `simper()` runs a similarity percentages analysis (SIMPER): it breaks
  the mean Bray-Curtis dissimilarity between each pair of groups down into
  each taxon's contribution, with its sd, ratio, the groups' mean
  abundances and the cumulative proportion, following ecopy's `simper()`,
  and gives each taxon a permutation p-value (Clarke 1993), adjusted within
  each pair by `p_adjust` (Holm's method by default). It gives the same
  results as `vegan::simper()` and `summary()` of it, as one data frame.
  Its p-values count permuted contributions at least as large as the
  observed one, as Clarke (1993) does; vegan counts only larger ones, which
  makes its p-values too small for small groups. As in the pairwise tests,
  each pair's permutations shuffle only its own samples (vegan shuffles all
  of them). It takes `groups` and `data` as `permanova()` does, `strata`,
  `norm` (`'none'` by default, as vegan) and sparse tables, which it never
  makes dense, and runs in parallel C: 300 samples and 3,000 taxa, three
  pairs with 99 permutations, take under 2 seconds on one thread, where
  `vegan::simper()` takes 27 seconds for one pair. Results are identical
  for any number of `cpus`. Its docs explain why SIMPER's ranking mixes up
  differences in mean and in variability (Warton et al. 2012).

* New `bioenv()` finds the subset of environmental variables whose distances
  among the samples best match the community distances (BIOENV; Clarke &
  Ainsworth 1993): for each subset size, the subset with the highest rank
  correlation with `dm`, as one data frame, with the best subset overall as
  attribute `best`. It follows scikit-bio's `bioenv()` and gives the same
  subsets and correlations as it and as `vegan::bioenv()`, with vegan's
  options: `method` (Spearman, Pearson or Kendall), `metric` (`'euclidean'`,
  `'mahalanobis'`, `'manhattan'` or `'gower'`), `max_vars` (vegan's
  `upto`) and `given` (vegan's `partial`, also with Kendall, which vegan
  refuses). Variables come from `data`, by `formula` or as its numeric
  columns, and samples are matched by name. Every subset is searched in
  parallel C, with no n-by-n matrix: 100 samples and 12 variables (4,095
  subsets) take 1.2 seconds on one thread, where `vegan::bioenv()` takes
  5.5. Results are identical for any number of `cpus`.

* New `varpart()` partitions the variation in a community among two to four
  sets of explanatory variables (Borcard et al. 1992): the fractions each
  set explains on its own, those they share, and the residuals, from the
  adjusted R-squared of (partial) RDAs of a count table (with `norm`, such
  as `"hellinger"`), db-RDAs of a `dist`, or, with `chisquare = TRUE`, CCAs
  with permutation-estimated adjusted R-squared (Peres-Neto et al. 2006). It
  gives the same fractions as `vegan::varpart()` with two, three and four
  sets, as one data frame named with vegan's letters (`[a]`, `[b]`, ...),
  with vegan's union and conditional fractions as attributes and a text
  Venn diagram when printed. Each set is a formula over `data`, a data
  frame, a matrix, or a vector, matched to the samples by name. For RDA and
  db-RDA, one pass over the data gives every union's R-squared, with no
  n-by-n matrix and sparse tables kept sparse: with 3,000 samples, 2,000
  taxa (5% non-zero) and three sets, it takes 0.13 seconds on 4 threads
  (Hellinger table or Bray-Curtis distances), where `vegan::varpart()`
  takes 3.5. Results are identical for any number of `cpus`.

* New `forward_select()` picks explanatory variables one at a time by
  adjusted R-squared, with a permutation test of each added variable and
  Blanchet et al.'s (2008) double stopping rule, for RDA of a count table
  (with `norm`) or db-RDA of a `dist`. It selects the same variables in the
  same order as `vegan::ordiR2step()`, with the same adjusted R-squared and
  F, and the same p-values given the same permutations. Candidates are the
  terms of a formula over `data` (factors enter whole, interactions after
  their main effects, `Condition()` terms are covariables). By default it
  first tests the global model and selects nothing if that is not
  significant; `max_vars` caps the number of variables, and `strata`
  restricts permutations. Every candidate's adjusted R-squared comes from
  one pass over the data, and only the chosen candidate of each step is
  tested: with 1,000 samples, 500 taxa, 11 candidates and 199
  permutations, it takes 0.17 seconds on 4 threads for the Hellinger table
  and 0.08 for Bray-Curtis distances, where `vegan::ordiR2step()` takes 58
  and 85. Results are identical for any number of `cpus`.

* New `alr()` and `ilr()`, the additive and isometric log-ratio transforms,
  complete the clr and rclr that `norm` already offers. `alr()` divides by a
  `denominator` feature that you choose (there is no default: scikit-bio
  uses the first feature and compositions the last). `ilr()` takes the
  Helmert basis that scikit-bio uses by default, a sequential binary
  partition (converted as scikit-bio's `sbp_basis()`), or any orthonormal
  basis, matched to the features by name. Both give the same coordinates as
  scikit-bio's `alr()` and `ilr()` (`alr()` to the last bit, `ilr()` to
  1e-15), take every table type with zeros handled by the shared
  `pseudocount` rule, and return a dense matrix in the input's orientation
  with `method`, `denominator` or `basis`, and `pseudocount` attributes,
  ready for `pca()`, `rda()`, `lda()` or `envfit()`. The default basis is
  applied by a running sum, without a dense product: 2,000 samples by
  5,000 taxa take 0.2 seconds, where `normalize(x, 'clr') %*% t(basis)`
  takes 30 with R's reference BLAS. Euclidean distances on `ilr()`
  coordinates are `aitchison()` distances, for any basis. Results are
  identical for any number of `cpus`.

* New `replace_zeros()` replaces the zeros of a count table with small
  positive values, so that log ratios (`norm = 'clr'`, `aitchison()`,
  `alr()`, `ilr()`) can be taken without a pseudocount: its result has no
  zeros, so they use none and give no warning. It runs any zero replacement
  method by name, as `count_ord()` does for ordinations, and
  `list_zero_replacements()` and `match_zero_replacement()` list and look
  up the methods. The first is `replace_zeros_mult()`, multiplicative
  replacement (Martín-Fernández et al. 2003): every zero becomes `delta`
  (one value, or one per sample; scikit-bio's `1 / D^2` by default) and
  each sample's non-zero proportions shrink by the same factor, keeping
  their ratios. Its results are bit-identical to scikit-bio's
  `multi_replace()`. Every method returns a dense matrix of proportions in
  the input's orientation, with `method`, its parameters and `n_replaced`
  attributes, and takes every table type: a sparse 2,000 by 5,000 table
  takes 0.2 seconds.

* New `log_ratio()` and `balance()` give one log ratio per sample between
  a `numerator` and a `denominator` group of features (names, indices or
  logical vectors), as a named vector that drops straight into metadata or
  a test. `log_ratio()` is the log of the ratio of summed abundances, as
  Qurro plots it (identical to Qurro's `qarcoal()`); it needs no
  pseudocount by default, and gives `NA`, with one warning, where a sum is
  0. `balance()` is the isometric log-ratio balance, the scaled log ratio
  of geometric means: bit-identical to the matching coordinate of `ilr()`
  with a sequential binary partition, and equal to scikit-bio's to 6e-16,
  without building the full basis; zeros follow the shared `pseudocount`
  rule over the whole table, and the output of `replace_zeros()` needs
  none. Both take every table type and keep sparse tables sparse: 2,000
  samples by 5,000 taxa take under 0.1 seconds.

* New `vlr()`, `rho_p()` and `phi_s()` measure proportionality between
  every pair of features, the compositionally valid alternative to
  correlating relative abundances (Lovell et al. 2015). `vlr()` is the
  variance of the log ratio, equal to scikit-bio's `pairwise_vlr()` to
  4e-15 but computed from differences rather than a covariance matrix, so
  nearly proportional pairs keep their precision. `rho_p()` is Erb &
  Notredame's proportionality correlation (propr's `"rho"`) and `phi_s()`
  the symmetric phi of Quinn et al. (propr's `"phs"`), against the clr by
  default or a `denominator` (one feature, an alr, or a set). Each returns
  a `dist` over the features, ready for `pcoa()` or `hclust()` (`rho_p()`,
  a similarity, for thresholding). They take every table type and keep it
  sparse, run in parallel, and give identical results for any `cpus`: all
  4.5 million pairs of 3,000 taxa in 500 samples take about a second on 4
  cores.

* New `qmp_scale()` and `qmp_rarefy()` turn relative abundances into
  absolute abundances when each sample's total microbial `load` was measured
  (flow cytometry, qPCR, a spike-in). `qmp_scale()` multiplies each sample's
  relative abundances by its load (Jian et al. 2020); `qmp_rarefy()` is
  Vandeputte et al.'s (2017) quantitative microbiome profiling, which first
  rarefies every sample to the same number of reads per unit of load. It
  follows the authors' `QMP.R` (BSD-3, Raes Lab), with the same sampling
  depths and targets, but uses ecodive's own rarefaction and `seed`. `load`
  is matched to the samples by name, or given as a column of `data`. Both
  keep the table's class and orientation, and sparse tables stay sparse:
  2,000 samples by 5,000 taxa take under 0.1 seconds.

* New `ancom()` runs ANCOM (Mandal et al. 2015), the analysis of
  composition of microbiomes: it tests every log ratio between two features
  for a difference between `groups` (`sig_test = 'anova'`, or `'kruskal'`),
  adjusts each feature's p-values (`p_adjust`, Holm by default), and
  declares the features whose W, the number of their significant log
  ratios, passes scikit-bio's cutoff rule. Its `W` and decisions are
  identical to scikit-bio's `ancom()` in every reference case, including
  its quirks (each feature's family includes a 0 for itself; numpy's
  cutoffs); the Kruskal-Wallis statistics are bit-identical to SciPy's and
  the F statistics more accurate (SciPy's cancel to 3e-12). It returns a
  per-feature table (`feature`, `W`, `reject`). Zeros follow the shared
  `pseudocount` rule, and the output of `replace_zeros()` needs none. The
  tests run in parallel C straight from the (sparse) table, identical for
  any `cpus`: all 4.5 million log ratios of 3,000 features in 200 samples
  take about 3 seconds on 4 cores (10 with `'kruskal'`).

* New `dirmult_ttest()` finds differentially abundant features the ALDEx2
  way (Fernandes et al. 2013, 2014), as scikit-bio's `dirmult_ttest()`:
  each sample's proportions are drawn `n_draws` times from their Dirichlet
  posterior (`pseudocount = 0.5` is the prior), every draw is
  clr-transformed, each feature is compared between the two sides of a
  `contrast` with Welch's t-test, and the results are averaged. It returns
  a per-feature table with the average `t`, `log2_fc` and `p_value`, the
  conservative interval of scikit-bio (the widest bounds over the draws),
  `p_adjusted` and `reject`. `contrast = c(treatment, reference)` names the
  groups; `NA` as the reference pools all other samples, which is also the
  default (the first group against the rest). Given the same draws, the
  results are scikit-bio's to about 1e-13. An optional scale model
  (`scale_sd`, following ALDEx3's `clr.sm()`; Nixon et al. 2025) adds
  uncertainty about the groups' total loads to every draw, so that a
  change in load isn't mistaken for changes in the features; it is off by
  default (`scale_sd = 0`, scikit-bio's procedure), and works best with
  more draws (ALDEx3 suggests `scale_sd = 0.5`). The draws come from a new Gamma
  sampler on PCG (Marsaglia & Tsang), with one random stream per sample
  and draw, so results are identical for any `cpus`: 500 samples by 5,000
  features, 128 draws, take about 5 seconds on 4 cores.

* New `ancombc()` and `ancombc2()` run ANCOM-BC and ANCOM-BC2 (Lin &
  Peddada 2020, 2024): each feature's log abundance is fitted with a
  linear model (`formula` over `data`, any `model.matrix()` formula,
  treatment contrasts), and each covariate's bias from the samples'
  differing sampling fractions is estimated by an EM algorithm and removed,
  giving bias-corrected log fold changes with standard errors, p-values
  (normal for ANCOM-BC, t for ANCOM-BC2), `p_adjusted` (`p_adjust` within
  each covariate) and `reject`. `ancombc2()` also corrects each sample for
  its sampling fraction and refits, widens the variances (`var_quantile`),
  and by default treats zeros as missing values, including ANCOMBC's
  handling of failed fits and missing reference levels; `ancombc()` uses
  the shared `pseudocount` rule. They follow scikit-bio's `ancombc()` and
  `ancombc2()` (themselves derived from the ANCOMBC package, whose
  Artistic-2.0 notice they carry) and match them to about 1e-8 in every
  reference case, and ANCOMBC to about 1e-5 where its EM converges. The
  result is a per-feature-and-covariate table that also carries the fit
  (coefficients, variances, covariance matrices, degrees of freedom), for
  post-hoc tests to come. The fits run in parallel C, identical for any
  `cpus` and for sparse and dense tables: 500 samples by 5,000 features
  take about 1 second.

* New post-hoc tests on an `ancombc()` or `ancombc2()` fit, for a factor
  with three or more levels, as scikit-bio's `ANCOMBCResult` methods:
  `ancombc_global()` (does the feature differ between any levels?),
  `ancombc_pairwise()` (every pair of levels, with mixed directional FDR),
  `ancombc_dunnett()` (every level against the reference, with a
  bootstrapped mdFDR screen) and `ancombc_trend()` (ordered patterns,
  increasing or decreasing by default, or your own constraint matrices,
  with a bootstrapped p-value). The deterministic statistics and p-values
  match scikit-bio's to about 1e-8, and the bootstrapped ones match exactly
  on scikit-bio's own null draws; ecodive draws its own on PCG (`seed`),
  one stream per feature, identical for any `cpus`. Comparisons are
  labelled with level names (`'severe - mild'`), and any categorical term
  of the fit can be tested, not only its `groups`. The trend test's
  constrained estimates come from an exact QP solver written for ecodive.
  New `structural_zeros()` finds features absent from a group (Kaul et al.
  2017), as scikit-bio's `struc_zero()`, returning a features-by-groups
  logical matrix.

* `rarefy(depth =)` now also takes one depth per sample (matched to the
  samples by name when named). `drop` and `warn` then compare each sample
  with its own depth. A single depth gives exactly the same results as
  before.

* New `normalize()` returns a count table after any of the transformations
  that `beta_div()`, `pca()` and `envfit()` apply with their `norm` argument
  (`'percent'`, `'chord'`, `'hellinger'`, `'binary'`, `'clr'`, `'rclr'`),
  with exactly the values they use. It keeps the table's orientation,
  dimnames and form: sparse Matrix and slam input stays sparse.
  `norm = 'clr'` gives a dense matrix, since every zero gets a value; the
  functions with a `norm` argument avoid that copy.

* New `norm = 'hellinger'` option, wherever `norm` is offered (beta
  diversity metrics, `pca()`, `pcoa_biplot()`, `envfit()` and
  `normalize()`): the square root of each sample's relative abundances.
  Euclidean distances of the result are Hellinger distances, and
  `pca(counts, norm = 'hellinger')` is the usual Hellinger PCA. Negative
  values are an error.

* `simper()` and `plsda()` now name `counts`, not `dm`, when `groups`,
  `strata` or the rows of `data` don't match the number of samples.

* Method names are matched more reliably by `match_metric()`,
  `match_dist_ord()`, `match_dist_test()` and `match_count_ord()`, and so by
  `beta_div()`, `alpha_div()`, `dist_ord()`, `dist_test()` and `count_ord()`:
  a name that is exactly a method's ID or alias, in any case, now selects it
  even when it also begins another method's name. For example, `"Jensen"`
  selects `jensen` rather than stopping as ambiguous with `jsd`
  (Jensen-Shannon divergence), and `dist_ord(dm, "Mantel")` now says to use
  `dist_test()`.

* New `pca()` runs a principal component analysis of a count table, after
  any of `beta_div()`'s normalizations (e.g. `norm = 'clr'` for an Aitchison
  PCA), optionally scaling features to unit variance (`scale = TRUE`). It
  gives the same scores, loadings, and variances as `prcomp()` on the
  normalized table (up to each axis's sign), but computes only the first few
  axes, in parallel C, straight from sparse input without making it dense.
  Its scores equal `pcoa()`'s on the matching Euclidean or Aitchison
  distances; its `loadings` attribute holds the features' weights. Results
  are identical for any number of `cpus`.

* New `rpca()` runs a robust Aitchison PCA (Martino et al. 2019) straight
  from a count table: a robust centered log-ratio transform that treats zeros
  as missing, OptSpace matrix completion, and a PCA of the completed matrix.
  It returns sample scores, feature loadings, and the fit's progress in the
  shared ordination format (with `eig` as variances, on `pca()`'s scale),
  and optionally gemelli's distances between samples as a `dist` object
  (`distances = TRUE`). A C port of `gemelli`'s `rpca()`, it gives the same
  results to within rounding error, but works only with the non-zero counts,
  so sparse tables stay sparse and nothing of size samples by features is
  formed. It takes every count-table class that `pca()` takes, and its
  results are identical for any number of `cpus`.

* New `rpca_phylo()` runs a phylogenetic robust Aitchison PCA (Martino et
  al. 2022), gemelli's `phylogenetic_rpca()`: `rpca()` of a table expanded
  to every node of a tree, each node holding the counts of the tips below it
  weighted by its branch length, so that loadings describe clades as well
  as tips. It cuts the tree down to the tips with counts as gemelli does
  (with its `min_depth` filter), returns the tree it used with the nodes
  named as gemelli names them, and optionally the node counts
  (`node_counts = TRUE`). A C port of gemelli 0.0.13 (with the tree code of
  scikit-bio and iow it calls), it gives gemelli's results to within
  rounding error; the node counts are summed in parallel from each sample's
  non-zero tips, without a dense table. It joins `count_ord()`
  (`count_ord(counts, "rpca_phylo", tree = tree)`), and its results are
  identical for any number of `cpus`.

* New `blas_threads()` gets or sets how many threads R's BLAS/LAPACK library
  uses, while R is running. It works with OpenBLAS, Intel MKL, BLIS, and
  FlexiBLAS, returns `NA` where the count can't be read (e.g. R's reference
  BLAS or Apple's Accelerate), and never raises an error, warning, or
  message. Setting is limited to what is safe: the count is never raised
  above the library's own, and OpenMP's thread count is left as it is. No
  other ecodive function calls it.

* Sparse counts with samples in rows (a `dgCMatrix` with `margin = 1`), or as
  triplets (`dgTMatrix`, slam's `simple_triplet_matrix`), are now imported
  with counting sorts in linear time, instead of a quicksort: 4 s became
  0.06 s for a 5000 x 4000 `dgCMatrix` with 2.5 million non-zero values, and
  26 s became 0.3 s for 9.4 million. Every alpha and beta diversity function
  benefits, and results are unchanged. Out-of-range indices in these inputs
  are now an error.

* Sparse counts: explicitly stored zeros are now ignored, and repeated
  entries of a `dgTMatrix` summed, as the Matrix package defines them, so
  every function gives the same results as for the dense matrix (before,
  `norm = 'binary'` counted a stored zero as present). Malformed sparse
  objects (wrong index types or lengths, bad column pointers) are now an
  error rather than read out of bounds, and tables with more than ~268
  million non-zero values no longer overflow an internal size.

* Threads now come from a pool that lasts the whole session, and jobs are
  shared out to whichever threads are free. On a machine busy with other
  processes, multithreaded calls are no longer much slower than
  single-threaded ones (`pcoa()` of 2,000 samples: 0.8 s with 6 threads,
  where restarting threads for every step took 16 s), and small jobs no
  longer start threads at all. Results are unchanged, and still the same
  for any `cpus`.

* Every function's `cpus` is now capped at `n_cpus()`, the cores available,
  so that asking for more threads than that costs nothing. Results are the
  same for any `cpus`. Without the parallelly package, `n_cpus()` now
  returns at most 2 under `R CMD check --as-cran`, as parallelly's does.

* Fixed `NA` in an integer or logical count table being read as
  -2147483648: every alpha and beta diversity function now treats it as
  `NA`, as it does in a double table. Before, `bray()` returned -1 for every
  pair with that sample.

* Fixed the `pairs` argument of the beta diversity functions, which counted
  the columns of `counts` as samples even when samples are rows
  (`margin = 1`, the default). Out-of-range pairs could crash R, and a table
  with one feature was rejected as having fewer than two samples.

* Fixed `faith()`, the UniFrac family and every other function taking a
  `tree`, which failed when samples were columns (`margin = 2`) and the tree
  had tips that were not in the table.



# ecodive 2.3.2

* Fixed `aitchison()`, and every beta metric called with `norm = 'clr'`,
  returning `NaN` (or, for `chebyshev()`, wrong values) for every sample
  pair. Two bugs were involved:

    - Data with no zeros (pseudocount 0) hit `log(0) * 0`, which is `NaN`.
    - The pseudocount was truncated to an integer in C, so values below 1,
      including the automatic `min(counts[counts > 0]) / 2`, became 0, and
      values like 1.5 were silently rounded down.

* `gower()` with `norm = 'clr'` now computes each feature's range, and the
  distance for features absent from both samples, from that feature's CLR
  values instead of treating absent features as 0. `gower()` with negative
  inputs (e.g. `norm = 'rclr'`) now also uses the full min-to-max range.



# ecodive 2.3.1

* `read_tree()` now handles quoted labels per the Newick standard:

    - `''` inside a quoted label is decoded as one literal apostrophe, so
      `'it''s'` is read as `it's`.
    - Square brackets inside quotes are literal characters rather than
      comments, so labels like `'o__[Saprospirales]'` are no longer truncated.
    - Newlines and whitespace inside quotes are preserved.
    - Comments (`[...]`) outside of quotes may contain apostrophes.

* `read_tree()` raises an informative error, with the position of the
  problem, for malformed input such as an unclosed quote, comment, or
  parenthesis. Previously these could return a garbled tree.

* `read_tree()`'s parser is now linear-time and non-recursive, making it
  about 4x faster on a typical 100,000-tip tree and much faster on very deep
  trees.



# ecodive 2.3.0

* Added `ecodive` JOSS citation information.
* New `robust_aitchison()` beta diversity function.
* Switch from `testthat` to dependency-free `tinytest` for test suite.
* `read_tree` fixes for edge cases.



# ecodive 2.2.6

* When using CLR normalization:

    - Warns when zeros are detected and `pseudocount` is `NULL`.
    - Then continues with `pseudocount <- min(counts[counts>0]) / 2`.

## Breaking Changes

* `norm` now defaults to `"none"` instead of `"percent"`.
* `norm` parameter was removed from some beta diversity functions.



# ecodive 2.2.4

## Breaking Changes

* The `rarefy()` function has been simplified:

    - `drop = FALSE` now returns samples with original counts when sum < depth.
    - `drop = TRUE` no longer drops zero-sum features - only zero-sum samples.
    - `depth` can no longer be a fraction (only positive integers now).
    - `n_samples` parameter has been removed.
    - `warn` parameter has been added.



# ecodive 2.2.3

* Coerce counts to double before rarefaction.



# ecodive 2.2.2

* Drop explicit zeros from sparse matrices after rarefaction.



# ecodive 2.2.1

* Corrected CRAN `rchk` issue.



# ecodive 2.2.0

* Dedicated rarefaction backends for dense and sparse matrices.
* Auto-conversion to col-major sparse matrix for distance/diversity calculations.



# ecodive 2.1.0

## Breaking Changes

* `rescale` parameter is superseded by `norm` parameter.


## Fixes

* `alpha_div()` and `beta_div()` will only forward applicable parameters.
* No longer crashes when `pairs = integer(0)`.
* `pseudocount = NULL` correctly selects smallest non-zero count.



# ecodive 2.0.0

## Breaking Changes

* Input matrix is now samples as rows.
* Removed `weighted` parameter.

## Added New Diversity Methods

* Updated Alpha Options: `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`

* Updated Beta Options: `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`



# ecodive 1.0.0

* This is the first release of ecodive.
