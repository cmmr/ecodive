# Changelog

## ecodive 2.4.0

- ecodive now covers a whole diversity analysis: besides alpha and beta
  diversity, it ordinates samples and tests what explains their
  distances (see below). New vignettes,
  [`vignette("ordination")`](https://cmmr.github.io/ecodive/articles/ordination.md)
  and
  [`vignette("stat_tests")`](https://cmmr.github.io/ecodive/articles/stat_tests.md),
  walk through both, and the README and
  [`vignette("ecodive")`](https://cmmr.github.io/ecodive/articles/ecodive.md)
  give an overview of everything the package does. The package title is
  now “Parallel and Memory-Efficient Ecological Diversity Analysis”.

- New [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  draw ecodive’s results with base graphics (only ‘graphics’ and
  ‘grDevices’, from R 3.6.0 on), taking their titles, axis labels,
  subtitles and legends from what each result holds, and computing
  nothing new:

  - Ordinations
    ([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
    [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
    [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
    [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
    [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md),
    …) draw two axes, labelled with their share of the variance, with
    the method’s fit in the subtitle (stress, R-squared, …).
    `plot(ord, "Species", iris)` colours the samples by a column of
    metadata, with a legend, by the same `groups` and `data` as the
    tests; supervised ordinations colour by their own groups. Loadings,
    constraint arrows, centroids and a Procrustes rotation’s target are
    drawn on top.
  - `plot(envfit(...), add = TRUE)` adds the fitted variables to an
    ordination plot.
  - Permutation tests show the permutation distribution of their
    statistic, with the observed value marked;
    [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) and
    [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
    show each group’s within-group distance or dispersion, and
    [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
    its Procrustes fit (each with `type = "null"` for the distribution).
  - [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md)
    gives a dot chart of the taxa that contribute most to each pair’s
    dissimilarity,
    [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
    a volcano, effect, MA or forest plot (see below),
    [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md) a dot
    chart of W,
    [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md) a
    Venn diagram of the fractions (two to four sets),
    [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md)
    the adjusted R-squared at each step, and
    [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md) the
    best correlation for each number of variables.
  - Distance matrices give a heatmap, in the samples’ order or sorted by
    `groups`, with a colour key: those of
    [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
    and every metric function, of
    [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
    [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
    [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) (over
    features), and
    [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md)’s `dm`.
    These now have class `c("ecodive_dist", "dist")`, so that
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) knows them;
    base R has no
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
    `dist` objects, and every function for `dist` objects works on them
    as before. Code that checks `identical(class(d), "dist")` or
    `class(d) == "dist"` should use `inherits(d, "dist")` instead.
  - [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
    and
    [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
    give a volcano, effect, MA or forest plot of one covariate’s log
    fold changes,
    [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
    and
    [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
    a heatmap of the log fold change of each feature in each comparison,
    and
    [`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
    and
    [`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
    a dot chart of their statistics.

- Every function that returns a distance matrix,
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  and its 38 metric functions,
  [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md), takes
  `clust`, to cluster the samples before returning it: one of
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)’s linkage
  methods (`"average"`, `"complete"`, `"single"`, `"ward.D"`,
  `"ward.D2"`, `"mcquitty"`, `"median"`, `"centroid"`), or `"pcoa"` to
  order them along the first axis of
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).
  Clustering **reorders the samples** of the result, with their labels
  (unnamed samples are labelled with their input positions); the
  matrix’s `order` attribute gives each one’s position in the input.
  ecodive’s ordinations and tests put the samples back in their input
  order before analysing them, so metadata lines up whether matched by
  name or by position; for other packages, reorder metadata with
  `attr(d, "order")`. `clust = NULL`, the default, leaves the matrix as
  it was, and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  then clusters it by average linkage for the heatmap; `clust = FALSE`
  (kept as `"none"`) keeps the input order in the plot too.

- New partition of beta diversity (Baselga 2010, 2012, 2013, 2017), as
  ‘betapart’ computes it:

  - Eight new beta metrics, run by
    [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
    like any other:
    [`sorensen_turnover()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md)
    (Simpson dissimilarity, `"betasim"`),
    [`sorensen_nestedness()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md),
    [`jaccard_turnover()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md)
    and
    [`jaccard_nestedness()`](https://cmmr.github.io/ecodive/reference/sorensen_turnover.md)
    of presence/absence; and
    [`bray_balanced()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md),
    [`bray_gradient()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md),
    [`soergel_balanced()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md)
    and
    [`soergel_gradient()`](https://cmmr.github.io/ecodive/reference/bray_balanced.md)
    (the Ruzicka components) of abundances, ported from ‘abdiv’ (MIT).
    Each pair of components adds up to
    [`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md),
    [`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md),
    [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) or
    [`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md).
  - New
    [`beta_partition()`](https://cmmr.github.io/ecodive/reference/beta_partition.md)
    returns all three distance matrices of a family
    (`index = "sorensen"`, `"jaccard"`, `"bray"` or `"ruzicka"`) from
    one pass over the pairs, as `betapart::beta.pair()` and
    `beta.pair.abund()`;
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
    pair’s turnover against its nestedness.
  - New
    [`beta_partition_multi()`](https://cmmr.github.io/ecodive/reference/beta_partition_multi.md)
    gives the multiple-site values of a whole set of samples, as
    `beta.multi()` and `beta.multi.abund()`, and with `n_sites`, of
    `n_draws` random subsets of that many samples (as `beta.sample()`,
    with ecodive’s own random numbers: `seed`);
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
    draws’ densities.
  - `beta_div(x, "Bray-Curtis")` (alias `"braycurtis"`, or `"bc"`) still
    selects [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md)
    now that other metric names begin with “Bray-Curtis”.

- Distance matrices, ordinations, permutation tests and
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) have
  a title, which
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) uses unless
  given `main`
  ([`?ecodive_display`](https://cmmr.github.io/ecodive/reference/ecodive_display.md)):
  the method’s name, then the title of the result it was computed from,
  such as “Principal Coordinates Analysis of Bray-Curtis Dissimilarity”.
  It is kept in one attribute, `display`, a list for ecodive’s own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html), beside a distance
  matrix’s `clust`; set `attr(x, "display")$title` to retitle a result.
  A distance matrix’s `method` attribute is its metric’s ID (such as
  `"bray"`), as [`stats::dist()`](https://rdrr.io/r/stats/dist.html)
  sets it. Results don’t record the arguments or earlier steps they were
  made with. Arithmetic on a distance matrix (such as `sqrt(d)`) drops
  its title and `method`, since they no longer describe its values.
  [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) no
  longer return their `denominator` and `pseudocount` as attributes.

- Ordinations, distance matrices and ecodive’s data frames (from
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  and their post-hoc tests,
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
  [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md),
  [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md) and
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md))
  print a preview, in the manner of the tibble package, rather than
  every value
  ([`?print_ecodive`](https://cmmr.github.io/ecodive/reference/print_ecodive.md)):
  a header with the object’s kind, size and title; the first 10 rows
  (all of up to 20), with each column’s type; the columns that fit the
  console; and a footer counting what was left out and naming the other
  attributes. A distance matrix shows the top-left corner of its lower
  triangle. `print(x, n = Inf)` shows every row, and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`unclass()`](https://rdrr.io/r/base/class.html) give the plain
  object. Every function’s examples print its result.

- New [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) runs
  a principal coordinates analysis (classical multidimensional scaling)
  on a `dist` object, such as from
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md).
  It gives the same coordinates and eigenvalues as
  [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html), `ape::pcoa()`,
  and `vegan::wcmdscale()` (up to each axis’s sign), but computes only
  the first few axes (`n_axes = 3` by default), in parallel C, without
  forming an n-by-n matrix. Its results are identical for any number of
  `cpus`. `correction = "lingoes"` or `"cailliez"` adds the smallest
  constant to the squared distances, or to the distances, that leaves no
  negative eigenvalue, giving the same constants, eigenvalues and
  coordinates as `ape::pcoa()`, `vegan::wcmdscale(add = )` and
  `cmdscale(add = TRUE)`, again without an n-by-n (or 2n-by-2n) matrix:
  with 2,000 samples, `"cailliez"` takes 20 seconds where
  `cmdscale(add = TRUE)` takes 2 minutes. Every result has `correction`
  and `constant` attributes.

- New
  [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md)
  projects the features of a table (taxa counts, or environmental
  variables) onto the axes of a
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) result,
  for drawing arrows on a PCoA plot. It computes the same scores as
  scikit-bio’s
  [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md),
  in parallel C, keeps sparse tables sparse, and takes every table type
  and `norm` that
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) does.
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) gives
  the same directions with permutation tests.

- Every ordination function takes the number of axes to return as
  `n_axes` (3 by default for the methods whose axes have eigenvalues:
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  and [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md);
  2 for [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md),
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md),
  [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md)). The
  name says what it counts; `k`, as
  [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) and vegan call
  it, often means a number of neighbours in these methods.

- New
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
  runs any distance-based ordination method by name, as
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  does for beta diversity metrics and
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
  does for statistical tests, with the method’s own arguments:
  `dist_ord(dm, "nmds")` is `nmds(dm)`.
  [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  and
  [`match_dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  list and look up the methods.
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md),
  with
  [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  and
  [`match_count_ord()`](https://cmmr.github.io/ecodive/reference/count_ords.md),
  does the same for the methods that take a count table rather than
  distances, [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md),
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md),
  [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) and
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md):
  `count_ord(counts, "pca")` is `pca(counts)`. Each wrapper points to
  the other when given the other’s method, and to
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
  when given a test. Every ordination method returns the same kind of
  result: a sample-by-axis matrix (of class `"ecodive_ord"`, which
  prints its attributes’ names rather than their contents) whose
  attributes start with `method`, `eig`, and `prop`. Every method
  rejects negative distances (beyond rounding error), and `seed` can be
  given to any of them through
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md).

- New [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) runs
  non-metric multidimensional scaling with the SMACOF algorithm of
  scikit-learn’s non-metric `MDS()`, ported to parallel C. It starts
  from the PCoA solution and from `n_starts` random configurations,
  keeps the one with the lowest Kruskal stress-1, and centers it,
  rotates it to principal axes, and scales it to the distances. Tied
  distances use Kruskal’s primary approach (as in `vegan::metaMDS()`),
  or scikit-learn’s secondary one with `ties = "secondary"`. The
  distances are sorted once and shared by all starts; nothing n-by-n is
  formed. Each start runs at most `max_iters` iterations, or stops when
  the stress improves by less than `tol`. Random starts come from a
  `seed` argument, and results are identical for any number of `cpus`.
  The result’s attributes are `stress`, `converged`, `n_iters`,
  `n_starts`, `best_start`, `n_best`, and `ties`.

- New [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) runs
  metric multidimensional scaling with the SMACOF algorithm of
  scikit-learn’s metric `MDS()`, fitting the distances themselves rather
  than their ranks. It shares
  [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md)’s engine,
  starts, arguments (without `ties`) and result format, and has no limit
  on the number of samples. `dist_ord(dm, "mds")` and `"metric MDS"`
  select it; `"classical MDS"`, `"pco"` and `"cmdscale"` select
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md).

- New [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md)
  runs an Isomap (isometric feature mapping) ordination: it links each
  sample to its `n_neighbors` nearest neighbors (or to those within a
  `radius`), replaces each distance with the shortest path through that
  graph, and runs
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) on those
  geodesic distances, which keeps long ecological gradients straight
  where [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)
  bends them into an arch. A C port of scikit-learn’s
  `Isomap(metric = "precomputed")`, it gives the same geodesic distances
  and embedding to rounding, breaks ties by sample order, and joins
  disconnected neighbor graphs at their closest samples, with a warning
  (`vegan::isomap()` stops instead; note that its `k` is the number of
  neighbors, while
  [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md)’s
  `n_axes` is the number of axes). The neighbor and shortest-path
  searches run in parallel without an n-by-n matrix, and results are
  identical for any number of `cpus`.

- New [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) runs
  t-distributed stochastic neighbor embedding (t-SNE) on a `dist`
  object, as a port of scikit-learn’s `TSNE` to parallel C.
  `gradient = "auto"` uses the exact gradient for up to 500 samples (or
  more than 3 dimensions) and the Barnes-Hut approximation (a quadtree
  in 2 dimensions, an octree in 3) above that. It finds nearest
  neighbors exactly from the distances and never forms an n-by-n matrix.
  It starts from the PCoA solution by default, sets the learning rate
  with scikit-learn’s rule (`learning_rate = "auto"`), stops early once
  the layout stops improving (at most `max_iters` iterations; the number
  run is the `n_iters` attribute), and lowers `perplexity` for small
  datasets. Its layouts match scikit-learn’s and `Rtsne::Rtsne()`’s in
  quality. Random starts come from a `seed` argument, and results are
  identical for any number of `cpus`.

- New [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md) runs
  a uniform manifold approximation and projection (UMAP) on a `dist`
  object, following `uwot::umap()` and its defaults: exact nearest
  neighbours, smooth-kNN calibration and fuzzy union (identical to
  uwot’s), a spectral start from the normalized graph Laplacian, and the
  layout optimization of uwot’s batch mode (`batch = TRUE`: per-epoch
  gradients, then Adam), which gives uwot’s coordinates from the same
  start and seeds. Every stage runs in parallel C, the neighbours
  straight from the dist triangle without an n-by-n matrix. Random
  numbers come from a `seed` argument, and results are identical for any
  number of `cpus`. `min_dist` defaults to the canonical 0.1
  (umap-learn’s), not uwot’s 0.01.

- New
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  runs a permutational multivariate analysis of variance on a `dist`
  object and sample metadata. Like ecodive’s other tests, it takes
  `groups`: one value per sample (categories or a numeric covariate), or
  a column name of `data`. For several variables, give `data` and a
  `formula` instead: factors, covariates, interactions, and
  `Condition()` terms. It gives the same table as `vegan::adonis2()`
  (Df, sums of squares, R2, F, and p-values), testing the model as a
  whole or each term sequentially, marginally, or by contrast, with
  optional `strata`. Given the same permutation matrix, its permuted F
  statistics and p-values are identical to `adonis2()`’s. It returns a
  list of class `"ecodive_test"`, like
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md): the
  model’s F, p-value, z-score (as `summary(vegan::permustats())`) and
  number of permutations come first, whatever `by` is, then its R2 and
  `adonis2()`’s table, as `table`. The sums of squares come straight
  from the dist triangle, without the n-by-n matrices `adonis2()` builds
  for each permutation, and permutations run in parallel. On one thread
  it is 25 to 50 times faster than `adonis2()` with 200 samples, and 100
  to 300 times faster with 1,000. Permutations come from a `seed`
  argument, and results are identical for any number of `cpus`. Like
  `adonis2()`, it also takes the distances on the left-hand side of a
  formula: `permanova(dm ~ treatment + ph, data)`.

- New [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) runs
  a multi-response permutation procedure on a `dist` object and each
  sample’s group, as `vegan::mrpp()`: the observed and expected weighted
  mean within-group distance (delta), the chance-corrected within-group
  agreement A, each group’s mean distance, and a permutation p-value,
  with group weights `'n'`, `'n-1'` or `'n(n-1)'` (vegan’s `weight.type`
  1 to 3) and optional `strata`. It also gives the standardized effect
  size z of `summary(vegan::permustats())`. Given the same permutation
  matrix, its permuted deltas, p-values and z are those of
  `vegan::mrpp()`. Each permutation costs one comparison per pair of
  samples, read straight from the dist triangle; with 500 samples and
  999 permutations it is about 60 times faster than `vegan::mrpp()` on
  one thread. Results are identical for any number of `cpus`. Like
  ecodive’s other permutation tests, it returns a list of class
  `"ecodive_test"` whose `statistic`, `p_value`, `z` and `n_perm` come
  first, and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) turns
  it into one row with the same columns for every test.

- New [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md)
  runs an analysis of similarities (ANOSIM) on a `dist` object and each
  sample’s group, given as a vector or as a column of sample metadata,
  with optional `strata`. Given the same permutation matrix, its R
  statistic, permuted R values and p-value match `vegan::anosim()`’s. It
  ranks the distances once, with a radix sort, and then each permutation
  needs only one pass over the ranks, without the n-by-n matrices
  `vegan::anosim()` builds; permutations run in parallel. With 300
  samples and 999 permutations it takes under a tenth of a second,
  against about seven for `vegan::anosim()`. Results are identical for
  any number of `cpus`.

- New
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
  tests the homogeneity of multivariate dispersions (PERMDISP), as
  `vegan::betadisper()` followed by `vegan::permutest()`: each sample’s
  distance to its group’s centroid (the default) or spatial median, each
  group’s mean distance, and the ANOVA F with a permutation p-value,
  plus optional pairwise t tests and small-group bias adjustment.
  Distances to centroids come from the distances within groups alone,
  with no eigendecomposition: they match `betadisper()`’s to rounding
  error, and are about 100 times faster with 1,000 samples.
  `betadisper()` defaults to spatial medians instead: pass
  `type = "median"` to compare with it.
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
  solves them exactly, where `betadisper()`’s optimizer stops at about
  1e-5. Given the same permutation matrix, its p-values are identical to
  `permutest()`’s. Permutations come from a `seed` argument, can be
  restricted to `strata`, and results are identical for any number of
  `cpus`.

- New [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)
  runs a Mantel test of the correlation between two distance matrices
  (Pearson, Spearman, or Kendall), and new
  [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md)
  a partial Mantel test, which removes the effect of a third (`given`).
  They give the same statistics as `vegan::mantel()` and
  `vegan::mantel.partial()`, and, given the same permutation matrix, the
  same permuted statistics and p-values. The second and third matrices
  may also be numeric vectors, one value per sample (such as age), and
  samples are matched by name. They add lower-tail and two-sided tests
  (`alternative`). Each permutation is one pass over the dist triangle,
  with no n-by-n matrix (Kendall’s tau-b by Knight’s merge sort rather
  than [`cor()`](https://rdrr.io/r/stats/cor.html)’s quadratic count),
  and permutations run in parallel. Permutations come from a `seed`
  argument, and results are identical for any number of `cpus`.

- New
  [`mantel_correlogram()`](https://cmmr.github.io/ecodive/reference/mantel_correlogram.md)
  draws a Mantel correlogram: how community similarity changes with
  (geographic or other) distance. It splits the pairs of samples into
  distance classes (Sturges’ rule by default, or `n_classes` or
  `breaks`) and tests each class with a Mantel test against the class’s
  model matrix, with p-values adjusted progressively (`p_adjust`,
  `progressive`). It gives the same classes, statistics, and p-values
  (raw and adjusted, given the same permutation matrix) as
  `vegan::mantel.correlog()`. `y` may be a distance matrix, coordinates,
  or columns of `data`. Every class’s statistic depends on a permutation
  only through a sum over the class, so each permutation is one pass
  over the pairs for all classes (for Kendall’s tau-b too), in parallel.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  correlogram.

- New
  [`mantel_regression()`](https://cmmr.github.io/ecodive/reference/mantel_regression.md)
  fits a multiple regression on distance matrices (MRM): one distance
  matrix (such as community dissimilarities) regressed on several others
  (geographic distance, differences in pH, …), with permutation p-values
  for each coefficient (pseudo-t) and for the whole model (F). Terms
  come as distance matrices, coordinates, one value per sample, or a
  formula of sample metadata, whose numbers become absolute differences,
  matrices Euclidean distances, and factors 0 (same) or 1 (different).
  `method = 'spearman'` regresses the ranks. It gives the same
  coefficients, R^2, F, and p-values (given the same permutations) as
  `ecodist::MRM()`. Each permutation is one pass over the pairs, in
  parallel, with `strata` for restricted permutations. `backward = TRUE`
  simplifies the model by backward elimination (Lichstein 2007): the
  least significant term is removed and the model refitted, with the
  same permutations, until every term left has a p-value of at most
  `alpha`; `elimination` records each step.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  standardized coefficients.

- New `dist_test(dm, test, groups, data, ...)` runs any of these tests
  by name, as
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md)
  does for ordination methods:
  `dist_test(dm, "mrpp", "treatment", meta)` is
  `mrpp(dm, "treatment", meta)`, and other arguments (such as
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md)’s
  `method`) pass through. It checks that the variable is of a kind the
  test takes: categorical groups for `anosim`, `mrpp` and `permdisp`,
  numeric values for `mantel` and `mantel_partial` (passed as their
  `y`), and either for `permanova`.
  [`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  and
  [`match_dist_test()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  list and look up the tests, with forgiving name matching (`"adonis2"`
  finds `permanova`, `"betadisper"` finds `permdisp`).

- [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
  [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) and
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)
  take `pairwise = TRUE` to also compare each pair of groups, as is
  usual after an overall test across several groups, and `p_adjust` to
  adjust the pairs’ p-values for multiple testing: any of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html),
  Holm’s method by default.
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
  [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md) and
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) rerun the
  test on each pair’s samples, as q2-diversity’s
  `beta_group_significance(pairwise=True)` and pairwiseAdonis’s
  `pairwise.adonis2()` do, with the pair’s statistics, z-score, sample
  size, and raw and adjusted p-values in a `pairwise` data frame. Given
  the same permutations, the pairs’ statistics and p-values are
  identical to vegan’s tests on each pair. Each pair’s random
  permutations come from `seed` and the pair, so results are identical
  for any number of `cpus`; a supplied `permutations` matrix is cut down
  to each pair’s samples. A pair that can’t be tested, such as two
  groups of one sample, gives `NA`s with a warning.
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md)’s
  pairwise t tests gain `p_adjusted` and `p_adjusted_perm` (each of its
  p-value columns adjusted across the pairs, which `vegan::permutest()`
  doesn’t do), and `z`, `n` and `n_perm` columns, so its pairwise table
  now has more columns than before.

- One rule aligns every per-sample argument (`groups`, `strata`, and the
  Mantel tests’ `y` and `given`) in
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md),
  [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md),
  [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md),
  [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md),
  [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md),
  [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md),
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md):
  with `data`, a column name or one value per row of `data`; without it,
  one value per sample, matched by name when named, else by position.
  The same checks and messages apply everywhere. As with
  `vegan::adonis2()`,
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  needs `groups` or a `formula` (use `formula = ~ .` to test every
  column of `data`). `permutations` is limited to 1e8.

- New [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)
  runs a distance-based redundancy analysis (db-RDA): the axes of a
  `dist` object that sample metadata explain, given as a formula such as
  `~ treatment + Condition(site)`. It gives the same eigenvalues,
  inertia, site scores, biplot scores and centroids as `vegan::dbrda()`
  (up to each axis’s sign; not `vegan::capscale()`, which differs for
  non-Euclidean distances), without forming an n-by-n matrix: the
  constrained axes come from a small eigenproblem, and unconstrained
  axes (`MDS1`, …) are only computed when the constraints give fewer
  than `n_axes` axes. Its results are identical for any number of
  `cpus`. As in `vegan`, `dbrda(dm ~ treatment, data)` works too, and a
  formula of only `Condition()` terms gives the unconstrained ordination
  after removing their effect.

- [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) takes
  the same `groups`, `data` and `formula` as
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md):
  one variable as `groups` (a vector or a column of `data`), or `data`
  and a `formula`, or `dm ~ terms`, so
  `dist_ord(dm, "dbrda", "treatment", meta)` works. Its results gain
  `r_squared` and `adj_r_squared`, as `vegan::RsquareAdj()`, including
  for partial models.

- New
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)
  runs a constrained analysis of principal coordinates. It gives the
  same results as `vegan::capscale()` (and
  `phyloseq::ordinate(method = "CAP")`), with any imaginary inertia from
  negative eigenvalues reported in `inertia`. Euclidean distances take
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)’s fast
  path; otherwise a threaded eigendecomposition makes it about five
  times faster than `vegan::capscale()` on one thread. Also available as
  `dist_ord(dm, "capscale", ...)`.

- New [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) runs a
  redundancy analysis of a count table: the axes of the (normalized)
  table that sample metadata explain (`RDA1`, …), then principal
  components of the rest (`PC1`, …). `norm = 'hellinger'` or `'chord'`
  gives transformation-based RDA (tb-RDA), and `Condition()` terms
  partial RDA. It takes `groups`, `data` and `formula` as
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) does
  (and `counts ~ terms`, as `vegan::rda()`), and with no model it is a
  PCA. It gives the same site scores, species scores (`loadings`),
  eigenvalues, inertia, R-squared, biplot scores and centroids as
  `vegan::rda()` (up to each axis’s sign), sharing
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)’s
  solver and
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)’s sparse
  table products, so the table is never made dense: with 3,000 samples
  and 2,000 taxa it takes half a second where `vegan::rda()` takes a
  minute. Results are identical for any number of `cpus`. Also available
  as `count_ord(counts, "rda", ...)`. Its permutation test is
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  of Euclidean distances of the normalized table.

- New [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) runs a
  canonical correspondence analysis (ter Braak 1986) of a table of
  abundances: the correspondence analysis axes, in the chi-square
  metric, that sample metadata explain (`CCA1`, …), then correspondence
  analysis axes of the rest (`CA1`, …). It takes `groups`, `data` and
  `formula` as
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) does,
  including `Condition()` for partial CCA and `counts ~ terms`, and with
  no model it is a correspondence analysis. It gives the same site
  scores, species scores (`loadings`), eigenvalues, inertia, R-squared,
  biplot scores and centroids as `vegan::cca()` (up to each axis’s
  sign), and the same permutation-based adjusted R-squared as
  `vegan::RsquareAdj()` for the same permutations (`permutations` and
  `seed`), partial CCA included. It shares
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md)’s solver,
  with the chi-square weights, and multiplies the table’s non-zero
  values only, so a sparse table is never made dense: with 3,000 samples
  and 2,000 taxa it takes half a second where `vegan::cca()` takes a
  minute. Results are identical for any number of `cpus`. Also available
  as `count_ord(counts, "cca", ...)`.

- New [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) runs a
  correspondence analysis of a table of abundances: the unconstrained
  ordination in the chi-square metric, with species scores (`loadings`).
  It gives the same site scores, species scores, eigenvalues and inertia
  as `vegan::ca()` (`vegan::cca()` without constraints),
  `ade4::dudi.coa()` and scikit-bio’s
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) (up to each
  axis’s sign), and is the same as
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) with no
  model, whose sparse solver it shares. Results are identical for any
  number of `cpus`. Also available as `count_ord(counts, "ca")`; `"CA"`
  selects it, and `dist_ord(dm, "ca")` points to
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md).

- New [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md)
  fits sample metadata (`data`) and feature tables (`counts`, such as
  taxa abundances) onto any ordination, for biplots: an arrow and
  r-squared for each numeric variable or feature, the centroid of each
  factor level, and permutation p-values. It gives the same arrows,
  r-squared and centroids as `vegan::envfit()`, and the same p-values
  for the same permutation matrix. Each permutation costs a few
  multiply-adds per non-zero value of the table, with no refitting, and
  sparse tables are never made dense, so thousands of taxa with 999
  permutations take well under a second (minutes in `vegan`). It fits
  every axis of the ordination by default (`vegan` fits the first two),
  and its results are identical for any number of `cpus`.

- New
  [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md)
  rotates, reflects, and optionally dilates one ordination to fit
  another of the same samples, and reports the Procrustes sum of
  squares, m^2: how similar, say, Bray-Curtis and UniFrac PCoAs (or 16S
  and metagenomic ones) of the same samples are. It gives the same
  rotation, translation, dilation, residuals and m^2 as
  `vegan::procrustes()` (and the same symmetric m^2 as SciPy’s
  `scipy.spatial.procrustes()`), but takes its arguments as
  `procrustes(ord, target)`: `vegan::procrustes(target, ord)`. vegan’s
  `scale` is called `dilate`. Inputs may be `dist` objects, ordinated
  with [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)
  first (`n_axes` axes), and samples are matched by name. The result is
  the rotated ordination, with the rest as attributes. New
  [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
  runs the Procrustes randomization test (PROTEST), giving the same
  statistic as `vegan::protest()` and the same p-values for the same
  permutation matrix, in parallel C, with `strata`, and identical
  results for any number of `cpus`.

- New [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md)
  runs a partial least squares discriminant analysis (PLS-DA): the
  components of a count table that best separate known groups, from PLS
  regression on a dummy matrix of the groups (ported from scikit-learn’s
  `PLSRegression`). It gives the same scores, weights, loadings and Y
  loadings as scikit-learn, and the same scores, explained variance and
  VIP as `mixOmics::plsda()` (up to each component’s sign), and also
  returns each feature’s VIP and each group’s centroid. Each component’s
  weights are computed directly rather than by NIPALS iteration, so they
  don’t depend on a convergence tolerance. It takes `groups` and `data`
  as
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  does, every `norm` (`'clr'` is the usual choice for microbiome
  counts), and sparse tables, which it never makes dense: 1,000 samples
  and 5,000 taxa take under a tenth of a second, where
  `mixOmics::plsda()` takes over a second. Results are identical for any
  number of `cpus`. `count_ord(counts, 'plsda', groups)` runs it too.
  Like every supervised method, it separates any groups it is given, so
  validate the separation (e.g. with
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)).

- New [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md)
  runs an orthogonal PLS discriminant analysis (OPLS-DA) of two groups:
  one predictive component (`p1`) that separates them and orthogonal
  components (`o1`, `o2`, …) for the rest of the table’s variation
  (ported from pyopls). It gives the same orthogonal scores, weights and
  loadings as pyopls, and the same scores, loadings, weights, R2X, R2Y
  and VIP as `ropls::opls()` (up to each component’s sign), and also
  returns each group’s centroid. It takes `groups`, `data`, every `norm`
  and every table type as
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) does,
  and never makes a sparse table dense: 1,000 samples, 5,000 taxa and
  two orthogonal components take about a tenth of a second, where
  `ropls::opls()` takes 8 (including its cross-validation). Results are
  identical for any number of `cpus`.
  `count_ord(counts, 'oplsda', groups)` runs it too. For more than two
  groups, use
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md).

- New [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) runs a
  linear discriminant analysis (Fisher’s discriminant axes, ported from
  scikit-learn’s `LinearDiscriminantAnalysis` with its `"svd"` solver).
  It gives the same scores, coefficients (`scaling`), group means and
  proportion of trace as
  [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) and its
  [`predict()`](https://rdrr.io/r/stats/predict.html) (up to each axis’s
  sign), with optional `prior` probabilities, and the same scores as
  scikit-learn up to a constant factor (scikit-learn divides the
  within-group covariance by n, not n - g). It takes `groups` and `data`
  as [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md)
  does, every `norm` and every table type, and stops with an error when
  the within-group covariance is singular, as it is for microbiome
  tables with more taxa than samples: run it on a few PCoA or PCA axes
  instead. 5,000 samples and 200 features take 0.2 seconds, where
  [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) takes 0.8.
  Results are identical for any number of `cpus`.
  `count_ord(counts, 'lda', groups)` runs it too.

- Ordinations of many axes are faster on several threads. The block
  Lanczos solver behind
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md),
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) and
  [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md)’s
  spectral start now splits its own steps on the basis (products with
  small matrices, dot products and, for more than four axes, block
  Gram-Schmidt) among `cpus`, each as far as its work is worth it, as
  well as the products with the data. On 4 threads,
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) of 50
  axes takes 1.4 s instead of 2.5 for 1,000 samples and 15 s instead of
  27 for 5,000, and
  [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) of 40
  axes for 2,000 samples 3.3 s instead of 5.6. Results stay identical
  for any number of `cpus`; for up to four axes they are unchanged, and
  for more they differ from before only by rounding.

- New [`capda()`](https://cmmr.github.io/ecodive/reference/capda.md)
  runs Anderson & Willis’s (2003) discriminant canonical analysis of
  principal coordinates (CAP, as in PRIMER-PERMANOVA+):
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) of the
  first m axes of
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md), for any
  distances, with m chosen by leave-one-out classification success (up
  to `max_pcoa_axes`, 50 by default; on ties, the smallest m), and a
  permutation test of the trace statistic (the sum of the squared
  canonical correlations), with `strata` or a supplied permutation
  matrix. It gives the same m, classification success, leave-one-out
  classes, squared canonical correlations and scores as
  `BiodiversityR::CAPdiscrim()`, whose leave-one-out posteriors
  underflow for samples far from every group, and which stops at a group
  of one sample
  ([`capda()`](https://cmmr.github.io/ecodive/reference/capda.md) counts
  such samples as misclassified). The success for every m considered and
  the leave-one-out table are returned. Results are identical for any
  number of `cpus`. `dist_ord(dm, 'capda', groups)` runs it too; `'CAP'`
  still selects
  [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md),
  which is db-RDA.
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md)’s error for
  a singular within-group covariance now suggests
  [`capda()`](https://cmmr.github.io/ecodive/reference/capda.md).

- New [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md)
  runs a similarity percentages analysis (SIMPER): it breaks the mean
  Bray-Curtis dissimilarity between each pair of groups down into each
  taxon’s contribution, with its sd, ratio, the groups’ mean abundances
  and the cumulative proportion, following ecopy’s
  [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md), and
  gives each taxon a permutation p-value (Clarke 1993), adjusted within
  each pair by `p_adjust` (Holm’s method by default). It gives the same
  results as `vegan::simper()` and
  [`summary()`](https://rdrr.io/r/base/summary.html) of it, as one data
  frame. Its p-values count permuted contributions at least as large as
  the observed one, as Clarke (1993) does; vegan counts only larger
  ones, which makes its p-values too small for small groups. As in the
  pairwise tests, each pair’s permutations shuffle only its own samples
  (vegan shuffles all of them). It takes `groups` and `data` as
  [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  does, `strata`, `norm` (`'none'` by default, as vegan) and sparse
  tables, which it never makes dense, and runs in parallel C: 300
  samples and 3,000 taxa, three pairs with 99 permutations, take under 2
  seconds on one thread, where `vegan::simper()` takes 27 seconds for
  one pair. Results are identical for any number of `cpus`. Its docs
  explain why SIMPER’s ranking mixes up differences in mean and in
  variability (Warton et al. 2012).

- New [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md)
  finds the subset of environmental variables whose distances among the
  samples best match the community distances (BIOENV; Clarke & Ainsworth
  1993): for each subset size, the subset with the highest rank
  correlation with `dm`, as one data frame, with the best subset overall
  as attribute `best`. It follows scikit-bio’s
  [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md) and
  gives the same subsets and correlations as it and as
  `vegan::bioenv()`, with vegan’s options: `method` (Spearman, Pearson
  or Kendall), `metric` (`'euclidean'`, `'mahalanobis'`, `'manhattan'`
  or `'gower'`), `max_vars` (vegan’s `upto`) and `given` (vegan’s
  `partial`, also with Kendall, which vegan refuses). Variables come
  from `data`, by `formula` or as its numeric columns, and samples are
  matched by name. Every subset is searched in parallel C, with no
  n-by-n matrix: 100 samples and 12 variables (4,095 subsets) take 1.2
  seconds on one thread, where `vegan::bioenv()` takes 5.5. Results are
  identical for any number of `cpus`.

- New [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md)
  partitions the variation in a community among two to four sets of
  explanatory variables (Borcard et al. 1992): the fractions each set
  explains on its own, those they share, and the residuals, from the
  adjusted R-squared of (partial) RDAs of a count table (with `norm`,
  such as `"hellinger"`), db-RDAs of a `dist`, or, with
  `chisquare = TRUE`, CCAs with permutation-estimated adjusted R-squared
  (Peres-Neto et al. 2006). It gives the same fractions as
  `vegan::varpart()` with two, three and four sets, as one data frame
  named with vegan’s letters (`[a]`, `[b]`, …), with vegan’s union and
  conditional fractions as attributes and a text Venn diagram when
  printed. Each set is a formula over `data`, a data frame, a matrix, or
  a vector, matched to the samples by name. For RDA and db-RDA, one pass
  over the data gives every union’s R-squared, with no n-by-n matrix and
  sparse tables kept sparse: with 3,000 samples, 2,000 taxa (5%
  non-zero) and three sets, it takes 0.13 seconds on 4 threads
  (Hellinger table or Bray-Curtis distances), where `vegan::varpart()`
  takes 3.5. Results are identical for any number of `cpus`.

- New
  [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md)
  picks explanatory variables one at a time by adjusted R-squared, with
  a permutation test of each added variable and Blanchet et al.’s (2008)
  double stopping rule, for RDA of a count table (with `norm`) or db-RDA
  of a `dist`. It selects the same variables in the same order as
  `vegan::ordiR2step()`, with the same adjusted R-squared and F, and the
  same p-values given the same permutations. Candidates are the terms of
  a formula over `data` (factors enter whole, interactions after their
  main effects, `Condition()` terms are covariables). By default it
  first tests the global model and selects nothing if that is not
  significant; `max_vars` caps the number of variables, and `strata`
  restricts permutations. Every candidate’s adjusted R-squared comes
  from one pass over the data, and only the chosen candidate of each
  step is tested: with 1,000 samples, 500 taxa, 11 candidates and 199
  permutations, it takes 0.17 seconds on 4 threads for the Hellinger
  table and 0.08 for Bray-Curtis distances, where `vegan::ordiR2step()`
  takes 58 and 85. Results are identical for any number of `cpus`.

- New
  [`principal_response_curves()`](https://cmmr.github.io/ecodive/reference/principal_response_curves.md)
  (PRC; van den Brink & ter Braak 1999) summarizes how treated
  communities diverge from a control over time: a partial RDA,
  `rda(counts ~ treatment:time + Condition(time))`, whose canonical
  coefficients give one curve per treatment (the control’s is 0) and
  whose species weights say how closely each taxon follows them. Its
  coefficients, weights and eigenvalues are those of `vegan::prc()` and
  [`summary()`](https://rdrr.io/r/base/summary.html) in each of vegan’s
  scalings (`scaling`, `"symmetric"` by default, as vegan’s
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)), up to each
  axis’s sign: here each curve’s largest coefficient is positive. It
  tests the first axis and all axes by permutation, with the same F and
  p-values as `anova(prc, first = TRUE)` and `anova(prc)` given the same
  permutations; `plots` permutes whole experimental units (each one’s
  series of samples kept together), and `strata` permutes within blocks.
  Both statistics of every permutation come from one product of the
  sparse table with the permuted design, in threaded C, the same for any
  `cpus`: for 480 samples and 2,000 taxa, both tests with 999 whole-plot
  permutations take 3 seconds on 4 threads, where `vegan::prc()` and its
  two [`anova()`](https://rdrr.io/r/stats/anova.html) calls take 12
  minutes. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  draws the PRC diagram, with the taxa’s weights on a scale at its
  right.

- New [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) and
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md), the
  additive and isometric log-ratio transforms, complete the clr and rclr
  that `norm` already offers.
  [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) divides by
  a `denominator` feature that you choose (there is no default:
  scikit-bio uses the first feature and compositions the last).
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) takes the
  Helmert basis that scikit-bio uses by default, a sequential binary
  partition (converted as scikit-bio’s `sbp_basis()`), or any
  orthonormal basis, matched to the features by name. Both give the same
  coordinates as scikit-bio’s
  [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) and
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md)
  ([`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) to the
  last bit, [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md)
  to 1e-15), take every table type with zeros handled by the shared
  `pseudocount` rule, and return a dense matrix in the input’s
  orientation with `method`, `denominator` or `basis`, and `pseudocount`
  attributes, ready for
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) or
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md). The
  default basis is applied by a running sum, without a dense product:
  2,000 samples by 5,000 taxa take 0.2 seconds, where
  `normalize(x, 'clr') %*% t(basis)` takes 30 with R’s reference BLAS.
  Euclidean distances on
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) coordinates
  are
  [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
  distances, for any basis. Results are identical for any number of
  `cpus`.

- New
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  replaces the zeros of a count table with small positive values, so
  that log ratios (`norm = 'clr'`,
  [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
  [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md)) can be
  taken without a pseudocount: its result has no zeros, so they use none
  and give no warning. It runs any zero replacement method by name, as
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
  does for ordinations, and
  [`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)
  and
  [`match_zero_replacement()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)
  list and look up the methods. The first is
  [`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md),
  multiplicative replacement (Martín-Fernández et al. 2003): every zero
  becomes `delta` (one value, or one per sample; scikit-bio’s `1 / D^2`
  by default) and each sample’s non-zero proportions shrink by the same
  factor, keeping their ratios. Its results are bit-identical to
  scikit-bio’s `multi_replace()`. Every method returns a dense matrix of
  proportions in the input’s orientation, with `method`, its parameters
  and `n_replaced` attributes, and takes every table type: a sparse
  2,000 by 5,000 table takes 0.2 seconds.

- New
  [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
  and [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md)
  give one log ratio per sample between a `numerator` and a
  `denominator` group of features (names, indices or logical vectors),
  as a named vector that drops straight into metadata or a test.
  [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
  is the log of the ratio of summed abundances, as Qurro plots it
  (identical to Qurro’s `qarcoal()`); it needs no pseudocount by
  default, and gives `NA`, with one warning, where a sum is

  0.  [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md)
      is the isometric log-ratio balance, the scaled log ratio of
      geometric means: bit-identical to the matching coordinate of
      [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) with a
      sequential binary partition, and equal to scikit-bio’s to 6e-16,
      without building the full basis; zeros follow the shared
      `pseudocount` rule over the whole table, and the output of
      [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
      needs none. Both take every table type and keep sparse tables
      sparse: 2,000 samples by 5,000 taxa take under 0.1 seconds.

- Qurro-style rank plots and log ratios (Fedarko et al. 2020; Morton et
  al. 2019’s reference frames):

  - `plot(fit, type = "rank")` draws every feature as a bar of its
    effect, sorted, for the results of
    [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
    [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
    [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
    [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
    [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
    [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
    and [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md),
    and `plot(ord, type = "rank")` the loadings on one axis of
    [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
    [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
    [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md),
    [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) and other
    ordinations with loadings. `numerator` and `denominator` highlight a
    log ratio’s features. `plot(ord, type = "n")` and other `type`s
    still reach
    [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html).
  - New `log_ratio_features(x, n)` picks the `n` top- and bottom-ranked
    features of such a result (Qurro’s “autoselection”), as `numerator`
    and `denominator` for
    [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
    or
    [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md).
    It selects the same features as Qurro 0.9.0’s own code in every
    reference case, ties included, and their
    [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
    matches Qurro’s `qarcoal()` to within one rounding of
    [`log()`](https://rdrr.io/r/base/Log.html).
  - [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
    and
    [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md)
    results have class `ecodive_log_ratio`, with a short
    [`print()`](https://rdrr.io/r/base/print.html) and a
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
    samples’ values by group (box plot with the samples over it, spread
    without the random number generator), against a numeric variable, or
    all together, leaving out and counting undefined values. They stay
    numeric vectors everywhere else; subsetting,
    [`sort()`](https://rdrr.io/r/base/sort.html) and arithmetic return
    plain named vectors.

- New [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) measure
  proportionality between every pair of features, the compositionally
  valid alternative to correlating relative abundances (Lovell et
  al. 2015). [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md)
  is the variance of the log ratio, equal to scikit-bio’s
  `pairwise_vlr()` to 4e-15 but computed from differences rather than a
  covariance matrix, so nearly proportional pairs keep their precision.
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) is Erb
  & Notredame’s proportionality correlation (propr’s `"rho"`) and
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) the
  symmetric phi of Quinn et al. (propr’s `"phs"`), against the clr by
  default or a `denominator` (one feature, an alr, or a set). Each
  returns a `dist` over the features, ready for
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) or
  [`hclust()`](https://rdrr.io/r/stats/hclust.html)
  ([`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), a
  similarity, for thresholding). They take every table type and keep it
  sparse, run in parallel, and give identical results for any `cpus`:
  all 4.5 million pairs of 3,000 taxa in 500 samples take about a second
  on 4 cores.

- New [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)
  estimates the correlations between taxa’s log absolute abundances from
  relative abundances with SparCC (Friedman & Alm 2012), the
  co-occurrence networks most often reported for microbiome data:
  Dirichlet draws of each sample’s fractions (so zeros need no
  pseudocount), basis variances under SparCC’s sparsity assumption, and
  the iterative exclusion of strongly correlated pairs, with each pair’s
  median over the draws. With the same settings it agrees with FastSpar
  and SpiecEasi to within the noise of the random draws (about 0.002 for
  counts in the thousands). `permutations` adds FastSpar’s p-values
  (resampled null datasets; never 0), parallel to the triangle. Like
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), it
  returns a `dist` over the features, which
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws on a
  diverging scale and `clust` clusters as `1 - r`. Draws run in
  parallel, with results identical for any `cpus`: 500 taxa in 200
  samples, 100 draws, take about 1.5 seconds on 4 cores. Its docs
  explain how it differs from
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md).

- New
  [`assoc_pairs()`](https://cmmr.github.io/ecodive/reference/assoc_pairs.md)
  lists the pairs of features whose correlation, from
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md) or
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) (or a
  correlation matrix), is at least a cutoff in absolute value, strongest
  first, as a data frame of `from`, `to` and `association`. With
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)’s
  p-values it adds `p_value` and `p_adjusted` (adjusted over all pairs
  by `p_adjust`, Benjamini-Hochberg by default), and `alpha` keeps only
  the significant pairs.

- New
  [`assoc_network()`](https://cmmr.github.io/ecodive/reference/assoc_network.md)
  turns associations between features
  ([`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), or a
  correlation matrix) into a network and analyses it as NetCoMi’s
  `netConstruct()` and `netAnalyze()` do: edges at a cutoff (and
  optionally significant), signed or unsigned edge lengths, degree,
  betweenness, closeness and eigenvector centralities of the largest
  connected component, hubs, fast greedy modularity clusters, and global
  properties (density, components, average path length, Barrat’s
  clustering coefficient, modularity, …). Written from the papers and
  NetCoMi’s documentation, it agrees with NetCoMi 1.2.0 to about 1e-13
  on 160 test networks. Shortest paths run in parallel in C, with
  results identical for any `cpus`: 3,000 features with 168,000 edges
  take about 4 seconds on 4 cores.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it on a
  circular layout.

- New
  [`assoc_network_compare()`](https://cmmr.github.io/ecodive/reference/assoc_network_compare.md)
  compares two groups’ networks as NetCoMi’s `netCompare()` does:
  differences in centralities and global properties, Jaccard indices of
  the most central features and of the hubs (with NetCoMi’s
  probabilities), and the adjusted Rand index of the clusterings with a
  permutation test. Given the count table and the groups, it tests every
  difference by permuting the samples’ groups and re-estimating both
  networks.

- New
  [`assoc_network_diff()`](https://cmmr.github.io/ecodive/reference/assoc_network_diff.md)
  tests every pair’s correlation for a difference between two groups
  with Fisher’s z, as NetCoMi’s `diffnet(diffMethod = "fisherTest")`.

- New
  [`cross_assoc()`](https://cmmr.github.io/ecodive/reference/cross_assoc.md)
  correlates every feature of one table (taxa) with every variable of
  another measured on the same samples (environmental factors,
  metabolites, clinical variables): Spearman, Pearson or Kendall’s
  tau-b, optionally partial (`given`), with p-values as SciPy’s
  `spearmanr()`, `pearsonr()` and `kendalltau()` compute them (or by
  permutation), adjusted over all pairs (Benjamini-Hochberg by default).
  Missing values in `y` are left out pair by pair, as `Hmisc::rcorr()`
  does. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
  the familiar heatmap, clustered, with stars on the significant cells,
  and [`as.matrix()`](https://rdrr.io/r/base/matrix.html) gives any
  column as a taxa x variables matrix. It agrees with SciPy,
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html), Hmisc, psych,
  ppcor and microeco to about 1e-14. Sparse tables stay sparse, and the
  correlations run in parallel in C, identical for any `cpus`.

- New
  [`feature_assoc()`](https://cmmr.github.io/ecodive/reference/feature_assoc.md)
  runs any association measure
  ([`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md),
  [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md),
  [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md),
  [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md)) by
  name, with
  [`list_feature_assocs()`](https://cmmr.github.io/ecodive/reference/feature_assocs.md)
  and
  [`match_feature_assoc()`](https://cmmr.github.io/ecodive/reference/feature_assocs.md).

- New [`qmp_scale()`](https://cmmr.github.io/ecodive/reference/qmp.md)
  and [`qmp_rarefy()`](https://cmmr.github.io/ecodive/reference/qmp.md)
  turn relative abundances into absolute abundances when each sample’s
  total microbial `load` was measured (flow cytometry, qPCR, a
  spike-in).
  [`qmp_scale()`](https://cmmr.github.io/ecodive/reference/qmp.md)
  multiplies each sample’s relative abundances by its load (Jian et
  al. 2020);
  [`qmp_rarefy()`](https://cmmr.github.io/ecodive/reference/qmp.md) is
  Vandeputte et al.’s (2017) quantitative microbiome profiling, which
  first rarefies every sample to the same number of reads per unit of
  load. It follows the authors’ `QMP.R` (BSD-3, Raes Lab), with the same
  sampling depths and targets, but uses ecodive’s own rarefaction and
  `seed`. `load` is matched to the samples by name, or given as a column
  of `data`. Both keep the table’s class and orientation, and sparse
  tables stay sparse: 2,000 samples by 5,000 taxa take under 0.1
  seconds.

- New
  [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md)
  clusters samples into community types as Arumugam et al. (2011) did:
  partitioning around medoids (PAM) on Jensen-Shannon distances (or any
  `dist`), with the number of types chosen by the Calinski-Harabasz
  index and each sample’s fit given by its silhouette. PAM, the
  silhouettes and the index are threaded C written from the published
  algorithms; medoids, assignments, objectives and silhouettes match
  [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html) and
  [`cluster::silhouette()`](https://rdrr.io/pkg/cluster/man/silhouette.html),
  and the index matches `clusterSim::index.G1()` as Arumugam et al.’s
  script calls it (`ch = "medoids"`), or is the distances’ pseudo-F
  (`ch = "centroids"`). The result, a data frame of class
  `ecodive_types` (see
  [`?ecodive_types`](https://cmmr.github.io/ecodive/reference/ecodive_types.md)),
  holds each sample’s type, as a factor to use as `groups` in any test,
  and as attributes each type’s medoid, size and mean composition (with
  its most abundant feature) and the index for every number of types
  tried. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
  the types’ mean compositions or the index;
  [`print()`](https://rdrr.io/r/base/print.html) shows the types, the
  selection and the first samples.

- New
  [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md)
  assigns vaginal microbiome samples to community state types with
  VALENCIA (France et al. 2020): each sample goes to the most similar of
  13 reference centroids (sub-CSTs I-A to V), by Yue and Clayton’s
  similarity, and its sub-CSTs are grouped into CSTs. It is a threaded C
  port of VALENCIA’s `Valencia.py` (MIT) that keeps sparse tables
  sparse; its similarities match `Valencia.py`‘s to 3e-15, with the same
  assignments. VALENCIA’s centroids ship as the dataset
  `valencia_centroids`, and revised ones can be given as `centroids`;
  `totals` gives each sample’s read count, as VALENCIA’s `read_count`
  column does. The result is an `ecodive_types` data frame with each
  sample’s sub-CST, CST, score and share of matched reads, and every
  similarity as an attribute. It warns when the centroids’ key taxa are
  missing from the table (VALENCIA matches taxa by exact name).
  `plot(x, "scores")` compares the samples’ scores with VALENCIA’s
  reference ones, as VALENCIA’s own plot does.

- New [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md)
  runs ANCOM (Mandal et al. 2015), the analysis of composition of
  microbiomes: it tests every log ratio between two features for a
  difference between `groups` (`sig_test = 'anova'`, or `'kruskal'`),
  adjusts each feature’s p-values (`p_adjust`, Holm by default), and
  declares the features whose W, the number of their significant log
  ratios, passes scikit-bio’s cutoff rule. Its `W` and decisions are
  identical to scikit-bio’s
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md) in
  every reference case, including its quirks (each feature’s family
  includes a 0 for itself; numpy’s cutoffs); the Kruskal-Wallis
  statistics are bit-identical to SciPy’s and the F statistics more
  accurate (SciPy’s cancel to 3e-12). It returns a per-feature table
  (`feature`, `W`, `reject`). Zeros follow the shared `pseudocount`
  rule, and the output of
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  needs none. The tests run in parallel C straight from the (sparse)
  table, identical for any `cpus`: all 4.5 million log ratios of 3,000
  features in 200 samples take about 3 seconds on 4 cores (10 with
  `'kruskal'`).

- New
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  finds differentially abundant features the ALDEx2 way (Fernandes et
  al. 2013, 2014), as scikit-bio’s
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md):
  each sample’s proportions are drawn `n_draws` times from their
  Dirichlet posterior (`pseudocount = 0.5` is the prior), every draw is
  clr-transformed, each feature is compared between the two sides of a
  `contrast` with Welch’s t-test, and the results are averaged. It
  returns a per-feature table with the average `t`, `log2_fc` and
  `p_value`, the conservative interval of scikit-bio (the widest bounds
  over the draws), `p_adjusted` and `reject`.
  `contrast = c(treatment, reference)` names the groups; `NA` as the
  reference pools all other samples, which is also the default (the
  first group against the rest). Given the same draws, the results are
  scikit-bio’s to about 1e-13. An optional scale model (`scale_sd`,
  following ALDEx3’s `clr.sm()`; Nixon et al. 2025) adds uncertainty
  about the groups’ total loads to every draw, so that a change in load
  isn’t mistaken for changes in the features; it is off by default
  (`scale_sd = 0`, scikit-bio’s procedure), and works best with more
  draws (ALDEx3 suggests `scale_sd = 0.5`). The draws come from a new
  Gamma sampler on PCG (Marsaglia & Tsang), with one random stream per
  sample and draw, so results are identical for any `cpus`: 500 samples
  by 5,000 features, 128 draws, take about 5 seconds on 4 cores.

- New [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  run ANCOM-BC and ANCOM-BC2 (Lin & Peddada 2020, 2024): each feature’s
  log abundance is fitted with a linear model (`formula` over `data`,
  any [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  formula, treatment contrasts), and each covariate’s bias from the
  samples’ differing sampling fractions is estimated by an EM algorithm
  and removed, giving bias-corrected log fold changes with standard
  errors, p-values (normal for ANCOM-BC, t for ANCOM-BC2), `p_adjusted`
  (`p_adjust` within each covariate) and `reject`.
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  also corrects each sample for its sampling fraction and refits, widens
  the variances (`var_quantile`), and by default treats zeros as missing
  values, including ANCOMBC’s handling of failed fits and missing
  reference levels;
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  uses the shared `pseudocount` rule. They follow scikit-bio’s
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  (themselves derived from the ANCOMBC package, whose Artistic-2.0
  notice they carry) and match them to about 1e-8 in every reference
  case, and ANCOMBC to about 1e-5 where its EM converges. The result is
  a per-feature-and-covariate table that also carries the fit
  (coefficients, variances, covariance matrices, degrees of freedom),
  for post-hoc tests to come. The fits run in parallel C, identical for
  any `cpus` and for sparse and dense tables: 500 samples by 5,000
  features take about 1 second.

- [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  gain ANCOMBC’s pseudocount sensitivity analysis, `sensitivity`: `TRUE`
  refits with pseudocounts 0.1, 0.5 and 1 added to every value (or give
  your own), and the table gains `passed_ss` (the decision is the same
  in every fit) and `reject_robust` (the decisions to use), as ANCOMBC’s
  `passed_ss` and `diff_robust`. The post-hoc tests repeat themselves on
  each fit and add the same columns, and the plots draw the robust
  decisions filled and the others hollow. It is on by default for
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  as in ANCOMBC (the main columns are still scikit-bio’s), with a
  warning when the table isn’t whole numbers (the default pseudocounts
  are on the scale of counts), and off by default for
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  for which neither reference offers it. The help page says plainly that
  without it the decisions can include many false discoveries.

- New post-hoc tests on an
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) or
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  fit, for a factor with three or more levels, as scikit-bio’s
  `ANCOMBCResult` methods:
  [`ancombc_global()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  (does the feature differ between any levels?),
  [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  (every pair of levels, with mixed directional FDR),
  [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  (every level against the reference, with a bootstrapped mdFDR screen)
  and
  [`ancombc_trend()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  (ordered patterns, increasing or decreasing by default, or your own
  constraint matrices, with a bootstrapped p-value). The deterministic
  statistics and p-values match scikit-bio’s to about 1e-8, and the
  bootstrapped ones match exactly on scikit-bio’s own null draws;
  ecodive draws its own on PCG (`seed`), one stream per feature,
  identical for any `cpus`. Comparisons are labelled with level names
  (`'severe - mild'`), and any categorical term of the fit can be
  tested, not only its `groups`. The trend test’s constrained estimates
  come from an exact QP solver written for ecodive. New
  [`structural_zeros()`](https://cmmr.github.io/ecodive/reference/structural_zeros.md)
  finds features absent from a group (Kaul et al. 2017), as scikit-bio’s
  `struc_zero()`, returning a features-by-groups logical matrix.

- New [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)
  runs any differential-abundance method by name, as
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
  runs ordinations:
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md),
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  and the new
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  and
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  ([`list_da_tests()`](https://cmmr.github.io/ecodive/reference/da_tests.md),
  [`match_da_test()`](https://cmmr.github.io/ecodive/reference/da_tests.md);
  `'Mann-Whitney'`, `'t-test'`, `'lm'` or `'ALDEx2'` work too). Its
  result is the method’s own, unchanged. It sets one convention for what
  to compare: group methods take `groups` (a column of `data`, or one
  value per sample); model methods take `formula`, with `groups` naming
  the term of interest, and `groups` alone is `formula = ~ groups`.
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  follow it: `formula` and `data` now default to `NULL`,
  `ancombc(counts, groups = grp)` fits `~ group`, and their `groups`
  need only be a term of `formula`; that it is categorical with three or
  more levels is now checked by the post-hoc tests that need it.
  [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md) and
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  results gain a `method` attribute, as
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md)’s
  has.

- New
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  and
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  test every feature of a normalized table with the classical tests that
  benchmarks found to control false discoveries (Wirbel et al. 2024;
  Pelto et al. 2025): Wilcoxon’s rank-sum test, Welch’s (or Student’s,
  `var_equal = TRUE`) t-test, and a linear model of any `formula`, each
  feature’s p-values adjusted with `p_adjust` (default `'holm'`, as the
  other DA methods). They give the results of
  [`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html),
  [`t.test()`](https://rdrr.io/r/stats/t.test.html) and `summary(lm())`
  / [`confint()`](https://rdrr.io/r/stats/confint.html) feature by
  feature, to about 1e-12, from a sparse table without a dense copy, in
  parallel C, identical for any `cpus`.
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  follows R 4.6.0’s
  [`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) on every R
  version: with ties and fewer than 50 values on each side, its p-value
  is exact given the ties (earlier R used the normal approximation
  there). `norm` chooses the values compared (`'clr'` by default;
  `'log'`, `'percent'`, `'rclr'`, whose zeros are left out, …). Their
  results print as tables and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) as volcano
  plots.

- New
  [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md)
  reruns any
  [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)
  method over the settings that have no right value, `pseudocounts`,
  zero replacements (`zeros`, through
  [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md))
  and
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)’s
  scale uncertainty (`scale_sds`), one at a time or `crossed`, and
  reports each feature’s decision in every run: one tidy table, a
  per-feature `summary` (`passed` when every run agrees, ANCOMBC’s
  `passed_ss` generalised, and `robust`), and two plots, tiles of the
  decisions and traces of the effects against a setting (for
  `scale_sds`, the display of ALDEx2’s `plotGamma()`). Every run is
  identical to the direct call with its setting, keeps the base run’s
  `seed`, and is the same for any `cpus`; ANCOM-BC’s pseudocount runs
  are the fits behind its `passed_ss`, reused when the base run made
  them.

- New [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md):
  LinDA (Zhou et al. 2022), linear models of each feature’s centred log2
  ratios with the compositional bias of every coefficient removed
  (estimated as the coefficients’ mode across features). Any `formula`,
  as
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md);
  winsorization of outliers (`outlier_pct`); zeros replaced by a
  pseudocount, by imputation proportional to sequencing depth, or
  adaptively (`zero_handling`, `corr_cut`); numeric covariates
  standardized (`standardize`); repeated measures with a random
  intercept, `~ group + (1 | subject)`, fitted by REML at the exact
  optimum, with Satterthwaite degrees of freedom (as lme4 and lmerTest,
  without them, and about 300 times faster than MicrobiomeStat’s mixed
  models). It matches MicrobiomeStat’s
  [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md) to
  about 1e-12 (fold changes, standard errors, the bias, p-values),
  except that its `'adaptive'` choice is applied (MicrobiomeStat 1.4
  reports it but adds the pseudocount regardless). Parallel C from a
  sparse table, the same for any `cpus`; one of
  [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)’s
  and
  [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md)’s
  methods, with a volcano
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- New
  [`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md):
  ALDEx3’s linear models (McGovern & Silverman 2026; Nixon et al. 2025),
  the regression form of
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md).
  Every Monte Carlo draw of the samples’ Dirichlet posteriors gets a
  scale model (`scale_model = 'clr'` or `'tss'` with uncertainty
  `scale_sd`, or measured or predicted loads, `load` with error
  `load_sd`), and each feature’s log2 abundance is regressed on any
  `formula`, with HC3 (default), HC0 or classical standard errors, or
  with a random intercept, `~ group + (1 | subject)` (REML and
  Satterthwaite, as
  [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md)).
  Estimates and p-values are averaged over the draws as ALDEx3’s
  `aldex()` and [`summary()`](https://rdrr.io/r/base/summary.html) do,
  and match ALDEx3 1.4.0 to about 1e-12 on the same draws. Its draws are
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)’s
  for the same `seed`, so with `scale_sd = 0` and two groups its
  estimate is
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)’s
  `log2_fc`. Parallel C, the same for any `cpus`; one of
  [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)’s
  methods (also as `"ALDEx3"`) and
  [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md)’s
  (`scale_sds`), with volcano, effect, MA, forest and rank plots.

- New
  [`maaslin3()`](https://cmmr.github.io/ecodive/reference/maaslin3.md):
  MaAsLin 3 (Nickols et al. 2026), ported from maaslin3 1.4.0. Each
  feature gets two models of any `formula`: its abundance where present
  (log2 proportions by default, zeros left out) and its prevalence (an
  augmented logistic regression, finite under separation). Abundance
  coefficients are tested against their median over the features (the
  median comparison, a compositional correction), a joint p-value
  combines the two models, and maaslin3’s warnings become per-model
  flags (`warn_prevalence` and its misfit checks). Repeated measures
  take a random intercept, `~ group + (1 | subject)`: REML with
  Satterthwaite degrees of freedom for abundance, a Laplace logistic
  mixed model (as lme4’s `glmer()`) for prevalence, both at their exact
  optimum, or the grouping as a fixed effect (`small_random_effects`).
  It matches maaslin3 to about 1e-10 with its median comparison’s draws
  given (the draws come from PCG streams, keyed by `seed`), reports
  every feature, and runs in parallel C from a sparse table, the same
  for any `cpus` (200 samples x 2,000 features in about a second on one
  core). One of
  [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)’s
  methods, with volcano, MA, forest and rank plots of either model.

- New
  [`sourcetracker()`](https://cmmr.github.io/ecodive/reference/sourcetracker.md)
  estimates what share of each sink sample’s community came from each of
  several source environments, and from an unknown source, with
  SourceTracker’s Gibbs sampler (Knights et al. 2011), ported from
  SourceTracker2. One `sources` vector gives each sample’s environment,
  with `NA` for the sinks; `loo = TRUE` tracks each environment against
  the others, and `features = TRUE` also returns each sink’s reads by
  source and feature. The sampler runs in threaded C, with one random
  stream per sink and restart, so results are reproducible and the same
  for any `cpus` (SourceTracker2’s parallel runs are neither). It
  returns a sinks x sources matrix of proportions (class
  `ecodive_source`) with their standard deviations, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  usual stacked bar chart.

- New
  [`decontam()`](https://cmmr.github.io/ecodive/reference/decontam.md)
  and
  [`decontam_not()`](https://cmmr.github.io/ecodive/reference/decontam.md)
  identify contaminant features with decontam’s frequency and prevalence
  tests (Davis et al. 2018), ported from decontam 1.32.0’s
  `isContaminant()` and `isNotContaminant()` (Artistic-2.0): a
  contaminant’s frequency falls as the sample’s DNA concentration
  (`conc`) rises, and it is more prevalent in the negative controls
  (`neg`). Every decontam method (`'frequency'`, `'prevalence'`,
  `'combined'`, `'minimum'`, `'either'`, `'both'`), batches with their
  three ways of combining scores, and decontam’s scores and calls to
  rounding. `conc`, `neg` and `batch` can be columns of `data`. Threaded
  C over a sparse table, the same for any `cpus`.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  score histogram (for choosing `threshold`) or each feature’s
  prevalence in the controls against the samples, and
  [`decontam_frequency_plot()`](https://cmmr.github.io/ecodive/reference/decontam_frequency_plot.md)
  a feature’s frequency against `conc` with both fitted models.

- New
  [`neutral_model()`](https://cmmr.github.io/ecodive/reference/neutral_model.md)
  fits Sloan et al.’s (2006) neutral community model as Burns et
  al. (2016) do: each feature’s occurrence frequency against its mean
  relative abundance, with the immigration probability `m` fitted by
  least squares (exactly, without an upper bound; a warning when it
  exceeds 1), its profile confidence interval and standard error,
  R-squared and RMSE, and the AIC and BIC of the neutral, binomial and
  Poisson models. Each feature is partitioned as above, within or below
  the model’s Wilson band (Wilson 1927). It follows neufit (Sieber et
  al. 2019, BSD-2); `pool` gives a separate source metacommunity, and
  `groups` (or a column of `data`) fits each group of samples on its
  own. It matches a 40-digit computation of the fit, Burns et al.’s code
  (via iCAMP’s `snm()` and tyRa’s `fit_sncm()`) to their optimizers’
  tolerances, and neufit’s predictions and bands to rounding. The
  per-feature sums run in threaded C over a sparse table, the same for
  any `cpus`. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  draws the fitted curve and its band over the features
  (`type = "binomial"` adds the binomial model’s curve).

- Longitudinal analyses, ported from q2-longitudinal 2026.10.0 (Bokulich
  et al. 2018, BSD-3), for repeated measures of the same subjects over
  time (`subjects` and numeric `states` given per sample or as columns
  of `data`):

  - New
    [`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md)
    gives each subject’s change in a variable from one state to the
    next, or from a `baseline` state, and
    [`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md)
    the distance between its two samples, from any `dist`. A vector of
    values, a column of `data`, or every feature of a count table (after
    `norm`) at once. Pairs follow q2’s rules (consecutive states over
    all samples; `replicates` = `"error"`, `"drop"` or a seeded
    `"random"` pick), and match q2’s differences and distances exactly.
  - New
    [`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md)
    computes q2’s volatility chart: each group’s mean at each state with
    Vega’s bootstrap confidence interval (on ecodive’s seeded
    generator), the global mean and the 2 and 3 SD control limits, and
    each subject’s trajectory.
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
    control chart.
  - New [`nmit()`](https://cmmr.github.io/ecodive/reference/nmit.md)
    gives the distances between subjects in how their features vary
    together over time (Zhang et al. 2017, NMIT): each subject’s
    Kendall, Pearson or Spearman correlations of its features, and the
    Frobenius or nuclear norm of their differences, in threaded C, the
    same for any `cpus`. It matches q2’s distances to about 1e-15;
    [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
    on the result is NMIT’s test.
  - New
    [`feature_volatility()`](https://cmmr.github.io/ecodive/reference/feature_volatility.md)
    runs q2’s feature-volatility pipeline: a
    [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
    regression of the state on the features, on a random 80% of the
    samples (`test_size`), with its accuracy on the rest, optionally
    after q2’s recursive feature elimination (scikit-learn’s RFECV;
    `rfe = TRUE`, `rfe_step`, `rfe_folds`); the `n_features` most
    important features with q2’s statistics (cumulative increase and
    decrease of the mean between states, net change, global mean,
    median, variance, SD and CV); and each one’s volatility chart
    numbers. Given q2’s `random_state`, the split, RFECV’s scores and
    selection, importances and predictions match q2’s exactly, and the
    statistics to about 1e-15.
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
    importances beside the chosen feature’s volatility chart.

- New
  [`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md)
  and
  [`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md)
  remove batch effects (study, sequencing run, extraction kit) from a
  table, keeping the effects of the covariates named in `formula`.
  [`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md)
  ports MMUPHin’s `adjust_batch()` (Ma et al. 2022; MMUPHin 2.0.0, MIT):
  zero-inflated empirical Bayes adjustment of log relative abundances,
  re-closed to each sample’s total, with zeros left as zeros
  (`zero_inflation = FALSE` fits every value).
  [`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md)
  ports sva 3.60.0’s `ComBat()` (Johnson et al. 2007; Artistic-2.0) for
  continuous values (or `norm = 'clr'` / `'log'`), with parametric or
  non-parametric priors, `mean_only` and `ref_batch`; its non-parametric
  priors are weighed on the log scale, so large batches don’t underflow
  to `NaN` as in sva. Both match their references to rounding (MMUPHin’s
  rounded counts exactly), in threaded C over a sparse table, the same
  for any `cpus`. `batch` can be a column of `data`. The result is the
  corrected table (dense, in the input’s orientation) of class
  `ecodive_batch`, with the batch effects before and after shrinkage as
  attributes; [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  draws MMUPHin’s diagnostics (the shrinkage of the batch effects, and
  each feature’s batch means before and after).

- New
  [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md)
  grows random forests (Breiman 2001) for classification (a factor,
  character or logical `y`) or regression (a numeric `y`), and gives
  q2-sample-classifier’s workflow (Bokulich et al. 2018):
  cross-validated predictions and class probabilities (stratified folds,
  or `folds`), accuracy (with q2’s baseline and ratio), ROC curves and
  AUC (with q2’s micro and macro averages), out-of-bag predictions and
  scores, and impurity and out-of-bag permutation importances. Vectors
  of `max_features` or `min_samples_leaf` make a grid, chosen inside
  each training set by out-of-bag error. The trees are an exact port of
  scikit-learn 1.9.1’s (BSD-3): given scikit-learn’s random draws, every
  node, importance and out-of-bag prediction matches, and q2’s
  cross-validated predictions, probabilities, importances and ROC curves
  match exactly. ecodive draws from PCG, one stream per tree, so a
  forest is the same for any `cpus`; threaded C over sparse or dense
  tables (float32, as scikit-learn).
  [`predict()`](https://rdrr.io/r/stats/predict.html) matches new
  samples’ features by name, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  importances, ROC curves, confusion matrix, predicted against observed
  values, or the grid’s out-of-bag errors.

- New
  [`moran_i()`](https://cmmr.github.io/ecodive/reference/moran_i.md):
  Moran’s I test of spatial autocorrelation, whether a variable (a
  taxon’s abundance, alpha diversity) is more alike at nearby sites,
  ported from PySAL esda 2.10.0’s `Moran` (BSD-3). `space` gives the
  sites’ positions as spatial weights (a matrix), geographic distances
  (a `dist`) or coordinates (both made into inverse-distance weights),
  and `style` row-standardises them (the default, as esda and spdep),
  makes them binary or variance-stabilising, or leaves them as given. It
  gives the expected value, the standard deviations and analytic
  p-values under normality and under randomisation, and a permutation
  test in threaded C over the packed weights, the same for any `cpus`. I
  and its moments match esda, spdep’s `moran.test()` and ape’s
  `Moran.I()` to about 1e-15, and the permuted statistics esda’s, given
  the same permutations.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  permutation distribution.

- New [`dbmem()`](https://cmmr.github.io/ecodive/reference/dbmem.md):
  distance-based Moran’s eigenvector maps (dbMEM, formerly PCNM; Borcard
  & Legendre 2002, Dray et al. 2006), uncorrelated spatial variables
  from broad to fine scales, to use as explanatory variables in
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md) or
  [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md).
  `space` gives the samples’ coordinates or geographic distances, as for
  [`moran_i()`](https://cmmr.github.io/ecodive/reference/moran_i.md);
  `threshold` is the truncation distance (by default the longest edge of
  the minimum spanning tree), and `autocor` keeps the MEMs with positive
  (the default), negative, non-zero or any Moran’s I, as adespatial’s
  `MEM.autocor`. The MEMs, their eigenvalues and Moran’s I (attributes
  `values` and `moran_i`) match adespatial’s
  [`dbmem()`](https://cmmr.github.io/ecodive/reference/dbmem.md) to
  about 1e-11 (up to sign), and its truncation vegan’s `pcnm()`. The
  eigenvectors come from a Householder reduction and the QL algorithm in
  threaded C, without the BLAS, the same for any `cpus`. The result is a
  matrix of class `ecodive_dbmem`;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
  MEM’s Moran’s I, or maps MEMs over the samples’ coordinates.

- New [`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md):
  LEfSe, LDA effect size biomarker discovery (Segata et al. 2011),
  ported from LEfSe 1.1.2: a Kruskal-Wallis test of each feature, the
  pairwise Wilcoxon tests between subgroups (`subgroups`,
  `min_subgroup`, `same_subgroups`, `one_against_one`), and the
  bootstrapped LDA score, with each sample rescaled to `total` (counts
  per million by default). Given LEfSe’s own jitter and subsamples, the
  scores match LEfSe’s to rounding, more features than samples included;
  the tests match R’s
  [`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html) and
  coin’s `wilcox_test()`. New `p_adjust` (default `'holm'`; LEfSe
  adjusts nothing, `'none'`). Threaded C over a sparse table, the
  bootstraps with one random stream each, so the same for any `cpus`.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws LEfSe’s
  bar chart of LDA scores. It joins
  [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)
  (`"lefse"`, alias `"lefser"`) and
  [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md)
  (zero replacements). The docs say why the plain tests are preferable.

- New
  [`indicator_species()`](https://cmmr.github.io/ecodive/reference/indicator_species.md):
  indicator species analysis (Dufrene & Legendre

  1997. with group combinations, De Caceres et al.’s (2010) multilevel
        pattern analysis. Each feature’s association with every group
        and combination of groups (`min_order`, `max_order`,
        `combinations`, by name or id) by the group-equalized indicator
        value (`index = 'indval_g'`, the default), the plain one
        (`'indval'`, with its specificity `a` and fidelity `b`), or the
        point-biserial correlation (`'r'`, `'r_g'`, with `negative`
        associations), the best combination, and a permutation test of
        it against the largest value over the combinations (`strata`, or
        a supplied matrix). Matches indicspecies’ `multipatt()` and
        `strassoc()` and labdsv’s `indval()` to rounding, with identical
        p-values given the same permutations; ties between combinations
        go the same way every time (multipatt breaks them at random).
        New `p_adjust` (default `'holm'`; multipatt adjusts nothing,
        `'none'`). Threaded C over a sparse table, permutations in fixed
        blocks, so the same for any `cpus`.
        [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
        the indicators’ values as bars, or their groups’ means as a
        heatmap. It joins
        [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md)
        (`"indicator_species"`, aliases `"indval"` and `"multipatt"`).

- New `norm = 'log'`, wherever `norm` is offered (beta diversity
  metrics, [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md),
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md),
  …): the log of each sample’s proportions after adding `pseudocount`
  (log-TSS), with the same pseudocount rule as `'clr'`, of which it is a
  shift by one constant per sample.

- `rarefy(depth =)` now also takes one depth per sample (matched to the
  samples by name when named). `drop` and `warn` then compare each
  sample with its own depth. A single depth gives exactly the same
  results as before.

- New
  [`beta_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md)
  and
  [`alpha_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md)
  rarefy a table many times (`times = 100` by default), compute any beta
  or alpha diversity metric from each rarefied table, and average the
  results: the mean (as `vegan::avgdist()` and mothur), the median, or
  for beta diversity the medoid (the matrix of the rarefaction closest
  to all the others, as in q2-boots). `sd = TRUE` adds each value’s
  standard deviation over the rarefactions. Iteration t is exactly
  `rarefy(counts, depth, seed = seed + t - 1)`, so results are
  reproducible and the same for any `cpus`; the mean needs memory for
  only a few distance matrices, and the median and medoid work through
  the pairs of samples in blocks. Titles say so: an ordination of the
  result is a “Principal Coordinates Analysis of Bray-Curtis
  Dissimilarity (rarefied 100×)”.

- New
  [`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md):
  alpha diversity against sequencing depth, one row per sample, metric
  and depth, with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Richness
  (`'observed'`) is computed exactly in closed form, as
  `vegan::rarefy()` and `rarecurve()` do (Hurlbert’s expected richness
  and Heck et al.’s standard deviation, without vegan’s rounding error);
  any other alpha metric, or several at once, by rarefying `times` times
  at each depth (as QIIME 2’s `alpha-rarefaction`, whose depth grid it
  uses), with the same seed at every depth, so that each rarefaction’s
  subsamples are nested. Every row also has the expected sample coverage
  (Chao & Jost 2012), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws depth,
  coverage-based or sample-completeness curves, by group, with a line at
  a planned rarefaction depth.

- New
  [`composition_barplot()`](https://cmmr.github.io/ecodive/reference/composition_barplot.md)
  draws stacked bars of relative abundance, one per sample and never
  averaged, ordered so they can be read: samples by their Aitchison PCoA
  axis 1 (or a clustering, `clust` as in
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)),
  within panels by `groups`; the most abundant features on the baseline.
  With `feature_groups` (e.g. each genus’s phylum) it shades the bars
  hierarchically, after microshades (Dahl et al. 2022): one hue per
  group, shades for its top members, the group’s own “Other” lightest.
  The subtitle says how much “Other” holds. `heatmap = TRUE` adds the
  drawn features’ centred clr values below, aligned with the bars. It
  returns the drawn proportions, colours and sample order invisibly;
  sparse tables stay sparse.

- New
  [`ternary_plot()`](https://cmmr.github.io/ecodive/reference/ternary_plot.md)
  draws three-part compositions in a ternary diagram, and returns the
  drawn coordinates invisibly. In samples mode (`parts`), each sample is
  a point at three features, or three sums of features (“phyla”; one
  part may be `NA`, for the rest). With `center = TRUE` the points are
  centred by perturbation (von Eynatten et al. 2002), and the gridlines
  become curves that still mark the original proportions.
  `ellipses = "confidence"` or `"prediction"` adds each group’s
  Hotelling region (F-based) from ilr coordinates, mapped back into the
  triangle. In features mode (`groups` with 3 levels), each feature is a
  point at its share of the three groups’ mean relative abundances,
  sized by abundance and coloured by the group that holds most of it.
  Sparse tables stay sparse.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an
  ordination takes `scaling` and `links`, for compositional biplots
  (Aitchison & Greenacre 2002). For
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) and
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md):
  `scaling = 'form'` (the default, as before), `'covariance'` (samples
  of unit variance and features at the loadings times `sqrt(eig)`, so
  that in a clr PCA the rays’ lengths approximate the clr values’
  standard deviations), or `'contribution'` (Greenacre 2013: the unit
  loadings times a stated round number, with the features contributing
  more than average highlighted). In a covariance biplot,
  `links = c("a", "b", "c")` joins the features’ tips with dashed links
  labelled with their length, which approximates (from below) the
  standard deviation of their log ratio; drawing them warns when the
  plane explains under half the variance, or when the PCA is not of clr
  values. For
  [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
  [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) and
  [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md),
  `scaling = 1`, `2` (the default, as stored) or `3` draws `vegan`’s
  scalings, matching `vegan::scores(scaling =)`.
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)’s result
  gains `norm` and `scale` attributes, which qualify its loadings.

- [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) takes
  a new last argument, `coverage`, to rarefy each sample to the smallest
  depth whose expected coverage reaches it (coverage-based rarefaction,
  Chao & Jost 2012), instead of to one `depth`.

- Eleven new alpha diversity metrics, ported from scikit-bio’s
  `skbio.diversity.alpha` (which QIIME 2 runs) and matching it to
  rounding: Pielou’s evenness
  [`pielou()`](https://cmmr.github.io/ecodive/reference/pielou.md)
  (QIIME 2’s `evenness_vector`),
  [`heip()`](https://cmmr.github.io/ecodive/reference/heip.md),
  [`simpson_evenness()`](https://cmmr.github.io/ecodive/reference/simpson_evenness.md)
  and
  [`mcintosh_evenness()`](https://cmmr.github.io/ecodive/reference/mcintosh_evenness.md);
  Good’s coverage
  [`goods_coverage()`](https://cmmr.github.io/ecodive/reference/goods_coverage.md);
  the Hill numbers
  [`hill()`](https://cmmr.github.io/ecodive/reference/hill.md) and the
  Rényi and Tsallis entropies
  [`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
  [`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md), of
  order `q` (default 2, as scikit-bio; `Inf` allowed); and Simpson’s and
  Strong’s dominance indices
  [`dominance()`](https://cmmr.github.io/ecodive/reference/dominance.md)
  and [`strong()`](https://cmmr.github.io/ecodive/reference/strong.md)
  and the Gini index
  [`gini()`](https://cmmr.github.io/ecodive/reference/gini.md). All run
  through
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
  which also accepts scikit-bio’s names (`'pielou_e'`, `'simpson_d'`,
  `'observed_features'`, `'enspie'`, …).

- [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  and
  [`alpha_div_rarefied()`](https://cmmr.github.io/ecodive/reference/rarefied_diversity.md)
  take the new `q` and `bias_corrected` arguments, before `cpus`.

- **Changed results:**
  [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md) is now
  the bias-corrected estimator, `S + F1 (F1 - 1) / (2 (F2 + 1))`, by
  default, as scikit-bio and QIIME 2: it is never `Inf` or `NaN`.
  `chao1(bias_corrected = FALSE)` gives the classic `S + F1^2 / (2 F2)`
  when a sample has singletons and doubletons, and the bias-corrected
  value otherwise, as scikit-bio.

- **Changed matching:** `alpha_div(x, 'gini')` is now the new Gini
  index; it used to match the Gini-Simpson index, which is `'simpson'`.
  Some short prefixes that matched one metric now match several, and
  give an error listing them (e.g. `'mc'`: `mcintosh` and
  `mcintosh_evenness`).

- Documentation:
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md)
  and
  [`robust_aitchison()`](https://cmmr.github.io/ecodive/reference/robust_aitchison.md)
  discuss sequencing depth (Schloss 2026) and show how to check an
  ordination against library size and a rarefied analysis;
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  recommends a `scale_sd` sensitivity check;
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md),
  [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md),
  [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) and
  [`philr()`](https://cmmr.github.io/ecodive/reference/philr.md) note
  that log-ratio transforms don’t improve machine learning;
  [`qmp_scale()`](https://cmmr.github.io/ecodive/reference/qmp.md)
  explains predicted loads. The reference to Nixon et al. (2025) now has
  the published title.

- [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc_pairwise()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  and
  [`ancombc_dunnett()`](https://cmmr.github.io/ecodive/reference/ancombc_global.md)
  results gain a `log_base` attribute (`exp(1)`: their `log_fc` is a
  natural log), and
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  results one of `2` (`log2_fc`), as
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  and
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  results have. Their plots label log fold changes from it.

- The plots of differential-abundance results take `type`: besides the
  volcano plot (the default), `"effect"`, each feature’s difference
  against its dispersion within the groups, with lines where the two are
  equal (Gloor et al. 2016; ALDEx2’s `"MW"` plot); `"ma"`, against its
  mean abundance (ALDEx2’s `"MA"`); and `"forest"`, the effects and
  intervals of the `n_features` most significant features. For
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  and
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  (all four types) and
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  (volcano and MA). The effect plot of a model result needs a two-level
  comparison (a covariate coded 0/1). Titles name the plot and the
  method (“Effect plot of ANCOM-BC2”), and subtitles the comparison.
  `type` is the second argument of the
  [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md),
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  and
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  methods, before `labels`, and the third of the
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
  [`da_lm()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  ones, after `covariate`.

- [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md)
  results gain three columns, for these plots: `pooled_sd`, the pooled
  within-group standard deviation, `mean`, the mean clr abundance (both
  averaged over the draws, in log2 units), and `effect`, each draw’s
  difference over its pooled standard deviation, averaged (Cohen’s d).
  `pooled_sd` and `effect` are ALDEx3’s `aldex.effect()` statistics and
  match them for the same draws.
  [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md) and
  [`ancombc2()`](https://cmmr.github.io/ecodive/reference/ancombc.md)
  results gain `resid_sd` and `mean` attributes (each feature’s residual
  standard deviation and mean clr), and
  [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)
  results a `mean` column, as
  [`da_ttest()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md)’s.

- [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) of a
  slam `simple_triplet_matrix` now draws exactly the same observations
  as the same table in any other class, for the same `seed`. Before,
  slam’s samples were seeded by their 1-based positions, so its draws
  differed from those of a `matrix` or `dgCMatrix`. Rarefying slam input
  therefore gives different (equally random) results than before.

- New
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
  returns a count table after any of the transformations that
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) and
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) apply
  with their `norm` argument (`'percent'`, `'chord'`, `'hellinger'`,
  `'binary'`, `'clr'`, `'rclr'`), with exactly the values they use. It
  keeps the table’s orientation, dimnames and form: sparse Matrix and
  slam input stays sparse. `norm = 'clr'` gives a dense matrix, since
  every zero gets a value; the functions with a `norm` argument avoid
  that copy.

- New `norm = 'hellinger'` option, wherever `norm` is offered (beta
  diversity metrics,
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
  [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md),
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) and
  [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)):
  the square root of each sample’s relative abundances. Euclidean
  distances of the result are Hellinger distances, and
  `pca(counts, norm = 'hellinger')` is the usual Hellinger PCA. Negative
  values are an error.

- [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md) and
  [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) now
  name `counts`, not `dm`, when `groups`, `strata` or the rows of `data`
  don’t match the number of samples.

- Method names are matched more reliably by
  [`match_metric()`](https://cmmr.github.io/ecodive/reference/metrics.md),
  [`match_dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ords.md),
  [`match_dist_test()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  and
  [`match_count_ord()`](https://cmmr.github.io/ecodive/reference/count_ords.md),
  and so by
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md),
  [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md),
  [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md),
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
  and
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md):
  a name that is exactly a method’s ID or alias, in any case, now
  selects it even when it also begins another method’s name. For
  example, `"Jensen"` selects `jensen` rather than stopping as ambiguous
  with `jsd` (Jensen-Shannon divergence), and `dist_ord(dm, "Mantel")`
  now says to use
  [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md).

- New [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) runs a
  principal component analysis of a count table, after any of
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)‘s
  normalizations (e.g. `norm = 'clr'` for an Aitchison PCA), optionally
  scaling features to unit variance (`scale = TRUE`). It gives the same
  scores, loadings, and variances as
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) on the normalized
  table (up to each axis’s sign), but computes only the first few axes,
  in parallel C, straight from sparse input without making it dense. Its
  scores equal
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md)’s on the
  matching Euclidean or Aitchison distances; its `loadings` attribute
  holds the features’ weights. Results are identical for any number of
  `cpus`.

- New [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) runs
  a robust Aitchison PCA (Martino et al. 2019) straight from a count
  table: a robust centered log-ratio transform that treats zeros as
  missing, OptSpace matrix completion, and a PCA of the completed
  matrix. It returns sample scores, feature loadings, and the fit’s
  progress in the shared ordination format (with `eig` as variances, on
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md)’s scale),
  and optionally gemelli’s distances between samples as a `dist` object
  (`distances = TRUE`). A C port of `gemelli`’s
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), it gives
  the same results to within rounding error, but works only with the
  non-zero counts, so sparse tables stay sparse and nothing of size
  samples by features is formed. It takes every count-table class that
  [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) takes, and
  its results are identical for any number of `cpus`.

- New
  [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md)
  runs a phylogenetic robust Aitchison PCA (Martino et al. 2022),
  gemelli’s `phylogenetic_rpca()`:
  [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) of a
  table expanded to every node of a tree, each node holding the counts
  of the tips below it weighted by its branch length, so that loadings
  describe clades as well as tips. It cuts the tree down to the tips
  with counts as gemelli does (with its `min_depth` filter), returns the
  tree it used with the nodes named as gemelli names them, and
  optionally the node counts (`node_counts = TRUE`). A C port of gemelli
  0.0.13 (with the tree code of scikit-bio and iow it calls), it gives
  gemelli’s results to within rounding error; the node counts are summed
  in parallel from each sample’s non-zero tips, without a dense table.
  It joins
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
  (`count_ord(counts, "rpca_phylo", tree = tree)`), and its results are
  identical for any number of `cpus`.

- New [`philr()`](https://cmmr.github.io/ecodive/reference/philr.md)
  computes the phylogenetic isometric log-ratio transform (PhILR;
  Silverman et al. 2017): one balance per internal node of a tree,
  between its two child clades, with philr’s part weights (`'uniform'`,
  `'gm.counts'`, `'enorm'`, `'anorm'` and their products with
  `'gm.counts'`, or a vector) and ILR weights (`'uniform'`, `'blw'`,
  `'blw.sqrt'`, `'mean.descendants'`, or a vector). Written from the
  paper and checked against philr 1.28.0’s outputs, it gives philr’s
  coordinates, in philr’s orientation and order, to within 1e-14 for
  every weighting. Tips that aren’t features are pruned from the tree,
  and polytomies (and an unrooted tree’s root) are resolved with
  zero-length branches as `ape::multi2di(random = FALSE)` does, with a
  warning; the tree used is returned. Each sample costs O(D), summing
  over clades from the tips up, with no basis matrix: 2,000 samples by
  5,000 taxa (sparse) take 0.2 seconds, where philr takes 160. Results
  are identical for any number of `cpus`. New
  [`philr_distance()`](https://cmmr.github.io/ecodive/reference/philr_distance.md),
  a phylogenetic beta metric in
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  (`beta_div(x, "philr", tree = tree)`), is the Euclidean distance
  between PhILR coordinates with the paper’s weights
  (`'enorm.x.gm.counts'` and `'blw.sqrt'`); with uniform weights, PhILR
  distances are Aitchison distances.

- New [`dpcoa()`](https://cmmr.github.io/ecodive/reference/dpcoa.md)
  runs a double principal coordinate analysis (DPCoA; Pavoine, Dufour &
  Chessel 2004), ade4’s
  [`dpcoa()`](https://cmmr.github.io/ecodive/reference/dpcoa.md) and
  phyloseq’s `ordinate(method = "DPCoA")`: samples and taxa in one
  space, the taxa placed by the square roots of their distances along a
  tree (or by any Euclidean dissimilarity among them, `taxa_dm`), each
  sample at the abundance-weighted centre of its taxa. It returns the
  sample coordinates, the taxa’s on the same axes (`taxa`), and each
  sample’s Rao quadratic entropy (`rao`). Written from the paper and
  checked against ade4 1.7-24’s and phyloseq 1.46.0’s outputs, it
  matches them to within 1e-12. With a tree, no taxa-by-taxa matrix is
  formed: block Lanczos runs on products computed from the sparse table
  with sums up and down the tree, so each costs time proportional to the
  non-zero counts plus the tree’s size (300 samples by 1,000 taxa: 0.02
  seconds, where ade4 takes 250; 5,000 by 20,000, sparse: 0.5 seconds).
  Tips that aren’t taxa are dropped, and polytomies and unrooted trees
  are accepted as they are (they don’t change path lengths). Taxa absent
  from every sample, which ade4 refuses, are kept with weight 0 and
  placed on the axes. It joins
  [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
  (`count_ord(counts, "dpcoa", tree = tree)`), and its results are
  identical for any number of `cpus`.

- New
  [`blas_threads()`](https://cmmr.github.io/ecodive/reference/blas_threads.md)
  gets or sets how many threads R’s BLAS/LAPACK library uses, while R is
  running. It works with OpenBLAS, Intel MKL, BLIS, and FlexiBLAS,
  returns `NA` where the count can’t be read (e.g. R’s reference BLAS or
  Apple’s Accelerate), and never raises an error, warning, or message.
  Setting is limited to what is safe: the count is never raised above
  the library’s own, and OpenMP’s thread count is left as it is. No
  other ecodive function calls it.

- Sparse counts with samples in rows (a `dgCMatrix` with `margin = 1`),
  or as triplets (`dgTMatrix`, slam’s `simple_triplet_matrix`), are now
  imported with counting sorts in linear time, instead of a quicksort: 4
  s became 0.06 s for a 5000 x 4000 `dgCMatrix` with 2.5 million
  non-zero values, and 26 s became 0.3 s for 9.4 million. Every alpha
  and beta diversity function benefits, and results are unchanged.
  Out-of-range indices in these inputs are now an error.

- Sparse counts: explicitly stored zeros are now ignored, and repeated
  entries of a `dgTMatrix` summed, as the Matrix package defines them,
  so every function gives the same results as for the dense matrix
  (before, `norm = 'binary'` counted a stored zero as present).
  Malformed sparse objects (wrong index types or lengths, bad column
  pointers) are now an error rather than read out of bounds, and tables
  with more than ~268 million non-zero values no longer overflow an
  internal size.

- Threads now come from a pool that lasts the whole session, and jobs
  are shared out to whichever threads are free. On a machine busy with
  other processes, multithreaded calls are no longer much slower than
  single-threaded ones
  ([`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) of 2,000
  samples: 0.8 s with 6 threads, where restarting threads for every step
  took 16 s), and small jobs no longer start threads at all. Results are
  unchanged, and still the same for any `cpus`.

- [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md) no
  longer uses the parallelly package, and ecodive no longer imports
  parallel. It reads the process’s CPU affinity, cgroup CPU quota (or
  Windows job object CPU cap) and core counts from the operating system,
  and follows the same job schedulers (Slurm, PBS, LSF, SGE, PJM,
  HTCondor), options, environment variables and
  `parallelly.availableCores.*` settings as
  `parallelly::availableCores()`, giving the same counts, in about 25
  microseconds rather than a few milliseconds. In a worker process (one
  forked by `mclapply()` and the like, or a worker of a socket cluster,
  mirai or clustermq), it returns 1. `options(ecodive.cpus = n)` or
  `R_ECODIVE_CPUS` caps ecodive’s threads.

- Every function’s `cpus` is now capped at the cores the process may
  use: its affinity, quota, scheduler allocation, `R CMD check`’s limit
  of 2 and `ecodive.cpus`. Settings such as `mc.cores` only lower the
  default,
  [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md), so
  that more threads can still be asked for, e.g. in a worker. Asking for
  more threads than the cap costs nothing; results are the same for any
  `cpus`.

- Fixed `NA` in an integer or logical count table being read as
  -2147483648: every alpha and beta diversity function now treats it as
  `NA`, as it does in a double table. Before,
  [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) returned
  -1 for every pair with that sample.

- Fixed the `pairs` argument of the beta diversity functions, which
  counted the columns of `counts` as samples even when samples are rows
  (`margin = 1`, the default). Out-of-range pairs could crash R, and a
  table with one feature was rejected as having fewer than two samples.

- Fixed [`faith()`](https://cmmr.github.io/ecodive/reference/faith.md),
  the UniFrac family and every other function taking a `tree`, which
  failed when samples were columns (`margin = 2`) and the tree had tips
  that were not in the table.

## ecodive 2.3.2

- Fixed
  [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md),
  and every beta metric called with `norm = 'clr'`, returning `NaN` (or,
  for
  [`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md),
  wrong values) for every sample pair. Two bugs were involved:

  - Data with no zeros (pseudocount 0) hit `log(0) * 0`, which is `NaN`.
  - The pseudocount was truncated to an integer in C, so values below 1,
    including the automatic `min(counts[counts > 0]) / 2`, became 0, and
    values like 1.5 were silently rounded down.

- [`gower()`](https://cmmr.github.io/ecodive/reference/gower.md) with
  `norm = 'clr'` now computes each feature’s range, and the distance for
  features absent from both samples, from that feature’s CLR values
  instead of treating absent features as 0.
  [`gower()`](https://cmmr.github.io/ecodive/reference/gower.md) with
  negative inputs (e.g. `norm = 'rclr'`) now also uses the full
  min-to-max range.

## ecodive 2.3.1

CRAN release: 2026-11-05

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  now handles quoted labels per the Newick standard:

  - `''` inside a quoted label is decoded as one literal apostrophe, so
    `'it''s'` is read as `it's`.
  - Square brackets inside quotes are literal characters rather than
    comments, so labels like `'o__[Saprospirales]'` are no longer
    truncated.
  - Newlines and whitespace inside quotes are preserved.
  - Comments (`[...]`) outside of quotes may contain apostrophes.

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  raises an informative error, with the position of the problem, for
  malformed input such as an unclosed quote, comment, or parenthesis.
  Previously these could return a garbled tree.

- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)’s
  parser is now linear-time and non-recursive, making it about 4x faster
  on a typical 100,000-tip tree and much faster on very deep trees.

## ecodive 2.3.0

- Added `ecodive` JOSS citation information.
- New
  [`robust_aitchison()`](https://cmmr.github.io/ecodive/reference/robust_aitchison.md)
  beta diversity function.
- Switch from `testthat` to dependency-free `tinytest` for test suite.
- `read_tree` fixes for edge cases.

## ecodive 2.2.6

CRAN release: 2026-04-14

- When using CLR normalization:

  - Warns when zeros are detected and `pseudocount` is `NULL`.
  - Then continues with `pseudocount <- min(counts[counts>0]) / 2`.

### Breaking Changes

- `norm` now defaults to `"none"` instead of `"percent"`.
- `norm` parameter was removed from some beta diversity functions.

## ecodive 2.2.4

### Breaking Changes

- The [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md)
  function has been simplified:

  - `drop = FALSE` now returns samples with original counts when sum \<
    depth.
  - `drop = TRUE` no longer drops zero-sum features - only zero-sum
    samples.
  - `depth` can no longer be a fraction (only positive integers now).
  - `n_samples` parameter has been removed.
  - `warn` parameter has been added.

## ecodive 2.2.3

- Coerce counts to double before rarefaction.

## ecodive 2.2.2

CRAN release: 2026-01-16

- Drop explicit zeros from sparse matrices after rarefaction.

## ecodive 2.2.1

CRAN release: 2025-12-12

- Corrected CRAN `rchk` issue.

## ecodive 2.2.0

CRAN release: 2025-12-03

- Dedicated rarefaction backends for dense and sparse matrices.
- Auto-conversion to col-major sparse matrix for distance/diversity
  calculations.

## ecodive 2.1.0

CRAN release: 2025-09-26

### Breaking Changes

- `rescale` parameter is superseded by `norm` parameter.

### Fixes

- [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  and
  [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md)
  will only forward applicable parameters.
- No longer crashes when `pairs = integer(0)`.
- `pseudocount = NULL` correctly selects smallest non-zero count.

## ecodive 2.0.0

CRAN release: 2025-09-21

### Breaking Changes

- Input matrix is now samples as rows.
- Removed `weighted` parameter.

### Added New Diversity Methods

- Updated Alpha Options:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`

- Updated Beta Options:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`

## ecodive 1.0.0

CRAN release: 2025-08-22

- This is the first release of ecodive.
