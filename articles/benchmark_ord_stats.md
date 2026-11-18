# Benchmarks: Ordination and Statistics

## Introduction

The [alpha and beta diversity
benchmarks](https://cmmr.github.io/ecodive/articles/benchmark.md)
compare ecodive’s diversity metrics with those of 15 other R packages.
This page does the same for the ordination methods and permutation
tests: each ecodive function is timed against the R functions it
replaces, on the same data, with ecodive running on 4 CPU threads.

Speed is worth nothing if the answers are wrong. So this page first sets
out what ecodive claims about its results, and the evidence behind each
claim:

1.  **Equivalence.** Every function gives the same result as an
    established implementation, to floating-point rounding wherever the
    method is deterministic.
2.  **Determinism.** The same call always returns the same result, bit
    for bit, on any number of threads, with or without
    [`set.seed()`](https://rdrr.io/r/base/Random.html).

Both claims are checked twice: by ecodive’s test suite, against stored
results from the reference packages, and again below on the benchmark
data itself.

  

## Equivalence

#### What “the same result” means

The methods fall into three groups, and “the same” means something
different for each.

- **Deterministic methods** (PCoA and its corrections, PCA, db-RDA,
  capscale, envfit, and the test statistics of PERMANOVA, PERMDISP,
  ANOSIM, MRPP and Mantel). There is one right answer, and ecodive
  returns it: scores, eigenvalues and statistics agree with the
  reference to about 1e-10 or better, relative to the largest value.
  Ordination axes may come out with the opposite sign, which carries no
  meaning, so signs are matched before comparing.

- **Permutation p-values.** A p-value depends on which random
  permutations were drawn. Given the *same* permutations as vegan
  (ecodive’s `permutations` argument accepts vegan’s permutation
  matrices), ecodive computes the same statistic for every one of them,
  and so the same p-value, exactly. With its own permutations, p-values
  differ from vegan’s only as much as two vegan runs with different
  seeds do.

- **Iterative methods with random starts or stochastic steps** (metric
  MDS, NMDS, t-SNE, UMAP). Each start converges to a different, equally
  valid layout, so coordinates can’t match any other package’s, or even
  vegan’s own from another seed. Here ecodive is checked two ways: each
  algorithm’s deterministic parts (affinities, neighbour graphs,
  gradients, steps from a shared starting layout) are compared with the
  reference implementation to rounding, and its final layouts are
  compared on quality measures (stress, KL divergence, neighbours kept),
  on which ecodive must do at least as well.

#### Evidence from the test suite

Each test below compares ecodive with results stored from the reference
package, so the tests need no reference package installed. The scripts
in
[data-raw/expected/](https://github.com/cmmr/ecodive/blob/main/data-raw/expected)
produce those results, and the tests in
[tests/tinytests/](https://github.com/cmmr/ecodive/blob/main/tests/tinytests)
run on GitHub CI on Linux, macOS and Windows. Tolerances are relative
differences.

| Function | Reference | What is compared | Tolerance | Test |
|:---|:---|:---|:---|:---|
| [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) | [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html), `ape::pcoa()`, `vegan::wcmdscale()` | Scores (up to sign), eigenvalues, proportions, total inertia | 1e-9 scores; 1e-10 the rest | [test-pcoa.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-pcoa.r) |
| `pcoa(correction =)` | `ape::pcoa()`, `vegan::wcmdscale()`, `cmdscale(add = TRUE)` | Lingoes and Cailliez constants, corrected eigenvalues, scores | 1e-10; 1e-8 scores | [test-pcoa_correct.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-pcoa_correct.r) |
| [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md) | scikit-bio 0.7.4 | Eigenvalues, coordinates, biplot scores | 1e-12; 1e-10 scores | [test-pcoa_biplot.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-pcoa_biplot.r) |
| [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) | scikit-learn 1.9.1 (SMACOF) | From the same start: iteration count (exact), stress at every 10th iteration, final layout. Quality: stress no worse than scikit-learn’s best | 1e-10 stress; 1e-8 layout | [test-mds.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-mds.r) |
| [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) | scikit-learn 1.9.1; `vegan::metaMDS()` | As [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md); and stress no worse than vegan’s | 1e-10 stress; 1e-8 layout | [test-nmds.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-nmds.r) |
| [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) | scikit-learn 1.9.1 | Neighbour graph (exact), geodesic distances, embedding, eigenvalues | 1e-12 geodesics; 1e-9 embedding | [test-isomap.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-isomap.r) |
| [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) | scikit-learn 1.9.1; `Rtsne::Rtsne()` | Affinities; exact and Barnes-Hut KL divergence and gradients; the first iterations’ steps. Quality: KL, trustworthiness and nearest-neighbour accuracy over 10 seeds | 1e-10 exact; 1e-6 where scikit-learn uses single precision | [test-tsne.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-tsne.r) |
| [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md) | `uwot::umap(batch = TRUE)` 0.2.5 | Nearest neighbours and fuzzy graph (exact), curve parameters, spectral start, optimised layout from the same start and seeds | 1e-14 graph; 0 (bit-identical) layout on Windows x86_64 | [test-umap.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-umap.r) |
| [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) | `vegan::dbrda()` | Site scores, eigenvalues, inertia, R², adjusted R², biplot scores, centroids | 1e-10 | [test-dbrda.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-dbrda.r) |
| [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md) | `vegan::capscale()` | As [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md), plus imaginary inertia | 1e-10 | [test-capscale.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-capscale.r) |
| [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) | [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) | Scores and loadings (up to sign), eigenvalues, proportions | 1e-10 | [test-pca.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-pca.r) |
| [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) | gemelli 0.0.13 (Python) | Every stage: rclr, starting subspaces, cost per iteration, completed matrix; final scores, loadings and distances | 1e-10 | [test-rpca.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-rpca.r) |
| [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) | `vegan::envfit()` | Arrows, r², centroids; p-values with vegan’s permutations | 1e-10; p-values identical | [test-envfit.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-envfit.r) |
| [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md) | `vegan::adonis2()` | The whole table (Df, SS, R², F); every permuted F; p-values | 1e-10; p-values identical | [test-permanova.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-permanova.r) |
| [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md) | `vegan::betadisper()`, `vegan::permutest()` | Distances to centres, dispersions, ANOVA table, every permuted F, p-values | 1e-10; p-values identical | [test-permdisp.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-permdisp.r) |
| [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md) | `vegan::anosim()` | R, every permuted R, mean ranks, p-values | 1e-14; p-values identical | [test-anosim.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-anosim.r) |
| [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) | `vegan::mrpp()` | δ, expected δ, A, every permuted δ, z, p-values | 1e-12; p-values identical | [test-mrpp.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-mrpp.r) |
| [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) | `vegan::mantel()`, `vegan::mantel.partial()` | r (Pearson, Spearman, Kendall; partial), every permuted r, p-values | 1e-12; p-values identical | [test-mantel.r](https://github.com/cmmr/ecodive/blob/main/tests/tinytests/test-mantel.r) |

#### Where ecodive deliberately differs

- **`permdisp(type = 'median')`.** vegan finds each group’s spatial
  median with a general-purpose optimiser that stops early, so its
  medians are only accurate to about 1e-5. ecodive solves for the exact
  median, and its results agree with an exact solution to 1e-10, and
  with vegan’s to about 1e-3. The default, `type = 'centroid'`, has no
  such difference.
- **[`umap()`](https://cmmr.github.io/ecodive/reference/umap.md)**
  always uses uwot’s batch optimiser (`batch = TRUE`), which can run in
  parallel and still give the same layout every time. uwot’s default
  optimiser can’t do both. Its `min_dist` default is 0.1, as in Python’s
  umap-learn (uwot’s is 0.01).
- **[`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md)**
  follows scikit-learn but computes in double precision throughout,
  where scikit-learn uses single precision for parts of the calculation.
  It starts from PCoA by default.
- **[`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) and
  [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md)** start
  from PCoA plus random starts and take relaxed (faster) SMACOF steps.
  They reproduce scikit-learn’s iterations exactly only when given its
  start and plain steps, which is how the tests check them.
- **[`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md)**
  reproduces `vegan::dbrda()`, and
  **[`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md)**
  reproduces `vegan::capscale()`. The two differ for non-Euclidean
  distances, in vegan as in ecodive.
- A few edge cases are handled differently: for example,
  [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) gives
  `NA` rather than 1 as the p-value of a constant variable, and
  [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) drops
  axes with non-positive eigenvalues instead of returning columns of
  zeros.

#### Evidence from the benchmark data

The same comparisons, repeated on the benchmark data set (1,006 samples;
see [Setup](#setup)). “Relative difference” is the largest absolute
difference divided by the largest absolute value in the reference
result. For the permutation tests, ecodive and vegan were both given the
same 999 permutations.

  

## Determinism

ecodive’s results depend only on the inputs and the `seed` argument.
They never depend on:

- **The number of threads.** Work is split into the same pieces whatever
  `cpus` is, and each piece is summed in a fixed order, so the
  floating-point rounding never changes.
- **R’s random number generator.** Methods with random parts take a
  `seed` argument (default `0`) and use their own generator. Calling the
  same function twice gives the same result, without
  [`set.seed()`](https://rdrr.io/r/base/Random.html), and ecodive never
  changes R’s random state.
- **Timing.** Threads never race to update shared values.

That matters for reproducible research: a p-value or a plot from a paper
can be regenerated exactly, on a laptop or on a 64-core server.

Many R packages are deterministic only after
[`set.seed()`](https://rdrr.io/r/base/Random.html), and some not even
then once they use more than one thread. Each function below was run
five times with its default arguments: twice in a row; twice more, each
after `set.seed(1)`; and once more after `set.seed(1)` with one thread
instead of four. A “yes” means the results were bit-for-bit identical.

DETERMINISM_NOTES

  

## Speed

#### Methodology

The `bench` R package timed each function and measured the memory it
allocated on R’s heap. All benchmarks ran on:

    CPU: 6-Core Intel i5-9600K @ 3.70GHz
    RAM: 64.0 GB
    OS:  Windows 11 Pro (64-bit, Build 26300)
    R:   4.6.1, with R's reference BLAS

ecodive used 4 threads (`cpus = 4`). Every other package that can use
several threads was given 4 as well: `Rtsne` (`num_threads`), `uwot`
(`n_threads`, `n_sgd_threads`), and vegan’s permutation tests and NMDS
(`parallel = 4`). On Windows, vegan’s `parallel` starts a cluster of R
processes on every call, which can make it slower than running serially,
so vegan’s tests were timed both ways.

A few things to keep in mind:

- Other packages’ eigendecompositions
  ([`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html),
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html), vegan’s constrained
  ordinations) call BLAS and LAPACK. With an optimised, multithreaded
  BLAS such as OpenBLAS or Intel MKL, they would run faster than with
  R’s reference BLAS used here. ecodive computes only the few axes that
  are returned, in its own threaded C, so its speed hardly depends on
  the BLAS.
- `bench` counts only memory allocated on R’s heap. That covers most of
  what these functions allocate, ecodive’s included, but not memory that
  C or Fortran code allocates directly.
- Where packages have different defaults, they were set to do comparable
  work: four starts for each NMDS and MDS that supports several, 999
  permutations for every test, and three axes for PCoA and PCA. The
  iterative methods then stop when their own convergence criteria are
  met.

#### Setup

The `gems` data set from the `rbiom` package has 1,006 stool samples
from children in four countries, with 767 OTUs, each sample’s diarrhea
status (case or control), age, and country. The benchmarks use its
Bray-Curtis distances, which every method takes as its input, except
PCA, which takes centered log-ratio transformed counts. The tests ask
whether country explains the distances, and the Mantel test correlates
them with differences in age. Distances and other inputs are computed
before the timed code.

Click to reveal the R code.

``` r

# Benchmarks and agreement checks for ecodive's ordination and statistical
# test functions, against the R packages they replace. Shown on the pkgdown
# site as vignettes/articles/benchmark_ord_stats.Rmd.
#
# Run from the package root:
#
#     Rscript vignettes/articles/benchmark_ord_stats/ord_stats.r
#
# Writes vignettes/articles/benchmark_ord_stats/ord_stats_results.rds (timings) and
# vignettes/articles/benchmark_ord_stats/ord_stats_agreement.rds (agreement with each reference), and
# vignettes/articles/figures/benchmarks_ord_stats.svg.


library(bench)
library(dplyr)
library(ggplot2)
library(ecodive)

CPUS <- 4L   # threads for ecodive, and for every reference that can use them


# GEMS dataset from rbiom: 1006 stool samples, 767 OTUs, and metadata.
gems   <- rbiom::gems
counts <- t(as.matrix(gems))                      # samples x OTUs
meta   <- as.data.frame(gems$metadata)
rownames(meta) <- meta$.sample
meta   <- meta[rownames(counts), c('diarrhea', 'age', 'country')]
n      <- nrow(counts)

# Inputs prepared outside the timed code.
dm     <- ecodive::bray(counts, cpus = CPUS)      # Bray-Curtis distances
age_dm <- dist(meta$age)                          # age differences
clr    <- log(counts + 1); clr <- clr - rowMeans(clr)  # CLR, pseudocount 1
ord    <- ecodive::pcoa(dm, n_axes = 2L, cpus = CPUS)  # an ordination for envfit

# One shared set of 999 permutations, so that ecodive and vegan test the very
# same permuted data sets and their p-values can be compared exactly.
set.seed(1)
perms <- permute::shuffleSet(n, 999)


# ---------------------------------------------------------------------------
# Agreement helpers
# ---------------------------------------------------------------------------

# Largest difference between two numeric objects, relative to the largest
# magnitude in the reference.
rel_diff <- function (x, ref) {
  x <- as.numeric(unlist(x)); ref <- as.numeric(unlist(ref))
  stopifnot(length(x) == length(ref))
  max(abs(x - ref)) / max(abs(ref))
}

# The same, for coordinates whose axes may come out with either sign.
rel_diff_signed <- function (x, ref) {
  x <- as.matrix(x); ref <- as.matrix(ref)
  for (j in seq_len(ncol(x)))
    if (sum(x[, j] * ref[, j]) < 0) x[, j] <- -x[, j]
  rel_diff(x, ref)
}

# Fraction of each sample's 15 nearest neighbours (in dm) that are also among
# its 15 nearest on the plot: how well an embedding keeps local structure.
dm_mtx <- as.matrix(dm)
nn_dm  <- apply(dm_mtx, 1L, function (d) order(d)[2:16])
nn_kept <- function (xy) {
  nn_xy <- apply(as.matrix(dist(xy)), 1L, function (d) order(d)[2:16])
  mean(vapply(seq_len(n), function (i) length(intersect(nn_dm[, i], nn_xy[, i])), 0)) / 15
}

# Correlation of two configurations after the best rotation, reflection and
# scaling (Procrustes; 1 = the same shape).
procrustes_r <- function (x, y) vegan::protest(x, y, permutations = 0)$t0

# Whether two p-values are exactly equal.
same_p <- function (a, b)
  sprintf('%s (%s)', if (identical(as.numeric(a), as.numeric(b))) 'identical' else 'DIFFERENT',
          paste(unique(c(a, b)), collapse = ' vs ') )

agreement <- list()
agree <- function (method, reference, quantity, value)
  agreement[[length(agreement) + 1L]] <<- data.frame(
    method = method, reference = reference, quantity = quantity,
    value = if (is.numeric(value)) signif(value, 3) else value )


# ---------------------------------------------------------------------------
# Timings
# ---------------------------------------------------------------------------

timed <- function (method, iterations, ...) {
  if (nzchar(Sys.getenv("QUICK"))) return (NULL)  # agreement checks only
  res <- bench::mark(..., iterations = iterations, check = FALSE, filter_gc = FALSE)
  res$method <- method
  res$package <- as.character(res$expression)
  print(res[, c('package', 'min', 'median', 'mem_alloc')] %>% arrange(median))
  res
}

results <- list()


# ---------------------------------------------------------------------------
# Principal Coordinates Analysis
# ---------------------------------------------------------------------------

results$pcoa <- timed(
  'PCoA', iterations = 10,
  'ecodive' = ecodive::pcoa(dm, n_axes = 3, cpus = CPUS),
  'stats'   = stats::cmdscale(dm, k = 3, eig = TRUE),
  'ape'     = ape::pcoa(dm),
  'vegan'   = vegan::wcmdscale(dm, k = 3, eig = TRUE),
  'labdsv'  = labdsv::pco(dm, k = 3),
  'ecodist' = ecodist::pco(dm) )

x  <- ecodive::pcoa(dm, n_axes = 3, cpus = CPUS)
cm <- stats::cmdscale(dm, k = 3, eig = TRUE)
ap <- ape::pcoa(dm)
agree('PCoA', 'stats::cmdscale', 'scores (3 axes)',     rel_diff_signed(x, cm$points))
agree('PCoA', 'stats::cmdscale', 'eigenvalues',         rel_diff(attr(x, 'eig'), cm$eig[1:3]))
agree('PCoA', 'ape::pcoa',       'scores (3 axes)',     rel_diff_signed(x, ap$vectors[, 1:3]))
agree('PCoA', 'ape::pcoa',       'variance explained',  rel_diff(attr(x, 'prop'), ap$values$Relative_eig[1:3]))


results$pcoa_cailliez <- timed(
  'PCoA (Cailliez)', iterations = 3,
  'ecodive' = ecodive::pcoa(dm, n_axes = 3, correction = 'cailliez', cpus = CPUS),
  'stats'   = stats::cmdscale(dm, k = 3, eig = TRUE, add = TRUE),
  'ape'     = ape::pcoa(dm, correction = 'cailliez'),
  'vegan'   = vegan::wcmdscale(dm, k = 3, eig = TRUE, add = 'cailliez') )

x  <- ecodive::pcoa(dm, n_axes = 3, correction = 'cailliez', cpus = CPUS)
vg <- vegan::wcmdscale(dm, k = 3, eig = TRUE, add = 'cailliez')
agree('PCoA (Cailliez)', 'vegan::wcmdscale', 'constant',        rel_diff(attr(x, 'constant'), vg$ac))
agree('PCoA (Cailliez)', 'vegan::wcmdscale', 'scores (3 axes)', rel_diff_signed(x, vg$points))
agree('PCoA (Cailliez)', 'vegan::wcmdscale', 'eigenvalues',     rel_diff(attr(x, 'eig'), vg$eig[1:3]))


# ---------------------------------------------------------------------------
# Permutation tests: same 999 permutations for timing (by count) and an exact
# comparison (by the shared `perms` matrix).
# ---------------------------------------------------------------------------

results$permanova <- timed(
  'PERMANOVA', iterations = 5,
  'ecodive'          = ecodive::permanova(dm, 'country', meta, permutations = 999, cpus = CPUS),
  'vegan'            = vegan::adonis2(dm ~ country, data = meta, permutations = 999),
  'vegan (parallel)' = vegan::adonis2(dm ~ country, data = meta, permutations = 999, parallel = CPUS) )

x  <- ecodive::permanova(dm, 'country', meta, permutations = perms, cpus = CPUS)
vg <- vegan::adonis2(dm ~ country, data = meta, permutations = perms)
agree('PERMANOVA', 'vegan::adonis2', 'pseudo-F',          rel_diff(x$statistic, vg$F[1]))
agree('PERMANOVA', 'vegan::adonis2', 'R-squared',         rel_diff(x$r2, vg$R2[1]))
agree('PERMANOVA', 'vegan::adonis2', 'permuted F (999)',  rel_diff(attr(x, 'f_perm'), attr(vg, 'F.perm')))
agree('PERMANOVA', 'vegan::adonis2', 'p-value',           same_p(x$p_value, vg$`Pr(>F)`[1]))


results$permdisp <- timed(
  'PERMDISP', iterations = 5,
  'ecodive'          = ecodive::permdisp(dm, 'country', meta, type = 'centroid', permutations = 999, cpus = CPUS),
  'vegan'            = vegan::permutest(vegan::betadisper(dm, meta$country, type = 'centroid'), permutations = 999),
  'vegan (parallel)' = vegan::permutest(vegan::betadisper(dm, meta$country, type = 'centroid'), permutations = 999, parallel = CPUS) )

x  <- ecodive::permdisp(dm, 'country', meta, type = 'centroid', permutations = perms, cpus = CPUS)
bd <- vegan::betadisper(dm, meta$country, type = 'centroid')
vg <- vegan::permutest(bd, permutations = perms)
agree('PERMDISP', 'vegan::betadisper', 'distances to centroid', rel_diff(x$distances, bd$distances))
agree('PERMDISP', 'vegan::permutest',  'F',                      rel_diff(x$statistic, vg$tab$F[1]))
agree('PERMDISP', 'vegan::permutest',  'permuted F (999)',       rel_diff(attr(x, 'f_perm'), vg$perm))
agree('PERMDISP', 'vegan::permutest',  'p-value',                same_p(x$p_value, vg$tab$`Pr(>F)`[1]))


results$anosim <- timed(
  'ANOSIM', iterations = 3,
  'ecodive'          = ecodive::anosim(dm, 'country', meta, permutations = 999, cpus = CPUS),
  'vegan'            = vegan::anosim(dm, meta$country, permutations = 999),
  'vegan (parallel)' = vegan::anosim(dm, meta$country, permutations = 999, parallel = CPUS) )

x  <- ecodive::anosim(dm, 'country', meta, permutations = perms, cpus = CPUS)
vg <- vegan::anosim(dm, meta$country, permutations = perms)
agree('ANOSIM', 'vegan::anosim', 'R',                rel_diff(x$statistic, vg$statistic))
agree('ANOSIM', 'vegan::anosim', 'permuted R (999)', rel_diff(attr(x, 'r_perm'), vg$perm))
agree('ANOSIM', 'vegan::anosim', 'p-value',          same_p(x$p_value, vg$signif))


results$mrpp <- timed(
  'MRPP', iterations = 3,
  'ecodive'          = ecodive::mrpp(dm, 'country', meta, permutations = 999, cpus = CPUS),
  'vegan'            = vegan::mrpp(dm, meta$country, permutations = 999),
  'vegan (parallel)' = vegan::mrpp(dm, meta$country, permutations = 999, parallel = CPUS) )

x  <- ecodive::mrpp(dm, 'country', meta, permutations = perms, cpus = CPUS)
vg <- vegan::mrpp(dm, meta$country, permutations = perms)
agree('MRPP', 'vegan::mrpp', 'delta',                rel_diff(x$statistic, vg$delta))
agree('MRPP', 'vegan::mrpp', 'A',                    rel_diff(x$A, vg$A))
agree('MRPP', 'vegan::mrpp', 'permuted delta (999)', rel_diff(attr(x, 'delta_perm'), vg$boot.deltas))
agree('MRPP', 'vegan::mrpp', 'p-value',              same_p(x$p_value, vg$Pvalue))


results$mantel <- timed(
  'Mantel', iterations = 5,
  'ecodive'          = ecodive::mantel(dm, age_dm, permutations = 999, cpus = CPUS),
  'vegan'            = vegan::mantel(dm, age_dm, permutations = 999),
  'vegan (parallel)' = vegan::mantel(dm, age_dm, permutations = 999, parallel = CPUS),
  'ade4'             = ade4::mantel.rtest(dm, age_dm, nrepet = 999),
  'ecodist'          = ecodist::mantel(dm ~ age_dm, nperm = 999, nboot = 0) )

x  <- ecodive::mantel(dm, age_dm, permutations = perms, cpus = CPUS)
vg <- vegan::mantel(dm, age_dm, permutations = perms)
agree('Mantel', 'vegan::mantel', 'r',                rel_diff(x$statistic, vg$statistic))
agree('Mantel', 'vegan::mantel', 'permuted r (999)', rel_diff(attr(x, 'r_perm'), vg$perm))
agree('Mantel', 'vegan::mantel', 'p-value',          same_p(x$p_value, vg$signif))


# ---------------------------------------------------------------------------
# Constrained ordination: country + age
# ---------------------------------------------------------------------------

results$dbrda <- timed(
  'db-RDA', iterations = 10,
  'ecodive' = ecodive::dbrda(dm, data = meta, formula = ~ country + age, cpus = CPUS),
  'vegan'   = vegan::dbrda(dm ~ country + age, data = meta) )

x  <- ecodive::dbrda(dm, data = meta, formula = ~ country + age, cpus = CPUS)
vg <- vegan::dbrda(dm ~ country + age, data = meta)
agree('db-RDA', 'vegan::dbrda', 'site scores (3 axes)',     rel_diff_signed(x, vegan::scores(vg, choices = 1:3, display = 'sites')))
agree('db-RDA', 'vegan::dbrda', 'constrained eigenvalues', rel_diff(attr(x, 'constrained_eig'), vg$CCA$eig))
agree('db-RDA', 'vegan::dbrda', 'R-squared',               rel_diff(attr(x, 'r_squared'), vegan::RsquareAdj(vg)$r.squared))


results$capscale <- timed(
  'capscale', iterations = 10,
  'ecodive' = ecodive::capscale(dm, data = meta, formula = ~ country + age, cpus = CPUS),
  'vegan'   = vegan::capscale(dm ~ country + age, data = meta) )

x  <- ecodive::capscale(dm, data = meta, formula = ~ country + age, cpus = CPUS)
vg <- vegan::capscale(dm ~ country + age, data = meta)
agree('capscale', 'vegan::capscale', 'site scores (3 axes)',     rel_diff_signed(x, vegan::scores(vg, choices = 1:3, display = 'sites')))
agree('capscale', 'vegan::capscale', 'constrained eigenvalues', rel_diff(attr(x, 'constrained_eig'), vg$CCA$eig))
agree('capscale', 'vegan::capscale', 'R-squared',               rel_diff(attr(x, 'r_squared'), vegan::RsquareAdj(vg)$r.squared))


# ---------------------------------------------------------------------------
# PCA of CLR-transformed counts. ecodive's time includes the transformation;
# the others are given `clr`, computed beforehand.
# ---------------------------------------------------------------------------

results$pca <- timed(
  'PCA', iterations = 10,
  'ecodive' = ecodive::pca(counts, norm = 'clr', pseudocount = 1, n_axes = 3, cpus = CPUS),
  'stats'   = stats::prcomp(clr, rank. = 3),
  'irlba'   = irlba::prcomp_irlba(clr, n = 3),
  'vegan'   = vegan::rda(clr),
  'ade4'    = ade4::dudi.pca(clr, scale = FALSE, scannf = FALSE, nf = 3) )

x  <- ecodive::pca(counts, norm = 'clr', pseudocount = 1, n_axes = 3, cpus = CPUS)
pr <- stats::prcomp(clr, rank. = 3)
agree('PCA', 'stats::prcomp', 'scores (3 axes)',   rel_diff_signed(x, pr$x))
agree('PCA', 'stats::prcomp', 'loadings (3 axes)', rel_diff_signed(attr(x, 'loadings'), pr$rotation))
agree('PCA', 'stats::prcomp', 'eigenvalues',       rel_diff(attr(x, 'eig'), pr$sdev[1:3]^2))


# ---------------------------------------------------------------------------
# envfit: fit diarrhea, age and country onto a 2-axis PCoA.
# ---------------------------------------------------------------------------

results$envfit <- timed(
  'envfit', iterations = 5,
  'ecodive' = ecodive::envfit(ord, meta, permutations = 999, cpus = CPUS),
  'vegan'   = vegan::envfit(ord, meta, permutations = 999) )

x  <- ecodive::envfit(ord, meta, permutations = perms, cpus = CPUS)
vg <- vegan::envfit(ord, meta, permutations = perms)
vg_r2 <- c(vg$vectors$r, vg$factors$r)
vg_p  <- c(vg$vectors$pvals, vg$factors$pvals)
x1    <- x[!duplicated(x$variable), ]
x1    <- x1[match(names(vg_r2), x1$variable), ]
agree('envfit', 'vegan::envfit', 'r-squared (3 variables)', rel_diff(x1$r2, vg_r2))
agree('envfit', 'vegan::envfit', 'arrow (age)', rel_diff(unlist(x['age', colnames(ord)]), vg$vectors$arrows['age', ]))
agree('envfit', 'vegan::envfit', 'centroids (6 levels)',
      rel_diff(as.matrix(x[!is.na(x$level), colnames(ord)]), vg$factors$centroids[rownames(x)[!is.na(x$level)], ]))
agree('envfit', 'vegan::envfit', 'p-values', same_p(x1$p_value, unname(vg_p)))


# ---------------------------------------------------------------------------
# Iterative ordinations. Each package runs four starts where it can (ecodive's
# default), so the work is comparable. Their outputs are compared by quality
# (stress, KL divergence, neighbours kept) and by shape (Procrustes), since
# different random starts give different, equally valid, coordinates.
# ---------------------------------------------------------------------------

results$nmds <- timed(
  'NMDS', iterations = 3,
  'ecodive' = ecodive::nmds(dm, n_axes = 2, n_starts = 4, seed = 1, cpus = CPUS),
  'vegan'   = vegan::metaMDS(dm, k = 2, try = 3, trymax = 3, parallel = CPUS, trace = 0),
  'MASS'    = MASS::isoMDS(dm, k = 2, trace = FALSE),
  'smacof'  = smacof::smacofSym(dm, ndim = 2, type = 'ordinal') )

x  <- ecodive::nmds(dm, n_axes = 2, n_starts = 4, seed = 1, cpus = CPUS)
vg <- vegan::metaMDS(dm, k = 2, try = 3, trymax = 3, parallel = CPUS, trace = 0)
ms <- MASS::isoMDS(dm, k = 2, trace = FALSE)
agree('NMDS', 'vegan::metaMDS', 'stress (ecodive vs vegan)', sprintf('%.4f vs %.4f', attr(x, 'stress'), vg$stress))
agree('NMDS', 'MASS::isoMDS',   'stress (ecodive vs MASS)',  sprintf('%.4f vs %.4f', attr(x, 'stress'), ms$stress / 100))
agree('NMDS', 'vegan::metaMDS', 'Procrustes correlation',    procrustes_r(vg$points, x))


results$mds <- timed(
  'Metric MDS', iterations = 3,
  'ecodive' = ecodive::mds(dm, n_axes = 2, n_starts = 4, seed = 1, cpus = CPUS),
  'smacof'  = smacof::smacofSym(dm, ndim = 2, type = 'ratio') )

x  <- ecodive::mds(dm, n_axes = 2, n_starts = 4, seed = 1, cpus = CPUS)
sm <- smacof::smacofSym(dm, ndim = 2, type = 'ratio')
agree('Metric MDS', 'smacof::smacofSym', 'stress-1 (ecodive vs smacof)', sprintf('%.4f vs %.4f', attr(x, 'stress'), sm$stress))
agree('Metric MDS', 'smacof::smacofSym', 'Procrustes correlation',       procrustes_r(sm$conf, x))


results$isomap <- timed(
  'Isomap', iterations = 5,
  'ecodive' = ecodive::isomap(dm, n_axes = 2, n_neighbors = 15, cpus = CPUS),
  'vegan'   = vegan::isomap(dm, ndim = 2, n_axes = 15) )

x  <- ecodive::isomap(dm, n_axes = 2, n_neighbors = 15, cpus = CPUS)
vg <- vegan::isomap(dm, ndim = 2, n_axes = 15)
agree('Isomap', 'vegan::isomap', 'scores (2 axes)', rel_diff_signed(x, vg$points))


results$tsne <- timed(
  't-SNE', iterations = 5,
  'ecodive' = ecodive::tsne(dm, n_axes = 2, seed = 1, cpus = CPUS),
  'Rtsne'   = Rtsne::Rtsne(dm, dims = 2, is_distance = TRUE, pca = FALSE, num_threads = CPUS) )

x  <- ecodive::tsne(dm, n_axes = 2, seed = 1, cpus = CPUS)
set.seed(1)
rt <- Rtsne::Rtsne(dm, dims = 2, is_distance = TRUE, pca = FALSE, num_threads = CPUS)
agree('t-SNE', 'Rtsne::Rtsne', 'KL divergence (ecodive vs Rtsne)', sprintf('%.3f vs %.3f', attr(x, 'kl_divergence'), tail(rt$itercosts, 1)))
agree('t-SNE', 'Rtsne::Rtsne', '15-NN kept (ecodive vs Rtsne)',    sprintf('%.3f vs %.3f', nn_kept(x), nn_kept(rt$Y)))


results$umap <- timed(
  'UMAP', iterations = 5,
  'ecodive' = ecodive::umap(dm, n_axes = 2, seed = 1, cpus = CPUS),
  'uwot'    = uwot::umap(dm, n_components = 2, min_dist = 0.1, batch = TRUE, n_threads = CPUS, n_sgd_threads = CPUS, seed = 1) )

x  <- ecodive::umap(dm, n_axes = 2, seed = 1, cpus = CPUS)
uw <- uwot::umap(dm, n_components = 2, min_dist = 0.1, batch = TRUE, n_threads = CPUS, n_sgd_threads = CPUS, seed = 1)
agree('UMAP', 'uwot::umap', '15-NN kept (ecodive vs uwot)', sprintf('%.3f vs %.3f', nn_kept(x), nn_kept(uw)))


# ---------------------------------------------------------------------------
# Determinism. Each function is run five times with default arguments:
#   1, 2  twice in a row, with no set.seed() call;
#   3, 4  twice after set.seed(1);
#   5     after set.seed(1) again, with 1 thread instead of 4.
# A runner takes the number of threads and returns the numbers to compare.
# Permutation tests use 99 permutations here, to save time.
# ---------------------------------------------------------------------------

determinism <- list()
det_check <- function (method, package, runner, threaded = TRUE) {
  a  <- runner(CPUS);  b  <- runner(CPUS)
  set.seed(1); s1 <- runner(CPUS)
  set.seed(1); s2 <- runner(CPUS)
  set.seed(1); t1 <- if (threaded) runner(1L) else NULL
  determinism[[length(determinism) + 1L]] <<- data.frame(
    method = method, package = package,
    rerun            = identical(a, b),
    rerun_same_seed  = identical(s1, s2),
    one_vs_4_threads = if (threaded) identical(t1, s1) else NA )
  print(tail(determinism, 1)[[1]])
}
num <- function (x) unname(as.numeric(unlist(x)))

# ecodive: every function, default arguments (seed = 0), cpus = 1 or 4.
det_check('PCoA',            'ecodive', function (th) ecodive::pcoa(dm, cpus = th))
det_check('PCoA (Cailliez)', 'ecodive', function (th) ecodive::pcoa(dm, correction = 'cailliez', cpus = th))
det_check('Metric MDS',      'ecodive', function (th) ecodive::mds(dm, cpus = th))
det_check('NMDS',            'ecodive', function (th) ecodive::nmds(dm, cpus = th))
det_check('Isomap',          'ecodive', function (th) ecodive::isomap(dm, n_neighbors = 15, cpus = th))
det_check('t-SNE',           'ecodive', function (th) ecodive::tsne(dm, cpus = th))
det_check('UMAP',            'ecodive', function (th) ecodive::umap(dm, cpus = th))
det_check('db-RDA',          'ecodive', function (th) ecodive::dbrda(dm, data = meta, formula = ~ country + age, cpus = th))
det_check('capscale',        'ecodive', function (th) ecodive::capscale(dm, data = meta, formula = ~ country + age, cpus = th))
det_check('PCA',             'ecodive', function (th) ecodive::pca(counts, norm = 'clr', pseudocount = 1, cpus = th))
det_check('envfit',          'ecodive', function (th) ecodive::envfit(ord, meta, cpus = th))
det_check('PERMANOVA',       'ecodive', function (th) ecodive::permanova(dm, 'country', meta, cpus = th))
det_check('PERMDISP',        'ecodive', function (th) ecodive::permdisp(dm, 'country', meta, cpus = th))
det_check('ANOSIM',          'ecodive', function (th) ecodive::anosim(dm, 'country', meta, cpus = th))
det_check('MRPP',            'ecodive', function (th) ecodive::mrpp(dm, 'country', meta, cpus = th))
det_check('Mantel',          'ecodive', function (th) ecodive::mantel(dm, age_dm, cpus = th))

# Other packages: the functions with random or threaded parts.
det_check('NMDS',      'vegan',         function (th) num(vegan::metaMDS(dm, k = 2, try = 3, trymax = 3, parallel = th, trace = 0)$points))
det_check('NMDS',      'MASS',          function (th) num(MASS::isoMDS(dm, k = 2, trace = FALSE)$points), threaded = FALSE)
det_check('t-SNE',     'Rtsne',         function (th) num(Rtsne::Rtsne(dm, is_distance = TRUE, pca = FALSE, num_threads = th)$Y))
det_check('UMAP',      'uwot',          function (th) num(uwot::umap(dm, min_dist = 0.1, n_threads = th, n_sgd_threads = th)))
det_check('UMAP',      'uwot (batch)',  function (th) num(uwot::umap(dm, min_dist = 0.1, n_threads = th, n_sgd_threads = th, batch = TRUE)))
det_check('PCA',       'irlba',         function (th) num(irlba::prcomp_irlba(clr, n = 3)$x), threaded = FALSE)
det_check('PERMANOVA', 'vegan',         function (th) num(attr(vegan::adonis2(dm ~ country, data = meta, permutations = 99, parallel = th), 'F.perm')))
det_check('PERMDISP',  'vegan',         function (th) num(vegan::permutest(vegan::betadisper(dm, meta$country, type = 'centroid'), permutations = 99, parallel = th)$perm))
det_check('ANOSIM',    'vegan',         function (th) num(vegan::anosim(dm, meta$country, permutations = 99, parallel = th)$perm))
det_check('MRPP',      'vegan',         function (th) num(vegan::mrpp(dm, meta$country, permutations = 99, parallel = th)$boot.deltas))
det_check('Mantel',    'vegan',         function (th) num(vegan::mantel(dm, age_dm, permutations = 99, parallel = th)$perm))
det_check('Mantel',    'ade4',          function (th) num(ade4::mantel.rtest(dm, age_dm, nrepet = 99)$sim), threaded = FALSE)
det_check('Mantel',    'ecodist',       function (th) num(ecodist::mantel(dm ~ age_dm, nperm = 99, nboot = 0)), threaded = FALSE)

determinism <- do.call(rbind, determinism)
print(determinism)


# ---------------------------------------------------------------------------
# Save and plot.
# ---------------------------------------------------------------------------

agreement <- do.call(rbind, agreement)
print(agreement, right = FALSE)

if (length(results)) {
timings <- bind_rows(results) %>%
  select(method, package, median, mem_alloc) %>%
  group_by(method) %>%
  mutate(speedup = as.numeric(median) / as.numeric(median[package == 'ecodive'])) %>%
  ungroup()
print(timings, n = Inf)

# Plain numbers (seconds, bytes), so the article can read them without bench.
saveRDS(
  with(timings, data.frame(
    method, package, seconds = as.numeric(median), bytes = as.numeric(mem_alloc), speedup )),
  'vignettes/articles/benchmark_ord_stats/ord_stats_results.rds', version = 2 )
saveRDS(agreement,     'vignettes/articles/benchmark_ord_stats/ord_stats_agreement.rds')
saveRDS(determinism,   'vignettes/articles/benchmark_ord_stats/ord_stats_determinism.rds')

writeLines(capture.output(sessionInfo()), 'vignettes/articles/benchmark_ord_stats/ord_stats_session.txt')
}


# ---------------------------------------------------------------------------
# Figure: median time of each package, per method.
# ---------------------------------------------------------------------------

if (length(results)) {

  order_methods <- c(
    'PCoA', 'PCoA (Cailliez)', 'Metric MDS', 'NMDS', 'Isomap', 't-SNE', 'UMAP',
    'db-RDA', 'capscale', 'PCA', 'envfit',
    'PERMANOVA', 'PERMDISP', 'ANOSIM', 'MRPP', 'Mantel' )

  fig_data <- timings %>%
    mutate(
      method  = factor(method, order_methods),
      seconds = as.numeric(median),
      is_eco  = package == 'ecodive',
      label   = ifelse(is_eco, 'ecodive', sprintf('%s  (%sx)', package, formatC(speedup, format = 'fg', digits = 2, big.mark = ','))) ) %>%
    group_by(method) %>%
    mutate(package = reorder(paste(method, package), -seconds)) %>%
    ungroup()

  fig <- ggplot(fig_data, aes(x = seconds, y = package)) +
    geom_segment(aes(x = 1e-3, xend = seconds, yend = package, colour = is_eco), linewidth = 0.6) +
    geom_point(aes(colour = is_eco), size = 2) +
    geom_text(aes(label = label, fontface = ifelse(is_eco, 2, 1)), hjust = -0.15, size = 3) +
    facet_wrap('method', ncol = 4, scales = 'free_y') +
    scale_y_discrete(labels = NULL) +
    scale_x_log10(
      labels = function (x) bench::as_bench_time(x),
      expand = expansion(mult = c(0.02, 0.9)) ) +
    scale_colour_manual(values = c(`TRUE` = '#D55E00', `FALSE` = 'grey45'), guide = 'none') +
    labs(
      x = sprintf('Median time (log scale); ecodive uses %d threads. In brackets: times slower than ecodive.', CPUS),
      y = NULL ) +
    theme_bw(base_size = 11) +
    theme(
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text         = element_text(face = 'bold', hjust = 0) )

  ggsave(
    'vignettes/articles/figures/benchmarks_ord_stats.svg', fig,
    width = 12, height = 9, device = svglite::svglite )
}
```

Click to reveal the R session information.

#### Results

SPEED_NOTES
