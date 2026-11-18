# Plan: benchmark ecodive's ordination and statistical tests

This file is the brief for a cloud Claude Code session. It runs the benchmarks,
finishes the pkgdown article, and commits the results. Read it all before
starting. Everything here was worked out in a local session on 2026-11-18,
which did not run the final benchmarks because the user's desktop was busy.


## 1. Goal

Publish a pkgdown article, `vignettes/articles/benchmark_ord_stats.Rmd`, that
does for ecodive's ordination and statistical-test functions what
`vignettes/articles/benchmark.Rmd` does for alpha and beta diversity. It must:

1. **Prove the results are correct.** "Correct" means identical, or very
   close, to the established implementations. Readers won't trust the speed
   numbers unless they trust the outputs. The article states the claims and
   shows two kinds of evidence:
   - the test suite: its reference fixtures and tolerances (the table is
     already written in the draft article);
   - the benchmark data: agreement numbers that the benchmark script computes.
2. **Prove the results are deterministic.** The same call gives bit-identical
   results on every run, with or without `set.seed()`, for any number of
   threads. The article lists, as measured, the other packages whose results
   are *not* deterministic.
3. **Compare speed.** ecodive runs with `cpus = 4`. It is timed against every
   R implementation of the same method that is available and works, each
   given 4 threads where it can use them.

Include every competitor. The user doesn't want to leave out any package
that might be as fast as ecodive, or faster.


## 2. Files that already exist (drafts, committed with this plan)

| File | State |
|:-----|:------|
| `vignettes/articles/benchmark_ord_stats/ord_stats.r` | Benchmark, agreement and determinism script. Its agreement and determinism parts have been tried (see §7). It **does not yet include** the packages added in §5.2, which you must add. |
| `vignettes/articles/benchmark_ord_stats.Rmd` | Draft article. It reads `vignettes/articles/benchmark_ord_stats/ord_stats_*.rds` and `ord_stats_session.txt` and embeds `vignettes/articles/benchmark_ord_stats/ord_stats.r` with `readLines()`. Two placeholders remain, `DETERMINISM_NOTES` and `SPEED_NOTES`, and the "Methodology" section still describes the user's Windows desktop. Update both once the results are in. |
| `vignettes/articles/benchmark.Rmd` | The existing alpha/beta benchmark page. Match its tone. |
| `data-raw/expected/*.r`, `tests/tinytests/test-*.r` | Fixture generators and tests. The draft article's test-suite table cites them. |

The script, this plan, the results and the article all live under
`vignettes/articles/`, which is in `.Rbuildignore` (so none of it goes into
the CRAN package) but is tracked by git (so pkgdown's CI can read the
results). Don't use the top-level `benchmark/` folder: it is git-ignored.
The script writes its outputs to `vignettes/articles/benchmark_ord_stats/`,
and is run from the package root.


## 3. Ground rules

- **Never add packages to ecodive's `DESCRIPTION`**, not even to Suggests.
  The benchmark packages are installed only in the cloud environment.
- **Build ecodive at `-O2`.** Install it from the repository with
  `R CMD INSTALL --preclean .`. Never time a `devtools::load_all()` build,
  which compiles at `-O0`. Delete any `src/*.o` and `src/*.so` files first.
- **Every text file must use LF line endings.**
- **Report results as measured.** If ecodive is slower at something, say
  so. In the local probe, `nmds()` with 4 starts took 8.9 s, against 2.1 s
  for one start of `vegan::monoMDS()`. If a check fails, such as a
  `DIFFERENT` p-value, stop and investigate before publishing. Don't hide it.
- **Commits.** The user's rule is: commit only to the default branch
  (`main`), and never push without asking. A cloud session may have to work
  on its own branch instead. If so, commit there, don't push to `main`, and
  say in your final report which branch holds the work, so that the user can
  merge it. Attribution lines go at the end of each commit message, as your
  system prompt says.


## 4. Environment

### 4.1 Machine

- Linux (Ubuntu 24.04 or similar), with **at least 8 vCPUs**, dedicated
  rather than burstable or shared-core, and at least 16 GB of RAM.
- Run the whole script pinned to 4 cores, so that ecodive and every
  reference package (including vegan's forked workers and BLAS threads)
  share the same 4 CPUs:

  ```sh
  taskset -c 0-3 Rscript vignettes/articles/benchmark_ord_stats/ord_stats.r > vignettes/articles/benchmark_ord_stats/ord_stats.log 2>&1
  ```

  The other cores absorb system noise. Record `lscpu`, `nproc`, `free -g`,
  the kernel version, and the cloud instance type in
  `vignettes/articles/benchmark_ord_stats/ord_stats_machine.txt`. The article's "Methodology" section
  must describe this machine, not the Windows desktop.
- Before starting, check that nothing else is using the CPU (`uptime`, `top`).

### 4.2 BLAS

`cmdscale()`, `prcomp()`, `ape::pcoa()`, vegan's eigendecompositions,
`irlba`, `rsvd` and others spend most of their time in BLAS and LAPACK.
Compare against the best case for them:

- **Main run: OpenBLAS, limited to 4 threads.** Install
  `libopenblas-openmp-dev` or `libopenblas0-pthread`, select it with
  `update-alternatives`, and set `OPENBLAS_NUM_THREADS=4` and
  `OMP_NUM_THREADS=4`. Confirm the choice with `sessionInfo()` (the BLAS and
  LAPACK lines) and with `ecodive::blas_threads()`.
- **Optional second run: R's reference BLAS.** Run only the BLAS-heavy
  methods (PCoA, PCoA with Cailliez, PCA, db-RDA, capscale) and add them to
  the article as a short "with reference BLAS" table. This shows how much of
  each competitor's speed comes from BLAS.

### 4.3 R and packages

R 4.6.x. Install with pak, which installs system requirements on Ubuntu:

```r
install.packages('pak')
pak::pak(c(
  'bench', 'dplyr', 'ggplot2', 'svglite', 'rbiom', 'permute', 'knitr', 'rmarkdown', 'pkgdown',
  'vegan', 'ape', 'labdsv', 'ecodist', 'ade4', 'MASS', 'smacof', 'fmds',
  'Rtsne', 'tsne', 'uwot', 'umap', 'manifoldsR',
  'irlba', 'RSpectra', 'bigstatsr', 'mlpack', 'pcaone', 'Rfast2', 'rsvd', 'FactoMineR',
  'LDM', 'GUniFrac', 'ncf',
  'bioc::scrapper', 'bioc::BiocNeighbors' ))
```

- `manifoldsR` is built from Rust source. Install a Rust toolchain first
  (`rustup`, the stable toolchain). If it still fails, leave it out and say
  so in the article.
- `scrapper` (Bioconductor) needs CMake.
- Then install ecodive from the checkout (see §3) and record
  `packageVersion()` for every package in the session file.


## 5. What to benchmark

### 5.1 Data and settings (already in the script)

- `rbiom::gems`: 1,006 samples × 767 OTUs, with metadata `diarrhea`, `age`
  and `country`.
- `dm` is `ecodive::bray(counts)`. There are no zero distances. `age_dm` is
  `dist(meta$age)`.
- `clr` is `log(counts + 1)`, centred by row. Precompute it for the other
  packages; ecodive's `pca()` does the transformation inside its timing,
  which is conservative.
- The tests use groups = `country`. The Mantel test uses `dm` against
  `age_dm`. 999 permutations throughout.
- One shared permutation matrix, `perms <- permute::shuffleSet(n, 999)` with
  `set.seed(1)`, is passed to both ecodive and vegan in the agreement checks,
  so their p-values can be compared exactly.
- Precompute any input conversion outside the timed code, such as
  `as.matrix(dm)`, kNN lists or a `bigstatsr` FBM, as the existing page
  does. Say so in the article whenever it favours the competitor.

### 5.2 Competitors per method

✓ = already in `ord_stats.r`. **+** = add it. The calls marked + were tried
in the local probe and work as written, unless a note says otherwise.

**PCoA (k = 3)**
- ✓ `stats::cmdscale`, `ape::pcoa`, `vegan::wcmdscale`, `labdsv::pco`, `ecodist::pco`
- + `fmds::pcoa(dmm, p = 3)`, where `dmm <- as.matrix(dm)`. C, single-threaded. Took 8.3 s locally.
- + `smacof::torgerson(dm, p = 3)`. It agrees with ecodive to 1.9e-14, so add that to the agreement table.
- Skip `ade4::dudi.pco`: it needs Euclidean input, and Bray-Curtis isn't Euclidean. Mention this in the article.

**PCoA with the Cailliez correction**
- ✓ `cmdscale(add = TRUE)`, `ape::pcoa(correction = 'cailliez')`, `vegan::wcmdscale(add = 'cailliez')`
- `fmds::pcoa(ac = TRUE)` exists, but it's unclear whether its additive
  constant is Cailliez's. Include it only if its constant matches vegan's
  `$ac`.

**Metric MDS (SMACOF)**
- ✓ `smacof::smacofSym(type = 'ratio')`
- + `fmds::fastmds(dmm, p = 2)`. It took 13 s locally. Its `n.stress` is
  squared: `sqrt(0.0978) = 0.3127`, the same as ecodive's stress-1 of
  0.3128. Use the square root in the agreement table.
- `fmds::fastermds(dmm, p = 2, z = cmdscale(dm, 2))` ran in 0.06 s but
  returned an object with no names. Investigate before including it.

**NMDS**
- ✓ `vegan::metaMDS(try = 3, trymax = 3, parallel = 4)`, which with the
  cmdscale start makes 4 starts, the same as ecodive's `n_starts = 4`;
  `MASS::isoMDS`; `smacof::smacofSym(type = 'ordinal')`
- + `fmds::fastordinalmds(dmm, p = 2, z = cmdscale(dm, 2))`. It fails
  without `z`, an fmds bug. It took 10.3 s locally, with stress-1
  `sqrt(n.stress) = 0.2236` (ecodive 0.2237, vegan 0.2236).
- Leave out `ecodist::nmds`: it is pure R and far too slow at n = 1006.
  Mention why.

**Isomap (n_neighbors = 15)**
- ✓ `vegan::isomap(dm, ndim = 2, k = 15)`. With the same k it agrees to
  1.4e-14. (k = 16 was wrong.) vegan is the only implementation that takes
  distances. dimRed, Rdimtools and RDRToolbox take features only; say so.

**t-SNE**
- ✓ `Rtsne::Rtsne(dm, is_distance = TRUE, pca = FALSE, num_threads = 4)`
- + `scrapper::runTsne(knn90, num.threads = 4)`. It takes a kNN list whose
  `index` and `distance` are **neighbours × samples** (transposed), with 90
  neighbours (3 × perplexity) taken from `dmm`. Build it outside the timing,
  and note that this favours scrapper. It took 3.8 s locally.
- + `tsne::tsne(dm, k = 2)`: pure R, about 150 s for 1000 iterations. Use
  `iterations = 1`.
- + `manifoldsR`, if it installs. It takes features or its own kNN object.
  Check its docs for a precomputed-kNN entry point (`rs_tsne_from_knn`),
  and leave it out if there's none.
- Quality columns, for every package: final KL divergence (where reported)
  and `nn_kept()`, the share of each sample's 15 nearest neighbours that the
  embedding keeps.

**UMAP (min_dist = 0.1, 15 neighbours)**
- ✓ `uwot::umap(dm, batch = TRUE, n_threads = 4, n_sgd_threads = 4)`
- + `uwot::umap(dm, min_dist = 0.1, n_threads = 4, n_sgd_threads = 4)`,
  uwot's default (non-batch) optimiser, which is faster but not
  reproducible. Time both.
- + `umap::umap(dmm, input = 'dist')`. It took 4.4 s locally, and its
  `min_dist` default is already 0.1.
- + `scrapper::runUmap(knn15, num.threads = 4, min.dist = 0.1)`, with the
  kNN list transposed as above. It took 2.0 s locally.
- + `manifoldsR::umap`, if possible (as for t-SNE).

**db-RDA and capscale (`~ country + age`)**
- ✓ `vegan::dbrda`, `vegan::capscale`. No other R package takes distances.
  `fastrda` (OpenMP) takes raw tables only, and MDMR and `dbstats::dblm`
  are other models. Mention these as not comparable.

**PCA (CLR, 3 components)**
- ✓ `stats::prcomp(rank. = 3)`, `irlba::prcomp_irlba`, `vegan::rda`,
  `ade4::dudi.pca`
- + `scrapper::runPca(t(clr), number = 3, num.threads = 4)`. It takes
  features × samples, and its scores are `t(res$components)`. It agrees with
  prcomp to 2e-10 and took 0.05 s locally, so it is likely the closest
  competitor.
- + `RSpectra::svds(clr, k = 3, opts = list(center = TRUE))`. 0.08 s.
- + `bigstatsr::big_randomSVD(fbm, bigstatsr::big_scale(center = TRUE, scale = FALSE), k = 3, ncores = 4)`,
  where `fbm <- bigstatsr::as_FBM(clr)` is built outside the timing. It
  took 11 s locally on Windows, probably in cluster start-up; on Linux, time
  `ncores = 4` and `ncores = 1`.
- + `Rfast2::pca(clr, k = 3, scale = FALSE, vectors = TRUE)`. The default
  `scale = TRUE` gives a correlation PCA. 2.4 s.
- + `rsvd::rpca(clr, k = 3, scale = FALSE)`. It is randomized and agrees
  with prcomp only to about 0.05; report that.
- + `FactoMineR::PCA(clr, scale.unit = FALSE, ncp = 3, graph = FALSE)`. 0.5 s.
- `pcaone::pcaone(clr_c, k = 3, method = 'alg1')` works, where `clr_c` is
  the column-centred `clr`. It is approximate (4e-5). The default method
  fails on a matrix this small ("window size…"), so include it with
  `method = 'alg1'` and say so.
- `mlpack::pca(input = clr, new_dimensionality = 3)` returned a 1006 × 0
  matrix. Investigate: try a transposed input, or reading the result as a
  list. If you can't fix it, leave it out with a note.

**envfit** (fit `diarrhea`, `age` and `country` onto a 2-axis PCoA)
- ✓ `vegan::envfit` (it has no `parallel` option)
- + `ecodist::vf(ord, meta['age'], nperm = 999)`. It fits numeric variables
  only, so label it "ecodist (age only)". It took 2.1 s locally.

**PERMANOVA (`~ country`)**
- ✓ `vegan::adonis2`, serial and with `parallel = 4`. On Linux, `parallel`
  forks, so it should help more than on Windows.
- + `GUniFrac::adonis3(dmm ~ country, data = meta, permutations = 999)`.
  The left-hand side must be a matrix, not a `dist`. 13 s.
- + `LDM::permanovaFL(counts ~ country, data = meta, dist.list = list(dmm), n.perm.max = 999, n.cores = 4)`
  failed locally with "subscript out of bounds" (n.cores = 1, Windows).
  Investigate its formula interface; its own examples use `otu.table ~ covariates`.
  It also stops early once `n.rej.stop` (default 100) permutations have
  been rejected, which doesn't happen here (p = 0.001). Note that its
  statistic is PERMANOVA-FL, a variant.

**PERMDISP, ANOSIM, MRPP**: vegan only (✓, serial and `parallel = 4`).
`Blossom::mrpp` was archived from CRAN in 2018; mention it.

**Mantel (`dm` against `age_dm`)**
- ✓ `vegan::mantel` (serial and `parallel = 4`), `ade4::mantel.rtest`, `ecodist::mantel(nboot = 0)`
- + `ape::mantel.test(dmm, as.matrix(age_dm), nperm = 999)`: 91 s.
- + `ncf::mantel.test(M1 = dmm, M2 = as.matrix(age_dm), resamp = 999, quiet = TRUE)`: 29 s.
- These two are slow, so use `iterations = 2` for this method.

**No R competitor, so not timed:** `pcoa_biplot()` (its reference is
scikit-bio) and `rpca()` (its reference is gemelli). The test-suite table
covers both.

### 5.3 Determinism (already in the script; extend it)

`det_check()` runs each function five times: twice with no `set.seed()`;
twice after `set.seed(1)`; and once after `set.seed(1)` with 1 thread instead
of 4. Add the new packages that have random or threaded parts:

- `scrapper` t-SNE, UMAP and PCA (`num.threads`)
- `umap::umap`
- `rsvd::rpca` and `pcaone`, which are randomized
- `bigstatsr::big_randomSVD` (`ncores`)
- `GUniFrac::adonis3`, `LDM::permanovaFL` (if fixed), `ape::mantel.test` and `ncf::mantel.test`
- `fmds::fastermds` (if included), which is stochastic

Use 99 permutations in the determinism checks. Leave out `tsne::tsne`,
which is too slow.

Local results, on Windows, for comparison:
- Every ecodive function was identical across reruns, seeds and thread counts.
- Not identical between runs without `set.seed()`: Rtsne, uwot (both
  modes), irlba, and vegan's adonis2, permutest, anosim, mrpp and mantel.
- uwot's default threaded optimiser differed **even with the same seed**,
  and between 1 and 4 threads.
- Identical without a seed: metaMDS (its best solution came from the
  deterministic cmdscale start), MASS, ade4 and ecodist Mantel. Their
  permutations don't come from R's random stream, which is worth checking
  and saying.

Re-measure all of this on Linux. Report only what you measure.


## 6. Steps

1. Set up the environment (§4). Write `vignettes/articles/benchmark_ord_stats/ord_stats_machine.txt`.
2. Add the §5.2 "+" competitors and the §5.3 determinism checks to
   `vignettes/articles/benchmark_ord_stats/ord_stats.r`. Add an agreement row for any new package whose
   output can be compared with ecodive's: torgerson, fmds stress, the
   scrapper and RSpectra PCA scores, and the adonis3 F and R².
3. **Dry run:** `QUICK=1 Rscript vignettes/articles/benchmark_ord_stats/ord_stats.r` skips the timings and
   runs only the agreement and determinism checks (about 10 minutes). Fix
   every error, and check that no agreement row says `DIFFERENT`.
4. **Full run**, pinned to 4 cores (§4.1), in the background, logging to
   `vignettes/articles/benchmark_ord_stats/ord_stats.log`. Expect about 1–2 hours; the slowest calls are
   ape's Mantel (91 s), tsne (150 s) and vegan's ANOSIM (47 s per run).
   It writes:
   - `vignettes/articles/benchmark_ord_stats/ord_stats_results.rds`, `_agreement.rds`, `_determinism.rds`, `_session.txt`
   - `vignettes/articles/figures/benchmarks_ord_stats.svg`
5. Optional: the reference-BLAS run (§4.2).
6. **Finish the article.**
   - Replace `DETERMINISM_NOTES` with a few sentences on what the table
     shows, naming the packages that weren't deterministic and how.
   - Replace `SPEED_NOTES` with the speed-up ranges per method, and where
     ecodive is not fastest.
   - Rewrite "Methodology" for the cloud machine and its BLAS.
   - Check every claim in the text against the measured tables.
   - Check that the figure is readable at about 16 methods with up to 10
     packages each; adjust `ncol` and the height in the script if not.
7. **Build and look at it:** run
   `pkgdown::build_article('articles/benchmark_ord_stats')`, then open the
   HTML and check:
   - the tables, the figure and the expandable code and session blocks;
   - that the links to the test files resolve
     (`https://github.com/cmmr/ecodive/blob/main/tests/tinytests/...`).

   Then add a link to the new page from `benchmark.Rmd`'s introduction, and
   from `vignettes/ordination.Rmd`'s "Speed and reproducibility" section.
8. Commit (§3): the script, the `.rds`, `.txt` and `.log` files, the SVG and
   the article. Report:
   - the headline numbers;
   - any competitor that was as fast as ecodive, or faster;
   - anything left out, and why;
   - the branch.


## 7. What the local dry run found (Windows, reference BLAS)

On the benchmark data, ecodive agrees with its reference to between 1e-16
and 7e-13, relative, for PCoA, Cailliez, PERMANOVA, PERMDISP, ANOSIM,
MRPP, Mantel, db-RDA, capscale, PCA, envfit and Isomap. Every p-value is
identical when given the same permutations.

The iterative methods compare as follows:

| Method | Quality measure | ecodive | Reference |
|:-------|:----------------|--------:|----------:|
| NMDS | stress | 0.2237 | vegan 0.2236, MASS 0.2250 |
| NMDS | Procrustes correlation with vegan | 0.986 | |
| Metric MDS | stress-1 | 0.3128 | smacof 0.3130 |
| Metric MDS | Procrustes correlation with smacof | 0.996 | |
| t-SNE | KL divergence | 0.972 | Rtsne 1.006 |
| t-SNE | 15-NN kept | 0.501 | Rtsne 0.499 |
| UMAP | 15-NN kept | 0.464 | uwot 0.467 |

Rough single-run times on Windows, before any benchmarking, as a sanity
check:

| Method | ecodive | Others |
|:-------|--------:|:-------|
| PCoA | 0.03 s | cmdscale 1.3 s, ape 3.0 s |
| Cailliez | 0.58 s | cmdscale 18 s |
| PERMANOVA | 0.03 s | adonis2 31 s, or 14 s with parallel = 4 |
| ANOSIM | 0.11 s | vegan 47 s |
| Mantel | 0.26 s | vegan 11 s, ade4 28 s, ecodist 3 s |
| NMDS | 8.9 s | monoMDS 2.1 s for one start, isoMDS 4.6 s |
| t-SNE | 0.63 s | Rtsne 1.1 s, scrapper 3.8 s |
| UMAP | 0.77 s | uwot 1.7 s, scrapper 2.0 s |
| PCA | 0.01 s | scrapper 0.05 s, RSpectra 0.08 s, irlba 0.14 s |
