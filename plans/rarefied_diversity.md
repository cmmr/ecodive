# `beta_div_rarefied()` and `alpha_div_rarefied()`: diversity over repeated rarefactions

**Priority:** 37. **Tier:** derived: loops over ecodive's own `rarefy()`,
`beta_div()` and `alpha_div()`. The averaging choices come from q2-boots
(BSD-3, read; see below); vegan's `avgdist()` (GPL-2) is output-only. Read
[README.md](README.md) first, and [compositional.md](compositional.md)
§ "The 2026-11-20 review".


## What it is

A single rarefaction throws reads away at random, so its distances carry
sampling noise. Rarefying many times and averaging the results keeps the
equal-depth control and removes most of that noise. This is mothur's
`dist.shared(subsample = T)`, vegan's `avgdist()` and QIIME 2's q2-boots
(Raspet et al. 2025). Schloss (2024, 2026) argues that rarefaction remains
the best control for uneven sequencing depth, and finds that it beats
robust Aitchison PCA. Rarefaction is also the most mentioned normalization
in microbiome papers, at about 3,600 a year (Europe PMC, research notes).
With these helpers, users can run both sides of that debate in one
package: `pcoa(beta_div_rarefied(counts, 'bray'))` next to `rpca(counts)`.


## Names

- README rule 2 (a companion is `<base>_<modifier>`, sorting with
  `beta_div()` / `alpha_div()`) and rule 5 (snake_case); neither name is in
  NAMESPACE. Rejected: `avgdist()` (vegan's; no alpha counterpart; doesn't
  sort with `beta_div()`), `rarefied_beta_div()` (doesn't sort), and
  `beta_div(times =)` (repetition is not an option of any metric in the
  references, which make it a separate function or action, and
  `beta_div()`'s formals are the metrics' shared interface via
  `mget(metric$params)`).
- **Same names elsewhere:** none known. Same purpose: `vegan::avgdist()`,
  q2-boots' `beta` / `alpha` / `core-metrics` actions, mothur's
  `dist.shared`.
- **No registry rows.** `METRICS` holds metrics (formulas over one table);
  these repeat any metric. They are plain functions that resolve `metric`
  with `match_metric()` (aliases and prefixes work). No wrapper passes
  `...` to them, so README's shadowing table doesn't apply.
- **Arguments:** `times` reuses `rarefy(times)` (README: reuse names with
  the same meaning) rather than an `n_<unit>` name. New, for
  parameter_names.md §1: `average` (how repeats are combined) and `sd`
  (also return the spread: an extra output, like `rpca(distances)`).


## Design

### Signatures

```r
beta_div_rarefied(
  counts, metric, depth = NULL, times = 100L, average = 'mean', sd = FALSE,
  margin = 1L, norm = 'none', pseudocount = NULL, power = 1.5, alpha = 0.5,
  tree = NULL, pairs = NULL, clust = NULL, seed = 0L, cpus = n_cpus() )

alpha_div_rarefied(
  counts, metric, depth = NULL, times = 100L, average = 'mean', sd = FALSE,
  cutoff = 10L, digits = 3L, tree = NULL, margin = 1L, seed = 0L,
  cpus = n_cpus() )
```

- `margin`, `norm`, `pseudocount`, `power`, `alpha`, `tree`, `pairs`,
  `clust`, `cutoff`, `digits` and `cpus` mean what they mean in
  `beta_div()` / `alpha_div()` (`@inheritParams`). `norm` is applied after
  each rarefaction. `alpha_div()`'s `norm` formal is passed to no metric,
  so it is left out.
- `depth`: one whole number ≥ 1. With `NULL`, it is `rarefy()`'s automatic
  depth, computed once from the totals. Move that code into an internal
  `auto_depth(sums)` and have `rarefy()` call it too. Per-sample depths
  (`qmp_rarefy()`'s use) are an error here.
- `times`: a whole number ≥ 1. Default 100, as `vegan::avgdist(iterations
  = 100)`. q2-boots has no default `n`.
- `average`: `'mean'` (default), `'median'`, or `'medoid'` (beta only).
  Matched like `norm`, case-insensitive with prefixes.
  - The default is the mean, as in `vegan::avgdist(meanfun = mean)` and
    mothur.
  - q2-boots offers `non-metric-mean`, `non-metric-median` and `medoid`
    for beta (default medoid), and `mean` / `median` for alpha (default
    median), all read from its source. The docs say how to get q2-boots'
    defaults.
- `seed`: iteration t (1, ..., `times`) is exactly `rarefy(counts, depth,
  seed = s_t, margin = margin, cpus = cpus, warn = FALSE)`. Here `s_t` is
  `rarefy(times =)`'s sequence `seed, seed + 1, ...`, wrapped as in
  `R/rarefy.r`. Iteration 1 is therefore `rarefy(counts, depth, seed)`.

### Steps

1. `validate_counts()`, which resolves `tree` from phyloseq/rbiom/TSE
   objects before rarefying (rarefied tables lose it). Then
   `assert_integer_counts()`, `match_metric(metric, div = 'beta' |
   'alpha')`, the `average` / `times` / `depth` checks, `validate_seed()`
   and `validate_cpus()`.
2. Compute the totals. Samples below `depth` are dropped, the same samples
   in every iteration, because `rarefy()` drops them. Record their names
   (or positions) in attribute `dropped`, and warn with their number. Fewer
   than 2 kept samples (1 for alpha) is an error. `pairs` indexes the
   distance matrix of the **kept** samples, as with
   `beta_div(rarefy(counts))`.
3. **Mean** (beta), with memory for at most three triangles:
   `acc <- acc + d_t` in iteration order, then `acc / times`. The result
   is bit-identical to `Reduce('+', dms) / times`. With `sd = TRUE`, add
   Welford's running mean and M2 (n − 1 denominator, as `stats::sd`).
   **Never hold all `times` matrices.**
4. **Median and medoid** (beta) need every iteration's value of a pair, so
   work in **blocks of pairs**:
   - Choose B pairs so that `B × times × 8` bytes stays under an internal
     budget (256 MB). Each block reruns the `times` rarefactions, which are
     exact under the same seeds, and computes only its pairs:
     `beta_div(r_t, pairs = block)`, giving column t of a `B × times`
     matrix.
   - **Median:** take each row's median with a small threaded C helper
     (`C_row_medians`: quickselect per row, `(a + b) / 2` for an even
     count, rows split independently of `cpus`). It must equal `apply(M,
     1, median)` to 1 ulp.
   - **Medoid** (q2-boots' `_medoid()`, read; BSD-3): accumulate the
     `times × times` matrix of summed squared differences between the
     iterations' condensed vectors, block by block (`dist(t(M))^2`, summed
     in block order). Take square roots, sum each row, and pick the first
     minimum (`np.argmin`). Then return that iteration's full matrix,
     recomputed from its seed, and record attribute `medoid` (its t).
   - Rarefying costs blocks × `times` runs. Say so in the docs, and say
     that the mean needs one pass.
   - If profiling shows that `beta_div(pairs =)` allocating a full-length
     result per call dominates the time, add an internal compact-output
     flag to `C_beta_div`. That is optional.
5. **Alpha:** keep the `n × times` matrix (small). The mean is the
   sequential sum divided by `times`; the median is `apply(M, 1, median)`;
   `sd` is `apply(M, 1, sd)`.
6. `clust` applies once, to the final matrix (`dist_clust()`).

The distances are computed with ecodive's own deterministic C, and the
averaging works in a fixed order, so results are identical for any `cpus`.
R's RNG is never touched.

### Output

- **Beta:** a `dist` of class `c('ecodive_dist', 'dist')` over the kept
  samples, with `beta_div()`'s `method` (the metric ID) and a
  `display$title` naming the rarefaction (see "Plotting"). Further
  attributes: `depth` (the depth used, chosen from the totals when
  `NULL`), `sd` (a `dist`, with `sd = TRUE`), `dropped`, and `medoid`.
- **Alpha:** a named numeric vector over the kept samples, as
  `alpha_div()`, with attributes `depth`, `dropped`, and `sd` (with
  `sd = TRUE`).

### Depth-sensitivity check (no new export)

`envfit()` already does it: `envfit(ord, data.frame(log_depth =
log(rowSums(counts))))` gives the r² of library size on any ordination,
with a permutation p-value, and `cor(ord, log(depth))` gives it per axis.
A new helper would duplicate these. Instead, as an optional last step,
put the recipe in this Rd page's examples and in the article (see Docs).
Revisit a `plot(ord, depth = )` overlay only if users ask for one.


## Reference implementations

- **ecodive itself** (exact): `rarefy()` + `beta_div()` / `alpha_div()`.
  The tests rebuild every average from those.
- **q2-boots** (github.com/caporaso-lab/q2-boots, main branch, 2026-11-20;
  BSD-3, "Copyright (c) 2023, QIIME 2"; file headers "Copyright (c) 2024,
  Caporaso Lab"). Read `_beta.py`, `_alpha.py` and `_resample.py`. They
  give the averaging options and the medoid's definition (Euclidean
  distance between condensed matrices, argmin of the sums). No code is
  ported, only those definitions, so no notice is needed; cite the paper.
  Its `replacement = True` mode (bootstrap) is not planned: `rarefy()`
  draws only without replacement.
- **vegan `avgdist()`** (2.7-8, GPL-2, output-only): mean (or median) of
  `iterations` rarefied `vegdist()`s via `rrarefy()`. Its draws come from
  R's RNG, so its results agree only within Monte Carlo error.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- **Beta:** `dist_result()` titles the matrix from the metric; then
  `set_display()` changes `display$title` to "<metric name> (rarefied
  <times>×)", e.g. "Bray-Curtis Dissimilarity (rarefied 100×)", so titles
  say so. The `sd` dist gets `title_of("SD", <that title>)`. `clust` sets
  `order` and `display$clust` as usual.
- **Alpha:** no `display`; it is a plain vector like `alpha_div()`'s.
- Downstream: `pcoa()` / `permanova()` of the beta result take its title
  through `ord_result(from =)` / `new_ecodive_test(from =)`, so theirs
  read "Principal Coordinates Analysis of Bray-Curtis Dissimilarity
  (rarefied 100×)". Test that.

## Tests (`tests/tinytests/test-rarefied_diversity.r`)

- `times = 1`: `identical()` (values and labels) to
  `beta_div(rarefy(counts, depth, seed), metric)`, and the same for alpha.
- **Mean:** `identical()` to `Reduce('+', lapply(rarefy(counts, depth,
  seed, times = k), beta_div, metric)) / k`, for Bray, Jaccard, Aitchison
  (with `pseudocount`) and weighted UniFrac (with `tree`). `sd` must match
  `apply(stack, 1, sd)` to 1e-12.
- **Median / medoid:** equal `apply(stack, 1, median)` (1e-15). The medoid
  is the iteration an R re-implementation of q2-boots' `_medoid()` picks,
  and its matrix is `identical()` to that iteration's. Force tiny blocks
  through the internal fitting function's block-budget argument, and check
  that multi-block results are `identical()` to single-block ones.
- **Statistics, with no reference package:**
  `alpha_div_rarefied(x, 'observed', depth, times = 2000)` agrees with the
  exact expected richness under rarefaction, Σᵢ (1 − C(N − Nᵢ, n) / C(N,
  n)), computed with `lchoose()`, to within 4 Monte Carlo SEs (from `sd`).
  The seeds are fixed, so the test is stable.
- **vegan fixture** (`data-raw/expected/rarefied_diversity.r`):
  `avgdist(x, sample = d, iterations = 5000)`. The ecodive mean with
  `times = 5000` lies within 4 combined SEs in every cell.
- `pairs`, `clust`, `norm` and `margin = 2`; dense, `dgCMatrix` and slam
  inputs; a phyloseq-like object with an embedded tree (where the suite
  can build one); a dropped sample (warning, `dropped`, labels).
- `identical()` for `cpus` 1, 2 and 4 (more at home), for every
  `average`.
- **Errors**, with messages: non-integer counts; `depth` ≤ 0, NA, or a
  vector; `times` 0, NA or 1.5; an unknown `average`; `'medoid'` for
  alpha; fewer than two samples left; an unknown metric.
- **Smoke:** one line each, `times = 3`.


## Docs

- One Rd page, `@name rarefied_diversity` (as `qmp`), title "Diversity
  Averaged over Repeated Rarefactions". It covers:
  - when to use it (uneven depth), the debate (McMurdie & Holmes 2014;
    Schloss 2024, 2026), the depth choice, and the trade-offs of `average`,
    including how to get q2-boots' defaults;
  - memory and time (one pass for the mean, blocks for median/medoid);
  - the depth-sensitivity recipe above;
  - "Same names in other packages" (`vegan::avgdist()`);
  - "Differences from q2-boots": no bootstrap mode, a different default
    `average`, PCG seeds.
- `pkgdown/_pkgdown.yml`: under "Utilities & Helpers", beside `rarefy`
  (compositional.md asks for them there).
- `R/rarefy.r`: its "Auto-Depth Selection" section says "the highest depth
  that retains at least 10%", but the code takes the **lowest** library
  size whose retained total (depth × samples kept) is ≥ 10% of all reads.
  Correct the wording, since both functions now share the rule.
- **Follow-ups from [compositional_docs.md](compositional_docs.md) (item
  34):**
  - in `rpca()`'s **Sequencing depth** section and in the compositional
    article's depth chunk, replace `robust_aitchison(rarefy(counts))` with
    `beta_div_rarefied(counts, 'robust_aitchison')`;
  - in `vignettes/articles/compositional.Rmd`, add a short "Rarefied or
    compositional?" comparison:
    `procrustes_test(pcoa(beta_div_rarefied(counts, 'bray')), rpca(counts))`.
- NEWS (2.4.0); NAMESPACE: two `export()` lines.


## References

- Raspet, I., Gehret, E., Herman, C., Meilander, J., Manley, A., Simard, A., Bolyen, E., & Caporaso, J. G. (2025). Facilitating bootstrapped and rarefaction-based microbiome diversity analysis with q2-boots. *F1000Research*, 14, 87. doi:10.12688/f1000research.156295.1
- Schloss, P. D. (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*, 9, e00354-23. doi:10.1128/msphere.00354-23
- Schloss, P. D. (2026). Rarefaction is better than robust Aitchison PCA and other compositional data analysis methods at controlling for uneven sequencing effort. *mSystems*, 11, e00930-26. doi:10.1128/msystems.00930-26
- McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology*, 10, e1003531. doi:10.1371/journal.pcbi.1003531
- Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. *Ecology*, 52, 577–586. doi:10.2307/1934145 (the expected-richness formula used in the tests)
