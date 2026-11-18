# Pairwise tests: `pairwise` and `p_adjust` on `permanova()`, `anosim()`, `mrpp()` and `permdisp()`

**Priority:** 2. **Tier:** port the approach from q2-diversity (BSD-3);
the statistics come from ecodive's own tests. Read [README.md](README.md)
first.


## What it is

After an overall test across several groups, test each pair of groups by
running the same test on just those two groups' samples. Microbiome
researchers use pairwise PERMANOVA in nearly every paper.


## Design (decided by the user, 2026-11-18)

- **New arguments on all four tests:** `pairwise = FALSE` and
  `p_adjust = 'holm'`.
  - `p_adjust` accepts any of `stats::p.adjust.methods` (`'holm'`,
    `'hochberg'`, `'hommel'`, `'bonferroni'`, `'BH'`, `'BY'`, `'fdr'`) or
    `'none'`, matched as `stats::p.adjust()` matches them.
  - Giving `p_adjust` explicitly without `pairwise = TRUE` is an error.
  - **The default is Holm for every test, `permdisp()` included.**
- **Output:**
  - A `pairwise` field holds a data frame with one row per pair, named
    like `"a-b"`. Its columns are the test's statistic (`F` and `R2` for
    PERMANOVA, `R` for ANOSIM, `delta` and `A` for MRPP), `n_perm`,
    `p_value` (raw) and **`p_adjusted`**. Raw values are never overwritten.
  - A `p_adjust` field records the method, and `print()` names it.
  - `NULL` when `pairwise = FALSE`.
  - The overall test (the core fields) isn't adjusted.
- **What gets adjusted:** each pairwise table's permutation p-values,
  across its pairs.
- **permdisp:** it already has `pairwise = TRUE`, with columns `t`, `Df`,
  `p_value` (from the t distribution) and `p_value_perm` (from the
  permutations), and adjusts nothing, as vegan's `permutest()` doesn't.
  Add `p_adjusted` and `p_adjusted_perm`, each adjusting its own column.
  This **changes `permdisp()`'s default output**, so it needs a NEWS
  entry.
- **How each pair is tested:**
  - Subset `dm` to the pair's samples and rerun the test (as q2-diversity
    does), with the same `permutations` count, `strata` (subset too) and
    `seed`.
  - For PERMANOVA with a `formula`, `pairwise` applies to one factor:
    require `groups` (or a single-factor formula), and error otherwise.
    Covariate-adjusted pairwise contrasts are out of scope for now.
  - Determinism: each pair's permutations come from the seed and the
    pair's index, never from the thread count. Pairs may run in parallel.
- **Through `dist_test()`:** `pairwise` and `p_adjust` pass via `...`, so
  `dist_test()` needs no change.
- **Pairs with too few samples:** e.g. a singleton group. Give `NA` for
  that pair, with a warning, rather than failing.


## Reference implementations

- **Port from:** q2-diversity, `beta_group_significance(pairwise=True)`,
  BSD-3, © 2016–2026 QIIME 2 development team,
  [LICENSE](https://raw.githubusercontent.com/qiime2/q2-diversity/dev/LICENSE).
  It subsets the distance matrix to each pair and reruns the test.
  - q2 uses Benjamini–Hochberg; ecodive defaults to Holm.
  - q2 handles one factor only.
- **Test statistics:** scikit-bio `permanova`, `anosim` and `permdisp`
  (BSD-3); ecodive's own tests already match vegan.
- **Output-only fixtures (GPL):**
  - pairwiseAdonis 0.4.1 (GPL-2|3): `pairwise.adonis()` and
    `pairwise.adonis2()`;
  - RVAideMemoire 0.9-83-12 (GPL-2): `pairwise.perm.manova()`;
  - vegan `permutest(betadisper, pairwise = TRUE)` for permdisp.

  With a shared permutation matrix (subset per pair), raw p-values must
  match exactly. Note that pairwiseAdonis draws new permutations per pair.


## Tests

- `p_adjusted == stats::p.adjust(p_value, method)` for every method, and
  equals `p_value` for `'none'`.
- Raw pairwise statistics and p-values match the fixtures:
  - for permanova, vegan `adonis2` on each subset (or pairwiseAdonis);
  - for anosim and mrpp, vegan on each subset;
  - for permdisp, vegan `permutest(pairwise = TRUE)`.

  Give the same permutation matrix, subset per pair.
- The overall results are unchanged when `pairwise = TRUE`.
- Invalid `p_adjust`, and `p_adjust` without `pairwise`, are errors with
  clear messages.
- A singleton group gives an `NA` row and a warning.
- Results are identical for `cpus` 1, 2 and 4.
- Add `pairwise = TRUE` to the existing smoke-test calls.


## Docs

- Document `pairwise` and `p_adjust` on each test, and add a "Pairwise
  comparisons" section to `vignettes/stat_tests.Rmd`.
- NEWS: the new arguments, and the change to `permdisp()`'s pairwise
  output.


## References

- Anderson, M. J. (2001). A new method for non-parametric multivariate
  analysis of variance. *Austral Ecology*, 26, 32–46.
- Holm, S. (1979). A simple sequentially rejective multiple test
  procedure. *Scandinavian Journal of Statistics*, 6, 65–70.


## Decisions made while implementing (2026-11-18)

Asked the user at the start, who chose the recommended option each time:

- **Pair seeds:** each pair's random permutations come from `seed` and the
  pair's index (its row in the table), not from `seed` alone. So pairs of
  equal size never share permutations, and a row can't be reproduced by
  running the test by hand on the subset with the same `seed` (the docs say
  so). `C_pair_seed()` in `src/permute.c` derives the pair's seed: one
  31-bit PCG draw from state `seed ^ 0x9E3779B97F4A7C15` and stream = the
  pair's index. The constant keeps these states apart from every test's own
  (a seed sign-extended to 64 bits). The internal `pair_seed()` exposes it,
  and the tests use it to rerun a pair by hand.
- **A supplied `permutations` matrix** is cut down to each pair by the
  *induced* permutation: each of the pair's samples goes to the first of the
  pair's samples after it on its permutation cycle (`induced_perms()` in
  `R/pairwise.r`). It is uniform when the full permutation is (checked by
  simulation: 60,000 draws, χ² p = 0.14) and keeps within-strata
  restrictions (tested). The fixture script writes the same thing as a plain
  loop, and the tests check that the two agree.
- **Pairs with a singleton group** are tested whenever the test can run on
  them (a 1-vs-5 PERMANOVA, MRPP or ANOSIM is reported). Only a pair the test
  fails on (an error, such as ANOSIM on two singletons) gets an all-`NA` row,
  with the warning "Pair a-b could not be tested: <the test's message>". A
  test that warns and returns `NA`s itself (PERMANOVA or MRPP on two
  singletons) has its warning passed on, prefixed "Pair a-b: ".
- **Extra columns `z` and `n`** (the pair's number of samples). The table's
  columns are: the statistics, `z`, `n`, `n_perm`, `p_value`, `p_adjusted`.
  So `F, R2, ...` for PERMANOVA, `R, ...` for ANOSIM, `delta, A, ...` for
  MRPP.

My own calls:

- **permdisp** got `z` (of `t` against `t_perm`, via `perm_z()`), `n` and
  `n_perm` too, for consistency. Its columns: `t, Df, z, n, n_perm,
  p_value, p_value_perm, p_adjusted, p_adjusted_perm`. Pairs with a
  singleton keep `NA` t statistics, as `permutest()` gives, but now warn:
  one warning per reason, listing the pairs ("Pairs a-d, b-d, c-d could not
  be tested: a group has only one sample, so it has no variance.").
  permdisp's pairs are still computed from the whole test (as `permutest()`),
  not rerun on subsets; its docs say so.
- **Argument order:** `pairwise, p_adjust` go just before `permutations` in
  all four (after `by` in permanova, `weights` in mrpp, `bias_adjust` in
  permdisp, `data` in anosim). No existing call passed `permutations` by
  position.
- **`p_adjust` matching:** `pmatch()` against `stats::p.adjust.methods`, as
  `p.adjust()`'s `match.arg()` does: exact names, or unique prefixes
  (`'bonf'`, `'hoc'`, `'f'`). An ambiguous prefix (`'h'`, `'B'`), `NA`, or a
  vector is an error that lists the choices. Giving `p_adjust` (even
  `'holm'`) with `pairwise = FALSE` is an error, as the plan says.
  `pairwise_method()` in `R/pairwise.r` does this for all four tests.
- **Pairs run one after another,** each test multithreaded over its
  permutations as usual, rather than pairs in parallel. Simpler, already
  deterministic for any `cpus`, and the per-pair work is the same.
- **Pair order and labels:** the factor's level order (sorted values for
  other vectors), as `utils::combn()`, labelled `"a-b"`; unused factor levels
  are dropped (no pairs for them). A factor's `NA` level is a group named
  "NA", as elsewhere.
- **PERMANOVA's pair** keeps the `R2` that `permanova()` reports on the
  subset, even when `F` is `NA` (two singletons give `R2 = 1`).
- **PERMANOVA's factor check** happens before any permutations run:
  `pairwise = TRUE` needs `groups` that are categorical, or a formula with
  one categorical term and no `Condition()`. Numeric groups, two terms, and
  `Condition()` are errors naming those options.
- **Subsetting `dm`** takes the pair's elements straight from the packed
  triangle (`dist_subset()`), with no n × n matrix.
- **`print()`:** a "Pairwise comparisons (p-values adjusted by Holm's
  method):" block after each test's own details ("Pairwise t tests (...)"
  for permdisp; "(p-values not adjusted)" for `'none'`).
- **Docs:** `pairwise` and `p_adjust` are shared `@param`s in
  `R/documentation.r` (permdisp overrides `pairwise`). A shared "Pairwise
  comparisons" section (`pairwise_section`, `@inheritSection`) explains the
  permutations, the adjustment and untestable pairs. A "Pairwise
  comparisons" section went into `vignettes/stat_tests.Rmd`, and a NEWS
  entry covers the new arguments and permdisp's extra columns.
- **References:**
  - vegan 2.6-4's `adonis2()`, `anosim()` and `mrpp()` (weights 1 and 3) on
    each pair's samples, with the induced permutation matrices, for Bray-Curtis
    and Euclidean distances, under free and within-block permutations
    (`data-raw/expected/pairwise.r`). Statistics agree to 1e-10, and
    p-values are identical.
  - RVAideMemoire 0.9-83-12's `pairwise.perm.manova()` (GPL; installed from
    CRAN into a scratch library; output only), with 9999 random
    permutations. Raw p-values agree within 0.015, and Holm-adjusted ones
    within 0.03.
  - pairwiseAdonis couldn't be used: it is on GitHub only, and attaching its
    repository was refused in this session. Its `pairwise.adonis2()` runs
    `adonis2()` on each pair's samples, which is what the vegan fixtures
    do.
  - permdisp's pairwise t tests were already checked against
    `permutest(pairwise = TRUE)`. Their adjusted columns are checked
    against `p.adjust()`.
- **Licensing:** no q2-diversity code was ported, only the approach (subset
  and rerun), so there is no `inst/COPYRIGHTS` stanza. `R/pairwise.r`'s
  header credits q2-diversity (BSD-3).
- **Fixture size:** `pairwise.rds` is about 48 KB. It holds the induced
  permutations of one case, for the cross-check against the fixture's
  loop.
- **Results:** the home suite (`NOT_CRAN=true`, `-O2` install, 4 cores) gave
  5580 results, all OK. CRAN mode gave 184, all OK. `R CMD check --no-manual
  --ignore-vignettes` gave only the known en_US.UTF-8 locale WARNING.
  `vignettes/stat_tests.Rmd` renders.
