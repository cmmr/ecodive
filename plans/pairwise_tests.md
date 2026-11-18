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
