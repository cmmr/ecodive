# `simper()`: similarity percentages

**Priority:** 8. **Tier:** port from ecopy (MIT), and add the permutation
p-values from Clarke (1993). Read [README.md](README.md) first.


## What it is

SIMPER (Clarke 1993) breaks down the average Bray-Curtis dissimilarity
between two groups into each taxon's contribution. It ranks the taxa that
"drive" the difference. In high demand in marine and community ecology;
vegan `simper()`.

**Caveat to document:** Warton, Wright & Wang (2012) show that SIMPER
confounds location effects with dispersion. Highly variable taxa rank high
even when their means don't differ. Say so in the docs, and suggest
`permdisp()` and per-taxon models as complements.


## Design

- **Signature:** `simper(counts, groups, data = NULL, margin = 1L, permutations = 999L, strata = NULL, seed = 0L, cpus = n_cpus())`.
  - Bray-Curtis only: SIMPER is defined for it. Use raw counts, as vegan
    does, which doesn't transform them; if users want relative abundances
    they can normalise first. Consider a `norm` argument for consistency;
    decide and document.
  - `groups` is categorical, handled as in the tests (`test_groups()`).
- **Result:** a plain function, not a wrapper method.
  - A data frame with one row per **pair of groups × taxon**, ordered
    within each pair by decreasing contribution.
  - Columns: `pair` ("a-b"), `taxon`, `average` (mean contribution to the
    pair's Bray-Curtis), `sd`, `ratio` (average / sd), `mean_a`, `mean_b`,
    `cumsum` (cumulative proportion), and `p_value` (the permutation
    p-value; `NA` when `permutations = 0`).
  - Class `c('ecodive_simper', 'data.frame')`, as `envfit()`'s result has
    a marker class. Optionally add `p_adjusted` with `p_adjust = 'holm'`,
    across the taxa within each pair, consistent with the pairwise tests.
    Raise this with the user: vegan doesn't adjust.
- **Masking:** `vegan::simper`.


## Reference implementations

- **Port from:** ecopy 0.1.2.2, `simper()` (MIT,
  [repo](https://github.com/Auerilas/ecopy)). It gives per pair the mean,
  sd, ratio, % and cumulative % of each taxon's Bray-Curtis contribution.
  - **No permutation p-values**: add them from Clarke (1993) and vegan's
    documentation. Permute the group labels and recompute each taxon's
    average contribution; p is the proportion of permuted averages ≥ the
    observed one.
  - ecopy normalises % to the sum of species means; vegan's `cusum` is
    relative to the total average dissimilarity. **Match vegan.**
  - ecopy breaks on pandas ≥ 2, so read its algorithm rather than running
    it.
- **Output-only fixtures (GPL):** vegan 2.7-x `simper()` and
  `summary(simper(...))`, with a supplied permutation matrix so that
  p-values compare exactly.
- **C:** for each pair of groups, sum each taxon's contribution
  |x_ik − x_jk| / Σ_k (x_ik + x_jk) over every between-group sample pair.
  Parallelise over sample pairs with fixed-order reduction per taxon.
  Permutations multiply that cost, so use the thread pool, and keep sparse
  tables sparse: a taxon absent from both samples contributes 0.


## Tests

- Fixtures from vegan, written by `data-raw/expected/simper.r`: dune with
  Management (4 groups), plus a random sparse table with 3 groups.
  Compare `average`, `sd`, `ratio`, `ava`/`avb`, `cusum` and `p` at 1e-10
  (p identical with shared permutations).
- The sum of `average` over taxa equals the mean between-group
  Bray-Curtis dissimilarity, computed with `bray()`.
- Sparse and dense input agree; results are identical for any `cpus`.
- Edge cases: two groups; a group of one sample; a taxon absent
  everywhere.
- A smoke-test line.


## Docs

- An Rd page with the Warton et al. caveat; add it to the reference index
  under Hypothesis Tests, or a new "Interpreting differences" group.
- NEWS.


## References

- Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in
  community structure. *Australian Journal of Ecology*, 18, 117–143.
- Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based
  multivariate analyses confound location and dispersion effects.
  *Methods in Ecology and Evolution*, 3, 89–101.


## Decisions made while implementing (2026-11-18)

Asked the user at the start, who chose the recommended option each time:

- **`norm`:** added, `norm = 'none'` by default (vegan uses the values as
  given), with `bray()`'s choices except `'clr'` and `'rclr'`, which give
  negative values that Bray-Curtis can't use (an error names the allowed
  ones). `pseudocount` is therefore not needed and not offered.
- **`p_adjust`:** added, `'holm'` by default, with a `p_adjusted` column:
  `stats::p.adjust()` over the taxa of each pair, so taxa whose `p_value` is
  `NA` aren't counted. Matched as in the pairwise tests (`pairwise_method()`).
  `p_value` stays unadjusted and matches vegan.
- **Permutations across pairs:** as the pairwise tests. Each pair shuffles
  only its own samples' labels, with seed `pair_seed(seed, k)` for pair k; a
  supplied matrix is cut down to each pair by `induced_perms()`; `strata` are
  cut down (and recoded) to the pair's samples. vegan instead shuffles all
  the samples' groups for every pair (its supplied matrix has a column per
  sample), so its p-values for a pair depend on the other groups' samples.
  The fixtures therefore run vegan on each pair's samples alone with the
  induced permutations, as `data-raw/expected/pairwise.r` does.

My own calls:

- **vegan's p-value rule, found from its outputs:** p = (1 + #{permuted
  average **>** observed}) / (1 + n_perm); permutations giving the same
  average (only swapping samples within groups) don't count. This is
  anti-conservative for small groups (2 vs 2: p = 1/(n+1) for the top taxon
  where the exact p is about 1/3). ecodive first matched it, and flagged it
  as an open question.
- **Changed by the user's choice (2026-11-18, after dbf0c2b): Clarke's
  "≥" rule.** p = (1 + #{permuted average ≥ observed}) / (1 + n_perm),
  with ties judged within a relative 1e-7 (both are sums in a fixed order,
  so true ties differ only by rounding). ecodive's p-values are therefore
  never smaller than vegan's on the same permutations, and equal for taxa
  that no permutation ties. The tests compare them with a direct R
  computation of the rule on the fixtures' permutations, and check that
  they are at least vegan's stored ones.
- **`p_value` is `NA` for taxa with `average` 0** (absent from both groups,
  or the same value in every sample), as vegan's outputs give.
- **Result:** one data frame, pairs in factor-level order ("a-b"), taxa by
  decreasing `average` (ties in column order: `order(decreasing = TRUE)`,
  which gave vegan's order in every fixture). Columns: `pair, taxon,
  average, sd, ratio, mean_a, mean_b, cumsum, p_value, p_adjusted`. `taxon`
  is the column number as text when `counts` has no column names. `ratio`
  is `NaN` for absent taxa (0/0) as in vegan; `sd` and `ratio` are `NA`
  when both groups have one sample. Attributes: `overall` (per pair, named),
  `n_perm`, `p_adjust`, `permutation` (`'free'`, `'within strata'`,
  `'supplied matrix'`). Class `c('ecodive_simper', 'data.frame')`, with no
  print method (as `envfit()`).
- **`mean_a`, `mean_b`** are means after `norm`. vegan's `ava`/`avb` follow
  its own pair order (order of first appearance), so the fixture script
  swaps them where vegan's pair is "b_a".
- **Empty samples are an error** ("Samples with no counts: ..."), as in
  ecopy: Bray-Curtis is 0/0 between two of them. NA, negative or infinite
  values are an error too.
- **C (`src/simper.c`):** the pair's samples are re-laid out taxon by taxon
  (CSC), and each taxon's sum runs over only the samples that have it
  (about nnz × n/2 operations; absent taxa cost nothing). The observed sums
  and the sd (two-pass: squared deviations over the pairs where the taxon is
  present, plus the zero pairs' `mean^2` times their count, so no
  cancellation) are split over blocks of 256 taxa; permutations over tasks
  of 16 permutations, each summed by one thread, with only integer counts
  combined. Results are identical for any `cpus`. 300 samples × 3,000
  taxa, three pairs, 99 permutations: 1.8 s on one thread, 0.8 s on four;
  `vegan::simper()` takes 27 s for one pair.
- **Argument order:** `counts, groups, data, margin, norm, p_adjust,
  permutations, strata, seed, cpus` (`norm` after `margin` as in
  `plsda()`; `p_adjust` before `permutations` as in the pairwise tests).
- **Licensing:** ecopy's sdist has no licence file; its `setup.py`,
  `PKG-INFO` and README say MIT. The standard MIT text, with "Copyright (c)
  2017 Nathan Lemoine" (the author, and the 0.1.2.2 release year; no
  copyright line exists in the source), is in `src/simper.c`'s header and a
  new `inst/COPYRIGHTS` stanza; Nathan Lemoine is `ctb`/`cph` in
  `Authors@R`.
- **Docs:** a new "Interpreting Differences" group in the pkgdown reference
  index; the Warton et al. (2012) caveat, a "Differences from
  vegan::simper()" list, and a "Same names in other packages" section;
  NEWS entry.
- **References:** vegan 2.7.8 and permute 0.10.0
  (`data-raw/expected/simper.r`): dune / Management with 999 free and 199
  within-`Use` block permutations, and a sparse 25 × 40 table with three
  unequal groups (an absent taxon, a constant one), raw and with
  `norm = 'percent'` against vegan on `decostand(, 'total')`. Every column
  agrees to about 5e-16 relative and the taxa come in the same order.
  p-values were identical under vegan's rule; under Clarke's rule they
  match a direct R computation instead (see above).
- **Open issue (not changed):** when `groups` has the wrong length, the
  shared `per_sample()` message says "`dm` has 4 samples" even for count
  tables (also in `plsda()`).
