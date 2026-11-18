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
