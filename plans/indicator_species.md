# `indicator_species()`: indicator species analysis (IndVal, multilevel pattern analysis)

**Priority:** 50. **Tier:** **clean room**. The reference, CRAN
**indicspecies 1.8.0** (`multipatt()`, `strassoc()`), is **GPL (>= 2)**
(CRAN `DESCRIPTION` verified 2026-11-20). labdsv's `indval()` is GPL too.
Work only from the papers (Dufrêne & Legendre 1997; De Cáceres & Legendre
2009; De Cáceres, Legendre & Moretti 2010; Tichý & Chytrý 2006),
indicspecies' **documentation** (`man/*.Rd`, the vignette
"IndicatorSpeciesAnalysis"; both read 2026-11-20) and its **outputs**. Never
open indicspecies' `R/` directory. No permissive implementation was found:
PyPI has no `indval` / `indicspecies` package, and the GitHub repositories
found (`IndicPiper`, `RedIndVal`, `SegFlora`) are GPL-3 wrappers or
scripts; mothur's `indicator` command is GPL-3. Read
[README.md](README.md) first.


## What it is

An indicator species (feature) is one that characterises a group of sites
(samples): found mostly there (**specificity**, A) and in most of its
samples (**fidelity**, B). Dufrêne & Legendre (1997) define IndVal = A × B
for single groups and test it by permutation. De Cáceres & Legendre (2009)
add the group-size-equalized IndVal.g and correlation indices (phi /
point-biserial r, and r.g); De Cáceres et al. (2010) extend the search to
**combinations of groups** ("multilevel pattern analysis"), so that a
feature shared by two of three habitats is found as an indicator of that
pair. indicspecies' `multipatt()` is the standard implementation, used in
microbiome papers as a per-feature association with sample groups
(soil, plant, water and gut studies).

Demand (retrieved 2026-11-20): Dufrêne & Legendre (1997) 2,642 OpenAlex
citations 2022–26; De Cáceres & Legendre (2009) 2,494; 1,487 Europe PMC
microbiome full-text mentions 2022–26 of "indicator species", indicspecies
or IndVal. The statistics are simple (group sums and counts per feature);
the work is pinning the formulas down from outputs and fitting the
permutation test to ecodive's machinery.


## Names

- **`indicator_species()`.** The method has no standard acronym: "ISA"
  is ambiguous, and `indval` is labdsv's function for the single-group
  IndVal only (rule 6: don't reuse a name whose meaning differs, since this
  function also runs r and r.g and group combinations). Spelled out (rules
  4 and 5), as `principal_response_curves`. `indval` and `multipatt` become
  aliases if it joins a registry (below). "Same names in other packages":
  none (`labdsv::indval()`, `indicspecies::multipatt()` differ in name).
- **`index`**: the association index, `'indval_g'` (default, multipatt's
  `"IndVal.g"`), `'indval'`, `'r'`, `'r_g'`; indicspecies' `func` values in
  snake case. Kept as an argument because every reference treats it as one
  (README: one function per variant only where references split them).
  Not `func` (ecodive's registries use `func` for a function) and not
  `method`. **Ask** if the user prefers `func`.
- **Registry: join `da_test()`?** It fits the shape: one categorical
  `groups`, a per-feature table with `p_value`, `p_adjusted`, `reject`,
  `alpha`, `method`; arguments clear the wrapper rule. Effect column `stat`,
  base `NA`, `zero_free_ok = FALSE` (B, the fidelity, counts zeros; a
  zero-replaced table makes every feature present everywhere). It is a test
  of association with groups rather than of differential abundance, and
  `da_sensitivity()`'s pseudocount grid has nothing to vary. Recommended:
  join (users list IndVal among DA tools, and `da_test()` comparisons across
  methods are useful), with aliases `indval`, `multipatt`; **ask the user**.


## Design

### Signature

```r
indicator_species(counts, groups, data = NULL, margin = 1L,
                  index = c('indval_g', 'indval', 'r', 'r_g'),
                  min_order = 1L, max_order = NULL, combinations = NULL,
                  negative = FALSE, alpha = 0.05, p_adjust = 'holm',
                  permutations = 999L, strata = NULL, seed = 0L,
                  cpus = n_cpus())
```

- `groups`, `data`: `test_group_values()` (as `simper()`; `ref =
  '`counts`'`).
- `min_order`, `max_order`: multipatt's `min.order`, `max.order` (sizes of
  group combinations considered; `max_order = 1` is `duleg = TRUE`, the
  original Dufrêne & Legendre analysis; `NULL` = the number of groups).
  multipatt's `duleg` is not kept (it equals `max_order = 1`).
- `combinations`: multipatt's `restcomb`, but as **names** (`'a'`,
  `'a+b'`) or the integer ids of the result's `combination_id`; `NULL` =
  all within the orders.
- `negative`: multipatt's `allow.negative` (r, r.g only).
- `permutations`, `strata`, `seed`: the shared permutation arguments;
  `perm_matrix()` and `src/permute.c` (`perm_source_init()`,
  `perm_source_get()`), so a supplied matrix makes p-values comparable
  exactly with indicspecies (`multipatt(permutations = )`), and `strata`
  restricts shuffles within blocks (indicspecies via `permute::how(blocks
  = )`, ✱ confirm the equivalence).
- `p_adjust`: over features (multipatt reports raw p-values; its vignette
  warns about multiple testing). ecodive's DA default `'holm'`; `'none'`
  reproduces indicspecies. **Ask**, with `simper()`'s precedent (it added
  `p_adjust = 'holm'` that vegan lacks).
- `alpha`: behind `reject` (`p_adjusted <= alpha`), as the DA methods;
  multipatt's `summary(alpha = 0.05)`.

### Indices (from the papers; ✱ = confirm from `strassoc()` / `multipatt()` outputs)

For feature \eqn{j}, abundances \eqn{a_i} over \eqn{n} samples, \eqn{K}
groups of sizes \eqn{n_k}, and a combination \eqn{C} (a set of groups,
\eqn{n_C = \sum_{k \in C} n_k} samples):

- **IndVal** (non-equalized; Dufrêne & Legendre 1997, De Cáceres &
  Legendre 2009 Appendix): \eqn{A = \sum_{i \in C} a_i / \sum_i a_i},
  \eqn{B = \#\{i \in C: a_i > 0\} / n_C}; `stat` \eqn{= \sqrt{AB}}
  (multipatt reports the square root; documented).
- **IndVal.g** (group-equalized): \eqn{\bar a_k} the mean in group \eqn{k};
  \eqn{A_g = \sum_{k \in C} \bar a_k / \sum_{k=1}^{K} \bar a_k}; \eqn{B} as
  above ✱ (whether B is also equalized for a combination, i.e. the mean of
  per-group fidelities, or pooled over \eqn{C}; De Cáceres et al. 2010 is
  the source; settle on outputs with unequal group sizes).
- **r** (phi for presence-absence, point-biserial for abundances): the
  Pearson correlation between \eqn{a} and the indicator of \eqn{C}.
- **r.g** (group-equalized, Tichý & Chytrý 2006 generalized by De Cáceres
  & Legendre 2009): the same correlation with every group given equal
  total weight, i.e. each sample of group \eqn{k} weighted \eqn{n/(K n_k)}
  ✱ (confirm this weighting reproduces `strassoc(func = "r.g")`).
- **Best combination**: the \eqn{C} with the largest index (with
  `negative`, the largest |r|, ties to the combination of fewer groups, as
  documented); ties otherwise ✱ (likely the first in combination order).
  When the best \eqn{C} is all groups, the p-value is `NA` (documented).
- **Combination order** ✱: multipatt numbers singletons first, then pairs,
  ... (the vignette's three-group example lists 7); confirm the order
  within a size (lexicographic by group order) from `$comb`, so
  `combination_id` matches multipatt's `index`.

### Permutation test (clean-room protocol first)

Permute the group labels (with `strata`), recompute, and count. Two
readings of "tests that best-matching pattern" are possible and must be
decided from outputs with a supplied permutation matrix before writing C:

- (a) the permuted statistic is the **same** combination's index (the
  observed best \eqn{C}), or
- (b) it is the **maximum over combinations** in each permutation (the
  multilevel analogue of labdsv's max-over-groups statistic).

Also confirm \eqn{p = (1 + \#\{\text{perm} \ge \text{obs}\})/(1 +
n_{perm})} and the tie tolerance (`perm_p()` uses vegan's
\eqn{\sqrt{\epsilon}}). Record the answers here. Crafted probe: three
groups, one feature, 20 hand-written permutations where (a) and (b) give
different counts.

### C and memory

- `src/indicator_species.c`: transpose the table by feature once (as
  `da_table_init()`, `norm = none`; never dense). Per feature, per
  permutation: group sums \eqn{\sum a}, sums of squares (r), counts present
  (B) by one pass over its non-zeros with the permuted group codes, then
  every combination's index from the \eqn{K} group totals (O(\eqn{2^K})
  per feature; combinations as bit masks; refuse \eqn{K > 20}).
- Threads over **blocks of permutations** (as `simper.c`): each task
  regenerates its permutations with `perm_source_get(ps, k, p)` and keeps
  integer exceedance counts per feature, so the combined counts are
  identical for any `cpus`. The observed statistics are permutation 0.
- Memory: O(D × threads) counts plus O(\eqn{2^K}) per thread; the
  feature × combination strengths are not kept.

### Result

`c('ecodive_indicator_species', 'data.frame')`, one row per feature:
`feature, combination` ("a+b", groups in factor order), `combination_id`
(multipatt's `index`), `n_groups` (the order), `stat`, `a`, `b` (IndVal
components, `NA` for r / r.g), `p_value, p_adjusted, reject`. Attributes:
`method = 'indicator_species'`, `index`, `groups` (levels), `alpha`,
`p_adjust`, `n_perm`, `permutation` (`'free'`, `'within strata'`,
`'supplied matrix'`, as `simper()`), `combinations` (the list
considered, id → name), `means` (group × feature means, for the plot).
Printing: `print_table_result`, header "Indicator species (IndVal.g)".


## Reference implementations

- **Output-only:** indicspecies 1.8.0 (`multipatt()`, `strassoc()`,
  `signassoc()`), labdsv (`indval()`, cross-check of `max_order = 1` with
  `index = 'indval'`).
- **Fixtures** (`data-raw/expected/indicator_species.r`): simulated
  tables, not indicspecies' `wetland` data (GPL package data must not be
  copied into ecodive's tests): 40 samples × 30 features in three unequal
  groups, abundances and presence-absence, with every `func`,
  `duleg = TRUE`, `max.order = 2`, `restcomb`, `allow.negative = TRUE`
  (r, r.g); a permutation matrix from `C_perm_indices()`'s convention
  passed as `permutations =`; `strassoc()` for every index on unequal group
  sizes (to pin the ✱ formulas); the probe above. Text or small `.rds`
  under `tests/tinytests/expected/`.
- **Clean-room record:** list in this plan, when implementing, each ✱ and
  the output that settled it (as linda's plan does).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `plot.ecodive_indicator_species(x, type = c('bar', 'heatmap'),
  n_features = 30, ...)`:
  - `'bar'` (default): horizontal bars of `stat` for the rejected
    features, grouped and coloured by `combination`, largest first within
    a combination.
  - `'heatmap'`: rejected features × groups, tiles shaded by the group
    mean (from `means`, each feature scaled to its maximum), with the
    feature's combination outlined.
- Title: "Indicator species (IndVal.g)" (index from the `index`
  attribute); the input is a count table, so no `display`. Reads `alpha`
  (subtitle: number rejected at `alpha`), `index`, `groups`, `means`.


## Tests

- `stat`, `a`, `b` equal `multipatt()` / `strassoc()` at 1e-12 for each
  index; best combinations identical; p-values identical with a supplied
  permutation matrix; `combination_id` equals multipatt's `index`.
- `max_order = 1`, `index = 'indval'` equals labdsv's IndVal (with its
  statistic rule, if (b) holds).
- Invariances: `index = 'indval_g'` is unchanged by duplicating every
  sample of one group (group means unchanged) ✱ (if B is pooled this
  holds; check).
- `identical()` for `cpus` 1, 2, 4 and any `seed` reused; dense = sparse;
  `margin = 2`.
- Edge cases: two groups; a group of one sample; a feature absent
  everywhere (`NA`); a feature present in every sample (best = all groups,
  p `NA`); `permutations = 0` (p `NA`); `combinations` naming an unknown
  group (error); `K > 20` (error); `negative = TRUE` with an IndVal index
  (error, as documented).
- Smoke line.


## Docs

- Rd page: A and B in words (specificity, fidelity) and formulas; IndVal
  vs IndVal.g (equalized for unequal group sizes; recommended by the
  vignette); r / r.g; group combinations; the permutation test and what it
  tests; multiple testing; that IndVal on relative abundances depends on
  normalization (A uses abundances, B only presence), so give a
  normalized table (`normalize()`), and rarefaction or depth differences
  change B.
- **"Differences from indicspecies"**: `p_adjust`; combinations named,
  not only numbered; the permutation streams (ecodive's PCG, same p-values
  with a supplied matrix); `duleg` folded into `max_order = 1`; no
  bootstrap confidence intervals (`strassoc(nboot.ci)`) and no
  `indicators()` (species combinations as indicators): later steps if
  asked.
- `parameter_names.md`: `index`, `min_order`, `max_order`,
  `combinations`, `negative` (new).
- NEWS; pkgdown (with `simper()` under "Interpreting Differences", or in
  "Differential abundance" if it joins `DA_TESTS`).


## References

- Dufrêne, M., & Legendre, P. (1997). Species assemblages and indicator
  species: the need for a flexible asymmetrical approach. *Ecological
  Monographs*, 67, 345–366.
  doi:10.1890/0012-9615(1997)067[0345:SAAIST]2.0.CO;2 (2,642 OpenAlex
  citations 2022–26, retrieved 2026-11-20).
- De Cáceres, M., & Legendre, P. (2009). Associations between species and
  groups of sites: indices and statistical inference. *Ecology*, 90,
  3566–3574. doi:10.1890/08-1823.1 (2,494 OpenAlex citations 2022–26).
- De Cáceres, M., Legendre, P., & Moretti, M. (2010). Improving indicator
  species analysis by combining groups of sites. *Oikos*, 119, 1674–1684.
  doi:10.1111/j.1600-0706.2010.18334.x
- Tichý, L., & Chytrý, M. (2006). Statistical determination of diagnostic
  species for site groups of unequal size. *Journal of Vegetation
  Science*, 17, 809–818. doi:10.1111/j.1654-1103.2006.tb02504.x


## Decisions made while implementing (2026-11-24)

**Asked the user at the start** (all three answered with the
recommendation):

1. **Join `da_test()`**: a `DA_TESTS` row, `indicator_species`, design
   `groups`, effect `stat`, base `NA`, `zero_free_ok = FALSE`,
   `signed = FALSE`; aliases `indval`, `multipatt`, `indicatorvalue`,
   `multilevelpatternanalysis` (checked against every registry: no
   collisions). pkgdown lists it under "Differential abundance".
2. **`index`**, not `func`. indicspecies' spellings (`'IndVal.g'`, `'r.g'`)
   are accepted too, matched ignoring case and non-letters.
3. **`p_adjust = 'holm'`** by default; `'none'` reproduces multipatt.

**Clean-room record: each ✱ and the output that settled it.** Probes are
in the session's scratch scripts; the decisive ones are re-run by
`data-raw/expected/indicator_species.r` and checked by the tests.

- **Formulas.** On 40 × 30 abundance and presence-absence tables in three
  unequal groups (8, 13, 19), multipatt's `$str` for all 7 (6)
  combinations equals: IndVal `sqrt(A B)` with `A = S_C / S`,
  `B = P_C / n_C`; IndVal.g with `A = sum_C m_k / sum_k m_k` and **B
  pooled over C** (`P_C / n_C`, not the mean of per-group fidelities); r =
  Pearson's correlation with the indicator of C; r.g = the weighted
  correlation with sample weights `1 / n_k` (each group weighs the same).
  Differences ≤ 1.4e-15 (exactly 0 for IndVal). So the plan's invariance
  test (duplicating one group's samples leaves IndVal.g unchanged) holds
  only for single groups (pooled B changes for a combination); the test
  uses `max_order = 1`.
- **labdsv's `indval()`** is `index = 'indval_g'` (not `'indval'`), single
  groups, squared (indicspecies' docs say so; the outputs agree to 1e-15).
  The plan's "max_order = 1, index = 'indval'" cross-check is therefore
  done with `'indval_g'`.
- **Permutation statistic: reading (b).** With a supplied permutation
  matrix (199 random permutations), multipatt's p-values equal the count
  of permutations whose **maximum over the considered combinations** is ≥
  the observed best value, with sample `i` taking group
  `cluster[perm[i]]`; reading (a) (the same combination) gave different
  counts. With `allow.negative = TRUE`, the maximum of |value|. Which
  combinations are maximized over: those considered (`max.order`,
  `restcomb`). The all-groups combination's IndVal is invariant under
  permutation, so including it or not changes nothing.
- **p-value formula and ties.** `(1 + #{perm ≥ obs}) / (1 + n_perm)`,
  compared **exactly**: a one-permutation probe whose permuted value is
  below the observed one by 1e-15 … 1e-2 (relative ~3e-17 … 3e-4) gives
  p = 1/2 every time. I first used a relative 1e-12 tolerance (ecodive's
  usual practice); it counted two permutations for one feature of the
  IndVal/duleg fixture whose values equal the observed one mathematically
  (from different group sums) but not in floating point, where multipatt
  doesn't count them. ecodive's IndVal arithmetic is bit-identical to
  multipatt's on these data, so exact comparison matches it; it is now
  exact. A permutation that leaves the best combination's samples in it
  gives bit-identical group sums (same summation order), so it always
  counts.
- **Best combination and ties.** multipatt breaks ties **at random**: the
  same data gave index 1 in one run and 4 in the next, and a probe with a
  perturbation of relative size e showed random choices for e ≤ 2e-5 and a
  fixed one for e ≥ 1e-4 (consistent with `max.col()`'s relative 1e-5
  tolerance; not confirmed from source). ecodive breaks ties
  deterministically: values within a relative 1e-12 of the largest are
  ties (so a combination and its complement, equal in |r| up to rounding,
  tie), and ties go to the first combination (fewest groups, then group
  order), which is multipatt's documented rule for `allow.negative`
  (confirmed: it always chose the smaller order). The fixture flags
  features whose two best values are within 1e-5 (`near_tie`: 4 of the
  920 feature rows over all cases, e.g. F10, the same in every sample,
  with single groups; complement pairs excluded for `allow.negative`);
  the tests compare their values but not their combinations or
  p-values.
- **Combination order and ids.** Singletons first, then pairs, ...,
  lexicographic within a size (`colnames(mp$comb)` for 3 and 4 groups:
  `1+2, 1+3, 1+4, 2+3, ...`). Character groups are sorted, factor groups
  keep their level order (multipatt's columns follow it). **multipatt's
  `index` is the position among the combinations considered**, not among
  all 2^K - 1: with `min.order = 2` or `restcomb = c(2, 5, 11, 15)` the
  first considered is index 1. So `combination_id` is the id among all
  combinations (what `restcomb` takes), equal to multipatt's `index` only
  when nothing before it is left out (`max_order` alone keeps them
  equal). The tests compare combination names.
- **`restcomb` with `max.order`** intersect (`restcomb = c(1, 2, 4, 7)`,
  `max.order = 2` considers 1, 2, 4): ecodive's `combinations` and the
  orders intersect too.
- **All groups.** For the correlations multipatt never considers the
  all-groups combination (undefined), even with `max.order = K` or
  `restcomb` naming it (it is dropped silently); ecodive caps `max_order`
  at K - 1 for them and drops it from `combinations` (an error only if
  nothing is left, or `min_order = K`). For the IndVal indices, a best
  combination of all groups gets p = `NA` in multipatt — **except with
  `restcomb`**, where multipatt tests it and gives p = 1 (21 features of
  the `restcomb_all` fixture). ecodive always gives `NA`; documented under
  "Differences", and the tests map multipatt's 1 to `NA` there.
- **Undefined indices.** A feature absent everywhere: multipatt gives
  `NA` index, stat, p (all `$str` NaN); a feature with the same value in
  every sample: r and r.g undefined (NaN). ecodive gives `NA` combination,
  id, stat and p. `strassoc()` returns 0 for these instead of NaN (the
  tests check that ours is `NA` exactly there). `multipatt(allow.negative
  = TRUE)` **stops with an error** on such a feature
  ("missing value where TRUE/FALSE needed"; no source text was printed),
  so its negative fixtures leave them out.
- **`strata`** = `permute::how(blocks = )`: ecodive's within-strata
  permutations (`C_perm_indices` with strata codes) keep every sample in
  its block, and given to multipatt as a matrix they give identical
  p-values (the `*/strata` fixtures). `strata` and a matrix together are
  an error, as in `simper()`.

**Other decisions**

- **Too many combinations.** Refused when more than 2^20 combinations
  would be considered (so 21 groups with the default `max_order`, as the
  plan's "K > 20" test asks), rather than refusing K > 20 outright: 21+
  groups work with `max_order = 2`. K > 30 is refused (bit masks).
- **`combinations` by name**: groups joined with `"+"` in any order
  (`'b+a'` is `'a+b'`), and group names may themselves contain `"+"` (a
  small backtracking parser). Ids are ranked/unranked arithmetically, so
  no list of all 2^K - 1 combinations is built.
- **No `norm` argument** (the plan's signature has none): values are used
  as given, and the docs say to normalize first. `zero_free_ok = FALSE`,
  so `da_sensitivity(zeros = )` refuses it (its message says the method
  "models counts", which is loose for IndVal; left as is, shared code).
- **Result**: columns as planned; `combination_id` and `n_groups` integer;
  `a`, `b` `NA` for r / r.g; `reject` is `FALSE` where `p_adjusted` is
  `NA`. The plan's `combinations` attribute is the considered ids named by
  combination. `means` is groups × features. No `display` (a count-table
  input); prints and plots are titled "Indicator species (IndVal.g)" etc.
- **`indicator_species_details()`** (internal) also returns every
  combination's value (features × combinations), for the tests'
  comparison with `$str`; the exported function never keeps it.
- **C**: `da_table_init()` (norm none) gives the table by feature; per
  feature and permutation one pass over its non-zeros (stored zeros
  skipped) for S_k, P_k and (r.g) the sum of squares around the feature's
  mean, then O(K) per combination. r uses the invariant total sum of
  squares, r.g the weighted variance `sum_k Q_k / n_k - K (mu - M/K)^2`
  from sums around the feature mean (no `Q - S^2/n` cancellation).
  Threads over blocks of 16 permutations; integer counts only. Observed
  values over blocks of 64 features.
- **Fixture data** are simulated (not indicspecies' `wetland`), the
  permutations are ecodive's own (`C_perm_indices`, seed 7) so the tests
  check random and supplied permutations against the same p-values.


## Outcome (from the implementing session's final report)

- **Commits:** 9737ada (`indicator_species()`, C, docs, fixtures, tests),
  then this one (ticks and outcome).
- **Files:** `R/indicator_species.r`, `src/indicator_species.c`,
  `data-raw/expected/indicator_species.r`,
  `tests/tinytests/expected/indicator_species.rds` (34 KB),
  `tests/tinytests/test-indicator_species.r` (512 tests); a `DA_TESTS`
  row, NAMESPACE, print label, NEWS, pkgdown, smoke line,
  `parameter_names.md`.
- **Reference agreement** (indicspecies 1.8.0, labdsv 2.3.1; 33 multipatt
  cases over 4 indices, abundance / presence-absence / relative tables,
  `duleg`, `max.order`, `min.order`, `restcomb`, `allow.negative`, strata,
  3 and 4 groups): every combination's value within 1.4e-15 (exactly 0
  for IndVal); best combinations identical except the 4 of 920 feature
  rows that multipatt breaks at random (near ties); p-values identical
  with ecodive's own random permutations and with the same permutations
  supplied as a matrix (and within strata). strassoc() for all indices
  and A, A.g, B within 1e-12 (it gives 0 where multipatt and ecodive give
  NaN/NA); labdsv::indval() = `'indval_g'`, single groups, squared, within
  1e-12.
- **Tests:** 30,149 at home (baseline 29,628), 434 in CRAN mode (432);
  `R CMD check --no-manual --ignore-vignettes`: only the existing
  installed-size NOTE. Identical results for `cpus` 1, 2, 4, dense,
  `Matrix`, `slam` and `margin = 2`.
- **Open:** bootstrap confidence intervals (`strassoc(nboot.ci)`) and
  `indicators()` (species combinations) not provided; `da_sensitivity()`'s
  refusal message for `zeros` says the method "models counts".
