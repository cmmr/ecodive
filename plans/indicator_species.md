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
