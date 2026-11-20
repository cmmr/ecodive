# `beta_partition()`, `beta_partition_multi()` and turnover / nestedness metrics: Baselga's partitioning of beta diversity

**Priority:** 52. **Tier:** formula + port (abdiv, MIT) for the pairwise
components; clean room for the multiple-site versions (betapart GPL,
output-only). Licences verified 2026-11-20:

- **betapart 1.6.1** (CRAN, published 2025-07-24): `License: GPL-2 |
  GPL-3 [expanded from: GPL (≥ 2)]` (CRAN page). Output-only: read its Rd
  pages and vignette, never its `R/`.
- **abdiv 0.2.0** (CRAN, 2020-01-20; github.com/kylebittinger/abdiv):
  `License: MIT + file LICENSE`; the LICENSE file reads "YEAR: 2020,
  COPYRIGHT HOLDER: Kyle Bittinger" (both checked). Its
  `R/beta_components.R` has `jaccard_turnover()`, `jaccard_nestedness()`,
  `sorenson_turnover()`, `sorenson_nestedness()`, `bray_curtis_balanced()`,
  `bray_curtis_gradient()`, `ruzicka_balanced()`, `ruzicka_gradient()` (and
  phylogenetic ones, see "Later"). It may be read and ported (pairwise only;
  plain R, one pair at a time). Not read while writing this plan.

The formulas are short and published (Baselga 2010, 2012, 2013, 2017), and
are given below. Read [README.md](README.md) first.


## What it is

Baselga (2010) splits a pairwise or multiple-site dissimilarity into
**turnover** (replacement of species by others) and **nestedness-resultant**
dissimilarity (one assemblage being a subset of the other). For two samples
with a shared species, b only in the first and c only in the second:

| Family | Total | Turnover | Nestedness |
|:-------|:------|:---------|:-----------|
| Sørensen (Baselga 2010) | β_sor = (b + c)/(2a + b + c) | β_sim = min(b, c)/(a + min(b, c)) (Simpson; Lennon et al. 2001) | β_sne = β_sor − β_sim = [(max − min)/(2a + b + c)] · [a/(a + min)] |
| Jaccard (Baselga 2012) | β_jac = (b + c)/(a + b + c) | β_jtu = 2 min(b, c)/(a + 2 min(b, c)) | β_jne = β_jac − β_jtu = [(max − min)/(a + b + c)] · [a/(a + 2 min)] |

The **abundance-based** analogue (Baselga 2013, 2017) uses, for abundances
x, y: A = Σ min(x_k, y_k), B = Σ max(x_k − y_k, 0), C = Σ max(y_k − x_k, 0):

| Family | Total | Balanced variation | Abundance gradient |
|:-------|:------|:-------------------|:-------------------|
| Bray-Curtis | d_BC = (B + C)/(2A + B + C) | d_BC-bal = min(B, C)/(A + min(B, C)) | d_BC-gra = d_BC − d_BC-bal |
| Ružička | d_RUZ = (B + C)/(A + B + C) | d_RUZ-bal = 2 min(B, C)/(A + 2 min(B, C)) | d_RUZ-gra = d_RUZ − d_RUZ-bal |

**Multiple-site** versions (Baselga 2010, 2012, 2017) summarize a whole set
of N sites in three numbers, with S_i the richness of site i, S_T the total
richness, and b_ij the number of species in i but not j:

- β_SOR = [Σ_{i<j} min(b_ij, b_ji) + Σ_{i<j} max(b_ij, b_ji)] /
  [2 (Σ_i S_i − S_T) + Σ min + Σ max]
- β_SIM = Σ min / [Σ_i S_i − S_T + Σ min];  β_SNE = β_SOR − β_SIM
- β_JAC = [Σ min + Σ max] / [Σ_i S_i − S_T + Σ min + Σ max]
- β_JTU = 2 Σ min / [Σ_i S_i − S_T + 2 Σ min];  β_JNE = β_JAC − β_JTU
- Abundance (Baselga 2017, betapart `beta.multi.abund()`): the same shapes
  with B_ij = Σ_k max(x_ik − x_jk, 0) for b_ij and
  A = Σ_k (Σ_i x_ik − max_i x_ik) for Σ S_i − S_T. ✱ Confirm A on betapart's
  outputs.

Because the multiple-site value depends on N, betapart's `beta.sample()`
averages it over random subsets of a fixed number of sites.

The partition is widely used (and debated: Carvalho et al. 2012 and Podani
& Schmera 2011 propose a different additive split, adespatial's
`beta.div.comp()`; the docs mention the alternative and the debate,
Baselga & Leprieur 2015).

**Demand** (retrieved 2026-11-20): Baselga (2010) 2,447 OpenAlex
citations 2022–26; betapart (Baselga & Orme 2012) 1,623; Europe PMC
microbiome full-text mentions 2022–26: "nestedness" AND "turnover" AND
"beta diversity" 101, "betapart" 98. Baseline: Bray-Curtis 28,581.


## Names and how it fits ecodive

**Each pairwise component is a beta metric** in `METRICS`, run by
`beta_div()` and usable everywhere a `dist` is (ordination, PERMANOVA,
Mantel), as the UniFrac variants are. Names follow rule 2 on the existing
base metric (`<base>_<modifier>`), as abdiv names them:

| id | Name (display) | weighted | alt_ids |
|:---|:---------------|:---------|:--------|
| `sorensen_turnover` | Simpson Dissimilarity (Sorensen Turnover) | NO | `simpsondissimilarity betasim` |
| `sorensen_nestedness` | Sorensen Nestedness | NO | `betasne` |
| `jaccard_turnover` | Jaccard Turnover | NO | `betajtu` |
| `jaccard_nestedness` | Jaccard Nestedness | NO | `betajne` |
| `bray_balanced` | Bray-Curtis Balanced Variation | YES | `bcbal braybal` |
| `bray_gradient` | Bray-Curtis Abundance Gradient | YES | `bcgra braygra` |
| `soergel_balanced` | Ruzicka Balanced Variation | YES | `ruzickabalanced ruzbal` |
| `soergel_gradient` | Ruzicka Abundance Gradient | YES | `ruzickagradient ruzgra` |

- `phylo NO`, `int_only NO`, `true_metric NO` for all (the components
  aren't metrics in general; ✱ check β_sim and β_jtu with a brute-force
  triangle search on random tables, and set YES only if proved).
- ecodive's `soergel()` **is** Ružička's dissimilarity (1 − Σ min / Σ max =
  vegan's quantitative Jaccard), so its components are `soergel_*`, with
  the Ružička names as aliases. **Ask the user** whether to also add
  `ruzicka` as an alias of `soergel` (recommended; collides with nothing
  today; check with the README snippet).
- `simpson` is an alpha metric (Gini-Simpson); `match_metric()` matches
  within the requested `div`, so the `simpsondissimilarity` alias is safe,
  but run the snippet.
- Note: `sorensen()`'s Rd says "2J/(A + B)", a similarity, while its C
  computes 1 − 2J/(A + B), the dissimilarity. Fix the Rd text in the same
  change (a doc bug found while planning).

**Two plain functions** give the whole partition at once:

- `beta_partition()`: the three pairwise matrices of one family (betapart's
  `beta.pair()` / `beta.pair.abund()`), computed in one pass.
- `beta_partition_multi()`: the multiple-site values (`beta.multi()`,
  `beta.multi.abund()`, and `beta.sample()` when `n_sites` is given).

`index` (new glossary row) selects the family, as betapart's
`index.family` does; it is an argument because betapart treats the family
as one (and the partition, not the metric, is the method). Not `family`,
which the glossary reserves for an outcome distribution (`coda_lasso()`).


## Design

### Signatures

```r
sorensen_turnover(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
# ... likewise sorensen_nestedness, jaccard_turnover, jaccard_nestedness

bray_balanced(counts, margin = 1L, norm = 'none', pseudocount = NULL,
              pairs = NULL, clust = NULL, cpus = n_cpus())
# ... likewise bray_gradient, soergel_balanced, soergel_gradient

beta_partition(counts, index = 'sorensen', margin = 1L, norm = 'none',
               pairs = NULL, cpus = n_cpus())

beta_partition_multi(counts, index = 'sorensen', margin = 1L,
                     norm = 'none', n_sites = NULL, n_draws = 100L,
                     seed = 0L, cpus = n_cpus())
```

- The incidence metrics take no `norm` (presence only), like `jaccard()`.
  The abundance metrics take `norm` like `bray()` (betapart uses the
  values as given; `norm = 'percent'` gives the relative-abundance
  version).
- `index`: `'sorensen'`, `'jaccard'` (incidence), `'bray'`, `'ruzicka'`
  (abundance; `'soergel'` accepted). `norm` is used only by the abundance
  indices (an error otherwise, as `simper()` refuses unusable norms).
- `beta_partition_multi()`: `n_sites` (betapart's `sites`) and `n_draws`
  (betapart's `samples`) resample `n_sites` sites without replacement
  `n_draws` times; `n_sites = NULL` computes the value for all sites only.

### C

Extend `src/beta_div.c`: the pairwise kernels need only (a, b, c) or
(A, B, C), which its macros already give.

- Incidence: `WITH_ABJ` gives A (richness of i), B (of j) and J (shared),
  so a = J, b = A − J, c = B − J. Four new `BDIV_*` codes with one-line
  formulas.
- Abundance: in `FOREACH_OTU`, accumulate Σ min(x, y) and the two totals;
  B = T_x − A, C = T_y − A. Four new codes.
- `beta_partition()` calls one kernel that writes **three** packed
  triangles in the same pass (a `C_beta_partition()` entry point with the
  `FOREACH_PAIR` macro, or `C_beta_div()` with an extra output), so the
  table is read once. Same pair split as every metric: identical for any
  `cpus`; `pairs` respected.
- Multiple-site: the sums over pairs of min(b_ij, b_ji) and max(...) come
  from the same per-pair (a, b, c), written to a packed triangle and then
  **summed serially in pair order** (no shared accumulator), plus
  per-feature passes for S_T and Σ_k max_i x_ik. `beta.sample()`-style
  draws: draw r picks its sites with PCG stream `(seed, r)` (Floyd's
  algorithm), computes its triangle over those sites only, and sums it;
  tasks are draws. Cost O(n_sites² S) per draw.
- Sparse tables stay sparse throughout.

### Results

- The eight metrics: `ecodive_dist`, via `dist_result()` (title from
  `METRICS`), as every metric.
- `beta_partition()`: a list of class `ecodive_beta_partition` with three
  `ecodive_dist` elements named `turnover`, `nestedness`, `total` (for the
  abundance indices: `balanced`, `gradient`, `total`), each a full metric
  result (e.g. `total` is identical to `sorensen()` / `jaccard()` / `bray()`
  / `soergel()`), plus attribute `index`. A short `print()` method (three
  lines with the means of each part), as list results write their own.
- `beta_partition_multi()`: class `ecodive_beta_partition_multi`, a list
  with `value` (named numeric: turnover, nestedness, total, for all sites)
  and, when `n_sites` is given, `draws` (an `n_draws` × 3 matrix) and
  `mean` (its column means); attributes `index`, `n_sites`, `n_draws`.
  Short `print()`.


## Reference implementations

- **Port:** abdiv 0.2.0 (MIT) for the pairwise component formulas,
  including its NaN conventions (its Rd: NaN when an input is empty or
  all zero). Cite it in the C file header and `inst/COPYRIGHTS` only if
  its code is followed closely; the formulas themselves are the papers'.
- **Output-only fixtures:** betapart 1.6.1 `beta.pair()`, `beta.multi()`
  (both families), `beta.pair.abund()`, `beta.multi.abund()` (both),
  and `beta.sample()` (statistically: it uses R's RNG), in
  `data-raw/expected/beta_partition.r`. Tables: `ex_counts`; a sparse
  40 × 300 random table with an empty sample, a sample nested in another,
  identical samples and a disjoint pair; and betapart's `ceram.s` /
  `ceram.n` example data (incidence; check that it is included in the
  package and small).
- **Cross-checks inside ecodive:** `total` equals the existing metric;
  turnover + nestedness = total exactly (by construction, but test the
  arithmetic order: compute nestedness as total − turnover, as betapart
  does).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- The eight metrics: the `ecodive_dist` heatmap, titled from `METRICS`.
- `plot(beta_partition(...))`: pairs as points of turnover (x) against
  nestedness (y), with dashed iso-lines of constant total (turnover +
  nestedness = c); `groups =` colours within-group pairs. Title: "Partition
  of Sorensen Dissimilarity" (from `index`).
- `plot(beta_partition_multi(...))` with draws: density or histogram of
  each component's draws (betapart's vignette plots these densities);
  without draws, a single stacked bar.


## Tests

- Fixtures at 1e-12 (the formulas are ratios of integers or sums).
- `turnover + nestedness == total`; `total` identical to the base metric;
  nested pair: turnover 0; disjoint pair: nestedness 0, turnover 1; equal
  richness: nestedness 0.
- Edge cases: both samples empty and one empty (match betapart / abdiv:
  NaN, ✱ confirm), a single site for `beta_partition_multi()` (error),
  `n_sites` > number of sites (error), `n_sites = 2`.
- `beta_div(x, 'betasim')` reaches `sorensen_turnover()`; aliases unique
  (README snippet in the test).
- `identical()` for `cpus` 1, 2, 4 (pairwise and resampled multi);
  dense = sparse; `pairs` respected.
- Smoke-test lines (one metric, both functions).


## Docs

The eight metrics in the beta metric list (`beta_div()`'s table and
`?beta_div`), a new `@family` "Partitioned metrics", and the two functions
under Beta Diversity in the reference index. Explain the concepts
(turnover vs nestedness; balanced variation vs abundance gradients), the
dependence of multiple-site values on the number of sites (hence
`n_sites`), the Baselga vs Podani/Carvalho debate, and "Differences from
betapart" (ecodive's RNG in `beta.sample()`-style draws; `norm`). Masking:
none (betapart uses dots; abdiv's `sorenson_*` are spelled differently).
NEWS.


## Later

- **Phylogenetic partition** (Leprieur et al. 2012; betapart
  `phylo.beta.pair()`; abdiv `unweighted_unifrac_turnover()`,
  `phylosor_turnover()`, ...; MIT): the same formulas with branch lengths
  (a = shared branch length, b and c unique branch length), on the UniFrac
  tree walk. Low extra effort once this is in; a separate small item.
- **Temporal** (`beta.temp()`) is just the pairwise values for matched
  samples: covered by `pairs`.


## References

- Baselga, A. (2010). Partitioning the turnover and nestedness components
  of beta diversity. *Global Ecology and Biogeography*, 19, 134–143.
  doi:10.1111/j.1466-8238.2009.00490.x
- Baselga, A. (2012). The relationship between species replacement,
  dissimilarity derived from nestedness, and nestedness. *Global Ecology
  and Biogeography*, 21, 1223–1232. doi:10.1111/j.1466-8238.2011.00756.x
- Baselga, A., & Orme, C. D. L. (2012). betapart: an R package for the
  study of beta diversity. *Methods in Ecology and Evolution*, 3, 808–812.
  doi:10.1111/j.2041-210X.2012.00224.x
- Baselga, A. (2013). Separating the two components of abundance-based
  dissimilarity: balanced changes in abundance vs. abundance gradients.
  *Methods in Ecology and Evolution*, 4, 552–557.
  doi:10.1111/2041-210X.12029
- Baselga, A. (2017). Partitioning abundance-based multiple-site
  dissimilarity into components: balanced variation in abundance and
  abundance gradients. *Methods in Ecology and Evolution*, 8, 799–808.
  doi:10.1111/2041-210X.12693
- Lennon, J. J., Koleff, P., Greenwood, J. J. D., & Gaston, K. J. (2001).
  The geographical structure of British bird distributions: diversity,
  spatial turnover and scale. *Journal of Animal Ecology*, 70, 966–979.
  doi:10.1046/j.0021-8790.2001.00563.x
- Carvalho, J. C., Cardoso, P., & Gomes, P. (2012). Determining the
  relative roles of species replacement and species richness differences
  in generating beta-diversity patterns. *Global Ecology and
  Biogeography*, 21, 760–771. doi:10.1111/j.1466-8238.2011.00694.x
- Podani, J., & Schmera, D. (2011). A new conceptual and methodological
  framework for exploring and explaining pattern in presence-absence data.
  *Oikos*, 120, 1625–1638. doi:10.1111/j.1600-0706.2011.19451.x
- Baselga, A., & Leprieur, F. (2015). Comparing methods to separate
  components of beta diversity. *Methods in Ecology and Evolution*, 6,
  1069–1079. doi:10.1111/2041-210X.12388
- Leprieur, F., Albouy, C., De Bortoli, J., et al. (2012). Quantifying
  phylogenetic beta diversity: distinguishing between "true" turnover of
  lineages and phylogenetic diversity gradients. *PLoS ONE*, 7, e42760.
  doi:10.1371/journal.pone.0042760
