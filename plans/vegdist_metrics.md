# `kulczynski()`, `mountford()`, `raup_crick()`, `binomial_deviance()`, `chao_jaccard()`, `cao()`, `mahalanobis_distance()`, `alt_gower()`, `chisq_distance()`: the vegdist metrics ecodive lacks

**Priority:** 79. **Tier:** formula; clean room against vegan (GPL-2,
output-only), with permissive references for some. Verified 2026-11-20:

- **vegan 2.7-8** (CRAN, published 2026-11-18): `License: GPL-2` (CRAN
  page). Work from `vegdist()`'s Rd (`man/vegdist.Rd`, read in full for
  this plan: its formula table and notes are quoted below in substance) and
  outputs. Never open `R/vegdist.R` or `src/vegdist.c`.
- **abdiv 0.2.0** (CRAN; `MIT + file LICENSE`, "COPYRIGHT HOLDER: Kyle
  Bittinger", 2020; checked): `weighted_kulczynski_second()` (vegan's
  quantitative Kulczynski), `binomial_deviance()` (Anderson & Millar 2004),
  `cy_dissimilarity()` (Cao; log10 and zeros as 0.1 per its Rd, so not
  vegan's scaling), `modified_mean_character_difference()` (its Rd points to
  vegan's `altGower`). May be read and ported.
- **SciPy** (BSD-3): `scipy.spatial.distance.mahalanobis` (given the inverse
  covariance) as a permissive cross-check; its `kulczynski1` is a binary
  similarity, not vegan's index.
- No permissive implementation was found for Mountford, Raup-Crick
  (analytic) or Chao's Jaccard.

Read [README.md](README.md) first.


## What ecodive has, and what is missing

ecodive's beta metrics (`METRICS`, `R/metrics.r`; kernels in
`src/beta_div.c`): aitchison, bhattacharyya, bray, canberra, chebyshev,
chord, clark, divergence, euclidean, gower, hamming, hellinger, horn,
jaccard (binary), jensen, jsd, lorentzian, manhattan, matusita, minkowski,
morisita, motyka, ochiai, psym_chisq, robust_aitchison, soergel, sorensen,
squared_chisq, squared_chord, squared_euclidean, topsoe, wave_hedges, the
five UniFracs and philr_distance.

`vegdist()`'s methods (vegan 2.7-8 Rd): manhattan, euclidean, canberra,
clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn,
mountford, raup, binomial, chao, cao, mahalanobis, chisq, chord,
hellinger, aitchison, robust.aitchison. Mapping:

| vegdist | ecodive |
|:--------|:--------|
| manhattan, euclidean, canberra, clark, bray, gower, morisita, horn, chord, hellinger, aitchison, robust.aitchison | same names (with `norm` for vegan's `binary = TRUE`) |
| jaccard (quantitative, 2B/(1 + B)) | **`soergel()`**: Ružička, 1 − Σ min / Σ max, which equals 2B/(1 + B). ecodive's `jaccard()` is vegan's `binary = TRUE`. Document the mapping in both Rd pages. |
| kulczynski, mountford, raup, binomial, chao, cao, mahalanobis, altGower | **missing**: this plan |
| chisq | **missing** (not ecodive's `squared_chisq` or `psym_chisq`, which are per-pair formulas; vegan's is Euclidean on chi-square-transformed data, the distance behind CA): this plan, see the question below |

`raup_crick()` here is vegdist's **analytic, presence/absence**
Raup-Crick (hypergeometric, equal species probabilities). It is **not**
vegan's `raupcrick()` (Chase et al. 2011's null-model version with
occupancy-weighted draws) and **not** Stegen's abundance-based RC_bray,
which is `raup_crick_bray()` in [assembly_processes.md](assembly_processes.md).
The Rd pages of all three say so.

**Demand** (retrieved 2026-11-20): Europe PMC microbiome full-text mentions
2022–26: "Kulczynski" 32 (the others are rarer). Baseline: Bray-Curtis
28,581. Value: completeness ("every `vegdist()` metric, threaded"), so
that users switching from vegan find what they used; low effort each.


## Names

Rule 6 and masking decide three names; the rest are vegan's.

| Function (METRICS id) | Display name | alt_ids | weighted | int_only | true_metric |
|:----------------------|:-------------|:--------|:---------|:---------|:------------|
| `kulczynski` | Kulczynski Dissimilarity | `kulczynskicody` | YES | NO | NO |
| `mountford` | Mountford Dissimilarity | | NO | NO | NO (Rd: non-metric) |
| `raup_crick` | Raup-Crick Dissimilarity | `raup raupcrick` | NO | NO | NO (Rd) |
| `binomial_deviance` | Binomial Deviance Dissimilarity | `binomial` | YES | NO | NO |
| `chao_jaccard` | Chao-Jaccard Dissimilarity | `chao` | YES | YES | NO |
| `cao` | Cao Dissimilarity (CYd) | `cyd cydissimilarity` | YES | YES | NO |
| `mahalanobis_distance` | Mahalanobis Distance | `mahalanobis` | YES | NO | YES |
| `alt_gower` | Alternative Gower Distance | `altgower modifiedmeancharacterdifference` | YES | NO | NO ✱ |
| `chisq_distance` | Chi-Square Distance | `chisq chisquare` | YES | NO | YES |

- `binomial_deviance`, not `binomial`: an exported `binomial()` would mask
  `stats::binomial()` (the GLM family) for every user. abdiv's name.
- `mahalanobis_distance`, not `mahalanobis`: would mask
  `stats::mahalanobis()` (different arguments). Pattern of
  `philr_distance()`.
- `chao_jaccard`, not `chao`: ecodive has the alpha metric `chao1`; vegan's
  index is Chao et al.'s Jaccard-type one (Chao's Sørensen type exists too,
  vegan's `chaodist()`).
- `chisq_distance`: next to `squared_chisq` and `psym_chisq`, a bare
  `chisq` would read as their square root, which it isn't. **Ask the user**
  whether to add it at all, and whether as this metric or as a `norm =
  'chisq'` choice for `euclidean()` (it needs column sums of the whole
  table, which `normalize()` doesn't do today). Recommended: the metric.
- `alt_gower`: snake_case of vegan's `altGower` (rule 5).
- `true_metric` ✱: set from a brute-force triangle search on random tables
  in the test; YES only if the literature proves it (Mahalanobis and
  chi-square are Euclidean, so YES).
- Check every alias with the README snippet (`chao` against `chao1` is in
  another `div`, but run it).


## Design

### Signatures

Like the existing metrics:

```r
kulczynski(counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus())
binomial_deviance(counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus())
alt_gower(counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus())
chisq_distance(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
mahalanobis_distance(counts, margin = 1L, norm = 'none', pseudocount = NULL, pairs = NULL, clust = NULL, cpus = n_cpus())
mountford(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
raup_crick(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
chao_jaccard(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
cao(counts, margin = 1L, pairs = NULL, clust = NULL, cpus = n_cpus())
```

- Presence/absence indices (`mountford`, `raup_crick`) take no `norm`, as
  `jaccard()`. Integer-count indices (`chao_jaccard`, `cao`) take no `norm`
  and call `assert_integer_counts()` (as `morisita()`): vegan's Rd says Cao
  "is intended for count (integer) data", and Chao's estimator counts
  singletons and doubletons. ✱ vegan doesn't refuse non-integers for Cao;
  **ask the user** whether ecodive should (recommended: yes for Chao,
  warn-only for Cao).
- `norm` where vegan users transform first (`decostand(x, "log")` then
  `altGower` is vegan's own example; Anderson et al. 2006). Note that
  vegan's `"log"` standardization is not ecodive's `norm = 'log'`
  (log-TSS); say so in `alt_gower()`'s Rd.

### Formulas (vegdist Rd, with x = x_ij, y = x_ik; sums over features)

- **Kulczynski:** 1 − ½ (Σ min(x, y) / Σ x + Σ min(x, y) / Σ y). Binary
  (`norm = 'binary'`): 1 − (J/A + J/B)/2, as the Rd states.
- **Binomial:** Σ_i [x log(x/n_i) + y log(y/n_i) − n_i log(½)] / n_i,
  n_i = x + y, over features with n_i > 0; 0 log 0 = 0. Binary: log(2)(A + B
  − 2J). No upper limit (Anderson & Millar 2004).
- **Cao (CYd):** (1/S) Σ_i [log(n_i/2) − (x log y + y log x)/n_i], natural
  logs, S = number of features present in either sample, n_i = x + y; zeros
  "replaced with arbitrary value 0.1" (Cao et al. 1997). ✱ Settle from
  outputs whether the 0.1 also enters n_i, and whether it replaces zeros
  before or after any scaling.
- **Mountford:** M = 1/α of Fisher's log-series, the positive root of
  exp(aM) + exp(bM) = 1 + exp((a + b − j)M) (a, b = richness of each
  sample, j = shared), solved by Newton's method in vegan; "if either a or
  b are equal to j ... the dissimilarity is 0"; "in the range 0 … log(2)".
  The Rd doesn't give the map from M to the returned dissimilarity, nor
  the j = 0 case. ✱ **Protocol first:** run vegan on crafted pairs over a
  grid of (a, b, j) and fit the returned values against candidates (M, M /
  log 2, 1 − M / log 2, with M from a high-precision root); also measure
  vegan's Newton tolerance (difference to the exact root). ecodive solves
  with a safeguarded Newton step (bracketed by bisection) to ~1e-14, so the
  match is to vegan's tolerance (probably ~1e-6 to 1e-10). The result
  depends only on (a, b, j): no table lookups inside the solve.
- **Raup-Crick (analytic):** "1 − prob(j), or based on the probability of
  observing at least j species in shared", from the hypergeometric
  distribution with S = **the number of features in the table** (so
  all-zero features change it; the Rd warns of this). ✱ The exact tail
  (1 − P(X ≥ j), 1 − P(X ≤ j), P(X ≤ j − 1) ...) is settled from outputs on
  a grid of (A, B, J, S), including j = 0 and identical samples (the Rd:
  "two communities with no shared species may have a dissimilarity slightly
  below one, and two identical communities may have dissimilarity slightly
  above zero"). In C: log-factorials tabulated once up to S (no `lgamma`
  in threads: glibc's sets a global); the hypergeometric pmf over the
  shorter tail, summed in a fixed order from the mode outward to keep
  precision, and the complementary tail used where the probability is near
  1 so `1 −` doesn't cancel. Match `phyper()` to 1e-10 relative.
- **Chao-Jaccard:** 1 − UV/(U + V − UV), U = C_j/N_j + (N_k − 1)/N_k ·
  a₁/(2a₂) · S_j/N_j, V symmetric (C_j: individuals in sample j of the
  species shared with k; N_j: total in j; a₁, a₂: shared species with one /
  two individuals in k; S_j: individuals in j of the shared species that
  have one individual in k; Chao et al. 2005). ✱ From outputs: a₂ = 0
  (Chao et al. 2005 use a₂ = 1 there, or the bias-corrected a₁(a₁ − 1)/2
  form), and whether U and V are capped at 1.
- **Mahalanobis:** "Euclidean distances of a matrix where columns are
  centred, have unit variance, and are uncorrelated" (Mardia et al. 1979):
  d² = (x − y)ᵀ Σ⁻¹ (x − y), Σ the sample covariance (n − 1). Compute
  without Σ⁻¹: Householder QR of the centred n × D table (as
  `src/feature_lm.c`; or Cholesky as `src/bioenv.c` does for its
  standardized distances), Y = √(n − 1) X_c R⁻¹, then the Euclidean kernel
  on Y. O(n D²) time, n × D dense memory; fine for D up to a few thousand.
  **Rank deficiency** (D ≥ n, or collinear features): vegan returns
  "nearly constant" noise and asks users to check; ecodive **errors** with
  "Mahalanobis distances need more samples than features and no
  collinear features; reduce the features first (e.g. pca())", using a
  relative tolerance on R's diagonal (✱ choose and document, e.g. 1e-9 ×
  max). A feature with zero variance is the same error.
- **altGower:** (1/NZ) Σ |x − y|, NZ = features not zero in both samples
  (Anderson et al. 2006); binary: (A + B − 2J)/(A + B − J). With `norm =
  'clr'` every feature is non-zero, so NZ = D; document.
- **Chi-square:** Euclidean distance after y_ij = √T · x_ij / (r_i √c_j)
  (row sums r, column sums c, total T; Legendre & Gallagher 2001;
  `decostand(x, "chi.square")`), over features with c_j > 0, so d²_jk = T
  Σ_i (x_ij/r_j − x_ik/r_k)² / c_i. ✱ Confirm the √T factor and the
  handling of empty rows on outputs.

### C

All in `src/beta_div.c` as new `BDIV_*` codes:

- Per-pair kernels on `FOREACH_OTU` / `WITH_ABJ`: kulczynski (Σ min and the
  two totals), binomial, cao (over the union of non-zero features), alt
  gower (Σ |x − y| and NZ), chao (shared-species sums and a₁, a₂),
  mountford and raup_crick (from A, B, J; Mountford's root solve per pair,
  Raup's tail sum per pair).
- Whole-table setups, as `gower_setup()` does for ranges: raup_crick
  (S = n_otus and the log-factorial table), chisq (column sums and T; then
  a weighted squared-difference kernel on proportions), mahalanobis
  (centring and QR in a setup step, serially or with the existing threaded
  helpers, then the Euclidean kernel on the dense Y; `C_beta_div()`'s
  ecomatrix is replaced by Y for that run).
- Everything inherits `FOREACH_PAIR`'s split: identical for any `cpus`,
  `pairs` respected, sparse tables never densified (except Mahalanobis'
  Y, which is dense by nature).

### Results

`ecodive_dist` via `dist_result()`, title from `METRICS`. vegan's
`maxdist` attribute is not added (ecodive's metrics don't carry it; the
Rd pages state each index's range instead).


## Reference implementations

- **Output-only fixtures:** vegan 2.7-8 `vegdist(x, m)` for every method
  above, with `binary = TRUE` where vegan defines it, in
  `data-raw/expected/vegdist_metrics.r`: `ex_counts`; vegan's `dune` and
  `varespec` (check they are small enough to store as text; else a
  subset); a sparse random 30 × 200 integer table with an empty sample,
  identical samples, disjoint samples and a nested pair; proportions (for
  the metrics that accept them); and the protocol grids for Mountford
  (a, b, j) and Raup-Crick (A, B, J, S), built as small binary tables.
- **Permissive cross-checks:** abdiv 0.2.0 for Kulczynski (= its
  `weighted_kulczynski_second`), binomial deviance, Cao (its log10 values
  × ln 10, if the 0.1 rule matches) and altGower; SciPy's `mahalanobis()`
  with `VI = inv(cov)` (`data-raw/expected/vegdist_metrics.py`).
- **Existing tests** (`tests/tinytests/test-beta_div.r`) keep reference
  values inline; add the new metrics to the same file or a new
  `test-vegdist_metrics.r` reading the stored fixtures.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots": every
metric is an `ecodive_dist`, drawn by `plot.ecodive_dist()` with its
`METRICS` display name. Nothing new.


## Tests

- Every metric against vegan at 1e-10 (Mountford: at the tolerance the
  protocol measured; record it in this plan).
- Identities: `kulczynski(norm = 'binary')` equals 1 − (J/A + J/B)/2 from
  plain R; `soergel(x)` equals `vegdist(x, "jaccard")` (add this fixture
  to document the mapping); `alt_gower(norm = 'binary')` equals vegan's
  binary form; `mahalanobis_distance()` equals `dist()` of the whitened
  table computed in plain R in the test (`chol(cov(x))`); `chisq_distance()`
  equals `dist()` of the plain-R chi-square transform.
- Raup-Crick: adding an all-zero feature changes the result (the
  documented dependence); a direct `phyper()` computation in the test
  matches.
- Errors: non-integers for `chao_jaccard()` (and `cao()` per the user's
  choice); Mahalanobis with D ≥ n; negative values for all.
- `identical()` for `cpus` 1, 2, 4; dense = sparse; `pairs` respected;
  `beta_div(x, 'binomial')`, `'chao'`, `'mahalanobis'`, `'raup'`,
  `'altGower'`, `'chisq'` reach the right functions.
- Smoke-test lines for each.


## Docs

One Rd page per metric (as now), each with its formula, base-R equivalent
(where short), range, and a "Differences from vegan" line when any ✱ choice
differs. Families: Kulczynski, binomial, Cao, Chao, altGower, Mahalanobis
and chi-square under "Abundance metrics"; Mountford and Raup-Crick under
"Presence/Absence metrics". Update `beta_div()`'s metric list and table.
In `jaccard()`'s and `soergel()`'s Rd, state the vegdist mapping. Masking
note: none (the names were chosen to avoid `stats`'). NEWS.


## References

- Kulczynski, S. (1927). Die Pflanzenassoziationen der Pieninen. *Bulletin
  International de l'Académie Polonaise des Sciences et des Lettres,
  Classe des Sciences Mathématiques et Naturelles, Série B*, Suppl. II,
  57–203.
- Faith, D. P., Minchin, P. R., & Belbin, L. (1987). Compositional
  dissimilarity as a robust measure of ecological distance. *Vegetatio*,
  69, 57–68. doi:10.1007/BF00038687 (quantitative Kulczynski)
- Mountford, M. D. (1962). An index of similarity and its application to
  classification problems. In P. W. Murphy (Ed.), *Progress in Soil
  Zoology* (pp. 43–50). Butterworths.
- Raup, D. M., & Crick, R. E. (1979). Measurement of faunal similarity in
  paleontology. *Journal of Paleontology*, 53, 1213–1227.
- Anderson, M. J., & Millar, R. B. (2004). Spatial variation and effects
  of habitat on temperate reef fish assemblages in northeastern New
  Zealand. *Journal of Experimental Marine Biology and Ecology*, 305,
  191–221. doi:10.1016/j.jembe.2003.12.011
- Chao, A., Chazdon, R. L., Colwell, R. K., & Shen, T.-J. (2005). A new
  statistical approach for assessing similarity of species composition
  with incidence and abundance data. *Ecology Letters*, 8, 148–159.
  doi:10.1111/j.1461-0248.2004.00707.x
- Cao, Y., Williams, W. P., & Bark, A. W. (1997). Similarity measure bias
  in river benthic Aufwuchs community analysis. *Water Environment
  Research*, 69, 95–106. doi:10.2175/106143097X125227
- Mahalanobis, P. C. (1936). On the generalized distance in statistics.
  *Proceedings of the National Institute of Sciences of India*, 2, 49–55.
- Mardia, K. V., Kent, J. T., & Bibby, J. M. (1979). *Multivariate
  Analysis*. Academic Press.
- Anderson, M. J., Ellingsen, K. E., & McArdle, B. H. (2006). Multivariate
  dispersion as a measure of beta diversity. *Ecology Letters*, 9,
  683–693. doi:10.1111/j.1461-0248.2006.00926.x
- Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
  transformations for ordination of species data. *Oecologia*, 129,
  271–280. doi:10.1007/s004420100716
- Chase, J. M., Kraft, N. J. B., Smith, K. G., Vellend, M., & Inouye, B. D.
  (2011). Using null models to disentangle variation in community
  dissimilarity from variation in α-diversity. *Ecosphere*, 2, art24.
  doi:10.1890/ES10-00117.1
