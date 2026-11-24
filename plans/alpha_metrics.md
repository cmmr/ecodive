# `pielou()`, `goods_coverage()`, `hill()` and other alpha metrics: ports from scikit-bio

**Priority:** 37. **Tier:** port (scikit-bio, `skbio.diversity.alpha`,
BSD-3; formulas only, each a few lines). Read [README.md](README.md) first.


## What it is

ecodive has 14 alpha metrics. scikit-bio (which QIIME 2's `diversity
alpha` runs) has 37. The missing ones that matter are short formulas over
one sample's counts, and they slot into the existing per-sample C switch
in `src/alpha_div.c` and the `METRICS` registry:

- **Pielou's evenness** \eqn{J = H / \ln S}: one of the four alpha outputs
  of QIIME 2's `core-metrics` (`evenness_vector`), with observed features,
  Shannon and Faith's PD, which ecodive already has. Pielou (1966): 3,163
  citations 2022–26 (OpenAlex); "Pielou" in 2,001 microbiome full texts
  (Europe PMC, 2022–26). Retrieved 2026-11-20.
- **Good's coverage** \eqn{1 - F_1 / N}: the standard check that sequencing
  depth was enough. Good (1953): 766 citations; "Good's coverage" in 1,236
  microbiome full texts (retrieved 2026-11-20).
- **Hill numbers** \eqn{{}^qD} of order \eqn{q}, and the Rényi and Tsallis
  entropies they map to: the unifying "effective number of species" of Hill
  (1973) and Jost (2006). Hill (1973): 2,727 citations; "Hill numbers" in
  319 microbiome full texts (retrieved 2026-11-20). They are also the
  quantity [hill_extrapolation.md](hill_extrapolation.md) interpolates and
  extrapolates.
- **Simpson's evenness, Simpson's dominance, Heip's and McIntosh's
  evenness, Strong's dominance, the Gini index:** cheap completions of
  scikit-bio's evenness and dominance groups, all offered by QIIME 2.

Baselines for scale: Bray-Curtis 28,581 microbiome full-text mentions;
LEfSe 19,755 (retrieved 2026-11-20).


## The metrics

All follow scikit-bio 0.7.5-dev (`skbio/diversity/alpha/_base.py`,
`_gini.py`; main branch read 2026-11-20). \eqn{x_i} are a sample's non-zero
counts, \eqn{N = \sum x_i}, \eqn{S} the number of non-zero features,
\eqn{p_i = x_i / N}, \eqn{H = -\sum p_i \ln p_i}, \eqn{D = \sum p_i^2}.

| id | Display name (`METRICS$name`) | Formula | `int_only` | skbio | Family |
|:---|:-----|:--------|:-----------|:------|:-------|
| `pielou` | Pielou's Evenness | \eqn{H / \ln S}; 1 when \eqn{S = 1} | NO | `pielou_e` | Evenness |
| `heip` | Heip's Evenness | \eqn{(e^H - 1)/(S - 1)}; 1 when \eqn{S = 1} | NO | `heip_e` | Evenness |
| `simpson_evenness` | Simpson's Evenness | \eqn{1 / (S D)} | NO | `simpson_e` | Evenness |
| `mcintosh_evenness` | McIntosh's Evenness | \eqn{\sqrt{\sum x_i^2} / \sqrt{(N - S + 1)^2 + S - 1}} | YES | `mcintosh_e` | Evenness |
| `goods_coverage` | Good's Coverage | \eqn{1 - F_1 / N} (\eqn{F_1}: singletons) | YES | `goods_coverage` | Coverage |
| `hill` | Hill Number | \eqn{(\sum p_i^q)^{1/(1-q)}}; \eqn{e^H} at \eqn{q = 1}; \eqn{1/\max p_i} at \eqn{q = \infty} | NO | `hill` | Diversity |
| `renyi` | Rényi Entropy | \eqn{\ln(\sum p_i^q)/(1 - q)}; \eqn{\ln S} at 0, \eqn{H} at 1, \eqn{-\ln \max p_i} at \eqn{\infty}; 0 when \eqn{S = 1} | NO | `renyi` | Diversity |
| `tsallis` | Tsallis Entropy | \eqn{(1 - \sum p_i^q)/(q - 1)}; \eqn{H} at 1, 0 at \eqn{\infty}; 0 when \eqn{S = 1} | NO | `tsallis` | Diversity |
| `dominance` | Simpson's Dominance Index | \eqn{D} | NO | `dominance`, `simpson_d` | Dominance |
| `strong` | Strong's Dominance Index | \eqn{\max_i (c_i / N - i / S)}, \eqn{c_i} the cumulative sum of the counts sorted in decreasing order | NO | `strong` | Dominance |
| `gini_index` | Gini Index | \eqn{\max(0, 1 - 2B)}, \eqn{B = \frac{1}{K}\sum_{j=1}^{K} L_j}, \eqn{L_j} the Lorenz curve (cumulative share) of **all \eqn{K} features of the table, zeros included**, sorted increasing ("rectangles") | NO | `gini_index` | Dominance |

`alt_ids` (letters only; checked with README's snippet): `pielou`:
`pieloue pielouj`; `heip`: `heipe`; `simpson_evenness`: `simpsone`;
`mcintosh_evenness`: `mcintoshe`; `goods_coverage`: `goods good`; `hill`:
`hillnumber`; `dominance`: `simpsond`; `gini_index`: none (see
**Names**). Add scikit-bio's other names for existing metrics too:
`observed`: `sobs observedfeatures`; `inv_simpson`: `enspie`; `berger`:
`bergerparkerd`; `mcintosh`: `mcintoshd`; `brillouin`: `brillouind`;
`fisher`: `fisheralpha`.

Edge cases, as scikit-bio: an empty sample is `NA` for every metric (the
existing `FOREACH_SAMPLE` already does this; scikit-bio gives `nan`).
\eqn{S = 1}: `pielou` and `heip` 1, `renyi` and `tsallis` 0 (scikit-bio's
explicit branches; vegan's `diversity()/log(specnumber())` gives `NaN` for
Pielou, so say so in the docs).


## Left out, and why

- `robbins` (\eqn{F_1/N}): exactly 1 − `goods_coverage`.
- `esty_ci`, `chao1_ci`, `lladser_ci`: intervals, two values per sample,
  which don't fit the one-value-per-sample interface of `alpha_div()`.
  Esty's interval could later be an extra output of `goods_coverage()`.
- `lladser_pe`: a random estimator drawn with NumPy's RNG, so it can't be
  matched; very rarely used.
- `michaelis_menten_fit`: a Powell optimisation over a rarefaction curve;
  superseded by [hill_extrapolation.md](hill_extrapolation.md)'s
  asymptotic estimators and rarely used.
- `kempton_taylor_q`: rarely used, two extra arguments (the quantiles), and
  scikit-bio's own docs flag its quantile rule as departing from Magurran's.
- `osd`, `singles`, `doubles`: counts, not diversity
  (`rowSums(counts == 1)`).
- `sobs`, `enspie`, `simpson_d`, `observed_features`: aliases of existing or
  new metrics (above).
- `phydiv` (generalised phylogenetic diversity: rooted/unrooted, abundance
  weighting `weight` in [0, 1]): a family of its own, with tree arguments to
  design; `faith()` covers the common case. A later item if asked for.
- A `base` for the logarithm (scikit-bio's `base` in `shannon`, `pielou_e`,
  `renyi`): ecodive's `shannon()` uses natural logs only; keep that, and say
  in the docs to divide by `log(b)`. Pielou's \eqn{J} doesn't depend on it.
- `finite = TRUE` (Simpson's index without replacement,
  \eqn{\sum x_i(x_i - 1)/(N(N-1))}): an option in scikit-bio only; the
  unbiased estimator appears in [hill_extrapolation.md](hill_extrapolation.md)
  as the \eqn{q = 2} asymptote. Revisit if asked.


## Names

- README rule 1: each metric is its usual short name, as the existing ones
  (`berger` for scikit-bio's `berger_parker_d`). Evenness versions of an
  existing index name are spelled out (rule 4): `simpson_evenness` and
  `mcintosh_evenness`, because `simpson` and `mcintosh` already mean the
  Gini-Simpson and McIntosh D indices. scikit-bio's `_e` / `_d` names are
  aliases.
- **`gini_index`, not `gini`** (**ask the user**): today `alpha_div(x,
  'gini')` prefix-matches the Gini-Simpson index (its display name's letters
  begin `gini`). An id `gini` would silently change what that call returns.
  With `gini_index`, `'gini'` becomes ambiguous (an error listing both),
  which is safe but breaks that prefix. The alternative is to keep `'gini'`
  meaning Gini-Simpson by adding it as `simpson`'s alias, and match the Gini
  index only as `gini_index`.
- `hill`, `renyi`, `tsallis` take the order **`q`** (new in
  parameter_names.md §1: "the order of a Hill number or of Rényi's or
  Tsallis's entropy, a single number ≥ 0, `Inf` allowed"). It is the
  standard symbol (Hill 1973, Jost 2006, Chao et al. 2014, iNEXT, hillR);
  scikit-bio calls it `order`, vegan `scales`.
- Same names in other packages: `vegan::renyi()` (a profile over `scales`,
  with `hill = TRUE` for Hill numbers); scikit-bio's `hill`, `renyi`,
  `tsallis`, `strong`, `dominance`. Say so in each Rd page.


## Design

### Signatures

```r
pielou(counts, margin = 1L, cpus = n_cpus())
heip(counts, margin = 1L, cpus = n_cpus())
simpson_evenness(counts, margin = 1L, cpus = n_cpus())
mcintosh_evenness(counts, margin = 1L, cpus = n_cpus())
goods_coverage(counts, margin = 1L, cpus = n_cpus())
hill(counts, q = 2, margin = 1L, cpus = n_cpus())
renyi(counts, q = 2, margin = 1L, cpus = n_cpus())
tsallis(counts, q = 2, margin = 1L, cpus = n_cpus())
dominance(counts, margin = 1L, cpus = n_cpus())
strong(counts, margin = 1L, cpus = n_cpus())
gini_index(counts, margin = 1L, cpus = n_cpus())

alpha_div(
  counts, metric, norm = 'percent', cutoff = 10L, digits = 3L,
  tree = NULL, margin = 1L, q = 2, cpus = n_cpus() )
```

- **`q` default (ask the user):** 2, scikit-bio's default for all three,
  so the fixtures compare defaults and "Differences from scikit-bio" stays
  short. The alternative is 1 (Jost 2006: the order that weighs features
  exactly by their frequency). Either default duplicates an existing metric
  (`inv_simpson` or `exp(shannon)`); the docs should tell users to give `q`.
- **`alpha_div()` gains `q`**, because it passes its formals to the metric
  with `mget(metric$params, environment())`. Put it **before `cpus`** (the
  glossary's "cpus always last") and after the released arguments, so no
  positional call to the CRAN version changes meaning. `alpha_div_rarefied()`
  gains `q` the same way (its `args` list in `R/rarefied.r`), and so do
  [rarefaction_curves.md](rarefaction_curves.md)'s functions.
- **`validate_q()`** in `R/validate.r` (`validate_args()` calls
  `validate_<arg>` for every formal): a single non-NA number ≥ 0; `Inf`
  allowed. Message: "`q` must be a single number of at least 0 (Inf is
  allowed)."
- ~~**Several metrics at once:** `alpha_div(counts, metric = c('shannon',
  'pielou', 'observed'))` returns an n × M matrix.~~ **Dropped** (the
  user, 2026-11-24, after it was first implemented in `9f1add2`): the
  return type must not depend on the input's length. `alpha_div()` and
  `alpha_div_rarefied()` take one metric and return a named vector; build a
  matrix with `sapply(metrics, alpha_div, counts = counts)`. See decision 4.

### C (`src/alpha_div.c`)

- New ids `ADIV_DOMINANCE` ... `ADIV_TSALLIS` (15 on, in both
  `R/alpha_div.r` and the C `#define`s), each a `FOREACH_SAMPLE` worker in
  the existing switch. Each sample is summed by one thread, in storage
  order, so results are identical for any `cpus`.
- `norm`: the proportion-based metrics are called with `norm = 'percent'`
  (as `shannon()`), so the C sees \eqn{p_i}; `goods_coverage` and
  `mcintosh_evenness` with `'none'` plus `assert_integer_counts()` (as
  `chao1()`), so the C sees counts.
- `q` reaches C through `sexp_extra` (a double), as `cutoff` and `digits`
  do. Branches exactly as scikit-bio: `q == 1`, `q == 0` (Rényi), `q ==
  Inf`, else the power sum. Compute \eqn{e^H} as `exp(-sum p log p)`, not
  scikit-bio's `prod(p^-p)`; the two differ by rounding only.
- `strong` and `gini_index` sort each sample's values: a per-thread buffer
  of the largest `nnz` (setup function, as `ace_setup()`), `qsort()` per
  sample. `gini_index` also needs \eqn{K} = `em->n_otus`, the number of
  features in the table: the zeros sort first and add nothing to the
  cumulative sums, so \eqn{B = \frac{1}{K}\sum_{j} L_j} over the non-zero
  values' positions. **Document** that the Gini index depends on how many
  (absent) features the table has, so filtering features changes it.
- No new `.Call` entry point; `C_alpha_div` keeps 6 arguments.

### Result

Unchanged shape: a named numeric vector, one value per sample, `NA` for an
empty sample; the matrix described above for several metrics. Like the
existing metrics, no attributes.


## Reference implementations

- **Port from:** scikit-bio 0.7.5-dev (github.com/scikit-bio/scikit-bio,
  main, read 2026-11-20), BSD-3, "Copyright (c) 2013--, scikit-bio
  development team" (`LICENSE.txt` read). The formulas above are each one
  line in `_base.py` / `_gini.py`; the edge-case branches (\eqn{S = 0, 1};
  \eqn{q = 0, 1, \infty}) are scikit-bio's. That is a port of behaviour,
  so put the BSD-3 notice in `src/alpha_div.c`'s header for the new
  section and add a `Files: src/alpha_div.c R/alpha_div.r` stanza to
  `inst/COPYRIGHTS` (the existing scikit-bio stanzas name other files).
  scikit-bio is already a `cph` in `Authors@R`: add "alpha diversity
  metrics" to its `comment`.
- **Fixtures:** `data-raw/expected/alpha_metrics.py` (scikit-bio from PyPI;
  record the version), writing text under `data-raw/expected/alpha_metrics/`:
  every new metric on `ex_counts` (written out by the `.r` side), on edge
  vectors (one feature; all singletons; a sample of zeros; ties for
  `strong`; a table with extra all-zero features for `gini_index`), and
  `hill` / `renyi` / `tsallis` at \eqn{q \in \{0, 0.5, 1, 2, 3, \infty\}}.
  `data-raw/expected/alpha_metrics.r` reads them into
  `tests/tinytests/expected/alpha_metrics.rds`.
- **Cross-checks (no fixture needed):** vegan's `renyi(x, scales, hill =
  TRUE)` and `diversity()` (GPL-2, 2.7-8; output-only) agree with `hill` /
  `renyi` and `pielou` where both are defined; hillR 0.5.2 (MIT;
  `hill_taxa()`) agrees with `hill`. Worth one comparison in the fixture
  script, not a dependency.


## Plotting

None: the results are plain vectors, as for the existing metrics.


## Tests (`tests/tinytests/test-alpha_div.r`, extended)

- Every new metric equals the scikit-bio fixture at 1e-12 (relative), dense
  input, `margin = 1` and `2`, `Matrix` and `slam` input.
- Identities, exactly or at 1e-14: `hill(q = 0) == observed`, `hill(q = 1)
  == exp(shannon)`, `hill(q = 2) == inv_simpson`, `hill(q = Inf) ==
  1/berger`, `renyi == log(hill)`, `tsallis(q = 2) == simpson`,
  `dominance == 1 - simpson`, `simpson_evenness == inv_simpson /
  observed`, `pielou == shannon / log(observed)`, `goods_coverage == 1 -
  rowSums(x == 1) / rowSums(x)`.
- Proportions give the same values as counts for the proportion-based
  metrics; `goods_coverage` and `mcintosh_evenness` reject proportions with
  the existing integer message.
- Edge cases: empty sample `NA`; \eqn{S = 1}; all singletons
  (`goods_coverage` 0); `gini_index` changes when all-zero features are
  added, as scikit-bio's does.
- `q`: `-1`, `NA`, `c(1, 2)`, `'a'` are errors with the message above;
  `Inf` works.
- `alpha_div()`: every new id and alias matches; `alpha_div(x, 'hill', q =
  0.5)` passes `q`; several metrics give the matrix whose columns are
  `identical()` to the single calls; a duplicated metric errors;
  `alpha_div_rarefied(x, 'hill', q = 1)` passes `q`. The `'gini'` behaviour
  the user chooses.
- `identical()` for `cpus` 1, 2 and 4.
- Smoke: one line per new export.


## Docs

- One Rd page per metric, as the existing ones (`@inherit documentation`,
  `@inherit adiv_percent_normalized` or `adiv_assert_integer`, a **Base R
  Equivalent** block, `@family`). New families **Evenness metrics** and
  **Coverage metrics**; `hill`, `renyi`, `tsallis` in Diversity, `dominance`,
  `strong`, `gini_index` in Dominance. Add the two subtitles to
  `pkgdown/_pkgdown.yml`'s Alpha Diversity section.
- `alpha_div()`: the metric list in `@param metric`, the "Integer Count
  Requirements" lists, `@param q`, and the several-metrics return value.
- `hill()`'s page: the Hill-number profile (\eqn{q = 0, 1, 2} are richness,
  \eqn{e^H} and inverse Simpson), and that observed Hill numbers depend on
  depth: see `rarefaction_curve()` and `hill_extrapolation()`.
- `goods_coverage()`'s page: it is the reference sample's coverage
  estimate of Good (1953); Chao & Jost (2012) give a less biased one, which
  `hill_extrapolation()` reports.
- "Differences from scikit-bio" (on `alpha_div()`'s page): natural logs only
  (no `base`); no `finite`; Gini by rectangles only; the `q` default if it
  differs.
- `vignettes/adiv.Rmd`: add the new metrics to its table; NEWS; NAMESPACE
  (11 `export()` lines); parameter_names.md §1: `q`.


## Open questions for the user

1. `gini_index` id and what `'gini'` should match (above).
2. `q` default: 2 (scikit-bio) or 1.
3. Not in this plan, noticed while reading: `chao1()` uses the classic
   \eqn{S + F_1^2/(2F_2)}, which is `Inf` when \eqn{F_2 = 0} (and `NaN` when
   \eqn{F_1 = F_2 = 0}); scikit-bio, vegan's `estimateR()` and iNEXT default
   to the bias-corrected \eqn{S + \frac{N-1}{N}\frac{F_1(F_1 - 1)}{2(F_2 +
   1)}}. Add `chao1(bias_corrected = )`? It would need a `bias_corrected`
   formal on `alpha_div()` too.


## References

- Pielou, E. C. (1966). The measurement of diversity in different types of biological collections. *Journal of Theoretical Biology*, 13, 131–144. doi:10.1016/0022-5193(66)90013-0
- Good, I. J. (1953). The population frequencies of species and the estimation of population parameters. *Biometrika*, 40, 237–264. doi:10.1093/biomet/40.3-4.237
- Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. *Ecology*, 54, 427–432. doi:10.2307/1934352
- Jost, L. (2006). Entropy and diversity. *Oikos*, 113, 363–375. doi:10.1111/j.2006.0030-1299.14714.x
- Rényi, A. (1961). On measures of entropy and information. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 547–561.
- Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics. *Journal of Statistical Physics*, 52, 479–487. doi:10.1007/BF01016429
- Heip, C. (1974). A new index measuring evenness. *Journal of the Marine Biological Association of the UK*, 54, 555–557. doi:10.1017/S0025315400022736
- Heip, C., & Engels, P. (1974). Comparing species diversity and evenness indices. *Journal of the Marine Biological Association of the UK*, 54, 559–563. doi:10.1017/S0025315400022748 (McIntosh's evenness)
- Smith, B., & Wilson, J. B. (1996). A consumer's guide to evenness indices. *Oikos*, 76, 70–82. doi:10.2307/3545749 (Simpson's evenness)
- Strong, W. L. (2002). Assessing species abundance unevenness within and between plant communities. *Community Ecology*, 3, 237–246. doi:10.1556/ComEc.3.2002.2.9
- Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688. doi:10.1038/163688a0
- Gini, C. (1912). *Variabilità e mutabilità*. Cuppini, Bologna.
- Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. *Ecology*, 93, 2533–2547. doi:10.1890/11-1952.1


## Decisions made while implementing

Asked upfront (2026-11-24), with the user's answers:

1. **Gini index id: `gini`**, and `alpha_div(x, 'gini')` now means the Gini
   index (it used to prefix-match the Gini-Simpson index). Display name
   "Gini Index", so `'gini_index'` (scikit-bio's name) matches it through
   the name's letters; the Gini-Simpson index stays `'simpson'` (or
   `'Gini-Simpson'`). NEWS says so under "Changed matching".
2. **`q` default 2**, as scikit-bio (so `hill()` equals `inv_simpson()` by
   default; the docs tell users to give `q`).
3. **`chao1(bias_corrected = TRUE)` by default** (open question 3), with a
   `bias_corrected` formal on `alpha_div()` and `alpha_div_rarefied()` too
   (after `q`, before `cpus`). The formula follows **scikit-bio exactly**:
   \eqn{S + F_1(F_1 - 1)/(2(F_2 + 1))}, *without* the \eqn{(N-1)/N} factor
   that this plan's question attributed to scikit-bio (scikit-bio 0.7.4's
   `_chao1.py` has none; vegan's `estimateR()` and iNEXT do). The docs say
   so. `bias_corrected = FALSE` is scikit-bio's too: the classic
   \eqn{S + F_1^2/(2F_2)} only when \eqn{F_1 > 0} and \eqn{F_2 > 0}, else the
   bias-corrected value, so it is never `Inf` or `NaN` (the old `chao1()`
   gave both). This changes released results; NEWS says so under "Changed
   results".
4. ~~**`alpha_div_rarefied()` takes several metrics**~~. First built (with
   several metrics in `alpha_div()` too) in `9f1add2`, then **reverted at
   the user's request** (2026-11-24): "I don't like mixing return object
   types based on input length. It can break things unexpectedly."
   `alpha_div()` and `alpha_div_rarefied()` again take exactly one metric
   (`match_metric()`'s "`metric` must be length 1" error otherwise) and
   always return a named vector. `q` and `bias_corrected` stay.

Made while implementing:

5. **Display name "Renyi Entropy"** (ASCII, as the registry's existing
   "Dice-Sorensen Dissimilarity"): with "Rényi", `list_metrics()` printed
   `R<U+00E9>nyi` in a C locale, and the letters-only matcher dropped the
   é (`'rnyientropy'`). The Rd pages and vignette spell it Rényi.
6. **`strong()` and `gini()` run on raw values** (`norm = 'none'`) and divide
   cumulative sums by the sample's total in C, in scikit-bio's order of
   operations; they still accept proportions (`int_only` NO). The
   proportion-based metrics use `norm = 'percent'`; `goods_coverage()` and
   `mcintosh_evenness()` use `'none'` plus `assert_integer_counts()`.
7. **Heip's evenness uses `exp(H)`** (as the plan said), not scikit-bio's
   `prod(p^-p)`: they differ by up to 1.4e-14 relative, the largest
   difference from scikit-bio.
8. **Gini docs:** by rectangles, one feature of \eqn{K} gives
   \eqn{1 - 2/K} (not \eqn{1 - 1/K}); an even sample gives 0 after the
   `max(0, .)` clamp. The page says the index depends on the table's
   absent features.
9. **Empty samples are `NA` for every metric**, as before; scikit-bio gives
   `nan`, except `0` for `chao1()` and `gini_index()`. Listed in
   `alpha_div()`'s "Differences from scikit-bio".
10. ~~**Several metrics**~~ (reverted, see 4): the internal
    `match_metrics()` and `multi_metric_matrix()` helpers are gone, and
    `alpha_div()`'s body is as before item 37. Its docs and `adiv.Rmd` show
    `sapply(metrics, alpha_div, counts = counts)` instead.
11. **`validate_q()`** coerces to double, so `q = 1L` works;
    `validate_bias_corrected()` takes only `TRUE` or `FALSE`.
12. **Aliases:** as planned, plus none for `gini` (its name gives
    `giniindex`), `strong` or `tsallis`. All checked with README's snippet
    against METRICS, ORDINATIONS, COUNT_ORDS, DIST_TESTS,
    ZERO_REPLACEMENTS, DA_TESTS and FEATURE_ASSOCS: no collisions. Short
    prefixes that matched one metric and now match several: `'mc'`
    (mcintosh, mcintosh_evenness), `'g'`, `'h'`, `'sim'`; new single
    matches: `'p'` pielou, `'r'` renyi, `'t'` tsallis, `'d'` dominance,
    `'e'` inv_simpson (via `enspie`).
13. **Exports** are in NAMESPACE's alpha-metric block (after `squares`),
    not at the end, and the NEWS entries follow the rarefied-diversity
    entry, to keep clear of the concurrent item-39 session's appends.
14. **Fixtures:** five tables (`ex_counts`; nine edge rows; the same with 12
    absent features; a 30 × 40 sparse random table; a 6 × 25 deep table
    with counts around 1e5) and 2,146 scikit-bio values, including the
    existing metrics under scikit-bio's names (a check of the aliases).
    vegan 2.7.8's `renyi()` (with and without `hill = TRUE`),
    `diversity() / log(specnumber())` and `1 - diversity('simpson')` are
    stored in the same `.rds` and tested. hillR was not run: it would pull
    a large dependency tree into the scratch library, and vegan's
    `renyi(hill = TRUE)` already checks the Hill numbers.
15. **Smoke tests** now require every alpha metric to be finite on the
    smoke table (chao1 no longer gives `NaN`), and call each new export
    once, plus `alpha_div_rarefied()` with `q`.
16. **Vignettes:** `adiv.Rmd` gains Evenness, Hill numbers, Coverage and
    "Several metrics at once" (an `sapply()` over `alpha_div()`) sections
    and the new formulas, and its
    `chao1()` outputs are the bias-corrected values; `ecodive.Rmd`'s and the
    README's `list_metrics()` outputs are regenerated.

**For items 38 and 48:** `alpha_div(counts, metric, q =, bias_corrected =)`
and `alpha_div_rarefied()` take **one** metric and return a named vector
(after the revert in decision 4). Several metrics mean one call each, e.g.
`sapply(metrics, alpha_div, counts = counts)`; a function that offers
several metrics should keep a fixed return type too. The internal
`validate_q()` and `validate_bias_corrected()` (in `R/validate.r`, reading
from the caller's environment) are there to reuse.


## Outcome (from the implementing session's final report)

**Item 37, done 2026-11-24.** Commits `9f1add2` (work), `4ea15be` (tick
and this outcome), and a follow-up that reverted several metrics at once
in `alpha_div()` and `alpha_div_rarefied()` at the user's request
(decision 4).

- Tests (`-O2` install, `LANG=C.UTF-8`): 22,736 at home (baseline 20,773),
  388 in CRAN mode (baseline 364), all passing. `R CMD check --no-manual
  --ignore-vignettes`: 1 NOTE (installed size, as before).
  After the revert (one metric per call): 22,714 at home and 387 in CRAN
  mode, all passing; `R CMD check` again 1 NOTE (installed size).
- scikit-bio 0.7.4: 2,146 values on five tables, every new metric, both
  `chao1()` forms, `hill` / `renyi` / `tsallis` at q = 0, 0.5, 1, 2, 3,
  Inf, and the existing metrics under scikit-bio's names. Largest relative
  difference 1.4e-14 (`heip`, `exp(H)` against `prod(p^-p)`); `chao1`,
  `goods_coverage`, `mcintosh_evenness`, `strong` and every q = 0 and
  q = Inf value are exactly equal. Dense, `margin = 2`, `Matrix` and
  `slam` input agree; results are identical for `cpus` 1, 2, 4 and 7.
- vegan 2.7.8 (output only): `renyi()` with and without `hill = TRUE`,
  Pielou's J and Simpson's dominance agree within 2.8e-14 (Rényi at
  q = 2); vegan's Pielou is `NaN` for one feature, ecodive's 1 (as
  scikit-bio).
- Speed (2,000 × 2,000 dense table, half zeros, 4 cores): 0.06–0.17 s
  per metric including reading the table, against 0.24 s for vegan's
  `renyi(hill = TRUE)` at one order and 0.35 s for its Pielou's J.
- Open: `chao1()` differs from vegan's `estimateR()` and iNEXT by their
  \eqn{(N-1)/N} factor (documented); `esty_ci`, `chao1_ci`, `phydiv` and
  `finite` remain left out, as planned.
