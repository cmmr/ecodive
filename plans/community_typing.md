# `enterotype()`, `valencia()`, `dmm()`: community types

**Priority:** 35 (`enterotype()`), 45 (`valencia()`), 66 (`dmm()`).
**Tier:**
- `enterotype()`: clean room. PAM is written from Kaufman & Rousseeuw;
  `cluster` and `clusterSim` are GPL (≥ 2), output-only.
- `valencia()`: port of VALENCIA. The code and reference centroids are
  **MIT**, verified 2026-11-20 from `LICENSE` at ravel-lab/VALENCIA commit
  `8559d45`: "Copyright (c) 2023 ravel-lab".
- `dmm()`: clean room from Holmes et al. (2012). Bioconductor
  DirichletMultinomial (LGPL-3) and mothur (GPL-3) are output-only.

Read [README.md](README.md) and [compositional.md](compositional.md) first.


## What it is, and demand

Community typing assigns each sample to one of a few discrete community
types. Citations below are from OpenAlex, retrieved 2026-11-19
(`research_notes/.../networks_typing.md`).

- **Enterotypes:** PAM on Jensen–Shannon distances, with the number of types
  chosen by the Calinski–Harabasz index. Arumugam et al. 2011 has 8,530
  citations, 3,182 of them in 2022–26. Costea et al. 2018, the consensus
  paper, has 1,277.
- **VALENCIA:** assigns vaginal samples to 13 reference sub-CSTs by
  nearest centroid. It has 499 citations, rising each year to 121 in 2026.
- **Dirichlet-multinomial mixtures (DMM):** a model-based clustering of
  counts, with the number of types chosen by the Laplace evidence. Holmes et
  al. 2012 has 976 citations, 506 of them in 2022–26; Ding & Schloss 2014
  has 975.

ecodive's `clust` argument only *orders* samples for heatmaps; none of
this exists yet (report, "Reconciliation").


## Names

- `enterotype()` names the procedure users ask for (rule 1). It is not
  `pam()`, which would mask `cluster::pam()` and name only one step. The PAM
  core is internal (`pam_fit()`); export it later only if asked.
- `valencia()` is the tool's own name (rule 1).
- `dmm()` is the field's acronym, used by Holmes et al., Ding & Schloss and
  mothur's documentation (rule 1). `dmn()` (DirichletMultinomial's function)
  is less recognised. `dirmult_mixture()` was considered, to sort with
  `dirmult_ttest()`, and rejected: that family is ALDEx2-style tests, not
  this model. **Masks** CRAN `dmm::dmm()` (dyadic mixed models; unrelated
  meaning; GPL). Say so under "Same names".
- No registry alias starts with `ent`, `val` or `dm` (checked 2026-11-20).
- **No wrapper.** The inputs differ (a table or distances, a table with
  VALENCIA's taxon names, raw counts), and README keeps table-returning
  methods as plain functions.
- New argument names (add to parameter_names.md): `n_types` and `centroids`.


## Shared output: `ecodive_types`

A `data.frame` of class `c('ecodive_<name>', 'ecodive_types',
'data.frame')`, **one row per sample**, row names = sample names:

| Column | Meaning |
|:-------|:--------|
| `sample` | sample name |
| `type` | factor: the assigned type (levels in type order) |
| method columns | `enterotype`: `silhouette`; `valencia`: `cst` (higher-order CST, factor), `score`; `dmm`: `probability` (posterior of the assigned type) |

Attributes:
- `types`: a data.frame, one row per type. Columns `type`, `n_samples`,
  `top_feature` (largest mean proportion; Arumugam's "driver"), plus
  per-method columns: `medoid` (enterotype), `weight` and `precision` (dmm:
  π and Σα), `reference_score_mean` / `_sd` (valencia).
- `profiles`: a types × features matrix. Mean proportions for `enterotype`;
  the Dirichlet means α/Σα for `dmm`; the centroids for `valencia`.
- `selection`: a data.frame, one row per candidate number of types, with
  `selected` marked. Columns `n_types`, `ch`, `silhouette` for enterotype;
  `n_types`, `laplace`, `aic`, `bic` for dmm; `NULL` for valencia.
- Method extras: `posterior` (dmm, n × K), `scores` (valencia, n × 13),
  `method`.

`res$type` can be used directly as `groups` in any test, matched by name.
`plot.ecodive_types(x, type = )`, base graphics, returning `x` invisibly:

- `'profiles'` (default): a heatmap or stacked bars of the top 10 features'
  profile per type.
- `'selection'`: the criterion against `n_types`, with the choice marked.
- `'scores'` (valencia): boxplots of `score` per sub-CST, with VALENCIA's
  reference mean ± sd. This ports `Valencia.py`'s plot: the values `ref_ave`
  and `ref_std` are MIT data in that file, and the colour scheme is
  `CST_color_scheme`.

A `print` method gives the counts per type and the selection.

---

## 1. `enterotype()` (priority 45)

```r
enterotype(x, n_types = 2:10, margin = 1L, cpus = n_cpus())
```

- `x`: counts (or proportions), giving `jensen(x, margin)` distances; or a
  `dist`, used as given, so users can try other distances. The paper says
  "**Jensen-Shannon distance**". That is ecodive's `jensen()` =
  √`jsd()` (one kernel, `BDIV_JSD`), not the divergence `jsd()`. The log
  base is irrelevant: PAM, CH and silhouette are all scale-invariant.
  Zeros need no pseudocount (0 log 0 = 0); some published scripts add 1e-6,
  so document that this changes results slightly. Users should aggregate to
  genus first (Arumugam et al. used genus abundances).
- `n_types`: the candidate numbers of types (a single value fixes it). Each
  must be between 2 and n − 1.

**Algorithm** (Arumugam et al. 2011, Methods: "Samples were clustered using
Jensen-Shannon distance and partitioning around medoid (PAM) clustering.
Optimal number of clusters was estimated using Calinski-Harabasz (CH)
index. We used the silhouette validation technique..."):

1. **PAM** (Kaufman & Rousseeuw 1990, ch. 2), clean room.
   - BUILD: the first medoid minimises \eqn{\sum_j d_{ij}}. Each later one
     maximises \eqn{\sum_j \max(D_j - d_{ij}, 0)}, with \eqn{D_j} the distance
     to the nearest chosen medoid.
   - SWAP: evaluate \eqn{\Delta T} for every (medoid, non-medoid) pair. Do
     the best swap if \eqn{\Delta T < 0}; repeat until none improves.
   - Ties go to the lowest index (check against `cluster::pam()`'s outputs).
   - Threaded over candidate swaps, with a fixed-order reduction to the
     minimum.
   - Types are numbered in the order of their medoids' positions. Report the
     BUILD and SWAP objectives.
2. **CH index** on distances: the pseudo-F,
   \eqn{T = \frac1n\sum_{i<j}d^2_{ij}},
   \eqn{W = \sum_k \frac1{n_k}\sum_{i<j\in k} d^2_{ij}},
   \eqn{CH = \frac{(T-W)/(k-1)}{W/(n-k)}}. For Euclidean d this is exactly
   Calinski & Harabasz (1974). It is PERMANOVA's F; reuse that code.
   `clusterSim::index.G1(..., centrotypes = "medoids")` may use a medoid
   variant. **Compare outputs**, and document any difference, or add the
   medoid form if the user wants parity.
3. **Silhouette** (Rousseeuw 1987): \eqn{s_i = (b_i - a_i)/\max(a_i,b_i)};
   0 for singletons. Report the per-sample value and the mean per k.
4. The chosen k is the CH maximum (the first on ties).

**Fixtures:** `cluster::pam(d, k)` (`id.med`, `clustering`, `objective`),
`cluster::silhouette()` and `clusterSim::index.G1()`, run on `jensen()` of a
genus table. Use a public table generated by data-raw, e.g. ecodive's own
`ex_counts`, or simulated three-cluster Dirichlet data. Medoids and
clustering must match exactly (when there are no ties), objectives and
silhouettes at 1e-12.

---

## 2. `valencia()` (priority 45; was 52 before the renumbering)

```r
valencia(counts, margin = 1L, centroids = valencia_centroids)
```

**Source read:** `Valencia.py` and `CST_centroids_012920.csv` at commit
`8559d45` (MIT; the README calls the file `CST_profiles_012920.csv`, but it
is `CST_centroids_...`). Don't read `Publication_materials/Cluster_testing/vegdist.R`:
it may derive from vegan (GPL), and its provenance is unclear.

**Algorithm (`Valencia.py`):**
1. Relative abundances: each taxon's count / the sample's `read_count`.
   VALENCIA takes `read_count` from an input column, which can exceed the
   taxa's sum. ecodive uses the row sum. Document this. **Ask the user**
   whether to add a per-sample totals argument.
2. Taxa are the **union** of the sample's and the centroids' columns,
   matched by exact name, with missing values 0. A sample taxon absent from
   the centroids adds \eqn{p^2} to the denominator.
3. The Yue–Clayton similarity to each of the 13 centroids:
   \eqn{\theta = \sum pq / (\sum (p-q)^2 + \sum pq)}. Compute it sparsely as
   \eqn{\sum pq/(\sum p^2 + \sum q^2 - \sum pq)}, which equals it to
   rounding (1e-12).
4. `type` (sub-CST) = the largest similarity, the first in VALENCIA's order
   on ties (pandas `idxmax`). The order is I-A, I-B, II, III-A, III-B, IV-A,
   IV-B, IV-C0–C4, V. `score` = that similarity.
5. `cst`: I-A/I-B → I; III-A/B → III; IV-C0–C4 → IV-C; the others as
   they are.
6. A sample with total 0 gives `NaN` scores and `NA` type (the Python gives
   NaN); warn.

**Data:** ship the 13 × 199 centroid matrix as an exported dataset
`valencia_centroids` (`data/valencia_centroids.rda`, built by
`data-raw/valencia_centroids.r` from the CSV at the pinned commit).
- The MIT notice goes in its Rd and in a `Files: data/valencia_centroids.rda
  R/valencia.r` stanza of `inst/COPYRIGHTS`. Add a `cph` person
  ("ravel-lab", VALENCIA) to `Authors@R`.
- Provenance: centroids of 13,160 profiles from 1,975 women (France et al.
  2020; the paper is CC-BY 4.0, and its article data are CC0). Only means are
  shipped, no individual data.
- `centroids` lets users pass revised centroids (e.g. with all *Prevotella*
  merged, which VALENCIA's README suggests).

**Taxon names** are VALENCIA's: species for key genera
(`Lactobacillus_crispatus`), `g_`/`f_` prefixes otherwise. Warn with the
list when any centroid's dominant taxon is absent from `counts`:
`Lactobacillus_crispatus`, `_iners`, `_gasseri`, `_jensenii`, `BVAB1`,
`Gardnerella_vaginalis`, `g_Streptococcus`, `g_Enterococcus`,
`g_Bifidobacterium`, `g_Staphylococcus`. Report the fraction of reads
matched in an attribute. Converting QIIME 2 names (`convert_qiime.py`, MIT)
is out of scope; mention it.

**Fixtures:** run `Valencia.py` (python3, pandas) on (a) 200 synthetic
samples drawn as Dirichlet around each centroid, plus unmatched taxa and an
all-zero row, and (b) a 50-sample subset of
`Publication_materials/Data_and_metadata/all_samples_taxonomic_composition_data.csv`
(MIT, same repository; keep the notice in data-raw). Similarities must
match at 1e-12; `subCST`, `CST` and `score` exactly.

---

## 3. `dmm()` (priority 55)

```r
dmm(counts, n_types = 1:7, margin = 1L, n_starts = 1L, max_iters = 1000L,
    tol = 1e-6, seed = 0L, cpus = n_cpus())
```

**Model** (Holmes, Harris & Quince 2012, eqs. 1–15, read from the paper's
equation figures). Sample i has counts \eqn{X_i} (total \eqn{J_i}); type k
has Dirichlet \eqn{\alpha_k = e^{\lambda_k}} (S features) and weight
\eqn{\pi_k}. The prior is iid \eqn{\mathrm{Gamma}(\eta, \nu)} (shape, rate)
on every \eqn{\alpha_{jk}}, which becomes \eqn{\propto \alpha^{\eta} e^{-\nu\alpha}}
in λ. **The values of η and ν are not in the text** (DirichletMultinomial
is believed to use η = ν = 0.1; unverified). Pin them by fitting the
fixtures under candidate values. Raw integer counts are required (warn
otherwise).

- **E-step** (eq. 12):
  \eqn{z_{ik} \propto \pi_k\, B(\alpha_k + X_i)/B(\alpha_k)}, on the log
  scale with `lgammafn`, and log-sum-exp per sample.
- **M-step:** \eqn{\pi_k = \frac1N\sum_i z_{ik}} (eq. 13). For each k,
  minimise
  \eqn{-\sum_i z_{ik}[\log B(\alpha_k+X_i) - \log B(\alpha_k)] +
  \nu\sum_j\alpha_{jk} - \eta\sum_j\lambda_{jk}}
  over \eqn{\lambda_k}, with gradient
  \eqn{-\alpha_{jk}\sum_i z_{ik}[\psi(\alpha_{jk}+X_{ij}) - \psi(\alpha_{jk})
  - \psi(A_k+J_i) + \psi(A_k)] + \nu\alpha_{jk} - \eta},
  where \eqn{A_k = \sum_j\alpha_{jk}}.
  - The paper used GSL's BFGS (GPL; don't read). Port **liblbfgs**
    instead (MIT, verified 2026-11-20 from `COPYING`: "Copyright (c) 1990
    Jorge Nocedal, (c) 2007-2010 Naoaki Okazaki"; github.com/chokkan/liblbfgs),
    with its notice in the C file and a `COPYRIGHTS` stanza.
  - Don't use R's `vmmin()`: it allocates through R and isn't
    thread-safe.
  - `digamma`/`trigamma` come from Rmath; check that they are thread-safe in
    the R versions supported.
- **Stop** when the change in the objective (eq. 11) is below `tol`
  (relative), or after `max_iters`.
- **Initialisation** isn't stated in the paper. DirichletMultinomial is
  believed to use soft k-means (unverified, and its source mustn't be
  read). ecodive uses k-means++ seeding on proportions (PCG) then Lloyd,
  with initial \eqn{\alpha_k} from the cluster's moment estimates, and
  repeats `n_starts` times, keeping the best posterior. It doesn't use
  `stats::kmeans()`, which uses R's RNG.
- **Model choice** (eqs. 14–15): Laplace,
  \eqn{-\log p(X|\mathcal H_K) \approx -\log p(X|\hat Q) - \log p(\hat Q) -
  \frac M2\log 2\pi + \frac12\log|H|}. H is block-diagonal per component in
  λ: S × S blocks with the paper's diagonal and off-diagonal second
  derivatives (the ψ, ψ₁ expressions of eq. 15, with
  \eqn{c_{ijk} = \alpha_{jk}+X_{ij}} and \eqn{C_{ik} = A_k + J_i}), plus the
  π terms \eqn{N/\pi_k^2}.
  - The paper's text doesn't define M, nor whether the multinomial constant
    \eqn{J_i!/\prod X_{ij}!} is included. Take both from DirichletMultinomial's
    `laplace()` / `AIC()` / `BIC()` outputs (constants that don't depend on K
    leave the choice unchanged, but parity needs them).
  - \eqn{\log|H|} comes from a Cholesky per block, O(K S³). Warn above S =
    2,000.
  - The chosen K minimises Laplace; AIC and BIC are reported too.
- Components are ordered by decreasing \eqn{\pi_k}, so labels are stable.
- **Threads:** the E-step is split over samples. The M-step objective and
  gradient are summed over samples in fixed blocks, then reduced in block
  order. Different `n_types` and starts are independent tasks. Results are
  identical for any `cpus`.

**Fixtures** (`data-raw/expected/dmm.r`): DirichletMultinomial
`dmn(count, k, seed)` for k = 1..4 on a simulated three-type table (D = 30,
n = 120) and on ecodive's `ex_counts` aggregated to 30 features. Save
`mixture()`, `fitted()`, `laplace()`, `AIC()`, `BIC()`.
- The optimisers differ, so match the fitted α and π at 1e-5 relative after
  ordering components, the posteriors at 1e-6 and Laplace at 1e-6.
- Assignments must match exactly when the posterior gap is > 1e-3.
- Run mothur's `get.communitytype` once as a second check, if available.

---

## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: `enterotype()` names its input, so it sets `display$title` with
  `set_display()` and `title_of("Enterotypes", <the distance's title>)`:
  "Enterotypes of Jensen-Shannon Distance" for counts (from the `jensen()`
  result it computed), or the given `dist`'s title. `valencia()` and
  `dmm()` take counts, which have no title, so their plots build the
  title from `method`: "VALENCIA Community State Types",
  "Dirichlet-Multinomial Mixtures".
- What the plots read, all stored when the result is made: `selection`
  (the per-`k` criterion values and the chosen `k`, for the
  model-selection plot), `profiles` (per-type mean composition, for the
  type profiles), `scores` (valencia), and the `type` factor's levels,
  a stable type order from which the colour key follows, so every plot of
  the same result agrees.
- The assignment can be passed as `groups` to tests and ordination plots:
  `res$type` is a plain `factor` with the type names as levels.

## Tests (all three)

- Fixtures as above.
- Hand-made cases:
  - PAM on a line of points with an obvious split;
  - CH equals the classical formula on Euclidean data;
  - silhouette against a direct base-R implementation;
  - VALENCIA on a sample identical to a centroid (score 1, that sub-CST);
  - a DMM on data drawn from a known 2-component mixture recovers it.
- Invariance: feature order; sample order (assignments permute with the
  samples); `dmm()` with the same `seed` gives the same result.
- `identical()` for `cpus` 1, 2 and 4 (PAM's swap search, DMM's E/M steps).
- Edge cases: n = 3; `n_types` out of range (error); duplicate samples
  (zero distances); an all-zero sample (error for `enterotype`/`dmm`, `NA`
  for `valencia`); a feature absent everywhere.
- `plot()` of each `type` runs (`pdf(NULL)`).
- Smoke-test lines for all three.


## Docs

- Three Rd pages, in the new pkgdown group "Compositional Data" →
  *Community types*. A short shared section explains the `ecodive_types`
  structure and how to use `type` as `groups`.
- Enterotypes: cite Costea et al. 2018's caution that enterotypes are a
  coarse, continuous stratification, and point to `dmm()` as the
  model-based alternative.
- VALENCIA: list the taxon naming rules.
- DMM: the cost (O(K S³) for the Laplace term).
- "Differences from ..." sections: cluster/clusterSim; Valencia.py's
  `read_count`; DirichletMultinomial's optimiser and initialisation.
- NEWS.


## References

- Arumugam, M., Raes, J., Pelletier, E., et al. (2011). Enterotypes of the
  human gut microbiome. *Nature*, 473, 174–180. doi:10.1038/nature09944
- Costea, P. I., Hildebrand, F., Arumugam, M., et al. (2018). Enterotypes in
  the landscape of gut microbial community composition. *Nature
  Microbiology*, 3, 8–16. doi:10.1038/s41564-017-0072-8
- Kaufman, L., & Rousseeuw, P. J. (1990). *Finding Groups in Data*. Wiley.
  doi:10.1002/9780470316801
- Calinski, T., & Harabasz, J. (1974). A dendrite method for cluster
  analysis. *Communications in Statistics*, 3, 1–27.
  doi:10.1080/03610927408827101
- Rousseeuw, P. J. (1987). Silhouettes. *J Comput Appl Math*, 20, 53–65.
  doi:10.1016/0377-0427(87)90125-7
- France, M. T., Ma, B., Gajer, P., et al. (2020). VALENCIA: a nearest
  centroid classification method for vaginal microbial communities based on
  composition. *Microbiome*, 8, 166. doi:10.1186/s40168-020-00934-6
- Yue, J. C., & Clayton, M. K. (2005). A similarity measure based on species
  proportions. *Commun Stat Theory Methods*, 34, 2123–2131.
  doi:10.1080/STA-200066418
- Ravel, J., et al. (2011). Vaginal microbiome of reproductive-age women.
  *PNAS*, 108, 4680–4687. doi:10.1073/pnas.1002611107
- Holmes, I., Harris, K., & Quince, C. (2012). Dirichlet multinomial
  mixtures: generative models for microbial metagenomics. *PLoS ONE*, 7,
  e30126. doi:10.1371/journal.pone.0030126
- Ding, T., & Schloss, P. D. (2014). Dynamics and associations of microbial
  community types across the human body. *Nature*, 509, 357–360.
  doi:10.1038/nature13178
- Liu, D. C., & Nocedal, J. (1989). On the limited memory BFGS method for
  large scale optimization. *Math Program*, 45, 503–528.
  doi:10.1007/BF01589116


## Decisions made while implementing (enterotype)

Item 35, 2026-11-20. `enterotype()` is in `R/enterotype.r` and `src/pam.c`;
the shared `ecodive_types` output (constructor, registry, `print()` and
`plot()`) is in `R/community_types.r`.

**Asked of the user up front:**
- **CH index: both, via an argument.** `ch = c('medoids', 'centroids')`,
  default `'medoids'` (parity with Arumugam et al.'s script, which calls
  `clusterSim::index.G1(t(data), cl, d = dist, centrotypes = "medoids")`);
  `'centroids'` is the distances' pseudo-F (PERMANOVA's F).
- **A `dist` input takes an optional `counts`**, used only for `profiles`
  and `top_feature` (and the medoid-form CH; below). Matched to the
  distances' labels by name, else by position. Without it, `profiles` is
  `NULL`, `top_feature` is `NA`, and `plot(type = 'profiles')` errors with
  a hint. Giving `counts` when `x` is already a count table is an error.
- Signature: `enterotype(x, n_types = 2:10, ch = c('medoids', 'centroids'),
  counts = NULL, margin = 1L, cpus = n_cpus())`.

**Settled from the references' outputs** (cluster 2.1.6, clusterSim
0.51.6, R 4.3.3; their source was never read):
- **BUILD's gain includes the candidate itself:** each later medoid
  maximises \eqn{\sum_j \max(D_j - d_{ij}, 0)} over every non-medoid
  \eqn{j}, *including* \eqn{j = i} (which gains \eqn{D_i}). Kaufman &
  Rousseeuw's wording ("another nonselected object j") reads as excluding
  it, but that version's BUILD objectives disagree with `pam()`'s (e.g.
  0.4018 vs 0.2824 on `ex_counts`, k = 2), and on varespec it then ends at
  different local optima. Including \eqn{i} matches every BUILD objective.
- **Ties** (pinned on 300 random tie-heavy integer grids, n = 10–30,
  k = 2–8, 2,100 fits):
  - BUILD: the **highest** index wins (first medoid and later ones). The
    plan's "lowest index" was a guess to check; the outputs say highest.
  - SWAP: the best \eqn{\Delta T}, the first in (lowest h, then lowest
    medoid) order; `h`-major and `i`-major orders were indistinguishable.
  - Assignment: a sample equally near two medoids joins the lower-indexed
    medoid (3,487/3,487 tied samples); **a medoid always keeps its own
    cluster**, even when another medoid is its duplicate (distance 0). The
    first version assigned it to the lower duplicate, leaving a cluster
    empty; 142 of the 2,100 fits caught that.
  - Result: medoids identical in 2,099/2,100 fits, partitions in
    2,100/2,100, objectives at 1e-12 in 2,099. The one miss (seed 235,
    k = 2): six samples tie for the first BUILD medoid, and `pam()` takes
    neither the highest nor the lowest; probably rounding in its own
    arithmetic. The plan only requires exact agreement without ties, so
    this is documented ("Differences from cluster and clusterSim") and
    not chased.
- **SWAP acceptance:** a swap is made when \eqn{\Delta T < -16
  \epsilon T} (relative), so rounding can't make it cycle; a cap of
  100 n swaps is a never-reached guard.
- **The medoid-form CH (`clusterSim`) is a mixed index.** With
  `centrotypes = "medoids"`, `index.G1(x, cl, d)` takes the medoids from
  `d` (each cluster's member with the smallest total within-cluster
  distance; \eqn{M} the same over all samples), but measures spread in
  `x`'s space: \eqn{T = \sum_i \lVert x_i - x_M\rVert^2},
  \eqn{W = \sum_i \lVert x_i - x_{m(i)}\rVert^2},
  \eqn{CH = ((T-W)/(k-1))/(W/(n-k))}. Found by fitting candidate formulas
  to its outputs: it returned NaN for `x` = zeros, changed with random `x`,
  and was unchanged when `x` was scaled. ecodive therefore uses the
  samples' **proportions** as `x` (what Arumugam et al. pass), from the
  count table or from `counts`; for a `dist` without `counts`, it uses the
  distances themselves for \eqn{e}, which equals `index.G1()` whenever the
  distances are Euclidean distances of its `x`. Matches to ≤ 1.1e-15
  relative.
- **The centroid form** is exactly `index.G1(x, cl)` for `d = dist(x)`
  (≤ 3.4e-15 on iris, k = 2–5), and equals `permanova()`'s F (tested). It
  is computed by a dedicated O(n²) loop in `src/pam.c` rather than through
  `C_permanova` (whose interface is built around permutations).
- **Numbering:** `cluster::pam()` numbers clusters by first appearance
  (1,000/1,000 random fits); ecodive numbers types by their medoids'
  positions, as the plan says. Documented.

**Other decisions:**
- **Type labels** are `"1"`, `"2"`, ... (factor levels in type order).
  `dmm()` can do the same; `valencia()`'s types are its sub-CST names.
- **Only a `dist` is read as distances;** a matrix is always a count
  table (a symmetric count table would otherwise be ambiguous).
- **`n_types`:** when omitted, `2:10` is cut to `2:min(10, n - 1)`; an
  explicit value outside 2..n−1 is an error. Fewer than 3 samples is an
  error. Values are sorted and de-duplicated.
- **Counts:** NA, negative values and a sample with total 0 are errors
  (the plan: error for an all-zero sample). No pseudocount.
- **Beyond the plan's columns:** `types$silhouette` (each type's mean
  silhouette, as `cluster`'s `clus.avg.widths`); `selection$objective`
  (PAM's mean distance to the nearest medoid, per k); attributes
  `objective = c(build =, swap =)` for the chosen k (the plan's "report
  the BUILD and SWAP objectives") and `ch` (which form made the `ch`
  column, a setting that qualifies a result, README rule 1). The C code
  also returns each sample's neighbouring cluster, not exposed.
- **Registry:** `COMMUNITY_TYPES` (id, display name, the `selection`
  criterion column, its label, and whether `max` or `min` wins) keyed by
  `method`, per README's "facts live in a registry". `print()` and the
  selection plot read it.
- **print():** a custom `print.ecodive_types()` (the plan's "counts per
  type and the selection"): the title; "n samples in k types, chosen by
  the <criterion> from <range>"; the `types` table; the `selection` table
  (when more than one k was tried); then the tibble-style preview of the
  samples, with the usual footer, all through `print_table_lines()`.
  Documented under `?print_ecodive`.
- **plot():** `type = c('profiles', 'selection')`. `'profiles'` draws
  stacked bars of the `n_features` (default 10) features with the
  largest mean proportion in any type, the rest pooled as "Other" (grey),
  bars labelled "type (n = …)", legend in a widened right margin, colours
  `plot_palette()`. `'selection'` plots the registry's criterion against
  `n_types`, circling the choice. `'scores'` is not added yet: valencia
  adds it to the `match.arg()` choices.
- **Threads:** BUILD (over candidates), SWAP (over non-medoids h, each
  computing every \eqn{T_{ih}} as an in-order sum over j) and silhouettes
  (over samples) use `run_parallel()`; the best candidate is chosen
  afterwards in index order. CH is a serial O(n²) loop. Results are
  `identical()` for cpus 1, 2, 4 (tests, n = 210; and all fixtures).
  Different `n_types` run one after another, each threaded.
- **Speed:** n = 2,000, 60 features, k = 2–6, from a `dist`: 1.07 s on one
  thread, 0.62 s on four, against 3.55 s for `cluster::pam()` +
  `silhouette()` + `index.G1()`.
- **Fixtures** (`data-raw/expected/enterotype.r`): `ex_counts` has only 4
  samples, so it is used for k = 2–3 only; vegan's `varespec` (24 samples)
  was added, with the simulated three-type Dirichlet table (60 × 30), a
  tie-heavy Manhattan grid with duplicated points, and iris (Euclidean,
  for the centroid form). Inputs are stored (8.9 KB); the tests recompute
  the distances.
- **Examples** use `iris[, 1:4]` (150 samples; types recover setosa vs the
  rest) because `ex_counts` is too small for `n_types = 2:10`.
- **Clean room:** no GPL source was read. The only GPL text seen was
  `cluster::pam()`'s R-level error message "Number of clusters 'k' must be
  in {1,2, .., n-1}; hence n >= 2", printed when a fixture script first
  asked for k = 4 of 4 samples; it was not used (ecodive's messages are
  its own).

**For `valencia()` and `dmm()`:** reuse `new_ecodive_types()`,
`types_profiles()`, `types_top_feature()`, `types_proportions()` and
`validate_n_types()` (its default-cutting rule), add a row to
`COMMUNITY_TYPES` (dmm: `criterion = 'laplace'`, `best = 'min'`;
valencia: no criterion, so `selection = NULL`, which `print()` and the
selection plot already handle), register `print`/`plot` through the
shared `ecodive_types` methods (no per-method S3 methods needed), and add
`'scores'` to `plot.ecodive_types()`. `print()` and `plot()` mark the
choice from `selection$selected`, never by recomputing a max or min, so
`dmm()`'s minimum-Laplace choice needs no change there (the `best` column
documents it).


## Outcome (from the implementing session's final report)

### `enterotype()`: item 35, done 2026-11-20

Commits `20f28e1` (work), `bdf2842` (tick); `582773b` (print locale fix).

- Tests: 18,326 at home, 337 in CRAN mode (after the locale fix).
- cluster 2.1.6 `pam()` / `silhouette()`, clusterSim 0.51.6: medoids and
  partitions identical; objectives 1.6e-16, silhouettes exact, CH 3.4e-15;
  2,099 of 2,100 tie-heavy fits identical (one six-way tie).
- Speed: n = 2,000, k = 2–6: 1.07 s (1 thread), 0.62 s (4); 3.55 s for
  cluster + clusterSim.
- Reusable for `valencia()` / `dmm()`: `new_ecodive_types()`,
  `COMMUNITY_TYPES`, shared `print` / `plot` (valencia needs a `'scores'`
  plot type).


## Decisions made while implementing (valencia)

Item 45, 2026-11-24. `valencia()` is in `R/valencia.r` and `src/valencia.c`;
its centroids are `data/valencia_centroids.rda` (built by
`data-raw/valencia_centroids.r`); the shared `ecodive_types` code in
`R/community_types.r` gained a registry row, a `'scores'` plot and a small
`print()` change.

**Asked of the user up front** (all four answered with the recommended
option):
- **A per-sample totals argument: yes.** `totals = NULL` (each sample's
  sum); otherwise one value per sample, matched by name when both are
  named, else by position, for exact parity with Valencia.py's
  `read_count`. Non-finite, negative, wrong-length, or 0 for a sample with
  counts is an error; a total below the sample's sum (beyond 1e-6
  relative, so rounded proportions with `totals = 1` pass) is allowed, as
  VALENCIA allows it, with a warning.
- **A sample sharing no taxa with any centroid** scores 0 against all 13;
  ecodive follows Valencia.py (pandas `idxmax`: the first, I-A, score 0),
  with a warning naming the samples.
- **The matched fraction is a per-sample column, `matched`** (the share of
  the sample's proportions in taxa that the centroids name; NA for a
  sample of total 0), not an attribute.
- **`type` always has all 13 levels** (the centroids' row names, in their
  order), and `types` always has one row per centroid (n_samples 0 where
  none), parallel to `profiles`; `cst` likewise has all its levels (I, II,
  III, IV-A, IV-B, IV-C, V). Tests already drop unused levels of `groups`
  (checked with `permanova()`).

**Other decisions:**
- **Signature:** `valencia(counts, centroids = valencia_centroids,
  totals = NULL, margin = 1L, cpus = n_cpus())`. `cpus` is added for
  `run_parallel()`; `margin` and `cpus` go last, as in the other count
  functions (the plan had `counts, margin, centroids`).
- **The default `centroids`** shows as `valencia_centroids` in the usage,
  but a lazy-loaded dataset isn't visible from the namespace, so the body
  takes `getExportedValue('ecodive', 'valencia_centroids')` when the
  argument is missing (the default promise is never forced). A user's own
  object of that name is used only when passed. `utils::globalVariables()`
  declares the name, for R CMD check's "no visible binding" NOTE.
- **`centroids`** may be a matrix (sub-CSTs in the row names, taxa in the
  column names, both unique; finite, non-negative; no all-zero centroid) or
  a data frame read from VALENCIA's file (its `sub_CST` column names the
  rows). Any number of centroids and any names are allowed: `cst` maps
  VALENCIA's names (I-A/I-B → I, III-A/B → III, IV-C0–C4 → IV-C) and keeps
  others; reference scores are `NA` and colours come from `plot_palette()`
  for other names.
- **Similarity:** \eqn{\sum pq/(\sum p^2 + \sum q^2 - \sum pq)} over each
  sample's non-zero taxa, with each centroid's \eqn{\sum q^2} precomputed.
  Its denominator is at least half of \eqn{\sum p^2 + \sum q^2}, so there
  is no cancellation; the largest difference from Valencia.py over all
  fixtures is 2.8e-15. A sample identical to a centroid (with `totals = 1`)
  scores exactly 1.
- **Assignment:** the largest similarity, first on ties, skipping NaN (as
  `idxmax` does); `score` is that similarity. A sample of total 0 has NaN
  scores, `NA` type and cst, NaN score, `NA` matched, and a warning
  (`print()` reports "k without a type").
- **Taxon-name warning:** instead of the plan's hard-coded list, the key
  taxa are each centroid's most abundant taxon where it makes up at least
  a quarter; for VALENCIA's centroids that is exactly the plan's 10 (IV-C0,
  an even community led by g_Corynebacterium_1 at 0.099, is left out). The
  warning lists the key taxa that `counts` has no feature for. Matching is
  by exact name only (no case folding or renaming); duplicated or missing
  feature names are errors.
- **A named vector** is one sample.
- **Result:** columns `sample`, `type`, `cst`, `score`, `matched`;
  attributes `types` (`type`, `n_samples`, `top_feature`, `cst`,
  `score_mean`, `reference_score_mean`, `reference_score_sd`), `profiles`
  (the centroids, as given), no `selection` (absent: `new_ecodive_types()`
  with `NULL`), `scores` (n × 13, sample and sub-CST names), `method`.
  `score_mean` (each sub-CST's mean score) is beyond the plan's columns.
  No `display` is set: counts have no title, so `print()` and `plot()` take
  "VALENCIA Community State Types" from `COMMUNITY_TYPES`.
- **Registry:** `COMMUNITY_TYPES` gained `valencia` with `criterion`,
  `criterion_label` and `best` all `NA`; `plot(type = 'selection')` then
  errors ("no criterion"), as it should.
- **print():** the header says "(k with samples)" when some types are
  empty, and "; k without a type" for NA types; unchanged for enterotype.
- **plot(type = 'scores'):** boxplots (no notches: R warns when notches
  pass the hinges, common with few samples) of `score` for each sub-CST with
  samples, filled with VALENCIA's `CST_color_scheme`, beside a diamond and
  bar for VALENCIA's reference mean ± sd (`ref_ave`, `ref_std`, stored in
  `types` when the result is made), y from 0 to 1, with a legend. Valencia.py
  fills its boxes black (`boxprops = dict(color = 'k')` with
  `patch_artist`); not copied. `'profiles'` shrinks the bar labels when
  there are more than 8 types (13 for valencia).
- **Data provenance and licence:** `data-raw/valencia_centroids.r`
  downloads the CSV at commit 8559d45 through jsDelivr and checks its MD5
  (40029ff8...), keeping the names as they are (`check.names = FALSE`:
  'BVAB2/BVAB3', 'g_Escherichia.Shigella'). The MIT notice is in
  `?valencia_centroids`, in `src/valencia.c` and `R/valencia.r`, in
  `data-raw/expected/valencia/LICENSE`, and in a `Files:` stanza of
  `inst/COPYRIGHTS` covering the C and R files, the `.rda` and the test
  fixture (which holds 50 publication samples). `Authors@R` gained
  "ravel-lab" (ctb, cph).
- **Fixtures** (`data-raw/expected/valencia.r` and `.py`): Valencia.py
  (pandas 2.2.3, numpy 1.26.4; pandas 3 breaks its positional `Series`
  indexing) on (a) 195 Dirichlet-multinomial samples around the centroids,
  with two unmatched taxa in some, read counts above the sum in a quarter,
  an all-zero sample, an unmatched-only sample and three mixtures; (b) 30 of
  them against centroids with I-B = I-A and III-A = III-B (exact ties, the
  first wins), passed as a data frame with `sub_CST`; (c) 50 samples of the
  publication data, about 4 per sub-CST. Similarities are written as hex
  floats. Valencia.py reproduces VALENCIA's published assignments for all
  50; its similarities differ from the published ones because the
  publication table names BVAB2/BVAB3 as g_Fastidiosipila (renamed,
  `valencia()` gives the published similarities to 6.5e-8). The
  publication case thus also covers a centroid taxon absent from `counts`
  and a sample taxon absent from the centroids.
- **Threads:** samples are split over threads (`run_parallel()`); each
  sample's sums run in a fixed order, so results are `identical()` for any
  `cpus` (tested 1–4; n = 200 and 234, above `run_parallel()`'s 100-task
  threshold).
- **Speed:** 100,000 samples × 196 taxa: 0.51 s on one thread from a
  dense matrix (mostly the dense-to-sparse import), 0.22 s on four from a
  dgCMatrix; Valencia.py takes 7.2 s for 200 samples.

