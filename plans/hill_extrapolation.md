# `hill_extrapolation()` and `hill_standardize()`: interpolated and extrapolated Hill numbers (iNEXT-style)

**Priority:** 64. **Tier:** **clean room**. iNEXT 3.0.2 (CRAN, 2025-07-30)
is **GPL ≥ 3** (`DESCRIPTION` read 2026-11-20). Work only from the papers
(Chao et al. 2014 main text and Appendix G; Chao & Jost 2012; Chao, Wang &
Jost 2013; Hsieh, Ma & Chao 2016), iNEXT's documentation (Rd reference
manual and vignette) and its **outputs**. Do not open iNEXT's `R/`
directory, nor the R code in the paper's Ecological Archives supplement
(M084-003, `Main_for_iNEXT.r`, which states no licence). Steps marked ✱
are not fully pinned down by the papers and must be settled from iNEXT's
outputs first. Read [README.md](README.md) first, and
[rarefaction_curves.md](rarefaction_curves.md), which shares the coverage
formulas.


## What it is

Observed diversity grows with sequencing depth, so samples of unequal depth
can't be compared directly. Chao et al. (2014) give, for the Hill numbers
of order \eqn{q = 0, 1, 2} (richness, \eqn{e^H}, inverse Simpson):

- **rarefaction** (interpolation) to any smaller depth \eqn{m < n}, by
  nearly unbiased closed forms;
- **extrapolation** to a larger depth \eqn{n + m^*}, up to about double
  the sample (reliable for \eqn{q = 0} only that far; for \eqn{q = 1, 2}
  further), toward each order's asymptotic estimate;
- **sample coverage** at every depth (Chao & Jost 2012), so samples can be
  compared at **equal completeness** rather than equal depth;
- **bootstrap confidence intervals** (Appendix G), unconditional on the
  reference sample.

iNEXT (Hsieh, Ma & Chao 2016) implements this and is widely cited: 3,689
citations 2022–26, Chao et al. (2014) 3,087 (OpenAlex); "iNEXT" in 263
microbiome full texts, "coverage-based rarefaction" / "sample coverage" in
375 (Europe PMC; all retrieved 2026-11-20). Microbiome use is held back by
cost: iNEXT is single-threaded R, slow at sequencing depths of 10⁴–10⁵
reads per sample over hundreds of samples. Threaded C with deterministic
bootstraps is the case for ecodive.

The per-sample standardized values (`hill_standardize()`) are what users
then test: they feed [alpha_test.md](alpha_test.md) directly.


## Permissive implementations: none usable (checked 2026-11-20)

- PyPI has no `inext`, `pyinext`, `iNEXT`, `hill-numbers` or similar
  package. scikit-bio has point Hill numbers only.
- GitHub `davidlevybooth/pynext` (2022): no licence, and its README says it
  is a "Python port of iNEXT.ind function from iNEXT", i.e. of GPL code.
  Excluded; only the README was read.
- GitHub `KilianMaire/hillrep` 0.4.0 (2026, MIT-labelled): its README says
  "Every estimator is a direct port of the corresponding iNEXT kernel", and
  its CHANGELOG "a port of iNEXT's general-q kernels (`Dq` / `Dq_2nd`)". A
  translation of GPL code can't be relicensed MIT, so treat it as GPL:
  **don't open its `src/`**, and don't use it for fixtures (iNEXT's own
  outputs serve). Only its README and CHANGELOG were read.
- So: clean room.


## Names

- **`hill_extrapolation()`** for the curves (iNEXT's `iNEXT()`) and
  **`hill_standardize()`** for diversity at a common depth or coverage
  (iNEXT's `estimateD()`). README rule 2: the base is `hill` (the Hill
  number metric of [alpha_metrics.md](alpha_metrics.md)), the modifiers
  say what is done to it, and the three sort together. `hill_extrapolation`
  names what is new compared with `rarefaction_curve()` (which also
  interpolates); `standardize` is the papers' word ("standardizing samples
  by completeness") and a verb like `rarefy`, `normalize`.
- Rejected: `inext()` (a package's brand, not the method's name; ecodive
  names methods), `hill_curve()` (too close to `rarefaction_curve()`),
  `estimate_d()`.
- Arguments: `q` (vector here; the glossary's `q`), `depths` (iNEXT's
  `size`), `max_depth` (`endpoint`), `n_depths` (`knots`; `n_<unit>`),
  `n_boot` (`nboot`), `conf_level` (`conf`), `coverage` / `depth` (the
  level of `estimateD()`, by `base`), `seed`, `margin`, `cpus`.
- Same names in other packages: none known. Same purpose: iNEXT's
  `iNEXT()`, `estimateD()`, `ChaoRichness()`, `ChaoShannon()`,
  `ChaoSimpson()`.


## Design

### Signatures

```r
hill_extrapolation(
  counts, q = c(0, 1, 2), depths = NULL, max_depth = NULL,
  n_depths = 40L, n_boot = 50L, conf_level = 0.95, margin = 1L,
  seed = 0L, cpus = n_cpus() )

hill_standardize(
  counts, q = c(0, 1, 2), coverage = NULL, depth = NULL, n_boot = 50L,
  conf_level = 0.95, margin = 1L, seed = 0L, cpus = n_cpus() )
```

- Each sample (row, with `margin = 1`) is one assemblage: iNEXT's
  `datatype = 'abundance'` on each column of a species × sites table.
  Integer counts only (`assert_integer_counts()`).
- `q`: values in \eqn{\{0, 1, 2\}} in the first version (Chao et al. 2014,
  Table 1); see **Later steps**.
- `max_depth`: `NULL` = twice each sample's size (iNEXT's `endpoint`); a
  number = the same endpoint for all. `depths`: explicit depths (iNEXT's
  `size`); else `n_depths` knots from 1 to the endpoint, about half on
  each side of \eqn{n} (✱ the exact knot rule, from iNEXT's `m` column for
  several \eqn{n} and endpoints), always including \eqn{n}.
- `n_boot = 0`: no intervals (iNEXT's `se = FALSE`, `nboot = 0`).
- `hill_standardize()`: give `coverage` (a level in (0, 1]) **or** `depth`
  (a whole number), not both. With neither, the coverage default of
  `estimateD(base = 'coverage')`: the smallest coverage among the samples
  extrapolated to twice their size (its docs). **Ask the user** whether the
  default should instead be `estimateD()`'s own default, `base = 'size'`;
  Chao & Jost (2012) recommend coverage.

### Estimators (abundance data; Chao et al. 2014)

For one sample: \eqn{n} reads, counts \eqn{X_i}, frequency counts
\eqn{f_k}, \eqn{S_{obs}}. All binomial ratios by running products, as in
[rarefaction_curves.md](rarefaction_curves.md) (shared helpers).
(Item 38 wrote them, in `src/rarefaction_curve.c`, declared in
`src/ecodive.h`: `rc_freq()` (frequency counts), `rc_absence()` (the
absence probabilities \eqn{\binom{N-k}{d}/\binom{N}{d}}),
`rc_coverage()` (\eqn{C_m}, \eqn{m < n}) and `rc_coverage_full()`
(Eq. 12 with Appendix G). Heck et al.'s variance there, `heck_variance()`,
avoids the cancellation of \eqn{a_{k+l} - a_k a_l} with a recurrence: see
its comment.)

- **Undetected richness** (Chao1, bias-adjusted; Appendix G):
  \eqn{\hat f_0 = \frac{n-1}{n}\frac{f_1^2}{2 f_2}} if \eqn{f_2 > 0}, else
  \eqn{\frac{n-1}{n}\frac{f_1(f_1 - 1)}{2}}.
- **Reference coverage** (Eq. 12; Appendix G for \eqn{f_2 = 0}):
  \eqn{\hat C(n) = 1 - \frac{f_1}{n}\frac{(n-1)f_1}{(n-1)f_1 + 2f_2}}, or
  with \eqn{(n-1)(f_1 - 1)/((n-1)(f_1-1) + 2)} when \eqn{f_2 = 0}.
- **Interpolation, \eqn{m < n}** (Table 1):
  - \eqn{q = 0}: \eqn{S_{obs} - \sum_i \binom{n - X_i}{m}/\binom{n}{m}}
    (= rarefaction_curves.md's expected richness; test the identity).
  - \eqn{q = 1}: \eqn{\exp\left[\sum_{k=1}^{m} -\frac{k}{m}\log\frac{k}{m}
    \hat f_k(m)\right]}, \eqn{\hat f_k(m) = \sum_{X_i \ge k}
    \binom{X_i}{k}\binom{n - X_i}{m - k}/\binom{n}{m}} (Eq. 8). Compute it
    per abundance class as the expectation of \eqn{-\frac{K}{m}\log
    \frac{K}{m}} over the hypergeometric \eqn{K}, with the pmf built from
    its mode outward by ratios and summed in a fixed order; cost O(\eqn{X})
    per class, O(\eqn{n}) per knot. Stop each tail once its terms fall
    below \eqn{10^{-17}} of the running sum (✱ check this keeps 1e-10
    agreement).
  - \eqn{q = 2}: \eqn{1 / \left[\frac{1}{m} + \frac{m-1}{m}\sum_i
    \frac{X_i(X_i - 1)}{n(n-1)}\right]} (Eq. 11b).
  - Coverage: \eqn{\hat C(m) = 1 - \sum_i \frac{X_i}{n}\binom{n -
    X_i}{m}/\binom{n-1}{m}}.
- **At \eqn{m = n}:** the observed Hill numbers (Eq. 7) and \eqn{\hat C(n)}.
- **Extrapolation, \eqn{n + m^*}** (Table 1):
  - \eqn{q = 0}: \eqn{S_{obs} + \hat f_0\left[1 - \left(1 - \frac{f_1}{n
    \hat f_0 + f_1}\right)^{m^*}\right]} (\eqn{S_{obs}} when \eqn{f_1 = 0}).
  - \eqn{q = 1}: \eqn{\exp\left[\frac{n}{n + m^*} H_{obs} + \frac{m^*}{n +
    m^*}\hat H\right]} (Eq. 10c), with \eqn{H_{obs} = -\sum \frac{X_i}{n}
    \log\frac{X_i}{n}}. ✱ iNEXT 3.x may use a later form; compare its
    extrapolated \eqn{q = 1} values with Eq. 10c first, and if they differ,
    with the discovery-rate extrapolation of Chao & Jost (2015). Record
    which.
  - \eqn{q = 2}: \eqn{1 / \left[\frac{1}{n + m^*} + \frac{n + m^* -
    1}{n + m^*}\sum_i \frac{X_i(X_i-1)}{n(n-1)}\right]} (Eq. 11c).
  - Coverage: \eqn{1 - \frac{f_1}{n}\left[\frac{(n-1)f_1}{(n-1)f_1 +
    2f_2}\right]^{m^* + 1}} (✱ the \eqn{f_2 = 0} variant, as for \eqn{\hat
    C(n)}).
- **Asymptotes** (iNEXT's `AsyEst`; Table 1 footnote):
  - \eqn{q = 0}: \eqn{S_{obs} + \hat f_0} (`ChaoRichness()`).
  - \eqn{q = 1}: \eqn{\exp(\hat H)}, Chao, Wang & Jost (2013), Eq. 10b:
    \eqn{\hat H = \sum_{k=1}^{n-1}\frac{1}{k}\sum_{1 \le X_i \le n - k}
    \frac{X_i}{n}\binom{n - X_i}{k}/\binom{n-1}{k} + \frac{f_1}{n}(1 -
    A)^{1-n}\left[-\log A - \sum_{r=1}^{n-1}\frac{1}{r}(1-A)^r\right]},
    \eqn{A = 2f_2/((n-1)f_1 + 2f_2)}; the second term is 0 when \eqn{f_1 =
    0}. ✱ \eqn{A} when \eqn{f_2 = 0} (likely \eqn{2/((n-1)(f_1-1)+2)}) and
    when \eqn{A = 1}. The first sum by abundance class with running
    products: O(\eqn{n K}).
  - \eqn{q = 2}: \eqn{n(n-1)/\sum_i X_i(X_i - 1)} (Good 1953;
    `ChaoSimpson(transform = TRUE)`).
- **Coverage-based curves** (iNEXT's `coverage_based`, and
  `hill_standardize(coverage =)`): for a level \eqn{C}, the depth with
  \eqn{\hat C(m) = C}. Extrapolation inverts in closed form, \eqn{m^* =
  \log(n(1 - C)/f_1)/\log B - 1}, \eqn{B = (n-1)f_1/((n-1)f_1 + 2f_2)}
  (non-integer). ✱ Interpolation: whether iNEXT takes a real \eqn{m}
  (root of a continuous extension) or interpolates between integer depths,
  and which coverage grid `coverage_based` uses; settle from outputs.

### Bootstrap (Appendix G, read)

Per sample:

1. The bootstrap assemblage: \eqn{S_{obs} + \lceil \hat f_0 \rceil}
   species. Detected species get \eqn{\hat p_i = \frac{X_i}{n}\left[1 -
   \hat\lambda\left(1 - \frac{X_i}{n}\right)^n\right]}, \eqn{\hat\lambda =
   (1 - \hat C(n)) / \sum_i \frac{X_i}{n}(1 - \frac{X_i}{n})^n}; each of
   the \eqn{\lceil \hat f_0 \rceil} undetected species gets \eqn{(1 - \hat
   C(n))/\lceil \hat f_0 \rceil}. (✱ \eqn{f_1 = 0}: \eqn{\hat C = 1},
   nothing undetected.)
2. Draw a multinomial sample of size \eqn{n} from it (✱ Appendix G says "a
   random sample of m individuals"; every estimator at every knot needs a
   reference sample, so read it as size \eqn{n} and confirm from iNEXT's
   interval widths), and recompute **every** output (each knot, each
   \eqn{q}, coverages, asymptotes) from it.
3. `n_boot` replicates; SE = their SD; interval = estimate ±
   \eqn{z_{(1+\text{conf})/2}} SE (Appendix G: ±1.96 SE). ✱ Whether iNEXT
   truncates the lower limit (at 0, or at \eqn{S_{obs}} for richness).
   Coverage-based intervals: ✱ whether each replicate re-inverts the
   coverage level (iNEXT 3.0: "confidence limits at fixed sample
   coverage").

**RNG and threads:** replicate \eqn{b} of sample \eqn{i} uses its own PCG
stream, `rng_seed(&r, seed, (uint64_t) i * n_boot + b)` (`rng_t` in
`src/rng_dist.c`), so work can be split over samples and replicates in any
way and give identical results for any `cpus`. The SD sums replicates in
replicate order. R's RNG is never touched; iNEXT's intervals agree only
within Monte Carlo error (a different random stream; README allows it).

**Multinomial draws:** conditional binomials, \eqn{O(S)} per replicate
rather than \eqn{O(n)}. Port the binomial sampler from NumPy's
`numpy/random/src/distributions/distributions.c` (`random_binomial`:
inversion for small \eqn{np}, BTPE otherwise; NumPy is BSD-3, "Copyright
(c) 2005-2026, NumPy Developers", `LICENSE.txt` read 2026-11-20; check the
file's own header when porting) onto `rng_unif()`. Notices: the BSD-3 text
in the new C file's header, an `inst/COPYRIGHTS` stanza, and "The NumPy
developers" as `cph` in `Authors@R`. (Fallback with no port: an alias
table and \eqn{n} categorical draws; slower for deep samples.)

### C (`src/hill_extrapolation.c`) and cost

- One task per (sample, replicate block); the point estimates are
  replicate 0. Per task: frequency counts, the running-product tables, and
  every knot's estimates.
- The \eqn{q = 1} interpolation dominates: O(\eqn{n}) per knot per
  replicate. 300 samples × 50,000 reads × 40 knots × 51 fits ≈ 3 × 10¹¹
  hypergeometric terms before tail cutting; profile, and say in the docs
  that `n_boot` and `n_depths` set the cost.
- Memory: O(\eqn{n + K}) per thread.

### Results

`hill_extrapolation()`: a list of class `ecodive_hill_extrapolation` (a
short `print()` method, as `print.ecodive_test()`), with data frames:

- `size`: `sample, q, depth, method` (`'rarefaction'`, `'observed'`,
  `'extrapolation'`), `estimate, lower, upper, coverage, coverage_lower,
  coverage_upper` (iNEXT's `size_based`: `m, Method, Order.q, qD,
  qD.LCL, qD.UCL, SC, SC.LCL, SC.UCL`).
- `coverage`: `sample, q, coverage, depth, method, estimate, lower, upper`
  (`coverage_based`).
- `asymptotic`: `sample, q, observed, estimate, se, lower, upper`
  (`AsyEst`).
- `info`: `sample, n, s_obs, coverage, f1, ..., f10` (`DataInfo`).

Attributes: `conf_level`, `n_boot`.

`hill_standardize()`: a numeric matrix, samples × orders, columns `q0`,
`q1`, `q2`, so it goes straight into `alpha_test()`. Attributes: `lower`,
`upper` (matrices, absent with `n_boot = 0`), `depth` and `coverage` (per
sample: the depth used and the coverage reached; one is the common level),
`method` (per sample: rarefaction, observed, extrapolation), `level` and
`base` (`'coverage'` or `'depth'`: what was held equal), `conf_level`.


## Reference implementations

- **Output-only:** iNEXT 3.0.2 (CRAN; GPL ≥ 3). `iNEXT(x, q = c(0, 1, 2),
  datatype = 'abundance', size, endpoint, knots, nboot, conf)`,
  `estimateD(x, q, base = 'size' | 'coverage', level)`, `ChaoRichness()`,
  `ChaoShannon()`, `ChaoSimpson()`, `DataInfo()`. Its docs (read): defaults
  `knots = 40`, `nboot = 50`, `conf = 0.95`, `endpoint` = double the
  reference size, about half the knots on each side of \eqn{n}; outputs as
  listed above.
- **Fixtures:** `data-raw/expected/hill_extrapolation.r` (iNEXT installed
  from CRAN into a scratch library): the black-box probes below, plus
  iNEXT's `spider` data (Chao et al. 2014's worked example; its Table 3
  frequency counts are also in the paper, so the paper's printed values
  are a second check: girdled \eqn{n = 168}, \eqn{S_{obs} = 26}, observed
  \eqn{q = 1, 2}: 12.06, 7.84; coverage 93%), three `ex_counts` samples,
  and one deep simulated sample. Store `size_based`, `coverage_based`,
  `AsyEst`, `DataInfo` as text under `data-raw/expected/hill_extrapolation/`;
  for intervals, `nboot = 2000` runs to compare SEs.
- **Black-box protocol (do first; record every answer in this plan):**
  (a) the knot rule, from `m` for \eqn{n} = 50, 168, 1000 and several
  `endpoint`s; (b) the \eqn{q = 1} extrapolation formula, from
  `size_based` at \eqn{m > n}; (c) \eqn{\hat H}, \eqn{\hat C} and
  \eqn{\hat f_0} with \eqn{f_2 = 0}, with \eqn{f_1 = 0}, and with all
  singletons; (d) the bootstrap sample size and any lower-limit truncation,
  from interval widths and limits on a tiny sample; (e) the
  coverage-based grid and interpolated \eqn{m}; (f) `estimateD()`'s
  default levels and `Method` labels.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots". The input
is a count table with no title, so no `display`.

```r
plot(x, type = 'size', q = NULL, samples = NULL, ...)
```

- `type`: `'size'` (iNEXT's type 1: estimate against depth),
  `'completeness'` (type 2: coverage against depth), `'coverage'` (type 3:
  estimate against coverage). Rarefaction solid, extrapolation dashed, the
  observed point marked; interval bands when `n_boot > 0`; one colour per
  sample, with a legend for up to about 10 samples (`samples =` picks a
  subset).
- `q`: which order (default the first); one order per plot.
- Titles: "Sample-size-based Rarefaction and Extrapolation",
  "Sample Completeness Curves", "Coverage-based Rarefaction and
  Extrapolation"; y label "Species richness" (\eqn{q = 0}), "Shannon
  diversity" (1), "Simpson diversity" (2), as Chao et al. (2014).
- Reads: the `size` / `coverage` tables and `conf_level` (band label).
  Nothing recomputed. Base graphics through `plot_do()`.


## Tests (`tests/tinytests/test-hill_extrapolation.r`)

- Point estimates (all knots, every \eqn{q}, coverage, asymptotes,
  `info`) equal iNEXT's fixtures at 1e-10; `hill_standardize()` equals
  `estimateD()` at 1e-10 for both bases and for levels inside and beyond
  the observed coverage.
- Identities: \eqn{q = 0} interpolation = `rarefaction_curve()`'s expected
  richness; at \eqn{m = n}, the observed values equal `observed()`,
  `exp(shannon())`, `inv_simpson()` (and `hill()` once it exists);
  coverage at \eqn{m < n} equals rarefaction_curves.md's.
- The spider example reproduces the paper's printed numbers.
- Intervals: with `n_boot = 2000`, each SE within 10% of iNEXT's at
  `nboot = 2000` (✱ widen if Monte Carlo error needs it; seeds fixed so
  the test is stable).
- Determinism: `identical()` for `cpus` 1, 2 and 4, and for `n_boot`
  replicate splits; the same `seed` gives identical intervals.
- The binomial sampler: mean and variance over 10⁶ draws for small and
  large \eqn{np}; multinomial sums equal \eqn{n}.
- Edge cases: \eqn{f_2 = 0}; \eqn{f_1 = 0} (coverage 1, no extrapolation
  gain); all singletons; one species; an empty sample (`NA` rows and a
  warning); `depths` beyond twice \eqn{n} (allowed; the docs warn
  extrapolation of richness is unreliable there).
- Errors: non-integer counts; `q` outside \{0, 1, 2\}; both `coverage` and
  `depth`; `coverage` outside (0, 1]; `n_boot` negative; `conf_level`.
- Plot: each `type` runs; title as above.
- Smoke: one line each, `n_boot = 2`.


## Docs

- Rd page `hill_extrapolation` (both functions + `plot()`): what
  interpolation, extrapolation and coverage mean; when to compare at equal
  depth versus equal coverage (Chao & Jost 2012); the reliability limits
  (richness to about twice the depth; Shannon and Simpson further); that
  `hill_standardize()` output feeds `alpha_test()`; cost and `n_boot`.
- "Differences from iNEXT": a different random stream (intervals differ
  within Monte Carlo error), abundance data only, \eqn{q \in \{0, 1, 2\}},
  `hill_standardize()`'s default base if it differs, the result shapes.
- "Same names in other packages": none.
- pkgdown: Alpha Diversity (a "Depth and coverage" subtitle, with
  `rarefaction_curve()`); `vignettes/adiv.Rmd`: an example; NEWS;
  NAMESPACE (`export()` ×2, `S3method(plot, ...)`, `S3method(print, ...)`).
- parameter_names.md §1: `n_boot`, `coverage`.


## Later steps

- Any \eqn{q \ge 0} (profiles): interpolation by Eq. 9a for any \eqn{q};
  extrapolation and asymptotes from Chao & Jost (2015).
- Incidence data (Table 2; iNEXT's `incidence_freq`): a group of samples
  as sampling units, for gamma diversity of a body site or treatment.
- The asymptotic estimators as alpha metrics (`METRICS` rows), e.g. Chao,
  Wang & Jost's Shannon estimator per sample.


## Dependencies

- `rng_t` / `rng_unif()` (`src/rng_dist.c`); the binomial port is new.
- Coverage and running-product helpers shared with
  [rarefaction_curves.md](rarefaction_curves.md).
- `hill()` from [alpha_metrics.md](alpha_metrics.md) for the identity
  tests (optional).
- [alpha_test.md](alpha_test.md) consumes `hill_standardize()`'s matrix.


## References

- Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. *Ecological Monographs*, 84, 45–67. doi:10.1890/13-0133.1 (main text read from the authors' copy; Appendix G read from Ecological Archives M084-003)
- Chao, A., & Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. *Ecology*, 93, 2533–2547. doi:10.1890/11-1952.1
- Chao, A., Wang, Y. T., & Jost, L. (2013). Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. *Methods in Ecology and Evolution*, 4, 1091–1100. doi:10.1111/2041-210X.12108
- Chao, A., & Jost, L. (2015). Estimating diversity and entropy profiles via discovery rates of new species. *Methods in Ecology and Evolution*, 6, 873–882. doi:10.1111/2041-210X.12349
- Colwell, R. K., Chao, A., Gotelli, N. J., Lin, S.-Y., Mao, C. X., Chazdon, R. L., & Longino, J. T. (2012). Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. *Journal of Plant Ecology*, 5, 3–21. doi:10.1093/jpe/rtr044
- Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). *Methods in Ecology and Evolution*, 7, 1451–1456. doi:10.1111/2041-210X.12613
- Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scandinavian Journal of Statistics*, 11, 265–270.
- Good, I. J. (1953). The population frequencies of species and the estimation of population parameters. *Biometrika*, 40, 237–264. doi:10.1093/biomet/40.3-4.237
- Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. *Ecology*, 54, 427–432. doi:10.2307/1934352
- Kachitvichyanukul, V., & Schmeiser, B. W. (1988). Binomial random variate generation. *Communications of the ACM*, 31, 216–222. doi:10.1145/42372.42381 (BTPE, the algorithm of NumPy's sampler)
