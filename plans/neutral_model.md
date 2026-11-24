# `neutral_model()`: Sloan's neutral community model

**Priority:** 51. **Tier:** port (neufit, BSD-2) + formula; Burns et
al.'s R code is output-only. References checked 2026-11-20:

| Reference | Licence (how verified) | Use |
|:----------|:-----------------------|:----|
| **neufit** (github.com/misieber/neufit, `neufit.py`, Michael Sieber 2018; last push 2025-11-13) | **BSD-2**: the licence text is in `neufit.py`'s header ("Copyright (C) 2018 Michael Sieber ... Redistribution and use ... permitted"); GitHub detects no LICENSE file | **port source** (read; it is short) |
| Burns et al. (2016) `sncm.fit()` (the paper's supplementary R code) | none stated: all rights reserved | output-only, if the user agrees to run it |
| tyRa 0.1.0 `fit_sncm()` (github.com/DanielSprockett/tyRa) | `DESCRIPTION` says "MIT + file LICENSE", but the LICENSE file is **missing** (404), and the function's docs say its author is Adam Burns and it was "originally published in Burns et al. (2016)" | **unclear provenance: don't port**; outputs only. Its `R/modeling.R` was read while writing this plan (before the provenance note was seen); nothing below depends on it that isn't in the papers or neufit, and the items it would have settled are marked ✱ for outputs |
| MicEco 0.10.0 `neutral.fit()` | GPL (≥ 3) (`DESCRIPTION`) | output-only; it fits m by **maximum likelihood**, not NLS (its Rd), so it is a different estimator |
| `minpack.lm::nlsLM`, `Hmisc::binconf` (what Burns' code calls) | GPL | not needed: ecodive needs only `stats::optimize()` and `stats::pbeta()` (base R, no new dependency) |

Read [README.md](README.md) first.


## What it is

Sloan et al. (2006) derived, for a neutral metacommunity with immigration,
the expected relationship between a taxon's **mean relative abundance in
the metacommunity** p and the **fraction of samples it is detected in**
(occurrence frequency). Under neutrality a taxon's relative abundance in a
local community is Beta(N m p, N m (1 − p))-distributed, where N is the
number of individuals per community and m the immigration probability, so
the probability of detection above the limit d = 1/N is

freq(p) = 1 − I_d(N m p, N m (1 − p)),

with I_d the regularized incomplete beta function (`pbeta(d, ...,
lower.tail = FALSE)`). Burns et al. (2016) popularized it for microbiomes:

- fit **m** by nonlinear least squares of observed occurrence on the curve,
  N fixed at the mean number of reads per sample;
- report **R²** (1 − SS_res / SS_tot), **Nm** (metacommunity size times
  immigration);
- draw **95% binomial (Wilson) confidence intervals** around the curve
  (n = number of samples) and partition taxa into **above** (more frequent
  than neutral: selected by the host / habitat, or dispersing well),
  **below** (less frequent: dispersal-limited or selected against) and
  **neutral**;
- compare the fit with a **binomial (random sampling) model** by AIC.

**Demand** (retrieved 2026-11-20): Sloan et al. (2006) 1,470 OpenAlex
citations 2022–26; Burns et al. (2016) 595; Europe PMC microbiome full-text
mentions 2022–26 of "neutral community model" / "Sloan neutral": 531.
Baseline: Bray-Curtis 28,581.

**Why ecodive helps.** Not speed (the fit is cheap) but availability: no
established CRAN package did it until the recent ncmR (not checked for
licence or method), so users copy Burns' supplementary script, tyRa or
MicEco (a different estimator). A documented,
tested function with a plot is the value.


## Names

`neutral_model()`. Rule 4: the acronyms (NCM, SNCM) are cryptic; "neutral
model" is how papers refer to it ("Sloan's neutral model"). The docs name
NCM, SNCM, MicEco's `neutral.fit()` and tyRa's `fit_sncm()` in "Same names
in other packages". New glossary rows: `pool` (a count table defining the
source metacommunity, Burns' `pool`), and `conf_level` (existing) for the
interval.


## Design

### Signature

```r
neutral_model(counts, pool = NULL, margin = 1L, conf_level = 0.95,
              cpus = n_cpus())
```

- `counts`: integer counts, **rarefied to equal depth** (Burns: "All samples
  must be rarefied to the same depth"; neufit rarefies itself). If depths
  differ, warn ("neutral_model() assumes equal depths; see rarefy()") and
  use N = the mean depth, as Burns' code does; don't rarefy silently.
- `pool`: optional count table (same features, by name) whose mean
  relative abundances define p; default the samples themselves.
- `conf_level`: the Wilson interval's level (Burns: 0.95).
- `cpus`: only for the per-feature sums over a large sparse table (C pass
  below); the fit is serial.

### Algorithm

1. **Per-feature summaries** (one pass over the stored values; a small C
   function `C_neutral_stats()` threaded over feature blocks, or
   `colSums`-style R on the sparse formats): N = mean sample total;
   p_i = mean over samples of count_i / N (equivalently total_i / (n N),
   neufit's form; ✱ Burns divides the mean count by N, the same thing);
   freq_i = fraction of samples with count_i > 0. With `pool`, p from
   `pool` and freq from `counts`. Drop taxa with p_i = 0 or freq_i = 0
   (Burns), and keep their names in an attribute `dropped`.
2. **Fit m** by least squares on the curve: minimize RSS(m) =
   Σ_i (freq_i − F(p_i; N, m))², F = `pbeta(1/N, N m p, N m (1 − p),
   lower.tail = FALSE)`. One parameter, so no general NLS is needed:
   - `stats::optimize()` (Brent) over log m in [log 1e-12, log 1e4] (Burns
     doesn't bound m; neufit bounds it to [0, 1]; see below) with
     `tol = 1e-12`, then a few Newton steps on dRSS/dm by central
     differences to polish to ~1e-12 relative.
   - **Multiple minima:** RSS(m) has been unimodal in practice, but check
     by evaluating RSS on a log grid (e.g. 50 points over the bracket)
     first and starting Brent in the best interval. If two interior minima
     are within 1e-8 relative RSS, warn.
   - `nlsLM` (Levenberg-Marquardt from m = 0.001) stops at relative
     tolerance ~1.5e-8 (`sqrt(.Machine$double.eps)`), so ecodive matches
     Burns' m to ~1e-7, and an exhaustive `optimize()` reference to 1e-10.
   - **m > 1:** the curve is defined for any m > 0, but m is a probability.
     Burns' fit can return m > 1 on near-saturated data; neufit caps at 1.
     ecodive follows Burns (no cap) and warns when m̂ > 1. ✱ Ask the user
     if they prefer neufit's bound.
3. **Interval for m:** the profile interval of the one-parameter NLS
   problem: the m where τ(m)² = (RSS(m) − RSS(m̂)) / s² equals the cutoff,
   s² = RSS(m̂)/(K − 1) (K taxa), solved by `uniroot()` on each side. ✱ The
   cutoff (`qt(1 − α/2, K − 1)²` or `qf(1 − α, 1, K − 1)`, the same; or a
   normal quantile) is settled from Burns' `confint()` output. R's
   `confint()` on an `nls` fit interpolates a sampled profile, so the
   reference's ends match ecodive's exact roots only to ~1e-3; the docs say
   so. Report `m_lower`, `m_upper` (Burns reports half-widths `m.ci`).
4. **Goodness of fit:** R² = 1 − Σ (freq − pred)² / Σ (freq − mean freq)²
   (Burns; can be negative); RMSE = sqrt(Σ (freq − pred)² / (K − 1)) (✱
   divisor from outputs).
5. **Wilson interval** (Wilson 1927) around each prediction, with x =
   pred × n and n = number of samples: centre (x + z²/2)/(n + z²), half-width
   z sqrt(x(n − x)/n + z²/4)/(n + z²), z = `qnorm(1 − (1 − conf_level)/2)`.
   neufit uses statsmodels' Wilson (BSD); Burns uses `Hmisc::binconf`. ✱
   Check the x = 0 and x = n edges against both.
6. **Partition:** `above` if freq > upper, `below` if freq < lower, else
   `neutral` (neufit and Burns agree; strict inequalities).
7. **Alternative models and AIC** (Burns et al. 2016, Methods): the
   binomial model freq = 1 − (1 − p)^N, i.e. `pbinom(1/N, N, p, lower.tail =
   FALSE)` (✱ form from outputs), and optionally the Poisson model
   1 − exp(−N p). For each model, AIC and BIC of a Gaussian likelihood on
   its residuals. The neutral model has parameters (m, σ); its maximum
   likelihood m is the least-squares m̂ (the Gaussian likelihood is maximal
   where RSS is), so log L = −K/2 (log(2π σ̂²) + 1), σ̂² = RSS/K, in closed
   form. ✱ For the binomial and Poisson models Burns fits a mean μ and σ to
   the residuals: log L at μ̂ = mean residual, σ̂² = variance (1/K). Burns
   obtains these with `stats4::mle()` (numerical BFGS), so expect agreement
   to ~1e-6, not 1e-10; ecodive's closed forms are exact. Settle the number
   of parameters counted (k = 2 for each, as `AIC(k = 2)` with two
   parameters) from outputs. Report ΔAIC = AIC_binomial − AIC_neutral.

All of this is plain R over K-length vectors (`pbeta`, `optimize`,
`uniroot` are base R). Only step 1 touches the table.

### Result

`c('ecodive_neutral_model', 'data.frame')`, one row per taxon kept, ordered
by `mean_abundance` (Burns and neufit sort by p): `feature, mean_abundance,
occurrence, predicted, lower, upper, partition` (factor: `above`,
`neutral`, `below`), and `predicted_binomial`.

Attributes (results, rule 1): `method = 'neutral_model'`, `m`, `m_lower`,
`m_upper`, `nm` (N m), `n_individuals` (N), `detection` (1/N),
`n_samples`, `r_squared`, `rmse`, `aic`, `bic`, `r_squared_binomial`,
`aic_binomial`, `bic_binomial` (and `_poisson` if kept), `conf_level`,
`dropped` (features without reads or without occurrences). Print with
`print_table_result` plus a header that shows m, Nm, R² and the counts per
partition (`print_df_label()`).


## Reference implementations

- **Port:** neufit (BSD-2): its model function (`beta_cdf()`: the
  difference of beta CDFs at 1 and 1/N), the p and occurrence definitions,
  the Wilson band and the above/below rule. Nothing else is needed from
  it; credit Michael Sieber in `inst/COPYRIGHTS` and `Authors@R`
  (`ctb`/`cph`) only if code is ported closely enough to need it (the
  formulas are the paper's; decide when writing).
- **Output-only fixtures** (`data-raw/expected/neutral_model.r` and
  `neutral_model.py`):
  - neufit on its `sim_data.csv` and on a rarefied `ex_counts`-like table:
    m, R², predictions, intervals (Python; lmfit's least squares, so m to
    ~1e-7).
  - Burns' `sncm.fit()` (if the user agrees to run the unlicensed
    supplementary code as a black box) or tyRa's `fit_sncm()` (outputs
    only): m, m.ci, R², RMSE, AIC / BIC of the neutral, binomial and
    Poisson models, predictions and intervals, on the same tables plus one
    with a `pool`.
  - MicEco `neutral.fit()` for the docs' comparison only (ML estimator,
    expected to differ slightly).
  - A high-precision reference for m: in the fixture script, plain R
    `optimize()` with `tol = 1e-15` from several brackets, stored.
- **Black-box protocol (first):** RMSE divisor; the binomial model's form;
  the AIC parameter count and whether μ is fitted for the neutral model too;
  the `confint()` cutoff; Wilson edges at 0 and n. Record answers here.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: "Neutral Community Model" (no input title: the input is a count
  table).
- `plot()`: occurrence frequency (y) against mean relative abundance (x,
  log scale), points coloured by `partition`; the fitted curve and the
  Wilson band drawn from `n_individuals`, `m`, `n_samples` and
  `conf_level` over a grid of p (a function of stored results, not of the
  data: allowed); subtitle "m = ..., Nm = ..., R² = ..." from attributes;
  legend with the counts per partition. `type = 'binomial'` adds the
  binomial model's curve.


## Tests

- m within 1e-10 relative of the high-precision `optimize()` reference and
  1e-6 of Burns / tyRa / neufit; R², predictions, Wilson bounds at 1e-10
  given the same m (feed the reference's m through an internal function);
  AIC / BIC at 1e-6 (numerical `mle()` in the reference).
- The closed-form Gaussian log-likelihood equals a direct `optim()` of the
  likelihood in the test at 1e-8.
- Simulated neutral data (Sloan's process: each sample a Dirichlet-
  multinomial with α = N m p, drawn in the fixture script) recovers m
  within its interval in most replicates (at home only).
- Edge cases: one sample (error: occurrence needs ≥ 2), all taxa in every
  sample (R² undefined → NA with a warning), unequal depths (warning),
  `pool` with missing features (error), proportions instead of counts
  (error: N must be a read count), a single taxon left after dropping.
- Dense and sparse inputs agree; `cpus` 1, 2, 4 identical.
- Smoke-test line.


## Docs

Rd page in a new "Community assembly" group (with
[assembly_processes.md](assembly_processes.md)'s functions). Explain the
model and its assumptions (equal depths, a single metacommunity, neutral
drift and immigration only), what m, Nm and the partitions mean and don't
mean (fit, not a test of neutrality; R² can be negative; above/below taxa
are candidates, not proof of selection), and "Differences from Burns et
al.'s code" (exact optimizer, closed-form AIC, profile interval by root
finding, no cap on m, ✱ unmatched items). Compare with MicEco's ML fit.
NEWS.


## References

- Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S., & Curtis,
  T. P. (2006). Quantifying the roles of immigration and chance in shaping
  prokaryote community structure. *Environmental Microbiology*, 8,
  732–740. doi:10.1111/j.1462-2920.2005.00956.x
- Burns, A. R., Stephens, W. Z., Stagaman, K., Wong, S., Rawls, J. F.,
  Guillemin, K., & Bohannan, B. J. M. (2016). Contribution of neutral
  processes to the assembly of gut microbial communities in the zebrafish
  over host development. *ISME Journal*, 10, 655–664.
  doi:10.1038/ismej.2015.142
- Sieber, M., Pita, L., Weiland-Bräuer, N., et al. (2019). Neutrality in
  the metaorganism. *PLoS Biology*, 17, e3000298.
  doi:10.1371/journal.pbio.3000298 (neufit)
- Wilson, E. B. (1927). Probable inference, the law of succession, and
  statistical inference. *Journal of the American Statistical
  Association*, 22, 209–212. doi:10.1080/01621459.1927.10502953
- Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*.
  Prentice-Hall. (The method of `stats::optimize()`.)


## Decisions made while implementing

Asked of the user at the start (2026-11-24):

- **m̂ > 1** (the ✱ point): **no cap**, as Burns et al.; a warning says
  when the fitted m exceeds 1 (the user's choice). neufit's bound at 1 is
  listed under "Differences".
- **Burns et al.'s reference:** the user chose to download the original
  supplementary script (Supplementary Code 1, `ismej2015142x4.txt`) as well
  as iCAMP and tyRa. It could not be fetched: nature.com now answers 404
  (ISME J moved to Oxford University Press), PMC serves a reCAPTCHA page,
  OUP answers 403, and the Wayback Machine's copy of the nature.com file
  was refused by the session's permission classifier (see "Outcome").
  Burns et al.'s code was therefore run through **iCAMP 1.9.1 `snm()`**
  (CRAN; its docs say it holds Burns et al.'s original code) and **tyRa
  0.1.0 `fit_sncm()`** (jsDelivr; credited to Adam Burns), both as black
  boxes: `sys.source()` of tyRa's `R/modeling.R` without viewing it, and
  only their documentation read.
- **Uncertainty of m:** the exact **profile interval** (`m_lower`,
  `m_upper`) **and** `m_se`, the delta-method standard error that neufit
  (lmfit) prints (the user's choice).
- **Extras (all three chosen):** the **Poisson model**
  (`predicted_poisson` and `*_poisson` attributes), **`plot(type =
  'binomial')`**, and **`groups` / `data`**: one fit per group, as
  subsetting `counts` would, returned in one table with a leading `group`
  column and attributes that become vectors named by group (`dropped`
  becomes a list).

Black-box protocol (settled from iCAMP's and tyRa's outputs on neufit's
`sim_data.csv`):

- **RMSE divisor:** K − 1 (`sqrt(SS / (K - 1))`), for all three models.
- **Binomial model:** `pbinom(1/N, N, p, lower.tail = FALSE)` = 1 − (1 −
  p)^N, computed as `-expm1(N * log1p(-p))` so that a non-integer mean
  depth works (matches `pbinom` to 1e-15). Poisson: 1 − e^(−Np).
- **AIC / BIC:** k = 2 for every model; BIC uses log K. The neutral model's
  likelihood has mean 0 and parameters (m, σ); the binomial's and
  Poisson's have a fitted mean and σ (variance with divisor K). iCAMP's
  `maxLL`, `binoLL`, `poisLL` are *negative* log-likelihoods. The closed
  forms match iCAMP's numerical `mle()` to 2e-6 (neutral) and 1e-10
  (binomial, Poisson), and a direct `optim()` in the test to 1e-8.
- **`confint()` cutoff:** F(1, K − 1) = t²(K − 1) (not the normal
  quantile): ecodive's exact ends match iCAMP's interpolated profile to
  ~1e-5. tyRa's `m.ci` is m − lower.
- **Wilson edges:** exactly 0 at x = 0 and 1 at x = n (Hmisc,
  statsmodels). ecodive computes the lower end as x² / (n (x + z²/2 +
  z√S)) and the upper as 1 − lower(n − x): the same without cancellation,
  matching both references to 1e-14 everywhere and exact at the edges.
- **Order of rows:** by mean abundance, ties in input order (`order(method
  = 'radix')`); the references order ties differently, so tests match rows
  by feature.

Decided while implementing:

- **Signature:** `neutral_model(counts, groups = NULL, data = NULL, pool =
  NULL, margin = 1L, conf_level = 0.95, cpus = n_cpus())`. `groups` and
  `data` come right after `counts`, as in every other ecodive function that
  takes them; `pool` (second in the plan's signature) follows.
- **`pool`:** p is each feature's share of **all** the pool's counts
  (iCAMP's definition). Burns et al.'s (tyRa) divides the feature's mean
  count in the pool by the N of `counts`, which is the same when the pool
  is rarefied to that depth (as in the fixture, which then matches tyRa).
  Matched by feature name (by position when either table has no feature
  names, with an error if the counts differ); a feature **present** in
  `counts` but missing from `pool` is an error (an absent one is simply
  dropped). Pool features not in `counts` still count in its total.
- **Dropping:** features with p = 0 or occurrence 0 are dropped (Burns,
  tyRa) and named in `dropped`; iCAMP keeps pool features with p = 0, so
  the pool case is compared with tyRa only.
- **N:** the unrounded mean depth (tyRa); iCAMP rounds it to a whole
  number, so the unequal-depth case is compared with tyRa only.
- **Empty samples** are left out with a warning (the plan didn't say).
- **Degenerate fits:** fewer than 2 features left → error (with `groups`,
  naming the group); every sample at least 2 is required per group; an
  exact fit (RSS = 0, e.g. every feature in every sample with predicted
  occurrence 1) → warning that m is not identifiable, interval `NA`; all
  occurrences equal → R² `NA` with a warning (the plan's case); AIC / BIC
  `NA` when the residual variance is 0; the profile interval's end is 0 or
  `Inf` when the cutoff is never reached on that side.
- **Partition margin:** above / below only by more than 1e-12. At
  occurrence 1 the band's upper end is 1 to within rounding, and the
  references' strict comparisons flipped on the last bit (iCAMP said
  "above" where ecodive's strict test said "neutral", and the reverse, on
  two fixture tables). Given a reference's m, ecodive's bands match its
  bands to 2e-16 and its partitions everywhere away from those edges.
- **Optimizer:** a 61-point grid over log m in [1e-12, 1e4]; Brent
  (`optimize(tol = 1e-12)`) around each grid minimum (a warning if two
  distinct minima are within 1e-8 relative RSS, and if the best is at the
  edge of the bracket); then `uniroot(tol = 1e-15)` on the gradient,
  whose per-feature derivatives are central differences with Richardson
  extrapolation. The profile ends are `uniroot()` roots in log m.
- **High-precision reference:** instead of R's `optimize(tol = 1e-15)`
  (the RSS is flat at its minimum, so `optimize()` can't place m better
  than ~1e-8), the Python fixture solves for the gradient's zero and the
  profile ends with **mpmath at 40 digits** (`betainc`). ecodive matches
  it to ~1e-13 in m and ~1e-16 at the interval ends, on five cases
  (neufit's `sim_data`, a Dirichlet-multinomial draw, a rarefied
  uneven-depth table, a near-saturated table with m̂ = 1.22, and a pool).
- **neufit's m** agrees only to ~1.5e-6 (lmfit's stopping rule; refitting
  neufit's model with `ftol = xtol = gtol = 1e-15` doesn't help, since
  lmfit's finite-difference Jacobian limits it). Given neufit's m,
  ecodive's predictions match to 1e-14 and its Wilson bands to 1e-14; R²
  matches to 1e-12 (it's flat in m). The test checks m at 1e-5.
- **The simulation test:** the plan's Dirichlet-multinomial simulation is
  not Sloan's model (a read count ≥ 1 is not x > 1/N): it biased m̂ by
  +15% and the interval covered m in 9 of 20 replicates. The test draws
  occurrences from the model itself (`rbeta`, x > 1/N); then m̂ is
  unbiased and the 95% interval covers m in 86–88% of 200 replicates (it
  assumes equal residual variances; occurrences are binomial). The docs say
  so; the test asks for ≥ 14 of 20.
- **MicEco's `neutral.fit()`** was not run; the docs mention its maximum
  likelihood estimator only qualitatively (iCAMP's `m.mle` showed the size
  of that difference: 1e-4 relative on `sim_data`).
- **C:** `C_neutral_stats()` (`src/neutral_model.c`) returns, per group
  and feature, the total count and the number of samples with the feature,
  in one pass over `da_table_init()`'s feature-major table, threaded over
  features with fixed-order sums (identical for any `cpus`; tested on a
  300 × 1000 table with 1, 2, 4 and 7 threads). Everything else is R over
  K-length vectors. 2,000 × 20,000 sparse (1.2M values) takes ~0.9 s, nearly
  all in `pbeta()`.
- **Licensing:** the R file carries neufit's BSD-2 notice (model function,
  p and occurrence definitions, Wilson band, above/below rule);
  `inst/COPYRIGHTS` has a stanza for it (and for neufit's `sim_data.csv`
  in the test fixture), and Michael Sieber is `ctb`/`cph` in `Authors@R`.
  No source text of Burns et al.'s code, iCAMP, tyRa or MicEco was read.
- **Display:** no `display` attribute (the input is a count table, so the
  title has no "of ..."): `print()` and `plot()` use "Neutral Community
  Model" (with groups, "Neutral Community Model: <group>", one panel per
  group, `group =` picks them). The print header shows m with its
  interval, Nm, R² and the partition counts, one line per group. Row
  subsets keep their attributes in R 4.3, so they print and plot (points
  of the subset over the full fit's curve).
- **pkgdown:** a new "Community Assembly" section after "Quality
  Control".
