# `linda()`: linear models for differential abundance (LinDA)

**Priority:** 38. **Tier:** **clean room**. LinDA's code is GPL: the
standalone `zhouhj1994/LinDA` is GPL-2 (per the review's notes) and the
maintained copy in CRAN **MicrobiomeStat 1.4 is GPL-3** (DESCRIPTION
verified 2026-11-20). Work only from the open-access paper (Zhou et al.
2022, CC-BY), MicrobiomeStat's `man/linda.Rd` (documentation, read), and
the function's **outputs**. Do not open `MicrobiomeStat/R/`. Read
[README.md](README.md) first, and [compositional.md](compositional.md).


## What it is

LinDA fits an ordinary linear model to each feature's **clr-transformed**
abundances and removes the compositional bias from every coefficient. If
\eqn{\log X_{is} = u_s\alpha_i + (1, c_s^\top)\beta_i + \varepsilon_{is}}
for absolute abundances, the clr data satisfy the same model with
\eqn{\alpha_i - \bar\alpha}, so OLS on clr values estimates
\eqn{\alpha_i - \bar\alpha}. When most features don't change, the
**mode** of those estimates is near \eqn{-\bar\alpha}; LinDA subtracts it
(paper, eqs. 4–7). It is the same family as ANCOM-BC (bias-corrected
log-linear regression; Wirbel 2024), much faster (100–1000× in the paper's
Table 1), and Yang & Chen (2023) found it "the only method that maintains
reasonable performance in the presence of strong compositional effects" in
correlated data. Pelto 2025: "acceptably sensitive with minor consistency
deficiencies". Uptake: 494 OpenAlex citations 2022–26 (8 → 66 → 97 → 153
→ 170 a year), retrieved 2026-11-19. compositional.md's 2026-11-20
revision overturned "not planned (GPL)": GPL rules out porting, not
clean-room work.


## Names

`linda()`: the method's published name in lowercase (rule 1; LinDA is
what users type, `MicrobiomeStat::linda()`). Registry row in `DA_TESTS`
([da_test.md](da_test.md)): design `formula`, effect `log2_fc`, base 2,
`zero_free_ok = TRUE`; no aliases needed. "Same names in other packages": MicrobiomeStat and LinDA
have `linda()` (same method).


## Design

### Signature

```r
linda(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
      zero_handling = c('adaptive', 'pseudocount', 'imputation'),
      pseudocount = 0.5, corr_cut = 0.1, outlier_pct = 0.03,
      standardize = TRUE, alpha = 0.05, p_adjust = 'BH', cpus = n_cpus())
```

- `formula`, `data`, `groups`: the model convention of
  [da_test.md](da_test.md) (`groups` = the term of interest, or alone
  `~ groups`); the design from the shared `da_design()` (ex
  `ancombc_design()`). LinDA's `meta.use` documentation says "numerical
  variables are scaled": `standardize = TRUE` centres and scales numeric
  covariates first (**verify** with outputs whether LinDA's reported
  coefficients are per SD).
- Counts vs proportions (LinDA's `feature.dat.type`): whole numbers →
  the count path; otherwise the proportion path (zeros → half the
  feature's smallest non-zero value). Documented; **ask the user** whether
  an explicit argument is preferred.
- `zero_handling`, `pseudocount`, `corr_cut`: LinDA's `zero.handling`,
  `pseudo.cnt`, `corr.cut`, with `'adaptive'` (its `adaptive = TRUE`
  default) choosing between the two. Not `zeros`, which
  [da_sensitivity.md](da_sensitivity.md) uses for zero-replacement methods.
- `outlier_pct`: winsorization (`is.winsor = TRUE`, `outlier.pct =
  0.03`); `0` turns it off.
- Filters (`prev.filter`, `mean.abund.filter`, `max.abund.filter`) are not
  carried over (as `ancombc()`): show the one-liners in the docs.
- `alpha`, `p_adjust`: LinDA's `alpha = 0.05`, `p.adj.method = "BH"`.

### Algorithm (paper and documentation; ✱ = to confirm from outputs)

With \eqn{Y} the n × m count table and \eqn{N_s} sample depths:

1. ✱ **Winsorize** each feature's high values at its `1 - outlier_pct`
   quantile (paper: "winsorization at quantile 0.97"; documentation:
   "replace outliers (high values)"). Confirm on what scale (counts or
   proportions), which quantile type, and whether it precedes zero handling,
   from the returned `feature.dat.use` table.
2. **Zeros** (count path):
   - `'pseudocount'`: add `pseudocount` (0.5) to every value.
   - `'imputation'`: a zero of feature \eqn{i} in sample \eqn{s} becomes
     \eqn{N_s / \max_{k: Y_{ik} = 0} N_k} (paper, Methods), so the deepest
     sample among the feature's zeros gets 1.
   - `'adaptive'` ✱: "the correlation p-value for any explanatory variable
     \eqn{\le} `corr_cut`" between \eqn{\log N_s} and the explanatory
     variables in `formula` → imputation, else pseudocount
     (documentation). Confirm the test (Pearson `cor.test`? `lm` F for
     factors?) and "any variable" vs "the variable of interest" (the paper
     says the latter) from outputs.
3. **clr in log2**: \eqn{W_{is} = \log_2 Y_{is} - \frac1m\sum_j \log_2
   Y_{js}} (the documentation: "centered log2-ratio"; outputs are
   `log2FoldChange`).
4. **OLS** per feature on the design \eqn{X} (intercept + covariates):
   \eqn{\tilde\beta_i = X^+ W_i}, \eqn{\hat\sigma_i^2 = \text{RSS}_i/(n -
   P)}, SE \eqn{= \sqrt{[(X^\top X)^{-1}]_{kk}\hat\sigma_i^2}} for each
   column \eqn{k} (paper eq. 5; `df` = "number of samples minus the number
   of explanatory variables (intercept included)").
5. **Bias** per non-intercept column \eqn{k}: \eqn{b_k = \widehat{\text{mode}}
   (\{\tilde\beta_{ik}\}_i)}, a kernel-density mode (paper eq. 6, any even
   kernel, bandwidth \eqn{h}). ✱ **The kernel, bandwidth and maximiser are
   not specified** in the paper or documentation. Candidates to test
   against LinDA's `bias` output: `density()` defaults (Gaussian,
   `bw.nrd0`, 512-point grid, `cut = 3`) and `x[which.max(y)]`; or the
   exact maximiser of the Gaussian KDE (Newton from the grid maximum).
   Paper's \eqn{\sqrt n} scaling cancels for a data-driven bandwidth.
   Settle this first; if only grid resolution separates two candidates,
   craft coefficient sets (designs with known OLS solutions) to tell them
   apart.
6. **Correct and test**: \eqn{\hat\beta_{ik} = \tilde\beta_{ik} - b_k},
   \eqn{T = \hat\beta/\text{SE}}, \eqn{p = 2F_{t, n-P}(-|T|)}, then
   `p_adjust` per covariate; `reject = p_adjusted <= alpha`.

### Sharing with `ancombc()`

Steps 3–4 are ANCOM-BC's initial fit without θ; use the shared OLS kernel
`src/feature_lm.c` ([da_test.md](da_test.md)) with the shared design.
LinDA has no missing values (zeros are filled), so no masked fits. The
mode is new: put the KDE-mode estimator in `src/feature_lm.c` too
(threaded over covariates; deterministic grid + refinement).

### Result

`data.frame` of class `c('ecodive_linda', 'data.frame')`, one row per
feature × non-intercept covariate (LinDA reports `variables` without the
intercept): `feature, covariate, log2_fc, se, t, dof, p_value,
p_adjusted, reject`. Attributes: `formula`, `groups`, `bias` (per
covariate, LinDA's `bias`), `zero_handling` (the one used, after
`'adaptive'`), `pseudocount`, `outlier_pct`, `alpha`, `p_adjust`,
`design`, `vcov` (feature × covariate × covariate, LinDA's `covariance`),
`mean_clr`, `resid_sd` (for [da_plots.md](da_plots.md)). LinDA's
`baseMean` ("2 to the power of the intercept coefficients (normalized by
one million)") is ambiguous; leave it out. `plot()` as
`plot.ecodive_ancombc()` (volcano, effect, ma, forest).

### Later: mixed effects

The paper's Remark 1 and `linda(formula = '~x + (1|id)')`: `lmer` per
feature on the clr data, the same mode correction, Satterthwaite df from
lmerTest. Needs a REML random-intercept fit in C; share it with
[maaslin3.md](maaslin3.md) and [dirmult_lm.md](dirmult_lm.md) (scikit-bio's
BSD `dirmult_lme()` random-intercept REML is a permissive reference).
Separate commit, after the fixed-effects version.


## Reference implementations

- **Output-only:** MicrobiomeStat 1.4 (CRAN, GPL-3) `linda()`; LinDA
  0.2.0 (GitHub, GPL-2 per the review's notes; not re-verified here).
- **Clean-room sources:** Zhou et al. (2022), open access (PMC9012043),
  Methods "Setup" to "Testing procedure"; `man/linda.Rd` of MicrobiomeStat
  1.4 (argument defaults and output definitions above).
- **Fixtures:** `data-raw/expected/linda.r` runs MicrobiomeStat on
  (a) two groups, (b) a numeric covariate, (c) a binary covariate plus two
  confounders (the paper's three settings), each with `adaptive = FALSE`
  and both zero handlings, and once with `adaptive = TRUE` where depth
  differs between groups; store `output`, `bias`, and `feature.dat.use`
  (to pin down steps 1–2) as text under `data-raw/expected/linda/`.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Attributes plots read: `alpha`, `log_base` (`2` if LinDA's
  log2-scale coefficients are kept, else `exp(1)`; say which), `bias`
  (the estimated mode per covariate, for an optional bias display), and
  `mean_clr`, `resid_sd` (computed when the result is made).
- `plot()`: [da_plots.md](da_plots.md)'s types (volcano, forest), titled
  from `method` via `DA_TESTS`, for example "Volcano plot of LinDA", with
  the covariate in the subtitle. The input is a count table, so no
  `display`.

## Tests

- Steps 1–4 exact at 1e-10 against `feature.dat.use` and the biased
  coefficients (`wald$beta`), once the ✱ items are settled; final
  `log2_fc`, SE, p at 1e-8 (the mode estimator's tolerance); `reject`
  identical.
- On a zero-free table with the proportion path (no pseudocount) and
  `outlier_pct = 0`, invariant to multiplying a sample by a constant (the clr cancels it); a
  pseudocount breaks this, so don't test it there.
- `identical()` for `cpus` 1, 2, 4; dense = sparse.
- Edge cases: n ≤ P (error), a single-level factor, NA in `data`,
  all-equal coefficients (mode = that value), D = 2, a feature of zeros.
- Smoke line.


## Docs

Rd page in *Differential abundance*: the model, the mode correction and
its assumption (most features unchanged), the adaptive zero handling,
units (log2), and a "Differences from LinDA" section (filters not carried
over; design coding as `ancombc()`; any ✱ choice that couldn't be matched
exactly). NEWS. New glossary entries: `zero_handling`, `corr_cut`,
`outlier_pct`.


## References

- Zhou, H., He, K., Chen, J., & Zhang, X. (2022). LinDA: linear models for
  differential abundance analysis of microbiome compositional data.
  *Genome Biology*, 23, 95. doi:10.1186/s13059-022-02655-5 (494 OpenAlex
  citations 2022–26, retrieved 2026-11-19).
- Yang, L., & Chen, J. (2023). Benchmarking differential abundance analysis
  methods for correlated microbiome sequencing data. *Briefings in
  Bioinformatics*, 24(1), bbac607. doi:10.1093/bib/bbac607.
- Wirbel, J., et al. (2024). *Genome Biology*, 25, 247.
  doi:10.1186/s13059-024-03390-9.
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130.


## Decisions made while implementing (item 38, 2026-11-20)

Asked of the user at the start:

- **`pseudocount` default: the shared rule (`NULL`)**, not LinDA's 0.5:
  nothing added to a table without zeros, else half the smallest non-zero
  (winsorized) value with the shared warning, which is 0.5 for counts with
  ones. An explicit `pseudocount = 0` is allowed without zeros and an error
  with them (as `ancombc()`). Fixtures pass `pseudocount = 0.5`.
- **Counts or proportions: detected** (all stored values whole numbers →
  counts), no argument.
- **`standardize = TRUE` default** (LinDA scales numeric covariates; see
  below), `FALSE` for per-unit coefficients.
- **The mode: LinDA's exact estimator**, matched to ~1e-15 (below).

Settled from outputs (MicrobiomeStat 1.4 and modeest 2.5.0, both GPL-3,
run only; neither's code was opened). No GPL source text reached this
session (MicrobiomeStat's own error messages, e.g. "the condition has
length > 1", quote no code).

- **The mode (step 5):** `wald$bias` equals `modeest::meanshift(b)` with
  its defaults to 1e-16 (found by testing candidates: `density()` grid
  maxima, exact Gaussian-KDE maxima by `optimize()` / `optim()` / `nlm()`,
  modeest's `parzen()`, `tsybakov()`; only `meanshift()` matched, and its
  default bandwidth is `bw.nrd0()`). From modeest's documentation
  (Rd: Gaussian kernel, start `par = shorth(x)`, `iter = 1000`,
  `tolerance = sqrt(eps)`) and outputs:
  - the iteration is \eqn{x \leftarrow \sum w b / \sum w}, \eqn{w =
    \phi((x - b)/h)}; it stops when \eqn{|x_{new} - x| < \sqrt\epsilon
    |x|} and returns **the previous** \eqn{x} (the returned value equals
    the iterate before the last step, in all 2,000 random cases, with the
    same iteration count as its `iterations` attribute); at the cap it
    returns the last iterate;
  - `shorth()`: the mean of the shortest window of `ceiling(n/2)` sorted
    values; among equally short windows, the one starting at
    `floor(mean(starts))` (matched modeest on 1,078 tied cases and 2,000
    untied ones; other tie rules failed).
  `linda_shorth()` in R, `bw.nrd0()` from stats, the iteration in C
  (`C_linda_mode`, long double sums). Agreement with modeest: 1e-14
  relative.
- **Winsorization (step 1):** for counts, proportions `y / N` (N: the
  original sample totals), capped per feature at the type-7 quantile
  `1 - outlier_pct`, multiplied back by N and **rounded** (`round()`, ties
  to even; `nearbyint()` in C). For proportions, the values as given are
  capped, not rounded. Only high values. Matched exactly (counts direct,
  other quantile types, no rounding all failed).
- **Zeros (step 2):** pseudocount added to every value; imputation as the
  paper, \eqn{M_s / \max_{k: Y_{ik} = 0} M_k} with **M the winsorized
  totals** (the original totals failed); proportions: half the feature's
  smallest non-zero (winsorized) value; a proportion feature with none is
  an error here (MicrobiomeStat gives `Inf` warnings, then fails).
- **Adaptive choice:** `summary(lm(log(M) ~ X))` t-tests of every design
  column but the intercept (M the winsorized totals); imputation if any
  p ≤ `corr_cut`. Matched MicrobiomeStat's reported choice in 300 of 300
  random cases (factor + numeric designs); per-column correlation tests,
  the F-test, Spearman and the unwinsorized totals did not. **MicrobiomeStat
  1.4 then adds the pseudocount whatever it chose** (its
  `feature.dat.use` and coefficients are those of the pseudocount even when
  it prints "Imputation approach is used."); `linda()` applies the choice,
  as documented and as the paper intends, and says so in "Differences from
  LinDA". The fixtures record MicrobiomeStat's reported choice and compare
  with its explicit run for that choice. `adaptive = FALSE` without an
  explicit `zero.handling` errors in MicrobiomeStat 1.4.
- **Standardization:** numeric (including integer) variables of `formula`
  are `scale()`d before `model.matrix()` (so `x:z` is the product of the
  scaled variables); logical and character ones are not. Added as
  `da_design(standardize = TRUE)`; a constant numeric variable is an error,
  and a non-finite design (e.g. `log(x)` of a standardized variable) is an
  error naming `standardize = FALSE` (LinDA can't take function terms at
  all: it reads parentheses as random effects).
- **clr and fit (steps 3–4):** log2, centred over all features per
  sample; OLS on the shared `flm_*` kernel (one QR of the design); `df` =
  n − rank (LinDA: n − P, the same for a full-rank design); p-values
  `2 pt(-|t|, df)`; `p_adjust` per covariate over non-NaN p-values
  (`da_adjust()`).

Made while implementing:

- **C (`src/linda.c`):** `C_linda_prep` (winsorized totals, each feature's
  smallest non-zero value and zero count), `C_linda_fit` (pass 1: per-sample
  sums of log2 values in blocks of 64 features, added in block order;
  pass 2: each feature's fit), `C_linda_mode`. The table is read by feature
  with `da_table_init(norm = none)`; nothing dense is kept. Results are
  identical for `cpus` 1, 2, 3, 4, 8 and for dense, `Matrix`, `slam`,
  `margin = 2` and integer input. 300 × 5000 counts: 0.47 s on one thread,
  0.28 s on four, against MicrobiomeStat's 12.7 s.
- **Result:** columns `feature, covariate, log2_fc, se, t, dof, p_value,
  p_adjusted, reject`, every design column but the intercept. Attributes
  `method`, `formula`, `groups`, `design` (standardized), `standardize`,
  `zero_handling` (what was done: `pseudocount`, `imputation`, or
  `half_min` for proportions), `pseudocount` (0 when none was added),
  `outlier_pct`, `alpha`, `p_adjust`, `log_base = 2`, `bias` (per reported
  covariate), `vcov` (D × P × P, intercept included, \eqn{s^2 (X'X)^{-1}}),
  and per feature `mean` (mean centred log2 value; named `mean`, not the
  plan's `mean_clr`, as `da_lm()` and `ancombc()` name it for item 35's
  MA plots) and `resid_sd`. `baseMean` left out, as planned.
- **Registry:** `LinDA, linda, formula, log2_fc, 2, TRUE`, no aliases (an
  alias starting `linear` would break `match_da_test('linear')` → `da_lm`).
- **Plots (after item 35 landed mid-task, at the user's request):**
  `plot.ecodive_linda()` calls `da_plot()` with linda's columns: volcano
  (default), effect (`resid_sd`, two-level covariates), MA (`mean`), and
  forest with t intervals `log2_fc ± qt((1 + conf_level)/2, dof) se`
  (`conf_level = 0.95`, as `plot.ecodive_ancombc()`'s Wald intervals).
  Default covariate from `da_term_covariate()`.
- **`da_sensitivity()`:** `pseudocounts` with `zero_handling =
  'imputation'` is an error (every run would be the same); with
  `'adaptive'` choosing imputation the runs are the same, as documented.
  The zeros runs give linda proportions without zeros (`half_min`, nothing
  to replace).
- **Edge cases:** at least two features; every sample with a positive
  total (counts); `~ 1` is an error; a rank-deficient design warns and its
  aliased columns are `NaN` (bias `NaN`); all-equal coefficients give
  themselves as the mode.

Reference results (`data-raw/expected/linda.r`, 20 cases on 7 tables:
two groups, a numeric covariate, a factor with confounders, an
interaction, a logical covariate, depth differing between groups, sparse
data, proportions, rare features whose coefficients tie, D = 2; five
adaptive cases, three choosing imputation and two the pseudocount):
log2 fold changes to 2.9e-15 absolute (2.8e-12 relative, on one near 0),
standard errors 8.4e-16, p-values 1.5e-14 and adjusted ones 1.3e-14
relative, bias 1.3e-15; `reject` and `dof` identical. (A first check on
other random tables gave up to 4e-11 relative, again on a fold change near
0.) Mode estimator: 42 vectors (ties, n = 2 to 500), shorth to 1e-13, mode
to 1e-14.

Not done (later): mixed models (`(1 | id)`), as planned.
