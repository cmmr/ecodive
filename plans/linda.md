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
