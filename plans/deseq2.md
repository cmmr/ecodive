# `deseq2()`: negative-binomial GLMs for differential abundance (DESeq2)

**Priority:** 61. **Tier:** port (**PyDESeq2 0.6.0**, MIT: `LICENSE` "MIT
License / Copyright (c) 2022 Owkin", PyPI 0.6.0 `license` field the same
text, `pyproject.toml` `license = { file = "LICENSE" }`; verified
2026-11-20; the repository now lives under scverse) + **clean room** for
what PyDESeq2 lacks (the likelihood-ratio test, R's `lowess()`, any
DESeq2-only behaviour). **DESeq2 1.52.0** (Bioc 3.23) is **LGPL (>= 3)**
(`DESCRIPTION` verified 2026-11-20): **output-only**, never open its `R/`
or `src/`; its Rd pages and vignette may be read. apeglm 1.34.0 is GPL-2
and ashr 2.2-63 GPL (>= 3) (both `DESCRIPTION`s verified): output-only.
Read [README.md](README.md) first, [da_test.md](da_test.md),
[normalization_factors.md](normalization_factors.md) and
[da_plots.md](da_plots.md).

> **Scope.** This is the largest DA item planned: four optimisation
> problems (size factors, dispersions, the trend, coefficients), outlier
> handling with a refit, two tests, independent filtering and, later,
> shrinkage. It is phased below; phase 1 alone is a complete, useful
> `deseq2()`. **Ask the user to approve the phasing before starting.**


## What it is

DESeq2 (Love, Huber & Anders 2014) fits, per feature \eqn{j}, a negative
binomial GLM \eqn{K_{ij} \sim NB(\mu_{ij}, \alpha_j)}, \eqn{\mu_{ij} =
s_i q_{ij}}, \eqn{\log q_{ij} = x_i^\top \beta_j}, with sample size factors
\eqn{s_i} and a dispersion \eqn{\alpha_j} shrunk towards a mean-dependent
trend shared by all features. Coefficients are tested by Wald (or
likelihood-ratio) tests, with Cook's-distance outlier handling and
independent filtering before Benjamini-Hochberg.

It is the most cited count model in genomics and widely used for microbiome
data since McMurdie & Holmes (2014) recommended it over rarefying (phyloseq
has `phyloseq_to_deseq2()`). Demand (retrieved 2026-11-20): Love et al.
(2014) 74,537 OpenAlex citations (all fields); **10,234 Europe PMC
microbiome full-text mentions 2022–26** (MaAsLin2: 2,232); McMurdie &
Holmes (2014) 1,338 citations.

**The docs must give the benchmark evidence.** Pelto et al. (2025, gut
datasets split into exploratory and validation halves): DESeq2 was among
the two most sensitive methods but among
the least consistent, and with permuted labels 10.5% of its unadjusted
p-values fell below 0.05 (nominal 5%). Wirbel et al. (2024) found that only
the classical tests, limma and fastANCOM controlled false discoveries;
DESeq2 (as the phyloseq vignette runs it) was not among them. Nearing et al.
(2022): intermediate. Weiss et al. (2017) and Hawinkel et al. (2019) also
report inflated FDR on microbiome data. As for LEfSe, the docs point to
`da_sensitivity()` and the plain tests.


## Names

- `deseq2()`: the method's name lowercased (rule 1). Not `deseq()`:
  DESeq (Anders & Huber 2010) is the earlier, different method. DESeq2's
  own entry point is `DESeq()`, PyDESeq2's `DeseqDataSet.deseq2()`; no R
  package exports `deseq2`. "Same names in other packages": none.
- **`DA_TESTS` row:** name `DESeq2`, id `deseq2`, design `formula`, effect
  `log2_fc`, base `2`, `zero_free_ok = FALSE` (it models integer counts).
  Aliases (letters only): `deseq`, `pydeseq` (check with README's
  snippet; none collides).
- Arguments clear the `da_test()` rule (none is a prefix of `counts`,
  `method`, `groups`, `data`, `formula`; `test` is fine there).


## Design

### Signature (phase 1, with later phases marked)

```r
deseq2(counts, formula = NULL, data = NULL, groups = NULL,
       contrast = NULL, margin = 1L, norm = 'rle_poscounts',
       fit_type = c('parametric', 'mean'), cooks = TRUE,
       refit_outliers = TRUE, min_replicates = 7L,
       independent_filter = TRUE, alpha = 0.1, p_adjust = 'BH',
       test = 'wald', reduced = NULL,                  # 'lrt': phase 2
       lfc_threshold = 0, alternative = 'two.sided',   # phase 3
       shrink = 'none',                                # 'apeglm': phase 3
       cpus = n_cpus())
```

- `formula`, `data`, `groups`: the model convention of
  [da_test.md](da_test.md); design from the shared `da_design()`
  (treatment contrasts, radix-sorted levels). DESeq2 also uses treatment
  contrasts with the first level as reference; PyDESeq2 (formulaic) the
  same.
- `contrast = NULL`: one row per feature and non-intercept design column
  (as `linda()`, `da_lm()`); `contrast = c(treatment, reference)` (the
  glossary's meaning: two levels of the `groups` factor): rows for that
  comparison only, from the contrast vector \eqn{c} (DESeq2's
  `results(contrast = c(var, trt, ref))`, PyDESeq2's `contrast`).
- **`norm`**: the size factors, by name, shared with
  [normalization_factors.md](normalization_factors.md):
  `'rle'` (DESeq2 `type = "ratio"`), `'rle_poscounts'` (`"poscounts"`), and
  any other size-factor norm of `size_factors()` (`'css'`, `'tmm'`,
  `'gmpr'`; DESeq2 users set these through `sizeFactors<-`). Default
  **`'rle_poscounts'`**: plain RLE fails whenever every feature has a zero,
  which is nearly every microbiome table (PyDESeq2 then switches to
  `"iterative"`; DESeq2 errors and suggests poscounts), and poscounts is
  what the DESeq2 vignette recommends for such data, derived from
  phyloseq's. **Ask.** `'iterate'` (DESeq2's `"iterate"`, PyDESeq2's
  `"iterative"`) is phase 3.
- `fit_type`: `'parametric'` (default) or `'mean'`. Not `'local'` (DESeq2
  fits it with locfit, GPL; no permissive reference) nor `'glmGamPoi'`.
- `cooks`, `refit_outliers`, `min_replicates`: DESeq2's `cooksCutoff`
  (on/off; the cutoff itself is fixed at \eqn{F_{p, n-p}(0.99)}),
  `minReplicatesForReplace = 7`.
- `alpha`: **0.1**, DESeq2's `results(alpha = 0.1)` (PyDESeq2 uses 0.05);
  it also sets the target of independent filtering, so it changes
  `p_adjusted`, not just `reject`. ecodive's other DA methods use 0.05.
  **Ask the user.** `p_adjust = 'BH'` (DESeq2's `pAdjustMethod`; the
  filtering is defined for it). Other methods are allowed without the
  filter's optimisation over BH (see below).
- `test`, `reduced` (phase 2): `'lrt'` with `reduced` a formula nested in
  `formula` (DESeq2's `test = "LRT", reduced = ~ 1`).
- `lfc_threshold`, `alternative` (phase 3): DESeq2's `lfcThreshold`
  (log2 units) and `altHypothesis` (`'greaterAbs'`, `'lessAbs'`,
  `'greater'`, `'less'`; ecodive spells them `'two.sided'` (=
  greaterAbs), `'equivalence'` (lessAbs), `'greater'`, `'less'`; ask).
- `shrink` (phase 3): `'apeglm'`-style LFC shrinkage (below).

Inputs: whole non-negative counts (DESeq2 and PyDESeq2 error otherwise);
`NA` an error. Features zero in every sample get `NA` results (both
references), not an error.

### Algorithm (port of PyDESeq2 0.6.0 `dds.py`, `ds.py`, `glm.py`, `dispersions.py`, `preprocessing.py`, `stats.py`, `grid_search.py`, `misc.py`)

1. **Size factors** \eqn{s_i} (`fit_size_factors()`, `deseq2_norm_*`):
   - `'rle'`: \eqn{\log g_j} = mean of \eqn{\log K_{ij}} over samples,
     features with any zero dropped; \eqn{s_i = \exp(\mathrm{median}_j(\log
     K_{ij} - \log g_j))}; no rescaling.
   - `'rle_poscounts'`: \eqn{\log g_j = \sum_{K_{ij} > 0} \log K_{ij} / n}
     (zeros add 0, divisor all \eqn{n}), features with any positive count;
     each sample's median over **its own positive** counts; then rescaled
     to geometric mean 1.
   These settle three of normalization_factors.md's "settle from outputs"
   points for RLE (median of log ratios; zeros leave each sample's median;
   rescaling for poscounts only). ✱ Confirm against DESeq2's
   `estimateSizeFactorsForMatrix()` outputs; if they agree, RLE in item 43
   becomes a **port from PyDESeq2** rather than clean room (update that
   plan's tier when implementing). The C lives in `src/size_factors.c`
   (item 43); whichever item lands first writes it.
2. **Initial means and rough dispersions** (`_fit_MoM_dispersions()`,
   `fit_rough_dispersions()`, `fit_moments_dispersions()`): normalized
   counts \eqn{K/s}; rough \eqn{\hat\alpha_{rde}} from a least-squares fit
   of the normalized counts on the design (fitted values floored at 1),
   moments \eqn{\hat\alpha_{mde} = (v - \bar{s^{-1}} m)/m^2}; start =
   \eqn{\min}, clipped to \eqn{[10^{-8}, \max(10, n)]}. Initial \eqn{\mu}:
   linear fit (`fit_lin_mu()`, floored at `min_mu = 0.5`) when the design
   has as many distinct rows as columns, else IRLS with the rough
   dispersions.
3. **Gene-wise dispersions** (`fit_alpha_mle()`): maximise the
   Cox-Reid-adjusted NB likelihood in \eqn{\log\alpha}: \eqn{\ell(\alpha) -
   \tfrac12 \log\det(X^\top W X)}, \eqn{W = \mu/(1 + \mu\alpha)}, within
   \eqn{[\log 10^{-8}, \log\max(10, n)]}, then clip. PyDESeq2 uses SciPy's
   L-BFGS-B (grid search on failure); DESeq2 a backtracking line search
   (`dispTol = 1e-6`, grid fallback). **ecodive solves the 1-D problem to
   full precision** (safeguarded Newton on \eqn{\log\alpha} with the
   analytic first derivative given in `dloss()` and a numerically checked
   second derivative, bracketed; grid search over 100 points as PyDESeq2's
   `grid_fit_alpha()` when Newton fails). Both references stop short of
   the optimum by their tolerances, so ecodive's estimates differ from
   either at about 1e-6–1e-4 relative, as linda's REML fits differ from
   lme4's; tests check that ecodive's objective is never worse.
4. **Trend** (`_fit_parametric_dispersion_trend()`): \eqn{\alpha_{tr}(\bar
   m) = a_0 + a_1/\bar m} (\eqn{\bar m} = mean normalized count), a
   **gamma GLM with identity link** of the gene-wise estimates on
   \eqn{1/\bar m}, refitted after dropping features with ratio to the fit
   \eqn{< 10^{-4}} or \eqn{\ge 15}, until \eqn{\sum (\log(a/a_{old}))^2 <
   10^{-6}}; non-convergence or a coefficient \eqn{\le 10^{-10}} falls back
   to `'mean'` with a warning (PyDESeq2 and DESeq2 both do). PyDESeq2
   minimises the gamma deviance with L-BFGS-B; DESeq2 calls `glm(family =
   Gamma(link = "identity"))`. ecodive: Fisher scoring to convergence
   (2 × 2), same optimum. The dropping rule is discrete: a feature at the
   boundary can flip between implementations; note it, test away from it.
   `'mean'`: the 0.1%-trimmed mean of the gene-wise estimates above
   \eqn{10^{-7}} (`_fit_mean_dispersion_trend()`; PyDESeq2's
   `trim_mean(proportiontocut = 0.001)`; ✱ DESeq2's trimming rule from
   outputs).
5. **Prior variance** (`fit_dispersion_prior()`): residuals \eqn{r_j =
   \log\hat\alpha_j - \log\alpha_{tr}} over features with \eqn{\hat\alpha
   \ge 10^{-6}}; \eqn{\sigma^2_r} = (MAD of \eqn{r}, scaled by
   \eqn{1/\Phi^{-1}(0.75)})²; prior variance \eqn{\max(\sigma^2_r -
   \psi_1((n-p)/2), 0.25)}. PyDESeq2 only warns when \eqn{n - p \le 3};
   ✱ check from DESeq2's documentation and outputs whether DESeq2 uses
   another estimate there (if so, follow DESeq2 and say so).
6. **MAP dispersions** (`fit_MAP_dispersions()`): step 3's objective plus
   the normal prior \eqn{-(\log\alpha - \log\alpha_{tr})^2/(2\sigma^2)},
   started at the trend. **Dispersion outliers**: \eqn{\log\hat\alpha >
   \log\alpha_{tr} + 2\sigma_r} keep the gene-wise estimate.
7. **Coefficients** (`irls_solver()`): IRLS for the NB GLM with log link,
   ridge \eqn{\lambda = 10^{-6}} on every coefficient, \eqn{\mu} floored at
   0.5 inside the iterations, stop when \eqn{|d - d_{old}|/(|d| + 0.1) <
   10^{-8}} (\eqn{d = -2\ell}); start from the QR solution of \eqn{\log(K/s
   + 0.1)}; any \eqn{|\beta| > 30} or 250 iterations → bounded
   quasi-Newton on the penalised likelihood (PyDESeq2: L-BFGS-B, then a
   grid for \eqn{p \le 2}). ✱ DESeq2's ridge is on the log2 scale
   (\eqn{\lambda/\log(2)^2} in natural units?): settle from outputs on a
   tiny-count feature where it matters; the difference is ~1e-6 otherwise.
   Returns \eqn{\beta}, the hat diagonals \eqn{h_{ii}} (for Cook's) and
   convergence. Reuse `src/feature_lm.c`'s QR for each weighted solve
   (\eqn{\sqrt{W}X}), which is more stable than PyDESeq2's normal
   equations; agreement stays at rounding when converged.
8. **Wald test** (`wald_test()`): \eqn{M = X^\top W X}, \eqn{H = (M +
   \Lambda)^{-1}}, SE of \eqn{c^\top\beta} = \eqn{\sqrt{c^\top H M H c}}
   (the sandwich with the ridge, as both references), \eqn{z = c^\top\beta
   / SE}, two-sided normal p. Reported in log2: `log2_fc` = \eqn{c^\top
   \beta/\log 2}, `lfc_se` = SE / \eqn{\log 2}.
9. **Cook's distances** (`calculate_cooks()`): \eqn{D_{ij} = r_{ij}^2/(p\,
   \tau_j) \cdot h_{ii}/(1 - h_{ii})^2} with Pearson residuals against a
   robust method-of-moments dispersion (`robust_method_of_moments_disp()`:
   trimmed variances per design cell, floor 0.04). Cutoff \eqn{F_{p,
   n-p}(0.99)}.
10. **Outlier replacement and refit** (`_replace_outliers()`,
    `_refit_without_outliers()`): in samples whose design cell has at least
    `min_replicates` samples, counts with \eqn{D > } cutoff are replaced by
    the feature's 20%-trimmed mean of normalized counts times \eqn{s_i}
    (truncated to an integer); affected features are refitted from step 3
    with the **original** trend and prior variance (steps 3, 6, 7); a
    feature left all-zero gets LFC 0, p 1.
11. **Cook's filtering** (`cooks_outlier()`): among samples in cells with
    \eqn{\ge 3} replicates, a feature whose largest \eqn{D} exceeds the
    cutoff (after replacement: the replaced samples' \eqn{D} set to 0) gets
    p `NA`, unless at least 3 samples have counts above that sample's.
12. **Independent filtering** (`_independent_filtering()`): 50 quantiles
    \eqn{\theta} of the base mean from the fraction of zero base means to
    0.95; at each, BH over features with base mean \eqn{\ge} the cutoff;
    count rejections at `alpha`; fit `lowess(theta, rejections, f = 1/5)`;
    choose the first \eqn{\theta} whose count exceeds the lowess maximum
    minus the RMS of residuals (\eqn{j = 0} when the maximum count \eqn{\le
    10}). Features below the chosen cutoff get `p_adjusted = NA`.
    **PyDESeq2's `lowess()` is not R's**: R's `lowess()` (Cleveland 1979;
    used by DESeq2 through genefilter) has `delta` interpolation, its own
    neighbour count and iteration rules. Write Cleveland's LOWESS from the
    1979 paper and `?lowess` (R's code is GPL: don't read it), match R's
    outputs on the 50-point inputs; ✱ record the rules settled. Only the
    chosen index \eqn{j} matters, so a near-tie can differ; test cases
    away from ties. With `p_adjust` other than BH, filter with that method
    the same way (documented as ecodive's extension).
13. `reject` = `p_adjusted <= alpha` (DESeq2's `padj < alpha`; ✱ the
    inequality, which only matters at exact ties).

**Phase 2: likelihood-ratio test** (clean room: DESeq2's `nbinomLRT` Rd,
Love et al. 2014 Methods, outputs). Fit the reduced design with the
**full model's** final dispersions; statistic \eqn{2(\ell_{full} -
\ell_{reduced})}, \eqn{\chi^2} with \eqn{p_{full} - p_{red}} df. ✱ Settle
from outputs: which dispersions the reduced fit uses, how Cook's filtering
and refitting apply, and what `log2_fc` / `lfc_se` DESeq2 reports with LRT
(the full model's chosen coefficient). Result: one row per feature, the
`covariate` column naming the terms dropped.

**Phase 3:**

- **`shrink = 'apeglm'`**: port PyDESeq2's `lfc_shrink()`,
  `_fit_prior_var()`, `nbinomFn()` / `nbinomGLM()` (MIT): a Cauchy prior
  of scale \eqn{\min(\sqrt{\hat A}, 1)} on one coefficient (\eqn{\hat A}
  from the MLEs and SEs, root-found in \eqn{[10^{-6}, 400]}), a wide normal
  (sd 15) on the others, MAP by quasi-Newton from PyDESeq2's start, SE from
  the inverse Hessian; p-values unchanged. It follows apeglm's model (Zhu
  et al. 2019) but is PyDESeq2's own code; match PyDESeq2 at optimiser
  tolerance and compare with `lfcShrink(type = "apeglm")` outputs (apeglm
  GPL: output-only), documenting differences. **Not** `type = "normal"`
  (DESeq2's original prior; not in PyDESeq2) nor `"ashr"` (GPL, no
  permissive code). Shrunk LFCs go in a `log2_fc_shrunk` column (the
  test stays on the MLE, as DESeq2 recommends).
- **`lfc_threshold`, `alternative`**: port `wald_test()`'s four
  alternatives.
- **`norm = 'iterate'`**: port `_fit_iterate_size_factors()` (alternating
  dispersion fits and a Powell minimisation of the trimmed NB likelihood);
  ✱ match DESeq2's `"iterate"` outputs only at optimiser tolerance.

Not planned: VST / rlog transforms (DESeq2's `vst()` exists in PyDESeq2
and could later be a `norm` choice: ask), `fitType = "local"`,
glmGamPoi, time-course `test = "LRT"` helpers, `collapseReplicates`,
interaction-contrast shortcuts beyond a contrast vector.

### C, threads and determinism

`src/deseq2.c`, threaded over blocks of features with `run_parallel()`;
each feature is independent once \eqn{s_i}, the trend and the prior are
known, so the pipeline is three threaded passes with serial steps between:

1. **Pass A** (threads): per feature, read its counts
   (`da_table_init(norm = none)` / `da_table_column()`: the table by
   feature, never dense), base mean, rough and moments dispersions, initial
   \eqn{\mu} (computed on the fly from \eqn{\beta}, never stored as
   \eqn{n \times D}), gene-wise dispersion.
2. **Serial**: trend (Fisher scoring over \eqn{D} features; sums in
   fixed-size feature blocks added in block order, so independent of
   `cpus`), MAD, prior variance.
3. **Pass B** (threads): MAP dispersion, IRLS, Wald statistics for every
   requested coefficient or contrast, Cook's distances reduced on the fly
   to what steps 10–11 need per feature (the max over qualifying samples,
   its sample, the count of larger counts, and the list of samples over the
   cutoff in replaceable cells), so no \eqn{n \times D} matrix is kept.
4. **Pass C** (threads, only features with replacements): steps 3, 6, 7 on
   the replaced counts.
5. **Main thread**: independent filtering, BH, `pnorm()`/`pchisq()`.

Per-feature work uses thread-local buffers (\eqn{O(np + p^2)}); all
reductions are per feature or in fixed block order, so results are
`identical()` for any `cpus`. No random numbers (phase 3's iterate and
shrinkage are deterministic too).

### Result

`c('ecodive_deseq2', 'data.frame')`, one row per feature × covariate (or
per feature for a `contrast` or the LRT): `feature, covariate, base_mean,
log2_fc, lfc_se, stat, p_value, p_adjusted, reject`. Per-feature
diagnostics as attributes (they don't repeat per covariate): `dispersion`
(final), `disp_gene`, `disp_trend`, `disp_map`, `disp_outlier`,
`cooks_outlier`, `refitted`, `converged`. Other attributes: `method =
'deseq2'`, `formula`, `groups`, `contrast`, `design`, `norm`,
`size_factors` (named by sample), `trend` (`c(a0 = , a1 = )` or the mean),
`fit_type` (the one used, after any fallback), `prior_var`,
`filter_threshold` (base-mean cutoff chosen) and `filter_theta`, `alpha`,
`p_adjust`, `test`, `log_base = 2`, `mean` (per feature
\eqn{\log_2} base mean for MA plots, named as `da_lm()`'s). Printing:
`print_table_result`, header "DESeq2" (+ covariate).


## Reference implementations

- **Port from:** PyDESeq2 0.6.0 (MIT), the files and functions named
  above. Headers of `src/deseq2.c` and `R/deseq2.r`: the MIT text,
  "Copyright (c) 2022 Owkin"; `inst/COPYRIGHTS` stanza; `Authors@R`:
  "Owkin" (`cph`), and the PyDESeq2 authors as `ctb` (Muzellec et al.).
- **Provenance note (tell the user):** PyDESeq2 describes itself as "a
  reimplementation from scratch rather than a port", yet several comments
  say "as in the R code" / "to reproduce DESeq2's behaviour", so its
  authors read DESeq2 (LGPL). ecodive relies on PyDESeq2's MIT licence as
  published, as it relies on scikit-bio's BSD-3 for code written from
  ANCOMBC; the difference is that LGPL, unlike Artistic-2.0, has no
  notice to carry. **Ask the user** whether that is acceptable.
- **Output-only:** DESeq2 1.52.0, apeglm 1.34.0; R's `lowess()`,
  `glm(Gamma)`.
- **Fixtures** (`data-raw/expected/deseq2.r` + `deseq2.py`): DESeq2 and
  PyDESeq2 in scratch libraries/environments on: (a) a 12-sample,
  2-group, 500-feature negative-binomial simulation without zeros
  (`norm = 'rle'`); (b) a sparse 40 × 800 microbiome-like table
  (`'rle_poscounts'`, three groups and a numeric covariate); (c) a table
  with planted outliers in a 7-replicate cell (refit) and in a 3-replicate
  cell (filter only); (d) `fit_type = 'mean'`; (e) LRT `~ g + x` vs `~ x`
  (phase 2); (f) apeglm shrinkage (phase 3). Store size factors, the three
  dispersion vectors, trend coefficients, prior variance, MLE
  coefficients, SEs, statistics, Cook's maxima, `replace`/`refit` flags,
  p-values, `padj` and the filter threshold
  (`metadata(res)$filterThreshold`), as text or small `.rds`.
- **Exact checks through injection:** an internal
  `deseq2_details(dispersions = , size_factors = )` fixes those inputs, so
  coefficients, SEs, Wald p, Cook's and filtering compare at 1e-8 with the
  references' own dispersions, isolating the optimiser differences of
  step 3.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `plot.ecodive_deseq2(x, type = c('volcano', 'ma', 'forest',
  'dispersion'), covariate = NULL, ...)` via `da_plot()`
  ([da_plots.md](da_plots.md)): volcano (default), MA (x: `mean`, the log2
  base mean; DESeq2's `plotMA()`), forest (Wald intervals `log2_fc ±
  qnorm((1 + conf_level)/2) lfc_se`, `conf_level = 0.95` as
  `plot.ecodive_ancombc()`). No effect plot (no within-group dispersion
  on the same scale).
- **`'dispersion'`**: DESeq2's `plotDispEsts()`: gene-wise (black), trend
  (red line), final (blue) dispersions against the base mean, log-log, with
  dispersion outliers circled. Reads `disp_gene`, `disp_trend`,
  `dispersion`, `disp_outlier`, `mean`.
- Title: "Volcano plot of DESeq2" etc. (`da_plot_title()`), covariate in
  the subtitle; the input is a count table, so no `display`. Reads
  `alpha`, `log_base`.


## Tests

- Size factors equal DESeq2's at 1e-12 (both norms); trend, prior variance
  and every injected-input quantity (coefficients, SEs, statistics,
  p-values, Cook's) at 1e-8 against both references; dispersions at 1e-4
  relative, with ecodive's objective \eqn{\le} the reference's at the
  reference's estimate (to 1e-12); final p-values at 1e-4 relative;
  `reject` identical except features within 1e-4 of `alpha` (listed).
- Independent filtering: the threshold identical to DESeq2's on the
  fixtures; `lowess()` equal to R's at 1e-12 on random 50-point inputs.
- `identical()` for `cpus` 1, 2, 4, 8; dense = `Matrix` = `slam`;
  `margin = 2`; `da_test(x, 'deseq2', ...)` identical to the direct call.
- Edge cases: a feature all zero (`NA` row, base mean 0); a feature
  nonzero in one sample; no replicates (n = p: error, as both
  references); every feature with a zero and `norm = 'rle'` (error naming
  `'rle_poscounts'`); non-integer counts (error); n − p ≤ 3 (warning);
  a cell with < 3 replicates (no Cook's filtering there); trend
  non-convergence (falls back to `'mean'` with a warning); rank-deficient
  design (error).
- Smoke line.


## Docs

Rd page in *Differential abundance*:

- The model, size factors, dispersion shrinkage, Wald test, Cook's and
  independent filtering, in a few paragraphs; units (log2); that size
  factors make it a relative (compositional) analysis like the others,
  with no compositional correction; poscounts and why.
- **Benchmarks** (above) and a pointer to `da_sensitivity()` and the plain
  tests; that `alpha` drives the independent filter.
- **"Differences from DESeq2 and PyDESeq2"**: dispersions and the trend
  solved to full precision (both references stop at optimiser
  tolerances, about 1e-6–1e-4 relative); the default `norm`
  (poscounts); `alpha = 0.1` as DESeq2 (PyDESeq2: 0.05); no `'local'`
  trend, glmGamPoi, VST / rlog, normal or ashr shrinkage; the filter with
  non-BH `p_adjust`; every covariate reported at once; and each ✱ settled
  in DESeq2's favour where the two references differ (lowess, LRT,
  small-df prior variance).
- `parameter_names.md`: `fit_type`, `cooks`, `refit_outliers`,
  `min_replicates`, `independent_filter`, `reduced`, `lfc_threshold`,
  `shrink` (new); `test` here selects Wald or LRT (as `ancom(sig_test)`
  selects its test).
- NEWS.


## Phasing and effort

| Phase | Content | Hardest part |
|:------|:--------|:-------------|
| 1 | Steps 1–13: size factors, dispersions (gene-wise, trend, prior, MAP), IRLS, Wald, Cook's, refit, independent filtering; `DA_TESTS` row; plots | The dispersion MLE with Cox-Reid (a log-determinant and its derivative per evaluation, bounds, fallbacks); the trend's discrete refits; Cook's refit bookkeeping; R's lowess |
| 2 | LRT | Settling DESeq2's LRT conventions from outputs only |
| 3 | apeglm-style shrinkage, `lfc_threshold` / `alternative`, `norm = 'iterate'` | Non-convex MAP with a Cauchy prior; Powell-based iterate size factors |

Overall effort: **high** (phase 1 alone is several times `linda()`).
Per-feature GLM work threads naturally and is where ecodive wins
(PyDESeq2 parallelises with joblib processes; DESeq2 with BiocParallel).


## References

- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of
  fold change and dispersion for RNA-seq data with DESeq2. *Genome
  Biology*, 15, 550. doi:10.1186/s13059-014-0550-8 (74,537 OpenAlex
  citations, all fields; 10,234 Europe PMC microbiome mentions 2022–26;
  retrieved 2026-11-20).
- Muzellec, B., Teleńczuk, M., Cabeli, V., & Andreux, M. (2023).
  PyDESeq2: a python package for bulk RNA-seq differential expression
  analysis. *Bioinformatics*, 39(9), btad547.
  doi:10.1093/bioinformatics/btad547
- Anders, S., & Huber, W. (2010). Differential expression analysis for
  sequence count data. *Genome Biology*, 11, R106.
  doi:10.1186/gb-2010-11-10-r106
- McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential
  expression analysis of multifactor RNA-Seq experiments with respect to
  biological variation. *Nucleic Acids Research*, 40, 4288–4297.
  doi:10.1093/nar/gks042 (Cox-Reid adjusted profile likelihood).
- Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior
  distributions for sequence count data: removing the noise and preserving
  large differences. *Bioinformatics*, 35, 2084–2092.
  doi:10.1093/bioinformatics/bty895
- Bourgon, R., Gentleman, R., & Huber, W. (2010). Independent filtering
  increases detection power for high-throughput experiments. *PNAS*, 107,
  9546–9551. doi:10.1073/pnas.0914005107
- Cleveland, W. S. (1979). Robust locally weighted regression and
  smoothing scatterplots. *Journal of the American Statistical
  Association*, 74, 829–836. doi:10.1080/01621459.1979.10481038
- McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why
  rarefying microbiome data is inadmissible. *PLoS Computational Biology*,
  10, e1003531. doi:10.1371/journal.pcbi.1003531 (1,338 citations).
- Weiss, S., et al. (2017). Normalization and microbial differential
  abundance strategies depend upon data characteristics. *Microbiome*, 5,
  27. doi:10.1186/s40168-017-0237-y
- Hawinkel, S., Mattiello, F., Bijnens, L., & Thas, O. (2019). A broken
  promise: microbiome differential abundance methods do not control the
  false discovery rate. *Briefings in Bioinformatics*, 20, 210–221.
  doi:10.1093/bib/bbx104
- Nearing, J. T., et al. (2022). *Nature Communications*, 13, 342.
  doi:10.1038/s41467-022-28034-z
- Wirbel, J., et al. (2024). *Genome Biology*, 25, 247.
  doi:10.1186/s13059-024-03390-9
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130
