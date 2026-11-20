# `rademu()`: fold-difference estimation with robust score tests (radEmu)

**Priority:** 67. **Tier:** port (radEmu 2.3.2.0, CRAN, **MIT**:
`DESCRIPTION` "MIT + file LICENSE", `LICENSE` "YEAR: 2023 / COPYRIGHT
HOLDER: radEmu authors"; verified 2026-11-20; pure R, no compiled code,
`Depends: MASS, Matrix, methods, magrittr, rlang`). Read
[README.md](README.md) first, and [compositional.md](compositional.md).


## What it is

radEmu (Clausen, Teichman & Willis 2026) estimates, for each taxon \eqn{j}
and covariate \eqn{k}, a **log fold-difference** relative to the typical
taxon, with no pseudocount and no reference taxon, and tests it with a
**robust score test** that keeps error rates controlled with few samples
or misspecified models. Its mean model is
\eqn{E[Y_{ij}] = \exp(z_i + x_i^\top\beta_j)} (Poisson-likelihood
estimating equations; \eqn{z_i} absorbs each sample's depth), identified by
a constraint on each row \eqn{\beta_{k\cdot}}: by default its **pseudo-Huber
smoothed median** is zero (`constraint_param = 0.1`). "An estimate of 1
... taxon X is exp(1) = 2.72 times more abundant ... when compared to
typical fold-differences in abundances of taxa" (`emuFit()` docs). It is
the Willis lab's successor to corncob (GPL).

**Why it is last (priority 67):**

- **Low uptake:** 15 OpenAlex / 3 Crossref citations, 6 Europe PMC
  mentions 2023–26 (retrieved 2026-11-19; OpenAlex includes the arXiv
  version).
- **High complexity:** Firth-penalised fitting by iterated data
  augmentation, a sparse \eqn{pJ \times pJ} information matrix (radEmu uses
  `Matrix`'s CHOLMOD), and one constrained null refit of the whole model
  **per score test** (\eqn{(p - 1)J} refits by default), with several
  stopping rules (`tolerance`, `tol_lik`, `tol_test_stat`,
  `tol_discrete`, `null_window`). Parity is therefore a tolerance
  exercise, not bit-matching.
- **Moving target:** a scalability follow-up ("Scalable fold-change
  estimation from microbial metagenomics", bioRxiv 2026,
  doi:10.64898/2026.10.08.757531; not read) may change the algorithm.
  Re-check before starting.
- Kept as future-proofing because it is MIT (portable), principled
  (no pseudocounts or reference taxa), and its robust score test answers
  the benchmarks' FDR concerns. The review (#12) rated it "High; future-
  proofing only".


## Names

`rademu()`: the package name lowercased (rule 1; `emuFit()` is its
function, but users know "radEmu"). `DA_TESTS` row: design `formula`,
effect `estimate`, base e, `zero_free_ok = FALSE` (zeros are information
here); aliases `rademu`, `emu`, `emufit`. "Same names in other packages":
radEmu's `emuFit()` (the original).


## Design

### Signature

```r
rademu(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
       subjects = NULL, penalize = TRUE, score_tests = TRUE,
       constraint_param = 0.1, tol = 1e-4, max_iters = 1000L,
       max_abs_coef = 250, conf_level = 0.95, alpha = 0.05,
       p_adjust = 'BH', cpus = n_cpus())
```

- Model convention of [da_test.md](da_test.md); design from `da_design()`.
- `subjects`: radEmu's `cluster` (cluster-robust scores for repeated
  measures); the glossary's name for per-sample subject ids (`ctf()`).
- `score_tests`: radEmu's `run_score_tests`; `TRUE` tests every
  non-intercept (k, j); a logical or a data frame of (covariate, feature)
  pairs (its `test_kj`) limits the work.
- `tol`, `max_iters`, `max_abs_coef`: its `tolerance`, `maxit`,
  `max_abs_B`. `conf_level` = 1 − its `alpha` for CIs. The null-fit
  controls (`maxit_null`, `tol_lik`, `tol_test_stat`, `tol_discrete`,
  `null_window`, `control_fn()`) keep radEmu's defaults internally.
- `alpha`, `p_adjust`: radEmu reports raw score p-values only; ecodive
  adds `p_adjusted` (BH) and `reject`, as for every DA method.

### Algorithm (read from radEmu 2.3.2.0 `R/`; to be re-read when porting)

1. **Checks** (`emuFit_check()`): taxa with all-zero counts are an error;
   **zero-comparison** parameters (for a factor with ≥ 3 levels, a level
   vs the reference where the taxon is zero in every sample of both,
   `zero_comparison_check()`) are flagged, and their p-values set to `NA`
   when at most `remove_zero_comparison_pvals = 0.01` (keep that rule).
2. **Estimation** (`emuFit_micro_penalized()`): start from \eqn{Y + 10^{-3}
   \bar Y}; then repeat: Firth augmentations from the current \eqn{B}
   (`get_augmentations()`: Kosmidis & Firth 2011, via the information of
   the equivalent Poisson model and its Cholesky factor), and an ML fit to
   \eqn{Y +} augmentations (`emuFit_micro()`: per-taxon Fisher scoring with
   backtracking line search, Armijo `c1 = 1e-4`, `max_step`, with column
   \eqn{j_{ref}} (the taxon with the most non-zero samples, first in ties,
   `get_j_ref()`) held at 0, \eqn{z_i = \log\sum_j Y_{ij} - \log\sum_j
   \exp(x_i^\top\beta_j)}); until max |ΔB| < `tol`. Finally each row is
   shifted so its constraint holds: \eqn{\beta_{k\cdot} \leftarrow
   \beta_{k\cdot} - g(\beta_{k\cdot})}.
3. **Constraint** \eqn{g} = `pseudohuber_median(x, d = 0.1)`: start at the
   median, iterate \eqn{y \leftarrow} weighted mean with weights
   \eqn{(1 + ((x - y)/d)^2)^{-1/2}} until |Δy| < 1e-8; its gradient
   (`dpseudohuber_median_dx()`) enters the constrained fits and scores.
4. **Robust Wald CIs** (`micro_wald()`): sandwich covariance from the
   information \eqn{I} and the empirical score covariance \eqn{D_y}
   (scores summed within `subjects` if given), with the constraint's
   gradient; \eqn{\hat\beta \pm z_{1-\alpha/2}\sqrt{v}}; Wald p
   \eqn{\chi^2_1}.
5. **Robust score test** per (k, j) (`score_test()`): refit under
   \eqn{H_0: \beta_{kj} = 0} subject to all constraints (null algorithm:
   `"discrete"` when the design is only categorical with an intercept and
   \eqn{J < 150}, else `"constraint_sandwich"`; augmented Lagrangian
   options `rho_init = 1`, `tau = 2`, `kappa = 0.8`), then
   \eqn{T = (U^\top I^{-1}H)^2 / (H^\top I^{-1} D_y I^{-1} H)\cdot
   \frac{n}{n-1}} (or \eqn{\frac{m}{m-1}} for \eqn{m} clusters; Guo et al.
   2005), \eqn{p = P(\chi^2_1 > T)} (`get_score_stat()`).
6. radEmu calls `rnorm()` only to classify the constraint function's type
   (no effect on results); ecodive needs no RNG.

### C

- Dense per-taxon blocks: the information has an arrow structure
  (\eqn{p \times p} blocks per taxon, coupled through the \eqn{n} sample
  effects), so solve by Schur complement on the \eqn{z} block instead of a
  sparse Cholesky ecodive doesn't have. Work this out against radEmu's
  `f_info()` / `macro_fisher_null()` before writing code.
- Parallel unit: **the score test** (each (k, j) null refit is
  independent, warm-started from \eqn{\hat B}); results gathered in test
  order, so identical for any `cpus`. Within the full fit, taxa updates
  parallel in fixed blocks.

### Result

`c('ecodive_rademu', 'data.frame')`, one row per feature × non-intercept
covariate: `feature, covariate, estimate` (natural-log fold-difference),
`ci_lower, ci_upper, score_stat, p_value` (score), `p_adjusted, reject,
zero_comparison`. Attributes: `formula`, `groups`, `coef` (the p × J
\eqn{B}, intercept included), `z` (sample effects), `converged` (per
test, radEmu's `score_test_hyperparams$converged`), `penalize`,
`constraint_param`, `subjects`. `plot()`: volcano / forest via
[da_plots.md](da_plots.md).


## Reference implementations

- **Port from:** radEmu 2.3.2.0 (R only). Header and `inst/COPYRIGHTS`
  stanza with the MIT notice; `Authors@R` `cph` "radEmu authors".
- **Fixtures:** `data-raw/expected/rademu.r`: `emuFit()` on a small table
  (n = 20, J = 15), a two-level factor (discrete null algorithm), and a
  numeric covariate (constraint_sandwich); with and without `cluster`;
  `penalize = FALSE` once. Store `coef`, `B`, `z_hat`,
  `score_test_hyperparams`.


## Tests

- `estimate` at 1e-6 (relative to `tol`), CIs at 1e-5, score statistics
  and p at 1e-4; `zero_comparison` identical.
- Each row of \eqn{B} satisfies its constraint to 1e-8.
- Rescaling a sample leaves `estimate` unchanged (absorbed by \eqn{z_i})
  to the fitting tolerance.
- `identical()` for `cpus` 1, 2, 4.
- Edge cases: an all-zero taxon (error), a zero-comparison level,
  intercept-only formula (message, no tests), J = 2, non-convergence
  (warning, `converged` attribute).
- Smoke line (no score tests, small table).


## Docs

Rd page in *Differential abundance*: the model and its "relative to the
typical taxon" reading; why score tests (the manuscript's recommendation);
cost (a refit per test; how to limit `score_tests`); and "Differences from
radEmu" (adjusted p-values added; design coding as `ancombc()`). Glossary:
`penalize`, `score_tests`, `constraint_param`, `max_abs_coef` (new). NEWS.


## References

- Clausen, D. S., Teichman, S., & Willis, A. D. (2026). Estimating ratios
  of means of multicategory data observed with sample and category
  perturbations. *Biometrika*, 113(2). doi:10.1093/biomet/asag009 (15
  OpenAlex / 3 Crossref, retrieved 2026-11-19; arXiv:2402.05231).
- Kosmidis, I., & Firth, D. (2011). Multinomial logit bias reduction via
  the Poisson log-linear model. *Biometrika*, 98, 755–759.
- Guo, X., Pan, W., Connett, J. E., Hannan, P. J., & French, S. A. (2005).
  Small-sample performance of the robust score test and its modifications
  in generalized estimating equations. *Statistics in Medicine*, 24,
  3479–3495. doi:10.1002/sim.2161.
- Martin, B. D., Witten, D., & Willis, A. D. (2020). Modeling microbial
  abundances and dysbiosis with beta-binomial regression. *Annals of
  Applied Statistics*, 14, 94–115. doi:10.1214/19-aoas1283 (corncob,
  GPL; context only).
