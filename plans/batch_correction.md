# Batch correction: `batch_mmuphin()`, `batch_combat()`, `batch_combat_seq()`, `batch_conqur()`

**Priority:** 53 (`batch_mmuphin()`, `batch_combat()`) and 81 (`batch_combat_seq()`, `batch_conqur()`); see "Order". **Tiers:**

| Function | Tier | Source and licence (checked 2026-11-20) |
|:---------|:-----|:----------------------------------------|
| `batch_mmuphin()` | **port** | MMUPHin 2.0.0 (Bioconductor 3.23 release), `adjust_batch()`: **MIT** (`DESCRIPTION` "MIT + file LICENSE"; `LICENSE` "YEAR: 2019 / COPYRIGHT HOLDER: Harvard School of Public Health"). The release tarball's `R/adjust_batch.R` and `R/helpers_adjust_batch.R` are identical to GitHub `master` (1.99.3, last commit 2c6cd74, 2026-01-31). |
| `batch_combat()` | **port** | sva 3.60.0 (Bioconductor 3.23), `ComBat()`: **Artistic-2.0** (`DESCRIPTION`). |
| `batch_combat_seq()` | **port + clean room** | sva 3.60.0 `ComBat_seq()` (Artistic-2.0) for the method; the edgeR functions it calls (edgeR 4.10.7, **GPL**, Bioconductor page) are output-only and must be written from edgeR's Rd pages and papers. |
| `batch_conqur()` | **clean room** | ConQuR 2.0 (GitHub wdl2459/ConQuR): **GPL-3.0** (`DESCRIPTION`, `LICENSE`). Only the paper (Nature Communications, CC-BY), the README, the vignette and the reference manual `ConQuR_2.0.pdf` may be read; its dependency quantreg 6.1 is **GPL (>= 2)**. |

Python ports are no way around this: pyComBat (`combat` 0.3.3 on PyPI) and
inmoose 0.9.1 are both **GPL-3.0-or-later** (PyPI classifiers / licence
expression, checked 2026-11-20). Read [README.md](README.md) first; the
Artistic-2.0 rules there (keep the notice, say how the port differs) apply
to `batch_combat()` and `batch_combat_seq()`.


## What it is

Multi-study (or multi-run, multi-extraction-kit) microbiome tables carry
**batch effects**: systematic per-batch shifts in feature abundance that
can swamp biology, and that make pooled analyses and classifiers learn the
batch. Batch correction estimates each batch's location (and scale)
effect per feature, **while preserving the effects of named covariates**
(disease, age), and removes it, returning a corrected table for the usual
downstream analyses (ordination, PERMANOVA, DA, prediction).

Demand (retrieved 2026-11-20; OpenAlex citations 2022–26, Europe PMC
microbiome full-text mentions 2022–26):

| Method | Citations | Mentions |
|:-------|:----------|:---------|
| MMUPHin (Ma et al. 2022) | 192 (8 in 2022 → 63 in 2026) | "MMUPHin": 134 |
| ConQuR (Ling et al. 2022) | 167 (→ 67 in 2026) | "ConQuR": 83 |
| ComBat-seq (Zhang et al. 2020) | 1,786 (all fields; mostly RNA-seq) | |
| any | | "batch effect" AND (ComBat OR MMUPHin OR ConQuR): 490 |

Wang & Lê Cao (2020) review the microbiome case: methods built for
microarrays and RNA-seq (ComBat, ComBat-seq) ignore zero inflation, which
is why MMUPHin and ConQuR exist.


## Order (recommended)

1. **`batch_mmuphin()` first.** MIT (readable, portable without the
   Artistic-2.0 obligations), the most mentioned microbiome-specific method,
   **deterministic** (closed-form per-feature least squares plus a short
   empirical-Bayes fixed-point iteration; no random numbers), and low to
   medium effort: about 600 lines of R to port, of which maybe 300 lines of
   C. It also yields ComBat's machinery.
2. **`batch_combat()`** right after: sva's parametric and non-parametric
   ComBat on a continuous table (log-scale metabolomics, clr values,
   log-transformed counts). It shares step 3–5's empirical-Bayes code with
   `batch_mmuphin()` (MMUPHin's `zero_inflation = FALSE` is "ComBat with
   parametric adjustment" on log2 TSS, with the differences listed below),
   so it is low effort.
3. **`batch_combat_seq()`**: Artistic-2.0, so its own code ports, but it
   rests on edgeR's dispersion estimation and NB GLM fitting, which are GPL
   and must be re-derived from documentation, then matched to edgeR's
   outputs. Medium-to-high effort and the main risk is exact parity of the
   dispersions. High citation count, but mostly outside microbiome work,
   and it is not designed for zero-inflated tables.
4. **`batch_conqur()`** last: the fastest-rising microbiome method, but GPL
   (clean room), and it needs a quantile-regression LP solver (Barrodale–
   Roberts simplex, Koenker & D'Orey 1987) and, for its penalised variants,
   a lasso logistic regression with CV and a quantile lasso. High effort.

**Deferred, with reasons:**
- **MMUPHin's meta-analysis of DA**: in MMUPHin 2.0 it is `maaslin_meta()`
  (it was `lm_meta()` on MaAsLin2 before 2.0; the user's note says
  `lm_meta`). It runs maaslin3 per batch, then a per-feature random-effects
  meta-analysis with `metafor::rma(method = 'REML')` (metafor 5.2-1, **GPL
  (>= 2)**: output-only). It should wait for `maaslin3()`
  ([maaslin3.md](maaslin3.md)) and become its own small plan,
  `maaslin3_meta()`: `maaslin3()` per batch (threaded over batches) plus a
  clean-room REML random-effects pooling (Viechtbauer 2005: a Fisher-scoring
  iteration on τ²), matched to `metafor::rma()` outputs. Low-to-medium
  effort once `maaslin3()` exists.
- **`discrete_discover()`** (per-batch `fpc::claraCBI` clustering with
  prediction strength; fpc GPL) and **`continuous_discover()`** (per-batch
  PCA loadings joined into a network and clustered with
  `igraph::cluster_optimal()`, which needs the GLPK integer-programming
  library): low uptake on their own, heavy machinery, and
  `enterotype()` / the planned `dmm()` already cover community typing.
  Not planned; revisit if users ask.
- Other methods (percentile normalisation, Gibbons et al. 2018; PLSDA-batch,
  Wang & Lê Cao 2023; RUV): not planned now.


## Names

- **A `batch_` family**, so the methods sort together, as `da_*` and
  `replace_zeros_*` do: `batch_mmuphin()`, `batch_combat()`,
  `batch_combat_seq()`, `batch_conqur()`.
- Why not the bare names (rule 1): `combat()` and `combat_seq()` would
  read as sva's functions, whose arguments and orientation differ (sva
  takes features × samples and `mod` matrices); `mmuphin()` is a package
  with four unrelated functions, not a method; MMUPHin's own
  `adjust_batch()` would mask it with a different signature.
- **A wrapper `batch_correct(counts, method, batch, formula, data, ...)`**
  with a `BATCH_CORRECTIONS` registry (`list_batch_corrections()`,
  `match_batch_correction()`), added with the second method, as
  `replace_zeros()` was. Aliases: `mmuphin`, `zicombat`; `combat`,
  `combatseq`. Members must not take `method` or any prefix of `method` or
  `counts` (add the row to README's wrapper table); `batch`, `formula` and
  `data` pass through, as `da_test()` passes `groups`. **Ask the user** to
  confirm the family prefix and the wrapper.
- **Same names in other packages:** none.


## Shared design

### Common arguments

```r
batch_mmuphin(counts, batch, formula = NULL, data = NULL, margin = 1L,
              zero_inflation = TRUE, pseudocount = NULL, tol = 1e-4,
              max_iters = 1000L, cpus = n_cpus())

batch_combat(counts, batch, formula = NULL, data = NULL, margin = 1L,
             norm = 'none', pseudocount = NULL, parametric = TRUE,
             mean_only = FALSE, ref_batch = NULL, cpus = n_cpus())

batch_combat_seq(counts, batch, formula = NULL, data = NULL, margin = 1L,
                 shrink = FALSE, shrink_disp = FALSE, n_shrink_features = NULL,
                 seed = 0L, cpus = n_cpus())

batch_conqur(counts, batch, ref_batch, formula = NULL, data = NULL,
             margin = 1L, logistic_lasso = FALSE, quantile_type = 'standard',
             simple_match = FALSE, lambda_quantile = '2p/n',
             interpolate = FALSE, delta = 0.4999,
             taus = seq(0.005, 0.995, by = 0.005), n_folds = 10L,
             seed = 0L, cpus = n_cpus())
```

- **`batch`:** one value per sample, or a column name of `data`
  (`per_sample()`, as `groups`); at least two batches, each with at least
  two samples (ComBat-seq's and MMUPHin's rule).
- **`formula`:** a one-sided formula over `data` naming the covariates
  whose effects are **kept** (MMUPHin's `covariates`, sva's `mod`,
  ComBat-seq's `group` + `covar_mod`, ConQuR's `covariates`). The design is
  `model.matrix()` of it without the intercept, bound to the batch
  indicators. Not `given`: the glossary's `given` is partialled *out*.
  Confounding with batch is an error, with MMUPHin's and sva's messages
  ("Covariates and batch are confounded!").
- **`tol`, `max_iters`:** MMUPHin's `conv` and `maxit` (glossary names).
  `ref_batch`: sva's `ref.batch` and ConQuR's `batch_ref` (`reference` is
  avoided as an argument name; parameter_names.md §5).
- **Output:** the corrected table, **same orientation and dimnames as the
  input** (`margin`), as a dense base matrix. Integer counts in give
  rounded counts out (MMUPHin's `back_transform_abd()`; ComBat-seq returns
  integers by construction). A sparse input is returned dense (ecodive has
  no Matrix dependency); MMUPHin and ComBat-seq keep zeros as zeros, so the
  docs show `Matrix::Matrix(x, sparse = TRUE)` to re-sparsify.
- **Result class:** `c('ecodive_batch', 'matrix', 'array')`, so `plot()`
  can draw the diagnostics; every ecodive function that takes `counts` must
  accept it (test it), and `[` drops the class and attributes as for any
  matrix. Attributes (README rule 1, results only): `method`;
  `batch_effects`, a data frame `feature`, `batch`, `gamma_hat`,
  `gamma_star`, `delta_hat`, `delta_star` (`NA` where a feature was not
  adjusted in a batch); `adjusted`, a logical per feature (MMUPHin's
  `ind_feature`); `pseudocount` (when one was chosen); and for the
  "means" plot, `batch_means` and `batch_means_adjusted` (features ×
  batches, mean relative abundance before and after) and
  `overall_means`. No `display`: count tables have no title.

### Threads and determinism

Every per-feature step (design subset, rank check, least squares,
standardisation, per-batch moments, adjustment) is independent: tasks are
blocks of features, written into per-feature slots. The empirical-Bayes
hyperparameters (per-batch means and variances over features) are reduced
in feature order on the main thread. The shrinkage iteration runs per
batch (one task per batch), sequentially over features inside, because its
stopping rule is a maximum over all features. No random numbers in
`batch_mmuphin()` or `batch_combat()`; results are `identical()` for any
`cpus`.


## `batch_mmuphin()` (port of MMUPHin 2.0.0 `adjust_batch()`)

MMUPHin takes features × samples; ecodive's `margin` handles orientation.
Steps, as read from `R/adjust_batch.R` and `R/helpers_adjust_batch.R`:

1. **Checks** (`check_feature_abd()`): no `NA`, no negatives; "proportions"
   if every value ≤ 1, "counts" if all whole, otherwise MMUPHin stops.
   ecodive also accepts other non-negative tables (not rounded on output);
   say so under "Differences".
2. **Transform:** `pseudocount = NULL` means half the smallest non-zero
   value of the **raw** table (`set_pseudo()`: 0.5 for counts); then
   per sample \eqn{\log_2} of TSS of (counts + pseudocount)
   (`normalize_features(TSS, pseudo_count)` then `transform_features(LOG)`).
3. **Which values enter** (`construct_ind()`): with `zero_inflation`, a
   cell whose raw value is 0 is left out of every fit and stays 0. A
   feature is adjusted in batch b only if, among its non-zero samples, at
   least two batches are present, the design (present batches + covariates)
   has full column rank (`qr()$rank`, tolerance 1e-7) and more rows than
   columns; batches present in fewer than 2 features are switched off;
   features adjusted in no batch are left unchanged. A batch adjusted in no
   feature is an error ("The following batch(es) either have no present
   features ...").
4. **Standardise per feature** (`standardize_feature()`): least squares on
   that feature's design; grand mean = mean over its samples of the batch
   part of the fit; pooled variance = `var()` (**n − 1**) of the residuals,
   1 if all-equal; \eqn{s = (y - \text{grand mean} - \text{covariate part})
   / \sqrt{\text{var}}}.
5. **Per-batch moments** (`fit_EB()`): \eqn{\hat\gamma_{fb}} = mean of
   \eqn{s} in batch b; \eqn{\hat\delta_{fb}} = **`sd()`** of \eqn{s} in
   batch b (1 if `NA` or 0). Priors per batch over features: \eqn{\bar\gamma},
   \eqn{\tau^2 = \mathrm{var}(\hat\gamma)}, and `aprior()` / `bprior()` of
   \eqn{\hat\delta} (method of moments, with a zero variance replaced by 1).
6. **Shrinkage** (`it_sol()`, per batch, over the features adjusted in it):
   iterate \eqn{\gamma^* = (\tau^2 n \hat\gamma + \delta^* \bar\gamma) /
   (\tau^2 n + \delta^*)}, \eqn{\delta^* = (\tfrac12 \sum (s - \gamma^*)^2 + b)
   / (n/2 + a - 1)} (a zero sum of squares is `NA`), from \eqn{\gamma^* =
   \hat\gamma}, \eqn{\delta^* = \hat\delta}, until \eqn{\max(|\Delta\gamma| /
   \gamma_{old}, |\Delta\delta| / \delta_{old}) \le} `tol` (note: the
   divisor is signed, so negative ratios never stop the loop early; keep
   this) or error after `max_iters`.
7. **Adjust** (`relocate_scale()`, `add_back_covariates()`): \eqn{(s -
   \gamma^*)/\sqrt{\delta^*}}, then × \eqn{\sqrt{\text{var}}} + the
   feature's standardisation mean.
8. **Back-transform** (`back_transform_abd()`): \eqn{2^{x}}, original
   zeros back to 0, TSS per sample, × the sample's original total, rounded
   if the input was counts.

**Reference quirks to decide with the user (recommend: replicate, and
document under "Differences from ComBat"):**
- step 5 uses the **standard deviation** as \eqn{\hat\delta} where ComBat
  uses the variance, so the inverse-gamma prior is fitted to SDs;
- step 4 uses `var()` (n − 1) where ComBat divides by n;
- the pseudocount is added before TSS, so non-zero values move slightly.

`zero_inflation = FALSE`: every cell enters every fit (MMUPHin's ComBat
mode); the quirks above still apply, which is why `batch_combat()` exists
separately.


## `batch_combat()` (port of sva 3.60.0 `ComBat()`)

For a continuous table (`norm = 'none'`: values used as given, e.g. log
metabolite intensities or clr values; `norm = 'clr'` or `'log'` computes
them first, with `pseudocount`). Port `ComBat()` and `helper.R`'s
`aprior`, `bprior`, `postmean`, `postvar`, `it.sol`, `int.eprior`:

- standardisation with `B.hat` from the full design, grand mean weighted
  by batch sizes (or the reference batch's coefficients with `ref_batch`),
  pooled variance with divisor n (or the reference batch's);
- \eqn{\hat\gamma} by least squares on the batch design, \eqn{\hat\delta} the
  per-batch **variance** (`rowVars`);
- `parametric = TRUE`: `it.sol()` with `conv = 1e-4` (fixed in sva);
  `FALSE`: `int.eprior()`, the non-parametric Monte-Carlo-free integral
  over the other features (O(features² × samples) per batch: threaded over
  features, each sum in feature order);
- `mean_only` (also forced when a batch has one sample, with sva's
  message);
- `ref_batch`: that batch's data are returned unchanged.
- Not ported: `prior.plots` (a `plot()` type instead), `BPPARAM`, `NA`
  handling (`Beta.NA`): `NA` is an error.
- File header: "C port of ComBat() from sva 3.60.0 (Artistic-2.0);
  restructured for ecodive, see plans/batch_correction.md", the Artistic
  text, and `inst/COPYRIGHTS`.


## `batch_combat_seq()` (port of sva 3.60.0 `ComBat_seq()`)

ComBat-seq fits a negative-binomial GLM per gene with batch and covariates,
removes the batch coefficients from the fitted means, pools the per-batch
dispersions, and maps each count to the batch-free NB distribution by
**quantile matching** (counts ≤ 1 unchanged; `1 + qnbinom(pnbinom(y - 1,
old), new)`; a tail probability within 1e-4 of 1 keeps the count). Read
from `R/ComBat_seq.R` and `R/helper_seq.R`:

1. Drop genes that are all-zero in any batch (returned unchanged).
2. Design = batch indicators + `model.matrix(~ group)` (if `full_mod` and
   `group` has ≥ 2 levels) + covariates, intercept columns dropped;
   ComBat-seq's confounding errors. In ecodive, `formula` names `group`
   and the covariates together (its column space is the same, so the
   batch coefficients are too); `formula = NULL` is ComBat-seq's null
   model. `full_mod` is not an argument.
3. Per batch: common dispersion (`edgeR::estimateGLMCommonDisp`, design
   = the covariate model on that batch's samples, or `NULL` when there are
   too few residual df); gene-wise dispersions
   (`estimateGLMTagwiseDisp(prior.df = 0)`), or the common one when df are
   too few.
4. `edgeR::glmFit(prior.count = 1e-4)` with offset log library size;
   \eqn{\alpha_g} = batch coefficients weighted by batch sizes; refit with
   offset + \eqn{\alpha_g}; \eqn{\hat\gamma} = batch coefficients,
   \eqn{\hat\mu} = fitted values.
5. `shrink = TRUE`: per batch, a leave-one-out likelihood-weighted average
   over other genes (`monte_carlo_int_NB`, O(G² n)); `n_shrink_features`
   (sva's `gene.subset.n`) samples that many genes with R's `sample()`:
   here PCG with `seed` (a seed-like difference only).
6. \eqn{\mu^* = \hat\mu e^{-\gamma^*}}, \eqn{\phi^*} = mean of the per-batch
   dispersions; quantile matching.

**The edgeR parts are clean room.** Write them from edgeR's Rd pages
(`estimateGLMCommonDisp`, `dispCoxReid`, `estimateGLMTagwiseDisp`,
`dispCoxReidInterpolateTagwise`, `adjustedProfileLik`,
`maximizeInterpolant`, `glmFit`, `addPriorCount`, `getOffset`) and
McCarthy, Chen & Smyth (2012) / Robinson & Smyth (2007): the Cox–Reid
adjusted profile likelihood, the grid of `grid.npts = 11` log2-dispersion
offsets over `grid.range = c(-6, 6)` around the common dispersion, the
interpolating spline's maximum, and how `prior.count` shrinks the
coefficients. Things the documentation may not settle (the spline type,
the scale and interval of the common-dispersion optimisation, the
`min.row.sum` filter, the GLM convergence rule) are settled **from edgeR's
outputs** on designed inputs, and recorded in this plan. Where R's own
`stats::optimize()` is the optimiser, call it from R (calling R is not
reading GPL code). `pnbinom()` / `qnbinom()` run in R, vectorised (R's API
is not thread-safe, as `src/ancom.c` notes).

Parity target: \eqn{\hat\gamma}, \eqn{\hat\mu} and the dispersions at about
1e-6 (iterative fits), and the adjusted counts **exactly** except where a
tail probability falls within that tolerance of a quantile boundary;
report how many cells differ in the fixture comparison.

Threads: per gene (dispersion profiles, GLM fits, quantile inputs); the
shrinkage integral per (batch, gene).


## `batch_conqur()` (clean room from Ling et al. 2022)

Defaults from the reference manual (`ConQuR_2.0.pdf`, usage of `ConQuR()`):
`logistic_lasso = F`, `quantile_type = "standard"` (also `"lasso"`,
`"composite"`), `simple_match = F`, `lambda_quantile = "2p/n"`,
`interplt = F`, `delta = 0.4999`, `taus = seq(0.005, 0.995, by = 0.005)`,
`num_core = 2`; `batch_ref` has no default. `ConQuR_libsize()` (library
size as a covariate) and `Tune_ConQuR()` exist too.

Per taxon, a two-part model with batch and covariates: a logistic
regression for presence (optionally lasso, λ by CV) and quantile
regressions of the count at each τ (optionally penalised); each sample's
observed count is located on its own fitted conditional distribution and
replaced by the value at the same percentile of the reference batch's
distribution (batch set to `ref_batch`, covariates kept); `simple_match`
skips the regressions and matches empirical quantiles. Settle from the
paper's Methods and supplement before coding, and record here: how a zero
is placed on the percentile scale (random within \eqn{[0, 1 - \hat\pi]},
or a fixed point?), the interpolation window (`interplt`, `delta`), the
quantile-crossing fix, and how CV folds are drawn (PCG with `seed`).
Needs: a Barrodale–Roberts simplex for quantile regression (Koenker &
D'Orey 1987; Koenker 2005), written from the papers; IRLS logistic
regression (exists in spirit in `maaslin3()`'s plan); coordinate-descent
lasso (shared with [coda_lasso.md](coda_lasso.md)). Threads per taxon.
Fixtures: ConQuR's outputs (GPL: run in a scratch library), compared on
PERMANOVA R² of batch before and after and per-taxon quantiles, since its
randomness (if any) is not ours. **Write a fuller plan for this method when
its turn comes**; this section records the scope only.


## Reference implementations and fixtures

- `data-raw/expected/batch_mmuphin.r`: `MMUPHin::adjust_batch()` (2.0.0)
  on MMUPHin's own `CRC_abd` / `CRC_meta` subset (5 studies; `batch =
  'studyID'`, `covariates = 'study_condition'`) and on a small simulated
  table with planted batch shifts, features absent from some batches,
  and a covariate; `zero_inflation` TRUE and FALSE; `control =
  list(verbose = FALSE, diagnostic_plot = NULL)`. Store the adjusted table
  and the internal \eqn{\hat\gamma}, \eqn{\gamma^*} (from
  `MMUPHin:::fit_EB()` / `fit_shrink()` called directly; reading MIT code to
  find them is fine). Match at 1e-10 before rounding (store the unrounded
  values from a proportions input too) and exactly after.
- `data-raw/expected/batch_combat.r`: `sva::ComBat()` parametric,
  non-parametric, `mean.only`, `ref.batch`, with a `mod`; 1e-10.
- `data-raw/expected/batch_combat_seq.r`: `sva::ComBat_seq()` default,
  with `group`, with `covar_mod`, `shrink = TRUE` (without
  `gene.subset.n`); plus edgeR's intermediate outputs
  (`estimateGLMCommonDisp`, `estimateGLMTagwiseDisp`, `glmFit`
  coefficients and fitted values) for the clean-room parts.
- ConQuR: see above.
- Fixtures under `tests/tinytests/expected/<name>.rds`, small.


## Plotting

`plot.ecodive_batch(x, type = c('shrinkage', 'means'), ...)`, MMUPHin's
diagnostic figure (`diagnostic_adjust_batch()`), in base graphics:

- `'shrinkage'`: \eqn{\hat\gamma} (x) against \eqn{\gamma^*} (y) per
  feature, coloured by batch (`plot_palette()`), with the 1:1 line: how much
  the empirical-Bayes step pulled each batch effect in.
- `'means'`: per feature and batch, the batch's mean relative abundance
  (y) against the feature's overall mean (x), before (black) and after
  (`PLOT_HIGHLIGHT`), with each feature's points joined by a line, and the
  1:1 line.
- Title: the method's display name from `BATCH_CORRECTIONS` ("MMUPHin
  Batch Correction", "ComBat", "ComBat-seq", "ConQuR"). Reads only the
  stored attributes. `batch_combat_seq()` stores \eqn{\hat\gamma} on the log
  scale of the NB means, so the shrinkage plot applies to it when `shrink =
  TRUE` (otherwise \eqn{\gamma^* = \hat\gamma}: say so in the subtitle).
- The real check of a correction is an ordination coloured by batch and
  PERMANOVA's R² for batch before and after: the docs show
  `permanova(bray(x), batch)` and `plot(pcoa(bray(x)), batch)` on both
  tables (MMUPHin's vignette does the same with vegan).


## Tests

- Fixtures above, per function.
- Invariants: a table with no batch effect (all batches drawn from one
  distribution) changes little (R² of batch near 0 before and after);
  permuting samples permutes the output; feature order doesn't matter;
  zeros stay zeros (`batch_mmuphin()`, `batch_combat_seq()`); the
  `ref_batch` is unchanged (`batch_combat()`).
- `identical()` for `cpus` 1, 2, 4; ComBat-seq with `n_shrink_features`
  identical for the same `seed`.
- The result passes through `bray()`, `pcoa()`, `permanova()`,
  `da_test()` unchanged in behaviour (class handled).
- Edge cases: one batch (error); a batch of one sample (error, or ComBat's
  `mean_only`); covariates confounded with batch (MMUPHin's / sva's
  errors); a feature present in one batch only (left unchanged); a sample
  of all zeros; `NA` (error).
- Smoke lines (small table, two batches).


## Docs

- A new pkgdown section "Batch correction" (with the wrapper's page).
- Each Rd: what is estimated and preserved; that covariates must be named
  in `formula` or their effects are removed if they differ between
  batches; that corrected tables are for exploration and pooled analyses,
  while DA across studies is better done with batch as a covariate (or
  meta-analysis) than on corrected values (MMUPHin's and Wang & Lê Cao's
  advice); compositional caveats (MMUPHin corrects log relative
  abundances and re-closes).
- "Differences from MMUPHin / sva / ConQuR": orientation (`margin`);
  `formula` instead of covariate names or `mod`; dense output; non-count
  tables accepted; PCG for ComBat-seq's subsampling; the replicated
  MMUPHin quirks; the clean-room edgeR parts and their measured agreement.
- `inst/COPYRIGHTS` and `Authors@R`: MMUPHin (MIT), sva (Artistic-2.0,
  with the modified-version statement).
- NEWS; parameter_names.md: `batch`, `ref_batch`, `zero_inflation`,
  `parametric`, `mean_only`, `shrink`, `shrink_disp`, `n_shrink_features`,
  ConQuR's arguments (`interpolate` for its `interplt`).


## References

- Ma, S., Shungin, D., Mallick, H., Schirmer, M., Nguyen, L. H., Kolde, R.,
  Franzosa, E., Vlamakis, H., Xavier, R., & Huttenhower, C. (2022).
  Population structure discovery in meta-analyzed microbial communities
  and inflammatory bowel disease using MMUPHin. *Genome Biology*, 23, 208.
  doi:10.1186/s13059-022-02753-4
- Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects
  in microarray expression data using empirical Bayes methods.
  *Biostatistics*, 8, 118–127. doi:10.1093/biostatistics/kxj037
- Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch
  effect adjustment for RNA-seq count data. *NAR Genomics and
  Bioinformatics*, 2, lqaa078. doi:10.1093/nargab/lqaa078
- Ling, W., Lu, J., Zhao, N., Lulla, A., Plantinga, A. M., Fu, W., Zhang,
  A., Liu, H., Song, H., Li, Z., Chen, J., Randolph, T. W., Koay, W. L. A.,
  White, J. R., Launer, L. J., Fodor, A. A., Meyer, K. A., & Wu, M. C.
  (2022). Batch effects removal for microbiome data via conditional
  quantile regression. *Nature Communications*, 13, 5418.
  doi:10.1038/s41467-022-33071-9
- McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential
  expression analysis of multifactor RNA-Seq experiments with respect to
  biological variation. *Nucleic Acids Research*, 40, 4288–4297.
  doi:10.1093/nar/gks042
- Robinson, M. D., & Smyth, G. K. (2007). Moderated statistical tests for
  assessing differences in tag abundance. *Bioinformatics*, 23,
  2881–2887. doi:10.1093/bioinformatics/btm453
- Cox, D. R., & Reid, N. (1987). Parameter orthogonality and approximate
  conditional inference. *JRSS B*, 49, 1–39.
  doi:10.1111/j.2517-6161.1987.tb01422.x
- Koenker, R., & D'Orey, V. (1987). Computing regression quantiles.
  *Applied Statistics*, 36, 383–393. doi:10.2307/2347802
- Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance
  estimators in the random-effects model. *Journal of Educational and
  Behavioral Statistics*, 30, 261–293. doi:10.3102/10769986030003261
- Wang, Y., & Lê Cao, K.-A. (2020). Managing batch effects in microbiome
  data. *Briefings in Bioinformatics*, 21, 1954–1970.
  doi:10.1093/bib/bbz105
- Gibbons, S. M., Duvallet, C., & Alm, E. J. (2018). Correcting for batch
  effects in case-control microbiome studies. *PLoS Computational
  Biology*, 14, e1006102. doi:10.1371/journal.pcbi.1006102


## Decisions made while implementing (item 53, 2026-11-24)

Asked of the user at the start:

- **Names:** the `batch_` family and the planned `batch_correct()`
  wrapper (item 81) with aliases `mmuphin` / `zicombat`, `combat`,
  `combatseq` are confirmed. The wrapper, `list_batch_corrections()` and
  `match_batch_correction()` are left to item 81; item 53 adds only an
  internal `BATCH_CORRECTIONS` table (`name`, `id`, `means`) in
  `R/batch_correct.r`, which `plot()` and `print()` read for the display
  name and the "means" axis label. Item 81 adds its rows, an `alt_ids`
  column and the wrapper. Neither function takes `method`, `m`, `c`, `co`
  or other prefixes of the wrapper's arguments.
- **MMUPHin's quirks** (standard deviation as \eqn{\hat\delta}, pooled
  variance over n − 1, pseudocount before TSS): replicated, and listed
  under "Differences from ComBat" in `?batch_mmuphin`.
- **`batch_combat()`'s iteration:** `tol = 1e-4` and `max_iters = 1000L`
  are arguments (the user's choice; sva fixes conv = 1e-4 and has no cap).
  More iterations, or NaN in the iteration, is an error.
- **Non-parametric priors:** the likelihood weights are computed on the
  log scale, relative to the largest, so large batches give finite values
  where sva's products underflow to 0 and give NaN; listed under
  "Differences from sva".

Settled while implementing:

- **Files:** `R/batch_correct.r` (ecodive's own: the registry, input
  checks `batch_input()`, `batch_design()` / `batch_covariates()`, the
  result `batch_result()`, `print()` and `plot()`), `R/batch_mmuphin.r` +
  `src/batch_mmuphin.c` (MMUPHin's MIT notice), `R/batch_combat.r` +
  `src/batch_combat.c` (sva's Artistic-2.0 notice and "how it differs").
  Item 81's methods can call `batch_input()`, `batch_design()`
  (`min_size` sets the smallest batch) and `batch_result()` (any
  per-feature × batch `gamma_hat` ... `delta_star`, means before and
  after) unchanged.
- **Reference versions:** MMUPHin 2.0.0 can't be installed into R 4.3.3
  (it imports maaslin3, which needs R ≥ 4.4). Its `R/adjust_batch.R`,
  `R/helpers_adjust_batch.R` and `R/helpers.R` differ from MMUPHin
  1.16.0's (Bioconductor 3.18) only in an author's e-mail address, so the
  fixture comes from 1.16.0, and `data-raw/expected/batch_mmuphin.r` also
  sources 2.0.0's own R files (`MMUPHIN_SRC`) and checks that they give
  identical results (they do). sva 3.60.0 installed from its source
  tarball into R 4.3.3 (dependencies from Bioconductor 3.18).
- **Least squares:** both ports use `src/feature_lm.c`'s pivoted QR
  (`flm_decompose()` / `flm_fit()`, R's `qr()` rank rule with tol 1e-7)
  instead of the references' normal equations: one QR of the whole design
  (ComBat; MMUPHin without zero inflation), or one per feature of its
  non-zero samples and present batches (MMUPHin with zero inflation, which
  also gives `construct_ind()`'s rank check). Agreement is ~1e-13 or
  better; very ill-conditioned covariates may agree less (documented).
- **Sums:** R's own algorithms where they matter for parity (`mean()` and
  `var()` with long double sums, `rowSums()` / `rowMeans()` in long
  double, BLAS's double sums in column order for the matrix products,
  `pow(2, x)` for `2^x`, `nearbyint()` for `round()`); every sum in
  feature or sample order, so results are identical for any `cpus`.
- **Threads:** features in blocks (values, fits, standardisation, batch
  moments, adjustment, back-transform), samples for the totals, one task
  per batch for the shrinkage iteration (its stopping rule is a maximum
  over features), and (batch, feature) tasks for the non-parametric
  integral.
- **Batches:** a factor's levels in order (unused ones dropped), else the
  sorted unique values (`sort(method = 'radix')`, locale-independent;
  MMUPHin's `as.factor()` sorts in the locale). `batch_mmuphin()` needs at
  least two samples per batch (the plan's rule); `batch_combat()` accepts
  a one-sample batch and, as sva, turns on `mean_only`, with a **warning**
  (sva: a message).
- **Covariates:** `formula` (one-sided; a character string or a two-sided
  formula's right-hand side also work) over `data`; characters and
  logicals become factors with sorted levels; the model matrix always has
  an intercept, which is dropped (so the batch indicators carry it). The
  column space, hence the result, doesn't depend on the contrasts. A
  categorical covariate with one value, covariates confounded with each
  other, and covariates confounded with batch are errors ("Covariates and
  batch are confounded: ..."); sva's dropping of an all-ones covariate
  column is not reproduced.
- **Inputs:** `batch_mmuphin()` takes any non-negative table; it rounds
  the result only for MMUPHin's "counts" (whole numbers, not all ≤ 1).
  `pseudocount = 0` is allowed with `zero_inflation = TRUE` (zeros never
  enter); otherwise it needs a table without zeros. `batch_combat()`
  takes any finite values with `norm = 'none'`, and non-negative ones with
  `norm = 'clr'` or `'log'` (natural logs, `normalize()`'s pseudocount
  rule and warning). `NA` is an error in both.
- **Errors for MMUPHin's stops** (all features single-batch, batches that
  can't be adjusted, the iteration cap, a shrunken scale of NA) carry
  ecodive's wording with the batch or feature named; MMUPHin's own "This
  shouldn't happen" checks are unreachable and not ported.
- **Silent where the references print:** MMUPHin's `verbose` messages and
  sva's progress messages and "Found %d genes with uniform expression"
  line are not printed; features ComBat leaves unchanged are marked
  `FALSE` in the `adjusted` attribute.
- **Result attributes:** `method`, `batch_effects` (data frame `feature`,
  `batch`, `gamma_hat`, `gamma_star`, `delta_hat`, `delta_star`; NA where
  not adjusted; batch-major rows, as `as.vector()` of MMUPHin's feature ×
  batch matrices), `adjusted`, then `pseudocount` (only when chosen from
  the data), `mean_only` (`batch_combat()`, as it can be switched on by
  the data), `n_iters` (per batch, when the parametric iteration ran),
  `batch_means`, `batch_means_adjusted` (feature × batch) and
  `overall_means`. `batch_mmuphin()`'s means are of relative abundances
  (TSS, an all-zero sample counting 0), `batch_combat()`'s of the values
  corrected (after `norm`). No `display`.
- **Printing:** `print.ecodive_batch()` previews the matrix like
  `print.ecodive_source()`, with the method's display name and how many
  features were adjusted; rows are called samples or features by matching
  the dimnames with the `adjusted` names.
- **Plotting:** `plot.ecodive_batch(type = c('shrinkage', 'means'))` in
  base graphics; the "means" plot draws each feature's batch means as a
  vertical line of points at its overall mean, before (black) and after
  (`PLOT_HIGHLIGHT`, offset 1% of the largest overall mean to the right,
  or 2% on a log x axis, as MMUPHin's figure offsets it); the 1:1 line in
  both. `mean_only` results say so in the subtitle.
- **Fixtures:** `tests/tinytests/expected/batch_mmuphin.rds` (159 KB, xz):
  MMUPHin's CRC data cut to 75 samples (15 per study) and the 40 most
  prevalent species, and a simulated count table (48 × 30, three batches,
  planted shifts, features absent from a batch, present in one batch only,
  and nearly absent), eight cases (zero inflation on/off, with and without
  covariates, a given pseudocount), with MMUPHin's internal estimates.
  `tests/tinytests/expected/batch_combat.rds` (77 KB): eleven cases of
  `ComBat()` (parametric and not, with and without `mod`, `mean.only`,
  `ref.batch`, a one-sample batch, features constant within a batch, and
  clr values).
