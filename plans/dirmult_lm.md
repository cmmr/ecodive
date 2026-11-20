# `dirmult_lm()`: linear models on Dirichlet draws, with scale models (ALDEx3)

**Priority:** 50. **Tier:** port (ALDEx3 1.4.0, CRAN, **MIT**:
`DESCRIPTION` "MIT + file LICENSE", `LICENSE` "YEAR: 2026 / COPYRIGHT
HOLDER: ALDEx3 authors"; verified 2026-11-20; already the source of
`dirmult_ttest()`'s scale model and credited in `src/dirmult.c` and
`inst/COPYRIGHTS`). Read [README.md](README.md) first,
[compositional.md](compositional.md), and [dirmult_ttest.md](dirmult_ttest.md).


## What it is

`dirmult_ttest()` compares two groups. ALDEx3 (McGovern & Silverman 2026)
fits **any linear model** to every Monte Carlo draw of the Dirichlet
posterior, after adding a **scale model** to each draw:
\eqn{\log_2 W = \log_2(\text{composition}) + \log_2(\text{scale})}. The
scale model states what is known about the samples' total loads and how
uncertain it is: the clr assumption with uncertainty (`clr.sm`), no
systematic difference (`tss.sm`), a measured or predicted load per sample
(`sample.sm`: flow cytometry, qPCR, spike-ins, or the MLP load predictor of
Nishijima et al. 2025), or a prior on a design coefficient
(`coefficient.sm`). Nixon et al. (2025) and Dos Santos et al. (2026) show
scale models restore FDR control; Nishijima et al. found microbial load
"the major determinant of gut microbiome variation". compositional.md's
2026-11-20 revision replaced the planned scikit-bio `dirmult_lme()` port
with this one.

Uptake is small (ALDEx3 / SR-MEM: 4 OpenAlex citations, 5–7 Europe PMC
mentions; Nixon 2025: 27 OpenAlex; retrieved 2026-11-19), but it extends
code ecodive already has and turns `load` (from `qmp_scale()`) into a
prior. Mixed models (SR-MEM) are a later step.


## Names

`dirmult_lm()`: the regression sibling of `dirmult_ttest()` (rule 2:
family prefix; compositional.md names it so). Not `aldex()` (rule 6, as
for `dirmult_ttest()`). `DA_TESTS` row: design `formula`, effect
`estimate`, base 2, `zero_free_ok = FALSE`; alias `aldex3` (**ask**, as
for `dirmult_ttest`'s `aldex2`).


## Design

### Signature

```r
dirmult_lm(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
           pseudocount = 0.5, scale_model = c('clr', 'tss'),
           scale_sd = 0.5, load = NULL, load_sd = NULL, n_draws = 128L,
           se_type = c('hc3', 'hc0', 'classical'), alpha = 0.05,
           p_adjust = 'BH', seed = 0L, cpus = n_cpus())
```

- Model convention of [da_test.md](da_test.md); design from `da_design()`.
  ALDEx3 takes any model matrix; ecodive's treatment coding and level order
  apply.
- `pseudocount`: the Dirichlet prior, 0.5, as ALDEx3 (`Y + 0.5`) and
  `dirmult_ttest()`.
- **`scale_model`**, not `scale` (glossary: `scale` standardises
  features). `'clr'` = `clr.sm()`, `'tss'` = `tss.sm()`.
- **`scale_sd`**: their `gamma`, log2 units. Default **0.5**, ALDEx3's
  `clr.sm(gamma = 0.5)`, where `dirmult_ttest()` defaults to 0 (scikit-bio
  parity). **Ask the user** which default to use for both.
- **`load`**, **`load_sd`**: `sample.sm()`. `load` (glossary, as
  `qmp_scale()`): per-sample total load, > 0, or a column of `data`;
  `load_sd`: the log2 SD of its error, one value or one per sample
  (required with `load`). With `load`, `scale_model`/`scale_sd` are
  ignored (an error if both given explicitly). Docs show predicted load
  from the MLP (Nishijima 2025; not ported: its licence is ambiguous, see
  compositional.md) passed this way with `load_sd` from its prediction
  error.
- `coefficient.sm()` (priors on design coefficients, `c.mu`/`c.sd`): a
  later option, if asked for.
- `n_draws`: 128, as `dirmult_ttest()`; ALDEx3 defaults to 2,000 (docs:
  use ≥ 1,000 with a scale model). **Ask the user.**
- `se_type`: ALDEx3's `test` (`"t.HC3"` default, `"t.HC0"`, `"t"`).
- `p_adjust`: ALDEx3's `p.adjust.method = "BH"`.

### Algorithm (read from ALDEx3 1.4.0 `R/aldex.R`, `aldex-internals.R`, `fflm.R`, `scale.R`, `dirichlet.R`)

For each draw \eqn{d}:

1. **Composition**: \eqn{G_{sf} \sim \text{Gamma}(x_{sf} + 0.5, 1)},
   \eqn{\log_2 \text{comp}_{sf} = (\log G_{sf} - \log\sum_g G_{sg})/\log 2}
   (`rLogDirichlet()`).
2. **Scale** \eqn{\log_2 S_{s}}:
   - clr: \eqn{-\text{mean}_f \log_2\text{comp}_{sf} + x_s^\top\Lambda},
     \eqn{\Lambda \sim N(0, \gamma^2 I_P)}, one per design column
     including the intercept;
   - tss: \eqn{x_s^\top\Lambda};
   - load: \eqn{N(\log_2 \text{load}_s, \text{load\_sd}_s^2)},
     independent per sample.
3. \eqn{Y = \log_2 W = \log_2\text{comp} + \log_2 S} (n × D).
4. **OLS** (`fflm()`): \eqn{\Theta = (X^\top X)^{-1}X^\top Y};
   \eqn{e = Y - X\Theta}; with \eqn{A = (X^\top X)^{-1}X^\top},
   \eqn{h = \text{diag}(XA)}:
   - classical: \eqn{se_p = \sqrt{[(X^\top X)^{-1}]_{pp}\,\sum_s e_s^2/(N-P)}};
   - HC0: \eqn{se_p = \sqrt{\sum_s A_{ps}^2 e_s^2}};
   - HC3: \eqn{se_p = \sqrt{\sum_s A_{ps}^2 e_s^2/(1-h_s)^2}}.

   \eqn{t = \Theta/se}; `p_upper = pt(t, N - P, lower.tail = FALSE)`;
   `p_lower = 1 - p_upper` (computed exactly so, for parity).
5. **Aggregate** (`aldex.pvals()`, two-sided): per draw
   \eqn{\ell = \min(1, 2p_{lower})}, \eqn{u = \min(1, 2p_{upper})};
   `p_value` = \eqn{\min(\overline{\ell}, \overline{u})} (means over
   draws); `p_adjusted` = \eqn{\min(\overline{\text{adj}(\ell)},
   \overline{\text{adj}(u)})}, where `adj` is `p.adjust(p_adjust)` across
   features, per covariate and draw. **Note:** ALDEx3 averages *adjusted*
   p-values; `dirmult_ttest()` adjusts *averaged* ones. Follow each
   reference; document the difference.
6. `estimate`, `se`: means over draws (`summary.aldex()`).

### Reuse of `dirmult_ttest()`'s C

- Factor phase A of `src/dirmult.c` into a shared draw routine: the same
  PCG stream per (seed, draw, table sample), `d * (n + 1) + s`, so for the
  same `seed` both functions see the **same Gamma draws** for a sample.
  `dirmult_ttest()` must stay bit-identical (its tests guard this).
- Scale normals from stream `d * (n + 1) + n` (as `dirmult_ttest()`'s one
  shift), drawn in design-column order (clr/tss) or sample order (load).
- Phase B: the OLS of [da_test.md](da_test.md)'s `src/feature_lm.c`
  (precompute \eqn{A}, \eqn{h}, \eqn{(X^\top X)^{-1}} once; blocks of
  features per thread). `pt()` on the main thread; BH per covariate per
  draw needs every feature's p: sort per covariate in threads (one task
  per covariate), then running sums in draw order on the main thread.
- Memory: one n × D draw, plus P × D running sums (×4: estimate, se,
  \eqn{\ell}, u, and their adjusted sums).

### Result

`c('ecodive_dirmult_lm', 'data.frame')`, one row per feature × design
column (intercept included, as `ancombc()`; ALDEx3's summary hides it by
default): `feature, covariate, estimate` (log2), `se, p_value,
p_adjusted, reject`. Attributes: `formula`, `groups`, `design`,
`scale_model`, `scale_sd`, `load_sd`, `n_draws`, `pseudocount`,
`se_type`, `p_adjust`, `alpha`, `seed`. `plot()`: volcano / forest
(estimate ± qt·se, labelled as approximate) / ma via
[da_plots.md](da_plots.md). No `overlap` statistic: ALDEx3's
`aitchison.mean()` is commented "Private function from ALDEx2" (GPL
provenance); don't port it.

### Later: mixed models (SR-MEM)

ALDEx3's `method = "lme4" / "nlme" / "blmm"` (McGovern & Silverman 2026).
`blmm` uses TMB/C++ (`src/blmm.cpp`, part of the MIT package, but TMB and
RcppEigen are dependencies ecodive can't take). Plan a C random-intercept
REML (shared with [linda.md](linda.md) and [maaslin3.md](maaslin3.md));
scikit-bio's BSD `dirmult_lme()` (`_lme.py`) is a permissive reference for
the closed-form random-intercept fit. Separate step.


## Reference implementations

- **Port from:** ALDEx3 1.4.0 as above. Extend the existing ALDEx3 stanza
  in `inst/COPYRIGHTS` and the header of the new C file.
- **Fixtures:** `data-raw/expected/dirmult_lm.r`: ALDEx3 accepts `scale`
  as an N × nsample matrix, and its Gamma draws come from R's `rgamma()` in
  a known order (`aldex.sampler()`: samples, then features, `nsample` at
  a time). So: `set.seed(1)`, run `aldex(Y, X, scale = <logScale matrix>,
  test = ...)`; `set.seed(1)` again and regenerate the same `rgamma()`
  calls to export the draws; feed both to an internal `dirmult_lm_details(
  gamma_draws =, scale_draws =)` (as `dirmult_ttest_details()`). Cover
  `se_type` × {two groups, three groups + numeric}, and `clr.sm`,
  `tss.sm`, `sample.sm` draws captured the same way.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- **`da` record** "Dirichlet-Multinomial Linear Model" via `call_step()`:
  `formula`, `scale_sd`, `load_sd`, the scale model (`"clr"`, `"tss"`,
  `"load"`), `se_type`, `n_draws`, `pseudocount`, `p_adjust`, `seed`.
- Attributes for plots: `alpha`, `log_base` (`2`, as ALDEx3), `conf_level`
  if intervals are returned; the input's `transform` record; and, when a
  measured or predicted `load` was used, its source in the record
  (`load_from`: `"measured"` or `"predicted"`), so titles can say so.
- `plot()`: [da_plots.md](da_plots.md)'s types per covariate.

## Tests

- With injected draws: `estimate`, `se`, `p_value`, `p_adjusted` match
  ALDEx3 at 1e-10 (it solves normal equations; ecodive uses QR).
- With `scale_sd = 0`, `scale_model = 'clr'`, `~ group` (two levels) and
  the same `seed`: `estimate` for the group column equals
  `dirmult_ttest()`'s `log2_fc` (contrast reversed as needed) at 1e-12;
  same draws, same mean difference.
- `identical()` for `cpus` 1, 2, 4; dense = sparse = `margin = 2`.
- Errors: N ≤ P, `load` ≤ 0 or `NA`, `load` without `load_sd`, `load`
  and `scale_sd` both given, invalid `se_type`; rank-deficient design.
- Smoke line.


## Docs

Rd page in *Differential abundance*: the scale model (with a table of
which to use when, following ALDEx3's own guidance), units (log2), the
p-value aggregation, the relation to `dirmult_ttest()`, a γ-sensitivity
pointer to [da_sensitivity.md](da_sensitivity.md), and "Same names in
other packages": ALDEx3's `aldex()` (the original). Glossary:
`scale_model`, `load_sd`, `se_type` (new). NEWS.


## References

- McGovern, K. C., & Silverman, J. D. (2026). Scale reliant mixed effects
  models enhance microbiome data analysis. *Microbiome*, 14, 117.
  doi:10.1186/s40168-026-02377-x (4 OpenAlex, retrieved 2026-11-19).
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). *Genome
  Biology*, 26, 139. doi:10.1186/s13059-025-03609-3.
- Dos Santos, S. J., et al. (2026). *PLoS Computational Biology*, 22(9),
  e1014728. doi:10.1371/journal.pcbi.1014728.
- Nishijima, S., et al. (2025). Fecal microbial load is a major determinant
  of gut microbiome variation and a confounder for disease associations.
  *Cell*, 188(1), 222–236. doi:10.1016/j.cell.2024.10.022.
- Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
  standard errors in the linear regression model. *The American
  Statistician*, 54, 217–224.
