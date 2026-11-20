# `fastancom()`: ANCOM from K marginal models

**Priority:** 54. **Tier:** **clean room**. fastANCOM 0.0.4
(github.com/ZRChao/fastANCOM, not on CRAN) is **GPL-3** (`DESCRIPTION`
verified 2026-11-20; the repository has no LICENSE file). Work only from
the paper, the package documentation (`man/fastANCOM.Rd`, README, the
author's RPubs manual; all read) and its **outputs**. Do not open its
`R/` directory. The paper (Bioinformatics application note) is not open
access and its publisher page was unreachable, so several steps below are
**unverified** (✱) and must be settled from outputs first. Read
[README.md](README.md) first, [compositional.md](compositional.md) and
[ancom.md](ancom.md).


## What it is

ANCOM tests all \eqn{K(K-1)/2} log ratios. fastANCOM "fits only K models
for log-transformed counts", "provides P-values ... and outputs log fold
changes" (abstract). Wirbel et al. (2024) found it the **only
microbiome-specific method** among those that "properly control false
discoveries at relatively high sensitivity" (with linear models, Wilcoxon,
t-test and limma), and it supports covariate adjustment (`Z`). Pelto et al.
(2025): "good consistency but low sensitivity". Uptake is small (22
Crossref, 17 Europe PMC citations; 12 Europe PMC mentions 2023–26;
retrieved 2026-11-19); it ranks here for its benchmark standing and low
cost.


## Names

`fastancom()`: the published name lowercased (rule 1), as ANCOM-BC is
`ancombc()` rather than `ancom_bc()`; it therefore sorts next to `ancom`.
`DA_TESTS` row: design `formula` (covariates allowed), effect `log2_fc`,
base 2 (✱), `zero_free_ok = TRUE`; alias `fastancom`. "Same names in other
packages": fastANCOM's `fastANCOM()` (the original).


## Design

### Signature

```r
fastancom(counts, formula = NULL, data = NULL, groups = NULL, margin = 1L,
          pseudocount = 0.5, alpha = 0.05, detect_rate = 0.7,
          ref_rate = 0.05, p_adjust = 'BH', cpus = n_cpus())
```

- fastANCOM's `x` (variable of interest) is the **term of interest**
  `groups` of [da_test.md](da_test.md)'s convention; its `Z` (covariates)
  are the other terms of `formula`. `groups` is required when `formula`
  has more than one term. A factor `x` with > 2 levels: fastANCOM has a
  `global` test and per-pair `results` (its docs; behaviour ✱). First
  version: two-level factors and numeric `x` only, an error otherwise;
  more levels in a later step once probed.
- `pseudocount = 0.5`: the Rd's `\usage` says `pseudo = 0.5` (its text
  says "default is 0.05"; the usage line is what runs). Added to every
  count ✱ (confirm: every value or zeros only).
- `alpha` = its `sig`; `detect_rate`, `ref_rate` = its `detect.rate`,
  `ref.rate` (snake case).
- Not carried over: `zero_cut = 0.9` (drops taxa with > 90% zeros),
  `lib_cut`, `struc_zero`, `outlier` (Cook's distance), `rand` (random
  effects): left to the user or later steps. Fixtures set `zero_cut = 1`.

### Algorithm (docs + reasoning; every ✱ to be confirmed from outputs)

1. \eqn{L = \log(Y + 0.5)}, base ✱ (the output is `log2FC`, so likely
   \eqn{\log_2}).
2. **Marginal models**: per taxon \eqn{k}, OLS of \eqn{L_k} on
   \eqn{[1, x, Z]} (the docs' "log-linear regression", "fits the marginal
   model"): \eqn{\hat\beta_k} (the coefficient of \eqn{x}), residuals
   \eqn{e_k}.
3. **Pairwise tests without pairwise fits**: with one design for all
   taxa, the log-ratio model's coefficient is \eqn{\hat\beta_k -
   \hat\beta_j} and its variance \eqn{c_{xx}\,\mathrm{var}(e_k - e_j)}
   (\eqn{c_{xx}} = the \eqn{x} entry of \eqn{(X^\top X)^{-1}}), so every
   pairwise Wald statistic follows from the \eqn{K} fits and the residual
   covariance \eqn{S = E^\top E/(n - P)}:
   \eqn{T_{kj} = (\hat\beta_k - \hat\beta_j)/\sqrt{c_{xx}(S_{kk} + S_{jj}
   - 2S_{kj})}}. This is exactly OLS on each log ratio, i.e. ANCOM with
   a t / F test. ✱ Reference distribution (normal, as "Wald" suggests, or
   t with \eqn{n - P} df) and whether p-values are adjusted per taxon
   before counting (ANCOM's Holm/BH row adjustment) or not.
4. **W**: \eqn{W_k = \#\{j: p_{kj} < \alpha\}} ("Reject.number"; the
   manual's example has 97–99 of 99). `reject` (its `REJECT`):
   \eqn{W_k > \text{detect\_rate}\,(K - 1)} ✱ (`>` or `≥`).
5. **Effect sizes** (`effect = TRUE`): the reference set \eqn{R} = the
   \eqn{\lceil \text{ref\_rate}\cdot K\rceil} ✱ taxa with the smallest W
   (ties ✱); `log2_fc`\eqn{_k = \hat\beta_k - \text{mean}_{j\in R}
   \hat\beta_j}; its SE \eqn{\sqrt{c_{xx}\,\mathrm{var}(e_k - \bar e_R)}}
   from \eqn{S}; Wald p ✱; q = `p.adjust(p_adjust)` ✱ (its output has
   `log2FC.qval`; BH likely).

**Black-box protocol (do first).** Run fastANCOM on crafted inputs and
compare its `joint` output (`beta`, `se`, `statistic`, `W`) and `final`
table with candidate formulas: (a) one taxon pair with known log-ratio
statistics decides normal vs t and adjusted vs raw; (b) W exactly at
\eqn{0.7(K-1)} decides `>`/`≥`; (c) `K` with \eqn{0.05K} non-integer and
tied W decides the reference-set rule; (d) a zero-free table with
`pseudo = 0` (if accepted) settles the log base and whether the pseudo
count is added everywhere. Record the answers in this plan.

### C and memory

- Marginal fits through `src/feature_lm.c` ([da_test.md](da_test.md)).
- \eqn{S} is K × K (5,000 taxa: 200 MB): don't store it. Process row
  blocks: each task forms \eqn{S_{b\cdot} = E_b^\top E/(n-P)} for its block
  of taxa (O(n K b)), counts its rows' W and keeps \eqn{S_{kR}}, then
  frees the block. Rows deal to threads as in `ancom.c`; results
  identical for any `cpus`. O(n K²) time, like one ANCOM pass of means,
  but no per-pair test over samples.
- The reference set needs all W first: two passes (W, then effects).

### Result

`c('ecodive_fastancom', 'data.frame')`, one row per taxon: `feature,
log2_fc, se, w` (integer), `p_value, p_adjusted, reject`. `reject` is the
W rule (fastANCOM's `REJECT`, as `ancom()`), not `p_adjusted <= alpha`;
the docs say how to filter on `p_adjusted` instead. Attributes: `formula`, `groups`, `alpha`,
`detect_rate`, `ref_rate`, `reference` (the reference taxa),
`pseudocount`, `p_adjust`. `plot()`: `plot_ranked()` of W (as `ancom()`)
by default, `type = 'volcano'` / `'forest'` via [da_plots.md](da_plots.md).


## Reference implementations

- **Output-only:** fastANCOM 0.0.4 (GitHub, GPL-3). Install from GitHub in
  the fixture generator only.
- **Fixtures:** `data-raw/expected/fastancom.r`: the black-box probes
  above; the README / manual example (`set.seed(123)`, 60 × 100 Poisson,
  first 10 taxa scaled by \eqn{2^b}); two groups with a numeric
  covariate `Z`; a numeric `x`. Store `joint` and `final` as text.
- **Cross-check inside ecodive:** for two groups, no covariates and ✱ t
  reference with no row adjustment, `w` should equal `ancom(sig_test =
  'anova', p_adjust = 'none')`'s W (F = t² with the same df): a strong
  consistency test once step 3 is pinned down.


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- **`da` record** "fastANCOM" via `call_step()`: `pseudocount`,
  `sig_level`/`alpha`, `p_adjust`, reference-set rule, `formula` or
  `groups`.
- Attributes for plots: `alpha`, `log_base` (`exp(1)`, check), the `W`
  cutoff (for the ANCOM-style dot chart), the input's `transform` record.
- `plot()`: the ANCOM dot chart and [da_plots.md](da_plots.md)'s volcano.

## Tests

- After the protocol: `w`, `reject` identical to fastANCOM; `log2_fc`,
  `se`, `p_value` at 1e-10.
- Invariant to multiplying a sample by a constant on a zero-free table
  without pseudocount (log ratios and reference differences cancel it).
- `identical()` for `cpus` 1, 2, 4; dense = sparse.
- Edge cases: K = 2; all W tied; `ref_rate` giving one taxon; n ≤ P;
  constant taxon (zero variance in S); NA in `data`; > 2-level `groups`
  (error, first version).
- Smoke line.


## Docs

Rd page in *Differential abundance*: what W and the reference set mean,
the equivalence with ANCOM's log-ratio tests, the benchmark evidence, and
"Differences from fastANCOM" (no filters, structural zeros, outlier
removal or random effects; design coding as `ancombc()`; any ✱ choice
that couldn't be matched exactly). Glossary: `detect_rate`, `ref_rate`
(new). NEWS.


## References

- Zhou, C., Wang, H., Zhao, H., & Wang, T. (2022). fastANCOM: a fast
  method for analysis of compositions of microbiomes. *Bioinformatics*,
  38(7), 2039–2041. doi:10.1093/bioinformatics/btac060 (22 Crossref,
  retrieved 2026-11-19).
- Mandal, S., et al. (2015). Analysis of composition of microbiomes: a
  novel method for studying microbial composition. *Microbial Ecology in
  Health and Disease*, 26, 27663.
- Wirbel, J., et al. (2024). *Genome Biology*, 25, 247.
  doi:10.1186/s13059-024-03390-9.
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130.
