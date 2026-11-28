# `fastancom()`: ANCOM from K marginal models

**Priority:** 80. **Tier:** **clean room**. fastANCOM 0.0.4
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Title: the plot type, " of " and the `DA_TESTS` name for `method`
  (`"fastancom"`), for example "Volcano plot of fastANCOM". The input is a
  count table, so no `display`.
- Attributes plots read: `alpha`, `log_base` (the base step 1 settles ✱),
  and `detect_rate`, from which the dot chart draws the `W` cutoff
  \eqn{\text{detect\_rate}\,(K - 1)}.
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


## Decisions made while implementing (item 80, 2026-11-27)

Asked of the user at the start:

- **`pseudocount`:** `NULL`, the shared rule (nothing added to a zero-free
  table, half the smallest non-zero value with a warning otherwise), not
  fastANCOM's 0.5 (the user's choice; the plan's signature said 0.5).
  `0` is allowed for a table without zeros, an error otherwise ("`counts`
  has zeros, and `pseudocount = 0` leaves their logs undefined.").
- **`p_adjust`:** `'holm'`, as ecodive's other DA methods (fastANCOM: BH).
- **More than two levels:** implemented now (the user's choice), as
  ecodive's own design, because fastANCOM 0.0.4's `global` and pairwise
  paths fail on every input tried (character, factor, integer 0/1/2 and
  1/2/3 `x`, named or not, `global` TRUE or FALSE, `zero_cut` 0.9 or 1:
  "contrasts can be applied only to factors with 2 or more levels", "0
  (non-NA) cases", "subscript out of bounds", "object 'res2' not found").
  Its internal function (`.fANCOM`, named in a warning) is not in the
  namespace, so it couldn't be observed without reading the source.
- **Result layout:** a `covariate` column (the user's choice), as `linda()`:
  one row per feature and comparison, feature-major.

Settled from fastANCOM's outputs (the black-box protocol; probes in the
session's scratch directory, not kept):

- **Log base and pseudocount:** \eqn{\log_2(Y + \text{pseudo})}, the
  pseudocount added to every value (zeros-only is off by 0.05 in β);
  `pseudo = 0` is accepted on a zero-free table. Proportions are used as
  given (plus 0.5). `x` is coded by `lm()`'s treatment contrasts (character
  levels sorted, a factor's kept, logical FALSE first); a numeric `x` is
  continuous. `Z` is a numeric vector or matrix (a data.frame errors).
- **Marginal model:** OLS of each feature on \eqn{[1, x, Z]}; `marginal$cov`
  is \eqn{c_{xx} E^\top E / (n - \mathrm{rank})} exactly.
- **Joint (pairwise) tests:** `joint$beta[k, j]` \eqn{= \beta_k - \beta_j},
  `se` \eqn{= \sqrt{C_{kk} + C_{jj} - 2 C_{kj}}}, `statistic` \eqn{= -|t|}.
  W counts \eqn{p < \alpha} with **no row adjustment**, p from the **t
  distribution** (not normal). The degrees of freedom are **\eqn{n - 2}
  whatever the covariates** (W fits only \eqn{n - 2} on a 300-feature table
  with 1, 2 or 3 numeric covariates, and a 3-level character `Z`); the
  final p-values too. So with covariates fastANCOM's variances use \eqn{n -
  p} and its tests \eqn{n - 2}: `fastancom()` uses \eqn{n - p} for both
  (documented). Without covariates the two agree.
- **Reference set:** the first `length(1:(ref.rate * K))` features of
  `order(W)` (stable: ties in feature order). R's `:` adds a fuzz, so 0.29
  × 100 = 28.999999999999996 gives 29, 19.95 gives 19, and anything below
  1 gives 1 (with the warning "Only take ... components as reference!"
  below 5, which `fastancom()` doesn't give). `log2FC` \eqn{= \beta_k -
  \mathrm{mean}_R \beta}.
- **Standard errors:** \eqn{\sqrt{C_{kk} + \sum_{R,R} C / m^2 - 2
  \sum_R C_{k\cdot} / m}} with \eqn{m} = `ref.rate * K` (at least 1), not
  the set's size: so they are the standard errors of `log2FC` only when
  `ref.rate * K` is a whole number (or below 1). Found from a feature alone
  in the reference set, whose SE grows as \eqn{(1 - 1/m)} with \eqn{m}.
  **Not copied:** `fastancom()` divides by the size for both, so that `se`
  is `log2_fc`'s standard error (documented); every fixture has a whole
  `ref.rate * K` or one below 1. p-values: \eqn{2 F_t(-|t|, n - 2)}; q-values:
  `p.adjust(, 'BH')` over the non-NaN p-values.
- **Decisions** (`REJECT`): with `effect = FALSE`, fastANCOM returns
  `detected` and `detect0.7` instead, which showed the rule:
  - `detect0.7` is \eqn{W \ge 0.7 K}, with 0.7 fixed whatever
    `detect.rate` (4,000 of 4,000 random tables);
  - `detected` is the **original ANCOM's empirical cutoff** (Mandal et al.
    2015, as its R code): \eqn{c_0 = \max W / K}, cutoffs \eqn{c_0 - (0.05,
    \ldots, 0.25)}, proportions with \eqn{W \ge K c_i}, τ = 0.02, \eqn{\nu}
    chosen among the **first four** cutoffs (scikit-bio's `ancom()`, which
    ecodive's `ancom()` follows, picks among the last four), detected = W
    at least the smallest \eqn{W \ge \nu K}; nothing for \eqn{K < 10} or
    \eqn{c_0 < 0.26} (boundary between 46/177 = 0.2599, none, and 13/50 =
    0.26, detected; 4,000 of 4,000);
  - `REJECT` \eqn{= (W \ge \text{detect.rate} \cdot K) \lor (\text{detected}
    \land W \ge K/2)} (5,500 of 5,500 random tables, `detect.rate` from 0.05
    to 0.95). Not \eqn{W > \text{detect.rate}\,(K - 1)} as the plan guessed.
  The rule is `fastancom_reject()`; the `cutoff` attribute is the smallest
  W it rejects, \eqn{\min(aK, \max(\text{up}, K/2))}.

Decided along the way:

- **Signature:** `fastancom(counts, formula = NULL, data = NULL, groups =
  NULL, margin = 1L, pseudocount = NULL, alpha = 0.05, detect_rate = 0.7,
  ref_rate = 0.05, p_adjust = 'holm', cpus = n_cpus())`. `ref_rate` in (0,
  0.5] (fastANCOM's docs: "between 0 and 0.5"); `detect_rate` and `alpha`
  strictly in (0, 1). `groups` is required when `formula` has more than one
  term. The model must have an intercept (`- 1` is an error), and a
  rank-deficient design is an error (`da_design()`'s warning is turned into
  one), as are \eqn{n \le p}, NA in `data` and single-level factors
  (`da_design()`'s checks).
- **Columns:** `feature, covariate, log2_fc, se, statistic, w, p_value,
  p_adjusted, reject`; `statistic` added (t, or F for a global test).
  Attributes: `method, formula, groups, design, alpha, detect_rate,
  ref_rate, p_adjust, pseudocount, log_base` (2), `dof, cutoff` (per
  comparison), `reference` (per comparison, feature names), `coef` (D × P,
  fastANCOM's marginal β in column `x`), `mean` (mean log2 abundance, for
  the MA plot).
- **Comparisons for a factor of L > 2 levels:** first the global test (row
  `covariate` = the term's name), then each level against the reference
  (`'statusb'`), then each other pair (`'statusc - statusb'`). Global: the
  F-test of the term's \eqn{q} coefficients on each log ratio, \eqn{F =
  d^\top A d / (q s^2)} with \eqn{A = (\Lambda^\top (X^\top X)^{-1}
  \Lambda)^{-1}}; without covariates its W is `ancom(sig_test = 'anova',
  p_adjust = 'none')`'s, as the two-group W is. Its effect row has
  `log2_fc` and `se` NaN, and the F-test of the feature against the
  reference set's mean coefficients as `statistic` / `p_value`. Pairs are
  contrasts in the full model (variance from all samples), not refits on
  the two groups' samples (which is probably what fastANCOM intended,
  judging from its "0 (non-NA) cases" error). Each comparison has its own W,
  reference set and decisions. A non-factor term of several columns
  (`poly(t, 2)`) gets the global test only.
- **C (`src/fastancom.c`):** three phases, each dealt to threads in turn:
  (1) per feature, `log2(x + pc)`, the fit with `feature_lm.c`'s QR (one
  decomposition, shared read-only) and its residuals (new `flm_resid()`,
  factored out of `flm_se()` with the same arithmetic), and each
  comparison's contrasts; (2) per row k of the pair triangle, \eqn{\lVert e_k
  - e_j \rVert^2} and every comparison's F, counted against the critical
  value `qf(1 - alpha, q, n - p)` (equivalent to \eqn{p < \alpha}) into
  per-thread integer counts summed afterwards; (3) per comparison, the
  reference set by a stable counting sort of W and its mean residuals (on
  the main thread, fixed order), then per feature the effect, SE and
  statistic. No \eqn{K \times K} matrix: memory is the n × K residuals.
  Results identical for any `cpus` (checked 1-8), dense, `Matrix`, `slam`
  and `margin = 2`.
- **Exact zero variance:** a feature whose log values are all equal gets
  coefficients (its value, 0, ...) and residuals of exactly 0, so two such
  features' log ratio has no test (0/0, not counted), as `ancom()`'s
  constant log ratios. A log ratio with zero residual and a non-zero
  difference has F = ∞ (counted).
- **No warning** for a small reference set (fastANCOM warns below 5).
- **Plots:** `plot.ecodive_fastancom(x, covariate, type = c('w', 'volcano',
  'ma', 'forest', 'rank'), ...)`: `'w'` is `plot_ranked()` (ancom's dot
  chart) titled "W of fastANCOM" with the `cutoff` line; the others are
  `da_plot()`'s ("Volcano plot of fastANCOM"). The default comparison is
  the global test when there is one. A global test has no fold changes:
  only `'w'`.
- **`da_sensitivity()`** works unchanged (the comparisons are its
  covariates); its docs list `fastancom`.
- **Fixtures** (`data-raw/expected/fastancom.r`): the manual's example, a
  numeric `x`, two groups with a numeric covariate (with
  `joint$statistic`, so the tests count W on \eqn{n - p}), zeros with other
  `sig`/`detect.rate`/`ref.rate`, `pseudo = 1`, `pseudo = 0`, proportions,
  8 features, `ref.rate = 0.29` (the `:` fuzz) and W ties; and 600 random
  tables' W, `detected` and `REJECT` for the rules.


## Outcome (from the implementing session's final report)

**Item 80, done 2026-11-28.** Commits `08c2201` (work), plus the tick commit.

- Tests: 37,518 at home (+2,035 from 35,483), 530 in CRAN mode (+4 from
  526). `R CMD check --no-manual --ignore-vignettes`: 1 NOTE (installed
  size, as before).
- fastANCOM 0.0.4, 10 cases (manual example, numeric `x`, covariate, zeros,
  `pseudo` 1 and 0, proportions, 8 features, `ref.rate = 0.29`, W ties): W
  and REJECT identical; marginal β, log2 FC, SE, statistic, p and q within
  1.1e-12 relative (p and q not compared for the covariate case, whose df
  differs by design; W there counted from its stored statistics on n − p).
  Decision rules identical on 600 random tables (5,500 in the probes).
  Without covariates W is `ancom(sig_test = 'anova', p_adjust = 'none')`'s,
  for two and three groups.
- Speed (4 cores, 3,000 features × 200 samples): 0.16 s (0.59 s on 1 core);
  fastANCOM 1.71 s (and three K × K matrices, 216 MB), `ancom()` 3.4 s.
- Open: fastANCOM's multi-group path couldn't be compared (it fails on
  every input); ecodive's global and pairwise tests are checked against an
  independent per-log-ratio `lm.fit()` computation instead.
