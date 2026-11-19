# Compositional methods: overview and shared conventions (2026-11-19)

Read [README.md](README.md) first. This file is the decision record for
the compositional-data methods, the counterpart of
[new_methods.md](new_methods.md) for the earlier batch. Each method's own
plan holds everything needed to build it; this file holds what they share.
**Where a method plan and this file differ, the method plan wins.**


## Why

Sequencing counts carry only relative information: a sample's total is
set by the sequencer, not by the community. Gloor et al. (2017) and many
since argue that analyses should therefore use **log ratios**, handle
**zeros** on purpose rather than by an arbitrary pseudocount, and avoid
correlations between relative abundances. ecodive already has
`norm = 'clr'` / `'rclr'`, `aitchison()`, `robust_aitchison()` and
`rpca()`, and its tests accept any distance. The plans below fill the
gaps around them.


## The plans

| Plan | Functions | Tier |
|:-----|:----------|:-----|
| [zero_replacement.md](zero_replacement.md) | `replace_zeros()` (wrapper), `replace_zeros_mult()`, `replace_zeros_czm()`, `replace_zeros_gbm()`, `replace_zeros_sq()`, `replace_zeros_bl()`, plus `list_zero_replacements()` / `match_zero_replacement()` | port (scikit-bio) + formula (zCompositions output-only) |
| [log_ratio_transforms.md](log_ratio_transforms.md) | `alr()`, `ilr()` | port (scikit-bio) |
| [log_ratio.md](log_ratio.md) | `log_ratio()`, `balance()` | port (Qurro, scikit-bio / gneiss) |
| [proportionality.md](proportionality.md) | `vlr()`, `rho_p()`, `phi_s()` | port (scikit-bio `vlr`) + formula (propr output-only) |
| [qmp.md](qmp.md) | `qmp_scale()`, `qmp_rarefy()`; per-sample `rarefy(depth =)` | formula over ecodive's `rarefy()` |
| [ancom.md](ancom.md) | `ancom()` | port (scikit-bio) |
| [dirmult_ttest.md](dirmult_ttest.md) | `dirmult_ttest()` | port (scikit-bio; scale model from ALDEx3, MIT) |
| [ancombc.md](ancombc.md) | `ancombc()`, `ancombc2()`; post-hoc `ancombc_global()`, `ancombc_pairwise()`, `ancombc_dunnett()`, `ancombc_trend()`; `structural_zeros()` | port (scikit-bio, BSD-3; derived from Artistic-2.0 ANCOMBC: both notices) |
| [philr.md](philr.md) | `philr()`, `philr_distance()` | port (basis, scikit-bio / gneiss) + formula (weights; philr output-only) |
| [rpca_phylo.md](rpca_phylo.md) | `rpca_phylo()` | port (gemelli) on ecodive's `rpca()` |
| [rpca_joint.md](rpca_joint.md) | `rpca_joint()` | port (gemelli) |
| [ctf.md](ctf.md) | `ctf()` | port (gemelli) |

[TODO.md](TODO.md) gives the order.

**Not planned** (reasons in the method plans that mention them):

- **SparCC:** the original's license is unknown (its Bitbucket repository
  is gone), FastSpar and SpiecEasi are GPL-3, and the MIT/BSD Python ports
  (dlegor/SparCC, fast_sparCC) say they follow the original closely, so
  their provenance is unclear. It could be written clean-room from
  Friedman & Alm (2012, CC-BY); not planned now, since `rho_p()` answers
  the same question without SparCC's sparsity assumption.
- **ALDEx2 itself, LinDA, MaAsLin2/3:** `dirmult_ttest()` covers the
  ALDEx2 approach; LinDA is GPL-3; MaAsLin2/3 are MIT and could be ported
  later, but they are general regression frameworks rather than
  compositional methods.
- **scikit-bio `dirmult_lme()`:** its default path is a closed-form
  random-intercept REML fit (`_lme.py`) that could be ported, but parity
  needs its grid and golden-section search and statsmodels-style standard
  errors reproduced exactly; revisit after `dirmult_ttest()`.
- **Compositional alpha diversity:** Shannon, Simpson and the Hill numbers
  are computed from proportions already, so they are scale invariant.
- **New compositional tests:** PERMANOVA, PERMDISP, ANOSIM, MRPP and the
  Mantel family on `aitchison()`, `robust_aitchison()`, `philr_distance()`
  or `rpca(distances = TRUE)` already are them.


## Names (README naming rules applied)

- **Families share a prefix** (rule 2), so they sort together:
  `replace_zeros_*`, `qmp_*`, `rpca_*`.
- **Standard acronyms stay** where the field uses them as function names
  (rule 1): `alr`, `ilr`, `philr`, `vlr`, `ancom`, `ancombc`, `ancombc2`,
  `ctf`. ANCOM-BC's post-hoc tests take its fitted result, so they are
  `ancombc_<test>` (rule 2); `structural_zeros()` is spelled out (rule 4).
- **Spelled out or qualified where an acronym would mislead** (rule 4):
  - `rho_p()`, not `rho()` (Spearman's ρ) — the ρp of Erb & Notredame
    (2016);
  - `phi_s()`, not `phi()` (the φ coefficient) or `phs()` (propr's code,
    which reads like a typo) — the symmetric φs of Quinn et al. (2017);
  - `philr_distance()`, not a second `philr()` — the transform keeps the
    package's name, the distance says what it returns;
  - `rpca_phylo()` and `rpca_joint()`, not `phylogenetic_rpca()` /
    `joint_rpca()` (gemelli's names), so they sort with `rpca()`;
  - `dirmult_ttest()`, scikit-bio's name, not `aldex()`: it resembles
    ALDEx2 but isn't it, and ALDEx2's name shouldn't be borrowed (rule 6).
- **Zero replacements keep zCompositions' method codes** (`gbm`, `sq`,
  `bl`, `czm`) because that's what users will have typed in `cmultRepl()`,
  with spelled-out aliases in the registry.
- **No `clr()` or `rclr()` functions.** They are `norm` choices, applied
  inside every function with a `norm` argument without making a dense
  copy, and `normalize()` returns them. A second route would be redundant.
  (`alr()`, `ilr()` and `philr()` can't be `norm` choices: they need
  extra inputs and change the number of features.)

### New argument names (add to [parameter_names.md](parameter_names.md) §1)

| Name | Meaning | Used by |
|:-----|:--------|:--------|
| `numerator`, `denominator` | features on the top and bottom of a log ratio: names, indices, or a logical vector | `log_ratio`, `balance`, `alr` (`denominator` only), `rho_p` / `phi_s` (`denominator`: `NULL` = the clr's geometric mean) |
| `basis` | an ILR basis, or a sequential binary partition | `ilr` |
| `load` | each sample's measured total microbial load | `qmp_scale`, `qmp_rarefy` |
| `contrast` | the samples to compare, `c(treatment, reference)`; `reference = NA` pools all other samples | `dirmult_ttest` |
| `part_weights`, `ilr_weights` | PhILR's taxon and balance weights (philr's `part.weights`, `ilr.weights`) | `philr` |
| `subjects`, `states` | per-sample subject ids and time points / conditions | `ctf` |
| `test_samples`, `n_test_samples` | held-out samples for cross-validation | `rpca_joint` |
| `max_power_iters` | a cap on power-method iterations (`max_<unit>` rule) | `ctf` |
| `fit` | a fitted result passed to a post-hoc test | `ancombc_global`, `ancombc_pairwise`, `ancombc_dunnett`, `ancombc_trend` |
| `n_draws` | an exact number of Monte Carlo draws (`n_<unit>` rule) | `dirmult_ttest` |
| `delta`, `frac` | the replacement value (proportion) and the fraction of the detection limit | `replace_zeros_mult` (`delta`), `replace_zeros_czm` (`frac`) |

`reference` is **not** used as an argument name: it would mean a feature
in `alr()` and a group in `dirmult_ttest()`. `denominator` and `contrast`
say which.


## Shared conventions

- **Logs are natural logs**, as in `normalize(norm = 'clr')`. Where a
  reference reports log2 (scikit-bio's `dirmult_ttest` fold changes), the
  plan says so and keeps the reference's base for that column only.
- **Zeros.** Every log-ratio function takes `pseudocount` with the shared
  rule in `validate_pseudocount()`: `NULL` means none if there are no
  zeros, else half the smallest non-zero value with a warning. The output
  of any `replace_zeros_*()` has no zeros, so it needs no pseudocount, and
  is the recommended input. (`log_ratio()` is the one exception, defaulting
  to `0`; see its plan.) `validate_pseudocount()` is currently keyed to
  `norm == 'clr'`; generalise it (e.g. a `force = TRUE` argument or a
  sibling `validate_log_pseudocount()`) rather than faking `norm`.
- **Orientation follows `margin`.** Transformed tables come back in the
  input's orientation, as `normalize()` does, so `margin` means the same
  on the way in and out.
- **`cpus` and determinism** as in README.md: random parts (`dirmult_ttest`,
  `qmp_rarefy`, the gemelli ports' initialisations) take `seed = 0L` and
  use PCG.


## Output structures (so results feed each other)

| Kind | Class | Functions | Feeds into |
|:-----|:------|:----------|:-----------|
| Zero-free table | dense base `matrix` of proportions, input orientation and dimnames; attributes `method` and the parameters used | `replace_zeros_*()` | any function taking `counts`, with `norm = 'clr'`, `aitchison()`, `alr()`, `ilr()`, `philr()`, `vlr()`, ... |
| Absolute-abundance table | the input's class (sparse stays sparse) | `qmp_scale()`, `qmp_rarefy()` | any function taking `counts` |
| Log-ratio coordinates | dense base `matrix`, input orientation; columns (or rows) are the new coordinates; attributes `method`, `basis` / `denominator`, `pseudocount` | `alr()`, `ilr()`, `philr()` | `pca()`, `rda()`, `euclidean()` / `beta_div(norm = 'none')`, `envfit()`, `lda()` ... (all accept negative values with `norm = 'none'`) |
| Per-sample values | named numeric vector; attributes `method`, `numerator`, `denominator` | `log_ratio()`, `balance()` | `groups` of any test, a column of `data`, `envfit()` |
| Sample distances | `dist` | `philr_distance()` | every `dm` function, `dist_ord()`, `dist_test()` |
| Feature associations | `dist` over **features** (packed lower triangle, as `beta_div()`), attribute `method` | `vlr()`, `rho_p()`, `phi_s()` | `pcoa()`/`hclust()` of features (`vlr`, `phi_s`); thresholding into networks (`rho_p`) |
| Ordinations | `ecodive_ord` from `ord_result()` (`method`, `eig`, `prop`, then extras such as `loadings`) | `rpca_phylo()`, `rpca_joint()`, `ctf()` | `envfit()`, `procrustes()`, plotting; `pcoa_biplot()`-style use of `loadings` |
| Per-feature results | `data.frame` of class `c('ecodive_<name>', 'data.frame')`, one row per feature (per feature × covariate for model-based methods, with a `covariate` column second); first column `feature`; then the method's statistics; then `p_value`, `p_adjusted` (where defined) and `reject` (logical) | `ancom()`, `dirmult_ttest()`, `ancombc()`, `ancombc2()` and their post-hoc tests | filtering, plotting; the shared columns let a future wrapper combine them |

`replace_zeros_*()` return proportions (each sample sums to 1), as
scikit-bio's `multi_replace()` and zCompositions' default `output =
'prop'` do. Every downstream log-ratio method is scale invariant, so
proportions lose nothing.


## Wrappers

- **New: `replace_zeros(counts, method, ...)`** with registry
  `ZERO_REPLACEMENTS` and `list_zero_replacements()` /
  `match_zero_replacement()`, shaped exactly like `count_ord()` (see
  [zero_replacement.md](zero_replacement.md)). Justified because five
  methods share one signature and one output, and the literature asks for
  sensitivity analyses across zero treatments, which is a loop over
  `list_zero_replacements()`. Its selector is `method`, like `dist_ord()`
  and `count_ord()`, so add a row to README.md's shadowing table:
  methods must not take `method` or a prefix of `method` / `counts`.
- **`count_ord()`:** `rpca_phylo()` joins (a count table plus `tree`,
  which passes through `...`). Add a `COUNT_ORDS` row.
- **`beta_div()`:** `philr_distance()` joins `METRICS` as a phylogenetic
  beta metric. `beta_div()` passes arguments with
  `mget(metric$params, environment())`, so every formal of a metric must
  also be a formal of `beta_div()`. See [philr.md](philr.md) for how to
  handle `part_weights` and `ilr_weights`.
- **No wrapper:**
  - `alr()`, `ilr()`, `philr()`: each needs a different extra input
    (`denominator`, `basis`, `tree`).
  - `log_ratio()`, `balance()`: two functions; numerator/denominator are
    their whole interface.
  - `vlr()`, `rho_p()`, `phi_s()`: three functions with one signature, so
    a wrapper is possible (`propr()` is one). Not yet: revisit if more
    feature-association measures arrive.
  - `ancom()`, `dirmult_ttest()`, `ancombc()`, `ancombc2()`: per-feature
    tables, which README.md keeps out of wrappers, and their inputs differ
    (`groups` vs `formula`). Their shared columns (above) leave room for a
    differential-abundance wrapper later if the user wants one.
  - `qmp_scale()`, `qmp_rarefy()`: need `load`.
  - `rpca_joint()` (several tables) and `ctf()` (subjects × times): two
    inputs or a different output; README.md keeps those as plain
    functions.


## Existing functions checked for splitting (2026-11-19)

The rule (one function per named method variant) was checked against the
functions these plans touch:

| Function / argument | Decision | Why |
|:--------------------|:---------|:----|
| `normalize(norm)` and every `norm` argument | **Keep** | `norm` is a preprocessing option shared by every count-table function, already in the "kept as arguments" table in TODO.md (`pca(norm)`). `normalize()` is its wrapper, as `beta_div()` is for metrics. |
| `pseudocount` | **Keep numeric** | Don't turn it into a selector (`pseudocount = 'gbm'`): that would put method modes inside an argument of ~30 functions. Zero replacement gets its own functions instead, whose output needs no pseudocount. |
| `aitchison()` / `robust_aitchison()` | Already split | |
| `rpca(distances)` | Keep | An extra output, already decided. |
| `rarefy(times)`, `rarefy(drop)` | Keep | Repetition and handling of short samples, not method variants. `rarefy(depth)` gains per-sample depths (see [qmp.md](qmp.md)), still one method. |

Nothing needs splitting.


## Docs

- A new pkgdown reference group, **"Compositional Data"**, after
  "Beta Diversity", with subtitles:
  - *Zeros*: `replace_zeros`, `replace_zeros_*`, `list_zero_replacements`,
    `match_zero_replacement`;
  - *Log-ratio transforms*: `alr`, `ilr`, `philr`, and a pointer to
    `normalize` (clr, rclr);
  - *Log ratios*: `log_ratio`, `balance`;
  - *Proportionality*: `vlr`, `rho_p`, `phi_s`;
  - *Differential abundance*: `ancom`, `ancombc`, `ancombc2`, the
    `ancombc_*` post-hoc tests, `structural_zeros`, `dirmult_ttest`;
  - *Absolute abundance*: `qmp_scale`, `qmp_rarefy`.
  `philr_distance` goes under Beta Diversity → Phylogenetic metrics; the
  `rpca_*` and `ctf` under Ordination.
- Consider a short vignette article, "Compositional analysis with
  ecodive", once several of these exist: zeros → transform → ordinate /
  test, and where QMP fits instead.
- **Same names in other packages** (each Rd page says so): `alr`, `ilr`
  (compositions, robCompositions), `philr` (philr: same meaning),
  `ancom` (scikit-bio, not R), `balance` (none known).


## References

- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
  (2017). Microbiome datasets are compositional: and this is not optional.
  *Frontiers in Microbiology*, 8, 2224.
- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
