# Compositional methods: overview and shared conventions (2026-11-19; revised 2026-11-20)

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

**Revision of 2026-11-20.** A literature review (2026-11-19; report in
`reports/Compositional methods for ecodive.md`, notes in
`research_notes/Compositional methods for ecodive/`, both kept out of git)
measured which compositional methods recent microbiome papers use and cite,
and which new methods appeared in 2023–2026. It overturned three decisions
below (SparCC, LinDA, MaAsLin3), added a second batch of plans (mostly
visualizations, differential-abundance helpers, networks and community
typing), and recorded two findings about methods ecodive already has. See
[The 2026-11-20 review](#the-2026-11-20-review) at the end.


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
| [coda_lasso.md](coda_lasso.md) | `coda_lasso()`, `coda_lasso_pairwise()` (added 2026-11-19) | port (zeroSum, MIT; coda4microbiome, MIT) |
| [philr.md](philr.md) | `philr()`, `philr_distance()` | port (basis, scikit-bio / gneiss) + formula (weights; philr output-only) |
| [rpca_phylo.md](rpca_phylo.md) | `rpca_phylo()` | port (gemelli) on ecodive's `rpca()` |
| [rpca_joint.md](rpca_joint.md) | `rpca_joint()` | port (gemelli) |
| [ctf.md](ctf.md) | `ctf()` | port (gemelli) |

**Second batch (2026-11-20).** Documentation and statistics:

| Plan | Functions | Tier |
|:-----|:----------|:-----|
| [compositional_docs.md](compositional_docs.md) | help-page and article changes (Schloss on RPCA and depth, scale-model guidance, ML caveats, predicted loads) | docs |
| [da_plots.md](da_plots.md) | `plot(type = 'effect' / 'ma' / 'forest')` for DA results; new dispersion and abundance columns in `dirmult_ttest()` | formula (CoDaSeq, MIT; ALDEx2 output-only) |
| [da_test.md](da_test.md) | `da_test()` wrapper (+ `list_da_tests()`, `match_da_test()`); `da_wilcoxon()`, `da_ttest()`, `da_lm()` | formula |
| [rarefied_diversity.md](rarefied_diversity.md) | `beta_div_rarefied()`, `alpha_div_rarefied()` | derived (on `rarefy()`) |
| [linda.md](linda.md) | `linda()` | clean room (LinDA GPL) |
| [sparcc.md](sparcc.md) | `sparcc()` | clean room (from the CC-BY paper) |
| [assoc_network.md](assoc_network.md) | `assoc_network()`, `assoc_network_compare()`, `assoc_network_diff()` | clean room (NetCoMi GPL) |
| [da_sensitivity.md](da_sensitivity.md) | `da_sensitivity()` | derived (on `da_test()`) |
| [maaslin3.md](maaslin3.md) | `maaslin3()` (**scope decision for the user**) | port (maaslin3, MIT) |
| [community_typing.md](community_typing.md) | `enterotype()`, `valencia()`, `dmm()` | clean room (PAM; DMM) + port (VALENCIA, MIT) |
| [dirmult_lm.md](dirmult_lm.md) | `dirmult_lm()` | port (ALDEx3, MIT) |
| [normalization_factors.md](normalization_factors.md) | new `norm` choices (`'css'`, `'tmm'`, `'rle'`, `'rle_poscounts'`, `'gmpr'`); `size_factors()` | port (CSS, Artistic-2.0) + clean room |
| [fastancom.md](fastancom.md) | `fastancom()` | clean room (GPL) |
| [spiec_easi.md](spiec_easi.md) | `spiec_easi_mb()`, `spiec_easi_glasso()` | clean room (SpiecEasi GPL) |
| [secom.md](secom.md) | `secom_pearson()`, `secom_distance()` | port (ANCOMBC, Artistic-2.0) |
| [rademu.md](rademu.md) | `rademu()` | port (radEmu, MIT) |

Visualizations:

| Plan | Functions | Tier |
|:-----|:----------|:-----|
| [biplot_scaling.md](biplot_scaling.md) | `plot(ord, scaling = )` for `pca()` / `rpca()` results, with links | formula |
| [log_ratio_plots.md](log_ratio_plots.md) | rank plots of DA effects and loadings (`plot(type = 'rank')`); `plot()` of `log_ratio()` / `balance()` results; `log_ratio_features()` | port (Qurro, BSD) |
| [ternary_plot.md](ternary_plot.md) | `ternary_plot()` | clean room |
| [coda_diagnostic_plots.md](coda_diagnostic_plots.md) | `variation_array()`, `zero_pattern()` | clean room |
| [composition_barplot.md](composition_barplot.md) | `composition_barplot()` | port (microshades idea, MIT) |
| [coda_dendrogram.md](coda_dendrogram.md) | `coda_dendrogram()` | clean room |
| [tree_plot.md](tree_plot.md) | `tree_plot()` | first principles |

The coda4microbiome-style signature plot is part of
[coda_lasso.md](coda_lasso.md) (its `plot()` method).

[TODO.md](TODO.md) gives the order.

**Not planned** (reasons in the method plans that mention them):

- **ALDEx2 itself:** `dirmult_ttest()` covers the ALDEx2 approach, and
  `dirmult_lm()` its regression form (via ALDEx3).
- **MaAsLin2:** log-TSS regression with no compositional correction; its
  successor MaAsLin3 corrects for compositionality and is planned instead
  (as a scope decision).
- **scikit-bio `dirmult_lme()`:** superseded by [dirmult_lm.md](dirmult_lm.md),
  which ports ALDEx3 (MIT), the source of ecodive's existing scale model;
  scikit-bio's closed-form random-intercept fit stays a reference for a
  later mixed-model step.
- **LOCOM, corncob, ZicoSeq, DACOMP:** GPL, with mixed benchmark standing;
  LinDA and fastANCOM cover the same ground.
- **Bayesian multinomial-logistic-normal models** (fido/pibble, PLN),
  **TEMPTED**, **compositional mediation** (CMM, LDM-med, CRAmed),
  **compositional deep learning**, **CoDaCoRe**, **BIRDMAn**: heavy
  machinery, GPL or unclear licences, or little uptake (see the review).
- **Porting the MLP microbial-load predictor** (Nishijima et al. 2025): an
  XGBoost model with an ambiguous "BSD license" string and no LICENSE
  file. Document how to pass its predicted load to `qmp_scale()` or a
  scale prior instead ([compositional_docs.md](compositional_docs.md)).
- **Heat trees, Sankey diagrams, bipartite "snowflake" plots:** need
  graph-layout engines and add little compositional value.
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
    ALDEx2 but isn't it, and ALDEx2's name shouldn't be borrowed (rule 6);
    `dirmult_lm()` is its regression sibling;
  - differential-abundance helpers share the `da_` prefix (`da_test()`,
    `da_wilcoxon()`, `da_ttest()`, `da_lm()`, `da_sensitivity()`). It is a
    prefix, so it doesn't clash with rule 3 (discriminant methods *end* in
    `da`). `da_ttest()` (Welch's test on one transformed table) and
    `dirmult_ttest()` (Welch's test on Monte Carlo draws) must each say how
    they differ from the other;
  - association networks are `assoc_network*()`, not `network()`, which
    would mask the network package's function;
  - plotting-only functions are `<thing>_plot()` (`ternary_plot()`,
    `tree_plot()`), or a noun with a `plot()` method when the object is
    worth keeping (`variation_array()`, `zero_pattern()`).
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
| `n_draws` | an exact number of Monte Carlo draws (`n_<unit>` rule) | `dirmult_ttest`, `dirmult_lm`, `sparcc` |
| `scaling` | how a biplot scales its scores and loadings (reserved for this in parameter_names.md) | `plot.ecodive_ord` |
| `threshold` | the association cutoff that makes an edge | `assoc_network` |
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
| Per-feature results | `data.frame` of class `c('ecodive_<name>', 'data.frame')`, one row per feature (per feature × covariate for model-based methods, with a `covariate` column second); first column `feature`; then the method's statistics; then `p_value`, `p_adjusted` (where defined) and `reject` (logical) | `ancom()`, `dirmult_ttest()`, `ancombc()`, `ancombc2()` and their post-hoc tests; planned `da_*()`, `linda()`, `maaslin3()`, `fastancom()`, `dirmult_lm()`, `rademu()` | filtering, plotting, `da_test()`, `da_sensitivity()` |
| Feature correlations | as feature associations (a similarity in a `dist` container, as `rho_p()`), with p-values where the method gives them | planned `sparcc()`, `spiec_easi_*()`, `secom_*()` | `assoc_network()` |
| Networks | see [assoc_network.md](assoc_network.md) | planned `assoc_network()` | `assoc_network_compare()`, `plot()` |
| Community types | see [community_typing.md](community_typing.md): one assignment per sample, plus per-type summaries | planned `enterotype()`, `valencia()`, `dmm()` | `groups` of any test; `plot()` |

**Plots.** Every result class gets a `plot()` method in base graphics that
returns its input invisibly (R/plot*.r). Alternative displays of the same
result are a `type` argument, as `plot.ecodive_mrpp(type = )` does.

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
- **Proposed (2026-11-20; needs the user's approval of the name):
  `feature_assoc(counts, method, ...)`** over `vlr()`, `rho_p()`,
  `phi_s()`, `sparcc()`, `spiec_easi_mb()`, `spiec_easi_glasso()`,
  `secom_pearson()` and `secom_distance()`: eight methods with one input and
  one output, and `assoc_network_compare()` must re-run whichever method
  made the networks on permuted samples. Not `assoc()`, which would mask
  `vcd::assoc()`. See [assoc_network.md](assoc_network.md). The shared
  output gains `similarity` and `features` flags in its record (set in
  [sparcc.md](sparcc.md)), replacing today's hard-coded `rho_p` checks.
- **New (2026-11-20): `da_test(counts, method, ...)`** with registry
  `DA_TESTS` and `list_da_tests()` / `match_da_test()`. The condition set
  here on 2026-11-19 ("if a third method joins") is met: `ancom()`,
  `dirmult_ttest()`, `ancombc()` and `ancombc2()` exist, and `linda()`,
  `maaslin3()`, `fastancom()`, `dirmult_lm()`, `rademu()` and the plain
  `da_*()` tests are planned. [da_test.md](da_test.md) settles how methods
  that take `groups` and methods that take `formula` share one wrapper.
  Add its row to README.md's shadowing table. `da_sensitivity()` runs any
  `da_test()` method over a grid of settings, building on
  `list_zero_replacements()` rather than beside it.
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
  - *Absolute abundance*: `qmp_scale`, `qmp_rarefy`;
  - *Associations and networks* (planned): `sparcc`, `spiec_easi_*`,
    `secom_*`, `assoc_network*`;
  - *Community types* (planned): `enterotype`, `valencia`, `dmm`;
  - *Plots* (planned): `ternary_plot`, `variation_array`, `zero_pattern`,
    `composition_barplot`, `coda_dendrogram`, `tree_plot`.
  `philr_distance` goes under Beta Diversity → Phylogenetic metrics; the
  `rpca_*` and `ctf` under Ordination; `beta_div_rarefied` and
  `alpha_div_rarefied` beside `rarefy`; `da_test` and the `da_*` helpers
  under *Differential abundance*.
- Consider a short vignette article, "Compositional analysis with
  ecodive", once several of these exist: zeros → transform → ordinate /
  test, and where QMP fits instead.
- **Same names in other packages** (each Rd page says so): `alr`, `ilr`
  (compositions, robCompositions), `philr` (philr: same meaning),
  `ancom` (scikit-bio, not R), `balance` (none known).


## The 2026-11-20 review

**What recent papers use** (OpenAlex citations 2022–26 / Europe PMC
microbiome full-text mentions 2022–26, retrieved 2026-11-19):

| Method | Citations | Mentions | In ecodive |
|:-------|----------:|---------:|:-----------|
| MaAsLin2 / MaAsLin3 | 2,850 / 246 | 2,284 / 155 | planned (`maaslin3()`) |
| ANCOM-BC / ANCOM-BC2 | 2,400 / 738 | 1,634 / 782 | yes |
| CLR transform | (Gloor 2017: 2,661) | ~3,000–3,460 | yes |
| SparCC | 1,848 | 1,363 | planned |
| SPIEC-EASI | 1,335 | 627 | planned |
| ALDEx2 | 1,038 + 685 | 1,178 | yes (`dirmult_ttest()`) |
| Aitchison distance | — | 1,036 | yes |
| NetCoMi | 651 | 327 | planned (`assoc_network*()`) |
| DMM | 506 | 328 | planned (`dmm()`) |
| RPCA ("robust Aitchison") | 501 | 247 | yes |
| LinDA | 494 | 122–664 (name collision) | planned |
| Ternary plots | — | 331 | planned |
| selbal / PhILR / Qurro / gneiss | 199 / 191 / 122 / 75 | 77 / 83 / 73 / 27 (all declining) | PhILR yes; Qurro-style plots planned |

Compositional practice is still a minority (Bray-Curtis is mentioned 28×
as often as Aitchison distance, LEfSe 8.6× as often as MaAsLin2), and CLR is
mentioned ~20× as often as ILR or ALR. Balance-selection tools are niche
and mostly declining. The fastest-growing methods (ANCOM-BC2, MaAsLin3,
LinDA) are refinements of linear models, and the sharpest recent papers
are benchmarks and critiques asking whether results survive changes of
settings. Hence the second batch's emphasis on displays that make results
legible (effect, MA, forest, rank, biplot scalings) and on sensitivity.

**Decisions of 2026-11-19 that the review changed:**

| Topic | Was | Now | Why |
|:------|:----|:----|:----|
| SparCC | not planned: "`rho_p()` answers the same question" | planned, clean room | That reasoning was wrong. ρp asks whether two taxa keep a constant ratio; SparCC estimates correlations between log *absolute* abundances, identifiable only under its sparsity assumption. Networks of the latter are what users report (SparCC: 1,848 citations vs propr's 159). The licence objection only rules out porting. |
| LinDA | not planned: GPL | planned, clean room | It is compositional by design (bias-corrected regression, ANCOM-BC's family), short, and the most robust method under strong compositional effects (Yang & Chen 2023). |
| MaAsLin3 | not planned: "general regression" | planned, as a scope decision | Right about MaAsLin2; MaAsLin3 corrects for compositionality and adds prevalence models, and its family has the highest uptake. Porting it widens ecodive into general per-feature regression: confirm with the user first. |
| A differential-abundance wrapper | later, "if a third method joins" | `da_test()` | Four methods exist. |
| `dirmult_lme()` | revisit after `dirmult_ttest()` | `dirmult_lm()` from ALDEx3 | `dirmult_ttest()` exists; ALDEx3 is MIT and the source of the scale model. |
| ANCOM-BC2 pseudocount sensitivity | left out (scikit-bio has none) | done upstream (`sensitivity`, on by default in `ancombc2()`) | The ANCOMBC developers default to it; Pelto et al. 2025 and a 2026 preprint report high false-discovery rates without it. |

**Findings about existing functions:**

- **Schloss (2026, mSystems, 10.1128/msystems.00930-26)** reports errors in
  the benchmark code that supported robust Aitchison PCA and finds that,
  corrected, it "was not insensitive to uneven sequencing effort and did
  not perform as well as rarefaction". `rpca()` and its siblings should
  cite it and the help pages shouldn't present RPCA as a depth correction
  ([compositional_docs.md](compositional_docs.md)); `beta_div_rarefied()`
  lets users run the rarefied alternative ([rarefied_diversity.md](rarefied_diversity.md)).
- Scale assumptions break FDR control in common tools (Dos Santos et al.
  2026; McGovern et al. 2023): recommend `scale_sd` sensitivity runs.
- CoDA transforms gave no advantage for random forests (Yerke et al.
  2024): don't promise predictive gains for clr, ILR or PhILR.

**Limits of the evidence:** OpenAlex rate limits (no citing-works
filter, so citations aren't restricted to microbiome papers); mentions are
not use; name collisions (LinDA, corncob, Songbird, gneiss, TEMPTED,
Gemelli); 2026 is a partial year; several findings rest on preprints.


## References

- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
  (2017). Microbiome datasets are compositional: and this is not optional.
  *Frontiers in Microbiology*, 8, 2224.
- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
  Chapman & Hall.
- Schloss, P. D. (2026). Rarefaction is better than robust Aitchison PCA
  and other compositional data analysis methods at controlling for uneven
  sequencing effort. *mSystems*. doi:10.1128/msystems.00930-26
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26, bbaf130.
  doi:10.1093/bib/bbaf130
- Wirbel, J., et al. (2024). *Genome Biology*.
  doi:10.1186/s13059-024-03390-9
