# `da_sensitivity()`: rerun a differential-abundance method over settings

**Priority:** 41. **Tier:** derived (a loop over `da_test()` and
`replace_zeros()`, with a tidy result and a `plot()` method). No reference
code is ported: ANCOMBC's sensitivity analysis is already in `ancombc()`
(Artistic-2.0 notices in place), and ALDEx2's `plotGamma` (GPL ≥ 3) is only
the idea for one display. Read [README.md](README.md) first, and
[compositional.md](compositional.md).


## What it is

A DA result can hinge on a setting that has no right value: the
pseudocount, the zero-replacement method, the assumed variation in total
load (`scale_sd`). The literature now treats checking this as part of the
analysis:

- **ANCOM-BC2.** The vignette "strongly recommend[s]" its pseudocount
  sensitivity analysis (pseudocounts 0.1, 0.5, 1; a taxon passes when its
  decision is the same in every run). Without it a 2026 preprint
  (Abdelmalek et al., 10.64898/2026.09.01.748571) saw FDR rise from about
  0.02 to about 0.82; Pelto et al. (2025) found ANCOM-BC2's results the
  least replicable of 14 methods. ecodive's `ancombc2()` already runs it by
  default (`sensitivity = TRUE`, commit 5856378).
- **Scale.** Nixon, Gloor & Silverman (2025; 27 OpenAlex citations,
  retrieved 2026-11-19) recast normalization as a scale assumption and
  recommend scale uncertainty; Dos Santos et al. (2026) show misspecified
  scale breaks FDR control and "implementing a scale model in ALDEx2 or
  ALDEx3 ameliorated unacceptably high FDRs"; McGovern et al. (2023) report
  positive predictive value as low as 9% under small scale errors.
  ALDEx2's `plotGamma` draws results against γ (= ecodive's `scale_sd`).
- **Zeros.** compositional.md designed `replace_zeros()` /
  `list_zero_replacements()` so that a zero-treatment sensitivity
  analysis is a loop; the review's reconciliation says to "build on the
  existing design, not beside it".

`da_sensitivity()` is that loop for any `da_test()` method, returning one
table and two plots (features × settings tiles; per-feature effect traces).


## Names

`da_sensitivity()`: the `da_` prefix of the DA helpers (compositional.md,
2026-11-20); the review's `sensitivity()` is too generic (and is a common
name in caret, epiR, ...). Not a registry member itself.


## Design

### Signature

```r
da_sensitivity(counts, method, ..., pseudocounts = NULL, zeros = NULL,
               scale_sds = NULL, crossed = FALSE)
```

- `counts`, `method`, `...`: as `da_test()`; `...` holds the method's
  arguments (`groups`, `data`, `formula`, `seed`, `cpus`, ...) and defines
  the **base run**.
- **The grids are plural** on purpose: arguments after `...` match only
  by exact name, so a singular `pseudocount` here would capture the
  method's own `pseudocount` and the base run could not set it. Plurals say
  "a vector of values to try", leave `pseudocount = 1` and `scale_sd = 0.5`
  free to pass through to the base run, and add no shadowing rule (no DA
  method may take an argument named exactly `pseudocounts`, `zeros`,
  `scale_sds` or `crossed`; add that to the wrapper-args test). **Ask the
  user**: the brief suggested singular names.
  - `pseudocounts`: `TRUE` = `c(0.1, 0.5, 1)` (ANCOMBC's), or positive
    values. For `dirmult_ttest()` / `dirmult_lm()` it varies the Dirichlet
    prior (documented: a different meaning, same name).
  - `zeros`: `TRUE` = `list_zero_replacements('id')`, or ids. Each run is
    `method(replace_zeros(counts, z), ...)`, whose output needs no
    pseudocount. Refused (error) for methods with `zero_free_ok = FALSE`
    in `DA_TESTS`: `dirmult_ttest`, `dirmult_lm` (the Dirichlet posterior
    needs counts), `maaslin3` (its prevalence model needs zeros), `rademu`.
  - `scale_sds`: `TRUE` = `c(0, 0.25, 0.5, 1)`, or non-negative values;
    only for methods with a `scale_sd` formal (else an error naming them).
  - At least one grid, else an error.
- `crossed = FALSE`: one setting varied at a time around the base; `TRUE`:
  the full factorial of the given grids.

### Runs

1. The base run: `da_test(counts, method, ...)`.
2. One run per setting, the base call plus the setting. For runs other
   than the base, `ancombc()` / `ancombc2()` get `sensitivity = FALSE` (or
   each run would refit three more times).
3. **Reuse, don't duplicate, ANCOM-BC's analysis**: for `ancombc` /
   `ancombc2` and `pseudocounts`, make one call with `sensitivity =
   pseudocounts` and take `attr(fit, 'sensitivity')$fits`; the numbers are
   then those behind `passed_ss` (an `identical()` test enforces it).
4. Random methods keep the base `seed` in every run, so differences
   reflect settings, not Monte Carlo noise.
5. Runs are sequential; each uses `cpus` internally. Results are
   identical for any `cpus` because each run is.

### Result

`data.frame` of class `c('ecodive_da_sensitivity', 'data.frame')`, long
format, one row per feature × [covariate ×] run:

- `feature`, `covariate` (model methods: every non-intercept covariate;
  the plot picks one), `run` (integer, 0 = base), `parameter`
  (`'base'`, `'pseudocount'`, `'zeros'`, `'scale_sd'`, or for crossed runs
  `'crossed'`), `setting` (a label, `"pseudocount = 0.5"`), `effect` (the
  registry's `effect` column; `NA` for `ancom`), `p_adjusted` (`NA` for
  `ancom`), `reject` (the method's own `reject`; for ANCOM-BC not
  `reject_robust`).
- Attributes: `method`, `effect_name`, `log_base`, `runs` (a data frame:
  `run`, `pseudocount`, `zeros`, `scale_sd`, `setting`), `fits` (the
  results, a list by run; documented as memory-heavy), and **`summary`**:
  per feature [× covariate] `n_runs`, `n_reject`, `prop_reject`,
  `passed` (the decision is the same in every run, ANCOMBC's `passed_ss`
  generalised), `robust` (`passed & base reject`), `sign_stable` (effect
  has one sign among runs where it is finite). The ANCOMBC2 vignette's
  "sensitivity score" (the proportion of pseudocount runs whose q-value
  exceeds alpha) is `1 - prop_reject` over the pseudocount runs; document
  the mapping.
- No `print()` or `summary()` method: the summary stays an attribute, as
  `ancom()` keeps its extras.

### `plot()`

```r
plot.ecodive_da_sensitivity(x, type = c('tiles', 'traces'),
                            covariate = NULL, parameter = NULL,
                            n_features = 30, ...)
```

- `'tiles'` (default): `image()` of features (rows: those rejected in any
  run, by `prop_reject` then |base effect|, top `n_features`) × runs
  (columns, grouped by `parameter` with vertical separators, base first).
  Cell colour: rejected and positive (red end of `hcl.colors(..., 'Blue-Red
  3')`), rejected and negative (blue end), not rejected (`grey90`), `NA`
  (white). Row labels in bold when `passed`. Margins as
  `ancombc_plot_heatmap()`.
- `'traces'`: for one `parameter` (default: the first numeric grid
  given), x = the setting's value (categorical positions for `zeros`), y =
  `effect`; a line per feature (the `n_features` rejected in the base run,
  else in any run), filled points where rejected, hollow otherwise; dotted
  `h = 0`. For `scale_sd` this is the `plotGamma`-style display, written
  from the idea only.
- `covariate`: as `plot.ecodive_ancombc()` (default: the term of
  interest). Returns `x` invisibly.


## Reference implementations

- ecodive's `ancombc_sensitivity()` (R/ancombc.r): the analysis it
  generalises, already carrying the ANCOMBC (Artistic-2.0) notices; no new
  porting.
- ANCOMBC2 vignette (Bioconductor; the sensitivity-score definition) and
  ALDEx2's manual (`plotGamma`, `aldex.senAnalysis`): documentation only.
- No fixtures from other packages: correctness is "each run equals the
  direct call".


## Step records and plotting

Follow [compositional.md](compositional.md) § "Step records and what plots read" (records: `R/steps.r`, `?ecodive_steps`).

- **`sensitivity` record** "Sensitivity Analysis" (`n_settings`,
  `crossed`, and which grids were given), plus the base method's `da`
  record, taken from the first fit, beside it.
- Each row of the tidy table names its setting in plain text (for example
  `"pseudocount 0.5"`, `"zeros: mult"`, `"scale_sd 1"`) so the tile plot's
  columns need no lookup; replacement settings take their names from the
  `ZERO_REPLACEMENTS` registry, or from the replaced table's `transform`
  record.
- `plot()` (`'tiles'`, `'traces'`): title `step_title("Sensitivity of",
  attr(x, 'da'))`; `alpha` and `log_base` come from the base method's
  attributes, stored on the result.

## Tests

- Every run's rows equal the direct call's (`identical()` on `effect`,
  `p_adjusted`, `reject`) for `ancom`, `ancombc2`, `dirmult_ttest`,
  `da_lm` (and each later method as it joins).
- ANCOM-BC pseudocount runs `identical()` to
  `attr(ancombc2(..., sensitivity = c(...)), 'sensitivity')$fits`, and
  `summary$passed` equals `passed_ss` when the grids are the same.
- `crossed = TRUE` gives the product of grid lengths runs (+ base);
  `zeros = TRUE` gives one run per `list_zero_replacements('id')`.
- `identical()` results for `cpus` 1, 2, 4.
- Errors: no grid; `scale_sds` for a method without `scale_sd`; `zeros`
  for `dirmult_ttest`; negative pseudocounts; unknown zero method.
- Both plot types draw on `pdf(NULL)`; `'traces'` for `ancom` errors (no
  effect). Smoke line (`da_sensitivity()` with `pseudocounts = TRUE` on a
  small `da_ttest`).


## Docs

Rd page in *Differential abundance*: why (the evidence above, briefly);
what each grid means per method (the pseudocount's three meanings: zero
patch, every-value pseudocount in ANCOM-BC, Dirichlet prior); how to read
both plots; that `passed` is ANCOMBC's `passed_ss` generalised. Add a
γ-sensitivity example to `dirmult_ttest()`'s docs pointing here. Glossary:
`pseudocounts`, `zeros`, `scale_sds`, `crossed` (new). NEWS.


## References

- Lin, H., & Peddada, S. D. (2024). *Nature Methods*, 21, 83–91.
  doi:10.1038/s41592-023-02092-7; ANCOMBC2 vignette (Bioconductor),
  "pseudo-count sensitivity analysis".
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating
  scale uncertainty in microbiome and gene expression analysis as an
  extension of normalization. *Genome Biology*, 26, 139.
  doi:10.1186/s13059-025-03609-3 (27 OpenAlex, retrieved 2026-11-19).
- Dos Santos, S. J., Murariu, A. C., Silverman, J. D., & Gloor, G. B.
  (2026). Not every gene is special: modelling scale controls the false
  discovery rate. *PLoS Computational Biology*, 22(9), e1014728.
  doi:10.1371/journal.pcbi.1014728 (2 OpenAlex, retrieved 2026-11-19).
- McGovern, K. C., Nixon, M. P., & Silverman, J. D. (2023). Addressing
  erroneous scale assumptions in microbe and gene set enrichment analysis.
  *PLoS Computational Biology*, 19(11), e1011659.
  doi:10.1371/journal.pcbi.1011659.
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130.
- Abdelmalek et al. (2026). A configuration-resolved benchmark of
  differential abundance analysis methods for human gut 16S rRNA microbiome
  data. *bioRxiv* (preprint, not peer reviewed).
  doi:10.64898/2026.09.01.748571.
