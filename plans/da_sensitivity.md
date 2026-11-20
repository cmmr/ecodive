# `da_sensitivity()`: rerun a differential-abundance method over settings

**Priority:** 34. **Tier:** derived (a loop over `da_test()` and
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
- Attributes: `method`, `effect_name`, `alpha`, `log_base`, `runs` (a
  data frame: `run`, `pseudocount`, `zeros`, `scale_sd`, `setting`),
  `fits` (the results, a list by run; documented as memory-heavy), and
  **`summary`**:
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


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- Each row of the tidy table names its setting in plain text (for example
  `"pseudocount 0.5"`, `"zeros: mult"`, `"scale_sd 1"`) so the tile plot's
  columns need no lookup; replacement settings take their names from the
  `ZERO_REPLACEMENTS` registry.
- `plot()` (`'tiles'`, `'traces'`): title "Sensitivity of " and the base
  method's display name from `DA_TESTS` (by the `method` attribute), for
  example "Sensitivity of ANCOM-BC2"; no `display` needed. `alpha` and
  `log_base` are copied from the base fit to top-level attributes of the
  result, for the effect axis label and the reject rule's legend.
- Bold row labels and filled points read `summary$passed` and the
  `reject` column; `plot()` recomputes nothing.

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


## Decisions made while implementing (2026-11-20)

Asked of the user at the start (all three answered with the
recommendation):

- **Grid names: plural** (`pseudocounts`, `zeros`, `scale_sds`, `crossed`),
  as the plan proposed. `test-wrapper_args.r` checks that
  `da_sensitivity()`'s arguments are `counts, method, ...` and these four,
  and that no DA method takes one of the four; README.md's shadowing table
  has a `da_sensitivity()` row.
- **`crossed = TRUE`: pseudocounts and zero replacements are alternatives**
  (one "zero treatment" factor), each crossed with every `scale_sds`; no run
  adds a pseudocount to a zero-replaced table. Runs are (number of
  pseudocounts + number of zero replacements) × number of scale_sds, plus
  the base (a grid not given counts as one level). With one grid, crossed
  and one-at-a-time give the same runs.
- **Printing: `print.ecodive_da_sensitivity <- print_table_result`** (the
  README's 2026-11-20 rule), header "Sensitivity of <registry name>"
  (`print_df_label()`). No `summary()` method: the summary stays an
  attribute, as planned.

Made while implementing:

- **Runs.** The base run is `da_test(counts, method, ...)`; every other run
  is the same call with its settings set by name in `...`. Names in `...`
  that partially match one of the method's arguments are completed first
  (`pseudo = 1` becomes `pseudocount = 1`), so a run's setting replaces the
  base's instead of being given twice. A run equal to the base (e.g.
  `pseudocount = 0.5` for `dirmult_ttest()`, whose default is 0.5) is run
  anyway, as ANCOMBC does.
- **Zero-replacement runs** give the method `replace_zeros(counts, z, margin
  = m)` with the base's `pseudocount` **dropped**, so the method's default
  applies; for every `zero_free_ok` method that default (`NULL`, or 0 for
  `ancombc2()`) adds nothing to a table without zeros (an explicit 0 is
  refused by the da_* tests' shared rule, so dropping it is the only form
  that works for all). Any input type is first converted as the methods
  convert it (`validate_counts()`, which also takes a phyloseq / rbiom /
  SummarizedExperiment object's orientation), and the run gets that
  `margin`.
- **ANCOM-BC reuse.** For `ancombc` / `ancombc2`, a pseudocount run (not
  crossed with anything) whose value is one of the base run's own
  `attr(fit, 'sensitivity')$pseudocounts` takes that stored fit; any other
  is fitted directly with `pseudocount = pc, sensitivity = FALSE`, which is
  the same call `ancombc_sensitivity()` makes. So with `ancombc2()`'s
  default (`sensitivity = TRUE`) and `pseudocounts = TRUE`, the pseudocount
  runs cost nothing. Every non-base ANCOM-BC run has `sensitivity = FALSE`.
- **Checks before running:** at least one grid; `pseudocounts` positive,
  finite, not repeated, and only for methods with a `pseudocount` formal;
  for methods with `norm`, refused unless the base norm (given or default)
  is `clr` or `log`, since every other norm ignores the pseudocount and all
  runs would be the same; `scale_sds` non-negative, only for methods with
  `scale_sd` (the error lists them, from the registry's functions); `zeros`
  `TRUE` or names matched with `match_zero_replacement()` (flexible
  matching, duplicates after matching refused), and refused for
  `zero_free_ok = FALSE` (the error lists the methods that can).
- **Warnings.** The base run's warnings are given as usual. The other runs'
  are collected and given once each after the loop, leaving out any the
  base run already gave (so a "Zeros detected" warning isn't repeated per
  run). An error in a run is re-raised as "The run with <setting> failed:
  <message>". `pseudocounts = TRUE` on a table that isn't whole numbers
  warns, as `ancombc2(sensitivity = TRUE)` does.
- **Setting labels** use the argument form: `"base"`, `"pseudocount =
  0.5"`, `"zeros = mult"` (the registry id), `"scale_sd = 1"`, and crossed
  runs joined by ", ". `parameter` is `'crossed'` only for runs that set
  two or more settings.
- **Table.** `covariate` is always a column (`NA` for the group methods),
  so every result has the same columns; intercept rows are dropped. Rows
  are run by run, each in the method's row order. `p_adjusted` and `effect`
  are `NA` for ancom.
- **Attributes:** `method`, `effect_name`, `alpha` (the base fit's),
  `log_base` (the base fit's `log_base` attribute when present, else from
  the registry's `base`: `'e'` → `exp(1)`, `'2'` → 2; so it is right
  before item 27 adds `log_base` to ancombc / dirmult results), `covariate`
  (model methods only: the first design column of `groups`'s term, else
  the first non-intercept column; `da_term_covariate()`, the rule of
  `plot.ecodive_da_lm()`), `runs`, `fits` and `summary`.
- **Summary** rows are the (feature, covariate) pairs in the base run's
  order. `n_runs` counts the runs with a row for the pair; `passed` is
  `n_reject` 0 or `n_runs`; `robust` is `passed` and rejected in the base
  run; `sign_stable` compares `sign()` of the finite effects (`NA` when none
  is finite, e.g. ancom).
- **Tiles plot:** cell codes drawn with `image()` (red / blue ends of
  `hcl.colors(11, 'Blue-Red 3')`, `grey90`, `grey35` for a rejection with
  no effect sign, i.e. ancom; white for no result), thin white lines
  between cells, black lines between parameter groups, the group names
  above, the setting's value below (the full setting for crossed runs).
  With no feature rejected in any run, the top `n_features` by base
  |effect| are shown (all grey) and the key says so. The key (colours,
  bold, counts, covariate, alpha) is an `mtext()` line under the labels;
  `sub =` replaces it.
- **Traces plot:** x on a log scale for pseudocounts (they span orders of
  magnitude; `log = ''` overrides), categorical for zeros; features are
  those rejected in the base run, else in any run, else all, largest base
  |effect| first; labels at the right end, bold when `passed`. With
  crossed runs, the other setting is held at its first value (named in the
  subtitle).
- **Title** "Sensitivity of " + the registry name, from
  `da_sens_plot_title()` (tested).
- **Docs:** cross-references from `da_test()`'s and `dirmult_ttest()`'s
  "See also", and a γ example in `dirmult_ttest()`'s examples, kept small
  because item 27 is editing those pages. Glossary row in
  parameter_names.md; NEWS; pkgdown (`da_sensitivity` under Differential
  abundance, its plot method under Plotting).

### Reference results

- **Each run equals the direct call** (`identical()` on the whole result,
  attributes included, so the table's `effect`, `p_adjusted` and `reject`
  too): ancom, ancombc, ancombc2, dirmult_ttest, da_wilcoxon, da_ttest and
  da_lm; base, pseudocount, zero-replacement and scale runs; dense,
  `Matrix` and transposed (`margin = 2`) input on a 30 × 200 table (135 of
  135 comparisons identical), plus the test file's 16 × 12 table.
- **ANCOM-BC2's `passed_ss`:** with pseudocounts 0.1, 0.5, 1, the
  pseudocount runs are `identical()` to `attr(ancombc2(..., sensitivity =
  c(0.1, 0.5, 1)), 'sensitivity')$fits` both when reused from the base run
  and when fitted directly (base `sensitivity = FALSE`), and
  `summary$passed` / `summary$robust` are identical to `passed_ss` /
  `reject_robust` (40 × 150 table, `~ grp` and `~ grp + age`: 20 base
  rejections, 19 robust, 1 not passed). Same for `ancombc()`.
- **Threads:** results `identical()` for `cpus` 1, 2 and 4 (ancom,
  ancombc2, dirmult_ttest, da_lm in the tests; every method on the larger
  table).
- No numerical differences anywhere: every comparison is exact.
