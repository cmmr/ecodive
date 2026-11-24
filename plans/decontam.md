# `decontam()`, `decontam_not()`: contaminant identification

**Priority:** 40. **Tier:** port (decontam 1.32.0, Bioconductor 3.23,
**Artistic-2.0**: `License: Artistic-2.0` on the Bioconductor 3.23 package
page and in the tarball's `DESCRIPTION`, which has no LICENSE file;
authors Benjamin Callahan (aut, cre), Nicole Marie Davis (aut), Felix G.M.
Ernst (ctb); verified 2026-11-20). Its R code (`R/decontam.R`, ~300 lines)
may be read and ported, keeping the Artistic-2.0 notice and stating how
the port differs (README § Licensing). It calls base R's `lm()`, `pf()`,
`prop.test()`, `fisher.test()` and `pchisq()`: those are R (GPL) and
**output-only**; their formulas are textbook and are confirmed from
outputs. Read [README.md](README.md) first.


## What it is

Reagent and laboratory DNA contaminates every sequencing run, and in
low-biomass samples it can be most of the reads (Salter et al. 2014;
Eisenhofer et al. 2019). decontam (Davis et al. 2018) scores each feature
with two signatures of contamination:

- **frequency**: a contaminant's relative abundance falls as the sample's
  total DNA concentration rises (it is a roughly constant amount of DNA
  diluted by more sample DNA), so its frequency goes as \eqn{1/C};
- **prevalence**: a contaminant is more prevalent in negative controls
  than in true samples.

The two can be combined, and are run per batch (sequencing run, extraction
plate). `isNotContaminant()` reverses the prevalence test's null for
low-biomass studies where most features may be contaminants.

The Fierer et al. (2025) consensus guidelines for low-biomass studies
(*Nature Microbiology*) ask for negative controls at every stage and for
reporting the contaminant-removal workflow; decontam is the tool most
papers then name. Demand (retrieved 2026-11-20): Davis et al. (2018) 3,160
OpenAlex citations 2022–26, rising from 478 in 2022 to 751 in 2026; 1,617
Europe PMC microbiome full-text mentions 2022–26. The method is simple
(per-feature two-sample statistics) and the port is small.


## Names

- **`decontam()`.** Rule 1: the short name users know. Methods sections
  say "contaminants were removed with decontam (prevalence method)"; the
  package exports `isContaminant()`, not `decontam()`, so nothing is
  masked and the name points users to the right reference. Spelled-out
  alternatives considered: `find_contaminants()` (clear, but sorts away
  from its companion and isn't what users search for) and `contaminants()`
  (a noun; reads as a dataset). Recommended: `decontam()`; **ask the
  user**.
- **`decontam_not()`** for `isNotContaminant()` (rule 2, `<base>_<modifier>`;
  one function per named variant, as decontam itself splits them). It
  mirrors "is not contaminant"; an alternative is
  `decontam_noncontaminant()`. **Ask.**
- Plain functions (no wrapper): not a DA method; no `DA_TESTS` row. pkgdown:
  a new "Quality control" group, or "Compositional Data" → first item.
- **`method`** keeps decontam's name and values (`'auto'`, `'frequency'`,
  `'prevalence'`, `'combined'`, `'minimum'`, `'either'`, `'both'`): the
  single reference treats them as one argument, and with no wrapper
  `method` is free (as `mantel(method)`).


## Design

### Signatures

```r
decontam(counts, conc = NULL, neg = NULL, data = NULL, margin = 1L,
         method = c('auto', 'frequency', 'prevalence', 'combined',
                    'minimum', 'either', 'both'),
         batch = NULL, batch_combine = c('minimum', 'product', 'fisher'),
         threshold = 0.1, norm = 'percent', cpus = n_cpus())

decontam_not(counts, neg, data = NULL, margin = 1L, threshold = 0.5,
             norm = 'percent', cpus = n_cpus())
```

- `conc`, `neg`, `batch`: one value per sample, or a column name of `data`,
  by the shared `per_sample()` / `match_rows()` rule (as `qmp_scale(load)`
  and `groups`). Matched by name, else by position.
- **`conc`** keeps decontam's name. The glossary's `load` (qmp) is a
  measured *total microbial load* used as a scale; decontam needs any
  quantity proportional to the input DNA (fluorometer reading of the PCR
  product, qPCR), and only its inverse relation matters. Reusing `load`
  would invite users to pass a load where decontam's docs ask for a DNA
  concentration. **Ask** whether to reuse `load` anyway.
- **`neg`**: logical (`TRUE` = negative control). `NA` samples are left
  out of both tests (decontam's documented intent; see "Differences").
- **`threshold`** keeps decontam's name and its **strict `<`** rule.
  ecodive's `alpha` would mislead: the frequency "p" is a heuristic score
  (an F ratio with equal degrees of freedom, below), not a p-value, and
  Davis et al. call both "scores". The glossary's `threshold` (an
  association cutoff) is a different quantity; note it in
  `parameter_names.md`. A length-2 `threshold` is the frequency and
  prevalence cutoffs for `'either'` / `'both'`.
- **`norm`** replaces decontam's `normalize = TRUE/FALSE`: `'percent'`
  (default; each sample divided by its total) or `'none'` (frequencies or
  equal-depth counts as given). Other norms are an error: both tests need
  non-negative relative abundances.
- `batch_combine`: decontam's `batch.combine`, snake case.
- Not carried over: `detailed = FALSE` (the result is always the table;
  `x$contaminant` is the logical vector), phyloseq `sample_data` lookups
  by name (ecodive takes `data`), `isContaminant(method = "chisq" /
  "fisher")` in the internal prevalence helper (not reachable from the
  exported function).

### Algorithm (ported from decontam 1.32.0 `R/decontam.R`)

0. Samples with zero total are removed with a warning naming how many
   (decontam's message). Normalize (`norm`). `method = 'auto'`:
   `'frequency'` with `conc` only, `'prevalence'` with `neg` only,
   `'combined'` with both. Frequency needs every `conc > 0` (all samples,
   controls included, as decontam checks). Batches: at least 2 samples
   each (error), warning at 4 or fewer.
1. **Frequency score** (`isContaminantFrequency()`), per batch, over the
   batch's non-control samples where the feature is present
   (\eqn{f > 0}), \eqn{k} of them: \eqn{y = \log f}, \eqn{l = \log c},
   \eqn{SS_0 = \sum (y - \bar y)^2} (intercept-only model: non-contaminant,
   frequency independent of \eqn{c}), \eqn{SS_1 = \sum ((y + l) -
   \overline{(y + l)})^2} (slope fixed at \eqn{-1}: contaminant),
   score \eqn{= F_{k-1, k-1}(SS_1/SS_0)} (lower tail, `pf()`). \eqn{k \le 1}
   → `NA`. In C, two passes over each feature's non-zero values (means,
   then sums of squares) give `lm()`'s residual sums to rounding.
2. **Prevalence score** (`isContaminantPrevalence()`), per batch, needing
   at least 2 samples with the feature and at least one control and one
   non-control: the 2 × 2 table (control / not) × (present / absent);
   - present in every sample of the batch → 0.5;
   - otherwise the one-sided "greater" (more prevalent in controls)
     two-proportion test of `prop.test()` with its continuity correction;
   - if that test would warn (✱ the chi-squared approximation warning:
     confirm from outputs that it is "an expected count below 5"), the
     Fisher mid-p instead: \eqn{P(X \ge x) - P(X = x)/2} under the
     hypergeometric with the table's margins (decontam's
     `fisher.pval()`: `p_greater - (p_greater + p_less - 1)/2`).
   ✱ Write `prop.test()`'s one-sided statistic from its documentation
   ("continuity correction is used only if it does not exceed the
   difference between sample and null proportions") and textbook formulas;
   confirm on a grid of tables (all margins up to 30 × 30, with and
   without the correction capping) that it matches to 1e-14. Vectorized
   over features in R from the C counts (`stats::phyper()`, `dhyper()`,
   `pnorm()`), not a loop.
3. **Combine batches** per feature: `'minimum'` (min over batches with a
   score; `NA` if none), `'product'` (product of the non-`NA` scores; **1
   when all are `NA`**, as decontam's `prod(na.rm = TRUE)`: keep the
   quirk, document it), `'fisher'` (`NA` → 0.5, then \eqn{P(\chi^2_{2B}
   > -2 \sum \log s)}).
4. **Overall score** `p`: the frequency or prevalence score;
   `'minimum'`: the smaller; `'combined'`: Fisher's method on the two,
   \eqn{P(\chi^2_4 > -2 \log(s_f s_p))}; `'either'`, `'both'`: `NA`.
5. **Call**: `contaminant = p < threshold` (`'either'`: \eqn{s_f < t_1}
   or \eqn{s_p < t_2}; `'both'`: and); `NA` → `FALSE`.
6. **`decontam_not()`** (`isNotContaminant()`): the prevalence score
   \eqn{s_p} as above, then \eqn{p = 1 - s_p}; `not_contaminant = p <
   threshold` (default 0.5); `NA` → `FALSE`. Prevalence only, as decontam.

### C and threads

`src/decontam.c`: one pass over the table by feature (transpose once, as
`da_table_init()` does; never dense), threaded over blocks of features.
Per feature and batch it returns \eqn{k}, \eqn{SS_0}, \eqn{SS_1} and the
four prevalence counts. Sums run over samples in sample order within a
feature, so results are identical for any `cpus`. `pf()`, the
two-proportion test, `phyper()` and `pchisq()` run vectorized in R.
Cost: O(nnz), trivial next to the I/O.

### Result

`c('ecodive_decontam', 'data.frame')`, one row per feature, decontam's
`detailed = TRUE` columns in snake case plus the counts behind the
prevalence plot: `feature, freq` (mean normalized abundance over all
samples, controls included, as decontam), `prev` (samples present),
`prev_neg`, `prev_pos` (present among controls / non-controls; `NA`
without `neg`), `p_freq, p_prev, p, contaminant`. `decontam_not()`:
`feature, freq, prev, prev_neg, prev_pos, p_prev` (\eqn{1 - s_p}), `p`,
`not_contaminant` (class `c('ecodive_decontam_not', 'data.frame')`).

Attributes: `method` (the function id, `'decontam'` or `'decontam_not'`),
`test` (the method used after `'auto'`), `threshold`, `batch_combine`,
`n_neg`, `n_pos` (control and non-control samples used), and, with
batches, `p_freq_batch` / `p_prev_batch` (batch × feature matrices, as
decontam's internal `p.freqs`). Printing: `print_table_result`, header
"decontam (prevalence)".

Docs show the removal one-liner: `counts[, !res$contaminant]` (and
`margin = 2`), since decontam leaves removal to the user.


## Reference implementations

- **Port from:** decontam 1.32.0 (Artistic-2.0). Header of `src/decontam.c`
  and `R/decontam.r`: the Artistic License 2.0 text, "Copyright (c)
  Benjamin Callahan and Nicole Marie Davis" (**verify** the holder named in
  the tarball; the DESCRIPTION has authors but no copyright line), and "C
  and R port of `isContaminant()`, `isNotContaminant()`,
  `isContaminantFrequency()`, `isContaminantPrevalence()` and
  `fish.combine()` from decontam 1.32.0; restructured for ecodive, see
  plans/decontam.md". `inst/COPYRIGHTS` stanza (format as ANCOMBC's);
  `Authors@R`: Benjamin Callahan and Nicole Marie Davis (`ctb`, `cph`).
- **Fixtures** (`data-raw/expected/decontam.r`): simulated tables only
  (decontam's `inst/extdata/st.rds` is Artistic-2.0 data; don't copy it
  into the tests): known contaminants (frequency \eqn{\propto 1/c}, more
  prevalent in controls), every `method`, two batches, all three
  `batch_combine`s, a length-2 `threshold`;
  `isNotContaminant()`; crafted 2 × 2 tables for the prop.test / Fisher
  branch.


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `plot.ecodive_decontam(x, type = c('scores', 'prevalence'))`:
  - `'scores'` (default): histogram of `p` (or of the score of the test
    run), with `threshold` as a dashed line. Davis et al. choose the
    threshold from this histogram's bimodality.
  - `'prevalence'`: `prev_neg / n_neg` against `prev_pos / n_pos`, one
    point per feature, contaminants highlighted (the decontam vignette's
    presence/absence scatter).
- Title: "Contaminant scores (decontam, prevalence)"; the input is a count
  table, so no `display`. Reads `threshold`, `test`, `n_neg`, `n_pos`,
  `contaminant`.
- decontam's `plot_frequency()` (a feature's frequency against `conc`, with
  both fitted lines) needs the per-sample values, which the result does not
  keep, and README forbids `plot()` from recomputing from data. **Ask**
  whether to add it as `plot(x, type = 'frequency', counts = , conc = ,
  features = )` (an exception to the rule), as a separate
  `decontam_frequency_plot()`, or not at all.


## Tests

- Every column equals decontam's at 1e-12 (scores) and `identical()`
  (calls, `prev`), for each `method`, `batch_combine`, `norm = 'none'`,
  and `decontam_not()`.
- The prevalence test's two branches, compared on the crafted tables with
  `prop.test()` and `fisher.test()` directly (base R; allowed in tests).
- `identical()` for `cpus` 1, 2, 4; dense = `Matrix` = `slam`;
  `margin = 2`.
- Edge cases: a zero-total sample (removed, warning); a feature present in
  one sample (`NA` scores, not a contaminant); all samples controls or none
  (prevalence `NA`); a batch of one sample (error) or four (warning);
  `conc` with a zero or `NA` (error); `NA` in `neg` (left out); `threshold`
  outside [0, 1]; length-2 `threshold` with `'frequency'` (error).
- Smoke lines for `decontam`, `decontam_not`.


## Docs

- Rd page (one for both): the two signatures of contamination; that the
  frequency score is a heuristic (equal-df F ratio), not a p-value; how to
  choose `threshold` (the score histogram; decontam's default 0.1, 0.5 for
  aggressive removal); controls should be extraction blanks (decontam's
  docs); a note on Fierer et al. (2025): controls at every stage, report
  the workflow, and statistical removal does not undo cross-contamination
  between samples.
- **"Differences from decontam"** (Artistic-2.0 §4): `NA` in `neg` leaves
  those samples out of the frequency test too (decontam's frequency branch
  subsets with `NA` and returns `NA` for every feature); `norm` instead of
  `normalize`; no phyloseq lookups; the result keeps the per-batch scores
  and the control counts; `detailed` dropped.
- `parameter_names.md`: `conc`, `neg`, `batch`, `batch_combine` (new);
  `threshold` here is a score cutoff (strict `<`).
- NEWS.


## References

- Davis, N. M., Proctor, D. M., Holmes, S. P., Relman, D. A., & Callahan,
  B. J. (2018). Simple statistical identification and removal of
  contaminant sequences in marker-gene and metagenomics data.
  *Microbiome*, 6, 226. doi:10.1186/s40168-018-0605-2 (3,160 OpenAlex
  citations 2022–26, 478 in 2022 to 751 in 2026; 1,617 Europe PMC
  microbiome mentions; retrieved 2026-11-20).
- Fierer, N., et al. (2025). Guidelines for preventing and reporting
  contamination in low-biomass microbiome studies. *Nature Microbiology*,
  10, 1570–1580. doi:10.1038/s41564-025-02035-2
- Salter, S. J., et al. (2014). Reagent and laboratory contamination can
  critically impact sequence-based microbiome analyses. *BMC Biology*, 12,
  87. doi:10.1186/s12915-014-0087-z
- Eisenhofer, R., et al. (2019). Contamination in low microbial biomass
  microbiome studies: issues and recommendations. *Trends in
  Microbiology*, 27, 105–117. doi:10.1016/j.tim.2018.11.003


## Decisions made while implementing (item 40, 2026-11-24)

Asked of the user at the start (all four answered with the
recommendation):

1. **Names: `decontam()` and `decontam_not()`.**
2. **`conc`**, not `load`.
3. **The frequency plot is a separate function,
   `decontam_frequency_plot(counts, features, conc, neg = NULL, data =
   NULL, margin = 1L, norm = 'percent', models = TRUE, log = TRUE, legend =
   'topright', ...)`,** so `plot()` keeps the no-recompute rule. It takes
   the data, draws one panel per feature (`par(mfrow = rev(n2mfrow(k)))`,
   restored on exit), fits the two lines over the true samples with the
   feature present (controls and `NA` `neg` left out, as the test), and
   returns the plotted values invisibly as a data frame. `features` selects
   by name, index or logical vector (`feature_set()`, as `numerator`).
4. **pkgdown: a new "Quality Control" section** (decontam,
   plot.ecodive_decontam, decontam_frequency_plot), after Source Tracking.

Made while implementing:

5. **Versions.** decontam 1.32.0 (Bioconductor 3.23) is pure R with
   `Depends: R (>= 3.4.1)`, so its source tarball installs on this R 4.3.3;
   the port and the committed fixtures come from it. 1.22.0 (Bioconductor
   3.18) was installed too: its `R/decontam.R` is **byte-identical** to
   1.32.0's, and `data-raw/expected/decontam.r` run with it writes
   byte-identical files (only `version.txt` differs). So the computations
   don't differ between the versions.
6. **Copyright holder:** neither the tarball nor the GitHub repository
   (read through jsDelivr) has a LICENSE file or a copyright line; the
   notices name the DESCRIPTION's authors, Benjamin Callahan and Nicole
   Marie Davis (`Authors@R`: `ctb`, `cph`). Felix G.M. Ernst (`ctb`, the
   S4 dispatch) isn't credited: none of that code is ported. The full
   Artistic-2.0 text is in the headers of `src/decontam.c` and
   `R/decontam.r` and in the `inst/COPYRIGHTS` stanza, with the "how it
   differs" paragraph in each.
7. **The prevalence score's formulas, confirmed from outputs** on every
   2 × 2 table with margins up to 30 × 30 (241,425 tables with a score):
   `prop.test()` warns **exactly** when an expected count `n * p` or
   `n * (1 - p)` (pooled `p`) is below 5; the one-sided statistic with the
   continuity correction `min(0.5, |p1 - p2| / (1/n1 + 1/n2))` matches to
   6.7e-15 (relative), and the Fisher mid-p is bit-identical. The mid-p is
   computed with decontam's own arithmetic, `pg - (pg + pl - 1) / 2` from
   the two `phyper()` tails (not `pg - dhyper()/2`), to keep its rounding.
   The tests repeat this on tables up to 10 × 10 through `decontam()`.
8. **`NA` in `neg`** leaves those samples out of both tests, including the
   prevalence test's "present in at least two samples" count (decontam
   counts them there: a feature in one control and one `NA` sample gets a
   score from decontam and `NA` here). They still count in `freq` and
   `prev`, as decontam's. Confirmed against decontam: its frequency test
   returns `NA` for every feature when `neg` has an `NA`.
9. **Per-sample codes in C instead of subsetting.** Empty samples (and, in
   `batch`, their batch code -1) are left out by a per-sample code, so the
   table is never copied or densified; decontam's warning text is kept
   ("Removed n samples with zero total counts (or frequency).").
10. **Batches.** Levels come from `factor()` of the kept samples (decontam's
    order, so `p_freq_batch` rows match its `p.freqs`), **unused levels
    dropped** (decontam's `table()` counts them and stops); an `NA` batch is
    an error; a batch with fewer than 2 samples is an error, 4 or fewer a
    warning, naming the batches. A batch with fewer than two true samples
    gets `NA` frequency scores, where decontam stops (`apply()` on a
    dropped dimension: "dim(X) must have a positive length").
11. **Kept quirks:** `batch_combine = 'product'` gives 1 when no batch has
    a score, `'fisher'` counts a missing score as 0.5, and both apply with
    a single batch too (as decontam); `'minimum'` and `'combined'`
    propagate `NA` (`pmin()` without `na.rm`).
12. **`neg` must be logical** (0/1 numbers are an error), and `conc` numeric
    and finite > 0 for every kept sample, controls included, as decontam
    checks (`NA` too). `conc` is only checked when the frequency test runs.
13. **`threshold`:** `[0, 1]`; for `'either'` / `'both'` one value is used
    for both without decontam's message, more than two is an error; the
    attribute is then length 2.
14. **Without `neg`:** `prev_neg`, `prev_pos` are `NA`, `n_neg = 0` and
    `n_pos` is every kept sample. `decontam_not()` has no `batch` (nor does
    `isNotContaminant()`), so no `batch_combine` attribute, and no `p_freq`
    column (decontam's is `1 - NA`).
15. **Print header:** "decontam (<test>)" and "decontam_not (prevalence)";
    no `display` (the input is a count table). Plot titles: "Contaminant
    scores (decontam, <test>)", "Non-contaminant scores (decontam_not,
    prevalence)", "Prevalence in controls and samples (...)". With
    `'either'` / `'both'` the scores plot overlays the frequency and
    prevalence histograms, each with its threshold. Plots take `legend`.
16. **Threads:** one task per feature (`f = thread_i; f += n_threads`),
    each writing only its own column of the batch × feature outputs, with
    sums over its samples in order: identical for any `cpus` (tested 1, 2,
    3, 4, 7 and on a 60 × 3000 table). `threads_for()` at ~30 ns per stored
    value.
17. **Fixtures:** simulated tables only (40 × 60 with 8 controls and 8
    contaminants; a 24 × 40 low-biomass table; 30 × 30 with three batches,
    one without controls; a 16 × 12 edge table with an empty sample and
    features in none, one and every sample). Every method × batch
    combination × with/without batches, one- and two-value thresholds, one
    input only (`'auto'`), `normalize = FALSE`; `isNotContaminant()` at two
    thresholds. Hex floats in `data-raw/expected/decontam/`,
    `tests/tinytests/expected/decontam.rds` (29 KB).


## Outcome (from the implementing session's final report)

**Item 40, done 2026-11-24.** Commits `5f01091` (work), `e23e787` (merge of
item 38's commits from `origin/main`: `NAMESPACE`, `R/print.r` and
`plans/parameter_names.md` conflicted, both sides kept, all additions),
and the tick and this outcome in the commit after it.

- Tests (`-O2` install, `LANG=C.UTF-8`): baseline at cbf9c7e 22,864 at
  home and 388 in CRAN mode; before the merge 26,343 and 394
  (`test-decontam.r` has 3,473); on the merged tree **26,600 at home, 396
  in CRAN mode**, all passing. `R CMD check --no-manual
  --ignore-vignettes`: 1 NOTE (installed size; 5.3 MB at the baseline
  too).
- decontam 1.32.0 (fixtures; 1.22.0 byte-identical): every method, batch
  combination, with and without batches, one- and two-value thresholds,
  `normalize = FALSE`, and `isNotContaminant()`: largest relative
  difference **1.4e-14** (scores, `freq`), tested at 1e-12; `prev` and
  every call identical. Prevalence formulas against `prop.test()` /
  `fisher.test()` on all 241,425 tables with margins up to 30 × 30:
  warning rule identical, 6.7e-15 (chi-squared branch), 0 (Fisher mid-p).
- Identical for `cpus` 1, 2, 3, 4, 7; dense, `dgCMatrix`, `dgTMatrix`,
  `dgeMatrix`, slam and `margin = 2` identical.
- Speed: 500 samples × 20,000 features (5% non-zero), 4 batches, both
  tests: 0.11 s on 1 thread, 0.09 s on 4 (the C pass is a small share; the
  rest is R's `pf()` and the prevalence formulas).
- Open: none known. decontam's `plot_condition()` (DNA concentration
  against an experimental condition) wasn't ported; it is a plain boxplot
  or scatter of metadata.
