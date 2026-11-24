# `lefse()`: LDA effect size (LEfSe) biomarker discovery

**Priority:** 41. **Tier:** port (LEfSe 1.1.2, SegataLab/lefse, **MIT**:
`license.txt` "Copyright (c) 2012, Nicola Segata and Curtis Huttenhower"
with the MIT text, `setup.py` `version='1.1.2'`; verified 2026-11-20 on
GitHub master). LEfSe's Python calls R through rpy2 for three steps, and
those R functions are **output-only**: `kruskal.test()` (R, GPL),
`coin::wilcox_test()` (coin, GPL-2) and `MASS::lda()` (MASS, GPL-2 | 3).
ecodive already reproduces the first and third (`ancom(sig_test =
'kruskal')`, `lda()` from scikit-learn's BSD-3 code), and the second is the
normal approximation of the rank-sum test that `da_wilcoxon()` already
computes. Cross-check: Bioconductor **lefser 1.22.0** (Bioc 3.23,
**Artistic-2.0**, `DESCRIPTION` verified 2026-11-20; its R code may be read
and ported, with the Artistic notices, but nothing from it is needed beyond
documenting how it differs). Read [README.md](README.md) first, and
[da_test.md](da_test.md).


## What it is

LEfSe (Segata et al. 2011) finds features that "explain" a difference
between classes in three steps:

1. a **Kruskal-Wallis** test of each feature across the classes (alpha
   0.05, no adjustment by default);
2. for features passing it, **pairwise Wilcoxon rank-sum tests** between
   the subclasses of different classes, all of which must be significant
   with the same direction ("biological consistency");
3. an **LDA effect size**: a linear discriminant analysis of the classes on
   all surviving features, repeated on 30 bootstrap subsamples; each
   feature's score is \eqn{\log_{10}(1 + m)}, with \eqn{m} the mean over
   bootstraps of half the sum of its between-class mean difference and its
   scaled discriminant coefficient. Features with |score| > 2 are reported.

Demand (retrieved 2026-11-20): Segata et al. (2011) 10,514 OpenAlex
citations 2022–26 and **19,755 Europe PMC microbiome full-text mentions
2022–26**, 8.9× MaAsLin2's 2,232 and close to Bray-Curtis (28,581): LEfSe
is the most used differential-abundance method in microbiome papers by a
wide margin.

**Its statistics are weak, and the docs must say so.** Nearing et al.
(2022, 38 datasets) found LEfSe among the methods calling the most features
significant, with the highest false discovery rates in their random-group
tests on filtered data (with edgeR), among the least consistent results
across datasets, and recommended avoiding it "without p-value correction".
**Correction to the brief:** Wirbel et al. (2024) and Pelto et al. (2025)
did **not** include LEfSe (neither paper's full text mentions it, checked
2026-11-20); they support the docs' advice to prefer the plain tests
(`da_wilcoxon()`, `da_ttest()`, `da_lm()`), not a LEfSe result. The docs
point users to `da_sensitivity()` and to those tests, as
[dca.md](dca.md) points to NMDS: ecodive offers LEfSe because users and
reviewers expect it, not because it is recommended.


## Names

- `lefse()`: the method's own lowercase name (rule 1); LEfSe is what users
  type and cite. No package exports `lefse()` in R (lefser has `lefser()`;
  LEfSe is a Python command, `lefse_run.py`). "Same names in other
  packages": none.
- The **LDA score** is not ecodive's `lda()`: README already warns that
  "LEfSe's 'LDA score' is an effect size from this method". `lda()`'s Rd
  sentence stays; `lefse()`'s Rd links to it.
- **`DA_TESTS` row** ([da_test.md](da_test.md)): name `LEfSe`, id `lefse`,
  design `groups`, effect `lda_score`, base `NA` (a log10 of a mean
  difference, not a log ratio; `da_sensitivity()` then reports
  `log_base = NA`), `zero_free_ok = TRUE` (it uses scaled abundances, and a
  zero-replaced table is a valid input). Aliases (letters only): `lefser`.
  **Not** `lineardiscriminantanalysiseffectsize`: an alias starting
  `linear` would break `match_da_test('linear')` (`da_lm`), as linda's
  plan found; and not `lda`, which is `COUNT_ORDS`'s. Check with README's
  alias snippet.
- It fits `da_test()`: one categorical grouping (`groups`), a `reject`
  column, `alpha`, a `method` attribute. Its own extra argument
  `subgroups` is not a prefix of any wrapper argument, and none of its
  arguments is a prefix of `counts`, `method`, `groups`, `data` or
  `formula` (`test-wrapper_args.r` checks).


## Design

### Signature

```r
lefse(counts, groups, data = NULL, subgroups = NULL, margin = 1L,
      total = 1e6, alpha = 0.05, wilcoxon = TRUE, wilcoxon_alpha = 0.05,
      lda_threshold = 2, n_boots = 30L, boot_fraction = 0.67,
      min_subgroup = 10L, one_against_one = FALSE, same_subgroups = FALSE,
      p_adjust = 'none', seed = 0L, cpus = n_cpus())
```

| ecodive | LEfSe (`lefse_run.py` / `lefse_format_input.py`) | Default |
|:--------|:--------------------------------------------------|:--------|
| `groups` | class (`-c`) | |
| `subgroups` | subclass (`-s`): a column of `data` or one value per sample | none: one subgroup per class |
| `total` | `-o` normalization value: each sample scaled to sum to `total`; `NULL` uses the values as given (LEfSe's default `-1`) | `1e6` (**ask**; see below) |
| `alpha` | `-a` (Kruskal-Wallis) | 0.05 |
| `wilcoxon` | `--wilc` | `TRUE` |
| `wilcoxon_alpha` | `-w` | 0.05 |
| `lda_threshold` | `-l` | 2 |
| `n_boots` | `-b` | 30 |
| `boot_fraction` | `-f` | 0.67 |
| `min_subgroup` | `--min_c` | 10 |
| `one_against_one` | `-y` (`multiclass_strat`) | `FALSE` (one against all) |
| `same_subgroups` | `-e` (`only_same_subcl`) | `FALSE` |
| `p_adjust` | none in LEfSe; lefser's `method` | `'none'` (**ask**) |

- `total = 1e6` is what the Galaxy/Huttenhower-lab instructions and lefser
  (`relativeAb()`, counts per million) use, and the LDA threshold of 2 is
  calibrated on that scale (the score is a log10 of a difference in these
  units). LEfSe's own command-line default is no scaling; say so.
  LEfSe sums only the top-level features of a hierarchical table
  (`numerical_values()`); ecodive's tables are flat, so the sample total is
  the sum of all features. Name `total`: no glossary word fits (`norm`
  takes names, not a target sum). **Ask** if the user prefers another.
- `p_adjust`: over features, on the Kruskal-Wallis p-values, before the
  `alpha` cut (lefser's behaviour; LEfSe has none). ecodive's other DA
  methods default to `'holm'`; `'none'` here reproduces what papers report,
  and Nearing et al. single out the missing correction. **Ask the user**
  which default; the docs recommend `'BH'` either way.
- Not carried over: `-s` (LEfSe's "correction" *widens* the Wilcoxon alpha
  to \eqn{1 - (1 - \alpha)^{k}} or \eqn{k\alpha}; document and skip), `-c`
  (Curtis's approach, marked BETA), `-r svm` (commented out in LEfSe
  1.1.2), `subject` (`-u`, used only for sorting), hierarchical clade
  expansion (`add_missing_levels()`: summing children into parents; users
  can supply clade tables), the cladogram.

### Algorithm (read from LEfSe 1.1.2 `lefse/lefse.py`, `lefse_run.py`, `lefse_format_input.py`)

Let \eqn{X} be the table scaled per sample to `total` (zero-total samples
stay 0, as LEfSe).

1. **Class means** (`get_class_means()`): each feature's mean in each
   class, over all samples. `enriched` = the class with the largest mean
   (first on ties); `log_max_mean` = \eqn{\log_{10}\max(1, \max_k
   \bar{x}_k)} (the second column of LEfSe's output).
2. **Kruskal-Wallis** (`test_kw_r()`): `kruskal.test(y ~ class)`,
   average ranks, tie correction, \eqn{\chi^2_{K-1}} p-value. Pass:
   `p_adjusted < alpha` (strict, as LEfSe). Reuse ancom.c's `kruskal_h()`
   (matches SciPy and R; refactor it out of `ancom.c` into a shared rank
   helper, e.g. `src/rank_tests.c`, used by both).
3. **Wilcoxon step** (`test_rep_wilcoxon_r()`), only for features passing 2
   and when `wilcoxon = TRUE`. Without `subgroups`, every class is one
   subgroup, so this is a two-class (or per-pair) rank-sum test. For each
   pair of classes \eqn{(a, b)} (`a < b` in sorted order), compare every
   subgroup of \eqn{a} with every subgroup of \eqn{b} (with
   `same_subgroups`, only subgroups with the same name; non-matching pairs
   count as passed):
   - if either subgroup has fewer than `min_subgroup` samples, only the
     **direction of the medians** is used (`med_comp`): the first
     comparison sets the direction; later ones must agree and differ
     (`sx != sy`);
   - otherwise the two-sided p-value of `coin::wilcox_test()`
     (asymptotic: the normal approximation with average ranks and the tie
     variance, **no continuity correction** ✱) must be `< 2 *
     wilcoxon_alpha` (LEfSe's doubling, kept), with the same median
     direction as the first comparison and `sx != sy`;
   - two subgroups whose values are all equal and identical fail.
   The pair differs when every comparison passes. Two classes: the feature
   passes if the pair differs. More classes: `one_against_one = TRUE`
   needs every pair to differ; `FALSE` needs one class to differ from all
   others. Port this control flow line by line; it has quirks (the first
   comparison under `med_comp` passes even with equal medians) that the
   docs list.
   ✱ Confirm from outputs that `coin::wilcox_test()`'s asymptotic p equals
   `wilcox.test(exact = FALSE, correct = FALSE)`'s; ecodive's normal
   branch in `src/da_tests.c` then serves (add a no-correction flag).
4. **Jitter** (start of `test_lda_r()`), for the LDA only: for each
   surviving feature and class, if the feature has at most
   \eqn{\max(n_c/2, 4)} distinct values in that class, replace each of
   its values \eqn{v} there by \eqn{|v + \varepsilon|},
   \eqn{\varepsilon \sim N(0, \max(0.05 v, 0.01)^2)} (sd from each value,
   not the vector: lefser's commented-out copy differs).
5. **Bootstraps** (`test_lda_r()`), \eqn{p} surviving features:
   - draw \eqn{\lfloor n \cdot \text{boot\_fraction}\rfloor} samples **with
     replacement**; redraw (up to 1000 times, then keep the last) while
     some class is missing, has fewer than \eqn{m = \max(1,
     \lfloor 0.5\, f^2 \min_k n_k\rfloor)} samples, or some feature has
     \eqn{\le m} distinct values in some class (\eqn{\le 1} when
     \eqn{m = 1}) (`contast_within_classes_or_few_per_class()`);
   - fit `MASS::lda(class ~ ., tol = 1e-10)` on the subsample; \eqn{w} =
     the first discriminant's coefficients (`scaling[, 1]`), \eqn{u =
     w/\|w\|}; effect \eqn{e = |\bar{s}_a - \bar{s}_b|} for the scores
     \eqn{s = X u} of each class pair; per feature \eqn{j}:
     \eqn{(|\bar{x}_{aj} - \bar{x}_{bj}| + |u_j| e)/2}, with the class means
     of the subsample (`z$means`), 0 for a class absent from it, and
     \eqn{|u_j| e = 0} when `NaN`. Multi-class: the same first discriminant
     for every pair.
6. **Score**: per feature and pair, the mean over bootstraps; \eqn{m_j} =
   the largest over pairs; `lda_score` \eqn{= \log_{10}(1 + m_j)} (always
   \eqn{\ge 0}: LEfSe's `copysign` sees a non-negative \eqn{m}). Reported
   when `lda_score > lda_threshold` (strict).
7. `reject` = passed 2, passed 3 (or `wilcoxon = FALSE`) and
   `lda_score > lda_threshold`.

### The LDA: reuse `lda()`, with MASS's rank-deficient path

LEfSe's LDA is `MASS::lda()`, which `lda()` already matches (scores,
`scaling`, means) for full-rank tables. Two gaps:

- **More features than samples.** LEfSe fits all surviving features at once
  on about \eqn{0.67 n} samples, so \eqn{p \ge n - g} is common. MASS keeps
  the within-group directions whose scaled singular value exceeds `tol`
  (documented: it "will reject variables and linear combinations of
  unit-variance variables whose variance is less than tol^2") and warns
  "variables are collinear"; scikit-learn's svd solver does the same with
  its own `tol`. ecodive's `lda()` stops instead (its "Singular tables"
  section). Add an internal path to `src/lda.c`: a wide-matrix SVD of the
  sphered residuals (QR of \eqn{Z^\top}, then the existing one-sided
  Jacobi), a rank cut at `tol`, and the scaling \eqn{S^{-1} V_r D_r^{-1}}.
  `lda()` itself keeps its error. ✱ Confirm MASS's exact rank rule (the
  scale `tol` applies to) from crafted outputs. With collinear features the
  direction of \eqn{u} depends on rounding (MASS's warning exists for that
  reason), so compare such cases at a loose tolerance only.
- **Threads.** `lda_fit()` uses `R_alloc()` and its own `run_parallel()`,
  so it can't run inside a worker. Write a worker-safe serial variant
  (caller-supplied buffers), and thread **over bootstraps**: each bootstrap
  draws from its own PCG stream (`seed`, bootstrap), sums are stored per
  bootstrap and averaged in bootstrap order, so results are identical for
  any `cpus`. Jitter draws use one stream per feature (`rng_norm()` in
  `src/rng_dist.c`). A feature constant within groups in one subsample (MASS
  errors, which would crash LEfSe) is left out of that bootstrap's fit, its
  contribution taken as its mean difference alone; listed in
  "Differences".

LEfSe's "LDA score" differs from the `eig` / `prop` of `lda()`: it is not
a discriminant ratio but an average of a raw mean difference (in units of
`total`) and a unit-norm coefficient times the class separation along the
axis. The first term dominates for abundant features, which is why the
score tracks log abundance (Nearing et al.'s and others' criticism).

### Result

`c('ecodive_lefse', 'data.frame')`, one row per feature (every feature,
not only reported ones): `feature, enriched, log_max_mean, kw_stat,
p_value` (Kruskal-Wallis), `p_adjusted, wilcoxon` (logical: passed step 3;
`NA` when step 2 failed or `wilcoxon = FALSE`), `lda_score` (`NA` when
not computed), `reject`. Attributes: `method = 'lefse'`, `groups` (the
classes, in order), `subgroups` (`NULL` or the subgroup names), `alpha`,
`wilcoxon_alpha`, `lda_threshold`, `n_boots`, `boot_fraction`, `total`,
`p_adjust`, `means` (class × feature means of the scaled table, for the
plot and for `enriched`). Printing: `print_table_result`, header "LEfSe:
a vs b" (add to `print_df_label()`).


## Reference implementations

- **Port from:** LEfSe 1.1.2 (MIT), functions named above. File header of
  `src/lefse.c` and `R/lefse.r`: the MIT text with "Copyright (c) 2012,
  Nicola Segata and Curtis Huttenhower"; `inst/COPYRIGHTS` stanza;
  `Authors@R`: Nicola Segata and Curtis Huttenhower (`ctb`, `cph`).
- **Output-only:** R's `kruskal.test()`, coin's `wilcox_test()`, MASS's
  `lda()`; never open their source. Their behaviour comes from their
  documentation and from outputs on crafted inputs.
- **Fixtures** (`data-raw/expected/lefse.py` + `lefse.r`):
  - Install LEfSe 1.1.2 (pip, with rpy2, R's coin and MASS) in a scratch
    environment. The `.py` generator imports `lefse.lefse` and wraps its
    `lrand` (Python's `random`) to **record** the jittered values and every
    bootstrap's sample indices, then calls `test_kw_r()`,
    `test_rep_wilcoxon_r()` and `test_lda_r()` directly (all features'
    LDA values, not only those above the threshold). Text outputs under
    `data-raw/expected/lefse/`.
  - Cases: two classes, 40 samples, no subgroups (Wilcoxon between the
    classes); two classes with two subgroups each, one subgroup under
    `min_subgroup` (median rule); three classes, both `one_against_one`
    settings; `same_subgroups = TRUE`; a feature with few distinct values
    (jitter); a case where features outnumber the bootstrap size
    (rank-deficient LDA).
  - `lefse.r`: lefser 1.22.0 on the two-class case (for "Differences"),
    and `wilcox.test(exact = FALSE, correct = FALSE)` vs
    `coin::wilcox_test()` on tied data (the ✱ in step 3).
  - An internal `lefse_details(jitter = , boots = )` injects the recorded
    draws so that the LDA scores compare exactly; with PCG draws they
    differ by a different seed only (README: that is fine).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- `plot.ecodive_lefse(x, n_features = Inf, ...)`: LEfSe's bar chart
  (`lefse_plot_res.py`, MIT): one horizontal bar per reported feature
  (`reject`), length `lda_score`, coloured by `enriched`; two classes:
  the first class's bars to the left (negative), sorted by score within
  class, as LEfSe draws them; more classes: grouped by class. x label
  "LDA score (log 10)". Base graphics, the palette of `plot_palette()`.
- Title: "LDA scores of LEfSe" (the `DA_TESTS` name via
  `da_method_name()`); subtitle the classes and the number reported. The
  input is a count table, so no `display`.
- Attributes read: `lda_threshold` (a dashed line), `groups` (legend
  order). No volcano: LEfSe has no effect with a p-value of the same
  feature set worth plotting together. Later step: a per-feature
  abundance strip chart (`lefse_plot_features.py`) from `means`.


## Tests

- Kruskal-Wallis statistics and p-values equal `kruskal.test()` (fixtures)
  at 1e-12; the Wilcoxon decisions and p-values equal LEfSe's; the set of
  features passing each step identical.
- With injected jitter and bootstraps: `lda_score` at 1e-10 for
  full-rank cases, 1e-6 for rank-deficient ones; `reject` identical.
- With PCG draws: `reject` agrees with LEfSe on clear-cut features (scores
  far from the threshold), and the scores are close (relative 0.05).
- `identical()` for `cpus` 1, 2, 4; dense = sparse; `margin = 2`.
- Invariant to scaling a sample by a constant (with `total`).
- Edge cases: a class with one sample; a feature absent everywhere; all
  features failing step 2 (empty `reject`, scores `NA`, no error); a
  single surviving feature; `n_boots = 1`; three or more classes; NA in
  `groups` or `subgroups`; invalid `boot_fraction`.
- Smoke line.


## Docs

Rd page in *Differential abundance*:

- The three steps and the score; that it compares relative abundances
  (scaled to `total`), not absolute ones.
- **"Use with care"**: Nearing et al.'s findings; that LEfSe reports no
  adjusted p-values by default; that the score grows with abundance; that
  bootstrapping makes scores (not the Kruskal-Wallis and Wilcoxon steps)
  depend on `seed`. Point to [da_wilcoxon()] / [da_lm()] (the benchmarks'
  choices) and to [da_sensitivity()] to see whether calls survive other
  settings.
- **"Differences from LEfSe"**: PCG instead of Python's Mersenne Twister
  (scores differ as a different seed would); flat tables (no clade
  expansion); `p_adjust`; `-s`, `-c`, `-r svm` and subjects not offered;
  a feature constant within a class in one subsample is left out of that
  fit instead of stopping.
- **"Differences from lefser"** (lefser 1.22.0, Artistic-2.0, read):
  lefser fits one LDA on all samples (no bootstrap, no jitter), uses
  `MASS::lda()`'s default `tol` (1e-4), gives *signed* scores from
  `(means[2] - means[1] + |coefficient|)/2` (so a feature richer in the
  first class gets a smaller magnitude than LEfSe's), allows two classes
  only, skips the Wilcoxon step without subclasses, compares subclasses
  without LEfSe's `min_c` median rule (requiring the z statistics' signs to
  agree instead), keeps `|score| >= threshold` (LEfSe: `>`), and offers
  `p.adjust` on both tests. ecodive follows LEfSe, which papers cite.
- `parameter_names.md`: `subgroups`, `total`, `wilcoxon`,
  `wilcoxon_alpha`, `lda_threshold`, `n_boots`, `boot_fraction`,
  `min_subgroup`, `one_against_one`, `same_subgroups` (new).
- NEWS; pkgdown "Differential abundance".


## References

- Segata, N., Izard, J., Waldron, L., Gevers, D., Miropolsky, L.,
  Garrett, W. S., & Huttenhower, C. (2011). Metagenomic biomarker discovery
  and explanation. *Genome Biology*, 12, R60.
  doi:10.1186/gb-2011-12-6-r60 (10,514 OpenAlex citations 2022–26;
  19,755 Europe PMC microbiome mentions 2022–26; retrieved 2026-11-20).
- Khleborodova, A., Gamboa-Tuz, S. D., Ramos, M., Segata, N., Waldron, L.,
  & Oh, S. (2024). lefser: implementation of metagenomic biomarker
  discovery tool, LEfSe, in R. *Bioinformatics*, 40(12), btae707.
  doi:10.1093/bioinformatics/btae707
- Nearing, J. T., et al. (2022). Microbiome differential abundance methods
  produce different results across 38 datasets. *Nature Communications*,
  13, 342. doi:10.1038/s41467-022-28034-z
- Wirbel, J., et al. (2024). *Genome Biology*, 25, 247.
  doi:10.1186/s13059-024-03390-9 (does not test LEfSe).
- Pelto, J., et al. (2025). *Briefings in Bioinformatics*, 26(2),
  bbaf130. doi:10.1093/bib/bbaf130 (does not test LEfSe).
- Fisher, R. A. (1936). The use of multiple measurements in taxonomic
  problems. *Annals of Eugenics*, 7, 179–188.
  doi:10.1111/j.1469-1809.1936.tb02137.x


## Decisions made while implementing (2026-11-24)

Asked of the user at the start:

- **`total = 1e6`, named `total`** (the recommendation): counts per
  million, as Galaxy and lefser; `NULL` uses the values as given (LEfSe's
  command-line default). The `total` attribute is `NA` for none.
- **`p_adjust = 'holm'`** (the user's choice, over the recommended
  `'none'`): consistent with ecodive's other DA methods, so the default
  result is *not* LEfSe's; `p_adjust = 'none'` reproduces it, and the docs
  recommend `'BH'`. Adjusted over the features with a p-value (`NaN` ones
  left out), before the strict `< alpha` cut.

Found from outputs (GPL code not read):

- **coin's asymptotic `wilcox_test()`** p-value equals `wilcox.test(exact =
  FALSE, correct = FALSE)` (300 random tied cases, 4.9e-15 relative): the
  normal approximation with the tie-corrected variance and no continuity
  correction, computed in R from C's rank sum and \eqn{\sum t^3 - t}.
- **`MASS::lda()`'s rank rule:** each variable is scaled by its pooled
  within-group sd with divisor **n − 1** (`var()` of the residuals), the
  SVD of the scaled residuals over \eqn{\sqrt{n - g}} keeps singular values
  **> tol** (crafted collinear pairs put the cut between the n − 1 and n − g
  predictions, and only n − 1 fits), and it stops ("appears to be constant
  within groups") when that sd is **< tol**. For p > n, my formula
  reproduced MASS's `scaling` to 1e-15 (ratio exactly 1).

Made while implementing:

- **C (`src/lefse.c`)**, four entry points: `C_lefse_kw` (class means and
  H, threaded over features, from a sparse table), `C_lefse_columns` (the
  dense rescaled columns of the features that reach the LDA),
  `C_lefse_wilcoxon` (per feature and pair of subgroups: rank sum, ties,
  numpy-style medians, and the both-constant-and-equal flag) and
  `C_lefse_lda` (jitter, then bootstraps threaded one per task).
- **Shared rank kernel:** ancom.c's sort and Kruskal-Wallis H moved to
  `src/rank_tests.c` (`rank_sort_vg()`, `rank_kruskal_h()`, plus
  `rank_sum_ties()`), `vg_t` to `ecodive.h`; ancom's code is otherwise
  unchanged (its 1,948 tests pass, bit-identical H).
- **The discriminant** is a worker-safe serial routine in `src/lefse.c`,
  not a change to `lda()` (which keeps its singular-table error): features
  with sd < 1e-10 left out (coefficient 0); for q ≤ m features a serial
  Householder QR then lda.c's one-sided Jacobi SVD (now exported as
  `lda_jacobi_svd()`); for q > m Jacobi directly on Z' (its columns become
  U D, giving the right singular vectors as Z'u / d). Rank cut d > 1e-10.
  Only the first discriminant is formed: Jacobi on B' (r × g), whose first
  column / its norm is B's first right singular vector. Degenerate
  subsamples (a class missing so g < 2, m − g < 1, rank 0, identical group
  means) give u = 0, so each pair contributes half its mean difference
  (LEfSe would crash or give NaN coefficients, which it zeroes).
- **Rescaling exactly as `numerical_values()`:** `value * (total / sum)`,
  the sum in feature order, zero-total samples 0, and its rounding of a
  feature with relative sd < 1e-10 to 6 decimals (`nearbyint`, half-even
  like Python's `round()`). The rescaled table is **bit-identical** to
  LEfSe's in every fixture.
- **Wilcoxon logic** ported line by line, vectorised over features: a pair
  of subgroups that both hold one and the same value fails; the first
  comparison passes under `min_subgroup` even with equal medians (direction
  "not less"); `p < 2 * wilcoxon_alpha`; one-against-all needs a class
  whose k − 1 pairs all differ.
- **`same_subgroups` matches subgroups by name** (documented difference):
  LEfSe compares the ends of the (renamed) subclass names after dropping
  `len(class)` characters, which also pairs some different names that occur
  in one class only. Without `subgroups` it has no effect (as in LEfSe,
  whose `<class>_subcl` names always match).
- **Orders:** classes and each class's subgroups by factor levels, else
  `sort(method = 'radix')` (LEfSe sorts by code point). `enriched` is the
  first class with the largest mean.
- **Random streams:** bootstrap b uses PCG stream (seed, b), with
  `pcg32_bounded()` draws; the jitter of feature f (its table index) uses
  stream (seed, 2^32 + f) and `rng_norm()`. The redraw check
  (`contast_within_classes_or_few_per_class()`) is ported exactly (up to
  1000 tries, then the last). Bootstrap sums are kept per bootstrap and
  averaged in order: identical for any `cpus` (tested 1, 2, 3, 4, 8, and
  the rank-deficient case on 1 and 4).
- **`run_parallel()` stays serial under 100 tasks**, so the bootstrap call
  passes `max(100, n_boots)` once `threads_for()` has judged the work
  (found when 4 threads gave no speed-up).
- **Result:** columns as the plan; `wilcoxon` is `NA` when not run;
  doubles' missing values are `NaN`; `enriched` and `log_max_mean` for
  every feature. Attributes as the plan plus `log_base = NA_real_` (the
  item-27 convention; the registry's `base` is `NA`). Print header
  "LEfSe: a vs b".
- **Checks:** `boot_fraction` in (0, 1], and floor(n × boot_fraction) must
  exceed the number of classes; `min_subgroup` ≥ 0; `lda_threshold` any
  finite number (negative: every feature that passed both tests is
  reported, as LEfSe's `-l` < 0 does, though here the scores are still
  computed); NA in `groups` or `subgroups` is an error.
- **`lefse_details()`** (internal) takes `jitter` and `boots` (1-based) to
  replace the draws, and returns the jittered values, the subsamples, each
  pair's averaged effects and every Wilcoxon p-value.
- **Plot:** LEfSe's layout (two classes: the first to the left, sorted by
  signed score; more: grouped by class), `plot_palette()` colours, dashed
  lines at ±`lda_threshold`, legend bottom right, "No features reported"
  when empty. Title "LDA scores of LEfSe" (from `da_method_name()`),
  subtitle "a vs b: k of D features reported".
- **Registry:** `LEfSe, lefse, groups, lda_score, NA, TRUE, lefser`. The
  alias collides with nothing; `match_da_test('l')` now lists three
  methods (it already listed two).
- **`da_sensitivity()`:** zero replacements only (no `pseudocount`
  argument, so `pseudocounts` errors with the usual message); added to its
  per-method `identical()` tests.
- **References used:** LEfSe 1.1.2 from PyPI (identical to jsDelivr's
  GitHub copy) with **rpy2 3.5.17** (rpy2 3.6 needs R ≥ 4.4). rpy2 3.5
  returns `None` for R's invisible results, which `test_lda_r()` reads
  (`mm <- z$means`), so the generator sets `robjects.r._invisible = False`;
  `lefse_format_input.py` also needs a subject row (`-u`) or it stops with
  `KeyError: 'subject'`. coin 1.4-3, MASS 7.3-60.0.1, R 4.3.3.
  **lefser 1.12.1** (Bioconductor 3.18, the newest that installs on R
  4.3.3), not 1.22.0: its code differs from the plan's description of
  1.22.0 (1.12.1 *does* bootstrap, 30 subsamples of 2/3, and jitters, with
  one sd per class column from its largest value), so "Differences from
  lefser" describes 1.12.1.
- **Fixtures** (`data-raw/expected/lefse.r` + `.py`, 9 cases, 56 KB rds):
  two classes (with an all-zero feature and a jittered one); two classes
  with two subclasses each, one under `min_c`; same subclass names with
  `-e 1`; three classes with `-y 0`, `-y 1` and `--wilc 0`; a wide case
  (40 features scored, 13 samples a subsample); one bootstrap at f = 0.5
  with `min_c = 5` (every pair tested); and deep samples whose rare
  features score below 2. The generator records each bootstrap's samples
  by wrapping `contast_within_classes_or_few_per_class()`, and the jitter
  from `test_lda_r()`'s in-place changes.
