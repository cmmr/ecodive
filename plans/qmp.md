# `qmp_scale()` and `qmp_rarefy()`: quantitative microbiome profiling

**Priority:** 18. **Tier:** port (`qmp_rarefy()`:
the authors' `QMP.R`, raeslab/QMP, BSD-3, © 2017 Raes Lab) and formula
(`qmp_scale()`), both a few lines over ecodive's own `rarefy()` C code. Read [README.md](README.md)
first, and [compositional.md](compositional.md) for the conventions the
compositional plans share.


## What it is

Sequencing gives relative abundances. Every compositional method
(log ratios, proportionality, ANCOM, ...) works around that, and can only
recover *ratios* between taxa. When a study has also **measured the total
microbial load** of each sample (flow cytometry cell counts, total 16S
qPCR, a spike-in, or DNA yield), the relative profile can be turned into an
estimate of absolute abundance instead. This is the other side of the
compositionality debate: if the denominator was measured, use it.

There are two named ways of doing it, so there are two functions (the
README's rule: one function per named method variant):

| Function | Method | What it does |
|:---------|:-------|:-------------|
| `qmp_scale()` | "Relative abundance × load" (Jian et al. 2020; most qPCR and spike-in studies; often also called QMP) | Divide each sample's counts by its total, multiply by its load. |
| `qmp_rarefy()` | Quantitative microbiome profiling (Vandeputte et al. 2017) | Rarefy each sample to the **same sampling depth relative to its load** (reads per cell), then scale to the load. |

`qmp_rarefy()` is Vandeputte's original: by rarefying each sample to an
equal number of reads per cell, every sample has the same detection limit
in absolute units, which plain scaling doesn't give.


## Names

- `qmp_<modifier>`, so that the pair sorts together (README naming rule 2;
  there is no bare `qmp()`, as there is no bare `unifrac()`).
- Aliases are for the docs only (no wrapper): "QMP", "quantitative
  microbiome profiling", "absolute abundance", "load scaling".
- Don't call either `normalize(norm = 'load')`: `norm` choices take no
  per-sample data and must work inside every function with a `norm`
  argument, and this one can't.


## Design

### Shared arguments

- `counts`, `margin`, `cpus`: as everywhere.
- **`load` (new glossary entry):** the total microbial load of each sample,
  in any absolute unit (cells per gram, 16S copies per µl, ...).
  - A numeric vector with one value per sample, or a column name of
    `data`, handled exactly as `groups` is in the tests: match by sample
    name when both have names (`per_sample()` / `match_rows()` in
    `R/design.r`), by position otherwise.
  - Finite and > 0; NA is an error naming the samples (don't silently drop
    samples).
- `data = NULL`: sample metadata, for `load` given as a column name.

### `qmp_scale(counts, load, data = NULL, margin = 1L, cpus = n_cpus())`

- Result: \eqn{x_{ij} / n_i \times L_i}, where \eqn{n_i} is sample i's
  total and \eqn{L_i} its load.
- **Output:** the same class, orientation and dimnames as `counts`,
  exactly as `normalize()` returns them (sparse stays sparse: zeros stay
  zero). Reuse `normalize()`'s output code: factor it into a helper
  rather than copying it.
- C: one pass, as `NORM_PERCENT` in `src/normalize.c`, with a per-sample
  multiplier. It could even be `normalize(counts, 'percent')` followed by a
  per-sample multiplication of the stored values in R; that's fine if it
  keeps sparse input sparse and doesn't copy a dense matrix twice.
- An empty sample (total 0) is an error naming it.

### `qmp_rarefy(counts, load, data = NULL, margin = 1L, seed = 0L, cpus = n_cpus())`

Vandeputte et al. (2017), Methods, "Quantitative microbiome profiling":

1. Each sample's **sampling depth** is \eqn{s_i = n_i / L_i} (reads per
   unit of load).
2. \eqn{s_{min} = \min_i s_i}.
3. Rarefy sample i, without replacement, to
   \eqn{t_i = \lfloor s_{min} L_i \rfloor} reads (\eqn{\le n_i} by
   construction; the sample with the smallest \eqn{s_i} keeps all its
   reads, up to the floor).
4. Scale the rarefied counts to the load:
   \eqn{q_{ij} = r_{ij} / t_i \times L_i}.

- **The authors' script** (`QMP.R`, BSD-3; read): rounds the
  (copy-number corrected) table **up** with `ceiling()` first; computes
  \eqn{s_i} and the targets \eqn{L_i s_{min}} as above; rarefies each
  sample with `phyloseq::rarefy_even_depth(rngseed = 711)` (phyloseq is
  AGPL: don't read it; ecodive's own `rarefy()` replaces it); and returns
  the rarefied counts **divided by their row sum** times \eqn{L_i} (so
  step 4 divides by the achieved total, \eqn{t_i}). Two things to settle
  from outputs: how a non-integer target is turned into a whole number of
  reads (floor is the natural reading; run the script on a table where it
  matters), and whether `qmp_rarefy()` should apply the `ceiling()` too
  (recommended: no, require integer counts and say how to reproduce the
  script with `ceiling()`).
- Counts must be integers (`assert_integer_counts()`, as `rarefy()`).
- **Reuse `src/rarefy.c`.** Generalise `C_rarefy` to take a per-sample
  target (a vector; a length-1 vector keeps today's behaviour) in place of
  the single static `target`. The per-sample PRNG seeding
  (`pcg32_srandom_r(&rng, seed, sam)`) already makes results independent
  of `cpus`.
  - **Also expose it** by letting `rarefy(depth =)` take one depth per
    sample (backward compatible: a single number works as now). `drop`
    and `warn` then compare each sample with its own depth. Note it in
    NEWS. This is the only change to an existing function.
- **Output:** the same class as `counts` (sparse stays sparse, as
  `rarefy()`), holding doubles (step 4 isn't integer).
- Attributes (on the base-matrix result; for sparse classes, document that
  they are attached to the object too): `sampling_depth` (the \eqn{s_i},
  named), `target` (the \eqn{t_i}).

### What the docs must say

- Neither method makes the data less noisy: the load measurement's error
  carries straight into every taxon (Vandeputte et al. 2017; Props et al.
  2017). Measure load in replicate.
- **16S copy-number correction** (dividing each taxon's counts by its 16S
  gene copy number) is a separate, earlier step in Vandeputte's pipeline.
  It isn't part of either function; show it in an example with `sweep()`.
- For a spike-in, the load is usually
  `spike_added * total_reads / spike_reads`; show it, and drop the spike-in
  feature before calling `qmp_*()`.
- The results are **absolute abundances, not compositions**: don't follow
  them with a CLR (which removes the load again). Bray-Curtis on QMP data
  is the usual next step, via `bray()` / `beta_div()` with `norm = 'none'`.
- Cross-link `normalize()`, `rarefy()` and the compositional functions.


## Wrapper

None. Two functions with a per-sample input don't fit `normalize()`'s
`norm` vocabulary (see "Names"), and a wrapper for two functions isn't
worth it.


## Reference implementations

- **`qmp_rarefy()`:** the authors' `QMP.R` (https://github.com/raeslab/QMP,
  BSD-3, © 2017 Raes Lab): port it and credit it in `inst/COPYRIGHTS` and
  `Authors@R` as README.md requires. Its rarefaction calls phyloseq (AGPL,
  output-only), so fixtures come from running the script; rarefaction is
  random, so compare only the deterministic parts
  exactly (\eqn{s_i}, \eqn{t_i}, each sample's total after rarefying) and
  the rest by distribution (the mean over many seeds approaches
  \eqn{x_{ij} / n_i \times L_i}).
- **`qmp_scale()`:** a formula; no reference needed. Fixtures can be
  computed in R in the generator script.


## Tests

- `qmp_scale()`: matches `sweep(counts / rowSums(counts), 1, load, '*')`
  at 1e-12 for every input class and both margins; sparse stays sparse;
  `load` by name, by position and as a `data` column all agree; row sums
  equal `load`.
- `qmp_rarefy()`: \eqn{t_i} and \eqn{s_i} exact; each sample's rarefied
  total equals \eqn{t_i}; each result's row sum equals `load` (to
  rounding); identical for `cpus` 1, 2 and 4 and for the same `seed`;
  differs for another `seed`; the sample with the smallest \eqn{s_i} is
  unchanged apart from scaling.
- `rarefy()` with a per-sample `depth` vector: same as rarefying each
  sample on its own with the same seed; a single depth still gives
  bit-identical results to the current release (store a fixture before
  changing the C).
- Errors: missing, NA, zero or negative `load`; length mismatch; names that
  don't match; an empty sample; non-integer counts for `qmp_rarefy()`.
- Smoke-test lines for both.


## Docs

- One Rd page for both functions (`@rdname qmp`), in a new pkgdown
  reference group "Compositional Data" (see
  [compositional.md](compositional.md)), subtitle "Absolute abundance".
- NEWS, including the `rarefy(depth =)` extension.


## References

- Vandeputte, D., Kathagen, G., D'hoe, K., et al. (2017). Quantitative
  microbiome profiling links gut community variation to microbial load.
  *Nature*, 551, 507–511. doi:10.1038/nature24460
- Jian, C., Luukkonen, P., Yki-Järvinen, H., Salonen, A., & Korpela, K.
  (2020). Quantitative PCR provides a simple and accessible method for
  quantitative microbiota profiling. *PLoS ONE*, 15, e0227285.
- Props, R., Kerckhof, F.-M., Rubbens, P., et al. (2017). Absolute
  quantification of microbial taxon abundances. *The ISME Journal*, 11,
  584–587.


## Decisions made while implementing (2026-11-19)

Asked of the user at the start (all answered with the recommendation):

- **Integer counts:** `qmp_rarefy()` requires whole numbers (as
  `rarefy()`); it does not apply the script's `ceiling()`. The docs show
  `qmp_rarefy(ceiling(x), load)` to do as the script does.
- **Non-integer targets:** floored, as the script's
  `phyloseq::rarefy_even_depth()` truncates them (confirmed from its
  output: 3.53 → 3, 7.41 → 7). But a value within 1e-9 (relative) of a whole
  number counts as that number, and targets are capped at \eqn{n_i}.
  Running the script showed why: \eqn{(n/L) \times L} for the
  smallest-depth sample came out as 101.99999… for \eqn{n = 102}, and the
  script rarefied it to 101 reads. In the fixtures, ecodive keeps 3 reads
  that the script loses (`expected/qmp.rds`, cases `random` and `hair`).
- **A target of 0 reads** (a load so small that \eqn{s_{min} L_i < 1}) is
  an error naming the samples (the script returns `NaN`s). So is an empty
  sample, in both functions.
- **Attributes:** `sampling_depth` and `target` (integer), named by sample,
  are attached to every output class. The docs say that subsetting or
  converting a Matrix/slam object drops them.
- **`rarefy(depth =)` with one depth per sample:** matched by name when both
  are named, else by position, by `per_sample()` (messages name `counts`).
  A length-1 `depth` keeps its old meaning even if named. `drop` and `warn`
  compare each sample with its own depth; the warning then reads "N samples
  have fewer observations than their `depth` and will be ...". No `data`
  argument was added to `rarefy()`.

Made while implementing:

- **`QMP.R` was fetched** (jsDelivr's GitHub mirror) and read; it is BSD-3.
  `R/qmp.r` ports its steps, so it carries the licence text, with a
  "modified port" note; `inst/COPYRIGHTS` has a `Files: R/qmp.r` stanza, and
  `Authors@R` credits "Raes Lab" (ctb, cph). A copy of the script and its
  licence is in `data-raw/expected/qmp/` for the fixture generator.
- **`C_rarefy`** now takes a vector of targets (length 1, or one per sample)
  in place of the static `target`; per-sample PCG seeding is unchanged. slam
  input indexes samples from 1 in that code, so its per-sample targets are
  shifted to match. A fixture of the pre-change outputs
  (`expected/rarefy.rds`, from `data-raw/expected/rarefy.r` run at 4600cb8)
  checks that a single depth is bit-identical for every class, both
  margins, and `cpus` 1, 2 and 4.
- **`qmp_scale()` runs in C** (`C_qmp_scale` in `src/normalize.c`, one pass
  per sample: `x / total * load`, as R evaluates `x / rowSums(x) * load`). It
  does *not* reuse `normalize(norm = 'percent')`, which leaves samples
  already summing to 1 (within 1e-6) as they are, so it would miss the
  1e-15 target. Its output goes through `ecomatrix_table()`, the code
  factored out of `C_normalize`, and then through `table_output()`, the R
  output code factored out of `normalize()`. Both functions return what
  `normalize()` returns: the same class, orientation and dimnames, with
  zeros not stored. `qmp_rarefy()` rarefies with `C_rarefy`, then scales
  with `C_qmp_scale`, dividing by the achieved total, which equals the
  target.
- **`counts`** may be any non-negative finite values for `qmp_scale()`
  (proportions give the same result as counts); negative values or NA are
  errors.
- **`load`** must be numeric (not a factor), finite and > 0, and errors name
  the samples. When `load` comes from `data`, messages name the column.
- **slam draws differ from the other classes'**, for `rarefy()` (already
  the case before this item) and so for `qmp_rarefy()`: slam's samples seed
  their PCG streams by 1-based index. That was left as it is, since changing
  it would change `rarefy()`'s results for slam input. Targets, depths and
  scaling agree for every class.
- **Reference testing:** the script was run on phyloseq 1.46.0 (output
  only). Its sampling depths and targets match exactly, apart from the
  rounding rescue above. Its scaling equals `qmp_scale()` of its own
  rarefied counts to 1e-15. With ecodive's rarefied counts substituted for
  phyloseq's, the script's output equals `qmp_rarefy()`'s with relative
  difference 0. Over 400 seeds, both means agree with each other and with
  \eqn{x / n \times L} (largest |z| 2.2 against the expectation, 2.1
  against the script, over 48 cells).
- **Test-environment note:** Debian's tinytest 1.4.1 exits a test file at
  any `exit_if_not()` call, even a passing one, so the Matrix/slam parts of
  several test files (and all of `test-qmp.r`) don't run with it. With
  tinytest 1.4.3 (CRAN) the suite gave 10,893 results at home and 246 in
  CRAN mode after this item; with 1.4.1, 10,213 at home (10,206 before).


## Outcome (from the implementing session's final report)

**Item 18, done 2026-11-19.** Commits `068453a` (work), `dea9e24` (tick).

- Tests (tinytest 1.4.3): 10,893 at home, 246 in CRAN mode.
- `QMP.R` (via jsDelivr; ported, Raes Lab credited) run with phyloseq 1.46.0:
  sampling depths identical; targets identical except 3 samples where the
  script loses a read to floating point (ecodive's rounding-safe floor keeps
  it); with ecodive's rarefied table substituted, its output equals
  `qmp_rarefy()` exactly; 400 seeds agree by distribution (|z| ≤ 2.2).
- `qmp_scale()` equals `x / rowSums(x) * load` exactly. `rarefy()` with one
  depth is bit-identical to before.
- Open: slam draws differed (fixed later in item 30, `963647f`).
