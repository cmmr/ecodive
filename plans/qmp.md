# `qmp_scale()` and `qmp_rarefy()`: quantitative microbiome profiling

**Priority:** see [TODO.md](TODO.md). **Tier:** formula (both are a few
lines over ecodive's own `rarefy()` C code). Read [README.md](README.md)
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

- **Check steps 3 and 4 against the authors' script** (raeslab's `QMP`
  R code, `rarefy_even_sampling_depth()`), whose license must be
  established before reading it (see "References"): whether it uses
  floor or round in step 3, and whether step 4 divides by \eqn{t_i} (as
  above) or by \eqn{s_{min}} (the two differ only by the rounding in
  step 3). Match it.
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

- **`qmp_rarefy()`:** the authors' R script, `rarefy_even_sampling_depth()`
  (GitHub, raeslab/QMP). **Before reading its source, check its license.**
  If it is MIT/BSD/Apache, port it and credit it as README.md requires.
  If it is GPL or unlicensed, it is output-only: write the function from
  the paper (the four steps above) and use the script only to make
  fixtures. Rarefaction is random, so compare only the deterministic parts
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
