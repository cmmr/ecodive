# Library-size normalizations as `norm` choices: CSS, TMM, RLE, GMPR; `size_factors()`

**Priority:** 43. **Tier:** port for CSS (metagenomeSeq 1.54.0,
Artistic-2.0: keep its notice) + clean room for TMM (edgeR, GPL ≥ 2), RLE
(DESeq2, LGPL ≥ 3) and GMPR (GUniFrac, GPL-3), from the papers, the
manuals, and outputs. Read [README.md](README.md) and
[compositional.md](compositional.md) first.


## What it is

Users coming from metagenomeSeq, edgeR, DESeq2 or GUniFrac expect their
library-size normalizations (CSS and TMM: ~120–145 microbiome mentions a
year each; GMPR 12–20; research notes). Each estimates one **size factor**
dₙ per sample, and the normalized table is xₙ / dₙ. Nixon et al. (2025)
read every such normalization as a point-mass *scale model*, so the same
factors can later centre `dirmult_ttest()`'s scale model (last section).


## Do they fit `norm`?

The rule (compositional.md, qmp.md): a `norm` choice takes no per-sample
input from outside the table and must work inside every function with a
`norm` argument.

- **CSS, TMM, RLE, GMPR: yes.** They need only the table, and every
  function with `norm` normalizes through one C function,
  `normalize(ecomatrix_t *em, ...)` (`src/normalize.c`), called from
  `pca.c` (`pca_table_init()`: `pca`, `rda`, `plsda`, `oplsda`, `dbrda`),
  `beta_div.c`, `envfit.c`, `pcoa_biplot.c`, `simper.c`, `alpha_div.c`, and
  via `C_normalize` (`normalize()`, `lda()`). The ecomatrix holds every
  sample, so a first pass computes the factors and a second divides.
- **New property, to document:** a sample's values depend on the other
  samples in the table (as `philr(part_weights = 'gm.counts')`'s already
  do), so adding or removing samples changes them. `beta_div(pairs =)`
  still uses the whole table (right). `envfit()` and `pcoa_biplot()`
  normalize samples they then skip (`map[j] < 0`, never NA-checked): pass
  `map` to the size-factor pass so only analysed samples set the factors.
  `alpha_div()`'s `norm` formal reaches no metric: nothing to do.
- **Wrench (Kumar et al. 2018): no.** `wrench(mat, condition, ...)` needs
  per-sample groups, and it is rare (≤ 7 mentions a year). Not planned;
  if wanted, port it (Artistic-2.0, Wrench 1.30.0) as a plain
  `size_factors_wrench(counts, groups, data = NULL, margin, cpus)` with
  `size_factors()`'s output.
- **Not planned:** edgeR's `TMMwsp` (the manual's description is too thin
  for clean-room work) and `upperquartile`; G-RLE/FTSS and IRS (licences
  unverified).


## Names

- **`norm` strings** (lowercase, like the existing ones): `'css'`,
  `'tmm'`, `'rle'` (DESeq2 `type = "ratio"`), `'rle_poscounts'` (DESeq2
  `type = "poscounts"`, for tables where every feature has a zero, which
  plain RLE can't handle) and `'gmpr'`.
- **`'rle'` vs `base::rle()`:** a string masks nothing, but say in the
  docs that it is relative log expression (Anders & Huber 2010), not
  run-length encoding (`?rle` finds base R's).
- **Prefixes:** `validate_norm()` uses `match.arg()`. `norm = 'r'` means
  `'rclr'` today and becomes an ambiguity error (NEWS). `'c'` is already
  ambiguous, `'p'` stays `'percent'`; exact matches win, so `'rle'` is
  never `'rle_poscounts'`.
- **`size_factors(counts, norm, margin = 1L, cpus = n_cpus())`**, exported.
  "Size factor" is the field's term for the per-sample divisor (Anders &
  Huber 2010; DESeq2 `sizeFactors()`; the GMPR paper); snake_case (rule
  5); not in NAMESPACE; no package has `size_factors` (BiocGenerics has
  `sizeFactors()`, edgeR `normLibSizes()`, scuttle `librarySizeFactors()`).
  Not `norm_factors()`: edgeR's `norm.factors` are a different quantity (a
  multiplier of the library size). **A function, not an attribute of
  `normalize()`**: the factors are wanted when no table is returned (to
  compare with DESeq2/edgeR, as offsets, as scale-model centres);
  attributes on Matrix/slam objects are lost on subsetting (qmp.md); and
  compositional.md keeps one route per thing. `normalize()`'s output is
  unchanged.


## Design

### Shared behaviour

Normalized value xₙf / dₙ: zeros stay zero, so sparse stays sparse and
`normalize()` keeps the input class. `pseudocount` is ignored, as for
every non-clr `norm`. Negative values or NA are errors ("cannot perform
TMM normalization", worded like the Hellinger check); non-integers are
allowed (the docs recommend raw counts). A sample whose dₙ is undefined
(empty, no usable features) is an error naming it and suggesting another
`norm` (a `norm` takes no arguments to fix it). dₙ is scaled so that the
normalized table is each reference's familiar output:

| `norm` | dₙ | normalized table equals | attributes of `size_factors()` |
|:-------|:---|:------------------------|:-------------------------------|
| `'css'` | normFactorₙ / 1000 | `metagenomeSeq::MRcounts(obj, norm = TRUE)` (`sl = 1000`) | `norm_factors` (metagenomeSeq's `normFactors`), `quantile` (p) |
| `'tmm'` | Nₙ fₙ / 10⁶ | `edgeR::cpm()` with TMM library sizes | `norm_factors` (edgeR's fₙ, product 1), `lib_size` (Nₙ), `reference` |
| `'rle'`, `'rle_poscounts'` | DESeq2's sₙ | `counts(dds, normalized = TRUE)` | — |
| `'gmpr'` | GMPR's sₙ | x / sₙ (the paper) | — |

`size_factors()` also takes `'percent'` (dₙ = Nₙ, the sample total; note
`normalize(norm = 'percent')` leaves samples already summing to 1 alone)
and `'none'` (1); other `norm`s error ("not a size-factor
normalization"). It returns a numeric vector named by sample, attribute
`norm`, and `normalize(x, n)` is `identical()` to dividing each sample by
`size_factors(x, n)`.

### `'css'`: port (metagenomeSeq `cumNormStatFast.R`, `cumNorm.R`, `cumNormMat.R`)

Read on GitHub master (1.49.1, Artistic-2.0); port from the Bioconductor
release **1.54.0** after checking these files are unchanged.

1. **Percentile p** (`cumNormStatFast(rel = 0.1)`). Sort each sample's
   non-zero values; L = the largest number of non-zero values in any
   sample. A sample with exactly one non-zero value is an error (as
   metagenomeSeq); an empty one is an error naming it (metagenomeSeq's
   `== 1` check misses it). ref[r], r = 1..L: the mean over samples of an
   L-row matrix holding each sample's values ascending at the bottom, 0
   above. q_i[r]: type-7 quantile of sample i's non-zero values at
   `seq(0, 1, length.out = L)`. d[r] = median over samples of |ref[r] −
   q_i[r]| (matrixStats `rowMedians`: mean of the two middle values for an
   even count). p = (first k with |d[k+1] − d[k]| / d[k+1] > 0.1) / L; if
   p ≤ 0.5, p = 0.5 (metagenomeSeq messages "Default value being used.";
   ecodive is silent and records p). If no k qualifies (identical samples,
   n = 1), metagenomeSeq fails (`if (NA <= 0.5)`): ecodive uses 0.5 and
   says so under "Differences".
2. **Factor** (`calcNormFactors(x, p)`): qₙ = type-7 p-quantile of sample
   n's **non-zero** values; normFactorₙ = Σ (xₙf − ε) over **all**
   features with xₙf − ε ≤ qₙ, ε = `.Machine$double.eps` (zeros add −ε
   each: keep the quirk).
3. **No L × n matrix.** Sort non-zeros per sample (threads over samples);
   ref[r] per row, summing samples in sample order; q_i[r] on the fly; row
   medians by quickselect (threads over rows). Results independent of
   `cpus`.
4. **Type-7** by its documented formula (Hyndman & Fan; `?quantile`):
   x₍ⱼ₎ + γ(x₍ⱼ₊₁₎ − x₍ⱼ₎), j = ⌊1 + (m − 1)p⌋. Don't read R's GPL source;
   matrixStats (Artistic-2.0) may be read. p is a discrete choice, so the
   fixtures must give **identical** p.

### `'tmm'`: clean room (Robinson & Oshlack 2010; edgeR 4.10.7 manual, `normLibSizes`)

Defaults from the manual, fixed: `logratioTrim = 0.3`, `sumTrim = 0.05`,
`doWeighting = TRUE`, `Acutoff = -1e10`. Nₙ = sample totals; features
zero in every sample are removed (no effect on f).

- **Reference r:** "the column whose count-per-million upper quartile is
  closest to the mean upper quartile" (manual); first on ties.
- **Sample k against r** (k = r: f = 1): features with xₖ > 0 and xᵣ > 0;
  M = log₂((xₖ/Nₖ)/(xᵣ/Nᵣ)), A = ½ log₂((xₖ/Nₖ)(xᵣ/Nᵣ)), A > Acutoff;
  trim 30% of M and 5% of A from each tail; log₂ fₖ = Σ wM / Σ w over the
  rest, w = 1 / [(Nₖ − xₖ)/(Nₖxₖ) + (Nᵣ − xᵣ)/(Nᵣxᵣ)] (the paper's
  delta-method variance).
- **Rescale** so ∏ f = 1 (manual).
- **Settle from crafted outputs:** the upper quartile's quantile type and
  that zeros count in it; how many features each trim removes (the paper
  gives fractions, not index rules) and how ties in M or A rank; f when no
  feature survives; edgeR NEWS's case of identical counts with unequal
  library sizes.
- Threads over samples, O(D log D) each.

### `'rle'`, `'rle_poscounts'`: clean room (Anders & Huber 2010; DESeq2 1.52.0 manual)

- `'rle'`: log gᶠ = mean over samples of log xₙf, features with any zero
  excluded; sₙ = exp(median over them of (log xₙf − log gᶠ)). No feature
  positive everywhere → error suggesting `'rle_poscounts'` or `'gmpr'`.
- `'rle_poscounts'`: log gᶠ = (Σ over positive xₙf of log xₙf) / n over all
  n samples ("the n-th root of the product of the non-zero counts"), for
  features with any positive count.
- **Settle from outputs:** for an even count, median of log ratios
  (geometric mean of the middle two) or of ratios; for `poscounts`,
  whether zeros leave each sample's median; whether factors are rescaled
  to geometric mean 1 (the manual says so when `geoMeans` is given).
- Feature log-sums in fixed sample blocks summed in block order
  (independent of `cpus`); medians with threads over samples.

### `'gmpr'`: clean room (Chen et al. 2018, CC-BY; GUniFrac 1.9 manual)

rᵢⱼ = median over shared features of xₖᵢ / xₖⱼ, computed in both
directions (for an even count, median(a/b) ≠ 1/median(b/a)); sᵢ =
(∏ⱼ rᵢⱼ)^(1/nᵢ) over the defined rᵢⱼ. Fixed manual defaults: a feature
counts only with ≥ `min_ct = 2` counts; a pair only with ≥ `intersect_no
= 4` shared features. **Settle from outputs:** whether `min_ct` applies to
each count; whether j = i is in the product (the paper's runs over
j = 1..n, rᵢᵢ = 1); the result for a sample with no qualifying pair;
rescaling. O(n² × shared) by merging sorted non-zero lists; threads over
i, logs summed in j order.

### C and R layout

- `src/norm_css.c`: the Artistic-2.0 text and copyright (from
  metagenomeSeq's DESCRIPTION authors; *verify* the exact notice in the
  1.54.0 tarball) and "C port of `cumNormStatFast()` and
  `calcNormFactors()` from metagenomeSeq 1.54.0; restructured for ecodive,
  see plans/normalization_factors.md". `src/size_factors.c` (MIT): TMM,
  RLE, GMPR.
- `normalize()` gains `const int *use` (NULL = all); codes `NORM_CSS 7`
  to `NORM_GMPR 11` are appended (`rpca.c` hard-codes 5); update every
  caller. `C_size_factors` serves `size_factors()` (`R/size_factors.r`).
- Update every list of `norm` names: `validate_norm()`, `NORM_*` and
  `validate_pseudocount()`'s negativity check (`R/validate.r`); messages
  in `R/normalize.r`; `norm_name()` (`R/dist_clust.r`); `varpart.r:355`;
  `simper()`'s message (the new norms are allowed: non-negative); the
  `@param norm` lists in `R/documentation.r`, `pca.r`, `pcoa_biplot.r`,
  `envfit.r`, `varpart.r`.
- `inst/COPYRIGHTS` stanza (format as ANCOMBC's); `Authors@R`:
  metagenomeSeq authors (ctb, cph).

### Centring a scale model later (not this item)

The normalization implies sample n's total Wₙ ∝ Σ_f xₙf / dₙ = Nₙ / dₙ,
i.e. μₙ = log₂ Nₙ − log₂ dₙ (log₂(1/fₙ) for TMM). In each Dirichlet draw,
log₂ absolute abundance = log₂ pₙf + μₙ, so the treatment-vs-reference
shift added to every clr difference is mean_T(ḡₙ + μₙ) − mean_R(ḡₙ + μₙ),
ḡₙ = mean_f log₂ pₙf, with `scale_sd` the spread around it; the clr
assumption μₙ = −ḡₙ gives today's 0. This is ALDEx3's `sample.sm(s.mu =
μ, s.sd)` (MIT, read). [dirmult_lm.md](dirmult_lm.md) owns the argument
(suggested `scale_mu`: per-sample log₂ values, or a `norm` name resolved
by `size_factors()`); a measured `load` gives μₙ = log₂ loadₙ.


## Fixtures (`data-raw/expected/normalization_factors.r`)

Install metagenomeSeq, edgeR, DESeq2 and GUniFrac in a scratch library;
tests never call them. Save `cumNormStatFast()`, `calcNormFactors()`,
`MRcounts(norm = TRUE)`; `normLibSizes(method = "TMM")` and `cpm()`;
`estimateSizeFactorsForMatrix(type = "ratio" / "poscounts")`; `GMPR()`
(features in rows). Tables: (a) 6 × 40, all features positive; (b) a
simulated 30 × 300 table, ~70% zeros; (c) crafted cases: ties in M and A,
even counts for every median, identical samples, a sample sharing < 4
features with all others, counts of 1 (`min_ct`), and shared-feature
counts 2–40 to pin down TMM's trim rule. Keep them small
(`saveRDS(version = 2)`).


## Plotting

Follows [README.md](README.md) § "Result attributes and plots".

- No plots. `normalize()` returns a plain table, as today: no `display`
  and no attributes saying how it was normalized. A metric or ordination
  run with `norm = 'tmm'` records nothing about the norm either.
- `size_factors()` returns the named vector with only its result
  attributes: `norm`, plus the computed values in the table above
  (`norm_factors`, `quantile`, `lib_size`, `reference`).

## Tests (`tests/tinytests/test-normalization_factors.r`)

- Factors and normalized tables match the fixtures at 1e-12 (relative);
  CSS's p identical.
- `normalize(x, n)` `identical()` to `sweep(x, 1, size_factors(x, n),
  '/')` for every input class and both margins; sparse stays sparse.
- Through every `norm` function: e.g. `bray(norm = 'tmm')` equals
  `bray(normalize(x, 'tmm'))`, likewise `pca`, `rda`, `lda`, `plsda`,
  `oplsda`, `simper`, `envfit(counts =)`, `pcoa_biplot`, `varpart`,
  `forward_select`. Adding an unmatched sample (even with NA) leaves
  `envfit()` and `pcoa_biplot()` unchanged.
- `identical()` for `cpus` 1, 2, 4.
- Edge cases: n = 1 and 2; an empty sample; one non-zero (CSS); no
  all-positive feature (`'rle'`); no qualifying GMPR pair; negatives;
  `norm = 'r'` now ambiguous; `size_factors(norm = 'clr')` error. Each
  error message checked.
- Smoke: `normalize(ex_counts, 'gmpr')`, `size_factors(ex_counts, 'css')`.


## Docs

`normalize()`'s `@param norm`: one bullet per method with its reference;
`@details`: values depend on every sample in the table; `'rle'` is not
`base::rle()`. `size_factors()` Rd: the table above (how to get edgeR's
`norm.factors` and metagenomeSeq's `normFactors`), "Differences from
metagenomeSeq" (Artistic-2.0 §4: the silent p fallback, the empty-sample
error), "Same names in other packages". pkgdown: "Compositional Data" →
*Absolute abundance*, beside `qmp_scale`. NEWS: the new norms and the
`'r'` ambiguity. (compositional_docs.md's machine-learning note stays
true: these are scalings, not log ratios.)


## References

- Paulson, J. N., Stine, O. C., Bravo, H. C., & Pop, M. (2013). Differential abundance analysis for microbial marker-gene surveys. *Nature Methods*, 10, 1200–1202. doi:10.1038/nmeth.2658
- Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11, R25. doi:10.1186/gb-2010-11-3-r25
- Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. *Genome Biology*, 11, R106. doi:10.1186/gb-2010-11-10-r106
- Chen, L., Reeve, J., Zhang, L., Huang, S., Wang, X., & Chen, J. (2018). GMPR: a robust normalization method for zero-inflated count data with application to microbiome sequencing data. *PeerJ*, 6, e4600. doi:10.7717/peerj.4600
- Kumar, M. S., Slud, E. V., Okrah, K., Hicks, S. C., Hannenhalli, S., & Corrada Bravo, H. (2018). Analysis and correction of compositional bias in sparse sequencing count data. *BMC Genomics*, 19, 799. doi:10.1186/s12864-018-5160-5
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization. *Genome Biology*, 26, 139. doi:10.1186/s13059-025-03609-3
- Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. *The American Statistician*, 50, 361–365. doi:10.1080/00031305.1996.10473566
