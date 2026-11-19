# To-do: new methods, by priority

Work top to bottom. Each item links to its plan; read [README.md](README.md)
first. Tick an item, and add the commit, when its method is merged with tests
and docs.

The order follows three rules the user set (2026-11-18):

1. Favour methods that are **simple**, that have an **obvious permissively
   licensed reference** to port, and that **microbiome researchers ask for
   now**.
2. **Clean-room** methods go to the bottom, with DPCoA first among them.
   "Clean room" means there is no permissive code to port, so a substantial
   algorithm must be written from the paper and checked against GPL outputs.
3. **Superseded** methods go last.

Methods marked *derived* have no permissive reference code, but are thin
layers over ecodive's own functions plus a published formula (varpart is
adjusted R² from partial RDAs, for example). They aren't clean-room in a
meaningful sense, so they sit just after the methods they build on.

| # | Do | Plan | Tier | Why here |
|--:|:---|:-----|:-----|:---------|
| 1 | [x] Split `mantel_partial()` out of `mantel(given =)` (3074a20) | [mantel_partial.md](mantel_partial.md) | housekeeping | It is a distinct named test (vegan's `mantel.partial()`), and the API is still unreleased, so splitting it now is cheap. |
| 2 | [x] Pairwise tests: `pairwise` and `p_adjust` on `permanova()`, `anosim()`, `mrpp()` and `permdisp()` (575da24) | [pairwise_tests.md](pairwise_tests.md) | port (q2-diversity, BSD-3) | In nearly every microbiome paper; simple, built on existing tests. |
| 3 | [x] `rda()`, including partial RDA (46db60b) | [rda.md](rda.md) | port (skbio, BSD-3) + ecodive's dbrda core | High demand (tb-RDA, phyloseq "RDA"). Several later items build on its core. |
| 4 | [x] `cca()` (e33c7e1) | [cca_ca.md](cca_ca.md) | port (skbio) | High demand in ecology; phyloseq "CCA". |
| 5 | [x] `ca()` (80c8a86) | [cca_ca.md](cca_ca.md) | port (skbio) | Cheap once `cca()`'s weighted core exists. |
| 6 | [x] `procrustes()` and `procrustes_test()` (a1c2a66) | [procrustes.md](procrustes.md) | port (SciPy, BSD-3) | In QIIME 2; simple. |
| 7 | [x] `plsda()` (b3c81a4) | [plsda_oplsda.md](plsda_oplsda.md) | port (scikit-learn, BSD-3) | Popular in microbiome work through mixOmics. |
| 8 | [x] `simper()` (dbf0c2b) | [simper.md](simper.md) | port (ecopy, MIT) + Clarke's permutation p-values | High demand in community ecology; simple. |
| 9 | [x] `bioenv()` (96381a6) | [bioenv.md](bioenv.md) | port (skbio) | In QIIME 2; simple. |
| 10 | [x] `lda()` (7323057) | [lda_capda.md](lda_capda.md) | port (scikit-learn) | Simple; `capda()` needs it. |
| 11 | [x] `oplsda()` (4c7e646) | [plsda_oplsda.md](plsda_oplsda.md) | port (pyopls, MIT) | Metabolomics demand; small on top of `plsda()`. |
| 12 | [x] `varpart()` (eb06b58) | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | High demand in ecology; formulas only. |
| 13 | [x] `forward_select()` (d0b6339) | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | Usually used alongside varpart. |
| 14 | [x] `alr()`, `ilr()` (51fd7d0) | [log_ratio_transforms.md](log_ratio_transforms.md) | port (scikit-bio, BSD-3) | Simplest compositional gap; `balance()` and `philr()` build on the ILR code. |
| 15 | [x] `replace_zeros()` wrapper (+ `list_zero_replacements()`, `match_zero_replacement()`) and `replace_zeros_mult()` (7d84f63) | [zero_replacement.md](zero_replacement.md) | port (scikit-bio, BSD-3) | Principled zeros for every log-ratio method; the wrapper is ready for the Bayesian methods (later item). |
| 16 | [x] `log_ratio()`, `balance()` (c56ea2a) | [log_ratio.md](log_ratio.md) | port (Qurro / scikit-bio, BSD-3) | Trivial; log ratios are what current compositional guidance recommends analysing. |
| 17 | [ ] `vlr()`, `rho_p()`, `phi_s()` | [proportionality.md](proportionality.md) | port (scikit-bio `vlr`) + formula (propr output-only) | High demand (propr); all-pairs C over features is ecodive's strength. |
| 18 | [ ] `qmp_scale()`, `qmp_rarefy()`; per-sample `rarefy(depth =)` | [qmp.md](qmp.md) | port (raeslab QMP.R, BSD-3) + formula | Simple, over `rarefy()`'s C; the absolute-abundance side of the debate. |
| 19 | [ ] `ancom()` | [ancom.md](ancom.md) | port (scikit-bio, BSD-3) | In QIIME 2 and widely used; D² independent tests suit threaded C. |
| 20 | [ ] `dirmult_ttest()` (then its optional scale model) | [dirmult_ttest.md](dirmult_ttest.md) | port (scikit-bio, BSD-3; scale model from ALDEx3, MIT) | ALDEx2-style effect sizes, high demand; needs a PCG Gamma sampler. |
| 21 | [ ] `rpca_phylo()` | [rpca_phylo.md](rpca_phylo.md) | port (gemelli, BSD-3) on `rpca()` | QIIME 2's phylogenetic RPCA; reuses `rpca()` and the UniFrac tree walk. Phylogeny-aware. |
| 22 | [ ] `ancombc()`, `ancombc2()` | [ancombc.md](ancombc.md) | port (scikit-bio, BSD-3, derived from Artistic-2.0 ANCOMBC: keep both notices) | Among the most used DA methods now; large (~1,000 lines of core). |
| 23 | [ ] `ancombc_global()`, `ancombc_pairwise()`, `ancombc_dunnett()`, `ancombc_trend()`, `structural_zeros()` | [ancombc.md](ancombc.md) | port (scikit-bio + ANCOMBC notices), as above | Post-hoc tests on a fitted `ancombc()` / `ancombc2()`; `ancombc_trend()` needs a small QP solver. |
| 24 | [ ] `capda()` | [lda_capda.md](lda_capda.md) | derived (`pcoa()` + `lda()`) | Anderson's CAP, popular in marine ecology. |
| 25 | [ ] `principal_response_curves()` | [principal_response_curves.md](principal_response_curves.md) | derived (partial `rda()`) | Experiments over time. |
| 26 | [ ] `mantel_correlogram()` | [mantel_correlogram.md](mantel_correlogram.md) | derived (`mantel()` on distance classes) | Distance decay, biogeography. |
| 27 | [ ] `moran_i()` | [moran_i.md](moran_i.md) | port (PySAL esda, BSD-3) | Simple; `dbmem()` needs it. |
| 28 | [ ] `diffusion_map()` | [diffusion_map.md](diffusion_map.md) | port (pydiffmap MIT / scanpy BSD-3) | Growing use; uses neighbour and eigensolver code ecodive already has. |
| 29 | [ ] `pacmap()` | [pacmap.md](pacmap.md) | port (PaCMAP, Apache-2.0) | Growing use; keep the Apache notices. |
| 30 | [ ] `rpca_joint()` | [rpca_joint.md](rpca_joint.md) | port (gemelli, BSD-3) | Multi-omics; a large port (~1,000 C lines) on `rpca()`'s OptSpace code. |
| 31 | [ ] `ctf()` | [ctf.md](ctf.md) | port (gemelli, BSD-3) | Longitudinal designs; the largest gemelli port, with reference quirks to settle with the user first. |
| 32 | [ ] `philr()`, `philr_distance()` | [philr.md](philr.md) | **clean room** (philr GPL-3; formulas from the open-access paper) | First among the clean-room methods: phylogeny-aware, high demand, and an O(D) formula per sample. |
| 33 | [ ] `dpcoa()` | [dpcoa.md](dpcoa.md) | **clean room** | Second among the clean-room methods, after `philr()`: phyloseq "DPCoA", phylogeny-aware. |
| 34 | [ ] `replace_zeros_czm()`, `replace_zeros_gbm()`, `replace_zeros_sq()`, `replace_zeros_bl()` | [zero_replacement.md](zero_replacement.md) | **clean room** (zCompositions GPL; the defining paper is closed access, so formulas are confirmed from outputs) | Simple once the formulas are pinned down; `cmultRepl()`'s GBM is the field's default. |
| 35 | [ ] `mantel_regression()` (MRM) | [mantel_regression.md](mantel_regression.md) | clean room | Simple algebra; biogeography. |
| 36 | [ ] `mirkat()` | [mirkat.md](mirkat.md) | clean room | In demand in microbiome methods work; needs Davies or Liu p-values written from the papers. |
| 37 | [ ] `splsda()` | [splsda.md](splsda.md) | clean room | mixOmics users; sparse PLS-DA is its own algorithm. |
| 38 | [ ] `coinertia()` | [coinertia.md](coinertia.md) | clean room (RV from hoggorm, BSD-2) | Multi-omics. |
| 39 | [ ] `dbmem()` | [dbmem.md](dbmem.md) | clean room | Spatial ecology. |
| 40 | [ ] `phate()` | [phate.md](phate.md) | clean room; the reference's license is non-commercial | Single-cell crossover; legal care needed even for fixtures. |
| 41 | [ ] `dca()` | [dca.md](dca.md) | clean room; **superseded** (criticised; NMDS preferred) | Kept for users who expect it, written to match vegan. |
| 42 | [ ] `largevis()` | [largevis.md](largevis.md) | port (uwot's BSD-2 headers); **superseded** by UMAP | Cheap, since `umap()` already uses these headers. |

**Compositional methods (added 2026-11-19).** The plans in
[compositional.md](compositional.md) were inserted by the same rules. The
permissive ports come first because they are simple and in high demand in
microbiome work now. `rho_p()`, `phi_s()` and `qmp_scale()` are formulas
on top of ported code, so they count as derived, like varpart. `philr()`
and the Bayesian zero replacements have no permissive code, so they go
with the clean-room methods. `philr()` leads that group as the simpler
of its two phylogeny-aware methods.

Dropped: `ordisurf()` (it would need a GAM written in C; the user dropped
it on 2026-11-18).


## Existing functions checked for splitting

2026-11-19: checked again for the compositional plans (`normalize(norm)`,
`pseudocount`, `rarefy()`). Nothing needs splitting; see
[compositional.md](compositional.md).

### 2026-11-18

The rule is one function per named method variant (see README).

**Split:**
- `mantel(given =)` becomes `mantel_partial()` (item 1).

**Kept as arguments,** because every reference treats them as options of
one method, not as separate methods:

| Argument | Why it stays an argument |
|:---------|:-------------------------|
| `pcoa(correction)` | ape, vegan and `cmdscale()` all treat Lingoes and Cailliez as options of PCoA. |
| `permdisp(type)` | vegan and Anderson's PERMDISP2 both treat centroid and median as options. |
| `permdisp(bias_adjust)` | An option in vegan too. |
| `mantel(method)` | The choice of correlation coefficient. |
| `mantel(alternative)` | A test option. |
| `nmds(ties)` | An option of SMACOF in scikit-learn and vegan. |
| `tsne(gradient)` | A computation route (exact or Barnes-Hut), not a different method. |
| `tsne(init)`, `umap(init)` | Starting layouts. |
| `pca(scale)` | `prcomp()`'s `scale.`. |
| `pca(norm)` | The transformation applied first, as in every function with `norm`. |
| `permanova(by)` | Which ANOVA table to compute. |
| `mrpp(weights)` | Group weights. |
| `rpca(distances)` | An extra output. |
| `envfit(data / counts)` | What is fitted, not how. |

**Already split:** `mds()` / `nmds()`, `dbrda()` / `capscale()`, and the
UniFrac family.
