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
| 4 | [ ] `cca()` | [cca_ca.md](cca_ca.md) | port (skbio) | High demand in ecology; phyloseq "CCA". |
| 5 | [ ] `ca()` | [cca_ca.md](cca_ca.md) | port (skbio) | Cheap once `cca()`'s weighted core exists. |
| 6 | [ ] `procrustes()` and `procrustes_test()` | [procrustes.md](procrustes.md) | port (SciPy, BSD-3) | In QIIME 2; simple. |
| 7 | [ ] `plsda()` | [plsda_oplsda.md](plsda_oplsda.md) | port (scikit-learn, BSD-3) | Popular in microbiome work through mixOmics. |
| 8 | [ ] `simper()` | [simper.md](simper.md) | port (ecopy, MIT) + Clarke's permutation p-values | High demand in community ecology; simple. |
| 9 | [ ] `bioenv()` | [bioenv.md](bioenv.md) | port (skbio) | In QIIME 2; simple. |
| 10 | [ ] `lda()` | [lda_capda.md](lda_capda.md) | port (scikit-learn) | Simple; `capda()` needs it. |
| 11 | [ ] `oplsda()` | [plsda_oplsda.md](plsda_oplsda.md) | port (pyopls, MIT) | Metabolomics demand; small on top of `plsda()`. |
| 12 | [ ] `varpart()` | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | High demand in ecology; formulas only. |
| 13 | [ ] `forward_select()` | [varpart_forward_select.md](varpart_forward_select.md) | derived (on `rda()` / `dbrda()`) | Usually used alongside varpart. |
| 14 | [ ] `capda()` | [lda_capda.md](lda_capda.md) | derived (`pcoa()` + `lda()`) | Anderson's CAP, popular in marine ecology. |
| 15 | [ ] `principal_response_curves()` | [principal_response_curves.md](principal_response_curves.md) | derived (partial `rda()`) | Experiments over time. |
| 16 | [ ] `mantel_correlogram()` | [mantel_correlogram.md](mantel_correlogram.md) | derived (`mantel()` on distance classes) | Distance decay, biogeography. |
| 17 | [ ] `moran_i()` | [moran_i.md](moran_i.md) | port (PySAL esda, BSD-3) | Simple; `dbmem()` needs it. |
| 18 | [ ] `diffusion_map()` | [diffusion_map.md](diffusion_map.md) | port (pydiffmap MIT / scanpy BSD-3) | Growing use; uses neighbour and eigensolver code ecodive already has. |
| 19 | [ ] `pacmap()` | [pacmap.md](pacmap.md) | port (PaCMAP, Apache-2.0) | Growing use; keep the Apache notices. |
| 20 | [ ] `dpcoa()` | [dpcoa.md](dpcoa.md) | **clean room** | Favoured among the clean-room methods: phyloseq "DPCoA", phylogeny-aware. |
| 21 | [ ] `mantel_regression()` (MRM) | [mantel_regression.md](mantel_regression.md) | clean room | Simple algebra; biogeography. |
| 22 | [ ] `mirkat()` | [mirkat.md](mirkat.md) | clean room | In demand in microbiome methods work; needs Davies or Liu p-values written from the papers. |
| 23 | [ ] `splsda()` | [splsda.md](splsda.md) | clean room | mixOmics users; sparse PLS-DA is its own algorithm. |
| 24 | [ ] `coinertia()` | [coinertia.md](coinertia.md) | clean room (RV from hoggorm, BSD-2) | Multi-omics. |
| 25 | [ ] `dbmem()` | [dbmem.md](dbmem.md) | clean room | Spatial ecology. |
| 26 | [ ] `phate()` | [phate.md](phate.md) | clean room; the reference's license is non-commercial | Single-cell crossover; legal care needed even for fixtures. |
| 27 | [ ] `dca()` | [dca.md](dca.md) | clean room; **superseded** (criticised; NMDS preferred) | Kept for users who expect it, written to match vegan. |
| 28 | [ ] `largevis()` | [largevis.md](largevis.md) | port (uwot's BSD-2 headers); **superseded** by UMAP | Cheap, since `umap()` already uses these headers. |

Dropped: `ordisurf()` (it would need a GAM written in C; the user dropped
it on 2026-11-18).


## Existing functions checked for splitting (2026-11-18)

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
