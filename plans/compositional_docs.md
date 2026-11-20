# Documentation: depth, scale, machine learning and predicted loads

**Priority:** 34. **Tier:** docs (roxygen comments, two articles, NEWS; no
new functions, no changed results). Read [README.md](README.md) first, and
[compositional.md](compositional.md) § "The 2026-11-20 review" for why.


## What it is

The 2026-11-19 literature review (`reports/Compositional methods for
ecodive.md`, `research_notes/Compositional methods for ecodive/`; out of
git) found four things that bear on existing functions, none of which
their help pages mention: Schloss's (2026) rebuttal of robust Aitchison
PCA (`rpca()`'s method) on uneven depth; FDR failures under misspecified
scale (`dirmult_ttest(scale_sd =)`); no machine-learning advantage for
CoDa transforms (Yerke et al. 2024); and predicted microbial loads
(Nishijima et al. 2025), which users will pass to `qmp_scale()`. Each page
should say so, cite it, and not overclaim in either direction.


## Ground rules

- Roxygen comments only; run `roxygen2::roxygenise()` once (README
  "Tooling"). Never edit `man/*.Rd`.
- **Use only the claims in this table**, quoted or paraphrased without
  going further. Never write that RPCA corrects or controls for depth, nor
  that rarefaction is required: present Schloss as a counterpoint with a
  check users can run. Label ecodive's own reasoning as inference.
- Name only functions that exist when this item is done.
  `beta_div_rarefied()` (item 37), `da_sensitivity()` (41) and
  `dirmult_lm()` (50) come later; their plans update these pages.
- Article chunks use only ecodive and base R. Describe a chunk's output
  only after knitting it.

| Claim (as it may be quoted) | Source | Checked |
|:----------------------------|:-------|:--------|
| Schloss "found numerous problems in the Python code that led to the support of robust Aitchison PCA"; corrected, "robust Aitchison PCA was not insensitive to uneven sequencing effort and did not perform as well as rarefaction"; "even using the benchmarking framework from the original analysis, rarefaction outperformed robust Aitchison PCA, other compositional data analysis methods, and other normalization methods" | Schloss 2026 | Europe PMC abstract, 2026-11-20 |
| Rarefaction "is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses" (the title) | Schloss 2024 | Crossref |
| In **differential set analysis** (enrichment of groups of taxa/genes), small errors in the scale assumption gave positive predictive values "as low as 9%" | McGovern et al. 2023 | abstract |
| Over 11 RNA-seq and other sequencing datasets, "scale misspecification results in poor control of the FDR by several commonly used tools"; "Implementing a scale model in ALDEx2 or ALDEx3 ameliorated unacceptably high FDRs; however, there was an inherent trade-off between satisfactory FDR control and high sensitivity- no tool offered both"; more scale uncertainty raised the smallest difference reported | Dos Santos et al. 2026 | abstract |
| Random forests, 65 metadata variables, 4 datasets: alr, clr and ilr (PhILR) "performed generally slightly worse or only as well as compositionally naïve transformations"; relative-abundance transforms did best "by a small but reliably statistically significant margin" | Yerke et al. 2024 | abstract |
| Load (cells/g) predicted "solely from relative abundance data"; "microbial load is the major determinant of gut microbiome variation"; trained on GALAXY/MicrobLiver (n = 1,894) and MetaCardis (n = 1,812) flow-cytometry loads; cross-validated Pearson r ≈ 0.67–0.68, cross-cohort r = 0.56 | Nishijima et al. 2025 | abstract; numbers from the research notes (paper PDF) |
| `MLP(input, profiler, training_data, output)`; `output = "load"` or `"qmp"` (= relative abundance × predicted load); profilers mOTUs 2.5/3, MetaPhlAn 3/4, DADA2 + RDP (genus 16S); training data `"galaxy"` or `"metacardis"` | MLP README | fetched 2026-11-20 |

Also fix `R/dirmult_ttest.r`'s Nixon et al. (2025) reference, which has
the preprint's title. Crossref: "Incorporating scale uncertainty in
microbiome and gene expression analysis as an extension of
normalization", *Genome Biology*, 26, 139.


## Changes

### 1. RPCA family (`R/rpca.r`, `R/rpca_phylo.r`, `R/beta_div.r`)

- **`rpca()`:** new `@section Sequencing depth:` after **Distances**:
  > `rpca()` is not a correction for uneven sequencing depth. The robust
  > clr divides each sample by the geometric mean of its own non-zero
  > counts, so scaling a sample's counts changes nothing; but a deeper
  > sample has more non-zero features, and which values are missing (the
  > zeros OptSpace completes) depends on depth. Schloss (2026) corrected
  > "numerous problems" in the benchmark that had supported robust
  > Aitchison PCA and found it "was not insensitive to uneven sequencing
  > effort and did not perform as well as rarefaction" (see also Schloss
  > 2024). To check a result, correlate its axes with log library size and
  > compare with a rarefied ordination:
  > `depth <- rowSums(counts); ord <- rpca(counts)`;
  > `cor(ord, log(depth[rownames(ord)]))` (per axis);
  > `envfit(ord, data.frame(log_depth = log(depth)))` (r², permutation p);
  > `pcoa(robust_aitchison(rarefy(counts)))`.

  Add both Schloss papers to `@references` after Martino et al. (2019).
- **`rpca_phylo()`:** `@inheritSection rpca Sequencing depth` and the two
  references, nothing more (Schloss tested `rpca()`'s method).
- **`robust_aitchison()`** ("specifically engineered for sparse
  datasets"): one `@details` sentence, an inference: "It uses the same
  robust clr as [rpca()], so the caveat about sequencing depth there
  applies too." Cite Schloss 2026. *Unverified:* whether Schloss tested the
  plain robust Aitchison distance; read the paper before writing more.
- **Planned `rpca_joint()` (30) and `ctf()` (31):** their Docs steps should
  add `@inheritSection rpca Sequencing depth` and both Schloss references;
  `rpca_joint()` should cite its now-published paper (Cordazzo Vargas,
  Martino, et al. 2026, *Mol. Syst. Biol.* 22, 1789–1808,
  doi:10.1038/s44320-026-00236-3), replacing its plan's "cite once
  published" note. This item does not edit those plans; whoever maintains
  them should copy this bullet in.

### 2. `dirmult_ttest()` (`R/dirmult_ttest.r`)

- `@param scale_sd`: append "Check how your calls depend on it (see
  **Choosing `scale_sd`**)."
- `@details`, after **Scale uncertainty**:
  > **Choosing `scale_sd`.** Report how the results depend on it: rerun
  > with `scale_sd = 0`, `0.5` and `1` (`n_draws = 1000` or more) and say
  > which features stay rejected. Wrong scale assumptions are costly: in
  > differential set analysis, McGovern et al. (2023) found positive
  > predictive values as low as 9% under small errors in the assumed scale;
  > over 11 RNA-seq and other sequencing datasets, Dos Santos et al. (2026)
  > found that misspecified scale broke false-discovery control in several
  > common tools, and that a scale model (ALDEx2, ALDEx3) reduced the
  > unacceptably high FDRs at a cost in sensitivity. Larger `scale_sd` raises the smallest change that
  > can be called. The package website's article "Compositional Analysis
  > of Microbiome Data" runs such a comparison.
- `@references`: correct Nixon et al.; add McGovern et al. (2023) and Dos
  Santos et al. (2026).

### 3. Machine learning (`R/documentation.r`, then four pages)

A shared block like `pseudocount_section`: `@name ml_section`,
`@keywords internal`, `@section Machine learning:`
> Log-ratio transforms (clr, rclr, alr, ilr, PhILR) serve interpretation
> and compositional geometry; don't expect them to improve prediction. With random forests, Yerke et al.
> (2024) found that alr, clr and ilr (PhILR) "performed generally slightly
> worse or only as well as compositionally naïve transformations", and
> relative abundances did best, by a small margin. Compare with
> `norm = 'percent'` before choosing a log-ratio predictor table.

`@inheritSection ml_section Machine learning` in `normalize()` (covers
`'clr'`/`'rclr'`), `alr()`, `ilr()` (`R/log_ratio_transforms.r`) and
`philr()`, each citing Yerke et al. None promises gains now, so nothing is
removed. (`philr()` cites Silverman et al.'s title "A phylogenetic
transform enhances analysis ..."; fine as a title, never paraphrase it as
a claim.)

### 4. Predicted loads (`R/qmp.r`)

New `@section Predicted loads:` after **Before and after**:
> Without a measured load, one can be predicted from a fecal profile with
> the MLP tool of Nishijima et al. (2025;
> \url{https://microbiome-tools.embl.de/mlp/}; an R package not on CRAN):
> `MLP(profile, profiler, training_data, "load")` gives cells per gram, to
> pass as `load` (`qmp_scale(counts, pred$load)`). MLP's own `"qmp"`
> output is the same calculation as `qmp_scale()`. But:
> * **A prediction is not a measurement.** It is computed from the same
>   relative profile, adds no independent information, has errors tied to
>   composition, and can't confirm a compositional result. Predicted and
>   measured loads correlated at r ≈ 0.67–0.68 within, and 0.56 across,
>   the training cohorts.
> * **Human stool only**, profiled with the supported tools (mOTUs,
>   MetaPhlAn, DADA2 + RDP), and trained on European cohorts.
> * That load is "the major determinant of gut microbiome variation" is a
>   reason to measure or adjust for load, not a validation of any one
>   prediction.

Cite Nishijima et al. *Unverified:* whether `MLP()` returns a data frame
whose row names are sample IDs (the README shows `load$load`); install it
in a scratch library and, if so, show `qmp_scale(counts, 'load', data =
pred)`. Don't copy MLP code or model files ("BSD license" string, no
LICENSE file; see compositional.md).

### 5. Article `vignettes/articles/compositional.Rmd`

- **Beta diversity**, after the `rpca` chunk and its "A clr PCA ..."
  paragraph (~line 337): a paragraph **Is it depth?** (Schloss, two
  sentences), then a chunk on the simulated data, which has `depth`
  (5,000–40,000 reads per sample):
  `cor(ord, log(depth[match(rownames(ord), rownames(counts))]))`;
  `ord_rar <- pcoa(robust_aitchison(rarefy(counts, seed = 1)), n_axes = 2)`;
  `procrustes_test(ord_rar, ord, seed = 1)` (signature `(ord, target,
  ...)`; check whether `rarefy()` drops samples here). Describe what it
  shows.
- **Differential abundance**, after "... could differ between your groups"
  (~line 608): one sentence each on Dos Santos et al. and McGovern et al.
- **Absolute abundance**, after the QMP introduction: one sentence on
  predicted loads and "a prediction is not a measurement". No chunk.
- **Checklist:** add "6. **Check depth sensitivity** of any ordination that
  doesn't rarefy: correlate its axes with log library size and compare with
  a rarefied analysis"; extend item 4 with "and, for `dirmult_ttest()`,
  more than one `scale_sd`".
- **References → Methods and reviews**, alphabetical, existing format:
  Dos Santos 2026, McGovern 2023, Nishijima 2025 (mark † unless open access
  is confirmed), Schloss 2024, Schloss 2026, Yerke 2024.

### 6. Article `vignettes/articles/bdiv_guide.Rmd`

Question 2 (line 41) calls robust Aitchison "the most statistically
rigorous approach" and "your primary choice". Keep it, and add: "Schloss
(2026) disputes this for uneven sequencing depth: in a corrected re-run of
the original benchmark, rarefaction did better than robust Aitchison PCA.
Check whether your ordination's axes follow library size (see `rpca()`'s
help), and consider a rarefied analysis alongside."

### 7. `NEWS.md` (2.4.0)

One bullet: "Documentation: `rpca()`, `rpca_phylo()` and
`robust_aitchison()` discuss sequencing depth (Schloss 2026);
`dirmult_ttest()` recommends a `scale_sd` sensitivity check; `normalize()`,
`alr()`, `ilr()` and `philr()` note that log-ratio transforms don't improve
machine learning; `qmp_scale()` explains predicted loads."


## `log_base` on existing DA results

This item also adds one top-level attribute that the DA plots
([da_plots.md](da_plots.md)) read: `log_base`, the base of the log fold
changes (`2` for `dirmult_ttest()`'s `log2_fc`, `exp(1)` for `ancombc()` and
`ancombc2()`'s `log_fc`), so axis labels say "log2 fold change" or "log fold
change" without guessing from column names. Their `plot()` methods label
axes from it, and a test checks one label of each. It is a result
attribute under [README.md](README.md) § "Result attributes and plots". The
retrofit of `transform`, `da` and `log_ratio` step records planned here on
2026-11-20 was dropped with the records themselves.


## Checks

`R CMD check` (examples, codoc); both articles knit
(`pkgdown::build_article()`); `grep -rn -i "depth" R/rpca*.r vignettes/`
shows no claim that RPCA handles depth. No new tests: nothing computed
changes. Keep new examples fast (`ex_counts`, or `iris` as `rpca()`'s do).


## References

- Schloss, P. D. (2026). Rarefaction is better than robust Aitchison PCA and other compositional data analysis methods at controlling for uneven sequencing effort. *mSystems*, 11, e00930-26. doi:10.1128/msystems.00930-26
- Schloss, P. D. (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere*, 9, e00354-23. doi:10.1128/msphere.00354-23
- Martino, C., et al. (2019). A novel sparse compositional technique reveals microbial perturbations. *mSystems*, 4, e00016-19. doi:10.1128/mSystems.00016-19
- Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization. *Genome Biology*, 26, 139. doi:10.1186/s13059-025-03609-3
- McGovern, K. C., Nixon, M. P., & Silverman, J. D. (2023). Addressing erroneous scale assumptions in microbe and gene set enrichment analysis. *PLOS Computational Biology*, 19, e1011659. doi:10.1371/journal.pcbi.1011659
- Dos Santos, S. J., Murariu, A. C., Silverman, J. D., & Gloor, G. B. (2026). Not every gene is special: modelling scale controls the false discovery rate when analysing high-throughput sequencing data. *PLOS Computational Biology*, 22, e1014728. doi:10.1371/journal.pcbi.1014728
- Yerke, A., Fry Brumit, D., & Fodor, A. A. (2024). Proportion-based normalizations outperform compositional data transformations in machine learning applications. *Microbiome*, 12, 45. doi:10.1186/s40168-023-01747-z
- Nishijima, S., et al. (2025). Fecal microbial load is a major determinant of gut microbiome variation and a confounder for disease associations. *Cell*, 188, 222–236. doi:10.1016/j.cell.2024.10.022
