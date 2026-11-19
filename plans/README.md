# Plans for new ecodive methods

Each file in this folder is a brief for one new method, or a small set of
related methods, written so that a fresh Claude session can implement it
without the conversation that produced it. Start with:

1. **[TODO.md](TODO.md)**: what to build next, in priority order.
2. **This file**: the rules every method follows.
3. **The method's own plan**, linked from TODO.md.

[new_methods.md](new_methods.md) is the original survey and
decision record, kept for reference. [compositional.md](compositional.md)
is the same for the compositional-data methods (zeros, log ratios,
proportionality, differential abundance, QMP and the gemelli ports). The per-method plans hold everything
from it that a method needs. **Where they differ** (function names, in
particular), **the per-method plans and this README win.**

`plans/` is in `.Rbuildignore`, so none of this goes into the CRAN package.


## What ecodive promises (non-negotiable)

- **Same results as the established implementation.** Match the reference
  to rounding (relative difference around 1e-10) wherever the method is
  deterministic. For methods with random starts or stochastic steps, match
  the deterministic parts exactly, and compare final results on quality
  measures (stress, KL divergence, neighbours kept). Differences that amount
  only to a different random seed are fine.
- **Less CPU time and memory:** multithreaded C (pthreads) through
  `run_parallel()` in `src/parallel.c`, which uses a persistent thread pool.
  Keep sparse tables sparse (`src/ecomatrix.c` reads every supported table
  type), and never form an n × n matrix when a packed triangle or a
  matrix-free product will do (see `src/pcoa.c` and `src/lanczos.c`).
- **Deterministic:** the same inputs and `seed` give bit-identical results
  on every run and for **any `cpus`**.
  - Split work into pieces that don't depend on the thread count.
  - Sum in a fixed order; no shared floating-point accumulators.
  - Random parts take `seed = 0L` and use PCG (`src/pcg_basic.c`), never
    R's RNG. ecodive never touches `.Random.seed`.
- **No new dependencies,** not even in Suggests. Code under `tests/` never
  calls a reference package.


## Licensing (read before porting anything)

- ecodive is MIT. Port code **only** from permissively licensed sources:
  MIT, BSD-2/3, Apache-2.0, Artistic-2.0 or public domain. Each method plan
  names the source and its license.
- **GPL / LGPL / AGPL / CeCILL / non-commercial code is output-only:** run
  it to make test fixtures, but don't read its source to write ecodive's.
  This covers vegan, ade4, adespatial, ecodist, BiodiversityR, mixOmics,
  MiRKAT, CompQuadForm, phyloseq, destiny, ropls, pairwiseAdonis,
  RVAideMemoire, graphtools and PHATE. Where a plan says "clean room", work
  only from the cited papers and the reference's *documentation and
  outputs*.
- When you port code:
  - put its copyright and license text in the ported C file's header (as
    `src/nmds.c` does for scikit-learn);
  - add a `Files:` stanza to `inst/COPYRIGHTS`;
  - add a `ctb`/`cph` person to `Authors@R` in `DESCRIPTION` if the source
    isn't already credited.
- Apache-2.0 code keeps the Apache text and notices (as
  `src/pcg_basic.c` does); it can't be relabelled MIT.
- **Artistic-2.0 code** (allowed since 2026-11-19) keeps the Artistic text
  and its copyright notice in each ported file; it can't be relabelled MIT.
  Its section 4 also asks a modified version to:
  - **document how it differs** from the original: say in the file header
    that it is a modified port (e.g. "C port of <function> from <package>
    <version>; restructured for ecodive, see the method's plan"), and list
    user-visible differences in the roxygen docs' "Differences from ..."
    section;
  - **bear a different name** and not stop the original being installed
    (section 4(b)): ecodive is a separate package, so this holds.
- **Code derived from Artistic-2.0 code** but distributed under a
  permissive licence (e.g. scikit-bio's `ancombc()`, written from ANCOMBC's
  source and released as BSD-3) carries **both** notices: the permissive
  one it was released under, and the Artistic-2.0 one of the code it came
  from, with the "how it differs" line above. Then the file complies
  whichever licence is held to apply.


## Where things go

| What | Where |
|:-----|:------|
| R function and roxygen docs | `R/<name>.r` (one file per method or family) |
| C code | `src/<name>.c`; register `.Call` entry points in `src/init.c`; shared declarations in `src/ecodive.h` |
| Exports | `NAMESPACE` **by hand**: it is not roxygen-managed |
| Rd files | Generated: never edit `man/*.Rd` by hand. Run `roxygen2::roxygenise()` once, after the roxygen comments are final, and commit whatever it writes (see "Tooling") |
| Reference index | `pkgdown/_pkgdown.yml` (the Ordination or Hypothesis Tests section) |
| News | `NEWS.md`, under the version being developed |
| Fixture generators | `data-raw/expected/<name>.r`, plus `<name>.py` for Python references, whose text outputs are committed under `data-raw/expected/<name>/` |
| Fixtures | `tests/tinytests/expected/<name>.rds`, written with `saveRDS(version = 2)` and kept small |
| Tests | `tests/tinytests/test-<name>.r`, plus one line per new export in `test-smoke.r` |


## Internals worth reusing

- **Ordination results:** `ord_result()` in `R/ord_utils.r` builds the
  sample × axis matrix of class `ecodive_ord`. Its attributes start with
  `method`, `eig` and `prop` (`NA` where undefined); method-specific
  attributes follow. `ord_flip()` fixes axis signs; `ord_k()` and
  `ord_keep()` check `n_axes` and drop non-positive axes.
- **Test results:** `new_ecodive_test()` in `R/ecodive_test.r` builds the
  object.
  - Core fields come first: `statistic` (named), `p_value`, `z`, `n_perm`.
  - Permuted statistics go in attribute `<stat>_perm`.
  - Classes are `c('ecodive_<name>', 'ecodive_test')`, with a
    `test_details()` method for printing.
  - `test_groups()` and `test_group_values()` handle `groups` and `data`
    the same way for every test.
  - `perm_p()` and `perm_z()` compute p-values and z-scores.
- **Permutations:** `perm_matrix()` in `R/design.r` turns `permutations`
  (a count or a matrix) into a matrix, and `src/permute.c` (`perm_source_*`)
  generates permutations reproducibly, with `strata`. A supplied matrix
  makes p-values comparable exactly with vegan's.
- **Metadata and formulas:**
  - `model_design()` handles `groups` / `data` / `formula`, including
    `Condition()`.
  - `match_rows()` and `per_sample()` align metadata with the samples.
  - `condition_basis()` and `design_basis()` build QR bases.
  - `perm_ss()` computes permuted sums of squares.
- **Constrained ordination:** `R/dbrda.r` (`constrained_ord()`,
  `dbrda_core()`, `dbrda_r_squared()`, `dbrda_centroids()`) and
  `src/dbrda.c`.
- **Eigen and SVD:** block Lanczos in `src/lanczos.c` (`lanczos()`,
  `lz_*`), used by pcoa, pca, dbrda and umap. `gower_mult()` is a
  matrix-free Gower product.
- **Neighbours:** `knn_from_dist()` in `src/knn.c`.
- **Count tables:** `validate_counts()`, `validate_norm()`,
  `validate_pseudocount()` and `validate_margin()` in `R/validate.r`, and
  `normalize()` in C (norms: none, percent, chord, hellinger, binary, clr,
  rclr).
- **Group kernels:** `src/group_kernel.h`, the within- and between-group
  distance sums that mrpp, anosim and permdisp share.
- **The uwot SGD port:** `src/umap.c`, which follows uwot's BSD-2 C++
  headers.


## Naming and wrappers

- **Function names follow the scheme in the next section.**
- **Methods with several named variants get one function each,** as the
  UniFrac family does (`unweighted_unifrac()`, `weighted_unifrac()`, ...).
  Use an argument only for options that every reference also treats as an
  argument (e.g. `pcoa(correction)`, `permdisp(type)`, `mantel(method)`).
- **Zero replacements** join `replace_zeros()`: add a row to
  `ZERO_REPLACEMENTS` (see [zero_replacement.md](zero_replacement.md)).
  **Ordinations of distances** join `dist_ord()`: add a row to `ORDINATIONS`
  in `R/dist_ord.r`. **Ordinations of a count table** join `count_ord()`:
  add a row to `COUNT_ORDS` in `R/count_ord.r`. **Tests of `dm` against
  one variable** join `dist_test()`: add a row to `DIST_TESTS` in
  `R/stat_tests.r`. Methods that return tables or take two inputs are plain
  functions, with no wrapper.
- **Aliases:**
  - Matching tries an exact name or alias first (`match_name()` in
    `R/utils.r`), then prefixes.
  - Keep `cap` for `capscale`.
  - Check every new alias against all the registries with the snippet
    below.
  - If a wrapper's matcher should point to another wrapper, follow how
    `match_dist_ord()` and `match_count_ord()` do it.
- **Argument names:** reuse existing ones where the meaning is the same
  (the full glossary is [parameter_names.md](parameter_names.md) §1):
  - `n_axes` (the number of axes returned), `axes` (**which** axes of an
    input ordination to use, as in `envfit()`), `cpus`, `seed`,
    `permutations`, `strata`;
  - `groups`, `data`, `formula`;
  - `given` (what is partialled out), `space` (sample positions:
    coordinates, geographic distances or spatial weights);
  - `margin`, `norm`, `pseudocount`;
  - `max_iters` (a cap), `n_iters` (attribute: the number run), `tol`,
    `n_starts`.
- **Arguments must pass cleanly through the wrappers.** R matches named
  arguments to a function's arguments before `...`, exactly or by
  **prefix**, so the wrappers' own arguments shadow a method's:

  | Wrapper | Its arguments | A method it runs must not take |
  |:--------|:--------------|:-------------------------------|
  | `dist_ord()` | `dm`, `method` | an argument named `method`, or any prefix of `method` or `dm` (`m`, `me`, `d`, ...) |
  | `count_ord()` | `counts`, `method` | an argument named `method`, or any prefix of `method` or `counts` (`m`, `c`, `co`, ...) |
  | `dist_test()` | `dm`, `test`, `groups`, `data` | an argument named `test`, or any prefix of `test`, `groups` or `data` (`t`, `g`, `gr`, `d`, `da`, ...) |
  | `replace_zeros()` | `counts`, `method` | an argument named `method`, or any prefix of `method` or `counts` (`m`, `c`, `co`, ...) |

  `groups` and `data` themselves are fine in `dist_test()`: it passes
  them through on purpose. Tests may take a `method` argument (as
  `mantel()` does, for the correlation), because `dist_test()` names its
  selector `test`. Plain functions with no wrapper are free of this
  constraint. `tests/tinytests/test-wrapper_args.r` enforces the rule for
  every registered method, so adding a registry row whose function breaks
  it fails the tests.
- **Pairwise tests** take `pairwise = FALSE` and `p_adjust = 'holm'`
  (see [pairwise_tests.md](pairwise_tests.md)).
- **Masking:** many names will mask vegan, ade4 or MASS functions when both
  packages are attached. Say so in a "Same names in other packages"
  section, as `dist_ord()`'s docs do.

## Function names: the scheme (set 2026-11-18)

Every planned name follows these rules, so that names group sensibly and
can't be mistaken for one another:

1. **A base method is its standard lowercase acronym or short name,** as
   R users know it from vegan and the literature: `pcoa`, `nmds`, `rda`,
   `cca`, `ca`, `dca`, `simper`, `bioenv`, `varpart`.
2. **A variant or companion of a base method is `<base>_<modifier>`,** so
   that a family sorts together in the index and in autocompletion:
   - `pcoa_biplot`;
   - `mantel_partial`, `mantel_correlogram`, `mantel_regression`;
   - `procrustes_test`.
3. **Every discriminant (supervised) method ends in `da`:** `lda`,
   `plsda`, `oplsda`, `splsda`, `capda`.
4. **Spell a name out when its acronym would be mistaken for something
   else:**
   - `principal_response_curves`, not `prc`, which looks like `prcomp`
     and would prefix-match pca's `prcomp` alias;
   - `mantel_regression`, not `mrm`, which looks like `mrpp`;
   - `moran_i`, not `moran`, which could be read as `dbmem`'s Moran's
     eigenvector maps;
   - `procrustes_test`, not `protest`, which reads as "pro test".

   The familiar acronym becomes an alias in the registry, where the
   function joins one, and is named in the docs.
5. **Use snake_case for multi-word names**, as the existing
   `pcoa_biplot`, `diffusion_map`, `forward_select` and `*_unifrac` do.
   Never use dots.
6. **Don't reuse vegan's name when vegan's meaning differs.** vegan's
   `capscale` is db-RDA, so Anderson's discriminant CAP is `capda`, and
   the `cap` alias stays with `capscale`.

**All functions**, existing and planned. ✓ means it exists.

| Family | Functions |
|:-------|:----------|
| Wrappers | ✓ `dist_ord`, ✓ `count_ord`, ✓ `dist_test`, ✓ `replace_zeros` (+ `list_*` / `match_*`) |
| Unconstrained, distances | ✓ `pcoa`, ✓ `mds`, ✓ `nmds`, ✓ `isomap`, ✓ `tsne`, ✓ `umap`, `diffusion_map`, `pacmap`, `phate`, `largevis` |
| Unconstrained, counts | ✓ `pca`, ✓ `rpca`, ✓ `rpca_phylo`, ✓ `ca`, `dca`, `dpcoa` |
| Unconstrained, several tables or time series | `rpca_joint`, `ctf` |
| Constrained | ✓ `dbrda`, ✓ `capscale`, ✓ `rda`, ✓ `cca`, `principal_response_curves` |
| Discriminant (`*da`) | ✓ `lda`, ✓ `plsda`, ✓ `oplsda`, `splsda`, ✓ `capda` |
| Fitting onto ordinations | ✓ `envfit`, ✓ `pcoa_biplot` |
| Comparing ordinations | ✓ `procrustes`, ✓ `procrustes_test`, `coinertia` |
| Model building | ✓ `varpart`, ✓ `forward_select`, ✓ `bioenv` |
| Tests of `dm` against one variable (`dist_test`) | ✓ `permanova`, ✓ `permdisp`, ✓ `anosim`, ✓ `mrpp`, ✓ `mantel`, ✓ `mantel_partial`, `mirkat` |
| Other distance tests and models | `mantel_correlogram`, `mantel_regression` |
| Composition | ✓ `simper` |
| Beta metrics (`beta_div`) | `philr_distance` |
| Zero replacement (`replace_zeros`) | ✓ `replace_zeros_mult`, `replace_zeros_czm`, `replace_zeros_gbm`, `replace_zeros_sq`, `replace_zeros_bl` |
| Log-ratio transforms | ✓ `alr`, ✓ `ilr`, `philr` (clr and rclr stay `norm` choices) |
| Log ratios | ✓ `log_ratio`, ✓ `balance` |
| Proportionality | ✓ `vlr`, ✓ `rho_p`, ✓ `phi_s` |
| Differential abundance | ✓ `ancom`, ✓ `ancombc`, ✓ `ancombc2`, ✓ `ancombc_global`, ✓ `ancombc_pairwise`, ✓ `ancombc_dunnett`, ✓ `ancombc_trend`, ✓ `structural_zeros`, ✓ `dirmult_ttest` |
| Absolute abundance | ✓ `qmp_scale`, ✓ `qmp_rarefy` |
| Spatial | `moran_i`, `dbmem` |

Names still worth a word of warning in the docs:
- `lda` also means latent Dirichlet allocation (topic models) and LEfSe's
  "LDA score" in microbiome work.
- `cca` also means canonical *correlation* analysis (scikit-learn, Julia).
- `rpca` also means robust or randomized PCA elsewhere.
- `mds` is SMACOF, not classical MDS (that's `pcoa`).


To check a new alias (run with ecodive installed):

```r
for (reg in c('METRICS', 'ORDINATIONS', 'COUNT_ORDS', 'DIST_TESTS', 'ZERO_REPLACEMENTS'))
  print(grep('^newalias', names(ecodive:::exact_names(reg)), value = TRUE))
```


## Tests

- **Two tiers** (`tests/tinytest.r`):
  - At home and on CI (`NOT_CRAN=true`): every test file.
  - On CRAN: only `test-smoke.r`, with at most 2 threads, no comparisons
    against stored or BLAS-computed references, and only class, shape and
    finiteness checks.
- Every method needs:
  - fixtures from the reference (matched at about 1e-10 where
    deterministic);
  - edge cases (n = 2 or 3, ties, duplicates, zero rows or columns, a
    singleton group, NA handling, invalid arguments with their messages);
  - a check that results are `identical()` for `cpus` 1, 2 and 4 (more at
    home);
  - a smoke-test line.
- Run the suite against an `-O2` install, not `load_all()`:

  ```sh
  rm -f src/*.o src/*.dll src/*.so
  R CMD INSTALL --preclean --library=<scratch lib> .
  cd tests && NOT_CRAN=true R_LIBS=<scratch lib> Rscript tinytest.r   # at home
  cd tests && R_LIBS=<scratch lib> Rscript tinytest.r                 # CRAN mode
  ```

- `devtools::test()` runs nothing: ecodive uses tinytest.


## Tooling gotchas

- **roxygen (set by the user, 2026-11-19):** `man/*.Rd` is generated from
  the roxygen comments in `R/`, so don't spend time on it.
  - Never edit an Rd file by hand; change the roxygen comment instead.
  - Run `roxygen2::roxygenise()` (or `devtools::document()`) **once**, when
    the roxygen comments are final, and commit everything it writes under
    `man/` (and `DESCRIPTION`'s `RoxygenNote`, if it changes), as written: no
    reverting, no hand-curation, even if a different
    roxygen2 version reformats files you didn't touch.
  - Don't inspect or diff the Rd output beyond that. `R CMD check` (codoc,
    examples) is the check that the docs are right.
  - `man/` stays committed because the GitHub workflows (R CMD check,
    pkgdown) use it as is and never run `document()`.
- **Stale object files:** `load_all()` and `roxygenise()` leave `-O0`
  `src/*.o` files, which a later `R CMD INSTALL` silently reuses (benchmarks
  then run about 10× slow). Delete them before installing.
- **Line endings:** LF only. Check with `git ls-files --eol`.
- **Windows `Rscript -e`:** an argument that starts with a newline crashes
  R 4.6.1. Write a script file instead.
- **Ignored `.r` files:** `vignettes/articles/.gitignore` ignores `*.R`,
  which on Windows matches `.r` too; `git add -f` is needed there.


## Git

- Commit straight to `main`; don't create branches. (Exception: a worktree
  session under `.claude/worktrees/` commits to its own `claude/*` branch,
  and the parent session merges it.)
- Commit freely, with pathspecs, because the user also edits and stages in
  this checkout. Never push without asking.
- End each commit message with the attribution line that your session's
  system prompt gives.
