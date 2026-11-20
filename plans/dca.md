# `dca()`: detrended correspondence analysis

**Priority:** 69. **Tier:** clean room, and **superseded**. vegan
`decorana()` (GPL-2) is output-only. **Build after `ca()`**
([cca_ca.md](cca_ca.md)). Read [README.md](README.md) first.


## Status and why it's near the bottom

DCA (Hill & Gauch 1980) "detrends" CA's second and later axes, to remove
the arch effect, and rescales axes to SD units of species turnover. Its
first-axis gradient length (in SD units) is a classic rule of thumb for
choosing between linear methods (PCA, RDA) and unimodal ones (CA, CCA).

But the method is widely criticised (Minchin 1987; Oksanen & Minchin
1997), and vegan's own documentation keeps downweighting and some
transforms "only for compliance". NMDS is now preferred. Some users still
expect it (phyloseq `ordinate("DCA")`, microeco), so implement it last,
and aim to match vegan.


## Licensing

- Hill's original DECORANA Fortran is public domain, according to vegan's
  `decorana.Rd`. But **no unmodified copy was found**.
- vegan's `src/decorana.f` includes Minchin's (1997) and Oksanen's (2001)
  changes, under GPL-2; cornpy (GPL) wraps it too.
- So: **no code to port**. Write it from Hill & Gauch (1980), Hill (1979)
  and Oksanen & Minchin (1997), and from vegan's *documentation* of its
  options. Never from `decorana.f`.
- "Faithful to Hill" and "matches vegan" are different targets: vegan's
  documentation says its preparatory steps "may differ somewhat from the
  original", and that bugs were fixed. **Aim to match vegan's outputs**,
  and document where Hill's original would differ.


## Design

- **Signature:** `dca(counts, margin = 1L, n_axes = 3L, n_segments = 26L, n_rescales = 4L, downweight = FALSE, cpus = n_cpus())`.
  - The defaults mirror vegan's `decorana(iweigh = 0, iresc = 4, ira = 0,
    mk = 26)`; `n_rescales` is vegan's `iresc`, an exact count of
    rescaling cycles.
  - `n_axes` defaults to 3, like every ordination whose axes have
    eigenvalues (see [parameter_names.md](parameter_names.md) §3.1).
    DCA traditionally reports 4 axes, and `n_axes = 4` still works.
  - `downweight` is vegan's `iweigh`. `ira = 1` (basic reciprocal
    averaging, without detrending) is just `ca()`, so leave that option
    out and point to `ca()`.
- **Result:** an `ecodive_ord` matrix with columns `DCA1..` (`DCA1..4`
  with `n_axes = 4`).
  - `eig` holds the "decorana values"; `prop` is `NA`, since DCA has no
    proper inertia decomposition.
  - Attributes `axis_lengths` (the gradient lengths in SD units, the most
    used output) and `loadings` (species scores).
- **Wrapper:** joins `COUNT_ORDS`, with aliases `dca`, `decorana` and
  `detrendedcorrespondenceanalysis`. Check `dca` against `dbrda`: no
  prefix clash.
- **Masking:** none (vegan's is `decorana`).


## Tests

- Fixtures from vegan 2.7-x `decorana()` on dune and varespec, with
  default and `iweigh = 1` settings.
  - Compare site and species scores (sign-matched), decorana values and
    axis lengths at 1e-6 or better.
  - Detrending by segments is iterative, so expect tolerances looser than
    1e-10. Document the achieved agreement.
- Results are identical for any `cpus`; a smoke-test line.


## Docs

- An Rd page with a "Prefer NMDS" note citing Minchin (1987).
- The method table in `vignettes/ordination.Rmd`, and NEWS.


## References

- Hill, M. O. (1979). *DECORANA — A FORTRAN program for detrended
  correspondence analysis and reciprocal averaging*. Cornell University.
- Hill, M. O., & Gauch, H. G. (1980). Detrended correspondence analysis:
  an improved ordination technique. *Vegetatio*, 42, 47–58.
- Minchin, P. R. (1987). An evaluation of the relative robustness of
  techniques for ecological ordination. *Vegetatio*, 69, 89–107.
- Oksanen, J., & Minchin, P. R. (1997). Instability of ordination results
  under changes in input data order: explanations and remedies. *Journal
  of Vegetation Science*, 8, 447–454.
