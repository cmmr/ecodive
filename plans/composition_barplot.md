# `composition_barplot()`: ordered, hierarchically shaded stacked bars

**Priority:** 61. **Tier:** port of an idea (microshades 1.13, MIT; its
colour logic was read, nothing is copied). Read [README.md](README.md)
first, and [compositional.md](compositional.md).


## What it is

Stacked bars of relative abundance are the most common microbiome figure,
and the most criticised. ecodive can offer a version whose design answers
the criticisms, plus a CoDa-consistent companion:

| Critique (source) | How this design answers it |
|:------------------|:---------------------------|
| Lengths of stacked segments are hard to compare across bars (Snowflake, Peeters et al. 2024) | Samples are **ordered by an Aitchison-geometry ordering** (PCoA axis 1 of `aitchison()` = clr PC1, or a clustering), so neighbouring bars are similar and gradients read as gradients; the most abundant taxon sits on the common baseline; the optional **clr heatmap** gives every taxon its own aligned row |
| "Other" hides rare taxa; colour is hard with many taxa (microshades, Dahl et al. 2022; Snowflake) | **Hierarchical shading** (microshades' idea): one hue per higher-level group, shades for its top members, the group's own "other" lightest; the subtitle says how much of the total "Other" holds |
| Averaging samples into group bars distorts (practitioner critiques: a MiSeq tutorial, a Biostars thread; not peer reviewed) | **Never averaged:** one bar per sample; groups are panels separated by a gap |
| Closure hides absolute load (Gloor et al. 2017, Fig. 1) | Docs say bars show proportions only; for loads use `qmp_scale()` and a plain `barplot()` |
| Stacked bars ignore phylogeny / DA (StackbarExtended, Cuisiniere & Santos 2024) | Out of scope; `tree_plot()` and the DA plots cover those |

No peer-reviewed paper critiques microbiome stacked bars on its own; the
critiques come from tool papers' motivations and tutorials (research note
Q9). The Rd should say so honestly.

**Uptake:** microshades 52 citations (2022); stacked bars are ubiquitous
but no count exists. Low effort; a better-ordered version of a figure users
make anyway. Not a compositional analysis, and the docs say so.


## Names

- `composition_barplot()`: draws immediately and returns the plotted values
  invisibly, as `graphics::barplot()` does (see
  [ternary_plot.md](ternary_plot.md) for the `_plot` rule); "composition"
  says what is stacked. No collisions in `NAMESPACE`; microbiome packages
  use `plot_bar()` (phyloseq), `plot_microshades()`, `comp_barplot()`
  (microViz, GPL): different names, no masking.
- Arguments from the glossary: `counts`, `groups`, `data`, `margin`,
  `pseudocount`, `clust` (the same values and meaning as in `beta_div()`:
  `DIST_CLUST`), `n_features`, `legend`, `col`, `cpus`. New:
  `feature_groups` (each feature's higher-level group, e.g. its phylum: the
  feature-side counterpart of `groups`), `n_feature_groups`, `heatmap`. Add
  them to [parameter_names.md](parameter_names.md) §1.


## Design

### Signature

```r
composition_barplot(counts, groups = NULL, data = NULL, feature_groups = NULL,
                    n_features = NULL, n_feature_groups = 5L, clust = 'pcoa',
                    heatmap = FALSE, margin = 1L, pseudocount = NULL,
                    legend = 'right', col = NULL, cpus = n_cpus(), ...)
```

- **Proportions:** each sample closed to 1 (`normalize(norm = 'percent')`
  / 100). Samples with total 0 are dropped with a warning.
- **Which features:**
  - without `feature_groups`: the `n_features` (default 12) with the
    largest mean proportion; the rest summed into "Other";
  - with `feature_groups` (a vector with one value per feature, matched by
    name, else by position; `NA` → "Unassigned"): the `n_feature_groups`
    groups with the largest total mean proportion, the others pooled as
    "Other"; within each kept group the `n_features` (default 4, as
    microshades' maximum) most abundant members, the rest of the group as
    "Other <group>". Ranking uses the mean proportion (microshades uses the
    sum, which ranks the same).
- **Sample order:** `clust = 'pcoa'` (default): order of PCoA axis 1 of
  `aitchison(counts, pseudocount = )`, i.e. the clr PC1; another
  `DIST_CLUST` method: `hclust` leaf order on the same distances; `'none'`
  or `FALSE`: input order. Use `dist_clust_order()`. With `groups`, samples
  are sorted by group, keeping that order within groups (the stable sort
  `plot.ecodive_dist()` uses), and a gap of one bar width separates groups,
  labelled beneath. `pseudocount` (shared rule) is used **only** for the
  ordering and the heatmap; the bars show raw proportions.
- **Stack order** (bottom to top): kept groups by decreasing abundance,
  within each its members by decreasing abundance then "Other <group>";
  "Other" on top. Without groups: features by decreasing abundance, "Other"
  on top.
- **Colours** (microshades' scheme, re-derived in HCL, not its hex table):
  - group \eqn{g} of \eqn{G}: hue \eqn{h_g} = `seq(15, 375, length.out = G
    + 1)[g]`; its members' shades `grDevices::hcl(h_g, c = 55, l =
    seq(35, 80, length.out = k + 1))`, darkest = most abundant, its "Other"
    the lightest; "Other" overall `'grey80'`;
  - without groups: `grDevices::hcl.colors(k, 'Dark 3')` (`plot_palette()`)
    plus `'grey80'`;
  - `col`: one colour per drawn segment, named or in stack order, overrides.
  - microshades also ships CVD-safe palettes; HCL ramps with distinct hues
    and a luminance ramp are reasonably CVD-robust; the docs suggest `col`
    for strict needs.
- **Drawing:** `graphics::barplot(P, space = s, border = NA, col = ...)`
  with `P` the segment × sample matrix in stack order, `s` 0.1 per bar (0
  above 60 samples) plus 1 at group starts; sample names on the axis when
  ≤ 60 samples (`las = 2`, `cex.axis = 0.6`), else none. Legend outside the
  plot on the right (`legend = 'right'`), below (`'bottom'`), or `NULL`;
  margin widened for it; entries top-down in stack order, grouped under
  bold group names. Subtitle: "Samples ordered by Aitchison PCoA axis 1;
  'Other' holds 23% of the total". `par()` / `layout()` restored on exit.
- **`heatmap = TRUE`:** `layout(matrix(1:2), heights = c(2, 1.4))`; below
  the bars, `image()` of the **centred clr** (clr value minus the feature's
  mean over samples) of the drawn features (not the "Other" segments), rows
  in legend order, columns at the same x positions (gaps kept as `NA`
  columns), `hcl.colors(65, 'Blue-Red 3')` symmetric about 0, feature names
  on the left, a small key on the right (`dist_plot_key()`). The clr uses
  all features (the pseudocount rule); the docs say it is the CoDa-consistent
  reading of the same samples.
- **Value:** invisibly, `list(proportions = P, col = , order = )`, with
  `P`'s rows the segments and columns the samples as drawn.


## Reference implementations

- **microshades 1.13** (github.com/KarstensLab/microshades; `LICENSE.md`:
  MIT, "Copyright (c) 2021 Erin Dahl"; `DESCRIPTION`: `MIT + file LICENSE`;
  verified 2026-11-20). Read: `R/color_mapping_functions.R`
  (`create_color_dfs()`: selected groups + "Other"; top
  `top_n_subgroups ≤ 4` per group by summed abundance; the subgroup "Other"
  forced to the lightest shade; "Other" group grey) and `R/palettes.R`
  (six 5-shade palettes, light to dark). Nothing is copied (the colours are
  re-derived), so no notice or `inst/COPYRIGHTS` stanza is needed; cite the
  paper. If the hex palettes are ever copied, add the MIT notice and a
  `Files:` stanza.
- The ordering and the clr heatmap are ecodive's own (`aitchison()`,
  `pcoa()`, `normalize()`).


## Tests (`test-composition_barplot.r`)

- Returned `proportions`: columns sum to 1 at 1e-15; "Other" = 1 − the
  drawn features; the top features are the largest means (ties by feature
  order); with `feature_groups`, groups and members chosen as specified,
  including `NA` → "Unassigned".
- Order: `clust = 'pcoa'` equals `order(pcoa(aitchison(x))[, 1])` (same
  sign convention as `dist_clust_order()`); `'none'` keeps input order;
  `groups` sorting is stable.
- Colours: one per segment; shades within a group darken with abundance;
  `col` by name.
- `pdf(NULL)`: returns invisibly; with `heatmap = TRUE`, `groups`,
  `legend = NULL`, 200 samples (no names drawn); `par('mfrow')` restored;
  `.Random.seed` untouched.
- Errors: `n_features = 0`; `feature_groups` of the wrong length; unknown
  `clust`; all-zero table.
- Smoke line.


## Docs

- Rd in "Plotting" → "Compositional displays". A short section "Why stacked
  bars mislead, and what this one does" with the table above (citations as
  given, and the note that the critiques are from tool papers and
  tutorials). Point to `variation_array()`, the covariance biplot and
  `ternary_plot()` for compositional readings. NEWS.


## References

- Dahl, E. M., Neer, E., Bowie, K. R., Leung, E. T., & Karstens, L. (2022).
  microshades: an R package for improving color accessibility and
  organization of microbiome data. *Microbiology Resource Announcements*,
  11, e00795-22. \doi{10.1128/mra.00795-22}
- Peeters, J., Bot, D. M., Rovelo Ruiz, G., & Aerts, J. (2024). Snowflake:
  visualizing microbiome abundance tables as multivariate bipartite graphs.
  *Frontiers in Bioinformatics*, 4, 1331043. \doi{10.3389/fbinf.2024.1331043}
- Cuisiniere, T., & Santos, M. M. (2024). StackbarExtended: a user-friendly
  stacked bar-plot representation incorporating phylogenetic information
  and microbiota differential abundance analysis. *F1000Research*.
  \doi{10.12688/f1000research.151662.1}
- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J.
  (2017). Microbiome datasets are compositional: and this is not optional.
  *Frontiers in Microbiology*, 8, 2224. \doi{10.3389/fmicb.2017.02224}
