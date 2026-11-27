# ecodive

`ecodive` is an R package for ecological diversity analysis, from a
count table to figures and p-values: alpha and beta diversity,
ordination, and permutation tests, all in parallelized and
memory-efficient C. It is designed to handle large datasets, such as
those common in microbiome research.

## Why ecodive?

Analyzing ecological diversity is often a computational bottleneck,
especially with large datasets. `ecodive` addresses this by providing:

- **High Performance:** `ecodive` is written in C and parallelized using
  pthreads, making it dramatically faster than other R packages.
  [Benchmarks](https://cmmr.github.io/ecodive/articles/benchmark.html)
  show it can be up to **14,000x faster** and use up to **200,000x less
  memory**.

- **Zero Dependencies:** The package has no external R dependencies,
  making it lightweight, stable, and easy to install. This also makes it
  an ideal and secure backend for other R packages.

- **Comprehensive Metrics:** It implements 14 alpha and 38 beta
  diversity metrics, including classic and phylogenetic-aware methods
  like **Faith’s PD**, the **complete UniFrac family** and the
  compositional **PhILR distance**.

- **Ordination and Statistics:** Ten ordination methods (PCoA, NMDS,
  UMAP, t-SNE, db-RDA, robust Aitchison PCA, and more), biplots, and
  five permutation tests (PERMANOVA, PERMDISP, ANOSIM, MRPP, and
  Mantel), so a whole analysis needs no other package.

- **Trusted Results:** Every method is tested against the established
  package it replaces, such as `vegan`, `ape`, `scikit-bio`,
  `scikit-learn`, `uwot`, and `gemelli`. Results are identical for any
  number of threads, and random methods take a `seed`.

- **Ease of Use:** The API is simple and integrates seamlessly with
  popular bioinformatics packages like
  [`phyloseq`](http://joey711.github.io/phyloseq/) and
  [`rbiom`](https://cmmr.github.io/rbiom/).

## What’s Included

| Task | Wrapper | Functions |
|:---|:---|:---|
| Alpha diversity | [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md) | [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md), [`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md), [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md), [`faith()`](https://cmmr.github.io/ecodive/reference/faith.md), [`pielou()`](https://cmmr.github.io/ecodive/reference/pielou.md), [`hill()`](https://cmmr.github.io/ecodive/reference/hill.md), and 19 more |
| Beta diversity | [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md) | [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md), [`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md), [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md), [`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md), and 42 more |
| Ordination | [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) | [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md), [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md), [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md), [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md), [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md), [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md), [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md), [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md), [`capda()`](https://cmmr.github.io/ecodive/reference/capda.md) |
| Ordination of counts | [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md) | [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md), [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md), [`rpca_phylo()`](https://cmmr.github.io/ecodive/reference/rpca_phylo.md), [`dpcoa()`](https://cmmr.github.io/ecodive/reference/dpcoa.md), [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md), [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md), [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md), [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md), [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md), [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md) |
| Biplots |  | [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md), [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) |
| Comparing ordinations |  | [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md), [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md) |
| Model building |  | [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md), [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md), [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md), [`principal_response_curves()`](https://cmmr.github.io/ecodive/reference/principal_response_curves.md), [`mantel_regression()`](https://cmmr.github.io/ecodive/reference/mantel_regression.md) |
| Statistical tests | [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md) | [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md), [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md), [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md), [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md), [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md), [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md) |
| Interpreting differences |  | [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md) |
| Spatial |  | [`moran_i()`](https://cmmr.github.io/ecodive/reference/moran_i.md), [`dbmem()`](https://cmmr.github.io/ecodive/reference/dbmem.md), [`mantel_correlogram()`](https://cmmr.github.io/ecodive/reference/mantel_correlogram.md) |
| Differential abundance | [`da_test()`](https://cmmr.github.io/ecodive/reference/da_test.md) | [`ancom()`](https://cmmr.github.io/ecodive/reference/ancom.md), [`ancombc()`](https://cmmr.github.io/ecodive/reference/ancombc.md), [`linda()`](https://cmmr.github.io/ecodive/reference/linda.md), [`maaslin3()`](https://cmmr.github.io/ecodive/reference/maaslin3.md), [`dirmult_ttest()`](https://cmmr.github.io/ecodive/reference/dirmult_ttest.md), [`dirmult_lm()`](https://cmmr.github.io/ecodive/reference/dirmult_lm.md), [`lefse()`](https://cmmr.github.io/ecodive/reference/lefse.md), [`indicator_species()`](https://cmmr.github.io/ecodive/reference/indicator_species.md), [`da_wilcoxon()`](https://cmmr.github.io/ecodive/reference/da_wilcoxon.md), [`da_sensitivity()`](https://cmmr.github.io/ecodive/reference/da_sensitivity.md) |
| Compositional data | [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md) | [`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md), [`replace_zeros_gbm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md), [`replace_zeros_sq()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md), [`replace_zeros_bl()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md), [`replace_zeros_czm()`](https://cmmr.github.io/ecodive/reference/replace_zeros_cmult.md), [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md), [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md), [`philr()`](https://cmmr.github.io/ecodive/reference/philr.md), [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md), [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md), [`log_ratio_features()`](https://cmmr.github.io/ecodive/reference/log_ratio_features.md), [`variation_array()`](https://cmmr.github.io/ecodive/reference/variation_array.md), [`zero_pattern()`](https://cmmr.github.io/ecodive/reference/zero_pattern.md) |
| Feature associations | [`feature_assoc()`](https://cmmr.github.io/ecodive/reference/feature_assoc.md) | [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md), [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md), [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md), [`sparcc()`](https://cmmr.github.io/ecodive/reference/sparcc.md), [`assoc_pairs()`](https://cmmr.github.io/ecodive/reference/assoc_pairs.md), [`assoc_network()`](https://cmmr.github.io/ecodive/reference/assoc_network.md), [`cross_assoc()`](https://cmmr.github.io/ecodive/reference/cross_assoc.md) |
| Absolute abundance |  | [`qmp_scale()`](https://cmmr.github.io/ecodive/reference/qmp.md), [`qmp_rarefy()`](https://cmmr.github.io/ecodive/reference/qmp.md) |
| Community types |  | [`enterotype()`](https://cmmr.github.io/ecodive/reference/enterotype.md), [`valencia()`](https://cmmr.github.io/ecodive/reference/valencia.md), [`dmm()`](https://cmmr.github.io/ecodive/reference/dmm.md) |
| Source tracking |  | [`sourcetracker()`](https://cmmr.github.io/ecodive/reference/sourcetracker.md) |
| Quality control |  | [`decontam()`](https://cmmr.github.io/ecodive/reference/decontam.md) |
| Community assembly |  | [`neutral_model()`](https://cmmr.github.io/ecodive/reference/neutral_model.md) |
| Longitudinal |  | [`first_differences()`](https://cmmr.github.io/ecodive/reference/first_differences.md), [`first_distances()`](https://cmmr.github.io/ecodive/reference/first_distances.md), [`volatility()`](https://cmmr.github.io/ecodive/reference/volatility.md), [`nmit()`](https://cmmr.github.io/ecodive/reference/nmit.md), [`feature_volatility()`](https://cmmr.github.io/ecodive/reference/feature_volatility.md) |
| Batch correction |  | [`batch_mmuphin()`](https://cmmr.github.io/ecodive/reference/batch_mmuphin.md), [`batch_combat()`](https://cmmr.github.io/ecodive/reference/batch_combat.md) |
| Prediction |  | [`random_forest()`](https://cmmr.github.io/ecodive/reference/random_forest.md) |
| Plotting |  | [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of distance matrices, ordinations, tests, and the other results; [`composition_barplot()`](https://cmmr.github.io/ecodive/reference/composition_barplot.md), [`ternary_plot()`](https://cmmr.github.io/ecodive/reference/ternary_plot.md) |
| Data preparation |  | [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md), [`rarefaction_curve()`](https://cmmr.github.io/ecodive/reference/rarefaction_curve.md), [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md), [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md) |

Each wrapper runs its functions by name, such as
`beta_div(counts, 'bray')` or `dist_ord(dm, 'nmds')`.

## Installation

The latest stable version can be installed from CRAN.

``` r

install.packages('ecodive')
```

The development version is available on GitHub.

``` r

install.packages('pak')
pak::pak('cmmr/ecodive')
```

## Usage

`ecodive` functions are straightforward to use. Here are a few examples.

### With `phyloseq` or `rbiom` objects

The easiest way to use `ecodive` is with a `phyloseq`, `rbiom`, or
`TreeSummarizedExperiment` object. These objects conveniently bundle the
count data and phylogenetic tree.

``` r

library(ecodive)
data(esophagus, package = 'phyloseq')
hmp50 <- rbiom::hmp50

# Calculate weighted UniFrac distance
w_unifrac <- weighted_unifrac(esophagus)
w_unifrac
#> # A dist: 3 samples (lower triangle)
#> # Weighted UniFrac
#>       B     C
#> C 0.105
#> D 0.14  0.142
#> # Attributes: method

# Calculate Faith's Phylogenetic Diversity
faith_pd <- faith(hmp50)
faith_pd[1:4]
#>   HMP01   HMP02   HMP03   HMP04 
#> 6.22296 8.59432 8.93375 9.86597
```

### With basic R objects

You can also provide the count data and phylogenetic tree as separate
objects.

The `ex_counts` and `ex_tree` objects are included with `ecodive`.

``` r

## Example Data ----------------------
counts <- rarefy(ex_counts)
counts[,1:4]
#>        Streptococcus Bacteroides Corynebacterium Haemophilus
#> Saliva           162           2               0         180
#> Gums             309           2               0          34
#> Nose               6           0             171           0
#> Stool              1         341               1           1


## Alpha Diversity -------------------
shannon(counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.35692121 1.10615349 0.07927797

faith(counts, tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    155    101    202


## Beta Diversity --------------------
bray(counts)
#> # A dist: 4 samples (lower triangle)
#> # Bray-Curtis Dissimilarity
#>       Saliva  Gums  Nose
#> Gums   0.426
#> Nose   0.98  0.983
#> Stool  0.988 0.988 0.991
#> # Attributes: method

weighted_unifrac(counts, tree = ex_tree)
#> # A dist: 4 samples (lower triangle)
#> # Weighted UniFrac
#>       Saliva Gums Nose
#> Gums    37
#> Nose    67.2   56
#> Stool  110    109  110
#> # Attributes: method
```

### Ordination

Ordinations place samples on a few axes for plotting. Most start from
beta diversity distances;
[`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md),
[`lda()`](https://cmmr.github.io/ecodive/reference/lda.md),
[`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) and
[`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md) start
from the counts. With only four samples, `ex_counts` is too small to
ordinate, so these examples use R’s built-in `iris` data: 150 flowers
from three species.

``` r

dm <- bray(iris[, 1:4])

ord <- pcoa(dm)
ord
#> # A matrix: 150 samples × 3 axes
#> # Principal Coordinates Analysis of Bray-Curtis Dissimilarity
#> # % of variance: 89.0, 9.3, 2.9
#>    PCoA1   PCoA2    PCoA3
#>  1 0.169 -0.0259 -0.00505
#>  2 0.171  0.0235  0.0308
#>  3 0.187 -0.0032  0.0322
#>  4 0.18   0.0113  0.0372
#>  5 0.176 -0.03   -0.00613
#>  6 0.136 -0.0398 -0.0431
#>  7 0.186 -0.0182  0.0156
#>  8 0.167 -0.0165  0.00287
#>  9 0.194  0.0406  0.0413
#> 10 0.171  0.011   0.0286
#> # … 140 more samples
#> # Use print(x, n = ...) to see more samples
#> # Attributes: method, eig, prop, trace, correction, constant

# Any method by name, coloured by a column of metadata.
plot(dist_ord(dm, 'pcoa'), 'Species', iris)
plot(dist_ord(dm, 'nmds'), 'Species', iris)
plot(dist_ord(dm, 'umap'), 'Species', iris)
```

![PCoA, NMDS and UMAP ordinations of the iris data, colored by
species](reference/figures/README-ordination.png)

### Statistical Tests

Permutation tests ask whether a variable, such as species or treatment,
explains the distances between samples.

``` r

permanova(dm, 'Species', iris)
#> Permutational multivariate analysis of variance (PERMANOVA) of Bray-Curtis Dissimilarity
#> 
#> F = 532.7, p = 0.001, z = 577.6
#> Permutations: 999 (free)
#> 
#> Model R2 = 0.8788. Model tested as a whole:
#> 
#>           Df SumOfSqs     R2     F Pr(>F)    
#> Model      2   2.3173 0.8788 532.7  0.001 ***
#> Residual 147   0.3197 0.1212                 
#> Total    149   2.6370 1.0000                 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Any test by name.
dist_test(dm, 'anosim', 'Species', iris)
#> Analysis of similarities (ANOSIM) of Bray-Curtis Dissimilarity
#> 
#> R = 0.8576, p = 0.001, z =  79.4
#> Permutations: 999 (free)
#> 
#> Mean rank of the distances, between and within groups:
#>            pairs mean_rank
#> Between     7500   7163.90
#> setosa      1225   2239.66
#> versicolor  1225   2457.85
#> virginica   1225   2418.15
```

## Available Methods

Use
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)
to browse the metrics available for calculating diversity,
[`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
for ordination methods, and
[`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
for statistical tests.

``` r

# Alpha Diversity
list_metrics('alpha', 'id')
#>  [1] "ace"               "berger"            "brillouin"        
#>  [4] "chao1"             "faith"             "fisher"           
#>  [7] "gini"              "simpson"           "goods_coverage"   
#> [10] "heip"              "hill"              "inv_simpson"      
#> [13] "margalef"          "mcintosh"          "mcintosh_evenness"
#> [16] "menhinick"         "observed"          "pielou"           
#> [19] "renyi"             "shannon"           "dominance"        
#> [22] "simpson_evenness"  "squares"           "strong"           
#> [25] "tsallis"          

# Beta Diversity
list_metrics('beta', 'id')
#>  [1] "aitchison"                 "bhattacharyya"            
#>  [3] "bray"                      "canberra"                 
#>  [5] "chebyshev"                 "chord"                    
#>  [7] "clark"                     "sorensen"                 
#>  [9] "divergence"                "euclidean"                
#> [11] "generalized_unifrac"       "gower"                    
#> [13] "hamming"                   "hellinger"                
#> [15] "horn"                      "jaccard"                  
#> [17] "jensen"                    "jsd"                      
#> [19] "lorentzian"                "manhattan"                
#> [21] "matusita"                  "minkowski"                
#> [23] "morisita"                  "motyka"                   
#> [25] "normalized_unifrac"        "ochiai"                   
#> [27] "philr_distance"            "psym_chisq"               
#> [29] "robust_aitchison"          "soergel"                  
#> [31] "squared_chisq"             "squared_chord"            
#> [33] "squared_euclidean"         "topsoe"                   
#> [35] "unweighted_unifrac"        "variance_adjusted_unifrac"
#> [37] "wave_hedges"               "weighted_unifrac"

# Ordination methods that take distances
list_dist_ords('id')
#> [1] "capda"    "capscale" "dbrda"    "isomap"   "mds"      "nmds"     "pcoa"    
#> [8] "tsne"     "umap"

# Statistical tests
list_dist_tests('id')
#> [1] "anosim"         "mantel"         "mantel_partial" "mrpp"          
#> [5] "permanova"      "permdisp"
```

## Documentation

The online manual for ecodive is available at
<https://cmmr.github.io/ecodive/>. It includes reference pages for each
function and these guides, which are also available offline with
[`vignette()`](https://rdrr.io/r/utils/vignette.html):

- [Introduction to
  ecodive](https://cmmr.github.io/ecodive/articles/ecodive.html) -
  [`vignette('ecodive')`](https://cmmr.github.io/ecodive/articles/ecodive.md)
- [Alpha Diversity](https://cmmr.github.io/ecodive/articles/adiv.html) -
  [`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)
- [Beta Diversity](https://cmmr.github.io/ecodive/articles/bdiv.html) -
  [`vignette('bdiv')`](https://cmmr.github.io/ecodive/articles/bdiv.md)
- [Ordination](https://cmmr.github.io/ecodive/articles/ordination.html) -
  [`vignette('ordination')`](https://cmmr.github.io/ecodive/articles/ordination.md)
- [Statistical
  Tests](https://cmmr.github.io/ecodive/articles/stat_tests.html) -
  [`vignette('stat_tests')`](https://cmmr.github.io/ecodive/articles/stat_tests.md)
- [Performance](https://cmmr.github.io/ecodive/articles/performance.html) -
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)

## Community guidelines

### Support

Bug reports, feature requests, and general questions can be submitted at
<https://github.com/cmmr/ecodive/issues>.

### Contributing

Pull requests are welcome. Please ensure contributed code is covered by
tests and documentation (add additional tests and documentation as
needed) and passes all automated tests.

New functions must leverage C and pthreads to minimize memory and CPU
time.

Please note that the ecodive project is released with a [Contributor
Code of Conduct](https://cmmr.github.io/ecodive/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Citation

If you use `ecodive` in your research, please cite our paper:

Smith et al., (2026). ecodive: Parallel and Memory-Efficient R Package
for Ecological Diversity Analysis. *Journal of Open Source Software*,
11(121), 9777, <https://doi.org/10.21105/joss.09777>

## License

`ecodive` is licensed under the [MIT
License](https://opensource.org/license/mit) © 2026 ecodive authors.

This package also contains third-party code:

- [**PCG Random Number Generation for C**](https://www.pcg-random.org):
  Licensed under the [Apache License, Version
  2.0](https://opensource.org/license/apache-2.0).
- [**scikit-learn**](https://scikit-learn.org) (SMACOF for MDS and NMDS,
  t-SNE, and Isomap): Licensed under the [BSD 3-Clause
  License](https://opensource.org/license/bsd-3-clause).
- [**scikit-bio**](https://scikit.bio) (PCoA biplot scores): Licensed
  under the [BSD 3-Clause
  License](https://opensource.org/license/bsd-3-clause).
- [**gemelli**](https://github.com/biocore/gemelli) (robust Aitchison
  PCA): Licensed under the [BSD 3-Clause
  License](https://opensource.org/license/bsd-3-clause).
- [**uwot**](https://github.com/jlmelville/uwot) (UMAP): Licensed under
  the [BSD 2-Clause
  License](https://opensource.org/license/bsd-2-clause).

Detailed copyright and licensing information for all components can be
found in the
[`inst/COPYRIGHTS`](https://github.com/cmmr/ecodive/blob/main/inst/COPYRIGHTS)
file.
