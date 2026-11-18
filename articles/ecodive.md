# Introduction to ecodive

## Introduction

Ecodive covers an ecological diversity analysis from start to finish:

- **Alpha diversity** metrics describe a single sample’s diversity.
- **Beta diversity** metrics indicate how different a pair of samples
  are from each other.
- **Ordinations** place samples on a few axes, so that similar samples
  plot near each other.
- **Statistical tests** check whether a variable, such as treatment or
  body site, explains the differences between samples.

In this guide, we’ll use the `ex_counts` dataset included with ecodive.
`ex_counts` is a feature table that enumerates how many times each
bacterial genera was observed on different body sites.

``` r

library(ecodive)

t(ex_counts)
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  793   22     1
#> Bacteroides            2    4    2   611
#> Corynebacterium        0    0  498     1
#> Haemophilus          180   87    2     1
#> Propionibacterium      1    1  251     0
#> Staphylococcus         0    1  236     1
```

In this example, the ‘features’ in our feature table are genera.
However, your own dataset can use whatever feature makes sense -
species, OTUs, ASVs, or even something completely unrelated to ecology.

## Alpha Diversity

Alpha diversity metrics describe how many different genera are present
in a sample. Depending on the metric, this can take into account the
number of unique genera (richness), how evenly the population is split
among genera (evenness), or how distantly related the genera are
(phylogenetic diversity).

The available alpha diversity metrics can be listed using
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md).

``` r

list_metrics('alpha')[,1:5]
#>                                        name          id phylo weighted int_only
#> 1  Abundance-based Coverage Estimator (ACE)         ace FALSE     TRUE     TRUE
#> 3                       Berger-Parker Index      berger FALSE     TRUE    FALSE
#> 6                           Brillouin Index   brillouin FALSE     TRUE     TRUE
#> 8                                     Chao1       chao1 FALSE     TRUE     TRUE
#> 15           Faith's Phylogenetic Diversity       faith  TRUE    FALSE    FALSE
#> 16                           Fisher's Alpha      fisher FALSE     TRUE     TRUE
#> 18                       Gini-Simpson Index     simpson FALSE     TRUE    FALSE
#> 23                    Inverse Simpson Index inv_simpson FALSE     TRUE    FALSE
#> 29                Margalef's Richness Index    margalef FALSE     TRUE     TRUE
#> 31                           McIntosh Index    mcintosh FALSE     TRUE     TRUE
#> 32               Menhinick's Richness Index   menhinick FALSE     TRUE     TRUE
#> 37                        Observed Features    observed FALSE    FALSE    FALSE
#> 40                  Shannon Diversity Index     shannon FALSE     TRUE    FALSE
#> 45               Squares Richness Estimator     squares FALSE     TRUE     TRUE
```

- Further reading:
  [`vignette('adiv')`](https://cmmr.github.io/ecodive/articles/adiv.md)

## Beta Diversity

Beta diversity metrics describe how different two samples are, based on
the genera observed in each. UniFrac metrics incorporate a phylogenetic
tree into this calculation.

The available beta diversity metrics can be listed using
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md).

``` r

list_metrics('beta')[,2:6]
#>                           id phylo weighted int_only true_metric
#> 2                  aitchison FALSE     TRUE    FALSE        TRUE
#> 4              bhattacharyya FALSE     TRUE    FALSE        TRUE
#> 5                       bray FALSE     TRUE    FALSE       FALSE
#> 7                   canberra FALSE     TRUE    FALSE        TRUE
#> 9                  chebyshev FALSE     TRUE    FALSE        TRUE
#> 10                     chord FALSE     TRUE    FALSE        TRUE
#> 11                     clark FALSE     TRUE    FALSE        TRUE
#> 12                  sorensen FALSE    FALSE    FALSE       FALSE
#> 13                divergence FALSE     TRUE    FALSE        TRUE
#> 14                 euclidean FALSE     TRUE    FALSE        TRUE
#> 17       generalized_unifrac  TRUE     TRUE    FALSE        TRUE
#> 19                     gower FALSE     TRUE    FALSE        TRUE
#> 20                   hamming FALSE    FALSE    FALSE        TRUE
#> 21                 hellinger FALSE     TRUE    FALSE        TRUE
#> 22                      horn FALSE     TRUE    FALSE       FALSE
#> 24                   jaccard FALSE    FALSE    FALSE        TRUE
#> 25                    jensen FALSE     TRUE    FALSE        TRUE
#> 26                       jsd FALSE     TRUE    FALSE        TRUE
#> 27                lorentzian FALSE     TRUE    FALSE       FALSE
#> 28                 manhattan FALSE     TRUE    FALSE        TRUE
#> 30                  matusita FALSE     TRUE    FALSE        TRUE
#> 33                 minkowski FALSE     TRUE    FALSE        TRUE
#> 34                  morisita FALSE     TRUE     TRUE       FALSE
#> 35                    motyka FALSE     TRUE    FALSE       FALSE
#> 36        normalized_unifrac  TRUE     TRUE    FALSE        TRUE
#> 38                    ochiai FALSE    FALSE    FALSE       FALSE
#> 39                psym_chisq FALSE     TRUE    FALSE       FALSE
#> 40          robust_aitchison FALSE     TRUE    FALSE        TRUE
#> 42                   soergel FALSE     TRUE    FALSE        TRUE
#> 43             squared_chisq FALSE     TRUE    FALSE       FALSE
#> 44             squared_chord FALSE     TRUE    FALSE       FALSE
#> 45         squared_euclidean FALSE     TRUE    FALSE       FALSE
#> 47                    topsoe FALSE     TRUE    FALSE        TRUE
#> 48        unweighted_unifrac  TRUE    FALSE    FALSE        TRUE
#> 49 variance_adjusted_unifrac  TRUE     TRUE    FALSE        TRUE
#> 50               wave_hedges FALSE     TRUE    FALSE       FALSE
#> 51          weighted_unifrac  TRUE     TRUE    FALSE        TRUE
```

## Example

### Rarefaction

The `ex_counts` feature table has 345 saliva observations, but nose has
1011 observations. This unequal sampling depth can cause systematic
biases. Specifically, rare genera will be observed more often in samples
with greater sampling depths, thereby artificially inflating the
observed richness.

The first step then is to rarefy `ex_counts` so that all samples have
the same number of observations. Rarefying randomly removes observations
from samples with more observations.

``` r

rowSums(ex_counts)
#> Saliva   Gums   Nose  Stool 
#>    345    886   1011    615 

counts <- rarefy(ex_counts)

rowSums(counts)
#> Saliva   Gums   Nose  Stool 
#>    345    345    345    345 

t(counts)
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  309    6     1
#> Bacteroides            2    2    0   341
#> Corynebacterium        0    0  171     1
#> Haemophilus          180   34    0     1
#> Propionibacterium      1    0   82     0
#> Staphylococcus         0    0   86     1
```

### Classic Metrics

These alpha and beta diversity metrics have been around for 50+ years
and don’t require a phylogenetic tree. The beta diversity functions can
take a `weighted = FALSE` argument to use only presence/absence
information instead of relative abundances.

``` r

## Alpha Diversity -------------------

shannon(counts)
#>     Saliva       Gums       Nose      Stool 
#> 0.74119910 0.35692121 1.10615349 0.07927797 


## Beta Diversity --------------------

bray(counts)
#>          Saliva      Gums      Nose
#> Gums  0.4260870                    
#> Nose  0.9797101 0.9826087          
#> Stool 0.9884058 0.9884058 0.9913043
```

### Phylogenetic Metrics

A phylogenetic tree enables alpha and beta diversity metrics to take
into account evolutionary relatedness between the observed genera,
generally giving higher diversity values for samples with more distantly
related genera. Faith (for alpha diversity) and UniFrac (for beta
diversity) are examples of phylogenetic metrics.

The `ex_tree` object included with ecodive provides the phylogenetic
tree for the genera in `ex_counts`. For your own datasets, you can use
ecodive’s
[`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
function to import a phylogenetic tree from a newick formatted string or
file.

``` r

## Alpha Diversity -------------------

faith(counts, tree = ex_tree)
#> Saliva   Gums   Nose  Stool 
#>    180    155    101    202 


## Beta Diversity --------------------

normalized_unifrac(counts, tree = ex_tree)
#>          Saliva      Gums      Nose
#> Gums  0.4328662                    
#> Nose  0.7928701 0.6767840          
#> Stool 0.9677535 0.9829736 0.9936121
```

### Distance Matrices

Beta diversity functions return a `dist` object. You can convert this to
a standard R matrix with the
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) function.

``` r

dm <- bray(counts)
dm
#>          Saliva      Gums      Nose
#> Gums  0.1428571                    
#> Nose  0.5000000 0.7142857          
#> Stool 0.3333333 0.2500000 0.3333333

mtx <- as.matrix(dm)
mtx
#>           Saliva      Gums      Nose     Stool
#> Saliva 0.0000000 0.1428571 0.5000000 0.3333333
#> Gums   0.1428571 0.0000000 0.7142857 0.2500000
#> Nose   0.5000000 0.7142857 0.0000000 0.3333333
#> Stool  0.3333333 0.2500000 0.3333333 0.0000000

mtx['Saliva', 'Nose']
#> [1] 0.5
```

## Ordination

An ordination turns a distance matrix into coordinates for a plot.
ecodive has eight methods that take distances, run by name with
[`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) or
directly by their own functions, and six
([`pca()`](https://cmmr.github.io/ecodive/reference/pca.md),
[`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md),
[`rda()`](https://cmmr.github.io/ecodive/reference/rda.md),
[`cca()`](https://cmmr.github.io/ecodive/reference/cca.md),
[`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) and
[`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md)) that
take the count table itself, run by name with
[`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md).

| Function | Method |
|:---|:---|
| [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) | Principal Coordinates Analysis (PCoA) |
| [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) | Metric Multidimensional Scaling |
| [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) | Non-metric Multidimensional Scaling (NMDS) |
| [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) | Isometric Feature Mapping (Isomap) |
| [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) | t-distributed Stochastic Neighbor Embedding |
| [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md) | Uniform Manifold Approximation and Projection |
| [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) | Distance-based Redundancy Analysis (db-RDA) |
| [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md) | Constrained Analysis of Principal Coordinates |
| [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) | Principal Component Analysis (PCA) |
| [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) | Robust Aitchison PCA |
| [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) | Redundancy Analysis (RDA) |
| [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) | Canonical Correspondence Analysis (CCA) |
| [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) | Correspondence Analysis (CA) |
| [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) | Partial Least Squares Discriminant Analysis |

[`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md)
and [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) add
arrows for the features or metadata that drive the axes.
[`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md)
and
[`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
compare two ordinations of the same samples.

Four samples are too few to ordinate, so this example uses R’s built-in
`iris` data: 150 flowers, four measurements each, and three species.
Every ordination returns a matrix with one row per sample and one column
per axis.

``` r

dm  <- bray(iris[, 1:4])
ord <- pcoa(dm)

head(ord, 3)
#>          PCoA1       PCoA2        PCoA3
#> [1,] 0.1692812 -0.02591990 -0.005053611
#> [2,] 0.1706354  0.02346821  0.030848760
#> [3,] 0.1870090 -0.00320247  0.032221661

# Percent of variance on each axis.
round(100 * attr(ord, 'prop'), 1)
#> [1] 89.0  9.3  2.9

plot(ord, col = iris$Species)
```

- Further reading:
  [`vignette('ordination')`](https://cmmr.github.io/ecodive/articles/ordination.md)

## Statistical Tests

ecodive’s permutation tests take the same distance matrix, and a
variable to test: a vector with one value per sample, or a column of a
data frame.

| Function | Test |
|:---|:---|
| [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md) | Permutational Multivariate Analysis of Variance |
| [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md) | Homogeneity of Multivariate Dispersions |
| [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md) | Analysis of Similarities |
| [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) | Multi-Response Permutation Procedure |
| [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) | Mantel Test |
| [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md) | Partial Mantel Test |

[`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
runs any of them by name.

``` r

permanova(dm, 'Species', iris)
#> Permutational multivariate analysis of variance (PERMANOVA)
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
```

- Further reading:
  [`vignette('stat_tests')`](https://cmmr.github.io/ecodive/articles/stat_tests.md)

## Performance with Large Datasets

`ecodive` is designed to be fast out-of-the-box. However, if you are
working with very large datasets (e.g., thousands of samples), you can
gain additional performance by following the recommendations in our
performance guide.

- Further reading:
  [`vignette('performance')`](https://cmmr.github.io/ecodive/articles/performance.md)
