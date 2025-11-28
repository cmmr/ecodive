# Alpha Diversity

## Input Matrix

Here we’ll use the `ex_counts` feature table included with ecodive. It
contains the number of observations of each bacterial genera in each
sample. In the text below, you can substitute the word ‘genera’ for the
feature of interest in your own data.

``` r
library(ecodive)

counts <- rarefy(ex_counts)
t(counts)
#>                   Saliva Gums Nose Stool
#> Streptococcus        162  309    6     1
#> Bacteroides            2    2    0   341
#> Corynebacterium        0    0  171     1
#> Haemophilus          180   34    0     1
#> Propionibacterium      1    0   82     0
#> Staphylococcus         0    0   86     1
```

## Alpha Diversity

Alpha diversity is a measure of diversity within a single sample.

Depending on the metric, it may measure **richness** and/or
**evenness**.

### Richness

Richness is how many genera are present in a sample. The simplest metric
is to count the non-zero genera. You can do this with base R’s
[`rowSums()`](https://rdrr.io/r/base/colSums.html) or with ecodive’s
[`observed()`](https://cmmr.github.io/ecodive/reference/adiv_functions.md).

``` r
rowSums(counts > 0)
#> Saliva   Gums   Nose  Stool 
#>      4      3      4      5 

observed(counts)
#> Saliva   Gums   Nose  Stool 
#>      4      3      4      5 
```

The Chao1 metric takes this a step further by including unobserved low
abundance genera, inferred using the number of times `counts == 1` vs
`counts == 2`.

``` r
# Infers 8 unobserved genera
chao1(c(1, 1, 1, 1, 2, 5, 5, 5))
#> [1] 16

# Infers less than 1 unobserved genera
chao1(c(1, 2, 2, 2, 2, 5, 5, 5))
#> [1] 8.125

# Datasets without 1s and 2s give Inf or NaN
chao1(counts)
#> Saliva   Gums   Nose  Stool 
#>    4.5    3.0    NaN    Inf 
```

### Evenness

Evenness is how equally distributed genera are within a sample. The
Simpson metric is a good measure of evenness.

``` r
# High Evenness
simpson(c(20, 20, 20, 20, 20))
#> [1] 0.8

# Low Evenness
simpson(c(100, 1, 1, 1, 1))
#> [1] 0.07507396

# Stool < Gums < Saliva < Nose
sort(simpson(counts))
#>      Stool       Gums     Saliva       Nose 
#> 0.02302037 0.18806133 0.50725478 0.63539593 
```

### Richness and Evenness

The Shannon diversity index weights both richness and evenness.

``` r
# Low richness, Low evenness
shannon(c(1, 1, 100))
#> [1] 0.1101001

# Low richness, High evenness
shannon(c(100, 100, 100))
#> [1] 1.098612

# High richness, Low evenness
shannon(1:100)
#> [1] 4.416898

# High richness, High evenness
shannon(rep(100, 100))
#> [1] 4.60517

# Stool < Gums < Saliva < Nose
sort(shannon(counts))
#>      Stool       Gums     Saliva       Nose 
#> 0.07927797 0.35692121 0.74119910 1.10615349 
```

### Phylogenetic Alpha Diversity

Faith’s phylogenetic diversity index incorporates a phylogenetic tree of
the genera in order to measure how many of the tree’s branches are
represented by each sample.

``` r
# ex_tree:
#
#       +----------44---------- Haemophilus
#   +-2-|
#   |   +----------------68---------------- Bacteroides  
#   |                      
#   |             +---18---- Streptococcus
#   |      +--12--|       
#   |      |      +--11-- Staphylococcus
#   +--11--|              
#          |      +-----24----- Corynebacterium
#          +--12--|
#                 +--13-- Propionibacterium


faith(c(Propionibacterium = 1, Corynebacterium = 1), tree = ex_tree)
#> [1] 60

faith(c(Propionibacterium = 1, Haemophilus = 1), tree = ex_tree)
#> [1] 82

# Nose < Gums < Saliva < Stool
sort(faith(counts, tree = ex_tree))
#>   Nose   Gums Saliva  Stool 
#>    101    155    180    202 
```
