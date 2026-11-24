# Alpha Diversity

## Input Matrix

We will use the `ex_counts` dataset included with ecodive. This feature
table contains counts of bacterial genera across various samples.

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

Alpha diversity measures diversity within a single sample. In `ecodive`,
metrics are grouped into categories based on the aspect of diversity
they quantify.

### Richness Metrics

Richness metrics estimate the number of distinct features (e.g., genera)
in a sample. The simplest metric,
[`observed()`](https://cmmr.github.io/ecodive/reference/observed.md),
counts features with non-zero abundance.

``` r

# Equivalent to rowSums(counts > 0)
observed(counts)
#> Saliva   Gums   Nose  Stool 
#>      4      3      4      5 
```

The **Chao1** estimator extends this by inferring the number of
unobserved, low-abundance features based on the ratio of singletons
(`counts == 1`) to doubletons (`counts == 2`).

``` r

# Infers 3 unobserved genera
chao1(c(1, 1, 1, 1, 2, 5, 5, 5))
#> [1] 11

# Infers no unobserved genera
chao1(c(1, 2, 2, 2, 2, 5, 5, 5))
#> [1] 8

chao1(counts)
#> Saliva   Gums   Nose  Stool 
#>      4      3      4     11 
```

By default,
[`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md) is the
bias-corrected estimator, as in scikit-bio and QIIME 2.
`chao1(bias_corrected = FALSE)` gives the classic formula, which is
larger when there are few doubletons.

### Diversity Metrics

Diversity metrics account for both richness and evenness (how equally
abundances are distributed).

**Simpson’s index** is often used as a measure of evenness, representing
the probability that two randomly selected individuals belong to
different species.

``` r

# High Evenness (0.8) vs Low Evenness (0.07)
simpson(c(20, 20, 20, 20, 20))
#> [1] 0.8
simpson(c(100, 1, 1, 1, 1))
#> [1] 0.07507396

# Stool < Gums < Saliva < Nose
sort(simpson(counts))
#>      Stool       Gums     Saliva       Nose 
#> 0.02302037 0.18806133 0.50725478 0.63539593 
```

The **Shannon diversity index** (entropy) is another common metric that
weights both richness and evenness.

``` r

# High richness, High evenness
shannon(rep(100, 100))
#> [1] 4.60517

# Stool < Gums < Saliva < Nose
sort(shannon(counts))
#>      Stool       Gums     Saliva       Nose 
#> 0.07927797 0.35692121 0.74119910 1.10615349 
```

### Dominance Metrics

Dominance metrics focus on the abundance of the most common species. The
**Berger-Parker** index is the proportional abundance of the single most
abundant feature.

``` r

# Stool is dominated by Bacteroides (341/345 counts -> ~0.99)
# Nose is more balanced; Corynebacterium is max (171/345 counts -> ~0.49)
sort(berger(counts))
#>      Nose    Saliva      Gums     Stool 
#> 0.4956522 0.5217391 0.8956522 0.9884058 
```

### Evenness Metrics

Evenness metrics measure how equally the observed features share the
sample, apart from how many there are. **Pielou’s evenness** is
Shannon’s entropy divided by its largest possible value, \ln n for n
observed features; QIIME 2’s `core-metrics` reports it.

``` r

# Perfectly even (1) vs dominated by one genus (0.13)
pielou(c(20, 20, 20, 20, 20))
#> [1] 1
pielou(c(100, 1, 1, 1, 1))
#> [1] 0.1344212

# Stool < Gums < Saliva < Nose
sort(pielou(counts))
#>      Stool       Gums     Saliva       Nose 
#> 0.04925817 0.32488369 0.53466213 0.79792108 
```

[`heip()`](https://cmmr.github.io/ecodive/reference/heip.md),
[`simpson_evenness()`](https://cmmr.github.io/ecodive/reference/simpson_evenness.md)
and
[`mcintosh_evenness()`](https://cmmr.github.io/ecodive/reference/mcintosh_evenness.md)
are other evenness indices.

### Hill Numbers

The **Hill numbers** ^qD are the “effective number of species”: how many
equally abundant features would give the same diversity. The order q
sets how much weight abundant features get, and unifies the common
indices: q = 0 is richness, q = 1 is the exponential of Shannon’s
entropy, q = 2 is the inverse Simpson index, and q = \infty is the
inverse Berger-Parker index.

``` r

# A diversity profile: q = 0, 1, 2 and Inf
sapply(c(0, 1, 2, Inf), function (q) hill(counts, q = q))
#>        [,1]     [,2]     [,3]     [,4]
#> Saliva    4 2.098450 2.029446 1.916667
#> Gums      3 1.428923 1.231620 1.116505
#> Nose      4 3.022709 2.742701 2.017544
#> Stool     5 1.082505 1.023563 1.011730
```

[`renyi()`](https://cmmr.github.io/ecodive/reference/renyi.md) and
[`tsallis()`](https://cmmr.github.io/ecodive/reference/tsallis.md) are
the entropies of the same orders.

### Coverage

**Good’s coverage**, 1 - F_1 / X_T, estimates the fraction of the
community’s individuals that belong to features seen in the sample:
values near 1 mean that deeper sequencing would find few new features.

``` r

goods_coverage(counts)
#>    Saliva      Gums      Nose     Stool 
#> 0.9971014 1.0000000 1.0000000 0.9884058 
```

### Several Metrics at Once

Give
[`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
several metric names to get a matrix with a column for each.

``` r

alpha_div(counts, c('observed', 'shannon', 'pielou', 'goods_coverage'))
#>        observed    shannon     pielou goods_coverage
#> Saliva        4 0.74119910 0.53466213      0.9971014
#> Gums          3 0.35692121 0.32488369      1.0000000
#> Nose          4 1.10615349 0.79792108      1.0000000
#> Stool         5 0.07927797 0.04925817      0.9884058
```

### Phylogenetic Metrics

Phylogenetic metrics use a phylogenetic tree to incorporate evolutionary
distance. **Faith’s Phylogenetic Diversity (PD)** calculates the total
branch length spanned by the features present in a sample.

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

## Formulas

Given:

- n : Number of features (e.g. species, OTUs, ASVs).
- X_i : Integer count of the i-th feature.
- X_T : Total of all counts (sequencing depth). X_T = \sum\_{i=1}^{n}
  X_i
- P_i : Proportional abundance of the i-th feature. P_i = X_i / X_T
- F_1 : Number of singletons (X_i = 1).
- F_2 : Number of doubletons (X_i = 2).
- H : Shannon’s entropy, -\sum\_{i=1}^{n} P_i \ln P_i.
- q : The order of a Hill number or entropy (`q`, default 2).

[TABLE]

### Abundance-based Coverage Estimator (ACE)

Given:

- r : Rare cutoff (features with \le r counts are considered rare).
- F\_{rare} : Number of rare features.
- F\_{abund} : Number of abundant features (\> r counts).
- X\_{rare} : Total counts belonging to rare features.
- C\_{ace} : Sample abundance coverage estimator.
- \gamma\_{ace}^2 : Estimated coefficient of variation.

C\_{ace} = 1 - \frac{F_1}{X\_{rare}}

\gamma\_{ace}^2 = \max\left\[\frac{F\_{rare}
\sum\_{i=1}^{r}i(i-1)F_i}{C\_{ace}X\_{rare}(X\_{rare} - 1)} - 1,
0\right\]

D\_{ace} = F\_{abund} + \frac{F\_{rare}}{C\_{ace}} +
\frac{F_1}{C\_{ace}}\gamma\_{ace}^2

### Faith’s Phylogenetic Diversity (Faith’s PD)

Given n branches with lengths L and a binary vector A indicating
presence (1) or absence (0) of descendants on each branch:

\sum\_{i = 1}^{n} L_i A_i
