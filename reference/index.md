# Package index

## Alpha Diversity

Calculate within-sample diversity using various indices.

- [`alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.md)
  : Alpha Diversity Wrapper Function

### Richness metrics

Estimators focused on the number of species (richness).

- [`ace()`](https://cmmr.github.io/ecodive/reference/ace.md) :
  Abundance-based Coverage Estimator (ACE)
- [`chao1()`](https://cmmr.github.io/ecodive/reference/chao1.md) : Chao1
  Richness Estimator
- [`margalef()`](https://cmmr.github.io/ecodive/reference/margalef.md) :
  Margalef's Richness Index
- [`menhinick()`](https://cmmr.github.io/ecodive/reference/menhinick.md)
  : Menhinick's Richness Index
- [`observed()`](https://cmmr.github.io/ecodive/reference/observed.md) :
  Observed Features
- [`squares()`](https://cmmr.github.io/ecodive/reference/squares.md) :
  Squares Richness Estimator

### Diversity metrics

Indices accounting for both richness and evenness.

- [`brillouin()`](https://cmmr.github.io/ecodive/reference/brillouin.md)
  : Brillouin Index
- [`fisher()`](https://cmmr.github.io/ecodive/reference/fisher.md) :
  Fisher's Alpha
- [`inv_simpson()`](https://cmmr.github.io/ecodive/reference/inv_simpson.md)
  : Inverse Simpson Index
- [`shannon()`](https://cmmr.github.io/ecodive/reference/shannon.md) :
  Shannon Diversity Index
- [`simpson()`](https://cmmr.github.io/ecodive/reference/simpson.md) :
  Gini-Simpson Index

### Dominance metrics

Indices focused on the abundance of the most dominant species.

- [`berger()`](https://cmmr.github.io/ecodive/reference/berger.md) :
  Berger-Parker Index
- [`mcintosh()`](https://cmmr.github.io/ecodive/reference/mcintosh.md) :
  McIntosh Index

### Phylogenetic metrics

Diversity based on phylogenetic branch lengths.

- [`faith()`](https://cmmr.github.io/ecodive/reference/faith.md) :
  Faith's Phylogenetic Diversity (PD)

## Beta Diversity

Calculate between-sample similarity or distance.

- [`beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.md) :
  Beta Diversity Wrapper Function

### Abundance metrics

Distance measures using quantitative abundance data.

- [`aitchison()`](https://cmmr.github.io/ecodive/reference/aitchison.md)
  : Aitchison distance
- [`bhattacharyya()`](https://cmmr.github.io/ecodive/reference/bhattacharyya.md)
  : Bhattacharyya distance
- [`bray()`](https://cmmr.github.io/ecodive/reference/bray.md) :
  Bray-Curtis dissimilarity
- [`canberra()`](https://cmmr.github.io/ecodive/reference/canberra.md) :
  Canberra distance
- [`chebyshev()`](https://cmmr.github.io/ecodive/reference/chebyshev.md)
  : Chebyshev distance
- [`chord()`](https://cmmr.github.io/ecodive/reference/chord.md) : Chord
  distance
- [`clark()`](https://cmmr.github.io/ecodive/reference/clark.md) :
  Clark's divergence distance
- [`divergence()`](https://cmmr.github.io/ecodive/reference/divergence.md)
  : Divergence
- [`euclidean()`](https://cmmr.github.io/ecodive/reference/euclidean.md)
  : Euclidean distance
- [`gower()`](https://cmmr.github.io/ecodive/reference/gower.md) : Gower
  distance
- [`hellinger()`](https://cmmr.github.io/ecodive/reference/hellinger.md)
  : Hellinger distance
- [`horn()`](https://cmmr.github.io/ecodive/reference/horn.md) :
  Horn-Morisita dissimilarity
- [`jensen()`](https://cmmr.github.io/ecodive/reference/jensen.md) :
  Jensen-Shannon distance
- [`jsd()`](https://cmmr.github.io/ecodive/reference/jsd.md) :
  Jensen-Shannon divergence (JSD)
- [`lorentzian()`](https://cmmr.github.io/ecodive/reference/lorentzian.md)
  : Lorentzian distance
- [`manhattan()`](https://cmmr.github.io/ecodive/reference/manhattan.md)
  : Manhattan distance
- [`matusita()`](https://cmmr.github.io/ecodive/reference/matusita.md) :
  Matusita distance
- [`minkowski()`](https://cmmr.github.io/ecodive/reference/minkowski.md)
  : Minkowski distance
- [`morisita()`](https://cmmr.github.io/ecodive/reference/morisita.md) :
  Morisita dissimilarity
- [`motyka()`](https://cmmr.github.io/ecodive/reference/motyka.md) :
  Motyka dissimilarity
- [`psym_chisq()`](https://cmmr.github.io/ecodive/reference/psym_chisq.md)
  : Probabilistic Symmetric Chi-Squared distance
- [`soergel()`](https://cmmr.github.io/ecodive/reference/soergel.md) :
  Soergel distance
- [`squared_chisq()`](https://cmmr.github.io/ecodive/reference/squared_chisq.md)
  : Squared Chi-Squared distance
- [`squared_chord()`](https://cmmr.github.io/ecodive/reference/squared_chord.md)
  : Squared Chord distance
- [`squared_euclidean()`](https://cmmr.github.io/ecodive/reference/squared_euclidean.md)
  : Squared Euclidean distance
- [`topsoe()`](https://cmmr.github.io/ecodive/reference/topsoe.md) :
  Topsoe distance
- [`wave_hedges()`](https://cmmr.github.io/ecodive/reference/wave_hedges.md)
  : Wave Hedges distance

### Presence/Absence metrics

Distance measures using binary (presence/absence) data.

- [`hamming()`](https://cmmr.github.io/ecodive/reference/hamming.md) :
  Hamming distance
- [`jaccard()`](https://cmmr.github.io/ecodive/reference/jaccard.md) :
  Jaccard distance
- [`ochiai()`](https://cmmr.github.io/ecodive/reference/ochiai.md) :
  Otsuka-Ochiai dissimilarity
- [`sorensen()`](https://cmmr.github.io/ecodive/reference/sorensen.md) :
  Dice-Sorensen dissimilarity

### Phylogenetic metrics

Distance measures based on phylogenetic branch lengths (UniFrac family).

- [`unweighted_unifrac()`](https://cmmr.github.io/ecodive/reference/unweighted_unifrac.md)
  : Unweighted UniFrac
- [`weighted_unifrac()`](https://cmmr.github.io/ecodive/reference/weighted_unifrac.md)
  : Weighted UniFrac
- [`normalized_unifrac()`](https://cmmr.github.io/ecodive/reference/normalized_unifrac.md)
  : Normalized Weighted UniFrac
- [`generalized_unifrac()`](https://cmmr.github.io/ecodive/reference/generalized_unifrac.md)
  : Generalized UniFrac (GUniFrac)
- [`variance_adjusted_unifrac()`](https://cmmr.github.io/ecodive/reference/variance_adjusted_unifrac.md)
  : Variance-Adjusted Weighted UniFrac

## Utilities & Helpers

Tools for programmatic access, data manipulation, and configuration.

- [`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)
  [`match_metric()`](https://cmmr.github.io/ecodive/reference/metrics.md)
  : Find and Browse Available Metrics
- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  : Read a newick formatted phylogenetic tree.
- [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) :
  Rarefy OTU counts.
- [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md) :
  Number of CPU Cores

## Datasets

Example data provided with the package.

- [`ex_counts`](https://cmmr.github.io/ecodive/reference/ex_counts.md) :
  Example counts matrix
- [`ex_tree`](https://cmmr.github.io/ecodive/reference/ex_tree.md) :
  Example phylogenetic tree
