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
- [`robust_aitchison()`](https://cmmr.github.io/ecodive/reference/robust_aitchison.md)
  : Robust Aitchison distance
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

## Compositional Data

Log ratios for relative abundances. The centered log ratios (clr and
rclr) are norm choices; see normalize().

### Zeros

Replace zeros before taking log ratios, instead of adding a pseudocount.
replace_zeros() runs any method by name.

- [`replace_zeros()`](https://cmmr.github.io/ecodive/reference/replace_zeros.md)
  : Zero Replacement
- [`replace_zeros_mult()`](https://cmmr.github.io/ecodive/reference/replace_zeros_mult.md)
  : Multiplicative Zero Replacement
- [`list_zero_replacements()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)
  [`match_zero_replacement()`](https://cmmr.github.io/ecodive/reference/zero_replacements.md)
  : Find and Browse Zero Replacement Methods

### Log-ratio transforms

Coordinates for analyses of the transformed features themselves. For
distances, use aitchison().

- [`alr()`](https://cmmr.github.io/ecodive/reference/alr.md) : Additive
  Log-Ratio Transform
- [`ilr()`](https://cmmr.github.io/ecodive/reference/ilr.md) : Isometric
  Log-Ratio Transform

### Log ratios

One log ratio per sample, between two groups of features, to use as a
variable in tests and models.

- [`log_ratio()`](https://cmmr.github.io/ecodive/reference/log_ratio.md)
  : Log Ratio of Summed Abundances
- [`balance()`](https://cmmr.github.io/ecodive/reference/balance.md) :
  Isometric Log-Ratio Balance

### Proportionality

Which features keep a constant ratio across samples, the compositionally
valid alternative to correlation. A dist over the features.

- [`vlr()`](https://cmmr.github.io/ecodive/reference/vlr.md) : Variance
  of Log Ratios
- [`rho_p()`](https://cmmr.github.io/ecodive/reference/rho_p.md) :
  Proportionality Correlation
- [`phi_s()`](https://cmmr.github.io/ecodive/reference/phi_s.md) :
  Symmetric Proportionality

## Ordination

Place samples on a few axes for plotting. dist_ord() runs the methods
that take beta diversity distances, and count_ord() the methods that
take the count table.

- [`dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ord.md) :
  Distance-based Ordination
- [`count_ord()`](https://cmmr.github.io/ecodive/reference/count_ord.md)
  : Count Table Ordination
- [`pcoa()`](https://cmmr.github.io/ecodive/reference/pcoa.md) :
  Principal Coordinates Analysis (PCoA)
- [`pcoa_biplot()`](https://cmmr.github.io/ecodive/reference/pcoa_biplot.md)
  : PCoA Biplot Scores
- [`mds()`](https://cmmr.github.io/ecodive/reference/mds.md) : Metric
  Multidimensional Scaling (MDS)
- [`nmds()`](https://cmmr.github.io/ecodive/reference/nmds.md) :
  Non-metric Multidimensional Scaling (NMDS)
- [`isomap()`](https://cmmr.github.io/ecodive/reference/isomap.md) :
  Isometric Feature Mapping (Isomap)
- [`tsne()`](https://cmmr.github.io/ecodive/reference/tsne.md) :
  t-distributed Stochastic Neighbor Embedding (t-SNE)
- [`umap()`](https://cmmr.github.io/ecodive/reference/umap.md) : Uniform
  Manifold Approximation and Projection (UMAP)
- [`dbrda()`](https://cmmr.github.io/ecodive/reference/dbrda.md) :
  Distance-based Redundancy Analysis (db-RDA)
- [`capscale()`](https://cmmr.github.io/ecodive/reference/capscale.md) :
  Constrained Analysis of Principal Coordinates (capscale)
- [`rda()`](https://cmmr.github.io/ecodive/reference/rda.md) :
  Redundancy Analysis (RDA)
- [`cca()`](https://cmmr.github.io/ecodive/reference/cca.md) : Canonical
  Correspondence Analysis (CCA)
- [`ca()`](https://cmmr.github.io/ecodive/reference/ca.md) :
  Correspondence Analysis (CA)
- [`pca()`](https://cmmr.github.io/ecodive/reference/pca.md) : Principal
  Component Analysis (PCA)
- [`rpca()`](https://cmmr.github.io/ecodive/reference/rpca.md) : Robust
  Aitchison PCA (RPCA)
- [`envfit()`](https://cmmr.github.io/ecodive/reference/envfit.md) : Fit
  Variables or Taxa onto an Ordination
- [`bioenv()`](https://cmmr.github.io/ecodive/reference/bioenv.md) :
  Best Subset of Environmental Variables (BIOENV)
- [`varpart()`](https://cmmr.github.io/ecodive/reference/varpart.md) :
  Variation Partitioning
- [`forward_select()`](https://cmmr.github.io/ecodive/reference/forward_select.md)
  : Forward Selection of Explanatory Variables
- [`procrustes()`](https://cmmr.github.io/ecodive/reference/procrustes.md)
  : Procrustes Rotation of One Ordination onto Another
- [`procrustes_test()`](https://cmmr.github.io/ecodive/reference/procrustes_test.md)
  : Procrustes Randomization Test (PROTEST)
- [`list_dist_ords()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  [`match_dist_ord()`](https://cmmr.github.io/ecodive/reference/dist_ords.md)
  : Find and Browse Distance-based Ordination Methods
- [`list_count_ords()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  [`match_count_ord()`](https://cmmr.github.io/ecodive/reference/count_ords.md)
  : Find and Browse Count Table Ordination Methods

### Supervised

Ordinations that look for the axes that separate known groups. They
always find some, so validate them with a permutation test such as
permanova().

- [`lda()`](https://cmmr.github.io/ecodive/reference/lda.md) : Linear
  Discriminant Analysis (LDA)
- [`plsda()`](https://cmmr.github.io/ecodive/reference/plsda.md) :
  Partial Least Squares Discriminant Analysis (PLS-DA)
- [`oplsda()`](https://cmmr.github.io/ecodive/reference/oplsda.md) :
  Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)

## Hypothesis Tests

Test whether sample metadata explain beta diversity distances.

- [`dist_test()`](https://cmmr.github.io/ecodive/reference/dist_test.md)
  : Distance Test Wrapper Function
- [`permanova()`](https://cmmr.github.io/ecodive/reference/permanova.md)
  : Permutational Multivariate Analysis of Variance (PERMANOVA)
- [`mantel()`](https://cmmr.github.io/ecodive/reference/mantel.md) :
  Mantel Test
- [`mantel_partial()`](https://cmmr.github.io/ecodive/reference/mantel_partial.md)
  : Partial Mantel Test
- [`mrpp()`](https://cmmr.github.io/ecodive/reference/mrpp.md) :
  Multi-Response Permutation Procedure (MRPP)
- [`anosim()`](https://cmmr.github.io/ecodive/reference/anosim.md) :
  Analysis of Similarities (ANOSIM)
- [`permdisp()`](https://cmmr.github.io/ecodive/reference/permdisp.md) :
  Homogeneity of Multivariate Dispersions (PERMDISP)
- [`list_dist_tests()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  [`match_dist_test()`](https://cmmr.github.io/ecodive/reference/dist_tests.md)
  : Find and Browse Distance Tests

## Interpreting Differences

Which taxa drive the differences between groups.

- [`simper()`](https://cmmr.github.io/ecodive/reference/simper.md) :
  Similarity Percentages (SIMPER)

## Utilities & Helpers

Tools for programmatic access, data manipulation, and configuration.

- [`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.md)
  [`match_metric()`](https://cmmr.github.io/ecodive/reference/metrics.md)
  : Find and Browse Available Metrics
- [`read_tree()`](https://cmmr.github.io/ecodive/reference/read_tree.md)
  : Read a newick formatted phylogenetic tree.
- [`rarefy()`](https://cmmr.github.io/ecodive/reference/rarefy.md) :
  Rarefy Observation Counts
- [`normalize()`](https://cmmr.github.io/ecodive/reference/normalize.md)
  : Normalize a Count Table
- [`n_cpus()`](https://cmmr.github.io/ecodive/reference/n_cpus.md) :
  Number of CPU Cores
- [`blas_threads()`](https://cmmr.github.io/ecodive/reference/blas_threads.md)
  : Threads of the BLAS Library

## Datasets

Example data provided with the package.

- [`ex_counts`](https://cmmr.github.io/ecodive/reference/ex_counts.md) :
  Example counts matrix
- [`ex_tree`](https://cmmr.github.io/ecodive/reference/ex_tree.md) :
  Example phylogenetic tree
