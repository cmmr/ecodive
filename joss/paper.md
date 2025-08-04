---
title: 'ecodive: Fast Implementations of Ecological Diversity Metrics in R'
authors:
  - name: Daniel P Smith
    orcid: 0000-0002-2479-2044
    corresponding: true
    affiliation: 1
  - name: Joseph F Petrosino
    orcid: 0000-0002-4046-6898
    affiliation: 1
affiliations:
 - name: Alkek Center for Metagenomics and Microbiome Research, Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, Texas, USA
   index: 1
   ror: 02pttbw34
tags:
  - R
  - microbiome
  - unifrac
  - bioinformatics
output: pdf_document
date: 04 August 2025
bibliography: paper.bib
---


![Ecodive package logo](../logo/ecodive.png){height="150pt"}


# Summary

In the context of ecology, diversity measures the distribution of different
species within in a community. This calculation may include the number of
species present, relative abundances, evolutionary relationships, or a
combination thereof. Alpha diversity metrics consider a single community in
isolation, whereas beta diversity metrics compute the dissimilarity between two
communities.

Applying diversity metrics to large collections of communities, for instance
thousands of gut microbiome samples, can offer insights into how specific
disease states may be predicted or diagnosed based on ecological "fingerprints".




# Statement of Need

Some diversity metrics, such as Faith's PD [@FaithPD] and UniFrac [@UniFrac],
require complex integration species counts with evolutionary distances.
Furthermore, processing thousands of communities is computationally intensive
and best implemented with parallel processing and compiled libraries. For these
reasons, the ecodive R package was developed to handle these challenges so that
R users don't have to.




## Related Works

There are currently five other R packages that can calculate alpha and beta
diversity metrics: abdiv [@abdiv], ampvis2 [@ampvis2], GUniFrac [@GUniFrac],
phyloseq [@phyloseq], and vegan [@vegan]. However, ecodive provides an
implementation which is both faster and more memory efficient.

The bench R package [@bench] was used to compare abdiv, ampvis2, ecodive,
GUniFrac, phyloseq, and vegan. The benchmarking runs are detailed in the
benchmark vignette, which is available from within R with
`vignette('benchmark')` and online at
<https://cmmr.github.io/ecodive/articles/benchmark.html>.


![UniFrac benchmarks. Ecodive is up to 15x to 2800x faster and uses 60x - 25000x less memory.](../man/figures/unifrac-benchmark.png)


![Classic beta diversity benchmarks. Ecodive is up to 23x to 160x faster and uses 0.8x to 640x less memory.](../man/figures/bdiv-benchmark.png)


![Alpha  diversity benchmarks. Ecodive is up to 10x to 40x faster and uses 5x to 25x less memory.](../man/figures/adiv-benchmark.png)




# Acknowledgements

# References
