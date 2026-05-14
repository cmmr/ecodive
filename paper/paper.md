---
title: 'ecodive: Parallel and Memory-Efficient R Package for Ecological Diversity Analysis'
authors:
- name: Daniel P Smith
  orcid: "0000-0002-2479-2044"
  corresponding: true
  affiliation: "1, 2"
- name: Sara J Javornik Cregeen
  orcid: "0009-0000-2698-6478"
  affiliation: "1, 2"
- name: Joseph F Petrosino
  orcid: "0000-0002-4046-6898"
  affiliation: "1, 2"
date: "29 January 2026"
output:
  word_document: default
  pdf_document: default
tags:
- R
- ecology
- bioinformatics
- microbiome
- diversity
- high-performance computing
affiliations:
- name: The Alkek Center for Metagenomics and Microbiome Research, Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX 77030, USA
  index: 1
- name: Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, TX, USA
  index: 2
  ror: 02pttbw34
bibliography: paper.bib
---

# Summary

Modern ecology relies on quantifying the biological complexity of communities across diverse scales, from microscopic human microbiomes to vast forest ecosystems. Researchers use diversity metrics to describe the variety within a site (alpha diversity) and the compositional differences between sites (beta diversity). `ecodive` is an R package built to compute these metrics with high efficiency. It provides a comprehensive, unified toolset that allows researchers to analyze large-scale datasets that were previously computationally prohibitive. By integrating optimized C-based algorithms with a parallel processing engine, `ecodive` enables rapid, high-throughput ecological analysis through a lightweight, portable framework.

# Statement of Need

A primary challenge in modern ecological analysis is the management of high-dimensional data. As sequencing technologies improve, datasets are growing to include thousands of samples and tens of thousands of unique taxa. Beta diversity calculations, which involve comparing every sample to every other sample, exhibit $O(n^2)$ complexity. This quadratic scaling creates a significant bottleneck; doubling a dataset's size quadruples the processing time and memory requirements, frequently exceeding the capacity of standard research workstations.

Furthermore, the software landscape for ecological metrics is fragmented. A researcher needing to calculate a specific set of indices - for example, Faith's Phylogenetic Diversity, Bray-Curtis dissimilarity, and UniFrac distances - often must install and manage multiple R packages (`picante`, `vegan`, `GUniFrac`), each with different dependencies, input formats, and performance limitations.

`ecodive` solves these problems by providing a centralized, high-performance library. It targets ecologists, microbiologists, and bioinformaticians who require a robust, dependency-free solution for diversity analysis. By consolidating 50 standard metrics into a single, optimized framework, it eliminates the need for "package hopping" and enables the analysis of massive datasets on standard hardware.

# State of the Field

While several R packages support diversity analysis, `ecodive` offers a unique scholarly contribution by bridging the gap between specialized high-performance tools and comprehensive ecological libraries. As shown in Table 1, the current landscape is defined by three primary trade-offs: analytical breadth, computational efficiency, and infrastructure overhead.

| R Package      | Alpha | Beta |      Uni     |      Par     |      Cmp     | Dep  | Citation               |
|:----------     |   ---:|  ---:|          ---:|          ---:|          ---:|  ---:|:-----------------------|
| `ecodive`      |    14 |   36 | $\checkmark$ | $\checkmark$ | $\checkmark$ |    0 | This work              |
| `abdiv`        |    17 |   48 | $\checkmark$ | $-$          | $-$          |    5 | @abdiv                 |
| `adiv`         |    42 |   58 | $-$          | $-$          | $-$          |   97 | @adiv                  |
| `ampvis2`      |     6 |   25 | $\checkmark$ | $-$          | $-$          |   75 | @ampvis2               |
| `ecodist`      |     0 |    9 | $-$          | $-$          | $\checkmark$ |   10 | @ecodist               |
| `entropart`    |    11 |    9 | $-$          | $-$          | $-$          |   82 | @entropart             |
| `GUniFrac`     |     0 |    4 | $\checkmark$ | $-$          | $\checkmark$ |   47 | @GUniFrac              |
| `labdsv`       |     2 |    7 | $-$          | $-$          | $\checkmark$ |    8 | @labdsv                |
| `OmicFlow`     |     6 |    9 | $\checkmark$ | $\checkmark$ | $\checkmark$ |   90 | @OmicFlow              |
| `parallelDist` |     0 |   26 | $-$          | $\checkmark$ | $\checkmark$ |    2 | @parallelDist          |
| `philentropy`  |     1 |   43 | $-$          | $\checkmark$ | $\checkmark$ |    4 | @philentropy           |
| `phyloregion`  |     9 |   10 | $\checkmark$ | $-$          | $-$          |   64 | @phyloregion           |
| `phyloseq`     |     7 |   43 | $\checkmark$ | $\checkmark$ | $-$          |   54 | @phyloseq              |
| `picante`      |    14 |    9 | $\checkmark$ | $-$          | $-$          |   11 | @picante               |
| `tabula`       |    16 |   12 | $-$          | $-$          | $-$          |    2 | @tabula                |
| `vegan`        |    17 |   52 | $-$          | $-$          | $\checkmark$ |    7 | @vegan                 |
Table: Comparison of community diversity metrics and architectural features across 16 R packages. Metrics were identified using a large language model protocol on package reference manuals to categorize metrics and filter for symmetry and ecological relevance. Alpha: Number of alpha diversity metrics; Beta: Number of beta diversity metrics; Uni: $\checkmark$ indicates UniFrac support; Par: Parallelized calculation of metrics; Cmp: Implementation uses compiled code (e.g. C/C++); Dep: Transitive count of hard R dependencies. Full methodology and LLM prompts are available in `ecodive`'s benchmarking documentation.

## 1. Analytical Breadth vs. Computational Scalability

The ecological software landscape often requires researchers to choose between a wide variety of metrics and the speed necessary for high-dimensional community data. As shown in Table 1, packages like `adiv` and `abdiv` offer significant analytical breadth, providing 100 and 65 total metrics, respectively. However, both lack parallelized execution (Par) and compiled backend implementations (Cmp), which limits their utility for large-scale datasets. `ecodive` overcomes this limitation by providing a comprehensive suite of 50 metrics - comparable in scope to the foundational `vegan` package (69 metrics) - while utilizing a parallelized C engine to ensure scalability.

## 2. Domain-Specific Functionality and UniFrac Support

General-purpose high-performance libraries often fail to address specific ecological needs, such as phylogenetic awareness. For instance, despite their high metric counts, neither `vegan` (69 metrics), `adiv` (100 metrics), nor `philentropy` (44 metrics) provide UniFrac support. While specialized tools like `picante` and `GUniFrac` do include these indices, they operate serially and carry significant dependency burdens. `ecodive` is unique in combining extensive general metrics with parallelized UniFrac support, functioning as a unified tool for both phylogenetic and non-phylogenetic analyses without sacrificing performance.

## 3. Infrastructure Stability and Dependency Management

A critical but often overlooked limitation in existing software is the "transitive burden" of recursive dependencies. Table 1 reveals significant infrastructure overhead across the field: `adiv` requires 97 hard dependencies, `OmicFlow` requires 90, and `phyloregion` requires 64. Such high counts increase maintenance complexity and supply chain vulnerability for research pipelines. `ecodive` addresses this by delivering its full feature set - including parallelization and compiled performance - with zero external R dependencies. This zero-dependency architecture ensures high community readiness and long-term stability for ecological research.

## The `ecodive` Contribution

`ecodive` offers a unique scholarly contribution by bridging these gaps into a single framework. It implements 50 symmetric metrics chosen specifically for their relevance to ecological ordination. Built on a parallelized C engine, the package matches or exceeds the performance of specialized tools while ensuring that results remain numerically identical to established packages. This performance is delivered within a zero-dependency architecture, which minimizes maintenance overhead and protects against the supply chain vulnerabilities common in complex bioinformatic pipelines.

# Software Design

The architecture of `ecodive` balances the user-friendly conventions of R with the raw performance of C. A critical design trade-off centered on data representation.

Most R users work with dense matrices where samples are rows and features are columns. Standard R functions like `dist()` expect this format. However, ecological matrices are typically 90-99% zeros (sparse). Storing them as dense matrices wastes gigabytes of RAM, and processing them row-by-row is cache-inefficient for many distance algorithms.

To address this, `ecodive` maintains the standard R interface (samples-as-rows) but fundamentally alters the backend data structure:

1.  **Transparent Conversion:** When a standard matrix is passed to `ecodive`, it is internally converted into a column-compressed sparse matrix (`dgCMatrix`) with samples transposed to columns. This incurs a one-time overhead but allows the C engine to skip zeros entirely and access memory in a cache-friendly, column-major pattern.
2.  **Transformation Bypass:** For extremely large datasets where the overhead of this transformation is non-trivial, users can manually provide data in the native `dgCMatrix` format (samples as columns). `ecodive` detects this optimized state and bypasses the transformation step, operating directly on the existing C pointers. This allows for "zero-copy" analysis of massive datasets.
3.  **Parallelization Strategy:** `ecodive` employs a direct implementation using the standard POSIX threads (`pthreads`) library, avoiding the memory duplication overhead of forking processes found in R's `parallel` package. This design enables fine-grained, dynamic load balancing, ensuring efficient execution even when calculating partial distance matrices.

# Research Impact Statement

`ecodive` has demonstrated immediate utility in high-dimensional microbiome studies. The core C algorithms in `ecodive` were originally developed for and deployed in the `rbiom` package [@rbiom]. As part of `rbiom`, these optimized metrics have already been utilized in diverse microbial ecology studies, including research on preterm infant microbiomes [@AhearnFord2025], dietary interventions [@DiMattia2025], and relationship satisfaction [@Cheng2023]. `ecodive` extracts these proven, high-performance components into a standalone, lightweight library to make them accessible to the broader R ecosystem without `rbiom`'s specific visualization and data structure dependencies.

To reflect realistic modern research conditions, comprehensive [benchmarking](https://cmmr.github.io/ecodive/articles/benchmark.html) was conducted on a workstation with a 6-core Intel i5-9600K (3.70GHz) CPU and 64GB RAM running Windows 11. While `ecodive`, `parallelDist`, `philentropy`, `phyloseq`, `OmicFlow`, and `ampvis2` leveraged all six available CPU cores, the results show that parallelization is not the sole determinant for performance. In several metrics, specialized serial implementations outperformed parallel-capable packages - for example, the serial `phyloregion` was significantly faster than the parallelized `ampvis2` for UniFrac calculations - indicating that efficient data handling and core algorithms are as critical as raw multi-threading.

![**Benchmarking results.** Execution time (x-axis) vs. peak memory usage (y-axis) for various diversity metrics across 16 R packages. `ecodive` (highlighted) consistently occupies the bottom-left quadrant, indicating high speed and low memory footprint. Note the log scale on both axes.](figures/fig1.svg)

* **Speed:** For the widely-used Unweighted UniFrac metric ($N=50$ samples), `ecodive` completed calculations in 6.9ms, compared to 2.6s for `picante` (~400x faster) and 320ms for `phyloseq` (~50x faster).
* **Scalability:** For standard Bray-Curtis dissimilarity ($N=1006$ samples), `ecodive` processed the matrix in 19ms, whereas `vegan` required 2.0s (~100x faster).
* **Memory:** `ecodive`'s sparse architecture reduced memory allocation for large operations from gigabytes (in `abdiv` or `tabula`) to megabytes, enabling analyses on laptops that previously required clusters.

`ecodive` is available for installation via CRAN and Conda-Forge, supported by comprehensive vignettes that guide users through metric selection and performance tuning to ensure immediate utility in diverse research environments.

# Example Usage

`ecodive` is designed for ease of use and integrates seamlessly with existing bioinformatics workflows, such as those using `phyloseq` objects. For example, calculating weighted UniFrac distances is straightforward:

``` r
data(esophagus, package = 'phyloseq')
ecodive::weighted_unifrac(esophagus)
#>           B         C
#> C 0.1050480          
#> D 0.1401124 0.1422409
```

# AI Usage Disclosure

Generative AI tools (Google Gemini) were used to assist in the drafting and revision of this manuscript and the generation of documentation. No AI tools were used to write the functional source code (R or C) of the software. All AI-generated text was critically reviewed, verified for accuracy, and edited by the authors.

# Acknowledgements

This study was supported by NIH/NIAD (Grant number U19 AI144297), and Baylor College of Medicine and Alkek Foundation Seed.

# References
