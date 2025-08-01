% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation.r
\name{documentation}
\alias{documentation}
\title{documentation}
\arguments{
\item{counts}{An OTU abundance matrix where each column is a sample, and
each row is an OTU. Any object coercible with \code{as.matrix()} can be
given here, as well as \code{phyloseq}, \code{rbiom}, \code{SummarizedExperiment},
and \code{TreeSummarizedExperiment} objects.}

\item{weighted}{If \code{TRUE}, the algorithm takes relative abundances into
account. If \code{FALSE}, only presence/absence is considered.}

\item{normalized}{For weighted UniFrac only, normalize distances by the
total branch length. Options: \code{TRUE} or \code{FALSE}.}

\item{alpha}{How much weight to give to relative abundances; a value
between 0 and 1, inclusive. Setting \code{alpha=1} is equivalent to
\code{normalized_unifrac()}.}

\item{tree}{A \code{phylo}-class object representing the phylogenetic tree for
the OTUs in \code{counts}. The OTU identifiers given by \code{colnames(counts)}
must be present in \code{tree}. Can be omitted if a tree is embedded with
the \code{counts} object or as \code{attr(counts, 'tree')}.}

\item{pairs}{Which combinations of samples should distances be
calculated for? The default value (\code{NULL}) calculates all-vs-all.
Provide a numeric or logical vector specifying positions in the
distance matrix to calculate. See examples.}

\item{cpus}{How many parallel processing threads should be used. The
default, \code{n_cpus()}, will use all logical CPU cores.}
}
\description{
documentation
}
\keyword{internal}
