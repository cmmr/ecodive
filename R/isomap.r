# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# isomap() follows scikit-learn's Isomap (sklearn/manifold/_isomap.py,
# version 1.9.1), and is distributed under its license as well:
#
#   BSD 3-Clause License
#
#   Copyright (c) 2007-2026 The scikit-learn developers.
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   * Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
#   * Neither the name of the copyright holder nor the names of its
#     contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#   THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
#   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#' Isometric Feature Mapping (Isomap)
#'
#' Ordinates samples from their distances, for plotting beta diversity, by
#' replacing each distance with the length of the shortest path between the
#' two samples through a graph that links each sample to its nearest
#' neighbors, then running principal coordinates analysis on those path
#' lengths. Recovers long ecological gradients that [pcoa()] bends into an
#' arch or horseshoe. A parallel C port of scikit-learn's `Isomap` with
#' `metric = "precomputed"`.
#'
#' @inheritParams documentation
#'
#' @param dm  A `dist` object, such as from [beta_div()] or [bray()], or a
#'        symmetric numeric matrix of distances. Must not contain `NA`s or
#'        negative values. Distances, and the geodesic distances summed from
#'        them, must not exceed \eqn{10^{150}}, since they are squared.
#'
#' @param n_axes   The number of axes to compute. Default: `3`. Reduced to
#'        \eqn{n - 1} for \eqn{n} samples.
#'
#' @param n_neighbors   How many nearest neighbors to link each sample to.
#'        Values larger than \eqn{n - 1} are reduced to \eqn{n - 1}.
#'        Default: `5`
#'
#' @param radius   Instead of `n_neighbors`, link each sample to every other
#'        sample within this distance (inclusive). A number of at least 0,
#'        or `NULL` (the default) to use `n_neighbors`. Give one or the
#'        other, not both.
#'
#' @return A numeric matrix with one row per sample and one column per axis,
#'   named `"Isomap1"`, `"Isomap2"`, etc., as described in "Ordination
#'   result" below. Its attributes are:
#'
#'   * `method`: `"isomap"`.
#'   * `eig`: The eigenvalue of each axis, on the same scale as [pcoa()]'s:
#'     of Gower's centered matrix of the squared geodesic distances.
#'   * `prop`: The proportion of the total variance of the geodesic
#'     distances that each axis explains, `eig / trace`.
#'   * `trace`: The total variance: the sum of all eigenvalues, including
#'     negative ones.
#'   * `n_neighbors`: The number of neighbors used, or `NA` with `radius`.
#'   * `radius`: The radius used, or `NA` with `n_neighbors`.
#'   * `components`: The number of connected components of the neighbor
#'     graph, before they were joined (see below). `1` if it was connected.
#'   * `membership`: Each sample's component, an integer vector named by
#'     sample. Components are numbered in order of their first sample.
#'
#' @details
#'
#' **Method**
#'
#' 1. **Neighbor graph.** Each sample is linked to its `n_neighbors` nearest
#'    other samples (or to every sample within `radius`), by an edge as long
#'    as the distance between them. The graph is undirected: a link from
#'    either end counts. Ties are broken by sample order, so that tied
#'    distances and duplicate samples (distance 0) give a well-defined graph.
#' 2. **Joining components.** If the graph falls apart into several
#'    connected components, as it often does for clustered data, every pair
#'    of components is joined by one edge between their two closest samples,
#'    with a warning. The joined graph's geodesics between clusters then run
#'    through those closest pairs. Raise `n_neighbors` to avoid this. If
#'    there are so many components that joining them would add more than 10
#'    times as many edges as the neighbor graph has (and more than 65,536),
#'    `isomap()` stops with an error instead (see "Speed and memory").
#' 3. **Geodesic distances.** The length of the shortest path through the
#'    graph between every pair of samples, by Dijkstra's algorithm, one
#'    search per sample (scikit-learn lets SciPy choose between Dijkstra
#'    and Floyd-Warshall; both give the same lengths, up to rounding).
#' 4. **Embedding.** Classical scaling of the geodesic distances, exactly as
#'    in [pcoa()]: the top `n_axes` eigenvectors of
#'    \eqn{B = -\frac{1}{2}JG^{(2)}J}, where \eqn{G^{(2)}} holds the squared
#'    geodesic distances and \eqn{J} is the centering matrix, scaled by the
#'    square roots of their eigenvalues. Axes whose eigenvalues aren't
#'    positive are dropped with a warning, as in `pcoa()`.
#'
#' With `n_neighbors` = \eqn{n - 1}, every pair of samples is linked, and
#' when `dm` obeys the triangle inequality (e.g. Euclidean distances), the
#' geodesic distances are `dm` itself, so `isomap()` gives the same result
#' as `pcoa()`. Small `n_neighbors` follow the data's shape more closely but
#' are more easily split into components and thrown off by noise.
#'
#'
#' **Speed and memory**
#'
#' Neighbors are found straight from `dm`, in parallel. The geodesic
#' distances are held as one more packed lower triangle the size of `dm`
#' (8 bytes per pair of samples: 10 GB for 50,000 samples), which is the
#' only large allocation; an \eqn{n \times n} matrix is never made. They
#' are kept in double precision, as `dm` is, since the eigenvectors need
#' it. Each search computes only the distances from its sample to later
#' samples, so each pair is computed once; searches run in parallel, `cpus`
#' at a time. The eigenvectors come from block Lanczos iteration, as in
#' `pcoa()`, which needs a few dozen vectors of length \eqn{n}. Besides
#' that, the graph takes about `44 * n * n_neighbors` bytes (or 32 bytes
#' per linked pair of samples with `radius`), and each thread 28 bytes per
#' sample.
#'
#' With \eqn{E} edges (about `n * n_neighbors` to twice that), the time is
#' about \eqn{n^2} to find neighbors, plus \eqn{n (n + E) \log n} for the
#' searches, plus a few dozen passes over the geodesic distances for the
#' eigenvectors, all divided among `cpus`.
#'
#' Joining \eqn{c > 1} components takes one more pass over `dm`, adds
#' \eqn{c (c - 1) / 2} edges (32 bytes each), and needs a table of 16 bytes
#' per pair of components for each thread that scans for the closest pairs
#' (one thread at least; more only while all their tables fit in 16 MB).
#' Every search may relax every joining edge, so a small `radius` or
#' `n_neighbors` that leaves many samples on their own would make memory
#' grow with \eqn{c^2} and time with \eqn{n c^2}. To keep both within a
#' small multiple of what a connected graph needs, `isomap()` refuses (with
#' an error saying how many components it found) when the joining edges
#' would number more than 10 times the neighbor graph's own edges, and more
#' than 65,536. That floor (about 3 MB) lets small datasets join freely: up
#' to 362 samples can all be on their own.
#'
#'
#' **Differences from scikit-learn**
#'
#' `isomap()` gives the same geodesic distances (`dist_matrix_`), up to
#' rounding, and the same embedding (`embedding_`) and eigenvalues as
#' `sklearn.manifold.Isomap(metric = "precomputed")`, except:
#'
#' * **Ties.** scikit-learn chooses among neighbors tied in distance by
#'   `numpy.argpartition()`, in an order that isn't specified, and so can
#'   link different ones among duplicate samples too. `isomap()` always
#'   prefers the earlier sample. When no tie straddles the
#'   `n_neighbors`-th neighbor, the two graphs are the same.
#' * **Axes with eigenvalues that aren't positive.** scikit-learn stops
#'   with an error when one of the top `n_axes` eigenvalues is substantially
#'   negative, and otherwise returns axes of zeros for them. `isomap()`
#'   drops them with a warning, as `pcoa()` does.
#' * **Rounding.** Each geodesic distance is summed along its path from the
#'   earlier of its two samples. SciPy may sum the same path from the other
#'   end, or with Floyd-Warshall's grouping, which can differ in the last
#'   bit or two. The eigenvectors come from block Lanczos rather than LAPACK
#'   or ARPACK. In the package's tests, geodesic distances agree to within
#'   \eqn{4 \times 10^{-16}} (relative), and coordinates to within
#'   \eqn{10^{-12}} of the largest one.
#' * **`n_neighbors` larger than \eqn{n - 1}.** scikit-learn stops with an
#'   error; `isomap()` uses \eqn{n - 1}, with a warning.
#' * **Very many components.** scikit-learn joins any number of them;
#'   `isomap()` refuses past the limit above.
#' * **Sign.** scikit-learn flips each axis so that its largest coordinate
#'   in absolute value is positive. `isomap()` does as [pcoa()] does: it
#'   makes positive the first coordinate whose absolute value is within a
#'   relative \eqn{10^{-8}} of the largest, so that rounding can't break a
#'   tie differently on different computers. The two only differ when an
#'   axis's largest absolute values are tied, or nearly so, with opposite
#'   signs.
#' * **Input and output.** Only `dist` objects or symmetric matrices of
#'   distances; not data matrices or sparse graphs. The diagonal of a matrix
#'   is ignored. `transform()` of new samples and `reconstruction_error()`
#'   are not offered. The components of the neighbor graph are returned as
#'   attributes.
#'
#'
#' **Coming from `vegan::isomap()`**
#'
#' * `vegan::isomap()`'s `k` is the number of neighbors; the number of axes
#'   (`vegan`'s `ndim`) is `isomap()`'s `n_axes`. Give the number of neighbors as
#'   `n_neighbors`.
#' * `vegan`'s `epsilon` corresponds to `radius`.
#' * When the neighbor graph falls apart into several components, `vegan`
#'   stops with an error (or, with `fragmentedOK = TRUE`, ordinates only
#'   the largest component). `isomap()` joins the components at their
#'   closest samples, with a warning, and returns every sample.
#'
#'
#' **Reproducibility**
#'
#' Given the same `dm` and arguments, results, including all attributes,
#' are bit-for-bit identical from run to run, in fresh R sessions, and for
#' any number of `cpus`. Each sample's neighbors, each joining edge, and
#' each search are computed by one thread in a fixed order of operations,
#' with comparisons that break ties by sample order, so tied distances,
#' equally short paths and duplicate samples can't make results depend on
#' which thread finishes first. The eigenvectors are as reproducible as
#' `pcoa()`'s. No random numbers are drawn from R, so `set.seed()` has no
#' effect and R's random number stream is left untouched. Results can
#' differ in the last bits between computers or compilers.
#'
#' `n_axes` larger than \eqn{n - 1} (the most axes possible) is reduced to
#' \eqn{n - 1}.
#'
#' @references
#' Tenenbaum, J. B., de Silva, V., & Langford, J. C. (2000). A global
#' geometric framework for nonlinear dimensionality reduction. *Science*,
#' 290(5500), 2319-2323. \doi{10.1126/science.290.5500.2319}
#'
#' Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
#' *Numerische Mathematik*, 1, 269-271. \doi{10.1007/BF01386390}
#'
#' Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python.
#' *Journal of Machine Learning Research*, 12, 2825-2830.
#' <https://jmlr.org/papers/v12/pedregosa11a.html>
#'
#' Gower, J. C. (1966). Some distance properties of latent root and vector
#' methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
#' \doi{10.1093/biomet/53.3-4.325}
#'
#' @inheritSection ordination_section Ordination result
#'
#' @seealso [pcoa()], [nmds()], [dist_ord()]
#'
#' @export
#' @examples
#'     dm  <- bray(iris[, 1:4])
#'     # Enough neighbours to link setosa to the other two species.
#'     ord <- isomap(dm, n_neighbors = 50)
#'     plot(ord, col = iris$Species)
#'     head(ord)
#'     str(ord)
#'
isomap <- function (dm, n_axes = 3L, n_neighbors = 5L, radius = NULL, cpus = n_cpus()) {

  validate_cpus()
  k <- ord_k(n_axes)

  by_radius <- !is.null(radius)

  if (by_radius) {

    # n_neighbors = NULL with a radius, as in scikit-learn, is fine.
    if (!missing(n_neighbors) && !is.null(n_neighbors))
      stop('Give `n_neighbors` or `radius`, not both.', call. = FALSE)

    # Not ord_arg(): an infinite radius is allowed.
    if (!is.numeric(radius) || length(radius) != 1 || is.na(radius) || radius < 0)
      stop('`radius` must be a number of at least 0, or NULL.', call. = FALSE)

    radius      <- as.double(radius)
    n_neighbors <- NA_integer_

  } else {

    explicit    <- !missing(n_neighbors)
    requested   <- n_neighbors
    n_neighbors <- ord_arg(
      n_neighbors, '`n_neighbors` must be a whole number of at least 1.',
      whole = TRUE, lo = 1 )
    radius      <- NA_real_
  }


  dv     <- dist_vector(dm)
  dm     <- dv$dm
  n      <- dv$n
  labels <- dv$labels

  k <- min(k, n - 1L)

  if (!by_radius && n_neighbors > n - 1L && n >= 2) {
    if (explicit)
      warning(sprintf(
        '`n_neighbors` = %s is more than the %d other samples; using %d.',
        format(requested), n - 1L, n - 1L ), call. = FALSE)
    n_neighbors <- n - 1L
  }

  if (n < 2) {

    # One sample or none: no axes.
    points     <- matrix(numeric(0), nrow = n, ncol = 0)
    eig        <- numeric(0)
    trace      <- 0
    components <- n
    membership <- seq_len(n)

    if (!by_radius) n_neighbors <- 0L

  } else {

    # Joining c components: one table of closest pairs (16 bytes per pair
    # of components) per scanning thread, with extra threads only while
    # all tables fit in 16 MB; refused if its c(c - 1) / 2 edges would be
    # more than 10 times the neighbor graph's own and more than 65,536.
    geo <- .Call(
      C_isomap, dm, n, if (by_radius) 0L else n_neighbors, radius, cpus,
      c(16 * 2^20, 10, 65536), FALSE )

    # dist_vector() caps each distance at 1e150, so that pcoa() can square
    # them, but a geodesic sums up to n - 1 of them. Hold the geodesics,
    # which go to pcoa's engine next, to the same limit.
    if (geo$max > 1e150)
      stop(
        '`dm` is too large: its geodesic distances (sums of distances) ',
        'exceed 1e150, too large to square; rescale `dm`.', call. = FALSE )

    components <- geo$components
    membership <- geo$membership

    if (components > 1)
      warning(sprintf(paste(
        'The neighbor graph has %d connected components; they were joined',
        'at their closest samples. Increase `%s` to avoid this.' ),
        components, if (by_radius) 'radius' else 'n_neighbors' ), call. = FALSE)

    res <- .Call(C_pcoa, geo$geo, n, k, cpus)
    rm(geo)

    if (!res$converged) # nocov start
      warning('Isomap eigenvectors did not converge; coordinates may be inexact.', call. = FALSE)
    # nocov end

    keep   <- ord_keep(res$eig, k)
    points <- res$points[, keep, drop = FALSE]
    eig    <- res$eig[keep]
    trace  <- res$trace
  }

  names(membership) <- labels

  ord_result(
    points, 'isomap', labels, prefix = 'Isomap',
    eig         = eig,
    prop        = if (trace > 0) eig / trace else eig,
    trace       = trace,
    n_neighbors = n_neighbors,
    radius      = radius,
    components  = components,
    membership  = membership )
}
