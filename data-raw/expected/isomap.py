#!/usr/bin/env python3
# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
"""
Reference fixtures for ecodive::isomap(), from scikit-learn's Isomap.

Usage (from the package root):

    python data-raw/expected/isomap.py data-raw/expected/isomap

Requires numpy, scipy and scikit-learn (fixtures were made with
scikit-learn 1.9.1, scipy 1.18.1 and numpy 2.5.5). Writes one text file
per dataset, <name>.txt. data-raw/expected/isomap.r converts them to
tests/tinytests/expected/isomap.rds, read by tests/tinytests/test-isomap.r.
There is no Python where the package is developed, so the committed text
files are the reference.

Each dataset's input distances are either given as coordinates plus a
metric ("euclidean" or "manhattan"), which the generator computes with the
same arithmetic as R's stats::dist(), so that R reproduces them bit for
bit, or listed in full ("dist").

The expected results are computed the way sklearn.manifold.Isomap's
_fit_transform() does with metric="precomputed", from scikit-learn's and
SciPy's own functions (connected_components, _fix_connected_components,
shortest_path(directed=False), KernelPCA), except for one step: the k
nearest neighbors are chosen by distance, then sample index. Isomap picks
among tied neighbors with numpy.argpartition(), in an unspecified order;
ecodive breaks ties by index. Isomap is also run end to end on every
dataset, and the line "sklearn_same" records whether its graph, geodesic
distances and embedding agree with the reference (always so when no tie
straddles a k-th neighbor).

File format: one item per line, "key value value ...", separated by single
spaces. Matrices are listed column by column.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import sklearn
from scipy.sparse.csgraph import connected_components, shortest_path
from sklearn.decomposition import KernelPCA
from sklearn.manifold import Isomap
from sklearn.neighbors import radius_neighbors_graph
from sklearn.utils.graph import _fix_connected_components


#--------------------------------------------------------------
# Distances with R's stats::dist() arithmetic: per pair, a
# left-to-right sum over the columns, then sqrt() for Euclidean.
#--------------------------------------------------------------
def r_dist(x, metric):
    n = x.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            s = 0.0
            for c in range(x.shape[1]):
                dev = float(x[i, c]) - float(x[j, c])
                s += dev * dev if metric == "euclidean" else abs(dev)
            d[i, j] = d[j, i] = np.sqrt(s) if metric == "euclidean" else s
    return d


def lower(d):
    """Packed lower triangle, column by column (R's dist order)."""
    n = d.shape[0]
    return np.concatenate([d[j + 1:, j] for j in range(n - 1)])


#--------------------------------------------------------------
# k nearest neighbors, by distance and then index, as a CSR
# graph like kneighbors_graph(mode="distance").
#--------------------------------------------------------------
def knn_graph(d, k):
    n = d.shape[0]
    rows, cols, vals = [], [], []
    for i in range(n):
        others = np.array([j for j in range(n) if j != i])
        order  = np.lexsort((others, d[i, others]))[:k]
        for j in others[order]:
            rows.append(i)
            cols.append(j)
            vals.append(d[i, j])
    return sp.csr_array((vals, (rows, cols)), shape=(n, n))


def reference(d, k_axes, n_neighbors=None, radius=None):
    """Isomap's _fit_transform() with the reference neighbor rule."""
    if n_neighbors is not None:
        nbg = knn_graph(d, n_neighbors)
    else:
        nbg = radius_neighbors_graph(d, radius=radius, metric="precomputed", mode="distance")
    ncomp, labels = connected_components(nbg)
    if ncomp > 1:
        # A CSR graph, as Isomap passes: assigning a 0 into a LIL matrix
        # would store nothing, dropping zero-length joining edges.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # SparseEfficiencyWarning
            nbg = _fix_connected_components(
                X=d, graph=sp.csr_array(nbg), n_connected_components=ncomp,
                component_labels=labels, mode="distance", metric="precomputed")
    geo = shortest_path(nbg, directed=False)
    kpca = KernelPCA(n_components=k_axes, kernel="precomputed")
    emb  = kpca.fit_transform(-0.5 * geo ** 2)
    return dict(graph=nbg, ncomp=ncomp, labels=labels, geo=geo,
                emb=emb, eig=kpca.eigenvalues_)


def edges(graph):
    """Undirected edges (i < j), 1-based, sorted, of a sparse graph."""
    g = sp.coo_array(graph)
    pairs = {(min(i, j), max(i, j)) for i, j in zip(g.row, g.col)}
    return sorted((i + 1, j + 1) for i, j in pairs)


def fmt(x):
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return repr(float(x))


def write(path, name, d, k_axes, n_neighbors=None, radius=None,
          coords=None, metric=None, geo=True, note=""):

    ref = reference(d, k_axes, n_neighbors, radius)

    # The same through sklearn's Isomap, end to end.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        iso = Isomap(n_neighbors=n_neighbors, radius=radius,
                     n_components=k_axes, metric="precomputed")
        emb = iso.fit_transform(d)

    same_geo = np.allclose(iso.dist_matrix_, ref["geo"], rtol=1e-13, atol=0)
    same_emb = np.allclose(emb, ref["emb"], rtol=0, atol=1e-10 * np.abs(emb).max())
    same_knn = True
    if n_neighbors is not None:
        same_knn = edges(knn_graph(d, n_neighbors)) == \
            edges(iso.nbrs_.kneighbors_graph(mode="distance"))
    same = same_knn and same_geo and same_emb

    n   = d.shape[0]
    out = []
    out.append("# isomap reference: scikit-learn %s, scipy %s, numpy %s"
               % (sklearn.__version__, scipy.__version__, np.__version__))
    if note:
        out.append("# " + note)
    out.append("name " + name)
    out.append("n %d" % n)
    out.append("k %d" % k_axes)
    out.append("n_neighbors " + ("NA" if n_neighbors is None else str(n_neighbors)))
    out.append("radius " + ("NA" if radius is None else fmt(radius)))
    if coords is not None:
        out.append("metric " + metric)
        out.append("ncol %d" % coords.shape[1])
        out.append("coords " + " ".join(fmt(v) for v in coords.T.ravel()))
    else:
        out.append("dist " + " ".join(fmt(v) for v in lower(d)))
    out.append("sklearn_same %d" % int(same))
    out.append("components %d" % ref["ncomp"])
    out.append("membership " + " ".join(str(v + 1) for v in ref["labels"]))
    e = edges(ref["graph"])
    out.append("edges " + " ".join("%d %d" % ij for ij in e))
    out.append("eig " + " ".join(fmt(v) for v in ref["eig"]))
    out.append("embedding " + " ".join(fmt(v) for v in ref["emb"].T.ravel()))
    if geo:
        out.append("geo " + " ".join(fmt(v) for v in lower(ref["geo"])))

    Path(path, "%s.txt" % name).write_text("\n".join(out) + "\n")
    print("%-14s n=%3d components=%2d sklearn_same=%d (knn %d, geo %d, emb %d)"
          % (name, n, ref["ncomp"], same, same_knn, same_geo, same_emb))


def swiss(rng, n, noise=0.0):
    t = 1.5 * np.pi * (1 + 2 * rng.random(n))
    h = 10 * rng.random(n)
    x = np.column_stack([t * np.cos(t), h, t * np.sin(t)])
    return np.round(x + noise * rng.standard_normal(x.shape), 3)


def main(outdir):

    Path(outdir).mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261117)

    # Swiss roll: a curled 2-D sheet in 3-D.
    x = swiss(rng, 40)
    write(outdir, "swiss", r_dist(x, "euclidean"), 3, n_neighbors=7,
          coords=x, metric="euclidean")

    # Bigger swiss roll, embedding only.
    x = swiss(rng, 150)
    write(outdir, "swiss150", r_dist(x, "euclidean"), 2, n_neighbors=10,
          coords=x, metric="euclidean", geo=False)

    # Radius neighbors, connected and not.
    x = swiss(rng, 36)
    d = r_dist(x, "euclidean")
    write(outdir, "radius", d, 2, radius=7.0, coords=x, metric="euclidean")
    write(outdir, "radius_split", d, 2, radius=2.0, coords=x, metric="euclidean")

    # Grid with Manhattan distances: integer distances, many ties.
    g = np.array([(i, j) for i in range(6) for j in range(6)], dtype=float)
    perm = rng.permutation(len(g))
    g = g[perm]
    write(outdir, "ties", r_dist(g, "manhattan"), 2, n_neighbors=3,
          coords=g, metric="manhattan", note="ties straddle the 3rd neighbor")

    # Clusters that the neighbor graph leaves apart.
    centers = np.array([[0, 0], [12, 1], [5, 14], [20, 18]])
    x = np.vstack([c + rng.standard_normal((9, 2)) for c in centers])
    x = np.round(x[rng.permutation(len(x))], 3)
    write(outdir, "clusters", r_dist(x, "euclidean"), 2, n_neighbors=3,
          coords=x, metric="euclidean")

    # Clusters on a grid: tied closest pairs between components.
    a = [(i, j) for i in range(3) for j in range(3)]
    b = [(i + 6, j) for i in range(3) for j in range(3)]
    c = [(i, j + 7) for i in range(3) for j in range(3)]
    x = np.array(a + b + c, dtype=float)
    x = x[rng.permutation(len(x))]
    write(outdir, "cluster_ties", r_dist(x, "manhattan"), 2, n_neighbors=2,
          coords=x, metric="manhattan", note="tied joining pairs")

    # Duplicate samples, including a group larger than n_neighbors + 1.
    x = np.round(rng.standard_normal((16, 2)) * 3, 3)
    x = np.vstack([x, x[[0, 0, 0, 0, 3, 3, 7]]])
    x = x[rng.permutation(len(x))]
    write(outdir, "duplicates", r_dist(x, "euclidean"), 2, n_neighbors=3,
          coords=x, metric="euclidean")

    # Communities along a long gradient (Bray-Curtis), listed in full.
    pos = np.sort(rng.random(30)) * 10
    opt = np.linspace(-1, 11, 16)
    cnt = np.round(60 * np.exp(-((pos[:, None] - opt[None, :]) ** 2) / 2))
    num = np.abs(cnt[:, None, :] - cnt[None, :, :]).sum(axis=2)
    den = cnt.sum(axis=1)[:, None] + cnt.sum(axis=1)[None, :]
    d   = np.round(num / den, 6)
    write(outdir, "gradient", d, 2, n_neighbors=4)

    # Components joined by an edge of length 0 (samples 3 and 4, 1-based).
    d = np.array([
        [0, 2, 0, 1, 1, 5],
        [2, 0, 1, 0, 5, 2],
        [0, 1, 0, 0, 5, 5],
        [1, 0, 0, 0, 5, 5],
        [1, 5, 5, 5, 0, 6],
        [5, 2, 5, 5, 6, 0]], dtype=float)
    write(outdir, "zero_join", d, 1, n_neighbors=1,
          note="the joining edge has length 0")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data-raw/expected/isomap")
