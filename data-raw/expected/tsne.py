"""Reference values for ecodive's tsne() from scikit-learn's t-SNE.

Writes the frozen scikit-learn reference values that data-raw/expected/tsne.r
stores in tests/tinytests/expected/tsne.rds for test-tsne.r, so that the
tests need neither Python nor scikit-learn. The committed outputs are the
source of truth: rerun this only to regenerate them.

Prerequisites: Python 3 with scikit-learn 1.9.1 (and its numpy and
scipy). From the repository root:

    python data-raw/expected/tsne.py
    Rscript data-raw/expected/tsne.r

Output: data-raw/expected/tsne/<dataset>.txt, one per dataset.
Each line is a name followed by whitespace-separated numbers.

Inputs are written first, rounded to 7 significant digits, and then
read back, so that both sides compute from exactly the same doubles.
Every reference value is computed from those rounded inputs.

  dist          Packed distances (scipy's condensed order, which is
                also R's dist order).
  counts        (bray only) The count table, row-major, and its
                dimensions as counts_dim. The distances are its
                Bray-Curtis dissimilarities.
  labels        Cluster labels, for the quality checks.
  perp          Perplexity for the affinity checks.
  P_dense       _joint_probabilities(D^2, perp), condensed.
  perp_nn       Perplexity for the sparse affinities.
  P_ptr/col/val _joint_probabilities_nn() on the 3 * perp_nn + 1
                nearest neighbors (CSR, 0-based), as TSNE builds it.
  Y<k>          A fixed embedding in k = 1, 2, 3 dimensions.
  exact<k>_kl, exact<k>_grad
                _kl_divergence(Y<k>, P_dense, dof = 1).
  bh<k>_kl, bh<k>_grad
                _kl_divergence_bh(Y<k>, P_sparse, dof = 1, angle = 0.5).
  bh2_t<theta>_grad
                The same, k = 2, for other angles.
  Y0            A starting embedding (k = 2) for the optimizer.
  traj<m>       Positions after m = 10 and 25 iterations of
                _gradient_descent() from Y0 with exact gradients: the
                early-exaggeration phase (P * 12, momentum 0.5,
                learning rate "auto").
  p2_start, p2_30
                A spread-out layout, and the positions after 30
                iterations of the second phase from it (P, momentum
                0.8, iterations 250 to 279), as _tsne() runs it.
  q_<init>_<stat>
                Quality of TSNE(metric = "precomputed", method =
                "exact", perplexity = qperp) over seeds 0-9:
                final KL divergence (kl), trustworthiness with 5
                neighbors (trust), and leave-one-out 5-nearest-
                neighbor label accuracy in the embedding (knn).
                <init> is "random" (scikit-learn's start for
                precomputed distances) or "pcoa" (classical scaling,
                scaled to standard deviation 1e-4, as ecodive does).
                Ten values each, by seed.
"""

import os

import numpy as np
import scipy.sparse as sp
import sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE, trustworthiness
from sklearn.manifold import _t_sne
from sklearn.neighbors import NearestNeighbors

assert sklearn.__version__ == "1.9.1", sklearn.__version__

OUT = os.path.join("data-raw", "expected", "tsne")


def r7(x):
    """Round to 7 significant digits, through the decimal string."""
    return np.array([float("%.7g" % v) for v in np.ravel(x)]).reshape(np.shape(x))


def fmt(v, digits):
    return " ".join(("%." + str(digits) + "g") % x for x in np.ravel(v))


class Writer:
    def __init__(self, name):
        self.path = os.path.join(OUT, "%s.txt" % name)
        self.lines = []

    def put(self, key, values, digits=12):
        self.lines.append(key + " " + fmt(values, digits))

    def close(self):
        with open(self.path, "w") as f:
            f.write("\n".join(self.lines) + "\n")
        print(self.path, os.path.getsize(self.path), "bytes")


def joint_nn(D, perp):
    n = D.shape[0]
    K = min(n - 1, int(3.0 * perp + 1))
    knn = NearestNeighbors(n_neighbors=K, metric="precomputed").fit(D)
    g = knn.kneighbors_graph(mode="distance")
    g.data **= 2
    return _t_sne._joint_probabilities_nn(g, perp, 0)


def classical_scaling(D, k):
    n = D.shape[0]
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (D**2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1][:k]
    Y = vecs[:, order] * np.sqrt(np.maximum(vals[order], 0))
    return Y / np.std(Y[:, 0]) * 1e-4


def knn_accuracy(Y, labels, m=5):
    D = squareform(pdist(Y))
    np.fill_diagonal(D, np.inf)
    hits = 0
    for i in range(len(labels)):
        nb = np.argsort(D[i], kind="stable")[:m]
        vals, counts = np.unique(labels[nb], return_counts=True)
        hits += vals[np.argmax(counts)] == labels[i]
    return hits / len(labels)


def dataset(name, d, labels, perp, perp_nn, qperp, counts=None):

    w = Writer(name)
    n = squareform(d).shape[0]

    w.put("dist", d, 7)
    if counts is not None:
        w.put("counts", counts, 7)
        w.put("counts_dim", counts.shape, 7)
    w.put("labels", labels, 7)
    w.put("perp", [perp])

    D = squareform(d)

    # Affinities.
    P = _t_sne._joint_probabilities(D**2, perp, 0)
    w.put("P_dense", P, 10)
    w.put("perp_nn", [perp_nn])
    Psp = sp.csr_matrix(joint_nn(D, perp_nn))
    Psp.sort_indices()
    w.put("P_ptr", Psp.indptr, 10)
    w.put("P_col", Psp.indices, 10)
    w.put("P_val", Psp.data, 10)

    # Objective and gradient at fixed embeddings.
    rng = np.random.default_rng(7)
    for k in (1, 2, 3):
        Y = r7(rng.standard_normal((n, k)) * 3)
        w.put("Y%d" % k, Y, 7)
        kl, g = _t_sne._kl_divergence(Y.ravel(), P, 1, n, k)
        w.put("exact%d_kl" % k, [kl])
        w.put("exact%d_grad" % k, g.reshape(n, k))
        kl, g = _t_sne._kl_divergence_bh(Y.ravel(), Psp, 1, n, k, angle=0.5)
        w.put("bh%d_kl" % k, [kl])
        w.put("bh%d_grad" % k, g.reshape(n, k), 9)
        if k == 2:
            for t in (0.0, 0.2, 0.8):
                _, g = _t_sne._kl_divergence_bh(Y.ravel(), Psp, 1, n, k, angle=t)
                w.put("bh2_t%g_grad" % t, g.reshape(n, k), 9)

    # Optimizer trajectory, exact gradients, double precision.
    Y0 = r7(rng.standard_normal((n, 2)) * 1e-4)
    w.put("Y0", Y0, 7)
    lr = max(n / 12.0 / 4, 50)
    for m in (10, 25):
        p, err, it = _t_sne._gradient_descent(
            _t_sne._kl_divergence, Y0.ravel().copy(), 0, m, n_iter_check=50,
            n_iter_without_progress=250, momentum=0.5, learning_rate=lr,
            min_gain=0.01, min_grad_norm=1e-7, args=[P * 12.0, 1, n, 2], kwargs={})
        w.put("traj%d" % m, p.reshape(n, 2))

    # The second phase (no exaggeration, momentum 0.8, gains and update
    # reset), from a spread-out layout, as _tsne() runs it after
    # iteration 250.
    Y2 = r7(rng.standard_normal((n, 2)) * 5)
    w.put("p2_start", Y2, 7)
    p, err, it = _t_sne._gradient_descent(
        _t_sne._kl_divergence, Y2.ravel().copy(), 250, 280, n_iter_check=50,
        n_iter_without_progress=300, momentum=0.8, learning_rate=lr,
        min_gain=0.01, min_grad_norm=1e-7, args=[P, 1, n, 2], kwargs={})
    w.put("p2_30", p.reshape(n, 2))

    # Quality over seeds.
    w.put("qperp", [qperp])
    starts = {"random": None, "pcoa": classical_scaling(D, 2)}
    for init, Yinit in starts.items():
        kl, tr, kn = [], [], []
        for seed in range(10):
            t = TSNE(
                n_components=2, perplexity=qperp, metric="precomputed",
                method="exact", random_state=seed,
                init="random" if Yinit is None else Yinit.copy())
            Y = t.fit_transform(D.copy())  # It squares D in place.
            kl.append(t.kl_divergence_)
            tr.append(trustworthiness(D, Y, n_neighbors=5, metric="precomputed"))
            kn.append(knn_accuracy(Y, labels))
            if Yinit is not None:
                break  # Deterministic: the seed is unused.
        w.put("q_%s_kl" % init, kl, 7)
        w.put("q_%s_trust" % init, tr, 7)
        w.put("q_%s_knn" % init, kn, 7)

    w.close()


def main():
    os.makedirs(OUT, exist_ok=True)

    # Three well-separated Gaussian clusters in 6 dimensions.
    rng = np.random.default_rng(2026)
    centers = rng.standard_normal((3, 6)) * 6
    labels = np.repeat(np.arange(3), 15)
    X = centers[labels] + rng.standard_normal((45, 6))
    d = r7(pdist(X))
    dataset("blobs", d, labels, perp=10.0, perp_nn=4.0, qperp=10.0)

    # Microbiome-like counts in two groups: negative binomial counts
    # around group-specific profiles. Bray-Curtis dissimilarities.
    rng = np.random.default_rng(11)
    n, p = 40, 25
    labels = np.repeat(np.arange(2), n // 2)
    base = rng.gamma(0.6, 1.0, size=(2, p))
    base[1] = base[1] * rng.gamma(2.0, 0.5, size=p)
    prof = base[labels] / base[labels].sum(axis=1, keepdims=True)
    depth = rng.integers(500, 3000, size=n)
    mu = prof * depth[:, None]
    counts = rng.negative_binomial(2, 2 / (2 + mu))
    counts[:, 0] += 1  # No empty samples.
    d = r7(pdist(counts, "braycurtis"))
    dataset("bray", d, labels, perp=5.0, perp_nn=3.0, qperp=8.0, counts=counts)


if __name__ == "__main__":
    main()
