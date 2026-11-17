"""Reference fixtures for ecodive's mds(), from scikit-learn's metric SMACOF.

Reads the dissimilarities that nmds.py wrote (run that first), and writes,
under data-raw/expected/mds/ (data-raw/expected/mds.r then converts them
to tests/tinytests/expected/mds.rds):

  mds-cases.csv    One row per case: its number (`id`), name,
                   dissimilarities, k, the number of iterations
                   scikit-learn ran, and the stress (as mds() measures it,
                   sqrt(raw stress / sum of squared dissimilarities)) of
                   scikit-learn's own best configurations from 8 random
                   starts (`skl_random`) and from classical MDS
                   (`skl_cmds`), for the quality comparison.

  mds-values.csv   Long format (id, field, index, value), with fields:
                     i  starting configuration (n x k, column-major)
                     s  raw stress after every 10th iteration, and after
                        the last (index = iteration)
                     x  configuration after the last iteration

Prerequisites: Python 3 with numpy and scikit-learn 1.9.1 (the version the
fixtures were made with). Run from the repository root:

    python data-raw/expected/nmds.py
    python data-raw/expected/mds.py

scikit-learn's _smacof_single() is called unmodified, with metric=True. To
record the stress of every iteration, it is run with max_iter = 1, 2, 3, ...:
the iterations are deterministic, so the run with max_iter = t ends with
iteration t's stress.
"""

import os

import numpy as np
import sklearn
from sklearn.manifold import ClassicalMDS
from sklearn.manifold._mds import _smacof_single, smacof

DIST = os.path.join("data-raw", "expected", "nmds")
OUT  = os.path.join("data-raw", "expected", "mds")
MAX_ITER = 300
EPS = 1e-6


def read_dist(name):
    with open(os.path.join(DIST, "nmds-dist-%s.txt" % name)) as f:
        v = np.array([float(s) for s in f.read().split()])
    n = int(round((1 + np.sqrt(1 + 8 * len(v))) / 2))
    D = np.zeros((n, n))
    i, j = np.triu_indices(n, 1)
    D[i, j] = v
    return D + D.T


def stress_mds(D, X):
    """Stress as mds() measures it: sqrt(raw stress / sum of squared
    dissimilarities)."""
    i, j = np.triu_indices(D.shape[0], 1)
    delta = D[i, j]
    d = np.sqrt(((X[i] - X[j]) ** 2).sum(-1))
    return np.sqrt(((d - delta) ** 2).sum() / (delta ** 2).sum())


STRESS = stress_mds


def fmt(x):
    return "%.17g" % x


def main():
    assert sklearn.__version__ == "1.9.1", sklearn.__version__
    rng = np.random.default_rng(20261118)

    cases = [
        ("bray30_k2", "bray30", 2),
        ("bray60_k3", "bray60", 3),
        ("ties50_k2", "ties50", 2),
        ("line40_k1", "line40", 1),
        ("dup12_k2", "dup12", 2),
    ]

    case_rows = ["id,case,dist,k,n_iter,skl_random,skl_cmds"]
    value_rows = ["id,field,index,value"]

    for cid, (case, dname, k) in enumerate(cases, 1):
        D = read_dist(dname)
        n = D.shape[0]
        init = np.array(["%.4f" % x for x in rng.uniform(0, 1, n * k)], dtype=float)
        init = init.reshape((n, k))

        X, _, n_iter = _smacof_single(
            D, metric=True, n_components=k, init=init.copy(),
            max_iter=MAX_ITER, eps=EPS, normalized_stress=False)

        trace = []
        for t in [t for t in range(1, n_iter + 1) if t % 10 == 0 or t == n_iter]:
            _, stress, it = _smacof_single(
                D, metric=True, n_components=k, init=init.copy(),
                max_iter=t, eps=EPS, normalized_stress=False)
            assert it == t
            trace.append((t, stress))

        skl = ["NA", "NA"]
        for idx, x in enumerate(init.T.ravel()):
            value_rows.append("%d,i,%d,%.4f" % (cid, idx + 1, x))
        for t, s in trace:
            value_rows.append("%d,s,%d,%s" % (cid, t, fmt(s)))
        for idx, x in enumerate(X.T.ravel()):
            value_rows.append("%d,x,%d,%s" % (cid, idx + 1, fmt(x)))

        if k >= 2:
            Xr, _ = smacof(
                D, metric=True, n_components=k, n_init=8, random_state=0,
                max_iter=MAX_ITER, eps=EPS)
            init_c = ClassicalMDS(metric="precomputed", n_components=k).fit_transform(D)
            Xc, _ = smacof(
                D, metric=True, n_components=k, init=init_c,
                max_iter=MAX_ITER, eps=EPS)
            skl = [fmt(STRESS(D, Xr)), fmt(STRESS(D, Xc))]

        case_rows.append("%d,%s,%s,%d,%d,%s,%s" % (cid, case, dname, k, n_iter, skl[0], skl[1]))

        print(case, "n_iter", n_iter)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "mds-cases.csv"), "w") as f:
        f.write("\n".join(case_rows) + "\n")
    with open(os.path.join(OUT, "mds-values.csv"), "w") as f:
        f.write("\n".join(value_rows) + "\n")


if __name__ == "__main__":
    main()
