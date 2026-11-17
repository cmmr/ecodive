"""Reference fixtures for ecodive's nmds(), from scikit-learn's non-metric SMACOF.

Writes, under data-raw/expected/nmds/ (data-raw/expected/nmds.r then
converts them to tests/tinytests/expected/nmds.rds):

  nmds-dist-<name>.txt   Dissimilarities, one per line, in the order of an R
                         `dist` object (the upper triangle, row by row).
                         Rounded to a few decimals, so that R and Python parse
                         them into identical doubles; rounding also makes ties.

  nmds-cases.csv         One row per case: its number (`id`), name,
                         dissimilarities, k, the number of iterations
                         scikit-learn ran, and the stress (as nmds()
                         measures it, with tied distances given equal
                         disparities) of scikit-learn's own best
                         configurations from 8 random starts (`skl_random`)
                         and from classical MDS (`skl_cmds`), for the
                         quality comparison.

  nmds-values.csv        Long format (id, field, index, value), with fields:
                           i  starting configuration (n x k, column-major)
                           s  raw stress after every 10th iteration, and
                              after the last (index = iteration)
                           x  configuration after the last iteration

Prerequisites: Python 3 with numpy and scikit-learn 1.9.1 (the version the
fixtures were made with). Run from the repository root:

    python data-raw/expected/nmds.py

scikit-learn's _smacof_single() is called unmodified. To record the stress of
every iteration, it is run with max_iter = 1, 2, 3, ...: the iterations are
deterministic, so the run with max_iter = t ends with iteration t's stress.
"""

import os

import numpy as np
import sklearn
from sklearn.manifold import ClassicalMDS
from sklearn.manifold._mds import _smacof_single, smacof

OUT = os.path.join("data-raw", "expected", "nmds")
MAX_ITER = 300
EPS = 1e-6


def condensed(D):
    """Upper triangle, row by row: the order of R's dist objects."""
    i, j = np.triu_indices(D.shape[0], 1)
    return D[i, j]


def square(v, n):
    D = np.zeros((n, n))
    i, j = np.triu_indices(n, 1)
    D[i, j] = v
    return D + D.T


def bray(counts):
    n = counts.shape[0]
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            x, y = counts[a], counts[b]
            D[a, b] = D[b, a] = np.abs(x - y).sum() / (x + y).sum()
    return D


def jaccard(pa):
    n = pa.shape[0]
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            x, y = pa[a] > 0, pa[b] > 0
            union = (x | y).sum()
            D[a, b] = D[b, a] = 1 - (x & y).sum() / union if union else 0
    return D


def datasets(rng):
    sets = {}

    # Random counts: few ties.
    for n in (30, 60):
        counts = rng.poisson(3, (n, 40)) * rng.binomial(1, 0.4, (n, 40))
        counts[:, 0] += 1
        sets["bray%d" % n] = (bray(counts), 6)

    # Presence/absence of 8 species: many tied dissimilarities, and
    # duplicate samples (zeros).
    pa = rng.binomial(1, 0.5, (50, 8))
    pa[:, 0] = 1
    pa[10] = pa[3]
    pa[20] = pa[3]
    sets["ties50"] = (jaccard(pa), 3)

    # A noisy line, for k = 1.
    t = np.sort(rng.uniform(0, 1, 40))
    pts = np.column_stack([t, rng.normal(0, 0.05, 40)])
    D = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    sets["line40"] = (D, 6)

    # Small, with duplicates and ties.
    counts = rng.poisson(2, (12, 6))
    counts[:, 0] += 1
    counts[7] = counts[2]
    sets["dup12"] = (bray(counts), 2)

    rounded = {}
    for name, (D, digits) in sets.items():
        v = np.round(condensed(D), digits)
        rounded[name] = (v, digits)
    return rounded


def stress_nmds(D, X):
    """Stress-1 as nmds(ties = 'secondary') measures it: the isotonic
    regression of the distances on the dissimilarities, tied
    dissimilarities sharing one disparity, zeros fixed at 0."""
    from sklearn.isotonic import isotonic_regression
    delta = condensed(D)
    d = condensed(np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)))
    dhat = np.zeros_like(d)
    nz = np.flatnonzero(delta != 0)
    o = nz[np.argsort(delta[nz], kind="stable")]
    vals, start, counts = np.unique(delta[o], return_index=True, return_counts=True)
    means = np.add.reduceat(d[o], start) / counts
    fit = isotonic_regression(means, sample_weight=counts)
    dhat[o] = np.repeat(fit, counts)
    return np.sqrt(((d - dhat) ** 2).sum() / (d ** 2).sum())


STRESS = stress_nmds


def fmt(x):
    return "%.17g" % x


def main():
    assert sklearn.__version__ == "1.9.1", sklearn.__version__
    rng = np.random.default_rng(20261117)
    os.makedirs(OUT, exist_ok=True)

    sets = datasets(rng)
    dists = {}
    for name, (v, digits) in sets.items():
        lines = ["%.*f" % (digits, x) for x in v]
        with open(os.path.join(OUT, "nmds-dist-%s.txt" % name), "w") as f:
            f.write("\n".join(lines) + "\n")
        # Parse the written text back, as R will.
        n = int(round((1 + np.sqrt(1 + 8 * len(v))) / 2))
        dists[name] = square(np.array([float(s) for s in lines]), n)

    cases = [
        ("bray30_k2", "bray30", 2),
        ("bray60_k2", "bray60", 2),
        ("bray60_k3", "bray60", 3),
        ("ties50_k2", "ties50", 2),
        ("ties50_k3", "ties50", 3),
        ("line40_k1", "line40", 1),
        ("dup12_k2", "dup12", 2),
    ]

    case_rows = ["id,case,dist,k,n_iter,skl_random,skl_cmds"]
    value_rows = ["id,field,index,value"]

    for cid, (case, dname, k) in enumerate(cases, 1):
        D = dists[dname]
        n = D.shape[0]
        init = np.array(["%.4f" % x for x in rng.uniform(0, 1, n * k)], dtype=float)
        init = init.reshape((n, k))

        X, _, n_iter = _smacof_single(
            D, metric=False, n_components=k, init=init.copy(),
            max_iter=MAX_ITER, eps=EPS, normalized_stress=False)

        trace = []
        for t in [t for t in range(1, n_iter + 1) if t % 10 == 0 or t == n_iter]:
            _, stress, it = _smacof_single(
                D, metric=False, n_components=k, init=init.copy(),
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
                D, metric=False, n_components=k, n_init=8, random_state=0,
                max_iter=MAX_ITER, eps=EPS)
            init_c = ClassicalMDS(metric="precomputed", n_components=k).fit_transform(D)
            Xc, _ = smacof(
                D, metric=False, n_components=k, init=init_c,
                max_iter=MAX_ITER, eps=EPS)
            skl = [fmt(STRESS(D, Xr)), fmt(STRESS(D, Xc))]

        case_rows.append("%d,%s,%s,%d,%d,%s,%s" % (cid, case, dname, k, n_iter, skl[0], skl[1]))

        print(case, "n_iter", n_iter)

    with open(os.path.join(OUT, "nmds-cases.csv"), "w") as f:
        f.write("\n".join(case_rows) + "\n")
    with open(os.path.join(OUT, "nmds-values.csv"), "w") as f:
        f.write("\n".join(value_rows) + "\n")


if __name__ == "__main__":
    main()
