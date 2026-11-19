"""Reference fixtures for ecodive's ancom().

Writes, under data-raw/expected/ancom/ (then data-raw/expected/ancom.r
converts them to tests/tinytests/expected/ancom.rds):

  table-<name>.csv  The input tables (samples x features): integer counts,
                    or hexadecimal floats (float.hex(), which R reads back
                    exactly) for real-valued tables.

  tables.csv        One row per table: name, pc (the pseudocount ancom()
                    uses by default: half the smallest non-zero value when
                    the table has zeros, else 0; hex float; scikit-bio is
                    run on x + pc), and groups (the samples' groups,
                    space-separated).

  results.csv       One row per run of scikit-bio's ancom(x + pc, groups):
                    table, sig_test ('anova' = "f_oneway", or 'kruskal'),
                    p_adjust (R's name; 'none' is p_adjust=None), alpha,
                    tau, theta (hex floats), W and Signif (space-separated,
                    Signif as 0/1).

  pairs.csv         Long format (table, sig_test, k, statistic, p_value):
                    for each pair of features i < j in the packed (row by
                    row) order, SciPy's statistic and p-value of the log
                    ratio log(x_i) - log(x_j), from _log_compare()'s inputs
                    and test, as hex floats.

  adjusted.csv      Long format (table, sig_test, p_adjust, i, j, value):
                    _adjust_pvalues() of _log_compare()'s p-value matrix,
                    row by row (axis = 1), as ancom() does, for every i != j
                    (1-based), as hex floats; tables of up to 15 features.

  cutoff.csv        scikit-bio's cutoff rule on crafted W: ancom() with
                    _log_compare() replaced by a p-value pattern with W[i]
                    values of 0.01 in row i (0.5 elsewhere) and
                    p_adjust=None, so that W is as chosen.
                    Columns D, tau, theta (hex floats), W and Signif
                    (space-separated, Signif as 0/1), and quirk (1 when
                    numpy.linspace's third cutoff, one rounding step above
                    0.15, changes Signif).

  dist.csv          SciPy's F and chi-squared upper tail probabilities
                    (special.fdtrc, special.chdtrc) on a grid: dist, df1,
                    df2 (0 for chi-squared), x, p.

Prerequisites: Python 3 with numpy, pandas, scipy and scikit-bio 0.7.4 (the
version the fixtures were made with). Run from the repository root:

    python data-raw/expected/ancom.py
    Rscript data-raw/expected/ancom.r
"""

import csv
import os
import warnings

import numpy as np
import pandas as pd
import scipy
import scipy.special
import scipy.stats
import skbio
from skbio.stats.composition import ancom
from skbio.stats.composition._ancom import _log_compare
from skbio.stats.composition._utils import _adjust_pvalues
from skbio.stats.composition.tests.test_ancom import AncomTests

assert skbio.__version__ == "0.7.4", skbio.__version__

warnings.simplefilter("ignore")   # SciPy's warnings about constant inputs

OUT = os.path.join("data-raw", "expected", "ancom")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261119)


# ---------------------------------------------------------------------------
# Tables: scikit-bio's own test tables, then ecodive's.
# ---------------------------------------------------------------------------

sk = AncomTests()
sk.setUp()

tables = {}
groups = {}

for i in (1, 2, 3, 4, 5, 6, 7, 9, 10):
    tables[f"skbio{i}"] = np.asarray(getattr(sk, f"table{i}"))
    groups[f"skbio{i}"] = [str(g) for g in getattr(sk, f"cats{i}")]

# The docstring example.
tables["doc"] = np.array([
    [12, 11, 10, 10, 10, 10, 10],
    [9,  11, 12, 10, 10, 10, 10],
    [1,  11, 10, 11, 10, 5,  9],
    [22, 21, 9,  10, 10, 10, 10],
    [20, 22, 10, 10, 13, 10, 10],
    [23, 21, 14, 10, 10, 10, 10]])
groups["doc"] = ["treatment"] * 3 + ["placebo"] * 3

# Three groups of unequal sizes, real-valued, some features shifted.
n3 = [7, 9, 8]
g3 = np.repeat([0, 1, 2], n3)
x = rng.lognormal(3, 1, size=(sum(n3), 15))
x[g3 == 1, :3] *= 4
x[g3 == 2, 3:5] *= 0.2
tables["three"] = x
groups["three"] = [["ctl", "lo", "hi"][g] for g in g3]

# Counts with zeros: ecodive's default pseudocount.
g2 = np.repeat([0, 1], [12, 12])
lam = np.tile(rng.uniform(0.5, 8, size=(1, 12)), (24, 1))
lam[g2 == 1, :2] *= 5
tables["zeros"] = rng.poisson(lam)
groups["zeros"] = [["A", "B"][g] for g in g2]

# Constant log ratios: features 0 and 1 identical (the same in every
# sample: NaN), feature 2 is 10 in one group and 20 in the other and
# feature 3 is 5 everywhere (constant within each group: F = Inf).
gc = np.repeat([0, 1], [5, 6])
x = rng.poisson(20, size=(11, 8)) + 1
x[:, 1] = x[:, 0]
x[:, 2] = np.where(gc == 0, 10, 20)
x[:, 3] = 5
tables["const"] = x
groups["const"] = [["a", "b"][g] for g in gc]

# Small integer counts: many ties for Kruskal-Wallis.
gt = np.repeat([0, 1, 2], [6, 6, 6])
x = rng.integers(1, 4, size=(18, 10))
x[gt == 2, 0] += 3
tables["ties"] = x
groups["ties"] = [["x", "y", "z"][g] for g in gt]

# Wide: more features, two groups, a block of shifted features.
gw = np.repeat([0, 1], [10, 12])
x = rng.lognormal(2, 0.7, size=(22, 30))
x[gw == 1, :6] *= 3
tables["wide"] = x
groups["wide"] = [["u", "v"][g] for g in gw]

# Singleton groups among others.
gs = np.array([0, 0, 0, 0, 1, 2, 2, 2])
x = rng.lognormal(2, 1, size=(8, 9))
tables["single"] = x
groups["single"] = [["p", "q", "r"][g] for g in gs]


def default_pc(tab):
    """ancom()'s pseudocount = NULL: half the smallest non-zero value."""
    return float(tab[tab > 0].min()) / 2 if (tab == 0).any() else 0.0


def hexf(v):
    return float(v).hex()


def write_table(name, tab):
    with open(os.path.join(OUT, f"table-{name}.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        if tab.dtype.kind == "f":
            w.writerows([[hexf(v) for v in row] for row in tab])
        else:
            w.writerows(tab.tolist())


TESTS = {"anova": scipy.stats.f_oneway, "kruskal": scipy.stats.kruskal}
ADJUST = {"none": None, "holm": "holm", "bonferroni": "bonferroni",
          "BH": "bh", "BY": "by"}

table_rows, result_rows, pair_rows, adj_rows = [], [], [], []

for name, tab in tables.items():
    write_table(name, tab)
    pc = default_pc(tab)
    table_rows.append([name, hexf(pc), " ".join(groups[name])])

    x = tab.astype(float) + pc
    grp = np.asarray(groups[name])
    uniq, labels = np.unique(grp, return_inverse=True)
    D = x.shape[1]
    ii, jj = np.triu_indices(D, k=1)

    for st, test in TESTS.items():

        # The statistics and p-values of each pair, as _log_compare() runs
        # the test.
        L = np.log(x)
        ratios = [L[labels == g][:, ii] - L[labels == g][:, jj]
                  for g in range(len(uniq))]
        stat, pv = test(*ratios)
        pmat = _log_compare(x, labels, len(uniq), test)
        assert np.array_equal(pmat[ii, jj], pv, equal_nan=True)
        for k in range(len(ii)):
            pair_rows.append([name, st, k + 1, hexf(stat[k]), hexf(pv[k])])

        for r_adj, sk_adj in ADJUST.items():
            adj = _adjust_pvalues(pmat.copy(), sk_adj, axis=1)
            for i in range(D if D <= 15 else 0):   # smaller tables only
                for j in range(D):
                    if i != j:
                        adj_rows.append([name, st, r_adj, i + 1, j + 1, hexf(adj[i, j])])

            params = [(0.05, 0.02, 0.1)]
            if name in ("skbio1", "three", "wide") and r_adj == "holm":
                params += [(0.1, 0.02, 0.1), (0.05, 0.1, 0.1),
                           (0.05, 0.02, 0.5), (0.01, 0.2, 0.3)]
            for alpha, tau, theta in params:
                res, _ = ancom(pd.DataFrame(x), pd.Series(grp), alpha=alpha,
                               tau=tau, theta=theta, p_adjust=sk_adj,
                               sig_test="f_oneway" if st == "anova" else st,
                               percentiles=[])
                result_rows.append([
                    name, st, r_adj, hexf(alpha), hexf(tau), hexf(theta),
                    " ".join(str(int(v)) for v in res["W"]),
                    " ".join(str(int(v)) for v in res["Signif"])])


def write(fname, header, rows):
    with open(os.path.join(OUT, fname), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


write("tables.csv", ["table", "pc", "groups"], table_rows)
write("results.csv", ["table", "sig_test", "p_adjust", "alpha", "tau",
                      "theta", "W", "signif"], result_rows)
write("pairs.csv", ["table", "sig_test", "k", "statistic", "p_value"], pair_rows)
write("adjusted.csv", ["table", "sig_test", "p_adjust", "i", "j", "value"], adj_rows)


# ---------------------------------------------------------------------------
# The cutoff rule on crafted W.
# ---------------------------------------------------------------------------

import skbio.stats.composition._ancom as _ancom_mod


def signif_with(W, D, tau, theta, steps):
    c_start = W.max() / D
    if c_start < theta:
        return np.zeros(D, dtype=bool)
    cutoff = c_start - steps
    prop_cut = (W[:, None] > D * cutoff).mean(axis=0)
    dels = np.abs(prop_cut - np.roll(prop_cut, -1))
    dels[-1] = 0
    if (dels[0] < tau) and (dels[1] < tau) and (dels[2] < tau):
        nu = cutoff[1]
    elif (dels[0] >= tau) and (dels[1] < tau) and (dels[2] < tau):
        nu = cutoff[2]
    elif (dels[1] >= tau) and (dels[2] < tau) and (dels[3] < tau):
        nu = cutoff[3]
    else:
        nu = cutoff[4]
    return W >= nu * D


cut_rows = []
orig_log_compare = _ancom_mod._log_compare
crng = np.random.default_rng(7)
try:
    for case in range(800):
        D = int(crng.choice([7, 13, 20, 40, 60, 100]))
        # W near the top (a few steps of D / 20 below the largest), and low.
        M = int(crng.integers(D // 4, D))
        top = M - crng.choice(np.arange(0, D // 4 + 1), size=D,
                              p=None if case % 2 else
                              np.ones(D // 4 + 1) / (D // 4 + 1))
        low = crng.integers(0, max(1, M // 3), size=D)
        W = np.where(crng.random(D) < crng.uniform(0.1, 0.9), top, low)
        W = np.clip(W, 0, D - 1)
        # A p-value pattern with W[i] values below alpha in row i, off the
        # diagonal (it need not be symmetric).
        pmat = np.full((D, D), 0.5)
        for i in range(D):
            others = np.delete(np.arange(D), i)
            pmat[i, crng.choice(others, size=W[i], replace=False)] = 0.01
        np.fill_diagonal(pmat, 0)
        _ancom_mod._log_compare = lambda *a, pm=pmat: pm.copy()
        tau = float(crng.choice([0.02, 0.05, 0.1, 0.3]))
        theta = float(crng.choice([0.1, 0.1, 0.3]))
        res, _ = ancom(pd.DataFrame(np.ones((4, D))), pd.Series([0, 0, 1, 1]),
                       tau=tau, theta=theta, p_adjust=None, percentiles=[])
        assert np.array_equal(res["W"].to_numpy(), W)
        sig = res["Signif"].to_numpy()
        plain = signif_with(W, D, tau, theta, np.array([0.05, 0.1, 0.15, 0.2, 0.25]))
        cut_rows.append([D, hexf(tau), hexf(theta),
                         " ".join(str(int(v)) for v in W),
                         " ".join(str(int(v)) for v in sig),
                         int(not np.array_equal(sig, plain))])
finally:
    _ancom_mod._log_compare = orig_log_compare

write("cutoff.csv", ["D", "tau", "theta", "W", "signif", "quirk"], cut_rows)
print("cutoff cases:", len(cut_rows), "where the quirk matters:",
      sum(r[-1] for r in cut_rows))


# ---------------------------------------------------------------------------
# The distributions on a grid.
# ---------------------------------------------------------------------------

dist_rows = []
xs = [0.0, 1e-8, 0.01, 0.3, 1.0, 2.5, 7.0, 15.0, 40.0, 120.0, 800.0]
for df1 in (1, 2, 3, 5):
    for df2 in (1, 2, 4, 9, 30, 118, 1000):
        for v in xs:
            dist_rows.append(["f", df1, df2, hexf(v),
                              hexf(scipy.special.fdtrc(df1, df2, v))])
    for v in xs:
        dist_rows.append(["chisq", df1, 0, hexf(v),
                          hexf(scipy.special.chdtrc(df1, v))])
write("dist.csv", ["dist", "df1", "df2", "x", "p"], dist_rows)

print("scikit-bio", skbio.__version__, "scipy", scipy.__version__,
      "numpy", np.__version__)
