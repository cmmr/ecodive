"""Reference fixtures for ecodive's dirmult_ttest().

scikit-bio's draws come from NumPy's PCG64 Gamma sampler, which ecodive
can't reproduce, so the deterministic parts are compared on the same
draws: each case's Gamma draws (NumPy, fixed seed) are written out, and
scikit-bio's dirmult_ttest() is run with its random generator replaced by
one that hands out those draws in turn. ecodive's tests feed the same
draws through dirmult_ttest_details(gamma_draws =).

Writes, under data-raw/expected/dirmult_ttest/ (then
data-raw/expected/dirmult_ttest.r converts them to
tests/tinytests/expected/dirmult_ttest.rds):

  table-<name>.csv  The input tables (samples x features), integer counts.

  cases.csv         One row per case: case, table, groups (the samples'
                    groups, space-separated), treatment and reference
                    (scikit-bio's arguments; empty for None), n_draws.

  gamma-<case>.csv  The case's Gamma draws, standard_gamma(table + 0.5):
                    columns draw, sample (0-based), then one hex float
                    (float.hex(), which R reads back exactly) per feature.

  results.csv       One row per run of scikit-bio's dirmult_ttest() on a
                    case's draws: case, p_adjust (R's name; 'none' is
                    p_adjust=None), then t, log2_fc, ci_lower, ci_upper,
                    p_value, p_adjusted (space-separated hex floats) and
                    reject (space-separated 0/1).

  perdraw.csv       Each draw's Welch statistics (scikit-bio's
                    _dirmult_draw() and _welch_draw_stats() on the draw,
                    then SciPy's stdtr() and stdtrit() as dirmult_ttest()
                    uses them): case, draw, feature (0-based), diff, se,
                    dof, t, p_tail (one-sided), lower, upper (natural logs,
                    95%), as hex floats.

  sampled.csv       scikit-bio's own run (its own draws, seed = 0, 128
                    draws) on the large-count table `deep`: table, groups,
                    treatment, reference, then the columns as in
                    results.csv, for the comparison by distribution.

  tdist.csv         SciPy's stdtr(df, x) and stdtrit(df, p) on a grid:
                    what ('stdtr' or 'stdtrit'), df, x (x or p), value.

Prerequisites: Python 3 with numpy, pandas, scipy and scikit-bio 0.7.4 (the
version the fixtures were made with). Run from the repository root:

    python data-raw/expected/dirmult_ttest.py
    Rscript data-raw/expected/dirmult_ttest.r
"""

import csv
import os
import warnings

import numpy as np
import pandas as pd
import scipy
import scipy.special
import skbio
import skbio.stats.composition._dirmult as dm
from skbio.stats.composition import dirmult_ttest

assert skbio.__version__ == "0.7.4", skbio.__version__

warnings.simplefilter("ignore")   # NumPy's warnings for a group of one

OUT = os.path.join("data-raw", "expected", "dirmult_ttest")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261119)
PC = 0.5


def hx(v):
    return " ".join(float(x).hex() for x in np.asarray(v, dtype=float))


# ---------------------------------------------------------------------------
# Tables.
# ---------------------------------------------------------------------------

tables = {}
groups = {}

# The docstring example (and test_dirmult_ttest_demo).
tables["doc"] = np.array([
    [20,  110, 100, 101, 100, 103, 104],
    [33,  110, 120, 100, 101, 100, 102],
    [12,  110, 100, 110, 100, 50,  90],
    [202, 201, 9,   10,  10,  11,  11],
    [200, 202, 10,  10,  13,  10,  10],
    [203, 201, 14,  10,  10,  13,  12]])
groups["doc"] = ["treatment"] * 3 + ["placebo"] * 3

# scikit-bio's small test table.
tables["small"] = np.array([
    [5, 3, 10], [8, 6, 15], [12, 9, 20], [15, 12, 25], [20, 15, 30]])
groups["small"] = ["Group1", "Group1", "Group2", "Group2", "Group2"]

# Three groups of unequal sizes, Poisson counts with zeros, some features
# shifted.
n3 = [5, 7, 6]
g3 = np.repeat([0, 1, 2], n3)
lam = np.tile(rng.lognormal(1.5, 1.2, size=(1, 12)), (sum(n3), 1))
lam[g3 == 1, :3] *= 4
lam[g3 == 2, 3:5] *= 0.25
tables["three"] = rng.poisson(lam)
groups["three"] = [["ctl", "lo", "hi"][g] for g in g3]

# A group of one sample: NaN statistics.
tables["single"] = rng.poisson(20, size=(6, 5))
groups["single"] = ["a", "b", "b", "b", "b", "b"]

# Wider and unbalanced, sparse.
gw = np.repeat([0, 1], [9, 14])
lam = np.tile(rng.lognormal(0.5, 1.5, size=(1, 25)), (23, 1))
lam[gw == 1, :5] *= 3
tables["wide"] = rng.poisson(lam)
groups["wide"] = [["x", "y"][g] for g in gw]

# Deep counts, for the comparison by distribution: little Dirichlet noise.
gd = np.repeat([0, 1], [10, 10])
p = rng.lognormal(0, 1, size=12)
p /= p.sum()
P = np.tile(p, (20, 1))
P[gd == 1, :3] *= 2.5
P[gd == 1, 3] *= 0.4
P /= P.sum(axis=1, keepdims=True)
tables["deep"] = np.vstack([rng.multinomial(200000, P[i]) for i in range(20)])
groups["deep"] = [["ctl", "trt"][g] for g in gd]

for nm, t in tables.items():
    np.savetxt(os.path.join(OUT, f"table-{nm}.csv"), t, fmt="%d", delimiter=",")


# ---------------------------------------------------------------------------
# Cases: a table, a contrast and its Gamma draws.
# ---------------------------------------------------------------------------

cases = [
    # case, table, treatment, reference, draws
    ("doc",         "doc",    "treatment", "placebo", 16),
    ("doc_default", "doc",    None,        None,      8),
    ("small",       "small",  "Group2",    "Group1",  8),
    ("three_pair",  "three",  "hi",        "lo",      8),
    ("three_pool",  "three",  "lo",        None,      8),
    ("three_dflt",  "three",  None,        None,      8),
    ("single",      "single", "a",         "b",       4),
    ("wide",        "wide",   "y",         "x",       12),
]

draws = {}
with open(os.path.join(OUT, "cases.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["case", "table", "groups", "treatment", "reference", "n_draws"])
    for case, tab, trt, ref, nd in cases:
        w.writerow([case, tab, " ".join(groups[tab]), trt or "", ref or "", nd])
        mat = tables[tab] + PC
        G = np.stack([rng.standard_gamma(mat) for _ in range(nd)])
        draws[case] = G
        with open(os.path.join(OUT, f"gamma-{case}.csv"), "w", newline="") as g:
            wg = csv.writer(g, lineterminator="\n")
            wg.writerow(["draw", "sample"] + [f"f{j}" for j in range(mat.shape[1])])
            for d in range(nd):
                for s in range(mat.shape[0]):
                    wg.writerow([d, s] + [float(v).hex() for v in G[d, s]])


class GivenDraws:
    """A stand-in for a NumPy Generator whose standard_gamma() returns the
    given draws in turn."""

    def __init__(self, G):
        self.G = G
        self.i = 0

    def standard_gamma(self, shape, size=None, out=None):
        g = self.G[self.i]
        self.i += 1
        assert g.shape == np.shape(shape)
        if out is None:
            return g.copy()
        out[...] = g
        return out


def run_given(case, tab, trt, ref, p_adjust):
    G = draws[case]
    orig = dm.get_rng
    dm.get_rng = lambda seed: GivenDraws(G)
    try:
        res = dirmult_ttest(
            pd.DataFrame(tables[tab]), pd.Series(groups[tab]), trt, ref,
            pseudocount=PC, draws=G.shape[0], p_adjust=p_adjust)
    finally:
        dm.get_rng = orig
    return res


def row_of(res):
    return [
        hx(res["T-statistic"]), hx(res["Log2(FC)"]), hx(res["CI(2.5)"]),
        hx(res["CI(97.5)"]), hx(res["pvalue"]), hx(res["qvalue"]),
        " ".join(str(int(v)) for v in res["Signif"])]


PADJ = [("holm", "holm"), ("bonferroni", "bonferroni"), ("BH", "bh"),
        ("BY", "by"), ("none", None)]

with open(os.path.join(OUT, "results.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["case", "p_adjust", "t", "log2_fc", "ci_lower", "ci_upper",
                "p_value", "p_adjusted", "reject"])
    for case, tab, trt, ref, nd in cases:
        for r_name, sk_name in PADJ:
            res = run_given(case, tab, trt, ref, sk_name)
            w.writerow([case, r_name] + row_of(res))


# Each draw's statistics, from scikit-bio's own pieces.
from scipy.special import stdtr, stdtrit
from skbio.stats.composition._utils import _check_grouping, _check_trt_ref_groups

with open(os.path.join(OUT, "perdraw.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["case", "draw", "feature", "diff", "se", "dof", "t",
                "p_tail", "lower", "upper"])
    for case, tab, trt, ref, nd in cases:
        mat = tables[tab] + PC
        grp, lab = _check_grouping(np.asarray(groups[tab]), mat)
        ti, ri = _check_trt_ref_groups(trt, ref, grp, lab)
        m = mat.shape[1]
        gen = GivenDraws(draws[case])
        for d in range(nd):
            z = dm._dirmult_draw(mat, gen)
            diff, se, dof, work = (np.empty(m) for _ in range(4))
            dm._welch_draw_stats(z[ti], z[ri], len(ti), len(ri), diff, se, dof, work)
            t = diff / se
            pt = stdtr(dof, -np.abs(t))
            q = stdtrit(dof, 0.975) * se
            lo, hi = diff - q, diff + q
            for j in range(m):
                w.writerow([case, d, j] + [float(v).hex() for v in
                           (diff[j], se[j], dof[j], t[j], pt[j], lo[j], hi[j])])


# scikit-bio with its own draws, on the deep table.
res = dirmult_ttest(
    pd.DataFrame(tables["deep"]), pd.Series(groups["deep"]), "trt", "ctl",
    seed=0)
with open(os.path.join(OUT, "sampled.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["table", "groups", "treatment", "reference", "t", "log2_fc",
                "ci_lower", "ci_upper", "p_value", "p_adjusted", "reject"])
    w.writerow(["deep", " ".join(groups["deep"]), "trt", "ctl"] + row_of(res))


# SciPy's t distribution on a grid.
with open(os.path.join(OUT, "tdist.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["what", "df", "x", "value"])
    dfs = np.concatenate([[1.0, 1.5, 2.0, 2.37, 3.9, 7.25, 15.5, 33.3, 120.7],
                          rng.uniform(1, 60, size=40)])
    for df in dfs:
        for x in (-40.0, -12.5, -6.0, -3.1, -2.0, -1.0, -0.37, -0.01, 0.0):
            w.writerow(["stdtr", df.hex(), x.hex(), float(stdtr(df, x)).hex()])
        for p in (0.95, 0.975, 0.995, 0.9995):
            w.writerow(["stdtrit", df.hex(), p.hex(), float(stdtrit(df, p)).hex()])

print("numpy", np.__version__, "scipy", scipy.__version__, "pandas", pd.__version__)
