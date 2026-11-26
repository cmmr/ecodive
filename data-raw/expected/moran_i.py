# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for moran_i() from PySAL esda.

Writes inputs and esda's Moran (esda 2.10.0, libpysal 4.15.1, both
BSD-3-Clause) of them to data-raw/expected/moran_i/. Then
data-raw/expected/moran_i.r folds them (with ape's and spdep's outputs)
into tests/tinytests/expected/moran_i.rds, so that the tests need no
Python. Run from the package root:

    python data-raw/expected/moran_i.py
    Rscript data-raw/expected/moran_i.r

Cases (each with its own spatial weights, written as a full matrix with a
zero diagonal, read by libpysal's full2W()):

  idw       30 random positions in a 10 x 10 square; inverse-distance
            weights 1 / d; x has a trend along the first coordinate.
  knn       25 random positions; inverse-distance weights to each
            sample's 4 nearest neighbours (not symmetric).
  island    12 samples with random symmetric weights, where sample 5 has
            none (an island).
  checker   A 4 x 4 grid with rook neighbours and alternating values:
            negative autocorrelation.
  small     4 samples with unequal raw weights.

For each case, and each of esda's transformations "r" (row-standardised),
"b" (binary), "v" (variance-stabilising) and "o" (as given), esda's Moran
is run with two_tailed=False and 199 permutations, and again with
two_tailed=True and none. The permutations are this script's own (written
1-based to perms-<case>.csv): numpy.random.permutation is replaced, while
esda runs, by a function that returns z[perm] for each of them in turn, so
that esda's own loop computes the permuted statistics (sims-<case>.csv)
and its p_sim from them, and ecodive can compute the same ones.

Every real number is written as a hexadecimal float (float.hex()).
"""

import csv
import os

import esda
import libpysal
import numpy as np
from libpysal.weights import full2W

assert esda.__version__ == "2.10.0", esda.__version__
assert libpysal.__version__ == "4.15.1", libpysal.__version__

OUT = os.path.join("data-raw", "expected", "moran_i")
os.makedirs(OUT, exist_ok=True)

N_PERM = 199
STYLES = ["r", "b", "v", "o"]


def hx(v):
    return float(v).hex()


def pairwise(xy):
    d = xy[:, None, :] - xy[None, :, :]
    return np.sqrt((d * d).sum(2))


rng = np.random.default_rng(20261125)
cases = {}

# idw
xy = rng.uniform(0, 10, (30, 2))
d = pairwise(xy)
w = np.zeros_like(d)
off = ~np.eye(30, dtype=bool)
w[off] = 1 / d[off]
x = xy[:, 0] + rng.normal(0, 2, 30)
cases["idw"] = (x, w, xy)

# knn
xy = rng.uniform(0, 10, (25, 2))
d = pairwise(xy)
np.fill_diagonal(d, np.inf)
w = np.zeros((25, 25))
for i in range(25):
    nn = np.argsort(d[i], kind="stable")[:4]
    w[i, nn] = 1 / d[i, nn]
x = 0.5 * xy[:, 1] + rng.normal(0, 1, 25)
cases["knn"] = (x, w, xy)

# island
w = rng.uniform(0, 1, (12, 12))
w = (w + w.T) / 2
w[w < 0.5] = 0
np.fill_diagonal(w, 0)
w[4, :] = 0
w[:, 4] = 0
x = rng.normal(5, 2, 12)
cases["island"] = (x, w, None)

# checker
g = np.arange(16)
r, c = g // 4, g % 4
w = ((np.abs(r[:, None] - r[None, :]) + np.abs(c[:, None] - c[None, :])) == 1).astype(float)
x = np.where((r + c) % 2 == 0, 1.0, -1.0) + rng.normal(0, 0.3, 16)
cases["checker"] = (x, w, np.column_stack([r, c]).astype(float))

# small
w = np.array([[0, 2, 1, 0.5], [2, 0, 3, 0], [1, 3, 0, 1], [0.5, 0, 1, 0]])
x = np.array([1.5, 2.0, 4.5, 3.0])
cases["small"] = (x, w, None)


def run_moran(x, w, style, perms, two_tailed):
    W = full2W(w)
    W.silence_warnings = True
    it = iter(perms)
    orig = np.random.permutation
    np.random.permutation = lambda z: z[next(it)]
    try:
        m = esda.Moran(x, W, transformation=style,
                       permutations=len(perms), two_tailed=two_tailed)
    finally:
        np.random.permutation = orig
    return m, W


stat_rows = []
for name, (x, w, xy) in cases.items():
    n = len(x)
    perms = [rng.permutation(n) for _ in range(N_PERM)]

    with open(os.path.join(OUT, f"x-{name}.csv"), "w", newline="") as fh:
        csv.writer(fh).writerows([[hx(v)] for v in x])
    with open(os.path.join(OUT, f"w-{name}.csv"), "w", newline="") as fh:
        csv.writer(fh).writerows([[hx(v) for v in row] for row in w])
    if xy is not None:
        with open(os.path.join(OUT, f"xy-{name}.csv"), "w", newline="") as fh:
            csv.writer(fh).writerows([[hx(v) for v in row] for row in xy])
    with open(os.path.join(OUT, f"perms-{name}.csv"), "w", newline="") as fh:
        csv.writer(fh).writerows([[int(v) + 1 for v in p] for p in perms])

    sims = []
    for style in STYLES:
        m, W = run_moran(x, w, style, perms, two_tailed=False)
        m2, _ = run_moran(x, w, style, [], two_tailed=True)
        assert m2.I == m.I
        stat_rows.append([
            name, style, n,
            hx(m.I), hx(m.EI), hx(m.VI_norm), hx(m.VI_rand),
            hx(m.z_norm), hx(m.z_rand),
            hx(m.p_norm), hx(m.p_rand), hx(m2.p_norm), hx(m2.p_rand),
            hx(m.p_sim), hx(m.EI_sim), hx(m.seI_sim), hx(m.z_sim),
            hx(W.s0), hx(W.s1), hx(W.s2)])
        sims.append([hx(v) for v in m.sim])

    with open(os.path.join(OUT, f"sims-{name}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(STYLES)
        wr.writerows(zip(*sims))

with open(os.path.join(OUT, "esda.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow([
        "case", "style", "n", "I", "EI", "VI_norm", "VI_rand", "z_norm", "z_rand",
        "p_norm", "p_rand", "p_norm2", "p_rand2", "p_sim", "EI_sim", "seI_sim",
        "z_sim", "s0", "s1", "s2"])
    wr.writerows(stat_rows)
