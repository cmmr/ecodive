"""Reference fixtures for ecodive's vlr().

Writes, under data-raw/expected/proportionality/ (then
data-raw/expected/proportionality.r converts them to
tests/tinytests/expected/proportionality.rds):

  table-<name>.csv  The input tables (samples x features): integer counts,
                    or hexadecimal floats (float.hex(), which R reads back
                    exactly) for real-valued tables.

  tables.csv        One row per table: name and pc, the pseudocount vlr()
                    uses by default (half the smallest non-zero value of the
                    whole table when it has zeros, else 0; hex float).

  values.csv        Long format (table, field, k, value), `value` a hex
                    float and `k` the 1-based position in the condensed
                    (packed lower triangle, column by column) vector over
                    the features. Fields:
                      pairwise  scikit-bio's pairwise_vlr(x + pc) (features
                                are its columns), .condensed_form()
                      vlr       scikit-bio's vlr(x_i + pc, x_j + pc) of each
                                pair, in the same order

Prerequisites: Python 3 with numpy and scikit-bio 0.7.4 (the version the
fixtures were made with). Run from the repository root:

    python data-raw/expected/proportionality.py
    Rscript data-raw/expected/proportionality.r
"""

import csv
import os

import numpy as np
import skbio
from skbio.stats.composition import pairwise_vlr, vlr

assert skbio.__version__ == "0.7.4", skbio.__version__

OUT = os.path.join("data-raw", "expected", "proportionality")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261121)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

base = rng.lognormal(3, 1, size=(20, 1))
nearprop = base * rng.uniform(1, 100, size=(1, 8)) * np.exp(
    rng.normal(0, 1e-4, size=(20, 8)))
nearprop[:, 6:] = rng.lognormal(3, 1, size=(20, 2))    # two unrelated

tables = {
    "skbio":    np.array([[1, 1, 2], [3, 5, 8], [13, 21, 55]]),
    "small":    rng.poisson(5, size=(12, 6)) + 1,
    "zeros":    rng.poisson(1.5, size=(30, 20)),
    "wide":     rng.poisson(8, size=(10, 60)) + 1,
    "prop":     rng.dirichlet(np.ones(25), size=15),
    "nearprop": nearprop,
}

for name, tab in tables.items():
    with open(os.path.join(OUT, f"table-{name}.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        if tab.dtype.kind == "f":
            w.writerows([[float(v).hex() for v in row] for row in tab])
        else:
            w.writerows(tab.tolist())


def default_pc(tab):
    """vlr()'s pseudocount = NULL: half the smallest non-zero value."""
    return float(tab[tab > 0].min()) / 2 if (tab == 0).any() else 0.0


rows = []
values = []

for name, tab in tables.items():
    pc = default_pc(tab)
    x = tab.astype(float) + pc
    D = x.shape[1]
    rows.append([name, float(pc).hex()])

    pw = pairwise_vlr(x).condensed_form()
    one = [vlr(x[:, i], x[:, j]) for i in range(D) for j in range(i + 1, D)]
    assert len(pw) == len(one) == D * (D - 1) // 2

    for field, vec in (("pairwise", pw), ("vlr", one)):
        for k, v in enumerate(vec):
            values.append([name, field, k + 1, float(v).hex()])


with open(os.path.join(OUT, "tables.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["table", "pc"])
    w.writerows(rows)

with open(os.path.join(OUT, "values.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["table", "field", "k", "value"])
    w.writerows(values)
