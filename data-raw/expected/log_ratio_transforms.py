"""Reference fixtures for ecodive's alr() and ilr(), from scikit-bio.

Writes, under data-raw/expected/log_ratio_transforms/ (then
data-raw/expected/log_ratio_transforms.r converts them to
tests/tinytests/expected/log_ratio_transforms.rds):

  counts-<name>.csv   Integer count tables (samples x features, with zeros),
                      one value per cell. The pseudocount is added in both
                      languages, so R and Python start from identical doubles.

  values.csv          Long format (case, field, row, col, value). `value` is
                      a hexadecimal float (float.hex()), which R reads back
                      exactly. Fields:
                        counts   the case's table (its name)
                        pc       the pseudocount added to every count
                        basis    an input basis or SBP (ilr cases)
                        result   scikit-bio's output, samples x coordinates

  cases.csv           One row per case: id, table, function, ref_idx (alr,
                      0-based), basis kind (default, sbp, orthonormal) and
                      the pseudocount.

  helmert-<D>.csv     _gram_schmidt_basis(D), hex floats, (D - 1) x D.
  sbp-<name>.csv      An SBP and sbp_basis() of it (hex floats).

Prerequisites: Python 3 with numpy and scikit-bio 0.7.4 (the version the
fixtures were made with). Run from the repository root:

    python data-raw/expected/log_ratio_transforms.py
    Rscript data-raw/expected/log_ratio_transforms.r
"""

import csv
import os

import numpy as np
import skbio
from skbio.stats.composition import alr, ilr, sbp_basis
from skbio.stats.composition._base import _gram_schmidt_basis

assert skbio.__version__ == "0.7.4", skbio.__version__

OUT = os.path.join("data-raw", "expected", "log_ratio_transforms")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261119)


def hexmat(a):
    return [[float(v).hex() for v in row] for row in np.atleast_2d(a)]


def write_hex(path, a):
    with open(path, "w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(hexmat(a))


def random_sbp(D):
    """A random sequential binary partition of D parts, (D - 1) x D."""
    rows = []

    def split(parts):
        if len(parts) < 2:
            return
        parts = list(rng.permutation(parts))
        k = int(rng.integers(1, len(parts)))
        row = np.zeros(D, dtype=int)
        row[parts[:k]] = 1
        row[parts[k:]] = -1
        rows.append(row)
        split(parts[:k])
        split(parts[k:])

    split(list(range(D)))
    return np.array(rows)


def random_basis(D):
    """A random orthonormal basis of clr space, (D - 1) x D."""
    A = rng.normal(size=(D, D - 1))
    A -= A.mean(axis=0)               # columns in clr space
    Q, _ = np.linalg.qr(A)            # orthonormal columns, still sum to 0
    return Q.T


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

tables = {
    "doc":    (np.array([[1, 3, 4, 2]]), 0.0),          # skbio's doc example (scaled)
    "small":  (rng.poisson(3, size=(12, 6)), 0.5),
    "medium": (rng.poisson(2, size=(30, 15)) * rng.integers(1, 50, size=15), 1.0),
    "wide":   (rng.poisson(5, size=(8, 60)), 0.5),
}

for name, (tab, _) in tables.items():
    with open(os.path.join(OUT, f"counts-{name}.csv"), "w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(tab.tolist())


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

cases = []
values = []


def add(table, fn, mat, ref_idx="", kind="", basis=None):
    cid = len(cases) + 1
    pc = tables[table][1]
    cases.append([cid, table, fn, ref_idx, kind, float(pc).hex()])
    if basis is not None:
        for i, row in enumerate(hexmat(basis)):
            for j, v in enumerate(row):
                values.append([cid, "basis", i + 1, j + 1, v])
    for i, row in enumerate(hexmat(mat)):
        for j, v in enumerate(row):
            values.append([cid, "result", i + 1, j + 1, v])


for name, (tab, pc) in tables.items():
    x = tab.astype(float) + pc
    D = x.shape[1]

    for ref in sorted({0, D // 2, D - 1}):
        add(name, "alr", alr(x, ref_idx=ref), ref_idx=ref)

    add(name, "ilr", ilr(x), kind="default")

    sbp = random_sbp(D)
    add(name, "ilr", ilr(x, sbp_basis(sbp)), kind="sbp", basis=sbp)

    if D <= 20:
        B = random_basis(D)
        add(name, "ilr", ilr(x, B), kind="orthonormal", basis=B)


with open(os.path.join(OUT, "cases.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "table", "fn", "ref_idx", "kind", "pc"])
    w.writerows(cases)

with open(os.path.join(OUT, "values.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "field", "row", "col", "value"])
    w.writerows(values)


# ---------------------------------------------------------------------------
# Bases
# ---------------------------------------------------------------------------

for D in (2, 3, 4, 7, 20, 60):
    write_hex(os.path.join(OUT, f"helmert-{D}.csv"), _gram_schmidt_basis(D))

sbps = {
    "doc": np.array([[1, 1, -1, -1, -1],
                     [1, -1, 0, 0, 0],
                     [0, 0, 1, -1, -1],
                     [0, 0, 0, 1, -1]]),
    "random20": random_sbp(20),
}
for name, sbp in sbps.items():
    with open(os.path.join(OUT, f"sbp-{name}.csv"), "w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(sbp.tolist())
    write_hex(os.path.join(OUT, f"sbp-{name}-basis.csv"), sbp_basis(sbp))

print("scikit-bio", skbio.__version__, "numpy", np.__version__, ":", len(cases), "cases")
