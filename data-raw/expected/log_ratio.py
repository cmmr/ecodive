"""Reference fixtures for ecodive's log_ratio() and balance().

Writes, under data-raw/expected/log_ratio/ (then data-raw/expected/log_ratio.r
converts them to tests/tinytests/expected/log_ratio.rds):

  table-<name>.csv  The input tables (samples x features): integer counts
                    with zeros, or hexadecimal floats (float.hex(), which R
                    reads back exactly) for the proportions table.

  cases.csv         One row per case: id, table, num and den (the 1-based
                    features of each side, separated by spaces), and pc,
                    the pseudocount balance() uses by default (half the
                    smallest non-zero value of the whole table when it has
                    zeros, else 0; hex float).

  values.csv        Long format (id, field, sample, value): `value` is a hex
                    float, or NA. Fields:
                      skbio     scikit-bio's ilr(x + pc, sbp_basis(sbp)),
                                first coordinate, where the SBP's first row
                                is the balance (+1 numerator, -1
                                denominator) and the other rows complete it
                      closed    sqrt(rs / (r + s)) * (mean(log x_N) -
                                mean(log x_D)), on x + pc, in numpy
                      qurro     Qurro's qarcoal() log ratio of the raw table
                                (no pseudocount): log(Num_Sum / Denom_Sum);
                                NA for the samples it drops (a zero sum)
                      lr_pc1    np.log(sum(x_N + 1) / sum(x_D + 1)), the
                                log ratio with pseudocount 1

Prerequisites: Python 3 with numpy, pandas, biom-format, scikit-bio 0.7.4 and
qurro 0.9.0 (the versions the fixtures were made with). qarcoal imports
qiime2 only for a type annotation (Metadata, used when `samples_to_use` is
given), so a stub module stands in for it. Run from the repository root:

    python data-raw/expected/log_ratio.py
    Rscript data-raw/expected/log_ratio.r
"""

import csv
import os
import sys
import types

import numpy as np
import pandas as pd
import biom
import skbio
from skbio.stats.composition import ilr, sbp_basis

# qarcoal's `from qiime2 import Metadata`: a stub, as qiime2 is not needed.
_q2 = types.ModuleType("qiime2")
_q2.Metadata = object
sys.modules.setdefault("qiime2", _q2)

import qurro
from qurro.qarcoal import qarcoal

assert skbio.__version__ == "0.7.4", skbio.__version__
assert qurro.__version__ == "0.9.0", qurro.__version__

OUT = os.path.join("data-raw", "expected", "log_ratio")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261120)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

tables = {
    "small":  rng.poisson(3, size=(12, 6)),
    "medium": rng.poisson(2, size=(30, 15)) * rng.integers(1, 50, size=15),
    "sparse": rng.poisson(0.6, size=(25, 40)),
    "wide":   rng.poisson(5, size=(8, 60)),
    "prop":   rng.dirichlet(np.ones(25), size=10),
}

for name, tab in tables.items():
    with open(os.path.join(OUT, f"table-{name}.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        if tab.dtype.kind == "f":
            w.writerows([[float(v).hex() for v in row] for row in tab])
        else:
            w.writerows(tab.tolist())


def default_pc(tab):
    """balance()'s pseudocount = NULL: half the smallest non-zero value."""
    return float(tab[tab > 0].min()) / 2 if (tab == 0).any() else 0.0


def splits(parts, D):
    """SBP rows splitting `parts` down to single parts."""
    rows = []
    if len(parts) < 2:
        return rows
    parts = list(rng.permutation(parts))
    k = int(rng.integers(1, len(parts)))
    row = np.zeros(D, dtype=int)
    row[parts[:k]] = 1
    row[parts[k:]] = -1
    rows.append(row)
    return rows + splits(parts[:k], D) + splits(parts[k:], D)


def full_sbp(num, den, D):
    """A (D - 1) x D SBP whose first row is num (+1) vs den (-1)."""
    first = np.zeros(D, dtype=int)
    first[num] = 1
    first[den] = -1
    rows = [first] + splits(num, D) + splits(den, D)
    rest = [j for j in range(D) if j not in set(num) | set(den)]
    if rest:
        top = np.zeros(D, dtype=int)
        top[num + den] = 1
        top[rest] = -1
        rows += [top] + splits(rest, D)
    sbp = np.array(rows)
    assert sbp.shape == (D - 1, D)
    return sbp


def qurro_log_ratio(tab, num, den):
    """Qurro's qarcoal() on the raw table; NaN where it drops a sample."""
    n, D = tab.shape
    samples = [f"S{i + 1}" for i in range(n)]
    features = [f"F{j + 1}" for j in range(D)]
    table = biom.Table(tab.T.astype(float), features, samples)
    taxon = ["other"] * D
    for j in num:
        taxon[j] = "num"
    for j in den:
        taxon[j] = "den"
    taxonomy = pd.DataFrame({"Taxon": taxon}, index=features)
    res = qarcoal(table, taxonomy, "num", "den")["log_ratio"]
    return np.array([res.get(s, np.nan) for s in samples])


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

cases = []
values = []


def hexv(v):
    return "NA" if np.isnan(v) else float(v).hex()


def add(table, num, den):
    tab = tables[table]
    n, D = tab.shape
    cid = len(cases) + 1
    pc = default_pc(tab)
    x = tab.astype(float) + pc
    r, s = len(num), len(den)

    sbp = full_sbp(num, den, D)
    sk = ilr(x, sbp_basis(sbp))[:, 0]

    lx = np.log(x)
    closed = np.sqrt(r * s / (r + s)) * (lx[:, num].mean(axis=1) - lx[:, den].mean(axis=1))

    if tab.dtype.kind == "f":
        qu = np.full(n, np.nan)   # Qurro reads count tables; skip
    else:
        qu = qurro_log_ratio(tab, num, den)

    t1 = tab.astype(float) + 1
    lr1 = np.log(t1[:, num].sum(axis=1) / t1[:, den].sum(axis=1))

    cases.append([
        cid, table, " ".join(str(j + 1) for j in num),
        " ".join(str(j + 1) for j in den), float(pc).hex()])
    for field, vec in (("skbio", sk), ("closed", closed), ("qurro", qu), ("lr_pc1", lr1)):
        for i, v in enumerate(vec):
            values.append([cid, field, i + 1, hexv(v)])


for name, tab in tables.items():
    D = tab.shape[1]
    perm = [int(j) for j in rng.permutation(D)]
    add(name, perm[:1], perm[1:2])                         # two features
    add(name, perm[:2], perm[2:5])                         # small groups
    k = int(rng.integers(1, D))
    add(name, sorted(perm[:k]), sorted(perm[k:]))          # every feature
    m = int(rng.integers(2, D))
    k = int(rng.integers(1, m))
    add(name, perm[:k], perm[k:m])                         # random sizes


with open(os.path.join(OUT, "cases.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "table", "num", "den", "pc"])
    w.writerows(cases)

with open(os.path.join(OUT, "values.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "field", "sample", "value"])
    w.writerows(values)

print("scikit-bio", skbio.__version__, "qurro", qurro.__version__,
      "numpy", np.__version__, ":", len(cases), "cases")
