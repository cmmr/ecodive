"""Reference fixtures for ecodive's replace_zeros_mult(), from scikit-bio.

Writes, under data-raw/expected/replace_zeros/ (then
data-raw/expected/replace_zeros.r converts them to
tests/tinytests/expected/replace_zeros.rds):

  cases.csv     One row per case: id, rows, cols, and delta ('' for the
                default, 'per_sample' for one value per sample).

  <id>-x.csv    The input table (samples x features), hexadecimal floats
                (float.hex(), read back exactly by as.numeric() in R).
  <id>-d.csv    The per-sample deltas, if any (hex floats, one per sample).
  <id>-r.csv    multi_replace() of it (hex floats, samples x features).

Cases cover each branch of numpy's pairwise summation (fewer than 8
features, up to 128, and above), integer counts and real-valued
proportions, a feature that is zero everywhere, a given delta, and one delta
per sample (an (n, 1) array, which multi_replace() broadcasts).

Prerequisites: Python 3 with numpy and scikit-bio 0.7.4 (the version the
fixtures were made with). Run from the repository root:

    python data-raw/expected/replace_zeros.py
    Rscript data-raw/expected/replace_zeros.r
"""

import csv
import os

import numpy as np
import skbio
from skbio.stats.composition import multi_replace

assert skbio.__version__ == "0.7.4", skbio.__version__

OUT = os.path.join("data-raw", "expected", "replace_zeros")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261119)


def write_hex(path, a):
    with open(path, "w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(
            [[float(v).hex() for v in row] for row in np.atleast_2d(a)]
        )


def counts(n, D, p_zero, scale=50):
    """Integer counts with about p_zero zeros; no sample is all zeros."""
    x = rng.poisson(rng.gamma(0.5, scale, size=(n, D))).astype(float)
    x[rng.random((n, D)) < p_zero] = 0
    for i in range(n):
        if not x[i].any():
            x[i, rng.integers(D)] = 1
    return x


cases = []


def case(id, x, delta=None):
    x = np.asarray(x, dtype=float)
    if delta is None:
        r = multi_replace(x)
        kind = ""
    elif np.ndim(delta) == 0:
        r = multi_replace(x, delta=float(delta))
        kind = float(delta).hex()
    else:
        d = np.asarray(delta, dtype=float).reshape(-1, 1)
        r = multi_replace(x, delta=d)
        kind = "per_sample"
        write_hex(os.path.join(OUT, f"{id}-d.csv"), d.T)
    r = np.atleast_2d(r)
    assert r.shape == x.shape
    write_hex(os.path.join(OUT, f"{id}-x.csv"), x)
    write_hex(os.path.join(OUT, f"{id}-r.csv"), r)
    cases.append([id, x.shape[0], x.shape[1], kind])


# scikit-bio's own examples (its docstring and tests).
case("doc", [[0.2, 0.4, 0.4, 0], [0, 0.5, 0.5, 0]])
case("plan", [[1, 2, 3, 0, 5], [1, 0, 0, 4, 5], [1, 2, 3, 4, 5]])
case("plan-delta", [[1, 2, 3, 0, 5], [1, 0, 0, 4, 5], [1, 2, 3, 4, 5]], 0.01)

# Every branch of the pairwise sum: D < 8, 8 <= D <= 128, D > 128 (with
# remainders).
for D in (2, 3, 7, 8, 9, 16, 23, 128, 129, 300, 1037):
    case(f"counts-{D}", counts(12 if D < 300 else 3, D, 0.4))

# Real-valued input (proportions, and arbitrary positive reals).
p = counts(10, 40, 0.5)
case("props", p / p.sum(axis=1, keepdims=True))
case("reals", counts(10, 40, 0.3) * rng.uniform(1e-3, 1e3, size=(10, 40)))

# A feature that is zero in every sample.
x = counts(8, 12, 0.3)
x[:, 5] = 0
case("zero-feature", x)

# A given delta, and one per sample.
x = counts(15, 30, 0.5)
case("delta", x, 1e-4)
case("per-sample", x, 0.5 / x.sum(axis=1))

# Samples without zeros, and one sample only.
case("no-zeros", counts(5, 10, 0.0) + 1)
case("one-sample", counts(1, 20, 0.3))

with open(os.path.join(OUT, "cases.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "rows", "cols", "delta"])
    w.writerows(cases)
