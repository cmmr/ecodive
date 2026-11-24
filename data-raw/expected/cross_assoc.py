# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for cross_assoc() from SciPy.

Reads the tables that data-raw/expected/cross_assoc.r writes to
data-raw/expected/cross_assoc/ (x.csv: the relative abundances that
cross_assoc() correlates, samples x taxa; y.csv: the environmental
variables, "NA" for missing values), both as hexadecimal floats
(float.hex(), which R reads back exactly), and writes scipy.csv: for every
taxon i and variable j (1-based), on the samples where y[, j] is not
missing, scipy.stats.spearmanr(), pearsonr() and
kendalltau(method='asymptotic'): each coefficient and its two-sided
p-value, as hexadecimal floats ("nan" where SciPy gives NaN).
data-raw/expected/cross_assoc.r then folds them into
tests/tinytests/expected/cross_assoc.rds.

Prerequisites: Python 3 with numpy and SciPy 1.18.1 (BSD-3-Clause), the
version the fixtures were made with. Run from the repository root:

    Rscript data-raw/expected/cross_assoc.r
    python data-raw/expected/cross_assoc.py
    Rscript data-raw/expected/cross_assoc.r
"""

import csv
import os
import sys
import warnings

import numpy as np
import scipy
from scipy import stats

assert scipy.__version__ == "1.18.1", scipy.__version__

DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "cross_assoc")


def read(name):
    with open(os.path.join(DIR, name), newline="") as f:
        rows = list(csv.reader(f))
    header, rows = rows[0], rows[1:]
    vals = [[np.nan if v == "NA" else float.fromhex(v) for v in row] for row in rows]
    return header, np.array(vals, dtype=float)


def hexf(v):
    v = float(v)
    return "nan" if np.isnan(v) else v.hex()


def main():
    _, x = read("x.csv")
    _, y = read("y.csv")
    out = [["i", "j", "n", "spearman", "spearman_p", "pearson", "pearson_p",
            "kendall", "kendall_p"]]
    for i in range(x.shape[1]):
        for j in range(y.shape[1]):
            ok = ~np.isnan(y[:, j])
            a, b = x[ok, i], y[ok, j]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = stats.spearmanr(a, b)
                p = stats.pearsonr(a, b)
                k = stats.kendalltau(a, b, method="asymptotic")
            out.append([i + 1, j + 1, int(ok.sum()),
                        hexf(s.statistic), hexf(s.pvalue),
                        hexf(p.statistic), hexf(p.pvalue),
                        hexf(k.statistic), hexf(k.pvalue)])
    with open(os.path.join(DIR, "scipy.csv"), "w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(out)


if __name__ == "__main__":
    main()
