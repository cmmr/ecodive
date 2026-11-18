# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for procrustes() from SciPy.

Reads the configurations that data-raw/expected/procrustes.r writes to
data-raw/expected/procrustes/ (<pair>_ord.csv, rotated onto
<pair>_target.csv) and writes scipy.spatial.procrustes()'s disparity, the
symmetric Procrustes sum of squares, for each pair to scipy_disparity.csv.
data-raw/expected/procrustes.r folds them into
tests/tinytests/expected/procrustes.rds. Generated with SciPy 1.18.1
(BSD-3-Clause):

    Rscript data-raw/expected/procrustes.r
    python data-raw/expected/procrustes.py
    Rscript data-raw/expected/procrustes.r

Every number is written with 17 significant digits.
"""

import os
import sys

import numpy as np
from scipy.spatial import procrustes


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "procrustes")


def main():
    rows = []
    for pair in ("YX", "RX", "NX"):
        ord_ = np.loadtxt(os.path.join(DIR, pair + "_ord.csv"), delimiter=",", skiprows=1, dtype=str)
        tgt = np.loadtxt(os.path.join(DIR, pair + "_target.csv"), delimiter=",", skiprows=1, dtype=str)
        ord_ = np.char.strip(ord_, '" ').astype(float)
        tgt = np.char.strip(tgt, '" ').astype(float)
        _, _, disparity = procrustes(tgt, ord_)
        rows.append("%s,%.17g" % (pair, disparity))
    with open(os.path.join(DIR, "scipy_disparity.csv"), "w", newline="\n") as f:
        f.write("pair,disparity\n" + "\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
