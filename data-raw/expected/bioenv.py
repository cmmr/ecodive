# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for bioenv() from scikit-bio.

Reads data-raw/expected/bioenv/dm.csv and env.csv (varespec's Bray-Curtis
dissimilarities and varechem, written by data-raw/expected/bioenv.r) and
writes scikit-bio's bioenv() of them, on its default options (Spearman,
Euclidean distances on z-scored variables), to skbio.csv in the same
directory. data-raw/expected/bioenv.r folds it into
tests/tinytests/expected/bioenv.rds, so that the tests need no Python.
Generated with scikit-bio 0.7.4 (BSD-3-Clause):

    Rscript data-raw/expected/bioenv.r      # writes the inputs
    python data-raw/expected/bioenv.py
    Rscript data-raw/expected/bioenv.r      # folds in the outputs

Every number is written with 17 significant digits.
"""

import os
import sys

import pandas as pd

from skbio import DistanceMatrix
from skbio.stats.distance import bioenv


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "bioenv")


def main():
    dm = pd.read_csv(os.path.join(DIR, "dm.csv"), index_col=0)
    dm.index = dm.index.astype(str)
    dm = DistanceMatrix(dm.values, ids=list(dm.index))
    env = pd.read_csv(os.path.join(DIR, "env.csv"), index_col=0)
    env.index = env.index.astype(str)
    res = bioenv(dm, env)
    res.to_csv(
        os.path.join(DIR, "skbio.csv"), float_format="%.17g", lineterminator="\n")


if __name__ == "__main__":
    main()
