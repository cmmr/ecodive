# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for ca() from scikit-bio.

Reads data-raw/expected/ca/dune.csv and varespec.csv (vegan's tables,
written by data-raw/expected/ca.r) and writes scikit-bio's correspondence
analysis of each, in scaling 2 (ca()'s), to CSV files in the same
directory. data-raw/expected/ca.r folds them into
tests/tinytests/expected/ca.rds, so that the tests need no Python.
Generated with scikit-bio 0.7.4 (BSD-3-Clause):

    Rscript data-raw/expected/ca.r      # writes the inputs
    python data-raw/expected/ca.py
    Rscript data-raw/expected/ca.r      # folds in the outputs

Every number is written with 17 significant digits.
"""

import os
import sys

import pandas as pd

from skbio.stats.ordination import ca


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "ca")


def write(df, name):
    df.to_csv(os.path.join(DIR, name), float_format="%.17g", lineterminator="\n")


def main():
    for name in ("dune", "varespec"):
        y = pd.read_csv(os.path.join(DIR, name + ".csv"), index_col=0)
        y.index = y.index.astype(str)
        res = ca(y.astype(float), scaling=2)
        write(res.samples, "skbio_%s_samples.csv" % name)
        write(res.features, "skbio_%s_features.csv" % name)
        write(res.eigvals.to_frame("eigvals"), "skbio_%s_eigvals.csv" % name)


if __name__ == "__main__":
    main()
