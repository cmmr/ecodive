# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for cca() from scikit-bio.

Reads data-raw/expected/cca/dune.csv and dune_env.csv (written by
data-raw/expected/cca.r: vegan's dune table and the dummy-coded
`~ Management + A1` design) and writes scikit-bio's CCA, in scaling 2
(cca()'s), to CSV files in the same directory. data-raw/expected/cca.r folds
them into tests/tinytests/expected/cca.rds, so that the tests need no Python.
Generated with scikit-bio 0.7.4 (BSD-3-Clause):

    Rscript data-raw/expected/cca.r      # writes the inputs
    python data-raw/expected/cca.py
    Rscript data-raw/expected/cca.r      # folds in the outputs

Every number is written with 17 significant digits.
"""

import os
import sys

import pandas as pd

from skbio.stats.ordination import cca


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "cca")


def write(df, name):
    df.to_csv(os.path.join(DIR, name), float_format="%.17g", lineterminator="\n")


def main():
    y = pd.read_csv(os.path.join(DIR, "dune.csv"), index_col=0)
    x = pd.read_csv(os.path.join(DIR, "dune_env.csv"), index_col=0)
    y.index = y.index.astype(str)
    x.index = x.index.astype(str)

    res = cca(y.astype(float), x, scaling=2)
    write(res.samples, "skbio_samples.csv")
    write(res.features, "skbio_features.csv")
    write(res.biplot_scores, "skbio_biplot.csv")
    write(res.eigvals.to_frame("eigvals"), "skbio_eigvals.csv")


if __name__ == "__main__":
    main()
