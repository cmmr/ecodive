# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for rda() from scikit-bio.

Reads data-raw/expected/rda/dune.csv and dune_env.csv (written by
data-raw/expected/rda.r: vegan's dune table and the dummy-coded
`~ Management + A1` design) and writes scikit-bio's RDA of the
Hellinger-transformed table, in scaling 2 (rda()'s), to CSV files in the same
directory. data-raw/expected/rda.r folds them into
tests/tinytests/expected/rda.rds, so that the tests need no Python.
Generated with scikit-bio 0.7.4 (BSD-3-Clause):

    Rscript data-raw/expected/rda.r      # writes the inputs
    python data-raw/expected/rda.py
    Rscript data-raw/expected/rda.r      # folds in the outputs

scikit-bio's "eigvals" are singular values s, not variances, and its
proportion_explained is s / sum(s), so neither is a share of variance;
test-rda.r checks rda()'s eig against s^2 / (n - 1) instead. Every number
is written with 17 significant digits.
"""

import os
import sys

import numpy as np
import pandas as pd

from skbio.stats.ordination import rda


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "rda")


def write(df, name):
    df.to_csv(os.path.join(DIR, name), float_format="%.17g", lineterminator="\n")


def main():
    y = pd.read_csv(os.path.join(DIR, "dune.csv"), index_col=0)
    x = pd.read_csv(os.path.join(DIR, "dune_env.csv"), index_col=0)
    y.index = y.index.astype(str)
    x.index = x.index.astype(str)

    v = y.to_numpy(dtype=float)
    h = np.sqrt(v / v.sum(axis=1, keepdims=True))
    h = pd.DataFrame(h, index=y.index, columns=y.columns)

    res = rda(h, x, scale_Y=False, scaling=2)
    write(res.samples, "skbio_samples.csv")
    write(res.features, "skbio_features.csv")
    write(res.biplot_scores, "skbio_biplot.csv")
    write(res.eigvals.to_frame("eigvals"), "skbio_eigvals.csv")


if __name__ == "__main__":
    main()
