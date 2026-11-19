# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for lda() from scikit-learn's LinearDiscriminantAnalysis.

Reads the cases that data-raw/expected/lda.r writes to
data-raw/expected/lda/ (cases.csv, and each case's normalized table
<name>_x.csv, groups <name>_groups.csv and priors <name>_prior.csv, in
sorted group order, as R's factor levels), fits
LinearDiscriminantAnalysis(solver = "svd", priors = prior), and writes the
transform() scores, scalings_ (the columns the scores use),
explained_variance_ratio_, means_ and xbar_ to <name>_sk_<attr>.csv in the
same directory, with 17 significant digits. data-raw/expected/lda.r folds
them into tests/tinytests/expected/lda.rds, so that the tests need no
Python.

Prerequisites: Python 3 with numpy, pandas and scikit-learn 1.9.1 (the
version the C code in src/lda.c was ported from):

    Rscript data-raw/expected/lda.r
    python data-raw/expected/lda.py
    Rscript data-raw/expected/lda.r
"""

import os
import sys

import numpy as np
import pandas as pd
import sklearn
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "lda")


def save(name, attr, value):
    path = os.path.join(DIR, f"{name}_sk_{attr}.csv")
    np.savetxt(path, np.atleast_2d(value) if np.ndim(value) < 2 else value,
               delimiter=",", fmt="%.17g")


def main():
    assert sklearn.__version__ == "1.9.1", sklearn.__version__

    cases = pd.read_csv(os.path.join(DIR, "cases.csv"))

    for name in cases["name"]:
        x = pd.read_csv(os.path.join(DIR, f"{name}_x.csv"), index_col=0).to_numpy(dtype=float)
        g = pd.read_csv(os.path.join(DIR, f"{name}_groups.csv"), index_col=0)["group"].astype(str)
        prior = pd.read_csv(os.path.join(DIR, f"{name}_prior.csv"))["prior"].to_numpy(dtype=float)
        assert list(np.unique(g)) == sorted(g.unique())

        m = LinearDiscriminantAnalysis(solver="svd", priors=prior).fit(x, g.to_numpy())
        scores = m.transform(x)
        save(name, "scores", scores)
        save(name, "scalings", m.scalings_[:, : scores.shape[1]])
        save(name, "evr", m.explained_variance_ratio_)
        save(name, "means", m.means_)
        save(name, "xbar", m.xbar_)


if __name__ == "__main__":
    main()
