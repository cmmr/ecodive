# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for plsda() from scikit-learn's PLSRegression.

Reads the cases that data-raw/expected/plsda.r writes to
data-raw/expected/plsda/ (cases.csv, and each case's normalized table
<name>_x.csv and groups <name>_groups.csv), fits
PLSRegression(n_components = k, scale = scale) to the table and the 0/1
dummy matrix of the groups (one column per group, in sorted order, as R's
factor levels), and writes x_scores_, x_weights_, x_loadings_ and
y_loadings_ to <name>_sk_<attr>.csv in the same directory, with 17
significant digits. data-raw/expected/plsda.r folds them into
tests/tinytests/expected/plsda.rds, so that the tests need no Python.

Each case is fitted with tol = 1e-28 and a large max_iter, so that NIPALS
converges fully (plsda() computes the converged weights directly). The
two-group case is also fitted with the default tol and max_iter
(<name>_sk_default_<attr>.csv): NIPALS then converges at once.

Prerequisites: Python 3 with numpy, pandas and scikit-learn 1.9.1 (the
version the C code in src/plsda.c was ported from):

    Rscript data-raw/expected/plsda.r
    python data-raw/expected/plsda.py
    Rscript data-raw/expected/plsda.r
"""

import os
import sys

import numpy as np
import pandas as pd
import sklearn
from sklearn.cross_decomposition import PLSRegression


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "plsda")


def write(name, suffix, model):
    for attr in ["x_scores_", "x_weights_", "x_loadings_", "y_loadings_"]:
        path = os.path.join(DIR, f"{name}_sk{suffix}_{attr.rstrip('_')}.csv")
        np.savetxt(path, getattr(model, attr), delimiter=",", fmt="%.17g")


def main():
    assert sklearn.__version__ == "1.9.1", sklearn.__version__

    cases = pd.read_csv(os.path.join(DIR, "cases.csv"))

    for _, case in cases.iterrows():
        name = case["name"]
        x = pd.read_csv(os.path.join(DIR, f"{name}_x.csv"), index_col=0)
        g = pd.read_csv(os.path.join(DIR, f"{name}_groups.csv"), index_col=0)["group"].astype(str)
        levels = sorted(g.unique())
        y = np.column_stack([(g == lev).to_numpy(dtype=float) for lev in levels])
        x = x.to_numpy(dtype=float)
        k = int(case["k"])
        scale = bool(case["scale"])

        m = PLSRegression(n_components=k, scale=scale, tol=1e-28, max_iter=100000).fit(x, y)
        assert max(m.n_iter_) < 100000, m.n_iter_
        write(name, "", m)

        if len(levels) == 2:
            m = PLSRegression(n_components=k, scale=scale).fit(x, y)
            write(name, "_default", m)


if __name__ == "__main__":
    main()
