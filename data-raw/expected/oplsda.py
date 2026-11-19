# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for oplsda() from pyopls's OPLS.

Reads the cases that data-raw/expected/oplsda.r writes to
data-raw/expected/oplsda/ (cases.csv, and each case's normalized table
<name>_x.csv and 0/1 group indicator <name>_groups.csv), and fits

  * pyopls's OPLS(n_components = k, scale = scale) to the table and the
    indicator, writing T_ortho_, W_ortho_ and P_ortho_; then
  * scikit-learn's PLSRegression(1, scale = False) to pyopls's filtered
    table (OPLS.transform()) and the indicator as pyopls centres and scales
    it, writing its x_scores_, x_weights_, x_loadings_, y_loadings_, and the
    R^2 of its fitted values (pyopls's OPLSValidator.r_squared_Y_).

to <name>_py_<attr>.csv in the same directory, with 17 significant digits.
(pyopls's README fits PLSRegression(1) with its default scale = True to the
filtered table, which rescales the filtered features; scale = False keeps
Trygg & Wold's predictive component, as ropls computes it.)
data-raw/expected/oplsda.r folds them into tests/tinytests/expected/oplsda.rds,
so that the tests need no Python.

Prerequisites: Python 3 with numpy, pandas, scikit-learn 1.9.1 and pyopls
20.3.post1 (the version the C code in src/oplsda.c was ported from):

    Rscript data-raw/expected/oplsda.r
    python data-raw/expected/oplsda.py
    Rscript data-raw/expected/oplsda.r
"""

import os
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
import sklearn
from pyopls import OPLS
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "oplsda")


def write(name, attr, value):
    path = os.path.join(DIR, f"{name}_py_{attr}.csv")
    value = np.asarray(value, dtype=float)
    if value.ndim < 2:
        value = value.reshape(-1, 1)
    np.savetxt(path, value, delimiter=",", fmt="%.17g")


def main():
    assert sklearn.__version__ == "1.9.1", sklearn.__version__
    assert version("pyopls") == "20.3.post1", version("pyopls")

    cases = pd.read_csv(os.path.join(DIR, "cases.csv"))

    for _, case in cases.iterrows():
        name = case["name"]
        x = pd.read_csv(os.path.join(DIR, f"{name}_x.csv"), index_col=0).to_numpy(dtype=float)
        y = pd.read_csv(os.path.join(DIR, f"{name}_groups.csv"), index_col=0)["group"].to_numpy(dtype=float)
        k = int(case["k"])
        scale = bool(case["scale"])

        opls = OPLS(n_components=k, scale=scale).fit(x, y)
        write(name, "t_ortho", opls.T_ortho_)
        write(name, "w_ortho", opls.W_ortho_)
        write(name, "p_ortho", opls.P_ortho_)

        z = opls.transform(x)
        yc = (y - opls.y_mean_) / opls.y_std_
        pls = PLSRegression(n_components=1, scale=False).fit(z, yc)
        write(name, "t_pred", pls.x_scores_)
        write(name, "w_pred", pls.x_weights_)
        write(name, "p_pred", pls.x_loadings_)
        write(name, "q_pred", pls.y_loadings_)
        write(name, "r2y", [r2_score(yc, pls.predict(z).ravel())])


if __name__ == "__main__":
    main()
