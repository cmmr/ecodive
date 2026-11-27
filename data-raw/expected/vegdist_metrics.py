"""SciPy's Mahalanobis distances, for ecodive's mahalanobis_distance().

Reads the tables that `Rscript data-raw/expected/vegdist_metrics.r inputs`
writes under data-raw/expected/vegdist_metrics/:

  maha-<name>.csv   Samples x features, no header.

and writes, for each:

  scipy-<name>.csv  The Mahalanobis distance of every pair of samples, in
                    the order of an R dist's lower triangle ((1, 2), (1, 3),
                    ..., (2, 3), ...), one per line, as repr() of a float
                    (which R reads back exactly).

Each distance is scipy.spatial.distance.mahalanobis(u, v, VI) with VI the
inverse of the features' sample covariance (numpy.cov, n - 1).

Prerequisites: Python 3 with numpy and SciPy. Run from the repository root:

    Rscript data-raw/expected/vegdist_metrics.r inputs
    python data-raw/expected/vegdist_metrics.py
"""

import glob
import os

import numpy as np
from scipy.spatial.distance import mahalanobis

DIR = os.path.join("data-raw", "expected", "vegdist_metrics")

for path in sorted(glob.glob(os.path.join(DIR, "maha-*.csv"))):
    name = os.path.basename(path)[len("maha-"):-len(".csv")]
    x = np.loadtxt(path, delimiter=",", ndmin=2)
    vi = np.linalg.inv(np.cov(x, rowvar=False))
    n = x.shape[0]
    out = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            out.append(repr(float(mahalanobis(x[i], x[j], vi))))
    with open(os.path.join(DIR, "scipy-" + name + ".csv"), "w") as f:
        f.write("\n".join(out) + "\n")
    print(name, x.shape, len(out))
