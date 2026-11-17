# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for pcoa() and pcoa_biplot() from scikit-bio.

Writes small CSV files to data-raw/expected/pcoa_biplot/ (or the directory
given as the first argument), which data-raw/expected/pcoa_biplot.r turns
into tests/tinytests/expected/pcoa_biplot.rds for test-pcoa_biplot.r, so
that the tests need no Python. Generated with scikit-bio 0.7.4 (BSD-3-Clause),
numpy 2.5.5, scipy 1.18.1, pandas 3.1.1:

    python data-raw/expected/pcoa_biplot.py
    Rscript data-raw/expected/pcoa_biplot.r

Inputs are written too, so that R reads exactly the distances and tables
scikit-bio was given. Every number is written with 17 significant digits.

scikit-bio has no Lingoes or Cailliez correction, so those are tested
against ape, vegan, and cmdscale() instead (data-raw/expected/pcoa_correct.r).
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa, pcoa_biplot


OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "pcoa_biplot")
AXES = 5


def write(df, name):
    df.to_csv(os.path.join(OUT, name), float_format="%.17g", lineterminator="\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(20260116)

    # Samples x features: sparse counts, with one all-zero feature and one
    # feature with the same count everywhere.
    n, p = 24, 12
    counts = rng.poisson(4, size=(n, p)) * rng.binomial(1, 0.5, size=(n, p))
    counts[:, 0] += 1
    counts[:, 10] = 0
    counts[:, 11] = 7
    samples = ["S%02d" % (i + 1) for i in range(n)]
    features = ["OTU%02d" % (j + 1) for j in range(p)]
    counts = pd.DataFrame(counts, index=samples, columns=features)

    # Continuous variables on different scales, one negative.
    env = pd.DataFrame(
        {
            "pH": rng.normal(6.5, 0.8, n),
            "temp": rng.normal(290, 5, n),
            "depth": rng.gamma(2, 3, n),
            "redox": rng.normal(-50, 40, n),
        },
        index=samples,
    )

    write(counts, "skbio_counts.csv")
    write(env, "skbio_env.csv")

    percent = counts.div(counts.sum(axis=1), axis=0)

    dists = {
        "bray": pdist(counts.values, "braycurtis"),
        "jaccard": pdist(counts.values > 0, "jaccard"),
        "euclidean": pdist(counts.values, "euclidean"),
    }

    for name, d in dists.items():

        dm = DistanceMatrix(squareform(d), ids=samples)
        write(pd.DataFrame(dm.data, index=samples, columns=samples), "skbio_%s_dist.csv" % name)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ord = pcoa(dm, dimensions=AXES)

        axes = ["PC%d" % (i + 1) for i in range(AXES)]
        write(
            pd.DataFrame(
                {
                    "eig": np.asarray(ord.eigvals),
                    "prop": np.asarray(ord.proportion_explained),
                },
                index=axes,
            ),
            "skbio_%s_eig.csv" % name,
        )
        write(ord.samples, "skbio_%s_coords.csv" % name)

        if name == "jaccard":
            continue

        for label, y in (("counts", counts), ("percent", percent), ("env", env)):
            scores = pcoa_biplot(ord, y).features
            write(scores, "skbio_%s_biplot_%s.csv" % (name, label))


if __name__ == "__main__":
    main()
