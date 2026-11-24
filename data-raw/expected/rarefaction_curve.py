# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for rarefaction_curve() from NumPy and pandas.

Reads the inputs that data-raw/expected/rarefaction_curve.r writes to
data-raw/expected/rarefaction_curve/ and writes:

* numpy.csv: q2-diversity's depth grid, np.linspace(min_depth, max_depth,
  num=steps, dtype=int) with min_depth = 1, for each (max_depth, steps)
  in grid_cases.csv (columns case, 1-based, and depth).
* pandas.csv: q2-diversity's per-sample summary of its iterations,
  pandas.DataFrame.describe() with its percentiles, applied to each row of
  iters.csv (values as hexadecimal floats); the 50% percentile is kept, as
  a hexadecimal float (float.hex(), which R reads back exactly).
* mpmath.csv: for each row of rows.csv (table, sample, depth), with the
  sample's counts from counts.csv, Hurlbert's (1971) expected richness S,
  Heck et al.'s (1975) standard deviation sd, and the expected coverage
  (Chao & Jost 2012; at the sample's total, Chao et al. 2014's Eq. 12, and
  their Appendix G when f2 = 0), from the published formulas, computed
  with mpmath at 50 significant digits (absence probabilities from
  loggamma()), as hexadecimal floats.
* versions.txt: the NumPy, pandas and mpmath versions.

Prerequisites: Python 3 with NumPy 2.5.5, pandas 3.1.1 (BSD-3-Clause) and
mpmath 1.4.1 (BSD-3-Clause), the versions the fixtures were made with. Run from the repository root:

    Rscript data-raw/expected/rarefaction_curve.r
    python data-raw/expected/rarefaction_curve.py
    Rscript data-raw/expected/rarefaction_curve.r
"""

import csv
import os
import sys

from collections import Counter, defaultdict

import mpmath
import numpy as np
import pandas as pd

DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "rarefaction_curve")

# q2-diversity's _compute_summary() percentiles.
PERC = [0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

with open(os.path.join(DIR, "grid_cases.csv"), newline="") as f:
    cases = [(int(float(r["max_depth"])), int(float(r["steps"]))) for r in csv.DictReader(f)]

with open(os.path.join(DIR, "numpy.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["case", "depth"])
    for i, (max_depth, steps) in enumerate(cases, start=1):
        for d in np.linspace(1, max_depth, num=steps, dtype=int):
            w.writerow([i, int(d)])

with open(os.path.join(DIR, "iters.csv"), newline="") as f:
    rows = list(csv.reader(f))
names = [r[0] for r in rows[1:]]
vals = np.array([[float.fromhex(v) for v in r[1:]] for r in rows[1:]])
df = pd.DataFrame(vals, index=names)
summary = df.apply(lambda r: r.describe(percentiles=PERC), axis=1)

with open(os.path.join(DIR, "pandas.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["row", "median"])
    for name, v in zip(names, summary["50%"]):
        w.writerow([name, float(v).hex()])

# High-precision closed forms.
mpmath.mp.dps = 50

counts = defaultdict(list)
with open(os.path.join(DIR, "counts.csv"), newline="") as f:
    for r in csv.DictReader(f):
        counts[(r["table"], r["sample"])].append(int(r["count"]))

with open(os.path.join(DIR, "rows.csv"), newline="") as f:
    rows = [(r["table"], r["sample"], int(r["depth"])) for r in csv.DictReader(f)]


def closed_forms(x, d):
    N = sum(x)
    f = Counter(x)
    ks = sorted(f)
    lgN = mpmath.loggamma(N + 1) - mpmath.loggamma(N - d + 1)
    cache = {}

    def a(k):  # C(N - k, d) / C(N, d)
        if k not in cache:
            cache[k] = (mpmath.exp(mpmath.loggamma(N - k + 1) - mpmath.loggamma(N - k - d + 1) - lgN)
                        if N - k >= d else mpmath.mpf(0))
        return cache[k]

    S = len(x) - sum(f[k] * a(k) for k in ks)
    var = mpmath.mpf(0)
    for i, k in enumerate(ks):
        var += f[k] * a(k) * (1 - a(k)) + f[k] * (f[k] - 1) * (a(2 * k) - a(k) ** 2)
        for l in ks[i + 1:]:
            var += 2 * f[k] * f[l] * (a(k + l) - a(k) * a(l))
    if d < N:
        cov = 1 - sum(mpmath.mpf(f[k]) * k * a(k) for k in ks) / (N - d)
    else:
        f1, f2 = f.get(1, 0), f.get(2, 0)
        if f1 == 0:
            cov = mpmath.mpf(1)
        elif f2 > 0:
            cov = 1 - mpmath.mpf(f1) / N * ((N - 1) * f1) / ((N - 1) * f1 + 2 * f2)
        else:
            cov = 1 - mpmath.mpf(f1) / N * ((N - 1) * (f1 - 1)) / ((N - 1) * (f1 - 1) + 2)
    return S, mpmath.sqrt(var) if var > 0 else mpmath.mpf(0), cov


with open(os.path.join(DIR, "mpmath.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["table", "sample", "depth", "S", "sd", "coverage"])
    for table, sample, d in rows:
        S, sd, cov = closed_forms(counts[(table, sample)], d)
        w.writerow([table, sample, d] + [float(v).hex() for v in (S, sd, cov)])

with open(os.path.join(DIR, "versions.txt"), "w") as f:
    f.write(np.__version__ + "\n" + pd.__version__ + "\n" + mpmath.__version__ + "\n")
