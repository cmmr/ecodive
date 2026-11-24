# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for first_differences(), first_distances() and nmit()
from q2-longitudinal 2026.10.0 (BSD-3, QIIME 2 development team).

q2-longitudinal's actions are thin wrappers around functions in its
q2_longitudinal/_utilities.py that take pandas objects: _first_differences()
(first-differences and first-distances) and _nmit() (nmit). This script
loads that one file from a checkout of the 2026.10.0 tag (set Q2L_DIR to
it) and calls those functions directly, so no QIIME 2 install is needed:

  * The package's __init__.py (which imports the QIIME 2 framework) is not
    run: an empty q2_longitudinal package object stands in for it, and
    _utilities.py is loaded as its submodule.
  * q2templates (QIIME 2's HTML templating, imported at the top of
    _utilities.py but used only by the visualizers) is an empty stand-in.
  * Its other imports are the real packages: numpy, pandas, scipy,
    statsmodels, seaborn, matplotlib, scikit-bio (DistanceMatrix), biom-format
    and patsy.

The public actions also do the following, which this script does in their
place: first_differences() with a table adds the feature as a metadata
column (_add_metric_to_metadata()); nmit() keeps the metadata rows of the
table's samples.

Reads data-raw/expected/longitudinal/{meta,rel,dm}.csv (written by
data-raw/expected/longitudinal.r; real numbers as hexadecimal floats) and
writes, to the same folder (real numbers as float.hex(), "nan" for NaN):

  fd-drop.csv, fd-baseline.csv, fd-baseline7.csv
      _first_differences() of shannon, replicate_handling = 'drop';
      consecutive states, baseline = 0, baseline = 7. Columns sample (q2's
      index: the later sample), difference, in q2's order.
  fd-table.csv
      The same ('drop', consecutive) for each feature of rel.csv: feature,
      sample, difference.
  fdist-drop.csv, fdist-baseline.csv
      _first_differences(distance_matrix = dm), 'drop'; consecutive and
      baseline = 0: sample, distance.
  errors.csv
      The messages of the errors q2 raises: replicate_handling = 'error';
      a baseline that is not a state; one state only.
  nmit-<corr>-<dist>.csv
      _nmit(rel, meta, 'subject', corr, dist) for corr in kendall, pearson,
      spearman and dist in fro, nuc: the distance matrix, with q2's ids.
  nmit-ids.csv
      Each of those ids (a subject's first sample) and its subject.
  versions.csv

Run from the repository root (see data-raw/expected/longitudinal.r):

    Q2L_DIR=/path/to/q2-longitudinal python -I data-raw/expected/longitudinal.py
"""

import csv
import importlib.machinery
import importlib.util
import os
import sys
import types
import warnings

import numpy as np
import pandas as pd
import scipy
import skbio

Q2L = os.environ["Q2L_DIR"]
DIR = os.path.join("data-raw", "expected", "longitudinal")


def load_utilities():
    pkg_dir = os.path.join(Q2L, "q2_longitudinal")
    pkg = types.ModuleType("q2_longitudinal")
    pkg.__path__ = [pkg_dir]
    pkg.__spec__ = importlib.machinery.ModuleSpec(
        "q2_longitudinal", None, is_package=True)
    pkg.__spec__.submodule_search_locations = [pkg_dir]
    sys.modules["q2_longitudinal"] = pkg
    sys.modules["q2templates"] = types.ModuleType("q2templates")
    spec = importlib.util.spec_from_file_location(
        "q2_longitudinal._utilities", os.path.join(pkg_dir, "_utilities.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def unhex(v):
    return np.nan if v in ("NA", "nan", "") else float.fromhex(v)


def hexf(v):
    v = float(v)
    return "nan" if np.isnan(v) else v.hex()


def read_meta():
    df = pd.read_csv(os.path.join(DIR, "meta.csv"), dtype=str, keep_default_na=False)
    df = df.set_index("sample")
    df.index.name = "id"
    df["state"] = df["state"].astype(float)
    df["shannon"] = [unhex(v) for v in df["shannon"]]
    return df


def read_matrix(name):
    df = pd.read_csv(os.path.join(DIR, name), dtype=str, keep_default_na=False).set_index("sample")
    return df.map(unhex).astype(float)


def write(name, header, rows):
    with open(os.path.join(DIR, name), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main():
    u = load_utilities()
    meta = read_meta()
    rel = read_matrix("rel.csv")
    dm = read_matrix("dm.csv")
    dm = skbio.DistanceMatrix(dm.values, ids=list(dm.index))

    def fd(md, **kw):
        return u._first_differences(md.copy(), "state", "subject", **kw)

    for name, base in [("drop", None), ("baseline", 0), ("baseline7", 7)]:
        s = fd(meta, metric="shannon", replicate_handling="drop",
               baseline=base)
        write("fd-%s.csv" % name, ["sample", "difference"],
              [[i, hexf(v)] for i, v in s.items()])

    rows = []
    for feat in rel.columns:
        md = pd.concat([meta, rel[[feat]]], axis=1, join="inner")
        s = fd(md, metric=feat, replicate_handling="drop")
        rows += [[feat, i, hexf(v)] for i, v in s.items()]
    write("fd-table.csv", ["feature", "sample", "difference"], rows)

    for name, base in [("drop", None), ("baseline", 0)]:
        s = fd(meta, metric=None, replicate_handling="drop", baseline=base,
               distance_matrix=dm)
        write("fdist-%s.csv" % name, ["sample", "distance"],
              [[i, hexf(v)] for i, v in s.items()])

    errors = []
    for case, md, kw in [
            ("replicate", meta,
             dict(metric="shannon", replicate_handling="error")),
            ("baseline", meta,
             dict(metric="shannon", replicate_handling="drop", baseline=5)),
            ("one_state", meta[meta["state"] == 3],
             dict(metric="shannon", replicate_handling="drop"))]:
        try:
            fd(md, **kw)
            msg = "(no error)"
        except Exception as e:  # noqa: BLE001
            msg = "%s: %s" % (type(e).__name__, e)
        errors.append([case, msg])
    try:
        u._validate_input_columns(meta[meta["state"] == 3], "subject", None,
                                  "state", "shannon")
        errors.append(["validate_one_state", "(no error)"])
    except Exception as e:  # noqa: BLE001
        errors.append(["validate_one_state", "%s: %s" % (type(e).__name__, e)])
    write("errors.csv", ["case", "message"], errors)

    ids = None
    for corr in ["kendall", "pearson", "spearman"]:
        for dist in ["fro", "nuc"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d = u._nmit(rel.copy(), meta.copy(), "subject",
                            corr_method=corr, dist_method=dist)
            ids = list(d.ids)
            write("nmit-%s-%s.csv" % (corr, dist), ["id"] + ids,
                  [[ids[i]] + [hexf(v) for v in d.data[i]]
                   for i in range(len(ids))])
    write("nmit-ids.csv", ["id", "subject"],
          [[i, meta.loc[i, "subject"]] for i in ids])

    import importlib.metadata as md
    write("versions.csv", ["package", "version"],
          [["python", sys.version.split()[0]], ["numpy", np.__version__],
           ["pandas", pd.__version__], ["scipy", scipy.__version__],
           ["scikit-bio", skbio.__version__],
           ["statsmodels", md.version("statsmodels")],
           ["q2-longitudinal", "2026.10.0"]])


if __name__ == "__main__":
    main()
