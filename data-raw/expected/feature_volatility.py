# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for feature_volatility() from q2-longitudinal 2026.10.0
and q2-sample-classifier 2026.10.0 (both BSD-3, QIIME 2 development team),
with scikit-learn 1.9.1.

q2-longitudinal's feature-volatility pipeline runs q2-sample-classifier's
regress-samples (split_table, then fit_regressor, then predict_regression and
the accuracy of the test samples' predictions), converts the table to
relative frequencies, and passes it with the importances to
plot-feature-volatility, which keeps the most important features
(_importance_filtering()) and summarizes each (_summarize_feature_stats()).
This script calls those same functions directly, so no QIIME 2 install is
needed:

  * q2_longitudinal/_utilities.py (Q2L_DIR, a checkout of its 2026.10.0
    tag) and q2_sample_classifier/utilities.py and visuals.py (Q2SC_DIR, a
    checkout of its 2026.10.0 tag) are loaded on their own, with empty
    package objects in place of the packages' __init__.py (which import the
    QIIME 2 framework).
  * qiime2.plugin (only get_available_cores() is imported) and q2templates
    (HTML templating, used only by the visualizers) are stand-ins, and so is
    qiime2.Metadata: an object whose to_dataframe() and to_series() return
    the state column, which is all _load_data() and _fit_estimator() call.
  * Their other imports are the real packages.

The pipeline's steps that this script does in their place: the table's
relative frequencies (feature-table relative-frequency: biom's
Table.norm(axis='sample')), and its filtering to the importances'
features.

Two runs with optimize_feature_selection = False (ecodive has no recursive
feature elimination), random_state 3 and 7, and one with the pipeline's own
optimize_feature_selection = True (random_state 3), for information.

Reads data-raw/expected/feature_volatility/{meta,counts}.csv and writes, to
the same folder, text files of "key<TAB>values" lines (values separated by
spaces; real numbers as float.hex()):

  case-rs<random_state>.txt   One run (see feature_volatility.r).
  rfe.txt                     The RFE run's importances.
  perm.txt                    numpy's RandomState(rs).permutation(n).
  versions.txt                Package versions.

Run from the repository root (see data-raw/expected/feature_volatility.r):

    Q2L_DIR=/path/to/q2-longitudinal Q2SC_DIR=/path/to/q2-sample-classifier \\
      python -I data-raw/expected/feature_volatility.py
"""

import importlib.machinery
import importlib.util
import os
import sys
import types
import warnings

import numpy as np
import pandas as pd
import sklearn

assert sklearn.__version__ == "1.9.1", sklearn.__version__
assert int(np.__version__.split(".")[0]) >= 2, np.__version__

DIR = os.path.join("data-raw", "expected", "feature_volatility")
Q2L = os.environ["Q2L_DIR"]
Q2SC = os.environ["Q2SC_DIR"]


def hexf(v):
    v = float(v)
    return "nan" if np.isnan(v) else v.hex()


def line(key, values):
    return key + "\t" + " ".join(values) + "\n"


def hexline(key, values):
    return line(key, [hexf(v) for v in values])


def load(pkg_name, pkg_dir, mods):
    pkg = types.ModuleType(pkg_name)
    pkg.__path__ = [pkg_dir]
    pkg.__spec__ = importlib.machinery.ModuleSpec(pkg_name, None, is_package=True)
    pkg.__spec__.submodule_search_locations = [pkg_dir]
    sys.modules[pkg_name] = pkg
    out = []
    for m in mods:
        spec = importlib.util.spec_from_file_location(
            pkg_name + "." + m, os.path.join(pkg_dir, m + ".py"))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        out.append(mod)
    return out


qiime2 = types.ModuleType("qiime2")
plugin = types.ModuleType("qiime2.plugin")
plugin.get_available_cores = lambda: 1
qiime2.plugin = plugin
sys.modules["qiime2"] = qiime2
sys.modules["qiime2.plugin"] = plugin
sys.modules["q2templates"] = types.ModuleType("q2templates")

import biom  # noqa: E402

lu, = load("q2_longitudinal", os.path.join(Q2L, "q2_longitudinal"), ["_utilities"])
vis, su = load("q2_sample_classifier", os.path.join(Q2SC, "q2_sample_classifier"),
               ["visuals", "utilities"])


class StateColumn:
    """What _load_data() and _fit_estimator() use of a MetadataColumn."""

    def __init__(self, series):
        self.series = series

    def to_series(self):
        return self.series.copy()

    def to_dataframe(self):
        return self.series.to_frame()


meta = pd.read_csv(os.path.join(DIR, "meta.csv"), dtype=str).set_index("sample")
meta.index.name = "id"
meta["state"] = meta["state"].astype(float)
counts = pd.read_csv(os.path.join(DIR, "counts.csv"), dtype=str).set_index("sample")
counts = counts.astype(float)
counts.index.name = "id"

samples = list(counts.index)
features = list(counts.columns)
table = biom.Table(counts.values.T, observation_ids=features, sample_ids=samples)
state = StateColumn(meta["state"])
state_md_col = meta[["state"]]

# feature-table relative-frequency.
rel = table.norm(axis="sample", inplace=False)
rel_df = rel.to_dataframe(dense=True).T


def regress(rs, rfe):
    """regress_samples(): split_table, fit_regressor, predict, accuracy."""
    X_train, X_test, y_train, y_test = su._prepare_training_data(
        table, state, "state", 0.2, rs, load_data=True, stratify=False,
        missing_samples="ignore")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est, imp = su._fit_estimator(
            X_train, state, "RandomForestRegressor", n_estimators=100,
            step=0.05, cv=5, random_state=rs, n_jobs=1,
            optimize_feature_selection=rfe, parameter_tuning=False,
            missing_samples="ignore", classification=False)
    pred = est.predict(su._extract_features(X_test))
    pred = pd.Series(pred, index=y_test.index, name="prediction")
    acc = vis._linear_regress(y_test, pred)
    return X_train, X_test, y_train, y_test, imp, pred, acc


def stats(imp, feature_count):
    """plot_feature_volatility()'s filtering and _summarize_feature_stats()."""
    # feature_volatility(): the relative table, filtered to the importances' features.
    tab = rel_df[[f for f in rel_df.columns if f in set(imp.index)]]
    tab, imps = lu._importance_filtering(tab, imp.copy(), None, feature_count)
    st = lu._summarize_feature_stats(tab, state_md_col)
    st = imps.join(st, how="inner")
    return st


def stats_lines(prefix, st):
    out = [line(prefix + "_feature", [str(i) for i in st.index]),
           hexline(prefix + "_importance", st["importance"].values)]
    for c in st.columns:
        if c != "importance":
            out.append(hexline(prefix + "_" + c, st[c].values))
    return out


for rs in [3, 7]:
    X_train, X_test, y_train, y_test, imp, pred, acc = regress(rs, False)
    assert list(X_train.ids()) == list(y_train.index)
    out = [line("random_state", [str(rs)]),
           line("train", list(X_train.ids())),
           line("test", list(X_test.ids())),
           line("importance_features", [str(i) for i in imp.index]),
           hexline("importance", imp["importance"].values),
           hexline("predicted", pred.values),
           hexline("accuracy", acc[["Mean squared error", "r-value", "r-squared",
                                     "Slope", "Intercept"]].values[0])]
    out += stats_lines("all", stats(imp, "all"))
    out += stats_lines("top10", stats(imp, 10))
    with open(os.path.join(DIR, "case-rs%d.txt" % rs), "w") as f:
        f.writelines(out)

_, _, _, _, imp, _, _ = regress(3, True)
with open(os.path.join(DIR, "rfe.txt"), "w") as f:
    f.write(line("importance_features", [str(i) for i in imp.index]))
    f.write(hexline("importance", imp["importance"].values))

with open(os.path.join(DIR, "perm.txt"), "w") as f:
    for rs, n in [(3, 7), (3, 78), (7, 78), (11, 1000), (2147483646, 300)]:
        p = np.random.RandomState(rs).permutation(n)
        f.write(line("perm_%d_%d" % (rs, n), [str(int(v)) for v in p]))

import scipy  # noqa: E402
with open(os.path.join(DIR, "versions.txt"), "w") as f:
    for m in [sklearn, np, scipy, pd, biom]:
        f.write(line(m.__name__, [m.__version__]))
    f.write(line("python", [sys.version.split()[0]]))
    f.write(line("q2_longitudinal", ["2026.10.0"]))
    f.write(line("q2_sample_classifier", ["2026.10.0"]))
