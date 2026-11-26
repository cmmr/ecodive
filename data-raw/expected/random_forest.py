# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for random_forest() from scikit-learn 1.9.1 and
q2-sample-classifier 2026.10.0 (both BSD-3).

Run after the first pass of data-raw/expected/random_forest.r, from the
package root:

  Q2SC_DIR=/path/to/q2-sample-classifier \\
    python -I data-raw/expected/random_forest.py

scikit-learn must be 1.9.1 (its class weights drive a weighted bootstrap
since 1.9) and numpy >= 2 (randint's default integer is int64 everywhere).

q2-sample-classifier's _fit_and_predict_cv() and its ROC helpers live in
q2_sample_classifier/utilities.py and visuals.py. This script loads those
two files from a checkout of the 2026.10.0 tag (Q2SC_DIR) and calls the
functions directly, so no QIIME 2 install is needed:

  * The package's __init__.py (which imports the QIIME 2 framework) is not
    run: an empty q2_sample_classifier package object stands in for it.
  * qiime2.plugin (only get_available_cores() is imported) and q2templates
    (HTML templating, used only by the visualizers) are stand-ins.
  * Its other imports are the real packages: scikit-learn, pandas, numpy,
    scipy, biom-format, seaborn and matplotlib.

Reads data-raw/expected/random_forest/{bin,three,reg}.csv and writes, to the
same folder, text files of "key<TAB>values" lines (values separated by
spaces; real numbers as float.hex()):

  cases.txt          The case names, one per line.
  case-<name>.txt    One fitted forest (see random_forest.r).
  rng.txt            numpy's legacy RandomState draws, our_rand_r(), and
                     np.sum() of some arrays.
  folds.txt          StratifiedKFold and KFold test folds (1-based).
  roc.txt            roc_curve() and auc() on stored scores.
  q2-cls.txt, q2-reg.txt   q2's nested cross-validation.
  dist.txt           OOB accuracy and importances over 20 random_states.
  versions.txt       Package versions.
"""

import importlib
import importlib.machinery
import importlib.util
import os
import sys
import types
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
import sklearn
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import roc_curve, auc, accuracy_score, mean_squared_error

assert sklearn.__version__ == "1.9.1", sklearn.__version__
assert int(np.__version__.split(".")[0]) >= 2, np.__version__

DIR = os.path.join("data-raw", "expected", "random_forest")
Q2SC = os.environ["Q2SC_DIR"]


def hexf(v):
    v = float(v)
    if np.isnan(v):
        return "nan"
    if np.isinf(v):
        return "inf" if v > 0 else "-inf"
    return v.hex()


def unhex(v):
    return np.nan if v in ("NA", "nan", "") else float.fromhex(v)


def line(key, values):
    return key + "\t" + " ".join(values) + "\n"


def hexline(key, values):
    return line(key, [hexf(v) for v in np.ravel(values)])


def intline(key, values):
    return line(key, [str(int(v)) for v in np.ravel(values)])


def read_table(name, numeric_y=False):
    df = pd.read_csv(os.path.join(DIR, name), dtype=str, keep_default_na=False)
    x = np.array([[unhex(v) for v in row] for row in df.iloc[:, 2:].values])
    y = np.array([unhex(v) for v in df["y"]]) if numeric_y else df["y"].values
    return df["sample"].values, list(df.columns[2:]), x, y


data = {
    "bin":   read_table("bin.csv"),
    "three": read_table("three.csv"),
    "reg":   read_table("reg.csv", numeric_y=True),
}


#----------------------------------------------------------------------
# Forests.
#----------------------------------------------------------------------

def fit_case(name, dname, sparse, params, random_state, n_trees, out):
    samples, feats, x, y = data[dname]
    reg = dname == "reg"
    kw = dict(params)
    if "class_weight" in kw and isinstance(kw["class_weight"], dict):
        pass
    est = (RandomForestRegressor if reg else RandomForestClassifier)(
        n_estimators=n_trees, random_state=random_state, oob_score=True,
        n_jobs=1, **kw)
    X = sp.csc_matrix(x) if sparse else x
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est.fit(X, y)

    k = 1 if reg else len(est.classes_)
    lines = [line("data", [dname]), line("sparse", [str(sparse)]),
             intline("k", [k]), intline("n_trees", [n_trees]),
             intline("random_state", [random_state])]
    for key, v in params.items():
        if isinstance(v, dict):
            vals = []
            for c, w in v.items():
                vals += [str(c), hexf(w)]
            lines.append(line("param_" + key, vals))
        else:
            lines.append(line("param_" + key, [str(v)]))

    offsets, cols = [0], {c: [] for c in [
        "left", "right", "feature", "threshold", "value", "impurity",
        "n_node_samples", "weighted_n_node_samples"]}
    for t in est.estimators_:
        tr = t.tree_
        offsets.append(offsets[-1] + tr.node_count)
        cols["left"].append(tr.children_left)
        cols["right"].append(tr.children_right)
        cols["feature"].append(tr.feature)
        cols["threshold"].append(tr.threshold)
        cols["value"].append(tr.value[:, 0, :].ravel())
        cols["impurity"].append(tr.impurity)
        cols["n_node_samples"].append(tr.n_node_samples)
        cols["weighted_n_node_samples"].append(tr.weighted_n_node_samples)
    lines.append(intline("offsets", offsets))
    for c in ["left", "right", "feature", "n_node_samples"]:
        lines.append(intline(c, np.concatenate(cols[c])))
    for c in ["threshold", "value", "impurity", "weighted_n_node_samples"]:
        lines.append(hexline(c, np.concatenate(cols[c])))

    n = x.shape[0]
    inbag = np.column_stack([np.bincount(s, minlength=n) for s in est.estimators_samples_])
    lines.append(intline("inbag", inbag.T.ravel()))   # n x trees, column-major
    lines.append(hexline("importance", est.feature_importances_))
    if reg:
        lines.append(hexline("pred", est.predict(X)))
        lines.append(hexline("oob", est.oob_prediction_))
    else:
        lines.append(hexline("pred", est.predict_proba(X).ravel()))
        lines.append(hexline("oob", est.oob_decision_function_.ravel()))
    lines.append(hexline("oob_score", [est.oob_score_]))

    with open(os.path.join(DIR, "case-" + name + ".txt"), "w") as f:
        f.writelines(lines)
    out.append(name)


cases = []
B, T, R = "bin", "three", "reg"
fit_case("bin_default",   B, False, {}, 7, 30, cases)
fit_case("bin_sparse",    B, True,  {}, 7, 30, cases)
fit_case("bin_entropy",   B, False, {"criterion": "entropy"}, 11, 20, cases)
fit_case("bin_balanced",  B, False, {"class_weight": "balanced"}, 3, 20, cases)
fit_case("bin_balsub",    B, False, {"class_weight": "balanced_subsample"}, 4, 20, cases)
fit_case("bin_dict",      B, False, {"class_weight": {"A": 1.0, "B": 3.0}}, 5, 20, cases)
fit_case("bin_maxsamp",   B, False, {"max_samples": 0.7}, 6, 20, cases)
fit_case("bin_leaf3",     B, False, {"min_samples_leaf": 3}, 8, 20, cases)
fit_case("bin_split5",    B, False, {"min_samples_split": 5}, 15, 20, cases)
fit_case("bin_depth3",    B, False, {"max_depth": 3}, 9, 20, cases)
fit_case("bin_mf1",       B, False, {"max_features": 1}, 10, 20, cases)
fit_case("bin_log2",      B, False, {"max_features": "log2"}, 12, 20, cases)
fit_case("bin_mf05",      B, False, {"max_features": 0.5}, 13, 20, cases)
fit_case("three_default", T, False, {}, 21, 25, cases)
fit_case("three_sparse",  T, True,  {}, 21, 25, cases)
fit_case("three_balanced_sparse", T, True,
         {"class_weight": "balanced", "max_samples": 0.8}, 22, 20, cases)
fit_case("three_balsub_sparse", T, True,
         {"class_weight": "balanced_subsample", "criterion": "entropy"}, 23, 20, cases)
fit_case("reg_default",   R, False, {}, 31, 20, cases)
fit_case("reg_sparse",    R, True,  {}, 31, 20, cases)
fit_case("reg_mf05",      R, False, {"max_features": 0.5}, 32, 20, cases)
fit_case("reg_leaf3",     R, False, {"min_samples_leaf": 3, "max_samples": 20}, 33, 20, cases)

with open(os.path.join(DIR, "cases.txt"), "w") as f:
    f.writelines(c + "\n" for c in cases)


#----------------------------------------------------------------------
# Draws.
#----------------------------------------------------------------------

def our_rand_r(state):
    # sklearn/utils/_random.pxd (BSD-3): a 32-bit xorshift.
    s = state[0] or 1
    s ^= (s << 13) & 0xFFFFFFFF
    s ^= s >> 17
    s ^= (s << 5) & 0xFFFFFFFF
    state[0] = s
    return s % (2147483647 + 1)


lines = [hexline("random_sample", np.random.RandomState(0).random_sample(1000))]
for n in [7, 60, 1000, 2147483647]:
    lines.append(hexline("randint_%d" % n, np.random.RandomState(1).randint(0, n, 1000)))
rs = np.random.RandomState(7)
lines.append(hexline("tree_seeds", [rs.randint(np.iinfo(np.int32).max) for _ in range(1000)]))
st = [12345]
lines.append(hexline("our_rand_r", [our_rand_r(st) for _ in range(1000)]))
g = np.random.default_rng(5)
for i, n in enumerate([1, 7, 8, 9, 100, 128, 129, 300, 1000]):
    a = g.random(n) * 10.0 ** g.integers(-3, 4, n)
    lines.append(hexline("np_sum_in_%d" % i, a))
    lines.append(hexline("np_sum_out_%d" % i, [np.sum(a)]))
with open(os.path.join(DIR, "rng.txt"), "w") as f:
    f.writelines(lines)


#----------------------------------------------------------------------
# Folds.
#----------------------------------------------------------------------

def test_folds(cv, x, y):
    out = np.zeros(len(y), dtype=int)
    for f, (_, test) in enumerate(cv.split(x, y)):
        out[test] = f + 1
    return out


with open(os.path.join(DIR, "folds.txt"), "w") as f:
    f.write(intline("three", test_folds(
        StratifiedKFold(5, shuffle=True, random_state=3), data["three"][2], data["three"][3])))
    f.write(intline("reg", test_folds(
        KFold(5, shuffle=True, random_state=3), data["reg"][2], data["reg"][3])))


#----------------------------------------------------------------------
# q2-sample-classifier.
#----------------------------------------------------------------------

def load_q2():
    pkg_dir = os.path.join(Q2SC, "q2_sample_classifier")
    pkg = types.ModuleType("q2_sample_classifier")
    pkg.__path__ = [pkg_dir]
    pkg.__spec__ = importlib.machinery.ModuleSpec(
        "q2_sample_classifier", None, is_package=True)
    pkg.__spec__.submodule_search_locations = [pkg_dir]
    sys.modules["q2_sample_classifier"] = pkg
    qiime2 = types.ModuleType("qiime2")
    plugin = types.ModuleType("qiime2.plugin")
    plugin.get_available_cores = lambda: 1
    qiime2.plugin = plugin
    sys.modules["qiime2"] = qiime2
    sys.modules["qiime2.plugin"] = plugin
    sys.modules["q2templates"] = types.ModuleType("q2templates")
    mods = {}
    for m in ["visuals", "utilities"]:
        spec = importlib.util.spec_from_file_location(
            "q2_sample_classifier." + m, os.path.join(pkg_dir, m + ".py"))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        mods[m] = mod
    return mods["utilities"], mods["visuals"]


import biom  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.feature_extraction import DictVectorizer  # noqa: E402

u, vis = load_q2()


def run_q2(dname, file, classification):
    samples, feats, x, y = data[dname]
    table = biom.Table(x.T, observation_ids=feats, sample_ids=list(samples))
    meta = pd.Series(y, index=list(samples), name="y")
    est = (RandomForestClassifier if classification else RandomForestRegressor)(
        n_jobs=1, n_estimators=20, random_state=5)
    est = Pipeline([("dv", DictVectorizer()), ("est", est)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scores, predictions, importances, tops, probabilities = u._fit_and_predict_cv(
            table, meta, est, None, 1,
            scoring=accuracy_score if classification else mean_squared_error,
            random_state=5, cv=5, stratify=classification,
            calc_feature_importance=True, parameter_tuning=False)
    cv = (StratifiedKFold if classification else KFold)(5, shuffle=True, random_state=5)
    features = u._extract_features(table)
    folds = test_folds(cv, features, meta)

    pred = predictions["prediction"].reindex(list(samples))
    lines = [intline("folds", folds), line("sample", list(samples))]
    if classification:
        lines.append(line("predictions", [str(v) for v in pred.values]))
        probs = probabilities.reindex(list(samples))
        classes = [str(c) for c in probs.columns]
        lines.append(line("classes", classes))
        lines.append(hexline("probabilities", probs.values.T.ravel()))  # column-major
        truth = meta.reindex(list(samples))
        bt = vis._binarize_labels(truth, probs.columns)
        fpr, tpr, roc_auc = vis._roc_per_class(bt, probs.values, probs.columns)
        fpr, tpr, roc_auc = vis._roc_micro_average(bt, probs.values, fpr, tpr, roc_auc)
        fpr, tpr, roc_auc = vis._roc_macro_average(fpr, tpr, roc_auc, probs.columns)
        for c in classes + ["micro", "macro"]:
            lines.append(hexline("roc_%s_fpr" % c, fpr[c]))
            lines.append(hexline("roc_%s_tpr" % c, tpr[c]))
            lines.append(hexline("roc_%s_auc" % c, [roc_auc[c]]))
    else:
        lines.append(hexline("predictions", pred.values))
    lines.append(line("importance_features", [str(i) for i in importances.index]))
    lines.append(hexline("importances", importances["importance"].values))
    lines.append(hexline("scores", scores))
    with open(os.path.join(DIR, file), "w") as f:
        f.writelines(lines)


run_q2("three", "q2-cls.txt", True)
run_q2("reg", "q2-reg.txt", False)


#----------------------------------------------------------------------
# ROC on stored scores.
#----------------------------------------------------------------------

g = np.random.default_rng(11)
truth = g.integers(0, 2, 40)
score = np.round(g.random(40) * 0.5 + truth * 0.3, 1)   # ties
lines = [hexline("bin_truth", truth), hexline("bin_score", score)]
fpr, tpr, thr = roc_curve(truth, score)
lines += [hexline("bin_fpr", fpr), hexline("bin_tpr", tpr),
          hexline("bin_thresholds", thr), hexline("bin_auc", [auc(fpr, tpr)])]
with open(os.path.join(DIR, "roc.txt"), "w") as f:
    f.writelines(lines)


#----------------------------------------------------------------------
# Distributional references: 500-tree forests on bin over 20 random_states.
#----------------------------------------------------------------------

samples, feats, x, y = data["bin"]
acc, imp = [], []
for s in range(20):
    est = RandomForestClassifier(n_estimators=500, random_state=s, oob_score=True, n_jobs=1)
    est.fit(x, y)
    acc.append(est.oob_score_)
    imp.append(est.feature_importances_)
with open(os.path.join(DIR, "dist.txt"), "w") as f:
    f.write(hexline("oob_accuracy", acc))
    f.write(hexline("importance", np.column_stack(imp).T.ravel()))   # p x 20, column-major

import scipy  # noqa: E402
import matplotlib  # noqa: E402
import seaborn  # noqa: E402
with open(os.path.join(DIR, "versions.txt"), "w") as f:
    for m in [sklearn, np, scipy, pd, biom, matplotlib, seaborn]:
        f.write(line(m.__name__, [m.__version__]))
    f.write(line("python", [sys.version.split()[0]]))
    f.write(line("q2_sample_classifier", ["2026.10.0"]))
