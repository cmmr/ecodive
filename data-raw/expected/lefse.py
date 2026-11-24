# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit

"""Reference values for lefse() from LEfSe 1.1.2.

Reads the cases that data-raw/expected/lefse.r writes to
data-raw/expected/lefse/ (cases.csv, and each input's counts
<input>_counts.csv and classes <input>_meta.csv, the samples sorted by class
and subclass as lefse_format_input.py sorts them), and runs LEfSe's own
functions on them as lefse_run.py does: numerical_values() (each sample
scaled to 1e6), get_class_means(), test_kw_r(), test_rep_wilcoxon_r() and
test_lda_r() (with LEfSe's tol, 1e-10). It records the random draws:
test_lda_r() jitters the values in place, and every call of
contast_within_classes_or_few_per_class() is recorded, so that each
bootstrap's samples are known. Writes, for each case <name>:

  <name>_norm.csv      the scaled table (samples x features)
  <name>_means.csv     each feature's class means (get_class_means())
  <name>_kw.csv        each feature's Kruskal-Wallis p-value
  <name>_wilcoxon.csv  the Wilcoxon step's result (1, 0, or NA if not run)
  <name>_jitter.csv    the jittered values of the features scored
  <name>_boots.csv     each bootstrap's samples (0-based), one per row
  <name>_lda.csv       test_lda_r()'s scores, every feature scored

with Python's repr() of every number (exact). data-raw/expected/lefse.r
folds them into tests/tinytests/expected/lefse.rds, so that the tests need
no Python.

Prerequisites: Python 3 with LEfSe 1.1.2 (pip install lefse==1.1.2), numpy
and rpy2 (3.5.17 was used, with R 4.3.3), and R with coin, MASS, survival, mvtnorm and
modeltools on R_LIBS. PYTHONHASHSEED=0 makes the set orders that LEfSe's
random draws follow reproducible:

    Rscript data-raw/expected/lefse.r
    PYTHONHASHSEED=0 python data-raw/expected/lefse.py
    Rscript data-raw/expected/lefse.r
"""

import csv
import os
import sys

import rpy2.robjects as robjects
from lefse import lefse as L
from lefse.lefse_format_input import get_class_slices, numerical_values, rename_same_subcl

# rpy2 3.5 returns None for R's invisible results (assignments such as
# `mm <- z$means`), which LEfSe 1.1.2's test_lda_r() reads; earlier rpy2
# returned them, as this setting does.
robjects.r._invisible = False


DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data-raw", "expected", "lefse")


def read_csv(path):
    with open(path) as f:
        return list(csv.reader(f))


def write_rows(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        if header is not None:
            w.writerow(header)
        for r in rows:
            w.writerow([repr(v) if isinstance(v, float) else v for v in r])


def run_case(case):

    name = case["name"]
    inp  = case["input"]

    counts = read_csv(os.path.join(DIR, inp + "_counts.csv"))
    meta   = read_csv(os.path.join(DIR, inp + "_meta.csv"))
    fnames = counts[0]
    rows   = counts[1:]
    n      = len(rows)

    feats = dict((fn, [r[j] for r in rows]) for j, fn in enumerate(fnames))
    feats = numerical_values(feats, 1000000.0)

    cls = {"class": [m[0] for m in meta[1:]]}
    sub = [m[1] for m in meta[1:]]
    if all(s == "" for s in sub):
        cls["subclass"] = [str(c) + "_subcl" for c in cls["class"]]
    else:
        cls["subclass"] = sub
    cls["subclass"] = rename_same_subcl(cls["class"], cls["subclass"])
    class_sl, subclass_sl, class_hierarchy = get_class_slices(
        list(zip(cls["class"], cls["subclass"], [None] * n)))

    write_rows(os.path.join(DIR, name + "_norm.csv"), fnames,
               [[feats[fn][i] for fn in fnames] for i in range(n)])

    L.init()
    kord, cls_means = L.get_class_means(class_sl, feats)
    write_rows(os.path.join(DIR, name + "_means.csv"), ["feature"] + list(kord),
               [[fn] + [float(v) for v in cls_means[fn]] for fn in fnames])

    alpha = 0.05
    wilc = int(case["wilc"])
    kw_rows, wil_rows = [], []
    for fn in list(feats.keys()):
        ok, pv = L.test_kw_r(cls, feats[fn], alpha, sorted(cls.keys()))
        kw_rows.append([fn, float(pv)])
        if not ok:
            del feats[fn]
            wil_rows.append([fn, "NA"])
            continue
        if not wilc:
            wil_rows.append([fn, "NA"])
            continue
        res = L.test_rep_wilcoxon_r(
            subclass_sl, class_hierarchy, feats[fn], 0.05, int(case["multiclass"]), 0, fn,
            int(case["min_c"]), int(case["same_subcl"]), False)
        wil_rows.append([fn, 1 if res else 0])
        if not res:
            del feats[fn]
    write_rows(os.path.join(DIR, name + "_kw.csv"), ["feature", "p"], kw_rows)
    write_rows(os.path.join(DIR, name + "_wilcoxon.csv"), ["feature", "pass"], wil_rows)

    # Record every subsample tried, and which the bootstraps kept.
    tries = []
    orig = L.contast_within_classes_or_few_per_class

    def recording(feats_, inds, min_cl, ncl):
        r = orig(feats_, inds, min_cl, ncl)
        tries.append((list(inds), r))
        return r

    L.contast_within_classes_or_few_per_class = recording
    scored = list(feats.keys())
    boots = int(case["n_boots"])
    lda = {}
    if scored:
        lda, _ = L.test_lda_r(cls, feats, class_sl, boots, float(case["f"]), 2.0, 0.0000000001, 3)
    L.contast_within_classes_or_few_per_class = orig

    kept, cur = [], 0
    for inds, r in tries:
        cur += 1
        if not r or cur == 1000:
            kept.append(inds)
            cur = 0
    assert len(kept) == (boots if scored else 0), (name, len(kept))

    write_rows(os.path.join(DIR, name + "_jitter.csv"), scored,
               [[feats[fn][i] for fn in scored] for i in range(n)])
    write_rows(os.path.join(DIR, name + "_boots.csv"), None, kept)
    write_rows(os.path.join(DIR, name + "_lda.csv"), ["feature", "score"],
               [[fn, float(lda[fn])] for fn in scored])


if __name__ == "__main__":
    for case in csv.DictReader(open(os.path.join(DIR, "cases.csv"))):
        run_case(case)
        print("done", case["name"])
