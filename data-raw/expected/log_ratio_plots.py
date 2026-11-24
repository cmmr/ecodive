"""Qurro 0.9.0's outputs for ecodive's log_ratio_features() and log_ratio().

Step 2 of 3 (see data-raw/expected/log_ratio_plots.r). Reads the tables and
rankings that step 1 wrote to data-raw/expected/log_ratio_plots/, and writes
there:

  qurro.csv   One row per case: id; features, the rank plot's features in
              the order of its dataset (Qurro's process_input()); top and
              bot, the features of Qurro's autoselection of n from each end
              of the ranking (feature_computation.js's filterFeatures() with
              "autoLiteralTop" and "autoLiteralBot", run under Node.js by
              log_ratio_plots.js), in the order Qurro returns them. IDs are
              separated by spaces.

  values.csv  Long format (id, sample, qarcoal, balance), hex floats or NA:
                qarcoal  qarcoal()'s log ratio of the summed counts of top
                         over bot, log(Num_Sum / Denom_Sum); NA for the
                         samples it drops (a zero sum)
                balance  the sample plot's value: display.js's sums over
                         top and bot, then computeBalance(), log(top) -
                         log(bot); NA where it is null

Prerequisites: Python 3 with numpy, pandas, biom-format and qurro 0.9.0, and
Node.js (any version with ES2015). qarcoal imports qiime2 only for a type
annotation, so a stub module stands in for it (as in log_ratio.py). Run from
the repository root:

    python data-raw/expected/log_ratio_plots.py
"""

import csv
import json
import os
import subprocess
import sys
import tempfile
import types

import numpy as np
import pandas as pd
import biom

_q2 = types.ModuleType("qiime2")
_q2.Metadata = object
sys.modules.setdefault("qiime2", _q2)

import qurro
from qurro.qarcoal import qarcoal
from qurro.generate import process_input

assert qurro.__version__ == "0.9.0", qurro.__version__

DIR = os.path.join("data-raw", "expected", "log_ratio_plots")
JS = os.path.join("data-raw", "expected", "log_ratio_plots.js")
QURRO_DIR = os.path.dirname(qurro.__file__)


def read_table(name):
    d = pd.read_csv(os.path.join(DIR, f"table-{name}.csv"), index_col=0)
    return d.astype(int)


def hexs(v):
    return "NA" if v is None or not np.isfinite(v) else float(v).hex()


tables = {nm: read_table(nm) for nm in ("syn", "da")}
rk = pd.read_csv(os.path.join(DIR, "rankings.csv"), dtype=str)
case_ids = list(dict.fromkeys(rk["id"]))

cases = []
for cid in case_ids:
    r = rk[rk["id"] == cid]
    tab = tables[r["table"].iloc[0]]
    feats = list(r["feature"])
    ranks = pd.DataFrame(
        {"r": [float.fromhex(v) for v in r["effect"]]}, index=pd.Index(feats)
    )

    # Qurro's input processing: the ranked features (and samples) its rank
    # plot keeps, in its dataset's order.
    bt = biom.Table(
        tab[feats].T.values.astype(float), feats, list(tab.index)
    )
    meta = pd.DataFrame({"x": ["a"] * tab.shape[0]}, index=list(tab.index))
    _, out_ranks, _, _, out_table = process_input(ranks, meta, bt)
    kept = list(out_ranks.index)

    counts = {f: {s: int(tab.at[s, f]) for s in tab.index} for f in kept}
    cases.append(
        {
            "id": cid,
            "n": int(r["n"].iloc[0]),
            "ranks": [{"id": f, "r": float(out_ranks.at[f, "r"])} for f in kept],
            "counts": counts,
            "samples": list(tab.index),
            "table": r["table"].iloc[0],
            "kept": kept,
        }
    )

with tempfile.TemporaryDirectory() as tmp:
    cin = os.path.join(tmp, "cases.json")
    cout = os.path.join(tmp, "out.json")
    with open(cin, "w") as fh:
        json.dump(
            [{k: c[k] for k in ("id", "n", "ranks", "counts", "samples")} for c in cases],
            fh,
        )
    subprocess.run(["node", JS, QURRO_DIR, cin, cout], check=True)
    with open(cout) as fh:
        js = {o["id"]: o for o in json.load(fh)}

with open(os.path.join(DIR, "qurro.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "features", "top", "bot"])
    for c in cases:
        o = js[c["id"]]
        w.writerow([c["id"], " ".join(c["kept"]), " ".join(o["top"]), " ".join(o["bot"])])

with open(os.path.join(DIR, "values.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "sample", "qarcoal", "balance"])
    for c in cases:
        o = js[c["id"]]
        tab = tables[c["table"]]
        feats = list(tab.columns)
        samples = list(tab.index)

        # qarcoal(): taxonomy strings mark the numerator and denominator.
        bt = biom.Table(tab.T.values.astype(float), feats, samples)
        tax = pd.DataFrame({"Taxon": ["other"] * len(feats)}, index=feats)
        tax.loc[o["top"], "Taxon"] = "num"
        tax.loc[o["bot"], "Taxon"] = "den"
        q = qarcoal(bt, tax, "num", "den")
        qv = q["log_ratio"] if "log_ratio" in q.columns else q.iloc[:, -1]

        for i, s in enumerate(samples):
            qa = float(qv[s]) if s in qv.index else None
            b = o["balance"][i]
            w.writerow([c["id"], s, hexs(qa), hexs(b)])

print(f"{len(cases)} cases written to {DIR}")
