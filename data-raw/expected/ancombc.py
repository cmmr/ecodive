"""Reference fixtures for ecodive's ancombc() and ancombc2().

Writes, under data-raw/expected/ancombc/ (then data-raw/expected/ancombc.r
converts them, with ANCOMBC's results, to tests/tinytests/expected/
ancombc.rds):

  table-<name>.csv  The input tables (samples x features): integer counts,
                    or hexadecimal floats (float.hex(), which R reads back
                    exactly).

  meta-<name>.csv   Their metadata: one column per variable; numbers as hex
                    floats. Column types are in tables.csv.

  tables.csv        name, types (space-separated name:type, type one of
                    'str' (sorted categories), 'cat' (a pandas Categorical,
                    levels in the order given by `levels`), 'num', 'bool'),
                    levels (name=a|b|c, space-separated, for 'cat').

  cases.csv         One row per run: case, table, method ('ancombc' or
                    'ancombc2'), formula (patsy), grouping ('' for None),
                    added (hex: added to the table before calling
                    ancombc(), as its docs ask; 0 for ancombc2),
                    pseudocount (hex: ancombc2's argument), var_quantile,
                    p_adjust (R's name), alpha, max_iter, tol (hex).

  results.csv       case, feature (1-based), covariate (patsy's name),
                    lfc, se, w, p, q (hex), signif (0/1): the primary table,
                    in its (feature, covariate) order.

  bias.csv          case, covariate (0-based), delta_em, delta_wls,
                    var_delta, n_iter: each covariate's _estimate_bias_em()
                    result, and its EM iterations (counted by wrapping
                    scipy's minimize()).

  fit.csv           case, feature, covariate (0-based), beta, var: the
                    fitted _beta_hat and _var_hat (FIT_CASES).

  vcov.csv          case, feature, i, j (0-based within the grouping's
                    columns), value: _vcov_hat (FIT_CASES).

  dof.csv           case, feature, dof (hex; ancombc2), estimable (the
                    feature's coefficients are estimable: 0/1).

  sf.csv            case, sample (0-based), theta: ANCOM-BC2's sampling
                    fractions (_sample_fractions()).

  Post-hoc tests, for the cases with a grouping (alpha and p_adjust the
  fit's):

  ph_global.csv     case, feature, w, p, q, signif: global_test().

  ph_pairs.csv      case, test ('pairwise' or 'dunnett'), feature,
                    comparison (scikit-bio's name), lfc, se, w, p, q,
                    signif: pairwise_test() and dunnett_test(bootstraps=20).

  ph_trend.csv      case, feature, w, p, q, signif: trend_test(bootstraps=20).

  ph_draws.csv      case, test, b, feature, j, value: for three cases
                    (DRAW_CASES), the null draws that dunnett_test() and
                    trend_test() used (scikit-bio's generator wrapped to
                    record them; trend: estimable features only), so that
                    ecodive can be given the same.

  ph_trend_big.csv  case, feature, p: trend_test(bootstraps=20000, seed=1)
                    on doc3, to compare ecodive's own draws in distribution.

  struc.csv         table, groups (space-separated), neg_lb, feature, then
                    one 0/1 per group: struc_zero().

Prerequisites: Python 3 with numpy, pandas, scipy, patsy and scikit-bio
0.7.4 (the version the fixtures were made with). Run from the repository
root:

    python data-raw/expected/ancombc.py
    Rscript data-raw/expected/ancombc.r
"""

import csv
import os
import warnings

import numpy as np
import pandas as pd
import scipy
import scipy.optimize
import skbio
import skbio.stats.composition._ancombc as AB
from skbio.stats.composition import ancombc, ancombc2

assert skbio.__version__ == "0.7.4", skbio.__version__

warnings.simplefilter("ignore")

OUT = os.path.join("data-raw", "expected", "ancombc")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261119)


def hx(x):
    return float(x).hex()


# ---------------------------------------------------------------------------
# Instrumentation: the EM's results and iterations, the sampling fractions.
# ---------------------------------------------------------------------------
CAPTURE = {}

_orig_em = AB._estimate_bias_em
_orig_sf = AB._sample_fractions


def _em(beta, var_hat, tol=1e-5, max_iter=100):
    # Count EM iterations: two minimize() calls per iteration.
    calls = [0]
    orig_min = AB.minimize

    def counting(*a, **k):
        calls[0] += 1
        return orig_min(*a, **k)

    AB.minimize = counting
    try:
        out = _orig_em(beta, var_hat, tol=tol, max_iter=max_iter)
    finally:
        AB.minimize = orig_min
    CAPTURE.setdefault("bias", []).append((*out, calls[0] // 2))
    return out


def _sf(*a, **k):
    out = _orig_sf(*a, **k)
    CAPTURE["sf"] = np.array(out)
    return out


AB._estimate_bias_em = _em
AB._sample_fractions = _sf


# ---------------------------------------------------------------------------
# Tables.
# ---------------------------------------------------------------------------
TABLES = {}


def add_table(name, x, meta, types, levels=None):
    TABLES[name] = (np.asarray(x), meta, types, levels or {})


# scikit-bio's ancombc() docstring examples.
x = np.array([[2, 0, 4, 7, 0, 0, 1], [1, 1, 0, 6, 5, 1, 10], [3, 0, 2, 9, 6, 0, 1],
              [0, 12, 1, 2, 0, 3, 2], [2, 2, 37, 0, 0, 7, 3], [10, 1, 0, 0, 4, 4, 3]])
add_table("doc1", x, pd.DataFrame({"status": ["healthy"] * 3 + ["sick"] * 3}),
          {"status": "str"})

x = np.array([[2, 0, 7, 0, 0, 2, 3, 2], [0, 2, 0, 1, 1, 2, 0, 0], [3, 1, 0, 9, 0, 1, 1, 0],
              [2, 0, 1, 8, 1, 0, 11, 46], [1, 0, 1, 1, 0, 0, 2, 2], [0, 3, 22, 1, 1, 1, 3, 0],
              [1, 7, 16, 1, 0, 0, 2, 2], [0, 5, 6, 1, 2, 1, 0, 1], [1, 7, 0, 2, 1, 0, 3, 2],
              [0, 6, 4, 2, 0, 2, 2, 1], [3, 13, 7, 0, 0, 0, 3, 9], [1, 8, 5, 1, 0, 0, 0, 0],
              [0, 5, 14, 1, 0, 1, 0, 1], [5, 26, 3, 2, 0, 3, 1, 3], [0, 18, 7, 0, 0, 3, 1, 0]])
add_table("doc2", x, pd.DataFrame({
    "status": ["mild"] * 5 + ["moderate"] * 5 + ["severe"] * 5,
    "age": [39., 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36]}),
    {"status": "str", "age": "num"})

# scikit-bio's ancombc2() docstring example.
x = np.array([[1, 10, 16, 0, 0, 12, 0, 7], [0, 10, 45, 2, 1, 25, 1, 9], [4, 8, 22, 1, 5, 14, 1, 5],
              [1, 11, 22, 1, 4, 14, 1, 3], [4, 7, 36, 0, 0, 19, 0, 9], [2, 30, 26, 0, 0, 12, 1, 0],
              [1, 36, 10, 1, 5, 15, 2, 9], [1, 25, 11, 2, 2, 8, 0, 5], [1, 36, 23, 0, 2, 11, 1, 2],
              [3, 60, 33, 3, 4, 14, 1, 6], [3, 57, 30, 3, 1, 3, 0, 4], [0, 65, 23, 0, 2, 2, 0, 8],
              [2, 46, 19, 0, 1, 2, 0, 1], [0, 41, 26, 0, 2, 0, 0, 3], [3, 46, 20, 0, 0, 3, 0, 0]])
add_table("doc3", x, pd.DataFrame({
    "status": ["mild"] * 5 + ["moderate"] * 5 + ["severe"] * 5,
    "age": [25., 48, 23, 35, 52, 51, 18, 29, 40, 29, 44, 39, 26, 37, 46]}),
    {"status": "str", "age": "num"})


def sim_table(n, D, groups, effects, zero_frac, depth, seed_shift=0):
    """Negative-binomial-ish counts with group effects and dropouts."""
    base = rng.normal(0, 1.5, D)
    sf = rng.normal(0, 0.6, n)
    mu = np.exp(base[None, :] + sf[:, None] + effects[groups])
    mu = mu / mu.sum(axis=1, keepdims=True) * depth[:, None]
    lam = rng.gamma(4.0, mu / 4.0)
    x = rng.poisson(lam)
    drop = rng.random((n, D)) < zero_frac
    x[drop] = 0
    return x


# A larger table: 3 sites, a numeric covariate and a binary one.
n, D = 24, 16
site = np.repeat(np.array(["gut", "oral", "skin"]), 8)
gi = np.repeat(np.arange(3), 8)
eff = np.zeros((3, D))
eff[1, :3] = 1.5
eff[2, 3:5] = -1.8
eff[2, 5:7] = 2.2
x = sim_table(n, D, gi, eff, 0.15, rng.integers(2000, 20000, n))
ph = np.round(rng.normal(6.5, 0.6, n), 2)
smoker = rng.random(n) < 0.4
add_table("sim", x, pd.DataFrame({"site": site, "ph": ph, "smoker": smoker}),
          {"site": "str", "ph": "num", "smoker": "bool"})

# Sparse: features absent from a whole group, one present in one group only,
# one observed in a single sample, and an all-but-one-zero sample.
n, D = 21, 12
trt = np.repeat(np.array(["ctrl", "dose_hi", "dose_lo"]), 7)
gi = np.repeat(np.arange(3), 7)
eff = np.zeros((3, D))
eff[1, :3] = 1.2
eff[2, 3:6] = -1.0
x = sim_table(n, D, gi, eff, 0.3, rng.integers(500, 5000, n))
x[x.sum(axis=1) == 0, 0] = 1
x[0:7, 7] = 0             # absent from the reference group (ctrl)
x[7:14, 7] = rng.integers(1, 9, 7)
x[0:7, 8] = 0             # only in dose_lo
x[7:14, 8] = 0
x[14:21, 8] = rng.integers(1, 9, 7)
x[:, 9] = 0
x[4, 9] = 7               # one sample
x[0:7, 10] = 0            # absent from ctrl and dose_hi, but one in dose_hi
x[7:14, 10] = 0
x[9, 10] = 3
x[14:21, 10] = rng.integers(1, 9, 7)
age = np.round(rng.uniform(20, 70, n))
add_table("sparse", x, pd.DataFrame({"trt": trt, "age": age}),
          {"trt": "str", "age": "num"})

# A factor with its own level order (a pandas Categorical), an interaction,
# and a covariate on a large scale (an ill-conditioned design).
n, D = 18, 10
dose = np.tile(np.array(["low", "mid", "high"]), 6)
sex = np.repeat(np.array(["f", "m"]), 9)
gi = np.tile(np.arange(3), 6)
eff = np.zeros((3, D))
eff[2, :3] = 1.4
eff[1, 3:5] = -1.1
x = sim_table(n, D, gi, eff, 0.2, rng.integers(1000, 8000, n))
days = np.round(rng.uniform(5000, 30000, n))
add_table("levels", x, pd.DataFrame({
    "dose": pd.Categorical(dose, categories=["low", "mid", "high"]),
    "sex": sex, "days": days}),
    {"dose": "cat", "sex": "str", "days": "num"}, {"dose": ["low", "mid", "high"]})

# Proportions with no zeros (a replace_zeros()-like table).
n, D = 12, 8
grp = np.repeat(np.array(["a", "b"]), 6)
gi = np.repeat(np.arange(2), 6)
eff = np.zeros((2, D))
eff[1, :2] = 1.0
x = sim_table(n, D, gi, eff, 0.0, rng.integers(3000, 6000, n)) + 0.5
x = x / x.sum(axis=1, keepdims=True)
add_table("prop", x, pd.DataFrame({"grp": grp}), {"grp": "str"})


# ---------------------------------------------------------------------------
# Cases.
# ---------------------------------------------------------------------------
CASES = []


def case(name, table, method, formula, grouping=None, added=0.0, pseudocount=0.0,
         var_quantile=0.05, p_adjust="holm", alpha=0.05, max_iter=100, tol=1e-5):
    CASES.append(dict(case=name, table=table, method=method, formula=formula,
                      grouping=grouping, added=added, pseudocount=pseudocount,
                      var_quantile=var_quantile, p_adjust=p_adjust, alpha=alpha,
                      max_iter=max_iter, tol=tol))


case("doc1", "doc1", "ancombc", "status", added=1.0)
case("doc2", "doc2", "ancombc", "status + age", "status", added=1.0)
case("doc3", "doc3", "ancombc2", "status + age", "status")
for pc in (0.5, 1.0):
    case(f"doc3_pc{pc}", "doc3", "ancombc2", "status + age", "status", pseudocount=pc)
case("doc3_bc", "doc3", "ancombc", "status + age", "status", added=0.5)

case("sim_bc", "sim", "ancombc", "site + ph + smoker", "site", added=0.5,
     p_adjust="BH", alpha=0.1)
case("sim_bc_it", "sim", "ancombc", "site + ph", "site", added=0.5, max_iter=3,
     tol=1e-3, p_adjust="bonferroni")
case("sim_bc2", "sim", "ancombc2", "site + ph + smoker", "site")
case("sim_bc2_pc", "sim", "ancombc2", "site + ph", "site", pseudocount=1.0,
     var_quantile=0.2, p_adjust="BY")
case("sim_bc2_vq0", "sim", "ancombc2", "site + ph", None, var_quantile=0.0,
     p_adjust="none")
case("sim_bc2_ia", "sim", "ancombc2", "site * smoker", "site")
case("sim_bc_ia", "sim", "ancombc", "site * ph", "site", added=0.5)
case("sim_bc2_noint", "sim", "ancombc2", "0 + site + ph", "site")

case("sparse_bc2", "sparse", "ancombc2", "trt + age", "trt")
case("sparse_bc2_one", "sparse", "ancombc2", "trt", "trt", var_quantile=0.1)
case("sparse_bc", "sparse", "ancombc", "trt + age", "trt", added=1.0)

case("levels_bc2", "levels", "ancombc2", "dose + sex + days", "dose")
case("levels_bc", "levels", "ancombc", "dose + sex + days", "dose", added=0.5)
case("levels_bc2_ia", "levels", "ancombc2", "dose * sex", "dose")

case("prop_bc", "prop", "ancombc", "grp")
case("prop_bc2", "prop", "ancombc2", "grp")


# ---------------------------------------------------------------------------
# Run and write.
# ---------------------------------------------------------------------------
SKBIO_ADJ = {"holm": "holm", "BH": "bh", "bonferroni": "bonferroni",
             "BY": "by", "none": None}


def write_table(name, x, meta, types, levels):
    with open(os.path.join(OUT, f"table-{name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        integer = np.issubdtype(x.dtype, np.integer)
        for row in x:
            w.writerow([str(int(v)) if integer else hx(v) for v in row])
    with open(os.path.join(OUT, f"meta-{name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(list(meta.columns))
        for i in range(meta.shape[0]):
            row = []
            for c in meta.columns:
                v = meta[c].iloc[i]
                t = types[c]
                if t == "num":
                    row.append(hx(v))
                elif t == "bool":
                    row.append("TRUE" if v else "FALSE")
                else:
                    row.append(str(v))
            w.writerow(row)


with open(os.path.join(OUT, "tables.csv"), "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["name", "types", "levels"])
    for name, (x, meta, types, levels) in TABLES.items():
        write_table(name, x, meta, types, levels)
        w.writerow([name, " ".join(f"{k}:{v}" for k, v in types.items()),
                    " ".join(f"{k}={'|'.join(v)}" for k, v in levels.items())])

f_cases = open(os.path.join(OUT, "cases.csv"), "w", newline="")
f_res = open(os.path.join(OUT, "results.csv"), "w", newline="")
f_bias = open(os.path.join(OUT, "bias.csv"), "w", newline="")
f_fit = open(os.path.join(OUT, "fit.csv"), "w", newline="")
f_vcov = open(os.path.join(OUT, "vcov.csv"), "w", newline="")
f_dof = open(os.path.join(OUT, "dof.csv"), "w", newline="")
f_sf = open(os.path.join(OUT, "sf.csv"), "w", newline="")
f_pg = open(os.path.join(OUT, "ph_global.csv"), "w", newline="")
f_pp = open(os.path.join(OUT, "ph_pairs.csv"), "w", newline="")
f_pt = open(os.path.join(OUT, "ph_trend.csv"), "w", newline="")
f_pd = open(os.path.join(OUT, "ph_draws.csv"), "w", newline="")
W = {k: csv.writer(v, lineterminator="\n") for k, v in dict(cases=f_cases, res=f_res, bias=f_bias, fit=f_fit,
                                       vcov=f_vcov, dof=f_dof, sf=f_sf, pg=f_pg, pp=f_pp,
                                       pt=f_pt, pd=f_pd).items()}
W["pg"].writerow(["case", "feature", "w", "p", "q", "signif"])
W["pp"].writerow(["case", "test", "feature", "comparison", "lfc", "se", "w", "p", "q", "signif"])
W["pt"].writerow(["case", "feature", "w", "p", "q", "signif"])
W["pd"].writerow(["case", "test", "b", "feature", "j", "value"])


# Cases whose fitted _beta_hat, _var_hat and _vcov_hat are written (the
# rest are covered by the primary table).
FIT_CASES = {"doc2", "doc3", "sparse_bc2", "sim_bc2_noint"}

# Cases whose null draws are recorded (the others' bootstrapped results are
# compared only in their deterministic parts).
DRAW_CASES = {"doc2", "doc3", "sparse_bc2"}


class Recorder:
    """A NumPy generator that records the draws it hands out."""

    def __init__(self, seed):
        self.rng = np.random.default_rng(seed)
        self.draws = []

    def standard_t(self, df, size):
        x = self.rng.standard_t(df, size=size)
        self.draws.append(np.array(x))
        return x

    def standard_normal(self, size):
        x = self.rng.standard_normal(size=size)
        self.draws.append(np.array(x))
        return x


def recorded(fun, seed):
    rec = Recorder(seed)
    orig = AB.get_rng
    AB.get_rng = lambda s: rec
    try:
        out = fun()
    finally:
        AB.get_rng = orig
    return out, rec.draws
W["cases"].writerow(["case", "table", "method", "formula", "grouping", "added",
                     "pseudocount", "var_quantile", "p_adjust", "alpha", "max_iter", "tol"])
W["res"].writerow(["case", "feature", "covariate", "lfc", "se", "w", "p", "q", "signif"])
W["bias"].writerow(["case", "covariate", "delta_em", "delta_wls", "var_delta", "n_iter"])
W["fit"].writerow(["case", "feature", "covariate", "beta", "var"])
W["vcov"].writerow(["case", "feature", "i", "j", "value"])
W["dof"].writerow(["case", "feature", "dof", "estimable"])
W["sf"].writerow(["case", "sample", "theta"])

for c in CASES:
    x, meta, types, levels = TABLES[c["table"]]
    CAPTURE.clear()
    xin = x + c["added"] if c["added"] else x
    if c["method"] == "ancombc":
        res = ancombc(xin, meta, c["formula"], grouping=c["grouping"],
                      max_iter=c["max_iter"], tol=c["tol"], alpha=c["alpha"],
                      p_adjust=SKBIO_ADJ[c["p_adjust"]])
    else:
        res = ancombc2(x, meta, c["formula"], grouping=c["grouping"],
                       pseudocount=c["pseudocount"], var_quantile=c["var_quantile"],
                       max_iter=c["max_iter"], tol=c["tol"], alpha=c["alpha"],
                       p_adjust=SKBIO_ADJ[c["p_adjust"]])

    W["cases"].writerow([c["case"], c["table"], c["method"], c["formula"],
                         c["grouping"] or "", hx(c["added"]), hx(c["pseudocount"]),
                         hx(c["var_quantile"]), c["p_adjust"], hx(c["alpha"]),
                         c["max_iter"], hx(c["tol"])])

    tab = res.result
    feats = list(tab.index.get_level_values(0))
    fidx = {f: i for i, f in enumerate(dict.fromkeys(feats))}
    for (f, cov), row in tab.iterrows():
        W["res"].writerow([c["case"], fidx[f] + 1, cov, hx(row["Log(FC)"]), hx(row["SE"]),
                           hx(row["W"]), hx(row["pvalue"]), hx(row["qvalue"]),
                           int(bool(row["Signif"]))])

    for i, b in enumerate(CAPTURE["bias"]):
        W["bias"].writerow([c["case"], i, hx(b[0]), hx(b[1]), hx(b[2]), b[3]])

    if "sf" in CAPTURE:
        for s, t in enumerate(CAPTURE["sf"]):
            W["sf"].writerow([c["case"], s, hx(t)])

    if c["grouping"]:
        beta, var = res._beta_hat, res._var_hat
        for f in range(beta.shape[0] if c["case"] in FIT_CASES else 0):
            for p in range(beta.shape[1]):
                W["fit"].writerow([c["case"], f, p, hx(beta[f, p]), hx(var[f, p])])
        vc = res._vcov_hat
        for f in range(vc.shape[0] if c["case"] in FIT_CASES else 0):
            for i in range(vc.shape[1]):
                for j in range(vc.shape[2]):
                    W["vcov"].writerow([c["case"], f, i, j, hx(vc[f, i, j])])
        dof = res._dof
        est = res._estimable
        for f in range(beta.shape[0]):
            d = np.nan if dof is None else (dof if np.isscalar(dof) else dof[f])
            e = 1 if est is None else int(np.all(est[f]))
            W["dof"].writerow([c["case"], f, hx(d), e])

        # Post-hoc tests.
        g = res.global_test()
        for f, (_, row) in enumerate(g.iterrows()):
            W["pg"].writerow([c["case"], f, hx(row["W"]), hx(row["pvalue"]),
                              hx(row["qvalue"]), int(bool(row["Signif"]))])

        def write_pairs(test, tab):
            fidx = {fid: i for i, fid in enumerate(dict.fromkeys(tab.index.get_level_values(0)))}
            for (fid, comp), row in tab.iterrows():
                W["pp"].writerow([c["case"], test, fidx[fid], comp, hx(row["Log(FC)"]),
                                  hx(row["SE"]), hx(row["W"]), hx(row["pvalue"]),
                                  hx(row["qvalue"]), int(bool(row["Signif"]))])

        write_pairs("pairwise", res.pairwise_test())

        keep = c["case"] in DRAW_CASES
        dn, draws = recorded(lambda: res.dunnett_test(bootstraps=20, seed=0), 11)
        write_pairs("dunnett", dn)
        for b, d in enumerate(draws if keep else []):
            for f in range(d.shape[0]):
                for j in range(d.shape[1]):
                    W["pd"].writerow([c["case"], "dunnett", b, f, j, hx(d[f, j])])

        tr, draws = recorded(lambda: res.trend_test(bootstraps=20, seed=0), 12)
        for f, (_, row) in enumerate(tr.iterrows()):
            W["pt"].writerow([c["case"], f, hx(row["W"]), hx(row["pvalue"]),
                              hx(row["qvalue"]), int(bool(row["Signif"]))])
        valid = res._posthoc_estimable()
        work = np.arange(beta.shape[0]) if valid is None else np.flatnonzero(valid)
        for b, d in enumerate(draws if keep else []):
            for r in range(d.shape[0]):
                for j in range(d.shape[1]):
                    W["pd"].writerow([c["case"], "trend", b, work[r], j, hx(d[r, j])])

        if c["case"] == "doc3":
            big = res.trend_test(bootstraps=20000, seed=1)
            with open(os.path.join(OUT, "ph_trend_big.csv"), "w", newline="") as fh:
                wb = csv.writer(fh, lineterminator="\n")
                wb.writerow(["case", "feature", "p"])
                for f, (_, row) in enumerate(big.iterrows()):
                    wb.writerow([c["case"], f, hx(row["pvalue"])])

for fh in (f_cases, f_res, f_bias, f_fit, f_vcov, f_dof, f_sf, f_pg, f_pp, f_pt, f_pd):
    fh.close()


# ---------------------------------------------------------------------------
# struc_zero(): its docstring example, and the sparse tables.
# ---------------------------------------------------------------------------
from skbio.stats.composition import struc_zero

x = np.array([[7, 1, 0, 11, 3, 1], [1, 1, 0, 13, 13, 0], [11, 5, 0, 1, 4, 1],
              [2, 2, 0, 16, 4, 0], [0, 1, 0, 0, 6, 0], [14, 8, 7, 9, 0, 5],
              [0, 7, 4, 1, 0, 26], [8, 1, 4, 28, 0, 10], [2, 2, 2, 4, 0, 5],
              [6, 4, 10, 1, 0, 9]])
add_table("struc", x, pd.DataFrame({"grouping": ["treatment"] * 5 + ["placebo"] * 5}),
          {"grouping": "str"})
with open(os.path.join(OUT, "tables.csv"), "a", newline="") as fh:
    csv.writer(fh, lineterminator="\n").writerow(["struc", "grouping:str", ""])
write_table("struc", x, TABLES["struc"][1], TABLES["struc"][2], {})

with open(os.path.join(OUT, "struc.csv"), "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["table", "grouping", "neg_lb", "groups", "feature", "zero"])
    for name, grouping in (("struc", "grouping"), ("doc2", "status"), ("doc3", "status"),
                           ("sparse", "trt"), ("levels", "dose"), ("sim", "site")):
        xx, meta = TABLES[name][0], TABLES[name][1]
        for neg_lb in (False, True):
            z = struc_zero(xx, meta, grouping, neg_lb=neg_lb)
            for f in range(z.shape[0]):
                w.writerow([name, grouping, int(neg_lb), " ".join(map(str, z.columns)), f,
                            " ".join(str(int(v)) for v in z.iloc[f].values)])

print("scikit-bio", skbio.__version__, "numpy", np.__version__, "scipy",
      scipy.__version__, "pandas", pd.__version__)
print(len(CASES), "cases")
