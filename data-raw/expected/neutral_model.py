"""Reference fixtures for ecodive's neutral_model(): neufit's outputs.

neufit (github.com/misieber/neufit, neufit.py, Michael Sieber 2018, BSD-2)
fits Sloan et al.'s (2006) neutral model with lmfit's least squares (m
bounded to [0, 1], N fixed at the read depth) and draws statsmodels'
Wilson intervals around the fit. This script runs neufit.py itself (with
runpy, as `python neufit.py <table>` would) on each table and writes what
it computed, so the fixtures are neufit's own outputs.

Tables (written here, under data-raw/expected/neutral_model/, as
features x samples tab-separated integer counts, neufit's input format):

  table-sim.tsv     neufit's own example, sim_data.csv (100 OTUs, 50
                    samples of 1,000 reads), copied from the checkout.
  table-dm.tsv      A neutral community drawn by Sloan's process (each
                    sample a Dirichlet-multinomial with alpha = N m p),
                    150 features, 40 samples of 2,000 reads, m = 0.02.
  table-uneven.tsv  The same process with unequal depths (500 to 1,500
                    reads): neufit rarefies it to the smallest depth with
                    NumPy's global RNG (seeded here). The rarefied table it
                    fitted is written as table-uneven_rarefied.tsv.
  table-sat.tsv     A near-saturated community (m = 1, 80 features, 25
                    samples of 1,000 reads, drawn with its own generator),
                    where the least-squares m is above 1 and neufit's
                    bound at 1 is active.
  table-pool.tsv    A source metacommunity for table-dm.tsv: 20 more
                    samples of 2,000 reads drawn from dm's metacommunity
                    with m = 0.5 (for neutral_model()'s `pool`; neufit has
                    no pool, so only the high-precision reference and the
                    R references use it).

Outputs (data-raw/expected/neutral_model/neufit-<table>.csv): one row per
feature in neufit's order (by mean abundance): feature, mean_abundance,
occurrence, predicted, lower, upper (hex floats); and neufit-stats.csv:
table, n_reads, n_samples, m, m_stderr, r_squared (neufit's own fit), and
m_tight (hex floats): neufit's model and lmfit Model refitted from neufit's
m with ftol = xtol = gtol = 1e-15, to show how much of the difference
from the exact least-squares m is lmfit's stopping rule.

hp-stats.csv: a high-precision reference (mpmath, 40 digits) for each
table (and dm with pool): the zero of the gradient of the residual sum of
squares in m (the least-squares m, unbounded), and the ends of its profile
interval, where (SS(m) - SS(m_hat)) / (SS(m_hat) / (K - 1)) equals the
0.95 quantile of F(1, K - 1) (that quantile from SciPy): table, K, m,
m_lower, m_upper (hex floats of the 40-digit values).

wilson.csv: statsmodels' proportion_confint(method='wilson') at the edges
(count 0 and n) and nearby, for conf levels 0.95 and 0.9.

Prerequisites: Python 3 with numpy, pandas, scipy, lmfit 1.3.4,
statsmodels 0.15.0, mpmath 1.4.1 and matplotlib, and a copy of neufit (the commit whose
neufit.py has the SHA-256 below) in the directory named by NEUFIT_DIR. Run
from the repository root:

    NEUFIT_DIR=/path/to/neufit python data-raw/expected/neutral_model.py
    Rscript data-raw/expected/neutral_model.r
"""

import csv
import hashlib
import io
import os
import runpy
import shutil
import sys
import tempfile
import contextlib

import matplotlib
matplotlib.use("Agg")

import lmfit
import mpmath as mp
import numpy as np
import scipy.stats
import pandas as pd
import statsmodels
from statsmodels.stats.proportion import proportion_confint

assert lmfit.__version__ == "1.3.4", lmfit.__version__
assert statsmodels.__version__ == "0.15.0", statsmodels.__version__
assert mp.__version__ == "1.4.1", mp.__version__

NEUFIT_DIR = os.environ["NEUFIT_DIR"]
NEUFIT = os.path.join(NEUFIT_DIR, "neufit.py")
SHA = "15c4bfefd568ddc72160db0a29b3d40c5b0179200b705d28d1d1d4e143a151b0"
with open(NEUFIT, "rb") as fh:
    sha = hashlib.sha256(fh.read()).hexdigest()
assert sha == SHA, sha

OUT = os.path.join("data-raw", "expected", "neutral_model")
os.makedirs(OUT, exist_ok=True)


def hx(v):
    return float(v).hex()


def write_table(name, counts, features, samples):
    """counts: features x samples."""
    df = pd.DataFrame(counts, index=features, columns=samples)
    df.index.name = ""
    path = os.path.join(OUT, f"table-{name}.tsv")
    df.to_csv(path, sep="\t")
    return path


def neutral_draw(rng, D, n, depths, m, shape, p=None):
    """Sloan's process: metacommunity p, each sample Dirichlet(N m p),
    then a multinomial draw of its depth. Returns (counts, p) when p is
    given, else the counts without all-zero features."""
    given = p is not None
    if not given:
        p = rng.gamma(shape, 1.0, D)
        p = p / p.sum()
    out = np.zeros((D, n), dtype=np.int64)
    for j in range(n):
        N = depths[j]
        x = rng.dirichlet(np.maximum(N * m * p, 1e-300))
        out[:, j] = rng.multinomial(N, x / x.sum())
    if given:
        return out
    keep = out.sum(1) > 0
    return out[keep]


# ---------------------------------------------------------------------------
# Tables.
rng = np.random.default_rng(20261124)

shutil.copyfile(os.path.join(NEUFIT_DIR, "sim_data.csv"), os.path.join(OUT, "table-sim.tsv"))

p_dm = rng.gamma(0.3, 1.0, 150)
p_dm = p_dm / p_dm.sum()
dm_all = neutral_draw(rng, 150, 40, [2000] * 40, 0.02, None, p_dm)
pool = neutral_draw(rng, 150, 20, [2000] * 20, 0.5, None, p_dm)
dm_keep = dm_all.sum(1) > 0
dm = dm_all[dm_keep]
dm_names = [f"F{i + 1}" for i in range(dm.shape[0])]
write_table("dm", dm, dm_names, [f"S{j + 1}" for j in range(dm.shape[1])])
# The pool has every feature of dm, and the features dm lacks, named P<i>.
pool_names = []
k = 0
for i in range(150):
    if dm_keep[i]:
        pool_names.append(dm_names[k])
        k += 1
    else:
        pool_names.append(f"P{i + 1}")
write_table("pool", pool, pool_names, [f"Q{j + 1}" for j in range(pool.shape[1])])

depths = rng.integers(500, 1501, 30)
un = neutral_draw(rng, 120, 30, depths, 0.05, 0.3)
write_table("uneven", un, [f"F{i + 1}" for i in range(un.shape[0])],
            [f"S{j + 1}" for j in range(un.shape[1])])

sat = neutral_draw(np.random.default_rng(20261124), 80, 25, [1000] * 25, 1.0, 0.2)
write_table("sat", sat, [f"F{i + 1}" for i in range(sat.shape[0])],
            [f"S{j + 1}" for j in range(sat.shape[1])])


# ---------------------------------------------------------------------------
# neufit on each table.
def run_neufit(name, seed=None):
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, f"{name}.tsv")
    shutil.copyfile(os.path.join(OUT, f"table-{name}.tsv"), path)
    if seed is not None:
        np.random.seed(seed)
    argv = sys.argv
    sys.argv = ["neufit.py", path]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            g = runpy.run_path(NEUFIT, run_name="__main__")
    finally:
        sys.argv = argv
    shutil.rmtree(tmp)
    return g


stats_rows = []
for name, seed in [("sim", None), ("dm", None), ("uneven", 20261124), ("sat", None)]:
    g = run_neufit(name, seed)
    of = g["occurr_freqs"]
    fit = g["beta_fit"]
    with open(os.path.join(OUT, f"neufit-{name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["feature", "mean_abundance", "occurrence", "predicted", "lower", "upper"])
        for feat, r in of.iterrows():
            w.writerow([feat, hx(r["mean_abundance"]), hx(r["occurrence"]),
                        hx(r["predicted_occurrence"]), hx(r["lower_conf_int"]),
                        hx(r["upper_conf_int"])])
    se = fit.params["m"].stderr
    # neufit's model, refitted from its m with tight tolerances.
    params = fit.params.copy()
    tight = g["beta_model"].fit(
        of["occurrence"], params, p=of["mean_abundance"],
        max_nfev=100000, fit_kws={"ftol": 1e-15, "xtol": 1e-15, "gtol": 1e-15})
    stats_rows.append([name, g["n_reads"], g["n_samples"], hx(fit.best_values["m"]),
                       hx(se) if se is not None else "NA", hx(g["r_square"]),
                       hx(tight.best_values["m"])])
    if name == "uneven":
        ab = g["abundances"]
        ab.index.name = ""
        ab.to_csv(os.path.join(OUT, "table-uneven_rarefied.tsv"), sep="\t")

with open(os.path.join(OUT, "neufit-stats.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["table", "n_reads", "n_samples", "m", "m_stderr", "r_squared", "m_tight"])
    w.writerows(stats_rows)


# ---------------------------------------------------------------------------
# High-precision reference (mpmath): the least-squares m and its profile
# interval, from exact p and occurrence (ratios of integers).
mp.mp.dps = 40


def read_counts(name):
    return pd.read_csv(os.path.join(OUT, f"table-{name}.tsv"), sep="\t", index_col=0)


def hp_fit(totals, occs, N, n, conf=0.95):
    # totals: the counts' (or the pool's) feature totals and their sum;
    # occs: the number of samples each feature is in.
    tot, grand = totals
    keep = [i for i in range(len(tot)) if tot[i] > 0 and occs[i] > 0]
    p = [mp.mpf(int(tot[i])) / mp.mpf(int(grand)) for i in keep]
    f = [mp.mpf(int(occs[i])) / n for i in keep]
    K = len(keep)
    Nm = mp.mpf(N)
    d = 1 / Nm

    def F(pi, m):
        return 1 - mp.betainc(Nm * m * pi, Nm * m * (1 - pi), 0, d, regularized=True)

    def rss(m):
        return mp.fsum((fi - F(pi, m)) ** 2 for pi, fi in zip(p, f))

    def grad(m):
        return mp.fsum((fi - F(pi, m)) * mp.diff(lambda t: F(pi, t), m) for pi, fi in zip(p, f))

    return p, f, K, rss, grad


def hp_solve(name, totals, occs, N, n, m0):
    p, f, K, rss, grad = hp_fit(totals, occs, N, n)
    m = mp.findroot(grad, (mp.mpf(m0) * (1 - mp.mpf("1e-4")), mp.mpf(m0) * (1 + mp.mpf("1e-4"))),
                    solver="secant", tol=mp.mpf(10) ** -60)
    r0 = rss(m)
    s2 = r0 / (K - 1)
    cut = mp.mpf(scipy.stats.f.ppf(0.95, 1, K - 1))
    tau = lambda t: (rss(t) - r0) / s2 - cut
    # Start from the delta-method half-width, then the secant method.
    J = mp.fsum(mp.diff(lambda t: 1 - mp.betainc(N * t * pi, N * t * (1 - pi), 0, mp.mpf(1) / N,
                                                 regularized=True), m) ** 2 for pi in p)
    hw = mp.sqrt(cut * s2 / J)
    lo = mp.findroot(tau, (m - hw, m - hw * mp.mpf("0.9")), solver="secant", tol=mp.mpf(10) ** -50)
    hi = mp.findroot(tau, (m + hw, m + hw * mp.mpf("1.1")), solver="secant", tol=mp.mpf(10) ** -50)
    return [name, K, hx(float(m)), hx(float(lo)), hx(float(hi)), hx(float(cut))]


m_start = {r[0]: float.fromhex(r[6]) for r in stats_rows}
hp_rows = []
for name in ["sim", "dm", "uneven_rarefied", "sat"]:
    tab = read_counts(name).to_numpy()
    n = tab.shape[1]
    N = tab.sum(0).mean()
    assert all(tab.sum(0) == N)
    tot = tab.sum(1)
    m0 = m_start["uneven" if name == "uneven_rarefied" else name]
    if name == "sat":
        m0 = 1.2   # neufit stops at its bound, 1
    hp_rows.append(hp_solve(name, (tot, tot.sum()), (tab > 0).sum(1), int(N), n, m0))

# dm with the pool: p from the pool's totals (all its features).
dm_tab = read_counts("dm")
pool_tab = read_counts("pool")
ptot = pool_tab.sum(1)
tot = np.array([ptot[f] for f in dm_tab.index])
hp_rows.append(hp_solve("dm_pool", (tot, ptot.sum()), (dm_tab.to_numpy() > 0).sum(1), 2000,
                        dm_tab.shape[1], m_start["dm"]))

with open(os.path.join(OUT, "hp-stats.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["table", "K", "m", "m_lower", "m_upper", "cutoff"])
    w.writerows(hp_rows)


# ---------------------------------------------------------------------------
# statsmodels' Wilson interval at and near the edges.
with open(os.path.join(OUT, "wilson.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["n", "conf_level", "count", "lower", "upper"])
    for n in [10, 50]:
        for cl in [0.95, 0.9]:
            for x in [0.0, 1e-12, 0.5, n / 3, n - 0.5, n - 1e-12, float(n)]:
                lo, hi = proportion_confint(x, n, alpha=1 - cl, method="wilson")
                w.writerow([n, cl, hx(x), hx(lo), hx(hi)])
