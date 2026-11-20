"""Reference fixtures for ecodive's sourcetracker().

SourceTracker2's Gibbs sampler draws from NumPy's legacy MT19937, which
ecodive can't reproduce, so its deterministic parts are compared exactly
and the sampler's results by distribution:

  * collapse_source_data(..., 'mean'): the pooled environments.
  * ConditionalProbability: the conditional probabilities of a read for
    given states of the sampler.
  * gibbs_sampler()'s schedule: which passes are followed by a draw. A
    ConditionalProbability that sends the one read of a sink to state p on
    pass p shows, in each draw's counts, the pass it was taken after.
  * gibbs(): its proportions over many seeds (np.random.seed(seed), one
    process), on simulated sinks that are known mixtures of deep sources
    and an unknown source, already rarefied (gibbs() doesn't rarefy).
  * Leave-one-out: gibbs(sources, None) fails. It picks the other sources
    with `sources[lambda x: x != source]`, which masks values instead of
    dropping the held-out row, so the sampler gets one source too many and
    collate_gibbs_results() can't fit its results. The loo runs therefore
    do what gibbs() means to: for each environment in turn,
    gibbs_sampler() with ConditionalProbability on the other environments
    (`sources.drop(index=source)`), then collate_gibbs_results(loo=True).

Writes, under data-raw/expected/sourcetracker/ (then
data-raw/expected/sourcetracker.r converts them to
tests/tinytests/expected/sourcetracker.rds):

  collapse-table.csv   Samples x features, integer counts.
  collapse-env.csv     One environment label per sample (sources only).
  collapse-result.csv  The collapsed table: env, then the counts.

  cp-source.csv        The E x F source table for the cases below.
  cp.csv               One row per case: alpha1, alpha2, beta, n, xi
                       (0-based), m_xiV, m_V, n_vnoti and the result (both
                       space-separated; the result as hex floats).

  schedule.csv         burnin, draws_per_restart, delay, then the passes
                       (1-based, space-separated) followed by a draw.

  sim-sources.csv      Environments x features: the (rarefied) sources.
  sim-sinks.csv        Sinks x features: the (rarefied) sinks.
  sim-truth.csv        Each sink's true mixing proportions (environments,
                       then Unknown).
  sim-<run>.csv        One row per seed and sink: run, seed, sink, then the
                       proportions and standard deviations (mpm and mps) as
                       hex floats, space-separated.
  sim-<run>-fas.csv    The feature tables of seed 1: sink, source, then the
                       counts (sums over the draws, as gibbs() returns them).

Prerequisites: Python 3 with numpy 2.5.5, pandas 2.3.3 (SourceTracker2 uses
DataFrame.applymap(), which pandas 3 removed) and scikit-bio 0.7.4, and
SourceTracker2's sources (github.com/biota/sourcetracker2, master,
2.0.1-dev; BSD-3), of which only sourcetracker/_sourcetracker.py is used.
Run from the repository root, giving the directory that holds the
`sourcetracker` package:

    python -I data-raw/expected/sourcetracker.py <sourcetracker2 dir>
    Rscript data-raw/expected/sourcetracker.r

The sampled runs take about half an hour on 4 cores.
"""

import csv
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(sys.argv[1]))
from sourcetracker._sourcetracker import (  # noqa: E402
    ConditionalProbability, collapse_source_data, collate_gibbs_results,
    gibbs, gibbs_sampler)

warnings.simplefilter("ignore")   # pandas' applymap() deprecation

OUT = os.path.abspath(os.path.join("data-raw", "expected", "sourcetracker"))
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20261120)


def hx(v):
    return " ".join(float(x).hex() for x in np.asarray(v, dtype=float).ravel())


def ints(v):
    return " ".join(str(int(x)) for x in np.asarray(v).ravel())


# ---------------------------------------------------------------------------
# collapse_source_data(): means that are not whole numbers, environment
# names that sort differently by locale ('C' < 'a' < 'b' in C locale), and
# a feature that is zero everywhere.
# ---------------------------------------------------------------------------

envs = ["b", "a", "C", "b", "a", "C", "b", "C", "a", "b"]
tab = rng.poisson(rng.lognormal(1.5, 1.2, size=(1, 9)), size=(10, 9))
tab[:, 4] = 0
tab[1, 2] = 7     # 'a' samples: a mean that truncates
ids = [f"s{i}" for i in range(10)]
ftab = pd.DataFrame(tab, index=ids, columns=[f"f{j}" for j in range(9)])
meta = pd.DataFrame({"Env": envs}, index=ids)
col = collapse_source_data(meta, ftab, ids, "Env", "mean")

np.savetxt(os.path.join(OUT, "collapse-table.csv"), tab, fmt="%d", delimiter=",")
with open(os.path.join(OUT, "collapse-env.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["env"])
    for e in envs:
        w.writerow([e])
with open(os.path.join(OUT, "collapse-result.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["env"] + list(col.columns))
    for e, row in col.iterrows():
        w.writerow([e] + [int(x) for x in row.values])


# ---------------------------------------------------------------------------
# ConditionalProbability.calculate_cp_slice() on given states.
# ---------------------------------------------------------------------------

src = rng.poisson(rng.lognormal(2, 1.5, size=(3, 1)) * rng.dirichlet(np.ones(12) * 0.5, size=3) * 40)
src[1, 5] = 0
np.savetxt(os.path.join(OUT, "cp-source.csv"), src, fmt="%d", delimiter=",")

with open(os.path.join(OUT, "cp.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["alpha1", "alpha2", "beta", "n", "xi", "m_xiV", "m_V", "n_vnoti", "result"])
    for a1, a2, b in [(0.001, 0.1, 10), (0.5, 0.001, 10), (0.01, 1, 0.5), (0, 0.1, 10)]:
        for _ in range(6):
            n = int(rng.integers(20, 2000))
            xi = int(rng.integers(0, 12))
            m_V = int(rng.integers(0, n))
            m_xiV = int(rng.integers(0, m_V + 1))
            c = rng.multinomial(n - 1, rng.dirichlet(np.ones(4)))
            cp = ConditionalProbability(a1, a2, b, src)
            cp.set_n(n)
            cp.precalculate()
            jp = cp.calculate_cp_slice(xi, m_xiV, m_V, c)
            w.writerow([repr(a1), repr(a2), repr(b), n, xi, m_xiV, m_V, ints(c), hx(jp)])


# ---------------------------------------------------------------------------
# gibbs_sampler()'s schedule. One read; on pass p the read can only go to
# state p (a one-hot probability vector), so each draw's counts say which
# pass it followed. np.random.shuffle() counts the passes.
# ---------------------------------------------------------------------------

class PassProbability:

    def __init__(self, n_states):
        self.V = n_states
        self.tau = 1
        self.pass_no = 0

    def set_n(self, n):
        pass

    def precalculate(self):
        pass

    def calculate_cp_slice(self, xi, m_xiV, m_V, n_vnoti):
        jp = np.zeros(self.V)
        jp[self.pass_no] = 1.0
        return jp


with open(os.path.join(OUT, "schedule.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["burnin", "draws_per_restart", "delay", "passes"])
    real_shuffle = np.random.shuffle
    for burnin, draws, delay in [(1, 1, 1), (100, 1, 1), (5, 3, 1), (5, 3, 4),
                                 (2, 4, 3), (10, 2, 10), (3, 5, 2)]:
        total = burnin + (draws - 1) * delay + 1
        cp = PassProbability(total + 1)

        def shuffle(x, cp=cp):
            cp.pass_no += 1
            real_shuffle(x)

        np.random.shuffle = shuffle
        try:
            np.random.seed(1)
            envcounts, _, _ = gibbs_sampler(np.array([1]), cp, 1, draws, burnin, delay)
        finally:
            np.random.shuffle = real_shuffle
        passes = [int(np.argmax(row)) for row in envcounts]
        w.writerow([burnin, draws, delay, ints(passes)])


# ---------------------------------------------------------------------------
# Simulated mixtures: 3 deep sources, an unknown source, 6 sinks.
# ---------------------------------------------------------------------------

F = 40
E = 3
dists = rng.dirichlet(np.ones(F) * 0.3, size=E + 1)    # the last is unknown
sources = np.vstack([rng.multinomial(2000, dists[e]) for e in range(E)])
truth = np.array([
    [0.5, 0.3, 0.2, 0.0],
    [0.7, 0.0, 0.0, 0.3],
    [0.2, 0.2, 0.2, 0.4],
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
    [0.1, 0.6, 0.3, 0.0]])
sinks = np.vstack([rng.multinomial(400, truth[i] @ dists) for i in range(len(truth))])

env_ids = ["soil", "water", "skin"]
sink_ids = [f"sink{i + 1}" for i in range(len(truth))]
feat_ids = [f"f{j}" for j in range(F)]
np.savetxt(os.path.join(OUT, "sim-sources.csv"), sources, fmt="%d", delimiter=",")
np.savetxt(os.path.join(OUT, "sim-sinks.csv"), sinks, fmt="%d", delimiter=",")
np.savetxt(os.path.join(OUT, "sim-truth.csv"), truth, delimiter=",")

src_df = pd.DataFrame(sources, index=env_ids, columns=feat_ids)
snk_df = pd.DataFrame(sinks, index=sink_ids, columns=feat_ids)

# run: (restarts, draws_per_restart, burnin, delay, loo, number of seeds)
RUNS = {
    "default": (10, 1, 100, 1, False, 40),
    "thinned": (4, 3, 20, 5, False, 40),
    "loo":     (10, 1, 100, 1, True, 40),
}


def one_run(job):
    run, seed = job
    restarts, draws, burnin, delay, loo, _ = RUNS[run]
    np.random.seed(seed)
    old = os.getcwd()
    os.chdir(TMP)   # cumulative_proportions() np.save()s to the working directory
    try:
        if not loo:
            mpm, mps, fas = gibbs(
                src_df, snk_df, alpha1=0.001, alpha2=0.1, beta=10,
                restarts=restarts, draws_per_restart=draws, burnin=burnin,
                delay=delay, jobs=1, create_feature_tables=True)
        else:
            res = []
            for source in src_df.index:
                cp = ConditionalProbability(0.001, 0.1, 10, src_df.drop(index=source).values)
                res.append(gibbs_sampler(
                    src_df.loc[source, :].values, cp, restarts, draws, burnin, delay))
            mpm, mps, fas = collate_gibbs_results(
                [i[0] for i in res], [i[1] for i in res], [i[2] for i in res],
                src_df.index, src_df.index, src_df.columns, True, loo=True)
    finally:
        os.chdir(old)
    return run, seed, mpm, mps, fas


TMP = os.path.join(OUT, "tmp")

if __name__ == "__main__":

    os.makedirs(TMP, exist_ok=True)
    jobs = [(run, seed) for run, cfg in RUNS.items() for seed in range(1, cfg[5] + 1)]
    # gibbs() starts a multiprocessing pool of its own, which the daemonic
    # workers of multiprocessing.Pool can't.
    with ProcessPoolExecutor(4) as p:
        results = list(p.map(one_run, jobs))
    for fn in os.listdir(TMP):
        os.remove(os.path.join(TMP, fn))
    os.rmdir(TMP)

    for run in RUNS:
        with open(os.path.join(OUT, f"sim-{run}.csv"), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["run", "seed", "sink", "proportions", "sd"])
            for r, seed, mpm, mps, fas in results:
                if r != run:
                    continue
                for sink in mpm.index:
                    w.writerow([run, seed, sink, hx(mpm.loc[sink].values), hx(mps.loc[sink].values)])
        r, seed, mpm, mps, fas = next(x for x in results if x[0] == run and x[1] == 1)
        with open(os.path.join(OUT, f"sim-{run}-fas.csv"), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["sink", "source"] + feat_ids)
            for sink, ft in zip(mpm.index, fas):
                for source, row in ft.iterrows():
                    w.writerow([sink, source] + [int(x) for x in row.values])
