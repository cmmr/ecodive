#!/usr/bin/env python3
# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
#
# Some lines repeat gemelli 0.0.13's optspace_helper() post-processing
# (BSD 3-Clause, Copyright (c) 2019, Cameron Martino); see inst/COPYRIGHTS.
"""
Reference results for ecodive::rpca(), from gemelli 0.0.13.

Writes the text fixtures data-raw/expected/rpca/<name>.txt, which
data-raw/expected/rpca.r converts into tests/tinytests/expected/rpca.rds for
tests/tinytests/test-rpca.r. There is no Python where the package is
developed, so the committed text files are the reference; this script is
only re-run to regenerate them.

Usage (from the package root):

    python -I data-raw/expected/rpca.py <gemelli-0.0.13 source dir> \
        data-raw/expected/rpca

Needs numpy, scipy, pandas, scikit-bio and biom-format (gemelli's own
requirements). gemelli's preprocessing module also imports `t2t` and `bp`,
which are only used for phylogenetic RPCA; they are stubbed out here if not
installed. gemelli itself is imported from its source directory and is not
modified.

For each table, the fixture holds:

  * The table and the arguments.
  * The end-to-end result of gemelli.rpca.rpca(): sample and feature
    loadings, eigenvalues (singular values), proportions explained, and the
    distance matrix.
  * For some tables, OptSpace's stages, from running
    gemelli.matrix_completion.MatrixCompletion on gemelli's matrix_rclr()
    of the table (as gemelli.rpca.optspace_helper() does), with
    gemelli.optspace's `svds` and `singular_values` wrapped to record their
    results: the rclr matrix, the initial factors (after OptSpace's
    scaling), the distortion after initialization and after each
    iteration, and OptSpace's final U, S, V. The final PCA repeats
    optspace_helper()'s post-processing lines with the same numpy/scipy
    calls.

gemelli's OptSpace starts from scipy's svds(), whose ARPACK start vector is
random. So that the fixtures are reproducible byte for byte, svds() is
called with a fixed seed here (gemelli.optspace.svds is replaced by a
wrapper that passes rng = numpy.random.default_rng(SEED)); the start only
changes the results at the level of ARPACK's tolerance.

The tables avoid robust clr values that come out exactly 0 for non-zero
counts, which gemelli (unlike ecodive) treats as missing.
"""

import sys
import types

import numpy as np
import scipy.linalg
from numpy.linalg import norm

# scikit-bio must be imported before any stubs.
import skbio  # noqa: F401
import skbio.diversity._util  # noqa: F401
import biom

for name in ['t2t', 't2t.nlevel', 'bp']:
    try:
        __import__(name)
    except ImportError:
        sys.modules[name] = types.ModuleType(name)
if not hasattr(sys.modules['t2t'], 'nlevel'):
    sys.modules['t2t'].nlevel = sys.modules['t2t.nlevel']
for attr in ['parse_newick', 'to_skbio_treenode']:
    if not hasattr(sys.modules['bp'], attr):
        setattr(sys.modules['bp'], attr, None)

sys.path.insert(0, sys.argv[1])

import gemelli  # noqa: E402
import gemelli.optspace as optspace  # noqa: E402
from gemelli.matrix_completion import MatrixCompletion  # noqa: E402
from gemelli.preprocessing import matrix_rclr  # noqa: E402
from gemelli.rpca import rpca as gemelli_rpca  # noqa: E402

assert gemelli.__version__ == '0.0.13', gemelli.__version__

OUT = sys.argv[2]

SEED = 20261117
_gemelli_svds = optspace.svds


def seeded_svds(obs, k, which='LM'):
    return _gemelli_svds(obs, k, which=which,
                         rng=np.random.default_rng(SEED))


optspace.svds = seeded_svds


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------

def fmt(x):
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if np.isnan(x):
        return 'NaN'
    return '%.17g' % x


def block(f, name, x):
    """A block: '@ name nrow ncol', then one line per row."""
    x = np.atleast_2d(np.asarray(x))
    f.write('@ %s %d %d\n' % (name, x.shape[0], x.shape[1]))
    for row in x:
        f.write(' '.join(fmt(v) for v in row) + '\n')


# ----------------------------------------------------------------------
# gemelli runs
# ----------------------------------------------------------------------

def end_to_end(counts, k, msc, mfc, mff):
    """gemelli.rpca.rpca() on a samples x features array."""
    n, p = counts.shape
    sids = ['S%d' % (i + 1) for i in range(n)]
    fids = ['F%d' % (j + 1) for j in range(p)]
    table = biom.Table(counts.T.astype(float), fids, sids)
    ord_res, dist_res = gemelli_rpca(
        table, n_components=k, min_sample_count=msc,
        min_feature_count=mfc, min_feature_frequency=mff)
    samples = ord_res.samples.iloc[:, :k]
    features = ord_res.features.iloc[:, :k]
    s_idx = [sids.index(s) + 1 for s in samples.index]
    f_idx = [fids.index(x) + 1 for x in features.index]
    dist = dist_res.filter(samples.index).data
    tri = dist[np.triu_indices(len(s_idx), 1)]  # R's dist order
    return dict(
        samples=np.array(s_idx), features=np.array(f_idx),
        u=samples.values, v=features.values,
        eig=ord_res.eigvals.values[:k],
        prop=ord_res.proportion_explained.values[:k],
        dist=tri)


def stages(counts, k):
    """OptSpace's stages on matrix_rclr(counts), recorded."""
    rec = dict(svds=[], sv=[])

    real_svds = optspace.svds
    real_sv = optspace.singular_values

    def svds(obs, k_, which):
        res = real_svds(obs, k_, which=which)
        rec['svds'].append(res)
        return res

    def singular_values(U, V, obs, mask):
        S = real_sv(U, V, obs, mask)
        err = norm(np.multiply(obs - U.dot(S).dot(V.T), mask), 'fro') / \
            np.sqrt(np.count_nonzero(mask))
        rec['sv'].append((U.copy(), V.copy(), S.copy(), err))
        return S

    optspace.svds = svds
    optspace.singular_values = singular_values
    try:
        rclr = matrix_rclr(counts.astype(float))
        opt = MatrixCompletion(n_components=k, max_iterations=5).fit(rclr)
    finally:
        optspace.svds = real_svds
        optspace.singular_values = real_sv

    n, m = rclr.shape
    U0, _, Vt0 = rec['svds'][0]
    u0 = U0 * np.sqrt(n)
    v0 = (Vt0 * np.sqrt(m)).T

    # optspace_helper()'s post-processing.
    X = opt.sample_weights @ opt.s @ opt.feature_weights.T
    X = X - X.mean(axis=0)
    X = X - X.mean(axis=1).reshape(-1, 1)
    u, s, v = scipy.linalg.svd(X)
    u = u[:, :k]
    v = v.T[:, :k]
    p = s**2 / np.sum(s**2)

    return dict(
        rclr=rclr, u0=u0, v0=v0,
        cost=np.array([e for (_, _, _, e) in rec['sv']]),
        U=opt.U, S=opt.s, V=opt.V,
        u=u, v=v, eig=s[:k], prop=p[:k], dist=opt.distance)


def penalty_stages(counts, k, step_size):
    """OptSpace with a small step_size, so that its Grassmann penalty
    terms (grassmann_manifold_one/two) are active. gemelli's rpca() always
    uses OptSpace's default step_size of 10000, where they are zero unless
    a row of U has a squared norm of at least 20000 * k."""
    rec = []
    real_sv = optspace.singular_values

    def singular_values(U, V, obs, mask):
        S = real_sv(U, V, obs, mask)
        err = norm(np.multiply(obs - U.dot(S).dot(V.T), mask), 'fro') / \
            np.sqrt(np.count_nonzero(mask))
        rec.append(err)
        return S

    optspace.singular_values = singular_values
    try:
        rclr = matrix_rclr(counts.astype(float))
        opt = optspace.OptSpace(n_components=k, max_iterations=5, tol=1e-5,
                                step_size=step_size)
        U, S, V = opt.solve(rclr)
    finally:
        optspace.singular_values = real_sv
    return dict(cost=np.array(rec), U=U, S=S, V=V)


# ----------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------

def clean(counts):
    """True if gemelli's missing values are exactly the zero counts."""
    with np.errstate(divide='ignore', invalid='ignore'):
        r = matrix_rclr(counts.astype(float))
    observed = np.isfinite(r) & (r != 0)
    return (np.all(observed == (counts > 0)) and
            np.all((counts > 0).sum(axis=1) >= 2) and
            np.all((counts > 0).sum(axis=0) >= 1))


def draw(make, seed):
    """make(rng) with the first seed from `seed` giving a clean table."""
    while True:
        counts = make(np.random.default_rng(seed))
        if clean(counts):
            return counts, seed
        seed += 1


def groups_table(rng, sizes, p, sd_group, sd_noise, depth, sd_depth,
                 dropout=0.0, sd_feature=1.0):
    """Counts from groups with different mean log-abundances."""
    base = rng.normal(0, sd_feature, p)
    rows = []
    for g, size in enumerate(sizes):
        mu = base + rng.normal(0, sd_group, p)
        for _ in range(size):
            logit = mu + rng.normal(0, sd_noise, p)
            prob = np.exp(logit - logit.max())
            prob /= prob.sum()
            d = depth * np.exp(rng.normal(0, sd_depth))
            x = rng.poisson(d * prob)
            if dropout > 0:
                x = x * (rng.random(p) >= dropout)
            rows.append(x)
    return np.array(rows, dtype=np.int64)


TABLES = [
    # name, maker, seed, k, filters, stages?
    ('groups', lambda rng: groups_table(
        rng, [7, 7, 6], 24, 1.5, 0.5, 2000, 0.8), 1, 3, (0, 0, 0), True),
    ('sparse', lambda rng: groups_table(
        rng, [12, 12], 60, 1.2, 0.6, 20000, 1.0, dropout=0.5,
        sd_feature=2.5), 100, 3, (0, 0, 0), True),
    ('filter', lambda rng: np.vstack([
        groups_table(rng, [9, 9], 25, 1.5, 0.5, 1000, 0.5, dropout=0.2),
        groups_table(rng, [2], 25, 1.5, 0.5, 30, 0.1)]),
        200, 2, (100, 20, 25), False),
    ('tall', lambda rng: groups_table(
        rng, [12, 12], 6, 1.0, 0.4, 500, 0.5), 300, 3, (0, 0, 0), False),
    ('wide', lambda rng: groups_table(
        rng, [3, 3], 40, 1.5, 0.5, 3000, 0.5, dropout=0.3), 400, 3,
        (0, 0, 0), True),
    ('khigh', lambda rng: groups_table(
        rng, [3, 3], 8, 1.5, 0.5, 1000, 0.5, dropout=0.2), 500, 5,
        (0, 0, 0), False),
    ('nozero', lambda rng: groups_table(
        rng, [5, 5], 8, 1.0, 0.3, 100000, 0.3) + 1, 600, 3, (0, 0, 0),
        False),
]


def main():
    for name, make, seed, k, (msc, mfc, mff), with_stages in TABLES:
        counts, seed = draw(make, seed)
        e2e = end_to_end(counts, k, msc, mfc, mff)
        if not with_stages:
            kept = counts[np.ix_(e2e['samples'] - 1, e2e['features'] - 1)]
            assert clean(kept), name

        path = '%s/%s.txt' % (OUT, name)
        with open(path, 'w') as f:
            f.write('# gemelli %s rpca() reference for ecodive; '
                    'generated by data-raw/expected/rpca.py\n'
                    % gemelli.__version__)
            f.write('# table seed %d, %d zeros of %d\n'
                    % (seed, np.sum(counts == 0), counts.size))
            block(f, 'counts', counts)
            block(f, 'k', [k])
            block(f, 'filters', [msc, mfc, mff])
            for key in ['samples', 'features', 'u', 'v', 'eig', 'prop',
                        'dist']:
                block(f, 'e2e_' + key, e2e[key] if e2e[key].ndim == 2
                      else [e2e[key]])

            if with_stages:
                st = stages(counts, k)
                for key in ['rclr', 'u0', 'v0', 'U', 'S', 'V', 'u', 'v']:
                    block(f, key, st[key])
                for key in ['cost', 'eig', 'prop']:
                    block(f, key, [st[key]])

                # The stage run repeats the end-to-end run, apart from
                # svds()'s random start.
                for key in ['eig', 'prop']:
                    assert np.allclose(st[key], e2e[key], rtol=1e-8), key
                for key in ['u', 'v']:
                    a, b = st[key], e2e[key]
                    b = b * np.sign(np.sum(a * b, axis=0))
                    assert np.allclose(a, b, atol=1e-8), key
            if name == 'groups':
                step = 2.0
                pen = penalty_stages(counts, k, step)
                block(f, 'pen_step', [step])
                for key in ['U', 'S', 'V']:
                    block(f, 'pen_' + key, pen[key])
                block(f, 'pen_cost', [pen['cost']])
        print(name, seed, counts.shape)


main()
