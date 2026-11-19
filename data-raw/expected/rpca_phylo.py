#!/usr/bin/env python3
# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
"""
Reference results for ecodive::rpca_phylo(), from gemelli 0.0.13's
phylogenetic_rpca(), with scikit-bio 0.7.4 and iow 1.0.8 (bp).

Writes the text fixtures data-raw/expected/rpca_phylo/<case>.txt, which
data-raw/expected/rpca_phylo.r converts into
tests/tinytests/expected/rpca_phylo.rds for
tests/tinytests/test-rpca_phylo.r. There is no Python where the package is
developed, so the committed text files are the reference; this script is
only re-run to regenerate them.

Usage (from the package root):

    python -I data-raw/expected/rpca_phylo.py <gemelli-0.0.13 source dir> \
        <gemelli's gemelli/q2/tests/data/tree.nwk> data-raw/expected/rpca_phylo

Needs numpy, scipy, pandas, biom-format, scikit-bio 0.7.4 and iow 1.0.8
(which provides `bp`). gemelli is imported from its source directory and is
not modified; its taxonomy module (`t2t`, only used with a taxonomy) is
stubbed out if not installed. gemelli's test tree.nwk is not in its sdist;
it is in gemelli's repository (gemelli/q2/tests/data/tree.nwk).

Each case gives gemelli a Newick tree and a table, and records:

  * The inputs: the Newick string, the counts (samples x tips, integers),
    the sample and feature ids, k, the filters and min_depth.
  * gemelli's node table, from the same calls as phylogenetic_rpca()
    (rpca_table_processing(), bp_read_phylogeny(), fast_unifrac()): the
    node ids (fids), each node's parent (as a 1-based column, 0 for the
    root), the branch lengths, the node counts, and matrix_rclr()'s robust
    clr matrix (NaN where missing).
  * The end-to-end result of phylogenetic_rpca(): the kept samples, sample
    and feature loadings, eigenvalues (singular values), proportions
    explained, and the distances.

Every real number is written as a hexadecimal float (float.hex()), which R
reads back exactly.

gemelli's OptSpace starts from scipy's svds(), whose ARPACK start vector is
random. So that the fixtures are reproducible byte for byte, svds() is
called with a fixed seed (as data-raw/expected/rpca.py does); the start only
changes the results at the level of ARPACK's tolerance.

The cases avoid robust clr values that come out exactly 0 for observed
entries, which gemelli (unlike ecodive) treats as missing: the script
checks that gemelli's missing entries are exactly ecodive's (count or
branch length 0, or a sample whose values are all equal).
"""

import io
import os
import sys
import tempfile
import types

import numpy as np

import skbio
import biom

for name in ['t2t', 't2t.nlevel']:
    try:
        __import__(name)
    except ImportError:
        sys.modules[name] = types.ModuleType(name)
if not hasattr(sys.modules['t2t'], 'nlevel'):
    sys.modules['t2t'].nlevel = sys.modules['t2t.nlevel']

sys.path.insert(0, sys.argv[1])

import importlib.metadata  # noqa: E402

import bp  # noqa: E402
import gemelli  # noqa: E402
import gemelli.optspace as optspace  # noqa: E402
from gemelli.preprocessing import (bp_read_phylogeny, fast_unifrac,  # noqa: E402
                                   matrix_rclr)
from gemelli.rpca import phylogenetic_rpca, rpca_table_processing  # noqa: E402

assert gemelli.__version__ == '0.0.13', gemelli.__version__
assert skbio.__version__ == '0.7.4', skbio.__version__
assert importlib.metadata.version('iow') == '1.0.8'

TREE_NWK = sys.argv[2]
OUT = sys.argv[3]
TMP = os.path.join(tempfile.mkdtemp(), 'tree.nwk')   # gemelli reads files

SEED = 20261119
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
    x = float(x)
    if np.isnan(x):
        return 'NaN'
    return x.hex()


def block(f, name, x):
    """A numeric block: '@ name nrow ncol', then one line per row."""
    x = np.atleast_2d(np.asarray(x))
    f.write('@ %s %d %d\n' % (name, x.shape[0], x.shape[1]))
    for row in x:
        f.write(' '.join(fmt(v) for v in row) + '\n')


def strings(f, name, x):
    """A block of strings: '$ name n', then one per line."""
    f.write('$ %s %d\n' % (name, len(x)))
    for s in x:
        assert '\n' not in s
        f.write(s + '\n')


# ----------------------------------------------------------------------
# Trees
# ----------------------------------------------------------------------

def rand_length(rng, p_zero=0.08, p_missing=0.03):
    u = rng.random()
    if u < p_zero:
        return ':0'
    if u < p_zero + p_missing:
        return ''
    return ':%.3f' % rng.uniform(0.02, 1.0)


def rand_tree(rng, tips, p_multi=0.25, p_unary=0.08, p_label=0.3):
    """A random Newick tree over `tips` (in order), with multifurcations,
    single-child nodes, zero and missing lengths, and some internal
    labels. Returns the string without the root's ';'."""
    count = [0]

    def build(ts):
        if len(ts) == 1:
            s = ts[0]
        else:
            nc = 2
            if len(ts) >= 3 and rng.random() < p_multi:
                nc = int(rng.integers(3, min(len(ts), 4) + 1))
            cuts = sorted(rng.choice(np.arange(1, len(ts)), nc - 1,
                                     replace=False))
            groups = np.split(np.array(ts, dtype=object), cuts)
            s = '(' + ','.join(build(list(g)) for g in groups) + ')'
            if rng.random() < p_label:
                count[0] += 1
                s += 'x%d' % count[0]
        if rng.random() < p_unary:
            s = '(' + s + rand_length(rng) + ')u%d' % count[0]
            count[0] += 1
        return s + rand_length(rng)

    ts = list(tips)
    nc = 3 if rng.random() < 0.5 else 2
    cuts = sorted(rng.choice(np.arange(1, len(ts)), nc - 1, replace=False))
    groups = np.split(np.array(ts, dtype=object), cuts)
    return '(' + ','.join(build(list(g)) for g in groups) + ')'


def gemelli_subtree(n_tips, rng):
    """gemelli's test tree.nwk, sheared (not collapsed) to `n_tips` of its
    tips, so that single-child nodes remain; as a Newick string."""
    with open(TREE_NWK) as fh:
        full = bp.parse_newick(fh.readline())
    names = sorted(n.name for n in
                   bp.to_skbio_treenode(full).tips())
    keep = set(rng.choice(names, n_tips, replace=False))
    sheared = bp.to_skbio_treenode(full.shear(keep))
    out = io.StringIO()
    sheared.write(out)
    return out.getvalue().strip()


# ----------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------

def groups_table(rng, sizes, p, sd_group, sd_noise, depth, sd_depth,
                 dropout=0.0, sd_feature=1.0):
    """Counts from groups with different mean log-abundances."""
    base = rng.normal(0, sd_feature, p)
    rows = []
    for size in sizes:
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


# ----------------------------------------------------------------------
# gemelli
# ----------------------------------------------------------------------

def run(case):
    counts, sids, fids_in = case['counts'], case['sids'], case['fids']
    with open(TMP, 'w') as fh:
        fh.write(case['newick'] + '\n')
    table = biom.Table(counts.T.astype(float), fids_in, sids)
    msc, mfc, mff = case['filters']
    k, md = case['k'], case['min_depth']

    # The stages, as phylogenetic_rpca() runs them.
    t = rpca_table_processing(table, msc, mfc, mff)
    tree = bp_read_phylogeny(t, TMP, md)
    cbn, tree_index, lengths, fids, _ = fast_unifrac(t, tree)
    assert all(tree_index['keep'].values())
    with np.errstate(divide='ignore', invalid='ignore'):
        rclr = matrix_rclr(cbn, branch_lengths=lengths)

    # Each column's parent column (1-based; 0 for the root).
    K = len(fids)
    parent = np.zeros(K, dtype=np.int64)
    for node_id, lo, hi in tree_index['child_index']:
        parent[lo:hi + 1] = node_id + 1
    assert parent[K - 1] == 0 and np.all(parent[:K - 1] > 0)

    # gemelli's missing entries must be exactly ecodive's.
    obs = (cbn > 0) & (lengths > 0)
    v = np.where(obs, np.log(np.where(obs, cbn, 1.0)) +
                 np.log(np.where(lengths > 0, lengths, 1.0)), np.nan)
    flat = np.array([(not np.any(o)) or np.max(r[o]) == np.min(r[o])
                     for r, o in zip(v, obs)])
    ours = obs & ~flat[:, None]
    gem = np.isfinite(rclr) & (rclr != 0)
    assert np.array_equal(ours, gem), case['name']

    # End to end.
    ord_res, dist_res, _, cbn_table, _ = phylogenetic_rpca(
        table, TMP, n_components=k, min_sample_count=msc,
        min_feature_count=mfc, min_feature_frequency=mff, min_depth=md)
    samples = list(ord_res.samples.index)
    assert samples == list(t.ids())
    assert list(ord_res.features.index) == fids
    assert np.array_equal(cbn_table.matrix_data.toarray().T, cbn)
    dist = dist_res.filter(samples).data
    tri = dist[np.triu_indices(len(samples), 1)]   # R's dist order

    return dict(
        samples=samples, fids=fids, parent=parent, lengths=lengths,
        node_counts=cbn, rclr=rclr, flat=flat.astype(np.int64),
        u=ord_res.samples.values[:, :k], v=ord_res.features.values[:, :k],
        eig=ord_res.eigvals.values[:k],
        prop=ord_res.proportion_explained.values[:k], dist=tri)


# ----------------------------------------------------------------------
# Cases
# ----------------------------------------------------------------------

def ids(prefix, n):
    return ['%s%d' % (prefix, i + 1) for i in range(n)]


def make_cases():
    cases = []

    # gemelli's own example (phylogenetic_rpca()'s docstring), whose
    # Newick root length (n1:1) gemelli drops.
    X = np.array([[9, 3, 0, 0], [9, 9, 0, 1], [0, 1, 4, 5],
                  [0, 0, 3, 4], [1, 0, 8, 9]]).T
    cases.append(dict(
        name='demo', counts=X, sids=['s1', 's2', 's3', 's4'],
        fids=['f1', 'f2', 'f3', 'f4', 'f5'],
        newick='(((f1:1,f2:1)n9:1,f3:1)n8:1,(f4:1,f5:1)n2:1)n1:1;',
        k=2, filters=(0, 0, 0), min_depth=0))

    # gemelli's test tree, cut to 150 tips, with a table on 100 of them.
    rng = np.random.default_rng(SEED)
    nwk = gemelli_subtree(150, rng)
    tips = [n.name for n in skbio.TreeNode.read([nwk]).tips()]
    used = list(rng.choice(tips, 100, replace=False))
    counts = groups_table(rng, [8, 8, 8], 100, 1.5, 0.5, 3000, 0.6,
                          dropout=0.4, sd_feature=1.5)
    cases.append(dict(
        name='gemelli', counts=counts, sids=ids('S', 24), fids=used,
        newick=nwk, k=3, filters=(0, 0, 0), min_depth=0))
    cases.append(dict(
        name='gemelli_filter', counts=counts, sids=ids('S', 24), fids=used,
        newick=nwk, k=4, filters=(1500, 30, 20), min_depth=3))

    # Random trees: multifurcations, single-child nodes, zero and missing
    # lengths, internal labels, and tips not in the table.
    rng = np.random.default_rng(SEED + 1)
    tips = ids('t', 40)
    nwk = rand_tree(rng, tips) + 'root:0.25;'
    used = list(rng.choice(tips, 32, replace=False))
    counts = groups_table(rng, [7, 7, 6], 32, 1.2, 0.5, 2000, 0.5,
                          dropout=0.3)
    cases.append(dict(
        name='random', counts=counts, sids=ids('S', 20), fids=used,
        newick=nwk, k=3, filters=(0, 0, 0), min_depth=0))
    cases.append(dict(
        name='random_depth', counts=counts, sids=ids('S', 20), fids=used,
        newick=nwk, k=3, filters=(0, 0, 0), min_depth=4))

    # A root left with a single child (all used tips in one clade, under a
    # chain of single-child nodes): gemelli merges them, adding the root
    # length once.
    rng = np.random.default_rng(SEED + 2)
    clade = rand_tree(rng, ids('a', 12), p_unary=0)
    nwk = ('(((' + clade + 'c1:0.4)c2:0.3)c3:0.2,' +
           rand_tree(rng, ids('b', 6), p_unary=0) + 'c4:0.7)root:0.6;')
    counts = groups_table(rng, [6, 6], 12, 1.5, 0.5, 1000, 0.5,
                          dropout=0.25)
    cases.append(dict(
        name='single_root', counts=counts, sids=ids('S', 12),
        fids=ids('a', 12), newick=nwk, k=3, filters=(0, 0, 0),
        min_depth=0))

    # A star tree of unit lengths: the tips plus an all-missing root.
    rng = np.random.default_rng(SEED + 3)
    counts = groups_table(rng, [5, 5], 8, 1.5, 0.4, 500, 0.4, dropout=0.2)
    cases.append(dict(
        name='star', counts=counts, sids=ids('S', 10), fids=ids('t', 8),
        newick='(' + ','.join('t%d:1' % (i + 1) for i in range(8)) + ');',
        k=3, filters=(0, 0, 0), min_depth=0))

    # Samples whose values are all equal (a single tip under unit
    # lengths): rows with no observed values in gemelli as here.
    rng = np.random.default_rng(SEED + 4)
    counts = groups_table(rng, [5, 5], 10, 1.5, 0.4, 800, 0.4, dropout=0.2)
    counts[0] = 0
    counts[0, 3] = 17
    counts[4] = 0
    counts[4, 7] = 5
    nwk = ('(((t1:1,t2:1):1,(t3:1,t4:1):1):1,((t5:1,t6:1):1,'
           '((t7:1,t8:1):1,(t9:1,t10:1):1):1):1);')
    cases.append(dict(
        name='flat', counts=counts, sids=ids('S', 10), fids=ids('t', 10),
        newick=nwk, k=3, filters=(0, 0, 0), min_depth=0))

    return cases


def main():
    for case in make_cases():
        res = run(case)
        path = '%s/%s.txt' % (OUT, case['name'])
        with open(path, 'w') as f:
            f.write('# gemelli %s phylogenetic_rpca() reference for '
                    'ecodive; generated by data-raw/expected/rpca_phylo.py\n'
                    % gemelli.__version__)
            strings(f, 'newick', [case['newick']])
            strings(f, 'sample_ids', case['sids'])
            strings(f, 'feature_ids', case['fids'])
            block(f, 'counts', case['counts'])
            block(f, 'k', [case['k']])
            block(f, 'filters', [list(case['filters'])])
            block(f, 'min_depth', [case['min_depth']])
            strings(f, 'samples', res['samples'])
            strings(f, 'fids', res['fids'])
            block(f, 'parent', [res['parent']])
            block(f, 'lengths', [res['lengths']])
            block(f, 'node_counts', res['node_counts'])
            block(f, 'flat', [res['flat']])
            block(f, 'rclr', res['rclr'])
            block(f, 'u', res['u'])
            block(f, 'v', res['v'])
            block(f, 'eig', [res['eig']])
            block(f, 'prop', [res['prop']])
            block(f, 'dist', [res['dist']])
        print(case['name'], case['counts'].shape, len(res['fids']),
              'nodes', int(res['flat'].sum()), 'flat')


main()
