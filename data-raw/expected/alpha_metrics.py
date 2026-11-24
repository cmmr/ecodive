"""Reference fixtures for ecodive's alpha diversity metrics.

Runs scikit-bio's skbio.diversity.alpha functions on every row (sample) of
the input tables that data-raw/expected/alpha_metrics.r writes, and writes
the results under data-raw/expected/alpha_metrics/ (which
data-raw/expected/alpha_metrics.r then folds into
tests/tinytests/expected/alpha_metrics.rds).

Input, from `Rscript data-raw/expected/alpha_metrics.r inputs`:

  table-<name>.csv  Integer counts, samples x features, no header.

Output:

  results.csv       One row per table, sample (1-based), metric and order:
                    table, sample, metric (ecodive's ID), q (the order for
                    hill, renyi and tsallis, else empty; 'inf' for
                    infinity), skbio (the scikit-bio call), value (a hex
                    float from float.hex(), which R reads back exactly, or
                    NaN / Inf).

Each metric is called on the sample's counts with its zeros, as scikit-bio
strips zeros itself; gini_index() keeps them, so that it sees every feature
of the table.

Prerequisites: Python 3 with numpy and scikit-bio 0.7.4 (the version the
fixtures were made with). Run from the repository root:

    Rscript data-raw/expected/alpha_metrics.r inputs
    python data-raw/expected/alpha_metrics.py
    Rscript data-raw/expected/alpha_metrics.r
"""

import csv
import glob
import math
import os
import warnings

import numpy as np
import skbio
from skbio.diversity import alpha

assert skbio.__version__ == '0.7.4', skbio.__version__

DIR = os.path.join('data-raw', 'expected', 'alpha_metrics')

ORDERS = [0, 0.5, 1, 2, 3, math.inf]

# ecodive ID -> scikit-bio function, for the metrics without an order.
METRICS = [
    ('pielou',            'pielou_e',        alpha.pielou_e),
    ('heip',              'heip_e',          alpha.heip_e),
    ('simpson_evenness',  'simpson_e',       alpha.simpson_e),
    ('mcintosh_evenness', 'mcintosh_e',      alpha.mcintosh_e),
    ('goods_coverage',    'goods_coverage',  alpha.goods_coverage),
    ('dominance',         'dominance',       alpha.dominance),
    ('strong',            'strong',          alpha.strong),
    ('gini',              'gini_index',      alpha.gini_index),
    ('chao1',             'chao1',           lambda x: alpha.chao1(x, bias_corrected=True)),
    ('chao1_classic',     'chao1(bias_corrected=False)',
                                             lambda x: alpha.chao1(x, bias_corrected=False)),
    # Existing ecodive metrics, as a cross-check of the aliases.
    ('observed',          'observed_features', alpha.observed_features),
    ('shannon',           'shannon',         alpha.shannon),
    ('simpson',           'simpson',         alpha.simpson),
    ('inv_simpson',       'inv_simpson',     alpha.inv_simpson),
    ('berger',            'berger_parker_d', alpha.berger_parker_d),
    ('mcintosh',          'mcintosh_d',      alpha.mcintosh_d),
    ('brillouin',         'brillouin_d',     alpha.brillouin_d),
    ('margalef',          'margalef',        alpha.margalef),
    ('menhinick',         'menhinick',       alpha.menhinick),
]

ORDERED = [
    ('hill',    alpha.hill),
    ('renyi',   alpha.renyi),
    ('tsallis', alpha.tsallis),
]


def fmt(v):
    v = float(v)
    if math.isnan(v):
        return 'NaN'
    if math.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    return v.hex()


def fmt_q(q):
    return 'inf' if math.isinf(q) else repr(q)


rows = []
for path in sorted(glob.glob(os.path.join(DIR, 'table-*.csv'))):
    name = os.path.basename(path)[len('table-'):-len('.csv')]
    tab  = np.loadtxt(path, delimiter=',', dtype=np.int64, ndmin=2)
    for i, x in enumerate(tab):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with np.errstate(all='ignore'):
                for mid, call, fun in METRICS:
                    rows.append([name, i + 1, mid, '', call, fmt(fun(x))])
                for mid, fun in ORDERED:
                    for q in ORDERS:
                        rows.append([
                            name, i + 1, mid, fmt_q(q),
                            '%s(order=%s)' % (mid, fmt_q(q)),
                            fmt(fun(x, order=q)) ])

with open(os.path.join(DIR, 'results.csv'), 'w', newline='') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['table', 'sample', 'metric', 'q', 'skbio', 'value'])
    w.writerows(rows)

print('Wrote %d results with scikit-bio %s.' % (len(rows), skbio.__version__))
