#!/usr/bin/env python3
# Copyright (c) 2026 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit
"""
Reference bases for ecodive::philr(), from scikit-bio 0.7.4's tree_basis()
(BSD-3): the unweighted ILR basis of a binary tree, in clr coordinates.

Writes the text fixtures data-raw/expected/philr/skbio-<case>.txt, which
data-raw/expected/philr.r folds into tests/tinytests/expected/philr.rds for
tests/tinytests/test-philr.r. There is no Python where the package is
developed, so the committed text files are the reference; this script is
only re-run to regenerate them.

Usage (from the package root):

    python -I data-raw/expected/philr.py data-raw/expected/philr

Each file is a series of blocks: "$ name n" then n lines of strings, or
"@ name nrow ncol" then one line per row of hexadecimal floats
(float.hex(), which R reads back exactly):

  newick   the tree, with every internal node named;
  tips     the tips in tree.tips() order: the basis's columns;
  nodes    the internal nodes in tree_basis()'s order: its rows;
  basis    tree_basis(tree), (D - 1) x D.

tree_basis() puts each node's *second* child in the numerator and orders
its rows level by level; philr() puts the first child in the numerator and
orders nodes by number (preorder). The test matches rows by node name and
flips the sign.
"""

import os
import sys

import numpy as np
import skbio
from skbio import TreeNode
from skbio.stats.composition import tree_basis

assert skbio.__version__ == '0.7.4', skbio.__version__

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)


def random_newick(n_tips, rng):
    """A random rooted binary tree with named tips (t1..) and nodes (n1..)."""
    nodes = ['t%d:%s' % (i + 1, repr(round(float(rng.uniform(0.05, 1)), 4)))
             for i in range(n_tips)]
    k = 0
    while len(nodes) > 1:
        i, j = sorted(rng.choice(len(nodes), 2, replace=False))
        b, a = nodes.pop(j), nodes.pop(i)
        k += 1
        nodes.append('(%s,%s)m%d:%s' % (a, b, k, repr(round(float(rng.uniform(0.05, 1)), 4))))
    root = nodes[0]
    return root[:root.rindex(':')] + ';'


def write_case(name, newick):
    tree = TreeNode.read([newick])
    basis, nodes = tree_basis(tree)
    basis = basis.toarray()
    tips = [t.name for t in tree.tips()]
    assert basis.shape == (len(tips) - 1, len(tips))
    with open(os.path.join(OUT, 'skbio-%s.txt' % name), 'w', newline='\n') as f:
        f.write('# scikit-bio %s tree_basis()\n' % skbio.__version__)
        f.write('$ newick 1\n%s\n' % newick)
        f.write('$ tips %d\n' % len(tips))
        f.write(''.join(t + '\n' for t in tips))
        f.write('$ nodes %d\n' % len(nodes))
        f.write(''.join(str(n) + '\n' for n in nodes))
        f.write('@ basis %d %d\n' % basis.shape)
        for row in basis:
            f.write(' '.join(float(v).hex() for v in row) + '\n')


write_case(
    'small',
    '(((a:0.3,b:0.2)ab:0.1,c:0.6)abc:0.2,((d:0.4,e:0.1)de:0.3,(f:0.2,(g:0.1,h:0.5)gh:0.2)fgh:0.1)defgh:0.4)root;')

write_case('random', random_newick(30, np.random.default_rng(20261119)))
