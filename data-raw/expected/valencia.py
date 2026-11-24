# Runs VALENCIA's Valencia.py on the inputs that data-raw/expected/valencia.r
# writes, and keeps its outputs: for each case, data-raw/expected/valencia/
# <case>-output.csv, with sampleID, the 13 similarities (<sub-CST>_sim, as
# hexadecimal floats, read back exactly), subCST, score (hexadecimal) and CST.
# Run from the package root:
#
#     python3 data-raw/expected/valencia.py <Valencia.py> <CST_centroids_012920.csv>
#
# both from VALENCIA at commit 8559d45 (the centroids file is the one that
# data-raw/valencia_centroids.r builds data/valencia_centroids.rda from).
#
# Valencia.py is VALENCIA's classifier (commit 8559d45 of
# github.com/ravel-lab/VALENCIA, MIT; Copyright (c) 2023 ravel-lab), run as
# its README says, with pandas 2.2.3 and numpy 1.26.4. (Valencia.py indexes a
# centroid's values by position, which pandas 3 no longer allows.)

import os
import subprocess
import sys
import tempfile

import pandas as pd

valencia = sys.argv[1]
d = os.path.join('data-raw', 'expected', 'valencia')

subs = ['I-A', 'I-B', 'II', 'III-A', 'III-B', 'IV-A', 'IV-B',
        'IV-C0', 'IV-C1', 'IV-C2', 'IV-C3', 'IV-C4', 'V']


def hexf(x):
    return float(x).hex()


with tempfile.TemporaryDirectory() as tmp:

    ref = {
        'synthetic':   sys.argv[2],
        'ties':        os.path.join(d, 'ties-centroids.csv'),
        'publication': sys.argv[2],
    }

    for case, centroids in ref.items():
        out = os.path.join(tmp, case)
        subprocess.run(
            [sys.executable, '-I', valencia, '-ref', centroids,
             '-i', os.path.join(d, case + '-input.csv'), '-o', out],
            check=True)
        res = pd.read_csv(out + '.csv', float_precision='round_trip')
        keep = pd.DataFrame({'sampleID': res['sampleID']})
        for s in subs:
            keep[s + '_sim'] = [hexf(x) for x in res[s + '_sim']]
        keep['subCST'] = res['subCST']
        keep['score'] = [hexf(x) for x in res['score']]
        keep['CST'] = res['CST']
        keep.to_csv(os.path.join(d, case + '-output.csv'), index=False)
