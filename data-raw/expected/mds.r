# Expected results for tests/tinytests/test-mds.r, from scikit-learn's
# metric SMACOF. Run from the package root:
#
#     Rscript data-raw/expected/mds.r
#
# Writes tests/tinytests/expected/mds.rds.
#
# The reference values themselves come from data-raw/expected/mds.py (Python
# with scikit-learn 1.9.1), whose text outputs are committed under
# data-raw/expected/mds/, from the dissimilarities under
# data-raw/expected/nmds/: rerun it only to change the cases. This script
# only converts those files. mds.rds has the same layout as nmds.rds (see
# data-raw/expected/nmds.r), with mds() measuring stress as
# sqrt(raw stress / sum of squared dissimilarities).

source(file.path('data-raw', 'expected', 'nmds.r'))

saveRDS(
  smacof_fixtures(file.path('data-raw', 'expected', 'mds'), 'mds'),
  file.path('tests', 'tinytests', 'expected', 'mds.rds'),
  version = 2 )
