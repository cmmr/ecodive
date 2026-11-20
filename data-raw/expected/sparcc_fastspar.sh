#!/bin/sh
# FastSpar 1.0.0 (GPL-3; output only) on the tables that sparcc_refs.r
# writes. Run from the package root, with FastSpar's binaries on the PATH
# (built from the jsDelivr mirror of github.com/scwatts/fastspar):
#
#     sh data-raw/expected/sparcc_fastspar.sh
#
# Writes FastSpar's median correlation and covariance tables, named
# fastspar_<table>_i<iterations>_x<exclude iterations>_{cor,cov}.tsv.
set -e
dir=data-raw/expected/sparcc

run () {   # table iterations exclusion_iterations
  fastspar --otu_table "$dir/$1.tsv" \
    --correlation "$dir/fastspar_$1_i$2_x$3_cor.tsv" \
    --covariance  "$dir/fastspar_$1_i$2_x$3_cov.tsv" \
    --iterations "$2" --exclude_iterations "$3" --threshold 0.1 \
    --threads 4 --seed 1 --yes > /dev/null 2>&1
}

run tight    100 10
run tight    100 17   # D - 3: FastSpar stops there for 20 features ("Exclusion
                      # iterations ending early; we need at least 4
                      # non-excluded pairs for SparCC")
run tight    100 20   # the same as 17
run moderate 1000 20
run d8       200 5    # D - 3
run d8       200 10   # the same as 5

# p-values for `moderate`: 99 bootstrap tables (fastspar_bootstrap), each
# run with the observed run's settings, then fastspar_pvalues. The
# bootstrap tables and their correlations are left in $work (not deleted).
work=${TMPDIR:-/tmp}/ecodive-sparcc-fastspar
mkdir -p "$work"
run moderate 20 10
fastspar_bootstrap --otu_table "$dir/moderate.tsv" --number 99 \
  --prefix "$work/boot" --seed 1 > /dev/null 2>&1
for f in "$work"/boot_*.tsv; do
  b=$(basename "$f")
  fastspar --otu_table "$f" --correlation "$work/cor_$b" --covariance "$work/cov_$b" \
    --iterations 20 --exclude_iterations 10 --threshold 0.1 --seed 1 --yes > /dev/null 2>&1
done
fastspar_pvalues --otu_table "$dir/moderate.tsv" \
  --correlation "$dir/fastspar_moderate_i20_x10_cor.tsv" \
  --prefix "$work/cor_boot_" --permutations 99 \
  --outfile "$dir/fastspar_moderate_i20_x10_pvalues.tsv" > /dev/null 2>&1
