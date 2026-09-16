#!/bin/bash
# run_b2.sh -- B2 n_loci>=7 (e3v5v6, founder-consistent) NJ backbone build.
# Same in-memory NJ path as run_b1.sh, adapted to B2. Needs a big-memory box
# (see SIZING below) -- not runnable on a normal workstation.
#
# SIZING: peak ~= 20*n^2 bytes. B2 = 238,111 kept + 1 synthroot = 238,112 tips
# -> ~1.06 TB peak. Est. runtime ~2-3 h. NOTHING ELSE may run on the box near
# the ceiling.
#
# Usage (from within this repo, on the big-memory box):
#   bash 1_build_nj_backbone/run_b2.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # <repo>/1_build_nj_backbone
REPO="$(dirname "$SCRIPT_DIR")"                              # <repo> (tree_building/)
cd "$REPO"
# If your R packages (Rcpp, RcppParallel, data.table, ape) live in a
# non-default library location, export R_LIBS before running this script --
# not needed for a normal `install.packages()` setup.

SIDE=B2
TSV="processed_data/e3v5v6.B2_tape_consensus.ge7_founderok.tsv.gz"
DEFINING="processed_data/e3v5v6.defining_sites.tsv"
OUT="results/B2/e3v5v6_nj_ge7.nwk"
DTT_DIR="$PWD/decenttree"
LOG="$PWD/logs/run_v6_B2_ge7_inmem.log"
NTHREADS=40                                 # set to the box's PHYSICAL core count

# ---- 1. preflight ---------------------------------------------------------
for f in "$TSV" "$DEFINING" "1_build_nj_backbone/run_full_distance.R"; do
    [ -f "$f" ] || { echo "FATAL: missing required file: $f"; exit 1; }
done
[ -d "$DTT_DIR" ] || { echo "FATAL: decenttree source not found at $DTT_DIR"; exit 1; }
grep -q 'v6 = "processed_data"' lib/paths.R \
    || { echo "FATAL: paths.R does not register DATA_VERSION=v6"; exit 1; }
Rscript -e 'library(Rcpp); library(RcppParallel)' \
    || { echo "FATAL: R packages not loadable -- install.packages(c('Rcpp','RcppParallel'))"; exit 1; }
n_in=$(( $(zcat < "$TSV" | wc -l) - 1 ))
echo "[run_b2] input cells: $n_in (expect 238111); target tips: $((n_in + 1))"
[ "$n_in" -eq 238111 ] || echo "WARNING: input row count $n_in != expected 238111"
command -v numactl >/dev/null || sudo apt-get install -y numactl || true
mkdir -p "$(dirname "$OUT")" logs

# ---- 2. launch the in-memory NJ, DETACHED ---------------------------------
NUMACTL=""; command -v numactl >/dev/null && NUMACTL="numactl --interleave=all"
echo "[run_b2] launching (detached) -> $OUT   log: $LOG"
setsid $NUMACTL env \
    DATA_VERSION=v6 OUTGROUP_MODE=synthroot SYNTHROOT_SIDE=B2 \
    QC_PASS_MODE=all \
    DEFINING_SITES_TSV="$DEFINING" \
    DTT_TSV="$TSV" \
    DTT_NJ_INMEM="$OUT" DTT_NJ_NDIGITS=0 \
    DTT_DIR="$DTT_DIR" \
    DTT_NTHREADS="$NTHREADS" R_LIBS="$R_LIBS" \
    Rscript 1_build_nj_backbone/run_full_distance.R \
    > "$LOG" 2>&1 < /dev/null &

echo "[run_b2] started. Monitor:  tail -f $LOG   and   free -g  (RSS must stay under the box's RAM)"
echo "[run_b2] when done, tree at: $REPO/$OUT  (guard enforces 238112 tips)"
