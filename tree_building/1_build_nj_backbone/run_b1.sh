#!/bin/bash
# run_b1.sh -- ONE-COMMAND B1 n_loci>=7 (e3v8, founder-consistent) NJ backbone
# build. Needs a big-memory box (see SIZING below) -- not runnable on a normal
# workstation.
#
# This is the exact in-memory NJ invocation used for the v8 build
# (run_full_distance.R's DTT_NJ_INMEM branch -> decenttree RapidNJ, no PHYLIP
# round-trip), run against the v8 ge7 founder-filtered B1 input (derived by
# 00_filter_highqual_consensus.R + 00a_reshape_v8_founder_genotypes.R).
#
# NOTE: The memory intense NJ tree building was run on a google cloud machine.
#
# SIZING (RAM is the hard constraint): peak ~= 20*n^2 bytes. B1 has 410,925 kept
# cells + 1 synthetic root = 410,926 tips -> ~3.1 TB peak. Est. runtime ~8-10 h
# (matrix build ~15 min, NJ the rest). NOTHING ELSE may run on the box near the
# ceiling.
#
# Usage (from within this repo, on the big-memory box):
#   bash 1_build_nj_backbone/run_b1.sh
# Launches DETACHED (survives SSH drop). Monitor with:
#   tail -f logs/run_v8_B1_ge7_inmem.log ; free -g
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # <repo>/1_build_nj_backbone
REPO="$(dirname "$SCRIPT_DIR")"                              # <repo> (tree_building/)
cd "$REPO"
# If your R packages (Rcpp, RcppParallel, data.table, ape) live in a
# non-default library location, export R_LIBS before running this script --
# not needed for a normal `install.packages()` setup.

SIDE=B1
TSV="processed_data/e3v8.B1_tape_consensus.ge7_founderok.tsv.gz"
DEFINING="../support_data/e3v8.defining_sites.tsv"
OUT="results/1-nj-backbone/nj_raw_B1.nwk"
DTT_DIR="$PWD/decenttree"
LOG="$PWD/logs/run_v8_B1_ge7_inmem.log"
NTHREADS=64                                 # set to the box's PHYSICAL core count

# ---- 1. preflight ---------------------------------------------------------
for f in "$TSV" "$DEFINING" "1_build_nj_backbone/run_full_distance.R"; do
    [ -f "$f" ] || { echo "FATAL: missing required file: $f"; exit 1; }
done
[ -d "$DTT_DIR" ] || { echo "FATAL: decenttree source not found at $DTT_DIR"; exit 1; }
# confirm DATA_VERSION=v8 is registered (paths.R must know v8, else it stop()s)
grep -q 'v8 = "processed_data"' lib/paths.R \
    || { echo "FATAL: paths.R does not register DATA_VERSION=v8"; exit 1; }
Rscript -e 'library(Rcpp); library(RcppParallel)' \
    || { echo "FATAL: R packages not loadable -- install.packages(c('Rcpp','RcppParallel'))"; exit 1; }
n_in=$(( $(zcat < "$TSV" | wc -l) - 1 ))
echo "[run_b1] input cells: $n_in (expect 410925); target tips: $((n_in + 1))"
[ "$n_in" -eq 410925 ] || echo "WARNING: input row count $n_in != expected 410925"
command -v numactl >/dev/null || sudo apt-get install -y numactl || true
mkdir -p "$(dirname "$OUT")" logs

# ---- 2. launch the in-memory NJ, DETACHED ---------------------------------
NUMACTL=""; command -v numactl >/dev/null && NUMACTL="numactl --interleave=all"
echo "[run_b1] launching (detached) -> $OUT   log: $LOG"
setsid $NUMACTL env \
    DATA_VERSION=v8 OUTGROUP_MODE=synthroot SYNTHROOT_SIDE=B1 \
    QC_PASS_MODE=all \
    DEFINING_SITES_TSV="$DEFINING" \
    DTT_TSV="$TSV" \
    DTT_NJ_INMEM="$OUT" DTT_NJ_NDIGITS=0 \
    DTT_DIR="$DTT_DIR" \
    DTT_NTHREADS="$NTHREADS" R_LIBS="$R_LIBS" \
    Rscript 1_build_nj_backbone/run_full_distance.R \
    > "$LOG" 2>&1 < /dev/null &

echo "[run_b1] started. Monitor:  tail -f $LOG   and   free -g  (RSS must stay under the box's RAM)"
echo "[run_b1] when done, the tree is at: $REPO/$OUT  (internal tip-count guard enforces 410926 tips)"
