#!/usr/bin/env bash
# run_side.sh -- DTT-ancestor placement onto the DATED per-side backbone, for one side.
#
#   ./run_side.sh B1
#   ./run_side.sh B2
#
# STEP 1  bestmatch.py       DTT-ancestor nearest-neighbour (genotype only)
# STEP 2  finalize_dated.py  d(C,P) pendant (matrix DTT units) -> /clock_rate -> days,
#                            grafted onto the dated backbone; writes the dated tree +
#                            place_qc.tsv (matched cell, attach score, pendant, days,
#                            validity).
#
# Run BOTH sides, then merge with 3_date_tree/merge_dated_subtrees.R (same script
# already used to merge the dated backbone -- dated_placed_{SIDE}.nwk is a dated tree
# just like that step's output):
#   LEFT_NWK=results/4-full-tree-placement/B1/dated_placed_B1.nwk \
#   RIGHT_NWK=results/4-full-tree-placement/B2/dated_placed_B2.nwk \
#   OUT_NWK=results/4-full-tree-placement/merged/merged_placed.nwk \
#     Rscript 3_date_tree/merge_dated_subtrees.R
set -euo pipefail
cd "$(dirname "$0")"
REPO="$(cd .. && pwd)"

SIDE="${1:?usage: run_side.sh B1|B2}"
export SIDE
export ALLCELLS="${ALLCELLS:-0}"        # 0 = pass_qc cells only (queries = pass_qc, <7 loci etc.)
export N_SITES="${N_SITES:-6}"

# ---- the dated backbone (3_date_tree's constrained-dating output) -------------
VARIANT="${VARIANT:-attempt1_minage_sourced_minB2h_l0.01}"
VDIR="$REPO/results/3-dated-tree/$SIDE/lsd2/nj99478/$VARIANT"
export BACKBONE_NWK="${BACKBONE_NWK:-$VDIR/time_tree.nwk}"
[ -f "$BACKBONE_NWK" ] || { echo "backbone not found: $BACKBONE_NWK" >&2; exit 1; }

# ---- clock rate: LSD2's own fitted rate for this side's dating run (see 3_date_tree,
# results/3-dated-tree/$SIDE/lsd2/nj99478/$VARIANT/dated -- "rate X, tMRCA" line).
# Hardcoded here (not re-parsed at run time) so the value is fixed and reviewable;
# override with env RATE if you re-date with a different backbone.
case "$SIDE" in
  B1) export RATE="${RATE:-0.2508}" ;;   # DTT/day, LSD2 iter3 fit
  B2) export RATE="${RATE:-0.3019}" ;;   # DTT/day, LSD2 iter3 fit
  *)  echo "unknown side: $SIDE" >&2; exit 1 ;;
esac

# ---- consensus (queries): the FULL per-side callset, not the ge7-founder-ok subset
# the backbone was built from -- placement's query universe is every cell.
export CONSENSUS_TSV="${CONSENSUS_TSV:-$REPO/processed_data/e3v5v6.${SIDE}_tape_consensus.tsv.gz}"

# ---- output ---------------------------------------------------------------------
OUTBASE="$REPO/results/4-full-tree-placement/${SIDE}"
mkdir -p "$OUTBASE"
export OUTDIR="${OUTDIR:-$OUTBASE/phase}"        # bestmatch checkpoints
export DATED_TREE_OUT="${DATED_TREE_OUT:-$OUTBASE/dated_placed_${SIDE}.nwk}"
export QC_OUT="${QC_OUT:-$OUTBASE/place_qc.tsv}"

# ---- self-documenting manifest of exactly what fed this run ---------------------
{
  echo "run          $(date '+%Y-%m-%dT%H:%M:%S%z')"
  echo "git_commit   $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo NA)"
  echo "side         $SIDE"
  echo "variant      $VARIANT"
  echo "backbone     $BACKBONE_NWK"
  echo "backbone_md5 $(md5 -q "$BACKBONE_NWK" 2>/dev/null || md5sum "$BACKBONE_NWK" | cut -d' ' -f1)"
  echo "consensus    $CONSENSUS_TSV"
  echo "allcells     $ALLCELLS"
  echo "rate         $RATE"
} > "$OUTBASE/manifest.txt"

echo ">>> [$SIDE] backbone : $BACKBONE_NWK"
echo ">>> [$SIDE] consensus: $CONSENSUS_TSV"
echo ">>> [$SIDE] rate     : $RATE DTT/day"
echo ">>> [$SIDE] outdir   : $OUTBASE"
echo ">>> STEP 1  bestmatch"
python3 bestmatch.py
echo ">>> STEP 2  finalize dated tree + QC"
python3 finalize_dated.py
echo ">>> [$SIDE] DONE -> $OUTBASE"
