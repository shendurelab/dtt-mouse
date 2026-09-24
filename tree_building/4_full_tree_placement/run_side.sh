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
#   LEFT_NWK=results/4-full-tree/B1/dated_placed_B1.nwk \
#   RIGHT_NWK=results/4-full-tree/B2/dated_placed_B2.nwk \
#   OUT_NWK=results/4-full-tree/merged_placed_check.nwk \
#     Rscript 3_date_tree/merge_dated_subtrees.R
set -euo pipefail
cd "$(dirname "$0")"
TB="$(cd .. && pwd)"          # tree_building/
REPO="$(cd ../.. && pwd)"     # repo root (support_data/ lives here)

SIDE="${1:?usage: run_side.sh B1|B2}"
export SIDE
export ALLCELLS="${ALLCELLS:-0}"        # 0 = pass_qc cells only (queries = pass_qc, <7 loci etc.)
export N_SITES="${N_SITES:-6}"

# ---- the dated backbone (3_date_tree's constrained, lineage-constrained min-age
# dating output for e3v8/ge7_founderok -- see tree_building/results/3-dated-tree) ---
VARIANT="${VARIANT:-minB2h_lineage_constrained}"
export BACKBONE_NWK="${BACKBONE_NWK:-$TB/results/3-dated-tree/perside_${SIDE}_${VARIANT}.nwk}"
[ -f "$BACKBONE_NWK" ] || { echo "backbone not found: $BACKBONE_NWK" >&2; exit 1; }

# ---- clock rate: LSD2's own fitted rate for the v8_ge7 constrained dating run of
# this side (mouse_sprint results/4-time-tree/v8_ge7/synthroot/$SIDE/lsd2/nj99478/
# attempt1_minage_sourced_minB2h_l0.01/dated -- "rate X, tMRCA" line, iter3/final;
# also restated verbatim in that repo's CURRENT_TREE.txt). dtt-mouse is v8-only, so
# these ARE the target rates (not a stale default); override with env RATE to test
# another value.
case "$SIDE" in
  B1) export RATE="${RATE:-0.251985}" ;;   # DTT/day, v8_ge7 constrained LSD2 fit
  B2) export RATE="${RATE:-0.303488}" ;;   # DTT/day, v8_ge7 constrained LSD2 fit
  *)  echo "unknown side: $SIDE" >&2; exit 1 ;;
esac

# ---- consensus (queries): founder-consistent cells of ANY n_loci (not just the
# ge7-founder-ok subset the backbone was built from -- placement's query universe
# is every QC-pass, founder-consistent cell). The manuscript Methods require
# founder-consistency "at both stages" (backbone AND placement queries), so this
# must NOT be the raw, unfiltered support_data consensus -- using that let 15,520
# founder-divergent cells get placed anyway (verified against the shipped tree).
# Derived by 1_build_nj_backbone/00_filter_highqual_consensus.R MIN_LOCI=4.
export CONSENSUS_TSV="${CONSENSUS_TSV:-$TB/processed_data/e3v8.${SIDE}_tape_consensus.founderok.tsv.gz}"
export FOUNDER_TABLE="${FOUNDER_TABLE:-$TB/processed_data/e3v8.supp_table1_founder_genotypes.wide.csv}"

# ---- output ---------------------------------------------------------------------
OUTBASE="$TB/results/4-full-tree/${SIDE}"
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
