#!/usr/bin/env bash
# =============================================================================
# run_constrained_dating.sh -- literature min-age ladder (population-size
# constraint) re-dating for per-blastomere trees, run AFTER run_dating.sh.
#
# For each side: build an LSD2 -d min-age datefile from a NODE RANKING (read
# off run_dating.sh's unconstrained dated tree) + a whole-embryo cell-count
# CEILING, then re-run LSD2 on the divergence tree (tree_nonneg.nwk) with -d,
# iterated ITERS times (re-ranking on the constrained tree each pass). Then
# merge B1+B2. This is the fix for the LTT-ceiling overshoot: see
# build_min_age_datefile.R for why plain LSD2 packs impossibly many early
# divisions into a sliver of calendar time near the root.
#
# REQUIRED env:
#   STAGEDIR           per-side dating namespace (contains {B1,B2}/lsd2/nj99478/...,
#                      from run_dating.sh)
#   RANK_TREE_B1       initial (unconstrained) dated tree used ONLY for the iter-1
#   RANK_TREE_B2         node ranking per side (branch lengths in days) -- run_dating.sh's
#                        time_tree.nwk. Paths relative to repo root.
#   MERGE_TOKEN        token in the merged filename (e.g. ge7)
# OPTIONAL env (defaults shown):
#   CEILING=processed_data/sample_matched_ceiling_sourced.csv
#   NULLBLEN_B1=0.01   LSD2 -l for the constrained dating (per side). Default 0.01
#   NULLBLEN_B2=0.01     (the data-justified value); MUST be tractable for the QP.
#   SIDE_FRAC_B1=0.5   fraction of the whole-embryo cell-count ceiling this side's
#   SIDE_FRAC_B2=0.5     tree is capped at (build_min_age_datefile.R's SIDE_FRAC arg).
#                        Default is an even split; set both to an asymmetric pair
#                        (e.g. 0.63/0.37) for a sensitivity run -- point STAGEDIR at
#                        a separate namespace so it doesn't overwrite the primary run.
#   VARIANT=attempt1_minage_sourced_minB2h_l0.01   output subdir + merged token part
#   ITERS=3  ROOT_AGE=1.5  TIPS_AGE=13.5  SEQLEN=66  MINB=0.083333
#   LSD_BIN=lsd2/src/lsd2
#
# Example (ranking from the minB2h unconstrained trees):
#   STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
#   RANK_TREE_B1=results/3-dated-tree/B1/lsd2/nj99478/minB2h/time_tree.nwk \
#   RANK_TREE_B2=results/3-dated-tree/B2/lsd2/nj99478/minB2h/time_tree.nwk \
#     bash 3_date_tree/run_constrained_dating.sh
#
# Safe to interrupt/rerun (LSD2 steps re-run; the final merge is the deliverable).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$REPO_DIR"

: "${STAGEDIR:?set STAGEDIR (per-side dating namespace)}"
: "${MERGE_TOKEN:?set MERGE_TOKEN (merged-filename token, e.g. ge7)}"
# SIDES: which side(s) to date this invocation (default both). Set SIDES=B1 to
# date only B1 (e.g. when B2's ranking isn't ready yet); the merge runs only once
# BOTH sides' constrained time trees exist. RANK_TREE_<side> is required only for
# the side(s) being dated.
SIDES="${SIDES:-B1 B2}"

LSD_BIN="${LSD_BIN:-$REPO_DIR/lsd2/src/lsd2}"
CEILING="${CEILING:-processed_data/sample_matched_ceiling_sourced.csv}"
ROOT_AGE="${ROOT_AGE:-1.5}"; TIPS_AGE="${TIPS_AGE:-13.5}"; SEQLEN="${SEQLEN:-66}"
ITERS="${ITERS:-3}"; MINB="${MINB:-0.083333}"
VARIANT="${VARIANT:-attempt1_minage_sourced_minB2h_l0.01}"
LOGDIR="$REPO_DIR/logs"; mkdir -p "$LOGDIR"

tag_for() { echo nj99478; }   # both v6 sides >50k tips -> nj99478
# Constrained dating nullblen (-l). Default 0.01 both sides (data-justified). An
# empty value -> -l omitted (LSD2 default). Must match what makes the QP tractable.
nullblen_for() { case "$1" in B1) echo "${NULLBLEN_B1-0.01}" ;; B2) echo "${NULLBLEN_B2-0.01}" ;; esac; }
rank_tree_for() { case "$1" in B1) echo "$RANK_TREE_B1" ;; B2) echo "$RANK_TREE_B2" ;; esac; }
side_frac_for() { case "$1" in B1) echo "${SIDE_FRAC_B1:-0.5}" ;; B2) echo "${SIDE_FRAC_B2:-0.5}" ;; esac; }

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---- 0. preflight -----------------------------------------------------------
[ -x "$LSD_BIN" ] || { echo "lsd2 not found at $LSD_BIN"; exit 1; }
[ -f "$CEILING" ] || { echo "ceiling CSV not found: $CEILING"; exit 1; }
for side in $SIDES; do
    nonneg="${STAGEDIR}/${side}/lsd2/$(tag_for "$side")/tree_nonneg.nwk"
    rank=$(rank_tree_for "$side")
    [ -f "$nonneg" ] || { echo "Missing divergence tree: $nonneg"; exit 1; }
    [ -n "$rank" ] && [ -f "$rank" ] || { echo "Missing ranking tree ($side): '${rank:-<unset RANK_TREE_$side>}'"; exit 1; }
done
Rscript -e 'library(ape)' >/dev/null 2>&1 || { echo "ape not loadable"; exit 1; }

# ---- 1. per-side: min-age ladder + constrained LSD2 (ITERS passes) -----------
date_side_constrained() {
    local side=$1
    local base="${STAGEDIR}/${side}/lsd2/$(tag_for "$side")"
    local nonneg="${base}/tree_nonneg.nwk"
    local vdir="${base}/${VARIANT}"
    mkdir -p "$vdir"; cp "$nonneg" "${vdir}/tree_nonneg.nwk"

    local nb; nb=$(nullblen_for "$side")

    local rank_tree; rank_tree=$(rank_tree_for "$side")   # iter-1 ranking = the unconstrained tree
    local iter
    for iter in $(seq 1 "$ITERS"); do
        local datefile="${vdir}/min_age.iter${iter}.datefile"
        log "[$side iter$iter] building min-age ladder from $rank_tree"
        ROOT_DAY="$ROOT_AGE" Rscript "${SCRIPT_DIR}/build_min_age_datefile.R" \
            "$rank_tree" "$CEILING" "$datefile" "$(side_frac_for "$side")" 2>&1 | tee "$LOGDIR/constr_${MERGE_TOKEN}_build_${side}_iter${iter}.log"

        log "[$side iter$iter] LSD2 -d $(basename "$datefile") (minB=$MINB, -l ${nb:-(default)})"
        LSD_BIN="$LSD_BIN" LSD_ROOT="$ROOT_AGE" LSD_TIPS="$TIPS_AGE" LSD_SEQLEN="$SEQLEN" \
        LSD_MINBLEN="$MINB" LSD_MINEXBLEN="$MINB" LSD_NULLBLEN="$nb" LSD_DATEFILE="$datefile" \
            "${SCRIPT_DIR}/run_lsd2.sh" "${vdir}/tree_nonneg.nwk" \
            2>&1 | tee "$LOGDIR/constr_${MERGE_TOKEN}_lsd2_${side}_iter${iter}.log"
        rank_tree="${vdir}/time_tree.nwk"                 # iters 2..N re-rank on the constrained tree
    done
}

for side in $SIDES; do date_side_constrained "$side" & done
wait

# ---- 2. merge B1 + B2 -- ONLY when BOTH sides were dated THIS invocation -----
# Guard against the single-side race: if SIDES=B1 (or B2) alone, the OTHER side's
# time_tree.nwk may exist but be an INTERMEDIATE iteration from a still-running
# run, which would produce a mixed-iteration merge. So only auto-merge when this
# invocation ran both sides. For a split run (one side at a time), do the merge
# explicitly once both sides are final.
both=0; echo "$SIDES" | grep -qw B1 && echo "$SIDES" | grep -qw B2 && both=1
B1_TT="${STAGEDIR}/B1/lsd2/$(tag_for B1)/${VARIANT}/time_tree.nwk"
B2_TT="${STAGEDIR}/B2/lsd2/$(tag_for B2)/${VARIANT}/time_tree.nwk"
if [ "$both" != 1 ] || [ ! -f "$B1_TT" ] || [ ! -f "$B2_TT" ]; then
    log "Merge skipped: this invocation must date BOTH sides (SIDES='$SIDES') AND both constrained time trees must exist (B1:$( [ -f "$B1_TT" ] && echo ok || echo MISSING) B2:$( [ -f "$B2_TT" ] && echo ok || echo MISSING)). For a split run, merge explicitly once both sides are final."
    exit 0
fi
MERGED="${STAGEDIR}/merged/merged_time_tree_${MERGE_TOKEN}_${VARIANT}.nwk"
mkdir -p "$(dirname "$MERGED")"
log "Merging B1 + B2 ($VARIANT) (MRCA_AGE=0, ROOT_AGE=$ROOT_AGE)..."
LEFT_NWK="${STAGEDIR}/B1/lsd2/$(tag_for B1)/${VARIANT}/time_tree.nwk" \
RIGHT_NWK="${STAGEDIR}/B2/lsd2/$(tag_for B2)/${VARIANT}/time_tree.nwk" \
OUT_NWK="$MERGED" MRCA_AGE=0 ROOT_AGE="$ROOT_AGE" \
    Rscript 3_date_tree/merge_dated_subtrees.R 2>&1 | tee "$LOGDIR/constr_${MERGE_TOKEN}_merge.log"

log "Done. Constrained merged time tree: $REPO_DIR/$MERGED"
