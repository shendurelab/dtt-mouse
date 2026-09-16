#!/bin/bash
# run_dating.sh -- per-blastomere tree dating: BAT neg-branch correction ->
# LSD2 dating -> merge into one whole-embryo time tree.
#
# Takes 2_root_tree's rooted, outgroup-dropped newick for each side (B1/B2)
# and produces one merged, whole-embryo time tree with branch lengths in
# days. B1/B2 are independent: BAT correction runs SERIALLY (memory-heavy in
# R); LSD2 dating runs B1 || B2 in PARALLEL (memory-light C++). Safe to
# interrupt/rerun: every step is skipped if its output already exists.
#
# REQUIRED env vars:
#   ROOTED_TREE_B1 / ROOTED_TREE_B2 : 2_root_tree's output for each side
#   STAGEDIR                        : output namespace (e.g. results/3-dated-tree)
#   MERGE_TOKEN                     : token in the merged filename (e.g. ge7)
# OPTIONAL env vars (defaults shown):
#   VARIANTS=minB2h        space-sep LSD2 min-branch variants (minB2h | noMinB)
#   NULLBLEN_B1=0.10       LSD2 -l for B1 (collapse near-zero NJ tie-break branches
#                          so the QP solver stays tractable on very large trees)
#   NULLBLEN_B2=(empty)    LSD2 -l for B2 (empty -> LSD2 default)
#   LSD_ROOT_AGE=1.5       per-side root date (days); blastomere founder = 2-cell
#   LSD_TIPS_AGE=13.5      tip date (days); E13.5
#   LSD_BIN=lsd2/src/lsd2
#
# Example:
#   ROOTED_TREE_B1=results/2-rooted-nj/B1/nj_rooted_ingroup.nwk \
#   ROOTED_TREE_B2=results/2-rooted-nj/B2/nj_rooted_ingroup.nwk \
#   STAGEDIR=results/3-dated-tree MERGE_TOKEN=ge7 \
#     bash 3_date_tree/run_dating.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$REPO_DIR"

: "${ROOTED_TREE_B1:?set ROOTED_TREE_B1 to 2_root_tree's B1 output}"
: "${ROOTED_TREE_B2:?set ROOTED_TREE_B2 to 2_root_tree's B2 output}"
: "${STAGEDIR:?set STAGEDIR to the output namespace dir}"
: "${MERGE_TOKEN:?set MERGE_TOKEN (merged-filename token, e.g. ge7)}"

LSD_BIN="${LSD_BIN:-$REPO_DIR/lsd2/src/lsd2}"
LOGDIR="$REPO_DIR/logs"; mkdir -p "$LOGDIR"

LSD_ROOT_AGE="${LSD_ROOT_AGE:-1.5}"
LSD_TIPS_AGE="${LSD_TIPS_AGE:-13.5}"
VARIANTS="${VARIANTS:-minB2h}"
minblen_for() { case "$1" in noMinB) echo "0";; minB2h) echo "0.083333";;  # 2h/24
    *) echo "Unknown dating variant: $1" >&2; exit 1;; esac; }
nullblen_for() { case "$1" in B1) echo "${NULLBLEN_B1:-0.10}";; B2) echo "${NULLBLEN_B2:-}";; esac; }
rooted_tree_for() { case "$1" in B1) echo "$ROOTED_TREE_B1";; B2) echo "$ROOTED_TREE_B2";; esac; }

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---- 0. preflight -------------------------------------------------------
for side in B1 B2; do
    f=$(rooted_tree_for "$side")
    [ -f "$f" ] || { echo "Missing rooted tree ($side): $f"; exit 1; }
done
[ -x "$LSD_BIN" ] || { echo "lsd2 not found at $LSD_BIN"; exit 1; }
Rscript -e 'library(ape); library(BAT)' \
    || { echo "R packages not loadable (ape, BAT)"; exit 1; }

# prep_tree.R tags its output subdir by tip count (>50k -> "nj99478").
tag_for_rooted_tree() {
    Rscript -e "n <- ape::Ntip(ape::read.tree('$1')); cat(if (n > 50000) 'nj99478' else sprintf('sub%dk', round(n / 1000)))"
}

# prep_side = BAT neg-branch correction (memory-heavy -> SERIAL).
prep_side() {
    local side=$1
    local rooted; rooted=$(rooted_tree_for "$side")
    local lsd_outdir="${STAGEDIR}/${side}/lsd2"
    local tag; tag=$(tag_for_rooted_tree "$rooted")
    local prepped="${lsd_outdir}/${tag}/tree_nonneg.nwk"
    if [ -f "$prepped" ]; then
        log "BAT-corrected tree ($side) already exists, skipping prep."
    else
        log "Prepping ($side)..."
        LSD_TREE="$rooted" LSD_OUTDIR="$lsd_outdir" \
            Rscript 3_date_tree/prep_tree.R 2>&1 | tee "$LOGDIR/prep_${MERGE_TOKEN}_${side}.log"
    fi
}

# lsd_side = LSD2 dating over the min-branch variants (memory-light -> PARALLEL).
lsd_side() {
    local side=$1
    local rooted; rooted=$(rooted_tree_for "$side")
    local lsd_outdir="${STAGEDIR}/${side}/lsd2"
    local tag; tag=$(tag_for_rooted_tree "$rooted")
    local prepped="${lsd_outdir}/${tag}/tree_nonneg.nwk"
    local variant
    for variant in $VARIANTS; do
        local vdir="${lsd_outdir}/${tag}/${variant}"
        local time_tree="${vdir}/time_tree.nwk"
        if [ -f "$time_tree" ]; then log "dated tree ($side, $variant) already exists, skipping."; continue; fi
        mkdir -p "$vdir"; cp "$prepped" "${vdir}/tree_nonneg.nwk"
        local minb; minb=$(minblen_for "$variant")
        local nullb; nullb=$(nullblen_for "$side")
        log "Dating ($side, $variant): LSD_ROOT=$LSD_ROOT_AGE, LSD_MINBLEN=$minb, LSD_NULLBLEN=${nullb:-(default)}"
        LSD_BIN="$LSD_BIN" LSD_ROOT="$LSD_ROOT_AGE" LSD_TIPS="$LSD_TIPS_AGE" \
        LSD_MINBLEN="$minb" LSD_MINEXBLEN="$minb" LSD_NULLBLEN="$nullb" \
            3_date_tree/run_lsd2.sh "${vdir}/tree_nonneg.nwk" \
            2>&1 | tee "$LOGDIR/lsd2_${MERGE_TOKEN}_${side}_${variant}.log"
    done
}

# ---- 1a. BAT prep, SERIAL (memory-heavy) --------------------------------
for side in B1 B2; do prep_side "$side"; done
# ---- 1b. LSD2 dating, B1 || B2 PARALLEL (memory-light) ------------------
for side in B1 B2; do lsd_side "$side" & done
wait

# ---- 2. merge B1 + B2 into one whole-embryo tree, per variant -----------
TAG_B1=$(tag_for_rooted_tree "$ROOTED_TREE_B1")
TAG_B2=$(tag_for_rooted_tree "$ROOTED_TREE_B2")
for variant in $VARIANTS; do
    MERGED="${STAGEDIR}/merged/merged_time_tree_${MERGE_TOKEN}_${variant}.nwk"
    if [ -f "$MERGED" ]; then log "merged tree ($variant) already exists, skipping."; continue; fi
    log "Merging B1 + B2 ($variant) (MRCA_AGE=0, ROOT_AGE=$LSD_ROOT_AGE)..."
    mkdir -p "$(dirname "$MERGED")"
    LEFT_NWK="${STAGEDIR}/B1/lsd2/${TAG_B1}/${variant}/time_tree.nwk" \
    RIGHT_NWK="${STAGEDIR}/B2/lsd2/${TAG_B2}/${variant}/time_tree.nwk" \
    OUT_NWK="$MERGED" MRCA_AGE=0 ROOT_AGE="$LSD_ROOT_AGE" \
        Rscript 3_date_tree/merge_dated_subtrees.R 2>&1 | tee "$LOGDIR/merge_${MERGE_TOKEN}_${variant}.log"
done

log "Done. Merged whole-embryo time trees (${MERGE_TOKEN}):"
for variant in $VARIANTS; do
    log "  $REPO_DIR/${STAGEDIR}/merged/merged_time_tree_${MERGE_TOKEN}_${variant}.nwk"
done
