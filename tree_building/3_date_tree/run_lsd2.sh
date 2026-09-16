#!/usr/bin/env bash
# run_lsd2.sh -- time-scale a rooted, non-negative-branch NJ tree with LSD2
# (least-squares dating). Takes prep_tree.R's cleaned newick and writes a
# dated (time-scaled) tree plus the LSD2 .result summary alongside it.
#
#   usage: 3_date_tree/run_lsd2.sh <tree_nonneg.nwk>
#
# Two calibration dates fix the timescale (root date, tips date); LSD2
# estimates the substitution rate and every internal node's date from those.
# Temporal constraints (node date <= descendant dates) are LSD2's default and
# left on. The input tree is already rooted on the founder-genotype ancestor,
# so we do NOT pass -r (which would re-root).
#
# Config via LSD_* environment variables (defaults below):
#   LSD_BIN     : lsd2 binary                (default lsd2/src/lsd2)
#   LSD_SEQLEN  : -s sequence length         (default 66 = 11 tapes x 6 sites)
#   LSD_ROOT    : -a root date               (default 1.5, blastomere founder)
#   LSD_TIPS    : -z tips date               (default 13.5, E13.5)
#   LSD_VAR     : -v variance mode           (default 1, variance from input branch lengths)
#   LSD_MINBLEN : -u minimum time-scaled branch length, in days (default 0, LSD2's own default)
#   LSD_MINEXBLEN : -U minimum EXTERNAL (tip) branch length, in days (default: same as LSD_MINBLEN)
#   LSD_NULLBLEN  : -l null-branch-length threshold for the QP solver (default: LSD2's own default)
#   LSD_DATEFILE  : -d node minimum-age constraints file (default: none -- unconstrained dating).
#                   Used by run_constrained_dating.sh to pass build_min_age_datefile.R's output.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TREE="${1:?usage: run_lsd2.sh <tree_nonneg.nwk>}"
BIN="${LSD_BIN:-lsd2/src/lsd2}"
SEQLEN="${LSD_SEQLEN:-66}"
ROOT="${LSD_ROOT:-1.5}"
TIPS="${LSD_TIPS:-13.5}"
VAR="${LSD_VAR:-1}"
MINBLEN="${LSD_MINBLEN:-0}"
MINEXBLEN="${LSD_MINEXBLEN:-$MINBLEN}"
NULLBLEN="${LSD_NULLBLEN:-}"
DATEFILE="${LSD_DATEFILE:-}"

OUT="$(dirname "$TREE")/dated"   # -> <tag>/dated.{nexus,date.nexus,nwk,...} + dated.result

lflag=(); [ -n "$NULLBLEN" ] && lflag=( -l "$NULLBLEN" )
dflag=(); [ -n "$DATEFILE" ] && dflag=( -d "$DATEFILE" )

echo "[lsd2] binary : $BIN"
"$BIN" -h 2>&1 | sed $'s/\x1b\\[[0-9;]*m//g' | head -1   # version banner
echo "[lsd2] input  : $TREE"
echo "[lsd2] dates  : root=$ROOT  tips=$TIPS   -s $SEQLEN  -v $VAR  -u $MINBLEN  -U $MINEXBLEN  -l ${NULLBLEN:-(default)}  -d ${DATEFILE:-(none)}"

"$BIN" -i "$TREE" -s "$SEQLEN" -a "$ROOT" -z "$TIPS" -v "$VAR" -u "$MINBLEN" -U "$MINEXBLEN" \
    ${dflag[@]+"${dflag[@]}"} ${lflag[@]+"${lflag[@]}"} -o "$OUT"

echo "[lsd2] done. outputs:"
ls -1 "${OUT}"* 2>/dev/null | sed 's/^/  /'
# LSD2 writes the run summary to the bare output-prefix file ($OUT), the day-scaled
# dated tree to $OUT.date.nexus, and a substitution-scaled tree to $OUT.nwk.
echo "[lsd2] --- ${OUT} (rate / tMRCA / objective) ---"
grep -iE "rate|tMRCA|objective function|Collapse" "${OUT}" 2>/dev/null || tail -20 "${OUT}"

# emit a clean day-unit newick + validate (tips at sampling day, root at root date).
LSD_TIPS="$TIPS" Rscript "$SCRIPT_DIR/time_tree.R" "$(dirname "$TREE")"
