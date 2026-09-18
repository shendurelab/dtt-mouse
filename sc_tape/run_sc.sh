#!/usr/bin/env bash
# Full single-cell circtape pipeline, embryo #3, from raw fastqs to the two
# per-blastomere tree-input matrices. Needs the raw data (see README); every
# stage after 00 is deterministic and re-runnable.
set -euo pipefail
cd "$(dirname "$0")"
: "${TAPE_RAW_DIR:?set TAPE_RAW_DIR to the directory holding the fastqs}"
: "${TAPE_DATA_DIR:=./data}"
export TAPE_RAW_DIR TAPE_DATA_DIR
PY=${PY:-python3}

# 1. overlap-merge the paired run (needs fastp); the single-end run needs nothing
$PY scripts/00_merge_pairs.py seq8_long

# 2. genotype: parse, UMI-collapse, de-noise, per-cell consensus, both runs in ONE pass
$PY scripts/sc_consensus_par.py seq8 seq8_long --out e3v8 --procs 8

# 3. blastomere routing -> e3v8.B1/B2_tape_consensus.tsv (+ set-aside, labels, report)
$PY scripts/blastomere_route.py "$TAPE_DATA_DIR/e3v8.cell_consensus.tsv" \
      --outdir "$TAPE_DATA_DIR" --out e3v8 \
      --min-n-loci 4 --max-doublet-loci 1 --min-mean-dom 0.95

# 4. reported summaries (need cell_metadata.v8.txt under TAPE_RAW_DIR)
$PY scripts/sc_recovery.py e3v8 || echo "skipped: needs cell_metadata.v8.txt"
