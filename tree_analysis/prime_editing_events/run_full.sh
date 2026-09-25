#!/bin/bash
# PROVENANCE ONLY -- this is the driver used for the original full-dataset run.
# It calls a sharding/consensus staging pipeline (sc1_shard.py, sc2p_consensus.py,
# sc4_tape.py, sc11_events.py) that is NOT part of this release; the analysis that
# the paper reports is 01_reparse_consensus_first.py .. 04_junction_analysis.py in
# this directory, which run from what ships in support_data/.
# Kept so the exact stages and ordering of the original run are on record.
#
# To re-run it you need that staging pipeline plus the raw reads (GEO GSE341627):
#   PY=/path/to/python  C8R_DIR=/path/to/staging  TAPE_RAW_DIR=/path/to/raw  ./run_full.sh
# Full-dataset run of the alternative (consensus-before-parse) pipeline.
set -e
PY="${PY:-python3}"
D="${C8R_DIR:?set C8R_DIR to the staging pipeline directory}"
OUT="${TAPE_RAW_DIR:?set TAPE_RAW_DIR}/c8sc_full"
cd "$D"
echo "[$(date +%H:%M:%S)] STAGE 1 shard (all cells, 256 shards)"
$PY sc1_shard.py --sub 1 --shards 256 --out "$OUT"
echo "[$(date +%H:%M:%S)] shards: $(du -sh $OUT | cut -f1)"
echo "[$(date +%H:%M:%S)] STAGE 2 per-UMI consensus (6 workers)"
$PY sc2p_consensus.py "$OUT" 8
echo "[$(date +%H:%M:%S)] STAGE 4 cross-UMI tape consensus"
$PY sc4_tape.py "$OUT"
echo "[$(date +%H:%M:%S)] DONE"
