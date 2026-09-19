#!/bin/bash
# Full-dataset run of the alternative (consensus-before-parse) pipeline.
set -e
PY=/Users/jay.shendure/Dropbox/claude/current/final_push/.venv/bin/python
D=/Users/jay.shendure/Dropbox/claude/current/final_push/comment8_reparse
OUT=/Users/shendure/tape_raw/c8sc_full
cd "$D"
echo "[$(date +%H:%M:%S)] STAGE 1 shard (all cells, 256 shards)"
$PY sc1_shard.py --sub 1 --shards 256 --out "$OUT"
echo "[$(date +%H:%M:%S)] shards: $(du -sh $OUT | cut -f1)"
echo "[$(date +%H:%M:%S)] STAGE 2 per-UMI consensus (6 workers)"
$PY sc2p_consensus.py "$OUT" 8
echo "[$(date +%H:%M:%S)] STAGE 4 cross-UMI tape consensus"
$PY sc4_tape.py "$OUT"
echo "[$(date +%H:%M:%S)] DONE"
