#!/usr/bin/env bash
# Run the full bulk TAPE pipeline for one embryo, stages 01 -> 06.
#
#   ./run_bulk.sh DTTz_3_S3          # full run (needs raw fastq in data/raw/)
#   ./run_bulk.sh DTTz_3_S3 --skip-ref   # reuse checked-in refs/, skip stage 01
#
# Stages 01 and 02 require the raw fastq under data/raw/ (or $TAPE_RAW_DIR); the
# raw data is not distributed with the code -- see README ("Raw data"). Stages
# 03-06 run from data/<tag>.patterns_full.pkl, which IS distributed, under
# example_data/. To run them without any raw data:
#     mkdir -p data && cp example_data/*.patterns_full.pkl data/
set -euo pipefail
cd "$(dirname "$0")"

TAG="${1:?usage: ./run_bulk.sh <embryo_tag> [--skip-ref]}"
SKIP_REF="${2:-}"
PY="${PYTHON:-python3}"

if [[ "$SKIP_REF" != "--skip-ref" ]]; then
  $PY scripts/01_derive_reference.py "$TAG"
else
  echo "[run] reusing checked-in refs/ (skipping stage 01)"
fi
$PY scripts/02_parse.py "$TAG"
$PY scripts/03_denoise.py "$TAG"
$PY scripts/04_depth_model.py "$TAG"
$PY scripts/05_integration_trees.py "$TAG"
$PY scripts/06_rarefaction.py "$TAG"
echo "[run] done: outputs in data/ and data/figures/"
