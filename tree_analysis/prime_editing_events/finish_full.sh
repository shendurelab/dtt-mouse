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
PY="${PY:-python3}"
D="${C8R_DIR:?set C8R_DIR to the staging pipeline directory}"
until grep -q "DONE" "$D/run_full.log"; do sleep 20; done
cd "$D"
echo "=== FULL DATASET ==="
$PY sc11_events.py ${DTT_TAPE_CALLS:?set DTT_TAPE_CALLS} \
    "$D/comment8_event_classes_FULL.csv"
echo; echo "=== cells / support / array length ==="
$PY - <<'PY'
import sys, os, gzip, statistics as st
from collections import Counter
sys.path.insert(0,os.environ["DTT_SC_TAPE"]); sys.path.insert(0,".")
os.environ.setdefault("TAPE_RAW_DIR",os.environ.get("TAPE_RAW_DIR",""))
import config
from tape.io import load_vocab
from sc11_events import walk
vocab=load_vocab(config.get_embryo("DTTz_3_S3").vocab_tsv)
cells=Counter(); nmol=[]; L=Counter(); exp=0; noterm=0; tot=0
with gzip.open("${DTT_TAPE_CALLS:?set DTT_TAPE_CALLS}","rt") as f:
    next(f)
    for line in f:
        p=line.rstrip("\n").split("\t"); tot+=1
        cells[p[0]]+=1; nmol.append(int(p[2]))
        prof,anch=walk(p[9].encode(),vocab)
        kinds=[k for k,_ in prof]
        if any(k not in ("unedited","edit_vocab") for k in kinds) or anch: continue
        seen=False; viol=False
        for k in kinds:
            if k=="unedited": seen=True
            elif k.startswith("edit") and seen: viol=True; break
        if viol: continue
        exp+=1
        if p[4]=="None": noterm+=1
        else: L[int(p[4])]+=1
print(f"tape calls {tot:,}   cells with >=1 call {len(cells):,}   tapes/cell median {st.median(list(cells.values())):.0f} mean {st.mean(list(cells.values())):.2f}")
print(f"molecules per call: median {st.median(nmol):.0f} mean {st.mean(nmol):.2f}")
print(f"expected calls {exp:,}   terminal located {sum(L.values()):,}   no terminal {noterm:,}")
t=sum(L.values())
print("array length: " + "  ".join(f"L{k}:{100*L[k]/t:.2f}%" for k in sorted(L) if 3<=k<=8))
PY
