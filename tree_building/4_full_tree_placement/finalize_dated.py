#!/usr/bin/env python3
"""STEP 2: build the DATED placed tree + placement QC.

Consumes STEP-1 (bestmatch.py) checkpoints and produces:
  * a dated per-side tree: backbone (days) + each placed query grafted with a
    time-scaled pendant  pendant_days = d(C,P) / clock_rate  (dated_tree.emit_dated);
  * place_qc.tsv: one row per placed query with the matched cell, attach score, the
    DTT pendant d(C,P) and residual d(T,P), the pendant in days, the irreversibility
    (is_valid_ancestor) flag, and whether the pendant was clamped to the anchor edge.

Reuses build_tree's checkpoint loader and parent resolver unchanged, and adds the
branch-length semantics (matrix-metric pendant, clock scaling), the validity
diagnostic, and the dated grafting.

Env (via config.py): SIDE, OUTDIR, CONSENSUS_TSV, BACKBONE_NWK (a DATED subtree, days),
N_SITES, ALLCELLS. Plus:
  RATE             -- per-side clock rate (clockrate.py)
  DATED_TREE_OUT   (default <HERE>/<SIDE>/dated_placed_<SIDE>.nwk)
  QC_OUT           (default <HERE>/<SIDE>/place_qc.tsv)
"""
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import clockrate
import dtt_lengths
import validity
from build_tree import load_bestmatch, resolve_parents
from dated_tree import emit_dated
from newick import tips_of


def main():
    t0 = time.time()
    here = os.path.dirname(os.path.abspath(__file__))
    outdir_default = os.path.join(here, C.SIDE)
    dated_out = os.environ.get(
        "DATED_TREE_OUT", os.path.join(outdir_default, f"dated_placed_{C.SIDE}.nwk"))
    qc_out = os.environ.get("QC_OUT", os.path.join(outdir_default, "place_qc.tsv"))
    os.makedirs(os.path.dirname(dated_out), exist_ok=True)
    os.makedirs(os.path.dirname(qc_out), exist_ok=True)

    cells_ids = np.load(f"{C.OUTDIR}/cellids_{C.SIDE}.npy", allow_pickle=True)
    isbb = np.load(f"{C.OUTDIR}/isbb_{C.SIDE}.npy")
    N = len(cells_ids)
    query_idx = np.where(~isbb)[0]

    gmatch, gscore, bbmatch, bbscore = load_bestmatch(C.OUTDIR, C.SIDE, N)
    parent, held_set = resolve_parents(query_idx, isbb, gmatch, bbmatch, bbscore)
    print(f"[finalize] {N:,} cells, {len(query_idx):,} queries, held {len(held_set):,}  "
          f"{time.time()-t0:.1f}s", file=sys.stderr)

    # genotypes aligned to the checkpoint order: loaded ONCE, shared by pendant + validity
    code, integ_cols = dtt_lengths._aligned_code(
        C.CONSENSUS_TSV, C.N_SITES, C.ALLCELLS, tips_of(C.BACKBONE_NWK), cells_ids)
    founder_len = dtt_lengths.load_founder_len(C.FOUNDER_TABLE, C.SIDE, integ_cols)

    pendant, residual, new_edits, shared_edits, n_shared = \
        dtt_lengths.pendant_components(code, parent, founder_len)   # d(C,P), d(T,P), + integer counts
    rate, rate_src = clockrate.get_rate(C.SIDE)
    pendant_days = pendant / rate
    print(f"[finalize] clock rate {C.SIDE}: {rate:.6f} DTT/day  (source: {rate_src})",
          file=sys.stderr)

    valid = validity.valid_flags(code, parent)

    nwk, placed, max_children, clamped = emit_dated(
        cells_ids, isbb, parent, held_set, query_idx, C.BACKBONE_NWK, pendant_days)
    with open(dated_out, "w") as fh:
        fh.write(nwk)

    # ---- placement QC table (one row per placed query) ----
    # attach_score = -d(C,T) at the chosen anchor (higher is always better).
    attach_score = np.where(parent == gmatch, gscore,
                            np.where(parent == bbmatch, bbscore, gscore))
    n_invalid = 0
    with open(qc_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cell_id", "matched_id", "match_is_query", "attach_score",
                    "n_shared_tapes", "shared_edits", "new_edits",
                    "d_C_P_dtt", "d_T_P_dtt", "pendant_days", "valid_attachment", "clamped"])
        for qi in query_idx:
            if qi in held_set or parent[qi] < 0:
                continue
            p = parent[qi]
            v = bool(valid[qi])
            n_invalid += (0 if v else 1)
            w.writerow([cells_ids[qi], cells_ids[p], (not bool(isbb[p])),
                        f"{float(attach_score[qi]):.4f}",
                        int(n_shared[qi]), int(shared_edits[qi]), int(new_edits[qi]),
                        f"{pendant[qi]:.4f}", f"{residual[qi]:.4f}",
                        f"{pendant_days[qi]:.4f}", v, (qi in clamped)])

    print(f"[finalize] {C.SIDE} DONE  ({time.time()-t0:.1f}s)")
    print(f"[finalize] SUMMARY: queries {len(query_idx):,}  placed {placed:,}  "
          f"held {len(held_set):,}  max-direct-children {max_children}  "
          f"invalid(irreversibility) {n_invalid:,}  clamped {len(clamped):,}")
    print(f"[finalize]   dated tree -> {dated_out}")
    print(f"[finalize]   place QC   -> {qc_out}")


if __name__ == "__main__":
    main()
