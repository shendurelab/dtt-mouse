#!/usr/bin/env python3
"""STEP 1: compute each query cell's DTT-closest genotype match.

Loads the consensus callset, encodes genotypes, builds the rarity-capped inverted
index over ALL cells (backbone anchors + queries), then finds, for every query, its
single DTT-closest cell overall and its DTT-closest backbone neighbour (dtt_match.py).
Results are checkpointed so finalize_dated.py can run -- and be re-run -- without
rescoring.

Queries = the non-backbone cells (backbone tips are fixed anchors, never moved).
The search space (anchors) is the FULL set: a query may attach to another query,
which is how genuine query-query relatives are kept.

Parallelism: the query list is sliced [START:END]; run several workers with
different START/END/OUTTAG and they each drop a bestmatch_<OUTTAG>.npz. The first
worker (OUTTAG=="0") also writes the shared cell ordering so all chunks agree on
global indices. finalize_dated.py globs all bestmatch_*.npz back together.

Env (all via config.py): SIDE, ALLCELLS, CAP, START, END, OUTTAG, OUTDIR,
CONSENSUS_TSV, BACKBONE_NWK. Plus TOPK (candidate-shortlist size, default 500).

Normally run via run_side.sh, which also sets CONSENSUS_TSV and BACKBONE_NWK (both
required -- no default). Direct run:
    SIDE=B1 CONSENSUS_TSV=... BACKBONE_NWK=... python3 bestmatch.py            # one worker
    SIDE=B1 CONSENSUS_TSV=... BACKBONE_NWK=... START=0 END=150000 OUTTAG=0 \\
        python3 bestmatch.py                                                  # chunk
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from genotypes import load_and_encode
from newick import tips_of
from similarity import token_freq, build_index
from dtt_match import best_matches_dtt


def main():
    os.makedirs(C.OUTDIR, exist_ok=True)

    backbone_tips = tips_of(C.BACKBONE_NWK)
    cells = load_and_encode(C.CONSENSUS_TSV, backbone_tips, C.N_SITES, C.ALLCELLS)

    # The canonical cell ordering is shared by every chunk; only worker 0 writes it.
    if C.OUTTAG == "0" and not os.path.exists(f"{C.OUTDIR}/cellids_{C.SIDE}.npy"):
        np.save(f"{C.OUTDIR}/cellids_{C.SIDE}.npy", cells.ids)
        np.save(f"{C.OUTDIR}/isbb_{C.SIDE}.npy", cells.isbb)
        np.save(f"{C.OUTDIR}/nloci_{C.SIDE}.npy", cells.nloci)

    freq = token_freq(cells.code, cells.G, C.N_SITES, C.TOK_STRIDE)
    postings = build_index(cells.code, cells.G, C.N_SITES, C.TOK_STRIDE, freq, C.CAP_FRAC)

    all_queries = np.where(~cells.isbb)[0]
    end = len(all_queries) if C.END < 0 else C.END
    q_slice = all_queries[C.START:end]

    m, s, mb, sb = best_matches_dtt(
        cells.code, cells.G, C.N_SITES, C.TOK_STRIDE, postings,
        q_slice, cells.isbb, cells.nloci,
        topk=int(os.environ.get("TOPK", "500")))

    out = f"{C.OUTDIR}/bestmatch_{C.OUTTAG}.npz"
    # q = global query index; m/s = best-overall match+score; mb/sb = best-backbone.
    np.savez(out, q=q_slice.astype(np.int64), m=m, s=s, mb=mb, sb=sb)
    print(f"[bestmatch] wrote {out}  ({len(q_slice):,} queries)", file=sys.stderr)


if __name__ == "__main__":
    main()
