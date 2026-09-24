#!/usr/bin/env python3
"""Dump one bestmatch.py checkpoint (OUTTAG=0) to a plain TSV for R to read.

Usage: dump_bestmatch_tsv.py OUTDIR SIDE OUT_TSV
Columns: cell_id, matched_id (best BACKBONE match; empty if unplaceable), mb_score.
"""
import sys

import numpy as np

outdir, side, out_tsv = sys.argv[1:4]
cellids = np.load(f"{outdir}/cellids_{side}.npy", allow_pickle=True)
z = np.load(f"{outdir}/bestmatch_0.npz")
q, mb, sb = z["q"], z["mb"], z["sb"]

with open(out_tsv, "w") as fh:
    fh.write("cell_id\tmatched_id\tmb_score\n")
    for qi, mi, s in zip(q, mb, sb):
        matched = cellids[mi] if mi >= 0 else ""
        fh.write(f"{cellids[qi]}\t{matched}\t{s:.6f}\n")
