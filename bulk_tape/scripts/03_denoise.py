#!/usr/bin/env python3
"""Stage 03 -- de-noise per-barcode patterns.

Usage:  python3 scripts/03_denoise.py DTTz_3_S3

Reads data/<tag>.patterns_full.pkl and applies the canonical funnel
(collapse -> chimera filter -> >=MIN_READS -> cross-barcode filter). Writes:
  data/<tag>.clean_patterns.pkl   {barcode: {pattern: count}} retained
  data/<tag>.denoise_funnel.tsv   per-barcode stage counts (QC)
Only barcodes with >= config.MIN_BC_READS total reads are analyzed.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import load_pickle, save_pickle
from tape.denoise import denoise_all


def main(tag: str):
    emb = config.get_embryo(tag)
    src = emb.out("patterns_full.pkl")
    if not src.exists():
        sys.exit(f"missing {src}; run stage 02 first.")
    full = load_pickle(src)
    bcs = {b: c for b, c in full.items() if sum(c.values()) >= config.MIN_BC_READS}
    print(f"[03] de-noising {emb.label}: {len(bcs)} barcodes >= {config.MIN_BC_READS} reads "
          f"({len(full)} total)")

    clean, funnels, _ = denoise_all(bcs)
    save_pickle(clean, emb.out("clean_patterns.pkl"))

    cols = ["raw", "collapsed", "chimera_filtered", "ge_min_reads", "cross_bc", "clean"]
    with open(emb.out("denoise_funnel.tsv"), "w") as f:
        f.write("tapebc\ttotal_reads\t" + "\t".join(cols) + "\n")
        for b in sorted(bcs, key=lambda b: -sum(bcs[b].values())):
            fn = funnels[b]
            f.write(f"{b}\t{sum(bcs[b].values())}\t" + "\t".join(str(fn[c]) for c in cols) + "\n")

    total_clean = sum(len(v) for v in clean.values())
    print(f"[03] {total_clean:,} clean lineages across {len(clean)} barcodes")
    print(f"[03] wrote {emb.out('clean_patterns.pkl')} + {emb.out('denoise_funnel.tsv')}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
