#!/usr/bin/env python3
"""Stage 00 (single-cell, paired only) -- merge R1/R2 mates.

Usage:  python3 scripts/00_merge_pairs.py seq3

Overlap-merges paired mates with fastp and writes amplicon-oriented reads
(merged where possible, R2 fallback otherwise) to data/<sample>.sc_reads.fastq.gz.
Single-end samples (e.g. seq1) skip this stage -- their R2 is used directly.
Requires fastp on PATH (see environment.yml) and paired-ordered R1/R2 fastqs.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.merge import merge_sample


def main(name: str):
    s = config.get_sc_sample(name)
    if s.layout != "paired":
        sys.exit(f"sample {name!r} is {s.layout}; no merge needed (use R2 directly).")
    for p in (s.r1_path, s.r2_path):
        if not p.exists():
            sys.exit(f"missing input: {p}")
    combined = s.out("sc_reads.fastq.gz")
    workdir = config.DATA_DIR / f"{name}.fastp"
    print(f"[00] merging {name}: {s.r1} + {s.r2}  (fastp)")
    stats = merge_sample(s.r1_path, s.r2_path, workdir, combined)
    print(f"[00] merged {stats['merged']:,} pairs ({stats['merge_rate']*100:.1f}%); "
          f"R2 fallback {stats['r2_fallback']:,}; total reads {stats['total']:,}")
    print(f"[00] wrote {combined}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
