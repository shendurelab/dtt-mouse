#!/usr/bin/env python3
"""Stage 01 -- derive the TAPE-BC whitelist + NNN vocabulary for an embryo.

Usage:  python3 scripts/01_derive_reference.py DTTz_3_S3

Samples config.SAMPLE_N reads from the embryo's fastq and writes
refs/<tag>.tapebc_whitelist.tsv and refs/<tag>.nnn_vocab.tsv. Run once per
library; the parse stage reuses these. Requires the raw fastq under RAW_DIR.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.reference import derive, write_whitelist, write_vocab


def main(tag: str):
    emb = config.get_embryo(tag)
    if not emb.fastq_path.exists():
        sys.exit(f"raw fastq not found: {emb.fastq_path}\n"
                 f"Restage it under {config.RAW_DIR} (see the README section 'Raw data').")
    print(f"[01] deriving references for {emb.label} ({tag}) from {emb.fastq_path.name}")
    whitelist, vocab, bc_counts, ins_peak = derive(emb.fastq_path)
    tot_bc = sum(bc_counts.values())
    write_whitelist(emb.whitelist_tsv, whitelist, tot_bc)
    write_vocab(emb.vocab_tsv, ins_peak, vocab)
    print(f"[01] {len(whitelist)} whitelist barcodes | {len(vocab)} vocab insertions "
          f"(per-(BC x position) peak >= {config.VOCAB_PEAK_MIN:,})")
    print(f"[01] wrote {emb.whitelist_tsv}")
    print(f"[01] wrote {emb.vocab_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
