#!/usr/bin/env python3
"""Per-barcode read-count rank-abundance data for the integration-barcode
supplementary panel (one facet per embryo): reads per TAPE-BC, rank-ordered,
annotated with whitelist membership + estimated copy number, plus the cutoff.

For each embryo we count extracted 12-bp barcodes over ALL reads, fold error
shadows onto their whitelist parent (post-Hamming correction, as the pipeline
does per read), and report the TOP_N most-abundant post-correction sequences
(the same number per embryo, for even facets) -- rank-ordered and annotated with
whitelist membership and estimated copy number (from stage 04).

Usage:  python3 scripts/barcode_rank_abundance.py DTTz_1_S1 DTTz_2_S2 DTTz_3_S3 DTTz_4_S4

Writes to tables/:
  suppfig1a_barcode_rank_abundance.csv   embryo, rank, tapebc, reads, freq, whitelisted, copies
  suppfig1a_cutoffs.csv                  embryo, n_whitelisted, total_bc_reads,
                                         freq_cutoff, reads_at_freq_cutoff, min_whitelisted_reads
"""
import sys, csv
from pathlib import Path
from collections import Counter
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from tape.io import iter_reads, load_whitelist
from tape.barcode import extract_bc, correct

TAGS = sys.argv[1:] or ["DTTz_1_S1", "DTTz_2_S2", "DTTz_3_S3", "DTTz_4_S4"]
TOP_N = 200              # report exactly this many post-correction sequences per embryo
K_FOLD = 5000            # fold the top-K raw barcodes (enough to capture all shadows + top N)
OUT = config.TABLES_DIR
OUT.mkdir(parents=True, exist_ok=True)

rows, cutoffs = [], []
for tag in TAGS:
    emb = config.get_embryo(tag)
    raw = Counter()
    for seq in iter_reads(emb.fastq_path):
        bc, _ = extract_bc(seq)
        if bc is not None and len(bc) == config.BC_LEN:
            raw[bc.decode()] += 1
    total = sum(raw.values())
    wl_bytes = load_whitelist(emb.whitelist_tsv)
    wl = set(w.decode() for w in wl_bytes)
    # Post-Hamming aggregation: fold each barcode onto its nearest whitelist member
    # within BC_CORRECT_MAX mismatches (the same correction the pipeline applies to
    # every read), so error shadows collapse into their real parent. Only the top
    # K_FOLD raw barcodes are folded (they contain every shadow of a real barcode
    # plus the whole top-N tail); barcodes not correcting to a whitelist member stay
    # as background. Reporting the top N regardless of count gives even facets.
    counts = Counter()
    for bc, n in raw.most_common(K_FOLD):
        corr = correct(bc.encode(), wl_bytes)
        counts[corr.decode() if corr is not None else bc] += n
    copies = {}
    dm = emb.out("depth_model.tsv")
    if dm.exists():
        for l in open(dm).read().splitlines()[1:]:
            c = l.split("\t"); copies[c[0]] = int(c[6])
    kept = counts.most_common(TOP_N)          # top N post-correction sequences (even facets)
    for rank, (bc, n) in enumerate(kept, 1):
        rows.append({"embryo": emb.label, "rank": rank, "tapebc": bc, "reads": n,
                     "freq": f"{n/total:.6f}", "whitelisted": int(bc in wl),
                     "copies": copies.get(bc, "")})
    wl_reads = [counts[b] for b in wl if b in counts]
    cutoffs.append({"embryo": emb.label, "tag": tag, "n_whitelisted": len(wl),
                    "total_bc_reads": total, "freq_cutoff": config.BC_MIN_FREQ,
                    "reads_at_freq_cutoff": int(round(config.BC_MIN_FREQ * total)),
                    "min_whitelisted_reads": min(wl_reads) if wl_reads else 0})
    print(f"[{tag}] {len(wl)} whitelisted / {len(counts):,} distinct barcodes; "
          f"cutoff ~{int(round(config.BC_MIN_FREQ*total)):,} reads (0.2%)")

with open(OUT / "suppfig1a_barcode_rank_abundance.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["embryo", "rank", "tapebc", "reads", "freq",
                                      "whitelisted", "copies"]); w.writeheader(); w.writerows(rows)
with open(OUT / "suppfig1a_cutoffs.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["embryo", "tag", "n_whitelisted", "total_bc_reads",
                                      "freq_cutoff", "reads_at_freq_cutoff",
                                      "min_whitelisted_reads"]); w.writeheader(); w.writerows(cutoffs)
print(f"wrote suppfig1a_barcode_rank_abundance.csv ({len(rows)} rows) + suppfig1a_cutoffs.csv")
