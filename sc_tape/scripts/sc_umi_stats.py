#!/usr/bin/env python3
"""Single-cell circTAPE PCR-duplication + UMI depth (SC results paragraph).

Read-level statistics that are NOT stored in the consensus outputs, so this
re-parses the circTAPE reads with the SAME front-end as sc_consensus.py
(header -> cell+UMI, whitelist-corrected barcode, walk_robust, vocab filter),
pools runs of one embryo, collapses UMIs to molecules exactly as the pipeline
does (tape.sc.collapse_umis), and reports:

  - reads kept, molecules (post-UMI-collapse), reads/UMI, PCR-duplication rate
  - molecules (UMIs) per cell: mean / median

Duplication rate = 1 - molecules/reads;  reads-per-UMI = reads/molecules.

Usage:  python3 scripts/sc_umi_stats.py seq1 seq3 --out e3 [--limit N]

Writes to tables/:
  <out>.sc_umi_stats.csv   metric, value  (one row per statistic)
"""
import sys, csv, gzip, statistics as st
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import load_whitelist, load_vocab
from tape.barcode import correct
from tape.parser import walk_robust, pattern_of
from tape.sc import parse_header, extract_bc_sc, collapse_umis


def _iter(path, limit):
    n = 0
    with gzip.open(str(path), "rb") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip(b"\n")
            fh.readline(); fh.readline()
            yield h.rstrip(b"\n"), s
            n += 1
            if limit and n >= limit:
                break


def main(sample_names, out, limit=None):
    samples = [config.get_sc_sample(n) for n in sample_names]
    emb = samples[0].embryo
    wl = load_whitelist(emb.whitelist_tsv)
    vocab = load_vocab(emb.vocab_tsv)

    # mols[(cell, bc)][umi][pattern] = reads   (pooled, identical to sc_consensus)
    mols = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    tot = kept = 0
    for s in samples:
        reads_path = (s.out("sc_reads.fastq.gz") if s.layout == "paired" else s.r2_path)
        if not reads_path.exists():
            sys.exit(f"missing input reads: {reads_path}")
        for header, seq in _iter(reads_path, limit):
            tot += 1
            ph = parse_header(header)
            if ph is None:
                continue
            cell, umi, sample = ph
            if s.sample_field and s.sample_field.encode() not in sample:
                continue
            bc, cs = extract_bc_sc(seq)
            if bc is None:
                continue
            corr = correct(bc, wl)
            if corr is None:
                continue
            state, e, ra = walk_robust(seq, cs)
            if e or not ra:
                continue
            if not all((x == "U") or (x[0] == "E" and x[1] in vocab) for x in state):
                continue
            mols[(cell.decode(), corr.decode())][umi][pattern_of(state)] += 1
            kept += 1
        print(f"[umi] {s.name}: read {tot:,} so far, kept {kept:,}")

    # collapse UMIs -> molecules; tally reads and molecules, and molecules per cell
    reads_kept = molecules = 0
    per_cell = defaultdict(int)
    for (cell, bc), umi_pat in mols.items():
        reads_kept += sum(sum(p.values()) for p in umi_pat.values())
        # raw UMI depth stats: keep EVERY UMI (min_reads_per_umi=1); the calling
        # pipeline applies the MIN_READS_PER_UMI floor, this report intentionally does not.
        m = int(sum(collapse_umis(umi_pat, min_reads_per_umi=1).values()))
        molecules += m
        per_cell[cell] += m
    pc = list(per_cell.values())

    reads_per_umi = reads_kept / molecules if molecules else 0.0
    dup = 1 - molecules / reads_kept if reads_kept else 0.0

    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    outp = config.TABLES_DIR / f"{out}.sc_umi_stats.csv"
    with open(outp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["reads_kept", reads_kept])
        w.writerow(["molecules_umis", molecules])
        w.writerow(["reads_per_umi", round(reads_per_umi, 2)])
        w.writerow(["pcr_duplication_rate", round(dup, 3)])
        w.writerow(["cells", len(pc)])
        w.writerow(["umis_per_cell_mean", round(st.mean(pc), 1)])
        w.writerow(["umis_per_cell_median", int(st.median(pc))])

    print(f"[umi] reads_kept={reads_kept:,} molecules={molecules:,} "
          f"reads/UMI={reads_per_umi:.2f} duplication={dup*100:.1f}%")
    print(f"[umi] UMIs/cell mean {st.mean(pc):.1f} median {int(st.median(pc))} "
          f"over {len(pc):,} cells")
    print(f"[umi] wrote {outp}")


if __name__ == "__main__":
    args = sys.argv[1:]
    limit = out = None
    if "--limit" in args:
        i = args.index("--limit"); limit = int(args[i + 1]); args = args[:i] + args[i + 2:]
    if "--out" in args:
        i = args.index("--out"); out = args[i + 1]; args = args[:i] + args[i + 2:]
    if not args:
        sys.exit(__doc__)
    main(args, out or "_".join(args), limit)
