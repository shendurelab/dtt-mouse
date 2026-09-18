#!/usr/bin/env python3
"""Single-cell TAPE: parse reads -> per-(cell, integration) consensus genotype.

Usage:  python3 scripts/sc_consensus.py seq1 seq3 --out e3 [--limit N]
        python3 scripts/sc_consensus.py seq1            [--limit N]   # single run

Multiple runs of the SAME embryo are POOLED (by cell, UMI, integration) into one
character matrix -- e.g. seq1 (single-end) + seq3 (paired, pre-merged) for embryo 3.

Input reads per run:
  paired run -> data/<run>.sc_reads.fastq.gz (from stage 00_merge_pairs)
  single run -> the run's R2 amplicon fastq directly

For each read: parse the cell id + UMI from the header, extract + whitelist-correct
the integration barcode (using the bulk-derived references), resolve the 6-site
array (shared walk_robust), and keep valid reads. Then, per (cell, integration):
UMIs are collapsed to molecules; the molecule genotypes are DE-NOISED exactly as
in bulk (error-collapse -> chimera removal -> >=2 support); and a consensus
genotype is called via the edit-chain model. Writes a cell x integration matrix + QC.

Outputs (in DATA_DIR):
  <out>.cell_consensus.tsv   cell x integration consensus genotypes (+ QC cols)
  <out>.sc_qc.tsv            per-cell QC summary
"""
import sys
import gzip
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import load_whitelist, load_vocab
from tape.barcode import correct
from tape.parser import walk_robust, pattern_of
from tape.denoise import denoise_barcode
from tape.sc import parse_header, extract_bc_sc, collapse_umis, consensus_locus, chain_to_str


def _iter_with_header(path, limit):
    """Yield (header_bytes, seq_bytes) from a (gzipped) fastq."""
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
    tags = {s.embryo_tag for s in samples}
    if len(tags) != 1:
        sys.exit(f"cannot pool runs from different embryos: {tags}")
    emb = samples[0].embryo
    wl = load_whitelist(emb.whitelist_tsv)          # references DERIVED FROM BULK
    vocab = load_vocab(emb.vocab_tsv)
    bc_order = [w.decode() for w in wl]
    print(f"[sc] pooling {[s.name for s in samples]} (embryo tag {emb.tag}); "
          f"{len(wl)} whitelist BCs, {len(vocab)} vocab -> out '{out}'")

    # mols[(cell, bc)][umi][pattern] = reads  (pooled across runs)
    mols = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    tot = kept = 0
    for s in samples:
        reads_path = (s.out("sc_reads.fastq.gz") if s.layout == "paired" else s.r2_path)
        if not reads_path.exists():
            hint = " (run stage 00_merge_pairs first)" if s.layout == "paired" else ""
            sys.exit(f"missing input reads: {reads_path}{hint}")
        stot = skept = 0
        for header, seq in _iter_with_header(reads_path, limit):
            stot += 1
            ph = parse_header(header)
            if ph is None:
                continue
            cell, umi, sample = ph
            if s.sample_field and s.sample_field.encode() not in sample:
                continue
            cell = cell.decode()                     # str key/output; umi stays bytes (Hamming)
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
            mols[(cell, corr.decode())][umi][pattern_of(state)] += 1
            skept += 1
        tot += stot; kept += skept
        print(f"[sc]   {s.name}: reads={stot:,} kept={skept:,} ({skept/max(1,stot)*100:.1f}%)")

    # per (cell, integration): UMI-collapse (>= MIN_READS_PER_UMI reads/UMI) -> bulk
    # de-noising -> edit-chain call, kept only if the dominant lineage clears
    # SC_MIN_MOLECULES and SC_DOMINANCE_MIN (else left missing). See tape/sc.py.
    cell_calls = defaultdict(dict)     # cell -> {bc: (genotype_str, branched, dominance, n_mol)}
    for (cell, bc), umi_pat in mols.items():
        mol_gts = collapse_umis(umi_pat)                    # PCR dups -> molecules (reads/UMI floor)
        mol_gts = denoise_barcode(mol_gts)                  # error-collapse + chimera + >=2 (as bulk)
        res = consensus_locus(mol_gts)                      # applies SC_MIN_MOLECULES + SC_DOMINANCE_MIN
        if res is None:
            continue
        chain, branched, dominance, n = res
        cell_calls[cell][bc] = (chain_to_str(chain), branched, dominance, n)

    mat_path = config.DATA_DIR / f"{out}.cell_consensus.tsv"
    qc_path = config.DATA_DIR / f"{out}.sc_qc.tsv"
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(mat_path, "w") as f:
        f.write("cell\t" + "\t".join(bc_order) + "\tn_loci\tn_doublet_loci\tmean_dominance\n")
        for cell in sorted(cell_calls):
            calls = cell_calls[cell]
            row = [calls[b][0] if b in calls else "" for b in bc_order]
            n_loci = len(calls); n_db = sum(1 for v in calls.values() if v[1])
            mdom = sum(v[2] for v in calls.values()) / max(1, n_loci)
            f.write(f"{cell}\t" + "\t".join(row) + f"\t{n_loci}\t{n_db}\t{mdom:.3f}\n")
    with open(qc_path, "w") as f:
        f.write("cell\tn_loci\tn_doublet_loci\tmean_dominance\ttotal_molecules\n")
        for cell in sorted(cell_calls):
            calls = cell_calls[cell]
            n_loci = len(calls); n_db = sum(1 for v in calls.values() if v[1])
            mdom = sum(v[2] for v in calls.values()) / max(1, n_loci)
            tmol = sum(v[3] for v in calls.values())
            f.write(f"{cell}\t{n_loci}\t{n_db}\t{mdom:.3f}\t{tmol}\n")

    n_cells = len(cell_calls)
    if n_cells:
        import statistics as st
        loci = [len(c) for c in cell_calls.values()]
        print(f"[sc] pooled reads={tot:,} kept={kept:,} | cells={n_cells:,}  "
              f"median loci/cell={int(st.median(loci))}  cells>=5 loci={sum(1 for l in loci if l>=5):,}")
    print(f"[sc] wrote {mat_path} + {qc_path}")


if __name__ == "__main__":
    args = sys.argv[1:]
    limit = None
    out = None
    if "--limit" in args:
        i = args.index("--limit"); limit = int(args[i + 1]); args = args[:i] + args[i + 2:]
    if "--out" in args:
        i = args.index("--out"); out = args[i + 1]; args = args[:i] + args[i + 2:]
    if not args:
        sys.exit(__doc__)
    if out is None:
        out = args[0] if len(args) == 1 else "_".join(args)
    main(args, out, limit)
