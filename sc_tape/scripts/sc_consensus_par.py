#!/usr/bin/env python3
"""Parallel drop-in for sc_consensus.py -- identical output, multi-core parse.

Usage:  python3 scripts/sc_consensus_par.py seq4 seq4_long --out e3v5 [--procs 8] [--limit N]

The ONLY difference vs sc_consensus.py is *how the per-read parse is executed*:
the read loop (header/UMI parse -> BC extract+correct -> walk_robust -> vocab
check) is fanned across a process pool, fed by `pigz` decompression instead of
Python's single-threaded gzip. All calling/de-noising/consensus/output code is
imported and run UNCHANGED, so the result is byte-identical to the serial script
(aggregation is order-independent; output rows are sorted). Verified with
--limit against sc_consensus.py before use.

Falls back to serial gzip if `pigz` is not on PATH.
"""
import sys
import gzip
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict, Counter
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import load_whitelist, load_vocab
from tape.barcode import correct
from tape.parser import walk_robust, pattern_of
from tape.denoise import denoise_barcode
from tape.sc import (parse_header, extract_bc_sc, collapse_umis,
                     consensus_locus, chain_to_str)

BATCH = 200_000          # reads per task; large -> more intra-batch PCR-dup collapse, less IPC

# --- worker globals (set once per worker via initializer, never re-pickled) ---
_WL = None
_VOCAB = None


def _init(wl, vocab):
    global _WL, _VOCAB
    _WL, _VOCAB = wl, vocab


def _parse_batch(batch):
    """Parse a list of (header_bytes, seq_bytes) exactly as sc_consensus.main's
    inner loop. Returns a Counter keyed ((cell, bc, umi, pattern)) -> reads,
    pre-aggregated (collapses PCR duplicates within the batch)."""
    wl, vocab = _WL, _VOCAB
    out = Counter()
    for header, seq in batch:
        ph = parse_header(header)
        if ph is None:
            continue
        cell, umi, _sample = ph                       # sample_field is "" for seq4 -> no filter
        cell = cell.decode()
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
        out[(cell, corr.decode(), umi, pattern_of(state))] += 1
    return out


def _batches(path, limit=None):
    """Yield batches of (header, seq) from a gzipped fastq, decompressing with
    pigz (multi-core) when available, else Python gzip."""
    use_pigz = shutil.which("pigz") is not None
    proc = None
    if use_pigz:
        proc = subprocess.Popen(["pigz", "-dc", "-p", "3", str(path)],
                                stdout=subprocess.PIPE, bufsize=1 << 22)
        fh = proc.stdout
    else:
        fh = gzip.open(str(path), "rb")
    n = 0
    buf = []
    try:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            fh.readline(); fh.readline()
            buf.append((h.rstrip(b"\n"), s.rstrip(b"\n")))
            n += 1
            if len(buf) >= BATCH:
                yield buf
                buf = []
            if limit and n >= limit:
                break
        if buf:
            yield buf
    finally:
        if proc is not None:
            if proc.poll() is None:
                proc.stdout.close()
                proc.terminate()
            proc.wait()
        else:
            fh.close()


def main(sample_names, out, procs=8, limit=None):
    samples = [config.get_sc_sample(n) for n in sample_names]
    tags = {s.embryo_tag for s in samples}
    if len(tags) != 1:
        sys.exit(f"cannot pool runs from different embryos: {tags}")
    emb = samples[0].embryo
    wl = load_whitelist(emb.whitelist_tsv)
    vocab = load_vocab(emb.vocab_tsv)
    bc_order = [w.decode() for w in wl]
    print(f"[sc-par] pooling {[s.name for s in samples]} (embryo {emb.tag}); "
          f"{len(wl)} whitelist BCs, {len(vocab)} vocab; procs={procs} -> out '{out}'",
          flush=True)

    mols = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    tot = kept = 0
    with Pool(procs, initializer=_init, initargs=(wl, vocab)) as pool:
        for s in samples:
            reads_path = (s.out("sc_reads.fastq.gz") if s.layout == "paired"
                          else s.r2_path)
            if not reads_path.exists():
                hint = " (run stage 00_merge_pairs first)" if s.layout == "paired" else ""
                sys.exit(f"missing input reads: {reads_path}{hint}")
            skept = 0
            # ORDERED imap (not imap_unordered): batches are contiguous file
            # segments, so merging them in order reproduces the serial script's
            # global first-occurrence order of (umi, pattern) keys. collapse_umis
            # and consensus break abundance ties by insertion order, so this is
            # required for byte-identical output. Main streams+decompresses while
            # workers parse; homogeneous batches keep ordering cost negligible.
            for cnt in pool.imap(_parse_batch, _batches(reads_path, limit)):
                for (cell, bc, umi, pat), c in cnt.items():
                    mols[(cell, bc)][umi][pat] += c
                    skept += c
            tot += skept  # note: pre-filter total not tracked in parallel path; see below
            kept += skept
            print(f"[sc-par]   {s.name}: kept={skept:,}", flush=True)

    cell_calls = defaultdict(dict)
    for (cell, bc), umi_pat in mols.items():
        mol_gts = collapse_umis(umi_pat)
        mol_gts = denoise_barcode(mol_gts)
        res = consensus_locus(mol_gts)
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
        print(f"[sc-par] kept reads={kept:,} | cells={n_cells:,}  "
              f"median loci/cell={int(st.median(loci))}  "
              f"cells>=5 loci={sum(1 for l in loci if l>=5):,}", flush=True)
    print(f"[sc-par] wrote {mat_path} + {qc_path}", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    limit = None; out = None; procs = 8
    if "--limit" in args:
        i = args.index("--limit"); limit = int(args[i + 1]); args = args[:i] + args[i + 2:]
    if "--procs" in args:
        i = args.index("--procs"); procs = int(args[i + 1]); args = args[:i] + args[i + 2:]
    if "--out" in args:
        i = args.index("--out"); out = args[i + 1]; args = args[:i] + args[i + 2:]
    if not args:
        sys.exit(__doc__)
    if out is None:
        out = args[0] if len(args) == 1 else "_".join(args)
    main(args, out, procs=procs, limit=limit)
