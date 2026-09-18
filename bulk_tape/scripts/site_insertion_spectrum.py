#!/usr/bin/env python3
"""Per-site insertion spectrum for a bulk embryo: how many distinct insertion
symbols exceed a frequency threshold at each site (and pooled over sites 1-2 vs
3-6). Answers 'N1 symbols >0.5% early, N2 late'.

Usage:  python3 scripts/site_insertion_spectrum.py DTTz_3_S3 [--min-freq 0.005]

Writes data/tables/<tag>.site_insertion_freq.csv (site, insertion, count, freq)
and prints the per-site-group counts above the threshold.
"""
import sys
from pathlib import Path
from collections import Counter, defaultdict
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from tape.io import iter_reads, load_whitelist
from tape.barcode import extract_bc, correct
from tape.parser import walk_robust

N = config.N_SITES


def main(tag, min_freq=0.005):
    emb = config.get_embryo(tag)
    wl = load_whitelist(emb.whitelist_tsv)
    site_ins = [Counter() for _ in range(N)]     # per site: insertion -> count
    for seq in iter_reads(emb.fastq_path):
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN:
            continue
        if correct(bc, wl) is None:
            continue
        state, err, _ = walk_robust(seq, cs)
        if err:
            continue
        for pos, s in enumerate(state):
            if s is not None and s != "U":
                site_ins[pos][s[1].decode()] += 1

    out = config.TABLES_DIR / f"{tag}.site_insertion_freq.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("site,insertion,count,freq_within_site\n")
        for pos in range(N):
            tot = sum(site_ins[pos].values()) or 1
            for ins, c in site_ins[pos].most_common():
                f.write(f"{pos+1},{ins},{c},{c/tot:.5f}\n")

    def group_count(sites):
        pooled = Counter()
        for pos in sites:
            pooled.update(site_ins[pos])
        tot = sum(pooled.values()) or 1
        return sum(1 for c in pooled.values() if c / tot >= min_freq), pooled, tot

    n12, _, t12 = group_count([0, 1])          # sites 1-2
    n36, _, t36 = group_count([2, 3, 4, 5])    # sites 3-6
    persite = [sum(1 for c in site_ins[p].values() if c / (sum(site_ins[p].values()) or 1) >= min_freq)
               for p in range(N)]
    print(f"[{tag}] insertions >= {min_freq*100:.1f}% frequency:")
    print(f"  per site (1..6): {persite}")
    print(f"  pooled sites 1-2: {n12}   pooled sites 3-6: {n36}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    args = sys.argv[1:]
    mf = 0.005
    if "--min-freq" in args:
        i = args.index("--min-freq"); mf = float(args[i + 1]); args = args[:i] + args[i + 2:]
    if len(args) != 1:
        sys.exit(__doc__)
    main(args[0], mf)
