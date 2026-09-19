#!/usr/bin/env python3
"""Clean array-length measurement: validate EVERY anchor gap from the array start
through the terminal, independent of where the parse happened to stop."""
import sys, os
from collections import Counter
sys.path.insert(0, "/Users/jay.shendure/Dropbox/claude/mouse_sprint/tape_pipeline")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("TAPE_RAW_DIR", "/Users/shendure/tape_raw")
import config
from config import ANCHOR, N_SITES
from tape.io import iter_reads, load_whitelist
from tape.barcode import extract_bc, correct
from tape.parser import find_terminal, _valid_gap
from c8r_reparse import walk_instrumented

def main(tag):
    emb = config.get_embryo(tag); wl = load_whitelist(emb.whitelist_tsv)
    raw = Counter(); clean = Counter(); clean_bc = Counter()
    stopmode = Counter(); n = 0
    for seq in iter_reads(emb.fastq_path):
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN: continue
        n += 1
        anchors = []
        i = seq.find(ANCHOR, cs)
        while i >= 0:
            anchors.append(i); i = seq.find(ANCHOR, i + 6)
        if len(anchors) < 2: continue
        tpos = find_terminal(seq, cs)
        if tpos < 0: continue
        ti = None
        for idx, a in enumerate(anchors):
            if abs(a - tpos) <= 1: ti = idx; break
        if ti is None or ti < 1: continue
        raw[ti] += 1
        # validate EVERY gap from anchors[0] to the terminal, and the 5' offset
        ok = (anchors[0] - cs == 5) and all(_valid_gap(anchors[k+1]-anchors[k]) for k in range(ti))
        if ok:
            clean[ti] += 1
            corr = correct(bc, wl)
            if corr is not None: clean_bc[(corr.decode(), ti)] += 1
            state, e, ra, d = walk_instrumented(seq, cs)
            stopmode[(ti, d["fwd_stop"])] += 1
    tr, tc = sum(raw.values()), sum(clean.values())
    print(f"BC-extractable {n:,};  terminal located {tr:,};  fully-validated spacing {tc:,} ({100.0*tc/tr:.1f}%)\n")
    print("  monomers   raw (term_idx)        fully-validated spacing")
    for L in sorted(set(raw) | set(clean)):
        if L > 12: continue
        print(f"   {L:>3}   {raw.get(L,0):>10,} {100.0*raw.get(L,0)/tr:6.2f}%   {clean.get(L,0):>10,} {100.0*clean.get(L,0)/tc:6.2f}%")
    print("\n  per TAPE-BC, fully-validated spacing (% of that BC's validated reads):")
    bcs = sorted({b for b, _ in clean_bc})
    Ls = [L for L in sorted({L for _, L in clean_bc}) if 3 <= L <= 8]
    print("   %-14s %10s  %s" % ("tape-BC", "n", "  ".join("L%d" % L for L in Ls)))
    for b in bcs:
        t = sum(v for (bb, _), v in clean_bc.items() if bb == b)
        print(f"   {b:<14} {t:>10,}  " + "  ".join("%5.1f%%" % (100.0*clean_bc.get((b,L),0)/t) for L in Ls))
    print("\n  where the forward walk stopped, for validated reads of each length:")
    for L in sorted({l for l, _ in stopmode}):
        if not (4 <= L <= 7): continue
        tt = sum(v for (l, _), v in stopmode.items() if l == L)
        top = sorted(((s, v) for (l, s), v in stopmode.items() if l == L), key=lambda x: -x[1])[:4]
        print(f"   L{L} n={tt:>9,}: " + "  ".join(f"{s}={100.0*v/tt:.0f}%" for s, v in top))

main(sys.argv[1] if len(sys.argv) > 1 else "DTTz_3_S3")
