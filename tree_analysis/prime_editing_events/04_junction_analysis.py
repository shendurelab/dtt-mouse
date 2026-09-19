#!/usr/bin/env python3
"""Characterise the junction-failure classes, GAP_SHORT, ordering violations,
and settle the cause of the site-5/6 duplication. Reviewer comment (8)."""
import sys, os
from collections import Counter
sys.path.insert(0, "/Users/jay.shendure/Dropbox/claude/mouse_sprint/tape_pipeline")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("TAPE_RAW_DIR", "/Users/shendure/tape_raw")
import config
from config import ANCHOR, LEADIN, N_SITES, GAP_MIN, GAP_MAX
from tape.io import iter_reads, load_whitelist, load_vocab
from tape.barcode import extract_bc, correct, hamming_le
from tape.parser import find_terminal, _valid_gap
from c8r_reparse import walk_instrumented

# canonical junction segment lengths: 4 = unedited (TGAT); 3bp insertion -> NNN+GGA+TGAT = 10;
# 4bp insertion (GATG) -> 11. Anything else implies an indel, not a substitution.
CANON = {4: "U", 10: "E3", 11: "E4"}


def main(tag):
    emb = config.get_embryo(tag)
    wl = load_whitelist(emb.whitelist_tsv); vocab = load_vocab(emb.vocab_tsv)
    notgat_len = Counter(); notgat_h1 = Counter()
    nogga_len = Counter(); nogga_h1 = Counter()
    gapshort = Counter()
    order_site = Counter(); order_len = Counter()
    off0 = Counter(); off0_res = Counter()
    termexact = Counter()

    for seq in iter_reads(emb.fastq_path):
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN: continue
        state, e, ra, d = walk_instrumented(seq, cs)

        # --- offset of the first detected anchor from the array start.
        # A correct parse puts monomer 1 at cstart + len(LEADIN) = +5.
        # A larger offset means the FIRST anchor was missed and the whole array is frame-shifted.
        a0 = seq.find(ANCHOR, cs)
        if a0 >= 0 and d["term_idx"] is not None:
            o = a0 - cs
            off0[(d["term_idx"], o if o <= 40 else 41)] += 1
            if ra:
                off0_res[(d["term_idx"], o if o <= 40 else 41)] += 1
                termexact[(d["term_idx"], seq.find(config.TERM, cs) >= 0)] += 1

        stop = d["fwd_stop"]
        if stop in ("JUNC_NO_TGAT", "JUNC_NO_GGA", "GAP_SHORT") and d.get("fwd_a0") is not None:
            pass
        # recompute the offending segment (walk_instrumented records site + detail only)
        if stop in ("JUNC_NO_TGAT", "JUNC_NO_GGA"):
            anchors = []
            i = seq.find(ANCHOR, cs)
            while i >= 0:
                anchors.append(i); i = seq.find(ANCHOR, i + 6)
            k = d["fwd_site"] - 1
            if k + 1 < len(anchors):
                seg = seq[anchors[k] + 10:anchors[k + 1]]
                L = len(seg)
                if stop == "JUNC_NO_TGAT":
                    notgat_len[(L, CANON.get(L, "odd"))] += 1
                    notgat_h1[(CANON.get(L, "odd"),
                               len(seg) >= 4 and hamming_le(seg[-4:], b"TGAT", 1) <= 1)] += 1
                else:
                    blk = seg[:-4]
                    nogga_len[(L, CANON.get(L, "odd"))] += 1
                    nogga_h1[(CANON.get(L, "odd"),
                              len(blk) >= 3 and hamming_le(blk[-3:], b"GGA", 1) <= 1)] += 1
        if stop == "GAP_SHORT":
            gapshort[d["fwd_detail"]] += 1
        if e:
            pat = tuple("U" if s == "U" else ("E" if s is not None else ".") for s in state)
            order_site[pat] += 1
            order_len[d["fwd_site"]] += 1

    def show(name, lenc, h1c, target):
        tt = sum(lenc.values())
        print(f"\n=== {name}: junction segment length (n = {tt:,}) ===")
        for (L, kind), n in sorted(lenc.items(), key=lambda x: -x[1])[:12]:
            print(f"   len {L:>3} ({kind:>3})  {n:>10,}  {100.0*n/tt:5.1f}%")
        print(f"   -- canonical length (substitution) vs odd length (indel):")
        can = sum(n for (L, k), n in lenc.items() if k != "odd")
        print(f"      canonical {can:,} ({100.0*can/tt:.1f}%)   odd {tt-can:,} ({100.0*(tt-can)/tt:.1f}%)")
        print(f"   -- among canonical-length, is {target.decode()} present at Hamming<=1?")
        for kind in ("U", "E3", "E4", "odd"):
            y, n = h1c.get((kind, True), 0), h1c.get((kind, False), 0)
            if y + n: print(f"      {kind:>3}: {y:>9,} of {y+n:>9,}  ({100.0*y/(y+n):5.1f}%)")

    show("JUNC_NO_TGAT", notgat_len, notgat_h1, b"TGAT")
    show("JUNC_NO_GGA", nogga_len, nogga_h1, b"GGA")

    print(f"\n=== GAP_SHORT gaps (n = {sum(gapshort.values()):,}) ===")
    for g, n in sorted(gapshort.items()): print(f"   gap {g:>3}  {n:>8,}")

    print(f"\n=== ORDERING VIOLATIONS (n = {sum(order_len.values()):,}) ===")
    print("   resolved state pattern (E=edit, U=unedited, .=unresolved):")
    for pat, n in order_site.most_common(10): print(f"     {''.join(pat)}  {n:>8,}")

    print("\n=== FIRST-ANCHOR OFFSET vs IMPLIED ARRAY LENGTH (fully-resolved reads) ===")
    print("   a correct parse has monomer 1 at cstart+5; +19/+25 means anchor 1 was missed")
    Ls = sorted({L for L, _ in off0_res if L is not None and 3 <= L <= 7})
    offs = sorted({o for _, o in off0_res})
    print("   %-14s %s" % ("array len", "  ".join(f"+{o}" for o in offs if o <= 30)))
    for L in Ls:
        tt = sum(n for (LL, _), n in off0_res.items() if LL == L)
        row = "  ".join("%5.1f%%" % (100.0*off0_res.get((L, o), 0)/tt) for o in offs if o <= 30)
        ex = termexact.get((L, True), 0); fz = termexact.get((L, False), 0)
        print(f"   L{L} n={tt:>9,}  {row}   | terminal exact-match: {100.0*ex/(ex+fz) if ex+fz else 0:.1f}%")


main(sys.argv[1] if len(sys.argv) > 1 else "DTTz_3_S3")
