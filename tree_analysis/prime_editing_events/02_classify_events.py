#!/usr/bin/env python3
"""Cross-tab forward-walk outcome against resolved_all, and classify every
insertion call in fully-resolved reads. Reviewer comment (8), bulk embryo #3."""
import sys, os
from collections import Counter
sys.path.insert(0, "/Users/jay.shendure/Dropbox/claude/mouse_sprint/tape_pipeline")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("TAPE_RAW_DIR", "/Users/shendure/tape_raw")
import config
from tape.io import iter_reads, load_whitelist, load_vocab
from tape.barcode import extract_bc, correct
from c8r_reparse import walk_instrumented

# 13-symbol vocabulary splits into the 8 written throughout the array and the 5
# confined to the 5'-most monomers (episomal epegRNA products; response to (6)).
EPISOMAL = {b"ACC", b"GGC", b"ACG", b"GAA", b"CCG"}

def main(tag):
    emb = config.get_embryo(tag)
    wl = load_whitelist(emb.whitelist_tsv); vocab = load_vocab(emb.vocab_tsv)
    persistent = {v for v in vocab} - EPISOMAL
    print(f"vocab {len(vocab)}: persistent {len(persistent)} {sorted(x.decode() for x in persistent)}")
    print(f"          episomal   {len(EPISOMAL & set(vocab))} {sorted(x.decode() for x in EPISOMAL)}")

    xt = Counter(); tot = bcok = 0
    calls = Counter(); lens = Counter(); nv = Counter(); readclass = Counter()
    for seq in iter_reads(emb.fastq_path):
        tot += 1
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN: continue
        bcok += 1
        corr = correct(bc, wl)
        state, e, ra, d = walk_instrumented(seq, cs)
        stop = "ORDER_VIOLATION" if e else (d["fwd_stop"] or "NONE")
        xt[(stop, ra)] += 1
        if not ra or corr is None: continue
        kinds = set()
        for s in state:
            if s == "U": continue
            nnn = s[1]
            k = "persistent8" if nnn in persistent else ("episomal5" if nnn in EPISOMAL else "nonvocab")
            calls[k] += 1; lens[(k, len(nnn))] += 1; kinds.add(k)
            if k == "nonvocab": nv[nnn.decode()] += 1
        readclass["nonvocab_any" if "nonvocab" in kinds else "all_vocab"] += 1

    print(f"\ntotal reads {tot:,}; BC-extractable {bcok:,}")
    print("\n=== FORWARD-WALK OUTCOME x ALL-SITES-RESOLVED (of BC-extractable) ===")
    print("  %-20s %12s %12s %12s  %s" % ("fwd_stop", "resolved", "NOT resolved", "total", "%resolved"))
    stops = sorted({s for s, _ in xt}, key=lambda s: -(xt.get((s,True),0)+xt.get((s,False),0)))
    for s in stops:
        r, n = xt.get((s,True),0), xt.get((s,False),0)
        print(f"  {s:<20} {r:>12,} {n:>12,} {r+n:>12,}  {(100.0*r/(r+n) if r+n else 0):6.1f}%")
    R = sum(v for (s,ok),v in xt.items() if ok)
    print(f"  {'TOTAL resolved':<20} {R:>12,}  ({100.0*R/bcok:.2f}% of BC-extractable)")

    tc = sum(calls.values())
    print(f"\n=== INSERTION CALLS in fully-resolved, whitelisted reads (n = {tc:,} calls) ===")
    for k in ("persistent8", "episomal5", "nonvocab"):
        print(f"  {k:<12} {calls[k]:>12,}  {100.0*calls[k]/tc:6.2f}%")
    print("\n  by insertion length:")
    for k in ("persistent8", "episomal5", "nonvocab"):
        row = {L: n for (kk, L), n in lens.items() if kk == k}
        print("   %-12s %s" % (k, "  ".join(f"{L}bp:{n:,}" for L, n in sorted(row.items()))))
    tr = sum(readclass.values())
    print(f"\n  reads with >=1 non-vocab call: {readclass['nonvocab_any']:,} of {tr:,} ({100.0*readclass['nonvocab_any']/tr:.2f}%)")
    print("\n  top non-vocabulary symbols:")
    for s, n in nv.most_common(15):
        print(f"   {s:<6} len {len(s)} {n:>10,}  {100.0*n/tc:5.2f}% of all calls")

main(sys.argv[1] if len(sys.argv) > 1 else "DTTz_3_S3")
