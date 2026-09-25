#!/usr/bin/env python3
"""Instrumented re-parse of bulk TAPE reads -- reviewer comment (8).

The production parser (mouse_sprint/tape_pipeline/tape/parser.py) resolves the
6-site array but records only WHETHER a read failed, never WHY: classify_junction
returns None for every unrecognised junction and walk_robust simply `break`s.
A deletion, a non-canonical junction, an out-of-range insertion and a read that
merely ran out of length all collapse into one bucket.

This fork mirrors walk_robust EXACTLY -- same constants, same control flow, same
break conditions -- but returns a reason code for every stop, plus the array
length implied by the anchor count. Nothing in the main pipeline is written to.

Validation: the funnel counts printed here must reproduce
data/<tag>.parse_summary.tsv exactly. If they don't, the fork has drifted.

Usage:  c8r_reparse.py <tag> [--limit N]
"""
from __future__ import annotations
import sys, os
from collections import Counter, defaultdict

PIPE = _os.path.join(_REPO, "sc_tape")
sys.path.insert(0, PIPE)
os.environ.setdefault("TAPE_RAW_DIR", os.environ.get("TAPE_RAW_DIR", ""))  # set TAPE_RAW_DIR; raw reads are on GEO GSE341627

import config
from config import ANCHOR, TERM, TERM_LONG, N_SITES, GAP_MIN, GAP_MAX
from tape.io import iter_reads, load_whitelist, load_vocab
from tape.barcode import extract_bc, correct, hamming_le
from tape.parser import find_terminal, _valid_gap

OUT = os.path.dirname(os.path.abspath(__file__))


def why_junction(seg: bytes):
    """classify_junction, but reporting the reason it fails.

    Returns (state, code, detail) where state is ("U",) | ("E", nnn) | None.
    """
    if not seg.endswith(b"TGAT"):
        return None, "JUNC_NO_TGAT", len(seg)
    block = seg[:-4]
    if block == b"":
        return ("U",), "OK_U", 0
    if block.endswith(b"GGA"):
        nnn = block[:-3]
        if 1 <= len(nnn) <= 6:
            return ("E", nnn), "OK_E", len(nnn)
        return None, "INS_LEN_OOR", len(nnn)      # insertion outside the 1-6 bp window
    return None, "JUNC_NO_GGA", len(block)


def walk_instrumented(seq: bytes, cstart: int):
    """Mirror of walk_robust with diagnostics. Returns (state, err, resolved_all, d)."""
    d = {"n_anchors": 0, "term_found": False, "term_idx": None,
         "fwd_stop": None, "fwd_site": None, "fwd_detail": None,
         "bwd_stop": None, "bwd_site": None, "bwd_detail": None,
         "order_violation": False, "readlen": len(seq),
         "fwd_a0": None, "gaps_clean": True}
    state = [None] * N_SITES

    anchors = []
    i = seq.find(ANCHOR, cstart)
    while i >= 0:
        anchors.append(i)
        i = seq.find(ANCHOR, i + 6)
    d["n_anchors"] = len(anchors)
    if len(anchors) < 2:
        d["fwd_stop"] = "NO_ANCHOR"
        return state, False, False, d

    tpos = find_terminal(seq, cstart)
    term_idx = None
    if tpos >= 0:
        for idx, a in enumerate(anchors):
            if abs(a - tpos) <= 1:
                term_idx = idx
                break
        if term_idx is None:
            anchors = sorted([a for a in anchors if a < tpos - 1] + [tpos])
            term_idx = anchors.index(tpos)
    d["term_found"] = tpos >= 0
    d["term_idx"] = term_idx

    # ---- FORWARD
    saw_u = False
    nfwd = min(len(anchors) - 1, N_SITES)
    for k in range(nfwd):
        gap = anchors[k + 1] - anchors[k]
        if not _valid_gap(gap):
            d["fwd_stop"] = "GAP_SHORT" if gap < GAP_MIN else "GAP_LONG"
            d["fwd_site"], d["fwd_detail"] = k + 1, gap
            d["fwd_a0"] = anchors[k]
            d["gaps_clean"] = False
            break
        st, code, detail = why_junction(seq[anchors[k] + 10:anchors[k + 1]])
        if st is None:
            d["fwd_stop"], d["fwd_site"], d["fwd_detail"] = code, k + 1, detail
            break
        if st[0] == "U":
            for j in range(k, N_SITES):
                state[j] = "U"
            saw_u = True
            d["fwd_stop"], d["fwd_site"] = "ORDERED_U", k + 1
            break
        state[k] = st
    else:
        if nfwd < N_SITES:
            d["fwd_stop"], d["fwd_site"] = "RAN_OUT_ANCHORS", nfwd + 1
        else:
            d["fwd_stop"] = "COMPLETE"

    # ---- BACKWARD from the terminal
    if not saw_u and term_idx is not None and term_idx >= 1:
        for j in range(N_SITES):
            s = N_SITES - 1 - j
            if state[s] is not None:
                d["bwd_stop"] = d["bwd_stop"] or "MET_FORWARD"
                break
            ai = term_idx - j
            if ai - 1 < 0 or ai >= len(anchors):
                d["bwd_stop"], d["bwd_site"] = "BWD_RANGE", s + 1
                break
            gap = anchors[ai] - anchors[ai - 1]
            if not _valid_gap(gap):
                d["bwd_stop"] = "BWD_GAP_SHORT" if gap < GAP_MIN else "BWD_GAP_LONG"
                d["bwd_site"], d["bwd_detail"] = s + 1, gap
                break
            st, code, detail = why_junction(seq[anchors[ai - 1] + 10:anchors[ai]])
            if st is None:
                d["bwd_stop"], d["bwd_site"], d["bwd_detail"] = "BWD_" + code, s + 1, detail
                break
            state[s] = "U" if st[0] == "U" else st

    seen_u = False
    for s in state:
        if s is None:
            continue
        if s == "U":
            seen_u = True
        elif seen_u:
            d["order_violation"] = True
            return state, True, False, d

    return state, False, all(x is not None for x in state), d


def main(tag: str, limit=None):
    emb = config.get_embryo(tag)
    wl = load_whitelist(emb.whitelist_tsv)
    vocab = load_vocab(emb.vocab_tsv)
    print(f"[c8r] {emb.label} ({tag}); {len(wl)} whitelist BCs, {len(vocab)} vocab symbols")
    print(f"[c8r] fastq: {emb.fastq_path}")

    tot = bc_ok = bc_wl = err = res_all = valid = 0
    fwd = Counter(); bwd = Counter(); fwd_site = Counter()
    gaps = Counter(); inslen = Counter(); nonvocab = Counter()
    arraylen = Counter(); arraylen_res = Counter(); arraylen_clean = Counter()
    term_missing = 0
    anchor_h1 = Counter(); gapfam = Counter(); bc_len = Counter(); bc_s56 = Counter()
    s56 = Counter()                      # (n_monomers, same/diff) for reads edited at 5 and 6
    nonvocab_reads = 0
    readlen_fail = Counter(); readlen_all = Counter()

    for seq in iter_reads(emb.fastq_path, limit):
        tot += 1
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != config.BC_LEN:
            continue
        bc_ok += 1
        corr = correct(bc, wl)
        if corr is not None:
            bc_wl += 1
        state, e, ra, d = walk_instrumented(seq, cs)

        readlen_all[len(seq) // 25 * 25] += 1
        if not d["term_found"]:
            term_missing += 1
        if d["term_idx"] is not None:
            arraylen[d["term_idx"]] += 1
            if d["gaps_clean"]:
                arraylen_clean[d["term_idx"]] += 1

        if e:
            err += 1
            fwd["ORDER_VIOLATION"] += 1
            continue

        fwd[d["fwd_stop"] or "NONE"] += 1
        if d["fwd_stop"] not in ("COMPLETE", "ORDERED_U"):
            fwd_site[(d["fwd_stop"], d["fwd_site"])] += 1
            readlen_fail[len(seq) // 25 * 25] += 1
            if d["fwd_stop"] in ("GAP_SHORT", "GAP_LONG"):
                gaps[d["fwd_detail"]] += 1
                if d["fwd_stop"] == "GAP_LONG":
                    g = d["fwd_detail"]
                    fam = {28: "UU", 34: "UE", 40: "EE", 41: "EE+1"}.get(g, "other")
                    gapfam[fam] += 1
                    a0 = d["fwd_a0"]
                    seg = seq[a0 + 6:a0 + g]          # between the two flanking anchors
                    hit = any(hamming_le(seg[i:i + 6], ANCHOR, 1) <= 1
                              for i in range(max(0, len(seg) - 5)))
                    anchor_h1[(fam, hit)] += 1
            if d["fwd_stop"] == "INS_LEN_OOR":
                inslen[d["fwd_detail"]] += 1
        if d["bwd_stop"]:
            bwd[d["bwd_stop"]] += 1

        if ra:
            res_all += 1
            if corr is not None and d["term_idx"] is not None:
                bc_len[(corr.decode(), d["term_idx"])] += 1
                a5, b6 = state[4], state[5]
                if a5 != "U" and b6 != "U":
                    bc_s56[(corr.decode(), d["term_idx"], a5[1] == b6[1])] += 1
            if d["term_idx"] is not None:
                arraylen_res[d["term_idx"]] += 1
            ok_vocab = True
            for s in state:
                if s != "U":
                    if s[1] in vocab:
                        pass
                    else:
                        ok_vocab = False
                        nonvocab[s[1].decode()] += 1
                    inslen[len(s[1])] += 0      # keep key space consistent
            if not ok_vocab:
                nonvocab_reads += 1
            if corr is not None and ok_vocab:
                valid += 1
            # site 5 vs 6 identity, by implied array length
            a, b = state[4], state[5]
            if a != "U" and b != "U":
                s56[(d["term_idx"], a[1] == b[1])] += 1

    # ---------------- report
    def pct(n): return 100.0 * n / tot if tot else 0.0
    print("\n=== FUNNEL (must match parse_summary.tsv) ===")
    for lab, n in [("total_reads", tot), ("bc_extractable", bc_ok), ("bc_in_whitelist", bc_wl),
                   ("ordering_errors", err), ("all_sites_resolved", res_all), ("valid_reads", valid)]:
        print(f"{lab:22s} {n:>10,}  {pct(n):6.2f}%")

    print("\n=== FORWARD-WALK OUTCOME (of BC-extractable reads) ===")
    for k, n in fwd.most_common():
        print(f"{k:22s} {n:>10,}  {100.0*n/bc_ok:6.2f}%")

    print("\n=== FAILURE SITE (non-terminating stops) ===")
    for (k, s), n in sorted(fwd_site.items(), key=lambda x: -x[1])[:25]:
        print(f"{k:22s} site {s}  {n:>10,}")

    print("\n=== ANCHOR GAPS AT FAILURE (bp; valid window %d-%d) ===" % (GAP_MIN, GAP_MAX))
    for g, n in sorted(gaps.items())[:40]:
        print(f"  gap {g:>4}  {n:>10,}")

    print("\n=== IMPLIED ARRAY LENGTH (# body monomers before terminal) ===")
    tota = sum(arraylen.values())
    for L, n in sorted(arraylen.items()):
        print(f"  {L:>3} monomers  {n:>10,}  {100.0*n/tota:6.3f}%   (resolved-all: {arraylen_res.get(L,0):,})")
    print(f"  terminal not found in {term_missing:,} reads ({pct(term_missing):.2f}% of all)")

    print("\n=== GAP_LONG: IS THE 'MISSING' ANCHOR PRESENT WITH ONE MISMATCH? ===")
    for fam in ("UU", "UE", "EE", "EE+1", "other"):
        h = anchor_h1.get((fam, True), 0); m = anchor_h1.get((fam, False), 0)
        if h + m:
            print(f"  gap family {fam:<5} n={h+m:>9,}  anchor recoverable at Hamming<=1: {h:>9,} ({100.0*h/(h+m):5.1f}%)")

    print("\n=== ARRAY LENGTH, reads with fully canonical anchor spacing ===")
    tc = sum(arraylen_clean.values())
    for L, n in sorted(arraylen_clean.items()):
        print(f"  {L:>3} monomers  {n:>10,}  {100.0*n/tc:6.3f}%")

    print("\n=== SITE 5 vs SITE 6 SYMBOL IDENTITY, by array length ===")
    lens = sorted({L for L, _ in s56}, key=lambda x: (x is None, x))
    for L in lens:
        same = s56.get((L, True), 0); diff = s56.get((L, False), 0)
        if same + diff:
            print(f"  {str(L):>4} monomers: same={same:>9,} diff={diff:>9,}  P(same)={same/(same+diff):.4f}")

    print("\n=== ARRAY LENGTH PER TAPE-BC (fully-resolved reads) ===")
    bcs = sorted({b for b, _ in bc_len})
    Ls = [L for L in sorted({L for _, L in bc_len}) if 3 <= L <= 8]
    print("  %-14s %10s  %s" % ("tape-BC", "total", "  ".join("L%d" % L for L in Ls)))
    for b in bcs:
        tt = sum(n for (bb, _), n in bc_len.items() if bb == b)
        cells = "  ".join("%5.1f%%" % (100.0 * bc_len.get((b, L), 0) / tt) for L in Ls)
        s5 = sum(n for (bb, L, sm), n in bc_s56.items() if bb == b and L == 5 and sm)
        d5 = sum(n for (bb, L, sm), n in bc_s56.items() if bb == b and L == 5 and not sm)
        print(f"  {b:<14} {tt:>10,}  {cells}   | L5 s5==s6: {s5:,} vs {d5:,}")

    print(f"\n=== NON-VOCABULARY INSERTIONS (in otherwise fully-resolved reads) ===")
    print(f"reads with >=1 non-vocab insertion: {nonvocab_reads:,} ({pct(nonvocab_reads):.2f}% of total)")
    print(f"distinct non-vocab symbols: {len(nonvocab):,}")
    for s, n in nonvocab.most_common(30):
        print(f"  {s:<8} len {len(s)}  {n:>10,}")

    print("\n=== OUT-OF-RANGE INSERTION LENGTHS AT FAILURE ===")
    for L, n in sorted(inslen.items()):
        if n: print(f"  {L:>3} bp  {n:>10,}")

    print("\n=== BACKWARD-PASS STOPS ===")
    for k, n in bwd.most_common():
        print(f"{k:24s} {n:>10,}")

    with open(os.path.join(OUT, f"{tag}.c8r_reasons.tsv"), "w") as f:
        f.write("category\tkey\tcount\n")
        for k, n in fwd.items(): f.write(f"fwd_stop\t{k}\t{n}\n")
        for (k, s), n in fwd_site.items(): f.write(f"fwd_stop_site\t{k}|{s}\t{n}\n")
        for g, n in gaps.items(): f.write(f"gap_at_failure\t{g}\t{n}\n")
        for L, n in arraylen.items(): f.write(f"array_len\t{L}\t{n}\n")
        for L, n in arraylen_res.items(): f.write(f"array_len_resolved\t{L}\t{n}\n")
        for L, n in arraylen_clean.items(): f.write(f"array_len_clean\t{L}\t{n}\n")
        for (fam, hit), n in anchor_h1.items(): f.write(f"gaplong_anchor_h1\t{fam}|{hit}\t{n}\n")
        for (L, sm), n in s56.items(): f.write(f"site56\t{L}|{'same' if sm else 'diff'}\t{n}\n")
        for s, n in nonvocab.items(): f.write(f"nonvocab_symbol\t{s}\t{n}\n")
        for k, n in bwd.items(): f.write(f"bwd_stop\t{k}\t{n}\n")
        for (b, L), n in bc_len.items(): f.write(f"bc_array_len\t{b}|{L}\t{n}\n")
        for (b, L, sm), n in bc_s56.items(): f.write(f"bc_site56\t{b}|{L}|{'same' if sm else 'diff'}\t{n}\n")
    print(f"\nwrote {tag}.c8r_reasons.tsv")


if __name__ == "__main__":
    a = sys.argv[1:]
    lim = None
    if "--limit" in a:
        i = a.index("--limit"); lim = int(a[i + 1]); a = a[:i] + a[i + 2:]
    main(a[0], lim)
