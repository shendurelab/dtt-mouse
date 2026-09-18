#!/usr/bin/env python3
"""
Validate the v4 blastomere routing against the recommendation's acceptance checks.

Inputs:
  --labels   <out>.routing_labels.tsv   (from blastomere_route.py)
  --treeclade FILE   two-column cell_id<TAB>{A|B} of the CURRENT merged-tree root
                     clades (the no-fix baseline). If a .nwk is given instead, its
                     root's two largest subclades are used.
  --meta     cell_metadata.v3.txt        (cell_id, major_trajectory, celltype)

Checks (from RECOMMENDATIONS_pipeline.md):
  1. Root-clade purity: do the tree's two root clades map cleanly onto B1/B2, and are
     the "misplaced" conflict cells now routed to the correct side (0 wrong-direction)?
  2. B1:B2 balance preserved vs the no-fix tree.
  3. Set-aside pile: small AND fate-flat (else keep as low-confidence, don't drop).
  4. Per-trajectory fate composition of B1 vs B2 preserved (not skewed by routing).
"""
import argparse, csv, sys, re
from collections import Counter, defaultdict


def load_labels(path):
    rows = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            rows[d["cell"]] = d
    return rows


def load_treeclade(path):
    if path.endswith(".nwk") or path.endswith(".newick"):
        return _parse_tree_clades(path)
    tc = {}
    with open(path) as f:
        for line in f:
            c, lab = line.rstrip("\n").split("\t")
            tc[c] = lab
    return tc


def _parse_tree_clades(path):
    sys.setrecursionlimit(50_000_000)
    s = open(path).read().strip()
    NUM = r'(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)'
    pos = 0
    def parse():
        nonlocal pos
        ch = []
        if s[pos] == '(':
            pos += 1
            while True:
                ch.append(parse())
                if s[pos] == ',': pos += 1; continue
                if s[pos] == ')': pos += 1; break
        m = re.match(r'([^(),:;]*)', s[pos:]); lab = m.group(1); pos += m.end()
        if pos < len(s) and s[pos] == ':':
            m = re.match(':'+NUM, s[pos:]); pos += m.end()
        return [lab, ch]
    root = parse()
    def tips(n):
        out, st = [], [n]
        while st:
            x = st.pop()
            (out if not x[1] else st).extend([x[0]] if not x[1] else x[1])
        return out
    kids = sorted(root[1], key=lambda c: -len(tips(c)))
    tc = {}
    for lab, k in (("A", kids[0]), ("B", kids[1])):
        for t in tips(k):
            if t: tc[t] = lab
    return tc


def load_meta(path):
    meta = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            meta[d["cell_id"]] = d.get("major_trajectory", "NA")
    return meta


def enrichment_table(cells, meta, ref_counts, ref_total, title):
    """Print fate composition of `cells` vs a reference distribution (enrichment)."""
    c = Counter(meta.get(x, "NA") for x in cells)
    tot = sum(c.values()) or 1
    print(f"\n{title}  (n={tot:,})")
    print(f"  {'trajectory':<28}{'this %':>8}{'ref %':>8}{'enrich':>8}")
    for k, _ in sorted(c.items(), key=lambda kv: -kv[1]):
        this = c[k] / tot * 100
        ref = ref_counts.get(k, 0) / ref_total * 100 if ref_total else 0
        en = this / ref if ref > 0 else float("inf")
        print(f"  {k[:27]:<28}{this:>7.1f}{ref:>8.1f}{en:>7.2f}x")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--treeclade")
    ap.add_argument("--meta")
    args = ap.parse_args()

    lab = load_labels(args.labels)
    cells = list(lab)
    b1 = [c for c in cells if lab[c]["blastomere"] == "B1"]
    b2 = [c for c in cells if lab[c]["blastomere"] == "B2"]
    aside = [c for c in cells if lab[c]["blastomere"] == ""]
    print("=" * 64)
    print("v4 routing validation")
    print("=" * 64)
    print(f"cells        : {len(cells):,}")
    print(f"routed B1    : {len(b1):,}")
    print(f"routed B2    : {len(b2):,}")
    print(f"set aside    : {len(aside):,} ({len(aside)/len(cells)*100:.2f}%)")
    print(f"B1:B2 ratio  : {len(b1)/(len(b1)+len(b2))*100:.1f} / "
          f"{len(b2)/(len(b1)+len(b2))*100:.1f}")

    # ---- 1. root-clade purity + misplaced fix ----
    if args.treeclade:
        tc = load_treeclade(args.treeclade)
        onA = {c for c in cells if tc.get(c) == "A"}
        onB = {c for c in cells if tc.get(c) == "B"}
        # orient: B1 <-> the clade it mostly contains
        b1s = set(b1); b2s = set(b2)
        a_b1 = len(b1s & onA); a_b2 = len(b2s & onA)
        B1clade, B2clade = ("A", "B") if a_b1 >= a_b2 else ("B", "A")
        cl = {"A": onA, "B": onB}
        print("\n--- check 1: root-clade purity / misplaced fix ---")
        print(f"tree clades matched : A={len(onA):,}  B={len(onB):,}")
        print(f"orientation         : B1<->tree{B1clade}, B2<->tree{B2clade}")
        mis_fixed = len(set(b1) & cl[B2clade])       # routed B1 but tree put in B2 clade
        wrong = len(set(b2) & cl[B1clade])            # routed B2 but tree put in B1 clade (should be ~0)
        pure_b1 = len(set(b1) & cl[B1clade])
        print(f"B1 in tree-{B1clade} (agree) : {pure_b1:,}")
        print(f"B1 in tree-{B2clade} (misplaced now fixed) : {mis_fixed:,}")
        print(f"B2 in tree-{B1clade} (wrong direction, want ~0) : {wrong:,}")
        aside_on_tree = sum(1 for c in aside if c in tc)
        print(f"tree cells set aside : {aside_on_tree:,}")

    # ---- 3/4. fate composition ----
    if args.meta:
        meta = load_meta(args.meta)
        allc = Counter(meta.get(c, "NA") for c in cells)
        tot = sum(allc.values())
        print("\n--- check 3: set-aside fate-flatness (enrichment vs all routed+aside) ---")
        # break aside by reason
        for reason in ("set_aside:no_founder", "set_aside:tie", "set_aside:multibleed"):
            sub = [c for c in aside if lab[c]["reason"] == reason]
            if sub:
                enrichment_table(sub, meta, allc, tot, f"set-aside [{reason}]")
        enrichment_table(aside, meta, allc, tot, "set-aside [ALL]")
        print("\n--- check 4: B1 vs B2 fate composition ---")
        enrichment_table(b1, meta, allc, tot, "routed B1")
        enrichment_table(b2, meta, allc, tot, "routed B2")


if __name__ == "__main__":
    main()
