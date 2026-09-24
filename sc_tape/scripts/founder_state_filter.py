#!/usr/bin/env python3
"""Founder-state filter -- the last cell filter before tree building.

The Methods state: "At both stages we further required that a cell carry no tape
genotype contradicting its blastomere's founder state, i.e. the near-fixed prefix
of edits written before the first cleavage."

This script derives that founder state per (integration, blastomere) and flags
every cell that contradicts it. For each integration it starts from the
blastomere-defining monomer-1 allele (tables/e3v8.defining_sites.tsv) and extends
one monomer at a time while the next monomer is edited and near-fixed (>= FIX,
default 0.99) among the cells that match the chain so far. A cell FAILS if any of
its recovered genotypes disagrees with its side's chain over the chain's length --
including an unedited monomer where the founder wrote an edit, since a founder
edit must be inherited by every descendant.

    python3 scripts/founder_state_filter.py <dir with e3v8.B{1,2}_tape_consensus.tsv[.gz]> \
            [--fix 0.99] [--defs tables/e3v8.defining_sites.tsv] [--out DIR]

Writes <side>.founder_state.tsv (cell_id, n_contradicting, founder_ok) per side and
prints the counts, split by the two tree-input bands (>=7 tapes, 4-6 tapes).

PROVENANCE -- READ THIS. The filter that produced the published trees was applied
upstream of the tree pipeline, which consumes files named
`*_tape_consensus.ge7_founderok.tsv.gz`; that code lives at
tree_building/1_build_nj_backbone/00_filter_highqual_consensus.R (+
00a_reshape_v8_founder_genotypes.R), NOT this script. This
script is a RECONSTRUCTION, validated against the published tip sets: on the v8
matrices it flags 4,195 of 4,195 (blastomere A) and 3,729 of 3,889 (B) of the
backbone-eligible cells that are absent from the 655,701-tip backbone, and 4,776
of 4,833 (A) and 2,622 of 2,741 (B) of the placement-eligible cells absent from
the 625,440 placed, with ZERO false positives in all four sets -- i.e. it never
flags a cell that is in a published tree. The 336 cells (2% of 15,658) it misses
are presumably a slightly looser near-fixed threshold than the 0.99 used here.
Use it to audit the published cell sets, not as a substitute for the original.
"""
import argparse
import collections
import csv
import gzip
import sys
from pathlib import Path

QC = {"n_loci", "n_doublet_loci", "mean_dominance", "pass_qc"}
ROOT = Path(__file__).resolve().parent.parent


def open_maybe_gz(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def find(d, side):
    for name in (f"e3v8.{side}_tape_consensus.tsv.gz", f"e3v8.{side}_tape_consensus.tsv",
                 f"sub.{side}_tape_consensus.tsv.gz", f"sub.{side}_tape_consensus.tsv",
                 f"{side}_tape_consensus.tsv.gz", f"{side}_tape_consensus.tsv"):
        p = Path(d) / name
        if p.exists():
            return p
    sys.exit(f"no {side} matrix in {d}")


def load(path):
    with open_maybe_gz(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        hdr = next(r)
        bcs = [c for c in hdr[1:] if c not in QC]
        idx = [hdr.index(b) for b in bcs]
        i_nl, i_pq = hdr.index("n_loci"), hdr.index("pass_qc")
        rows = [(row[0], int(row[i_nl]), int(row[i_pq]), tuple(row[j] for j in idx))
                for row in r]
    return bcs, rows


def derive_chains(bcs, rows, seed, fix, max_depth=6):
    """Extend each integration's founder chain while the next monomer is near-fixed."""
    chain = {b: [seed[b]] for b in bcs if b in seed}
    for depth in range(1, max_depth):
        grew = False
        for b in list(chain):
            ch = chain[b]
            if len(ch) != depth:
                continue
            j = bcs.index(b)
            cnt = collections.Counter()
            for _c, _nl, pq, toks in rows:
                if not pq:
                    continue
                v = toks[j]
                if v == "NA":
                    continue
                g = v.split("|")
                if g[:depth] == ch:
                    cnt[g[depth]] += 1
            t = sum(cnt.values())
            if not t:
                continue
            a, n = cnt.most_common(1)[0]
            if a != "U" and n / t >= fix:
                chain[b] = ch + [a]
                grew = True
        if not grew:
            break
    return chain


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix_dir")
    ap.add_argument("--defs", default=str(ROOT / "tables" / "e3v8.defining_sites.tsv"))
    ap.add_argument("--fix", type=float, default=0.99,
                    help="near-fixed fraction required to extend a founder chain [0.99]")
    ap.add_argument("--out", default=None, help="write per-cell flags here [matrix_dir]")
    a = ap.parse_args()
    out = Path(a.out or a.matrix_dir)

    defs = collections.defaultdict(dict)
    with open(a.defs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            defs[row["integration"]][row["blastomere"]] = row["allele"]

    for side in ("B1", "B2"):
        bcs, rows = load(find(a.matrix_dir, side))
        seed = {b: d[side] for b, d in defs.items() if side in d}
        chain = derive_chains(bcs, rows, seed, a.fix)
        print(f"--- {side}: founder chains (>= {a.fix} fixed at each extension)")
        for b, ch in sorted(chain.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            print(f"      {b}  {'-'.join(ch)}")
        tally = collections.Counter()
        path = out / f"{side}.founder_state.tsv"
        with open(path, "w") as fh:
            fh.write("cell_id\tn_contradicting\tfounder_ok\n")
            for cell, nl, pq, toks in rows:
                n = 0
                for b, v in zip(bcs, toks):
                    if v == "NA" or b not in chain:
                        continue
                    ch = chain[b]
                    if v.split("|")[:len(ch)] != ch:
                        n += 1
                ok = int(n == 0)
                fh.write(f"{cell}\t{n}\t{ok}\n")
                if pq:
                    band = "ge7" if nl >= 7 else ("4to6" if nl >= 4 else "lt4")
                    tally[(band, "total")] += 1
                    tally[(band, "fail")] += 0 if ok else 1
        for band in ("ge7", "4to6"):
            tot, fail = tally[(band, "total")], tally[(band, "fail")]
            if tot:
                print(f"   pass_qc & {band:>4} tapes: {tot:>8,} eligible, "
                      f"{fail:>6,} contradict the founder state ({100*fail/tot:.2f}%)")
        print(f"   wrote {path}")


if __name__ == "__main__":
    main()
