#!/usr/bin/env python3
"""Upper bound on the parser frame-shift leak into single-cell genotype calls.

See tests/test_known_issue_frame_shift.py for the defect itself. A 5'-shifted
call loses monomer 1, whose symbol is the blastomere founder allele, so calls
whose monomer-1 symbol is NEITHER side's founder allele cap the leak. The bound
is an over-estimate: it also contains ambient bleed and every other monomer-1
discordance. It is blind where the true monomer-2 symbol equals the side's own
founder allele (the same symbol written twice), and that fraction is measured
here and used to correct the bound.

Usage:
    python3 scripts/frame_shift_bound.py <dir-with-e3v8.B1/B2_tape_consensus.tsv[.gz]> \
                                         [--defs tables/e3v8.defining_sites.tsv]

Published values (full callset): 10,783 of 8,311,724 testable calls = 0.130%;
17.2% of shifts invisible -> corrected bound 0.157%.
"""
import argparse
import collections
import csv
import gzip
import sys
from pathlib import Path

QC = {"n_loci", "n_doublet_loci", "mean_dominance", "pass_qc"}


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def find(d, side):
    for name in (f"e3v8.{side}_tape_consensus.tsv.gz", f"e3v8.{side}_tape_consensus.tsv",
                 f"{side}_tape_consensus.tsv.gz", f"{side}_tape_consensus.tsv"):
        p = Path(d) / name
        if p.exists():
            return p
    sys.exit(f"no {side} matrix found in {d}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix_dir")
    ap.add_argument("--defs", default=str(Path(__file__).resolve().parent.parent
                                          / "tables" / "e3v8.defining_sites.tsv"))
    a = ap.parse_args()

    defs = collections.defaultdict(dict)
    with open(a.defs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            defs[row["integration"]][row["blastomere"]] = row["allele"]

    m2 = collections.defaultdict(collections.Counter)   # (bc,side) -> true monomer-2 usage
    n_ed = collections.Counter()
    n_off = collections.Counter()
    n_all = n_u1 = 0
    for side in ("B1", "B2"):
        with open_maybe_gz(find(a.matrix_dir, side)) as fh:
            r = csv.reader(fh, delimiter="\t")
            hdr = next(r)
            bcs = [c for c in hdr[1:] if c not in QC]
            idx = [hdr.index(b) for b in bcs]
            for row in r:
                for b, j in zip(bcs, idx):
                    v = row[j]
                    if v == "NA":
                        continue
                    g = v.split("|")
                    n_all += 1
                    if g[0] == "U":
                        n_u1 += 1
                        continue
                    own = defs[b].get(side)
                    if own is None:            # no defining allele on this side
                        continue
                    n_ed[(b, side)] += 1
                    if g[0] == own:
                        m2[(b, side)][g[1]] += 1
                    else:
                        n_off[(b, side)] += 1

    E, O = sum(n_ed.values()), sum(n_off.values())
    if not E:
        sys.exit("no testable calls found -- wrong matrices or defining-site file?")
    print(f"non-NA calls                     {n_all:>12,}")
    print(f"  monomer-1 unedited             {n_u1:>12,}  ({100*n_u1/n_all:.3f}%)")
    print(f"testable (side has a defining allele) {E:>7,}")
    print(f"off-founder monomer-1            {O:>12,}  ({100*O/E:.3f}%)")
    print(f"  per integration x side range   "
          f"{100*min(n_off[k]/n_ed[k] for k in n_ed):.3f}% - "
          f"{100*max(n_off[k]/n_ed[k] for k in n_ed):.3f}%")
    num = den = 0.0
    for k in n_ed:
        b, side = k
        t = sum(m2[k].values())
        if t:
            num += (m2[k][defs[b][side]] / t) * n_ed[k]
            den += n_ed[k]
    p = num / den
    print(f"shift invisible (true m2 == own founder allele) {100*p:.2f}% of calls")
    print(f"CORRECTED UPPER BOUND            {100*O/E/(1-p):.3f}% of calls")
    print("(upper bound only: also contains ambient bleed and any other monomer-1 "
          "discordance)")


if __name__ == "__main__":
    main()
