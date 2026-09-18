#!/usr/bin/env python3
"""KNOWN ISSUE (open, not fixed): backward-pass frame shift, and why single cell
is largely protected from it.

`tape/parser.walk_robust` resolves the monomer array twice -- forward from
anchors[0] and backward from the terminal landmark. The two number the sites from
opposite ends of the anchor list, so they agree only when that list holds exactly
N_SITES entries ahead of the terminal. If the FIRST anchor fails to match (one
mismatch in the 6-bp ANCHOR does it), the forward pass runs to completion one
site out of frame and the backward pass fills the vacated slot from the junction
the forward pass just consumed: monomer 1's symbol is lost and the last symbol
appears twice. A partially edited array shifts the same way with no duplicate.

The parser is shared VERBATIM with the bulk pipeline, so the defect is identical;
what differs is what the single-cell funnel then does with a shifted READ:

  1. UMI consensus       a shifted read must be the dominant read of its molecule
  2. >=2 molecules       tape/denoise.py drops any genotype with one molecule, so
                         the same shift must recur in two independent molecules of
                         the same (cell, integration)
  3. >=3 molecules and   a shifted genotype branches at depth 0 rather than being a
     >=0.9 dominance     prefix, so when it does not win it lowers dominance and the
                         locus is left MISSING rather than miscalled

Measured upper bound in the published callset: a 5'-shifted call loses monomer 1,
whose symbol is the blastomere founder allele (>=96.9% of recovered calls on a
side carry that side's allele). Calls whose monomer-1 symbol is NEITHER side's
founder allele therefore cap the leak: 10,783 of 8,311,724 testable calls =
0.130%, or ~0.16% after correcting for the 17.2% of cases in which a shift would
be invisible because the true monomer-2 symbol equals the founder allele. That is
an upper bound -- it also contains ambient bleed and every other monomer-1
discordance. Compare a ~0.5% read-level excess in bulk.
`scripts/frame_shift_bound.py` re-derives it from the per-blastomere matrices.

This file PINS THE CURRENT (WRONG) parser behaviour so the defect stays visible
and any future fix has a repro, and asserts the three suppression steps above.

    python3 tests/test_known_issue_frame_shift.py
"""
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config                                                      # noqa: E402
from tape.denoise import denoise_barcode                           # noqa: E402
from tape.parser import pattern_of, walk_robust                    # noqa: E402
from tape.sc import (chain_to_str, collapse_umis, consensus_locus,  # noqa: E402
                     extract_bc_sc)

BC = b"CGGGGAATTGTA"                            # a real embryo-#3 integration barcode
SYMS = [b"ACT", b"GCC", b"AAG", b"CCC", b"ACG", b"GAA"]
BROKEN_ANCHOR = b"GGTGTGCACG"                   # GGTGAGCACG with one substitution
FAILURES = []


def sc_read(n_edited=6, break_anchor=None, with_terminal=True):
    """Synthesise one single-cell amplicon read (starts at the barcode)."""
    array = b""
    for i in range(config.N_SITES):
        body = BROKEN_ANCHOR if i == break_anchor else config.MONO
        junction = (SYMS[i] + b"GGA" + b"TGAT") if i < n_edited else b"TGAT"
        array += body + junction
    if with_terminal:
        array += config.TERM_LONG
    return (BC + config.BC_RIGHT_SPACER + config.BC_RIGHT + config.LEADIN
            + array + b"CCCCCCCCCC")


def parse(read):
    _bc, cstart = extract_bc_sc(read)
    state, err, resolved_all = walk_robust(read, cstart)
    genotype = tuple("?" if s is None else ("U" if s == "U" else s[1].decode())
                     for s in state)
    return genotype, err, resolved_all


def check(label, got, want, should_be=None):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {label}")
    print(f"        current  {got}")
    if should_be is not None:
        print(f"        correct  {should_be}   <-- what a fix must produce")
    if not ok:
        print(f"        expected {want}")
        FAILURES.append(label)


TRUTH = tuple(s.decode() for s in SYMS)
SHIFTED = ("GCC", "AAG", "CCC", "ACG", "GAA", "GAA")
print(f"\nground truth for a fully edited array: {TRUTH}\n")

print("control -- intact read parses correctly")
check("intact, fully edited", parse(sc_read())[0], TRUTH)

print("\nTHE DEFECT -- first anchor lost, read reaches the terminal")
got, err, resolved = parse(sc_read(break_anchor=0))
check("monomer 1 dropped, monomer 6 duplicated", got, SHIFTED, should_be=TRUTH)
check("...and the corrupted read passes the ordering check", err, False)
check("...and passes all-sites-resolved, so it is counted VALID", resolved, True)

print("\nsilent variant -- partially edited array shifts with no duplicate to betray it")
check("only monomers 1-4 edited, first anchor lost",
      parse(sc_read(n_edited=4, break_anchor=0))[0],
      ("GCC", "AAG", "CCC", "U", "U", "U"),
      should_be=("ACT", "GCC", "AAG", "CCC", "U", "U"))

print("\ncases that are already SAFE (the read is discarded, not corrupted)")
for label, kwargs in [("interior anchor lost (a3)", dict(break_anchor=2)),
                      ("first anchor lost, no terminal landmark",
                       dict(break_anchor=0, with_terminal=False))]:
    _got, _err, resolved = parse(sc_read(**kwargs))
    check(f"{label} -> not all monomers resolved", resolved, False)

print("\nWHY SINGLE CELL IS PROTECTED -- the molecule funnel, on the shifted genotype")

# step 1: the shifted read has to win its own molecule's consensus
check("a shifted read that is a minority within its UMI is voted out",
      collapse_umis({b"AAAAAAAA": {TRUTH: 3, SHIFTED: 1}}), Counter({TRUTH: 1}))

# step 2: one shifted molecule is removed before consensus is even attempted
check("a single shifted molecule is dropped by the >=2-molecule floor",
      denoise_barcode(Counter({TRUTH: 8, SHIFTED: 1})), {TRUTH: 8})

# step 3: two shifted molecules survive de-noising but only lower dominance
kept = denoise_barcode(Counter({TRUTH: 20, SHIFTED: 2}))
check("two shifted molecules do survive de-noising (not a 1-symbol neighbour)",
      set(kept), {TRUTH, SHIFTED})
res = consensus_locus(Counter(kept))
check("...at 20:2 (dominance 0.91) the call is still the true genotype, flagged branched",
      (chain_to_str(res[0]), res[1]), ("-".join(TRUTH), True))
check("...but at 6:2 (dominance 0.75) the locus is left MISSING, not miscalled",
      consensus_locus(Counter({TRUTH: 6, SHIFTED: 2})), None)
check("a shifted genotype is CALLED only if it dominates the locus outright",
      chain_to_str(consensus_locus(Counter({SHIFTED: 9}))[0]), "-".join(SHIFTED))

print()
if FAILURES:
    print(f"{len(FAILURES)} check(s) did not match the pinned behaviour -- the parser or "
          f"the funnel changed. If that was a deliberate fix, update this file.")
    sys.exit(1)
print("known issue reproduced exactly as documented (still unfixed), and the "
      "single-cell suppression steps behave as described")
