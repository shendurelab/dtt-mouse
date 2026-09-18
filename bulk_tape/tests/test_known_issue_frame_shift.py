#!/usr/bin/env python3
"""KNOWN ISSUE (open, not fixed): backward-pass symbol duplication.

`tape/parser.walk_robust` resolves the monomer array twice -- forward from
anchors[0] and backward from the terminal landmark. The two number the sites
from opposite ends of the anchor list, so they agree only when that list holds
exactly N_SITES entries ahead of the terminal.

If the FIRST anchor fails to match (one mismatch in the 6-bp ANCHOR does it),
every remaining gap is still valid, so the forward pass runs to completion one
site out of frame and leaves the last slot empty; the backward pass then fills
that slot from the same a6->TERM junction the forward pass just consumed. Site
1's symbol is lost and the last symbol appears twice.

This file PINS THE CURRENT (WRONG) BEHAVIOUR so the defect stays visible and any
future fix has a repro to work against. It passes today by asserting what the
parser does now, not what it should do; the `should_be` values are the target.

    python3 tests/test_known_issue_frame_shift.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config                                    # noqa: E402
from tape.barcode import extract_bc              # noqa: E402
from tape.parser import walk_robust, pattern_of  # noqa: E402

BC = b"CGGGGAATTGTA"                             # a real embryo-#3 integration barcode
SYMS = [b"ACT", b"GCC", b"AAG", b"CCC", b"ACG", b"GAA"]
BROKEN_ANCHOR = b"GGTGTGCACG"                    # GGTGAGCACG with one substitution
FAILURES = []


def build_read(n_edited=6, break_anchor=None, with_terminal=True):
    """Synthesise one amplicon read. `break_anchor` is a 0-based monomer index."""
    array = b""
    for i in range(config.N_SITES):
        body = BROKEN_ANCHOR if i == break_anchor else config.MONO
        junction = (SYMS[i] + b"GGA" + b"TGAT") if i < n_edited else b"TGAT"
        array += body + junction
    if with_terminal:
        array += config.TERM_LONG
    return (b"AAAA" + config.BC_LEFT + BC + config.BC_RIGHT_SPACER
            + config.BC_RIGHT + config.LEADIN + array + b"CCCCCCCCCC")


def parse(read):
    _bc, cstart = extract_bc(read)
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
print(f"\nground truth for a fully edited array: {TRUTH}\n")

print("control -- intact read parses correctly")
check("intact, fully edited", parse(build_read())[0], TRUTH)

print("\nTHE DEFECT -- first anchor lost, read reaches the terminal")
got, err, resolved = parse(build_read(break_anchor=0))
check("site 1 dropped, site 6 duplicated",
      got, ("GCC", "AAG", "CCC", "ACG", "GAA", "GAA"), should_be=TRUTH)
check("...and the corrupted read passes the ordering check", err, False)
check("...and passes all-sites-resolved, so it is counted VALID", resolved, True)

print("\nsilent variant -- partially edited array shifts with no duplicate to betray it")
check("only sites 1-4 edited, first anchor lost",
      parse(build_read(n_edited=4, break_anchor=0))[0],
      ("GCC", "AAG", "CCC", "U", "U", "U"),
      should_be=("ACT", "GCC", "AAG", "CCC", "U", "U"))

print("\ncases that are already SAFE (the read is discarded, not corrupted)")
for label, kwargs in [("interior anchor lost (a3)", dict(break_anchor=2)),
                      ("first anchor lost, no terminal landmark",
                       dict(break_anchor=0, with_terminal=False))]:
    _got, _err, resolved = parse(build_read(**kwargs))
    check(f"{label} -> not all sites resolved", resolved, False)

print()
if FAILURES:
    print(f"{len(FAILURES)} check(s) did not match the pinned behaviour -- the parser "
          f"changed. If that was a deliberate fix, update this file.")
    sys.exit(1)
print("known issue reproduced exactly as documented (still unfixed)")
