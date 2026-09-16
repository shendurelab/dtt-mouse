#!/usr/bin/env python3
"""Hand-worked checks for the matrix-metric pendant d(C,P) (dtt_lengths.pendants_from_code).

Encoding (genotypes.py): -1 missing tape, 0 unedited site, >=1 a specific edit token.
Verifies: d(C,P) counts C's private edits ONLY over tapes shared with T; additivity
d(C,P)+d(T,P) == d(C,T); the missing-tape and identical-genotype cases.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from dtt_lengths import pendants_from_code


def _code(*cells):
    """Stack per-cell (G,S) int lists into an (N,G,S) int16 array."""
    return np.array(cells, dtype=np.int16)


def test_extra_edit_on_shared_tape():
    # C = T plus one further edit (e3) on tape0; identical tape1.
    T = [[1, 2, 0], [1, 0, 0]]
    C = [[1, 2, 3], [1, 0, 0]]
    code = _code(C, T)                       # idx0=C (query), idx1=T (parent)
    parent = np.array([1, -1])
    pend, resid = pendants_from_code(code, parent)
    # tape0: depth_C 3, prefix 2 -> 1 ; tape1: depth_C 1, prefix 1 -> 0 ; /ns=2
    assert abs(pend[0] - 0.5) < 1e-12, pend[0]
    assert abs(resid[0] - 0.0) < 1e-12, resid[0]
    # additivity vs matrix d(C,T): tape0 3+2-2*2=1, tape1 1+1-2*1=0 -> mean 0.5
    assert abs((pend[0] + resid[0]) - 0.5) < 1e-12


def test_private_edit_on_nonshared_tape_ignored():
    # C has edits on tape1 that T is MISSING -> excluded from the shared-tape average.
    T = [[1, 0, 0], [-1, -1, -1]]
    C = [[1, 0, 0], [1, 2, 0]]
    code = _code(C, T)
    parent = np.array([1, -1])
    pend, resid = pendants_from_code(code, parent)
    assert abs(pend[0] - 0.0) < 1e-12, pend[0]   # only shared tape0 counts, where C==T
    assert abs(resid[0] - 0.0) < 1e-12


def test_identical_genotypes_zero_pendant():
    G = [[1, 2, 0], [1, 0, 0]]
    code = _code(G, G)
    parent = np.array([1, -1])
    pend, resid = pendants_from_code(code, parent)
    assert pend[0] == 0.0 and resid[0] == 0.0


def test_divergent_edit_counts_as_private():
    # both edited site1 but differently: shared prefix ends at site0; each carries a
    # private edit below the ancestor -> d(C,P)=d(T,P)=0.5 (one private edit / 1 tape... here /1? ns=1 tape)
    T = [[1, 2, 0]]
    C = [[1, 5, 0]]
    code = _code(C, T)
    parent = np.array([1, -1])
    pend, resid = pendants_from_code(code, parent)
    # tape0: depth_C 2, depth_T 2, prefix 1 -> d(C,P)=1, d(T,P)=1 ; ns=1
    assert abs(pend[0] - 1.0) < 1e-12, pend[0]
    assert abs(resid[0] - 1.0) < 1e-12, resid[0]
    assert abs((pend[0] + resid[0]) - 2.0) < 1e-12   # d(C,T)=2+2-2*1=2


def run():
    for fn in (test_extra_edit_on_shared_tape, test_private_edit_on_nonshared_tape_ignored,
               test_identical_genotypes_zero_pendant, test_divergent_edit_counts_as_private):
        fn(); print(f"  ok  {fn.__name__}")
    print("test_dtt_pendant: ALL PASS")


if __name__ == "__main__":
    run()
