#!/usr/bin/env python3
"""Checks for the DTT-ancestor nearest-neighbour search (dtt_match.best_matches_dtt),
built on the real inverted index (similarity.token_freq / build_index) -- the one
piece of core placement logic exercised only indirectly (via the R cross-check and
the end-to-end smoke test) elsewhere in this test suite.

Encoding (genotypes.py): -1 missing tape, 0 unedited site, >=1 a specific edit token.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from similarity import token_freq, build_index
from dtt_match import best_matches_dtt

STRIDE = 8192
N_SITES = 3


def _code(*cells):
    return np.array(cells, dtype=np.int16)


def _match(code, isbb, nloci, query_idx, G=1, cap_frac=1.0, topk=500):
    freq = token_freq(code, G, N_SITES, STRIDE)
    postings = build_index(code, G, N_SITES, STRIDE, freq, cap_frac)
    return best_matches_dtt(code, G, N_SITES, STRIDE, postings,
                            np.array(query_idx), np.array(isbb), np.array(nloci),
                            topk=topk)


def test_finds_correct_closest_match():
    # idx0 T_close: shares Q's whole 2-edit prefix -> d=1 ; idx1 T_far: shares only
    # 1 edit -> d=2. Q must attach to T_close, not T_far.
    T_close = [[1, 2, 0]]
    T_far = [[1, 0, 0]]
    Q = [[1, 2, 3]]
    code = _code(T_close, T_far, Q)
    isbb = [True, True, False]
    nloci = [1, 1, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[2])
    assert m[0] == 0, m[0]                # closest = T_close (idx0)
    assert abs(-s[0] - 1.0) < 1e-9, s[0]   # d(Q,T_close) = 3+2-2*2 = 1
    assert mb[0] == 0                      # backbone-restricted match agrees (only bb is idx0/1)
    print("  ok  test_finds_correct_closest_match")


def test_tie_break_prefers_more_shared_tapes():
    # A shares BOTH tapes with Q (d=0, nsh=2); B shares only tape0 (d=0, nsh=1).
    # Distance ties at 0 -- more shared tapes (more evidence) must win: A over B.
    A = [[1, 0, 0], [1, 0, 0]]
    B = [[1, 0, 0], [-1, -1, -1]]
    Q = [[1, 0, 0], [1, 0, 0]]
    code = _code(A, B, Q)
    isbb = [True, True, False]
    nloci = [1, 1, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[2], G=2)
    assert m[0] == 0, m[0]                 # A (idx0), not B (idx1)
    assert abs(-s[0]) < 1e-9
    print("  ok  test_tie_break_prefers_more_shared_tapes")


def test_tie_break_prefers_higher_nloci():
    # C and D are genotype-identical to Q (both d=0, both nsh=2) -- tie-break falls
    # to n_loci (best-quality anchor): D (nloci=10) must win over C (nloci=5).
    C = [[1, 0, 0], [1, 0, 0]]
    D = [[1, 0, 0], [1, 0, 0]]
    Q = [[1, 0, 0], [1, 0, 0]]
    code = _code(C, D, Q)
    isbb = [True, True, False]
    nloci = [5, 10, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[2], G=2)
    assert m[0] == 1, m[0]                 # D (idx1), not C (idx0)
    print("  ok  test_tie_break_prefers_higher_nloci")


def test_tie_break_prefers_lower_index():
    # E and F tie on distance AND n_loci -- last tie-break criterion is lowest
    # global index: E (idx0) must win over F (idx1).
    E = [[1, 0, 0], [1, 0, 0]]
    F = [[1, 0, 0], [1, 0, 0]]
    Q = [[1, 0, 0], [1, 0, 0]]
    code = _code(E, F, Q)
    isbb = [True, True, False]
    nloci = [7, 7, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[2], G=2)
    assert m[0] == 0, m[0]
    print("  ok  test_tie_break_prefers_lower_index")


def test_best_backbone_match_restricted_to_backbone():
    # Q's overall-closest cell is another query (P, d=0); its closest BACKBONE cell
    # is farther (T, d=1). m must be P; mb must be T, not P.
    T = [[1, 2, 0]]     # backbone, d(Q,T) = 3+2-2*2 = 1
    P = [[1, 2, 3]]     # query anchor, genotype-identical to Q -> d(Q,P) = 0
    Q = [[1, 2, 3]]
    code = _code(T, P, Q)
    isbb = [True, False, False]
    nloci = [1, 1, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[2])
    assert m[0] == 1, m[0]                 # overall closest = P (idx1), a query
    assert mb[0] == 0, mb[0]               # backbone-restricted closest = T (idx0)
    assert abs(-sb[0] - 1.0) < 1e-9, sb[0]
    print("  ok  test_best_backbone_match_restricted_to_backbone")


def test_no_shared_token_returns_unmatched():
    # Q shares no edit token with anything (all-unedited) -> no candidate at all.
    T = [[1, 2, 0]]
    Q = [[0, 0, 0]]
    code = _code(T, Q)
    isbb = [True, False]
    nloci = [1, 1]
    m, s, mb, sb = _match(code, isbb, nloci, query_idx=[1])
    assert m[0] == -1, m[0]
    assert mb[0] == -1, mb[0]
    print("  ok  test_no_shared_token_returns_unmatched")


def run():
    for fn in (test_finds_correct_closest_match, test_tie_break_prefers_more_shared_tapes,
               test_tie_break_prefers_higher_nloci, test_tie_break_prefers_lower_index,
               test_best_backbone_match_restricted_to_backbone,
               test_no_shared_token_returns_unmatched):
        fn()
    print("test_dtt_match: ALL PASS")


if __name__ == "__main__":
    run()
