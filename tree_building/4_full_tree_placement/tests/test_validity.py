#!/usr/bin/env python3
"""Checks for the irreversibility (is_valid_ancestor) diagnostic (validity.valid_flags).

invalid <=> some (tape,site) where T edited (>=1), C resolved (!= -1), C != T.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from validity import valid_flags


def _code(*cells):
    return np.array(cells, dtype=np.int16)


def test_divergent_edit_invalid():
    T = [[1, 2, 0]]
    C = [[1, 5, 0]]                       # C resolves site1 differently -> irreversible
    v = valid_flags(_code(C, T), np.array([1, -1]))
    assert v[0] == False


def test_missing_edit_invalid():
    T = [[1, 2, 0]]
    C = [[1, 0, 0]]                       # C lacks T's site1 edit (observed unedited) -> invalid
    v = valid_flags(_code(C, T), np.array([1, -1]))
    assert v[0] == False


def test_missing_tape_is_wildcard_valid():
    T = [[1, 2, 0], [1, 0, 0]]
    C = [[1, 2, 0], [-1, -1, -1]]         # C missing tape1 entirely -> "?" wildcard, valid
    v = valid_flags(_code(C, T), np.array([1, -1]))
    assert v[0] == True


def test_descendant_valid():
    T = [[1, 2, 0]]
    C = [[1, 2, 3]]                       # C = T + further edit -> valid descendant
    v = valid_flags(_code(C, T), np.array([1, -1]))
    assert v[0] == True


def run():
    for fn in (test_divergent_edit_invalid, test_missing_edit_invalid,
               test_missing_tape_is_wildcard_valid, test_descendant_valid):
        fn(); print(f"  ok  {fn.__name__}")
    print("test_validity: ALL PASS")


if __name__ == "__main__":
    run()
