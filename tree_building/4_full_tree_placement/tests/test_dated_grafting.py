#!/usr/bin/env python3
"""Checks for the dated caterpillar emitter (dated_tree.emit_dated).

Backbone is a tiny DATED tree (days, ultrametric to present=3.0). Verifies:
  * time-scaled pendants place a query's divergence node at present - pendant_days;
  * the grafted tree stays ultrametric (every named tip at the present) with no
    negative branches;
  * genotype-identical attachers (pendant 0) form a POLYTOMY -- a node with >=3 members
    and NO 0-length INTERNAL edge (terminal 0-length branches are allowed);
  * a pendant exceeding the anchor's edge budget is CLAMPED and reported.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from dated_tree import emit_dated
from newick import parse_structure

PRESENT = 3.0
# ultrametric dated backbone: A,B under a node at depth 2.0; C direct; all tips at 3.0
BACKBONE = "((A:1.0,B:1.0):2.0,C:3.0);"


def _write_backbone(tmp):
    p = os.path.join(tmp, "bb.nwk")
    open(p, "w").write(BACKBONE)
    return p


def _depths(nwk):
    """Return (name->depth-from-root) for named tips and the parse arrays."""
    parent, name, kids, blen = parse_structure(nwk)
    depth = [0.0] * len(parent)
    # nodes are created parent-before-child in parse_structure, so a forward pass works
    for v in range(1, len(parent)):
        b = float(blen[v]) if blen[v] else 0.0
        depth[v] = depth[parent[v]] + b
    tipdepth = {name[v]: depth[v] for v in range(len(parent)) if not kids[v] and name[v]}
    return tipdepth, parent, name, kids, blen, depth


def _emit(tmp, cells_ids, isbb, parent, pendant_days):
    bb = _write_backbone(tmp)
    query_idx = np.where(~np.array(isbb))[0]
    return emit_dated(np.array(cells_ids, dtype=object), np.array(isbb),
                      np.array(parent), set(), query_idx, bb, np.array(pendant_days, float))


def test_cherry_time_scaling_and_ultrametric(tmp_path=None):
    tmp = str(tmp_path) if tmp_path is not None else _mktmp()
    # Q (idx3) attaches to A (idx0) with pendant 0.5 days -> divergence node at 3.0-0.5=2.5
    cells = ["A", "B", "C", "Q"]
    isbb = [True, True, True, False]
    parent = [-1, -1, -1, 0]
    pend = [0.0, 0.0, 0.0, 0.5]
    nwk, placed, maxch, clamped = _emit(tmp, cells, isbb, parent, pend)
    assert placed == 1 and len(clamped) == 0
    tipdepth, *_ = _depths(nwk)
    for t in ("A", "B", "C", "Q"):
        assert abs(tipdepth[t] - PRESENT) < 1e-6, (t, tipdepth[t])   # ultrametric to present
    # Q's terminal branch == its pendant (0.5): Q sits 0.5 below its divergence node
    _, parent_a, name_a, kids_a, blen_a, depth_a = _depths(nwk)
    qv = [v for v in range(len(name_a)) if name_a[v] == "Q"][0]
    assert abs(float(blen_a[qv]) - 0.5) < 1e-6
    assert abs(depth_a[parent_a[qv]] - (PRESENT - 0.5)) < 1e-6   # node age = present - pendant
    print("  ok  test_cherry_time_scaling_and_ultrametric")


def test_identical_group_is_polytomy_no_internal_zero(tmp_path=None):
    tmp = str(tmp_path) if tmp_path is not None else _mktmp()
    # Q1 distinct (0.5); Q2,Q3 identical to A (pendant 0) -> polytomy with A, no internal 0-edge
    cells = ["A", "B", "C", "Q1", "Q2", "Q3"]
    isbb = [True, True, True, False, False, False]
    parent = [-1, -1, -1, 0, 0, 0]
    pend = [0, 0, 0, 0.5, 0.0, 0.0]
    nwk, placed, maxch, clamped = _emit(tmp, cells, isbb, parent, pend)
    assert placed == 3 and len(clamped) == 0
    tipdepth, parent_a, name_a, kids_a, blen_a, depth_a = _depths(nwk)
    for t in ("A", "B", "C", "Q1", "Q2", "Q3"):
        assert abs(tipdepth[t] - PRESENT) < 1e-6, (t, tipdepth[t])
    # no NEGATIVE branches
    for v in range(1, len(blen_a)):
        assert (float(blen_a[v]) if blen_a[v] else 0.0) >= -1e-9
    # NO 0-length INTERNAL edge (internal node = has children); terminal 0-edges allowed
    for v in range(1, len(blen_a)):
        b = float(blen_a[v]) if blen_a[v] else 0.0
        if kids_a[v]:                                   # internal node
            assert b > 1e-9, f"internal node {v} has 0-length edge"
    # A, Q2, Q3 share ONE parent node with >=3 members (the polytomy)
    def par_of(nm): return parent_a[[v for v in range(len(name_a)) if name_a[v] == nm][0]]
    pa, p2, p3 = par_of("A"), par_of("Q2"), par_of("Q3")
    assert pa == p2 == p3, (pa, p2, p3)
    assert len(kids_a[pa]) >= 3
    print("  ok  test_identical_group_is_polytomy_no_internal_zero")


def test_clamp_when_pendant_exceeds_edge(tmp_path=None):
    tmp = str(tmp_path) if tmp_path is not None else _mktmp()
    # Q attaches to A but pendant 5.0 > A's terminal edge (1.0 day) -> clamp to the edge
    cells = ["A", "B", "C", "Q"]
    isbb = [True, True, True, False]
    parent = [-1, -1, -1, 0]
    pend = [0, 0, 0, 5.0]
    nwk, placed, maxch, clamped = _emit(tmp, cells, isbb, parent, pend)
    assert 3 in clamped, clamped                        # Q (idx3) recorded as clamped
    tipdepth, parent_a, name_a, kids_a, blen_a, depth_a = _depths(nwk)
    for t in ("A", "B", "C", "Q"):
        assert abs(tipdepth[t] - PRESENT) < 1e-6, (t, tipdepth[t])
    # Q's divergence node clamped to A's parent depth (2.0): edge above it == 0
    qv = [v for v in range(len(name_a)) if name_a[v] == "Q"][0]
    assert abs(depth_a[parent_a[qv]] - 2.0) < 1e-6      # not present-5.0 (=negative), clamped to 2.0
    print("  ok  test_clamp_when_pendant_exceeds_edge")


def _mktmp():
    import tempfile
    return tempfile.mkdtemp(prefix="dated_test_")


def run():
    test_cherry_time_scaling_and_ultrametric()
    test_identical_group_is_polytomy_no_internal_zero()
    test_clamp_when_pendant_exceeds_edge()
    print("test_dated_grafting: ALL PASS")


if __name__ == "__main__":
    run()
