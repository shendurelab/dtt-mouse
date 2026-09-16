#!/usr/bin/env python3
"""Emit the placed tree onto the DATED backbone, with time-scaled pendants.

The backbone here is a DATED per-side subtree (branch lengths in DAYS, ultrametric to
the present, e.g. 13.5). Every cell -- backbone or placed -- is a real cell sampled at
the present, so a placed query C must land AT the present with its own branch equal to
its time since diverging from its anchor:  pendant_days(C) = d(C,P) / clock_rate.

Grafting an anchor A (a backbone tip, or a placed query hosting further queries) that
has attached children q_1..q_k (pendant_days p_1..p_k): we cannot keep A's original
tip position and hang the children below it -- that would push them past the present.
Instead we spend A's edge budget E (its dated branch length) placing the children as a
caterpillar ordered by pendant DEPTH, so each child and A itself reach the present:

    order children by capped pendant DESC (c_1 >= ... >= c_k), c_i = min(p_i, E)
    node N_1 at depth-from-present c_1 (oldest divergence), ..., N_k at c_k (youngest)
      edge(parent -> N_1) = E - c_1
      edge(N_i -> N_{i+1}) = c_i - c_{i+1}
      edge(N_i -> child_i subtree) = c_i        (child at present; recurse with budget c_i)
      edge(N_k -> A's own leaf copy) = c_k       (A at present)

so child_i's branch back to its divergence node is exactly its pendant c_i, A stays at
the present, and all inserted node ages are monotone. CLAMP: c_i = min(p_i, E) keeps a
node from predating A's parent; a cell whose pendant exceeds the available edge is
capped (attached at A's parent depth) and recorded. k == 1 is the ordinary cherry.

POLYTOMY on tied pendant: children sharing a pendant depth carry no signal to order
them, so they share ONE node (a soft polytomy) rather than a 0-length-edge ladder. In
particular children with pendant 0 -- genotype-IDENTICAL to the anchor -- become direct
siblings of the anchor's own leaf copy: a polytomy of identical cells at the present.
Only terminal branches of truly identical cells stay length 0; no 0-length internal
edge is ever emitted.

Query->query chains recurse: a child that itself hosts queries is expanded with budget
= its own capped pendant. The whole walk is iterative (these are caterpillar trees with
hundreds of thousands of tips; recursion would overflow).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from newick import parse_structure

CLAMP_EPS = 1e-9


def _build_attach_kids(N, parent, held_set, query_idx):
    """attach_kids[c] = placed query cells whose resolved parent is c."""
    attach_kids = [[] for _ in range(N)]
    placed = 0
    for qi in query_idx:
        if qi in held_set:
            continue
        p = parent[qi]
        if p < 0:
            continue
        attach_kids[p].append(int(qi))
        placed += 1
    return attach_kids, placed


def _expand(anchor, budget, attach_kids, pendant_days, cells_ids, clamped):
    """Return the Newick fragment for `anchor`'s expansion (top branch included).

    A leaf (no attached queries) returns `label:budget`; otherwise the pendant-ordered
    caterpillar described in the module docstring. Iterative post-order.
    """
    def make_frame(node, bud):
        ks = attach_kids[node]
        caps = []
        for c in ks:
            p = float(pendant_days[c])
            cap = p if p <= bud else bud
            if p > bud + CLAMP_EPS:
                clamped.add(c)
            caps.append((c, cap))
        caps.sort(key=lambda x: -x[1])                       # deepest divergence first
        return {"node": node, "budget": bud, "caps": caps,
                "frags": [None] * len(caps), "i": 0}

    stack = [make_frame(anchor, budget)]
    ret = None
    while stack:
        fr = stack[-1]
        if not fr["caps"]:                                   # leaf: label:budget
            frag = f"{cells_ids[fr['node']]}:{fr['budget']:.6f}"
            stack.pop()
            if stack:
                par = stack[-1]
                par["frags"][par["i"] - 1] = frag
            else:
                ret = frag
            continue
        if fr["i"] < len(fr["caps"]):                        # descend into next child
            c, cap = fr["caps"][fr["i"]]
            fr["i"] += 1
            stack.append(make_frame(c, cap))
            continue
        # all children resolved -> assemble from the inside out, GROUPING children
        # that share a pendant depth into ONE node (a soft polytomy). This collapses
        # what would be 0-length internal edges; in particular children with pendant 0
        # (genotype-identical to the anchor) become direct siblings of the anchor's own
        # leaf copy -- a polytomy of identical cells at the present, not a resolved
        # ladder. caps is sorted DESC, so equal values are already adjacent.
        node, budget_, caps, frags = fr["node"], fr["budget"], fr["caps"], fr["frags"]
        label = cells_ids[node]
        levels = []                                          # [value, [frag indices]] deepest-first
        for i, (_c, cap) in enumerate(caps):
            key = round(cap, 9)
            if levels and levels[-1][0] == key:
                levels[-1][1].append(i)
            else:
                levels.append([key, [i]])
        # innermost node (smallest pendant): its own-depth children + the anchor leaf
        v_inner = levels[-1][0]
        members = [frags[i] for i in levels[-1][1]] + [f"{label}:{v_inner:.6f}"]
        s = "(" + ",".join(members) + ")"
        for L in range(len(levels) - 2, -1, -1):             # wrap out to larger pendants
            edge = levels[L][0] - levels[L + 1][0]
            members = [frags[i] for i in levels[L][1]] + [s + f":{edge:.6f}"]
            s = "(" + ",".join(members) + ")"
        s = s + f":{budget_ - levels[0][0]:.6f}"             # outermost node's edge to the parent
        stack.pop()
        if stack:
            par = stack[-1]
            par["frags"][par["i"] - 1] = s
        else:
            ret = s
    return ret


def emit_dated(cells_ids, isbb, parent, held_set, query_idx, backbone_nwk, pendant_days):
    """Grow the dated backbone with time-scaled query attachments and emit Newick.

    Returns (newick_string, placed, max_children, clamped_set). Backbone internal edges
    keep their verbatim (day) lengths; each backbone leaf hosting queries is replaced by
    its `_expand` caterpillar spending that leaf's dated terminal edge as the budget.
    clamped_set is the set of cell indices whose pendant exceeded the available edge.
    """
    N = len(cells_ids)
    attach_kids, placed = _build_attach_kids(N, parent, held_set, query_idx)
    max_children = max((len(k) for k in attach_kids), default=0)
    clamped = set()

    bp, bname, bkids, bblen = parse_structure(open(backbone_nwk).read().strip())
    name2g = {cells_ids[i]: i for i in range(N) if isbb[i]}

    out = []
    stack = [("open", 0)]
    while stack:
        kind, v = stack.pop()
        if kind == "close":
            out.append(")")
            if bblen[v]:
                out.append(":" + bblen[v])
            continue
        if kind == "sep":
            out.append(",")
            continue
        # open node v
        if not bkids[v]:                                     # backbone leaf
            g = name2g.get(bname[v])
            if g is not None and attach_kids[g]:
                budget = float(bblen[v]) if bblen[v] else 0.0
                out.append(_expand(g, budget, attach_kids, pendant_days,
                                   cells_ids, clamped))
            else:
                out.append(bname[v])
                if bblen[v]:
                    out.append(":" + bblen[v])
            continue
        out.append("(")                                     # backbone internal
        ch = bkids[v]
        stack.append(("close", v))
        for i in range(len(ch) - 1, -1, -1):
            stack.append(("open", ch[i]))
            if i > 0:
                stack.append(("sep", None))
    out.append(";")
    return "".join(out), placed, max_children, clamped
