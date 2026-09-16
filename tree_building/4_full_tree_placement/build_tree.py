#!/usr/bin/env python3
"""Checkpoint loading and parent resolution, shared by finalize_dated.py.

The rule for growing the tree is deliberately simple (two fancier variants --
bootstrap re-routing and edge/LCA attachment -- were both tried and both
underperformed this one):

    parent(query) = its single DTT-closest cell (backbone cell or another query)

Consequences and the bookkeeping they require:
  * A query whose closest cell is another query chains query -> query -> ... ->
    backbone. These chains nest naturally and are what preserve genuine
    query-query relatives.
  * Backbone cells are the roots of these chains (they are fixed, parent = none).
  * A query whose closest-cell search found nothing on the backbone side but has a
    general closest cell still attaches; a query with a backbone match but no
    better overall attaches to the backbone.
  * Mutual-nearest-cell cycles among query-only cells (A's closest is B, B's is A,
    neither on the backbone) would never reach a root. We break each cycle by
    redirecting its best-backbone-scoring member to its backbone match, so the
    whole cycle grounds out; if that member has no backbone match anywhere either,
    the whole cycle is held.

Not every query ends up in the tree: resolve_parents() returns a `held_set` of
query indices that could not be grounded to a backbone cell -- a query with NO
indexed token at all (no shared edit with anything), or one whose chain runs into
an ungroundable cycle (above). finalize_dated.py drops held cells from the output
tree and does not write them a place_qc.tsv row; the run's stderr summary
(`held N,NNN`) is the only record of how many were dropped and is worth checking.
"""
import glob
import sys

import numpy as np


def load_bestmatch(outdir, side, N):
    """Reassemble per-cell best-match arrays from all bestmatch_*.npz chunks.

    Returns (gmatch, gscore, bbmatch, bbscore), each length N indexed by global cell
    index (-1 / 0 where a cell was never a query, i.e. backbone cells).
    """
    gmatch = np.full(N, -1, np.int64)
    gscore = np.zeros(N, np.float32)
    bbmatch = np.full(N, -1, np.int64)
    bbscore = np.zeros(N, np.float32)
    files = sorted(glob.glob(f"{outdir}/bestmatch_*.npz"))
    if not files:
        sys.exit(f"[build_tree] no bestmatch_*.npz in {outdir}; run bestmatch.py first")
    for f in files:
        z = np.load(f)
        q = z["q"]
        gmatch[q] = z["m"]
        gscore[q] = z["s"]
        bbmatch[q] = z["mb"]
        bbscore[q] = z["sb"]
    return gmatch, gscore, bbmatch, bbscore


def resolve_parents(query_idx, isbb, gmatch, bbmatch, bbscore):
    """Assign each query a parent and make sure every query reaches a backbone cell.

    Returns (parent, held_set):
      parent[i]  parent cell index (-1 for backbone roots / held cells)
      held_set   set of query indices that could not be grounded and are dropped
    """
    N = len(isbb)
    parent = np.full(N, -1, np.int64)
    held = []
    # initial parent = closest cell, falling back to the backbone match, else hold.
    for qi in query_idx:
        if gmatch[qi] < 0:
            if bbmatch[qi] >= 0:
                parent[qi] = bbmatch[qi]
            else:
                held.append(qi)
        else:
            parent[qi] = gmatch[qi]

    # Walk each query's parent chain; guarantee it terminates at a backbone cell.
    reaches = np.zeros(N, np.int8)      # 1 = this cell's chain reaches the backbone
    reaches[isbb] = 1
    held_set = set(held)
    for qi in query_idx:
        if qi in held_set or reaches[qi]:
            continue
        path = []
        cur = qi
        while True:
            if reaches[cur]:                        # joined a known-good chain
                for n in path:
                    reaches[n] = 1
                break
            if cur in held_set:                     # chain runs into a held cell
                for n in path:
                    parent[n] = parent[cur] if parent[cur] >= 0 else bbmatch[n]
                break
            if cur in path:                         # cycle detected
                cyc = path[path.index(cur):]
                best = max(cyc, key=lambda n: bbscore[n])   # ground its best-bb member
                parent[best] = bbmatch[best] if bbmatch[best] >= 0 else best
                if bbmatch[best] < 0:               # no backbone anchor anywhere -> hold all
                    for n in cyc:
                        held_set.add(n)
                    break
                for n in path:                      # restart resolution from qi
                    reaches[n] = 0
                cur = qi
                path = []
                continue
            path.append(cur)
            cur = parent[cur]
            if cur < 0:                             # ran off the end with no anchor
                for n in path:
                    reaches[n] = 0
                held_set.add(path[0])
                break
    return parent, held_set
