#!/usr/bin/env python3
"""DTT-ANCESTOR placement scoring: pick each query's sibling tip by MINIMISING the
project's own cell-cell DTT distance d(C,T) -- the same distance
(1_build_nj_backbone/dtt_distance.R) that builds the NJ matrix, the backbone, and the
grafted pendant now also chooses the attachment. (There is no "k" here -- each cell
attaches to its single closest tip -- so we call it the DTT-ANCESTOR, not a kNN.)

Two-stage search for tractability (an exact all-cell argmin is far too costly):
(1) SHORTLIST -- using the inverted index, count how many indexed edit tokens each
cell shares with the query (a bincount, unweighted), and keep the TOPK cells with
the most shared tokens; the DTT-closest cell shares the query's whole edited prefix,
so it shares the most tokens and is in the shortlist. (2) RANK -- compute the EXACT
d(C,T) over the shortlist and take the argmin. Stage 1 is a candidate filter only;
the ranking is pure DTT distance.

Tie-break among cells at the minimum distance (DTT distance is coarse, so ties are
common): most shared tapes (most evidence) -> highest n_loci (best-quality anchor) ->
lowest global index. Returns (m, s, mb, sb): m/s the best-overall match+score (s = -d,
so "higher is better"), mb/sb the best-backbone match+score.
"""
import sys
import time

import numpy as np

from similarity import token_id


def _dtt_over(cand, qcode, qdepth, qrecov, depth, recov):
    """Exact d(C,T) (mean over shared tapes) for each cell in `cand`. Returns (d, nsh)."""
    cc = _CODE[cand]                                            # (nc, G, S)
    shared = qrecov[None, :] & recov[cand]                      # (nc, G) tapes in both
    pref = np.cumprod((cc >= 1) & (cc == qcode[None]), axis=2).sum(axis=2)   # (nc, G) shared prefix
    per = qdepth[None, :] + depth[cand] - 2 * pref             # (nc, G) per-tape distance
    per = np.where(shared, per, 0)
    nsh = shared.sum(axis=1)
    d = np.where(nsh > 0, per.sum(axis=1) / np.maximum(nsh, 1), np.inf)
    return d, nsh


_CODE = None   # module-global handle to the genotype array (set in best_matches_dtt)


def _rank_min(d, nsh, nloci_c, cand):
    """Index into `cand` of the min-distance cell, tie-broken by
    (shared tapes desc, n_loci desc, global index asc). None if all inf."""
    dmin = d.min()
    if not np.isfinite(dmin):
        return None, dmin
    tie = np.where(d == dmin)[0]
    if len(tie) == 1:
        return int(tie[0]), dmin
    order = np.lexsort((cand[tie], -nloci_c[tie], -nsh[tie]))
    return int(tie[order[0]]), dmin


def best_matches_dtt(code, G, n_sites, stride, postings, query_idx,
                     isbb, nloci, topk=500):
    """For each query, its DTT-closest cell overall and its DTT-closest backbone tip.

    Returns arrays aligned to query_idx: m/s best-overall (s = -d), mb/sb best-backbone.
    """
    global _CODE
    _CODE = code
    t0 = time.time()
    N = code.shape[0]
    nq = len(query_idx)
    m = np.full(nq, -1, np.int64); s = np.zeros(nq, np.float32)
    mb = np.full(nq, -1, np.int64); sb = np.zeros(nq, np.float32)

    depth = (code >= 1).sum(axis=2)           # (N, G) edited sites per tape
    recov = code[:, :, 0] != -1               # (N, G) tape recovered (present)

    for r, qi in enumerate(query_idx):
        qcode = code[qi]
        # --- stage 1: shortlist by shared indexed-token count (bincount) ---
        idx_lists = []
        for g in range(G):
            for si in range(n_sites):
                a = int(qcode[g, si])
                if a < 1:
                    continue
                arr = postings.get(token_id(g, si, a, n_sites, stride))
                if arr is not None:
                    idx_lists.append(arr)
        if not idx_lists:
            continue
        cnt = np.bincount(np.concatenate(idx_lists), minlength=N)
        cnt[qi] = 0
        cand = np.nonzero(cnt)[0]
        if len(cand) == 0:
            continue
        if len(cand) > topk:
            keep = np.argpartition(cnt[cand], -topk)[-topk:]     # top-K by shared-token count
            cand = cand[keep]

        # --- stage 2: exact d(C,T) over the shortlist ---
        d, nsh = _dtt_over(cand, qcode, depth[qi], recov[qi], depth, recov)
        nloci_c = nloci[cand]
        best_local, dmin = _rank_min(d, nsh, nloci_c, cand)
        if best_local is None:
            continue
        m[r] = int(cand[best_local]); s[r] = -float(dmin)

        bbmask = isbb[cand]
        if bbmask.any():
            bl, dbmin = _rank_min(np.where(bbmask, d, np.inf), nsh, nloci_c, cand)
            if bl is not None:
                mb[r] = int(cand[bl]); sb[r] = -float(dbmin)

        if r % 20000 == 0:
            print(f"[dtt-match] scored {r:,}/{nq:,}  {time.time()-t0:.1f}s", file=sys.stderr)

    print(f"[dtt-match] scored {nq:,} queries (topk={topk})  {time.time()-t0:.1f}s",
          file=sys.stderr)
    return m, s, mb, sb
