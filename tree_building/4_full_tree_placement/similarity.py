#!/usr/bin/env python3
"""Build the rarity-capped inverted index that dtt_match.py's shortlist stage
searches: `postings[token] -> array of cell indices carrying that token`.

We never materialise an N x N matrix. To shortlist candidates for a query we walk
only its own tokens and gather the cells sharing each one via `postings`.

Rarity cap
----------
Tokens carried by more than `cap` cells are dropped from the index entirely — they
are homoplastic / early-shared edits that link unrelated lineages and would bloat
every query's candidate shortlist without adding lineage signal.
"""
import sys
import time
from collections import defaultdict

import numpy as np


def token_id(g, s, allele_id, n_sites, stride):
    """Pack (integration g, site s, allele_id) into one int token id."""
    return (g * n_sites + s) * stride + allele_id


def token_freq(code, G, n_sites, stride):
    """Token frequencies over all cells: dict token_id -> count.

    Only informative alleles (code>=1) are counted.
    """
    freq = defaultdict(int)
    for g in range(G):
        for s in range(n_sites):
            col = code[:, g, s]
            vals, counts = np.unique(col[col >= 1], return_counts=True)
            for a, n in zip(vals.tolist(), counts.tolist()):
                freq[token_id(g, s, a, n_sites, stride)] += n
    return freq


def build_index(code, G, n_sites, stride, freq, cap_frac):
    """Build the inverted index over anchor cells.

    Returns postings: token -> np.int64 array of cell indices carrying that token
    (only tokens with freq <= cap_frac*N are indexed).
    """
    t0 = time.time()
    N = code.shape[0]
    cap = int(cap_frac * N)
    postings = defaultdict(list)
    for i in range(N):
        for g in range(G):
            for s in range(n_sites):
                a = int(code[i, g, s])
                if a < 1:
                    continue
                t = token_id(g, s, a, n_sites, stride)
                if freq[t] <= cap:
                    postings[t].append(i)
    postings = {t: np.array(v, dtype=np.int64) for t, v in postings.items()}
    print(f"[similarity] index built over {N:,} anchors, "
          f"{len(postings):,} tokens  {time.time()-t0:.1f}s", file=sys.stderr)
    return postings
