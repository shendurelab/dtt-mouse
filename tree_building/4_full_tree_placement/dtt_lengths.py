#!/usr/bin/env python3
"""Pendant length for a placed cell, in the backbone's OWN distance units.

We attach query C as sibling of the cell T it matched, under their inferred common
ancestor P (the shared EDITED PREFIX per tape, the SciPhy sequential-editing MRCA).
The pendant we assign C is the distance d(C, P) measured with the SAME metric that
built the NJ distance matrix / backbone — 1_build_nj_backbone/dtt_distance.R
(`dtt_distance` / its bit-identical kernel dtt_distance_scale.cpp):

    per shared tape:  tape_distance = depth1 + depth2 - 2 * shared_edit_prefix
    cell distance   = MEAN of that over tapes recovered in BOTH cells.

For P = shared-prefix ancestor of (C, T), depth_P (per shared tape) = shared_edit_prefix,
so per shared tape  d(C,P) = depth_C - prefix  and  d(T,P) = depth_T - prefix, and
    d(C,P) + d(T,P) = d(C,T)      (additive through P; verified to 1e-15).
Because it is the matrix metric, d(C,P) is in the backbone's edge units EXACTLY
(verified: on backbone cherries, matrix d(C,T) == tree edge_a+edge_b, ratio 1.000).
That is what makes the downstream  pendant_days = d(C,P) / clock_rate  valid with the
per-side LSD2 rate (subst-per-day on the SAME scale) and NO unit fudge factor.

KEY DIFFERENCE from a naive per-tape pendant: the private-edit numerator is summed
ONLY over tapes recovered in BOTH C and T (the shared tapes the matrix averages over),
not over all of C's tapes. Averaging C's own-tape edits over the shared-tape count
inflated the pendant ~3x and put it on the wrong scale for the clock rate.

Model per integration ("tape", a length-N_SITES token vector; encoding from genotypes.py):
    code == -1  tape/site missing ("?"/"None")   -- whole integration absent
    code ==  0  unedited ("ETY")
    code >=  1  a specific edit token

Returned by compute_pendants(): (pendant, residual) float arrays indexed by global
cell index --  pendant[qi] = d(C,P) (C's new edge, DTT units), residual[qi] = d(T,P)
(the anchor-side residual, for additivity checks; not written into the DTT tree).
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genotypes import load_and_encode

# Supp. Table 1's blastomere label -> our SIDE (B1/B2) column prefix.
_FOUNDER_COL = {"B1": "A", "B2": "B"}


def load_founder_len(founder_csv, side, integ_cols):
    """Per-tape founder-prefix length to exclude from the shared-edit-prefix metric.

    Returns an (G,) int array aligned to `integ_cols` (same g-axis as `code`):
    founder_len[g] = number of leading sites of that integration fixed (>=90%
    consensus) within `side`'s blastomere at its founding -- i.e. shared by every
    cell of the side by construction, so not evidence of closer-than-average
    relatedness between a query and its anchor. Integrations absent from the table
    (or "(unedited)") get 0.
    """
    col = f"{_FOUNDER_COL[side]}_n_edits"
    n_edits = {}
    with open(founder_csv) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            n_edits[row["integration"]] = int(row[col])
    return np.fromiter((n_edits.get(g, 0) for g in integ_cols), dtype=np.int64,
                       count=len(integ_cols))


def _aligned_code(consensus_tsv, n_sites, allcells, backbone_tips, cellids_order):
    """Load + encode genotypes and reorder rows to match `cellids_order` exactly.

    The best-match checkpoints fix a global cell order (cellids_<SIDE>.npy); we index
    the freshly-encoded genotypes into that same order so `code[i]` is the genotype of
    the cell at global index i.
    """
    cells = load_and_encode(consensus_tsv, backbone_tips, n_sites, allcells)
    pos = {cid: i for i, cid in enumerate(cells.ids)}
    missing = [c for c in cellids_order if c not in pos]
    if missing:
        sys.exit(f"[dtt] {len(missing):,} checkpoint cells missing from consensus "
                 f"(e.g. {missing[:3]}); genotype source does not match the checkpoint")
    order = np.fromiter((pos[c] for c in cellids_order), dtype=np.int64,
                        count=len(cellids_order))
    return cells.code[order], cells.integ_cols


def pendant_components(code, parent, founder_len=None):
    """Return per-cell pendant components given aligned genotype `code` and parents.

    Arrays indexed by global cell index (0 for backbone/held cells):
      pendant      = d(C, P) = new_edits / n_shared  (C's new edge, matrix DTT units)
      residual     = d(T, P) = anchor-side residual  (for the additivity check)
      new_edits    = INTEGER count of edits C acquired on its own branch below P,
                     summed over tapes shared with T   ("new edits per query")
      shared_edits = INTEGER count of edits C and T agree on BEYOND each tape's
                     blastomere-founder prefix (the shared edited prefix, minus the
                     founder_len[g] leading sites every cell of the side already
                     shares by construction), summed over shared tapes
                     ("overlapping edits over shared tapes")
      n_shared     = INTEGER number of tapes recovered in BOTH C and T

    `founder_len`, if given, is a (G,) int array (see dtt_lengths.load_founder_len):
    the per-tape blastomere-founder prefix length to exclude from `shared_edits` only
    -- it does NOT affect pendant/residual/new_edits, which are the matrix-metric
    dating quantities and already net the shared prefix out of both branches symmetrically.
    Fully vectorised over all placed queries.
    """
    N = code.shape[0]
    pendant = np.zeros(N, dtype=np.float64)
    residual = np.zeros(N, dtype=np.float64)
    new_edits = np.zeros(N, dtype=np.int64)
    shared_edits = np.zeros(N, dtype=np.int64)
    n_shared = np.zeros(N, dtype=np.int64)

    qi = np.where(parent >= 0)[0]                      # every placed cell (has a parent)
    if len(qi) == 0:
        return pendant, residual, new_edits, shared_edits, n_shared

    cQ = code[qi]                                      # (nq, G, S)
    cP = code[parent[qi]]                              # (nq, G, S)

    # tapes recovered (present, not all-missing) in each cell; a tape is all -1 iff
    # absent, so site 0 decides. SHARED tapes = recovered in both -> what the matrix
    # metric averages over.
    q_obs = cQ[:, :, 0] != -1                          # (nq, G)
    p_obs = cP[:, :, 0] != -1
    shared = q_obs & p_obs                             # (nq, G)
    shared3 = shared[:, :, None]                       # (nq, G, 1)

    # shared EDITED PREFIX per tape: True from site 0 until the first non-matching edit
    # (a differing edit or an unedited/missing site). Naturally 0 on non-shared tapes.
    match_edit = (cQ >= 1) & (cQ == cP)
    prefix = np.cumprod(match_edit, axis=2).astype(bool)
    np_prefix = (prefix & shared3).sum(axis=(1, 2))    # overlapping (shared) edits

    # shared_edits (panel-facing metric only): drop each tape's founder-prefix length
    # from its match-run before summing, so a founder edit -- shared by every cell of
    # the blastomere by construction -- isn't reported as evidence of a closer anchor.
    if founder_len is not None:
        prefix_len = prefix.sum(axis=2)                        # (nq, G) per-tape run length
        informative_len = np.maximum(prefix_len - founder_len[None, :], 0)
        shared_edits_reported = (informative_len * shared).sum(axis=1)
    else:
        shared_edits_reported = np_prefix

    # private edits summed ONLY over shared tapes (this is the matrix-metric numerator)
    depth_Q_shared = ((cQ >= 1) & shared3).sum(axis=(1, 2))
    depth_T_shared = ((cP >= 1) & shared3).sum(axis=(1, 2))
    branch_Q = np.maximum(depth_Q_shared - np_prefix, 0)   # new edits C acquired below P
    branch_T = np.maximum(depth_T_shared - np_prefix, 0)

    nshared = shared.sum(axis=1)
    ns = np.maximum(nshared, 1)                        # denominator, floored at 1

    pendant[qi] = branch_Q / ns                        # d(C, P), backbone DTT units
    residual[qi] = branch_T / ns                       # d(T, P)
    new_edits[qi] = branch_Q
    shared_edits[qi] = shared_edits_reported
    n_shared[qi] = nshared
    print(f"[dtt] d(C,P) pendants for {len(qi):,} placed cells: "
          f"mean {pendant[qi].mean():.3f}  median {np.median(pendant[qi]):.3f}  "
          f"max {pendant[qi].max():.3f}  (matrix DTT units); "
          f"new_edits median {int(np.median(branch_Q))}  "
          f"shared_edits median {int(np.median(shared_edits_reported))}"
          + ("" if founder_len is None else " (founder-prefix excluded)"),
          file=sys.stderr)
    return pendant, residual, new_edits, shared_edits, n_shared


def pendants_from_code(code, parent):
    """(pendant, residual) only -- thin wrapper over pendant_components(), used by
    the test suite (test_dtt_pendant.py, test_pendant_matches_R.py)."""
    pendant, residual, _ne, _se, _ns = pendant_components(code, parent)
    return pendant, residual


def compute_pendants(consensus_tsv, n_sites, allcells, backbone_tips,
                     cellids_order, parent):
    """Convenience wrapper: load + align genotypes, then pendants_from_code()."""
    code, _integ_cols = _aligned_code(consensus_tsv, n_sites, allcells, backbone_tips,
                                      cellids_order)
    return pendants_from_code(code, parent)
