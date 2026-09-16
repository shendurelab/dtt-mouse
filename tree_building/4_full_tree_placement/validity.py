#!/usr/bin/env python3
"""Irreversibility (is_valid_ancestor) DIAGNOSTIC for DTT-ancestor attachments.

Edits are irreversible, so a placed cell C cannot descend from (or share an ancestor
consistent with) a matched cell T that carries an edit C resolves differently. The
DTT-ancestor search has no such check built in (it minimises distance, not validity),
so it can occasionally propose an attachment that violates this. We do NOT reject
those attachments here (the topology stays as matched); we FLAG them, so the QC table
records how many attachments violate irreversibility.

invalid  <=>  any (tape, site) with  T edited (code_T >= 1)
                               AND  C resolved (code_C != -1, i.e. tape observed)
                               AND  C differs (code_C != code_T)
A tape missing in C (code_C == -1 across it) never triggers -- "?" is a wildcard.
Encoding is genotypes.py's: -1 missing, 0 unedited, >=1 a specific edit token.
"""
import numpy as np


def valid_flags(code, parent):
    """Return a bool array indexed by global cell index.

    valid[i] is True for placed cell i (parent >= 0) iff its attachment to code[parent]
    respects irreversibility; True by default for non-placed cells (they have no edge).
    """
    N = code.shape[0]
    valid = np.ones(N, dtype=bool)
    qi = np.where(parent >= 0)[0]
    if len(qi) == 0:
        return valid
    cQ = code[qi]                                    # (nq, G, S)
    cP = code[parent[qi]]
    invalid_site = (cP >= 1) & (cQ != -1) & (cQ != cP)
    valid[qi] = ~invalid_site.any(axis=(1, 2))
    return valid
