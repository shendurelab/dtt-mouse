#!/usr/bin/env python3
"""Load a consensus callset and encode each cell's genotype into a dense integer
array of edit-token ids.

The genotype of a cell is its set of TAPE edits across all integrations. We model
each cell as a fixed-shape array `code[cell, integration, site]` where each entry is:

    -1  the integration is absent from this cell (no chain reported)
     0  the integration is present but this site is uninformative (unedited /
        dropout / ambiguous: "", "NA", "U", "-")
   >=1  a small integer id for the specific edit allele seen at (integration, site)

Allele ids are assigned per (integration, site) from a vocabulary built over the
whole callset, so id 3 at (integration 2, site 0) is the same physical edit for
every cell. Only code>=1 entries carry lineage signal; the placement search works
entirely from those.

Returned `Cells` bundles everything the downstream steps need and nothing else.
"""
import csv
import gzip
import sys
import time
from dataclasses import dataclass

import numpy as np

from config import QC_COLS


@dataclass
class Cells:
    ids: np.ndarray          # (N,) object array of cell_id strings, the canonical order
    code: np.ndarray         # (N, G, N_SITES) int16, encoded genotypes (see module docstring)
    isbb: np.ndarray         # (N,) bool: True if this cell is a fixed backbone anchor
    nloci: np.ndarray        # (N,) int16: reported informative loci (for QC/reporting)
    G: int                   # number of integration columns
    n_sites: int             # sites per chain (== config.N_SITES)
    integ_cols: list         # (G,) integration barcode names, in `code`'s g-axis order


def _informative(site_val):
    """True if a per-site token carries a real edit (not a dropout/ambiguity marker)."""
    return site_val not in ("", "NA", "U", "-")


def load_and_encode(consensus_tsv, backbone_tips, n_sites, allcells):
    """Read the gzipped consensus TSV and encode it.

    Parameters
    ----------
    consensus_tsv : path to the gzipped, header-bearing consensus table.
    backbone_tips : set[str] of cell_ids that are fixed anchors (from the backbone).
    n_sites       : edit sites per integration chain (config.N_SITES).
    allcells      : if True keep every cell with n_loci>=1; else keep only pass_qc==1.

    Returns a `Cells`. Cell order is the file order of the kept rows and is the
    canonical global index used everywhere downstream (best-match npz, tree build).
    """
    t0 = time.time()

    # --- pass 1: read the rows we keep, holding the raw chain strings ---
    rows = []                      # (cell_id, [chain per integration], n_loci)
    with gzip.open(consensus_tsv, "rt") as fh:
        r = csv.reader(fh, delimiter="\t")
        hdr = next(r)
        integ_cols = [c for c in hdr[1:] if c not in QC_COLS]
        colpos = {c: hdr.index(c) for c in integ_cols}
        i_pass = hdr.index("pass_qc")
        i_nloci = hdr.index("n_loci")
        for row in r:
            nl = int(row[i_nloci])
            keep = (nl >= 1) if allcells else (int(row[i_pass]) == 1)
            if keep:
                rows.append((row[0], [row[colpos[c]] for c in integ_cols], nl))

    G = len(integ_cols)
    N = len(rows)
    ids = np.array([r[0] for r in rows], dtype=object)
    nloci = np.array([r[2] for r in rows], dtype=np.int16)
    isbb = np.array([r[0] in backbone_tips for r in rows], dtype=bool)
    print(f"[genotypes] kept {N:,} cells (backbone {int(isbb.sum()):,}); "
          f"{G} integrations  {time.time()-t0:.1f}s", file=sys.stderr)

    # --- pass 2a: build the per-(integration, site) allele vocabulary ---
    # vocab[(g, s)] = {allele_string: small_int_id starting at 1}
    vocab = {}
    for _cid, chains, _nl in rows:
        for g, chain in enumerate(chains):
            if chain and chain != "NA":
                for s, allele in enumerate(chain.split("|")):
                    if s >= n_sites:
                        break
                    if _informative(allele):
                        d = vocab.setdefault((g, s), {})
                        if allele not in d:
                            d[allele] = len(d) + 1

    # --- pass 2b: encode into the dense int array ---
    code = np.full((N, G, n_sites), -1, dtype=np.int16)   # default: integration absent
    for i, (_cid, chains, _nl) in enumerate(rows):
        for g, chain in enumerate(chains):
            if chain and chain != "NA":
                code[i, g, :] = 0                          # present -> sites default uninformative
                for s, allele in enumerate(chain.split("|")):
                    if s >= n_sites:
                        break
                    if _informative(allele):
                        code[i, g, s] = vocab[(g, s)][allele]
    print(f"[genotypes] encoded {time.time()-t0:.1f}s", file=sys.stderr)

    return Cells(ids=ids, code=code, isbb=isbb, nloci=nloci, G=G, n_sites=n_sites,
                integ_cols=integ_cols)
