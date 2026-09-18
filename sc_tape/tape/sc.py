"""Single-cell TAPE front-end.

Differs from bulk only at the ends of the shared parser:
  - reads carry a cell id + UMI in the header (parse_header);
  - reads begin at the integration barcode rather than after a GGCATG flank
    (extract_bc_sc);
  - genotypes are resolved per (cell, integration, UMI), PCR duplicates are
    collapsed to molecules by UMI, and a single consensus genotype is called per
    (cell, integration) via the shared edit-chain model.

Everything in between -- classify_junction, walk_robust, the vocabulary, the
de-noising funnel -- is the same code the bulk pipeline uses.
"""
from __future__ import annotations
from collections import Counter, defaultdict

from config import (BC_RIGHT, BC_RIGHT_SPACER, BC_LEN, N_SITES, UMI_MERGE_DIST,
                    MIN_READS_PER_UMI, SC_MIN_MOLECULES, SC_DOMINANCE_MIN)
from tape.barcode import hamming_le
from tape.editchain import build_trie, dominant_lineage, edit_chain


def parse_header(header: bytes):
    """@<cell_id>,PCR-<well>,<UMI>,<sample>  ->  (cell_id, umi, sample) or None."""
    h = header[1:] if header[:1] == b"@" else header
    parts = h.split(b",")
    if len(parts) < 4:
        return None
    return parts[0], parts[2], parts[3]


def extract_bc_sc(seq: bytes):
    """Extract the integration barcode from a single-cell read.

    SC reads begin at the barcode: [12-bp BC][spacer G]CTCTGG[array...]. Anchor on
    CTCTGG, take the barcode immediately 5' of it (dropping the constant spacer).
    Returns (bc, cstart) where cstart is the array start (past CTCTGG), or (None,None).
    """
    spacer = len(BC_RIGHT_SPACER)
    need = BC_LEN + spacer
    p = seq.find(BC_RIGHT)
    if p < need:                                   # not enough room for BC + spacer
        return None, None
    bc = seq[p - need:p - spacer]                  # 12-bp barcode (spacer stripped)
    if len(bc) != BC_LEN:
        return None, None
    return bc, p + len(BC_RIGHT)


def collapse_umis(umi_patterns: dict, min_reads_per_umi: int = MIN_READS_PER_UMI) -> Counter:
    """Collapse PCR duplicates to molecules within one (cell, integration).

    umi_patterns: {umi: {pattern: reads}}. UMIs within UMI_MERGE_DIST of a more
    abundant UMI are merged into it (directional). Each surviving UMI is one
    molecule whose genotype is its dominant pattern. Returns Counter{pattern: n_molecules}.

    A merged UMI with fewer than `min_reads_per_umi` reads is dropped (1-read UMIs
    are enriched for PCR/ambient error). Pass min_reads_per_umi=1 to keep every UMI
    (e.g. for raw per-UMI read statistics).
    """
    # total reads per umi, most abundant first
    order = sorted(umi_patterns, key=lambda u: -sum(umi_patterns[u].values()))
    merged = {}                                    # umi -> combined {pattern: reads}
    reps = []                                      # surviving representative umis
    for u in order:
        hit = None
        if UMI_MERGE_DIST > 0:
            for r in reps:
                if len(r) == len(u) and hamming_le(u, r, UMI_MERGE_DIST) <= UMI_MERGE_DIST:
                    hit = r
                    break
        if hit is None:
            reps.append(u)
            merged[u] = dict(umi_patterns[u])
        else:
            for p, c in umi_patterns[u].items():
                merged[hit][p] = merged[hit].get(p, 0) + c
    mol_genotypes = Counter()
    for u in reps:
        if sum(merged[u].values()) < min_reads_per_umi:          # (#1) reads-per-UMI floor
            continue
        best = max(merged[u].items(), key=lambda kv: kv[1])[0]   # molecule's dominant pattern
        mol_genotypes[best] += 1
    return mol_genotypes


def consensus_locus(mol_genotypes: Counter,
                    min_molecules: int = SC_MIN_MOLECULES,
                    dom_min: float = SC_DOMINANCE_MIN):
    """Call the consensus genotype for one (cell, integration) from its molecule
    genotypes, using the edit-chain model.

    The dominant lineage is the deepest dominant tip; its support = its own UMIs plus
    all prefix (3' dropout) UMIs on its path, and `dominance` = that support / all
    molecules (only branching lineages -- a different edit at the same site -- count
    against it). The locus is CALLED only when the dominant lineage is both:
      (#2) supported by >= `min_molecules` molecules, and
      (#3) dominant with ratio >= `dom_min`;
    otherwise it is left MISSING (return None).

    Returns (chain, branched, dominance, n_support) where n_support is the molecules
    on the called lineage.
    """
    if sum(mol_genotypes.values()) == 0:
        return None
    root = build_trie(mol_genotypes)
    chain, n_support, dominance, branched = dominant_lineage(root)
    if n_support < min_molecules:       # (#2) under-supported
        return None
    if dominance < dom_min:             # (#3) not clearly dominant (doublet / mixture)
        return None
    return chain, branched, dominance, int(round(n_support))


def chain_to_str(chain) -> str:
    """Render an edit chain as a genotype string, e.g. GCC-AAG-CCC (or '-' if none)."""
    return "-".join(chain) if chain else "-"
