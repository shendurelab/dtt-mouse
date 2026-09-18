"""De-noising: the single canonical implementation of the bulk filtering funnel.

The original workspace had this logic copy-pasted (and subtly diverging) across
four scripts. It lives here once; both the bulk stages and the single-cell
pipeline import it. Order matters and mirrors the validated pipeline:

    raw patterns
      -> collapse      : fold 1-symbol-error neighbors into abundant parents
      -> dechimera     : drop UCHIME-style prefix-A/suffix-B PCR chimeras
      -> >= MIN_READS  : drop patterns below the read-support floor
      -> cross-barcode : drop deep patterns whose max-read "home" is another barcode
                         (TAPE-BC-swap recombinants)

Every function takes/returns {pattern_tuple: count} dicts. Patterns are 6-tuples
of insertion strings / 'U'. `count` is reads (bulk) or molecules (single cell).
"""
from __future__ import annotations
from collections import defaultdict

from config import N_SITES, COLLAPSE_RATIO, MIN_READS, CROSS_BC_DEPTH


def neighbors(pattern):
    """Yield all patterns one substitution away at an EDITED site (never mutate a
    'U' -- an unedited site carries no sequence to misread)."""
    for i in range(N_SITES):
        s = pattern[i]
        if s == "U":
            continue
        for j in range(len(s)):
            for ch in "ACGT":
                if ch != s[j]:
                    yield pattern[:i] + (s[:j] + ch + s[j + 1:],) + pattern[i + 1:]


def collapse(counts: dict, ratio: float = COLLAPSE_RATIO) -> dict:
    """Error-collapse: fold each pattern into its most abundant 1-symbol neighbor
    when that neighbor is >= `ratio`x its count (and strictly larger). Processes
    rarest-first so error chains collapse toward their true parent."""
    cur = dict(counts)
    for pat, _ in sorted(counts.items(), key=lambda x: x[1]):
        c = cur.get(pat, 0)
        if c == 0:
            continue
        best = None
        thresh = c * ratio
        for q in neighbors(pat):
            cn = cur.get(q, 0)
            if cn >= thresh and cn > c:
                thresh = cn
                best = q
        if best is not None:
            cur[best] = cur.get(best, 0) + c
            del cur[pat]
    return cur


def dechimera(counts: dict) -> dict:
    """UCHIME-style chimera removal. Walking most-abundant-first, drop a
    FULLY-EDITED pattern that can be explained as prefix(A)+suffix(B) of two
    already-seen more-abundant patterns. Partially-edited patterns (containing a
    'U') are exempt: their unedited tails are legitimately shared by many real
    lineages, so they must not be flagged as chimeric."""
    items = sorted(counts.items(), key=lambda x: -x[1])
    pre = [set() for _ in range(N_SITES)]
    suf = [set() for _ in range(N_SITES)]
    keep = {}
    for pat, c in items:
        chim = ("U" not in pat) and any(
            pat[:k] in pre[k] and pat[k:] in suf[k] for k in range(1, N_SITES))
        if not chim:
            keep[pat] = c
        for k in range(1, N_SITES):
            pre[k].add(pat[:k])
            suf[k].add(pat[k:])
    return keep


def depth(pattern) -> int:
    """Number of edited sites in a pattern."""
    return sum(1 for s in pattern if s != "U")


def denoise_barcode(counts: dict, ratio: float = COLLAPSE_RATIO,
                    min_reads: int = MIN_READS) -> dict:
    """collapse -> dechimera -> drop < min_reads, for a single barcode."""
    clean = dechimera(collapse(counts, ratio))
    return {p: c for p, c in clean.items() if c >= min_reads}


def denoise_all(per_bc: dict, ratio: float = COLLAPSE_RATIO,
                min_reads: int = MIN_READS, cross_depth: int = CROSS_BC_DEPTH):
    """Run the full funnel across all barcodes and apply the cross-barcode filter.

    Args:
      per_bc: {barcode: {pattern: count}} of raw patterns.
    Returns:
      (clean, funnels) where
        clean[bc]  = {pattern: count} of retained patterns (cross-BC removed)
        funnels[bc]= dict of stage counts for a QC funnel table
    """
    # per-barcode collapse + dechimera + read floor
    stage = {}
    funnels = {}
    for bc, counts in per_bc.items():
        collapsed = collapse(counts, ratio)
        dechim = dechimera(collapsed)
        ge = {p: c for p, c in dechim.items() if c >= min_reads}
        stage[bc] = ge
        funnels[bc] = {"raw": len(counts), "collapsed": len(collapsed),
                       "chimera_filtered": len(dechim), "ge_min_reads": len(ge)}

    # cross-barcode: a pattern's "home" is the barcode where it has the most reads
    glob = defaultdict(dict)
    for bc, counts in stage.items():
        for p, c in counts.items():
            glob[p][bc] = c
    home = {p: max(d, key=d.get) for p, d in glob.items()}

    def is_cross(p, bc):
        return depth(p) >= cross_depth and home[p] != bc

    clean = {}
    for bc, counts in stage.items():
        kept = {p: c for p, c in counts.items() if not is_cross(p, bc)}
        clean[bc] = kept
        funnels[bc]["cross_bc"] = len(counts) - len(kept)
        funnels[bc]["clean"] = len(kept)
    return clean, funnels, home
