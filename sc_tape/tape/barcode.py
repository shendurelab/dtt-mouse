"""TAPE integration-barcode handling: extraction, Hamming distance, correction."""
from __future__ import annotations

from config import BC_LEFT, BC_RIGHT, BC_RIGHT_SPACER, BC_CORRECT_MAX


def hamming_le(a: bytes, b: bytes, k: int) -> int:
    """Hamming distance between equal-length prefixes of a and b, early-exiting
    once it exceeds k. A return value > k means "more than k" (not exact)."""
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return d
    return d


def extract_bc(seq: bytes):
    """Extract the TAPE-BC = bytes between BC_LEFT (GGCATG) and the next BC_RIGHT (CTCTGG).

    The barcode is followed by a constant BC_RIGHT_SPACER ('G') before CTCTGG; it is
    stripped so the returned barcode is the true 12-bp sequence. cstart (the array
    start) is the index just past CTCTGG and is unaffected by the strip.

    Returns (bc_bytes, cstart), or (None, None) if no plausible barcode is found.
    The length window [5, 20] is deliberately loose; downstream code enforces the
    exact length (config.BC_LEN = 12).
    """
    l = seq.find(BC_LEFT)
    if l < 0:
        return None, None
    bcs = l + len(BC_LEFT)
    r = seq.find(BC_RIGHT, bcs)
    if r < 0 or not (5 <= r - bcs <= 20):
        return None, None
    bc = seq[bcs:r]
    if BC_RIGHT_SPACER and bc.endswith(BC_RIGHT_SPACER):
        bc = bc[:-len(BC_RIGHT_SPACER)]     # drop the constant spacer -> true barcode
    return bc, r + len(BC_RIGHT)


def correct(bc: bytes, whitelist: list[bytes], maxd: int = BC_CORRECT_MAX):
    """Correct an observed barcode to its unique nearest whitelist member.

    Returns the whitelist barcode iff there is a single nearest member within
    Hamming `maxd` (i.e. the second-nearest is strictly farther). Ambiguous or
    too-distant barcodes return None -- this is what rejects Hamming-1 "error
    shadows" of abundant barcodes rather than merging them.
    """
    best = None
    bd = maxd + 1        # best distance
    sd = bd              # second-best distance
    for w in whitelist:
        d = hamming_le(bc, w, maxd)
        if d < bd:
            sd = bd
            bd = d
            best = w
        elif d < sd:
            sd = d
    return best if (best is not None and bd <= maxd and sd > bd) else None
