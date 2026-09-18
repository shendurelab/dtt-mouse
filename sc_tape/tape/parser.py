"""Robust TAPE array parser.

`walk_robust` resolves each of the 6 monomer sites as independently as possible:
  - FORWARD from the 5' anchor (high-quality end) with ordered-U inference for
    the lightly-edited majority;
  - BACKWARD from the 3' terminal landmark for heavily-edited reads where the
    forward read-through truncates before reaching the deep sites.

This is the piece SHARED VERBATIM between the bulk and single-cell pipelines --
the read body downstream of the barcode is identical in both assays, so only
barcode extraction (see barcode.py) differs between them.

KNOWN ISSUE (open; affects both assays) -- see the block above the backward pass
in `walk_robust`: the two passes index the anchor list from opposite ends, so a
read that loses its FIRST anchor is parsed one site out of frame by the forward
pass and correctly registered by the backward pass. The two then read the same
junction into adjacent slots, duplicating the last symbol and dropping site 1.
Repro + prevalence estimate: tests/test_known_issue_frame_shift.py.
"""
from __future__ import annotations

from config import (ANCHOR, TERM, TERM_LONG, N_SITES, GAP_MIN, GAP_MAX)
from tape.barcode import hamming_le

# Site-state values returned by the parser:
#   "U"            -- site is unedited
#   ("E", b"NNN")  -- site is edited with insertion NNN
#   None           -- site could not be resolved


def find_terminal(seq: bytes, start: int) -> int:
    """Locate the 3' terminal monomer (end of the array). Exact match on the
    10-bp landmark first; else a fuzzy scan for the 18-bp landmark within 2
    mismatches. Returns the index, or -1 if not found."""
    p = seq.find(TERM, start)
    if p >= 0:
        return p
    nl = len(TERM_LONG)
    end = len(seq) - nl + 1
    for i in range(start, end):
        if hamming_le(seq[i:i + nl], TERM_LONG, 2) <= 2:
            return i
    return -1


def _valid_gap(g: int) -> bool:
    """True if the spacing between two consecutive anchors is a valid monomer step."""
    return GAP_MIN <= g <= GAP_MAX


def classify_junction(seg: bytes):
    """Classify the bytes between two monomer bodies.

    Returns ("U",) for an unedited junction (trailing 'TGAT' with nothing before),
    ("E", nnn) for an edit (insertion NNN + the 'GGA' key + 'TGAT'), or None on
    a corrupt/unrecognized junction.
    """
    if not seg.endswith(b"TGAT"):
        return None
    block = seg[:-4]
    if block == b"":
        return ("U",)
    if block.endswith(b"GGA"):
        nnn = block[:-3]
        if 1 <= len(nnn) <= 6:
            return ("E", nnn)
    return None


def walk_robust(seq: bytes, cstart: int):
    """Resolve the 6-site array from a read.

    Returns (state, err, resolved_all):
      state        -- list of length N_SITES; each entry "U" | ("E", nnn) | None
      err          -- True if an ordering violation was detected (edit after a U)
      resolved_all -- True if all 6 sites were resolved with no error
    """
    state = [None] * N_SITES

    # All body/terminal anchor positions from the array start.
    anchors = []
    i = seq.find(ANCHOR, cstart)
    while i >= 0:
        anchors.append(i)
        i = seq.find(ANCHOR, i + 6)
    if len(anchors) < 2:
        return state, False, False

    # Locate the terminal anchor (3' bound). If the fuzzy landmark is not already
    # an exact anchor, splice its position in.
    tpos = find_terminal(seq, cstart)
    term_idx = None
    if tpos >= 0:
        for idx, a in enumerate(anchors):
            if abs(a - tpos) <= 1:
                term_idx = idx
                break
        if term_idx is None:
            anchors = sorted([a for a in anchors if a < tpos - 1] + [tpos])
            term_idx = anchors.index(tpos)

    # ---- FORWARD: body1 = anchors[0]; junction k spans anchors[k]..anchors[k+1] = site k+1
    saw_u = False
    for k in range(min(len(anchors) - 1, N_SITES)):
        if not _valid_gap(anchors[k + 1] - anchors[k]):
            break                                   # missing/garbled anchor
        st = classify_junction(seq[anchors[k] + 10:anchors[k + 1]])
        if st is None:
            break
        if st[0] == "U":
            for j in range(k, N_SITES):             # ordered: everything after is unedited
                state[j] = "U"
            saw_u = True
            break
        state[k] = st

    # ---- BACKWARD from the terminal: fill remaining high sites for deep reads
    #
    # KNOWN ISSUE (open -- flagged, not fixed). The forward pass numbers sites
    # from anchors[0]; this pass numbers them back from `term_idx`. The two
    # agree only when the anchor list holds exactly N_SITES entries ahead of the
    # terminal. If the FIRST anchor fails to match (one mismatch in the 6-bp
    # ANCHOR is enough), every remaining gap is still valid, so the forward pass
    # runs to completion one site out of frame -- writing site 2's symbol into
    # slot 0, ..., site 6's into slot 4 -- and leaves slot 5 empty. This pass
    # then fills slot 5 from the *same* a6->TERM junction the forward pass just
    # consumed. Result: site 1's symbol is lost and the last symbol appears
    # twice, e.g. (ACT,GCC,AAG,CCC,ACG,GAA) -> (GCC,AAG,CCC,ACG,GAA,GAA).
    #
    # The read survives every downstream check: no ordering violation, all six
    # sites resolved, every symbol in the vocabulary. De-noising cannot catch it
    # either -- a frame shift is not a one-symbol neighbour (collapse), not a
    # prefix/suffix join (dechimera), and it keeps its own integration barcode
    # (cross-BC filter).
    #
    # A missing INTERIOR anchor is safe: it leaves an over-long gap that both
    # passes reject (_valid_gap), so the read is dropped as unresolved. A
    # truncated read with no terminal landmark is also safe -- this pass never
    # runs, so the shifted read fails all-sites-resolved and is discarded. It is
    # specifically the first anchor, on a read that does reach the terminal.
    #
    # A partially-edited array shifts the same way with no duplicate to betray
    # it: the forward U-fill sets saw_u, this pass never runs, and the genotype
    # is silently shifted and one site too shallow.
    #
    # Prevalence on embryo #3: the 5'-shift-plus-duplicate signature covers
    # ~2.6% of fully-edited reads against a ~2.0% chance-match background from
    # the mirrored (parser-impossible) 3'-shift, i.e. an excess of roughly 0.5%,
    # concentrated in genotypes at a few percent of their parent's depth. A
    # lower bound: the partially-edited variant above leaves no signature.
    #
    # Fix direction when revisited: register both passes against the same frame
    # -- require term_idx == N_SITES before trusting the backward walk, or make
    # the forward pass anchor on the CTCTGG/lead-in offset rather than on
    # anchors[0], so a missing a1 leaves a detectable hole instead of a shift.
    if not saw_u and term_idx is not None and term_idx >= 1:
        for j in range(N_SITES):                    # j=0 -> site6, j=1 -> site5, ...
            s = N_SITES - 1 - j
            if state[s] is not None:
                break
            ai = term_idx - j                       # anchor just 3' of site (s+1)
            if ai - 1 < 0 or ai >= len(anchors):
                break
            if not _valid_gap(anchors[ai] - anchors[ai - 1]):
                break
            st = classify_junction(seq[anchors[ai - 1] + 10:anchors[ai]])
            if st is None:
                break
            state[s] = "U" if st[0] == "U" else st

    # ---- ordering-consistency check: no edit may follow an unedited site
    seen_u = False
    for s in state:
        if s is None:
            continue
        if s == "U":
            seen_u = True
        elif seen_u:
            return state, True, False               # ordering violation

    return state, False, all(x is not None for x in state)


def pattern_of(state) -> tuple:
    """Convert a resolved 6-site state into a hashable pattern tuple of strings:
    'U' for unedited, or the decoded insertion (e.g. 'GCC') for an edit."""
    return tuple("U" if s == "U" else s[1].decode() for s in state)
