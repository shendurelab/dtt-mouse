"""Derive the per-embryo TAPE-BC whitelist and insertion vocabulary (BULK only).

Run once per bulk library on a sample of reads; the resulting whitelist + vocab
tsvs are then reused by every downstream stage AND by the single-cell pipeline
for the same embryo (bulk has deep, clean per-integration coverage, so it is the
right place to define the reference sets).

Vocabulary rule: an insertion is admitted iff at some single (TAPE-BC x site)
combination it is seen >= VOCAB_PEAK_MIN times. Real edits are clonal and pile up
at a specific integration+position, so genuine insertions peak far above error
(on embryo 3 the lowest admitted insertion peaks at 72,667 against 5,446 for the
highest rejected one, a 13.3x gap), while 1-bp error "shadows" of abundant
insertions stay below the threshold. See config.VOCAB_PEAK_MIN for the per-embryo
figures.
"""
from __future__ import annotations
from collections import Counter, defaultdict

from config import (SAMPLE_N, BC_LEN, BC_MIN_FREQ, BC_MIN_SEP, VOCAB_PEAK_MIN)
from tape.io import iter_reads
from tape.barcode import extract_bc, correct, hamming_le
from tape.parser import walk_robust


def derive(fastq_path, sample_n: int = SAMPLE_N, peak_min: int = VOCAB_PEAK_MIN):
    """Derive references.

    The whitelist is frequency-based, so it is built from a SAMPLE_N-read sample
    (fractions are stable). The vocabulary uses an ABSOLUTE per-(BC x position)
    peak count, which scales with sequencing depth, so it is measured over ALL
    reads -- otherwise a cold (lightly edited) library's peaks fall below the
    fixed threshold on a small sample.

    Returns (whitelist, vocab, bc_counts, ins_peak):
      whitelist : list of (barcode_bytes, sample_count), rank-ordered
      vocab     : set of admitted insertion byte-strings
      bc_counts : Counter of all extracted barcodes
      ins_peak  : {insertion_bytes: max count at any single (BC, position)}
    """
    # ---- pass 1: barcode whitelist (frequency + Hamming separation), sampled ----
    bc_counts = Counter()
    for seq in iter_reads(fastq_path, sample_n):
        bc, _ = extract_bc(seq)
        if bc is not None and len(bc) == BC_LEN:
            bc_counts[bc] += 1
    tot_bc = sum(bc_counts.values())
    whitelist = []
    for bc, c in bc_counts.most_common():
        if tot_bc == 0 or c / tot_bc < BC_MIN_FREQ:
            break
        if all(hamming_le(bc, w, BC_MIN_SEP - 1) >= BC_MIN_SEP for w, _ in whitelist):
            whitelist.append((bc, c))
    wl = [w for w, _ in whitelist]

    # ---- pass 2: per-(BC, position) insertion counts over ALL reads -> peak ----
    cell = defaultdict(Counter)          # (bc_bytes, pos) -> Counter(insertion_bytes)
    for seq in iter_reads(fastq_path):
        bc, cs = extract_bc(seq)
        if bc is None or len(bc) != BC_LEN:
            continue
        corr = correct(bc, wl)
        if corr is None:
            continue
        state, err, _ = walk_robust(seq, cs)
        if err:
            continue
        for pos, s in enumerate(state):
            if s is not None and s != "U":
                cell[(corr, pos)][s[1]] += 1

    ins_peak = {}
    for counts in cell.values():
        for ins, k in counts.items():
            if k > ins_peak.get(ins, 0):
                ins_peak[ins] = k
    vocab = {ins for ins, pk in ins_peak.items() if pk >= peak_min}
    return whitelist, vocab, bc_counts, ins_peak


def write_whitelist(path, whitelist, tot_bc: int, corrected: Counter | None = None):
    with open(path, "w") as f:
        f.write("rank\ttapebc\tsample_count\tsample_freq\tfull_corrected_count\n")
        for i, (bc, c) in enumerate(whitelist, 1):
            cc = corrected.get(bc, 0) if corrected else 0
            f.write(f"{i}\t{bc.decode()}\t{c}\t{c / max(1, tot_bc):.5f}\t{cc}\n")


def write_vocab(path, ins_peak: dict, vocab: set):
    """Write every observed insertion with its per-(BC,position) peak count and
    an in_vocab flag (peak >= VOCAB_PEAK_MIN)."""
    with open(path, "w") as f:
        f.write("nnn\tlen\tin_vocab\tbcpos_peak\n")
        for nnn, pk in sorted(ins_peak.items(), key=lambda kv: -kv[1]):
            f.write(f"{nnn.decode()}\t{len(nnn)}\t{int(nnn in vocab)}\t{pk}\n")
