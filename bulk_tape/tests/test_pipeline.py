#!/usr/bin/env python3
"""Reproducibility check for the bulk TAPE pipeline.

Runs stages 03 (de-noise), 04 (depth model + copy number), 05 (edit-chain
folding) and 06 (rarefaction) from the committed per-embryo `patterns_full`
inputs in `example_data/`, and asserts the published values. No raw fastq and no
network needed; the whole thing takes a few seconds.

    python3 tests/test_pipeline.py

Stages 01-02 (reference derivation and read parsing) need the raw fastqs, which
are not distributed with the code -- see README ("Raw data"). Their outputs are
what `refs/` and `example_data/` hold, and the reference files are checked
against their recorded contents here.
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

WORK = Path(tempfile.mkdtemp(prefix="bulk_tape_test_"))
os.environ["TAPE_DATA_DIR"] = str(WORK)

import config                                    # noqa: E402  (after TAPE_DATA_DIR)
from tape.io import load_pickle, load_whitelist, load_vocab   # noqa: E402
from tape.denoise import denoise_all             # noqa: E402
from tape.depth import estimate_lambda, genomic_equivalents, cells_per_pattern  # noqa: E402
from tape.copynumber import estimate_copies      # noqa: E402
from tape.editchain import build_trie, fold_dropout, tips     # noqa: E402
from tape.rarefaction import chao1               # noqa: E402

FAILURES = []


def check(label, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {label}: {got}" + ("" if ok else f"  (expected {want})"))
    if not ok:
        FAILURES.append(f"{label}: got {got}, expected {want}")


def load_clean(tag):
    """Stage 03 on the committed patterns_full input."""
    full = load_pickle(ROOT / "example_data" / f"{tag}.patterns_full.pkl")
    kept = {b: c for b, c in full.items() if sum(c.values()) >= config.MIN_BC_READS}
    clean, _funnels, _home = denoise_all(kept)
    return kept, clean


# --------------------------------------------------------------------------
# References (stage 01 output, checked in)
# --------------------------------------------------------------------------
print("\nreferences (refs/, derived from bulk by stage 01)")
for tag, n_bc, n_vocab in [("DTTz_2_S2", 17, 11), ("DTTz_3_S3", 11, 13), ("DTTz_6_S6", 5, 1)]:
    emb = config.get_embryo(tag)
    check(f"{tag} whitelist barcodes", len(load_whitelist(emb.whitelist_tsv)), n_bc)
    check(f"{tag} vocabulary insertions", len(load_vocab(emb.vocab_tsv)), n_vocab)

# the vocabulary threshold sits in a clean gap -- re-derived from refs/ so the
# figures quoted in config.py / tape/reference.py / README.md cannot go stale
for tag, n_adm, lo_adm, hi_rej in [("DTTz_2_S2", 11, 14589, 2350),
                                   ("DTTz_3_S3", 13, 72667, 5446),
                                   ("DTTz_6_S6", 1, 28813, 2997)]:
    rows = [l.rstrip("\n").split("\t")
            for l in open(config.get_embryo(tag).vocab_tsv)][1:]
    admitted = [int(r[3]) for r in rows if r[2] == "1"]
    rejected = [int(r[3]) for r in rows if r[2] == "0"]
    check(f"{tag} admitted insertions", len(admitted), n_adm)
    check(f"{tag} lowest admitted peak", min(admitted), lo_adm)
    check(f"{tag} highest rejected peak", max(rejected), hi_rej)
    check(f"{tag} threshold falls inside the gap",
          max(rejected) < config.VOCAB_PEAK_MIN <= min(admitted), True)

# every whitelisted barcode matches the NNNNNAANNNNN design
for tag in ("DTTz_2_S2", "DTTz_3_S3", "DTTz_6_S6"):
    wl = load_whitelist(config.get_embryo(tag).whitelist_tsv)
    check(f"{tag} barcodes are 12 bp with the AA spacer at 6-7",
          all(len(b) == config.BC_LEN and b[5:7] == b"AA" for b in wl), True)

# --------------------------------------------------------------------------
# Stage 03: de-noising funnel
# --------------------------------------------------------------------------
print("\nstage 03 -- de-noising")
kept3, clean3 = load_clean("DTTz_3_S3")
check("embryo 3 barcodes >= MIN_BC_READS", len(kept3), 11)
check("embryo 3 clean lineages", sum(len(v) for v in clean3.values()), 4989)

kept2, clean2 = load_clean("DTTz_2_S2")
check("embryo 2 barcodes >= MIN_BC_READS", len(kept2), 17)
check("embryo 2 clean lineages", sum(len(v) for v in clean2.values()), 1562)

kept6, clean6 = load_clean("DTTz_6_S6")
check("embryo 6 barcodes >= MIN_BC_READS", len(kept6), 5)
check("embryo 6 clean lineages", sum(len(v) for v in clean6.values()), 12)

# --------------------------------------------------------------------------
# Stage 04: depth model + copy number
# --------------------------------------------------------------------------
print("\nstage 04 -- depth model and copy number")
ge3 = {b: genomic_equivalents(c.values(), estimate_lambda(c.values())) for b, c in clean3.items()}
check("embryo 3 genomic-equivalent range", (min(ge3.values()), max(ge3.values())), (179, 366))

copies3, unit3, tot3 = estimate_copies({b: sum(c.values()) for b, c in clean3.items()})
check("embryo 3 total integrations", tot3, 11)
check("embryo 3 multi-copy barcodes", sum(1 for c in copies3.values() if c > 1), 0)

copies2, unit2, tot2 = estimate_copies({b: sum(c.values()) for b, c in clean2.items()})
check("embryo 2 total integrations", tot2, 35)
check("embryo 2 copy numbers > 1",
      sorted((b, c) for b, c in copies2.items() if c > 1),
      [("CGAAGAAAAGTG", 8), ("TAGGAAAAGTGG", 6), ("TGGCGAAAAAAC", 6), ("TTGGTAACCCCA", 2)])

# --------------------------------------------------------------------------
# Stage 05: edit-chain dropout folding
# --------------------------------------------------------------------------
print("\nstage 05 -- edit-chain trie and dropout folding")
folded = 0
for counts in clean3.values():
    root = build_trie(counts)
    fold_dropout(root, config.DROPOUT_FOLD_RATIO)
    folded += sum(1 for _ in tips(root))
check("embryo 3 lineages after folding (from 4,989)", folded, 4542)

# --------------------------------------------------------------------------
# Stage 06: rarefaction
# --------------------------------------------------------------------------
print("\nstage 06 -- rarefaction (Chao1)")
n_bc_curves = 0
for b, counts in clean3.items():
    cells = cells_per_pattern(counts, estimate_lambda(counts.values()))
    S, S_chao, f1, f2, f0 = chao1(list(cells.values()))
    if S > 0 and S_chao >= S:
        n_bc_curves += 1
check("embryo 3 integrations with a valid Chao1 curve", n_bc_curves, 11)

shutil.rmtree(WORK, ignore_errors=True)
print()
if FAILURES:
    print(f"{len(FAILURES)} CHECK(S) FAILED")
    sys.exit(1)
print("all checks passed")
