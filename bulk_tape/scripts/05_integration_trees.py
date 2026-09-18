#!/usr/bin/env python3
"""Stage 05 -- per-integration lineage trees (the endpoint).

Usage:  python3 scripts/05_integration_trees.py DTTz_3_S3 [--no-fold] [--cells]

Reads data/<tag>.clean_patterns.pkl and renders one lineage trie per TAPE
integration (barcode). By default it applies edit-chain dropout folding
(config.DROPOUT_FOLD_RATIO) so 3'-truncated reads are attributed to their deeper
parent lineage -- the single-cell-consistent handling. Flags:
  --no-fold   keep every observed depth as a distinct tip (raw pattern-tuple view)
  --cells     weight tips by estimated cells (round(reads/lambda)) instead of reads
Writes data/figures/<tag>.integration_trees.png
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import config
from tape.io import load_pickle
from tape.editchain import build_trie, fold_dropout, tips
from tape.depth import estimate_lambda, cells_per_pattern
from tape.tree import render_all, _pad


def _fold_counts(counts: dict, ratio: float) -> dict:
    """Return {padded-pattern: weight} after dropout folding a barcode's lineages."""
    root = build_trie(counts)
    fold_dropout(root, ratio)
    return {_pad(chain): w for chain, w in tips(root)}


def main(tag: str, fold: bool = True, cells: bool = False):
    emb = config.get_embryo(tag)
    src = emb.out("clean_patterns.pkl")
    if not src.exists():
        sys.exit(f"missing {src}; run stage 03 first.")
    clean = load_pickle(src)

    per_bc = {}
    for bc, counts in clean.items():
        c = counts
        if cells:
            lam = estimate_lambda(counts.values())
            c = cells_per_pattern(counts, lam)
            if not c:
                continue
        if fold:
            c = _fold_counts(c, config.DROPOUT_FOLD_RATIO)
        per_bc[bc] = c

    config.FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = config.FIG_DIR / f"{tag}.integration_trees.png"
    weight = "cells" if cells else "reads"
    subtitle = "edit-chain folded" if fold else "raw pattern tuples"
    render_all(per_bc, out, weight_name=weight,
               title=f"{emb.label} -- per-integration TAPE lineage trees ({subtitle})")
    n_lineages = sum(len(v) for v in per_bc.values())
    print(f"[05] rendered {len(per_bc)} integrations, {n_lineages:,} lineages "
          f"(weight={weight}, fold={fold})")
    print(f"[05] wrote {out}")


if __name__ == "__main__":
    args = sys.argv[1:]
    fold = "--no-fold" not in args
    cells = "--cells" in args
    args = [a for a in args if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    main(args[0], fold=fold, cells=cells)
