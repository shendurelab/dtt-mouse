#!/usr/bin/env python3
"""Export per-integration bulk TAPE lineage genotypes for tree building.

For each embryo-3 TAPE-BC integration, emit the set of distinct LINEAGE
GENOTYPES so a colleague can build a lineage tree for EACH TAPE IN ISOLATION
from the bulk amplicon data (no single-cell data used here).

A lineage genotype is the pipeline's canonical cell-weighted, dropout-folded
edit-chain tip (identical to `scripts/05_integration_trees.py --cells`):

  clean_patterns.pkl {tapebc: {6-site pattern: reads}}   (stage 03)
    -> lambda = single-genomic-equivalent read depth per barcode (tape/depth.py)
    -> n_cells(pattern) = round(reads / lambda), keep >= 1 cell
    -> edit-chain trie + fold_dropout (config.DROPOUT_FOLD_RATIO): 3'-truncated
       reads are attributed to their deeper parent lineage
    -> each surviving trie tip = one distinct lineage genotype (weight = n_cells)

Sites are written as six '|'-separated tokens with 'U' for unedited, matching the
tokenization consumed by the tree-building code (src/1-dist-mat).

Usage:  python3 scripts/export_bulk_lineage_genotypes.py [DTTz_3_S3]
Writes: tables/<tag>.bulk_lineage_genotypes.csv
        tables/<tag>.bulk_lineage_genotypes_summary.csv
"""
import sys, csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config
from tape.io import load_pickle
from tape.editchain import build_trie, fold_dropout, tips, edit_chain
from tape.depth import estimate_lambda, cells_per_pattern

N = config.N_SITES


def sites_token(chain) -> str:
    """Edit chain -> six '|'-separated site tokens, 'U' for unedited."""
    sites = list(chain) + ["U"] * (N - len(chain))
    return "|".join(sites[:N])


def main(tag: str = "DTTz_3_S3"):
    emb = config.get_embryo(tag)
    src = emb.out("clean_patterns.pkl")
    if not src.exists():
        sys.exit(f"missing {src}; run stages 01-03 first.")
    clean = load_pickle(src)

    tables = ROOT / "tables"
    tables.mkdir(exist_ok=True)
    out_long = tables / f"{tag}.bulk_lineage_genotypes.csv"
    out_sum = tables / f"{tag}.bulk_lineage_genotypes_summary.csv"

    long_rows, summary_rows = [], []
    bc_order = sorted(clean, key=lambda b: -sum(clean[b].values()))

    for bc in bc_order:
        counts = clean[bc]
        lam = estimate_lambda(counts.values())

        # read-weighted fold -> reads per surviving genotype (reference column)
        rroot = build_trie(counts)
        fold_dropout(rroot, config.DROPOUT_FOLD_RATIO)
        reads_by_chain = {ch: w for ch, w in tips(rroot)}
        raw_reads = {}
        for pat, w in counts.items():
            raw_reads[edit_chain(pat)] = raw_reads.get(edit_chain(pat), 0) + w

        # canonical cell-weighted, dropout-folded genotypes (== 05 --cells)
        cc = cells_per_pattern(counts, lam)
        croot = build_trie(cc)
        fold_dropout(croot, config.DROPOUT_FOLD_RATIO)
        geno = sorted(tips(croot), key=lambda cw: (-cw[1], cw[0]))

        for i, (chain, ncells) in enumerate(geno, 1):
            reads = reads_by_chain.get(chain, raw_reads.get(chain, 0))
            long_rows.append([
                bc, f"{bc}_g{i:03d}", sites_token(chain),
                len(chain), int(round(ncells)), int(round(reads)),
            ])
        summary_rows.append([
            bc, round(lam, 1), int(sum(counts.values())),
            len(geno), int(round(sum(w for _, w in geno))),
        ])

    with out_long.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tapebc", "genotype_id", "genotype",
                    "edit_depth", "n_cells", "n_reads"])
        w.writerows(long_rows)
    with out_sum.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tapebc", "lambda", "total_reads", "n_genotypes", "total_cells"])
        w.writerows(summary_rows)

    rng = (min(r[3] for r in summary_rows), max(r[3] for r in summary_rows))
    print(f"[export] {len(summary_rows)} integrations, {len(long_rows)} genotypes "
          f"(per-tapebc range {rng[0]}-{rng[1]})")
    print(f"[export] wrote {out_long.relative_to(ROOT)}")
    print(f"[export] wrote {out_sum.relative_to(ROOT)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "DTTz_3_S3")
