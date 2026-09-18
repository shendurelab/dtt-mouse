#!/usr/bin/env python3
"""Single-cell circTAPE recovery per integration (Fig 2 / SC results paragraph).

Consumes the committed single-cell outputs (no read re-parsing):
  <tag>.cell_consensus.tsv  cell x integration consensus matrix (sc_consensus.py)
  <tag>.sc_qc.tsv           per-cell QC incl. total_molecules
  <CELL_METADATA>           transcriptome cell table -> recovery denominator

For each integration barcode it counts how many cells carry a consensus genotype
and expresses that as a fraction of ALL transcriptome cells (the metadata row
count, e.g. 213,069) -- the denominator used in the manuscript. Also reports the
molecules-per-cell distribution (mean/median) over cells with >=1 consensus.

Usage:  python3 scripts/sc_recovery.py            # tag defaults to "e3"
        python3 scripts/sc_recovery.py <tag>

Writes to tables/:
  <tag>.sc_recovery_by_integration.csv  integration, n_cells_recovered, recovery_pct
                                        (+ header: denominator, cells-with-consensus,
                                         molecules/cell mean+median)
"""
import sys, csv, statistics as st
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

QC_COLS = {"n_loci", "n_doublet_loci", "mean_dominance"}


def load_metadata_cells(path):
    """Cell ids from the transcriptome metadata (col 0), skipping the header."""
    with open(path) as f:
        next(f)
        return {ln.split("\t", 1)[0] for ln in f if ln.strip()}


def main(tag: str):
    mat = config.DATA_DIR / f"{tag}.cell_consensus.tsv"
    qc = config.DATA_DIR / f"{tag}.sc_qc.tsv"
    meta = config.RAW_DIR / config.CELL_METADATA
    for p in (mat, qc, meta):
        if not p.exists():
            sys.exit(f"missing {p}")

    meta_cells = load_metadata_cells(meta)
    denom = len(meta_cells)

    with open(mat) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    hdr, data = rows[0], rows[1:]
    bc_cols = [c for c in hdr[1:] if c not in QC_COLS]
    idx = {c: hdr.index(c) for c in bc_cols}
    cell_i = 0

    recovered = {b: 0 for b in bc_cols}
    n_mat = 0
    n_outside = 0
    for r in data:
        n_mat += 1
        if r[cell_i] not in meta_cells:
            n_outside += 1
        for b in bc_cols:
            # "recovered" = a consensus genotype was called at this integration,
            # INCLUDING cells whose recorder was unedited there ("-"): that is a
            # real consensus observation, not a recovery failure. Only "" (no
            # consensus at all) is a miss.
            if r[idx[b]] != "":
                recovered[b] += 1

    # molecules per cell, from QC. Averaged over ALL transcriptome cells (the 213,069
    # denominator): cells with no recovered genotype contribute 0, matching the
    # recovery denominator. (Over recovered cells only it is slightly higher.)
    with open(qc) as f:
        qrows = list(csv.reader(f, delimiter="\t"))
    tmi = qrows[0].index("total_molecules")
    mol_recov = [int(r[tmi]) for r in qrows[1:]]
    mol = mol_recov + [0] * (denom - len(mol_recov))

    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out = config.TABLES_DIR / f"{tag}.sc_recovery_by_integration.csv"
    ordered = sorted(bc_cols, key=lambda b: -recovered[b])
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"# circTAPE recovery per integration; denominator = {denom:,} "
                    f"transcriptome cells ({config.CELL_METADATA})"])
        w.writerow([f"# cells with >=1 consensus = {n_mat:,} ({n_mat/denom*100:.1f}%); "
                    f"cells in matrix not in metadata = {n_outside}"])
        w.writerow([f"# molecules/cell over all {denom:,} cells (0 if none recovered): "
                    f"mean {st.mean(mol):.1f}, median {int(st.median(mol))}; "
                    f"over recovered cells only: mean {st.mean(mol_recov):.1f}, "
                    f"median {int(st.median(mol_recov))}"])
        w.writerow(["integration", "n_cells_recovered", "recovery_pct"])
        for b in ordered:
            w.writerow([b, recovered[b], round(recovered[b] / denom * 100, 1)])

    pcts = sorted((recovered[b] / denom * 100 for b in bc_cols), reverse=True)
    good = [p for p in pcts if p >= 50]
    print(f"[sc_recovery] denom={denom:,}  cells-with-consensus={n_mat:,}  "
          f"outside-metadata={n_outside}")
    print(f"[sc_recovery] molecules/cell over all {denom:,}: mean {st.mean(mol):.1f} "
          f"median {int(st.median(mol))} (recovered-only mean {st.mean(mol_recov):.1f})")
    print(f"[sc_recovery] {len(good)}/{len(pcts)} integrations >=50% "
          f"(range {min(good):.0f}-{max(good):.0f}%); "
          f"low: {', '.join(f'{p:.1f}%' for p in pcts if p < 50)}")
    print(f"[sc_recovery] wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "e3")
