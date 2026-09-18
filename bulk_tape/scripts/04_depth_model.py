#!/usr/bin/env python3
"""Stage 04 -- read-depth model: estimate lambda and genomic equivalents per barcode.

Usage:  python3 scripts/04_depth_model.py DTTz_3_S3

Reads data/<tag>.clean_patterns.pkl, estimates the single-genomic-equivalent
read depth lambda (log-mode) per barcode, and reports genomic equivalents
(cells) = sum(round(reads/lambda)). Writes:
  data/<tag>.depth_model.tsv           per-barcode lambda + cell estimates
  data/figures/<tag>.depth_model.png   per-barcode log read-depth histograms
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from tape.io import load_pickle
from tape.depth import estimate_lambda, genomic_equivalents
from tape.copynumber import estimate_copies


def main(tag: str):
    emb = config.get_embryo(tag)
    src = emb.out("clean_patterns.pkl")
    if not src.exists():
        sys.exit(f"missing {src}; run stage 03 first.")
    clean = load_pickle(src)
    bcs = sorted(clean, key=lambda b: -sum(clean[b].values()))
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)

    ncol = 4
    nrow = (len(bcs) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.3 * nrow))
    axes = np.atleast_1d(axes).ravel()

    rows = []
    for k, b in enumerate(bcs):
        reads = np.array(sorted(clean[b].values(), reverse=True), float)
        lam = estimate_lambda(reads)
        N = genomic_equivalents(reads, lam)
        rows.append((b, len(reads), int(reads.sum()), lam, float(np.median(reads)), N))
        ax = axes[k]
        ax.hist(np.log10(reads), bins=30, color="#0C6B74", alpha=.85)
        ax.axvline(np.log10(lam), color="#D1495B", lw=1.5, label=f"lambda={lam:.0f}")
        for m in (2, 3, 4):
            ax.axvline(np.log10(lam * m), color="#D1495B", ls=":", lw=.8)
        ax.set_title(f"{b}\nN~{N:,} cells", fontsize=8)
        ax.set_xlabel("log10 reads/lineage")
        ax.set_ylabel("# lineages")
        ax.legend(fontsize=6)
    for j in range(len(bcs), len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"{emb.label} -- read-depth per clean lineage; lambda = single genomic "
                 f"equivalent (red), multiplets at 2x/3x/4x (dotted)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(config.FIG_DIR / f"{tag}.depth_model.png"), dpi=130)
    plt.close(fig)

    # ---- genomic copy number per integration barcode ----
    # per-BC total reads = copy-number readout (see tape/copynumber.py); cells =
    # genomic_equivalents / copies (a multi-copy barcode over-counts cells).
    per_bc_reads = {b: sum(clean[b].values()) for b in bcs}
    copies, unit, total_integrations = estimate_copies(per_bc_reads)

    with open(emb.out("depth_model.tsv"), "w") as f:
        f.write("tapebc\tn_lineages\treads\tlambda\tmedian_reads\tgenomic_equivalents\t"
                "copies\tcells\n")
        for b, npat, tot, lam, med, N in rows:
            cp = copies[b]
            f.write(f"{b}\t{npat}\t{tot}\t{lam:.1f}\t{med:.0f}\t{N}\t{cp}\t{int(round(N/cp))}\n")

    Ns = [r[5] for r in rows]
    n_multi = sum(1 for c in copies.values() if c > 1)
    print(f"[04] genomic equivalents/barcode: mean {np.mean(Ns):.0f}, range {min(Ns)}-{max(Ns)}")
    print(f"[04] copy number: single-copy unit ~{unit:,.0f} reads; "
          f"{n_multi} multi-copy barcode(s); "
          f"total integrations = {total_integrations} across {len(bcs)} valid barcodes")
    if n_multi:
        for b in sorted(bcs, key=lambda b: -copies[b]):
            if copies[b] > 1:
                print(f"        {b}: {copies[b]} copies ({per_bc_reads[b]/unit:.1f}x unit)")
    print(f"[04] wrote {emb.out('depth_model.tsv')} + figures/{tag}.depth_model.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
