#!/usr/bin/env python3
"""Stage 06 -- lineage-diversity rarefaction per integration.

Usage:  python3 scripts/06_rarefaction.py DTTz_3_S3

For each integration barcode, converts clean per-genotype read counts to cells
(round(reads/lambda)) and asks how the number of distinct lineage genotypes grows
with the number of cells sampled: Hurlbert interpolation up to the observed
sample size and Chao1 extrapolation beyond it (tape/rarefaction.py). Writes:
  tables/<tag>.rarefaction_curve.csv     tapebc, cells_sampled, expected_lineages, region
  tables/<tag>.rarefaction_summary.csv   tapebc, n_cells, S_obs, chao1, pct_captured, f1, f2
  data/figures/<tag>.rarefaction.png     accumulation curves (one per integration)
"""
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from tape.io import load_pickle
from tape.depth import estimate_lambda, cells_per_pattern
from tape.rarefaction import accumulation, chao1


def main(tag: str):
    emb = config.get_embryo(tag)
    src = emb.out("clean_patterns.pkl")
    if not src.exists():
        sys.exit(f"missing {src}; run stage 03 first.")
    clean = load_pickle(src)
    bcs = sorted(clean, key=lambda b: -sum(clean[b].values()))
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)

    curve_rows, summary_rows = [], []
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for b in bcs:
        reads = list(clean[b].values())
        lam = estimate_lambda(reads)
        counts = list(cells_per_pattern(clean[b], lam).values())   # cells per genotype
        if sum(counts) <= 0:
            continue
        mgrid, curve, n = accumulation(counts, config.RAREFACTION_EXTRAP,
                                       config.RAREFACTION_INTERP_PTS, config.RAREFACTION_EXTRAP_PTS)
        S_obs, S_chao, f1, f2, f0 = chao1(counts)
        pct = 100.0 * S_obs / S_chao if S_chao else 0.0
        summary_rows.append({"tapebc": b, "n_cells": n, "S_obs": S_obs,
                             "chao1": round(S_chao, 1), "pct_captured": round(pct, 1),
                             "f1": f1, "f2": f2})
        for m, s in zip(mgrid, curve):
            curve_rows.append({"tapebc": b, "cells_sampled": int(m),
                               "expected_lineages": round(float(s), 2),
                               "region": "observed" if m <= n else "extrapolated"})
        obs = mgrid <= n
        ax.plot(mgrid[obs], curve[obs], "-", color="#0C6B74", lw=1, alpha=.6)
        ax.plot(mgrid[~obs], curve[~obs], "--", color="#0C6B74", lw=1, alpha=.4)
        ax.plot(n, S_obs, "o", color="#16242E", ms=4, zorder=5)

    ax.set_xlabel("cells (genomic equivalents) sampled")
    ax.set_ylabel("expected # distinct lineage genotypes")
    ax.set_title(f"{emb.label} -- lineage-genotype rarefaction per integration\n"
                 f"(solid = observed, dashed = Chao1 extrapolation; dot = observed sample)",
                 fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(str(config.FIG_DIR / f"{tag}.rarefaction.png"), dpi=140)
    plt.close(fig)

    with open(config.TABLES_DIR / f"{tag}.rarefaction_curve.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["tapebc", "cells_sampled", "expected_lineages", "region"])
        w.writeheader(); w.writerows(curve_rows)
    with open(config.TABLES_DIR / f"{tag}.rarefaction_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["tapebc", "n_cells", "S_obs", "chao1",
                                          "pct_captured", "f1", "f2"])
        w.writeheader(); w.writerows(summary_rows)

    pj = np.mean([r["pct_captured"] for r in summary_rows]) if summary_rows else 0
    print(f"[06] {len(summary_rows)} integrations; mean fraction of estimated (Chao1) "
          f"richness captured = {pj:.0f}%")
    print(f"[06] wrote tables/{tag}.rarefaction_curve.csv + rarefaction_summary.csv "
          f"+ figures/{tag}.rarefaction.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
