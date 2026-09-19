#!/usr/bin/env python3
"""STEP 5 -- redraw of the A-vs-B replication scatter (fig. S7C), fixing two things in h3's version.

1. THREE CLASSES instead of two:
     - the 34 captured divisions (pooled FDR<1% AND pooled fold>=3, obs>=3)
     - 15 pairings significant at FDR<1% in BOTH halves that never reach 3-fold in both
     - the remaining 357 pairings testable in both halves
   (Note the second class is exactly "FDR in both, not 3-fold in both": none of those 15 reaches
   3-fold in both halves, so the class is disjoint from the captured set by construction.)

2. AXIS LIMITS THAT SHOW THE DATA. h3 set the lower limit to 0.55 while folds run down to 0.10,
   so 254 of the 406 plotted pairings fell outside the axes and were silently dropped -- while the
   legend quoted a correlation computed over all 406. Limits now span the full range.

The 3-fold guides are per-half references; the >=3-fold criterion is applied to the POOLED test, so
captured points may legitimately sit below a guide (5 of 34 do, and are labelled).

Reads figSX_captured_divisions_AB_replication.csv (shipped source data; no recomputation).
Writes figS7C_AB_replication_draft.png.
"""
import csv, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = f"{HERE}/figSX_captured_divisions_AB_replication.csv"
OUT  = f"{HERE}/figS7C_AB_replication_draft.png"
FDR, FOLD = 0.01, 3.0
LABEL_BELOW_GUIDE = True

INK, INK2, MUTED, GRID = "#1f2328", "#4a5158", "#8b9198", "#d7dade"
# validated categorical palette (light surface): all six checks pass
#   scripts/validate_palette.js "#c8553d,#e0a13c,#2f6db3" --mode light
RED, AMBER, BLUE = "#c8553d", "#e0a13c", "#2f6db3"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

rows = [r for r in csv.DictReader(open(SRC)) if r["tested_in_both"] == "1"]
f = lambda r, k: float(r[k])
xa = np.array([f(r, "fold_A") for r in rows]); xb = np.array([f(r, "fold_B") for r in rows])
cap = np.array([r["captured"] == "1" for r in rows])
sig_both = np.array([f(r, "q_A") < FDR and f(r, "q_B") < FDR for r in rows])
fold_both = np.array([f(r, "fold_A") >= FOLD and f(r, "fold_B") >= FOLD for r in rows])
c2 = (~cap) & sig_both & (~fold_both)          # FDR in both halves, not 3-fold in both
c3 = (~cap) & (~c2)
assert cap.sum() + c2.sum() + c3.sum() == len(rows)
rho = spearmanr(xa, xb); rlog = pearsonr(np.log2(xa), np.log2(xb))
print(f"tested in both {len(rows)}: captured {cap.sum()}, FDR-both-only {c2.sum()}, other {c3.sum()}")
print(f"Spearman rho={rho.statistic:.3f} p={rho.pvalue:.2e}; Pearson r(log2)={rlog.statistic:.3f}")
print(f"captured below a 3-fold guide: {int((cap & ~fold_both).sum())}")

fig, ax = plt.subplots(figsize=(6.6, 6.6))
lo = min(xa.min(), xb.min()) * 0.8; hi = max(xa.max(), xb.max()) * 1.35
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
# the region satisfying >=3-fold in BOTH halves, as a light field rather than bare guides
ax.add_patch(plt.Rectangle((FOLD, FOLD), hi - FOLD, hi - FOLD, facecolor="#f2f4f6",
                           edgecolor="none", zorder=0))
ax.plot([lo, hi], [lo, hi], ls="--", lw=0.9, color=MUTED, zorder=1)
ax.axhline(FOLD, color=GRID, lw=0.9, zorder=1); ax.axvline(FOLD, color=GRID, lw=0.9, zorder=1)
ax.axhline(1, color=GRID, lw=0.6, ls=":", zorder=1); ax.axvline(1, color=GRID, lw=0.6, ls=":", zorder=1)
# size + ring are the secondary encoding, so identity is never colour-alone
ax.scatter(xa[c3], xb[c3], s=16, c=BLUE, alpha=0.35, linewidths=0, zorder=2,
           label=f"tested in both halves ({int(c3.sum())})")
ax.scatter(xa[c2], xb[c2], s=42, c=AMBER, edgecolors="white", linewidths=1.1, zorder=3,
           label=f"FDR<1% in both halves, <3-fold in one ({int(c2.sum())})")
ax.scatter(xa[cap], xb[cap], s=62, c=RED, edgecolors="white", linewidths=1.4, zorder=4,
           label=f"captured division ({int(cap.sum())})")
# The 5 captured pairings that do not clear 3-fold in BOTH halves are drawn like the rest -- they
# are captured on the POOLED test, and the guides at 3 are per-half references. Listed here so the
# legend can say so without a fourth marker style.
if LABEL_BELOW_GUIDE:
    print("  captured but <3-fold in one half (plotted as ordinary captured points):")
    for i in np.flatnonzero(cap & ~fold_both):
        print(f"    A={xa[i]:>6.2f} B={xb[i]:>6.2f}  {rows[i]['progenitor']} -> {rows[i]['postmitotic']}")
ax.text(FOLD * 1.12, hi / 1.06, "≥3-fold in both halves", fontsize=7.2, color=MUTED, va="top")
ax.set_xlabel("fold enrichment, blastomere A (B1)")
ax.set_ylabel("fold enrichment, blastomere B (B2)")
ax.set_title("Independent replication across separately reconstructed half-embryos\n"
             f"Spearman ρ = {rho.statistic:.2f} (n = {len(rows)})   "
             f"Pearson r(log₂) = {rlog.statistic:.2f}", fontsize=10, color=INK, pad=10)
ax.legend(fontsize=7.6, frameon=False, loc="upper left")
ax.grid(alpha=0.18, which="both", lw=0.5); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT, dpi=170)
print("wrote", os.path.basename(OUT))
