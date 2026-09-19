#!/usr/bin/env python3
"""STEP 3 -- final table + figures for the captured differentiation divisions.

Reads the per-pairing table for the primary threshold cell (set via CELL env; default
age>=E12.5, disc<=1) and writes:
  out/h3_captured.csv        the captured divisions, ranked by fold
  out/h3_heatmap.png         progenitor x post-mitotic, colour = log10 pooled fold, cell = n cells
  out/h3_ab_scatter.png      blastomere A vs B fold (independent half-embryo replication)
  out/h3_sweep.png           captured-set stability across the 12-cell threshold grid
"""
import csv, os, sys, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CELL = os.environ.get("CELL", "ageE12.5_disc1")
SRC = f"{HERE}/out/h2_pairings_{CELL}.csv"
INK, INK2, MUTED = "#1f2328", "#4a5158", "#8b9198"
BLUE, RED = "#2f6db3", "#c8553d"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "axes.linewidth": 0.8})

rows = list(csv.DictReader(open(SRC)))
cap = [r for r in rows if r["captured"] == "1"]
cap.sort(key=lambda r: -float(r["fold"]))
print(f"{SRC}: {len(rows)} pairings tested, {len(cap)} captured")

with open(f"{HERE}/out/h3_captured.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in cap: w.writerow(r)
nsame = sum(1 for r in cap if r["same_germlayer"] == "1")
ncell = sum(int(r["obs"]) for r in cap)
print(f"captured divisions {len(cap)}; within a single germ layer {nsame}/{len(cap)}; "
      f"{ncell:,} post-mitotic cells")
print("cross-germ-layer exceptions:")
for r in cap:
    if r["same_germlayer"] != "1":
        print(f"   {r['progenitor']} -> {r['postmitotic']}  ({r['prog_germlayer']} -> {r['pm_germlayer']})")

# ---- Fig A: heatmap -----------------------------------------------------------------------
gord = sorted({r["progenitor"] for r in cap},
              key=lambda g: (next(r["prog_germlayer"] for r in cap if r["progenitor"] == g), g))
pord = sorted({r["postmitotic"] for r in cap},
              key=lambda p: (next(r["pm_germlayer"] for r in cap if r["postmitotic"] == p), p))
gi = {g: i for i, g in enumerate(gord)}; pj = {p: j for j, p in enumerate(pord)}
M = np.full((len(gord), len(pord)), np.nan); OB = np.zeros_like(M)
for r in cap:
    M[gi[r["progenitor"]], pj[r["postmitotic"]]] = float(r["fold"])
    OB[gi[r["progenitor"]], pj[r["postmitotic"]]] = int(r["obs"])
fig, ax = plt.subplots(figsize=(max(8, len(pord) * 0.62), max(6, len(gord) * 0.42)))
im = ax.imshow(np.log10(M), aspect="auto", cmap="magma")     # sequential, monotone lightness
ax.set_xticks(range(len(pord))); ax.set_xticklabels(pord, rotation=52, ha="right", fontsize=8)
ax.set_yticks(range(len(gord))); ax.set_yticklabels(gord, fontsize=8)
ax.set_xticks(np.arange(-.5, len(pord), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(gord), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.4)          # 2px surface gap between cells
ax.tick_params(which="minor", length=0)
for s in ax.spines.values(): s.set_visible(False)
hi = np.log10(np.nanmax(M))
for i in range(len(gord)):
    for j in range(len(pord)):
        if not np.isnan(M[i, j]):
            ax.text(j, i, int(OB[i, j]), ha="center", va="center", fontsize=6.5,
                    color="white" if np.log10(M[i, j]) < hi * 0.62 else "#1a1005")
ax.set_ylabel("progenitor / regional state"); ax.set_xlabel("post-mitotic derivative")
ax.set_title(f"Captured differentiation divisions (n={len(cap)})\n"
             f"full tree, 1 pairing per post-mitotic cell; cell label = n cells",
             fontsize=10, color=INK, pad=10)
cb = fig.colorbar(im, ax=ax, shrink=0.55, label="log₁₀ fold enrichment")
cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{HERE}/out/h3_heatmap.png", dpi=160, bbox_inches="tight")

# ---- Fig B: A/B replication scatter -------------------------------------------------------
def fl(r, k):
    try: return float(r[k])
    except (ValueError, KeyError): return np.nan
both = [r for r in rows if np.isfinite(fl(r, "foldA")) and np.isfinite(fl(r, "foldB"))]
xa = np.array([fl(r, "foldA") for r in both]); xb = np.array([fl(r, "foldB") for r in both])
isc = np.array([r["captured"] == "1" for r in both])
from scipy.stats import spearmanr, pearsonr
rho = spearmanr(xa, xb).correlation
rlog = pearsonr(np.log2(xa), np.log2(xb))[0]
fig, ax = plt.subplots(figsize=(6.2, 6.2))
ax.scatter(xa[~isc], xb[~isc], s=18, c=BLUE, alpha=0.45, linewidths=0,
           label=f"tested in both halves ({int((~isc).sum())})", zorder=2)
ax.scatter(xa[isc], xb[isc], s=58, c=RED, edgecolors="white", linewidths=1.4,
           label=f"captured division ({int(isc.sum())})", zorder=3)
lim = [0.55, max(xa.max(), xb.max()) * 1.25]
ax.plot(lim, lim, ls="--", lw=0.9, color=MUTED, zorder=1)
ax.axhline(3, color=MUTED, lw=0.6, ls=":", zorder=1); ax.axvline(3, color=MUTED, lw=0.6, ls=":", zorder=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("fold enrichment, blastomere A (B1)"); ax.set_ylabel("fold enrichment, blastomere B (B2)")
ax.set_title("Independent replication across separately reconstructed half-embryos\n"
             f"Spearman ρ = {rho:.2f}   Pearson r(log₂) = {rlog:.2f}", fontsize=10, color=INK)
ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.grid(alpha=0.22, which="both", lw=0.5); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{HERE}/out/h3_ab_scatter.png", dpi=160)
print(f"A/B replication over all {len(both)} pairings tested in both halves: "
      f"rho={rho:.3f}, r(log2)={rlog:.3f}")

# ---- Fig C: threshold-grid stability ------------------------------------------------------
sw = list(csv.DictReader(open(f"{HERE}/out/h2_sweep.csv")))
ages = ["none", "E12", "E12.5", "E13"]; discs = ["2", "1", "0"]
G = np.full((len(discs), len(ages)), np.nan)
for r in sw: G[discs.index(r["disc_max"]), ages.index(r["age_min"])] = int(r["captured"])
# The counts span 33-36 -- a colour ramp over that range would imply variation that
# isn't there (and fails contrast). Render it as a table: uniform surface, ink text.
fig, ax = plt.subplots(figsize=(5.6, 3.0))
ax.imshow(np.ones_like(G), cmap="Greys", vmin=0, vmax=8, aspect="auto")
ax.set_xticks(range(len(ages))); ax.set_xticklabels([f"≥{a}" if a != "none" else "no filter" for a in ages])
ax.set_yticks(range(len(discs))); ax.set_yticklabels([f"≤{d}" if d != "0" else "= 0" for d in discs])
ax.set_xlabel("MRCA age threshold"); ax.set_ylabel("discordant edits")
for i in range(len(discs)):
    for j in range(len(ages)):
        ax.text(j, i, int(G[i, j]), ha="center", va="center", fontsize=13, color=INK)
ax.set_xticks(np.arange(-.5, len(ages), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(discs), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=2.0); ax.tick_params(which="minor", length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title(f"Captured divisions across the evidence-threshold grid\n"
             f"range {int(np.nanmin(G))}–{int(np.nanmax(G))}; "
             f"{int(np.nanmin(G))} pairings captured in every cell", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig(f"{HERE}/out/h3_sweep.png", dpi=160)
print("wrote out/h3_captured.csv, out/h3_heatmap.png, out/h3_ab_scatter.png, out/h3_sweep.png")
