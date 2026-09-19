#!/usr/bin/env python3
"""STEP 18 -- fig. S8C: K=15 vs K=200 co-occurrence enrichment, blastomere A and B side by side.

c2_blocks.py drew this as ONE panel using the mean of A and B. This draws the two independently
reconstructed half-embryos as separate panels, each with its own enrichments on both axes, while the
colour still encodes the joint call (significant in BOTH halves) that defines the reported counts.
So the panels show the raw agreement between halves and the colour shows what was called.

Significance, as in c1/c2: z = (obs - exp)/sqrt(exp) >= 3 and obs >= 5, in both blastomeres, at that K.

Palette note: the "both K" category was #ffd400, which sits outside the usable lightness band and
reaches only 1.39:1 against a white surface -- nearly invisible for the 8 points that matter most.
Replaced with #a3195b; validated with the dataviz palette checker (all checks pass).

Writes out/c18_s8c_k15_vs_k200_AB.png and .csv.
"""
import csv, os, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Z_THR, MIN_OBS = 3.0, 5
INK, INK2, MUTED = "#1f2328", "#4a5158", "#8b9198"
GREY, ORANGE, BLUE, CRIMSON = "#c9ced3", "#ff7a1a", "#1f6feb", "#a3195b"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

def load(bl, K):
    d = {}
    for r in csv.DictReader(open(f"{HERE}/out/cooccur_{bl}_K{K}_final.csv")):
        if r["celltype_1"] == r["celltype_2"]: continue
        d[tuple(sorted((r["celltype_1"], r["celltype_2"])))] = (
            float(r["obs"]), float(r["exp"]), float(r["log2_enr"]))
    return d
D = {(bl, K): load(bl, K) for bl in ("B1", "B2") for K in (15, 200)}
pairs = sorted(set(D[("B1", 15)]) & set(D[("B2", 15)]) & set(D[("B1", 200)]) & set(D[("B2", 200)]))
def sig(bl, K, k):
    o, e, _ = D[(bl, K)][k]
    return (o - e) / math.sqrt(max(e, 1e-9)) >= Z_THR and o >= MIN_OBS
s15 = np.array([sig("B1", 15, k) and sig("B2", 15, k) for k in pairs])
s200 = np.array([sig("B1", 200, k) and sig("B2", 200, k) for k in pairs])
cat = np.where(s15 & s200, 3, np.where(s15, 1, np.where(s200, 2, 0)))
print(f"{len(pairs)} off-diagonal pairs: K=15 {int(s15.sum())}, K=200 {int(s200.sum())}, "
      f"both {int((s15&s200).sum())}, K=15 only {int((cat==1).sum())}, "
      f"K=200 only {int((cat==2).sum())}, neither {int((cat==0).sum())}")

SPEC = [(0, GREY, 12, None, f"not enriched in both halves ({int((cat==0).sum())})"),
        (1, ORANGE, 40, "white", f"K = 15 only ({int((cat==1).sum())})"),
        (2, BLUE, 40, "white", f"K = 200 only ({int((cat==2).sum())})"),
        (3, CRIMSON, 56, "white", f"both K ({int((cat==3).sum())})")]
fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.9), sharex=True, sharey=True)
allv = [D[(bl, K)][k][2] for bl in ("B1", "B2") for K in (15, 200) for k in pairs]
lim = [min(allv) - 0.3, max(allv) + 0.3]
for ax, bl, nm in zip(axes, ("B1", "B2"), ("A (B1)", "B (B2)")):
    x = np.array([D[(bl, 200)][k][2] for k in pairs])
    y = np.array([D[(bl, 15)][k][2] for k in pairs])
    ax.plot(lim, lim, ls="--", lw=0.9, color=MUTED, zorder=1)
    ax.axhline(0, color=MUTED, lw=0.5, zorder=1); ax.axvline(0, color=MUTED, lw=0.5, zorder=1)
    for c, col, s, ec, lab in SPEC:
        m = cat == c
        ax.scatter(x[m], y[m], s=s, c=col, edgecolor=ec, lw=0.7 if ec else 0,
                   zorder=2 + c, label=lab if bl == "B1" else None)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("log₂ enrichment, K = 200  (ancient clades)")
    ax.set_title(f"blastomere {nm}", fontsize=11, color=INK, pad=6)
    ax.grid(alpha=.22, lw=.5); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("log₂ enrichment, K = 15  (recent clades)")
axes[0].legend(fontsize=8, frameon=False, loc="upper left")
fig.suptitle("Cell-type pair co-occurrence at two clade scales, in each independently "
             "reconstructed half-embryo", fontsize=11.5, color=INK, y=0.99)
fig.tight_layout(); fig.savefig(f"{HERE}/out/c18_s8c_k15_vs_k200_AB.png", dpi=160, bbox_inches="tight")
NAME = {0: "neither", 1: "K15_only", 2: "K200_only", 3: "both_K"}
with open(f"{HERE}/out/c18_s8c_k15_vs_k200_AB.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["celltype_1", "celltype_2", "category",
                                    "log2_enr_A_K15", "log2_enr_A_K200",
                                    "log2_enr_B_K15", "log2_enr_B_K200"])
    for i, k in enumerate(pairs):
        w.writerow([k[0], k[1], NAME[cat[i]],
                    f"{D[('B1',15)][k][2]:.4f}", f"{D[('B1',200)][k][2]:.4f}",
                    f"{D[('B2',15)][k][2]:.4f}", f"{D[('B2',200)][k][2]:.4f}"])
print("wrote c18_s8c_k15_vs_k200_AB.png / .csv")
