#!/usr/bin/env python3
"""STEP 8 -- final figures for the time sweep: the carpet (all rows labelled) and a dendrogram of
cell-type relationships whose branch points carry real developmental dates.

COUPLING DEPTH. For each cell-type pair, the earliest time slice at which the pair is still
significantly co-enriched (z>=3, obs>=5, in both blastomeres). Pairs never significant at any
slice are censored at E13.5, just past the last slice.

WHY THE HIERARCHY CARRIES NO TIME AXIS. Coupling depth is not a divergence time in the
phylogenetic sense: co-occurrence measures shared lineage restriction, so two fates still
co-enriched within clades founded at E8.5 are concentrated in a restricted pool of E8.5 lineages
and are the MOST tightly related, not the least. Small depth therefore means tightly coupled, and
because average linkage merges small distances first, merge heights increase toward the root: a
scaled axis would place the root at the latest date and the deepest couplings at the tips. No
orientation fixes that -- putting early time on the left forces the root to the right, and putting
the root on the left forces time to run backwards. The tree is therefore drawn UNSCALED (uniform
rank spacing) with each resolved group annotated by its actual date. Detectability is also bounded
at early times: no pair is significant before E8.5, because clades founded then are large enough
that observed ~ expected under the null.

Average linkage over the coupling-depth matrix; the tree is built on the 74 cell types that couple
with at least one partner. The 8 that never couple are reported as an unresolved set.

Outputs c8_carpet.png, c8_dendrogram.png, c8_tree_merges.csv.
"""
import numpy as np, csv, json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram, cophenet, fcluster, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CK = _os.path.join(_REPO, "support_data")
Z_THR, MIN_OBS, CENSOR = 3.0, 5, 13.5
INK, INK2, MUTED = "#1f2328", "#4a5158", "#8b9198"
DIV = LinearSegmentedColormap.from_list("so", [
    (0.00, "#4dd2ff"), (0.25, "#1f6feb"), (0.50, "#000000"), (0.75, "#ff7a1a"), (1.00, "#ffd400")])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

d = np.load(f"{HERE}/out/c7_timesweep.npz", allow_pickle=True)
types = list(d["types"]); T = len(types)
TT = [float(t) for t in d["TIMES"] if f"B1|{float(t)}|enr" in d and f"B2|{float(t)}|enr" in d]
def sg(bl, t):
    o, x = d[f"{bl}|{t}|obs"], d[f"{bl}|{t}|exp"]
    return ((o - x) / np.sqrt(np.maximum(x, 1e-9)) >= Z_THR) & (o >= MIN_OBS)
SIG = {t: sg("B1", t) & sg("B2", t) for t in TT}
ENR = {t: (d[f"B1|{t}|enr"] + d[f"B2|{t}|enr"]) / 2 for t in TT}
iu = np.triu_indices(T, 1)

def depth_matrix(sigfn):
    D = np.full((T, T), np.nan)
    for t in TT:
        S = sigfn(t); D[S & np.isnan(D)] = t
    return D
Dboth = depth_matrix(lambda t: SIG[t])
DA = depth_matrix(lambda t: sg("B1", t)); DB = depth_matrix(lambda t: sg("B2", t))

gl = {}
for fn, okf in [(f"{CK}/germ_layer_map_validated.csv", lambda r: r["status"] in ("data-backed", "tree-resolved")),
                (f"{CK}/germ_layer_map_v6.csv", lambda r: r["confidence"] == "clear")]:
    for r in csv.DictReader(open(fn)):
        if r["celltype"] not in gl and r["germ_layer"] and okf(r): gl[r["celltype"]] = r["germ_layer"]

# Within a pair label, list the cell type that appears EARLIER in the Qiu et al. developmental
# graph first (depth from Oocyte; see c11_qiu_depth.py). Guidepost only -- it affects no statistic.
# Unmapped types and ties fall back to alphabetical, so the order is deterministic.
qd = {}
for r in csv.DictReader(open(f"{HERE}/out/qiu_celltype_depth.csv")):
    if r["qiu_depth_from_oocyte"] != "": qd[r["celltype"]] = int(r["qiu_depth_from_oocyte"])
def ordpair(x, y):
    """(earlier, later) by prior-graph depth, ties alphabetical."""
    kx = (qd.get(x, 10**6), x); ky = (qd.get(y, 10**6), y)
    return (x, y) if kx <= ky else (y, x)

# ---------------- carpet, every row labelled -------------------------------------------------
pairs = [(a, b) for a, b in zip(*iu) if not np.isnan(Dboth[a, b])]
curve = np.array([[ENR[t][a, b] for t in TT] for a, b in pairs])
sigm = np.array([[SIG[t][a, b] for t in TT] for a, b in pairs])
dv = np.array([Dboth[a, b] for a, b in pairs])
o1 = np.lexsort((-curve.max(1), dv))
lo, hi = np.percentile(curve, [2, 98])
norm = TwoSlopeNorm(vcenter=0.0, vmin=min(lo, -.1), vmax=max(hi, .1))
nrow = len(pairs)
fig, ax = plt.subplots(figsize=(13.0, max(11.0, 0.145 * nrow + 2.2)))
im = ax.imshow(curve[o1], cmap=DIV, norm=norm, aspect="auto", interpolation="nearest")
xt = [i for i, t in enumerate(TT) if abs(t - round(t)) < 1e-6]
ax.set_xticks(xt); ax.set_xticklabels([f"E{TT[i]:.0f}" for i in xt], fontsize=9)
ax.set_xlabel("developmental time of the clade ancestor  (earlier to the left)", labelpad=8)
lab = ["{}  |  {}".format(*ordpair(types[pairs[i][0]], types[pairs[i][1]])) for i in o1]
ax.set_yticks(range(nrow)); ax.set_yticklabels(lab, fontsize=5.6)
ax.tick_params(axis="y", length=0, pad=2)
ys, xs = np.where(sigm[o1])
ax.plot(xs, ys, ls="", marker=".", ms=1.6, color="white", alpha=.6)
ax.set_ylim(nrow - .5, -.5)
ax.set_title(f"Cell-type coupling across developmental time  (n = {nrow} pairs)\n"
             "clades defined by the age of their dated ancestor; rows ordered by coupling depth; "
             "white dots = significant in both blastomeres", fontsize=11, pad=12)
for sp in ax.spines.values(): sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=.30, pad=.015, label="log₂ enrichment (mean of A, B)")
cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{HERE}/out/c8_carpet.png", dpi=155, bbox_inches="tight")
print(f"carpet: {nrow} rows, all labelled")

# ---------------- dendrogram with dated branch points ----------------------------------------
keep = sorted({a for a, b in pairs} | {b for a, b in pairs})
kn = [types[i] for i in keep]
never = [t for t in types if t not in kn]
def tree(D):
    M = D[np.ix_(keep, keep)].copy()
    M[np.isnan(M)] = CENSOR
    np.fill_diagonal(M, 0.0); M = (M + M.T) / 2
    return linkage(squareform(M, checks=False), method="average"), M
Zb, Mb = tree(Dboth); Za, _ = tree(DA); Zc, _ = tree(DB)
cb_, _ = cophenet(Za, squareform(tree(DA)[1], checks=False))
coA = cophenet(Za); coB = cophenet(Zc)
rho = spearmanr(coA, coB).correlation; rp = pearsonr(coA, coB)[0]
print(f"\ncell types in the tree: {len(kn)}; never coupled with anything: {len(never)}")
for t in never: print(f"   unresolved: {t}")
print(f"\nA-vs-B cophenetic reproducibility of the dated tree: Spearman rho = {rho:.3f}, "
      f"Pearson r = {rp:.3f}")

# NO QUANTITATIVE TIME AXIS. Merge heights increase toward the root, so the root necessarily sits
# at the LATEST date and the deepest couplings at the tips -- the inverse of a phylogeny. That is a
# property of the metric (small coupling depth = tightly coupled, and average linkage merges small
# distances first), not of the drawing, so no orientation makes a scaled axis read correctly:
# putting early on the left forces the root to the right, and putting the root on the left forces
# time to run backwards. Instead the topology is drawn unscaled -- merge heights replaced by their
# rank, so spacing is uniform -- and each resolved group is annotated with its actual date. The
# root is then a conventional left-hand node and the dates are exact.
Zt = Zb.copy()
Zt[:, 2] = np.arange(1, len(Zt) + 1)                  # rank spacing; real dates go on as text
rank_to_age = {i + 1: float(Zb[i, 2]) for i in range(len(Zb))}
dn = dendrogram(Zt, labels=kn, orientation="left", no_plot=True)
XLO, XHI = len(Zt) + 1.0, -0.6                        # left = root, right = tips
fig, ax = plt.subplots(figsize=(11.2, 0.176 * len(kn) + 2.2))
for xs, ys in zip(dn["dcoord"], dn["icoord"]):
    ax.plot(xs, ys, color="#7f868d", lw=0.9, solid_capstyle="round", zorder=2)
# date every resolved group that closes at or before E11.5
for xs, ys in zip(dn["dcoord"], dn["icoord"]):
    r = int(round(xs[1]))
    age = rank_to_age.get(r)
    if age is None or age > 11.5: continue
    ax.annotate(f"E{age:.2f}", xy=(xs[1], (ys[1] + ys[2]) / 2), xytext=(-3, 0),
                textcoords="offset points", fontsize=5.6, color="#8a3324", ha="right",
                va="center", zorder=4)
lv = dn["ivl"]
ax.set_ylim(0, 10 * len(lv))
ax.set_yticks([5 + 10 * i for i in range(len(lv))]); ax.set_yticklabels(lv, fontsize=6.4)
ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")   # labels where the leaf stems end
ax.set_xlim(XLO, XHI)
ax.set_xticks([])
ax.set_xlabel("branch points labelled with the developmental time at which the group's coupling "
              "is still detectable (horizontal distance not to scale)", labelpad=10, fontsize=9)
ax.set_title("Cell types grouped by the depth at which their lineage coupling remains detectable\n"
             "dates are coupling depths in E-days, not divergence times; groups closing earliest "
             "stay co-restricted deepest", fontsize=10.5, pad=12)
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(axis="y", length=0, pad=3)
# group names in the left margin, aligned to the vertical centre of each group's leaves
GRP = [("Haematopoiesis", ["Microglia", "Primitive erythroid cells", "Hematopoietic stem cells (Cd34+)"]),
       ("Retina / eye field", ["Retinal progenitor cells", "Naive retinal progenitor cells"]),
       ("Gut tube", ["Gut", "Midgut/Hindgut epithelial cells"]),
       ("Kidney", ["Nephron progenitors", "Ureteric bud"]),
       ("Olfactory", ["Olfactory epithelial cells", "Olfactory sensory neurons"]),
       ("Otic", ["Otic sensory neurons", "Otic epithelial cells"]),
       ("Neural crest / PNS", ["Neural crest (PNS neurons)", "Melanocyte cells"]),
       ("Keratinocytes", ["Pre-epidermal keratinocytes", "Branchial arch epithelium"])]
for nm, anchors in GRP:
    ys = [5 + 10 * lv.index(a) for a in anchors if a in lv]
    if not ys: continue
    ax.annotate(nm, xy=(XLO + 0.06, float(np.mean(ys))), fontsize=6.8, color=INK,
                va="center", ha="left", annotation_clip=False, style="italic")
fig.tight_layout(); fig.savefig(f"{HERE}/out/c8_dendrogram.png", dpi=155, bbox_inches="tight")

# merge table
n = len(kn); lab2 = {i: kn[i] for i in range(n)}
with open(f"{HERE}/out/c8_tree_merges.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["merge", "age_E_day", "n_members", "members"])
    for i, (x, y, h, cnt) in enumerate(Zb):
        x, y = int(x), int(y)
        mem = (lab2[x].split(" + ") if x < n else lab2[x].split(" + ")) + \
              (lab2[y].split(" + ") if y < n else lab2[y].split(" + "))
        lab2[n + i] = " + ".join(mem)
        w.writerow([i + 1, round(float(h), 3), int(cnt), "; ".join(mem)])
print("\nfirst 14 merges (deepest couplings first):")
for i, (x, y, h, cnt) in enumerate(Zb[:14]):
    print(f"   E{h:.2f}  n={int(cnt)}  {lab2[n+i][:110]}")
print("\nwrote c8_carpet.png, c8_dendrogram.png, c8_tree_merges.csv")
